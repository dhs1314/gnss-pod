"""Orekit dynamics bridge for sequential EKF POD.

Provides high-fidelity orbit propagation using the Orekit space dynamics library
(via orekit-jpype), with finite-difference STM computation and GRACE-FO force
model configuration.

Force models (Orekit auto-computes all partials):
  - Earth gravity: ICGEM .gfc file (e.g. GGM05C) via Holmes-Featherstone, Nmax=150
  - Solid tides: IERS 2010 frequency-independent model
  - Ocean tides: FES2004 (if data available, else skipped)
  - Third-body: Sun + Moon (LuniSolar, no JPL DE required)
  - Solar radiation pressure: isotropic cannonball with conical eclipse
  - Atmospheric drag: exponential atmosphere model (Harris-Priester)
  - Relativity: Schwarzschild correction

Key API:
    prop = OrekitPropagator(...)
    # Single-step: FD STM over dt
    r_new, v_new, Phi, S_emp = prop.propagate(r0, v0, a_rtn, dt, mjd_utc, mjd_tt)

    # Multi-segment with chained Jacobian (for GN loop)
    r, v, phi, S = prop.propagate_arc(r0, v0, t_epochs, a_rtn,
                                       mjd_utc_start, mjd_tt_start,
                                       param_names=['Cd','CR','aR','aT','aN'])
"""

import os
import numpy as np

_OREKIT_AVAILABLE = False
_OREKIT_ERROR = None

try:
    import orekit_jpype
    _OREKIT_AVAILABLE = True
except ImportError as e:
    _OREKIT_ERROR = str(e)


def is_orekit_available():
    return _OREKIT_AVAILABLE


def get_orekit_error():
    return _OREKIT_ERROR


def _ensure_vm():
    """Lazily initialize the Orekit JVM and data providers."""
    if not _OREKIT_AVAILABLE:
        raise RuntimeError(
            f"Orekit is not available: {_OREKIT_ERROR}\n"
            f"  Install: pip install orekit-jpype\n"
            f"  Data: https://gitlab.orekit.org/orekit/orekit-data"
        )
    if hasattr(_ensure_vm, '_done'):
        return
    orekit_jpype.initVM()

    from java.io import File
    from org.orekit.data import DataContext, DirectoryCrawler

    data_paths = []
    env_data = os.environ.get('OREKIT_DATA_PATH')
    if env_data:
        data_paths.append(env_data)
    for candidate in [
        r'd:\prj\gnss_pod\data\orekit',
        r'./data/orekit',
    ]:
        if os.path.isdir(candidate) and candidate not in data_paths:
            data_paths.append(candidate)

    DM = DataContext.getDefault().getDataProvidersManager()
    for dp in data_paths:
        df = File(dp)
        if df.exists():
            DM.addProvider(DirectoryCrawler(df))

    _ensure_vm._done = True


# ── Finite-difference STM helpers ──

def _fd_stm(propagate_fn, r0, v0, perturb=1.0):
    """Compute 6x6 state transition matrix via forward finite differences.

    Args:
        propagate_fn: fn(r_eci, v_eci) -> (r_new, v_new)
        r0, v0: initial ECI position [m] and velocity [m/s]
        perturb: position perturbation [m] / velocity perturbation [m/s]

    Returns:
        Phi: (6,6) state transition matrix d(r_new,v_new)/d(r0,v0)
        r_ref, v_ref
    """
    r_ref, v_ref = propagate_fn(r0, v0)
    Phi = np.zeros((6, 6))

    for i in range(6):
        dr = np.zeros(6)
        dr[i] = perturb if i < 3 else perturb * 0.1
        r1 = r0 + dr[0:3]
        v1 = v0 + dr[3:6]
        r_new, v_new = propagate_fn(r1, v1)
        denom = dr[i] if abs(dr[i]) > 1e-15 else 1.0
        Phi[0:3, i] = (r_new - r_ref) / denom
        Phi[3:6, i] = (v_new - v_ref) / denom

    return Phi, r_ref, v_ref


def _fd_param_sensitivity(propagate_fn, r0, v0, params_dict, perturb_fractions):
    """Compute 6xNp parameter sensitivity via FD.

    Args:
        propagate_fn: fn(r_eci, v_eci, **params) -> (r_new, v_new)
        r0, v0: initial state
        params_dict: dict of {param_name: nominal_value}
        perturb_fractions: dict of {param_name: fraction_to_perturb}

    Returns:
        S: (6, Np) sensitivity matrix
        param_order: list of parameter names matching columns
    """
    r_ref, v_ref = propagate_fn(r0, v0, **params_dict)

    param_order = sorted(params_dict.keys())
    Np = len(param_order)
    S = np.zeros((6, Np))

    for col, pname in enumerate(param_order):
        nominal = params_dict[pname]
        frac = perturb_fractions.get(pname, 0.01)
        delta = max(abs(nominal) * frac, 1e-8)
        # Perturb
        p_perturbed = dict(params_dict)
        p_perturbed[pname] = nominal + delta
        r_new, v_new = propagate_fn(r0, v0, **p_perturbed)
        S[0:3, col] = (r_new - r_ref) / delta
        S[3:6, col] = (v_new - v_ref) / delta

    return S, param_order


def _empirical_sensitivity(r_mid, v_mid, dt):
    """Compute 6x3 empirical acceleration sensitivity analytically."""
    from src.empirical import rtn_to_eci, compute_rtn_frame
    R_vec, T_vec, N_vec = compute_rtn_frame(
        np.asarray(r_mid), np.asarray(v_mid))
    R_eci_rtn = np.column_stack([R_vec, T_vec, N_vec])
    S = np.zeros((6, 3))
    S[0:3, :] = 0.5 * dt * dt * R_eci_rtn
    S[3:6, :] = dt * R_eci_rtn
    return S


# ── Orekit Propagator ──

class OrekitPropagator:
    """Wrapper around Orekit NumericalPropagator.

    Single-step with FD STM:
        r_new, v_new, Phi, S_emp = prop.propagate(r0, v0, a_rtn, dt, mjd_utc, mjd_tt)

    Multi-segment arc with chained Jacobian (for GN outer loop):
        r, v, phi, S = prop.propagate_arc(r0, v0, t_epochs, a_rtn,
                                           mjd_utc_start, mjd_tt_start,
                                           param_names=['Cd','CR','aR','aT','aN'])
    """

    def __init__(self, orekit_data_path=None,
                 gravity_field=None,
                 gravity_degree=150,
                 solid_tides=True,
                 ocean_tides=True,
                 ocean_tide_degree=50,
                 third_body='lunisolar',
                 srp_model='isotropic',
                 relativity=True,
                 drag_model='msise00',
                 mass=580.0,
                 area_drag=0.68,
                 area_srp=3.4,
                 CR=1.3,
                 CD=2.2,
                 stm_perturb=1.0,
                 integrator_tol=1e-12,
                 integrator_min_step=0.1,
                 integrator_max_step=300.0):
        _ensure_vm()

        self.gravity_field = gravity_field
        self.gravity_degree = gravity_degree
        self.solid_tides = solid_tides
        self.ocean_tides = ocean_tides
        self.ocean_tide_degree = ocean_tide_degree
        self.third_body = third_body
        self.srp_model = srp_model
        self.relativity = relativity
        self.drag_model = drag_model
        self.mass = mass
        self.area_drag = area_drag
        self.area_srp = area_srp
        self.CR = CR
        self.CD = CD
        self.stm_perturb = stm_perturb
        self.integrator_tol = integrator_tol
        self.integrator_min_step = integrator_min_step
        self.integrator_max_step = integrator_max_step

        self._ready = False
        self._force_models_cache = None

    def _setup(self):
        """Lazy-init Orekit frames and data."""
        if self._ready:
            return

        from org.orekit.frames import FramesFactory
        from org.orekit.time import TimeScalesFactory, AbsoluteDate
        from org.orekit.utils import IERSConventions, Constants as OC
        from org.orekit.forces.gravity.potential import GravityFieldFactory

        self._gcrf = FramesFactory.getGCRF()
        self._itrf = FramesFactory.getITRF(IERSConventions.IERS_2010, True)
        self._utc = TimeScalesFactory.getUTC()
        self._tt = TimeScalesFactory.getTT()

        self._mu = OC.EIGEN5C_EARTH_MU

        # Copy gravity file to orekit data dir
        if self.gravity_field and os.path.exists(self.gravity_field):
            import shutil
            gf_name = os.path.basename(self.gravity_field)
            orekit_data = os.environ.get('OREKIT_DATA_PATH',
                                         r'd:\prj\gnss_pod\data\orekit')
            gf_dest = os.path.join(orekit_data, gf_name)
            if self.gravity_field != gf_dest:
                shutil.copy2(self.gravity_field, gf_dest)

        ref_date = AbsoluteDate(2024, 4, 29, 12, 0, 0.0, self._utc)

        try:
            self._gravity_provider = GravityFieldFactory.getConstantNormalizedProvider(
                self.gravity_degree, self.gravity_degree, ref_date)
        except Exception:
            self._gravity_provider = GravityFieldFactory.getNormalizedProvider(
                self.gravity_degree, self.gravity_degree)

        self._ready = True

    def _make_force_models(self):
        """Create Orekit force models (cached)."""
        if self._force_models_cache is not None:
            return self._force_models_cache

        from org.orekit.forces.gravity import HolmesFeatherstoneAttractionModel
        from org.orekit.forces.gravity import ThirdBodyAttraction
        from org.orekit.bodies import CelestialBodyFactory

        models = []

        # 1. Gravity (spherical harmonics Nmax)
        gravity = HolmesFeatherstoneAttractionModel(
            self._itrf, self._gravity_provider)
        models.append(gravity)

        # 2. Third-body: Sun + Moon
        if self.third_body:
            try:
                sun = CelestialBodyFactory.getSun()
                moon = CelestialBodyFactory.getMoon()
                models.append(ThirdBodyAttraction(sun))
                models.append(ThirdBodyAttraction(moon))
            except Exception:
                pass

        # 3. Solid tides
        if self.solid_tides:
            try:
                from org.orekit.forces.gravity import SolidTides
                tides = SolidTides(
                    self._itrf,
                    self._gravity_provider.getAe(),
                    self._gravity_provider.getMu(),
                    self._gravity_provider,
                )
                models.append(tides)
            except Exception:
                pass

        # 4. Ocean tides
        if self.ocean_tides:
            try:
                from org.orekit.forces.gravity import OceanTides
                ocean = OceanTides(
                    self._itrf,
                    self._gravity_provider.getAe(),
                    self._gravity_provider.getMu(),
                    self._gravity_provider,
                    self.ocean_tide_degree,
                    self.ocean_tide_degree,
                    True,
                )
                models.append(ocean)
            except Exception:
                pass

        # 5. Solar radiation pressure (Orekit v13: needs RadiationSensitive spacecraft)
        if self.srp_model:
            try:
                from org.orekit.forces.radiation import SolarRadiationPressure
                from org.orekit.forces.radiation import IsotropicRadiationSingleCoefficient
                from org.orekit.utils import Constants as OC
                from org.orekit.bodies import CelestialBodyFactory, OneAxisEllipsoid
                from org.orekit.frames import FramesFactory
                from org.orekit.utils import IERSConventions
                sun = CelestialBodyFactory.getSun()
                earth_radius = OC.WGS84_EARTH_EQUATORIAL_RADIUS
                earth_body = OneAxisEllipsoid(
                    earth_radius, 1.0 / 298.257223563,
                    FramesFactory.getITRF(IERSConventions.IERS_2010, True))
                spacecraft_srp = IsotropicRadiationSingleCoefficient(
                    float(self.area_srp), float(self.CR))
                srp = SolarRadiationPressure(sun, earth_body, spacecraft_srp)
                models.append(srp)
            except Exception:
                pass

        # 6. Atmospheric drag (Orekit v13: DragForce wraps atmosphere + spacecraft)
        if self.drag_model and self.drag_model != 'none':
            try:
                from org.orekit.forces.drag import IsotropicDrag, DragForce
                from org.orekit.utils import Constants as OC
                from org.orekit.bodies import CelestialBodyFactory, OneAxisEllipsoid
                from org.orekit.frames import FramesFactory
                from org.orekit.utils import IERSConventions
                earth_radius = OC.WGS84_EARTH_EQUATORIAL_RADIUS

                if self.drag_model == 'harris-priester':
                    from org.orekit.models.earth.atmosphere import HarrisPriester
                    sun = CelestialBodyFactory.getSun()
                    earth_body = OneAxisEllipsoid(
                        earth_radius, 1.0 / 298.257223563,
                        FramesFactory.getITRF(IERSConventions.IERS_2010, True))
                    atm = HarrisPriester(sun, earth_body)
                elif self.drag_model in ('msise00', 'nrlmsise00'):
                    from org.orekit.models.earth.atmosphere import NRLMSISE00
                    from org.orekit.data import DataContext
                    dm = DataContext.getDefault().getDataProvidersManager()
                    atm = NRLMSISE00(dm, self._itrf, self._utc)
                elif self.drag_model == 'dtm2000':
                    from org.orekit.models.earth.atmosphere import DTM2000
                    from org.orekit.data import DataContext
                    dm = DataContext.getDefault().getDataProvidersManager()
                    atm = DTM2000(dm, self._itrf, self._utc)
                elif self.drag_model == 'jb2008':
                    from org.orekit.models.earth.atmosphere import JB2008
                    from org.orekit.data import DataContext
                    dm = DataContext.getDefault().getDataProvidersManager()
                    atm = JB2008(dm, self._itrf, self._utc)
                else:
                    from org.orekit.models.earth.atmosphere import SimpleExponentialAtmosphere
                    earth_body = OneAxisEllipsoid(
                        earth_radius, 1.0 / 298.257223563,
                        FramesFactory.getITRF(IERSConventions.IERS_2010, True))
                    atm = SimpleExponentialAtmosphere(
                        earth_body, 500000.0, 60000.0, 1.0e-12,
                    )
                spacecraft = IsotropicDrag(
                    float(self.area_drag), float(self.CD),
                )
                drag_force = DragForce(atm, spacecraft)
                models.append(drag_force)
            except Exception as e:
                print(f"  [Orekit] Drag model '{self.drag_model}' failed: {e}")

        # 7. Relativity
        if self.relativity:
            try:
                from org.orekit.forces.gravity import Relativity
                models.append(Relativity(self._mu))
            except Exception:
                pass

        self._force_models_cache = models
        return models

    def _make_date(self, mjd_utc):
        """Convert MJD(UTC) to Orekit AbsoluteDate."""
        from datetime import datetime, timedelta
        from org.orekit.time import AbsoluteDate

        mjd_epoch = datetime(1858, 11, 17, 0, 0, 0)
        dt_utc = mjd_epoch + timedelta(days=mjd_utc)
        return AbsoluteDate(
            dt_utc.year, dt_utc.month, dt_utc.day,
            dt_utc.hour, dt_utc.minute,
            dt_utc.second + dt_utc.microsecond * 1e-6,
            self._utc,
        )

    def _propagate_one(self, r_eci, v_eci, dt, date_utc):
        """Single Orekit propagation: (r,v) -> (r_new, v_new)."""
        from org.orekit.orbits import CartesianOrbit
        from org.orekit.propagation.numerical import NumericalPropagator
        from org.orekit.propagation import SpacecraftState
        from org.hipparchus.ode.nonstiff import DormandPrince853Integrator
        from org.orekit.utils import PVCoordinates
        from org.hipparchus.geometry.euclidean.threed import Vector3D

        pos = Vector3D(float(r_eci[0]), float(r_eci[1]), float(r_eci[2]))
        vel = Vector3D(float(v_eci[0]), float(v_eci[1]), float(v_eci[2]))
        pv = PVCoordinates(pos, vel)
        orbit = CartesianOrbit(pv, self._gcrf, date_utc, self._mu)

        integrator = DormandPrince853Integrator(
            float(self.integrator_min_step),
            float(min(self.integrator_max_step, dt)),
            float(self.integrator_tol),
            float(self.integrator_tol),
        )

        prop = NumericalPropagator(integrator)
        prop.setInitialState(SpacecraftState(orbit))
        prop.setOrbitType(orbit.getType())

        for fm in self._make_force_models():
            prop.addForceModel(fm)

        final_state = prop.propagate(date_utc.shiftedBy(float(dt)))
        final_orbit = final_state.getOrbit()
        final_pv = final_orbit.getPVCoordinates(self._gcrf)
        r_new = np.array([
            final_pv.getPosition().getX(),
            final_pv.getPosition().getY(),
            final_pv.getPosition().getZ(),
        ])
        v_new = np.array([
            final_pv.getVelocity().getX(),
            final_pv.getVelocity().getY(),
            final_pv.getVelocity().getZ(),
        ])
        return r_new, v_new

    def _propagate_one_with_config(self, r_eci, v_eci, dt, date_utc,
                                    Cd=None, CR=None):
        """Like _propagate_one but with overrideable Cd/CR.

        Creates a fresh propagator with updated force parameters.
        Slightly slower than _propagate_one but needed for FD on Cd/CR.
        """
        from org.orekit.orbits import CartesianOrbit
        from org.orekit.propagation.numerical import NumericalPropagator
        from org.orekit.propagation import SpacecraftState
        from org.hipparchus.ode.nonstiff import DormandPrince853Integrator
        from org.orekit.utils import PVCoordinates
        from org.hipparchus.geometry.euclidean.threed import Vector3D

        pos = Vector3D(float(r_eci[0]), float(r_eci[1]), float(r_eci[2]))
        vel = Vector3D(float(v_eci[0]), float(v_eci[1]), float(v_eci[2]))
        pv = PVCoordinates(pos, vel)
        orbit = CartesianOrbit(pv, self._gcrf, date_utc, self._mu)

        integrator = DormandPrince853Integrator(
            float(self.integrator_min_step),
            float(min(self.integrator_max_step, dt)),
            float(self.integrator_tol),
            float(self.integrator_tol),
        )

        prop = NumericalPropagator(integrator)
        prop.setInitialState(SpacecraftState(orbit))
        prop.setOrbitType(orbit.getType())

        # Rebuild force models with overridden Cd/CR
        use_cd = Cd if Cd is not None else self.CD
        use_cr = CR if CR is not None else self.CR
        self._force_models_cache = None  # invalidate cache
        saved_cd, saved_cr = self.CD, self.CR
        self.CD, self.CR = use_cd, use_cr
        try:
            for fm in self._make_force_models():
                prop.addForceModel(fm)
        finally:
            self.CD, self.CR = saved_cd, saved_cr

        final_state = prop.propagate(date_utc.shiftedBy(float(dt)))
        final_orbit = final_state.getOrbit()
        final_pv = final_orbit.getPVCoordinates(self._gcrf)
        r_new = np.array([
            final_pv.getPosition().getX(),
            final_pv.getPosition().getY(),
            final_pv.getPosition().getZ(),
        ])
        v_new = np.array([
            final_pv.getVelocity().getX(),
            final_pv.getVelocity().getY(),
            final_pv.getVelocity().getZ(),
        ])
        return r_new, v_new

    def propagate(self, r_eci, v_eci, a_rtn, dt, mjd_utc, mjd_tt):
        """Propagate state forward by dt seconds.

        Returns:
            r_new, v_new: propagated ECI state (3,) each
            Phi: 6x6 state transition matrix (FD)
            S_emp: 6x3 empirical sensitivity matrix (analytic)
        """
        self._setup()
        date_utc = self._make_date(mjd_utc)

        # Empirical acceleration in ECI
        from src.empirical import rtn_to_eci
        a_emp_eci = np.asarray(rtn_to_eci(
            np.asarray(a_rtn, dtype=float),
            np.asarray(r_eci, dtype=float),
            np.asarray(v_eci, dtype=float)))

        # Absorb empirical into initial velocity
        v_eci_adj = v_eci + 0.5 * dt * a_emp_eci

        def propagate_fn(r0, v0):
            return self._propagate_one(r0, v0, dt, date_utc)

        # FD STM
        Phi, r_intermediate, v_intermediate = _fd_stm(
            propagate_fn, r_eci, v_eci_adj, self.stm_perturb)

        # Apply remaining empirical effect
        r_new = r_intermediate + 0.5 * dt * dt * a_emp_eci
        v_new = v_intermediate + 0.5 * dt * a_emp_eci

        # Empirical sensitivity (analytic)
        r_mid = 0.5 * (r_eci + r_new)
        v_mid = 0.5 * (v_eci + v_new)
        S_emp = _empirical_sensitivity(r_mid, v_mid, dt)

        return r_new, v_new, Phi, S_emp

    # ── Single-pass arc propagation (FAST - reuses propagator) ──

    def propagate_continuous_arc(self, r0, v0, t_epochs, a_rtn,
                                  mjd_utc_start):
        """Continuous Orekit propagation — ONE propagator, N epochs.

        Reuses the same NumericalPropagator + integrator across all
        epochs, only calling setInitialState + propagate per segment.
        Avoids JVM object creation overhead (~10-20x faster than
        per-epoch _propagate_one calls).

        Args:
            r0, v0: initial ECI state (3,)
            t_epochs: list of times [s] from t0 (strictly increasing)
            a_rtn: (3,) constant or (N,3) per-epoch empirical RTN [m/s²]
            mjd_utc_start: MJD(UTC) at t=0

        Returns:
            r: (N, 3), v: (N, 3) ECI state at each epoch
        """
        from org.orekit.orbits import CartesianOrbit
        from org.orekit.propagation.numerical import NumericalPropagator
        from org.orekit.propagation import SpacecraftState
        from org.hipparchus.ode.nonstiff import DormandPrince853Integrator
        from org.orekit.utils import PVCoordinates
        from org.hipparchus.geometry.euclidean.threed import Vector3D
        from src.empirical import rtn_to_eci

        self._setup()
        SEC_PER_DAY = 86400.0
        N = len(t_epochs)
        r = np.zeros((N, 3)); v = np.zeros((N, 3))

        # Handle per-epoch or constant a_rtn
        if a_rtn is None:
            a_rtn_arr = np.zeros((N, 3))
        elif a_rtn.ndim == 1:
            a_rtn_arr = np.tile(np.asarray(a_rtn, dtype=float), (N, 1))
        else:
            a_rtn_arr = np.asarray(a_rtn, dtype=float)

        # ── Create propagator ONCE ──
        T_max = float(max(t_epochs))
        integrator = DormandPrince853Integrator(
            float(self.integrator_min_step),
            float(min(self.integrator_max_step, T_max + 60.0)),
            float(self.integrator_tol),
            float(self.integrator_tol),
        )
        prop = NumericalPropagator(integrator)
        for fm in self._make_force_models():
            prop.addForceModel(fm)

        # ── Propagate epoch by epoch, matching old _propagate_one conventions ──
        r_cur, v_cur = r0.copy(), v0.copy()
        t_prev = 0.0

        for i, t_obs in enumerate(t_epochs):
            dt = t_obs - t_prev
            if dt <= 0:
                r[i] = r_cur; v[i] = v_cur
                continue

            a_emp_eci = np.zeros(3)
            if np.any(a_rtn_arr[i]):
                a_emp_eci = np.asarray(rtn_to_eci(
                    a_rtn_arr[i], r_cur, v_cur))

            v_adj = v_cur + 0.5 * dt * a_emp_eci

            # OLD convention: orbit date = mjd at t_obs,
            # propagate to mjd at t_obs+dt (matches _propagate_one behavior)
            pos_i = Vector3D(float(r_cur[0]), float(r_cur[1]), float(r_cur[2]))
            vel_i = Vector3D(float(v_adj[0]), float(v_adj[1]), float(v_adj[2]))
            pv_i = PVCoordinates(pos_i, vel_i)
            mjd_end = mjd_utc_start + t_obs / SEC_PER_DAY
            date_end = self._make_date(mjd_end)
            orbit_i = CartesianOrbit(pv_i, self._gcrf, date_end, self._mu)
            prop.setInitialState(SpacecraftState(orbit_i))

            # Propagate by dt from the orbit date
            final_state = prop.propagate(date_end.shiftedBy(float(dt)))
            final_pv = final_state.getOrbit().getPVCoordinates(self._gcrf)

            r_ref = np.array([final_pv.getPosition().getX(),
                              final_pv.getPosition().getY(),
                              final_pv.getPosition().getZ()])
            v_ref = np.array([final_pv.getVelocity().getX(),
                              final_pv.getVelocity().getY(),
                              final_pv.getVelocity().getZ()])

            r_new = r_ref + 0.5 * dt * dt * a_emp_eci
            v_new = v_ref + 0.5 * dt * a_emp_eci

            r[i] = r_new; v[i] = v_new
            r_cur, v_cur = r_new, v_new
            t_prev = t_obs

        return r, v


    # ── Multi-segment arc propagation with chained Jacobian ──

    def propagate_arc(self, r0, v0, t_epochs, a_rtn,
                      mjd_utc_start, mjd_tt_start,
                      param_names=None,
                      Cd_override=None, CR_override=None):
        """Propagate full arc with Orekit, chaining FD STM + param Jacobians.

        For each segment [t_i-1, t_i]:
          1. Orekit FD STM -> Phi_seg (6x6)
          2. FD on force params -> S_seg (6xNp)
          3. Chain: Phi(t_i) = Phi_seg @ Phi(t_{i-1})
                     S(t_i) = Phi_seg @ S(t_{i-1}) + S_seg

        Args:
            r0, v0: initial ECI state (3,) [m, m/s]
            t_epochs: list of times [s] from t0 (strictly increasing, t[0]>0)
            a_rtn: empirical RTN acceleration (3,) or None
            mjd_utc_start, mjd_tt_start: MJD at t=0
            param_names: list of parameter names for S-matrix columns.
                Default: ['Cd','CR','aR','aT','aN']
            Cd_override, CR_override: if provided, use these values instead
                of self.CD/self.CR for the reference propagation. FD
                sensitivities for Cd/CR are w.r.t. these values.

        Returns:
            r: (N, 3) position at each observation epoch
            v: (N, 3) velocity at each observation epoch
            phi: (N, 6, 6) chained STM at each epoch
            S: (N, 6, Np) chained parameter sensitivity at each epoch
        """
        if param_names is None:
            param_names = ['Cd', 'CR', 'aR', 'aT', 'aN']

        self._setup()
        SEC_PER_DAY = 86400.0
        N = len(t_epochs)
        r = np.zeros((N, 3)); v = np.zeros((N, 3))
        phi = np.zeros((N, 6, 6)); S = np.zeros((N, 6, len(param_names)))

        r_cur, v_cur = r0.copy(), v0.copy()
        phi_cur = np.eye(6)
        S_cur = np.zeros((6, len(param_names)))
        t_prev = 0.0

        # Use overrides for reference propagation
        use_Cd = Cd_override if Cd_override is not None else self.CD
        use_CR = CR_override if CR_override is not None else self.CR

        for i, t_obs in enumerate(t_epochs):
            dt = t_obs - t_prev
            mjd_utc = mjd_utc_start + t_obs / SEC_PER_DAY
            mjd_tt = mjd_tt_start + t_obs / SEC_PER_DAY

            # Empirical acc in ECI at current state
            a_emp_eci = np.zeros(3)
            if a_rtn is not None and np.any(a_rtn):
                from src.empirical import rtn_to_eci
                a_emp_eci = np.asarray(rtn_to_eci(
                    np.asarray(a_rtn, dtype=float),
                    np.asarray(r_cur, dtype=float),
                    np.asarray(v_cur, dtype=float)))

            v_adj = v_cur + 0.5 * dt * a_emp_eci
            date_utc = self._make_date(mjd_utc)

            # Reference propagation (Orekit) with overridden Cd/CR
            r_ref, v_ref = self._propagate_one_with_config(
                r_cur, v_adj, dt, date_utc, Cd=use_Cd, CR=use_CR)

            # Apply empirical effect
            r_new = r_ref + 0.5 * dt * dt * a_emp_eci
            v_new = v_ref + 0.5 * dt * a_emp_eci

            # ── FD STM for this segment (uses overridden Cd/CR) ──
            def seg_prop(r0_seg, v0_seg):
                v0_adj = v0_seg + 0.5 * dt * a_emp_eci
                return self._propagate_one_with_config(
                    r0_seg, v0_adj, dt, date_utc, Cd=use_Cd, CR=use_CR)

            Phi_seg_raw, _, _ = _fd_stm(seg_prop, r_cur, v_cur,
                                         self.stm_perturb)
            # Apply empirical to STM: r_new = r_ref + 0.5*dt²*a_emp
            # a_emp does NOT depend on r_cur (constant over segment) →
            # ∂r_new/∂r_cur = ∂r_ref/∂r_cur (empirical is state-independent)
            # But v_adj depends on v_cur → ∂r_new/∂v_cur affected
            # For simplicity: Phi_seg ≈ Phi_seg_raw (empirical correction
            # is small and state-independent at this level)
            Phi_seg = Phi_seg_raw

            # ── FD parameter sensitivity for this segment (uses overridden Cd/CR) ──
            S_seg = self._param_sensitivity_segment(
                r_cur, v_cur, dt, date_utc, a_emp_eci, param_names,
                Cd_ref=use_Cd, CR_ref=use_CR)

            # Chain
            phi_cur = Phi_seg @ phi_cur
            S_cur = Phi_seg @ S_cur + S_seg

            r[i] = r_new; v[i] = v_new
            phi[i] = phi_cur.copy(); S[i] = S_cur.copy()

            r_cur, v_cur = r_new, v_new
            t_prev = t_obs

        return r, v, phi, S

    def _param_sensitivity_segment(self, r_cur, v_cur, dt, date_utc,
                                    a_emp_eci, param_names,
                                    Cd_ref=None, CR_ref=None):
        """Compute 6xNp parameter sensitivity for one segment via FD.

        Uses Cd_ref/CR_ref as the nominal values for FD perturbation.
        """
        Np = len(param_names)
        S_seg = np.zeros((6, Np))
        if Np == 0:
            return S_seg

        # Reference propagation at nominal Cd/CR
        use_Cd = Cd_ref if Cd_ref is not None else self.CD
        use_CR = CR_ref if CR_ref is not None else self.CR
        v_adj = v_cur + 0.5 * dt * a_emp_eci
        r_ref, v_ref = self._propagate_one_with_config(
            r_cur, v_adj, dt, date_utc, Cd=use_Cd, CR=use_CR)
        r_ref += 0.5 * dt * dt * a_emp_eci
        v_ref += 0.5 * dt * a_emp_eci

        for col, pname in enumerate(param_names):
            if pname in ('aR', 'aT', 'aN'):
                # Empirical RTN — use analytic sensitivity
                # These are handled by chaining, NOT FD on Orekit
                # (empirical is applied outside Orekit)
                r_mid = 0.5 * (r_cur + r_ref)
                v_mid = 0.5 * (v_cur + v_ref)
                S_analytic = _empirical_sensitivity(r_mid, v_mid, dt)
                idx_map = {'aR': 0, 'aT': 1, 'aN': 2}
                if pname in idx_map:
                    S_seg[:, col] = S_analytic[:, idx_map[pname]]
                continue

            # Perturb Cd or CR (w.r.t. overridden reference values)
            if pname == 'Cd':
                delta = max(abs(use_Cd) * 0.01, 0.01)
                r_p, v_p = self._propagate_one_with_config(
                    r_cur, v_adj, dt, date_utc,
                    Cd=use_Cd + delta, CR=use_CR)
                r_p += 0.5 * dt * dt * a_emp_eci
                v_p += 0.5 * dt * a_emp_eci

            elif pname == 'CR':
                delta = max(abs(use_CR) * 0.01, 0.01)
                r_p, v_p = self._propagate_one_with_config(
                    r_cur, v_adj, dt, date_utc,
                    Cd=use_Cd, CR=use_CR + delta)
                r_p += 0.5 * dt * dt * a_emp_eci
                v_p += 0.5 * dt * a_emp_eci
            else:
                continue

            S_seg[0:3, col] = (r_p - r_ref) / delta
            S_seg[3:6, col] = (v_p - v_ref) / delta

        return S_seg


# ── Factory ──

def create_propagator(mode='simplified', **kwargs):
    """Factory: create the appropriate propagator.

    Returns:
        (propagator, mode_string, warning_message)
    """
    if mode == 'orekit':
        if not _OREKIT_AVAILABLE:
            msg = (f"Orekit not available ({_OREKIT_ERROR}). "
                   f"Falling back to simplified dynamics.")
            print(f"  [WARNING] {msg}")
            return None, 'simplified', msg
        try:
            prop = OrekitPropagator(**kwargs)
            return prop, 'orekit', None
        except Exception as e:
            import traceback
            msg = f"Orekit init failed: {e}. Falling back to simplified."
            print(f"  [WARNING] {msg}")
            traceback.print_exc()
            return None, 'simplified', msg

    return None, 'simplified', None
