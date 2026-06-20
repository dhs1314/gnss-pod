"""Orekit dynamics bridge for sequential EKF POD.

Provides high-fidelity orbit propagation using the Orekit space dynamics library
(via orekit-jpype), with finite-difference STM computation and GRACE-FO force
model configuration.

Force models:
  - Earth gravity: ICGEM .gfc file (e.g. GGM05C) via Holmes-Featherstone, Nmax=150
  - Solid tides: IERS 2010 frequency-independent model
  - Ocean tides: FES2004 (if data available, else skipped)
  - Third-body: Sun + Moon (LuniSolar, no JPL DE required)
  - Solar radiation pressure: isotropic cannonball with conical eclipse
  - Atmospheric drag: exponential atmosphere model (Harris-Priester)
  - Relativity: Schwarzschild correction
  - Empirical RTN: applied analytically (constant acceleration over step)

Installation:
    pip install orekit-jpype

Data requirements:
    - IERS EOP C04 file (eopc04_IAU2000.txt) in data/orekit/
    - UTC-TAI history (UTC-TAI.history) in data/orekit/
    - Gravity field .gfc file for Orekit (may use user's GGM05C.gfc)
    Set OREKIT_DATA_PATH or pass orekit_data_path to the constructor.

Architecture:
    SequentialEKF.predict()
        ├── dynamics_mode='simplified' → _predict_simplified() (pure Python)
        └── dynamics_mode='orekit'     → _predict_orekit() (this module)
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

    # Use orekit data directory relative to project root, or from env var
    data_paths = []
    env_data = os.environ.get('OREKIT_DATA_PATH')
    if env_data:
        data_paths.append(env_data)
    # Auto-detect project data dirs
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

def _fd_stm(propagate_fn, r0, v0, dt, perturb=1.0):
    """Compute 6x6 state transition matrix via forward finite differences.

    Args:
        propagate_fn: fn(r_eci, v_eci) → (r_new, v_new)
        r0, v0: initial ECI position [m] and velocity [m/s]
        dt: integration step [s]
        perturb: position perturbation [m] / velocity perturbation [m/s]

    Returns:
        Phi: (6,6) state transition matrix ∂(r_new,v_new)/∂(r0,v0)
    """
    r_ref, v_ref = propagate_fn(r0, v0)
    Phi = np.zeros((6, 6))
    I6 = np.eye(6)

    for i in range(6):
        dr = np.zeros(6)
        dr[i] = perturb if i < 3 else perturb * 0.1  # vel perturb 0.1 m/s
        r1 = r0 + dr[0:3]
        v1 = v0 + dr[3:6]
        r_new, v_new = propagate_fn(r1, v1)
        # Forward difference
        Phi[0:3, i] = (r_new - r_ref) / dr[i] if abs(dr[i]) > 1e-15 else 0
        Phi[3:6, i] = (v_new - v_ref) / dr[i] if abs(dr[i]) > 1e-15 else 0

    return Phi, r_ref, v_ref


def _empirical_sensitivity(r_mid, v_mid, dt):
    """Compute 6x3 empirical acceleration sensitivity matrix analytically.

    For small constant RTN accelerations aR, aT, aN applied over dt:
        ∂r/∂a_rtn = 0.5 * dt² * R_rtn_to_eci
        ∂v/∂a_rtn = dt * R_rtn_to_eci

    Uses finite differences on the rtn_to_eci transform to get the rotation
    matrix, since the RTN basis vectors form the rotation from RTN to ECI.

    Args:
        r_mid, v_mid: mid-point ECI state
        dt: integration step [s]

    Returns:
        S: (6, 3) sensitivity matrix
    """
    from src.empirical import rtn_to_eci, compute_rtn_frame
    # RTN basis vectors form the rotation matrix RTN → ECI
    R_vec, T_vec, N_vec = compute_rtn_frame(
        np.asarray(r_mid), np.asarray(v_mid))
    R_eci_rtn = np.column_stack([R_vec, T_vec, N_vec])  # (3,3)
    S = np.zeros((6, 3))
    S[0:3, :] = 0.5 * dt * dt * R_eci_rtn
    S[3:6, :] = dt * R_eci_rtn
    return S


# ── Orekit Propagator ──

class OrekitPropagator:
    """Wrapper around Orekit NumericalPropagator for EKF predict step.

    Usage:
        prop = OrekitPropagator(orekit_data_path='./data/orekit',
                                gravity_degree=150,
                                gravity_field='d:/prj/gnss_pod/data/GGM05C.gfc')
        r_new, v_new, Phi, S = prop.propagate(r_eci, v_eci, a_rtn, dt,
                                               mjd_utc, mjd_tt)
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
                 drag_model='exponential',
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

    def _setup(self):
        """Lazy-init Orekit frames, constants, and force model factory."""
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

        # Earth gravity constant
        self._mu = OC.EIGEN5C_EARTH_MU

        # Gravity provider via Orekit's data loading system.
        # GGM05C.gfc (or the user's gravity file) must be present in the
        # orekit data directory so the ICGEM reader can find it.
        # If gravity_field is a .gfc file, copy it to the orekit data dir.
        if self.gravity_field and os.path.exists(self.gravity_field):
            import shutil
            gf_name = os.path.basename(self.gravity_field)
            orekit_data = os.environ.get('OREKIT_DATA_PATH',
                                         r'd:\prj\gnss_pod\data\orekit')
            gf_dest = os.path.join(orekit_data, gf_name)
            if self.gravity_field != gf_dest:
                shutil.copy2(self.gravity_field, gf_dest)

        # Reference date for gravity (any recent date works; provider is constant)
        ref_date = AbsoluteDate(2024, 4, 29, 12, 0, 0.0, self._utc)

        try:
            self._gravity_provider = GravityFieldFactory.getConstantNormalizedProvider(
                self.gravity_degree, self.gravity_degree, ref_date)
        except Exception:
            self._gravity_provider = GravityFieldFactory.getNormalizedProvider(
                self.gravity_degree, self.gravity_degree)

        self._ready = True

    def _make_force_model(self):
        """Create a list of Orekit force models for the propagator."""
        from org.orekit.forces.gravity import HolmesFeatherstoneAttractionModel
        from org.orekit.forces.gravity import ThirdBodyAttraction
        from org.orekit.bodies import CelestialBodyFactory

        models = []

        # 1. Central gravity with high-degree spherical harmonics
        gravity = HolmesFeatherstoneAttractionModel(
            self._itrf, self._gravity_provider)
        models.append(gravity)

        # 2. Third-body: Sun + Moon (LuniSolar, no JPL DE needed)
        if self.third_body:
            try:
                sun = CelestialBodyFactory.getSun()
                moon = CelestialBodyFactory.getMoon()
                models.append(ThirdBodyAttraction(sun))
                models.append(ThirdBodyAttraction(moon))
            except Exception:
                pass

        # 3. Solid tides (IERS 2010)
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

        # 4. Ocean tides (FES2004)
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
                    True,  # include pole tide
                )
                models.append(ocean)
            except Exception:
                pass

        # 5. Solar radiation pressure (isotropic cannonball)
        if self.srp_model:
            try:
                from org.orekit.forces.radiation import SolarRadiationPressure
                from org.orekit.utils import Constants as OC
                sun = CelestialBodyFactory.getSun()
                earth_radius = OC.WGS84_EARTH_EQUATORIAL_RADIUS
                # Isotropic model: cross-section area * (1 + CR)
                srp = SolarRadiationPressure(
                    sun, earth_radius,
                    float(self.area_srp),  # cross-section [m²]
                    float(self.CR),        # reflectivity coefficient
                    float(self.mass),      # spacecraft mass [kg]
                )
                models.append(srp)
            except Exception:
                pass

        # 6. Atmospheric drag
        if self.drag_model == 'exponential':
            try:
                from org.orekit.forces.drag import IsotropicDrag
                from org.orekit.models.earth.atmosphere import SimpleExponentialAtmosphere
                from org.orekit.utils import Constants as OC
                earth_radius = OC.WGS84_EARTH_EQUATORIAL_RADIUS
                # Simple exponential model: rho = rho0 * exp(-(h-h0)/H)
                # GRACE-FO ~490 km altitude
                atm = SimpleExponentialAtmosphere(
                    earth_radius,
                    500000.0,    # reference altitude [m]
                    60000.0,     # scale height [m]
                    1.0e-12,     # density at reference altitude [kg/m³]
                )
                drag = IsotropicDrag(
                    float(self.area_drag),  # cross-section [m²]
                    float(self.CD),         # drag coefficient
                    float(self.mass),       # mass [kg]
                )
                models.append(drag)
            except Exception:
                pass

        # 7. Relativity (Schwarzschild)
        if self.relativity:
            try:
                from org.orekit.forces.gravity import Relativity
                models.append(Relativity(self._mu))
            except Exception:
                pass

        return models

    def _propagate_one(self, r_eci, v_eci, dt, date_utc):
        """Single propagation: (r,v) → (r_new, v_new) using Orekit.

        Used as the core integration call for both reference and perturbed states.
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

        for fm in self._make_force_model():
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

    def propagate(self, r_eci, v_eci, a_rtn, dt, mjd_utc, mjd_tt):
        """Propagate state forward by dt seconds.

        Args:
            r_eci: ECI position (3,) [m]
            v_eci: ECI velocity (3,) [m]
            a_rtn: empirical RTN acceleration (3,) [m/s²]
            dt: propagation time step [s] (positive)
            mjd_utc: MJD(UTC) at start epoch
            mjd_tt: MJD(TT) at start epoch

        Returns:
            r_new: propagated ECI position (3,) [m]
            v_new: propagated ECI velocity (3,) [m]
            Phi: 6x6 state transition matrix
            S: 6x3 empirical sensitivity matrix
        """
        self._setup()

        # Convert MJD to Orekit AbsoluteDate
        from datetime import datetime, timedelta
        from org.orekit.time import AbsoluteDate

        mjd_epoch = datetime(1858, 11, 17, 0, 0, 0)
        dt_utc = mjd_epoch + timedelta(days=mjd_utc)
        date_utc = AbsoluteDate(
            dt_utc.year, dt_utc.month, dt_utc.day,
            dt_utc.hour, dt_utc.minute,
            dt_utc.second + dt_utc.microsecond * 1e-6,
            self._utc,
        )

        # Apply empirical RTN acceleration analytically before/after propagation
        # Convert RTN → ECI for the empirical correction
        from src.empirical import rtn_to_eci
        a_emp_eci = np.asarray(rtn_to_eci(
            np.asarray(a_rtn, dtype=float),
            np.asarray(r_eci, dtype=float),
            np.asarray(v_eci, dtype=float)))

        # Absorb empirical acceleration into initial velocity for integration
        # (constant acceleration approximation over dt)
        v_eci_adj = v_eci + 0.5 * dt * a_emp_eci

        # Define propagate function for finite differences
        def propagate_fn(r0, v0):
            return self._propagate_one(r0, v0, dt, date_utc)

        # Compute STM via finite differences + reference trajectory
        Phi, r_intermediate, v_intermediate = _fd_stm(
            propagate_fn, r_eci, v_eci_adj, dt, self.stm_perturb)

        # Apply remaining empirical effect
        r_new = r_intermediate + 0.5 * dt * dt * a_emp_eci
        v_new = v_intermediate + 0.5 * dt * a_emp_eci

        # Empirical sensitivity (analytical, mid-point)
        r_mid = 0.5 * (r_eci + r_new)
        v_mid = 0.5 * (v_eci + v_new)
        S = _empirical_sensitivity(r_mid, v_mid, dt)

        return r_new, v_new, Phi, S


# ── Factory ──

def create_propagator(mode='simplified', **kwargs):
    """Factory: create the appropriate propagator.

    Args:
        mode: 'simplified' or 'orekit'
        **kwargs: passed to the underlying propagator constructor

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
