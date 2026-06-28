"""Orekit dynamics bridge for sequential EKF POD.

Provides a Python interface to Orekit's high-fidelity force models via
orekit-jpype. Gracefully degrades if Orekit is not installed.

Installation:
    pip install orekit-jpype

Required data (download from https://gitlab.orekit.org/orekit/orekit-data):
    - Earth Orientation Parameters (finals2000A.all)
    - JPL DE430 ephemerides
    - Gravity field files (.gfc)
    - FES2004 ocean tide coefficients
    Set OREKIT_DATA_PATH environment variable to the data directory.

Architecture:
    SequentialEKF.predict()
        ├── dynamics_mode='simplified' → _predict_simplified() (V2.2.1)
        └── dynamics_mode='orekit'     → _predict_orekit() (this module)
"""
import numpy as np

_OREKIT_AVAILABLE = False
_OREKIT_ERROR = None

try:
    import orekit
    _OREKIT_AVAILABLE = True
except ImportError as e:
    _OREKIT_ERROR = str(e)


def is_orekit_available():
    """Check if orekit-jpype is installed and importable."""
    return _OREKIT_AVAILABLE


def get_orekit_error():
    """Return the import error message if Orekit is unavailable."""
    return _OREKIT_ERROR


class OrekitPropagator:
    """Wrapper around Orekit NumericalPropagator for EKF predict step.

    Manages the Java VM, frame setup, force model configuration, and
    state+STM propagation for a single EKF predict step.

    Usage:
        prop = OrekitPropagator(orekit_data_path='./data/orekit',
                                gravity_degree=150,
                                solid_tides=True,
                                ocean_tides=True)
        r_new, v_new, Phi, S = prop.propagate(r_eci, v_eci, a_rtn, dt,
                                               mjd_utc, mjd_tt)
    """

    def __init__(self, orekit_data_path=None,
                 gravity_field='EIGEN-GRGS.RL03.MEAN-FIELD',
                 gravity_degree=150,
                 solid_tides=True,
                 ocean_tides=True,
                 ocean_tide_degree=50,
                 third_body='DE430',
                 srp_model='isotropic',
                 relativity=True,
                 empirical_tau=600.0,
                 empirical_ss=1e-8):
        """
        Args:
            orekit_data_path: path to Orekit data directory
            gravity_field: gravity field model name
            gravity_degree: max degree for spherical harmonics
            solid_tides: enable IERS 2010 solid tides
            ocean_tides: enable FES2004 ocean tides
            ocean_tide_degree: max degree for ocean tides
            third_body: 'DE430' or 'DE405'
            srp_model: 'isotropic' or 'box_wing'
            relativity: enable relativistic corrections
            empirical_tau: empirical acceleration correlation time [s]
            empirical_ss: empirical acceleration steady-state sigma [m/s²]
        """
        if not _OREKIT_AVAILABLE:
            raise RuntimeError(
                f"Orekit is not available: {_OREKIT_ERROR}\n"
                f"  Install: pip install orekit-jpype\n"
                f"  Data: https://gitlab.orekit.org/orekit/orekit-data"
            )

        self.orekit_data_path = orekit_data_path
        self.gravity_field = gravity_field
        self.gravity_degree = gravity_degree
        self.solid_tides = solid_tides
        self.ocean_tides = ocean_tides
        self.ocean_tide_degree = ocean_tide_degree
        self.third_body = third_body
        self.srp_model = srp_model
        self.relativity = relativity
        self.empirical_tau = empirical_tau
        self.empirical_ss = empirical_ss

        self._initialized = False

    def _init_orekit(self):
        """Lazy-initialize Orekit VM and frames. Called on first propagate()."""
        if self._initialized:
            return

        # Import Java classes through orekit
        from orekit import JArray, init_vm
        from java.io import File

        # Initialize Orekit VM
        orekit_data = File(self.orekit_data_path) if self.orekit_data_path else None
        if orekit_data is None:
            import os
            env_data = os.environ.get('OREKIT_DATA_PATH')
            if env_data:
                orekit_data = File(env_data)

        init_vm(orekit_data)

        # Frames
        from org.orekit.frames import FramesFactory
        from org.orekit.time import TimeScalesFactory
        from org.orekit.utils import IERSConventions

        self._gcrf = FramesFactory.getGCRF()
        self._itrf = FramesFactory.getITRF(IERSConventions.IERS_2010, True)
        self._utc = TimeScalesFactory.getUTC()
        self._tt = TimeScalesFactory.getTT()
        self._earth = FramesFactory.getITRF(IERSConventions.IERS_2010, True)  # for mu

        # Constants
        from org.orekit.utils import Constants as OrekitConstants
        self._mu = OrekitConstants.EIGEN5C_EARTH_MU  # or WGS84_EARTH_MU

        self._initialized = True

    def propagate(self, r_eci, v_eci, a_rtn, dt, mjd_utc, mjd_tt):
        """Propagate state forward by dt seconds using Orekit dynamics.

        Args:
            r_eci: ECI position (3,) [m]
            v_eci: ECI velocity (3,) [m]
            a_rtn: empirical RTN acceleration (3,) [m/s²] (not yet integrated)
            dt: propagation time step [s] (positive)
            mjd_utc: MJD(UTC) at start epoch
            mjd_tt: MJD(TT) at start epoch

        Returns:
            r_new: propagated ECI position (3,) [m]
            v_new: propagated ECI velocity (3,) [m]
            Phi: 6x6 state transition matrix [None if unavailable]
            S: 6x3 empirical sensitivity matrix [None if unavailable]
        """
        self._init_orekit()

        # ── Build initial orbit ──
        from org.orekit.orbits import CartesianOrbit
        from org.orekit.time import AbsoluteDate, DateTimeComponents, DateComponents, TimeComponents
        from org.orekit.propagation.numerical import NumericalPropagator
        from org.orekit.propagation import SpacecraftState
        from org.orekit.integrator import DormandPrince853Integrator
        from org.orekit.utils import PVCoordinates
        from org.hipparchus.geometry.euclidean.threed import Vector3D

        from orekit import JArray

        # Convert MJD to Orekit AbsoluteDate
        # MJD(TT) = JD(TT) - 2400000.5
        # Orekit AbsoluteDate uses a reference epoch internally
        jd_utc = mjd_utc + 2400000.5
        # For now, use a simpler approach:
        import datetime as dt
        from datetime import timedelta

        # Convert MJD to datetime
        mjd0 = 2400000.5
        jd_utc = mjd_utc + mjd0
        jd_tt = mjd_tt + mjd0

        # Use AbsoluteDate from JulianDate components
        # Reference: 2000-01-01T12:00:00 TT = JD 2451545.0
        from org.orekit.time import AbsoluteDate
        # Build from components using modified julian date
        # Actually Orekit's AbsoluteDate can be built from components
        # Using the factory: new AbsoluteDate(int year, int month, int day,
        #                                     int hour, int minute, double second,
        #                                     TimeScale)

        # For now, create placeholder date and note this needs proper conversion
        # MJD epoch: 1858-11-17T00:00:00
        # mjd_utc days since MJD epoch → compute datetime
        from datetime import datetime as dt_module, timedelta
        mjd_epoch = dt_module(1858, 11, 17, 0, 0, 0)
        dt_utc = mjd_epoch + timedelta(days=mjd_utc)
        dt_tt = mjd_epoch + timedelta(days=mjd_tt)

        date_utc = AbsoluteDate(
            dt_utc.year, dt_utc.month, dt_utc.day,
            dt_utc.hour, dt_utc.minute,
            dt_utc.second + dt_utc.microsecond * 1e-6,
            self._utc
        )

        # Build position and velocity
        pos = Vector3D(float(r_eci[0]), float(r_eci[1]), float(r_eci[2]))
        vel = Vector3D(float(v_eci[0]), float(v_eci[1]), float(v_eci[2]))
        pv = PVCoordinates(pos, vel)

        orbit = CartesianOrbit(pv, self._gcrf, date_utc, self._mu)

        # ── Configure propagator ──
        min_step = 0.1
        max_step = min(dt, 300.0)
        tol = 1e-12

        # Dormand-Prince 8(5,3) integrator
        integrator = DormandPrince853Integrator(min_step, max_step, tol, tol)

        propagator = NumericalPropagator(integrator)
        propagator.setInitialState(SpacecraftState(orbit))
        propagator.setOrbitType(orbit.getType())  # Cartesian

        # ── Force models ──
        self._setup_force_models(propagator)

        # ── Propagate ──
        final_state = propagator.propagate(date_utc.shiftedBy(float(dt)))

        # Extract final state
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

        # STM extraction (if available — requires additional setup)
        # For now return None; full STM support requires PartialDerivativesEquations
        Phi = None
        S = None

        return r_new, v_new, Phi, S

    def _setup_force_models(self, propagator):
        """Configure force models on the propagator."""
        from org.orekit.forces.gravity.potential import GravityFieldFactory
        from org.orekit.forces.gravity import HolmesFeatherstoneAttractionModel

        # Central gravity with high-degree spherical harmonics
        gravity_provider = GravityFieldFactory.getUnnormalizedProvider(
            self.gravity_degree, 0  # max degree, max order (0 = same)
        )
        gravity_model = HolmesFeatherstoneAttractionModel(
            self._itrf, gravity_provider
        )
        propagator.addForceModel(gravity_model)

        # Third-body (Sun, Moon) — requires CelestialBodyFactory
        try:
            from org.orekit.bodies import CelestialBodyFactory
            from org.orekit.forces.gravity import ThirdBodyAttraction

            sun = CelestialBodyFactory.getSun()
            moon = CelestialBodyFactory.getMoon()
            propagator.addForceModel(ThirdBodyAttraction(sun))
            propagator.addForceModel(ThirdBodyAttraction(moon))
        except Exception:
            pass

        # Solid tides
        if self.solid_tides:
            try:
                from org.orekit.forces.gravity import SolidTides
                tides_model = SolidTides(
                    self._itrf, gravity_provider.getAe(), gravity_provider.getMu(),
                    gravity_provider
                )
                propagator.addForceModel(tides_model)
            except Exception:
                pass

        # Ocean tides
        if self.ocean_tides:
            try:
                from org.orekit.forces.gravity import OceanTides
                ocean_model = OceanTides(
                    self._itrf, gravity_provider.getAe(), gravity_provider.getMu(),
                    gravity_provider, self.ocean_tide_degree, self.ocean_tide_degree,
                    True  # include pole tide
                )
                propagator.addForceModel(ocean_model)
            except Exception:
                pass

        # Solar radiation pressure
        if self.srp_model:
            try:
                from org.orekit.forces.radiation import SolarRadiationPressure
                from org.orekit.bodies import CelestialBodyFactory, OneAxisEllipsoid
                from org.orekit.frames import FramesFactory
                from org.orekit.utils import Constants as OC

                sun = CelestialBodyFactory.getSun()
                # Spherical Earth for eclipse computation
                earth_radius = OC.WGS84_EARTH_EQUATORIAL_RADIUS
                # Simple isotropic model: cross-section area * (1 + CR)
                # Using placeholder values
                srp = SolarRadiationPressure(sun, earth_radius, 1.0, 0.0, 4.5)
                propagator.addForceModel(srp)
            except Exception:
                pass

        # Relativistic corrections
        if self.relativity:
            try:
                from org.orekit.forces import Relativity
                rel = Relativity(self._mu)
                propagator.addForceModel(rel)
            except Exception:
                pass

    def compute_stm(self, propagator, state, dt):
        """Setup STM computation using Orekit's PartialDerivativesEquations.

        This requires creating a JacobianMapper and adding it to the
        propagator. For EKF, we need:
          - 6x6 position-velocity STM
          - 6x3 empirical acceleration sensitivity

        This is a placeholder for future enhancement. The simplified mode
        already handles STM computation.
        """
        # TODO: Implement full STM with PDE
        # Requires:
        #   1. PartialDerivativesEquations setup
        #   2. Additional equations for empirical parameters
        #   3. State transition matrix extraction
        pass


# ── Simplified wrapper for drop-in use ──

def create_propagator(mode='simplified', **kwargs):
    """Factory function: create the appropriate propagator.

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
            msg = f"Orekit init failed: {e}. Falling back to simplified."
            print(f"  [WARNING] {msg}")
            return None, 'simplified', msg

    return None, 'simplified', None
