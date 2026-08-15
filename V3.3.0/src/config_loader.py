"""Configuration loader for YAML-based pipeline config (V3.2.1).
Loads ALL parameters from YAML, returns a unified cfg dict.
CLI args override YAML defaults for dates/hours/satellite/gps mode.
"""
import yaml
from pathlib import Path


def load_config(config_path):
    """Load YAML config and return full dict."""
    with open(config_path, 'r', encoding='utf-8') as f:
        cfg = yaml.safe_load(f)
    return cfg


def apply_config_to_args(cfg, args):
    """Apply YAML config to argparse namespace. CLI args take precedence over YAML when explicitly set."""
    # Evaluation dates & arcs (CLI overrides YAML)
    if 'evaluation' in cfg:
        ev = cfg['evaluation']
        if ev.get('dates'):
            args.dates = ','.join(str(d) for d in ev['dates'])
        if ev.get('arc_hours'):
            args.hours = ','.join(str(h) for h in ev['arc_hours'])
        if ev.get('min_coverage') is not None:
            args.min_coverage = float(ev['min_coverage'])

    # Satellite (CLI --grace-id overrides YAML if not default 'C')
    if 'satellite' in cfg:
        sat = cfg['satellite']
        args.grace_id = sat.get('id', args.grace_id)
        if sat.get('mission'):
            args.mission = sat['mission']

    # GPS products
    if 'gps_products' in cfg:
        gp = cfg['gps_products']
        if gp.get('mode') == 'broadcast':
            args.broadcast = True

    # Code orbit
    if 'observations' in cfg:
        obs = cfg['observations']
        if obs.get('code_arc_hours', 0) > 0:
            args.code_arc_hours = float(obs['code_arc_hours'])

    return args


def get_cfg_param(cfg, *keys, default=None):
    """Get a nested config parameter with fallback default.
    e.g. get_cfg_param(cfg, 'ekf', 'sigma_phase', default=0.20)
    """
    d = cfg
    for k in keys[:-1]:
        if not isinstance(d, dict):
            return default
        d = d.get(k, {})
    return d.get(keys[-1], default) if isinstance(d, dict) else default


def build_pipeline_config(cfg):
    """Extract all pipeline parameters from YAML config into a flat dict."""
    return {
        # ── QC ──
        'qc_min_sv_coverage':    get_cfg_param(cfg, 'qc', 'min_sv_coverage', default=0.40),
        'qc_snr_l1_min':         get_cfg_param(cfg, 'qc', 'snr_l1_min', default=30),
        'qc_snr_l2_min':         get_cfg_param(cfg, 'qc', 'snr_l2_min', default=25),
        'qc_mp_threshold':       get_cfg_param(cfg, 'qc', 'mp_threshold', default=3.5),
        'qc_mw_std_corrupted':   get_cfg_param(cfg, 'qc', 'mw_std_corrupted', default=10.0),
        'qc_mw_std_noisy':       get_cfg_param(cfg, 'qc', 'mw_std_noisy', default=0.30),
        'qc_mw_std_unstable':    get_cfg_param(cfg, 'qc', 'mw_std_unstable', default=0.50),
        'qc_wl_fix_residual':    get_cfg_param(cfg, 'qc', 'wl_fix_residual', default=0.35),
        'qc_nl_fix_residual':    get_cfg_param(cfg, 'qc', 'nl_fix_residual', default=0.30),
        'qc_robust_reweight':    get_cfg_param(cfg, 'qc', 'robust_reweight', default=True),
        'qc_min_wl_epochs':      get_cfg_param(cfg, 'qc', 'min_wl_epochs', default=5),
        'qc_max_rejected_sv_ratio': get_cfg_param(cfg, 'qc', 'max_rejected_sv_ratio', default=0.80),

        # ── EKF ──
        'ekf_sigma_phase':       get_cfg_param(cfg, 'ekf', 'sigma_phase', default=0.20),
        'ekf_sigma_code':        get_cfg_param(cfg, 'ekf', 'sigma_code', default=0.30),
        'ekf_sigma_acc':         get_cfg_param(cfg, 'ekf', 'sigma_acc_process', default=1e-3),
        'ekf_tau_emp':           get_cfg_param(cfg, 'ekf', 'tau_emp', default=600.0),
        'ekf_sigma_emp_ss':      get_cfg_param(cfg, 'ekf', 'sigma_emp_ss', default=1e-8),
        'ekf_sigma_zwd_rw':      get_cfg_param(cfg, 'ekf', 'sigma_zwd_rw', default=1e-9),
        'ekf_chi2_short':        get_cfg_param(cfg, 'ekf', 'chi2_short', default=25),
        'ekf_chi2_long':         get_cfg_param(cfg, 'ekf', 'chi2_long', default=100),
        'ekf_el_min_deg':        get_cfg_param(cfg, 'ekf', 'el_min_deg', default=5.0),
        'ekf_clock_rw_short':    get_cfg_param(cfg, 'ekf', 'clock_rw_short', default=0.0004),
        'ekf_clock_rw_long':     get_cfg_param(cfg, 'ekf', 'clock_rw_long', default=0.001),
        'ekf_elev_exp_phase':    get_cfg_param(cfg, 'ekf', 'elev_exp_phase', default=1.0),
        'ekf_elev_exp_code_short': get_cfg_param(cfg, 'ekf', 'elev_exp_code_short', default=1.0),
        'ekf_elev_exp_code_long': get_cfg_param(cfg, 'ekf', 'elev_exp_code_long', default=0.70),
        'ekf_ar_min_epochs':     get_cfg_param(cfg, 'ekf', 'ar_min_epochs', default=6),
        'ekf_mw_max_epochs':     get_cfg_param(cfg, 'ekf', 'mw_max_epochs', default=200),
        'ekf_predict_only':     get_cfg_param(cfg, 'ekf', 'predict_only', default=False),
        'ekf_geo_from_orbit':   get_cfg_param(cfg, 'ekf', 'geo_from_orbit', default=False),
        'ekf_use_windup':        get_cfg_param(cfg, 'ekf', 'use_phase_windup', default=True),
        'ekf_use_relativity':    get_cfg_param(cfg, 'ekf', 'use_relativity', default=True),
        'ekf_use_cycle_slip':    get_cfg_param(cfg, 'ekf', 'use_cycle_slip', default=True),
        'ekf_dynamics_mode':     get_cfg_param(cfg, 'ekf', 'dynamics_mode', default='simplified'),

        # ── GN ──
        'gn_max_iter':           get_cfg_param(cfg, 'gn_loop', 'max_iter', default=6),
        'gn_damping':            get_cfg_param(cfg, 'gn_loop', 'damping', default=0.5),
        'gn_prior_r0':           get_cfg_param(cfg, 'gn_loop', 'prior_r0', default=1.0),
        'gn_prior_v0':           get_cfg_param(cfg, 'gn_loop', 'prior_v0', default=0.01),
        'gn_prior_emp':          get_cfg_param(cfg, 'gn_loop', 'prior_emp', default=1e-7),

        # ── Gravity ──
        'gravity_nmax':          get_cfg_param(cfg, 'gravity', 'nmax', default=150),

        # ── Dynamics ──
        'dynamics_mode':         get_cfg_param(cfg, 'dynamics', 'dynamics_mode', default='simplified'),
        'srp_model':             get_cfg_param(cfg, 'dynamics', 'srp_model', default='isotropic'),
        'drag_model':            get_cfg_param(cfg, 'dynamics', 'drag_model', default='exponential'),
    }
