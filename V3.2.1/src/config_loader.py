"""Configuration loader for YAML-based pipeline config (V3.2.0 container)."""
import yaml
from pathlib import Path


def load_config(config_path):
    """Load YAML config and resolve paths / expand patterns."""
    with open(config_path, 'r', encoding='utf-8') as f:
        cfg = yaml.safe_load(f)
    return cfg


def apply_config_to_args(cfg, args):
    """Apply YAML config to argparse namespace, CLI args take precedence."""
    # Evaluation dates & arcs
    if 'evaluation' in cfg:
        ev = cfg['evaluation']
        if ev.get('dates') and not getattr(args, '_dates_set', False):
            args.dates = ','.join(str(d) for d in ev['dates'])
        if ev.get('arc_hours') and not getattr(args, '_hours_set', False):
            args.hours = ','.join(str(h) for h in ev['arc_hours'])
        if ev.get('min_coverage') is not None:
            setattr(args, 'min_coverage', ev['min_coverage'])

    # Satellite
    if 'satellite' in cfg:
        sat = cfg['satellite']
        if sat.get('id') and args.grace_id == 'C':
            args.grace_id = sat['id']

    # GPS products
    if 'gps_products' in cfg:
        gp = cfg['gps_products']
        if gp.get('mode') == 'broadcast':
            args.broadcast = True

    # Code orbit
    if 'observations' in cfg:
        obs = cfg['observations']
        if obs.get('code_arc_hours', 0) > 0:
            setattr(args, 'code_arc_hours', obs['code_arc_hours'])
        if not obs.get('gnv1b_pattern'):
            args.skip_code_orbit = False  # force self-consistent mode

    return args
