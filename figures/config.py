"""
Configuration: parameter names, ranges, output schema, time presets.
All angular parameters use pi-normalized representation (×π = radians).
Gamma is also a scan parameter (was hardcoded 1.0).
"""
GAMMA = 1.0                      # default when not scanning/fixed

PARAM_RANGES = {
    'Gamma':    (0.01, 30.0),
    'Delta':    (-30.0, 30.0),
    'theta_a1': (0.0,   2.0),
    'theta_a2': (0.0,   2.0),
    'theta_b2': (0.0,   2.0),
    'theta_b3': (0.0,   2.0),
    'phi1':     (0.0,   2.0),
    'phi2':     (0.0,   2.0),
    'phi3':     (0.0,   2.0),
}

PARAM_NAMES = ['Gamma', 'Delta', 'theta_a1', 'theta_a2', 'theta_b2',
               'theta_b3', 'phi1', 'phi2', 'phi3']

CHANNELS = ['T1', 'R1', 'T2', 'R2', 'T3', 'R3']

OUTPUT_COLS = PARAM_NAMES + ['score', 'channel'] + CHANNELS

# ── Time presets: (scan_fraction, de_pop, de_iter_per_seed) ──
# scan_fraction: portion of total time given to scan
# de_pop, de_iter_per_seed: tuned for the optimize phase
TIME_PRESETS = {
    600:      (0.20, 20,  200),   # 10 min
    3600:     (0.15, 30,  600),   # 1 h
    18000:    (0.10, 50,  2000),  # 5 h
    36000:    (0.08, 80,  4000),  # 10 h
}

def preset_for(seconds: float) -> tuple:
    """Return (scan_fraction, de_pop, de_iter_per_seed) for given time budget."""
    best = None
    for t, vals in sorted(TIME_PRESETS.items()):
        if seconds <= t:
            return vals
        best = vals
    return best  # largest preset
