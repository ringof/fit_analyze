#!/usr/bin/env python3
"""
chung_fit.py — CdA/Crr regression from accumulated ride DATA lines.

Reads structured DATA lines produced by fit_analyze.py and solves for
aerodynamic drag area (CdA) using the Chung method power balance:

    P = 0.5 × ρ × CdA × (v + v_hw)² × v  +  Crr × m × g × v

This is linear in CdA and Crr, so ordinary least squares (OLS) gives
a direct solution — no iterative solver or external library needed.

Usage:
    python3 chung_fit.py <data_file>
    python3 chung_fit.py <data_file> --solve-crr
    grep "DATA:" *_summary.txt | python3 chung_fit.py -
"""

import argparse
import math
import re
import sys

# ── constants ────────────────────────────────────────────────────────────────

RIDER_MASS  = 65.9      # kg
BIKE_MASS   = 13.6      # kg (30 lbs)
TOTAL_MASS  = RIDER_MASS + BIKE_MASS  # 79.5 kg
G           = 9.80665   # m/s²
DEFAULT_RHO = 1.225     # kg/m³ — sea level ISA
DEFAULT_CRR = 0.006     # Schwalbe Marathon GreenGuard at ~70 psi
MPH_TO_MS   = 0.44704   # 1 mph in m/s

# Sanity bounds (warn, not reject)
CDA_LO, CDA_HI = 0.20, 0.45   # m²
CRR_LO, CRR_HI = 0.002, 0.012

# t-distribution critical values for 95% two-tailed CI (α/2 = 0.025)
_T_TABLE = {
     1: 12.706,  2:  4.303,  3:  3.182,  4:  2.776,  5:  2.571,
     6:  2.447,  7:  2.365,  8:  2.306,  9:  2.262, 10:  2.228,
    11:  2.201, 12:  2.179, 13:  2.160, 14:  2.145, 15:  2.131,
    16:  2.120, 17:  2.110, 18:  2.101, 19:  2.093, 20:  2.086,
    25:  2.060, 30:  2.042, 40:  2.021, 50:  2.009, 60:  2.000,
    80:  1.990, 100: 1.984, 120: 1.980,
}

def t_critical(df):
    """95% two-tailed t critical value for given degrees of freedom."""
    if df < 1:
        return float('inf')
    if df in _T_TABLE:
        return _T_TABLE[df]
    if df > 120:
        return 1.960  # normal approximation
    # interpolate between nearest keys
    keys = sorted(_T_TABLE.keys())
    for i in range(len(keys) - 1):
        if keys[i] < df < keys[i + 1]:
            lo, hi = keys[i], keys[i + 1]
            frac = (df - lo) / (hi - lo)
            return _T_TABLE[lo] + frac * (_T_TABLE[hi] - _T_TABLE[lo])
    return 1.960


# ── parsing ──────────────────────────────────────────────────────────────────

def parse_data_line(line):
    """Parse a DATA line into a dict, or return None if not a DATA line.

    Handles lines like:
      DATA: 2026-04-20 | pwr=109.9W | spd=6.215m/s | NP=117W | wind=8mph_W
            | hdg=150deg | hw=-4.0mph | rho=1.225 | Crr=0.0060 | EI=0.776 | CLEAN

    The Crr field may be absent in older summaries (defaults to DEFAULT_CRR).
    Lines may be prefixed with a filename from grep output.
    """
    # strip grep-style filename prefix  (path/to/file.txt:  DATA: ...)
    m = re.search(r'DATA:\s*(\S+)', line)
    if not m:
        return None

    date = m.group(1)

    def field(name, pattern=r'[-+]?\d+\.?\d*'):
        f = re.search(name + r'=(' + pattern + r')', line)
        return float(f.group(1)) if f else None

    pwr = field('pwr')
    spd = field('spd')        # m/s
    np_val = field('NP')
    hw  = field('hw')         # mph
    rho = field('rho')
    crr = field('Crr')
    ei  = field('EI')
    hdg = field('hdg')

    # wind string (e.g. "8mph_W")
    wm = re.search(r'wind=(\d+)mph_(\w+)', line)
    wind_spd = int(wm.group(1)) if wm else None
    wind_dir = wm.group(2) if wm else None

    # status is the last pipe-delimited token
    status = line.rstrip().rsplit('|', 1)[-1].strip()

    if pwr is None or spd is None or hw is None:
        return None  # incomplete DATA line

    return {
        'date':     date,
        'pwr':      pwr,
        'spd':      spd,       # m/s
        'np':       np_val,
        'wind_spd': wind_spd,
        'wind_dir': wind_dir,
        'hdg':      hdg,
        'hw':       hw,        # mph (positive = headwind)
        'rho':      rho if rho is not None else DEFAULT_RHO,
        'crr':      crr if crr is not None else DEFAULT_CRR,
        'ei':       ei,
        'status':   status,
    }


def load_data(filepath, include_stops=False):
    """Read DATA lines from a file (or stdin if filepath is '-').

    Returns list of parsed ride dicts.  Skips NO-STRAND lines.
    By default only CLEAN rides; --include-stops adds STOPS-CHECK.
    """
    allowed = {'CLEAN'}
    if include_stops:
        allowed.add('STOPS-CHECK')

    if filepath == '-':
        lines = sys.stdin.readlines()
        source = 'stdin'
    else:
        with open(filepath) as f:
            lines = f.readlines()
        source = filepath

    rides = []
    skipped = {'no_strand': 0, 'stops': 0}

    for line in lines:
        d = parse_data_line(line)
        if d is None:
            continue
        if d['status'] == 'NO-STRAND':
            skipped['no_strand'] += 1
            continue
        if d['status'] not in allowed:
            skipped['stops'] += 1
            continue
        rides.append(d)

    return rides, skipped, source


# ── regression ───────────────────────────────────────────────────────────────

def build_regression_terms(rides, mass, rho_override=None):
    """Compute per-ride regression terms from the power balance equation.

    For each ride:
        A_i = 0.5 × ρ × |v_air| × v_air × v    (aero coefficient of CdA)
        R_i = m × g × v                           (rolling coefficient of Crr)
        P_i = measured average power

    v_air = v + v_hw (rider ground speed + headwind component).
    Uses signed form |v_air|×v_air instead of v_air² to correctly handle
    the rare case where tailwind exceeds rider speed (assistive drag).

    Returns (A, R, P, dates) as parallel lists.
    """
    A, R, P, dates = [], [], [], []
    for r in rides:
        v = r['spd']                              # m/s
        v_hw = r['hw'] * MPH_TO_MS                # headwind in m/s
        rho = rho_override if rho_override else r['rho']

        v_air = v + v_hw                          # relative airspeed
        a_i = 0.5 * rho * abs(v_air) * v_air * v # signed aero term
        r_i = mass * G * v                        # rolling term
        p_i = r['pwr']                            # measured power

        A.append(a_i)
        R.append(r_i)
        P.append(p_i)
        dates.append(r['date'])

    return A, R, P, dates


def ols_cda_only(A, R, P, crr):
    """Single-parameter OLS: solve for CdA with known Crr.

    P_adj = P - Crr × R   (subtract rolling resistance)
    CdA = Σ(A·P_adj) / Σ(A²)
    """
    n = len(A)
    P_adj = [P[i] - crr * R[i] for i in range(n)]

    sum_A2   = sum(a * a for a in A)
    sum_A_Pa = sum(A[i] * P_adj[i] for i in range(n))

    if sum_A2 == 0:
        raise ValueError("All aero terms are zero — no speed/wind variation")

    cda = sum_A_Pa / sum_A2

    # predicted power and residuals
    pred = [cda * A[i] + crr * R[i] for i in range(n)]
    resid = [P[i] - pred[i] for i in range(n)]

    # R² (relative to mean of P)
    p_mean = sum(P) / n
    ss_res = sum(r * r for r in resid)
    ss_tot = sum((P[i] - p_mean) ** 2 for i in range(n))
    r_sq = 1 - ss_res / ss_tot if ss_tot > 0 else float('nan')

    # residual standard error, standard error of CdA, 95% CI
    df = n - 1
    rse = math.sqrt(ss_res / df) if df > 0 else float('nan')
    se_cda = rse / math.sqrt(sum_A2) if sum_A2 > 0 else float('nan')
    t_crit = t_critical(df)
    ci_lo = cda - t_crit * se_cda
    ci_hi = cda + t_crit * se_cda

    return {
        'cda': cda, 'se_cda': se_cda, 'ci': (ci_lo, ci_hi),
        'r_sq': r_sq, 'rse': rse, 'pred': pred, 'resid': resid,
        'df': df, 'n': n,
    }


def ols_cda_crr(A, R, P):
    """Two-parameter OLS: solve for CdA and Crr simultaneously.

    Normal equations (2×2):
        [Σ(A²)    Σ(A·R)] [CdA]   [Σ(A·P)]
        [Σ(A·R)   Σ(R²) ] [Crr] = [Σ(R·P)]
    """
    n = len(A)

    sum_A2  = sum(A[i] * A[i] for i in range(n))
    sum_R2  = sum(R[i] * R[i] for i in range(n))
    sum_AR  = sum(A[i] * R[i] for i in range(n))
    sum_AP  = sum(A[i] * P[i] for i in range(n))
    sum_RP  = sum(R[i] * P[i] for i in range(n))

    det = sum_A2 * sum_R2 - sum_AR * sum_AR
    if abs(det) < 1e-20:
        raise ValueError(
            "Singular matrix — rides have too little speed/wind diversity.\n"
            "Need more rides with varied conditions for two-parameter fit."
        )

    cda = (sum_R2 * sum_AP - sum_AR * sum_RP) / det
    crr = (sum_A2 * sum_RP - sum_AR * sum_AP) / det

    # predicted power and residuals
    pred = [cda * A[i] + crr * R[i] for i in range(n)]
    resid = [P[i] - pred[i] for i in range(n)]

    # R²
    p_mean = sum(P) / n
    ss_res = sum(r * r for r in resid)
    ss_tot = sum((P[i] - p_mean) ** 2 for i in range(n))
    r_sq = 1 - ss_res / ss_tot if ss_tot > 0 else float('nan')

    # residual SE, parameter standard errors from (X'X)^-1
    df = n - 2
    rse = math.sqrt(ss_res / df) if df > 0 else float('nan')

    # (X'X)^-1 diagonal elements
    inv_diag_cda = sum_R2 / det
    inv_diag_crr = sum_A2 / det

    se_cda = rse * math.sqrt(inv_diag_cda) if inv_diag_cda >= 0 else float('nan')
    se_crr = rse * math.sqrt(inv_diag_crr) if inv_diag_crr >= 0 else float('nan')

    t_crit = t_critical(df)
    ci_cda = (cda - t_crit * se_cda, cda + t_crit * se_cda)
    ci_crr = (crr - t_crit * se_crr, crr + t_crit * se_crr)

    return {
        'cda': cda, 'se_cda': se_cda, 'ci_cda': ci_cda,
        'crr': crr, 'se_crr': se_crr, 'ci_crr': ci_crr,
        'r_sq': r_sq, 'rse': rse, 'pred': pred, 'resid': resid,
        'df': df, 'n': n,
    }


# ── output ───────────────────────────────────────────────────────────────────

def format_results(result, rides, dates, mode, source, skipped, crr_used,
                    mass):
    """Format regression output matching fit_analyze.py visual style."""
    lines = []
    w = lines.append

    w("=" * 62)
    w("CHUNG METHOD REGRESSION — CdA from ride data")
    w("=" * 62)

    # ── input section
    n = result['n']
    w("")
    w("── INPUT " + "─" * 53)

    w(f"  Data file:       {source}")

    skip_parts = []
    if skipped['stops']:
        skip_parts.append(f"{skipped['stops']} STOPS-CHECK skipped")
    if skipped['no_strand']:
        skip_parts.append(f"{skipped['no_strand']} NO-STRAND skipped")
    skip_note = f"  ({', '.join(skip_parts)})" if skip_parts else ""
    w(f"  Rides loaded:    {n}{skip_note}")

    if mode == 'cda_only':
        w(f"  Mode:            CdA-only (Crr fixed at {crr_used:.4f})")
    else:
        w(f"  Mode:            CdA + Crr (two-parameter)")
    if mass == TOTAL_MASS:
        w(f"  Total mass:      {mass:.1f} kg"
          f" (rider {RIDER_MASS:.1f} + bike {BIKE_MASS:.1f})")
    else:
        w(f"  Total mass:      {mass:.1f} kg (override)")

    # ── result section
    w("")
    w("── RESULT " + "─" * 52)

    cda = result['cda']
    w(f"  CdA              {cda:.3f} m²")

    if mode == 'cda_only':
        ci = result['ci']
        w(f"  95% CI           [{ci[0]:.3f}, {ci[1]:.3f}]")
        w(f"  Std error        {result['se_cda']:.3f} m²")
    else:
        ci_cda = result['ci_cda']
        w(f"  95% CI           [{ci_cda[0]:.3f}, {ci_cda[1]:.3f}]")
        w(f"  Std error        {result['se_cda']:.3f} m²")
        w("")
        crr = result['crr']
        w(f"  Crr              {crr:.4f}")
        ci_crr = result['ci_crr']
        w(f"  95% CI           [{ci_crr[0]:.4f}, {ci_crr[1]:.4f}]")
        w(f"  Std error        {result['se_crr']:.4f}")
        # compare to manufacturer value
        delta = crr - DEFAULT_CRR
        w(f"  vs manufacturer  {DEFAULT_CRR:.4f} (delta {delta:+.4f})")

    w("")
    r_sq = result['r_sq']
    w(f"  R²               {r_sq:.3f}" if not math.isnan(r_sq) else
      f"  R²               n/a (insufficient data)")
    w(f"  Residual SE      {result['rse']:.1f} W" if not math.isnan(result['rse']) else
      f"  Residual SE      n/a")

    # ── sanity warnings
    warnings = []
    if cda < CDA_LO or cda > CDA_HI:
        warnings.append(
            f"  ** CdA {cda:.3f} outside expected range"
            f" [{CDA_LO:.2f}, {CDA_HI:.2f}] — check data quality")
    if mode == 'cda_crr':
        crr = result['crr']
        if crr < CRR_LO or crr > CRR_HI:
            warnings.append(
                f"  ** Crr {crr:.4f} outside expected range"
                f" [{CRR_LO:.3f}, {CRR_HI:.3f}] — check data quality")
    if cda < 0:
        warnings.append("  ** CdA is NEGATIVE — bad data, remove outlier rides")
    if mode == 'cda_crr' and result['crr'] < 0:
        warnings.append("  ** Crr is NEGATIVE — bad data, remove outlier rides")

    if warnings:
        w("")
        for wn in warnings:
            w(wn)

    # ── per-ride residuals
    rse = result['rse']
    w("")
    w("── PER-RIDE RESIDUALS " + "─" * 41)
    w("  Date         Pwr   Spd    HW     Pred   Resid")
    w("               (W)   (mph)  (mph)  (W)    (W)")

    for i, r in enumerate(rides):
        spd_mph = r['spd'] / MPH_TO_MS
        pred_w = result['pred'][i]
        res_w = result['resid'][i]
        flag = "  !!" if (not math.isnan(rse) and rse > 0
                          and abs(res_w) > 2 * rse) else ""
        w(f"  {dates[i]}  {r['pwr']:5.1f}  {spd_mph:5.1f}  {r['hw']:5.1f}"
          f"  {pred_w:5.1f}  {res_w:+5.1f}{flag}")

    if not math.isnan(rse) and rse > 0:
        w("  !! = residual > 2 SE — investigate ride conditions")

    # ── reference table
    w("")
    w("── REFERENCE " + "─" * 49)
    w("  0.24-0.28  drops / aggressive road position")
    w("  0.28-0.32  hoods position, typical road cycling")
    w("  0.32-0.40  relaxed / upright position")

    # ── note for small samples
    if n <= 4:
        w("")
        w(f"  Note: only {n} ride{'s' if n != 1 else ''} — confidence intervals"
          f" are wide.")
        w("  More rides with varied wind/speed will tighten the estimate.")

    w("")
    return '\n'.join(lines)


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='CdA/Crr regression from fit_analyze DATA lines (Chung method)')
    parser.add_argument('data_file',
        help='File with DATA lines, or - for stdin')
    parser.add_argument('--include-stops', action='store_true',
        help='Include STOPS-CHECK rides (default: CLEAN only)')
    parser.add_argument('--solve-crr', action='store_true',
        help='Two-parameter mode: solve for both CdA and Crr')
    parser.add_argument('--crr', type=float, default=None,
        help=f'Override Crr for CdA-only mode (default: {DEFAULT_CRR})')
    parser.add_argument('--mass', type=float, default=None,
        help=f'Override total mass in kg (default: {TOTAL_MASS:.1f})')
    parser.add_argument('--rho', type=float, default=None,
        help='Override air density for all rides (default: per-ride from data)')
    args = parser.parse_args()

    crr_used = args.crr if args.crr is not None else DEFAULT_CRR
    mass = args.mass if args.mass is not None else TOTAL_MASS

    # ── load data
    try:
        rides, skipped, source = load_data(args.data_file, args.include_stops)
    except FileNotFoundError:
        print(f"Error: file not found: {args.data_file}", file=sys.stderr)
        sys.exit(1)

    if len(rides) == 0:
        print("Error: no usable DATA lines found.", file=sys.stderr)
        if skipped['no_strand']:
            print(f"  ({skipped['no_strand']} NO-STRAND lines skipped)",
                  file=sys.stderr)
        if skipped['stops']:
            print(f"  ({skipped['stops']} STOPS-CHECK lines skipped"
                  f" — try --include-stops)", file=sys.stderr)
        sys.exit(1)

    if len(rides) < 2:
        print("Error: need at least 2 rides for regression.", file=sys.stderr)
        print("  Accumulate more rides with varied wind/speed conditions.",
              file=sys.stderr)
        sys.exit(1)

    if args.solve_crr and len(rides) < 3:
        print("Error: need at least 3 rides for two-parameter fit (CdA + Crr).",
              file=sys.stderr)
        print(f"  Have {len(rides)} rides. Use CdA-only mode or add more data.",
              file=sys.stderr)
        sys.exit(1)

    # ── build regression terms
    A, R, P, dates = build_regression_terms(rides, mass, args.rho)

    # ── solve
    mode = 'cda_crr' if args.solve_crr else 'cda_only'
    try:
        if mode == 'cda_only':
            result = ols_cda_only(A, R, P, crr_used)
        else:
            result = ols_cda_crr(A, R, P)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    # ── check for negative coefficients
    if result['cda'] < 0:
        print("Error: regression produced negative CdA — data is inconsistent.",
              file=sys.stderr)
        print("  Rides in dataset:", file=sys.stderr)
        for i, r in enumerate(rides):
            print(f"    {dates[i]}  pwr={r['pwr']:.1f}W  spd={r['spd']:.3f}m/s"
                  f"  hw={r['hw']:.1f}mph", file=sys.stderr)
        print("  Try removing outlier rides.", file=sys.stderr)
        sys.exit(1)

    if mode == 'cda_crr' and result['crr'] < 0:
        print("Warning: regression produced negative Crr — likely insufficient"
              " data diversity.", file=sys.stderr)
        print("  Falling back to CdA-only mode with manufacturer Crr.",
              file=sys.stderr)
        mode = 'cda_only'
        result = ols_cda_only(A, R, P, DEFAULT_CRR)
        crr_used = DEFAULT_CRR

    # ── output
    report = format_results(result, rides, dates, mode, source, skipped, crr_used,
                            mass)
    print(report)


if __name__ == '__main__':
    main()
