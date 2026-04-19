# CLAUDE.md

## Project Overview

**fit_analyze** — CLI tool that processes Garmin FIT files for cycling performance analysis on the Silver Strand route (San Diego). Outputs text summaries, gnuplot visualizations, and structured Chung method data for aerodynamic regression.

## Quick Reference

```bash
# Run analysis
python3 fit_analyze.py <file.fit> [--ftp 200] [--rest-hr 72] [--max-hr 189] [--plot]

# Override wind (skips KNZY auto-fetch)
python3 fit_analyze.py <file.fit> --wind-speed 12 --wind-dir SE

# Debug weather fetch
python3 knzy_test.py --date 2026-04-17 --time-utc 11:23
```

## Dependencies

- **Python 3.8+**
- **fitparse** (`pip install fitparse`) — FIT binary parsing
- **gnuplot** (optional) — visualization when `--plot` is used

No other external dependencies. Standard library only otherwise.

## Architecture

Single-script design (`fit_analyze.py`, ~800 lines). Key functions:

| Function | Purpose |
|----------|---------|
| `analyze()` | Core pipeline: parse FIT, compute metrics, detect Strand section |
| `despike_conservative()` | Cross-channel de-spike for metrics (cadence validates power, vice versa) |
| `despike()` | Aggressive de-spike for plot (speed-aware + Hampel filter) |
| `fetch_knzy_observations()` | Scrape NWS weather station HTML |
| `find_best_knzy_obs()` | Match observation to ride time (±120 min) |
| `normalized_power()` | Coggan NP from 30s rolling average |
| `headwind_mph()` | Vector wind component along route heading |
| `plot_strand()` | Generate gnuplot visualization |

Data flow: FIT file → fitparse → `analyze()` → text summary + plot data → file output

## Key Configuration (top of fit_analyze.py)

- `RIDER_MASS = 65.9` kg — update when weight changes
- `DEFAULT_FTP = 200` W — functional threshold power
- `DEFAULT_CDA = 0.26` m² — drag area estimate
- `LA_TZ` — uses `zoneinfo.ZoneInfo("America/Los_Angeles")` for auto DST; falls back to UTC-7 on Python < 3.9
- Strand GPS bounding boxes define section detection (~32.678N to ~32.593N)

## Coding Conventions

- Pure Python, minimal dependencies
- Constants at module top level (ALL_CAPS)
- Functions use snake_case
- Inline comments explain domain-specific formulas (cycling physics, meteorology)
- Graceful degradation: optional features (gnuplot, KNZY) fail with warnings, not crashes
- Active riding detection: `cadence > 10 rpm OR power > 20W`

## Output Files

Generated alongside the input FIT file:
- `<basename>_summary.txt` — text metrics report
- `<basename>_strand.png` — gnuplot chart (only with `--plot`)

## Known Limitations

- KNZY HTML scraping is fragile to page structure changes
- Timezone falls back to fixed PDT (UTC-7) on Python < 3.9 (no zoneinfo)
- GPS section detection depends on signal quality
- CdA is a single estimate, not regression-fitted (future: `chung_fit.py`)

## Git Notes

- `.gitignore` excludes: `*.fit`, `*_summary.txt`, `__pycache__/`, `fit_files/`
- FIT files are binary and large — kept local only
- Sample data lives in `fit_files/` (not tracked)
