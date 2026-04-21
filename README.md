# fit_analyze

Cycling FIT file preprocessor for performance tracking on the Silver Strand route,
San Diego. Produces a concise text summary suitable for pasting into a chat or log,
and accumulates structured data for aerodynamic (CdA/Crr) regression via the
Chung method.

## Files

| File | Purpose |
|------|---------|
| `fit_analyze.py` | Main analysis script |
| `chung_fit.py` | CdA/Crr regression from accumulated ride data |
| `knzy_test.py` | Standalone KNZY weather fetch debugger |

## Dependencies

```bash
pip install fitparse
```

Python 3.8+ standard library otherwise (no numpy, no pandas).

## Usage

### Basic — auto-fetches wind from KNZY

```bash
python3 fit_analyze.py ride.fit
```

### With manual wind override

```bash
python3 fit_analyze.py ride.fit --wind-speed 12 --wind-dir SE
```

### Full options

```bash
python3 fit_analyze.py ride.fit \
  --ftp 200 \
  --rest-hr 72 \
  --max-hr 189 \
  --wind-speed 12 \
  --wind-dir SE \
  --route-heading 150
```

A `.txt` summary is saved alongside the `.fit` file automatically.

## Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `fit_file` | required | Path to `.fit` file |
| `--ftp` | 200 | Functional threshold power (W) |
| `--rest-hr` | 72 | Resting heart rate (bpm) |
| `--max-hr` | 189 | Maximum heart rate (bpm) |
| `--wind-speed` | auto | KNZY wind speed (mph). If omitted, fetched automatically. |
| `--wind-dir` | auto | KNZY wind direction (e.g. SE, NW, E) — direction wind blows FROM |
| `--route-heading` | 150 | Travel direction in degrees (Strand southbound ≈ 150°) |

## What it produces

### Full-ride metrics
- Duration, distance, speed, ascent, TSS
- Time in motion vs coasting
- Heart rate zones (Karvonen method)
- Cadence distribution (active periods only — see Active Detection below)
- Power distribution, normalized power (Garmin session value), intensity factor, W/kg
- Mechanical kJ and caloric estimate (see kJ note below)

### Silver Strand section (miles 3.8–11.0, geolocation-anchored)
The open-water section is detected by GPS bounding boxes rather than distance,
so it works correctly regardless of ferry wait duration or minor route variation.
This section is the only part of the route with unobstructed wind fetch, making
it the appropriate domain for aerodynamic analysis and the most consistent basis
for ride-to-ride comparison.

Reported separately:
- Avg power, speed, HR, cadence median
- Section normalized power (computed from per-second data — see NP note below)
- Efficiency Index (NP / avg active HR)
- Stop detection: flags any stop >30s that may affect Chung data quality

### Wind conditions (KNZY / NAS North Island)
Attempts to auto-fetch the KNZY observation table from the National Weather Service.
Selects the observation closest to (and preceding) the ride start time within a
±2 hour window. Calculates headwind component along the route heading and an
estimated aerodynamic power penalty.

KNZY fetching is HTML-scrape based and therefore fragile. If the fetch fails
(network error, page structure change), the script continues without wind data
and reports the failure reason. Treat wind context as useful but approximate.

### Chung method data block
A structured single-line DATA record per ride, using Strand section values:

```
DATA: 2026-04-17 | pwr=122.1W | spd=5.452m/s | NP=138W | wind=12mph_SE | hdg=150deg | hw=11.6mph | rho=1.225 | EI=0.895 | CLEAN
```

Accumulate these lines across rides and feed them to `chung_fit.py` for
CdA/Crr regression (see [CdA regression](#cda-regression-chung_fitpy) below).

## Configuration constants

Edit directly in `fit_analyze.py`:

```python
RIDER_MASS      = 65.9    # kg — update as weight changes (slow parameter)
DEFAULT_FTP     = 200     # W — update after FTP test
DEFAULT_MAX_HR  = 189     # bpm — observed maximum
DEFAULT_REST_HR = 72      # bpm — morning resting average
DEFAULT_CDA     = 0.26    # m² — back-calculated from one ride; refine via Chung
DEFAULT_CRR     = 0.006   # Schwalbe Marathon GreenGuard at ~70 psi
BIKE_MASS       = 13.6    # kg (30 lbs)

# Active riding detection thresholds
ACTIVE_CAD_MIN  = 10      # rpm
ACTIVE_PWR_MIN  = 20      # W — fallback when cadence sensor drops out

# Silver Strand GPS bounding boxes (±0.005° ≈ 500m)
STRAND_ENTRY_LAT =  32.678
STRAND_ENTRY_LON = -117.175
STRAND_EXIT_LAT  =  32.593
STRAND_EXIT_LON  = -117.122
STRAND_BOX_DEG   =  0.005
```

`RIDER_MASS` is not a command-line argument — it changes slowly enough that
editing the constant directly is appropriate.

## Notes on specific metrics

### Normalized power (NP)

Two NP values appear in the output with different sources:

**Full-ride NP** is taken directly from the Garmin session summary
(`normalized_power` field in the FIT file). This is Garmin's own calculation
over the complete activity.

**Strand section NP** is computed from per-second power data using the standard
Coggan definition: 30-second rolling mean of power raised to the 4th power,
then the mean of those values raised to 1/4. This is necessary because the
Garmin session value covers the full ride, not the section. The two values
will differ slightly due to scope and floating-point rounding. The docstring
in `normalized_power()` documents the implementation.

### Efficiency Index (EI)

`NP / avg active HR` in W/bpm. A derived metric relating power output to
cardiovascular cost. On a controlled course with consistent conditions, EI
rising over successive rides reflects cardiovascular adaptation. Wind affects
HR independently of power, so headwind conditions will produce lower EI values
than calm conditions at equivalent effort.

Full-ride EI uses Garmin session NP. Strand section EI uses the section-computed
NP. The Strand value is more consistent ride to ride because it excludes
sheltered urban sections with different HR dynamics.

EI is one signal among several. It should be interpreted alongside power, HR
zones, cadence, and wind conditions rather than in isolation.

### De-spiking (plot visualization only)

Power and cadence traces contain near-zero values from light pedaling,
momentary recovery, and natural force variation through the stroke cycle.
Analysis of the data shows most of these are real riding behavior, not sensor
dropouts — only ~30 per ride are isolated single-second zeros consistent with
a sampling rate mismatch between the Edge 500's 1Hz recording and the
PowerTap hub's 4Hz transmission. These are too few to affect calculated metrics.

**All metrics use raw uncorrected data.** This avoids systematic power
inflation that would bias the Chung CdA regression.

**Plot only:** `despike()` applies aggressive filtering (speed-aware dropout
replacement + Hampel filter) for visual clarity. The plot subtitle notes
the filter and correction count.

### Active detection

A record is classified as active if `cadence > 10 rpm OR power > 20W`. The
OR condition handles brief cadence sensor dropouts during genuine pedaling,
which would otherwise be misclassified as coasting and distort every
"pedaling only" metric. All statistics labelled "active" or "pedaling only"
use this definition.

### kJ and caloric estimates

`Mechanical work (kJ)` comes from the FIT file's `total_work` session field,
which accumulates power meter output integrated over time.

The caloric estimate uses the approximation that mechanical kJ ≈ net exercise
kcal. This holds because cycling mechanical efficiency is approximately 24%
and the unit conversions cancel to within a few percent. It is an approximation.

The Garmin caloric estimate uses heart rate to infer metabolic rate and includes
resting metabolic contribution. As cardiovascular fitness improves and HR drops
for the same power output, the Garmin estimate will fall even though actual
mechanical work is unchanged. The kJ figure is not affected by this.

### CdA (coefficient of drag area)

Set to 0.26 m², back-calculated from one ride's power, speed, and wind data.
This is a rough estimate. Use `chung_fit.py` with accumulated ride data to
produce a better-constrained value. The aero penalty figure in the output is
approximate context, not a precise measurement.

## Wind source

[KNZY — NAS North Island, Coronado](https://forecast.weather.gov/data/obhistory/KNZY.html)

Located near the Coronado ferry landing, close to the start of the route.
Observations are published at :52 past each hour. The script selects the
observation with the smallest time offset within a −120 to +30 minute window
relative to ride start (UTC converted to Pacific local time via `zoneinfo`).

The KNZY page is fetched by HTTP with HTML parsing. Page structure changes
will break the parser. Use `knzy_test.py` to diagnose fetch issues.

## Debugging KNZY fetch

```bash
python3 knzy_test.py --date 2026-04-17 --time-utc 11:23
```

Shows the full fetch, HTML structure diagnosis, all parsed observations,
and time-matching logic. Useful when auto-fetch produces unexpected results.

## Known limitations

- KNZY fetching is HTML-scrape based and will break if NWS changes page structure.
- Timezone handling uses `zoneinfo` (Python 3.9+) for automatic PDT/PST. Falls back to fixed UTC−7 on older Python.
- CdA regression needs 5+ rides with varied wind for tight confidence intervals.
- Section detection depends on GPS quality; poor GPS near the entry box may
  cause the section to be missed or mis-anchored.
- Active detection (cad > 10 or pwr > 20W) is a heuristic and may misclassify
  some records during unusual conditions.

## CdA regression (chung_fit.py)

Companion script that solves for aerodynamic drag area (CdA) from accumulated
DATA lines using the Chung method power balance equation:

```
P = 0.5 × ρ × CdA × (v + v_hw)² × v  +  Crr × m × g × v
```

This is linear in CdA and Crr, so ordinary least squares gives a direct solution.
No external libraries required.

### Usage

```bash
# CdA-only mode (Crr fixed at manufacturer value 0.006)
grep "DATA:" fit_files/*_summary.txt | python3 chung_fit.py -

# From a data file
python3 chung_fit.py ride_log.txt

# Two-parameter mode (solve for both CdA and Crr)
python3 chung_fit.py ride_log.txt --solve-crr

# Include STOPS-CHECK rides (default: CLEAN only)
python3 chung_fit.py ride_log.txt --include-stops

# Override Crr or total mass
python3 chung_fit.py ride_log.txt --crr 0.005
python3 chung_fit.py ride_log.txt --mass 82.0
```

Needs 2+ rides for CdA-only, 3+ for `--solve-crr`. Confidence intervals tighten
significantly with 5+ rides across varied wind/speed conditions. Uses average power
(not NP) — the power balance describes steady-state equilibrium.

## Future work

- Temperature correction to air density (RHO) from KNZY temperature field
- Configurable route profiles beyond Silver Strand
- Remove Python 3.9+ requirement for zoneinfo (backports exist but add a dependency)
