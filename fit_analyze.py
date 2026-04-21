#!/usr/bin/env python3
"""
fit_analyze.py — Cycling FIT file preprocessor
Usage: python3 fit_analyze.py <file.fit> [--ftp 200] [--rest-hr 72] [--max-hr 189]
       [--wind-speed 12] [--wind-dir SE] [--route-heading 150] [--plot]

If --wind-speed/--wind-dir are omitted, the script automatically fetches
the KNZY (NAS North Island) observation closest to the ride start time.

Silver Strand open section is detected by geolocation, not distance,
so it works correctly regardless of ferry wait or route variation.
"""

import sys
import argparse
import statistics
import math
import re
import calendar
import urllib.request
from datetime import datetime, timezone, timedelta

try:
    from fitparse import FitFile
except ImportError:
    print("ERROR: fitparse not installed. Run: pip install fitparse")
    sys.exit(1)

# ── profile defaults ─────────────────────────────────────────────────────────
DEFAULT_FTP     = 200
DEFAULT_REST_HR = 72
DEFAULT_MAX_HR  = 189
RIDER_MASS      = 65.9   # kg

# ── physics constants ────────────────────────────────────────────────────────
RHO         = 1.225  # air density kg/m³ at sea level
DEFAULT_CDA = 0.26   # m² — back-calculated from 2026-04-17 ride, refine via Chung
DEFAULT_CRR = 0.006  # Schwalbe Marathon GreenGuard at ~70 psi (from manufacturer data)
BIKE_MASS   = 13.6   # kg (30 lbs)

# ── Silver Strand geolocation anchors ────────────────────────────────────────
# ±0.005° ≈ 500m boxes — wide enough for GPS scatter and road position variation
STRAND_ENTRY_LAT =  32.678
STRAND_ENTRY_LON = -117.175
STRAND_EXIT_LAT  =  32.593
STRAND_EXIT_LON  = -117.122
STRAND_BOX_DEG   =  0.005

# Valid heading range for Strand southbound travel
STRAND_HEADING_MIN = 115
STRAND_HEADING_MAX = 175

# Flag stops within Strand section longer than this
STOP_FLAG_SECS = 30

# Active riding detection — either cadence or power above threshold
# Handles sensor dropouts where one signal is briefly missing
ACTIVE_CAD_MIN = 10   # rpm
ACTIVE_PWR_MIN = 20   # W — fallback when cadence absent

# ── KNZY weather station ─────────────────────────────────────────────────────
KNZY_URL    = "https://forecast.weather.gov/data/obhistory/KNZY.html"
try:
    from zoneinfo import ZoneInfo
    LA_TZ = ZoneInfo("America/Los_Angeles")
except ImportError:
    LA_TZ = None  # Python < 3.9 fallback: use fixed UTC-7
PDT_OFFSET  = timedelta(hours=-7)   # fallback when zoneinfo unavailable

# ── wind direction lookup ─────────────────────────────────────────────────────
WIND_DIR_DEG = {
    'N':0,'NNE':22.5,'NE':45,'ENE':67.5,
    'E':90,'ESE':112.5,'SE':135,'SSE':157.5,
    'S':180,'SSW':202.5,'SW':225,'WSW':247.5,
    'W':270,'WNW':292.5,'NW':315,'NNW':337.5,
    'CALM':None,'VRB':None,'NONE':None
}


# ── utility functions ─────────────────────────────────────────────────────────
def despike(series, spd_list, dropout_thresh=5, speed_thresh=2.0, lookback=15):
    """
    Two-pass de-spike filter for power/cadence sensor dropouts.

    Pass 1 — dropout replacement: near-zero values while rider is moving
    (speed > speed_thresh) are replaced with the median of clean neighbors
    within ±lookback samples. Catches clustered multi-sample dropout runs
    that defeat pure Hampel (consecutive zeros corrupt the window median).

    Pass 2 — Hampel filter: standard sliding-window median + MAD on the
    cleaned series catches remaining isolated outliers. negative_only=True
    preserves sprint peaks.

    Returns (filtered_list, spike_mask).
    """
    n = len(series)
    cleaned = list(series)
    spikes = [False] * n

    # Pass 1: replace obvious dropouts (near-zero while moving)
    for i in range(n):
        if series[i] > dropout_thresh:
            continue
        if spd_list[i] < speed_thresh:
            continue  # genuinely stopped/slow — not a dropout
        neighbors = []
        for j in range(max(0, i - lookback), min(n, i + lookback + 1)):
            if j != i and series[j] > dropout_thresh:
                neighbors.append(series[j])
        if len(neighbors) >= 3:
            med = statistics.median(neighbors)
            if med > dropout_thresh * 4:  # neighbors substantially higher
                cleaned[i] = med
                spikes[i] = True

    # Pass 2: Hampel filter on cleaned series (negative_only)
    half = 3  # window=7: ±3 samples
    threshold = 3.0
    for i in range(n):
        lo = max(0, i - half)
        hi = min(n, i + half + 1)
        w = cleaned[lo:hi]
        med = statistics.median(w)
        mad = statistics.median([abs(x - med) for x in w])
        sigma = 1.4826 * mad  # scale MAD to approximate std dev
        if sigma == 0:
            continue
        if cleaned[i] < med - threshold * sigma:
            cleaned[i] = med
            spikes[i] = True

    return cleaned, spikes



def utc_to_local(dt_utc):
    """Convert UTC datetime to America/Los_Angeles, with fixed UTC-7 fallback."""
    if LA_TZ is not None:
        return dt_utc.astimezone(LA_TZ)
    return dt_utc + PDT_OFFSET


def sc_to_deg(sc):
    return sc * (180.0 / 2**31)


def fmt_duration(seconds):
    h = int(seconds // 3600)
    m = int((seconds % 3600) // 60)
    s = int(seconds % 60)
    if h:
        return f"{h}h {m:02d}m {s:02d}s"
    return f"{m}m {s:02d}s"


def in_box(lat, lon, center_lat, center_lon, radius=STRAND_BOX_DEG):
    return abs(lat - center_lat) <= radius and abs(lon - center_lon) <= radius


def bearing(lat1, lon1, lat2, lon2):
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    return math.degrees(math.atan2(dlon, dlat)) % 360


def headwind_mph(wind_spd, wind_dir_str, route_hdg):
    """
    Headwind component in mph. Positive = against rider.
    wind_dir_str: direction wind blows FROM (KNZY meteorological convention).
    """
    wind_from = WIND_DIR_DEG.get(wind_dir_str.upper())
    if wind_spd <= 0 or wind_from is None:
        return 0.0
    wind_to = (wind_from + 180) % 360
    angle = math.radians((wind_to - route_hdg) % 360)
    return -wind_spd * math.cos(angle)


def aero_penalty_watts(hw_mph_val, rider_spd_mph, cda=DEFAULT_CDA):
    hw_ms = hw_mph_val * 0.44704
    v_ms  = rider_spd_mph * 0.44704
    p_calm = 0.5 * RHO * cda * v_ms**2 * v_ms
    p_wind = 0.5 * RHO * cda * (v_ms + hw_ms)**2 * v_ms
    return p_wind - p_calm


def compute_zones(rest_hr, max_hr):
    hrr = max_hr - rest_hr
    return [
        ("Z1 Recovery",  rest_hr + 0.50*hrr, rest_hr + 0.60*hrr),
        ("Z2 Aerobic",   rest_hr + 0.60*hrr, rest_hr + 0.70*hrr),
        ("Z3 Tempo",     rest_hr + 0.70*hrr, rest_hr + 0.80*hrr),
        ("Z4 Threshold", rest_hr + 0.80*hrr, rest_hr + 0.90*hrr),
        ("Z5 VO2max",    rest_hr + 0.90*hrr, max_hr),
    ]


def normalized_power(power_series, window=30):
    """
    Compute normalized power per Coggan (2003).
    Standard definition: 30-second rolling mean of power^4, then mean of
    those values raised to 1/4. Applied to a per-second power series.
    Used here for section-level NP (Strand) where Garmin session NP is
    unavailable. Garmin session NP is used for full-ride display.
    """
    if len(power_series) < window:
        return statistics.mean(power_series) if power_series else 0
    rolling4 = []
    for i in range(window, len(power_series) + 1):
        w = power_series[i-window:i]
        rolling4.append((sum(p**4 for p in w) / window) ** 0.25)
    return statistics.mean(rolling4) if rolling4 else 0


# ── KNZY fetch and parse ──────────────────────────────────────────────────────
def parse_wind_str(wind_str):
    """
    Parse KNZY wind field: 'SE 12', 'Calm', 'NW 10 G 18', 'Vrbl 5' etc.
    Returns (direction_str, speed_mph).
    """
    wind_str = wind_str.strip()
    if not wind_str or wind_str.lower() == 'calm':
        return 'CALM', 0.0
    parts = wind_str.split()
    if len(parts) >= 2:
        direction = parts[0].upper()
        if direction == 'VRBL':
            direction = 'VRB'
        try:
            speed = float(parts[1])
            return direction, speed
        except ValueError:
            pass
    return 'CALM', 0.0


def get_cells(row_html):
    cells = re.findall(r'<t[dh][^>]*>(.*?)</t[dh]>', row_html, re.DOTALL | re.IGNORECASE)
    return [re.sub(r'<[^>]+>', '', c).strip() for c in cells]


def fetch_knzy_observations():
    """
    Fetch KNZY page and return list of observations:
    [{'day': int, 'hour': int, 'minute': int, 'wind_dir': str, 'wind_spd': float,
      'temp_f': str, 'raw': str}, ...]
    Returns None if fetch fails.
    """
    try:
        req = urllib.request.Request(
            KNZY_URL,
            headers={'User-Agent': 'Mozilla/5.0 (compatible; fit_analyze/1.0)'}
        )
        with urllib.request.urlopen(req, timeout=10) as resp:
            html = resp.read().decode('utf-8', errors='replace')
    except Exception as e:
        return None, str(e)

    rows = re.findall(r'<tr[^>]*>(.*?)</tr>', html, re.DOTALL | re.IGNORECASE)
    observations = []
    for row in rows:
        cells = get_cells(row)
        if len(cells) < 6:
            continue
        try:
            day = int(cells[0])
            time_parts = cells[1].split(':')
            if len(time_parts) != 2:
                continue
            hour, minute = int(time_parts[0]), int(time_parts[1])
            wind_dir, wind_spd = parse_wind_str(cells[2])
            temp_f = cells[6] if len(cells) > 6 else ''
            observations.append({
                'day': day,
                'hour': hour,
                'minute': minute,
                'wind_dir': wind_dir,
                'wind_spd': wind_spd,
                'temp_f': temp_f,
                'raw': f"{cells[0]} {cells[1]} | {' '.join(cells[2].split())} | {temp_f}°F",
            })
        except (ValueError, IndexError):
            continue

    return observations, None


def find_best_knzy_obs(observations, ride_start_utc):
    """
    Find the KNZY observation immediately before or closest to ride start.
    ride_start_utc: datetime in UTC.
    KNZY times are local (Pacific). Observations are at :52 past each hour.
    Returns the best matching observation dict, or None.
    """
    if not observations:
        return None

    # Convert ride start to local Pacific time
    ride_local = utc_to_local(ride_start_utc)

    # Build full datetimes for observations using ride date as anchor.
    # KNZY page only provides day-of-month, so we infer month/year from
    # the ride date. Handle month rollover: if obs day is much larger than
    # ride day (e.g., obs=31, ride=1), the obs is from the prior month.
    best = None
    best_delta = None

    for obs in observations:
        # Anchor to ride's month/year
        obs_year  = ride_local.year
        obs_month = ride_local.month

        if obs['day'] - ride_local.day > 15:
            # Obs day is much larger — must be prior month
            obs_month -= 1
            if obs_month < 1:
                obs_month = 12
                obs_year -= 1
        elif ride_local.day - obs['day'] > 15:
            # Ride day is much larger — obs is from next month
            obs_month += 1
            if obs_month > 12:
                obs_month = 1
                obs_year += 1

        # Clamp obs day to valid range for the resolved month
        max_day = calendar.monthrange(obs_year, obs_month)[1]
        obs_day = obs['day']
        if obs_day > max_day:
            print(f"  DEBUG: KNZY obs day {obs_day} clamped to {max_day} for {obs_year}-{obs_month:02d}")
            obs_day = max_day

        obs_dt = ride_local.replace(year=obs_year, month=obs_month,
                                    day=obs_day, hour=obs['hour'],
                                    minute=obs['minute'], second=0,
                                    microsecond=0)
        delta_minutes = (obs_dt - ride_local).total_seconds() / 60

        # Prefer the closest preceding observation (delta <= 0), up to 2 hours before.
        # Fall back to future observations (up to +30 min) only if none preceding.
        if -120 <= delta_minutes <= 30:
            if best_delta is None:
                best = obs
                best_delta = delta_minutes
            elif delta_minutes <= 0 and best_delta > 0:
                # Preceding beats any future observation
                best = obs
                best_delta = delta_minutes
            elif delta_minutes <= 0 and best_delta <= 0:
                # Both preceding — pick closest to ride start
                if delta_minutes > best_delta:
                    best = obs
                    best_delta = delta_minutes
            elif delta_minutes > 0 and best_delta > 0:
                # Both future (fallback) — pick closest
                if delta_minutes < best_delta:
                    best = obs
                    best_delta = delta_minutes

    return best


# ── main analysis ─────────────────────────────────────────────────────────────
def analyze(fit_path, ftp, rest_hr, max_hr,
            wind_spd_arg=None, wind_dir_arg=None, route_hdg=150):

    fitfile = FitFile(fit_path)

    # ── session summary ───────────────────────────────────────────────────────
    session = {}
    for msg in fitfile.get_messages('session'):
        for f in msg:
            if f.value is not None:
                session[f.name] = f.value

    ride_start_utc = session.get('start_time')
    if ride_start_utc and ride_start_utc.tzinfo is None:
        # FIT datetimes are UTC
        ride_start_utc = ride_start_utc.replace(tzinfo=timezone.utc)

    # ── per-second records ────────────────────────────────────────────────────
    records = []
    for msg in fitfile.get_messages('record'):
        data = {f.name: f.value for f in msg if f.value is not None}
        if 'timestamp' not in data:
            continue
        records.append({
            'ts':  data['timestamp'],
            'hr':  data.get('heart_rate', 0),
            'cad': data.get('cadence', 0),
            'pwr': data.get('power', 0),
            'spd': data.get('speed', 0),
            'dist':data.get('distance', 0),
            'lat': sc_to_deg(data['position_lat'])  if 'position_lat'  in data else None,
            'lon': sc_to_deg(data['position_long']) if 'position_long' in data else None,
            'alt': data.get('enhanced_altitude', data.get('altitude', 0)),
        })

    if not records:
        print("ERROR: No records found in FIT file.")
        sys.exit(1)

    total_secs = len(records)

    # ── full-ride active subset ───────────────────────────────────────────────
    active      = [r for r in records if r['cad'] > ACTIVE_CAD_MIN or r['pwr'] > ACTIVE_PWR_MIN]
    hr_all      = [r['hr']  for r in records if r['hr']  > 0]
    hr_active   = [r['hr']  for r in active  if r['hr']  > 0]
    pwr_active  = [r['pwr'] for r in active  if r['pwr'] > 0]
    cad_active  = [r['cad'] for r in active  if r['cad'] > 0]
    spd_active  = [r['spd'] for r in active  if r['spd'] > 0]
    active_secs = len(active)
    coast_secs  = total_secs - active_secs

    avg_hr_active  = statistics.mean(hr_active)  if hr_active  else 0
    avg_pwr_active = statistics.mean(pwr_active) if pwr_active else 0
    avg_spd_ms     = statistics.mean(spd_active) if spd_active else 0
    avg_spd_mph    = avg_spd_ms * 2.237
    # EI uses Garmin session NP for full ride — consistent with displayed NP
    # Custom normalized_power() used only for Strand section where session NP unavailable
    w_per_kg       = avg_pwr_active / RIDER_MASS

    duration       = session.get('total_timer_time', total_secs)
    distance_m     = session.get('total_distance', 0)
    distance_mi    = distance_m / 1609.34
    total_work_j   = session.get('total_work', 0)
    total_work_kj  = total_work_j / 1000
    norm_power_ses = session.get('normalized_power', 0)
    # EI uses Garmin session NP for full ride — consistent with displayed NP
    # Custom normalized_power() used only for Strand section where session NP unavailable
    ei_full        = norm_power_ses / avg_hr_active if avg_hr_active and norm_power_ses else 0
    garmin_kcal    = session.get('total_calories', 0)
    total_ascent_m = session.get('total_ascent', 0)
    tss            = session.get('training_stress_score', None)
    stored_ftp     = session.get('threshold_power', None)
    intensity_fac  = norm_power_ses / ftp if ftp else None

    # ── KNZY wind fetch ───────────────────────────────────────────────────────
    knzy_obs       = None
    knzy_err       = None
    knzy_auto      = False
    wind_spd       = wind_spd_arg if wind_spd_arg is not None else 0.0
    wind_dir       = wind_dir_arg if wind_dir_arg is not None else 'CALM'

    # Auto-fetch if no wind args provided
    if wind_spd_arg is None and ride_start_utc is not None:
        observations, knzy_err = fetch_knzy_observations()
        if observations:
            knzy_obs = find_best_knzy_obs(observations, ride_start_utc)
            if knzy_obs:
                wind_spd   = knzy_obs['wind_spd']
                wind_dir   = knzy_obs['wind_dir']
                knzy_auto  = True

    has_wind = wind_spd > 0 and wind_dir.upper() not in ('CALM', 'VRB', 'NONE')
    hw_val   = headwind_mph(wind_spd, wind_dir, route_hdg) if has_wind else 0.0
    aero_pen = aero_penalty_watts(hw_val, avg_spd_mph) if has_wind else 0.0

    if has_wind:
        if hw_val > 1:
            wind_desc   = f"{hw_val:.1f} mph headwind component"
            wind_effect = "suppressing speed / elevating HR for given power"
        elif hw_val < -1:
            wind_desc   = f"{abs(hw_val):.1f} mph tailwind component"
            wind_effect = "boosting speed / reducing HR for given power"
        else:
            wind_desc   = "near-crosswind — minimal forward effect"
            wind_effect = "minor impact on speed or HR"

    # ── Silver Strand section — geolocation detection ─────────────────────────
    strand_start_idx = None
    for i, r in enumerate(records):
        if r['lat'] is None:
            continue
        if in_box(r['lat'], r['lon'], STRAND_ENTRY_LAT, STRAND_ENTRY_LON):
            for j in range(i+1, min(i+20, len(records))):
                if records[j]['lat'] is not None:
                    hdg = bearing(r['lat'], r['lon'],
                                  records[j]['lat'], records[j]['lon'])
                    if STRAND_HEADING_MIN <= hdg <= STRAND_HEADING_MAX:
                        strand_start_idx = i
                    break
            if strand_start_idx is not None:
                break

    strand_end_idx = None
    for i in range(len(records)-1, -1, -1):
        r = records[i]
        if r['lat'] is None:
            continue
        if in_box(r['lat'], r['lon'], STRAND_EXIT_LAT, STRAND_EXIT_LON):
            strand_end_idx = i
            break

    strand_found = (strand_start_idx is not None and
                    strand_end_idx   is not None and
                    strand_end_idx   > strand_start_idx)

    if strand_found:
        sr         = records[strand_start_idx:strand_end_idx+1]
        sr_active  = [r for r in sr if r['cad'] > ACTIVE_CAD_MIN or r['pwr'] > ACTIVE_PWR_MIN]
        sr_pwr     = [r['pwr'] for r in sr_active if r['pwr'] > 0]
        sr_hr      = [r['hr']  for r in sr_active if r['hr']  > 0]
        sr_cad     = [r['cad'] for r in sr_active if r['cad'] > 0]
        sr_spd     = [r['spd'] for r in sr_active if r['spd'] > 0]

        sr_avg_pwr     = statistics.mean(sr_pwr)   if sr_pwr  else 0
        sr_avg_hr      = statistics.mean(sr_hr)    if sr_hr   else 0
        sr_avg_cad_med = statistics.median(sr_cad) if sr_cad  else 0
        sr_avg_spd_ms  = statistics.mean(sr_spd)   if sr_spd  else 0
        sr_avg_spd_mph = sr_avg_spd_ms * 2.237
        sr_np          = normalized_power([r['pwr'] for r in sr_active])
        sr_ei          = sr_np / sr_avg_hr if sr_avg_hr else 0
        sr_active_pct  = len(sr_active) / len(sr) * 100 if sr else 0
        sr_start_mi    = sr[0]['dist'] / 1609.34
        sr_end_mi      = sr[-1]['dist'] / 1609.34
        sr_dur_mins    = len(sr) / 60

        strand_stops = []
        in_stop = False
        stop_start = 0
        for si, r in enumerate(sr):
            if r['cad'] <= ACTIVE_CAD_MIN and r['pwr'] <= ACTIVE_PWR_MIN and not in_stop:
                in_stop = True
                stop_start = si
            elif (r['cad'] > ACTIVE_CAD_MIN or r['pwr'] > ACTIVE_PWR_MIN) and in_stop:
                dur = si - stop_start
                if dur > STOP_FLAG_SECS:
                    strand_stops.append((sr[stop_start]['dist']/1609.34, dur))
                in_stop = False

        chung_clean = "CLEAN" if not strand_stops else "STOPS-CHECK"
    else:
        sr_avg_pwr = sr_avg_hr = sr_avg_cad_med = 0
        sr_avg_spd_ms = sr_avg_spd_mph = sr_np = sr_ei = sr_active_pct = 0
        sr_start_mi = sr_end_mi = sr_dur_mins = 0
        strand_stops = []
        chung_clean  = "NO-STRAND"

    # ── zones ─────────────────────────────────────────────────────────────────
    zones = compute_zones(rest_hr, max_hr)
    zone_counts = []
    for name, lo, hi in zones:
        count = sum(1 for h in hr_all if lo <= h < hi)
        zone_counts.append((name, int(lo), int(hi), count))
    above_max = sum(1 for h in hr_all if h >= max_hr)

    cad_buckets = [(0,60),(60,70),(70,80),(80,90),(90,100),(100,120),(120,999)]
    pwr_buckets = [(0,50),(50,100),(100,150),(150,200),(200,250),(250,999)]

    # ── output ────────────────────────────────────────────────────────────────
    out = []
    def p(s=""): out.append(s)

    p("=" * 62)
    p(f"RIDE SUMMARY  {session.get('start_time','')}")
    p("=" * 62)

    p("\n── OVERVIEW ────────────────────────────────────────────────")
    p(f"  Duration (active):   {fmt_duration(duration)}")
    p(f"  Distance:            {distance_mi:.2f} miles  ({distance_m/1000:.2f} km)")
    p(f"  Avg speed:           {avg_spd_mph:.1f} mph")
    p(f"  Total ascent:        {total_ascent_m} m  ({total_ascent_m*3.281:.0f} ft)")
    if tss:
        p(f"  TSS:                 {tss:.1f}")

    p("\n── TIME IN MOTION ──────────────────────────────────────────")
    p(f"  Pedaling:            {fmt_duration(active_secs)}  ({active_secs/total_secs*100:.1f}%)")
    p(f"  Coasting/stopped:    {fmt_duration(coast_secs)}  ({coast_secs/total_secs*100:.1f}%)")

    p("\n── HEART RATE ──────────────────────────────────────────────")
    if hr_all:
        p(f"  Avg (whole ride):    {statistics.mean(hr_all):.1f} bpm")
    p(f"  Avg (pedaling only): {avg_hr_active:.1f} bpm")
    p(f"  Max recorded:        {max(hr_all) if hr_all else 0} bpm  (profile max: {max_hr})")
    p(f"\n  Zone distribution (whole ride, {max_hr} max / {rest_hr} rest):")
    for name, lo, hi, count in zone_counts:
        pct = count/total_secs*100
        bar = "█" * int(pct/2)
        p(f"    {name:20s} [{lo:3d}-{hi:3d}] {pct:5.1f}%  {fmt_duration(count):>10s}  {bar}")
    if above_max:
        p(f"    {'Above max':20s} [>{max_hr:3d}] {above_max/total_secs*100:5.1f}%  {fmt_duration(above_max):>10s}")

    p("\n── CADENCE (pedaling only, excl. coasting zeros) ───────────")
    if cad_active:
        p(f"  Mean:    {statistics.mean(cad_active):.1f} rpm")
        p(f"  Median:  {statistics.median(cad_active):.0f} rpm")
        p(f"  Max:     {max(cad_active)} rpm")
    p(f"\n  Distribution:")
    for lo, hi in cad_buckets:
        count = sum(1 for c in cad_active if lo <= c < hi)
        pct = count/active_secs*100 if active_secs else 0
        bar = "█" * int(pct/2)
        p(f"    {lo}-{hi if hi<999 else '+'} rpm:  {pct:5.1f}%  {fmt_duration(count):>10s}  {bar}")

    p("\n── POWER ───────────────────────────────────────────────────")
    p(f"  Avg (whole ride):    {statistics.mean([r['pwr'] for r in records]):.0f} W  (incl. coasting)")
    p(f"  Avg (pedaling only): {avg_pwr_active:.0f} W")
    p(f"  Normalized power:    {norm_power_ses} W  (session)")
    p(f"  Max:                 {max(r['pwr'] for r in records)} W")
    p(f"  FTP (device):        {stored_ftp} W  |  FTP (analysis): {ftp} W")
    if intensity_fac:
        p(f"  Intensity Factor:    {intensity_fac:.3f}")
    p(f"  W/kg (active avg):   {w_per_kg:.2f}")
    p(f"\n  Distribution (pedaling only):")
    for lo, hi in pwr_buckets:
        count = sum(1 for pw in pwr_active if lo <= pw < hi)
        pct = count/active_secs*100 if active_secs else 0
        bar = "█" * int(pct/2)
        p(f"    {lo}-{hi if hi<999 else '+'} W:    {pct:5.1f}%  {fmt_duration(count):>10s}  {bar}")

    p("\n── ENERGY & CALORIES ───────────────────────────────────────")
    p(f"  Mechanical work:     {total_work_kj:.1f} kJ")
    p(f"  Power-based kcal:    ~{total_work_kj:.0f} kcal  (kJ≈kcal approximation, see README)")
    p(f"  Garmin estimate:     {garmin_kcal} kcal  (HR-based, incl. BMR)")

    p("\n── SILVER STRAND SECTION (geolocation-anchored) ────────────")
    if strand_found:
        p(f"  Detected:            mile {sr_start_mi:.1f} to {sr_end_mi:.1f}  ({sr_dur_mins:.0f} min)")
        p(f"  Active pedaling:     {sr_active_pct:.1f}% of section")
        p(f"  Avg power (active):  {sr_avg_pwr:.0f} W")
        p(f"  Avg speed:           {sr_avg_spd_mph:.1f} mph")
        p(f"  Avg HR:              {sr_avg_hr:.1f} bpm")
        p(f"  Cadence median:      {sr_avg_cad_med:.0f} rpm")
        p(f"  Normalized power:    {sr_np:.0f} W")
        sr_rr_watts = DEFAULT_CRR * sr_avg_spd_ms * (RIDER_MASS + BIKE_MASS) * 9.81
        sr_aero_watts = 0.5 * RHO * DEFAULT_CDA * sr_avg_spd_ms**3
        p(f"  Efficiency Index:    {sr_ei:.3f} W/bpm  ← primary trend metric")
        p(f"  Est. rolling loss:   {sr_rr_watts:.0f} W  (Crr={DEFAULT_CRR})")
        p(f"  Est. aero drag:      {sr_aero_watts:.0f} W  (CdA={DEFAULT_CDA})")
        if strand_stops:
            p(f"  !! Stops >30s:       {len(strand_stops)} detected")
            for sdist, sdur in strand_stops:
                p(f"     Mile {sdist:.2f}: {sdur}s stop — may affect Chung data")
        else:
            p(f"  Stops >30s:          none — clean for Chung regression")
    else:
        p(f"  Not detected — GPS entry/exit boxes not matched.")
        p(f"  Entry box: {STRAND_ENTRY_LAT}N ±{STRAND_BOX_DEG}°, {STRAND_ENTRY_LON}W ±{STRAND_BOX_DEG}°")
        p(f"  Exit box:  {STRAND_EXIT_LAT}N ±{STRAND_BOX_DEG}°, {STRAND_EXIT_LON}W ±{STRAND_BOX_DEG}°")

    p("\n── WIND CONDITIONS (KNZY / NAS North Island) ───────────────")
    if knzy_auto and knzy_obs:
        ride_local = utc_to_local(ride_start_utc).strftime('%H:%M %Z') if ride_start_utc else '?'
        p(f"  Auto-fetched from KNZY")
        p(f"  Ride start:          {ride_local}")
        p(f"  Observation used:    Day {knzy_obs['day']} {knzy_obs['hour']:02d}:{knzy_obs['minute']:02d} PDT")
        p(f"  Raw observation:     {knzy_obs['raw']}")
        if wind_spd > 0:
            p(f"  Parsed:              {wind_spd:.0f} mph from {wind_dir}")
        else:
            p(f"  Parsed:              Calm")
        p(f"  To override: --wind-speed {wind_spd:.0f} --wind-dir {wind_dir}")
    elif knzy_err and wind_spd_arg is None:
        p(f"  KNZY fetch failed: {knzy_err}")
        p(f"  Provide manually: --wind-speed N --wind-dir DIR")
        p(f"  Source: forecast.weather.gov/data/obhistory/KNZY.html")
    elif wind_spd_arg is not None:
        p(f"  Manually provided:   {wind_spd:.0f} mph from {wind_dir.upper()}")

    if has_wind:
        p(f"  Route heading:       {route_hdg}°")
        p(f"  Effect on route:     {wind_desc}")
        p(f"  Implication:         {wind_effect}")
        p(f"  Est. aero penalty:   {aero_pen:+.0f} W vs calm at {avg_spd_mph:.1f} mph")
        p(f"  Note: Power meter reads actual mechanical output.")
        p(f"        EI naturally captures wind — use EI for trend comparison.")
    else:
        p(f"  Conditions:          Calm or not recorded")

    p("\n── EFFICIENCY INDEX (primary fitness trend metric) ─────────")
    p(f"  Full ride EI:        {ei_full:.3f} W/bpm")
    if strand_found:
        p(f"  Strand section EI:   {sr_ei:.3f} W/bpm  ← use this for trend")
        p(f"  (Strand EI strips sheltered sections — more comparable ride to ride)")
    if has_wind and hw_val > 2:
        p(f"  Note: {hw_val:.1f} mph headwind elevated HR — calm days will show higher EI.")
    if intensity_fac:
        p(f"  NP/FTP:              {intensity_fac:.3f}")
    p(f"  W/kg (active):       {w_per_kg:.2f}")

    p("\n── CHUNG METHOD DATA (accumulate across rides) ──────────────")
    p("  # Paste DATA lines into a ride log for CdA/Crr regression.")
    p("  # Uses Strand section values — open fetch, consistent heading.")
    p("  # Fields: date | pwr | spd | NP | wind | hdg | hw | rho | Crr | EI | status")
    chung_date = str(session.get('start_time', '')).split(' ')[0]
    chung_pwr  = sr_avg_pwr    if strand_found else avg_pwr_active
    chung_spd  = sr_avg_spd_ms if strand_found else avg_spd_ms
    chung_np   = sr_np         if strand_found else norm_power_ses
    chung_ei   = sr_ei         if strand_found else ei_full
    p(f"  DATA: {chung_date}"
      f" | pwr={chung_pwr:.1f}W"
      f" | spd={chung_spd:.3f}m/s"
      f" | NP={chung_np:.0f}W"
      f" | wind={wind_spd:.0f}mph_{wind_dir.upper()}"
      f" | hdg={route_hdg}deg"
      f" | hw={hw_val:.1f}mph"
      f" | rho={RHO:.3f}"
      f" | Crr={DEFAULT_CRR:.4f}"
      f" | EI={chung_ei:.3f}"
      f" | {chung_clean}")

    p("\n" + "=" * 62)
    p("Paste output into chat for quick discussion.")
    p("Drop raw .fit file only for deep dives.")
    p("=" * 62)

    plot_data = {
        'records':          records,
        'strand_start_idx': strand_start_idx if strand_found else None,
        'strand_end_idx':   strand_end_idx   if strand_found else None,
        'wind_desc':        wind_desc if has_wind else '',
        'strand_stats': {
            'active_pct':  sr_active_pct,
            'ei':          sr_ei,
            'stops':       len(strand_stops),
        } if strand_found else None,
    }
    return "\n".join(out), plot_data



def plot_strand(records, strand_start_idx, strand_end_idx, wind_desc="", fit_path="", strand_stats=None):
    """
    Plot Silver Strand section data using gnuplot.
    Produces both a PNG file and an interactive window.
    Requires gnuplot to be installed (checked gracefully).
    """
    import shutil
    import subprocess
    import tempfile
    import os

    if not shutil.which('gnuplot'):
        print("WARNING: gnuplot not found. Install gnuplot to enable plotting.")
        print("         Ubuntu/Debian: sudo apt install gnuplot")
        print("         conda:         conda install -c conda-forge gnuplot")
        return

    sr = records[strand_start_idx:strand_end_idx+1]
    if not sr:
        print("No Strand section data to plot.")
        return

    # Determine output PNG path alongside the FIT file
    if fit_path:
        base = fit_path.replace('.fit', '_strand.png')
    else:
        base = 'strand_plot.png'

    # De-spike power and cadence (sensor dropouts)
    raw_pwr = [r['pwr'] for r in sr]
    raw_cad = [r['cad'] for r in sr]
    spd_list = [r['spd'] for r in sr]
    filt_pwr, pwr_spikes = despike(raw_pwr, spd_list)
    filt_cad, cad_spikes = despike(raw_cad, spd_list)
    n_pwr_spikes = sum(pwr_spikes)
    n_cad_spikes = sum(cad_spikes)
    print(f"Despiked: {n_pwr_spikes} power dropouts, {n_cad_spikes} cadence dropouts corrected")

    # Active mask uses filtered data so dropouts don't create false inactive gaps
    active_mask = [filt_cad[i] > ACTIVE_CAD_MIN or filt_pwr[i] > ACTIVE_PWR_MIN
                   for i in range(len(sr))]

    # Compute averages from filtered active values for reference lines
    pwr_vals  = [filt_pwr[i] for i in range(len(sr)) if active_mask[i] and filt_pwr[i] > 0]
    hr_vals   = [r['hr']  for r, a in zip(sr, active_mask) if a and r['hr']  > 0]
    spd_vals  = [r['spd'] * 2.237 for r, a in zip(sr, active_mask) if a and r['spd'] > 0]
    avg_pwr   = statistics.mean(pwr_vals)  if pwr_vals  else 0
    avg_hr    = statistics.mean(hr_vals)   if hr_vals   else 0
    avg_spd   = statistics.mean(spd_vals)  if spd_vals  else 0

    # Write data file: time_min pwr hr cad spd active filt_pwr filt_cad pwr_spike cad_spike
    tmpdir = tempfile.mkdtemp()
    data_file = os.path.join(tmpdir, 'strand.dat')
    with open(data_file, 'w') as f:
        f.write("# time_min  pwr  hr  cad  spd_mph  active  filt_pwr  filt_cad  pwr_spike  cad_spike\n")
        for i, r in enumerate(sr):
            t   = i / 60.0
            pwr = r['pwr']
            hr  = r['hr']  if r['hr']  > 0 else 0
            cad = r['cad'] if r['cad'] > 0 else 0
            spd = r['spd'] * 2.237
            act = 1 if active_mask[i] else 0
            fp  = filt_pwr[i]
            fc  = filt_cad[i]
            ps  = 1 if pwr_spikes[i] else 0
            cs  = 1 if cad_spikes[i] else 0
            f.write(f"{t:.4f}  {pwr}  {hr}  {cad}  {spd:.2f}  {act}  {fp:.0f}  {fc:.0f}  {ps}  {cs}\n")

    dist_start = sr[0]['dist']  / 1609.34
    dist_end   = sr[-1]['dist'] / 1609.34
    dur_mins   = len(sr) / 60.0

    # Build entry timestamp from Strand start record
    entry_ts = sr[0].get('ts')
    if entry_ts:
        entry_local = utc_to_local(entry_ts.replace(tzinfo=timezone.utc) if entry_ts.tzinfo is None else entry_ts)
        date_str = entry_local.strftime('%Y-%m-%d %H:%M %Z')
    else:
        date_str = ''

    # Title line 1: date/time and route info
    title_line1 = (f"{date_str}  —  "
                   f"Silver Strand  mile {dist_start:.1f}-{dist_end:.1f}"
                   f"  {dur_mins:.0f} min"
                   + (f"  |  {wind_desc}" if wind_desc else ""))

    # Title line 2: key stats
    stats_parts = []
    if strand_stats:
        stats_parts.append(f"Active {strand_stats['active_pct']:.0f}%")
        stats_parts.append(f"EI {strand_stats['ei']:.3f}")
        if strand_stats['stops'] > 0:
            stats_parts.append(f"{strand_stats['stops']} stops >30s")
        else:
            stats_parts.append("no stops >30s")
    stats_parts.append(f"plot filtered ({n_pwr_spikes}pwr {n_cad_spikes}cad)")
    title_line2 = "  |  ".join(stats_parts)

    # HR zone boundaries for reference lines
    zones = compute_zones(DEFAULT_REST_HR, DEFAULT_MAX_HR)
    z_boundaries = [int(lo) for _, lo, _ in zones] + [int(zones[-1][2])]

    # gnuplot script — renders to PNG then reopens as window
    gp_script = f"""
set terminal pngcairo size 1200,900 font "Sans,9"
set output "{base}"
set multiplot layout 4,1 title "{title_line1}\\n{title_line2}" font "Sans,10"
set xrange [0:{dur_mins:.2f}]
set grid ytics lc rgb "#cccccc" lw 0.5
set key top right font "Sans,8"
set lmargin 10
set rmargin 4

# Panel 1: Power (filtered, with dropout markers at corrected value)
set ylabel "Power (W)"
set yrange [0:*]
plot "{data_file}" using 1:($6==1 ? $7 : 1/0) with lines lc rgb "#2196F3" lw 1 title "Power", \
     {avg_pwr:.0f} with lines lc rgb "#2196F3" lw 1 dt 2 title "avg {avg_pwr:.0f}W", \
     "{data_file}" using 1:($6==0 ? $7 : 1/0) with lines lc rgb "#cccccc" lw 1 notitle

# Panel 2: Heart rate with zone lines
set ylabel "HR (bpm)"
set yrange [{z_boundaries[0]-5}:{z_boundaries[-1]+5}]
"""
    for i, bnd in enumerate(z_boundaries):
        if i < 5:
            arrow = ('set arrow from 0,' + str(bnd) + ' to ' +
                     f'{dur_mins:.2f},' + str(bnd) +
                     ' nohead lc rgb "#999999" lw 0.5 dt 2 front' + '\n')
            gp_script += arrow
    gp_script += f"""
plot "{data_file}" using 1:($3 > 0 ? $3 : 1/0) with lines lc rgb "#E53935" lw 1 title "HR", \
     {avg_hr:.0f} with lines lc rgb "#E53935" lw 1 dt 2 title "avg {avg_hr:.0f} bpm"
unset arrow

# Panel 3: Cadence (filtered, with dropout markers at corrected value)
set ylabel "Cadence (rpm)"
set yrange [0:*]
plot "{data_file}" using 1:($6==1 ? $8 : 1/0) with lines lc rgb "#43A047" lw 1 title "Cadence", \
     90 with lines lc rgb "#888888" lw 1 dt 3 title "90 rpm target", \
     "{data_file}" using 1:($6==0 ? $8 : 1/0) with lines lc rgb "#cccccc" lw 1 notitle

# Panel 4: Speed
set ylabel "Speed (mph)"
set xlabel "Time in section (minutes)"
set yrange [0:*]
plot "{data_file}" using 1:($6==1 ? $5 : 1/0) with lines lc rgb "#FB8C00" lw 1 title "Speed", \
     {avg_spd:.1f} with lines lc rgb "#FB8C00" lw 1 dt 2 title "avg {avg_spd:.1f} mph", \
     "{data_file}" using 1:($6==0 ? $5 : 1/0) with lines lc rgb "#cccccc" lw 1 notitle

unset multiplot
"""

    gp_file = os.path.join(tmpdir, 'strand.gp')
    with open(gp_file, 'w') as f:
        f.write(gp_script)

    # Render PNG
    try:
        result = subprocess.run(['gnuplot', gp_file],
                                capture_output=True, text=True)
        if result.returncode != 0:
            print(f"gnuplot error: {result.stderr}")
            return
        print(f"Plot saved to: {base}")
    except Exception as e:
        print(f"gnuplot failed: {e}")
        return

    # Open interactive window using wxt terminal
    png_header = 'set terminal pngcairo size 1200,900 font "Sans,9"' + '\n' + f'set output "{base}"'
    wxt_header = 'set terminal wxt size 1200,900 font "Sans,9" persist'
    wxt_script = gp_script.replace(png_header, wxt_header)
    wxt_file = os.path.join(tmpdir, 'strand_wxt.gp')
    with open(wxt_file, 'w') as f:
        f.write(wxt_script)

    try:
        subprocess.Popen(['gnuplot', wxt_file])
    except Exception as e:
        print(f"Interactive window failed: {e}")
        print(f"PNG is available at: {base}")


def main():
    parser = argparse.ArgumentParser(
        description="FIT file preprocessor — Silver Strand geolocation-anchored"
    )
    parser.add_argument("fit_file")
    parser.add_argument("--ftp",           type=int,   default=DEFAULT_FTP)
    parser.add_argument("--rest-hr",       type=int,   default=DEFAULT_REST_HR)
    parser.add_argument("--max-hr",        type=int,   default=DEFAULT_MAX_HR)
    parser.add_argument("--wind-speed",    type=float, default=None,
                        help="KNZY wind speed mph. If omitted, auto-fetched from KNZY.")
    parser.add_argument("--wind-dir",      type=str,   default=None,
                        help="KNZY wind direction e.g. SE, E, NW (wind FROM).")
    parser.add_argument("--plot", action="store_true",
                        help="Plot Strand section data via gnuplot (PNG + interactive window).")
    parser.add_argument("--route-heading", type=int,   default=150,
                        help="Travel direction in degrees (Strand southbound ~150).")
    args = parser.parse_args()

    result, plot_data = analyze(
        args.fit_file, args.ftp, args.rest_hr, args.max_hr,
        args.wind_speed, args.wind_dir, args.route_heading
    )
    print(result)

    if args.plot and plot_data:
        if plot_data['strand_start_idx'] is None or plot_data['strand_end_idx'] is None:
            print("WARNING: --plot requires Strand section detection. No plot generated.")
        else:
            wind_desc = plot_data.get('wind_desc', '')
            plot_strand(plot_data['records'],
                        plot_data['strand_start_idx'],
                        plot_data['strand_end_idx'],
                        wind_desc=wind_desc,
                        fit_path=args.fit_file,
                        strand_stats=plot_data.get('strand_stats'))

    out_path = args.fit_file.replace('.fit', '_summary.txt')
    try:
        with open(out_path, 'w') as f:
            f.write(result)
        print(f"\n[Summary saved to: {out_path}]")
    except Exception:
        pass



if __name__ == "__main__":
    main()
