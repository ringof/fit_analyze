#!/usr/bin/env python3
"""
knzy_test.py — Standalone KNZY weather fetch debugger
Usage: python3 knzy_test.py [--date 2026-04-17] [--time-utc 11:23]
"""

import urllib.request
import re
import sys
import argparse
from datetime import datetime, timezone, timedelta

KNZY_URL   = "https://forecast.weather.gov/data/obhistory/KNZY.html"
PDT_OFFSET = timedelta(hours=-7)


def get_cells(row_html):
    """Extract text from both TD and TH cells, case-insensitive."""
    cells = re.findall(r'<t[dh][^>]*>(.*?)</t[dh]>', row_html,
                       re.DOTALL | re.IGNORECASE)
    return [re.sub(r'<[^>]+>', '', c).strip() for c in cells]


def parse_wind(wind_str):
    wind_str = wind_str.strip()
    if not wind_str or wind_str.lower() in ('calm', 'vrbl'):
        return ('CALM' if 'calm' in wind_str.lower() else 'VRB'), 0.0
    parts = wind_str.split()
    if len(parts) >= 2:
        direction = parts[0].upper()
        if direction == 'VRBL':
            direction = 'VRB'
        try:
            return direction, float(parts[1])
        except ValueError:
            pass
    return 'CALM', 0.0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--date',     default=None)
    parser.add_argument('--time-utc', default=None)
    args = parser.parse_args()

    print(f"Fetching: {KNZY_URL}")
    try:
        req = urllib.request.Request(
            KNZY_URL,
            headers={'User-Agent': 'Mozilla/5.0 (compatible; knzy_test/1.0)'}
        )
        with urllib.request.urlopen(req, timeout=10) as resp:
            html = resp.read().decode('utf-8', errors='replace')
        print(f"Fetch OK — {len(html)} bytes")
    except Exception as e:
        print(f"FETCH FAILED: {e}")
        sys.exit(1)

    # ── diagnose HTML structure ───────────────────────────────────────────────
    print()
    print("=== HTML STRUCTURE DIAGNOSIS ===")

    # Count various tag types
    for tag in ['<tr', '<TR', '<td', '<TD', '<th', '<TH', '<table', '<TABLE']:
        count = html.lower().count(tag.lower())
        print(f"  {tag!r:12s}: {count}")

    # Show a sample of what's around the data table
    # Look for the word 'Calm' or a wind direction which appears in obs rows
    idx = html.find('Calm')
    if idx == -1:
        idx = html.find(' SE ')
    if idx == -1:
        idx = html.find(' NW ')
    if idx >= 0:
        print(f"\n=== HTML AROUND FIRST WIND DATA (index {idx}) ===")
        print(html[max(0,idx-300):idx+300])
    else:
        print("\nNo 'Calm' or wind direction found in HTML — showing middle section")
        mid = len(html)//2
        print(html[mid:mid+500])

    # ── try row extraction with case-insensitive pattern ─────────────────────
    print()
    rows_ci = re.findall(r'<tr[^>]*>(.*?)</tr>', html,
                         re.DOTALL | re.IGNORECASE)
    print(f"=== ROWS (case-insensitive <tr>): {len(rows_ci)} ===")

    # Show first 10 rows with their cell counts
    for i, row in enumerate(rows_ci[:10]):
        cells = get_cells(row)
        preview = ' | '.join(cells[:5]) if cells else '(empty)'
        print(f"  Row {i:2d}: {len(cells)} cells — {preview}")

    # ── parse observations ────────────────────────────────────────────────────
    print()
    observations = []
    for row in rows_ci:
        cells = get_cells(row)
        if len(cells) < 3:
            continue
        try:
            day = int(cells[0])
            if not (1 <= day <= 31):
                continue
            time_parts = cells[1].split(':')
            if len(time_parts) != 2:
                continue
            hour, minute = int(time_parts[0]), int(time_parts[1])
            wind_dir, wind_spd = parse_wind(cells[2])
            temp_f = cells[6] if len(cells) > 6 else '?'
            observations.append({
                'day': day, 'hour': hour, 'minute': minute,
                'wind_dir': wind_dir, 'wind_spd': wind_spd,
                'temp_f': temp_f,
            })
        except (ValueError, IndexError):
            continue

    print(f"=== PARSED OBSERVATIONS: {len(observations)} ===")
    if observations:
        print(f"{'Day':>4} {'Time':>6} {'Dir':>6} {'Spd':>6} {'Temp':>5}")
        print("-" * 35)
        for obs in observations:
            print(f"  {obs['day']:2d}  {obs['hour']:02d}:{obs['minute']:02d}"
                  f"  {obs['wind_dir']:>6}  {obs['wind_spd']:>5.0f}mph"
                  f"  {obs['temp_f']:>4}°F")
    else:
        print("Still no observations — need to inspect HTML structure further")
        # Show all rows with >3 cells
        print("\nRows with >3 cells:")
        for i, row in enumerate(rows_ci):
            cells = get_cells(row)
            if len(cells) > 3:
                print(f"  Row {i}: {cells[:6]}")

    # ── time matching ─────────────────────────────────────────────────────────
    if observations and args.date and args.time_utc:
        print()
        print("=== TIME MATCHING ===")
        ride_utc = datetime.strptime(
            f"{args.date} {args.time_utc}:00", "%Y-%m-%d %H:%M:%S"
        ).replace(tzinfo=timezone.utc)
        ride_pdt = ride_utc + PDT_OFFSET
        print(f"Ride UTC: {ride_utc}  →  PDT: {ride_pdt.strftime('%H:%M day %d')}")

        best = None
        best_delta = None
        for obs in observations:
            day_diff = obs['day'] - ride_pdt.day
            delta = (day_diff * 1440 +
                     (obs['hour'] - ride_pdt.hour) * 60 +
                     (obs['minute'] - ride_pdt.minute))
            in_window = -120 <= delta <= 30
            marker = ""
            if in_window and (best_delta is None or abs(delta) < abs(best_delta)):
                best = obs
                best_delta = delta
                marker = " ← BEST"
            if abs(delta) < 180:  # only show nearby
                print(f"  Day {obs['day']:2d} {obs['hour']:02d}:{obs['minute']:02d}"
                      f" | delta={delta:+4d}m | {'IN' if in_window else 'out'}"
                      f" | {obs['wind_dir']} {obs['wind_spd']:.0f}mph{marker}")

        print()
        if best:
            print(f"RESULT: {best['wind_dir']} {best['wind_spd']:.0f} mph")
        else:
            print("RESULT: No match found")


if __name__ == "__main__":
    main()
