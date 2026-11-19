"""Forecast GTA VI (PS5) trophy earn rates on Exophase, one week after launch.

Entry for make.fun competition gta6-trophy-rarity-forecast-01.

Method (per trophy, all my own design):
  * Work on the logit scale, logit(p) = ln(p / (1 - p)), with rates clamped to
    0.01..99.99 percent before taking the logit, so 0% and 100% stay finite.
  * Fit a straight line of logit(rate) against ln(1 + hours since launch,
    2026-11-19T00:00:00Z) by weighted least squares over every read in the
    case, weighting recent reads more (exponential decay, 12-hour half-life).
    Exophase refreshes dilute rates as the tracked-player base grows, and
    that dilution shows up as a roughly linear drift on this scale, so the
    line extrapolates it to the target instead of freezing the last read.
  * Spread: the fit's weighted residual scatter, combined in quadrature with
    a random-walk term that widens with the horizon (0.5 logit units per
    square-root day) and a refresh-jump allowance (0.5 logit units) for the
    unobserved Exophase refreshes between the case cut and the target. The
    launch-night evidence justifies the jump term: the first refresh cut
    most rates roughly in half in logit terms while the history was flat.
  * Quantiles come from the standard normal at the case's levels, mapped
    back through the logistic function and clamped to 0.01..99.99.

Pure standard library, no network, no randomness: two runs on the same
input write byte-identical output.
"""

import json
import math
import os
from datetime import datetime, timezone
from pathlib import Path

INPUT_DIR = Path(os.environ.get("MAKEFUN_INPUT_DIR", "/input"))
OUTPUT_DIR = Path(os.environ.get("MAKEFUN_OUTPUT_DIR", "/output"))

LAUNCH = datetime(2026, 11, 19, 0, 0, 0, tzinfo=timezone.utc)

# Recency weighting: a read loses half its weight every HALF_LIFE_HOURS.
HALF_LIFE_HOURS = 12.0
# Random-walk spread, logit units per square-root day of horizon.
RW_SIGMA = 0.5
# Allowance for unobserved Exophase refresh jumps, logit units.
JUMP_SIGMA = 0.5
# Floor on the residual scatter so a perfectly flat history still carries
# honest uncertainty.
RESID_FLOOR = 0.05
# Safety valve on how far the extrapolated median may move from the last
# observed rate, in logit units.
MAX_LOGIT_SHIFT = 4.0

LO = 0.01
HI = 99.99


def parse_time(text):
    """Parse an ISO-8601 UTC timestamp (accepts a trailing Z)."""
    return datetime.fromisoformat(text.replace("Z", "+00:00"))


def clamp_pct(rate):
    if rate is None:
        return None
    return min(HI, max(LO, float(rate)))


def logit_of(pct):
    p = clamp_pct(pct) / 100.0
    return math.log(p / (1.0 - p))


def pct_of(logit_value):
    if logit_value >= 0:
        p = 1.0 / (1.0 + math.exp(-logit_value))
    else:
        e = math.exp(logit_value)
        p = e / (1.0 + e)
    return min(HI, max(LO, p * 100.0))


def norm_ppf(level):
    """Inverse standard normal CDF (Acklam's rational approximation)."""
    a = [-3.969683028665376e+01, 2.209460984245205e+02,
         -2.759285984469534e+02, 1.383577518672690e+02,
         -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02,
         -1.556989798598866e+02, 6.680131188771972e+01,
         -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01,
         -2.400758277161838e+00, -2.549732539343734e+00,
         4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01,
         2.445134137142996e+00, 3.754408661907416e+00]
    p_low = 0.02425
    p_high = 1.0 - p_low
    if level < p_low:
        q = math.sqrt(-2.0 * math.log(level))
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / \
               ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0)
    if level <= p_high:
        q = level - 0.5
        r = q * q
        return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / \
               (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0)
    q = math.sqrt(-2.0 * math.log(1.0 - level))
    return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / \
            ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0)


def round4(value):
    """Round to 4 decimals, half up, returning an int for whole numbers."""
    from decimal import Decimal, ROUND_HALF_UP
    rounded = Decimal(repr(value)).quantize(Decimal("0.0001"), rounding=ROUND_HALF_UP)
    if rounded == rounded.to_integral_value():
        return int(rounded)
    return float(rounded)


def weighted_line(xs, ys, ws):
    """Weighted least squares y ~ a + b*x. Returns (a, b) or None if degenerate."""
    sw = sum(ws)
    if sw <= 0:
        return None
    mx = sum(w * x for w, x in zip(ws, xs)) / sw
    my = sum(w * y for w, y in zip(ws, ys)) / sw
    sxx = sum(w * (x - mx) ** 2 for w, x in zip(ws, xs))
    if sxx <= 1e-12:
        return None
    sxy = sum(w * (x - mx) * (y - my) for w, x, y in zip(ws, xs, ys))
    b = sxy / sxx
    return (my - b * mx, b)


def forecast_trophy(reads, target_time, horizon_days, levels):
    """reads: list of (time, rate-or-None) oldest first for one trophy."""
    usable = [(t, r) for (t, r) in reads if r is not None]
    if not usable:
        centre = logit_of(LO)
        sigma = math.sqrt(RESID_FLOOR ** 2 + (RW_SIGMA * math.sqrt(horizon_days)) ** 2 + JUMP_SIGMA ** 2)
        return [round4(pct_of(centre + norm_ppf(lv) * sigma)) for lv in levels]

    last_time, last_rate = usable[-1]
    last_logit = logit_of(last_rate)

    xs = [math.log1p(max(0.0, (t - LAUNCH).total_seconds() / 3600.0)) for (t, _) in usable]
    ys = [logit_of(r) for (_, r) in usable]
    ws = [0.5 ** ((last_time - t).total_seconds() / 3600.0 / HALF_LIFE_HOURS) for (t, _) in usable]

    line = weighted_line(xs, ys, ws) if len(set(xs)) >= 2 else None
    if line is None:
        centre = last_logit
        resid_sd = RESID_FLOOR
    else:
        a, b = line
        x_target = math.log1p(max(0.0, (target_time - LAUNCH).total_seconds() / 3600.0))
        centre = a + b * x_target
        # Keep the extrapolation honest: never move further than the valve
        # allows from where the trophy was last seen.
        centre = min(last_logit + MAX_LOGIT_SHIFT, max(last_logit - MAX_LOGIT_SHIFT, centre))
        centre = min(logit_of(HI), max(logit_of(LO), centre))
        sw = sum(ws)
        resid_var = sum(w * (y - (a + b * x)) ** 2 for w, x, y in zip(ws, xs, ys)) / sw
        resid_sd = math.sqrt(max(resid_var, 0.0))

    sigma = math.sqrt(max(resid_sd, RESID_FLOOR) ** 2
                      + (RW_SIGMA * math.sqrt(horizon_days)) ** 2
                      + JUMP_SIGMA ** 2)
    quantiles = [round4(pct_of(centre + norm_ppf(lv) * sigma)) for lv in levels]
    # Monotonicity is guaranteed by construction (levels ascending, sigma > 0),
    # but enforce it explicitly so rounding can never break it.
    for i in range(1, len(quantiles)):
        if quantiles[i] < quantiles[i - 1]:
            quantiles[i] = quantiles[i - 1]
    return quantiles


def forecast_case(case):
    target_time = parse_time(case["target"])
    horizon_days = float(case["horizon_minutes"]) / 1440.0
    levels = [float(lv) for lv in case["levels"]]
    trophies = case["trophies"]
    history = case.get("history") or []

    rows = []
    for index, trophy in enumerate(trophies):
        reads = []
        for read in history:
            rates = read.get("rates") or []
            rate = rates[index] if index < len(rates) else None
            reads.append((parse_time(read["at"]), rate))
        rows.append({
            "id": trophy["id"],
            "q": forecast_trophy(reads, target_time, horizon_days, levels),
        })
    return {"trophies": rows}


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for path in sorted(INPUT_DIR.glob("*.in")):
        case = json.loads(path.read_text(encoding="utf-8"))
        forecast = forecast_case(case)
        out_path = OUTPUT_DIR / (path.stem + ".out")
        out_path.write_text(json.dumps(forecast, separators=(",", ":")) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
