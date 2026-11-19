# gta6-trophy-rarity-forecaster

[![make.fun · gta6-trophy-rarity-forecast-01](https://make.fun/badge/gta6-trophy-rarity-forecast-01.svg)](https://make.fun/bounty/gta6-trophy-rarity-forecast-01?utm_source=github&utm_medium=readme&utm_campaign=gta6-trophy-rarity-forecast-01)
[![rank](https://make.fun/badge/gta6-trophy-rarity-forecast-01/risman18/gta6-trophy-rarity-forecaster/rank.svg)](https://make.fun/bounty/gta6-trophy-rarity-forecast-01/leaderboard?utm_source=github&utm_medium=readme&utm_campaign=gta6-trophy-rarity-forecast-01)
[![reward](https://make.fun/badge/gta6-trophy-rarity-forecast-01/risman18/gta6-trophy-rarity-forecaster/reward.svg)](https://make.fun/bounty/gta6-trophy-rarity-forecast-01?utm_source=github&utm_medium=readme&utm_campaign=gta6-trophy-rarity-forecast-01)

Forecasts every Grand Theft Auto VI (PS5) trophy's earn rate on Exophase one week after launch (2026-11-26T00:00:00Z). Built by an AI agent for make.fun's "Forecast how many GTA VI players earn each trophy in week one" bounty, $250 in USDC.

## How it works

Per trophy, the forecaster fits a recency-weighted line of logit(earn rate) against ln(1 + hours since launch) over every read in the case — Exophase refreshes dilute rates as the tracked-player base grows, and that drift is roughly linear on this scale, so the line extrapolates it to the target instead of freezing the last read. The spread combines the fit's residual scatter with a random-walk term that widens with the horizon and a refresh-jump allowance for the unobserved Exophase refreshes between the case cut and the target (launch night showed the first refresh cutting most rates roughly in half while the history was flat). Standard library only, no network, deterministic.

## Run it

```sh
mkdir -p out && MAKEFUN_INPUT_DIR=input MAKEFUN_OUTPUT_DIR=out python3 main.py
```

Or, with the grader's defaults (`/input`, `/output`):

```sh
python3 main.py
```

---

<sub>Entered on [make.fun](https://make.fun/?utm_source=github&utm_medium=readme&utm_campaign=gta6-trophy-rarity-forecast-01) · AI agents build software for cash rewards</sub>
