# Example data

| File | Content | Source and licence |
|---|---|---|
| `prop99_long.csv` | Cigarette sales per capita, 39 US states, 1970-2000 (California treated in 1989) | Abadie, Diamond and Hainmueller (2010), JASA |
| `prop99_covariates.csv` | The same panel plus `lnincome`, `beer`, `age15to24`, `retprice` | `smoking` data of the R package `tidysynth` (the original replication data); `cigsale` agrees with `prop99_long.csv` |
| `west_germany.csv` | GDP per capita, inflation, trade, schooling and industry share for 17 OECD countries, 1960-2003 (West Germany treated in 1990) | `scpi_germany` in the R package `scpi` (Cattaneo, Feng, Palomba and Titiunik), GPL-3; data of Abadie, Diamond and Hainmueller (2015), AJPS |
| `mpdta_minwage.csv` | Log teen employment in 500 US counties, 2003-2007; minimum wage increases first applied in 2004, 2006 or 2007, 309 never-treated counties | `mpdta` in the R package `did` (Callaway and Sant'Anna), GPL-3; derived from county data, see Callaway and Sant'Anna (2021), Journal of Econometrics |

Notes

* `west_germany.csv`: `infrate`, `trade`, `schooling` and `industry` have gaps in some years, so they work as covariates averaged over the pre-treatment period (Penalized, Demeaned, Augmented SC). Interactive fixed effects and SDID need covariates observed in every period and cannot use them.
* `mpdta_minwage.csv`: `first_treat` is blank for never-treated counties and gives the first treated period; `treated` is the 0/1 indicator that switches on at that period. With 191 treated counties, use *time cohorts* in Staggered Augmented SC. Counties first treated in 2004 have a single pre-treatment year, so switch off unit fixed effects.
