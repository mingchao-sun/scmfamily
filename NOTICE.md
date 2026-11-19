# Third-party code and data

scmfamily is distributed under the GNU General Public License, version 3 (see `LICENSE`).
Parts of it are ports of other open-source software and are used under their own licences.
The notices below must be kept with every copy of the module.

## synthdid (R/sdid_engine.R)

`R/sdid_engine.R` is a port of the R package *synthdid*
(<https://github.com/synth-inference/synthdid>) by Dmitry Arkhangelsky, Susan Athey,
David A. Hirshberg, Guido W. Imbens and Stefan Wager. The package is licensed
"GPL (>= 2) | BSD_3_clause + file LICENSE"; this module uses it under the GPL. The BSD 3-clause
notice of the package reads:

> Copyright 2019, Stanford University
>
> Redistribution and use in source and binary forms, with or without modification, are permitted
> provided that the following conditions are met:
>
> 1. Redistributions of source code must retain the above copyright notice, this list of conditions
> and the following disclaimer.
>
> 2. Redistributions in binary form must reproduce the above copyright notice, this list of
> conditions and the following disclaimer in the documentation and/or other materials provided with
> the distribution.
>
> 3. Neither the name of the copyright holder nor the names of its contributors may be used to
> endorse or promote products derived from this software without specific prior written permission.
>
> THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
> IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
> FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
> CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
> DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
> DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
> IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
> OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

## augsynth (R/augsynth_engine.R, R/multisynth_engine.R)

`R/augsynth_engine.R` and `R/multisynth_engine.R` are ports of the R package *augsynth*
(<https://github.com/ebenmichael/augsynth>) by Eli Ben-Michael, licensed under the MIT licence:

> MIT License
>
> Copyright (c) 2018 Elijahu Ben-Michael
>
> Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
> associated documentation files (the "Software"), to deal in the Software without restriction,
> including without limitation the rights to use, copy, modify, merge, publish, distribute,
> sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
> furnished to do so, subject to the following conditions:
>
> The above copyright notice and this permission notice shall be included in all copies or
> substantial portions of the Software.
>
> THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
> NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
> NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
> DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
> OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## Classic synthetic control (R/scm_engine.R, R/scm_plots.R, R/scm_help_text.R, R/scm.b.R)

Copied from the author's standalone module *synthcontrol* (GPL-3), which is unchanged. Its author (the author of this module) releases the copy in this module under GPL (>= 2).

## Example data

| File | Source | Licence of the source |
|---|---|---|
| `prop99_long.csv` | Abadie, Diamond and Hainmueller (2010), California tobacco sales | replication data of the paper |
| `prop99_covariates.csv` | `smoking` dataset of the R package *tidysynth* (the same replication data with covariates) | tidysynth: MIT + file LICENSE |
| `west_germany.csv` | `scpi_germany` of the R package *scpi* (Cattaneo, Feng, Palomba and Titiunik); data of Abadie, Diamond and Hainmueller (2015) | scpi: GPL-3 |
| `mpdta_minwage.csv` | `mpdta` of the R package *did* (Callaway and Sant'Anna) | did: GPL-3 |

## R packages the module depends on

jmvcore, R6, ggplot2, quadprog (GPL >= 2), osqp (Apache 2.0), Matrix, glmnet (GPL-2 only, used by the
Elastic-Net analysis; this is why the module is licensed GPL (>= 2) and not GPL-3). *osqp* (Apache 2.0) is
installed as a separate R package and called through its R interface, not copied into the module.
