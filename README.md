
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![minimal R
version](https://img.shields.io/badge/R%3E%3D-4.2.3-brightgreen.svg)](https://cran.r-project.org/)
[![Licence](https://img.shields.io/github/license/mashape/apistatus.svg)](http://choosealicense.com/licenses/mit/)

## Research compendium for selection coefficients of malaria drug resitance

This is a working R compendium (think R package but for reproducible
analysis). A good overview on research compendiums, see the [R for
Reproducible Research](https://annakrystalli.me/rrresearch/index.html)
course.

### Installation

    git clone https://github.com/bailey-lab/selmar.git
    cd selmar
    open selmar.Rproj

Next, if `renv` has been used in this repository (look out for
`renv.lock`) then use `renv::restore` to set up package dependencies.
Otherwise `devtools::install_dev_deps()` will install all required
packages, as specified in the Imports in DESCRIPTION.

### Overview

The structure within analysis is as follows:

    analysis/
        |
        ├── 01_xxxxx /           # analysis scripts used for generating figures
        |
        ├── figures/              # location of figures produced by the analysis scripts
        |
        ├── data/
        │   ├── DO-NOT-EDIT-ANY-FILES-IN-HERE-BY-HAND
        │   ├── raw_data/       # data obtained from elsewhere
        │   └── derived_data/   # data generated during the analysis

### Compendium DOI:

<https://doi.org/X/X>

The files at the URL above will generate the results as found in the
publication.

### The R package

This repository is organized as an R package. There are no/negligable R
functions exported in this package - the majority of the R code is in
the analysis directory. The R package structure is here to help manage
dependencies, to take advantage of continuous integration, and so we can
keep file and data management simple. For any R packages that are used
frequently in this repository, they are documented in `R/` and are used
in the analysis folder using `devtools::load_all()`.

To download the package source as you see it on GitHub, for offline
browsing, use this line at the shell prompt (assuming you have Git
installed on your computer):

``` r
git clone https://github.com/bailey-lab/selmar.git
```

Once the download is complete, open the `selmar.Rproj` in RStudio to
begin working with the package and compendium files. We will endeavour
to keep all package dependencies required listed in the DESCRIPTION.

<!-- To add this once all the analysis is done -->

In addition we use `renv` to track package dependencies for
reproducibility. Please use `renv::restore` to restore the state of the
project and see <https://rstudio.github.io/renv/articles/renv.html> for
more information.

### To-Dos
Week of Monday 05/08
- [ ] Remove y-log and check figures Uganda
- [ ] Remove y-log and check figures SE Asia
- [ ] Finalize the predictions for SE Asia
- [ ] Finalize the predictions for Uganda
- [ ] Finalize the metrics of predictions for SE Asia
- [ ] Add forecasting to methods
- [ ] Add forecasting to conclusion
- [ ] Add forecasting to discussion
- [ ] Update discussion
- [ ] Update Figure 2
- [ ] Finalize literature review
- [x] Finalize supplement

Week of Monday 05/08
- [ ] Polishing Manuscript
- [ ] Supplement figures
- [ ] Take Manuscript offline and add all references
- [ ] Format Code for GitHub
- [ ] Finalize ReadMe on GitHub

### Licenses

Code: [MIT](http://opensource.org/licenses/MIT) year: 2023, copyright
holder: OJ Watson

Data: [CC-0](http://creativecommons.org/publicdomain/zero/1.0/)
attribution requested in reuse
