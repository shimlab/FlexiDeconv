# FlexiDeconv: flexible cell type deconvolution of the Spatial Transcriptomics

FlexiDeconv is a semi-supervised cell type deconvolution method of the Spatial Transcriptomics that infers cell type proportions across pixels by flexibly leveraging the reference, which consists of cell type gene expression profiles typically derived from scRNA-seq data. The advantage of FlexiDeconv lies within the usage of reference as prior information, while actively estimating cell type gene expression profiles from the spatial data. FlexiDeconv allows users to assign weights to the provided cell type gene expression profiles, reflecting the extent to which reference contributes to estimating gene expression relative to the ST data.

This package implements FlexiDeconv, please report any bugs to [issues](https://github.com/shimlab/FlexiDeconv/issues).

# License

All source code and software in this repository is free software; you can redistribute it and/or modify it under the terms of the [GNU General Public License](https://www.gnu.org/licenses/gpl-3.0.html) as published by the [Free Software Foundation](https://www.fsf.org/); either version 3 of the License, or (at your option) any later version. See the [LICENSE](LICENSE) file for the full text of the license.

# Installation

Currently, `FlexiDeconv` is only available in R, via `remotes`. It can be installed as follows:

``` r
require(remotes)
remotes::install_github('shimlab/FlexiDeconv')
```

Note that the R version must be \>= 4.3.0

Installation should be completed within a few minutes.

## Vignette

If you wish to view the vignette, please follow the following step.

1.  Install `FlexiDeconv` while enabling `build_vignettes`.

``` r
require(remotes)
remotes::install_github('shimlab/FlexiDeconv', build_vignettes = TRUE)
```

2.  Load the package.

``` r
library(FlexiDeconv)
```

3.  View the vignette (or also view [here](https://shimlab.github.io/FlexiDeconv/FlexiDeconv_Vignette.html)).

``` r
vignette("FlexiDeconv_Vignette")
```

## Other Vignette

- [Interpreting results] (https://shimlab.github.io/FlexiDeconv/Inference_Vignette.html)
