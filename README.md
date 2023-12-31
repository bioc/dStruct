# dStruct: method for identifying differential reactive regions from RNA structurome profiling data.

dStruct is a statistical method for identifying regions that display altered reactivity patterns between two groups of samples. It can perform *de novo* discovery or identify regions from a list provided by the user. The latter case is called *guided discovery*. dStruct is compatible with a diverse range of structure profiling technologies, accounts for biological variation and controls the false discovery rate in a multiple testing context.

## Getting started

To start with, download and install the latest versions of [R](https://cran.r-project.org/) and [RStudio](https://www.rstudio.com/products/rstudio/).

**To install via Bioconductor, please visit [http://www.bioconductor.org/packages/release/bioc/html/dStruct.html](http://www.bioconductor.org/packages/release/bioc/html/dStruct.html)**

dStruct can be installed directly from source. First, install `devtools` by executing the following in RStudio.

`install.packages("devtools")`

Next, in RStudio run the following.

    devtools::install_github("dataMaster-Kris/dStruct")

This should install the package in R. Check by executing the following command.

`library(dStruct)`

***

**Constructive feedback is most welcome! If you have issues to report, or feature requests, or any questions, please file issues or initiate discussion on this GitHub page instead of emailing. This will be beneficial for other users as well. Thanks for your interest!**

***

## Usage

dStruct takes reactivities of groups of samples, say A and B, for all transcripts under consideration. The reactivities for each transcript must be a data frame object with columns labeled as A1, A2, ... and B1, B2, ... The numerals in column names indicate sample number. Let us say that `reactivity` represents one such data frame, that there are 3 samples of each group and that the user needs to search for a minimum length of 11 nt. Note that this is the default value for this parameter. For other parameters that can be specified, refer the manual. Unless specified, `dStruct` assumes that the user intends to run the analysis with default settings. _De novo_ discovery of differential regions can be done by executing the following command.

`dStruct(reactivity, reps_A = 3, reps_B = 3, min_length = 11)`

Reactivities for all transcripts can be stored together in a list object, say `rlist`, with one data frame for each transcript. All the list elements must have unique names. In this case,  _de novo_ discovery can be done simultaneously for all transcripts.

`dStructome(rlist, reps_A = 3, reps_B = 3, min_length = 11)`

Users can specify the number of cores to be used by `dStructome` for parallel processing of transcripts.

For guided discovery, data frame with reactivities should contain only information for regions of interest. There should be a separate data frame for each region. Let one such data frame be `reac_region1`. If there are say, 3 samples of each group, guided discovery can be performed as follows.

`dStructGuided(reac_region1, reps_A = 3, reps_B = 3)`

All the data frames can also be combined in a single list object, say `rlist`. Then, guided discovery could be done in parallel for all regions using the following command.

`dStructome(rlist, reps_A = 3, reps_B = 3, method = "guided")`

For other options available in dStruct package, refer the dStruct manual or email the contributors.

## Citation

Choudhary, K., Lai, YH., Tran, E. J., Aviran, S. [dStruct: identifying differentially reactive regions from RNA structurome profiling data.](https://doi.org/10.1186/s13059-019-1641-3) Genome Biology, 20, 40 (2019). 
