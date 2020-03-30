# DE-using-R

#### What's this I found?

This is a basic walkthrough on how to perform a differential gene expression analysis (**DE**).

A in detail description of steps taken and commentary on results obtained you can fined [here](https://hudogriz2.shinyapps.io/report/) (it needs some time for the server to build the report).

It uses packages from *Bioconductor*, mostly `affyPLM` and `limma`. In the file `p_install.R` are installation commands contained for all the packages required to complete the analysis, and obtain fancy plots.

#### How can I do it myself?

To get started you'll need to either clone the repository or download the files. Important are the `run_me`, `feret_script.R` and `feret_functions.R`. 

`run_me` which is a _bash_ script, downloads the data (not included here), creates the components for the “targets” file and starts the `feret_script.R` script. from here the magic happens. `feret_script.R` contains the whole process and `ferets_functions.R` contains the functions used in the process.

`report.Rmd` contains the code for the report you can find live on the above link.

##### Have fun!
