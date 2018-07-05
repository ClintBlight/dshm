# dshm
Density surface Hurdle modelling

This package contains functions to fit density surface models using a Hurdle approach. It also contains useful function to split transect lines into segments, buffer segments, and check for overlapping segments.

To install the dshm package:

1) Be sure you have the devtools and countreg packages installed. Devtools is available on CRAN and it can be installed running install.packages("devtools"). Contreg can be download as .tar file here in the /local_repository folder and installed in R studio (go to Packages->Intall->Install Archive File->Search for the downloaded .tar file). You can also install countreg by running install.packages("countreg", repos="http://R-Forge.R-project.org").

2) Install dshm by typing devtools::install_github("FilippoFranchini/dshm")

Note: Parallelization is currently available only for macs. Make sure you have the doMC package installed.

Enjoy dshm
