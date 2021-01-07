SYSTEM REQUIREMENTS:
The latest version of R (>3.6.3) should be installed from https://cran.r-project.org
and the latest version of JAGS (>4.3.0) should be installed from https://sourceforge.net/projects/mcmc-jags/files. 

INSTALLATION GUIDE:
Open the Script_NMI.R script in your R session and adapt the parameters of the analysis (lines 36-38).
Parameters used in Broennimann et al. 2020 are 
grain=10 # resolution of the climatic data. 
envelope="kde" # choise of the envelop methode
level=99 # level of inclusion of rare climatic conditions for kde and mve.
The calculation of NMI for all species should take about 1.6 hours with CPU 2.6 GHz / RAM 8 Go.
The Bayesian Mixed Models should take about 15 minutes with CPU 2.6 GHz / RAM 8 Go.

DEMO:
The example given in Fig 1C of Broennimann et al. 2020 for Alces alces can be performed by running lines 102-190.
The calculation time should be about 1 minute. 

INSTRUCTIONS FOR USE:
The script can easily be adapted to run on your own data. At lines 42-52 you can choose to import a different shapefile representing your native range(s), 
a different set of introduction data (provided that the file contains a column with success or failure of the introductions) or different climatic data