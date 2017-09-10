# openpopscr
Fits Jolly-Seber spatial capture-recapture by maximum likelihood.

## Download
Download the latest stable release: https://github.com/r-glennie/openpopscr/files/1290198/openpopscr_1.0.0.tar.gz

## Install
In R, if downloaded package is saved at "openpopscr_1.0.0.tar.gz" then install with command 
```
install.packages("openpopscr_1.0.0.tar.gz", repos = NULL)
```
In R studio, you can select the "Packages" tab, click "Install", select "Install from: Package Archive File", and select file "openpopscr_1.0.0.tar.gz". 

The package requires you have a C compiler installed on your system. Windows users may need to install R-tools for this reason. It is assumed Linux and Mac users have a compiler installed. 


## Help 
A vignette describing how to use the package is included: 
```
library(openpopscr)
vignette("openpopscr") 
```
