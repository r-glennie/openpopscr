# openpopscr
Fits open population spatial capture-recapture with
activity centre movement by maximum likelihood.

## Install 

In R, the latest release can be install using the <code>devtools</code> package
with the command: 

```
devtools::install_github("r-glennie/openpopscr", build = TRUE, build_opts = c("--no-resave-data", "--no-manual"))
```

The package requires you have a <code>C</code> compiler installed on your system. 
Windows users may need to install <code>R-tools</code> for this reason. 
It is assumed Linux and Mac users have a compiler installed. 

## Use

The package consists of different <em>classes</em>. Each class describes a 
particular kind of <em>object</em>. Objects have properties and functions that
you can use. In this package, there are two main types of classes: models and
data. The <code>ScrData</code> class is used to create data objects; other
classes such as <code>ScrModel</code>, <code>CjsModel</code>, <code>ScrTransientModel</code>, are used to create model objects that you can fit to your data and 
obtain inference from. 

There is a vignette for each class. To see a list of all vignettes, use the
command: 

```
browseVignettes("openpopscr")
```
To learn how to use the package, I recommend reading the <code>ScrData</code>
and <code>ScrModel</code> vignettes first. The rest are very similar to these
vignettes, making it easier to use other classes once the basics are familiar.

Also, the folder inst/examples/ in the repository contains code showing how to fit
each type of model to simulated data. 

## Issues 

This package is under development. I'd be grateful for anyone who has issues
with the package to open an issue or contact me. 


