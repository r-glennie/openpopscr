# openpopscr
Fits open population spatial capture-recapture with
activity centre movement by maximum likelihood.

## Install 

In R, the latest release can be install using the <code>devtools</code> package
with the command: 

```
devtools::install_github("r-glennie/openpopscr@v1.1.0", build_vignettes = TRUE)
```

The package requires you have a <code>C</code> compiler installed on your system. 
Windows users may need to install <code>R-tools</code> for this reason. 
It is assumed Linux and Mac users have a compiler installed. 
For parallel processing and faster performance, install <code>openMP</code>. 

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
browseVignettes()
```
To learn how to use the package, I recommend reading the <code>ScrData</code>
and <code>ScrModel</code> vignettes first. The rest are very similar to these
vignettes, making it easier to use other classes once the basics are familiar. 

## Issues 

This package is under development. I'd be grateful for anyone who has issues
with the package to open an issue or contact me. 


