# Convert accession numbers to taxonomy

<!--
[![Build Status](https://travis-ci.org/sherrillmix/taxonomizr.svg?branch=master)](https://travis-ci.org/sherrillmix/taxonomizr)
[![codecov](https://codecov.io/gh/sherrillmix/taxonomizr/branch/master/graph/badge.svg)](https://codecov.io/gh/sherrillmix/taxonomizr)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/taxonomizr)](https://cran.r-project.org/package=taxonomizr)
-->

## Introduction

`antiprofr` provides some simple functions to read in and display measurements of antibody profiling in serum.

The major functions are:
 * `readAnti`: read in antibody profiling data
 * `convertStackedToMat`: average stacked data to a single value per patient-antigen
 * `plotAnti`: plot the raw data from positive controls and patient serum
 * `plotHeat`: plot the average patient-antigen values

And a simple use case might look like (see below for more details):


```r
library(antir)
raw<-read.csv('antibodyData.csv',stringsAsFactors=FALSE,header=FALSE)
stacked<-readAnti(raw,p24Cut=10)
odMat<-convertStackedToMat(stacked)
par(mar=c(6,10,1,3))
plotHeat(odMat,main='',filterLess=200,filterMore=50000,scaleMain='Dilution reaching OD450(p24=10pg)')
```

## Installation
The package is available from github, use the [<code>devtools</code>](https://github.com/hadley/devtools) library and run:

```r
devtools::install_github("sherrillmix/antir")
```

To use the library, load it in R:

```r
library(antir)
```

## Usage

### Data format

This package is designed for antibody profiling data where each row given the OD for various dilutions of a patient serum with a given antigen (or positive or negative control). The data should be in a 15 column .csv where the first column gives the patient ID (blanks are filled down), the second column gives the antigen and columns 3-14 give OD values for given dilutions of the patient serum for that antigen.

### Reading in data
In most cases, the user will be reading in data from a .csv file (if an Excel file then Save As .csv in Excel or OpenOffice):


```r
raw<-read.csv('antibodyData.csv',stringsAsFactors=FALSE,header=FALSE)
```

Here we'll use data cached from the package:


```r
raw<-antibodyData
head(raw)
```

```
##        V1               V2     V3     V4     V5     V6     V7     V8
## 1 L016-01               BB 0.1018 0.0740 0.0663 0.0588 2.0930 1.0313
## 2                    HSV-1 0.1192 0.0861 0.0781 0.0662 0.1205 0.0819
## 3                    CMV-N 0.2488 0.1784 0.0901 0.0743 0.1951 0.1198
## 4                    CMV-C 0.2174 0.1201 0.0874 0.0720 0.2575 0.1087
## 5                    HSV-2 0.4526 0.2213 0.1303 0.0945 0.4012 0.2289
## 6         Negative Control 0.1906 0.0928 0.0832 0.0690 0.1350 0.0877
##       V9    V10    V11    V12    V13    V14
## 1 0.5778 0.2599 2.0792 1.0701 0.5203 0.2457
## 2 0.0785 0.0630 0.1156 0.0811 0.0674 0.0630
## 3 0.0928 0.0712 0.2090 0.1226 0.0896 0.0765
## 4 0.0903 0.0701 0.2009 0.1171 0.0867 0.0718
## 5 0.1247 0.0933 0.3473 0.2005 0.1288 0.0856
## 6 0.0720 0.0623 0.1425 0.0845 0.0717 0.0640
```


