# Version 2.5.3
## Major Changes
* Added lambda statistic

## Bug Fix
* Corrected log to be log(x+1) and changed CPM(x+1) to CPM(x)

# Version 2.5.2
## Major Changes
* Added edgeR as a normalization method
* Added edgeR as a differential expression analysis method

# Version 2.5.1
## Major Changes
* Added TB data example

## Minor Changes
* Changed sapply to vapply in nb check

# Version 2.2.5
## Minor Changes
* Added sva batch correction method (for unknown variation correction)

# Version 2.2.4
## Major Changes
* Updated variation ratio analysis to be log transformed
* Added umap exploratory option

# Version 2.2.3
## Major Changes
* Added limma as a batch correction
* Added umap plot option
* Added ellipses to pca plot

## Minor Changes
* Set p-val plot x scale to always be 0 to 1
* Removed Intercept from p-val violin plots

# Version 2.2.2
## Minor Changes
* Updated SE object upload to allow assays of any name (no longer require one
  assay to be called "counts")
  
# Version 2.2.1
## Minor Changes
* Added pval information to the DESeq2 binomial evaluation

# Version 2.1.6
## Bug Fixes
* Corrected code for proper division of less than or 20+ samples

#Version 2.1.5
## bug Fixes
* Coerce variable of interest to a factor

# Version 2.1.4
## Major Changes
* Added negative binomial check for 20+ samples to DESeq2

# Version 2.1.3
## Major Changes
* Added negative binomial check for less than 20 samples to DESeq2

# Version 2.1.2
## Major Changes
* Added Variation Ratio Statistic to the explained variation tab

## Minor Changes
* Removed extra "Samples" column from example data
* Uploaded bladder example data batch variable as a factor

# Version 2.1.1
## Bug Fixes
* Updated imports to include shinyjs
* Updated imports to remove dendextend which is no longer utilized
* Corrected typos in Intro vignette

# Version 2.0.0

## Bug Fixes
* Fixed various errors for the Bioconductor build

## Major Changes
* Created interactive Shiny interface
* Allow user upload of data
* Allow user download of Batch corrected and/or normalized data
* Added Example Data Functionality

## Minor Changes
* Added a `NEWS.md` file to track changes to the package.

## Changes made to dendrogram
* Added `dendrogram_color_palette.R` for coloring dendrogram
* Updated `dendrogram.R` allowing batch & category to plot together
