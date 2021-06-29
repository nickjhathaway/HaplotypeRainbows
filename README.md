# Haplotype Rainbow

This is a collection of tools in R to create haplotype "rainbows" with a myriad input as opposed to SNP barcodes which only have binary input. 

# Install  

Install the most recent release using devtools  

```r
devtools::install_github('nickjhathaway/rainbowHaplotypes')
```
Or the developmental branch 

```r
devtools::install_github('nickjhathaway/rainbowHaplotypes@develop')
```

# Input 

The tools need a minimum of 4 columns 

1.  Sample name column
2.  A loci name column 
3.  A population id for the haplotype for that loci
4.  A within sample relative abundance of the haplotype


These tools were developed to work with [SeekDeep](https://github.com/bailey-lab/SeekDeep) data output and so it assumes the default for the columns above are as follows:


1.  **s_Sample** - Sample name column
2.  **p_name** - A loci name column 
3.  **h_popUID** - A population id for the haplotype for that loci
4.  **c_AveragedFrac** - A within sample relative abundance of the haplotype

But you can specify column names when using the functions 

## Prep
There are two main steps for creating the haplotypes, first is to process the input data to create new data table to create the figures from. The library relies heavily on the [tidyverse](https://www.tidyverse.org/) library.  

Below uses example data set from *Plasmodium falciparum* lab isolates and mixtures  
 
```r
library(rainbowHaplotypes)

# load example data 
data("pfisolateExample") 

#prep data 
pfIsosHeomeV1_prep = prepForRainbow(pfIsosHeomeV1)
 
#prep data when columns have other than default names 

pfIsosHeomeV1_otherName_prep = prepForRainbow(pfIsosHeomeV1_otherName,
                                              sampleCol = Sample, 
                                              targetCol = loci, 
                                              popUIDCol = ID, 
                                              relAbundCol = freq)

```

## Plotting  

Below creates a ggplot object from the prep data, this object can either be plotted or further manipulated as needed first  

```r
genRainbowHapPlotObj(pfIsosHeomeV1_prep)

genRainbowHapPlotObj(pfIsosHeomeV1_otherName_prep,
                                              sampleCol = Sample, 
                                              targetCol = loci, 
                                              popUIDCol = ID, 
                                              relAbundCol = freq) 

```

This creates a plot where with samples on the y-axis and targets/loci on the x-axis. The within sample frequencies will be taken into account and will adjust the bars accordingly. The colors have meaning in each color, e.g. the same color within a column is the same haplotype, but colors across columns don't relate to each other. 

The package was developed so the colors denoting each major haplotype slightly in hue in each column/loci which ends up creating a repeating "rainbow" across (with default period of 11). 

![example]()

### Manipulating plotting  

By default the plotting 

```r

```
