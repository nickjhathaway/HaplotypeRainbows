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

```r

```
