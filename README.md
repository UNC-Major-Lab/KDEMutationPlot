# KDEMutationPlot
Creates a figure using kernel density estimation to show the mutational landscape of a protein.

You might need to install a few packages for the script to work: Biostrings, ggplot2, zoo, and gridExtra.

```R
install.packages(c("ggplot2","zoo","gridExtra"))
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")
```

Example output (before manually adjusting label positions):

![Image of example output] (https://github.com/UNC-Major-Lab/KDEMutationPlot/blob/master/example_data/TYK2_KDE.png)
