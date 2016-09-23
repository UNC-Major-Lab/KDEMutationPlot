# KDEMutationPlot
Creates a figure using kernel density estimation to show the mutational landscape of a protein.

You might need to install a few packages for the script to work: Biostrings, ggplot2, zoo, and gridExtra.

```R
install.packages(c("ggplot2","zoo","gridExtra"))
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")
```

You can run the script via RStudio after updating the hardcoded arguments, or you can run it from the command line using:

```bash
Rscript mutation_kde.R example_data/mutations.txt example_data/aa.txt example_data/domain.txt example_data/P29597.fasta example_data/test.pdf
```

Example output (before manually adjusting label positions):

![Image of example output] (https://github.com/UNC-Major-Lab/KDEMutationPlot/blob/master/example_data/TYK2_KDE.png)
