library(Biostrings)
library(ggplot2)
library(zoo)
library(gridExtra)

# Get command line arguments, or hardcode them
#args <- commandArgs(trailingOnly = TRUE)
args <- c("example_data/mutations.txt", "example_data/aa.txt", "example_data/domain.txt", "example_data/P29597.fasta", "example_data/TYK2_KDE.pdf")

# Theme that removes all axes, backgrounds, etc
theme_empty <- theme_bw()
theme_empty$line <- element_blank()
theme_empty$rect <- element_blank()
theme_empty$strip.text <- element_blank()
theme_empty$axis.text <- element_blank()
theme_empty$plot.title <- element_blank()
theme_empty$axis.title <- element_blank()
theme_empty$legend.position <- "none"
theme_empty$plot.margin <- structure(c(0, 0, 0, 0), unit = "lines", valid.unit = 3L, class = "unit")

# Read data files. See samples for file format
mutation.data <- read.table(args[1], sep="\t", header=T)  # mutations.txt is the sample input file
highlight.data <- read.table(args[2], sep="\t", header=T) # aa.txt is the sample input file
domain.data <- read.table(args[3], sep="\t", header=T)    # domain.txt is the sample input file
seq.data <- readAAStringSet(args[4], "FASTA")             # P29597.fasta is the sample input file

# Need to know full sequence length to set x-axis range
seq.length <- length(seq.data[[1]])

# Create new PDF file
pdf(file=args[5],width=7.0,height=3.5)                  # TYK2_KDE.pdf is the sample output file

##############################################################################
########## First Panel Begin #################################################
##############################################################################

# Only using missense mutations for kernel density estimation (KDE)
mutation.data.missense <- subset(mutation.data, Mutation.Type=="Substitution - Missense")

# Need to know maximum value of KDE to place vertical lines and labels for highlighted residues
y.max <- max(density(mutation.data.missense$Residue,adjust=.05)$y)

# Plot kernel density estimation and highlighted residues
# Adjust is to .05 because it seems to be an aesthetically pleasing value. There is no science behind the choice
p1 <- ( ggplot(data=mutation.data.missense, aes(x=Residue)) 
      + stat_density(adjust=.05, fill='#EEF1F6', color="black") 
      + geom_segment(data=highlight.data,aes(x=AA,y=0,xend=AA,yend=y.max), linetype = "dashed")
      + geom_text(data=highlight.data,aes(label=label, x=AA, y=y.max*1.1),size=4)
      + scale_x_continuous(limits=c(1,seq.length))
      + scale_y_continuous(limits=c(0,y.max*1.2))
      + theme_empty
      )

##############################################################################
########### First Panel End ##################################################
##############################################################################





##############################################################################
######### Second Panel Begin #################################################
##############################################################################

# Sorting the mutations by mutation type first so that each type will be contiguous
mutation.data <- mutation.data[with(mutation.data, order(Mutation.Type)), ]
# Count what number duplicate each residue is so they can be plotted at different y values
mutation.data$level <- rollapply(mutation.data$Residue, width=nrow(mutation.data), partial=T, align="right", FUN = function(x) sum(x==x[length(x)]))

# Adjusting margins to bring the graphs closer
theme.p2 <- theme_empty
theme.p2$plot.margin <- structure(c(-0.9, 0, 0, 0), unit = "lines", valid.unit = 3L, class = "unit")
theme.p2$legend.position <- "bottom"

# Plot individual mutations
p2 <- ( ggplot(data=mutation.data,aes(x=Residue,y=-level, color=Mutation.Type))
      + geom_point(shape="|",size=3) 
      + scale_x_continuous(limits=c(1,seq.length))
      + scale_y_continuous(limits=c(-(max(mutation.data$level+1)),0))
      + scale_color_brewer(palette = "Set1")
      + theme.p2
      )

##############################################################################
########## Second Panel End ##################################################
##############################################################################







##############################################################################
########## Third Panel Begin #################################################
##############################################################################

# Adding a fake domain that spans the length of the protein
full.length.domain <- data.frame(Domain.Name="",Start=1,End=seq.length,height=2)
domain.data <- cbind(domain.data,data.frame(height=4))
domain.data <- rbind(full.length.domain,domain.data)

# Adjusting margins to bring this graph above the legend of the previous graphs
theme.p3 <- theme_empty
theme.p3$plot.margin <- structure(c(-6, 0, 6, 0), unit = "lines", valid.unit = 3L, class = "unit")

# Plot the domain cartoon
p3 <- ( ggplot(data=domain.data, aes(xmin = Start, xmax = End, ymin = -height/2, ymax = height/2,fill=Domain.Name))
      + geom_rect(color="black")
      + geom_text(aes(x=Start+(End-Start)/2, y=0, label=Domain.Name), size=3)
      + scale_x_continuous(limits=c(1,seq.length))
      + scale_y_continuous(limits=c(-13,13))
      + scale_fill_brewer(palette = "Set2")
      + theme.p3
      )

##############################################################################
########### Third Panel End ##################################################
##############################################################################

# Arrange the graphs vertically
grid.arrange(p1, p2, p3, ncol=1)

# Close the PDF file
dev.off()