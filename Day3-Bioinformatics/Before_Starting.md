Because some of these packages may take long to download and isntall, its advisable to install them before starting the class
Just copy and paste the code chunks into R Studio and wait

# Install CRAN packages
install.packages(c(
  "tidyverse",
  "ggpubr",
  "pheatmap",
  "data.table"
))

# Install Bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install Bioconductor Packages
BiocManager::install(c(
  "Biobase",
  "GEOquery",
  "DESeq2",
  "org.Hs.eg.db",
  "ashr",
  "clusterProfiler",
  "msigdbr",
  "fgsea"
))
