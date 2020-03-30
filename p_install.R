if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("affyPLM")
BiocManager::install("vsn")
BiocManager::install("limma")
BiocManager::install("canine2.db")

# there are some coplications with Vennerable
# if not working try to instal the following
BiocManager::install("RBGL")
BiocManager::install("graph")
install.packages("reshape")
install.packages(
  "Vennerable",
  repos="http://r-forge.r-project.org",
  type="source"
  )

install.packages("gplots")
install.packages("MASS")
install.packages('shiny')
install.packages("magick")
