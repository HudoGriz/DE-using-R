setwd("/home/zorya/Documents/bioinfo_sel_ascpects/ferets/")

library(affyPLM)
library(vsn)
library(limma)
library(canine2.db)
library(Vennerable)
library(gplots)


source("ferets_functions.R")

## prepering data
targets_f <- read.table(
  "sample_descripton.txt",
  sep = ","
)
targets_f <- cbind(
  targets_f,
  read.table("sample_names.txt")
)
names(targets_f) <- c(
  "tissue",
  "treatment",
  "time",
  "repetition",
  "samples"
)

targets_f$treatment <- gsub("-", "_", targets_f$treatment)
targets_f$samples <- as.character(targets_f$samples)

write.csv(
  targets_f,
  "targets.csv",
  row.names = FALSE
  )

ab_f <- ReadAffy(
  filenames = paste0(
    "samples/", targets_f$samples, ".CEL.gz"
    )
  )

abcdf_f <- cdfName(ab_f)
cat("CDF:", abcdf_f, "\n")

# create chip images & pseudo-images & histogram
do.plot.all(
  ab = ab_f,
  targ = targets_f
)

ab.gc_f <- bg.adjust.gcrma(
  ab_f,
  affinity.source = "reference",
  type = "affinities",
  GSB.adjust = T,
  fast = F,
  optical.correct = F
)
# vsn between array normalization
ab.vsn_f <- justvsn(ab.gc_f)

# affter corection comparison
do.plot.all(
  targ = targets_f,
  ab = ab.gc_f,
  p = "corection/background"
)

do.plot.all(
  targ = targets_f,
  ab = ab.vsn_f,
  p = "corection/between",
  grf = FALSE,
  grm = FALSE
)

## probe summarization - obtaining probe set gene level estimate
## create data sets
sets <- list(ab_f, ab.gc_f, ab.vsn_f)
names <- c("raw", "background", "crossprobe")
data.sets <- do.data.sets(
  sets = sets,
  names = names
)

##################
## DE
print("Starting with DE")
tissues <- c("Lung", "Blood")
methods <- c("ls", "robust")
adj.method <- "BH"
coefs <- c("xaction", "ifn", "virus")
for (set.n in seq_along(data.sets)) {
  set <- data.sets[set.n]
  e <- set[[1]]$e
  pset_f <- set[[1]]$pset_f
  set.name <- names(set)
  print(set.name)
  # save(file = "data.sets.Rdata", data.sets)
  top.g.l <- list()
  all.genes <- c()
  for (tis in tissues) {
    if (tis == "Blood") {
      days <- 2
    } else {
      days <- unique(targets_f$time)[-1]
    }

    for (day in days) {
      print(paste0("processing: ", tis, " day-", day))

      d_targets <- targets_f[
        (targets_f$time == day | targets_f$time == 0)
        & targets_f$tissue == tis, ]

      ## differentialy expression analysis with limma
      print("differentialy expression analysis with limma")
      levs <- sort(unique(d_targets$treatment))
      L. <- factor(d_targets$treatment, levels = levs)

      ## experiment design
      print("experiment design")
      design <- model.matrix(~ -1 + L.)
      rownames(design) <- L.

      ## contrast
      print("contrast")
      cont.matrix <-
        makeContrasts(
          base = L.uninfected,
          virus = L.SARS_CoV - L.uninfected,
          ifn = L.IFN_a2b - L.uninfected,
          xaction = L.SARS_CoV - L.uninfected - (L.IFN_a2b - L.uninfected),
          levels = design
        )

      for (method in methods) {
        ## linear model fitting
        print(paste0("linear model fitting, method: ", method))
        sel <- paste0(d_targets$samples, ".CEL.gz")
        fit <- lmFit(e[, colnames(e) %in% sel], design, method = method)
        fit2 <- contrasts.fit(fit, cont.matrix)
        fit2.eb <- eBayes(fit2)

        tt_l <- list()

        for (c in coefs) {
          print(paste("toptable for coef:", c))
          tt <- get.tt(
            fit = fit2.eb,
            adj.method = adj.method,
            co = c,
            so = "P",
            b = TRUE,
            g = TRUE
          )

          n_sig.genes_f <- sum(tt$adj.P.Val < 0.05, na.rm = TRUE)
          topgenes <- tt[tt[, "adj.P.Val"] < 0.05, ]
          topgenes <- na.omit(topgenes)
          topups <- topgenes[topgenes[, "logFC"] > 1, ]
          topdowns <- topgenes[topgenes[, "logFC"] < -1, ]
          top.g.l[[paste0(tis, ".day_", day)]] <- topgenes

          all.genes <- c(all.genes, rownames(topgenes))
          tt_l[[c]] <- tt
        }

        ## graphs
        do.DGEplots()
      }
    }
  }
}
print("Process compleated!")