sub <- function(tab, trash = 0.05) {
  rownames(tab[tab$adj.P.Val < trash, ])
}

get.anotation <- function(ref = canine2SYMBOL, tt) {
  anotation <- toTable(ref)

  po <- match(row.names(tt), anotation$probe_id)
  g.names <- anotation[, 2][po]

  return(g.names)
}

get.tt <- function(fit,
                   adj.method,
                   co = "xaction",
                   so = "P",
                   b = TRUE,
                   g = TRUE,
                   ref = canine2SYMBOL) {
  tt <- topTable(
    fit = fit,
    adjust = adj.method,
    coef = co,
    sort.by = "P",
    number = Inf
  )

  if (b == TRUE) {
    tt$P <- 1 / (1 + exp(-tt$B))
  }

  if (g == TRUE) {
    tt$gene <- get.anotation(ref = ref, tt)
  }

  return(tt)
}

get.plot.raw <- function(targ,
                         ab,
                         p = "raw") {
  for (i in 1:length(targ$samples))
  {
    name <- paste(
      "./images/",
      p,
      "/chip_images/",
      i, ".jpg",
      sep = ""
      )
    jpeg(name)
    image(ab[, i], main = targ$samples[i])
    dev.off()
  }
}

get.plot.fit <- function(targ,
                         ab,
                         p = "raw",
                         data.fited = FALSE) {
  if (data.fited) {
    pset_f <- ab
  } else {
    pset_f <- fitPLM(ab, background = F, normalize = F)
  }

  for (i in 1:length(targ$samples))
  {
    name <- paste(
      "./images/",
      p,
      "/chip_pseudo_images/",
      i, ".jpg",
      sep = ""
      )
    jpeg(name)
    image(pset_f, which = i, type = "resids",
          main = targ$samples[i], add.legend=TRUE)
    dev.off()
  }
}

get.plot.hist <- function(targ,
                          ab,
                          p = "raw",
                          t = "") {
  color <- colours()[1:length(targ$samples)]
  name <- paste0(
    "./images/",
    p,
    "/histogram",
    t, ".jpg"
    )
  jpeg(name)
  hist(
    ab,
    lwd = 2,
    which = "both",
    col = color,
    ylab = "Density",
    xlab = "Log2 intensities",
    main = paste0("Histogram of ", p, " data")
  )
  dev.off()
}

get.plot.ma <- function(targ,
                        ab,
                        p = "raw") {
  for (i in 1:length(targ$samples))
  {
    name <- paste(
      "./images/",
      p,
      "/MAplot/",
      i, ".jpg",
      sep = ""
      )
    jpeg(name)
    MAplot(ab, which = i)
    dev.off()
  }
}

do.plot.all <- function(
                        targ,
                        ab,
                        gr = TRUE,
                        grf = TRUE,
                        grh = TRUE,
                        grm = TRUE,
                        p = "raw",
                        data.fited = FALSE) {
  if (gr) {
    get.plot.raw(
      targ,
      ab,
      p = p
    )
    print("get.plot.raw DONE!")
  }
  if (grf) {
    get.plot.fit(
      targ,
      ab,
      p = p,
      data.fited
    )
    print("get.plot.fit DONE!")
  }
  if (grh) {
    get.plot.hist(
      targ,
      ab,
      p = p
    )
    print("get.plot.his DONE!")
  }
  if (grm) {
    get.plot.ma(
      targ,
      ab,
      p = p
    )
    print("get.plot.ma DONE!")
  }
}

do.data.sets <- function(
                         sets,
                         names) {
  data.sets <- list()
  for (d in 1:length(sets)) {
    pset_f <- fitPLM(sets[[d]], background = F, normalize = F)
    e <- coefs(pset_f)
    x <- list(
      pset_f = pset_f,
      e = e
    )
    data.sets[[names[d]]] <- x
  }
  return(data.sets)
}

get.pname <- function(type,
                      format = "jpg",
                      tis. = tis,
                      day. = day,
                      set.name. = set.name,
                      method. = method) {
  name <- paste0(
    "./images/DGEresults/",
    set.name., "/",
    method., "/",
    type, "/",
    tis., ".day-", day., ".", format
  )
}

###
do.DGEplots <- function() {
  ## graphs
  # heat map
  print("heat map")
  eh <- e[rownames(topgenes), sel]
  row.names(eh) <- get.anotation(tt = eh)
  titl <- paste0(tis, "s ", "Day-", day)

  col <- apply(d_targets, 1, function(x) {
    paste0(x[2], ".rep", x[4])
  })
  colnames(eh) <- col

  name <- get.pname(
    type = "heatmap",
    format = "png"
  )
  png(name)
  heatmap.2(
    eh,
    col = redgreen(75),
    scale = "row",
    key = TRUE,
    symkey = TRUE,
    density.info = "none",
    trace = "none",
    cexRow = 0.8,
    cexCol = 0.6,
    main = titl
  )
  dev.off()

  # volcano
  print("volcano")
  name <- get.pname(type = "volcano")
  jpeg(name)
  volcanoplot(
    fit2.eb,
    coef = 3,
    highlight = 10,
    main = titl
    )
  dev.off()

  # Venn
  print("Venn")
  tt_v <- list(
    ifn = sub(tt_l$ifn),
    virus = sub(tt_l$virus)
  )

  vv <- Venn(tt_v,)
  X11()
  plot(
    vv,
    doWeights = FALSE
  )
  grid.text(titl, y=0.9, gp=gpar(col="red", cex=2))
  savePlot(
    filename =
      get.pname(
        type = "venn",
        format = "png"
      )
  )
  dev.off()

  print("Finished!")
}
