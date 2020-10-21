generate_marker_genes <- function(method,data){
  switch(method,
         cellassign = markergene.cellassign(data),
         seurat = markergene.seurat(data),
         stop("No such marker gene generating method")
  )
}

markergene.cellassign <- function(data){
  
  require(org.Hs.eg.db)
  require(edgeR)
  ensembl_map <- dplyr::transmute(as_tibble(rowData(data)),ENSEMBLID=EnsembleId,SYMBOL=geneName)
  #> 'select()' returned 1:1 mapping between keys and columns
  gene_annotations <- ensembl_map %>%
    dplyr::rename(GeneID=ENSEMBLID,
                  Symbol=SYMBOL)
  cell_lines <- unique(colData(data)$label)
  dge <- DGEList(counts = counts(data), 
                 group = colData(data)$label, 
                 genes = gene_annotations, 
                 remove.zeros = TRUE)
  genes_to_keep <- rowSums(cpm(dge$counts) > 0.5) >= 2
  dge_filt <- dge[genes_to_keep,]
  dge_filt <- calcNormFactors(dge_filt, method="TMM")
  design <- model.matrix(~ 0+dge_filt$samples$group)
  colnames(design) <- levels(dge_filt$samples$group)
  v <- voom(dge_filt, design)
  fit <- lmFit(v, design)
  args <- map(1:ncol(combn(cell_lines,2)), ~{str_glue("{combn(cell_lines,2)[,.][[1]]} - {combn(cell_lines,2)[,.][[2]]}")})
  args$levels <- design
  contrast.matrix <- do.call('makeContrasts',args)

  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  tt <- topTable(fit2, n=Inf)
  tt_sig <- tt %>%
    dplyr::filter(adj.P.Val < 0.05)
  ##########
  
  
  
  lfc_table <- tt_sig[,c("H2228...H1975", "HCC827...H1975")]
  lfc_table <- lfc_table %>%
    dplyr::mutate(H1975=0,
                  H2228=H2228...H1975,
                  HCC827=HCC827...H1975) %>%
    dplyr::select(H1975, H2228, HCC827)
  rownames(lfc_table) <- tt_sig$GeneID
  lfc_table <- as.matrix(lfc_table)
  lfc_table <- lfc_table - rowMins(lfc_table)
  lfc_table <- as.data.frame(lfc_table)
  binarize <- function(x, threshold) {
    x[x <= threshold] <- -Inf
    x[x > -Inf] <- 1
    x[x == -Inf] <- 0
    return(x)
  }
  # Find the biggest difference
  maxdiffs <- apply(lfc_table, 1, function(x) max(diff(sort(x))))
  
  #
  thres_vals <- apply(lfc_table, 1, function(x) sort(x)[which.max(diff(sort(x)))])
  expr_mat_thres <- plyr::rbind.fill(lapply(1:nrow(lfc_table), function(i) {
    binarize(lfc_table[i,], thres_vals[i])
  }))
  rownames(expr_mat_thres) <- rownames(lfc_table)
  marker_gene_mat <- expr_mat_thres[(maxdiffs >= quantile(maxdiffs, c(.99))) 
                                    & (thres_vals <= log(2)),] %>%
    as.matrix
  suppressMessages({
    symbols <- plyr::mapvalues(
      rownames(marker_gene_mat),
      from = gene_annotations$GeneID,
      to = gene_annotations$Symbol
    )
  })
  
  is_na <- is.na(symbols)
  
  marker_gene_mat <- marker_gene_mat[!is_na,]
  rownames(marker_gene_mat) <- symbols[!is_na]
}

markergene.seurat <- function(data){
  
}