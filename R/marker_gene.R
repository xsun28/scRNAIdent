source('R/markergene_config.R')
markergene.generate_marker_genes <- function(method,data,file_name){
  switch(method,
         cellassign = markergene.cellassign(data,file_name),
         seurat = markergene.seurat(data,file_name),
         stop("No such marker gene generating method")
  )
}
###using cell assign to find marker gene
markergene.cellassign <- function(data,file_name){
  
  require(org.Hs.eg.db)
  require(edgeR)
  marker_config <- markergene.config.cellassign
  rowData(data)$geneName <- rownames(data)
  ensembl_map <- dplyr::transmute(as_tibble(rowData(data)),SYMBOL=geneName)
  #> 'select()' returned 1:1 mapping between keys and columns
  gene_annotations <- ensembl_map %>%
    dplyr::rename(Symbol=SYMBOL)
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
  base_cell_type <- if(purrr::is_null(marker_config$base_cell_type)) unique(colData(data)$label)[[1]] else marker_config$base_cell_type
  args <- purrr::map(1:ncol(combn(cell_lines,2)), ~{ if(combn(cell_lines,2)[,.][[1]]==base_cell_type)
                                                    return(str_glue("{combn(cell_lines,2)[,.][[2]]} - {combn(cell_lines,2)[,.][[1]]}"))
                                              str_glue("{combn(cell_lines,2)[,.][[1]]} - {combn(cell_lines,2)[,.][[2]]}")
                                              })
  args$levels <- design
  contrast.matrix <- do.call('makeContrasts',args)

  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  tt <- topTable(fit2, n=Inf)
  tt_sig <- tt %>%
    dplyr::filter(adj.P.Val < 0.05)
  ##########
  
  
  #####how to determine the base cell type?????

  diff_baseline_celltypes <- unlist(purrr::map(cell_lines[which(cell_lines!=base_cell_type)],~{str_glue("{.}...{base_cell_type}")}))
  lfc_table <- tt_sig[,diff_baseline_celltypes]
  colnames(lfc_table) <- purrr::map(diff_baseline_celltypes,~{str_split(.,"\\.\\.\\.")[[1]][1]})
  lfc_table[[base_cell_type]] <- 0
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
  write_rds(marker_gene_mat,str_glue("{file_name}.RDS"))
  marker_gene_mat
}

###using seurat to find marker gene
markergene.seurat <- function(data,file_name){
  require(Seurat)
  cnts <- counts(data)
  seuset <- CreateSeuratObject(cnts, project='suerat_marker_gene')
  seuset <- NormalizeData(object = seuset)
  Labels <- colData(data)$label
  Idents(seuset) <- Labels
  marker_config <- markergene.config.seurat
  DEgenes <- FindAllMarkers(seuset, only.pos = marker_config$only_pos, min.pct = marker_config$min_pct, logfc.threshold = marker_config$logfc_threshold)
  Markers <- matrix(nrow = marker_config$marker_gene_num,ncol = length(unique(Labels)))
  colnames(Markers) <- unique(Labels)
  for (i in unique(Labels)){
    TempList <- DEgenes$gene[((DEgenes$cluster == i) & (DEgenes$avg_logFC > 0))]
    MarkerGenes <- DEgenes$p_val_adj[DEgenes$cluster == i]
    print(MarkerGenes[1:marker_config$marker_gene_num])
    if (length(TempList) >= marker_config$marker_gene_num){
      Markers[,i] <- TempList[1:marker_config$marker_gene_num]
    }
    else{
      if(length(TempList) > 0){
        Markers[c(1:length(TempList)),i] <- TempList
      }
    }
  }
  selected_marker_genes <- unique(unlist(as.list(Markers)))
  marker_gene_mat <- matrix(0, length(selected_marker_genes), length(unique(Labels)))
  rownames(marker_gene_mat) <- selected_marker_genes
  colnames(marker_gene_mat) <- unique(Labels)
  for(type in colnames(Markers)){
    type_marker_genes <- Markers[,type]
    marker_gene_mat[,type][type_marker_genes] <- 1
  }
  write_rds(marker_gene_mat,str_glue("{file_name}.RDS"))
  marker_gene_mat
}