#' Infer ADAR activity with the dataset-specific regulon
#'
#' This function is used to get ADAR activity scores starting from an expression
#' matrix with genes as rows and samples as columns, and the gene regulation
#' network obtained by the expression matrix.
#'
#' @usage yourSignature(mat, net, genes, organism = "", file_name = "")
#' @param mat Normalized expression matrix used by ARACNe to compute the
#' network
#' @param net Network obtained from ARACNe
#' @param genes List of genes of interest
#' @param organism "hsapiens" or "mmusculus"
#' @param file_name How to nominate the output file
#' @return Matrix with the activity scores of the genes in the reg_list
#' @examples
#'
#' # Load the example data
#' data(example_network)
#' data(example_matrix)
#' data(example_genes)
#'
#' # Infer ADAR activity
#' res <- yourSignature(example_matrix, example_network, example_genes,
#' organism = "hsapiens", file_name = "interactors_scores")
#'
#' @import utils
#' @importFrom igraph graph_from_adjacency_matrix as_data_frame
#' @importFrom dplyr filter
#' @importFrom viper aracne2regulon viper
#' @export

yourSignature <- function(mat, net, genes, organism = "", file_name = "") {

  g  <- graph_from_adjacency_matrix(net, weighted = TRUE, mode = "undirected")
  dataf <- as_data_frame(g)

  if (organism == "hsapiens") {
    genes_selected <- c(regulators_list, genes)
  }
  else {genes_selected <- c(regulators_list_Mm, genes)}

  df_new <- data.frame()
  for (g in genes_selected) {
    df1 <- filter(dataf, dataf[, 1] == g)
    df2 <- filter(dataf, dataf[, 2] == g)
    df2 <- df2[, c(2, 1, 3)]
    colnames(df2) <- c("from", "to", "weight")
    df3 <- rbind(df1, df2)
    df_new <- rbind(df_new, df3)
  }

  df_list <- split(df_new, df_new$from)
  df_list <- lapply(df_list, function(dataf) {dataf[, -1]})

  file_adj_path <- file.path(getwd(), paste0(file_name, ".adj"))

  if (file.exists(file_adj_path)) {
    file.remove(file_adj_path)
  }

  i <- 1
  while (i < length(df_list)) {
    a <- c(names(df_list[i]), paste(df_list[[i]]$to, df_list[[i]]$weight))
    write.table(x = rbind(a), file = file_adj_path, row.names = FALSE,
                col.names = FALSE, quote = FALSE, append = TRUE, sep = " ")
    i <- i + 1
  }

  adj_file <- readLines(file_adj_path)
  modified_file <- gsub(" ", "\t", adj_file)
  writeLines(modified_file, file_adj_path)

  regulon <- aracne2regulon(file_adj_path, mat, verbose = FALSE)
  res <- viper(mat, regulon[grep("ADAR", ignore.case = TRUE, names(regulon))],
               verbose = FALSE, minsize = 0)

  res[which(tolower(rownames(res)) == "adar"), ] <-
    (res[which(tolower(rownames(res)) == "adar"), ] -
       min(res[which(tolower(rownames(res)) == "adar"), ])) /
    (max(res[which(tolower(rownames(res)) == "adar"), ]) -
       min(res[which(tolower(rownames(res)) == "adar"), ]))

  res[which(tolower(rownames(res)) == "adarb1"), ] <-
    (res[which(tolower(rownames(res)) == "adarb1"), ] -
       min(res[which(tolower(rownames(res)) == "adarb1"), ])) /
    (max(res[which(tolower(rownames(res)) == "adarb1"), ]) -
       min(res[which(tolower(rownames(res)) == "adarb1"), ]))

  res

}

#' Compute ADAR activity with the cancer-related signature
#'
#' The function returns ADAR activity scores for each sample/spot/cell depending
#' on the starting expression matrix given as input. It employs a signature
#' derived from the union of two datasets of A549 cell line samples, which is
#' ideal to measure ADAR activity in a cancer-related context.
#'
#' @usage tumoralSignature(mat, object)
#' @param mat Either a numeric expression matrix with genes as rows and samples
#' as columns or an ExpressionSet from a Seurat object, SingleCellExperiment or
#' SpatialExperiment
#' @param object Character string indicating the type of object, either bulk or
#' eset
#' @return Matrix with ADAR activity scores or eset object updated
#' @examples
#'
#' # Load the regulon
#' data(tumoral_regulon)
#'
#' # Load the example expression matrix
#' data(example_matrix)
#'
#' # Compute ADAR activity
#' res <- tumoralSignature(example_matrix, object = "bulk")
#'
#' @importFrom viper viper
#' @importFrom GSVA ssgseaParam gsva
#' @importFrom GSEABase GeneSetCollection GeneSet
#' @importFrom SummarizedExperiment assay
#' @export

tumoralSignature <- function(mat, object = c("bulk", "eset")) {

  if (object == "bulk") {
    regulon <- tumoral_regulon
    res <- viper(mat, regulon, verbose = FALSE, minsize = 0)
    res["ADAR", ] <- (res["ADAR", ] - min(res["ADAR", ])) /
      (max(res["ADAR", ]) - min(res["ADAR", ]))
    res["ADARB1", ] <- (res["ADARB1", ] - min(res["ADARB1", ])) /
      (max(res["ADARB1", ]) - min(res["ADARB1", ]))

    return(res)
  }

  else {

    if (class(mat) == "Seurat") {
      counts <- mat@assays$RNA$data
    }
    else {counts <- assay(mat, "logcounts")}

    tumoral_adar1_sign_genes <- names(tumoral_regulon$ADAR$tfmode)
    tumoral_adar1_gs <- GeneSet(tumoral_adar1_sign_genes,
                                setName = "tumoral_adar1")
    tumoral_adar2_sign_genes <- names(tumoral_regulon$ADARB1$tfmode)
    tumoral_adar2_gs <- GeneSet(tumoral_adar2_sign_genes,
                                setName = "tumoral_adar2")

    all_sets <- GeneSetCollection(list(tumoral_adar1_gs, tumoral_adar2_gs))
    ssgsea_obj <- ssgseaParam(counts, geneSets = all_sets)
    ssgsea <- gsva(ssgsea_obj)
    ssgsea[1, ] <- (ssgsea[1, ] - min(ssgsea[1, ])) / (max(ssgsea[1, ]) -
                                                         min(ssgsea[1, ]))
    ssgsea[2, ] <- (ssgsea[2, ] - min(ssgsea[2, ])) / (max(ssgsea[2, ]) -
                                                         min(ssgsea[2, ]))

    mat[["adar1_activity"]] <- ssgsea[1, ]
    mat[["adar2_activity"]] <- ssgsea[2, ]

    return(mat)
  }

}

#' Compute ADAR activity with the neuronal signature
#'
#' The function returns ADAR activity scores for each sample/spot/cell depending
#' on the starting expression matrix given as input. It employs a signature
#' derived from a dataset of hESC h9 cell line samples, which is ideal to
#' measure ADAR activity in a human neuronal context.
#'
#' @usage neuronalSignature(mat, object)
#' @param mat Either a numeric expression matrix with genes as rows and samples
#' as columns or an ExpressionSet from a Seurat object, SingleCellExperiment or
#' SpatialExperiment
#' @param object Character string indicating the type of object, either bulk or
#' eset
#' @return Matrix with ADAR activity scores
#' @examples
#'
#' # Load the regulon
#' data(neuronal_regulon)
#'
#' # Load the example expression matrix
#' data(example_matrix)
#'
#' # Compute ADAR activity
#' res <- neuronalSignature(example_matrix, object = "bulk")
#'
#' @importFrom viper viper
#' @importFrom GSVA ssgseaParam gsva
#' @importFrom GSEABase GeneSetCollection GeneSet
#' @importFrom SummarizedExperiment assay
#' @export

neuronalSignature <- function(mat, object = c("bulk", "eset")) {

  if (object == "bulk") {
    regulon <- neuronal_regulon
    res <- viper(mat, regulon, verbose = FALSE, minsize = 0)
    res["ADAR", ] <- (res["ADAR", ] - min(res["ADAR", ])) /
      (max(res["ADAR", ]) - min(res["ADAR", ]))
    res["ADARB1", ] <- (res["ADARB1", ] - min(res["ADARB1", ])) /
      (max(res["ADARB1", ]) - min(res["ADARB1", ]))
    res
  }

  else {

    if (class(mat) == "Seurat") {
      counts <- mat@assays$RNA$data
    }
    else {counts <- assay(mat, "logcounts")}

    neuronal_adar1_sign_genes <- names(neuronal_regulon$ADAR$tfmode)
    neuronal_adar1_gs <- GeneSet(neuronal_adar1_sign_genes,
                                 setName = "neuronal_adar1")

    ssgsea_obj <- ssgseaParam(counts, geneSets = neuronal_adar1_gs)
    ssgsea <- gsva(ssgsea_obj)
    ssgsea[1, ] <- (ssgsea[1, ] - min(ssgsea[1, ])) / (max(ssgsea[1, ]) -
                                                         min(ssgsea[1, ]))
    mat[["adar1_activity"]] <- ssgsea[1, ]

    return(mat)

  }

}

#' Compute ADAR activity with the mouse-specific neuronal signature
#'
#' The function returns ADAR activity scores for each sample/spot/cell depending
#' on the starting expression matrix given as input. It employs a signature
#' derived from a dataset of neurons and whole brain of mice, which is
#' ideal to measure ADAR activity in a mouse neuronal context.
#'
#' @usage neuronalMmSignature(mat, object)
#' @param mat Either a numeric expression matrix with genes as rows and samples
#' as columns or an ExpressionSet from a Seurat object, SingleCellExperiment or
#' SpatialExperiment
#' @param object Character string indicating the type of object, either bulk or
#' eset
#' @return Matrix with ADAR activity scores
#' @examples
#'
#' # Load the regulon
#' data(mouse_neuronal_regulon)
#'
#' # Load the example expression matrix
#' data(example_matrix_Mm)
#'
#' # Compute ADAR activity
#' res <- neuronalMmSignature(example_matrix_Mm, object = "bulk")
#'
#' @importFrom viper viper
#' @importFrom GSVA ssgseaParam gsva
#' @importFrom GSEABase GeneSetCollection GeneSet
#' @importFrom SummarizedExperiment assay
#' @export

neuronalMmSignature <- function(mat, object = c("bulk", "eset")) {

  if (object == "bulk") {
    regulon <- mouse_neuronal_regulon
    res <- viper(mat, regulon, verbose = FALSE, minsize = 0)
    res["Adar", ] <- (res["Adar", ] - min(res["Adar", ])) /
      (max(res["Adar", ]) - min(res["Adar", ]))
    res["Adarb1", ] <- (res["Adarb1", ] - min(res["Adarb1", ])) /
      (max(res["Adarb1", ]) - min(res["Adarb1", ]))
    res["Adarb2", ] <- (res["Adarb2", ] - min(res["Adarb2", ])) /
      (max(res["Adarb2", ]) - min(res["Adarb2", ]))
    res
  }

  else {

    if (class(mat) == "Seurat") {
      counts <- mat@assays$RNA$data
    }
    else {counts <- assay(mat, "logcounts")}

    Mm_adar1_sign_genes <- names(mouse_neuronal_regulon$Adar$tfmode)
    Mm_adar1_gs <- GeneSet(Mm_adar1_sign_genes,
                                setName = "neuronalMm_adar1")
    Mm_adar2_sign_genes <- names(mouse_neuronal_regulon$Adarb1$tfmode)
    Mm_adar2_gs <- GeneSet(Mm_adar2_sign_genes,
                                setName = "neuronalMm_adar2")

    all_sets <- GeneSetCollection(list(Mm_adar1_gs, Mm_adar2_gs))
    ssgsea_obj <- ssgseaParam(counts, geneSets = all_sets)
    ssgsea <- gsva(ssgsea_obj)
    ssgsea[1, ] <- (ssgsea[1, ] - min(ssgsea[1, ])) / (max(ssgsea[1, ]) -
                                                         min(ssgsea[1, ]))
    ssgsea[2, ] <- (ssgsea[2, ] - min(ssgsea[2, ])) / (max(ssgsea[2, ]) -
                                                         min(ssgsea[2, ]))
    mat[["adar1_activity"]] <- ssgsea[1, ]
    mat[["adar2_activity"]] <- ssgsea[2, ]

    return(mat)

  }

}
