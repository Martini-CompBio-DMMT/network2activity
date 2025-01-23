#' Get the data frame from the network
#'
#' This function is used to get a data frame from the file obtained by the
#' network reconstruction with ARACNe.
#'
#' @usage getDf(net)
#' @param net network obtained from ARACNe
#' @return data frame
#' @examples
#'
#' # Load the example network
#' data(network)
#'
#' # Get the data frame from the network obtained by ARACNe
#' net_df <- getDf(network)
#'
#' @importFrom igraph graph_from_adjacency_matrix as_data_frame
#' @export

getDf <- function(net) {

  g  <- graph_from_adjacency_matrix(net, weighted=TRUE, mode="undirected")
  dataf <- as_data_frame(g)

  dataf

}


#' Get a list of genes with their interactors
#'
#' This function extracts gets the list of regulators, and it adds the genes
#' for which we want to compute the activity with VIPER.
#' Then it filters the data frame obtained by the network to consider only the
#' regulators of the gene of interest and, for each regulator, it gets a list of
#' their interactors and the corresponding scores.
#'
#' @usage getInteractorsList(dataf, reg_list)
#' @param dataf data frame obtained by the network
#' @param reg_list list with the regulator genes
#' @return list
#' @examples
#'
#' # Load the list of regulators
#' data(regulators_list)
#'
#' # Load the example data frame derived from the network
#' data(network_df)
#'
#' # Get the list of the regulators with their interactors
#' interactors <- getInteractorsList(network_df, regulators_list)
#'
#' @importFrom dplyr filter
#' @export

getInteractorsList <- function(dataf, reg_list) {

  genes_selected <- reg_list

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

  df_list

}

#' Get the adjacency matrix file as requested by VIPER
#'
#' This function is used to obtain the adjacency matrix file in the format
#' requested by VIPER to obtain the regulon with the "aracne2regulon" function.
#'
#' @usage writeAdjFile(df_list, file_name = "")
#' @param df_list list of genes with their interactors
#' @param file_name how to nominate the file
#' @return .adj file
#' @examples
#'
#' # Load the list of regulators
#' data(regulators_list)
#'
#' # Load the example data frame derived from the network
#' data(network_df)
#'
#' # Get the list of the regulators with their interactors
#' interactors <- getInteractorsList(network_df, regulators_list)
#'
#' writeAdjFile(interactors, file_name = "interactors")
#'
#' @import utils
#' @export

writeAdjFile <- function(df_list, file_name = "") {

  file_txt <- file.path(getwd(), paste0(file_name, ".txt"))
  file_adj <- file.path(getwd(), paste0(file_name, ".adj"))

  if (file.exists(file_txt)) {
    file.remove(file_txt)
  }

  if (file.exists(file_adj)) {
    file.remove(file_adj)
  }

  i = 1
  file(file_txt, "w")
  while (i < length(df_list)) {
    a <- c(names(df_list[i]), paste(df_list[[i]]$to, df_list[[i]]$weight))
    b <- write.table(x = rbind(a), file = file_txt, row.names = F,
                     col.names = F, quote = F, append = T)
    i = i+1
  }

  txt_file <- readLines(file_txt)
  adj_file <- file_adj
  modified_file <- gsub(" ", "\t", txt_file)
  writeLines(modified_file, adj_file)

  file_adj

}

#' Infer ADAR activity with the dataset-specific signature
#'
#' This function is used to infer the protein activity with VIPER using the
#' signature computed from the original dataset.
#'
#' @usage inferenceWithYourSignature(adjfile, mat)
#' @param adjfile file .adj with the interactors and their scores
#' @param mat expression matrix
#' @return VIPER results
#' @examples
#'
#' # Load the list of regulators
#' data(regulators_list)
#'
#' # Load the example data frame derived from the network
#' data(network_df)
#'
#' # Load the example expression matrix
#' data(vsd_matrix)
#'
#' # Get the list of the regulators with their interactors
#' interactors <- getInteractorsList(network_df, regulators_list)
#'
#' # Get the file in the format requested by VIPER
#' interactors_adj <- writeAdjFile(interactors, file_name = "interactors")
#'
#' # Infer protein activity with VIPER
#' inferenceWithYourSignature(interactors_adj, vsd_matrix)
#'
#' @import viper
#' @export

inferenceWithYourSignature <- function(adjfile, mat) {

  regul <- aracne2regulon(adjfile, mat, verbose = FALSE)
  vpres <- viper(mat, regul, verbose = FALSE)
  vpres

}

#' Infer ADAR activity with the A549 cells signature
#'
#' This function is used to infer the protein activity with VIPER using the
#' signature computed from the A549 cell line.
#'
#' @usage inferenceWithA549Signature(mat)
#' @param mat expression matrix
#' @return VIPER results
#' @examples
#'
#' # Load the regulon
#' data(regul_A549_cells)
#'
#' # Load the example expression matrix
#' data(vsd_matrix)
#'
#' # Infer protein activity with VIPER
#' inferenceWithA549Signature(vsd_matrix)
#'
#' @import viper
#' @export

inferenceWithA549Signature <- function(mat) {

  regul <- regul_A549_cells
  vpres <- viper(mat, regul, verbose = FALSE, minsize = 5)
  vpres

}

#' Infer ADAR activity with the neurons-specific signature
#'
#' This function is used to infer the protein activity with VIPER using the
#' signature computed from hESC h9 neurons.
#'
#' @usage inferenceWithNeuronsSignature(mat)
#' @param mat expression matrix
#' @return VIPER results
#' @examples
#'
#' # Load the regulon
#' data(regul_hESC_neurons)
#'
#' # Load the example expression matrix
#' data(vsd_matrix)
#'
#' # Infer protein activity with VIPER
#' inferenceWithNeuronsSignature(vsd_matrix)
#'
#' @import viper
#' @export

inferenceWithNeuronsSignature <- function(mat) {

  regul <- regul_hESC_neurons
  vpres <- viper(mat, regul, verbose = FALSE, minsize = 5)
  vpres

}

#' Infer ADAR activity with the interferon-specific signature
#'
#' This function is used to infer the protein activity with VIPER using the
#' signature computed from HEK293T cells treated with interferon.
#'
#' @usage inferenceWithInterferonSignature(mat)
#' @param mat expression matrix
#' @return VIPER results
#' @examples
#'
#' # Load the regulon
#' data(regul_interferon_specific)
#'
#' # Load the example expression matrix
#' data(vsd_matrix)
#'
#' # Infer protein activity with VIPER
#' inferenceWithInterferonSignature(vsd_matrix)
#'
#' @import viper
#' @export

inferenceWithInterferonSignature <- function(mat) {

  regul <- regul_interferon_specific
  vpres <- viper(mat, regul, verbose = FALSE, minsize = 5)
  vpres

}
