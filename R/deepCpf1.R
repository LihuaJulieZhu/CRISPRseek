#' DeepCpf1 Algorithm for predicting CRISPR-Cpf1 gRNA Efficacy
#'
#' @description DeepCpf1 algorithm from https://doi.org/10.1038/nbt.4061, which
#' takes in 34 bp target sequences with/without chromatin 
#' accessibility information and returns predicted CRISPR-Cpf1 gRNA efficacy 
#' for each input sequence.
#'
#' @param in_df A dataframe formatted with three columns:
#' \itemize{
#' \item{gRNA} - gRNA Name
#' \item{sequence} - 34 bp target sequences as specified by 
#' http://deepcrispr.info/, i.e., 4bp before the 5' PAM, 4bp PAM, 20bp gRNA,
#' and 6bp after 3' of gRNA.
#' \item{chrom} - Optional binary variable indicating chromatin accessibility 
#' information with 1 indicating accessible and 0 not accessible.
#' }
#' @return a dataframe with two added columns:
#' \itemize{
#' \item{Seq-DeepCpf1 Score} - does NOT take into account chromatin 
#' accessibility information
#' \item{DeepCpf1 Score} - DOES take into account chromatin accessibility 
#' information if accessibility information is provided
#' }
#'
#' @details Having chromatin accessibility information will aid in the accuracy
#' of the scores, but one can still get accurate scoring with only the 
#' 34 bp target sequences.
#'
#' @author Paul Scemama and Lihua Julie Zhu
#' 
#' @references Kim et al., Deep learning improves prediction of CRISPR–Cpf1 
#' guide RNA activityNat Biotechnol 36, 239–241 (2018). 
#' https://doi.org/10.1038/nbt.4061
#' 
#' @import keras
#' @importFrom mltools one_hot
#' @importFrom dplyr select ends_with mutate_if
#' @importFrom data.table as.data.table
#' @importFrom stats predict
#'
#' @export
#'
#' @examples
#' library(keras)
#' library(mltools)
#' library(dplyr)
#' library(data.table)
#'
#' use_implementation("tensorflow")
#'
#' in_df = data.frame('gRNA' = paste("gRNA", 1:4), sequence =
#' c('GTTATTTGAGCAATGCCACTTAATAAACATGTAA',
#'  'TGACTTTGAATGGAGTCGTGAGCGCAAGAACGCT',
#'  'GTTATTTGAGCAATGCCACTTAATAAACATGTAA',
#'  'TGACTTTGAATGGAGTCGTGAGCGCAAGAACGCT'),
#'  chrom = c(0,1, 0, 1))
#'
#'  deepCpf1(in_df)


deepCpf1 <- function(in_df){
  # use tensorflow implementation
  use_implementation('tensorflow')

  ##### PREPROCESSING #####
  in_df$sequence <- toupper(in_df$sequence)
  SEQ <- as.data.table(strsplit(in_df$sequence,"")) %>%
    mutate_if(is.character, as.factor) %>%
    one_hot
  CA = in_df$chrom * 100

  SEQ <- array(c(t(select(as.data.frame(SEQ), ends_with("A"))),
                 t(select(as.data.frame(SEQ), ends_with("C"))),
                 t(select(as.data.frame(SEQ), ends_with("G"))),
                 t(select(as.data.frame(SEQ), ends_with("T")))),
               dim = c(nrow(in_df), 34L, 4L))
  
  CA <- array(as.numeric(unlist(CA)), dim = c(length(in_df$chrom), 1))

  #### MAKE MODELS

  #### SEQUENCE DEEP CPF1 (Disregards Chromatin Accessibility)

  # Define Input Tensor
  seq_input <- layer_input(shape = c(34,4))

  # Create output tensor
  seq_predictions <- seq_input %>%
    layer_conv_1d(filters = 80, kernel_size = 5, activation = 'relu') %>%
    layer_average_pooling_1d(pool_size = 2) %>%
    layer_flatten() %>%
    layer_dropout(rate = 0.3) %>%
    layer_dense(units = 80, activation = 'relu') %>%
    layer_dropout(rate = 0.3) %>%
    layer_dense(units = 40, activation = 'relu') %>%
    layer_dropout(rate = 0.3) %>%
    layer_dense(units = 40, activation = 'relu') %>%
    layer_dropout(rate = 0.3) %>%
    layer_dense(units = 1, activation = 'linear')

  # create model
  seq_model <- keras_model(inputs = seq_input, outputs = seq_predictions)

  ## get path for external data (weights for model)
  seq_weights_filepath <- system.file("extdata/DeepCpf1", 
                                     "real_weights/Seq_deepCpf1_weights.h5", 
                                     package = "CRISPRseek")

  # Load weights for the seq-model
  load_model_weights_hdf5(seq_model, seq_weights_filepath)

  # predict on test data
  seq_model.score = predict(seq_model, SEQ, batch_size = 50, verbose = 0)

  #### DeepCpf1 - REGARDS Chromatin Accessibility

  # define input tensors
  seq_input_portion <- layer_input(shape = c(34,4))
  ca_input_portion <- layer_input(shape = c(1,1))

  # create output tensor for sequential portion
  seq_output <- seq_input_portion %>%
    layer_conv_1d(filters = 80, kernel_size = 5, activation = 'relu',
                  input_shape = c(34,4)) %>%
    layer_average_pooling_1d(pool_size = 2) %>%
    layer_flatten() %>%
    layer_dropout(rate = 0.3) %>%
    layer_dense(units = 80, activation = 'relu') %>%
    layer_dropout(rate = 0.3) %>%
    layer_dense(units = 40, activation = 'relu') %>%
    layer_dropout(rate = 0.3) %>%
    layer_dense(units = 40, activation = 'relu')

  # create output tensor for ca portion
  ca_output <- ca_input_portion %>%
    layer_dense(units = 40, activation = 'relu')

  # create output tensor for final portion
  final_output <- layer_multiply(inputs = 
        c(seq_output, ca_output)) %>%  # piece-wise multiplication
    layer_dropout(rate = 0.3) %>%
    layer_dense(1, activation = 'linear')

  # create model
  seqCA_model <- keras_model(inputs = c(seq_input_portion, ca_input_portion), 
                             outputs = final_output)


  ## get path for external data (weights for model)
  CA_weights_filepath <-
    system.file("extdata/DeepCpf1", "real_weights/DeepCpf1_weights.h5", 
                package = "CRISPRseek")


  # load  weights
  load_model_weights_hdf5(seqCA_model, CA_weights_filepath)

  # predict
  seqCA_model.score = predict(seqCA_model, 
                              x = list(SEQ,CA), 
                              batch_size = 50, 
                              verbose = 0) * 3


  
  # add columns with scores to the input dataframe
  in_df$seq.DeepScore = seq_model.score
  in_df$DeepScore = seqCA_model.score

  in_df
}
