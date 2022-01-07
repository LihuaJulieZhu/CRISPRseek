#' DeepCpf1 Algorithm for predicting CRISPR-Cpf1 gRNA Efficacy
#'
#' @description DeepCpf1 algorithm from https://doi.org/10.1038/nbt.4061, which
#' takes in 34 bp target sequences with/without chromatin
#' accessibility information and returns predicted CRISPR-Cpf1 gRNA efficacy
#' for each input sequence.
#'
#' @param extendedSequence Sequences containing gRNA plus PAM plus flanking
#' sequences. Each sequence should be 34 bp long as specified by
#' http://deepcrispr.info/, i.e., 4bp before the 5' PAM, 4bp PAM, 20bp gRNA,
#' and 6bp after 3' of gRNA.
#' @param chrom_acc Optional binary variable indicating chromatin accessibility
#' information with 1 indicating accessible and 0 not accessible.

#' @return a numeric vector with prediced CRISPR-Cpf1 gRNA efficacy
#' taking into account chromatin accessibility
#' information if accessibility information is provided
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
#' @importFrom dplyr select ends_with mutate_if slice
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
#' extendedSequence <- c('GTTATTTGAGCAATGCCACTTAATAAACATGTAA',
#'  'TGACTTTGAATGGAGTCGTGAGCGCAAGAACGCT',
#'  'GTTATTTGAGCAATGCCACTTAATAAACATGTAA',
#'  'TGACTTTGAATGGAGTCGTGAGCGCAAGAACGCT')
#' chrom_acc <- c(0,1, 0, 1)
#'
#' if (interactive()) {
#'  deepCpf1(extendedSequence = extendedSequence, chrom_acc = chrom_acc)
#' }


deepCpf1 <- function(extendedSequence, chrom_acc){
  # use tensorflow implementation
  use_implementation('tensorflow')
  if (missing(extendedSequence))
     stop("extendedSequence is required for predicting efficacy using DeepCpf1 algorithm!")
  len <- unlist(lapply(extendedSequence, nchar))
  if (sum(len == 34, na.rm = TRUE) == 0)
  {
      warning("None of the extendedSequences has length of 34 which is required for DeepCpf1 algorithm!")
      return(rep(NA, length(extendedSequence)))
  }
  ##### PREPROCESSING #####
  sequence <- toupper(extendedSequence)

  if (!missing(chrom_acc))
  {
      in_df <- as.data.frame(cbind(sequence = sequence, chrom_acc = chrom_acc, len = len))
      in_df <- subset(in_df, in_df[, 3] == 34)
      CA = as.numeric(in_df[,2]) * 100
      CA <- array(as.numeric(CA), dim = c(length(CA), 1))
  }
  else
  {
      in_df <- as.data.frame(cbind(sequence = sequence, len = len))
      in_df <- subset(in_df, in_df[, 2] ==  34)
  }

  SEQ <- as.data.table(strsplit(paste0("ACGT", in_df[,1]),"")) %>%
    mutate_if(is.character, as.factor) %>%
    one_hot(dropUnusedLevels = FALSE, dropCols = FALSE) %>%
    slice(-c(1:4))


  SEQ <- array(c(t(select(as.data.frame(SEQ), ends_with("A"))),
                 t(select(as.data.frame(SEQ), ends_with("C"))),
                 t(select(as.data.frame(SEQ), ends_with("G"))),
                 t(select(as.data.frame(SEQ), ends_with("T")))),
               dim = c(nrow(in_df), 34L, 4L))


  #### MAKE MODELS
  # Define Input Tensor
  if (missing(chrom_acc))
  {
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
     effi = predict(seq_model, SEQ, batch_size = 50, verbose = 0)
  }
  else
  {
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
  effi = predict(seqCA_model,
                              x = list(SEQ,CA),
                              batch_size = 50,
                              verbose = 0) * 3


  }
   if (missing(chrom_acc))
   {
     temp <- cbind(extendedSequence = in_df[,1], effi = effi)
     as.numeric(temp[match(sequence, temp[,1]),2]) /100
   }
   else
   {
      temp <- cbind(extendedSequence = paste0(in_df[,1],in_df[,2]), effi = effi)
      as.numeric(temp[match(paste0(sequence,chrom_acc), temp[,1]),2]) /100
   }
}
