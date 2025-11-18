#' @title
#' DNA Sequence Encoding Function
#'
#' @description
#' This function encodes DNA strings into a data frame of factor columns, where each column represents a nucleotide position.
#'
#' @param dna_strings A character vector of DNA strings, all of the same length.
#' @return A data frame where each column is a factor (levels: "A", "T", "C", "G") representing a nucleotide position in the input DNA strings.
#' @examples
#' library(m6APrediction)
#' @import stats
dna_encoding <- function(dna_strings){
  nn <- nchar( dna_strings[1] )
  seq_m <- matrix( unlist( strsplit(dna_strings, "") ), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}

#' Multiple Samples Prediction for m6A Sites
#'
#' This function predicts m6A status and probability for multiple samples using a trained machine learning model.
#'
#' @param ml_fit A trained machine learning model (e.g., from randomForest).
#' @param feature_df A data frame containing features for multiple samples. Required columns: "gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction", "evolutionary_conservation", "DNA_5mer".
#' @param positive_threshold A numeric value (default 0.5) to determine if a prediction is "Positive" (if probability exceeds this threshold).
#' @return The input feature data frame with two additional columns: "predicted_m6A_prob" (predicted probability of being m6A positive) and "predicted_m6A_status" (factor with levels "Negative", "Positive").
#' @examples
#'  ml_fit=readRDS(system.file("extdata", "rf_fit.rds", package="m6APrediction"))
#'  feature_df=read.csv(system.file("extdata", "m6A_input_example.csv",package="m6APrediction"))
#'  prediction_multiple(ml_fit, feature_df)
#' @import randomForest
#' @export
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length",
                  "distance_to_junction", "evolutionary_conservation", "DNA_5mer")
                %in% colnames(feature_df)))

  feature_df$RNA_type   <- factor(feature_df$RNA_type,
                                  levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region,
                                  levels = c("CDS", "intron", "3'UTR", "5'UTR"))

  seq_df <- dna_encoding(feature_df$DNA_5mer)

  input_df <- cbind(feature_df[, c("gc_content", "RNA_type", "RNA_region",
                                   "exon_length", "distance_to_junction",
                                   "evolutionary_conservation")],
                    seq_df)

  prob <- predict(ml_fit, newdata = input_df, type = "prob")[, "Positive"]
  status <- ifelse(prob > positive_threshold, "Positive", "Negative")

  feature_df$predicted_m6A_prob <- prob
  feature_df$predicted_m6A_status <- factor(status, levels = c("Negative", "Positive"))

  return(feature_df)
}

#' Single Sample Prediction for m6A Sites
#'
#' This function takes individual feature values and u#'m6A probability and status for a single sample. It returns the predicted values as a named vector.
#'
#' @param ml_fit A trained machine learning model (e.g., from randomForest).
#' @param gc_content A numeric value for the GC content.
#' @param RNA_type A character string specifying the RNA type (one of "mRNA", "lincRNA", "lncRNA", "pseudogene").
#' @param RNA_region A character string specifying the RNA region (one of "CDS", "intron", "3'UTR", "5'UTR").
#' @param exon_length A numeric value for the exon length.
#' @param distance_to_junction A numeric value for the distance to junction.
#' @param evolutionary_conservation A numeric value for evolutionary conservation.
#' @param DNA_5mer A character string representing the 5-mer DNA sequence.
#' @param positive_threshold A numeric value (default 0.5) to determine if a prediction is "Positive".
#' @return A named vector with two elements: "predicted_m6A_prob" (predicted probability) and "predicted_m6A_status" (either "Positive" or "Negative").
#' @examples
#' # Assume ml_fit is a trained model
#' # prediction_single(ml_fit, gc_content = 0.6, RNA_type = "mRNA", RNA_region = "CDS",
#' #                   exon_length = 12, distance_to_junction = 100,
#' #                   evolutionary_conservation = 0.8, DNA_5mer = "ATCGT")
#' @import randomForest
#' @export
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length, distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){

  feature_df <- data.frame(
    gc_content, RNA_type, RNA_region, exon_length,
    distance_to_junction, evolutionary_conservation, DNA_5mer,
    stringsAsFactors = FALSE
  )

  feature_df$RNA_type   <- factor(feature_df$RNA_type,
                                  levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region,
                                  levels = c("CDS", "intron", "3'UTR", "5'UTR"))

  seq_df <- dna_encoding(feature_df$DNA_5mer)
  seq_df[] <- lapply(seq_df, function(x) factor(x, levels = c("A","T","C","G")))

  input_df <- cbind(
    feature_df[, c("gc_content","RNA_type","RNA_region",
                   "exon_length","distance_to_junction","evolutionary_conservation")],
    seq_df
  )

  prob <- predict(ml_fit, newdata = input_df, type = "prob")[, "Positive"]
  status <- ifelse(prob > positive_threshold, "Positive", "Negative")

  returned_vector <- c("predicted_m6A_prob" = prob,
                       "predicted_m6A_status" = status)
  return(returned_vector)
}

