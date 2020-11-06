#' Define cells that are significantly expressing mRNA using negative control information
#' @param input input data.table (wide-format) with marker expression values and cellID
#' @param cellID column name for unique cell identifier
#' @param mRNA_channels character vector with mRNA channel names of interest
#' @param negative_control channel name of negative control
#' @param threshold (optional) numeric threshold for adjusted p-value cut-off - (default = 0.05)
#' @param return_calc_metrics (optional) returns additional data as data.table (p-values, adjusted p-values...), data.table and SCE object are stored in a list (default = FALSE)
#' @return returns SCE object with new colData entires for each mRNA_channels entry (1 = significant expression, 0 = non-significant)
#' @export

compute_difference <- function(input_sce,
                               cellID = NULL,
                               assay_name = NULL,
                               threshold = 0.05,
                               mRNA_channels,
                               negative_control,
                               return_calc_metrics = FALSE){
  # input validity control
  .checkInputmRNADetection(input_sce, cellID, assay_name, threshold, mRNA_channels, negative_control)

  # create data.table with assay values
  dat_wide <- as.data.table(t(assay(input_sce, assay_name)))

  # add cellID
  dat_wide[,id := input_sce[[cellID]]]

  # create table for binary output
  col_to_keep = c("id", mRNA_channels)
  binary_output = dat_wide[,..col_to_keep]

  # create empty list for output metrics
  output_list <- list()

  # loop through all channels of interest and calculate difference between channel and negative control (and compute p-values)
  for (i in mRNA_channels) {

    # calculate difference between signal and negative control and compute stats (pvalue, adjusted)
    diff = as.data.table(dat_wide[[as.character(i)]] - dat_wide[[as.character(negative_control)]])
    diff[, scaled_diff := scale(V1)]
    diff[, mean_expresson_chemokine := dat_wide[[as.character(i)]]]
    diff[, mean_DapB := dat_wide[[negative_control]]]
    diff[, p_value := pnorm(-abs(scaled_diff), mean = 0 , sd = 1)]
    diff[, adj_p_value := p.adjust(as.numeric(p_value), method = "BH")]
    colnames(diff) = c("diff", "scaled_diff", "mean_chemokine", "mean_negative_control", "p_value", "padj")

    # assign binary class (p-value below threshold, scaled difference bigger than 0)
    binary_output[,as.character(i) := ifelse(diff$padj <= threshold &
                                               diff$scaled_diff > 0,
                                          1, 0)]

    # add info to SCE object (new colData entry)
    colData(input_sce)$cur_channel <- binary_output[[as.character(i)]]

    # overwrite column name with channel name
    names(colData(input_sce))[length(names(colData(input_sce)))] <- i

    # append diff matrix to output list if return_calc_metrics is TRUE
    if(return_calc_metrics == TRUE){
      output_list[[i]] <- diff
    }
  }

  # return options
  if(return_calc_metrics == FALSE){
    return(input_sce)
  } else{
    output_list[["output_sce"]] <- input_sce
    return(output_list)
  }
}
