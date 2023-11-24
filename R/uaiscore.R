






#' The function uaiscore is designed for predictive modeling of gene expression data.
#'
#' @param eset Required matrix or data frame containing the gene expression data to be used for prediction. If not a matrix, it will be coerced into one.
#' @param pdata An optional data frame containing sample data to be used during gene expression prediction. Defaults to NULL.
#' @param id_pdata An optional character vector indicating the name of the sample ID column to be used in pdata data frame. Defaults to "ID".
#' @param scale An optional logical value indicating whether to scale the gene expression data. Default value is FALSE.
#' @param check_eset An optional logical value indicating whether to remove outlier genes from the gene expression data. Default value is FALSE.
#'
#' @return
#' @export
#' @author Dongqiang Zeng
#' @examples
uaiscore<-function(eset, pdata = NULL, id_pdata = "ID", scale = FALSE, check_eset = FALSE){

  if(!is.matrix(eset)) eset<-as.matrix(eset)
  ###########################################

  ############################################
  if(scale){
    cat(crayon::green(">>>-- Scaling data...\n"))
    eset<-t(scale(t(eset)))
  }
  #############################################
  if(check_eset){
    cat(crayon::green(">>>-- Removing outlier genes...\n"))
    genes<- rownames(eset)
    genes<- IOBR::feature_manipulation(data = eset, feature = genes, is_matrix = TRUE, print_result = T)
    eset<-eset[rownames(eset)%in%genes, ]
  }
  #############################################
  data("uai_feas", package = "uaiscore")
  ###########################################
  ###########################################
  data("uai_model", package = "uaiscore")
  ############################################
  # https://github.com/r-lib/crayon
  cat(crayon::green(">>>-- Predicting new data with uaiscore model...\n"))
  res<-      predict_cox_model(sur_model    = uai_model,
                                  eset_new     = eset,
                                  pdata_new    = pdata,
                                  feas         = NULL,
                                  feas_genes   = uai_feas,
                                  prefix       = c("\\-","_"),
                                  id_pdata_new = id_pdata)

  cat(crayon::green(">>>-- DONE! \n"))

  return(res)

}
