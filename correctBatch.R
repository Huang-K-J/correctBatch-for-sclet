combine_VariableFeatures<- function (all.batch, include_fields = NULL,equiweight = TRUE, ncells = NULL) 
{ 
  sub.dec <- vector("list", length(all.batch))
  for (i in seq_along(all.batch)) {
    sub.dec[[i]] <- rowData(all.batch[[i]])
  }
  collected <- scuttle:::.unpackLists(sub.dec)
  if (is.null(ncells)) {
    ncells <- rep(10L, length(collected))
  }
  if (is.null(include_fields)) {
    include_fields <- exclude_some_fields(collected, c("ID", "Symbol"))
  }
  combine_blocked(collected,geometric = FALSE,equiweight = equiweight, 
                  ncells = ncells,  ave.fields = include_fields)
}


exclude_some_fields <- function (collected,exclude) 
{
  all.numerics <- list()
  for (i in seq_along(collected)) {
    x <- collected[[i]]
    all.numerics[[i]] <- colnames(x)[vapply(x, is.numeric, 
                                            TRUE)]
  }
  setdiff(Reduce(intersect, all.numerics), exclude)
}



combine_blocked <- function (blocks, geometric = FALSE, equiweight=TRUE, ncells=ncells,
                             ave.fields =c("mean", "variance", "variance.expected", "variance.standardized"))
{ 
  valid = ncells >= 2L
  weights = ncells
  if (length(blocks) == 1L) {
    return(blocks[[1]])
  }
  rn <- unique(lapply(blocks, rownames))
  if (length(rn) != 1L) {
    stop("gene should be the same")
  }
  if (equiweight) {
    weights <- rep(1, length(blocks))
  }
  else if (is.null(weights)) {
    stop("'weights must be specified")
  }
  original <- blocks
  if (!any(valid)) {
    stop("no entry of 'blocks' has positive weights")
  }
  blocks <- blocks[valid]
  weights <- weights[valid]
  combined <- list()
  for (i in ave.fields) {
    extracted <- lapply(blocks, "[[", i = i)
    extracted <- mapply("*", extracted, weights, SIMPLIFY = FALSE, 
                        USE.NAMES = FALSE)
    averaged <- Reduce("+", extracted)/sum(weights)
    combined[[i]] <- averaged
  }
  output <- DataFrame(combined, row.names = rn[[1]])
  output$per.block <- do.call(DataFrame, c(lapply(original, 
                                                  I), list(check.names = FALSE)))
  output
} 


BatchRemover <- function (..., batch = NULL, restrict = NULL, correct.all = TRUE, 
          assay.type = "counts", PARAM = FastMnnParam(), multi.norm.args = list(), 
           model.var.args = list(),hvg.args=5000) 
{ 
  all.batches <- scuttle:::.unpackLists(...)
  all.genes <- lapply(all.batches, rownames)
  universe <- Reduce(intersect, all.genes)
  if (length(universe) == 0) {
    stop("no genes remaining in the intersection")
  }
  non.same <- FALSE
  for (i in seq_along(all.batches)) {
    if (!identical(universe, rownames(all.batches[[i]]))) {
      non.same <- TRUE
      all.batches[[i]] <- all.batches[[i]][universe, , 
                                           drop = FALSE]
    }
  }
  
    all.batches <- do.call(multiBatchNorm, c(all.batches, list(batch = batch, 
                 assay.type = assay.type, preserve.single = TRUE), multi.norm.args))
   
    dec <- combine_VariableFeatures(all.batches,include_fields = NULL,equiweight = TRUE, ncells = NULL) 
    i <- order(dec$variance.standardized, decreasing = TRUE)
    dec <- dec[i, ]
    hvgs <- rownames(dec)[1:hvg.args]
    corrected <- batchelor:::batchCorrect(all.batches, batch = batch, restrict = restrict, 
                              correct.all = correct.all, subset.row = hvgs, PARAM = PARAM)
    list( hvgs = hvgs, corrected = corrected)
  } 

















