#' Removes the effect of the cell-cycle
#'
#' \code{ccRemover} returns a data matrix with the effects of the cell-cycle
#' removed.
#'
#' Implements the algorithm described in Barron, M. & Li, J.
#' "Identifying and removing the cell-cycle effect from scRNA-Sequencing data"
#' (2016), Scientific Reports. This function takes a normalized,
#' log-transformed and centered matrix of scRNA-seq  expression data
#' and a list of genes which are known to be related to the cell-cycle effect.
#' It then captures the main sources of variation in the data and determines
#' which of these are related to the cell-cycle before removing those that are.
#' Please see the original manuscript for further details.
#'
#' @param dat A list containing a data frame , \code{x}, that contains gene expression
#' measurements with each column representing a sample and each row
#' representing a gene and a logical vector, \code{if_cc}, that indicates
#' which of the genes/rows are related to the cell-cycle or factor of interest.
#'
#' It is recommended that the elements of x are log-transformed and centered
#' for each gene. For example if \code{x} contains TPM measurements then we
#' suggest the following two-steps:
#' Step 1: \code{dat$x <- log(dat$x + 1)}
#' Step 2: \code{dat$x} - rowMeans(dat$x)
#' ccRemover requires that the samples have been properly normalized for
#' sequencing depth and we recommend doing so prior to applying the above steps.
#'
#' The \code{if_cc} vector must be the same length as the number of rows in
#' \code{x} and have elements equal to \code{TRUE} for genes which are related
#' to the cell-cycle and and elements equal to \code{FALSE} for genes which
#' are unrelated to the cell-cycle.
#' @param cutoff The significance cutoff for identifying sources of variation
#' related to the cell-cycle. The default value is 3, which roughly corresponds
#' to a p-value of 0.01.
#' @param max_it The maximum number of iterations for the algorithm. The
#' default value is 4.
#' @param nboot The number of bootstrap repititions to be carried out on each
#' iteration to determine the significance of each component.
#' @param ntop The number of components considered tested at each iteration as
#' cell-cycle effects. The default value if 10
#' @param bar Whether to display a progress bar or not. The progress bar will
#' not work in R-markdown enviornments so this option may be turned off. The
#' default value is \code{TRUE}.
#'
#' @return A data matrix with the effects of the cell-cycle removed.
#'
#' @examples
#' set.seed(10)
#' # Load in example data
#' data(t.cell_data)
#' head(t.cell_data[,1:5])
#' # Center data and select small sample for example
#' t_cell_data_cen <- t(scale(t(t.cell_data[,1:20]), center=TRUE, scale=FALSE))
#' # Extract gene names
#' gene_names <- rownames(t_cell_data_cen)
#' # Determine which genes are annotated to the cell-cycle
#' cell_cycle_gene_indices <- gene_indexer(gene_names,
#' species = "mouse", name_type = "symbol")
#' # Create "if_cc" vector
#' if_cc <- rep(FALSE,nrow(t_cell_data_cen))
#' if_cc[cell_cycle_gene_indices] <- TRUE
#' # Move data into list
#' dat <- list(x=t_cell_data_cen, if_cc=if_cc)
#' # Run ccRemover
#' \dontrun{
#'  xhat <- ccRemover(dat, cutoff = 3, max_it = 4, nboot = 200, ntop = 10)
#' }
#' # Run ccRemover with reduced bootstrap repetitions for example only
#' xhat <- ccRemover(dat, cutoff = 3, max_it = 4, nboot = 20, ntop = 10)
#' head(xhat[,1:5])
#' # Run ccRemover with more compoents considered
#' \dontrun{
#' xhat <- ccRemover(dat, cutoff = 3, max_it = 4, nboot = 200, ntop = 15)
#'  }
#'
ccRemover <- function(dat, cutoff=3, max_it=4, nboot=200, ntop=10, bar=TRUE, cores=48,  rank.=6000)
  #Li: one can add more options to stats::prcomp, e.g. rank.=3, for choosing different PC nr to keep.
{
  ## check arguments
  if (!is.list(dat)) stop("dat has to be a list!")
  if (is.null(dat$x)) stop("dat has to have an element x!")
  if (is.null(dat$if_cc)) stop("dat has to have an element if.cc!")

  if (!is.vector(dat$if_cc)) stop("dat$if_cc has to be a vector!")
  if (!is.logical(dat$if_cc)) stop("dat$if_cc has to be a vector of
                                   TRUE or FALSE!")

  if (!is.matrix(dat$x) & !is.data.frame(dat$x)) stop("dat$x has to
                                                      be a matrix or data frame")
  if (nrow(dat$x) != length(dat$if_cc)) stop("The number of rows of dat$x has
                                             to be the same as the length of
                                             dat$if.cc!")

  if (sum(dat$x < 0) == 0) warning("All the values of dat$x are
                                   nonnegative--did you perform log(x) or
                                   log(x+1) transformation? This transformation
                                   is ususally necessary.")
  if (sum(rowSums(dat$x) > 1 | rowSums(dat$x) < -1) > 0) warning("It is highly
                                                                 recommended
                                                                 that the rows
                                                                 of dat$x are
                                                                 centered.")
  if (sum(dat$if_cc) < 1) stop("No genes are cell-cycle genes! Please check
                               your dat$if_cc vector")
  if(sum(!dat$if_cc) < 1) stop("All genes are cell-cycle genes! Please check
                               your dat$if_cc vector")
  percentage_cc <- sum(dat$if_cc)/length(dat$if_cc)
  cat(percentage_cc, " of genes are cell-cycle genes")
  if(percentage_cc <= 0.01) warning("Few genes are cell-cycle genes. This may
                                    lead to unreliable estimation of the
                                    cell-cycle effect. Are you sure your
                                    dat$if_cc vector is correct?")
  if(percentage_cc >= 0.99) warning("Almost all genes are cell-cycle genes.
                                    This may lead to unreliable estimation of
                                    the cell-cycle effect. Are you sure your
                                    dat$if_cc vector is correct?")
  ## begin calculation
  xy <- dat$x[dat$if_cc, ]
  xn <- dat$x[!dat$if_cc, ]

  for (i in 1 : max_it)
  {
    cat("\nIteration ", i, "...", fill=TRUE)

    ## calculate the test statistic and the t statistic
    res_boot <- bootstrap_diff(xy=xy, xn=xn, nboot=nboot, cores=cores, bar, rank.=rank.)
    cat("The bootstrap results on the top", ntop, "components are:")
    print(res_boot[1 : ntop, ])
#          xn_load   xy_load   diff_load t_load_boot #t_load_boot=diff_load/sd
#   PC1  9.681281 11.164245  1.48296444  3.74402433
#   PC2  6.932106  6.221538 -0.71056800 -0.89587019
#   PC3  5.043708  5.091521  0.04781354  0.08173771
#   PC4  4.212430  5.938865  1.72643456  5.74015903
#   PC5  3.019651  2.733409 -0.28624251 -0.45537545
#   PC6  2.875084  3.217810  0.34272615  0.42149293
#   PC7  2.682083  2.371485 -0.31059809 -1.60642301
#   PC8  2.585086  2.634081  0.04899532  0.11869654
#   PC9  2.465026  2.381195 -0.08383108 -0.76575403
#   PC10 2.406141  2.021187 -0.38495483 -2.60754835
    ## decide which components to remove
    which_cc <- which(res_boot$t_load_boot[1 : ntop] >= cutoff)
      # 1, 4
    ## when there is nothing to remove, end the iteration
    if (length(which_cc) == 0)
    {
      cat("No more cell-cycle effect is detected.", fill=TRUE)
      break
    }

    ## remove the cell-cycle effect
    cat("The follow components are removed:", which_cc, fill=TRUE)
    xn_pca <- prcomp(xn, scale.=FALSE) #Li remove stats::
    xn_hat <- xn # row: genes; col: cells
    xy_hat <- xy
    for (i in 1 : length(which_cc))
    {
      xn_hat <- xn_hat - (xn %*% xn_pca$rotation[, which_cc[i]]) %*%
        #xn(row-gene, col-cell) %*% xn_pca$rotation[, which_cc[i]] (row-cell, col-PC), output: row, genes; col, PCs.
        t(xn_pca$rotation[, which_cc[i]])#row: PCs, col: cells
        #The output of the second half of the above equation: row-genes, col-cells
        #The second half of the above equation: first load cc effect (cc dominating PCs) onto genes,
                #Then load the cc effect onto cells
        #The entire equation above: the original gene expression subtracted by the cc effect loading on each gene/cell.
      xy_hat <- xy_hat - (xy %*% xn_pca$rotation[, which_cc[i]]) %*%
        t(xn_pca$rotation[, which_cc[i]])
    }
    xn <- xn_hat
    xy <- xy_hat
  }

  xhat <- dat$x
  xhat[dat$if_cc, ] <- xy
  xhat[!dat$if_cc, ] <- xn

  return(xhat)
}


#' Calculates the difference in the loading score for cell-cycle and control
#' genes
#'
#' This function is only used internally inside ccRemover. The function
#' calcualtes the average load difference on the cell-cycle and control genes.
#' Bootstrap resampling is then used to provide a score for each component.
#' Please see the original manuscript for the mathematical details.
#'
#' @param xy The data for the genes which are annotated to the cell-cycle, i.e.
#' those genes for which "if_cc" is \code{TRUE}.
#' @param xn The data for the genes which are not annotated to the cell-cycle,
#' control genes, genes for which "if_cc" is \code{FALSE}
#' @param nboot The number of bootstrap repititions to be carried out on each
#' iteration to determine the significance of each component.
#' @param bar Whether to display a progress bar or not. The progress bar will
#' not work in R-markdown enviornments so this option may be turned off. The
#' default value is \code{TRUE}.
#'
#'
#' @return A data frame containing the loadings for each component on the
#' cell-cycle and control genes as well as the difference between the loadings
#' and the bootstrapped statistic for each loading.
#'
bootstrap_diff <- function(xy, xn, nboot=200, cores=48, bar=TRUE, rank.=6000)
{
  res0 <- get_diff(xy, xn, rank.=rank.)
  val <- length(res0$xn_load) #min(ncol(xn), nrow(xn))
  #diff_load_boot <- matrix(NA, nrow=val, ncol=nboot)
  cat("Bootstrapping...")
  if (bar == TRUE){
    pb <- utils::txtProgressBar(min = 1, max = nboot, style = 3)
  }
  
  
  
  
  bootRun <- function(i){
    if (bar == TRUE){
      utils::setTxtProgressBar(pb, i)
    }
    res <- get_diff(xy[sample(1 : nrow(xy), nrow(xy), replace=TRUE), ], 
                    #replace=T make some genes disappear, while others overrepresented.
                    xn[sample(1 : nrow(xn), nrow(xn), replace=TRUE), ], rank.=rank.)
                    #Bug: the random sampling may produce cells with identical gene expression profile,
                        #Leading to a 
    #diff_load_boot[, i] <- res$diff_load 
    res$diff_load
    #From each bootstrap run, we will get $nPC diff_load values
    #All together, we will get a matrix: nrow= nPC, ncol=nboot
    if(bar == TRUE){
    close(pb)}
  }
  

  diff_load_boot <- mclapply(1:nboot, bootRun, mc.cores = cores)
  
  
  #for (i in 1 : nboot)
  #{
  #  if (bar == TRUE){
  #    utils::setTxtProgressBar(pb, i)
  #  }
  #  res <- get_diff(xy[sample(1 : nrow(xy), nrow(xy), replace=TRUE), ], 
  #                  #replace=T make some genes disappear, while others overrepresented.
  #                  xn[sample(1 : nrow(xn), nrow(xn), replace=TRUE), ], ...)
  #                  #Bug: the random sampling may produce cells with identical gene expression profile,
  #                      #Leading to a 
  #  diff_load_boot[, i] <- res$diff_load 
  #  #From each bootstrap run, we will get $nPC diff_load values
  #  #All together, we will get a matrix: nrow= nPC, ncol=nboot
  #}
  
  sd1 <- mclapply(diff_load_boot, stats::sd, mc.cores = cores) %>% Reduce(c,.)
  #sd1 <- apply(diff_load_boot, 1, stats::sd)#1: row/PC #the diff_load sd at each PC among all bootstrap runs
  
  res0$t_load_boot <- res0$diff_load / sd1 #Bootstraping here are for getting sd
  #res0$diff_load: the original diff_load, without bootstrapping, a vector of nPC long
  

  return(res0)
}


#' Calculates the average load difference between the cell-cycle genes and
#' control genes for each component.
#'
#' This is an interal function for use by "bootstrap_diff" only.
#'
#' @param xy The data for the genes which are annotated to the cell-cycle, i.e.
#' those genes for which "if_cc" is \code{TRUE}.
#' @param xn The data for the genes which are not annotated to the cell-cycle,
#' control genes, genes for which "if_cc" is \code{FALSE}
#'
#' @return A data frame containing the loadings for each component on the
#' cell-cycle and control genes.

get_diff <- function(xy, xn, rank.)
{
  xn_pca <- prcomp(xn, scale.=FALSE) #Li removed stats::
  #xn_pca %>% names()
  #[1] "sdev"     "rotation" "center"   "scale"    "x"  
  #sdev: diagonal matrix?
  #rotation: the variance of each cell explained at each PC; row-cell, col-PC
  #x: the loading of the genes onto the new reduced space; row-gene, col-PC
  #Be aware that the PCs are variance summary of cells not genes
  #xn_pca$x: row-genes, col-PCs
  #rank. optionally, a number specifying the maximal rank, 
  #i.e., maximal number of principal components to be used. 
  #Can be set as alternative or in addition to tol, 
  #useful notably when the desired rank is considerably smaller than the dimensions of the matrix.
  xn_proj <- xn %*% (xn_pca$rotation[, 1:rank.]) #output: row: genes, col: PCs
  xy_proj <- xy %*% (xn_pca$rotation[, 1:rank.])

  xn_load <- sqrt(colMeans(xn_proj ^ 2))
  xy_load <- sqrt(colMeans(xy_proj ^ 2))

  return(data.frame(xn_load=xn_load, xy_load=xy_load,
                    diff_load=xy_load-xn_load))
  #diff_load is a vector of nPC long
}
