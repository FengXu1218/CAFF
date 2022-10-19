#' build_pseudo_order
#' This function is used to build a quasi-timing sequence using
#' a number of different methods.
#'
#' @param data Expression profile data of cell cycle genes
#' @param reverse Whether to change the direction of pseudo order
#' @param method Which method to use
#'
#' @return
#' @export
#'
#' @examples
build_pseudo_order <- function(data, reverse = FALSE,
                               method = c("Default", "tricycle", "reCAT")){
  if (method == "Default") {
    ## default method ===============
    pca <- prcomp(t(data))
    score <- pca$x
    sdev <- pca$sdev
    if (reverse == T) {
      score[,2] <- -score[,2]
    }

    # convert rectangular coordinates to polar coordinates:
    THETA <- apply(score[,1:2], 1, cart2pol)[1,]
    names(THETA) <- colnames(data)
    RHO <- apply(score[,1:2], 1, cart2pol)[2,]

    # determine the pseudo order according to the theta value:
    pseudo_order <- order(THETA)
    pseudo_order_rank <- rank(THETA)
    # names(pseudo_order) <- colnames(data)

    return(list(pseudo_order = pseudo_order,
                pseudo_order_rank = pseudo_order_rank,
                score = score,
                sdev = sdev,
                THETA = THETA))
  } else if (method == "tricycle"){
    ## tricycle method ==============
    library(tricycle)
    tricycle_order <- estimate_cycle_position(data,
                                              gname.type = "SYMBOL",
                                              species = "human")

    tricycle_order_rank <- rank(tricycle_order)
    return(tricycle_order_rank)
  } else if (method == "reCAT"){
    ## reCAT method ==============
    setwd("reCAT-master/R")
    source("get_ordIndex.R")
    reCAT_order <- get_ordIndex(t(data), 10)
    setwd("../../")
    return(reCAT_order)
  } else {
    message("The method must be one of 'Default', 'tricycle',or 'reCAT'.")
  }
}

#' cor_scatter_plot
#'
#'
#' This function is used to plot a scatter plot that maximizes correlation
#'
#' @param index_data a data frame.
#'
#' @return
#' @export
#'
#' @examples
cor_scatter_plot <- function(index_data){
  p_list <- list()
  library(Hmisc)

  t <- 1
  for (i in 1:(ncol(index_data)-1)) {
    for (j in 2:ncol(index_data)) {
      if (j > i) {
        res_out <- rings2aligment(index_data[,j], index_data[,i])
        res <- rcorr(index_data[,j][res_out$output_posi], index_data[,i])
        p_value <- signif(res$P[1,2], 2)
        cor_value <- round(res$r[1,2], 2)

        index_data_new <- index_data[,c(j,i)]
        index_data_new[,1] <- index_data_new[res_out$output_posi,1]
        colnames(index_data_new) <- c("x", "y")

        p <- ggplot(index_data_new, aes(x = x, y = y))+
          geom_point(color = "#988d7b")+
          geom_smooth(method = "lm", formula = y ~ x,
                      fill = "#b2e7fa", color = "#00aeef", alpha = 0.8)+
          theme_bw()+
          xlab(colnames(index_data)[j])+
          ylab(colnames(index_data)[i])+
          theme(
            panel.grid = element_blank(),
            axis.title = element_text(face = "bold.italic"),
            plot.title = element_text(hjust = 0.5, size = 10)
          )+
          labs(title = paste0("r =", cor_value, ", q = ", p_value))
        p_list[[t]] <- p
        t <- t+1
      }
    }
  }

  return(p_list)
}


#' rings2aligment
#'
#'
#' This function is used to align the ring structure
#' Adjust the order position of the elements in input_vector,
#' Maximizes the Spearman correlation coefficients of input_vector and reference_vector
#'
#' @param input_vector a vector
#' @param reference_vector another vector
#'
#' @return output_posi, output_coef
#' @export
#'
#' @examples
#' v1=c(0.6,0.7,0.8,0.1,0.2,0.3,0.4,0.5)
#' v2=c(1:8)
#' res <- rings2aligment(v1, v2)
#' res
#' $output_posi
#' [1] 4 5 6 7 8 1 2 3
#'
#' $output_coef
#' [1] 1
#' v1[res$output_posi]

rings2aligment <- function(input_vector, reference_vector){
  len_vector <- length(input_vector)

  temp_order <- matrix(NA, nrow = len_vector, ncol = len_vector)
  for (i in 1:len_vector) {
    temp_order[i:len_vector,i] <- 1:(len_vector-i+1)
    if (i>1) {
      temp_order[1:i-1,i]=(len_vector-i+2):len_vector;
    }
  }

  temp_coef <- c()
  for (i in 1:len_vector) {
    temp_coef[i] <- cor(input_vector[temp_order[,i]],
                        reference_vector, method ='spearman')
  }
  a <- max(abs(temp_coef));
  b <- which.max(abs(temp_coef))

  output_coef <- temp_coef[b]
  output_posi <- temp_order[,b]

  return(list(output_posi = output_posi, output_coef = output_coef))
}
