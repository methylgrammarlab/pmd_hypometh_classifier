rm(list = ls())
setwd(
  "C:\\Users\\liorf\\OneDrive\\Documents\\University\\year 3\\Project\\proj_scwgbs\\classifier"
)


library(ggplot2)
library(dplyr)

####################################################################################################
#                                         functions                                                #
####################################################################################################

calculate.means <- function(seqs, scores, nuc) {
  masked.scores <-
    replace(scores, which(seqs != nuc, arr.ind = T), NA)
  
  means <- colMeans(masked.scores, na.rm = T)
  return(means)
}

create.mean.df <- function(seq.path, score.path) {
  seqs <- read.csv(seq.path)
  scores <- read.csv(score.path)
  df <-
    data.frame(
      a = calculate.means(seqs, scores, 1),
      c = calculate.means(seqs, scores, 2),
      g = calculate.means(seqs, scores, 3),
      t = calculate.means(seqs, scores, 4)
    )
  df$ID <- seq.int(nrow(df))
  return(df)
}

plot.means <- function(means) {
  # create plot
  p <- ggplot(means, aes(x = ID, y = c)) +
    geom_line() +
    geom_point() +
    geom_smooth(method = "gam") +
    theme_minimal() +
    theme(# panel.grid.major.x = element_blank(),
      # axis.ticks.x = element_blank(),
      axis.text.x = element_blank()) +
    ylab("mean attribution score") +
    xlab("position")
  
  # save plot
  out.path = 'nn_c_attribution.png'
  ggsave(out.path)
}

####################################################################################################
#                                            main                                                  #
####################################################################################################

seq.path <- "small_seq.csv"
score.path <- "small_score.csv"
means.df <- create.mean.df(seq.path, score.path)
plot.means(means.df)
