rm(list = ls())
setwd(
  "C:\\Users\\liorf\\OneDrive\\Documents\\University\\year 3\\Project\\proj_scwgbs\\classifier"
)


library(ggplot2)
library(dplyr)
library(reshape2)


####################################################################################################
#                                         functions                                                #
####################################################################################################

calculate.means <- function(seqs, scores, nuc) {
  masked.scores <-
    replace(scores, which(seqs != nuc, arr.ind = TRUE), NA)

  means <- colMeans(masked.scores, na.rm = TRUE)
  return(means)
}

create.mean.df <- function(seq.path, score.path) {
  seqs <- read.csv(seq.path, header = FALSE)
  scores <- read.csv(score.path, header = FALSE)
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

combine.mean.df <- function(prone.means.df, resistant.means.df) {
  prone.means.df$type <- "prone"
  resistant.means.df$type <- "resistant"
  
  means.df <- rbind(prone.means.df, resistant.means.df)
  
  means.df$dist <- sapply(means.df$ID, function(x) ifelse(x <= 75, x - 75, x - 76))
  return(means.df)
}

plot.means <- function(means) {
  # create plot
  p <- ggplot(means, aes(x = dist, y = c, color = type)) +
    geom_line() +
    theme_minimal() +
    scale_x_continuous(n.breaks = 20) +
    ylab("mean C attribution score") +
    xlab("distance from CpG site")

  # save plot
  out.path = 'bian_nn_c_attribution.png'
  ggsave(out.path)
}

plot.all.means <- function(means) {
  # melt data
  means.long <-
    melt(
      means,
      id.vars = c("ID", "type", "dist"),
      measure.vars = c("a", "c", "g", "t"),
      variable.name = "nuc",
      value.name = "attr"
    )
  
  # create plot
  p <- ggplot(means.long, aes(x = dist, y = attr, color = type)) +
    geom_line() +
    facet_grid(nuc ~ ., scales = "free_y") +
  theme_minimal() +
    scale_x_continuous(n.breaks = 20) +
    ylab("mean attribution score") +
    xlab("distance from CpG site")
  
  # save plot
  out.path = 'bian_nn_attribution_free.png'
  ggsave(out.path)
}

####################################################################################################
#                                            main                                                  #
####################################################################################################

prone.seq.path <- "attribute_data\\bian_prone_seq.csv"
prone.score.path <- "attribute_data\\bian_prone_gradients.csv"
prone.means.df <- create.mean.df(prone.seq.path, prone.score.path)

resistant.seq.path <- "attribute_data\\bian_resistant_seq.csv"
resistant.score.path <- "attribute_data\\bian_resist_gradients.csv"
resistant.means.df <- create.mean.df(resistant.seq.path, resistant.score.path)

mean.df <- combine.mean.df(prone.means.df, resistant.means.df)

plot.means(mean.df)
plot.all.means(mean.df)

mean.df + 3
