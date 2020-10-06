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

create.mean.df <- function(dataset, hypo, medians) {
  base.path <- "orig_meth_above_0.5\\gradients\\"
  seq.path <- paste0(base.path, dataset, "_", hypo, "_seq.csv")
  score.path <- paste0(base.path, dataset, "_", hypo, "_gradients.csv")
  seqs <- read.csv(seq.path, header = FALSE)
  scores <- read.csv(score.path, header = FALSE)
  df.median <- medians[dataset, hypo]
  scores <- replace(scores, which(scores <= df.median, arr.ind = TRUE), 0)
  scores <- replace(scores, which(scores > df.median, arr.ind = TRUE), 1)
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

plot.means <- function(means, dataset) {
  means[which(((means$dist > -5) & (means$dist < 5))), c("a", "c", "g", "t")] = NaN
  # create plot
  p <- ggplot(means, aes(x = dist, y = c, color = type)) +
    geom_line() +
    theme_minimal() +
    theme(text = element_text(size=30)) +
    scale_x_continuous(
      breaks = c(-70, -50, -30, -10, 10, 30, 50, 70)) +
    ylab("mean C attribution score") +
    xlab("distance from CpG site")

  # save plot
  # p
  out.path = paste0("orig_meth_above_0.5\\graphs\\", dataset, '_nn_c_attribution.png')
  out.path = paste0("C:\\Users\\liorf\\OneDrive\\Documents\\University\\year 3\\Project\\proj_scwgbs\\for_slides\\pdfs\\attribution_plots\\", dataset, '_nn_c_attribution.pdf')
  ggsave(out.path, width = 12.05, height = 8)
  
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
  p
  # out.path = 'zhou_nn_attribution_free.png'
  # ggsave(out.path)
}

####################################################################################################
#                                            main                                                  #
####################################################################################################

medians <-
  data.frame(prone = c(zhou = 0.006320351, bian = 0.0062969915),
             resistant = c(zhou = 0.0007945345, bian = 0.001437973))

dataset <- "bian"

prone.means.df <- create.mean.df(dataset, "prone", medians)

resistant.means.df <- create.mean.df(dataset, "resistant", medians)

mean.df <- combine.mean.df(prone.means.df, resistant.means.df)

plot.means(mean.df, dataset)
# plot.all.means(mean.df)

