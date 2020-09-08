rm(list = ls())
setwd(
  "C:\\Users\\liorf\\OneDrive\\Documents\\University\\year 3\\Project\\proj_scwgbs\\classifier"
)


library(ggplot2)
library(dplyr)
library(reshape2)
library(hrbrthemes)
library(zoo)
library(scales)
library(RColorBrewer)



####################################################################################################
#                                         functions                                                #
####################################################################################################

##############################################
#               data parsing                 #
##############################################


read.data <- function(dataset, hypo, k, dmers) {
  path <- paste0("attribute_data\\", dataset, "_", hypo, "_", k, "mer_all_seq.csv")
  seqs <- read.csv(path, header = FALSE)
  seqs <- data.frame(lapply(seqs, factor, levels=dmers))
  return(seqs)
}

create.counts <- function(seqs) {
  data.counts <- as.data.frame(sapply(X = seqs, FUN = table))
  # data.counts$nuc <- rownames(data.counts)
  return(data.counts)
}

melt.counts <- function(data.counts, seqs) {
  # Create frequencies and smooth frequencies
  data.freqs <- data.counts[, 1:ncol(data.counts) - 1] / nrow(seqs)
  data.freqs.smooth <- t(data.freqs)
  data.freqs.smooth <- rollmean(data.freqs.smooth, 3)
  data.freqs.smooth <- data.frame(t(data.freqs.smooth))
  
  # add nucleotide column
  data.counts$nuc <- rownames(data.counts)
  data.freqs$nuc <- rownames(data.freqs)
  data.freqs.smooth$nuc <- rownames(data.freqs.smooth)
  
  # melt data
  data.count.melted <- melt(id.vars = "nuc", data.counts, variable.name = "ID", value.name = "count")
  data.freqs.melted <- melt(id.vars = "nuc", data.freqs, variable.name = "ID", value.name = "freq")
  data.freqs.smooth.melted <- melt(id.vars = "nuc", data.freqs.smooth, variable.name = "ID", value.name = "smooth.freq")
  
  data.count.melted$ID <- as.numeric(data.count.melted$ID)
  data.freqs.melted$ID <- as.numeric(data.freqs.melted$ID)
  data.freqs.smooth.melted$ID <- as.numeric(data.freqs.smooth.melted$ID)
  
  #combine data
  data.merged <- merge(data.count.melted, data.freqs.melted, all=TRUE) %>%
    merge(data.freqs.smooth.melted, all = TRUE)
  
  
  
  # data.counts.melted <- melt(id.vars = "nuc", data.counts, variable.name = "ID", value.name = "count")
  # data.counts.melted$freq <- data.counts.melted$count/nrow(seqs)
  
  # data.counts.melted$location <- as.numeric(data.counts.melted$ID)
  data.merged$location <- sapply(data.merged$ID, function(x) ifelse(x <= 75, x - 75, x - 76))
  return(data.merged)
}

create.df <- function(dataset, hypo, k, dmers) {
  data <- read.data(dataset, hypo, k, dmers)
  data.counts <- create.counts(data)
  data.melted <- melt.counts(data.counts, data)
  data.melted$hypo <- hypo
  return(data.melted)
}

combine.df <- function(prone.df, resistant.df) {
  df <- rbind(prone.df, resistant.df)
  return(df)
}

##############################################
#                  plots                     #
##############################################

plot.single.frequencies <- function(freq.df, nucleotide, dataset, k) {
  freq.df$grp <- sapply(freq.df$location, function(x) ifelse(x < -3, "low", ifelse(x > 3, "high", "mid")))
  freq.df$grp <- with(freq.df, paste0(grp, hypo))
  freq.df$further.loc <- freq.df$location
  # freq.df[which(freq.df$location == 0), c("count", "freq")] = NaN
  freq.df[which((freq.df$location >= -4) & (freq.df$location <= 4)), c("further.loc", "freq", "smooth.freq")] = NaN
  # create plot
  p <- freq.df %>%
    filter(nuc == nucleotide) %>%
    ggplot(aes(x = further.loc, y = smooth.freq, color = hypo)) +
    geom_point() + 
    geom_line(aes(group = grp)) +
    # geom_smooth(method = "loess") +
    # geom_smooth(span = 0.04, se = F, aes(group = groups)) + 
    theme_minimal() +
    theme(text = element_text(size=20)) +
    # scale_x_continuous(n.breaks = 20) +
    scale_x_continuous(
      breaks = c(-70, -60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70)) +
    ylab(paste(nucleotide, "frequency")) +
    xlab("distance from CpG site") +
    ylim(0, NA)
  
  # save plot
  # p
  out.path = paste0("frequency_plots\\", dataset, '_nn_', nucleotide, '_', k, 'mer_all_frequencies.png')
  if (k == 2) {
    h <- 7
  } else {
    h <- 8
  }
  ggsave(out.path, width = 12.05, height = h)
}

plot.all.frequencies <- function(freq.df, dataset, k) {
  # freq.df[which(freq.df$location == 0), c("count", "freq")] = NaN
  freq.df$grp <- sapply(freq.df$location, function(x) ifelse(x < -3, "low", ifelse(x > 3, "high", "mid")))
  freq.df$grp <- with(freq.df, paste0(grp, hypo))
  freq.df$further.loc <- freq.df$location
  freq.df[which((freq.df$location >= -4) & (freq.df$location <= 4)), c("further.loc", "freq", "smooth.freq")] = NaN
  
  # create plot
  p <- ggplot(freq.df, aes(x = further.loc, y = smooth.freq, color = hypo)) +
    geom_point() + 
    geom_line(aes(group = grp)) +
    facet_grid(nuc ~ .) + #, scales = "free_y") +
    theme_minimal() +
    theme(text = element_text(size=20)) +
    # scale_x_continuous(n.breaks = 20) +
    scale_x_continuous(
      breaks = c(-70, -60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70)) +
    ylab("frequency") +
    xlab("distance from CpG site")+
    ylim(0, NA)
  
  # save plot
  # p
  out.path = paste0("frequency_plots\\", dataset, '_nn_', k, 'mer_all_frequencies.png')
  ggsave(out.path)
}

plot.heatmap <- function(freq.df, dataset, k) {
  # create plot
  p <- ggplot(freq.df, aes(x = as.numeric(ID), y = nuc, fill = freq)) +
    geom_tile() +
    theme_minimal() +
    scale_x_continuous(n.breaks = 20) +
    ylab("nucleotide") +
    xlab("distance from CpG site") +
    scale_fill_viridis_c(direction = -1) +
    facet_grid(hypo ~ .)
  
  # save plot
  # p
  out.path = paste0("frequency_plots\\", dataset, '_nn_', k, 'mer_all_frequencies_heatmap.png')
  ggsave(out.path)
}

plot.single.heatmap <- function(freq.df, dataset, k, type) {
  # create plot
  p <- freq.df %>%
    filter(hypo == "prone") %>%
    ggplot(aes(x = as.numeric(ID), y = nuc, fill = freq)) +
    geom_tile() +
    theme_minimal() +
    scale_x_continuous(n.breaks = 20) +
    ylab("nucleotide") +
    xlab("distance from CpG site") +
    scale_fill_viridis_c(direction = -1)# +

  # save plot
  # p
  out.path = paste0("frequency_plots\\", dataset, '_nn_', k, 'mer_all_', type, '_frequencies_heatmap.png')
  ggsave(out.path)
}

plot.log.heatmap <- function(prone.df, resistant.df, dataset, k, w, h, order = NULL) {
  freq.df <- prone.df
  freq.df$smooth.freq <- log2(freq.df$smooth.freq / resistant.df$smooth.freq)
  freq.df[which((freq.df$smooth.freq == Inf) | (freq.df$smooth.freq == -Inf) | (freq.df$smooth.freq < -50)), c("freq", "smooth.freq")] = NaN
  
  #clustering
  long.freq.df <- freq.df[, c("nuc", "location", "ID", "smooth.freq")]
  wide.freq.df <- dcast(freq.df, nuc~ID, value.var = "smooth.freq")
  if (is.null(order)) {
    ord <- hclust( dist(wide.freq.df, method = "euclidean"))$order
  } else {
    ord <- order
  }
  long.freq.df$nuc <- factor(long.freq.df$nuc, levels = wide.freq.df$nuc[ord])
  
  nocg.freq.df <- long.freq.df %>%
    filter(!stringr::str_detect(nuc, "CG"))
  nocg.freq.df[which((nocg.freq.df$location >= -5) & (nocg.freq.df$location <= 5)), "smooth.freq"] = NaN
  
  
  # nocg.min <- min(nocg.freq.df$smooth.freq, na.rm = TRUE)
  # nocg.max <- max(nocg.freq.df$smooth.freq, na.rm = TRUE)
  # all.min <- min(long.freq.df$smooth.freq, na.rm = TRUE)
  # # all.min <- -14.506370	
  # # all.min <- -12.713923
  # all.max <- max(long.freq.df$smooth.freq, na.rm = TRUE)
  # 
  # # 11 colors
  # colour.values <- c(
  #   all.min,
  #   (all.min + nocg.min) / 2,
  #   nocg.min,
  #   (2/3) * nocg.min,
  #   (1/3) * nocg.min,
  #   0,
  #   (1/3) * nocg.max,
  #   (2/3) * nocg.max,
  #   nocg.max,
  #   (nocg.max + all.max) / 2,
  #   all.max
  # )
  # # colours <- rev(brewer.pal(11, "Spectral"))
  # # colours[6] <- "#FFFFFF"
  # 
  # colours <- c("#4dac26", muted("blue"), "white", muted("red"), "orange")
  # 
  # # # 5 colors
  # # colour.values <- c(
  # #   all.min,
  # #   nocg.min,
  # #   0,
  # #   nocg.max,
  # #   all.max
  # # )
  # # colours <- rev(brewer.pal(5, "Spectral"))
  # # colours[3] <- "#FFFFFF"
  nocg.min.95 <-
    quantile(nocg.freq.df$smooth.freq[nocg.freq.df$smooth.freq < 0],
             na.rm = T,
             probs = c(0.05))[[1]]
  nocg.max95 <-
    quantile(nocg.freq.df$smooth.freq[nocg.freq.df$smooth.freq > 0],
             na.rm = T,
             probs = c(0.95))[[1]]
  all.min.95 <-
    quantile(long.freq.df$smooth.freq[long.freq.df$smooth.freq < 0],
             na.rm = T,
             probs = c(0.03))[[1]]
  all.max95 <-
    quantile(long.freq.df$smooth.freq[long.freq.df$smooth.freq > 0],
             na.rm = T,
             probs = c(0.95))[[1]]
  all.min <- min(long.freq.df$smooth.freq, na.rm = TRUE)
  all.max <- max(long.freq.df$smooth.freq, na.rm = TRUE)
  
  colour.values <- c(all.min, 0, all.max)
  # nocg <- c(-1, 1) * max(abs(c(nocg.min.95, nocg.max95)))
  # all <- c(-1, 1) * min(max(abs(c(all.min.95, all.max95))), all.max)
  nocg <- c(nocg.min.95, nocg.max95)
  all <- c(all.min.95, all.max95)
  
  
  # final colors (7)
  colours <- c("#4dac26", "#4dac26", muted("blue"), "white", muted("red"), "orange", "orange")
  colour.values <- sort(c(colour.values, nocg, all))
  
  
  # create plot
  p <- ggplot(long.freq.df, aes(x = location, y = nuc, fill = smooth.freq)) +
    geom_tile() +
    theme_minimal() +
    theme(text = element_text(size=20)) +
    scale_x_continuous(
        breaks = c(-70, -60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70)) +
    # scale_fill_gradientn(colours = c(muted("blue"), muted("blue"), "white", muted("red"), muted("red")),
    #                      values = rescale(c(all.min, nocg.min, 0, nocg.max, all.max)),
    #                      breaks = c(nocg.min, 0, nocg.max),
    #                      na.value = "#d3d3d3",
    #                      labels = scales::number_format(accuracy = 0.1)) +
    # guides(fill = guide_colourbar(barwidth = 0.8, barheight = 10)) +
      guides(fill = guide_colourbar(barwidth = 0.5, barheight = 21)) +
      # scale_fill_distiller(palette="Spectral", limits = c(-1, 1) * max(abs(long.freq.df$smooth.freq), na.rm = T)) +
      scale_fill_gradientn(colours = colours,
                           values = rescale(colour.values)
                           # breaks = colour.values,
                           # labels = scales::number_format(accuracy = 0.1),
                           # limits = c(all.min, all.max)
                           ) +
    ylab("nucleotide") +
    xlab("distance from CpG site") +
    labs(fill =  expression(paste(
      log[2], bgroup("(", frac(scriptstyle(prone_freq), scriptstyle(resistant_freq)), ")")
    ))) +
    ggtitle(paste0(dataset, " ", k, "mer"))

  # save plot
  # p
  out.path = paste0("frequency_plots\\heat\\", dataset, '_nn_', k, 'mer_all_frequencies_log_heatmap_GBWRO_nsym.png')
  ggsave(out.path, width = w, height = h)
  return(ord)
}

# plot.log.heatmap(zhou2.prone.df, zhou2.resistant.df, dataset, k)

plot.log.heatmap.central <- function(prone.df, resistant.df, dataset, k, h) {
  freq.df <- prone.df
  freq.df$smooth.freq <- log2(freq.df$smooth.freq / resistant.df$smooth.freq)
  freq.df[which((freq.df$smooth.freq == Inf) | (freq.df$smooth.freq == -Inf)), c("freq", "smooth.freq")] = NaN
  
  freq.df[which(((freq.df$location <= -5) | (freq.df$location >= 5)) & (!stringr::str_detect(freq.df$nuc, "CG"))), "smooth.freq"] = NaN
  
  
  #clustering
  long.freq.df <- freq.df[, c("nuc", "location", "ID", "smooth.freq")]
  wide.freq.df <- dcast(freq.df, nuc~ID, value.var = "smooth.freq")
  ord <- hclust( dist(wide.freq.df, method = "euclidean"))$order
  long.freq.df$nuc <- factor(long.freq.df$nuc, levels = wide.freq.df$nuc[ord])

    all.min <- min(long.freq.df$smooth.freq, na.rm = TRUE)
  all.max <- max(long.freq.df$smooth.freq, na.rm = TRUE)
  
  # create plot
  p <- ggplot(long.freq.df, aes(x = location, y = nuc, fill = smooth.freq)) +
    geom_tile() +
    theme_minimal() +
    theme(text = element_text(size=20)) +
    scale_x_continuous(
      breaks = c(-70, -60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70)) +
    # scale_fill_gradient2(na.value = "#d3d3d3") +
    scale_fill_gradientn(colours = c(muted("blue"), "white", muted("red")),
                         values = rescale(c(all.min, 0, all.max)),
                         # breaks = c(all.min, nocg.min, 0, nocg.max, all.max),
                         na.value = "#d3d3d3",
                         # labels = scales::number_format(accuracy = 0.1)
    ) +
    guides(fill = guide_colourbar(barwidth = 0.5, barheight = 35)) +
    ylab("nucleotide") +
    xlab("distance from CpG site") +
    labs(fill =  expression(paste(
      log[2], bgroup("(", frac(scriptstyle(prone_freq), scriptstyle(resistant_freq)), ")")
    ))) 
  
  # save plot
  # p
  out.path = paste0("frequency_plots\\", dataset, '_nn_', k, 'mer_all_frequencies_log_heatmap_central.png')
  ggsave(out.path, height = h)
}

plot.histogram <- function(prone.df, resistant.df, dataset, k, w, h) {
  freq.df <- prone.df
  freq.df$smooth.freq <- log2(freq.df$smooth.freq / resistant.df$smooth.freq)
  freq.df[which((freq.df$smooth.freq == Inf) | (freq.df$smooth.freq == -Inf)), c("freq", "smooth.freq")] = NaN
  
  #clustering
  long.freq.df <- freq.df[, c("nuc", "location", "ID", "smooth.freq")]

  nocg.freq.df <- long.freq.df %>%
    filter(!stringr::str_detect(nuc, "CG"))
  nocg.freq.df[which((nocg.freq.df$location >= -5) & (nocg.freq.df$location <= 5)), "smooth.freq"] = NaN
  
  
  nocg.min.95 <-
    quantile(nocg.freq.df$smooth.freq[nocg.freq.df$smooth.freq < 0],
             na.rm = T,
             probs = c(0.05))[[1]]
  nocg.max95 <-
    quantile(nocg.freq.df$smooth.freq[nocg.freq.df$smooth.freq > 0],
             na.rm = T,
             probs = c(0.95))[[1]]
  all.min.95 <-
    quantile(long.freq.df$smooth.freq[long.freq.df$smooth.freq < 0],
             na.rm = T,
             probs = c(0.03))[[1]]
  all.max95 <-
    quantile(long.freq.df$smooth.freq[long.freq.df$smooth.freq > 0],
             na.rm = T,
             probs = c(0.98))[[1]]
  all.min <- min(long.freq.df$smooth.freq, na.rm = TRUE)
  all.max <- max(long.freq.df$smooth.freq, na.rm = TRUE)
  
  # nocg.min <- min(nocg.freq.df$smooth.freq, na.rm = TRUE)
  # nocg.max <- max(nocg.freq.df$smooth.freq, na.rm = TRUE)
  # all.min <- min(long.freq.df$smooth.freq, na.rm = TRUE)
  # all.max <- max(long.freq.df$smooth.freq, na.rm = TRUE)
  
  if (k == 3) {
    if (dataset == "bian") {
      all.min <- -12.713923
    }
    else {
      all.min <- -14.506370
    }
  }
  
  # colour.values <- c(
  #   all.min,
  #   nocg.min,
  #   0,
  #   nocg.max,
  #   all.max
  # )
  
  colour.values <- c(all.min, 0, all.max)
  nocg <- c(nocg.min.95, nocg.max95)
  all <- c(all.min.95, all.max95)
  # nocg <- c(-1, 1) * max(abs(c(nocg.min.95, nocg.max95)))
  # all <- c(-1, 1) * min(max(abs(c(all.min.95, all.max95))), all.max)
  
  long.freq.df$class <- "nocg"
  long.freq.df[which(stringr::str_detect(long.freq.df$nuc, "CG")), "class"] = "CG"
  long.freq.df[which((nocg.freq.df$location >= -5) & (nocg.freq.df$location <= 5)), "class"] = "central"
  
  
  ggplot(long.freq.df, aes(x = smooth.freq, after_stat(density), fill = class)) +
    geom_histogram(bins = 500, na.rm = T) +
    geom_vline(xintercept = colour.values, linetype = "dashed", size = 0.4) +
    geom_vline(xintercept = nocg, linetype = "dashed", size = 0.4, color = "red") +
    geom_vline(xintercept = all, linetype = "dashed", size = 0.4, color = "blue") +
    ggtitle(paste0(dataset, " ", k, "mer"))
  
  out.path = paste0("frequency_plots\\heat\\", dataset, '_nn_', k, 'mer_all_frequencies_histogram_nsym.png')
  ggsave(out.path)
  # return(long.freq.df)
}

####################################################################################################
#                                            main                                                  #
####################################################################################################

main <- function(dataset, k, dmers) {
  prone.df <- create.df(dataset, 'prone', k, dmers[[k]])
  resistant.df <- create.df(dataset, 'resistant', k, dmers[[k]])
  
  combined.df <- combine.df(prone.df, resistant.df)
  
  if (k == 1) {
    plot.single.frequencies(combined.df, 'C', dataset, k)
    plot.all.frequencies(combined.df, dataset, k)
    # plot.log.heatmap(prone.df, resistant.df, dataset, k, 4)
  }
  
  if (k == 2) {
    # plot.single.frequencies(combined.df, 'CG', dataset, k)
    plot.log.heatmap(prone.df, resistant.df, dataset, k, 13, 6.8)
    # plot.log.heatmap.central(prone.df, resistant.df, dataset, k, NA)
  }
  
  if (k == 3) {
    plot.log.heatmap(prone.df, resistant.df, dataset, k, 8, 15)
    # plot.log.heatmap.central(prone.df, resistant.df, dataset, k, 15)
  }
}

dmers <- list(c("A", "C", "G", "T"), 
              c("AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT"), 
              c("AAA","AAC","AAG","AAT","ACA","ACC","ACG","ACT","AGA","AGC","AGG","AGT","ATA","ATC","ATG","ATT","CAA","CAC","CAG","CAT","CCA","CCC","CCG","CCT","CGA","CGC","CGG","CGT","CTA","CTC","CTG","CTT","GAA","GAC","GAG","GAT","GCA","GCC","GCG","GCT","GGA","GGC","GGG","GGT","GTA","GTC","GTG","GTT","TAA","TAC","TAG","TAT","TCA","TCC","TCG","TCT","TGA","TGC","TGG","TGT","TTA","TTC","TTG","TTT"))

#################################################################################
dataset <- "bian"

k <- 2
bian.2.prone.df <- create.df(dataset, 'prone', k,  dmers[[k]])
bian.2.resistant.df <- create.df(dataset, 'resistant', k, dmers[[k]])

k <- 3
bian.3.prone.df <- create.df(dataset, 'prone', k,  dmers[[k]])
bian.3.resistant.df <- create.df(dataset, 'resistant', k, dmers[[k]])

dataset <- "zhou"

k <- 2
zhou.2.prone.df <- create.df(dataset, 'prone', k,  dmers[[k]])
zhou.2.resistant.df <- create.df(dataset, 'resistant', k, dmers[[k]])

k <- 3
zhou.3.prone.df <- create.df(dataset, 'prone', k,  dmers[[k]])
zhou.3.resistant.df <- create.df(dataset, 'resistant', k, dmers[[k]])

dataset <- "zhou"

k <- 2
plot.histogram(zhou.2.prone.df, zhou.2.resistant.df, dataset, k, 13, 6.8)
zhou.2.ord <- plot.log.heatmap(prone.df = zhou.2.prone.df,
                               resistant.df = zhou.2.resistant.df,
                               dataset = dataset, 
                               k = k, 
                               w = 13, 
                               h = 6.8)
k <- 3
plot.histogram(zhou.3.prone.df, zhou.3.resistant.df, dataset, k, 13, 6.8)
zhou.3.ord <- plot.log.heatmap(prone.df = zhou.3.prone.df,
                               resistant.df = zhou.3.resistant.df,
                               dataset = dataset, 
                               k = k, 
                               w = 8, 
                               h = 15)

dataset <- "bian"

k <- 2
plot.histogram(bian.2.prone.df, bian.2.resistant.df, dataset, k, 13, 6.8)
plot.log.heatmap(prone.df = bian.2.prone.df,
                 resistant.df = bian.2.resistant.df,
                 dataset = dataset, 
                 k = k, 
                 w = 13, 
                 h = 6.8,
                 order = zhou.2.ord)
k <- 3
plot.histogram(bian.3.prone.df, bian.3.resistant.df, dataset, k, 13, 6.8)
plot.log.heatmap(prone.df = bian.3.prone.df,
                 resistant.df = bian.3.resistant.df,
                 dataset = dataset, 
                 k = k, 
                 w = 8, 
                 h = 15,
                 order = zhou.3.ord)


print('starting bian')
# main('bian', 1, ten.perc, dmers)
# main('bian', 2, ten.perc, dmers)
# main('bian', 3, ten.perc, dmers)
# 
# print('starting zhou')
# # main('zhou', 1, ten.perc, dmers)
# main('zhou', 2, ten.perc, dmers)
# main('zhou', 3, ten.perc, dmers)

# print('done!')


