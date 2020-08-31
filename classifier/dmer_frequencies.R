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


####################################################################################################
#                                         functions                                                #
####################################################################################################

##############################################
#               data parsing                 #
##############################################


read.data <- function(dataset, hypo, k, perc, dmers) {
  path <- paste0("attribute_data\\", dataset, "_", hypo, "_", k, "mer_all_seq.csv")
  seqs <- read.csv(path, header = FALSE)
  # seqs <- seqs[1:perc,]
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

create.df<- function(dataset, hypo, k, perc, dmers) {
  data <- read.data(dataset, hypo, k, perc, dmers)
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
  ggsave(out.path)
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

plot.log.heatmap <- function(prone.df, resistant.df, dataset, k, h) {
  freq.df <- prone.df
  freq.df$smooth.freq <- log2(freq.df$smooth.freq / resistant.df$smooth.freq)
  freq.df[which((freq.df$smooth.freq == Inf) | (freq.df$smooth.freq == -Inf)), c("freq", "smooth.freq")] = NaN
  
  #clustering
  long.freq.df <- freq.df[, c("nuc", "location", "ID", "smooth.freq")]
  wide.freq.df <- dcast(freq.df, nuc~ID, value.var = "smooth.freq")
  ord <- hclust( dist(wide.freq.df, method = "euclidean"))$order
  long.freq.df$nuc <- factor(long.freq.df$nuc, levels = wide.freq.df$nuc[ord])
  
  nocg.freq.df <- long.freq.df %>%
    filter(!stringr::str_detect(nuc, "CG"))
  nocg.freq.df[which((nocg.freq.df$location >= -5) & (nocg.freq.df$location <= 5)), "smooth.freq"] = NaN
  
  
  nocg.min <- min(nocg.freq.df$smooth.freq, na.rm = TRUE)
  nocg.max <- max(nocg.freq.df$smooth.freq, na.rm = TRUE)
  all.min <- min(long.freq.df$smooth.freq, na.rm = TRUE)
  all.max <- max(long.freq.df$smooth.freq, na.rm = TRUE)
  
  # create plot
  p <- ggplot(long.freq.df, aes(x = location, y = nuc, fill = smooth.freq)) +
    geom_tile() +
    theme_minimal() +
    theme(text = element_text(size=20)) +
    scale_x_continuous(
        breaks = c(-70, -60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70)) +
    scale_fill_gradientn(colours = c(muted("blue"), muted("blue"), "white", muted("red"), muted("red")),
                         values = rescale(c(all.min, nocg.min, 0, nocg.max, all.max)),
                         breaks = c(nocg.min, 0, nocg.max),
                         na.value = "#d3d3d3", 
                         labels = scales::number_format(accuracy = 0.1)) +
    guides(fill = guide_colourbar(barwidth = 0.8, barheight = 13)) +
    ylab("nucleotide") +
    xlab("distance from CpG site") +
    labs(fill =  expression(paste(
      log[2], bgroup("(", frac(scriptstyle(prone_freq), scriptstyle(resistant_freq)), ")")
    ))) 
  
  # save plot
  p
  # out.path = paste0("frequency_plots\\", dataset, '_nn_', k, 'mer_all_frequencies_log_heatmap_capped_colours.png')
  # ggsave(out.path, height = h)
}

plot.log.heatmap(zhou2.prone.df, zhou2.resistant.df, dataset, k)

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



####################################################################################################
#                                            main                                                  #
####################################################################################################

main <- function(dataset, k, ten.perc, dmers) {
  prone.df <- create.df(dataset, 'prone', k, ten.perc[["bian"]], dmers[[k]])
  resistant.df <- create.df(dataset, 'resistant', k, ten.perc[["bian"]], dmers[[k]])
  
  combined.df <- combine.df(prone.df, resistant.df)
  
  if (k == 1) {
    plot.single.frequencies(combined.df, 'C', dataset, k)
    plot.all.frequencies(combined.df, dataset, k)
    plot.log.heatmap(prone.df, resistant.df, dataset, k, 4)
  }
  
  if (k == 2) {
    # plot.single.frequencies(combined.df, 'CG', dataset, k)
    # plot.log.heatmap(prone.df, resistant.df, dataset, k, NA)
    plot.log.heatmap.central(prone.df, resistant.df, dataset, k, NA)
  }
  
  if (k == 3) {
    plot.log.heatmap(prone.df, resistant.df, dataset, k, 15)
    plot.log.heatmap.central(prone.df, resistant.df, dataset, k, 15)
  }
}

ten.perc <- c(bian = 64779, zhou = 105933)
dmers <- list(c("A", "C", "G", "T"), 
              c("AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT"), 
              c("AAA","AAC","AAG","AAT","ACA","ACC","ACG","ACT","AGA","AGC","AGG","AGT","ATA","ATC","ATG","ATT","CAA","CAC","CAG","CAT","CCA","CCC","CCG","CCT","CGA","CGC","CGG","CGT","CTA","CTC","CTG","CTT","GAA","GAC","GAG","GAT","GCA","GCC","GCG","GCT","GGA","GGC","GGG","GGT","GTA","GTC","GTG","GTT","TAA","TAC","TAG","TAT","TCA","TCC","TCG","TCT","TGA","TGC","TGG","TGT","TTA","TTC","TTG","TTT"))


# print('starting bian')
# main('bian', 1, ten.perc, dmers)
# main('bian', 2, ten.perc, dmers)
# main('bian', 3, ten.perc, dmers)
# 
# print('starting zhou')
# main('zhou', 1, ten.perc, dmers)
# main('zhou', 2, ten.perc, dmers)
# main('zhou', 3, ten.perc, dmers)
# 
# print('done!')



# k <- 2
k <- 3
dataset <- 'bian'

zhou2.prone.df <- create.df(dataset, 'prone', k, ten.perc[["bian"]], dmers[[k]])
zhou2.resistant.df <- create.df(dataset, 'resistant', k, ten.perc[["bian"]], dmers[[k]])
# 
# #---------------------------
# freq.df <- zhou2.prone.df
# freq.df$smooth.freq <- log2(freq.df$smooth.freq / zhou2.resistant.df$smooth.freq)
# freq.df[which((freq.df$smooth.freq == Inf) | (freq.df$smooth.freq == -Inf)), c("freq", "smooth.freq")] = NaN
# 
# #clustering
# long.freq.df <- freq.df[, c("nuc", "location", "ID", "smooth.freq")]
# wide.freq.df <- dcast(freq.df, nuc~ID, value.var = "smooth.freq")
# ord <- hclust( dist(wide.freq.df, method = "euclidean"))$order
# #---------------------------
# 
# 
# 
# 
# k <- 2
# h <- NA
# # k <- 3
# # h <- 15
# dataset <- 'bian'
# 
# bian2.prone.df <- create.df(dataset, 'prone', k, ten.perc[["bian"]], dmers[[k]])
# bian2.resistant.df <- create.df(dataset, 'resistant', k, ten.perc[["bian"]], dmers[[k]])
# 
# #------------------
# freq.df <- bian2.prone.df
# freq.df$smooth.freq <- log2(freq.df$smooth.freq / bian2.resistant.df$smooth.freq)
# freq.df[which((freq.df$smooth.freq == Inf) | (freq.df$smooth.freq == -Inf)), c("freq", "smooth.freq")] = NaN
# 
# #clustering
# long.freq.df <- freq.df[, c("nuc", "location", "ID", "smooth.freq")]
# wide.freq.df <- dcast(freq.df, nuc~ID, value.var = "smooth.freq")
# # ord <- hclust( dist(wide.freq.df, method = "euclidean"))$order
# long.freq.df$nuc <- factor(long.freq.df$nuc, levels = wide.freq.df$nuc[ord])
# 
# nocg.freq.df <- long.freq.df %>%
#   filter(!stringr::str_detect(nuc, "CG"))
# nocg.freq.df[which((nocg.freq.df$location >= -5) & (nocg.freq.df$location <= 5)), "smooth.freq"] = NaN
# 
# 
# nocg.min <- min(nocg.freq.df$smooth.freq, na.rm = TRUE)
# nocg.max <- max(nocg.freq.df$smooth.freq, na.rm = TRUE)
# all.min <- min(long.freq.df$smooth.freq, na.rm = TRUE)
# all.max <- max(long.freq.df$smooth.freq, na.rm = TRUE)
# 
# # create plot
# p <- ggplot(long.freq.df, aes(x = location, y = nuc, fill = smooth.freq)) +
#   geom_tile() +
#   theme_minimal() +
#   theme(text = element_text(size=20)) +
#   scale_x_continuous(
#     breaks = c(-70, -60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70)) +
#   scale_fill_gradientn(colours = c(muted("blue"), muted("blue"), "white", muted("red"), muted("red")),
#                        values = rescale(c(all.min, nocg.min, 0, nocg.max, all.max)),
#                        breaks = c(all.min, nocg.min, 0, nocg.max, all.max),
#                        na.value = "#d3d3d3", 
#                        labels = scales::number_format(accuracy = 0.1)) +
#   guides(fill = guide_colourbar(barwidth = 0.5, barheight = 13)) +
#   ylab("nucleotide") +
#   xlab("distance from CpG site") +
#   labs(fill =  expression(paste(
#     log[2], bgroup("(", frac(scriptstyle(prone_freq), scriptstyle(resistant_freq)), ")")
#   ))) 
# 
# # save plot
# # p
# out.path = paste0("frequency_plots\\", dataset, '_nn_', k, 'mer_all_frequencies_log_heatmap_capped_colours_zhou_order.png')
# ggsave(out.path, height = h)
# #------------------------------