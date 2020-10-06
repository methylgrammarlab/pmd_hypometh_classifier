rm(list = ls())
setwd(
  "C:\\Users\\liorf\\OneDrive\\Documents\\University\\year 3\\Project\\proj_scwgbs\\classifier"
)


library(ggplot2)
library(dplyr)
library(reshape2)
library(zoo)
library(scales)
library(RColorBrewer)



####################################################################################################
#                                         functions                                                #
####################################################################################################

##############################################
#               data parsing                 #
##############################################


read.data <- function(dataset, hypo, k, dmers, cg.density) {
  path <- paste0("orig_meth_above_0.5\\kmer\\density_", cg.density, "\\", dataset, "_", hypo, "_", k, "mer_all_seq.csv")
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

create.df <- function(dataset, hypo, k, dmers, cg.density) {
  data <- read.data(dataset, hypo, k, dmers, cg.density)
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



get.heatmap.order <- function(prone.df, resistant.df, dataset) {
  freq.df <- prone.df
  freq.df$smooth.freq <- log2(freq.df$smooth.freq / resistant.df$smooth.freq)
  freq.df[which((freq.df$smooth.freq == Inf) | (freq.df$smooth.freq == -Inf) | (freq.df$smooth.freq < -50)), c("freq", "smooth.freq")] = NaN
  
  #clustering
  long.freq.df <- freq.df[, c("nuc", "location", "ID", "smooth.freq")]
  wide.freq.df <- dcast(freq.df, nuc~ID, value.var = "smooth.freq")
  ord <- hclust( dist(wide.freq.df, method = "euclidean"))$order
  return(ord)    
}

plot.log.heatmap <- function(prone.df, resistant.df, ord, dataset, cg.density, k, w, h) {
  freq.df <- prone.df
  freq.df$smooth.freq <- log2(freq.df$smooth.freq / resistant.df$smooth.freq)
  freq.df[which((freq.df$smooth.freq == Inf) | (freq.df$smooth.freq == -Inf)), c("freq", "smooth.freq")] = NaN
  
  # ordering
  long.freq.df <- freq.df[, c("nuc", "location", "ID", "smooth.freq")]
  wide.freq.df <- dcast(freq.df, nuc~ID, value.var = "smooth.freq")
  long.freq.df$nuc <- factor(long.freq.df$nuc, levels = wide.freq.df$nuc[ord])
  
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
             probs = c(0.95))[[1]]
  all.min <- min(long.freq.df$smooth.freq, na.rm = TRUE)
  all.max <- max(long.freq.df$smooth.freq, na.rm = TRUE)
  
  colour.values <- c(all.min, 0, all.max)
  nocg <- c(nocg.min.95, nocg.max95)
  all <- c(all.min.95, all.max95)
  
  
  # final colors (7)
  colours <- c("orange", "orange", muted("red"), "white", muted("blue"), "#4dac26", "#4dac26")
  colour.values <- sort(c(colour.values, nocg, all))
  
  
  # create plot
  p <- ggplot(long.freq.df, aes(x = location, y = nuc, fill = smooth.freq)) +
    geom_tile() +
    theme_minimal() +
    theme(text = element_text(size=30)) +
    scale_x_continuous(
      breaks = c(-70, -50, -30, -10, 10, 30, 50, 70)) +
    guides(fill = guide_colourbar(barwidth = 0.5, barheight = 21)) +
      scale_fill_gradientn(colours = colours,
                           values = rescale(colour.values)
                           ) +
    ylab("nucleotide") +
    xlab("distance from CpG site") +
    labs(fill =  expression(paste(
      log[2], bgroup("(", frac(scriptstyle(prone_freq), scriptstyle(resistant_freq)), ")")
    ))) #+
    # ggtitle(paste0(dataset, " ", k, "mer"))

  # save plot
  # p
  out.path = paste0("orig_meth_above_0.5\\graphs\\", dataset, '_nn_cg_density_', cg.density, '_', k, 'mer_all_frequencies_log_heatmap_final.png')
  # out.path = paste0("C:\\Users\\liorf\\OneDrive\\Documents\\University\\year 3\\Project\\proj_scwgbs\\for_slides\\pdfs\\frequency_heatmap\\", dataset, '_nn_cg_density_', cg.density, '_', k, 'mer_all_frequencies_log_heatmap_final.pdf')
  ggsave(out.path, width = w, height = h)
}




####################################################################################################
#                                            main                                                  #
####################################################################################################

main <- function(dataset, k, dmers) {
  prone.df <- create.df(dataset, 'prone', k, dmers[[k]], cg.density)
  resistant.df <- create.df(dataset, 'resistant', k, dmers[[k]], cg.density)
  
  combined.df <- combine.df(prone.df, resistant.df)
  
  if (k == 1) {
    plot.single.frequencies(combined.df, 'C', dataset, k)
    plot.all.frequencies(combined.df, dataset, k)
    # plot.log.heatmap(prone.df, resistant.df, dataset, k, 4)
  }
  
  if (k == 2) {
    # plot.single.frequencies(combined.df, 'CG', dataset, k)
    plot.log.heatmap(prone.df, resistant.df, dataset, cg.density, k, 13, 6.8)
    # plot.log.heatmap.central(prone.df, resistant.df, dataset, k, NA)
  }
  
  if (k == 3) {
    plot.log.heatmap(prone.df, resistant.df, dataset, cg.density, k, 8, 15)
    # plot.log.heatmap.central(prone.df, resistant.df, dataset, k, 15)
  }
}

dmers <- list(c("A", "C", "G", "T"), 
              c("AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT"), 
              c("AAA","AAC","AAG","AAT","ACA","ACC","ACG","ACT","AGA","AGC","AGG","AGT","ATA","ATC","ATG","ATT","CAA","CAC","CAG","CAT","CCA","CCC","CCG","CCT","CGA","CGC","CGG","CGT","CTA","CTC","CTG","CTT","GAA","GAC","GAG","GAT","GCA","GCC","GCG","GCT","GGA","GGC","GGG","GGT","GTA","GTC","GTG","GTT","TAA","TAC","TAG","TAT","TCA","TCC","TCG","TCT","TGA","TGC","TGG","TGT","TTA","TTC","TTG","TTT"))

cg.density <- "below"

#################################################################################
dataset <- "bian"

k <- 2
bian.2.prone.df <- create.df(dataset, 'prone', k,  dmers[[k]], cg.density)
bian.2.resistant.df <- create.df(dataset, 'resistant', k, dmers[[k]], cg.density)

k <- 3
bian.3.prone.df <- create.df(dataset, 'prone', k,  dmers[[k]], cg.density)
bian.3.resistant.df <- create.df(dataset, 'resistant', k, dmers[[k]], cg.density)

dataset <- "zhou"

k <- 2
zhou.2.prone.df <- create.df(dataset, 'prone', k,  dmers[[k]], cg.density)
zhou.2.resistant.df <- create.df(dataset, 'resistant', k, dmers[[k]], cg.density)

k <- 3
zhou.3.prone.df <- create.df(dataset, 'prone', k,  dmers[[k]], cg.density)
zhou.3.resistant.df <- create.df(dataset, 'resistant', k, dmers[[k]], cg.density)

# -------------------------------------------------------------------
zhou.2.ord <- get.heatmap.order(prone.df = zhou.2.prone.df, resistant.df = zhou.2.resistant.df)
zhou.3.ord <- get.heatmap.order(prone.df = zhou.3.prone.df, resistant.df = zhou.3.resistant.df)

dataset <- "zhou"

k <- 2
# plot.histogram(zhou.2.prone.df, zhou.2.resistant.df, dataset, k, 13, 6.8)
plot.log.heatmap(zhou.2.prone.df, zhou.2.resistant.df, zhou.2.ord, dataset, cg.density, k, 13, 6.8)

k <- 3
# plot.histogram(zhou.3.prone.df, zhou.3.resistant.df, dataset, k, 13, 6.8)
plot.log.heatmap(zhou.3.prone.df, zhou.3.resistant.df, zhou.3.ord, dataset, cg.density, k, 10, 18)

dataset <- "bian"
k <- 2
# plot.histogram(bian.2.prone.df, bian.2.resistant.df, dataset, k, 13, 6.8)
plot.log.heatmap(bian.2.prone.df, bian.2.resistant.df, zhou.2.ord, dataset, cg.density, k, 13, 6.8)
k <- 3
# plot.histogram(bian.3.prone.df, bian.3.resistant.df, dataset, k, 13, 6.8)
plot.log.heatmap(bian.3.prone.df, bian.3.resistant.df, zhou.3.ord, dataset, cg.density, k, 10, 18)


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

