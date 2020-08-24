rm(list = ls())
setwd(
  "C:\\Users\\liorf\\OneDrive\\Documents\\University\\year 3\\Project\\proj_scwgbs\\classifier"
)


library(ggplot2)
library(dplyr)
library(reshape2)
library(hrbrthemes)


####################################################################################################
#                                         functions                                                #
####################################################################################################

read.data <- function(dataset, hypo, k, perc, dmers) {
  path <- paste0("attribute_data\\", dataset, "_", hypo, "_", k, "mer_seq.csv")
  seqs <- read.csv(path, header = FALSE) #, as.is = FALSE)
  seqs <- seqs[1:perc,]
  seqs <- data.frame(lapply(seqs, factor, levels=dmers))
  return(seqs)
}

create.counts <- function(seqs) {
  data.counts <- as.data.frame(sapply(X = seqs, FUN = table))
  data.counts$nuc <- rownames(data.counts)
  return(data.counts)
}

melt.counts <- function(data.counts, seqs) {
  data.counts.melted <- melt(id.vars = "nuc", data.counts, variable.name = "ID", value.name = "count")
  data.counts.melted$freq <- data.counts.melted$count/length(rownames(seqs))
  
  data.counts.melted$location <- as.numeric(data.counts.melted$ID)
  data.counts.melted$location <- sapply(data.counts.melted$location, function(x) ifelse(x <= 75, x - 75, x - 76))
  return(data.counts.melted)
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


plot.single.frequencies <- function(freq.df, nucleotide, dataset, k) {
  freq.df[which(freq.df$location == 0), c("count", "freq")] = NaN
  # create plot
  p <- freq.df %>%
    filter(nuc == nucleotide) %>%
    ggplot(aes(x = location, y = freq, color = hypo)) +
    geom_line() +
    theme_minimal() +
    scale_x_continuous(n.breaks = 20) +
    ylab(paste(nucleotide, "frequency")) +
    xlab("distance from CpG site")
  
  # save plot
  
  out.path = paste0("frequency_plots\\", dataset, '_nn_', nucleotide, '_', k, 'mer_frequencies.png')
  ggsave(out.path)
}

plot.all.frequencies <- function(freq.df, dataset, k) {
  freq.df[which(freq.df$location == 0), c("count", "freq")] = NaN
  # create plot
  p <- ggplot(freq.df, aes(x = location, y = freq, color = hypo)) +
    geom_line() +
    facet_grid(nuc ~ .) + #, scales = "free_y") +
    theme_minimal() +
    scale_x_continuous(n.breaks = 20) +
    ylab("frequency") +
    xlab("distance from CpG site")
  
  # save plot
  # p
  out.path = paste0("frequency_plots\\", dataset, '_nn_', k, 'mer_frequencies.png')
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
  out.path = paste0("frequency_plots\\", dataset, '_nn_', k, 'mer_frequencies_heatmap.png')
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
  out.path = paste0("frequency_plots\\", dataset, '_nn_', k, 'mer_', type, '_frequencies_heatmap.png')
  ggsave(out.path)
}

main <- function(dataset, k, ten.perc, dmers) {
  prone.df <- create.df(dataset, 'prone', k, ten.perc[["bian"]], dmers[[k]])
  resistant.df <- create.df(dataset, 'resistant', k, ten.perc[["bian"]], dmers[[k]])
  
  combined.df <- combine.df(prone.df, resistant.df)
  
  if (k == 1) {
    plot.single.frequencies(combined.df, 'C', dataset, k)
    plot.all.frequencies(combined.df, dataset, k)
  }
  
  if (k == 2) {
    plot.single.frequencies(combined.df, 'CG', dataset, k)
    plot.heatmap(combined.df, dataset, k)
  }
  
  if (k == 3) {
    plot.single.heatmap(combined.df, dataset, k, "prone")
    plot.single.heatmap(combined.df, dataset, k, "resistant")
  }
}

####################################################################################################
#                                            main                                                  #
####################################################################################################


ten.perc <- c(bian = 64779, zhou = 105933)
dmers <- list(c("A", "C", "G", "T"), 
              c("AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT"), 
              c("AAA","AAC","AAG","AAT","ACA","ACC","ACG","ACT","AGA","AGC","AGG","AGT","ATA","ATC","ATG","ATT","CAA","CAC","CAG","CAT","CCA","CCC","CCG","CCT","CGA","CGC","CGG","CGT","CTA","CTC","CTG","CTT","GAA","GAC","GAG","GAT","GCA","GCC","GCG","GCT","GGA","GGC","GGG","GGT","GTA","GTC","GTG","GTT","TAA","TAC","TAG","TAT","TCA","TCC","TCG","TCT","TGA","TGC","TGG","TGT","TTA","TTC","TTG","TTT"))

print('starting bian')
main('bian', 1, ten.perc, dmers)
main('bian', 2, ten.perc, dmers)
main('bian', 3, ten.perc, dmers)

print('starting zhou')
main('zhou', 1, ten.perc, dmers)
main('zhou', 2, ten.perc, dmers)
main('zhou', 3, ten.perc, dmers)

print('done!')