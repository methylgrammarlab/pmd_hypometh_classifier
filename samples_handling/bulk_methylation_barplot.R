rm(list=ls())
setwd("C:\\Users\\liorf\\OneDrive\\Documents\\University\\year 3\\Project\\proj_scwgbs\\samples_handling")


library(ggplot2)
library(dplyr)


####################################################################################################
#                                         functions                                                #
####################################################################################################

read.data <- function(methylation.path) {
  data <- read.csv(methylation.path, as.is=F)
  return(data)
}


create.plot <- function(data) {
  # create plot
  p <- ggplot(data, aes(x = reorder(name, methylation), y = methylation, fill = tumor)) +
    geom_bar(stat = 'identity', width = 1) +
    theme_minimal() +
    theme(
      panel.grid.major.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      text = element_text(size=30)
    ) +
    scale_fill_manual(values = c(YES = "#F380A3", NO = "#6FBE44"), labels = c(YES = "tumor", NO = "normal cell")) +
    ylim(0, 1) +
    ylab("mean methyaltion") +
    xlab("tumor cells (self sorted)")
  
  # save plot
  out.path = paste0('my_files\\final_graphs\\bulk', '_cell_methylation_by_', 'subtype','.png')
  ggsave(out.path)
}

####################################################################################################
#                                            main                                                  #
####################################################################################################

methylation.path <- "bulk_mean_meth_update.csv"

data <- read.data(methylation.path)
create.plot(data)

print('done!')


