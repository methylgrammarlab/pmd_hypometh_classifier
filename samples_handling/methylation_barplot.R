rm(list=ls())
setwd("C:\\Users\\liorf\\OneDrive\\Documents\\University\\year 3\\Project\\proj_scwgbs\\samples_handling")


library(ggplot2)
library(dplyr)


####################################################################################################
#                                         functions                                                #
####################################################################################################

read.data <- function(methylation.path) {
  data <- read.csv(methylation.path, as.is=F)
  
  # update factor order
  data$lesion = factor(data$lesion, levels = c("NC", "PT", "LN", "MO", "ML", "MP"))
  data$region = factor(data$region, levels = c("NC", "PT1", "PT2", "PT3", "PT4",
                                               "PT5", "PT6", "LN1", "LN2", "LN3",
                                               "LN4", "LN5", "ML1", "ML2", "ML3", 
                                               "ML4", "ML5", "ML6", "MP1", "MP2", 
                                               "MP3", "MP4", "MP5", "MO1", "MO2", 
                                               "MO3", "MO4", "MO5", "MO6"))
  data$sublineage = factor(data$sublineage, levels = c("NC", "A0", "A1", "A2", "A3",
                                                       "A4", "A5",
                                                       "A6", "A7", "A8", "A9", "B",
                                                       "B0", "B1", "B2", "B3", "C0", 
                                                       "C1", "C2", "C3", "C4", "C5", 
                                                       "undefined"))
  return(data)
}


create.plot <- function(data, patient.name, color.by, my_colour, custom_labels) {
  patient <- data %>%
    filter(patient == patient.name) %>%
    droplevels()
  
  num.tumor = nrow(patient[patient$lesion != "NC", ])
  perc = 20
  
  # create plot
  p <- ggplot(patient, aes(x = reorder(X, mean), y = mean, fill = !!as.symbol(color.by))) +
    geom_bar(stat = 'identity', width = 1) +
    annotate(
      "rect",
      xmin = 0,
      xmax = (num.tumor / 100) * perc,
      ymin = 0,
      ymax = 1,
      color = 'black',
      fill = 'grey',
      alpha = 0.1,
      size = 1
    ) +
    annotate(
      "text", 
      x = ((num.tumor / 100) * perc) / 2, 
      y = 0.97, 
      label=paste0("low ", perc, "%")
    ) +
    annotate(
      "rect",
      xmin = num.tumor - (num.tumor / 100) * perc,
      xmax = num.tumor,
      ymin = 0,
      ymax = 1,
      color = 'black',
      fill = 'grey',
      alpha = 0.1,
      size = 1
    ) +
    annotate(
      "text", 
      x = (num.tumor - (num.tumor / 100) * perc) + (num.tumor - (num.tumor - (num.tumor / 100) * perc)) / 2, 
      y = 0.97, 
      label=paste0("high ", perc, "%")
    ) +
    theme_minimal() +
    theme(
      panel.grid.major.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank()
    ) +
    scale_fill_manual(values = my_colour[[color.by]], labels = custom_labels) +
    ylim(0, 1) +
    ylab("mean methyaltion") +
    xlab("tumor cells (self sorted)") +
    ggtitle(patient.name)
  
  # save plot
  out.path = paste0(patient.name, '_cell_methylation_by_', color.by, '.png')
  ggsave(out.path)
}

####################################################################################################
#                                            main                                                  #
####################################################################################################

my_colour = list(
  lesion = c(NC = "#6FBE44", PT = "#F380A3", LN = "#834D40", MO = "#D1A349", ML = "#1D489E", 
             MP = "#815BA6"),
  region = c(NC = "#6FBE44", PT1 = "#F57E32", PT2 = "#3C8A45", PT3 = "#e4373b", PT4 = "#53ABDA",
             PT5 = "#EBE94F", PT6 = "#EFB0BC", LN1 = "#B6B6BA", LN2 = "#B28A8E", LN3 = "#866c6f",
             LN4 = "#67532B", LN5 = "#514321", ML1 = "#94D6E4", ML2 = "#4872B7", ML3 = "#1C49A0", 
             ML4 = "#333463", ML5 = "#464B7D", ML6 = "#3C3E69", MP1 = "#CFB0D3", MP2 = "#BC85BB", 
             MP3 = "#8254A2", MP4 = "#842F8D", MP5 = "#632E65", MO1 = "#E7CF9F", MO2 = "#E2C68E", 
             MO3 = "#E2C68E", MO4 = "#D9B36B", MO5 = "#D5AB5B", MO6 = "#D1A349"),
  sublineage = c(NC = "#6FBE44", A0 = "#F8B4C0", A1 = "#E0E346", A2 = "#6F51A1", A3 = "#F89C31",
                 A4 = "#EF292A", A5 = "#A45AA4",
                 A6 = "#993232", A7 = "#2256A6", A8 = "#BC84BA", A9 = "#ED3095", B = "#3B86C6",
                 B0 = "#3B86C6", B1 = "#66CAD4", B2 = "#6D8841", B3 = "#28898A", C0 = "#E7CFA0", 
                 C1 = "#DDBD7E", C2 = "#D1A34A", C3 = "#B89979", C4 = "#AA845D", C5 = "#8C6E4A", 
                 undefined = "#BBBABC")
)
custom_labels = c(NC = "NC: Normal Cell", PT = "PT: Primary Tumor", LN = "LN: Lymph Node Metastasis", MO = "MO: Omental Metastasis",
                  ML = "ML: Liver Metastasis", MP = "MP: Post-treatment Liver Metastasis")

methylation.path <- "avg_data_all_mean_coverage.csv"

data <- read.data(methylation.path)
create.plot(data, 'CRC13', 'sublineage', my_colour, custom_labels)


print('done!')
