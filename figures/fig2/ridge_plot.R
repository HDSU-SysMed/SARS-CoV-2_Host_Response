##########################################
##
## Ridge Plot | Distribution of LFC
##
##########################################

library(xlsx)
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(ggridges)

data.l2fc <- read.csv('Data_l2fc_bgcorr_Host_Exp3_interpolation_201217.txt', sep='\t')#[,c(-1)]

data.l2fc <- data.l2fc %>% pivot_longer(!geneid, names_to = c("stat_label"))

check_stat_label <- function(x){
  if(startsWith(x, "X")){
    return("mean")
  } else {
    return("std")
  }
}

data.l2fc <- data.l2fc %>% mutate(., 
                     label=map(stat_label, check_stat_label),
                     t=map(stat_label, function(x) gsub("h", "", 
                                                  gsub("sem.", "", 
                                                  gsub("X", "", x)))
                     ))

data.l2fc$t <- as.numeric(data.l2fc$t)
data.l2fc <- data.l2fc[, !names(data.l2fc) %in% c("stat_label")]

##########################################

trace(grDevices:::png, quote({
  if (missing(type) && missing(antialias)) {
    type <- "cairo-png"
    antialias <- "subpixel"
  }
}), print = FALSE)

##########################################

data.l2fc %>% subset(label == 'mean') %>% group_by(t) %>% summarise(max=max(value), min=min(value))


require(scales)
library(gridExtra)
library(ggbeeswarm)


data.l2fc.toplot <- data.l2fc %>% subset(t != 0) %>% subset(label == 'mean')

qt_cut <- 0.005

vert_lines_left <- data.l2fc.toplot %>% group_by(t) %>%
                            summarize(x=quantile(value, qt_cut),
                                   xend=x,
                                   y=mean(t), 
                                   yend=1.5*y
                                   )

vert_lines_right <- data.l2fc.toplot %>% group_by(t) %>%
  summarize(x=quantile(value, 1-qt_cut),
            xend=x,
            y=mean(t), 
            yend=1.5*y
  )

theme_set(ggthemes::theme_base())

#----------------------------------------------------------#


qt_cut <- 0.0001

vert_lines_down <- data.l2fc.toplot %>% group_by(t) %>%
  summarize(ymin=quantile(value, qt_cut),
            ymax=quantile(value, qt_cut),
            )

vert_lines_top <- data.l2fc.toplot %>% group_by(t) %>%
  summarize(ymin=quantile(value, 1-qt_cut),
            ymax=quantile(value, 1-qt_cut),
  )

movement_plot <- ggplot(data.l2fc.toplot, aes(x=factor(t), y=value, fill=factor(1))) + 
  geom_hline(yintercept = c(-1, 1), linetype="dashed") +
  geom_violin(trim=FALSE, size=0.3) +
  geom_errorbar(data = vert_lines_top, 
                aes(x=factor(t), ymin = ymin, ymax = ymax), 
                width = 0.2,
                size=1.5,
                color="#a2142f",
                inherit.aes = FALSE) +
  geom_errorbar(data = vert_lines_down, 
                aes(x=factor(t), ymin = ymin, ymax = ymax), 
                width = 0.2,
                size=1.5,
                color="#0072bd",
                inherit.aes = FALSE) +
  stat_summary(
    fun.data = "mean_sdl",  fun.args = list(mult = 1), 
    geom = "pointrange", color = "black", size=0.1
  ) +
  scale_fill_manual(values=c('grey')) + 
  scale_y_continuous(position = "right") + 
  labs(x = "Hours (post infection)", y = "Log2-Fold Change") + 
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 10, color="black", margin = margin(t = 5)),
        axis.title.y = element_text(size = 10, color="black", margin = margin(r = 5)),
        axis.text.x=element_text(size = 10, color="black", margin = margin(t = 5)),
        axis.text.y=element_text(size = 10, color="black", margin = margin(r = 5)),
        plot.margin = unit(c(0,0,0.5,0), "cm"),
        panel.background=element_blank(),
        panel.border = element_rect(color="black", size=0.7),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank()
        )

movement_plot


ggsave("C:\\Users\\lukas\\Projects\\SARS_COV_2_Run3\\manuscript\\plots\\Fig2_A.png", p,
       width = 5,
       height = 2,
       dpi = 300)