library(ggplot2)
library(patchwork)

layout <- "
##AAAAAAA
CCBBBBBBB
CCBBBBBBB
CCBBBBBBB
CCBBBBBBB
CCBBBBBBB
CCBBBBBBB
"

heatmap_plot_legend <- heatmap_plot + cowplot::draw_plot(legend, 80, .71, 0, 120)

patched_plot <- movement_plot + heatmap_plot_legend + dendrogram_plot + plot_layout(design = layout)

ggsave("C:\\Users\\lukas\\Projects\\SARS_COV_2_Run3\\manuscript\\plots\\Fig2.png", patched_plot,
       width = 10,
       height = 10,
       dpi = 300)