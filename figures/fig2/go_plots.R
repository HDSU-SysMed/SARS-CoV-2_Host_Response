#----------------------------------------#
#
#   GO Plot | Time plot over GO Terms
#
#----------------------------------------#

library(tidyverse)
library(ggplot2)

data.l2fc.go <- read.csv('Data_l2fc_bgcorr_Host_Exp3_interpolation_fine_0to48h_201217_GO_Terms.csv')[, -c(1)]

# Get list of top 15 GO terms 

data.l2fc.go.convert <- data.l2fc.go %>% mutate(P.DE_LOG = as.numeric(-log10(P.DE)))

#----------------------------------------#

data.l2fc.go.convert <- na.omit(data.l2fc.go.convert)
data.l2fc.go.convert <- data.l2fc.go.convert %>% arrange(t, desc(P.DE_LOG))

#----------------------------------------#

top.terms <- (data.l2fc.go.convert %>%
                group_by(t) %>% 
                filter(row_number() <= 15))$Term
top.terms <- unique(top.terms)

data.l2fc.go.subset <- data.l2fc.go.convert %>% subset(Term %in% top.terms)

#----------------------------------------#

library(GO.db)

data.l2fc.go.subset.go.terms <- unique(data.l2fc.go.subset$Term)

all.go.terms <- as.list(GOTERM)
all.go.terms.names <- as.data.frame(lapply(all.go.terms, function(x) x@Term))
all.go.terms.names <- as.data.frame(t(all.go.terms.names))
all.go.terms.names <- all.go.terms.names %>% rownames_to_column()
colnames(all.go.terms.names) <- c("GO_ID", "GO_Term")

data.l2fc.go.subset.go.terms <- all.go.terms.names %>% subset(GO_Term %in% data.l2fc.go.subset.go.terms)

#----------------------------------------#

# data.l2fc.go.subset.go.terms.filtered.revigo <- read.csv("go/data.l2fc.go.subset.go.terms.filtered.revigo.csv")[,-c(1)]
# data.l2fc.go.subset.go.terms.filtered.revigo <- data.l2fc.go.subset.go.terms.filtered.revigo %>% subset(eliminated == 0)

#data.l2fc.go.subset.revigo <- data.l2fc.go.subset[data.l2fc.go.subset$Term %in% data.l2fc.go.subset.go.terms.filtered.revigo$description,]

data.l2fc.go.subset <- data.l2fc.go.subset.go.terms
data.l2fc.go.subset <- na.omit(data.l2fc.go.subset)

#----------------------------------------#

x <- data.l2fc.go.subset %>%
  pivot_wider(id_cols = "Term",
              names_from = "t", 
              values_from = "P.DE_LOG"
  ) %>% column_to_rownames(var = "Term") 

xx <- as.data.frame(lapply(x, function(e) as.data.frame(e)))
rownames(xx) <- rownames(x)
colnames(xx) <- colnames(x)

x[is.na(x)] <- 0

dd.row <- as.dendrogram(hclust(d = dist(x = x)))

dd.order <- order.dendrogram(dd.row)

data.l2fc.go.subset.ordered <- data.l2fc.go.subset 


data.l2fc.go.subset.ordered$Term <- factor(x = data.l2fc.go.subset.ordered$Term ,
                                           levels = rownames(x)[dd.order], 
                                           ordered = TRUE)
data.l2fc.go.subset.ordered$t <- as.numeric(data.l2fc.go.subset.ordered$t)


#----------------------------------------#

library(dendextend)
library(wesanderson)

theme_set(ggthemes::theme_base())

dendro_data_row <- ggdendro::dendro_data(dd.row, type = "rectangle")

dendrogram_plot <- ggplot() + 
  geom_segment(data=dendro_data_row$segments, aes(x=x, y=y, xend=xend, yend=yend), size=0.35) + 
  #geom_text(data=dendro_data_row$labels, aes(x=x, y=y-5, label=label, hjust=0), family="Calibri", size=3) +
  coord_flip() + 
  scale_y_reverse(limits=c(100, 0), expand=c(0, 0)) +
  scale_x_discrete(expand=c(0.005, 0)) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        plot.margin = unit(c(0,0,0,0), "cm"))

dendrogram_plot

#----------------------------------------#







library(RColorBrewer)
library(extrafont)

#colors <- c('#ffffff', '#ffda9e', '#f9c566', '#f3b941')
colors <- c('#FFFFFF','#F2C6B2','#E68C66','#D95319')

#colors <- c('#ffffff', '#f7b799', '#e27b62', '#cd4f45', '#b41c2d')

data.l2fc.go.subset.ordered <- data.l2fc.go.subset.ordered %>% mutate(P.DE_LOG_transform = sqrt(P.DE_LOG))

theme_set(ggthemes::theme_base())

# Heatmap upregulated genes
heatmap_plot <- ggplot(data = data.l2fc.go.subset.ordered, aes(t, Term, fill=P.DE_LOG)) +
  geom_tile() +
  xlab("Hours (post infection)") + 
  ylab("") + 
  scale_x_continuous(expand=c(0, 0), limits = c(0.5, 47.5)) +
  scale_y_discrete(position="right") + 
  scale_fill_gradientn(
    colors = colors
  ) +
  theme(axis.line = element_blank(),
        axis.title.x = element_text(size = 10, color="black", margin = margin(t = 5)),
        axis.title.y = element_text(size = 10, color="black", margin = margin(r = 5)),
        axis.text.x=element_text(size = 10, color="black", margin = margin(t = 5)),
        axis.text.y=element_text(size = 10, color="black", margin = margin(r = 5)),
        text = element_text(size = 11),
        legend.position = "none",
        plot.margin = unit(c(0,0,0,0), "cm"),
        panel.background=element_blank(),
        panel.border = element_rect(size=0.7),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())


heatmap_plot


#------------------------------------------------#

# Heatmap upregulated genes
legend_plot_extract <- ggplot(data = data.l2fc.go.subset.ordered, aes(t, Term, fill=P.DE_LOG)) +
  geom_tile() +
  xlab("Hours post infection") + 
  ylab("") + 
  theme_bw() + 
  theme(
    axis.line = element_line(colour = "black"),
    panel.background = element_blank(),
    text = element_text(size = 10)
  ) + 
  scale_x_continuous(expand=c(0, 0), limits = c(0.5, 47.5)) +
  scale_y_discrete(position="right") + 
  scale_fill_gradientn(
    colors = colors
  ) + labs(fill = "-log10 p-value") + 
  guides(fill = guide_colourbar(barheight = unit( 0.1 , "in" ),
                                    ticks.colour = "black",
                                    ticks.linewidth = 1,
                                    frame.colour = "black",
                                    direction = "horizontal", 
                                    title.hjust=0.5,
                                    title.position = "top",
                                    frame.linewidth = 1)) 

legend_plot_extract

legend <- cowplot::get_legend(legend_plot_extract)

cowplot::ggdraw(legend)

#------------------------------------------------#

# Assemble final plot

final_plot_F <- cowplot::insert_yaxis_grob(final_plot, dendro_row, width = grid::unit(2, "cm"), position = "left")

final_plot_F <- cowplot::ggdraw(final_plot_F) + cowplot::draw_plot(legend, .65, .71, .5, .5)

final_plot_F

#------------------------------------------------#

ggsave(filename = "C:\\Users\\lukas\\Projects\\SARS_COV_2_Run3\\manuscript\\plots\\GO_Terms_Fig.svg", 
       width = 10,
       height = 8,
       scale = 1,
       dpi=700)