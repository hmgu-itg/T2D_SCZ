firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  return(x)
}

setwd("C:/Users/ana.arruda/OneDrive - Helmholtz Zentrum München/Projects/SCZ_T2D/")

my_colors <- RColorBrewer::brewer.pal(9, "RdPu")[c(5,7,9)]

################### ConsensusPathDb #####################
dt <- data.table::fread("pathway_analysis/consensus_likely_level5.csv")

dt <- dt[, term_name:=firstup(term_name)]
dt <- dt[order(`q-value`)]
lvl=5

ggplot2::ggplot(dt[c(1:7,10,14)], aes(x=-log(`p-value`), y=reorder(term_name, -log(`p-value`), sum))) +
  geom_col(fill=my_colors[1], alpha=0.75) +
  theme_bw() + 
  ylab("") +
  theme(axis.text = element_text(colour = "black", size = 15),
        axis.title = element_text(colour = "black", size = 15),
        legend.position = "none")+
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 60)) +
  xlab("-log(p-value)")

ggplot2::ggsave(paste0("paper/figures/pathway_analysis_likely_", lvl, ".jpg"), width=8, height=2.7, dpi=600)
