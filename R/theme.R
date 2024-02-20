# Author:       Stephen-John Sammut
# Description:  ggplot theme for manuscript
#=========================================================================

library (ggplot2)
library (ggpubr)

theme_manuscript <- function(base_size=12, base_family="arial") {
  library(grid)
  (ggthemes::theme_foundation(base_size=base_size)
    + theme(plot.title = element_text(face = "bold",hjust = 0.5,size = base_size),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title.y = element_text(angle=90,vjust =2,size = base_size),
            axis.title.x = element_text(vjust = -0.2,size = base_size),
            axis.text = element_text(size = base_size), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.5, "cm"),
            legend.spacing = unit(0, "cm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(size=base_size)
    ))
}

stat_p_eq     <- stat_poly_eq(aes(label = paste(..rr.label..)), color="#2F4858",
                              label.x.npc = "left", label.y.npc = 0.99,
                              formula = y~x, parse = TRUE, size = 1.8)
stat_f_glance <- stat_fit_glance(method = "lm", label.y = 0.86, label.x = 0.05,
                                 method.args = list(formula = y~x),
                                 aes(label = sprintf('italic(P)~"="~%.30f',
                                                     as.numeric(as.character(
                                                       formatC(stat(p.value),
                                                               format = "e", 
                                                               digits = 1))))), 
                                 size = 1.8,parse = TRUE)
  
cols.StemCladePrivate <- c("#4281A4","#F17105","#60A561")

cols.organ<- c("#33a02c","#fb9a99","#e31a1c","#ff7f00","#1f78b4","#a6cee3","#6a3d9a")
names(cols.organ)<-c("Bone","Breast","Liver","Lymph node","Lung/pleura","Brain/meninges","Pericardium")

cols.cloneClass<-c("#EEC643","#79B4A9", "#C06E52","#D5CAD6")
names(cols.cloneClass)<-c("C","D","A","B")



