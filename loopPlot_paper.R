# install.packages("cli")
# install.packages("devtools")
# devtools::install_github("matherealize/looplot")
library("devtools")
library("looplot")
library(dplyr)
res_selct$na_percnt<- format(res_selct$na_percnt,digits = 2)
res_selct$na_percnt[res_selct$na_percnt==51]<-50

resdf_w<-select(res_selct,na_percnt,rho0,r,m,p,Effloss)
names(resdf_w)[4] <- "Subject"
names(resdf_w)[5] <- "Period"

#change the order of cac
resdf_w$r= factor(resdf_w$r, levels=c('1','0.95','0.8'))
  
#pdf(file="looplot_paper2.pdf",width = 15, height = 9)
add_name = "G:\\Shared drives\\Ehsan PhD work\\Paper_1\\plots\\"
png(paste0(add_name,"looplot_paper1",".png"), width = 1500, height = 1000)

#Sample Size\nicc=0.05, cac_HH=1, cac_BE or cac_DTD=0.8
nested_loop_plot(resdf = resdf_w, 
                 x = "na_percnt", steps = c("r","rho0"),
                 grid_rows = "Subject", grid_cols = "Period",
                 #draw = "add_steps", 
                 steps_y_base = -15, steps_y_height = 10, steps_y_shift =12,
                 x_name = "Removal percentage of cluster-period cells (%)", y_name = "Efficiency loss (%)",
                 spu_x_shift = 60,
                 # connect_spus = TRUE,
                 point_shapes= c(15,19),
                 colors = scales::brewer_pal(palette = "Set1"),
                 line_linetypes = c(2,3),
                 steps_color = "dimgrey",
                 steps_values_annotate = TRUE, steps_annotation_size = 7,
                 steps_annotation_nudge = 0.3,
                 steps_annotation_color= "dimgrey",
                 hline_intercept = c(0,100), 
                 y_expand_add = c(20, NULL), 
                 #legend_name = "",
                 steps_names=c("CAC","ICC"),
                 sizes =2,base_size = 15,
                 post_processing = list(
                   add_custom_theme = list(
                     axis.text.x = element_text(angle = -90, 
                                                vjust = 0.5, 
                                                size = 19),     
                      # axis.text = element_text(size = 15,face='bold'),
                      axis.title.x = element_text(size = 20,face='bold'),
                      axis.title.y = element_text(size = 20,face='bold'),
                      legend.title= element_blank(), # element_text(size=15,face='bold'),
                      legend.text = element_blank(),  # element_text(size=18,face='bold'),
                      legend.position="none",
                      strip.text.y =element_text(size = 20,face='bold'),
                      strip.text.x =element_text(size = 20,face='bold')
                   )))


dev.off()

