# install.packages("rsconnect")
# library("rsconnect")

library("shiny")
library("ggplot2")
library("reshape2")
library("plyr")
library("swCRTdesign")
library("matrixcalc")
library("scales")
library("tidyverse")
library("shinythemes")
library("Matrix")
library("plotly")
library("RColorBrewer")
library("gganimate")
library("gifski")
library("tidyr")
library("animation")
library("htmlwidgets")
library("plotly")
library("rsconnect")
require("gridExtra")

#setwd("~/Google Drive/Shared drives/Ehsan PhD work/Codes/")
setwd("G:\\Shared drives\\Ehsan PhD work\\Codes\\Git\\Iterative-Removals-SW\\IterativeInfContent")
source("CRTVarAdj_func.R", local=TRUE)
source("IterRemoval_func.R", local=TRUE)
#source("ICcell_appfunc.R", local=TRUE)
# # Functions for generating design matrices
SWdesmat <- function(Tp) {
  Xsw <- matrix(data=0, ncol = Tp, nrow = (Tp-1))
  for(i in 1:(Tp-1)) {
    Xsw[i,(i+1):Tp] <- 1
  }
  return(Xsw)
}

# #######

###############make all of them the same shape! last design (information content matrix) should be looking empty for all
###############e.g. m=20 rho0=0.05 r=0.95 type=1 T=7
###############Add this adjustment to codes (think about it)
################NO NEED TO ADJUST AS THE FIGURE IS SHOWN NA FOR EVERY SITUATION

# Xdlist
# dlist
# (K+1-k,T+1-t)
m=50
rho0=.01
r=0.8
# #lbl="Hussey and Hughes model"
# #lbl="Exponential decay model"
# #1 decay 0 HH
type=1
Tp=5
# # 0 removing without any consideration, 1 removing only one pair at each step
cpnum=1
effsize=0.20

#m <- res_m[which(res_m$na_percnt)==median(res_m$na_percnt), ]
#I couldn't understand why I used this line below in the next graph, I need to check it later
# linetype=~as.factor(r)
# 
# 
# cac_type <- unique(res_m$r)
# linetypes <- setNames(c('solid', 'dot', 'dash', 'longdash', 'dashdot'), cac_type)
# 
# plot = function(data, m, T,rho0) {
##relative variance plots accross cac
# plot_ly(data = res_m, x = ~na_percnt, y = ~Effloss, type="scatter",linetype=~as.factor(r),linetypes=linetypes,
#         mode="lines", hoverinfo="text", hoverlabel=list(bordercolor=NULL, font=list(size=14)),
#        # line=list(dash = "dash",color="#FF6666"),
#         text=~paste("RPercent:", round(na_percnt,2),"%", "<br>Effloss:", round(Effloss,2),"%",
#         "<br>Power:", format(round(power,4)*100,2),"%")) %>%
#         layout(xaxis=list(title="Removal percentage of cluster-period cells (%)", titlefont=list(size=18), showline=TRUE,
#         tickmode="auto", tickfont=list(size=16), nticks=11, range=c(0,100),ticks="inside",
#         mirror=TRUE, showgrid=FALSE),
#         yaxis=list(title="Efficiency loss (%)", titlefont=list(size=18), tickfont=list(size=16),
#                     mirror=TRUE, showline=TRUE,range=c(0,100)),legend = list(title=list(text='<b> cac </b>')))#%>%
#         # title=list(text=paste0("\n","m=",m,","," ","T=",T,","," ","icc","=",rho0,","," ","effsize=",effsize),
#         #                      y =0.95,font=list(size = 15))) %>%
#          #add_segments(y = 0, yend = 50, x = 50, xend = 50,line = list(dash = "dash",color="#00BFFF"), showlegend=FALSE,hoverinfo='skip')
#          # title=list(text=paste(if (type == 1) {paste0("plot"," ",l," ","(","Discrete time decay",")")
#          # } else if (type == 0){paste0("plot"," ",l," ","(","Exchangeable",")")},
#          # paste0("\n","m=",m,","," ","T=",T,","," ","icc","=",rho0,","," ","effsize=",effsize))
#          #                      ,y =0.95,font=list(size = 15)))
# 
# }
# 
# 
# p=plot(res_m,m,T,rho0)
# save_image(p, file="Eff_loss_plot.png",width = 1150, height = 550)

#single line efficiency loss or power
# plot = function(data, m, T,rho0) {
# 
#   plot_ly(data = res_m, x = ~na_percnt, y = ~power, type="scatter",name = "power", mode="lines",
#           hoverinfo="text", hoverlabel=list(bordercolor=NULL, font=list(size=14)),line=list(dash ="dash",color="#00CC00", width = 6),
#           text=~paste("RPercent:", round(na_percnt,2),"%", "<br>Effloss:", round(Effloss,2),"%",
#                       "<br>Power:", round(power,2),"%")) %>%
#   layout(xaxis=list(title="Removal percentage of cluster-period cells (%)", titlefont=list(size=22,face="bold"), showline=TRUE,
#                                 tickmode="auto", tickfont=list(size=16), nticks=11, range=c(0,100),ticks="inside",
#                               mirror=TRUE, showgrid=FALSE),
#         yaxis=list(title="Power (%)", titlefont=list(size=20,face="bold"), tickfont=list(size=16),
#                    mirror=TRUE, showline=TRUE,range=c(0,100)),legend = list(title=list(text='<b> Type </b>')))%>%
#    add_segments(y = 0, yend = 50, x = 50, xend = 50,line = list(dash = "dash",color="#00BFFF"), showlegend=FALSE,hoverinfo='skip')
# }
# 
# p=plot(res_m,m,T,rho0)
# #save_image(p, file="efflossPos.png",width = 1150, height = 550)
# save_image(p, file="pwr.png",width = 1150, height = 550)
# 
# ay <- list(
#   tickfont = list(color = "#00CC00"),
#   overlaying = "y",
#   side = "right",
#   title = "<b>Power (%)</b>")
# 
# ###add power
# plot = function(data, m, T,rho0) {
#   #relative variance plots accross cac
#   plot_ly(data = res_m, x = ~na_percnt, y = ~Effloss, type="scatter",name = "Effloss", mode="lines",
#          hoverinfo="text", hoverlabel=list(bordercolor=NULL, font=list(size=14)),line=list(dash ="dash",color="#FF6666"),
#           text=~paste("RPercent:", round(na_percnt,2),"%", "<br>Effloss:", round(Effloss,2),"%",
#                       "<br>Power:", round(power,2),"%")) #%>%
#   # add_trace(data = res_m, x = ~na_percnt, y = ~power, type="scatter",name = "power",yaxis = "y2",mode="lines",
#   #           line=list(dash ="dot",color="#00CC00")) %>%
#   # layout(yaxis2 = ay,xaxis=list(title="Removal percentage of cluster-period cells (%)", titlefont=list(size=18), showline=TRUE,
#   #                               tickmode="auto", tickfont=list(size=16), nticks=11, range=c(0,100),ticks="inside",
#   #                               mirror=TRUE, showgrid=FALSE),
#   #        yaxis=list(title="Efficiency loss (%)", titlefont=list(size=18), tickfont=list(size=16),
#   #                   mirror=TRUE, showline=TRUE,range=c(0,100)),legend = list(title=list(text='<b> Type </b>')))%>%
#   #   # title=list(text=paste0("\n","m=",m,","," ","T=",T,","," ","icc","=",rho0,","," ","effsize=",effsize),
#   #   #                      y =0.95,font=list(size = 15))) %>%
#   #   add_segments(y = 0, yend = 50, x = 50, xend = 50,line = list(dash = "dash",color="#00BFFF"), showlegend=FALSE,hoverinfo='skip')
#     # title=list(text=paste(if (type == 1) {paste0("plot"," ",l," ","(","Discrete time decay",")")
#     # } else if (type == 0){paste0("plot"," ",l," ","(","Exchangeable",")")},
#     # paste0("\n","m=",m,","," ","T=",T,","," ","icc","=",rho0,","," ","effsize=",effsize))
#     #                      ,y =0.95,font=list(size = 15)))
# }
# 
# p=plot(res_m,m,T,rho0)
# save_image(p, file="efflossPow_paper1.png",width = 1150, height = 550)
# 
# ###Paper, Fig5 & with some modification Fig 7
# ####two separate icc (0.15, 0.14) to do this generate two datasets

# add_name = "G:\\Shared drives\\Ehsan PhD work\\Paper_1\\plots\\Fig 7"
# 
# 
# #efficiency loss and power in one figure (Fig 5)
# fig1<- plot_ly(data = res_m, x = ~na_percnt, y = ~Effloss, type="scatter",name = "Discrete-time Decay", mode="lines",
#         hoverinfo="text", hoverlabel=list(bordercolor=NULL, font=list(size=14)),line=list(dash ="dot",color="#8b0000"),
#         #,line=list(dash = "dash",color="#FF4500")#013220 darkgreen
#         text=~paste("RPercent:", round(na_percnt,2),"%", "<br>Effloss:", round(Effloss,2),"%",
#                     "<br>Power:", round(power,2),"%")) %>%
#   layout(xaxis=list(title="Removal percentage of cluster-period cells (%)", titlefont=list(size=18), showline=TRUE,
#                     tickmode="auto", tickfont=list(size=16), nticks=11, range=c(0,100),ticks="inside",
#                     mirror=TRUE, showgrid=FALSE),
#          yaxis=list(title="Efficiency loss (%)", titlefont=list(size=18), tickfont=list(size=16),
#                     mirror=TRUE, showline=TRUE,range=c(0,100)),legend = list(title=list(text='<b>Type</b>')))%>%
#   # title=list(text=paste0("\n","m=",m,","," ","T=",T,","," ","icc","=",rho0,","," ","effsize=",effsize),
#   #                      y =0.95,font=list(size = 15))) %>%
# add_segments(y = 100, yend = 100, x = 0, xend = 100,line = list(dash = "solid",color="black"), showlegend=FALSE,hoverinfo='skip')%>%
# # title=list(text=paste(if (type == 1) {paste0("plot"," ",l," ","(","Discrete time decay",")")
# # } else if (type == 0){paste0("plot"," ",l," ","(","Exchangeable",")")},
# # paste0("\n","m=",m,","," ","T=",T,","," ","icc","=",rho0,","," ","effsize=",effsize))
# #                      ,y =0.95,font=list(size = 15)))
#  add_trace(data = res_m1, x = ~na_percnt, y = ~Effloss, type="scatter",name = "Exchangeable",mode="lines",
#            line=list(dash ="solid",color="#8b0000"))
# #save_image(p1, file="effloss_paper1_fig6.png",width = 1150, height = 550)
# 
# fig2<-plot_ly(data = res_m, x = ~na_percnt, y = ~power, type="scatter",name = "Discrete-time Decay", mode="lines",
#           hoverinfo="text", hoverlabel=list(bordercolor=NULL, font=list(size=14)),line=list(dash ="dot",color="#8b0000"),
#           #,line=list(dash = "dash",color="#FF4500")#013220 darkgreen
#           text=~paste("RPercent:", round(na_percnt,2),"%", "<br>Effloss:", round(Effloss,2),"%",
#                       "<br>Power:", round(power,2),"%"),showlegend=F) %>%
#     layout(xaxis=list(title="Removal percentage of cluster-period cells (%)", titlefont=list(size=18), showline=TRUE,
#                       tickmode="auto", tickfont=list(size=16), nticks=11, range=c(0,100),ticks="inside",
#                       mirror=TRUE, showgrid=FALSE),
#            yaxis=list(title="Power (%)", titlefont=list(size=18), tickfont=list(size=16),
#                       mirror=TRUE, showline=TRUE,range=c(0,100)),legend = list(title=list(text='<b>Type</b>'))) %>%
#     add_trace(data = res_m1, x = ~na_percnt, y = ~power, type="scatter",name = "Exchangeable",mode="lines",
#               line=list(dash ="solid",color="#8b0000"),showlegend=F)
# 
# fig <- subplot(fig1, fig2, shareX=TRUE, titleY = TRUE,nrows=2) %>%
#   layout(title = '',showlegend = TRUE)
# 
# save_image(fig, file="efflossPow_paper1_fig5.png",width = 1150, height = 550)



l=0
cac = c(0.8,0.95)
#cac = 0.95
#tp = c(0,1)

#LOOPLOT
res_selct <- data.frame(iter=integer(),
                    variance=integer(),
                    power=integer(),
                    r=integer(),
                    RVariance=integer(),
                    Effloss=integer(),
                    na_percnt=integer(),m=integer(),p=integer(),rho0=integer())
p <-list ()
for (m in c(10,100)){
  for (T in c(5,10)){
    for (rho0 in c(0.01,0.05,0.15)){

      l=l+1
      
#HTML Efficiency loss plots for different settings 
      #datasets for fig 5
res_m <- data.frame(iter=integer(),
                      variance=integer(),
                      power=integer(),
                      r=integer(),
                      RVariance=integer(),
                      Effloss=integer(),
                      na_percnt=integer(),m=integer(),p=integer(),rho0=integer())

for (type in c(0,1)){

  for (r in  if (type ==0){r = 1}  else {r = cac})
    {
 
#LOOPLOT       
res_sb <- data.frame(iter=integer(),
                         variance=integer(),
                         power=integer(),
                         r=integer(),
                         RVariance=integer(),
                         Effloss=integer(),
                         na_percnt=integer(),m=integer(),p=integer(),rho0=integer())
    
         #Put parameters here
         melted_varmatexcl<- melt(IterRemoval(Tp,m, rho0, r, type,cpnum)[[1]])
         
         #varmat_excl$value<-round(varmat_excl$value, 4)
         melted_desmatexcl<-  melt(IterRemoval(Tp,m, rho0, r, type,cpnum)[[2]])
         names(melted_desmatexcl)[names(melted_desmatexcl)=="value"] <- "Xdvalue"
         melted_varmatexcl_t<- merge(melted_varmatexcl, melted_desmatexcl,all.y=TRUE,by = c('Var1','Var2','L1'))
         melted_varmatexcl_t$value<-round(melted_varmatexcl_t$value, 4)

        #color_palette <-colorRampPalette(c( "yellow", "red"))(length(table(varmat_excl)))

        pal <- colorRampPalette(brewer.pal(8, "YlOrRd"))(length(unique(melted_varmatexcl_t$value))-1)

        names(melted_varmatexcl_t)[names(melted_varmatexcl_t)=="Var1"] <- "Sequence"
        names(melted_varmatexcl_t)[names(melted_varmatexcl_t)=="Var2"] <- "Period"
        names(melted_varmatexcl_t)[names(melted_varmatexcl_t)=="L1"] <- "iter"
        
        Xdes <- SWdesmat(Tp)
        varmatall<- IterRemoval(Tp,m, rho0, r, type,cpnum)[[3]]
          
        Tp <- ncol(Xdes)
        K <- nrow(Xdes)

        ###Need power values for removal plot####
        iter=1:length(varmatall)
        df=as.data.frame(varmatall)

        pow <- function(vars, effsize, siglevel=0.05){
          z <- qnorm(siglevel/2)
          pow <- pnorm(z + sqrt(1/vars)*effsize)
          return(pow)
        }

        # Calculate power for a set of variances, a given effect size and sig level
        powdf <- function(df, effsize, siglevel=0.05){
          powvals <- apply(df, MARGIN=2, pow, effsize, siglevel)
          powdf <- data.frame(iter, df$varmatall,powvals*100)
          colnames(powdf) <- c("iter","variance","power")
          return(powdf)
        }

        res <- powdf(df,effsize)
        res$r <- r
   
        res <- cbind(res,res$variance[1]/res$variance,(1-(res$variance[1]/res$variance))*100)
        colnames(res) <- c("iter","variance","power","r","Rvariance","Effloss")
        
        #count NA Xdvalue
        #res_na <- melted_varmatexcl_t %>% group_by(iter) %>% summarise(na_percnt = (sum(is.na(value))/(T*K))*100)
        res_na <- melted_varmatexcl_t %>% group_by(iter) %>% summarise(na_percnt = (sum(is.na(Xdvalue))/(T*K))*100)
        res_2 <- merge(res, res_na,by = c('iter'))
        
        res_2$m <- m
        res_2$p <- Tp
        res_2$rho0 <- rho0
        
        res_m <- rbind(res_m,res_2)
        #select for looplot
        res_sb<-res_2[c(which.min(abs(20-res_2$na_percnt)),which.min(abs(50-res_2$na_percnt)),which.min(abs(80-res_2$na_percnt))),]
        res_selct <-rbind(res_selct,res_sb)

        #end of loop for type and cac and type
  }
}
      #save plots

     #p[[l]]= plot(res_m,m,T,rho0)
      #print(i)
      #htmlwidgets::saveWidget(p,paste0("G:/Shared drives/Ehsan PhD work/Outputs/RV/RV_","m",m,"_","T",T,"_","icc"," ",rho0,".html"))
      #end of loop for desired settings

    }
  }
}

# p

#heatmap plot
# 
#     melted_varmatexcl_t <- merge(res, melted_varmatexcl_t, by = "iter", all = TRUE)
# 
# 
#     p<-ggplot(melted_varmatexcl_t,aes(Period,Sequence,frame=iter))+
#       geom_tile(aes(fill=factor(value)),colour = "grey50") +
#       scale_y_reverse(breaks=c(1:K)) +
#       scale_x_continuous(breaks=c(1:T)) +
#       theme(panel.grid.minor = element_blank()) +
#       theme(legend.position="none")+
#       # geom_label(data = melted_varmatexcl_t,aes(label= round(value,4)),fontface = "bold",
#       #          colour = "white",size = 4) +
#       geom_text(aes(label= Xdvalue),color = "black",size = 7, check_overlap = T)+
#       #scale_fill_distiller
#       scale_fill_manual(values = pal,breaks=levels(melted_varmatexcl_t$value)[seq(90, 150, by=5)],na.value="gray")
# 
#     p1<-ggplotly(p) %>%
#       animation_opts(frame = 500,transition =0,redraw = TRUE) %>%
#       animation_slider(currentvalue = list(prefix = "Iter: ", font = list(color="darkblue"))) %>%
#       partial_bundle(local = FALSE)%>%
#       animation_button(visible = T)%>%
#     animation_slider(visible =T)
#     #error seems ok: this is for geom_text
# p1
# 
# htmlwidgets::saveWidget(as_widget(p1), "index.html")
#  
# 
#  
# ####PRESNTATION
# ####PAPER1
#plot_list =list()
# ####PRESNTATION
# #add_name = "G:\\Shared drives\\Ehsan PhD work\\Publications\\43rd ISCB\\Figures_2\\"
# ####PAPER1
#add_name = "G:\\Shared drives\\Ehsan PhD work\\Paper_1\\plots\\Fig 3\\"
###Confirmation report
#melted_varmatexcl_t <- merge(res, melted_varmatexcl_t, by = "iter", all = TRUE)
# add_name = "G:\\Shared drives\\Ehsan PhD work\\confirmation\\Fig 2"
# png(paste0(add_name,"_t",".png"), width = 2000, height = 1000)
# titl=c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(l)")
# 
# for (i in unique(melted_varmatexcl_t$iter) ) {
# 
#    # png(paste0(add_name,i,".png"), width = 500, height = 500)
# 
#     p<-ggplot(subset(melted_varmatexcl_t,iter==i),aes(Period,Sequence))+
#     geom_tile(aes(fill=factor(value)),colour = "grey50") +
#     geom_text(x=0.5,y=-0.35,hjust=0,aes(label=paste0("Power:",format(round(power[[1]],2),2),"%")),
#               size=4)+
#     scale_y_reverse(breaks=c(1:K)) +
#     scale_x_continuous(breaks=c(1:T)) +
#       labs(title = titl[[i]],caption =  paste("iteration:",i))+
#       theme(plot.title = element_text(hjust = 0,face="bold", colour="Black", size = 18),
#             plot.caption = element_text(size = 14,face = "bold",hjust = 1),
#             axis.title.x = element_text(face="bold", colour="Black", size = 18),
#             axis.title.y = element_text(face="bold", colour="Black", size =18),
#             axis.text.x = element_text(face="bold", color="Black",
#                                        size=14),
#             axis.text.y = element_text(face="bold", color="Black",
#                                        size=14)) +
#     theme(panel.grid.minor = element_blank()) +
#     theme(legend.position="none")+
#     # geom_label(data = melted_varmatexcl_t,aes(label= round(value,4)),fontface = "bold",
#     #          colour = "white",size = 4) +
#     geom_text(aes(label= value),color = "black",size = 7, check_overlap = T)+
#     #scale_fill_distiller
#     scale_fill_manual(values = pal,breaks=levels(factor(melted_varmatexcl_t$value)),na.value="gray")
# 
#     plot_list[[i]] = p
#     print(plot_list[[i]])
#     dev.off()
# }
# # 
# # 
# grid.arrange(plot_list[[1]], plot_list[[2]] , plot_list[[3]], plot_list[[4]],
#              plot_list[[5]], plot_list[[6]] , plot_list[[7]], plot_list[[8]], ncol=4)



####PAPER1
# plot_list =list()
# # # Fig6 & Appendix B 
# add_name = "G:\\Shared drives\\Ehsan PhD work\\Paper_1\\plots\\Fig 7\\"
# # CONFIRMATION SEMINAR
# #add_name = "G:\\Shared drives\\Ehsan PhD work\\Confirmation\\materials\\ex3\\"
# 
# #png(paste0(add_name,"_t7",".png"), width = 1000, height = 750)
# titl1=paste("(",letters, ")",sep="")
# titl2=paste("(","a",letters, ")",sep="")
# titl=c(titl1,titl2)
# 
# for (i in unique(melted_varmatexcl_t$iter) ) {
# 
#    #png(paste0(add_name,i,".png"), width = 500, height = 500)
# 
#   p<-ggplot(subset(melted_varmatexcl_t,iter==i),aes(Period,Sequence))+
#     geom_tile(aes(fill=factor(value)),colour = "grey50") +
#     scale_y_reverse(breaks=c(1:K)) +
#     scale_x_continuous(breaks=c(1:T)) +
#     labs(title = paste0(titl[[i]]))+
#     theme(plot.title = element_text(hjust = 0,face="bold", colour="Black", size = 18),
#           axis.title.x = element_text(face="bold", colour="Black", size = 18),
#           axis.title.y = element_text(face="bold", colour="Black", size =18),
#           axis.text.x = element_text(face="bold", color="Black",
#                                      size=14),
#           axis.text.y = element_text(face="bold", color="Black",
#                                      size=14)) +
#     theme(panel.grid.minor = element_blank()) +
#     theme(legend.position="none")+
#     # geom_label(data = melted_varmatexcl_t,aes(label= round(value,4)),fontface = "bold",
#     #          colour = "white",size = 4) +
#     geom_text(aes(label=value),color = "black",size = 4, check_overlap = T)+
#     #scale_fill_distiller
#     scale_fill_manual(values = pal,breaks=levels(factor(melted_varmatexcl_t$value)),na.value="gray")
# 
#   plot_list[[i]] = p
#   print(plot_list[[i]])
#   #dev.off()
# }
# 
# png(paste0(add_name,"_B1",".png"), width = 2000, height = 500)
# grid.arrange(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],ncol=4)
# dev.off()
# 
# png(paste0(add_name,"_B2",".png"), width = 2000, height = 500)
# grid.arrange(plot_list[[5]],plot_list[[6]],plot_list[[7]],plot_list[[8]],ncol=4)
# dev.off()
# 
# png(paste0(add_name,"_B3",".png"), width = 2000, height = 500)
# grid.arrange(plot_list[[9]],plot_list[[10]],plot_list[[11]],plot_list[[12]],ncol=4)
# dev.off()
# 
# png(paste0(add_name,"_B4",".png"), width = 2000, height = 500)
# grid.arrange(plot_list[[13]],plot_list[[14]],plot_list[[15]],plot_list[[16]],ncol=4)
# dev.off()
# 
# png(paste0(add_name,"_B5",".png"), width = 2000, height = 500)
# grid.arrange(plot_list[[17]],plot_list[[18]],plot_list[[19]],plot_list[[20]],ncol=4)
# dev.off()
# 
# png(paste0(add_name,"_B6",".png"), width = 2000, height = 500)
# grid.arrange(plot_list[[21]],plot_list[[22]],plot_list[[23]],plot_list[[24]],ncol=4)
# dev.off()
# 
# png(paste0(add_name,"_B7",".png"), width = 2000, height = 500)
# grid.arrange(plot_list[[25]],plot_list[[26]],plot_list[[27]],plot_list[[28]],ncol=4)
# dev.off()
# 
# png(paste0(add_name,"_B8",".png"), width = 2000, height = 500)
# grid.arrange(plot_list[[29]],plot_list[[30]],plot_list[[31]],plot_list[[32]],ncol=4)
# dev.off()
# 
# png(paste0(add_name,"_B9",".png"), width = 2000, height = 500)
# grid.arrange(plot_list[[33]],plot_list[[34]],plot_list[[35]],plot_list[[36]],ncol=4)
# dev.off()
# 
# png(paste0(add_name,"_B10",".png"), width = 2000, height = 500)
# grid.arrange(plot_list[[37]],plot_list[[38]],plot_list[[39]],plot_list[[40]],ncol=4)
# dev.off()
# 
# png(paste0(add_name,"_B11",".png"), width = 2000, height = 500)
# grid.arrange(plot_list[[41]],plot_list[[42]],plot_list[[43]],plot_list[[44]],ncol=4)
# dev.off()
# 
#trial schemtaic
# color_palette <-colorRampPalette(c( "white", "grey"))(2)
# # 
# # add_name = "G:\\Shared drives\\Ehsan PhD work\\Ehsan\\Presentations\\ViCBiostat_July14\\schematic\\"
# # png(paste0(add_name,2,".png"), width = 500, height = 500)
# plot_list[[i+1]]=ggplot(data =subset(melted_varmatexcl_t,iter==1), aes(x=Period, y=Sequence, fill = factor(Xdvalue))) +
#   geom_tile( colour = "grey50") +
#   scale_y_reverse(breaks=c(1:K)) +
#   scale_x_continuous(breaks=c(1:T)) +
#   theme(plot.title = element_text(hjust = 0,face="bold", colour="Black", size = 18),
#         axis.title.x = element_text(face="bold", colour="Black", size = 18),
#         axis.title.y = element_text(face="bold", colour="Black", size =18),
#         axis.text.x = element_text(face="bold", color="Black",
#                                    size=14),
#         axis.text.y = element_text(face="bold", color="Black",
#                                    size=14)) +
#   theme(panel.grid.minor = element_blank()) +
#   geom_text(aes(Period, Sequence, label = Xdvalue), color = "black", size = 10) +
#   scale_fill_manual(values = color_palette) +  theme(legend.position="none")

# plot_list[[i+2]]=ggplot(data =subset(melted_varmatexcl_t,iter==24), aes(x=Period, y=Sequence, fill = factor(Xdvalue))) +
#   geom_tile( colour = "grey50") +
#   scale_y_reverse(breaks=c(1:K)) +
#   scale_x_continuous(breaks=c(1:T)) +
#   theme(plot.title = element_text(hjust = 0,face="bold", colour="Black", size = 18),
#         axis.title.x = element_text(face="bold", colour="Black", size = 18),
#         axis.title.y = element_text(face="bold", colour="Black", size =18),
#         axis.text.x = element_text(face="bold", color="Black",
#                                    size=14),
#         axis.text.y = element_text(face="bold", color="Black",
#                                    size=14)) +
#   theme(panel.grid.minor = element_blank()) +
#   geom_text(aes(Period, Sequence, label = Xdvalue), color = "black", size = 10) +
#   scale_fill_manual(values = color_palette) +  theme(legend.position="none")
# 
# grid.arrange(plot_list[[i+1]],plot_list[[i+2]],plot_list[[1]], plot_list[[24]], ncol=2)
# dev.off()


# 
#trial schemtaic
# color_palette <-colorRampPalette(c( "white", "grey"))(2)
# 
# #add_name = "G:\\Shared drives\\Ehsan PhD work\\Ehsan\\Presentations\\ViCBiostat_July14\\schematic\\"
# add_name = "G:\\Shared drives\\Ehsan PhD work\\Confirmation\\materials\\"
# png(paste0(add_name,"ts1",".png"), width = 500, height = 500)
# ggplot(data =subset(melted_varmatexcl_t,iter==1), aes(x=Period, y=Sequence, fill = factor(Xdvalue))) +
#   geom_tile( colour = "grey50") +
#   scale_y_reverse(breaks=c(1:K)) +
#   scale_x_continuous(breaks=c(1:T)) +
#   theme(plot.title = element_text(hjust = 0,face="bold", colour="Black", size = 18),
#         axis.title.x = element_text(face="bold", colour="Black", size = 18),
#         axis.title.y = element_text(face="bold", colour="Black", size =18),
#         axis.text.x = element_text(face="bold", color="Black",
#                                    size=14),
#         axis.text.y = element_text(face="bold", color="Black",
#                                    size=14)) +
#   geom_text(aes(Period, Sequence, label = Xdvalue), color = "black", size = 10) +
#   scale_fill_manual(values = color_palette) +  theme(legend.position="none")
# dev.off()
# 
# 
# # #relative variance plot
# # p <- plot_ly(res, height=500, width=800, x=~iter, y=~Rvariance, name="RVariance", type="scatter",
# #              mode="lines", hoverinfo="text", hoverlabel=list(bordercolor=NULL, font=list(size=16)),
# #              text=~paste("Iteration:", iter, "<br>RVariance:", round(Rvariance, 3)),
# #              line=list(color="#F8766D", width=4, dash="dash"))%>%
# #   layout(xaxis=list(title="Iteration", titlefont=list(size=18), showline=TRUE,
# #                     tickmode="auto", tickfont=list(size=16), nticks=6, ticks="inside",
# #                     mirror=TRUE, showgrid=FALSE),
# #          yaxis=list(title="Relative Variance", titlefont=list(size=18), tickfont=list(size=16),
# #                     mirror=TRUE, showline=TRUE),
# #          # title=list(text=paste("m=",input$m,",","T=",input$T,",","rho0=",input$rho0,",","r=",
# #          #input$r,",","effsize=",
# #          #                       input$effsize,"\n","Removing one pair at each step=",cstus), y =1),
# #          legend=list(orientation="h", xanchor="center", yanchor="bottom", x=0.5, y=-0.5, font=list(size=16)),
# #          margin=list(l=100, r=40))
# # print(p)
# # 


