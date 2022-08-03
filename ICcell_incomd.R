install.packages("rsconnect")
library("rsconnect")


#setwd("~/Google Drive/Shared drives/Ehsan PhD work/Codes/")
#setwd("G:\\Shared drives\\Ehsan PhD work\\Codes\\Git\\Iterative-Removals-SW")
source("CRTVarAdj_func.R", local=TRUE)
source("ICcell_appfunc.R", local=TRUE)
# Functions for generating design matrices
SWdesmat <- function(T) {
  Xsw <- matrix(data=0, ncol = T, nrow = (T-1))
  for(i in 1:(T-1)) {
    Xsw[i,(i+1):T] <- 1
  }
  return(Xsw)
}
#######

###############make all of them the same shape! last design (information content matrix) should be looking empty for all
###############e.g. m=20 rho0=0.05 r=0.95 type=1 T=7
###############Add this adjustment to codes (think about it)
################NO NEED TO ADJUST AS THE FIGURE IS SHOWN NA FOR EVERY SITUATION

# Xdlist
# dlist
# (K+1-k,T+1-t)
m=90
rho0=0.15
r=0.95
#lbl="Hussey and Hughes model"
#lbl="Exponential decay model"
#1 decay 0 HH
type=1
T=5
# 0 removing without any consideration, 1 removing only one pair at each step
cnum=1
effsize=0.36
#
# pal <- c("#FF4500", "#00FF00", "#00BFFF")
# pal <- setNames(pal, c("0.95", "0.8","0.5"))
# # 
# m <- res_m[which(res_m$na_percnt)==median(res_m$na_percnt), ]
#linetype=~as.factor(r)

cac_type <- unique(res_m$r)
linetypes <- setNames(c('solid', 'dot', 'dash', 'longdash', 'dashdot'), cac_type)

plot = function(data, m, T,rho0) {
  #relative variance plots accross cac
plot_ly(data = res_m, x = ~na_percnt, y = ~Effloss, type="scatter",linetype=~as.factor(r),linetypes=linetypes,
        mode="lines", hoverinfo="text", hoverlabel=list(bordercolor=NULL, font=list(size=14)),
        line=list(dash = "dash",color="#FF6666"),
        text=~paste("RPercent:", round(na_percnt,2),"%", "<br>Effloss:", round(Effloss,2),"%",
        "<br>Power:", format(round(power,4)*100,2),"%")) %>%
        layout(xaxis=list(title="Removal percentage of cluster-period cells (%)", titlefont=list(size=18), showline=TRUE,
        tickmode="auto", tickfont=list(size=16), nticks=11, range=c(0,100),ticks="inside",
        mirror=TRUE, showgrid=FALSE),
        yaxis=list(title="Efficiency loss (%)", titlefont=list(size=18), tickfont=list(size=16),
                    mirror=TRUE, showline=TRUE,range=c(0,100)),legend = list(title=list(text='<b> cac </b>')))%>% 
        # title=list(text=paste0("\n","m=",m,","," ","T=",T,","," ","icc","=",rho0,","," ","effsize=",effsize),
        #                      y =0.95,font=list(size = 15))) %>% 
         add_segments(y = 0, yend = 50, x = 50, xend = 50,line = list(dash = "dash",color="#00BFFF"), showlegend=FALSE,hoverinfo='skip')
         # title=list(text=paste(if (type == 1) {paste0("plot"," ",l," ","(","Discrete time decay",")")
         # } else if (type == 0){paste0("plot"," ",l," ","(","Exchangeable",")")},
         # paste0("\n","m=",m,","," ","T=",T,","," ","icc","=",rho0,","," ","effsize=",effsize))
         #                      ,y =0.95,font=list(size = 15)))
  
}


p=plot(res_m,m,T,rho0)
save_image(p, file="Eff_loss_plot.png",width = 1150, height = 550)

#single line
plot = function(data, m, T,rho0) {
  #relative variance plots accross cac
  plot_ly(data = res_m, x = ~na_percnt, y = ~Effloss, type="scatter",name = "Effloss", mode="lines",
          hoverinfo="text", hoverlabel=list(bordercolor=NULL, font=list(size=14)),line=list(dash ="dash",color="#FF6666", width = 6),
          text=~paste("RPercent:", round(na_percnt,2),"%", "<br>Effloss:", round(Effloss,2),"%",
                      "<br>Power:", round(power,2),"%")) %>%
  layout(xaxis=list(title="Removal percentage of cluster-period cells (%)", titlefont=list(size=18), showline=TRUE,
                                tickmode="auto", tickfont=list(size=16), nticks=11, range=c(0,100),ticks="inside",
                              mirror=TRUE, showgrid=FALSE),
        yaxis=list(title="Efficiency loss (%)", titlefont=list(size=18), tickfont=list(size=16),
                   mirror=TRUE, showline=TRUE,range=c(0,100)),legend = list(title=list(text='<b> Type </b>')))# %>% 
   # add_segments(y = 0, yend = 50, x = 50, xend = 50,line = list(dash = "dash",color="#00BFFF"), showlegend=FALSE,hoverinfo='skip')
}

p=plot(res_m,m,T,rho0)
save_image(p, file="efflossPos.png",width = 1150, height = 550)

ay <- list(
  tickfont = list(color = "#00CC00"),
  overlaying = "y",
  side = "right",
  title = "<b>Power (%)</b>")

###add power
plot = function(data, m, T,rho0) {
  #relative variance plots accross cac
  plot_ly(data = res_m, x = ~na_percnt, y = ~Effloss, type="scatter",name = "Effloss", mode="lines",
         hoverinfo="text", hoverlabel=list(bordercolor=NULL, font=list(size=14)),line=list(dash ="dash",color="#FF6666"),
          text=~paste("RPercent:", round(na_percnt,2),"%", "<br>Effloss:", round(Effloss,2),"%",
                      "<br>Power:", round(power,2),"%")) #%>%
  # add_trace(data = res_m, x = ~na_percnt, y = ~power, type="scatter",name = "power",yaxis = "y2",mode="lines",
  #           line=list(dash ="dot",color="#00CC00")) %>% 
  # layout(yaxis2 = ay,xaxis=list(title="Removal percentage of cluster-period cells (%)", titlefont=list(size=18), showline=TRUE,
  #                               tickmode="auto", tickfont=list(size=16), nticks=11, range=c(0,100),ticks="inside",
  #                               mirror=TRUE, showgrid=FALSE),
  #        yaxis=list(title="Efficiency loss (%)", titlefont=list(size=18), tickfont=list(size=16),
  #                   mirror=TRUE, showline=TRUE,range=c(0,100)),legend = list(title=list(text='<b> Type </b>')))%>% 
  #   # title=list(text=paste0("\n","m=",m,","," ","T=",T,","," ","icc","=",rho0,","," ","effsize=",effsize),
  #   #                      y =0.95,font=list(size = 15))) %>% 
  #   add_segments(y = 0, yend = 50, x = 50, xend = 50,line = list(dash = "dash",color="#00BFFF"), showlegend=FALSE,hoverinfo='skip')  
    # title=list(text=paste(if (type == 1) {paste0("plot"," ",l," ","(","Discrete time decay",")")
    # } else if (type == 0){paste0("plot"," ",l," ","(","Exchangeable",")")},
    # paste0("\n","m=",m,","," ","T=",T,","," ","icc","=",rho0,","," ","effsize=",effsize))
    #                      ,y =0.95,font=list(size = 15)))
}

p=plot(res_m,m,T,rho0)
save_image(p, file="efflossPow_paper1.png",width = 1150, height = 550)

###Paper
####two separate icc (0.15, 0.14)
plot = function(data, m, T,rho0) {
plot_ly(data = res_m, x = ~na_percnt, y = ~Effloss, type="scatter",name = "Discrete-time Decay", mode="lines",
        hoverinfo="text", hoverlabel=list(bordercolor=NULL, font=list(size=14)),line=list(dash ="dot",color="#FF6666"),
        #,line=list(dash = "dash",color="#FF4500")
        text=~paste("RPercent:", round(na_percnt,2),"%", "<br>Effloss:", round(Effloss,2),"%",
                    "<br>Power:", round(power,2),"%")) %>%
  layout(xaxis=list(title="Removal percentage of cluster-period cells (%)", titlefont=list(size=18), showline=TRUE,
                    tickmode="auto", tickfont=list(size=16), nticks=11, range=c(0,100),ticks="inside",
                    mirror=TRUE, showgrid=FALSE),
         yaxis=list(title="Efficiency loss (%)", titlefont=list(size=18), tickfont=list(size=16),
                    mirror=TRUE, showline=TRUE,range=c(0,100)),legend = list(title=list(text='<b>Type</b>')))%>% 
  # title=list(text=paste0("\n","m=",m,","," ","T=",T,","," ","icc","=",rho0,","," ","effsize=",effsize),
  #                      y =0.95,font=list(size = 15))) %>% 
  add_segments(y = 0, yend = 50, x = 50, xend = 50,line = list(dash = "dash",color="#00BFFF"), showlegend=FALSE,hoverinfo='skip')%>%
# title=list(text=paste(if (type == 1) {paste0("plot"," ",l," ","(","Discrete time decay",")")
# } else if (type == 0){paste0("plot"," ",l," ","(","Exchangeable",")")},
# paste0("\n","m=",m,","," ","T=",T,","," ","icc","=",rho0,","," ","effsize=",effsize))
#                      ,y =0.95,font=list(size = 15)))
  add_trace(data = res_m1, x = ~na_percnt, y = ~Effloss, type="scatter",name = "Exchangeable",mode="lines",
            line=list(dash ="solid",color="#B266FF"))

}


p=plot(data,m,T,rho0)
save_image(p, file="effloss_paper1_fig5.png",width = 1150, height = 550)

l=0
#cac = c(0.95,0.8,0.5,0.2)
cac = 0.95
#tp = c(0,1)


#p <-list ()
# c(20,50,100)
# c(0.01,0.05,0.1)
# c(5,7,8,10)
for (m in 90){
  for (T in 5){
    for (rho0 in 0.15){
      
      l=l+1

res_m <- data.frame(iter=integer(),
                      variance=integer(),
                      power=integer(),
                      r=integer(),
                      RVariance=integer(),
                      Effloss=integer())


for (type in c(0,1)){

  for (r in  if (type ==0){r = 1}  else {r = cac})
    {

#r = cac
        
        K=T-1
        Xdes <- SWdesmat(T)
        varmatall <- c()
        varmatall<- CRTVarGeneralAdj(Xdes,m,rho0,r,type)
        varmat_excl<-matrix(data=NA, nrow=nrow(Xdes), ncol=ncol(Xdes))
        
        
        #needs to be optimised in repeating function
        IC_func2 = function(Xdes,varmatall){
          for(i in 1:nrow(Xdes)){
            for (j in 1:ncol(Xdes)){
              if(is.na(Xdes[i,j])==TRUE | is.na(Xdes[K-i+1,T-j+1])==TRUE){
                varmat_excl[i,j] <- NA
                varmat_excl[K-i+1,T-j+1] <- NA
              }
              else if(is.na(Xdes[i,j])==FALSE & is.na(Xdes[K-i+1,T-j+1])==FALSE) {
                Xdesij <- Xdes
                Xdesij[i,j] <- NA
                Xdesij[K-i+1,T-j+1] <- NA
                varmat_excl[i,j] <- CRTVarGeneralAdj(Xdesij,m,rho0,r,type)/varmatall
                varmat_excl[K-i+1,T-j+1] <- varmat_excl[i,j]
              }
            }
          }
          return(varmat_excl)
        }
        mval <- list()
        Xdlist <- list()
        dlist <- list()
        Xdlist[[1]] <- Xdes
        dlist [[1]]<- IC_func2(Xdlist[[1]],varmatall)
        varmatall[1] <- varmatall
        
        # output <- tryCatch({expr},
        #                    warning = function(w) {warning-handler-code},
        #                    error = function(e) {error-handler-code})
        # 
        
        #remove all low-infomration content cells and updating.
        for (i in 2:((T*K)/2)){
        #lapply(2:((T*K)/2))
          mval[[i-1]] <- tryCatch(which(IC_func2(Xdlist[[i-1]],varmatall[i-1])==min(IC_func2(Xdlist[[i-1]],varmatall[i-1]),na.rm = TRUE), arr.ind = TRUE), warning=function(w) NA)
        
          if (is.na(mval[[i-1]][1])) {
            dlist[[i-1]]<- NULL
            break
          }
          
          Xdlist[[i]]=Xdlist[[i-1]]
          
          
            if (cnum==0){
          # loop is repeating, it might be better to adjust code
          tryCatch(for (j in 1:dim(mval[[i-1]])[1]){
            Xdlist[[i]][mval[[i-1]][[j]],mval[[i-1]][[dim(mval[[i-1]])[1]+j]]]<- NA
            Xdlist[[i]][K+1-mval[[i-1]][[j]],T+1-mval[[i-1]][[dim(mval[[i-1]])[1]+j]]]<- NA
          }, error=function(e) NA)
            }
          #remove the smallest cluster and period, and removing the corresponding pair. 
          else if (cnum==1){
              # if (dim(mval[[i-1]])[1]==1) {
              #   Xdlist[[i]][mval[[i-1]][[1]],mval[[i-1]][[dim(mval[[i-1]])[1]+1]]]<- NA
              #   Xdlist[[i]][K+1-mval[[i-1]][[1]],T+1-mval[[i-1]][[dim(mval[[i-1]])[1]+1]]]<- NA
              #   }else if (dim(mval[[i-1]])[1]==2 & (mval[[i-1]][[2]]==K+1-mval[[i-1]][[1]]) & (mval[[i-1]][[dim(mval[[i-1]])[1]+2]]==T+1-mval[[i-1]][[dim(mval[[i-1]])[1]+1]])){
              #     Xdlist[[i]][mval[[i-1]][[1]],mval[[i-1]][[dim(mval[[i-1]])[1]+1]]]<- NA
              #     Xdlist[[i]][K+1-mval[[i-1]][[1]],T+1-mval[[i-1]][[dim(mval[[i-1]])[1]+1]]]<- NA
              #   }
            # else if (dim(mval[[i-1]])[1]>=2 & (mval[[i-1]][[2]]!=K+1-mval[[i-1]][[1]] | mval[[i-1]][[dim(mval[[i-1]])[1]+2]]!=T+1-mval[[i-1]][[dim(mval[[i-1]])[1]+1]])) {
                    mval[[i-1]] <- mval[[i-1]][order(mval[[i-1]][,1],mval[[i-1]][,2]),]
                    Xdlist[[i]][mval[[i-1]][[1]],mval[[i-1]][[dim(mval[[i-1]])[1]+1]]]<- NA
                    Xdlist[[i]][K+1-mval[[i-1]][[1]],T+1-mval[[i-1]][[dim(mval[[i-1]])[1]+1]]]<- NA
                  }
           # }
          #varmatall[i] <- round(CRTVarGeneralAdj(Xdlist[[i]],m,rho0,r,type),10)
          varmatall[i] <-tryCatch(CRTVarGeneralAdj(Xdlist[[i]],m,rho0,r,type),error=function(e) NA)
          if (is.na(varmatall[i])) {
            Xdlist[[i]]<- NULL
            break
          }
          dlist[[i]] = IC_func2(Xdlist[[i]],varmatall[i])
        }
        
        
         # melted_varmatexcl_t <- data.frame( Var1=integer(),
         #                                   Var2=integer(),
         #                                   value=integer(),
         #                                   iter=integer())

        # for (i in 1:length(dlist)){
        #   varmat_excl<-round(dlist[[i]], 4)
        #   melted_varmatexcl <- melt(varmat_excl)
        #   melted_desmatexcl <- melt(Xdlist)
        #   melted_varmatexcl$iter<- i
        #   melted_varmatexcl_t <- rbind(melted_varmatexcl, melted_varmatexcl_t)
        # }
        #  
         
         melted_varmatexcl<- melt(dlist)
         #varmat_excl$value<-round(varmat_excl$value, 4)
         melted_desmatexcl<- melt(Xdlist)
         names(melted_desmatexcl)[names(melted_desmatexcl)=="value"] <- "Xdvalue"
         melted_varmatexcl_t<- merge(melted_varmatexcl, melted_desmatexcl,all.y=TRUE,by = c('Var1','Var2','L1'))
         melted_varmatexcl_t$value<-round(melted_varmatexcl_t$value, 4)

        #color_palette <-colorRampPalette(c( "yellow", "red"))(length(table(varmat_excl)))

        pal <- colorRampPalette(brewer.pal(8, "YlOrRd"))(length(unique(melted_varmatexcl_t$value))-1)

        names(melted_varmatexcl_t)[names(melted_varmatexcl_t)=="Var1"] <- "Sequence"
        names(melted_varmatexcl_t)[names(melted_varmatexcl_t)=="Var2"] <- "Period"
        names(melted_varmatexcl_t)[names(melted_varmatexcl_t)=="L1"] <- "iter"

        T <- ncol(Xdes)
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
        
        res_m <- rbind(res_m,res_2)

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

    melted_varmatexcl_t <- merge(res, melted_varmatexcl_t, by = "iter", all = TRUE)
    
    
    p<-ggplot(melted_varmatexcl_t,aes(Period,Sequence,frame=iter))+
      geom_tile(aes(fill=factor(value)),colour = "grey50") +
      scale_y_reverse(breaks=c(1:K)) +
      scale_x_continuous(breaks=c(1:T)) +
      theme(panel.grid.minor = element_blank()) +
      theme(legend.position="none")+
      # geom_label(data = melted_varmatexcl_t,aes(label= round(value,4)),fontface = "bold",
      #          colour = "white",size = 4) +
      geom_text(aes(label= Xdvalue),color = "black",size = 7, check_overlap = T)+
      #scale_fill_distiller
      scale_fill_manual(values = pal,breaks=levels(melted_varmatexcl_t$value)[seq(90, 150, by=5)],na.value="gray")
      
    p1<-ggplotly(p) %>% 
      animation_opts(frame = 500,transition =0,redraw = TRUE) %>%  
      animation_slider(currentvalue = list(prefix = "Iter: ", font = list(color="darkblue"))) %>%
      partial_bundle(local = FALSE)%>% 
      animation_button(visible = T)%>%
    animation_slider(visible =T)
    #error seems ok: this is for geom_text
p1

htmlwidgets::saveWidget(as_widget(p1), "index.html")
 

 
####PRESNTATION

plot_list =list()
add_name = "G:\\Shared drives\\Ehsan PhD work\\Publications\\43rd ISCB\\Figures_2\\"
 
for (i in unique(melted_varmatexcl_t$iter)) {
  
    png(paste0(add_name,i,".png"), width = 500, height = 500) 
  
    p<-ggplot(subset(melted_varmatexcl_t,iter==i),aes(Period,Sequence))+
    geom_tile(aes(fill=factor(value)),colour = "grey50") +
    scale_y_reverse(breaks=c(1:K)) +
    scale_x_continuous(breaks=c(1:T)) +
    theme(panel.grid.minor = element_blank()) +
    theme(legend.position="none")+
    # geom_label(data = melted_varmatexcl_t,aes(label= round(value,4)),fontface = "bold",
    #          colour = "white",size = 4) +
    geom_text(aes(label= Xdvalue),color = "black",size = 8, check_overlap = T)+
    #scale_fill_distiller
    scale_fill_manual(values = pal,breaks=levels(factor(melted_varmatexcl_t$value)),na.value="gray") 
  
    plot_list[[i]] = p
    print(plot_list[[i]])
    dev.off()
}

#trial schemtaic
color_palette <-colorRampPalette(c( "white", "grey"))(2)

add_name = "G:\\Shared drives\\Ehsan PhD work\\Ehsan\\Presentations\\ViCBiostat_July14\\schematic\\"
png(paste0(add_name,2,".png"), width = 500, height = 500) 
ggplot(data =subset(melted_varmatexcl_t,iter==1), aes(x=Period, y=Sequence, fill = factor(Xdvalue))) + 
  geom_tile( colour = "grey50") +
  scale_y_reverse(breaks=c(1:K)) +
  scale_x_continuous(breaks=c(1:T)) +
  theme(panel.grid.minor = element_blank()) +           
  geom_text(aes(Period, Sequence, label = Xdvalue), color = "black", size = 10) +
  scale_fill_manual(values = color_palette) +  theme(legend.position="none")
dev.off()


# #relative variance plot
# p <- plot_ly(res, height=500, width=800, x=~iter, y=~Rvariance, name="RVariance", type="scatter",
#              mode="lines", hoverinfo="text", hoverlabel=list(bordercolor=NULL, font=list(size=16)),
#              text=~paste("Iteration:", iter, "<br>RVariance:", round(Rvariance, 3)),
#              line=list(color="#F8766D", width=4, dash="dash"))%>%
#   layout(xaxis=list(title="Iteration", titlefont=list(size=18), showline=TRUE,
#                     tickmode="auto", tickfont=list(size=16), nticks=6, ticks="inside",
#                     mirror=TRUE, showgrid=FALSE),
#          yaxis=list(title="Relative Variance", titlefont=list(size=18), tickfont=list(size=16),
#                     mirror=TRUE, showline=TRUE),
#          # title=list(text=paste("m=",input$m,",","T=",input$T,",","rho0=",input$rho0,",","r=",
#          #input$r,",","effsize=",
#          #                       input$effsize,"\n","Removing one pair at each step=",cstus), y =1),
#          legend=list(orientation="h", xanchor="center", yanchor="bottom", x=0.5, y=-0.5, font=list(size=16)),
#          margin=list(l=100, r=40))
# print(p)
# 


