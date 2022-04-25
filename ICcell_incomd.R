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
#1.048 2.24 2.69

# #(K+1-k,T+1-t)
# #m=100
# #rho0=0.05
# #r=0.95
# #lbl="Hussey and Hughes model"
# #lbl="Exponential decay model"
# #1 decay 0 HH
# type=1
# #T=5
# # 0 removing without any consideration, 1 removing only one pair at each step
# #cnum=0
# effsize=0.2

pal <- c("#FF4500", "#00FF00", "#00BFFF")
pal <- setNames(pal, c("0.95", "0.8", "0.5"))

plot = function(data, m, T,rho0) {
  #relative variance plots accross cac
plot_ly(data = res_m, x = ~iter, y = ~Rvariance,  type="scatter",linetype=~as.factor(r), colors = pal,
        mode="lines", hoverinfo="text", hoverlabel=list(bordercolor=NULL, font=list(size=14)),
        text=~paste("Iteration:", iter, "<br>RVariance:", round(Rvariance,3),
        "<br>Power:", format(round(power,4)*100,2),"%")) %>%
        layout(xaxis=list(title="Iteration", titlefont=list(size=18), showline=TRUE,
        tickmode="auto", tickfont=list(size=16), nticks=6, ticks="inside",
        mirror=TRUE, showgrid=FALSE),
        yaxis=list(title="Relative variance", titlefont=list(size=18), tickfont=list(size=16),
                    mirror=TRUE, showline=TRUE),
         title=list(text=paste("plot",l,"\n","m=",m,",","T=",T,",","rho0=",rho0,",",
                               "effsize=",effsize), y =0.98))
  
}

#p <-list ()

l=0

for (m in c(50,100)){
  for (T in c(5,7,8,10)){
    for (rho0 in c(0.01,0.05,0.1)){
      l=l+1
      
      res_m <- data.frame(iter=integer(),
                          variance=integer(),
                          power=integer(),
                          r=integer(),
                          RVariance=integer())
      
      for (r in c(0.95,0.8,0.5)){
        
        K=T-1
        Xdes <- SWdesmat(T)
        varmatall <- c()
        varmatall<- CRTVarGeneralAdj(Xdes,m,rho0,r,type)
        varmat_excl<-matrix(data=NA, nrow=nrow(Xdes), ncol=ncol(Xdes))
        
        
        #needs to be optimised in repeating function
        rep_func2 = function(Xdes,varmatall){
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
        dlist [[1]]<- rep_func2(Xdlist[[1]],varmatall)
        varmatall[1] <- varmatall
        
        #remove all low-infomration content cells and updating.
        for (i in 2:((T*K)/2)){
          mval[[i-1]] <- which(rep_func2(Xdlist[[i-1]],varmatall[i-1])==min(rep_func2(Xdlist[[i-1]],varmatall[i-1]),na.rm = TRUE), arr.ind = TRUE)
          Xdlist[[i]]=Xdlist[[i-1]]
          
          # if (cnum==0){
          # needs to be optimised in repeating loop
          for (j in 1:dim(mval[[i-1]])[1]){
            Xdlist[[i]][mval[[i-1]][[j]],mval[[i-1]][[dim(mval[[i-1]])[1]+j]]]<- NA
            Xdlist[[i]][K+1-mval[[i-1]][[j]],T+1-mval[[i-1]][[dim(mval[[i-1]])[1]+j]]]<- NA
          }
          # pattern for odd and even periods are like below
          if (sum(colSums(!is.na(Xdlist[[i]])))==4 | sum(colSums(!is.na(Xdlist[[i]])))==2) {
            if (sum(colSums(!is.na(dlist[[i-1]])))==2){
              dlist[[i-1]]<- NULL
              varmatall <- varmatall[-(i-1)] 
            }
            break
          }
          #}
          # else if (cnum==1){
          #     if (dim(mval[[i-1]])[1]==1) {
          #       Xdlist[[i]][mval[[i-1]][[1]],mval[[i-1]][[dim(mval[[i-1]])[1]+1]]]<- NA
          #       Xdlist[[i]][K+1-mval[[i-1]][[1]],T+1-mval[[i-1]][[dim(mval[[i-1]])[1]+1]]]<- NA
          #       }else if (dim(mval[[i-1]])[1]==2 & (mval[[i-1]][[2]]==K+1-mval[[i-1]][[1]]) & (mval[[i-1]][[dim(mval[[i-1]])[1]+2]]==T+1-mval[[i-1]][[dim(mval[[i-1]])[1]+1]])){
          #         Xdlist[[i]][mval[[i-1]][[1]],mval[[i-1]][[dim(mval[[i-1]])[1]+1]]]<- NA
          #         Xdlist[[i]][K+1-mval[[i-1]][[1]],T+1-mval[[i-1]][[dim(mval[[i-1]])[1]+1]]]<- NA
          #       }
          #    else if (dim(mval[[i-1]])[1]>=2 & (mval[[i-1]][[2]]!=K+1-mval[[i-1]][[1]] | mval[[i-1]][[dim(mval[[i-1]])[1]+2]]!=T+1-mval[[i-1]][[dim(mval[[i-1]])[1]+1]])) {
          #           mval[[i-1]] <- mval[[i-1]][order(mval[[i-1]][,1],mval[[i-1]][,2]),]
          #           Xdlist[[i]][mval[[i-1]][[1]],mval[[i-1]][[dim(mval[[i-1]])[1]+1]]]<- NA
          #           Xdlist[[i]][K+1-mval[[i-1]][[1]],T+1-mval[[i-1]][[dim(mval[[i-1]])[1]+1]]]<- NA
          #         }
          #       #pattern for odd and even periods are like below
          #       if ((T %% 2 != 0 & sum(colSums(!is.na(Xdlist[[i]])))==4) | (T %% 2 == 0 & sum(colSums(!is.na(Xdlist[[i]])))==2)) {
          #         break
          #       }
          #   }
          #varmatall[i] <- round(CRTVarGeneralAdj(Xdlist[[i]],m,rho0,r,type),10)
          varmatall[i] <- CRTVarGeneralAdj(Xdlist[[i]],m,rho0,r,type)
          dlist[[i]] = rep_func2(Xdlist[[i]],varmatall[i])
        }
        
        
        melted_varmatexcl_t <- data.frame( Var1=integer(),
                                           Var2=integer(),
                                           value=integer(),
                                           iter=integer())
        
        for (i in 1:length(dlist)){
          varmat_excl<-round(dlist[[i]], 4)
          melted_varmatexcl <- melt(varmat_excl)
          melted_varmatexcl$iter<- i
          melted_varmatexcl_t <- rbind(melted_varmatexcl, melted_varmatexcl_t)
        }
        
        #color_palette <-colorRampPalette(c( "yellow", "red"))(length(table(varmat_excl)))
        
        pal <- colorRampPalette(brewer.pal(8, "YlOrRd"))(length(unique(melted_varmatexcl_t$value))-1)
        
        names(melted_varmatexcl_t)[names(melted_varmatexcl_t)=="Var1"] <- "Sequence"
        names(melted_varmatexcl_t)[names(melted_varmatexcl_t)=="Var2"] <- "Period"
        
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
          powdf <- data.frame(iter, df$varmatall,powvals)
          colnames(powdf) <- c("iter","variance","power")
          return(powdf)
        }
        
        res <- powdf(df,effsize)
        res$r <- r
        
        res <- cbind(res,res$variance/res$variance[1])
        colnames(res) <- c("iter","variance","power","r","Rvariance")
        
        res_m <- rbind(res_m,res)
        
        #end of loop for cac
      }    
      #save plots
   
      p[[l]]= plot(res_m,m,T,rho0)
      #print(i)
      #htmlwidgets::saveWidget(p,paste0("G:/Shared drives/Ehsan PhD work/Outputs/RV/RV_","m",m,"_","T",T,"_","icc"," ",rho0,".html"))
      #end of loop for desired settings
    }
  }
}
#p

# combineWidgets(list = p, title = "")
# htmlwidgets::saveWidget(p, "test.html")
# saveWidget(p, "p1.html", selfcontained = F, libdir = "lib")

#paste0("G:/Shared drives/Ehsan PhD work/Outputs/RV/RV_",m,"_",T,"_",rho0,".pdf")
#heatmap plot 

# melted_varmatexcl_t <- merge(res, melted_varmatexcl_t, by = "iter", all = TRUE)
# 
#   p<-ggplot(melted_varmatexcl_t,aes(x=Period, y=Sequence,fill=factor(value),frame=iter))+
#     geom_tile(colour = "grey50") +
#     scale_y_reverse(breaks=c(1:K)) +
#     scale_x_continuous(breaks=c(1:T)) +
#     theme(panel.grid.minor = element_blank()) +
#     geom_text(x=0.9,y=-0.4,hjust=0,aes(label=paste0("Power:",format(round(power,4)*100,2),"%")),
#               size=4,fontface="")+ 
#     theme(legend.position="none")+
#     geom_label(data = melted_varmatexcl_t,aes(label= round(value,4),fontface = "bold"),
#     colour = "white",size = 4) +
#     scale_fill_manual(values = pal, breaks=levels(melted_varmatexcl_t$value)[seq(90, 150, by=5)],
#                       na.value="gray")
#     
#     p1<- ggplotly(p) %>% 
#       animation_opts(frame = 500, transition = 0,redraw = TRUE) %>% 
#       animation_slider(currentvalue = list(prefix = "Iter: ", font = list(color="orange")))

# htmlwidgets::saveWidget(as_widget(p1), "index.html")
#relative variance plot
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



