df<-read.table('c:/Users/zhoud2/Dropbox/POLYGENIC_LARGE_EFFECT/simulation/simulation.txt',header = T,stringsAsFactors = F)


#rm NA all cases are carriers
df<-df[!is.na(df[,10]),]

summary<-df[1,1:6]

i=1
for (h2_p in c(0.1,0.2,0.3,0.4,0.5)){
  for (N_causal_SNPs in c(100)){
    for (vr in c(0.01,0.02,0.03,0.05,0.1)){
      for (freq_r in c(0.005,0.01,0.02)){
        for (k in c(0.02,0.05,0.1)){
          for (N_samples in c(1000,2000,3000,5000,10000)){
            print(i)
            h2_p_pos<-which(df$h2_p==h2_p)
            N_causal_SNPs_pos<-which(df$N_causal_SNPs==N_causal_SNPs)
            vr_pos<-which(df$vr==vr)
            freq_r_pos<-which(df$freq_r==freq_r)
            k_pos<-which(df$k==k)
            N_samples_pos<-which(df$N_samples==N_samples)
            
            tmp_pos<-Reduce(intersect,list(h2_p_pos,N_causal_SNPs_pos,vr_pos,freq_r_pos,k_pos,N_samples_pos))
            
            df_tmp<-df[tmp_pos,]
            summary[i,]<-c(df_tmp[1,1:6])
            
            #power
            summary[i,7]<-length(which(df_tmp$p_wilcox<0.05))/nrow(df_tmp) 
            
            #penetrance
            summary[i,8]<-mean(df_tmp$penetrance)
            
            #mean_polygenetic_in_cases
            summary[i,9]<-mean(df_tmp$mean_prs_cases)
            
            #h2_poly
            summary[i,10]<-mean(df_tmp$h2_poly_est)
            
            i=i+1
            
          }
        }
      }
    }
  }
}

colnames(summary)[7:10]<-c('power','penetrance','mean_prs_cases','h2_poly_est')


write.csv(summary,'c:/Users/zhoud2/Dropbox/POLYGENIC_LARGE_EFFECT/simulation/summary.csv',quote = F,row.names = F)







#---plot---power---
library(gcookbook)
library(ggplot2)
library(RColorBrewer)
library(grid)


summary<-read.csv('c:/Users/zhoud2/Dropbox/POLYGENIC_LARGE_EFFECT/simulation/summary.csv',head=T,stringsAsFactors = F)
colnames(summary)<-c('h2_polygenic','N_causal_SNPs','h2_rare','freq_r','prevalence','N_samples','power','penetrance','mean_prs_cases','h2_poly')


#-------------------------------------------------
#h2_polygenic
#h2_rare
#freq_r=0.003       
#prevalence=0.05
#N_samples


for (freq_r in c(0.005,0.01,0.02)){
  for(prevalance in c(0.02,0.05,0.10)){
    
    par(mfcol=c(1,1))
    
    df_plot<-summary[summary$freq_r==freq_r & summary$prevalence==prevalance,]
    
    png(paste0('c:/Users/zhoud2/Dropbox/POLYGENIC_LARGE_EFFECT/simulation/plots/freqr',freq_r,'_k',prevalance,'.png'),width = 1800,height = 1000,res = 150)
    
    p1<-ggplot(df_plot[df_plot$h2_rare==0.01,], aes(x=factor(h2_polygenic), y=power, colour=factor(N_samples),group=N_samples)) + geom_line(size=1.5,position=position_dodge(0.2)) + geom_point(position=position_dodge(0.2), size=3,shape=17) + ggtitle(paste0('k=',prevalance,' freq_r=',freq_r,', h2_rare=0.01')) + ylim(0, 1)
    
    p2<-ggplot(df_plot[df_plot$h2_rare==0.02,], aes(x=factor(h2_polygenic), y=power, colour=factor(N_samples),group=N_samples)) + geom_line(size=1.5,position=position_dodge(0.2)) + geom_point(position=position_dodge(0.2), size=3,shape=17) + ggtitle(paste0('k=',prevalance,' freq_r=',freq_r,', h2_rare=0.02')) + ylim(0, 1)
    
    p3<-ggplot(df_plot[df_plot$h2_rare==0.03,], aes(x=factor(h2_polygenic), y=power, colour=factor(N_samples),group=N_samples)) + geom_line(size=1.5,position=position_dodge(0.2)) + geom_point(position=position_dodge(0.2), size=3,shape=17) + ggtitle(paste0('k=',prevalance,' freq_r=',freq_r,', h2_rare=0.03')) + ylim(0, 1)
    
    p4<-ggplot(df_plot[df_plot$h2_rare==0.05,], aes(x=factor(h2_polygenic), y=power, colour=factor(N_samples),group=N_samples)) + geom_line(size=1.5,position=position_dodge(0.2)) + geom_point(position=position_dodge(0.2), size=3,shape=17) + ggtitle(paste0('k=',prevalance,' freq_r=',freq_r,', h2_rare=0.05')) + ylim(0, 1)
    
    p5<-ggplot(df_plot[df_plot$h2_rare==0.10,], aes(x=factor(h2_polygenic), y=power, colour=factor(N_samples),group=N_samples)) + geom_line(size=1.5,position=position_dodge(0.2)) + geom_point(position=position_dodge(0.2), size=3,shape=17) + ggtitle(paste0('k=',prevalance,' freq_r=',freq_r,', h2_rare=0.10')) + ylim(0, 1)
    
    multiplot(p1, p2, p3, p4, p5, cols=3)
    dev.off()
  }
}







multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols),byrow = T)
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}






#plot penetrance
#install.packages("vioplot")
library("vioplot")
library(RColorBrewer)
display.brewer.all£¨type = "seq"£©

#barplot(rep(1,6),col=brewer.pal(9, "YlOrRd")[3:7])

summary<-read.csv('c:/Users/Dan/Dropbox/POLYGENIC_LARGE_EFFECT/simulation/summary.csv',head=T,stringsAsFactors = F)

# png('c:/Users/Dan/Dropbox/POLYGENIC_LARGE_EFFECT/simulation/plots/h2_polygenic_penetrance.png',height = 600,width = 800,res=150)
# vioplot(summary[summary$h2_p==0.1,8],summary[summary$h2_p==0.2,8],summary[summary$h2_p==0.3,8],summary[summary$h2_p==0.4,8],summary[summary$h2_p==0.5,8],names=c('0.1','0.2','0.3','0.4','0.5'),xlab='h2_polygenic',ylab='penetrance',main='',col=brewer.pal(9, "YlOrRd")[3:7])
# dev.off()

png('c:/Users/Dan/Dropbox/POLYGENIC_LARGE_EFFECT/simulation/plots/h2_polygenic_penetrance_mean.png',height = 800,width = 800,res=150)
plot(summary$mean_prs_cases,summary$penetrance)
abline(lm(summary$penetrance ~ summary$mean_prs_cases),lty=2,col='red')
dev.off()


#plot h2_poly est true

summary<-read.table('c:/Users/Dan/Dropbox/POLYGENIC_LARGE_EFFECT/simulation/simulation.txt',head=T,stringsAsFactors = F)

png('c:/Users/Dan/Dropbox/POLYGENIC_LARGE_EFFECT/simulation/plots/h2_poly_true_est.png',height = 800,width = 800,res=150)
vioplot(summary[summary$h2_p==0.1,13],summary[summary$h2_p==0.2,13],summary[summary$h2_p==0.3,13],summary[summary$h2_p==0.4,13],summary[summary$h2_p==0.5,13],names=c('0.1','0.2','0.3','0.4','0.5'),xlab='h2_poly_true',ylab='h2_poly_estimated',main='',col=brewer.pal(9, "YlOrRd")[3:7])
dev.off()























