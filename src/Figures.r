#--------------------utility------------------------------

library(gcookbook)
library(ggplot2)
library(RColorBrewer)
library(grid)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL,by_row=T) {
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
                     ncol = cols, nrow = ceiling(numPlots/cols),byrow = by_row)
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


#----inf probit---

i=1

for(sample in c('case')){ #,'control'
  
  for(disease_model in c('probit')){  #,'logit'
    
    single_var<-list()
    
    for(genetic_model in c('a')){ #,'b','c'
      
      if(genetic_model=='a'){
        genetic_model_label='Polygenic'
      }else if(genetic_model=='b'){
        genetic_model_label='Negative selection'
      }else{
        genetic_model_label='LD-adjusted kinship'
      }
      
      summary<-read.table(paste0('~/../Dropbox/POLYGENIC_LARGE_EFFECT/simulation/results/result_',disease_model,'_',genetic_model,'.txt'),head=T,stringsAsFactors = F)
      summary<-summary[,c(colnames(summary)[1:7],paste0(sample,'_inverse_power'))]
      colnames(summary)<-c('h2_poly','N_causal_SNPs','h2_rare','freq_rare','N_samples','prevalence','pi0','inverse_utility')
      
      #--h2_poly--
      
      #h2_poly=c(0.1,0.2,0.3,0.4,0.5)[3]
      h2_rare=c(0.01,0.02,0.03,0.05,0.10)[3]
      freq_rare=c(1e-4,0.001,0.005,0.01,0.05)[3]
      prevalence=c(0.001,0.005,0.01,0.02,0.05)[3]
      #N_samples=c(500,1000,2000,3000,5000,10000)
      
      df_plot<-summary[summary$h2_rare==h2_rare & summary$freq_rare==freq_rare & summary$prevalence==prevalence,]
      
      single_var[[paste0(sample,'_',genetic_model,'_h2_poly')]]<-ggplot(df_plot, aes(x=factor(h2_poly), y=inverse_utility, colour=factor(N_samples),group=N_samples)) + 
        geom_line(size=1.5) +#,position=position_dodge(0.2)
        geom_point(size=3,shape=17) + #position=position_dodge(0.2), 
        ggtitle(paste0('freq_LEV(f)=',freq_rare,', h2_LEV=',h2_rare)) + 
        ylim(0,1) +
        xlab('h2_cSEV') + 
        ylab('utility') +
        labs(col = "sample size",tag = letters[i])+
        theme(legend.margin=margin(0,0,0,0),
              legend.box.margin=margin(-10,0,-10,-10),
              plot.title = element_text(size = 10),
              panel.border = element_blank(),
              panel.background = element_blank(),
              panel.grid =element_blank(),
              axis.line = element_line(colour = "black"))
      i=i+1
      
      
      #--h2_rare--
      h2_poly=c(0.1,0.2,0.3,0.4,0.5)[3]
      #h2_rare=c(0.01,0.02,0.03,0.05,0.10)[3]
      freq_rare=c(1e-4,0.001,0.005,0.01,0.05)[3]
      prevalence=c(0.001,0.005,0.01,0.02,0.05)[3]
      #N_samples=c(500,1000,2000,3000,5000,10000)
      
      df_plot<-summary[summary$h2_poly==h2_poly & summary$freq_rare==freq_rare & summary$prevalence==prevalence,]
      
      single_var[[paste0(sample,'_',genetic_model,'_h2_rare')]]<-ggplot(df_plot, aes(x=factor(h2_rare), y=inverse_utility, colour=factor(N_samples),group=N_samples)) + 
        geom_line(size=1.5) + 
        geom_point(size=3,shape=17) + 
        ggtitle(paste0('freq_LEV(f)=',freq_rare,', h2_cSEV=',h2_poly)) + 
        ylim(0,1) +
        xlab('h2_LEV') + 
        ylab('utility') +
        labs(col = "sample size",tag = letters[i])+
        theme(legend.margin=margin(0,0,0,0),
              legend.box.margin=margin(-10,0,-10,-10),
              plot.title = element_text(size = 10),
              panel.border = element_blank(),
              panel.background = element_blank(),
              panel.grid =element_blank(),
              axis.line = element_line(colour = "black"))
      i=i+1
      
      #--freq_rare--
      
      h2_poly=c(0.1,0.2,0.3,0.4,0.5)[3]
      h2_rare=c(0.01,0.02,0.03,0.05,0.10)[3]
      #freq_rare=c(1e-4,0.001,0.005,0.01,0.05)[3]
      prevalence=c(0.001,0.005,0.01,0.02,0.05)[3]
      #N_samples=c(500,1000,2000,3000,5000,10000)
      
      df_plot<-summary[summary$h2_rare==h2_rare & summary$prevalence==prevalence & summary$h2_poly==h2_poly,]
      
      single_var[[paste0(sample,'_',genetic_model,'_rare_freq')]]<-ggplot(df_plot, aes(x=factor(freq_rare), y=inverse_utility, colour=factor(N_samples),group=N_samples)) + 
        geom_line(size=1.5) + 
        geom_point(size=3,shape=17) + 
        ggtitle(paste0('h2_LEV=',h2_rare,', h2_cSEV=',h2_poly)) + 
        ylim(0,1) +
        xlab('freq_LEV (f)') + 
        ylab('utility') +
        labs(col = "sample size",tag = letters[i])+
        theme(legend.margin=margin(0,0,0,0),
              legend.box.margin=margin(-10,0,-10,-10),
              plot.title = element_text(size = 10),
              panel.border = element_blank(),
              panel.background = element_blank(),
              panel.grid =element_blank(),
              axis.line = element_line(colour = "black"))
      i=i+1
      
      
      #--prevalence--
      
      h2_poly=c(0.1,0.2,0.3,0.4,0.5)[3]
      h2_rare=c(0.01,0.02,0.03,0.05,0.10)[3]
      freq_rare=c(1e-4,0.001,0.005,0.01,0.05)[3]
      #prevalence=c(0.001,0.005,0.01,0.02,0.05)[3]
      #N_samples=c(500,1000,2000,3000,5000,10000)
      
      df_plot<-summary[summary$h2_rare==h2_rare & summary$freq_rare==freq_rare & summary$h2_poly==h2_poly,]
      
      single_var[[paste0(sample,'_',genetic_model,'_prevalence')]]<-ggplot(df_plot, aes(x=factor(prevalence), y=inverse_utility, colour=factor(N_samples),group=N_samples)) + 
        geom_line(size=1.5) + 
        geom_point(size=3,shape=17) + 
        ggtitle(paste0('h2_LEV=',h2_rare,' freq_LEV(f)=',freq_rare,', h2_cSEV=',h2_poly)) + 
        ylim(0,1) +
        xlab('prevalence (K)') + 
        ylab('utility') +
        labs(col = "sample size",tag = letters[i])+
        theme(legend.margin=margin(0,0,0,0),
              legend.box.margin=margin(-10,0,-10,-10),
              plot.title = element_text(size = 10),
              panel.border = element_blank(),
              panel.background = element_blank(),
              panel.grid =element_blank(),
              axis.line = element_line(colour = "black"))
      i=i+1
      
    }
  }
}

par(mfcol=c(2,2))
png(paste0('~/../Dropbox/POLYGENIC_LARGE_EFFECT/simulation/plots/utility/case_',disease_model,'_inf.png'),width = 2300,height = 2000,res = 275)
multiplot(single_var$case_a_h2_poly,single_var$case_a_h2_rare,single_var$case_a_rare_freq,single_var$case_a_prevalence,cols=2,by_row = T)
dev.off()



#----negative+ld probit---

i=1

for(sample in c('case')){ #,'control'
  
  for(disease_model in c('probit')){  #,'logit'
    
    single_var<-list()
    
    for(genetic_model in c('b','c')){ #,'b','c'
      
      if(genetic_model=='a'){
        genetic_model_label='Polygenic'
      }else if(genetic_model=='b'){
        genetic_model_label='Negative selection'
      }else{
        genetic_model_label='LD-adjusted kinship'
      }
      
      summary<-read.table(paste0('~/../Dropbox/POLYGENIC_LARGE_EFFECT/simulation/results/result_',disease_model,'_',genetic_model,'.txt'),head=T,stringsAsFactors = F)
      summary<-summary[,c(colnames(summary)[1:7],paste0(sample,'_inverse_power'))]
      colnames(summary)<-c('h2_poly','N_causal_SNPs','h2_rare','freq_rare','N_samples','prevalence','pi0','inverse_utility')
      
      #--h2_poly--
      
      #h2_poly=c(0.1,0.2,0.3,0.4,0.5)[3]
      h2_rare=c(0.01,0.02,0.03,0.05,0.10)[3]
      freq_rare=c(1e-4,0.001,0.005,0.01,0.05)[3]
      prevalence=c(0.001,0.005,0.01,0.02,0.05)[3]
      #N_samples=c(500,1000,2000,3000,5000,10000)
      
      df_plot<-summary[summary$h2_rare==h2_rare & summary$freq_rare==freq_rare & summary$prevalence==prevalence,]
      
      single_var[[paste0(sample,'_',genetic_model,'_h2_poly')]]<-ggplot(df_plot, aes(x=factor(h2_poly), y=inverse_utility, colour=factor(N_samples),group=N_samples)) + 
        geom_line(size=1.5) + 
        geom_point(size=3,shape=17) + 
        ggtitle(paste0(genetic_model_label,'\nfreq_LEV(f)=',freq_rare,', h2_LEV=',h2_rare)) + 
        ylim(0,1) +
        xlab('h2_cSEV') + 
        ylab('utility') +
        labs(col = "sample size",tag = letters[i])+
        theme(legend.margin=margin(0,0,0,0),
              legend.box.margin=margin(-10,0,-10,-10),
              plot.title = element_text(size = 10),
              panel.border = element_blank(),
              panel.background = element_blank(),
              panel.grid =element_blank(),
              axis.line = element_line(colour = "black"))
      i=i+1
      
      
      #--h2_rare--
      h2_poly=c(0.1,0.2,0.3,0.4,0.5)[3]
      #h2_rare=c(0.01,0.02,0.03,0.05,0.10)[3]
      freq_rare=c(1e-4,0.001,0.005,0.01,0.05)[3]
      prevalence=c(0.001,0.005,0.01,0.02,0.05)[3]
      #N_samples=c(500,1000,2000,3000,5000,10000)
      
      df_plot<-summary[summary$h2_poly==h2_poly & summary$freq_rare==freq_rare & summary$prevalence==prevalence,]
      
      single_var[[paste0(sample,'_',genetic_model,'_h2_rare')]]<-ggplot(df_plot, aes(x=factor(h2_rare), y=inverse_utility, colour=factor(N_samples),group=N_samples)) + 
        geom_line(size=1.5) + 
        geom_point(size=3,shape=17) + 
        ggtitle(paste0(genetic_model_label,'\nfreq_LEV(f)=',freq_rare,', h2_cSEV=',h2_poly)) + 
        ylim(0,1) +
        xlab('h2_LEV') + 
        ylab('utility') +
        labs(col = "sample size",tag = letters[i])+
        theme(legend.margin=margin(0,0,0,0),
              legend.box.margin=margin(-10,0,-10,-10),
              plot.title = element_text(size = 10),
              panel.border = element_blank(),
              panel.background = element_blank(),
              panel.grid =element_blank(),
              axis.line = element_line(colour = "black"))
      i=i+1
      
      #--freq_rare--
      
      h2_poly=c(0.1,0.2,0.3,0.4,0.5)[3]
      h2_rare=c(0.01,0.02,0.03,0.05,0.10)[3]
      #freq_rare=c(1e-4,0.001,0.005,0.01,0.05)[3]
      prevalence=c(0.001,0.005,0.01,0.02,0.05)[3]
      #N_samples=c(500,1000,2000,3000,5000,10000)
      
      df_plot<-summary[summary$h2_rare==h2_rare & summary$prevalence==prevalence & summary$h2_poly==h2_poly,]
      
      single_var[[paste0(sample,'_',genetic_model,'_rare_freq')]]<-ggplot(df_plot, aes(x=factor(freq_rare), y=inverse_utility, colour=factor(N_samples),group=N_samples)) + 
        geom_line(size=1.5) + 
        geom_point(size=3,shape=17) + 
        ggtitle(paste0(genetic_model_label,'\nh2_LEV=',h2_rare,', h2_cSEV=',h2_poly)) + 
        ylim(0,1) +
        xlab('freq_LEV (f)') + 
        ylab('utility') +
        labs(col = "sample size",tag = letters[i])+
        theme(legend.margin=margin(0,0,0,0),
              legend.box.margin=margin(-10,0,-10,-10),
              plot.title = element_text(size = 10),
              panel.border = element_blank(),
              panel.background = element_blank(),
              panel.grid =element_blank(),
              axis.line = element_line(colour = "black"))
      i=i+1
      
      
      #--prevalence--
      
      h2_poly=c(0.1,0.2,0.3,0.4,0.5)[3]
      h2_rare=c(0.01,0.02,0.03,0.05,0.10)[3]
      freq_rare=c(1e-4,0.001,0.005,0.01,0.05)[3]
      #prevalence=c(0.001,0.005,0.01,0.02,0.05)[3]
      #N_samples=c(500,1000,2000,3000,5000,10000)
      
      df_plot<-summary[summary$h2_rare==h2_rare & summary$freq_rare==freq_rare & summary$h2_poly==h2_poly,]
      
      single_var[[paste0(sample,'_',genetic_model,'_prevalence')]]<-ggplot(df_plot, aes(x=factor(prevalence), y=inverse_utility, colour=factor(N_samples),group=N_samples)) + 
        geom_line(size=1.5) + 
        geom_point(size=3,shape=17) + 
        ggtitle(paste0(genetic_model_label,'\nh2_LEV=',h2_rare,' freq_LEV(f)=',freq_rare,', h2_cSEV=',h2_poly)) + 
        ylim(0,1) +
        xlab('prevalence (K)') + 
        ylab('utility') +
        labs(col = "sample size",tag = letters[i])+
        theme(legend.margin=margin(0,0,0,0),
              legend.box.margin=margin(-10,0,-10,-10),
              plot.title = element_text(size = 10),
              panel.border = element_blank(),
              panel.background = element_blank(),
              panel.grid =element_blank(),
              axis.line = element_line(colour = "black"))
      i=i+1
      
    }
  }
}

par(mfcol=c(2,4))
png(paste0('~/../Dropbox/POLYGENIC_LARGE_EFFECT/simulation/plots/utility/case_',disease_model,'_negative+ld.png'),width = 2800,height = 1200,res = 180)
multiplot(single_var$case_b_h2_poly,
          single_var$case_b_h2_rare,
          single_var$case_b_rare_freq,
          single_var$case_b_prevalence,
          single_var$case_c_h2_poly,
          single_var$case_c_h2_rare,
          single_var$case_c_rare_freq,
          single_var$case_c_prevalence,
          cols=4,by_row = T)
dev.off()

#----negative+ld probit---

i=1

for(sample in c('case')){ #,'control'
  
  for(disease_model in c('logit')){  #,'logit'
    
    single_var<-list()
    
    for(genetic_model in c('a','b','c')){ #,'b','c'
      
      if(genetic_model=='a'){
        genetic_model_label='Polygenic'
      }else if(genetic_model=='b'){
        genetic_model_label='Negative selection'
      }else{
        genetic_model_label='LD-adjusted kinship'
      }
      
      summary<-read.table(paste0('~/../Dropbox/POLYGENIC_LARGE_EFFECT/simulation/results/result_',disease_model,'_',genetic_model,'.txt'),head=T,stringsAsFactors = F)
      summary<-summary[,c(colnames(summary)[1:7],paste0(sample,'_inverse_power'))]
      colnames(summary)<-c('h2_poly','N_causal_SNPs','h2_rare','freq_rare','N_samples','prevalence','pi0','inverse_utility')
      
      #--h2_poly--
      
      #h2_poly=c(0.1,0.2,0.3,0.4,0.5)[3]
      h2_rare=c(0.01,0.02,0.03,0.05,0.10)[3]
      freq_rare=c(1e-4,0.001,0.005,0.01,0.05)[3]
      prevalence=c(0.001,0.005,0.01,0.02,0.05)[3]
      #N_samples=c(500,1000,2000,3000,5000,10000)
      
      df_plot<-summary[summary$h2_rare==h2_rare & summary$freq_rare==freq_rare & summary$prevalence==prevalence,]
      
      single_var[[paste0(sample,'_',genetic_model,'_h2_poly')]]<-ggplot(df_plot, aes(x=factor(h2_poly), y=inverse_utility, colour=factor(N_samples),group=N_samples)) + 
        geom_line(size=1.5) + 
        geom_point(size=3,shape=17) + 
        ggtitle(paste0(genetic_model_label,'\nfreq_LEV(f)=',freq_rare,', h2_LEV=',h2_rare)) + 
        ylim(0,1) +
        xlab('h2_cSEV') + 
        ylab('utility') +
        labs(col = "sample size",tag = letters[i])+
        theme(legend.margin=margin(0,0,0,0),
              legend.box.margin=margin(-10,0,-10,-10),
              plot.title = element_text(size = 10),
              panel.border = element_blank(),
              panel.background = element_blank(),
              panel.grid =element_blank(),
              axis.line = element_line(colour = "black"))
      i=i+1
      
      #--h2_rare--
      h2_poly=c(0.1,0.2,0.3,0.4,0.5)[3]
      #h2_rare=c(0.01,0.02,0.03,0.05,0.10)[3]
      freq_rare=c(1e-4,0.001,0.005,0.01,0.05)[3]
      prevalence=c(0.001,0.005,0.01,0.02,0.05)[3]
      #N_samples=c(500,1000,2000,3000,5000,10000)
      
      df_plot<-summary[summary$h2_poly==h2_poly & summary$freq_rare==freq_rare & summary$prevalence==prevalence,]
      
      single_var[[paste0(sample,'_',genetic_model,'_h2_rare')]]<-ggplot(df_plot, aes(x=factor(h2_rare), y=inverse_utility, colour=factor(N_samples),group=N_samples)) + 
        geom_line(size=1.5) + 
        geom_point(size=3,shape=17) + 
        ggtitle(paste0(genetic_model_label,'\nfreq_LEV(f)=',freq_rare,', h2_cSEV=',h2_poly)) + 
        ylim(0,1) +
        xlab('h2_LEV') + 
        ylab('utility') +
        labs(col = "sample size",tag = letters[i])+
        theme(legend.margin=margin(0,0,0,0),
              legend.box.margin=margin(-10,0,-10,-10),
              plot.title = element_text(size = 10),
              panel.border = element_blank(),
              panel.background = element_blank(),
              panel.grid =element_blank(),
              axis.line = element_line(colour = "black"))
      i=i+1
      
      #--freq_rare--
      
      h2_poly=c(0.1,0.2,0.3,0.4,0.5)[3]
      h2_rare=c(0.01,0.02,0.03,0.05,0.10)[3]
      #freq_rare=c(1e-4,0.001,0.005,0.01,0.05)[3]
      prevalence=c(0.001,0.005,0.01,0.02,0.05)[3]
      #N_samples=c(500,1000,2000,3000,5000,10000)
      
      df_plot<-summary[summary$h2_rare==h2_rare & summary$prevalence==prevalence & summary$h2_poly==h2_poly,]
      
      single_var[[paste0(sample,'_',genetic_model,'_rare_freq')]]<-ggplot(df_plot, aes(x=factor(freq_rare), y=inverse_utility, colour=factor(N_samples),group=N_samples)) + 
        geom_line(size=1.5) + 
        geom_point(size=3,shape=17) + 
        ggtitle(paste0(genetic_model_label,'\nh2_LEV=',h2_rare,', h2_cSEV=',h2_poly)) + 
        ylim(0,1) +
        xlab('freq_LEV (f)') + 
        ylab('utility') +
        labs(col = "sample size",tag = letters[i])+
        theme(legend.margin=margin(0,0,0,0),
              legend.box.margin=margin(-10,0,-10,-10),
              plot.title = element_text(size = 10),
              panel.border = element_blank(),
              panel.background = element_blank(),
              panel.grid =element_blank(),
              axis.line = element_line(colour = "black"))
      i=i+1
      
      #--prevalence--
      
      h2_poly=c(0.1,0.2,0.3,0.4,0.5)[3]
      h2_rare=c(0.01,0.02,0.03,0.05,0.10)[3]
      freq_rare=c(1e-4,0.001,0.005,0.01,0.05)[3]
      #prevalence=c(0.001,0.005,0.01,0.02,0.05)[3]
      #N_samples=c(500,1000,2000,3000,5000,10000)
      
      df_plot<-summary[summary$h2_rare==h2_rare & summary$freq_rare==freq_rare & summary$h2_poly==h2_poly,]
      
      single_var[[paste0(sample,'_',genetic_model,'_prevalence')]]<-ggplot(df_plot, aes(x=factor(prevalence), y=inverse_utility, colour=factor(N_samples),group=N_samples)) + 
        geom_line(size=1.5) + 
        geom_point(size=3,shape=17) + 
        ggtitle(paste0(genetic_model_label,'\nh2_LEV=',h2_rare,' freq_LEV(f)=',freq_rare,', h2_cSEV=',h2_poly)) + 
        ylim(0,1) +
        xlab('prevalence (K)') + 
        ylab('utility') +
        labs(col = "sample size",tag = letters[i])+
        theme(legend.margin=margin(0,0,0,0),
              legend.box.margin=margin(-10,0,-10,-10),
              plot.title = element_text(size = 10),
              panel.border = element_blank(),
              panel.background = element_blank(),
              panel.grid =element_blank(),
              axis.line = element_line(colour = "black"))
      i=i+1
      
    }
  }
}

par(mfcol=c(3,4))
png(paste0('~/../Dropbox/POLYGENIC_LARGE_EFFECT/simulation/plots/utility/case_',disease_model,'_all.png'),width = 2800,height = 1800,res = 200)
multiplot(single_var$case_a_h2_poly,
          single_var$case_a_h2_rare,
          single_var$case_a_rare_freq,
          single_var$case_a_prevalence,
          single_var$case_b_h2_poly,
          single_var$case_b_h2_rare,
          single_var$case_b_rare_freq,
          single_var$case_b_prevalence,
          single_var$case_c_h2_poly,
          single_var$case_c_h2_rare,
          single_var$case_c_rare_freq,
          single_var$case_c_prevalence,
          cols=4,by_row = T)
dev.off()




#--------------------power------------------------------

library(gcookbook)
library(ggplot2)
library(RColorBrewer)
library(grid)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL,by_row=T) {
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
                     ncol = cols, nrow = ceiling(numPlots/cols),byrow = by_row)
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




power_plot<-function(genetic_model,freq_rare,h2_rare,prevalence,h2_poly){
  
  if(genetic_model=='a'){
    genetic_model_label='Polygenic'
  }else if(genetic_model=='b'){
    genetic_model_label='Negative selection'
  }else{
    genetic_model_label='LD-adjusted kinship'
  }
  
  summary<-read.table(paste0('~/../Dropbox/POLYGENIC_LARGE_EFFECT/simulation/results/power/result_',disease_model,'_',genetic_model,'_pi0_',level,'.txt'),head=T,stringsAsFactors = F)
  summary<-summary[,c('h2_p','N_causal_SNPs','vr','freq_r','N_samples','k','pi0',paste0(sample,'_inverse_power'))]
  colnames(summary)<-c('h2_poly','N_causal_SNPs','h2_rare','freq_rare','N_samples','prevalence','pi0','inverse_power')
  
  df_plot<-summary
  
  sub_plot<-ggplot(df_plot, aes(x=factor(pi0), y=inverse_power, colour=factor(N_samples),group=N_samples)) +
    geom_line(size=1.5) + #,position=position_dodge(0.2)
    geom_point(size=3,shape=17) + #position=position_dodge(0.2), 
    ggtitle(paste0(genetic_model_label,'\nfreq_LEV(f)=',freq_rare,', h2_LEV=',h2_rare,',\nh2_cSEV=',h2_poly,', prevalence(K)=',prevalence)) +
    ylim(0,1) +
    xlab('pi0') +
    ylab('power') +
    labs(col = "sample size",tag = letters[i])+
    theme(legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-10,0,-10,-10),
          plot.title = element_text(size = 10),
          panel.border = element_blank(),
          panel.background = element_blank(),
          panel.grid =element_blank(),
          axis.line = element_line(colour = "black"))
  i<<-i+1
  
  return(sub_plot)
}


#--probit--inf---

i=1

for (sample in c('case')){
  for(disease_model in c('probit')){
    plot_list<-list()
    for(genetic_model in c('a')){
      #very low
      level='very_low'
      #plot_list[[paste0(genetic_model,'_',level)]]<-power_plot(genetic_model=genetic_model,freq_rare=0.001,h2_rare=0.02,prevalence=0.005,h2_poly=0.2)
      
      #low
      level='low'
      plot_list[[paste0(genetic_model,'_',level)]]<-power_plot(genetic_model=genetic_model,freq_rare=0.005,h2_rare=0.03,prevalence=0.01,h2_poly=0.3)
      
      #median
      level='median'
      #plot_list[[paste0(genetic_model,'_',level)]]<-power_plot(genetic_model=genetic_model,freq_rare=0.01,h2_rare=0.05,prevalence=0.02,h2_poly=0.4)
      
      #high
      level='high'
      plot_list[[paste0(genetic_model,'_',level)]]<-power_plot(genetic_model=genetic_model,freq_rare=0.05,h2_rare=0.1,prevalence=0.05,h2_poly=0.5)
    }
  }
}

par(mfcol=c(1,2))
png(paste0('~/../Dropbox/POLYGENIC_LARGE_EFFECT/simulation/plots/power/',sample,'_power_',disease_model,'_inf.png'),width = 1500,height = 700,res = 175)
power_result<-multiplot(plot_list$a_low,
                        plot_list$a_high,
                        cols=2,by_row = F)
print(power_result)
dev.off()

#--probit--negative ld kinship---

i=1

for (sample in c('case')){
  for(disease_model in c('probit')){
    plot_list<-list()
    for(genetic_model in c('b','c')){
      #very low
      level='very_low'
      #plot_list[[paste0(genetic_model,'_',level)]]<-power_plot(genetic_model=genetic_model,freq_rare=0.001,h2_rare=0.02,prevalence=0.005,h2_poly=0.2)
      
      #low
      level='low'
      plot_list[[paste0(genetic_model,'_',level)]]<-power_plot(genetic_model=genetic_model,freq_rare=0.005,h2_rare=0.03,prevalence=0.01,h2_poly=0.3)
      
      #median
      level='median'
      #plot_list[[paste0(genetic_model,'_',level)]]<-power_plot(genetic_model=genetic_model,freq_rare=0.01,h2_rare=0.05,prevalence=0.02,h2_poly=0.4)
      
      #high
      level='high'
      plot_list[[paste0(genetic_model,'_',level)]]<-power_plot(genetic_model=genetic_model,freq_rare=0.05,h2_rare=0.1,prevalence=0.05,h2_poly=0.5)
    }
  }
}

par(mfcol=c(2,2))
png(paste0('~/../Dropbox/POLYGENIC_LARGE_EFFECT/simulation/plots/power/',sample,'_power_',disease_model,'_negative+ld.png'),width = 1500,height = 1300,res = 175)
power_result<-multiplot(plot_list$b_low,
                        plot_list$c_low,
                        plot_list$b_high,
                        plot_list$c_high,
                        cols=2,by_row = F)
print(power_result)
dev.off()



#--logit---

i=1

for (sample in c('case')){
  for(disease_model in c('logit')){
    plot_list<-list()
    for(genetic_model in c('a','b','c')){
      #very low
      level='very_low'
      #plot_list[[paste0(genetic_model,'_',level)]]<-power_plot(genetic_model=genetic_model,freq_rare=0.001,h2_rare=0.02,prevalence=0.005,h2_poly=0.2)
      
      #low
      level='low'
      plot_list[[paste0(genetic_model,'_',level)]]<-power_plot(genetic_model=genetic_model,freq_rare=0.005,h2_rare=0.03,prevalence=0.01,h2_poly=0.3)
      
      #median
      level='median'
      #plot_list[[paste0(genetic_model,'_',level)]]<-power_plot(genetic_model=genetic_model,freq_rare=0.01,h2_rare=0.05,prevalence=0.02,h2_poly=0.4)
      
      #high
      level='high'
      plot_list[[paste0(genetic_model,'_',level)]]<-power_plot(genetic_model=genetic_model,freq_rare=0.05,h2_rare=0.1,prevalence=0.05,h2_poly=0.5)
    }
  }
}

par(mfcol=c(3,2))
png(paste0('~/../Dropbox/POLYGENIC_LARGE_EFFECT/simulation/plots/power/',sample,'_power_',disease_model,'_negative+ld.png'),width = 1500,height = 2000,res = 175)
power_result<-multiplot(plot_list$a_low,
                        plot_list$b_low,
                        plot_list$c_low,
                        plot_list$a_high,
                        plot_list$b_high,
                        plot_list$c_high,
                        cols=2,by_row = F)
print(power_result)
dev.off()











#--------------------OR------------------------------


library(gcookbook)
library(ggplot2)
library(RColorBrewer)
library(grid)
library(reshape2)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL,by_row=T) {
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
                     ncol = cols, nrow = ceiling(numPlots/cols),byrow = by_row)
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


#------------------probit--logit----------------------

for(disease_model in c('probit','logit')){
  
  for(genetic_model in c('a','b','c')){
    
    i=1
    
    single_var<-list()
    
    if(genetic_model=='a'){
      genetic_model_label='Infinitesimal architecture'
    }else if(genetic_model=='b'){
      genetic_model_label='Negative selection'
    }else{
      genetic_model_label='LD-adjusted kinship'
    }
    
    summary<-read.table(paste0('~/../Dropbox/POLYGENIC_LARGE_EFFECT/simulation/results/OR/result_',disease_model,'_',genetic_model,'.txt'),head=T,stringsAsFactors = F)
    colnames(summary)[1:7]<-c('h2_poly','N_causal_SNPs','h2_rare','freq_rare','N_samples','prevalence','pi0')
    
    summary<-summary[which(summary$prevalence %in% c(0.005,0.01,0.02)),]
    
    
    #--h2_poly--
    
    #h2_poly=c(0.1,0.2,0.3,0.4,0.5)[3]
    h2_rare=c(0.01,0.02,0.03,0.05,0.10)[3]
    freq_rare=c(0.001,0.005,0.01)[2]
    prevalence=c(0.005,0.01,0.02)[2]
    N_samples=c(500,1000,2000,3000,5000,10000)[6]
    
    df_tmp<-summary[summary$h2_rare==h2_rare & summary$freq_rare==freq_rare & summary$prevalence==prevalence,]
    
    df_plot_beta<-melt(df_tmp,measure.vars = c('A_top_0.01_beta','A_top_0.05_beta','A_top_0.1_beta','R_beta'),id.vars = 'h2_poly')
    df_plot_se<-melt(df_tmp,measure.vars = c('A_top_0.01_se','A_top_0.05_se','A_top_0.1_se','R_se'),id.vars = 'h2_poly')
    n=5
    group=c(rep('cSEV top 1%',n),rep('cSEV top 5%',n),rep('cSEV top 10%',n),rep('LEV',n))
    df_plot<-data.frame(parameter=df_plot_beta[,1],group=group,beta=df_plot_beta$value,se=df_plot_se$value)
    df_plot$parameter<-paste0(' ',df_plot$parameter)
    df_plot$OR=2.718^df_plot$beta
    df_plot$conf.low=2.718^(df_plot$beta-1.96*df_plot$se)
    df_plot$conf.high=2.718^(df_plot$beta+1.96*df_plot$se)
    
    dodger = position_dodge(width = 0.3)
    
    df_plot$group = factor(df_plot$group,levels(df_plot$group)[c(4,1,3,2)]) 
    
    p_h2_poly<-ggplot(df_plot, aes(y = OR, x = parameter, colour = group)) +
      geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
                      position = dodger,
                      size = 0.8) +
      geom_hline(yintercept = 1.0, linetype = "dotted", size = 1) +
      scale_y_log10(breaks = c(0.5,1, 2, 5,10,20,50,100),
                    minor_breaks = NULL) +
      labs(y = "odds ratio (OR)", x = "h2_cSEV") +
      coord_flip(ylim = c(0.5, 100)) +
      theme_bw()+
      labs(tag = letters[i])+
      theme(legend.title = element_blank(),
            legend.margin=margin(0,0,0,0),
            legend.box.margin=margin(-10,0,-10,-10))
    i=i+1
    
    
    
    #--h2_rare--
    
    h2_poly=c(0.1,0.2,0.3,0.4,0.5)[3]
    #h2_rare=c(0.01,0.02,0.03,0.05,0.10)
    freq_rare=c(0.001,0.005,0.01)[2]
    prevalence=c(0.005,0.01,0.02)[2]
    N_samples=c(500,1000,2000,3000,5000,10000)[6]
    
    df_tmp<-summary[summary$h2_poly==h2_poly & summary$freq_rare==freq_rare & summary$prevalence==prevalence,]  ######
    
    df_plot_beta<-melt(df_tmp,measure.vars = c('A_top_0.01_beta','A_top_0.05_beta','A_top_0.1_beta','R_beta'),id.vars = 'h2_rare') #####
    df_plot_se<-melt(df_tmp,measure.vars = c('A_top_0.01_se','A_top_0.05_se','A_top_0.1_se','R_se'),id.vars = 'h2_rare') ########
    n=5 ########
    group=c(rep('cSEV top 1%',n),rep('cSEV top 5%',n),rep('cSEV top 10%',n),rep('LEV',n))
    df_plot<-data.frame(parameter=df_plot_beta[,1],group=group,beta=df_plot_beta$value,se=df_plot_se$value)
    df_plot$parameter<-paste0(' ',df_plot$parameter)
    df_plot$OR=2.718^df_plot$beta
    df_plot$conf.low=2.718^(df_plot$beta-1.96*df_plot$se)
    df_plot$conf.high=2.718^(df_plot$beta+1.96*df_plot$se)
    
    dodger = position_dodge(width = 0.3)
    
    df_plot$group = factor(df_plot$group,levels(df_plot$group)[c(4,1,3,2)]) 
    
    p_h2_rare<-ggplot(df_plot, aes(y = OR, x = parameter, colour = group)) +
      geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
                      position = dodger,
                      size = 0.8) +
      geom_hline(yintercept = 1.0, linetype = "dotted", size = 1) +
      scale_y_log10(breaks = c(1, 2, 5,10,20,50,100,1000),
                    minor_breaks = NULL) +
      labs(y = "odds ratio (OR)", x = "h2_LEV") +
      coord_flip(ylim = c(0.5, 5000)) +
      theme_bw()+
      labs(tag = letters[i])+
      theme(legend.title = element_blank(),
            legend.margin=margin(0,0,0,0),
            legend.box.margin=margin(-10,0,-10,-10))
    i=i+1
    
    
    #--freq--
    h2_poly=c(0.1,0.2,0.3,0.4,0.5)[3]
    h2_rare=c(0.01,0.02,0.03,0.05,0.10)[3]
    #freq_rare=c(0.001,0.005,0.01)[2]
    prevalence=c(0.005,0.01,0.02)[2]
    N_samples=c(500,1000,2000,3000,5000,10000)[6]
    
    df_tmp<-summary[summary$h2_poly==h2_poly & summary$h2_rare==h2_rare & summary$prevalence==prevalence,]  ######
    
    df_plot_beta<-melt(df_tmp,measure.vars = c('A_top_0.01_beta','A_top_0.05_beta','A_top_0.1_beta','R_beta'),id.vars = 'freq_rare') #####
    df_plot_se<-melt(df_tmp,measure.vars = c('A_top_0.01_se','A_top_0.05_se','A_top_0.1_se','R_se'),id.vars = 'freq_rare') ########
    n=3 ########
    group=c(rep('cSEV top 1%',n),rep('cSEV top 5%',n),rep('cSEV top 10%',n),rep('LEV',n))
    df_plot<-data.frame(parameter=df_plot_beta[,1],group=group,beta=df_plot_beta$value,se=df_plot_se$value)
    df_plot$parameter<-paste0(' ',df_plot$parameter)
    df_plot$OR=2.718^df_plot$beta
    df_plot$conf.low=2.718^(df_plot$beta-1.96*df_plot$se)
    df_plot$conf.high=2.718^(df_plot$beta+1.96*df_plot$se)
    
    dodger = position_dodge(width = 0.3)
    
    df_plot$group = factor(df_plot$group,levels(df_plot$group)[c(4,1,3,2)]) 
    
    p_freq_rare<-ggplot(df_plot, aes(y = OR, x = parameter, colour = group)) +
      geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
                      position = dodger,
                      size = 0.8) +
      geom_hline(yintercept = 1.0, linetype = "dotted", size = 1) +
      scale_y_log10(breaks = c(1, 2, 5,10,20,50,100,1000,10000),
                    minor_breaks = NULL) +
      labs(y = "odds ratio (OR)", x = "freq_LEV (f)") +
      coord_flip(ylim = c(0.5, 10000)) +
      theme_bw()+
      labs(tag = letters[i])+
      theme(legend.title = element_blank(),
            legend.margin=margin(0,0,0,0),
            legend.box.margin=margin(-10,0,-10,-10))
    i=i+1
    
    #--k--
    h2_poly=c(0.1,0.2,0.3,0.4,0.5)[3]
    h2_rare=c(0.01,0.02,0.03,0.05,0.10)[3]
    freq_rare=c(0.001,0.005,0.01)[2]
    #prevalence=c(0.005,0.01,0.02)[2]
    N_samples=c(500,1000,2000,3000,5000,10000)[6]
    
    df_tmp<-summary[summary$h2_poly==h2_poly & summary$h2_rare==h2_rare & summary$freq_rare==freq_rare,]  ######
    
    df_plot_beta<-melt(df_tmp,measure.vars = c('A_top_0.01_beta','A_top_0.05_beta','A_top_0.1_beta','R_beta'),id.vars = 'prevalence') #####
    df_plot_se<-melt(df_tmp,measure.vars = c('A_top_0.01_se','A_top_0.05_se','A_top_0.1_se','R_se'),id.vars = 'prevalence') ########
    n=3 ########
    group=c(rep('cSEV top 1%',n),rep('cSEV top 5%',n),rep('cSEV top 10%',n),rep('LEV',n))
    df_plot<-data.frame(parameter=df_plot_beta[,1],group=group,beta=df_plot_beta$value,se=df_plot_se$value)
    df_plot$parameter<-paste0(' ',df_plot$parameter)
    df_plot$OR=2.718^df_plot$beta
    df_plot$conf.low=2.718^(df_plot$beta-1.96*df_plot$se)
    df_plot$conf.high=2.718^(df_plot$beta+1.96*df_plot$se)
    
    dodger = position_dodge(width = 0.3)
    
    df_plot$group = factor(df_plot$group,levels(df_plot$group)[c(4,1,3,2)]) 
    
    p_prevalence<-ggplot(df_plot, aes(y = OR, x = parameter, colour = group)) +
      geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
                      position = dodger,
                      size = 0.8) +
      geom_hline(yintercept = 1.0, linetype = "dotted", size = 1) +
      scale_y_log10(breaks = c(1, 2, 5,10,20,50,100),
                    minor_breaks = NULL) +
      labs(y = "odds ratio (OR)", x = "prevalence (K)") +
      coord_flip(ylim = c(0.5, 100)) +
      theme_bw()+
      labs(tag = letters[i])+
      theme(legend.title = element_blank(),
            legend.margin=margin(0,0,0,0),
            legend.box.margin=margin(-10,0,-10,-10))
    i=i+1
    
    
    #---plot---
    #par(mfcol=c(2,2))
    png(paste0('~/../Dropbox/POLYGENIC_LARGE_EFFECT/simulation/plots/OR/',disease_model,'_',genetic_model_label,'.png'),width = 1500,height = 1200,res = 175)
    p<-multiplot(p_h2_poly,p_h2_rare,p_freq_rare,p_prevalence,cols=2,by_row = T)
    print(p)
    dev.off()
  }
}






















