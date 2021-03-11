library(ggplot2)

h2_poly=0.5 #h2 Poly
h2_levs=0.1 #h2 LEV (total)
n_levs=5  #number of LEVs
freq_snp=0.10 #freq of LEVs
h2_levs_each=h2_levs/n_levs #h2 lev for each 
k=0.01 #prevalence
sample_size=10000 #sample size

#generate the LEV burden for each LEV
R_single_func<-function(){
  dosage=rbinom(n = 1, size = 2, prob = freq_snp)
  effect=(h2_levs_each/2/freq_snp/(1-freq_snp))^0.5
  noise=rnorm(n=1,mean=0,sd=0.002) #make a bit noise to sightly shuffle the dots up and down
  return(effect*dosage+noise)
}

#sum up the LEV burden to get the R
R_sum_func<-function(){
  sum(replicate(n_levs,R_single_func()))
}

#LEV
set.seed(2020)
R<-replicate(sample_size,R_sum_func()) 
#Poly
A<-rnorm(sample_size, mean = 0, sd=h2_poly^0.5) 

df<-data.frame(A=A,R=R)
df$AR=df$A+df$R

h2_e=1-h2_poly-h2_levs 

#Liability
set.seed(2021)
df$L=df$AR+rnorm(nrow(df),0,h2_e^0.5)

#threshold
threshold=as.numeric(quantile(df$L,(1-k)))

p_hist<-ggplot(df, aes(x=AR)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666",bw=0.3) +
  geom_vline(aes(xintercept=threshold),color="red", linetype="dashed", size=1)+
  xlab('Disease liability')+
  ylab('Density')+
  theme(legend.position= "none",
        panel.background = element_blank(),
        panel.grid =element_blank(),
        axis.line = element_line(colour = "black"))+
  
  labs(tag='a')


df$color=ifelse((df$L)>threshold,'case','control')

case<-df[df$color=='case',]
control<-df[df$color=='control',]

p_scatter<-ggplot(df, aes(x=A, y=R)) + 
  geom_point(data = df, aes(x=A, y=R,color=color),size=0.3)+
  scale_color_manual(values=c('case'='red','control'='black'))+
  
  geom_point(data = case, aes(x=A, y=R),color='red',size=0.3) + 
  #geom_smooth(method = "lm", data = case,color='red',size=0.1)+  
  
  xlab('Polygenic Burden (PB)')+
  ylab('Genetic burden due to \nLarge-Effect Variants (LEVs)')+
  
  theme(legend.position="top",
        legend.title=element_blank(),
        panel.background = element_blank(),
        panel.grid =element_blank(),
        axis.line = element_line(colour = "black"))+
  
  labs(tag='b')


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



#---plot---
pdf(paste0('~/../Desktop/Figure 1.pdf'),width = 7,height = 4)
multiplot(p_hist,p_scatter,cols=2,by_row = T)
dev.off()


