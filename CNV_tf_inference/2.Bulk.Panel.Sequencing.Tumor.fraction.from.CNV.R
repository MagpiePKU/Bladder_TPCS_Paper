
sid = commandArgs(T)[1]
cnr_file <- commandArgs(T)[2]


library(ggplot2)
library(dplyr)
library(reshape2)
library(tidyr)
library(readr)
library(GenomicRanges)
library(bezier)
library(ggforce)
library(RColorBrewer)
library(ggthemes)
library(cowplot)
library(mixtools)
library(mclust)

get_fraction_df = function(df){
  df$stat = 0
  if(min(abs(df$index)) == 0){
    df[df$index == 0,]$stat = 2
    df[df$index != 0,]$stat = 3
  }else{
    df$stat = df$index / min(abs(df$index))
    df[abs(df$index) == min(abs(df$index)),]$stat = 2
  }
  
  
  info = (df[df$stat == 2,]$mu - df[df$stat == 2,]$diff_between_steps)[1]
  info1 = (df[df$stat == 2,]$mu - df[df$stat == 2,]$diff_between_steps/2)[1]
  df[df$mu <= info1,]$stat = 1
  if (nrow(df[df$mu <= info, ]) > 0){
    df[df$mu <= info, ]$stat = 0
  }
  
  cn1 = (df[df$stat == 2,]$diff_between_steps)[1] / 2
  df$cn = round(df$diff_between_steps / cn1)
  
  df[df$stat <2,]$cn = 0
  
  df$cn_sum = cumsum(df$cn)
  df$cn_state_mod = df$cn_sum
  if (nrow(df[df$cn_sum <2 & df$stat != df$cn_sum,]) > 0){
    df[df$cn_sum <2 & df$stat != df$cn_sum,]$cn_state_mod = df[df$cn_sum <2 & df$stat != df$cn_sum,]$stat
  }
  return(df)
}


cnrdat = read.delim(cnr_file)
X <- data.matrix(cnrdat$log2)

selected_k_list <- list()
for (i in seq(1:10)){
  mod <- Mclust(X)
  selected_k_list[[i]] = mod$G
  
}
#names(selected_k_list) = 'k'
selected_k_df = t(bind_cols(selected_k_list))
selected_k_df_freq = as.data.frame(table(selected_k_df))
selected_k =as.numeric(as.character(selected_k_df_freq[max(selected_k_df_freq$Freq) == selected_k_df_freq$Freq,]$selected_k_df))

mod <- Mclust(X,G=selected_k)
df = as.data.frame(mod$parameters$mean)
colnames(df) <- 'mu2'
df <- arrange(df,mu2)

df$mu = 2^(df$mu2 + 1)
df$last_step_mu <- c(2,df$mu[1:(nrow(df)-1)])
df$diff_between_steps = df$mu - df$last_step_mu 
df$index = round(df$mu,1) - 2 

df = get_fraction_df(df)

lm(df,formula =mu ~ cn_state_mod) -> linear_model
linear_summary = summary(linear_model)
coefficients = linear_summary$coefficients[2]
r2 = linear_summary$adj.r.squared
b = linear_summary$coefficients[1]
b_info = 1 - b / 2

if(r2<0.85 & !is.na(r2)){
  df = df[-1,]
  df = get_fraction_df(df)
  lm(df,formula =mu ~ cn_state_mod) -> linear_model
  linear_summary = summary(linear_model)
  coefficients = linear_summary$coefficients[2]
  r2 = linear_summary$adj.r.squared
  b = linear_summary$coefficients[1]
  b_info = 1 - b / 2
}

result_output = data.frame(sample = sid,coefficients = coefficients,r2 =r2,cn_stat_mod = paste(df$cn_state_mod,collapse = ","),k=selected_k,Intercept=b,Intercept_info = b_info)

output_file = paste0("./result/",sid,".txt")
output_file_df = paste0("./result/",sid,"_mu_df.info")
write.table(result_output,output_file,quote=F,sep = "\t",row.names = FALSE)
write.table(df,output_file_df,quote=F,sep = "\t",row.names = FALSE)


pdf(paste0("./result/",sid,".pdf"),width=5,height=3.5)
ggplot(df,aes(cn_state_mod,mu)) + geom_point() + geom_smooth(method = 'lm') + theme_bw() + xlab('Absolute Copy Number') + ylab('Observed CNV') + ggtitle(sid,subtitle = paste0('Tumor Fraction by CNV = ',round(linear_model$coefficients[2]*100,digits = 3),'%'))
dev.off()

outputjson = jsonlite::toJSON(result_output,pretty=T)
sink(file=paste0("./result/",sid,'.cnv_fraction.json'))
cat(outputjson)
sink()


#ggplot(df,aes(cn_state_mod,mu)) + geom_point() + geom_smooth(method = 'lm') + theme_bw() + xlab('Absolute Copy Number') + ylab('Observed CNV') + ggtitle(sid,subtitle = paste0('Tumor Fraction by CNV = ',round(linear_model$coefficients[2]*100,digits = 3),'%'))



