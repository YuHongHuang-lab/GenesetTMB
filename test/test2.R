nettop_feature2_TCGA = as.data.frame(nettop_feature2_TCGA)
colnames(nettop_feature2_TCGA) = c('alphen', 'cluster_coef')
nettop_feature2_TCGA$pathid = str_remove(rownames(nettop_feature2_TCGA), 'HALLMARK_')
ggplot(nettop_feature2_TCGA, aes(alphen, cluster_coef))+geom_point(alpha=0.8, color='#E64B35FF')+
  theme_cowplot()+theme(legend.title = element_blank(),
                        legend.position = 'none')+
  geom_vline(xintercept = 0.05, linetype='dotdash', color='black')+
  geom_hline(yintercept = 0.05, linetype='dotdash', color='black')
