mapping_df <- data.frame(sample=c("cat_372-1","cat_372-2","cat_372-3","cat_372-4"),
                         input_reads=c(31601662,32027326,35651767,29949567),
                         unique_mapped=c(28089846,28234031,30565908,25657698),
                         multiple_mapped=c(1497164,1459920,1581943,1178593),
                         too_many_mapped=c(4201,4127,4580,4074))

mapping_df$unmapped <- mapping_df$input_reads-mapping_df$unique_mapped-
  mapping_df$multiple_mapped-mapping_df$too_many_mapped

mapping_long <- gather(mapping_df,condition,count,unique_mapped:unmapped,factor_key=TRUE)

mapping_long$percent <- mapping_long$count/mapping_long$input_reads*100

mapping_long$percent_lab <- paste0(sprintf("%.2f",
                                           mapping_long$count/mapping_long$input_reads*100),
                                   "%")

ggplot(mapping_long,aes(x=sample,y=percent,fill=condition))+
  geom_col()+
  #scale_y_continuous(labels = scales::percent)+
  geom_text(aes(label=percent_lab), position=position_stack(vjust=0.5), size=4)
