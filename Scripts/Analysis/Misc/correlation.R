cluster_names_all_10x <- c("TEM", "TTE", "TN", "Cyto_TEM", "PRE_EX", "TCM")
names_to_check = c("ref_HPCA_all_10x", "ref_BE_all_10x", "ref_DICE_all_10x", "ref_NH_all_10x", "ref_Monaco_all_10x")



temp = as.data.frame(table(all_10x@meta.data$ref_HPCA_all_10x, all_10x$Cluster)) %>% 
  pivot_wider(names_from = Var2, values_from = Freq)
temp2 = temp$Var1
temp = temp[,-1]
rownames(temp) = temp2
colnames(temp) = cluster_names_all_10x

list_o = list()





for (i in names(names_to_check)) {
  temp = as.data.frame(table(all_10x@meta.data[[`i`]], all_10x$Cluster)) %>% 
    pivot_wider(names_from = Var2, values_from = Freq)
  temp2 = temp$Var1
  temp = temp[,-1]
  rownames(temp) = temp2
  colnames(temp) = cluster_names_all_10x
  temp = t(temp)
  list_o[[length(list_o) +1]] <- temp
}



temp2 = cor(temp)




temp2 = temp %>% pivot_wider(names_from = Var2, values_from = Freq)
