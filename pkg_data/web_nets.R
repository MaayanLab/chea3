
net = chea3::libs2network()
net = chea3::toTFNet(net)
net = plyr::ddply(net,plyr::.(lib),function(sub){
  min = min(sub$weight)
  max = max(sub$weight)
  sub$weight = scale(sub$weight, center = min, scale = max-min)
  return(sub)
})
net[is.na(net$weight),"weight"] = 0.5
net1 = plyr::ddply(net,plyr::.(source,target),function(sub){
  return(data.frame(
    source = unique(sub$source),
    target = unique(sub$target),
    weight_sum = sum(sub$weight, na.rm = T),
    libs = paste(sub$lib,collapse = ","),
    num_libs = length(unique(sub$lib)),
    stringsAsFactors = F
  ))
})

net1$score = net1$weight_sum*net1$num_libs

net2 = net1[net1$num_libs>=2,]
n = 0.02
#take top 10% of targets for each source
net3 = plyr::ddply(net2,plyr::.(source),function(sub){
  sub$rank = rank(sub$weight_sum, ties.method = "random")
  sub = sub[order(sub$rank, decreasing = T),][1:ceiling((n*nrow(sub))),]
  return(sub)
})

setwd("/users/alexandrakeenan/Desktop/")
write.table(data.frame(source = net3$source,
  target = net3$target, stringsAsFactors = F),"tempnet.tsv",sep = "\t",quote = F, col.names = T, row.names = F)
