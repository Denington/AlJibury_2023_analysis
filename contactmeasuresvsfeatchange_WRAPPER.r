mults <- T

#
cOff.U <- 0.08
cOff.O <- 0.05

source('~/Documents/Form_function_scripts/contactmeasures_vs_featurechange.r')

outstats$cOff.O <- cOff.O ; outstats$cOff.U <- cOff.U ; outtab$cOff.O <- cOff.O ; outtab$cOff.U <- cOff.U

outS1 <- outstats
outT1 <- outtab
#
cOff.U <- 0.6
cOff.O <- 0.05

source('~/Documents/Form_function_scripts/contactmeasures_vs_featurechange.r')
outstats$cOff.O <- cOff.O ; outstats$cOff.U <- cOff.U ; outtab$cOff.O <- cOff.O ; outtab$cOff.U <- cOff.U

outS1 <- rbind(outS1,outstats)
outT1 <- rbind(outT1,outtab)


#
cOff.O <- 0.001
for(i in c(0.0001,0.001,0.01,0.05,0.075,0.1,0.25)){
  print(i)
  cOff.U <- i
  source('~/Documents/Form_function_scripts/contactmeasures_vs_featurechange.r')
  outstats$cOff.O <- cOff.O ; outstats$cOff.U <- cOff.U ; outtab$cOff.O <- cOff.O ; outtab$cOff.U <- cOff.U
  outS1 <- rbind(outS1,outstats)
  outT1 <- rbind(outT1,outtab)
  
}

#
cOff.O <- 0.01
for(i in c(0.0001,0.001,0.01,0.05,0.075,0.1,0.25)){
  print(i)
  cOff.U <- i
  source('~/Documents/Form_function_scripts/contactmeasures_vs_featurechange.r')
  outstats$cOff.O <- cOff.O ; outstats$cOff.U <- cOff.U ; outtab$cOff.O <- cOff.O ; outtab$cOff.U <- cOff.U
  outS1 <- rbind(outS1,outstats)
  outT1 <- rbind(outT1,outtab)
  
}

#
cOff.O <- 0.05
for(i in c(0.0001,0.001,0.01,0.05,0.075,0.1,0.25)){
  print(i)
  cOff.U <- i
  source('~/Documents/Form_function_scripts/contactmeasures_vs_featurechange.r')
  outstats$cOff.O <- cOff.O ; outstats$cOff.U <- cOff.U ; outtab$cOff.O <- cOff.O ; outtab$cOff.U <- cOff.U
  outS1 <- rbind(outS1,outstats)
  outT1 <- rbind(outT1,outtab)
  
}

#
cOff.O <- 0.1
for(i in c(0.0001,0.001,0.01,0.05,0.075,0.1,0.25)){
  print(i)
  cOff.U <- i
  source('~/Documents/Form_function_scripts/contactmeasures_vs_featurechange.r')
  outstats$cOff.O <- cOff.O ; outstats$cOff.U <- cOff.U ; outtab$cOff.O <- cOff.O ; outtab$cOff.U <- cOff.U
  outS1 <- rbind(outS1,outstats)
  outT1 <- rbind(outT1,outtab)
  
}

#
cOff.O <- 0.25
for(i in c(0.0001,0.001,0.01,0.05,0.075,0.1,0.25)){
  print(i)
  cOff.U <- i
  source('~/Documents/Form_function_scripts/contactmeasures_vs_featurechange.r')
  outstats$cOff.O <- cOff.O ; outstats$cOff.U <- cOff.U ; outtab$cOff.O <- cOff.O ; outtab$cOff.U <- cOff.U
  outS1 <- rbind(outS1,outstats)
  outT1 <- rbind(outT1,outtab)
  
}

#####

outS1 <- outS1[outS1$cOff.U <= outS1$cOff.O,]
outT1 <- outT1[outT1$cOff.U <= outT1$cOff.O,]

write.table(outS1,'~/Desktop/feature_3Dchange_stats.tsv',sep='\t',
            quote = F,row.names = F)

write.table(outT1,'~/Desktop/feature_3Dchange_counts.tsv',sep='\t',
            quote = F,row.names = F)
