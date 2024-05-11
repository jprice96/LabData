setwd("/Users/Jack/Dropbox/Phylo")

# Load in data
orf.matrix = read.csv(file = "ORFidentitymatrix_hattci_5Dec21.csv")
orf.matrix = orf.matrix[,-c(1)]
orf.list = read.delim(file = "ORFresults_hattci_15Nov21.txt", header = TRUE, sep = " ")

# Create species/strain categories
species.list = data.frame(species = unlist(strsplit(as.character(orf.list$species), split = "_"))[seq(1, length(orf.list$species)*3, 3)],
                          strain = unlist(strsplit(as.character(orf.list$species), split = "_"))[seq(2, length(orf.list$species)*3, 3)])

# Convert to longform
longform.orf = data.frame(
                          species1 = rep(species.list$species, times = length(orf.matrix[,1])),
                          strain1 = rep(species.list$strain, times = length(orf.matrix[,1])),
                          orf1 = rep(c(1:length(orf.matrix[,1])), times = length(orf.matrix[,1])),
                          identity = unlist(orf.matrix, use.names = FALSE),
                          species2 = rep(species.list$species, each = length(orf.matrix[1,])),
                          strain2 = rep(species.list$strain, each = length(orf.matrix[1,])),
                          orf2 = rep(c(1:length(orf.matrix[,1])), each = length(orf.matrix[1,]))
                          )

longform.orf = longform.orf[complete.cases(longform.orf),]

save(longform.orf, file = "longform_orf_table_15Nov.csv")
load(file = "longform_orf_table.csv")


paralog.orf = longform.orf[which(longform.orf$identity > 0.75),]
common.orf = table(paralog.orf$orf1)
common.orf = common.orf[which(common.orf>10)]

orf.identity = paralog.orf
orf.identity$orfid = NA

for(i in 1:length(orf.identity$orfid)){
  if(is.na(orf.identity$orfid[i])==TRUE){
    orf.identity$orfid[i] = i
  }
  orf.identity$orfid[which(orf.identity$orf1 == orf.identity$orf1[i])] = orf.identity$orfid[i]
  orf.identity$orfid[which(orf.identity$orf2 == orf.identity$orf1[i])] = orf.identity$orfid[i]
  orf.identity$orfid[which(orf.identity$orf1 == orf.identity$orf2[i])] = orf.identity$orfid[i]
  orf.identity$orfid[which(orf.identity$orf2 == orf.identity$orf2[i])] = orf.identity$orfid[i]
}

orf.identity$orfid = as.numeric(as.factor(orf.identity$orfid))

common.orf = table(orf.identity$orfid)
common.orf = common.orf[which(common.orf>=10)]

common.orf = orf.identity[which(orf.identity$orfid %in% names(common.orf)),]
common.orf = common.orf[order(common.orf$orfid),]
common.orf = common.orf[,-c(4:7)]
common.orf = common.orf[!duplicated(common.orf$orf1),]

recount.orf = table(common.orf$orfid)
recount.orf = recount.orf[which(recount.orf>=10)]
common.orf = common.orf[which(common.orf$orfid %in% names(recount.orf)),]
common.orf$orfid = as.numeric(as.factor(common.orf$orfid))

## common.orf is now a data.frame indicating each orf and it's paralog id where the paralog appears at least 10 times (incl dupl within strains)
## n = 139 paralogs

full.data = common.orf
full.data$salt = NA
full.data$habt = NA
full.data$temp = NA
full.data$blum = NA
full.data$symb = NA
        #Af As Aw Sk Va Va Vc Vc Vc Vc Vf Vf Vh Vj Vm Vp Vr Vs Vt Vv
salt = c(0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0)  # 0=normal, 1=high salinity
habt = c(0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1)  # 0=marine, 1=brackish
temp = c(0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0)  # 0=mesophilic, 1=psycrophilic
blum = c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1)  # 0=normal, 1=bioluminescent
symb = c(1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0)  # 0=normal, 1=symbiotic

for(i in 1:length(unique(common.orf$species1))){
  full.data$salt[which(full.data$species1 == levels(full.data$species1)[i])] = salt[i]
  full.data$habt[which(full.data$species1 == levels(full.data$species1)[i])] = habt[i]
  full.data$temp[which(full.data$species1 == levels(full.data$species1)[i])] = temp[i]
  full.data$blum[which(full.data$species1 == levels(full.data$species1)[i])] = blum[i]
  full.data$symb[which(full.data$species1 == levels(full.data$species1)[i])] = symb[i]
}

save(full.data, file = "Paralogs_wEco.csv")
load(file = "Paralogs_wEco.csv")

salt.data = NA
habt.data = NA
temp.data = NA
blum.data = NA
symb.data = NA

for(i in 1:length(unique(full.data$orfid))){
  store.data = full.data[which(full.data$orfid == unique(full.data$orfid)[i]),]
  a = length(unique(store.data$species1[which(store.data$salt==1)]))
  b = length(unique(store.data$species1[which(store.data$salt==0)]))
  c = length(unique(full.data$species1[which(full.data$salt==1)])) - a
  d = length(unique(full.data$species1[which(full.data$salt==0)])) - b
  store.matrix = matrix(data = c(a,b,c,d), nrow=2)
  salt.data[i] = unlist(chisq.test(store.matrix)[3])
  
  a = length(unique(store.data$species1[which(store.data$habt==1)]))
  b = length(unique(store.data$species1[which(store.data$habt==0)]))
  c = length(unique(full.data$species1[which(full.data$habt==1)])) - a
  d = length(unique(full.data$species1[which(full.data$habt==0)])) - b
  store.matrix = matrix(data = c(a,b,c,d), nrow=2)
  habt.data[i] = unlist(chisq.test(store.matrix)[3])
  
  a = length(unique(store.data$species1[which(store.data$temp==1)]))
  b = length(unique(store.data$species1[which(store.data$temp==0)]))
  c = length(unique(full.data$species1[which(full.data$temp==1)])) - a
  d = length(unique(full.data$species1[which(full.data$temp==0)])) - b
  store.matrix = matrix(data = c(a,b,c,d), nrow=2)
  temp.data[i] = unlist(chisq.test(store.matrix)[3])
  
  a = length(unique(store.data$species1[which(store.data$blum==1)]))
  b = length(unique(store.data$species1[which(store.data$blum==0)]))
  c = length(unique(full.data$species1[which(full.data$blum==1)])) - a
  d = length(unique(full.data$species1[which(full.data$blum==0)])) - b
  store.matrix = matrix(data = c(a,b,c,d), nrow=2)
  blum.data[i] = unlist(chisq.test(store.matrix)[3])
  
  a = length(unique(store.data$species1[which(store.data$symb==1)]))
  b = length(unique(store.data$species1[which(store.data$symb==0)]))
  c = length(unique(full.data$species1[which(full.data$symb==1)])) - a
  d = length(unique(full.data$species1[which(full.data$symb==0)])) - b
  store.matrix = matrix(data = c(a,b,c,d), nrow=2)
  symb.data[i] = unlist(chisq.test(store.matrix)[3])
}

which(habt.data<=0.05)
#Salt Significant ORFID: 21 (99% are VALGI)
#Habt Significant ORFID: 128 (could be something)
#Temp Significant ORFID: 50, 60, 61, 88 (each could be something)
#Blum Significant ORFID: 3 (could be something), 8 (could be something)
#Symb Significant ORFID: 3 (could be something), 8 (could be something), 21 (99% are VALGI)

length(full.data[which(full.data$orfid==50),7])
#3 = AFISC ES114 orf-26
#8 = AFISC ES114 orf-17
#50 = VALGI VB1833 orf-1226
#60 = VALGI VB1833 orf-1205
#61 = VALGI VB1833 orf-1206
#88 = VANGU 879116 orf-1460
#128 = VANGU 425 orf-1363

#21 = VALGI K01m1 orf-587



### ARE PHYLO.DIST AND ENVIRONMENT CORRELATED?
library(ecodist)

dist.matrix = read.csv(file = "PhyloDist_Core23.csv")
dist.matrix = as.matrix(dist.matrix[,-c(1)])
rownames(dist.matrix) = colnames(dist.matrix)
dist.matrix = dist.matrix[4:23, 4:23]   # Remove the photobacterium which aren't in full data
species.list = row.names(dist.matrix)
dist.matrix = as.dist(dist.matrix)

# Create blank matrix template
env.matrix = matrix(0, ncol = length(species.list), nrow = length(species.list))
colnames(env.matrix) = species.list
rownames(env.matrix) = species.list

# Create a disimilarity matrix for each environment
salt.matrix = env.matrix
for(i in 1:length(species.list)){
  for(j in 1:length(species.list)){
    species1 = tail(full.data$salt[which(full.data$species1 == species.list[i])],1)
    species2 = tail(full.data$salt[which(full.data$species1 == species.list[j])],1)
    salt.matrix[i,j] = abs(species1 - species2)
  }
}
salt.matrix = as.dist(salt.matrix)
mantel.salt = mantel(salt.matrix ~ dist.matrix, nperm = 100)

habt.matrix = env.matrix
for(i in 1:length(species.list)){
  for(j in 1:length(species.list)){
    species1 = tail(full.data$habt[which(full.data$species1 == species.list[i])],1)
    species2 = tail(full.data$habt[which(full.data$species1 == species.list[j])],1)
    habt.matrix[i,j] = abs(species1 - species2)
  }
}
habt.matrix = as.dist(habt.matrix)
mantel.habt = mantel(habt.matrix ~ dist.matrix, nperm = 100)

temp.matrix = env.matrix
for(i in 1:length(species.list)){
  for(j in 1:length(species.list)){
    species1 = tail(full.data$temp[which(full.data$species1 == species.list[i])],1)
    species2 = tail(full.data$temp[which(full.data$species1 == species.list[j])],1)
    temp.matrix[i,j] = abs(species1 - species2)
  }
}
temp.matrix = as.dist(temp.matrix)
mantel.temp = mantel(temp.matrix ~ dist.matrix, nperm = 100)

blum.matrix = env.matrix
for(i in 1:length(species.list)){
  for(j in 1:length(species.list)){
    species1 = tail(full.data$blum[which(full.data$species1 == species.list[i])],1)
    species2 = tail(full.data$blum[which(full.data$species1 == species.list[j])],1)
    blum.matrix[i,j] = abs(species1 - species2)
  }
}
blum.matrix = as.dist(blum.matrix)
mantel.blum = mantel(blum.matrix ~ dist.matrix, nperm = 100)

symb.matrix = env.matrix
for(i in 1:length(species.list)){
  for(j in 1:length(species.list)){
    species1 = tail(full.data$symb[which(full.data$species1 == species.list[i])],1)
    species2 = tail(full.data$symb[which(full.data$species1 == species.list[j])],1)
    symb.matrix[i,j] = abs(species1 - species2)
  }
}
symb.matrix = as.dist(symb.matrix)
mantel.symb = mantel(symb.matrix ~ dist.matrix, nperm = 100)




### PHYLOGENETIC LOGISTIC REGRESSION FOR A BINARY RESPONSE
library(ape)
library(phytools)
library(phylolm)

#test/example code
set.seed(123456)
tre = rtree(50)
x = rTrait(n=1, phy=tre)
X = cbind(rep(1,50),x)
y = rbinTrait(n=1, phy=tre, beta=c(-1,0.5), alpha=1, X=X)
set.seed(1234567)
x = rTrait(n=1, phy=tre)
X = cbind(rep(1,50),x)
x = rbinTrait(n=1, phy=tre, beta=c(-1,0.5), alpha=1, X=X)
dat = data.frame(trait01 = y, predictor = x)
fit = phyloglm(trait01~predictor, phy=tre, data=dat, boot=100, method = "logistic_IG10")
summary(fit)
coef(fit)
vcov(fit)


dat = data.frame(trait01 = y, predictor = sample(y, length(y)))
# my assumption here is to "ignore" (Intercept) significance and just care about predictor

#With my data

tre = read.tree(file = "integraseTree_Jan1622.newick")
plotTree(tre, offset=1)


orf.presence = data.frame(species = tre$tip.label, strain = NA)
for(i in 1:length(orf.presence$species)){
  orf.presence$strain[i] = strsplit(orf.presence$species[i], split = "_")[[1]][2]
  orf.presence$species[i] = strsplit(orf.presence$species[i], split = "_")[[1]][1]
  orf.presence[i,3:(length(common.orf.short[1,]))] = NA
  if(is.na(match(orf.presence$strain[i], common.orf.short$strain))==FALSE){
    orf.presence[i,3:(length(common.orf.short[1,]))] = common.orf.short[match(orf.presence$strain[i], common.orf.short$strain)[1], 3:(length(common.orf.short[1,]))]
    for(j in 3:length(common.orf.short[1,])){
      if(orf.presence[i,j]>0){orf.presence[i,j]=1}
    }
  }
}

missing.tips = tre$tip.label[which(is.na(orf.presence[,3])==TRUE)]
tre = drop.tip(phy = tre, tip = missing.tips)
orf.presence = orf.presence[complete.cases(orf.presence[,3]),]

orf.environment = orf.presence[,1:2]
for(i in 1:length(orf.environment$strain)){
  orf.environment$salt[i] = full.data$salt[which(full.data$strain1==orf.environment$strain[i])[1]]
  orf.environment$habt[i] = full.data$habt[which(full.data$strain1==orf.environment$strain[i])[1]]
  orf.environment$temp[i] = full.data$temp[which(full.data$strain1==orf.environment$strain[i])[1]]
  orf.environment$blum[i] = full.data$blum[which(full.data$strain1==orf.environment$strain[i])[1]]
  orf.environment$symb[i] = full.data$symb[which(full.data$strain1==orf.environment$strain[i])[1]]
}


results.matrix = matrix(data = NA, nrow = length(orf.presence[1,])-2, ncol = 5)
corr.matrix = results.matrix
for(i in 1:length(results.matrix[,1])){
  data.orf = orf.presence[,i+2]
  names(data.orf) = tre$tip.label
  for(j in 1:length(results.matrix[1,])){
    data.env = orf.environment[,j+2]
    names(data.env) = tre$tip.label
    
    dat = data.frame(orf = data.orf, env = data.env)
    fit = suppressWarnings(try(phyloglm(orf~env, phy = tre, data = dat, boot = 100, method = "logistic_IG10"), silent = TRUE))
    
    
    
    if(length(fit)>1){
      corr.matrix[i,j] = fit$logLik
      
      fit = summary(fit)
      results.matrix[i,j] = fit$coefficients[12]
    }
  }
}

write.csv(x = results.matrix, file = "phyloglm_results_3Jan.csv")


# independent variables: environment, phylo.dist
# dependent variables: ORF presence

seed = 1234
data.orf = sample(x, size = length(data.orf), replace = TRUE)
data.env = sample(y, size = length(data.env), replace = TRUE)



### Plot a phylo tree with data labels
setwd("/Users/Jack/Dropbox/Phylo")
library(ape)
library(phytools)
library(phylolm)

#Convert common.orf to shortform data.table
common.orf.short = data.frame(species = NA, strain = unique(common.orf$strain1))
add.columns = as.character(unique(common.orf$orfid))
common.orf.short[,add.columns] = 0
for(i in 1:length(common.orf.short$species)){
  common.orf.short$species[i] = as.character(common.orf$species1[which(common.orf$strain1 == common.orf.short$strain[i])])[1]
  add.data = table(common.orf$orfid[which(common.orf$strain1==common.orf.short$strain[i])])
  common.orf.short[i,names(add.data)] = add.data
  }

jpeg(file = "TestHeat.jpeg", width = 500, height = 2000)
colheat = c("grey90", "wheat", "gold", "orange", "darkred", "purple")
plot(x = c(0, length(unique(common.orf$orfid))), y = c(0, length(common.orf.short$strain)),
     type = "n", xaxt = "n", yaxt = "n", ylab = "Species", xlab = "ORFid", bty="n")
for(i in 1:length(unique(common.orf$orfid))){
  for(j in 1:length(common.orf.short$strain)){
    heat = 1
    if(common.orf.short[j,i+2]==1){heat = 2}
    if(common.orf.short[j,i+2]==2){heat = 3}
    if(common.orf.short[j,i+2]==3){heat = 4}
    if(common.orf.short[j,i+2]==4){heat = 5}
    if(common.orf.short[j,i+2]>=5){heat = 6}
    rect(xleft = i-1, ybottom = j-1, xright = i, ytop = j, col = colheat[heat], border = NA)
  }
}
dev.off()

strain.tree = read.tree(file = "integraseTree_Jan1622.newick")
strain.tree$edge.length = NULL

jpeg(file = "TestTree.jpeg", width = 500, height = 2000)
plotTree(tree = strain.tree, fsize = 1, lwd = 1)
dev.off()

strain.order = rev(strain.tree$tip.label)
orf.store = common.orf.short

common.orf.short = data.frame(species = NA, strain = strain.order)
add.columns = as.character(unique(common.orf$orfid))
common.orf.short[,add.columns] = (-1)
for(i in 1:length(common.orf.short$strain)){
  common.orf.short$species[i] = strsplit(common.orf.short$strain[i], split = "_")[[1]][1]
  common.orf.short$strain[i] = strsplit(common.orf.short$strain[i], split = "_")[[1]][2]
  if(is.na(match(common.orf.short$strain[i], orf.store$strain))==FALSE){
  common.orf.short[i,-c(1:2)] = orf.store[which(orf.store$strain==common.orf.short$strain[i]), -c(1:2)]
  }
}



species.tree = read.tree(file = "Core_Genome_23.newick")
species.tree$edge.length = NULL

jpeg(file = "TestTree23.jpeg", width = 500, height = 2000)
plotTree(tree = species.tree, fsize = 4, lwd = 1)
dev.off()

eco.data = data.frame(species = species.tree$tip.label, salt = (-1), habt = (-1), temp = (-1), blum = (-1), symb = (-1))
for(i in 1:length(eco.data$species)){
  add.data = match(eco.data$species[i], full.data$species1)
  if(is.na(add.data)==FALSE){
    eco.data[i,-c(1)] = full.data[add.data[1],5:9] 
  }
}


col.salt = c("grey80", "pink", "darkgreen")
jpeg(file = "TestBubble.jpeg", width = 500, height = 2000)
plot(x = c(0,1), y = c(0, 2*length(eco.data$species)), type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "S H T B Sy", ylab = "")
for(i in 1:length(eco.data$species)){
  points(x=0.1, y=i*2, pch = 21, bg = col.salt[eco.data$salt[i]+2], col = col.salt[eco.data$salt[i]+2], cex = 2)
  points(x=0.15, y=i*2, pch = 21, bg = col.salt[eco.data$habt[i]+2], col = col.salt[eco.data$habt[i]+2], cex = 2)
  points(x=0.2, y=i*2, pch = 21, bg = col.salt[eco.data$temp[i]+2], col = col.salt[eco.data$temp[i]+2], cex = 2)
  points(x=0.25, y=i*2, pch = 21, bg = col.salt[eco.data$blum[i]+2], col = col.salt[eco.data$blum[i]+2], cex = 2)
  points(x=0.3, y=i*2, pch = 21, bg = col.salt[eco.data$symb[i]+2], col = col.salt[eco.data$symb[i]+2], cex = 2)
}
dev.off()




integrase.tree = read.tree(file = "integraseTree_Jan1622.newick")



unique.species = species.list[match(unique(species.list$strain), species.list$strain),]
unique.species = as.character(interaction(unique.species, sep = "_"))

            # 
            # shared.orf = matrix(0, nrow = length(common.orf), ncol = length(unique.species))
            # shared.orf.species = names(common.orf)
            # for(i in 1:length(shared.orf.species)){
            #   shared.orf.species[i] = paste0(c(as.character(paralog.orf$species1[match(shared.orf.species[i], paralog.orf$orf1)]),
            #                                    as.character(paralog.orf$strain1[match(shared.orf.species[i], paralog.orf$orf1)]),
            #                                    as.character(shared.orf.species[i])),
            #                                  collapse = "_")
            # }
            # 
            # rownames(shared.orf) = shared.orf.species
            # colnames(shared.orf) = unique.species
            # 
            # for(i in 1:length(shared.orf.species)){
            #   orf.match = paralog.orf[which(paralog.orf$orf1 == names(common.orf)[i]),]
            #   orf.match = as.character(interaction(orf.match[,5:6], sep = "_"))
            #   for(j in 1:length(unique.species)){
            #     shared.orf[i,j] = sum(orf.match==unique.species[j])
            #   }
            #   orf.origin = paste0(strsplit(rownames(shared.orf)[i], split = "_")[[1]][1:2], collapse = "_")
            #   orf.origin = match(orf.origin, colnames(shared.orf))
            #   shared.orf[i,orf.origin] = shared.orf[i,orf.origin]+1
            # }
            # 
            # write.csv(shared.orf, file = "shared_orf_matrix.csv")

count.matrix = matrix(0, ncol = length(unique.species), nrow = length(unique.species))
rownames(count.matrix) = unique.species
colnames(count.matrix) = unique.species
full.count.matrix = count.matrix

common.orf.df = subset(paralog.orf, paralog.orf$orf1 %in% as.numeric(names(common.orf)))

for(i in 1:length(unique.species)){
  species1 = unlist(strsplit(unique.species[i], split = "_"))
  store.set = common.orf.df[which(common.orf.df$species1 == species1[1]),]
  store.set = store.set[which(store.set$strain1 == species1[2]),]
  for(j in 1:length(unique.species)){
    species2 = unlist(strsplit(unique.species[j], split = "_"))
    store.set2 = store.set[which(store.set$species2 == species2[1]),]
    store.set2 = store.set2[which(store.set2$strain2 == species2[2]),]
    count.matrix[i,j] = length(store.set2$identity)
  }
  count.matrix[i,i] = length(store.set$identity)
}

write.csv(count.matrix, file = "commonORFcount.csv")



for(i in 1:length(unique.species)){
  species1 = unlist(strsplit(unique.species[i], split = "_"))
  store.set = paralog.orf[which(paralog.orf$species1 == species1[1]),]
  store.set = store.set[which(store.set$strain1 == species1[2]),]
  for(j in 1:length(unique.species)){
    species2 = unlist(strsplit(unique.species[j], split = "_"))
    store.set2 = store.set[which(store.set$species2 == species2[1]),]
    store.set2 = store.set2[which(store.set2$strain2 == species2[2]),]
    full.count.matrix[i,j] = length(store.set2$identity)
  }
  full.count.matrix[i,i] = length(store.set$identity)
}

write.csv(count.matrix, file = "fullORFcount.csv")
