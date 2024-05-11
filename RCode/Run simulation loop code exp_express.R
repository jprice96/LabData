setwd("/Users/Jack/Dropbox/RCode/SimFunctions")
source("SimLoop_function_exponent_express.R")
source("Generic_Functions.R")
library(foreach)
library(doParallel)

n.cores = detectCores() - 1
cl = makeCluster(n.cores, type = "PSOCK", outfile = "")
registerDoParallel(cl)

set.seed(12345)
start.time = Sys.time()

#Starting population size, this can be provided as a single value, or as a data.frame with multiple starting integron configs
#start.pop = 10^5
start.pop = data.frame(10^5)
colnames(start.pop) = c("1-2-3-4")

# Integron and evolutionary variables
selection = seq(1.5, 3.5, 0.5)  #Wild-type death rate, relative to birth rate of 1
partial.express = 0.3       #Value from 0 to 1, is the base in the exponential expression gradient (position is exponent)

cassette.excision = 10^((-1)*seq(3, 5, 0.25))
cassette.reinsertion = 0.9

integron.length = 4
adaptive.cassette = "4"

num.replicates = 1000

for(c in 2:length(selection)){
  for(a in 4:length(cassette.excision)){
    for(b in 1:length(cassette.reinsertion)){
      for(g in 1:length(integron.length)){
        for(h in 1: length(adaptive.cassette)){
          data.files = list("Parameters" = list(), "Population" = list())
          integron = seq(1, integron.length[g], 1)
          integron = as.character(integron)
          
          cassette.shuffle = NA
          cassette.loss = NA
          cassette.shuffle = cassette.excision[[a]]*cassette.reinsertion[[b]]
          cassette.loss = cassette.excision[[a]]*(1-cassette.reinsertion[[b]])
          
          data.files$Population = foreach(i=1:num.replicates) %dopar% {
            multi.sim.fun(1, start.pop, selection[c], partial.express, cassette.shuffle, cassette.loss, integron, integron.length[g], adaptive.cassette[h])
          }

          for(i in 1:num.replicates){
            data.files$Parameters[[i]] = data.frame("Excision" = cassette.excision[a], "Reinsertion" = cassette.reinsertion[b], "Selection" = selection[c], "PartialExpression" = partial.express,
                                                    "k" = integron.length[g], "AdaptiveCassette" = paste(adaptive.cassette[[h]], sep="", collapse="-"), "Replicate" = i)
            }
          
          print(data.files$Parameters[length(data.files$Parameters)])
          print(Sys.time())
          
          setwd("/Users/Jack/Dropbox/RCode/SimFunctions/SimData")
          save(data.files, file = paste(c("17-09-21", a, b, c), collapse = "_"))
        }
      }
    }
  }
}

end.time = Sys.time()
time.taken = end.time - start.time
print(time.taken)

save(data.files, file="18-08-21_Test.RData")
