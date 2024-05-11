# selection x partial expression results in
###     death.rate = selection - (selection*(partial.express)^position) + partial.express
###     where selection is the death.rate of wild-type in relation to birth.rate (1)
###     and partial.express is a value between 0 and 1, where 1/partial.express is the cumulative position multiplier (so 0.5 means pos 2 is 2x more that pos 1, pos 3 is 4x more etc)

simulation_loop <- function(pop, start.pop, selection, partial.express, 
                            cassette.shuffle, cassette.loss, integron, integron.length, adaptive.cassette){
  #Set cycle time (in generations)
  cyc.time = 15

  #Create the population
  if(class(pop)=="numeric"){
    population = data.frame(time = 0, pop.size = pop, pop)
    names(population)[3] = paste(integron, sep="", collapse="-")
  }
  
  if(class(pop)=="data.frame"){
    population = data.frame(time = 0, pop.size = sum(pop[1,]), pop)
    names(population)[3:length(population[1,])] = names(pop)
  }
  
  #Starting values/placeholders
  t = 1
  extend = 0
  
  #Run the population until 5 "generations" have passed (average time unit, not literal)
  while(tail(population$time, 1)<cyc.time){
    
    #Create a table providing the current expected number of events per generation
    gen.rates = as.data.frame(matrix(ncol=5, nrow=length(names(population))-2 ))
    names(gen.rates) = c("Genotype", "Birth", "Death", "Shuffle", "Loss")
    gen.rates$Genotype = names(population)[-(1:2)]
    
    #Calculate expected births and deaths first
    for(i in 1:length(gen.rates$Genotype)){
      #Find the location of the adaptive cassette within the integron
      location = NA
      location = match(adaptive.cassette, unlist(strsplit(gen.rates$Genotype[i], split="-")))
      if(is.na(location)==TRUE){location = integron.length*10}
     
      #Assign a rate based off the cassette location and the current number of individuals
      gen.rates$Birth[i] = unlist(population[t, (i+2)] * 1.0 * (1 - (population[t, (i+2)]/start.pop)))
      gen.rates$Death[i] = population[t, (i+2)] * (selection - (selection*(partial.express^(location-1))))
    }
     
    #Calculate the tau period
    tau = NA
    
    alpha = unlist(c(gen.rates$Birth, gen.rates$Death))
    alpha = alpha[alpha!=0]
    alpha.sum = sum(alpha)
    alpha.mean = alpha.sum / length(alpha)
    alpha.var = var(alpha)
    
    tau[1] = (0.01 * alpha.sum)/alpha.mean
    tau[2] = ((0.01^2) * (alpha.sum^2))/alpha.var
    if(sum(tau, na.rm=TRUE)==0){tau = 0}
    tau = min(tau, na.rm=TRUE)

    #Correction for large tau
    if(extend>0){
      tau = 1/alpha.sum
      extend = extend -1
    }
    
    if(tau<(1/alpha.sum)){
      extend = 100
      tau = 1/alpha.sum
    }
    
    if(population$pop.size[t] == 0){
      tau = 5 - population[t, 1]
    }
    
    if(tau == 0){tau = 15 - population[t,1]}
    # if(tau == 0){print("Error: Tau equals 0")}
    
    #Set next time point
    population[t+1, 1] = population[t, 1] + tau
    
    #Introduce stochasticity of birth/death events
    events = gen.rates
    for(i in 1:length(events$Genotype)){
      events$Birth[i] = rpois(1, events$Birth[i] * tau)
      events$Death[i] = rpois(1, events$Death[i] * tau)
    }
    
    #Apply birth/death events
    for(i in 1:length(events$Genotype)){
      population[t+1, i+2] = population[t, i+2] + events$Birth[i] - events$Death[i]
    }
    
    #Population clean-up
    for(i in 3:length(population[1,])){
      if(population[t+1, i]<0){
        population[t+1, i] = 0
      }
    }
    
    #Calculate expected number of integrase events
    for(i in 1:length(gen.rates$Genotype)){ #Shuffle events
      #Find the location of the adaptive cassette within the integron
      location = NA
      location = match(adaptive.cassette, unlist(strsplit(gen.rates$Genotype[i], split="-")))
      #Allow for NA processing
      if(is.na(location)==TRUE){location = integron.length*10}
        
      #Assign a rate based off the cassette location and the current number of individuals
      #First set all as nonadaptive to start
      gen.rates$Shuffle[i] = events$Birth[i] * cassette.shuffle
      #Lastly assign any first position phenotypes
      for(j in 1:length(adaptive.cassette)){
        if(location[j]==1){
          gen.rates$Shuffle[i] = 0
        }
      }
    }
    #Assign integrase events within time tau
    for(i in 1:length(events$Genotype)){
      events$Shuffle[i] = round(rpois(1, gen.rates$Shuffle[i]))
      #Loss events
      gen.rates$Loss[i] = events$Shuffle[i] * cassette.loss
      events$Loss[i] = round(rpois(1, gen.rates$Loss[i]))
    }
    
    
    #Apply integrase events
    for(i in 1:length(events$Genotype)){
      
      #Apply shuffle events
      if(events$Shuffle[i]>0){
        for(j in 1:events$Shuffle[i]){
          #Excise a randomly chosen cassette and reinsert it at the front of the integron
          old.integron = unlist(strsplit(events$Genotype[i], split="-"))
          excised.region = sample(old.integron, 1)
          new.integron = old.integron[-match(excised.region, old.integron)]
          new.integron = c(excised.region, new.integron)
          new.integron = paste(new.integron, sep="", collapse="-")
          
          #Remove previous birth result
          old.genotype = match(events$Genotype[i], names(population))
          population[t+1, old.genotype] = population[t+1, old.genotype] - 1
          
          #Is the new genotype pre-existing in the population?
          check.genotype = match(new.integron, names(population))
          if(is.na(check.genotype) == FALSE){
            population[t+1, check.genotype] = population[t+1, check.genotype] + 1
          }
          
          #If not, create a new genotype column
          if(is.na(check.genotype) == TRUE){
            population[t+1, length(population[1,])+1] = 1
            names(population)[length(population)] = new.integron
          }
        }
      }
      
      #Apply loss events
      if(events$Loss[i]>0){
        for(j in 1:events$Loss[i]){
          #Excise a randomly chosen cassette
          old.integron = unlist(strsplit(events$Genotype[i], split="-"))
          excised.region = sample(old.integron, 1)
          new.integron = old.integron[-match(excised.region, old.integron)]
          new.integron = paste(new.integron, sep="", collapse="-")
          
          #Remove previous birth result
          old.genotype = match(events$Genotype[i], names(population))
          population[t+1, old.genotype] = population[t+1, old.genotype] - 1
          
          #Is the new integron a NULL genotype?
          if(is.na(new.integron) == TRUE){
            new.integron = "0"
          }
          
          #Is the new genotype pre-existing in the population?
          check.genotype = match(new.integron, names(population))
          if(is.na(check.genotype) == FALSE){
            population[t+1, check.genotype] = population[t+1, check.genotype] + 1
          }
          
          #If not, create a new genotype column
          if(is.na(check.genotype) == TRUE){
            population[t+1, length(population[1,])+1] = 1
            names(population)[length(population)] = new.integron
          }
        }
      }
    }
    
    #Population clean-up
    for(i in 3:length(population[1,])){
      if(population[t+1, i]<0){
        population[t+1, i] = 0
      }
    }
    
    #Total up population
    population$pop.size[t+1] = sum(population[t+1, 3:length(population[1,])])
    if(population$pop.size[t+1] == 0){population$time[t+1] == 15}
    if(population$pop.size[t+1] >= (start.pop*2)){population$time[t+1] == 15}
    
    
    #Advance to next time step
    t = t+1
  }
  
  
  return(population)
}


#Function allowing multiple successive simulations on the same population

multi.sim.fun = function(sim.num,
                         start.pop, selection, partial.express,
                         cassette.shuffle, cassette.loss, integron, integron.length, adaptive.cassette){
  
  if(length(selection)==1){selection = rep(selection, sim.num)}
  if(length(partial.express)==1){partial.express = rep(partial.express, sim.num)}
  if(length(cassette.shuffle)==1){cassette.shuffle = rep(cassette.shuffle, sim.num)}
  if(length(cassette.loss)==1){cassette.loss = rep(cassette.loss, sim.num)}
  if(length(adaptive.cassette)==1){adaptive.cassette = rep(adaptive.cassette, sim.num)}
  
  pop = start.pop
  start.pop = sum(start.pop)
  
  sim.result = simulation_loop(pop, start.pop, selection[1], partial.express[1],
                               cassette.shuffle[1], cassette.loss[1], integron, integron.length, adaptive.cassette[1])
  
  if(sim.num>1){
    for(i in 2:sim.num){
      final.pop = sim.result[length(sim.result[,1]),-c(1:2)]
      sim.result = dplyr::bind_rows(sim.result, simulation_loop(final.pop, start.pop, selection[i], partial.express[i],
                                                                cassette.shuffle[i], cassette.loss[i], integron, integron.length, adaptive.cassette[i]))
    }
  }
  
  return(sim.result)
  
}











