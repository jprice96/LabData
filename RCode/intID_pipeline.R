##### Script for the detection of integrons and their features in bacterial chromosomes

### Packages needed
library(sys)
library(seqinr)
# library(RSelenium)
set.seed(1234)
### Input
# Designate a folder containing one or more .fasta files of bacterial chromosomes for analysis (remember to keep this folder up to date)
file.names = list.files(path = "/home/jack/Desktop/HMMER/StrainSequences/", pattern = "fasta")
file.names = sample(file.names, 20, replace = FALSE)

# Create a blank data frame for locating integrases called "integrase.results"
species.names = unlist(file.names)
for(i in 1:length(species.names)){
  species.names[i] = unlist(strsplit(species.names[i], split="[.]"))[1]
}
integrase.results = data.frame(species = NA, int.start = NA, int.end = NA, sequence = NA)

### Integrase Detection
# Change working directory to the location of the HMMER needed objects
setwd("/home/jack/Desktop/HMMER/StrainSequences/")

# Call on terminal to use HMMER to search for possible integrases
for(i in 1:length(file.names)){
  # Use HMMER function 'nhmmer' to search through chromosome files using the HMM profile 'integrase.hmm' (built independently using HMMER)
  hmm.out = exec_internal(cmd = "nhmmer", args = c("integrase_fresh.hmm", file.names[i]))
  hmm.out = as_text(hmm.out$stdout)
  
  # Identify the top match (if any) and note the start and end location of the proposed integrase in integrase.results
  add.integrase = data.frame(species = NA, int.start = NA, int.end = NA, sequence = NA)
  
  top.match = unlist(strsplit(hmm.out[(match("Scores for complete hits:", hmm.out)+3):(match("Annotation for each hit  (and alignments):", hmm.out)-1)], split = "     "))
  for(x in 1:length(top.match)){
    top.match[x] = strsplit(unlist(top.match[x]), split = " ")
    top.match[[x]] = top.match[[x]][-which(top.match[[x]]=="")]
    if(length(top.match[[x]])==0){top.match[x] = NA}
  }
  top.match = top.match[!is.na(top.match)]
  # For if no matches have been found
  if(top.match[[1]][1] == "[No"){top.match = NA}
  # For other text messages
  suppressWarnings(
  for(x in 1:length(top.match)){
    if(is.na(as.numeric(top.match[[x]][1]))==TRUE){ top.match[x] = NA}
  })
  top.match = top.match[!is.na(top.match)]
  
  if(length(top.match)>0){
    for(x in 1:length(top.match)){
      top.match[[x]] = top.match[[x]][1:6]
    }
    top.match = as.data.frame(t(matrix(unlist(top.match), nrow = 6)))
    colnames(top.match) = c("Evalue", "score", "bias", "sequence", "start", "end")
    add.integrase = data.frame(species = species.names[i],
                               int.start = as.numeric(top.match$start),
                               int.end = as.numeric(top.match$end),
                               sequence = NA)
    
    # Extract actual sequence of proposed integrase from fasta file, storing it in integrase.results
    for(j in 1:length(add.integrase$int.start)){
      if(is.na(add.integrase$int.end[j])==FALSE){
        genome.data = read.fasta(file = file.names[i], as.string = FALSE, strip.desc = TRUE)
        add.integrase$sequence[j] = paste(genome.data[[1]][add.integrase$int.start[j]:add.integrase$int.end[j]], collapse="")
      }
    }
    
    integrase.results = rbind(integrase.results, add.integrase)
  }
}

# Remove any fasta files which failed to produce an integrase from the list of queries for following analyses
integrase.results = integrase.results[-1,]

# Remove any proposed integrase genes below a certain size threshold (likely indicating a false positive)
for(i in 1:length(integrase.results$species)){
  if(abs(integrase.results$int.start[i] - integrase.results$int.end[i]) < 850){
    integrase.results$sequence[i] = NA
  }
}
integrase.results = integrase.results[which(is.na(integrase.results$sequence)!=TRUE),]
species.names = unique(integrase.results$species)

### attC Site Detection
# Save current output for transfer to Linux
setwd("/home/jack/Dropbox/Phylo")
write.table(integrase.results[,1:3], file = "integrase_storage_Jan16.txt")

# The following code needs to be run in a Linux terminal
setwd("/home/jack/Dropbox/Phylo")
integrase.storage = read.table(file = "integrase_storage_Jan16.txt") #Need to change file path to match Linux path

# Determine direction of integrase, will assume that integron runs in opposing direction
integrase.storage$direction = 1
for(i in 1:length(integrase.storage$direction)){
  if(integrase.storage$int.start[i]<integrase.storage$int.end[i]){
    integrase.storage$direction[i] = (-1)
  }
}
# integrase.storage$direction = integrase.storage$direction*(-1)

# HattCI attempt
setwd("/home/jack/Desktop/HattCI") # Need to add genome seq to Dropbox and update path

attc.data = data.frame(species = NA, pos_beg = NA, pos_end = NA, Vscore = NA)

for(i in 1:length(integrase.storage$species)){
  integron.seq = read.fasta(file = paste(c("StrainSequences/", integrase.storage$species[i], ".fasta"), sep = "", collapse = ""))
  slide = integrase.storage$direction[i]
  finished = FALSE
  # Script will continue to extend length of integrase until it cannot find an attC site
  while(finished == FALSE){
    seq.window = integron.seq[[1]][(integrase.storage$int.start[i]+slide):(integrase.storage$int.start[i]+slide+(5000*integrase.storage$direction[i]))]
    write.fasta(sequences = seq.window, names = integrase.storage$species[i], file.out = "hattci_input.fasta")
    
    test.attc = exec_internal(cmd = "HattCI-master/hattci.out", args = c("hattci_input.fasta", "hattci_output.txt"))
    test.attc = readLines("hattci_output.txt")
    if(length(test.attc)==1){
      finished = TRUE
      print(c(integrase.storage$species[i], "ended"))
    }
    if(length(test.attc)>1){
      test.attc = readLines("outHattCI.fasta")
      test.attc = test.attc[seq(1, length(test.attc), 3)]
      
      add.attc = data.frame(species = rep(paste0(integrase.storage$species[i], "-", i), length(test.attc)), pos_beg = NA, pos_end = NA, Vscore = NA)
      for(j in 1:length(test.attc)){
        add.data = unlist(strsplit(test.attc[j], split = "_"))
        add.attc$pos_beg[j] = (as.numeric(add.data[4])*integrase.storage$direction[i] + slide)
        add.attc$pos_end[j] = (as.numeric(unlist(strsplit(add.data[5], split = " "))[1])*integrase.storage$direction[i] + slide)
        add.attc$Vscore[j] = unlist(strsplit(add.data[5], split = " "))[3]
      }
      attc.data = rbind(attc.data, add.attc)
      
      slide = slide + (4000*integrase.storage$direction[i])
      print(c(integrase.storage$species[i], "slide"))
      
      if(file.exists("outHattCI.fasta")==TRUE){
        file.remove("outHattCI.fasta")
      }
    }
  }
}
setwd("/home/jack/Dropbox/Phylo")
write.table(attc.data, file = "attc_storage_hattci_16Jan.txt")


# Use a sliding window approach to have FOLDALIGN detect attachment sites
setwd("/home/jack/Desktop/FOLDALIGN/foldalign.2.5.3") # Need to add genome seq to Dropbox and update path

attc.data = data.frame(species = NA, pos_beg = NA, pos_end = NA, pvalue = NA)

for(i in 8:length(integrase.storage$species)){
  integron.seq = read.fasta(file = paste(c("StrainSequences/", integrase.storage$species[i], ".fasta"), sep = "", collapse = ""))
  slide = integrase.storage$direction[i]
  finished = FALSE
  # Script will continue to extend length of integrase until it cannot find an attC site
  while(finished == FALSE){
    seq.window = integron.seq[[1]][(integrase.storage$int.start[i]+slide):(integrase.storage$int.start[i]+slide+(5000*integrase.storage$direction[i]))]
    write.fasta(sequences = seq.window, names = integrase.storage$species[i], file.out = "foldalign_input.fasta")
    
    # Use FOLDALIGN to compare sequence window with a reference attC
    test.attc = exec_internal(cmd = "bin/foldalign", args = c("-plot_score", "-ID", "'TestSequences'", "reference_attc.fasta", "foldalign_input.fasta"))
    test.attc2 = rawToChar(test.attc$stdout)
    writeLines(test.attc2, "testOut.txt")
    test.attc3 = exec_internal(cmd = "bin/locateHits", args = c("testOut.txt"))
    test.attc4 = rawToChar(test.attc3$stdout)
    writeLines(test.attc4, "testHits.txt")
    
    # Extract table of results and then take only pvalues of 0.01 and below
    header.names = c("seq1_id", "seq1_beg", "seq1_end", "species", "pos_beg", "pos_end", "foldscore", "pvalue", "rank")
    foldalign.output = read.table(file = "testHits.txt")
    colnames(foldalign.output) = header.names
    foldalign.output$pos_beg = foldalign.output$pos_beg + slide
    foldalign.output$pos_end = foldalign.output$pos_end + slide
    
    putative_attc = foldalign.output[which(foldalign.output$pvalue<=0.5),c(4,5,6,8)]
    attc.data = rbind(attc.data, putative_attc)
    
    # Check if complete
    if(length(putative_attc$species>0)){
      slide = slide + 4000
      print(c(integrase.storage$species[i], "slide"))}
    if(length(putative_attc$species)==0){
      finished = TRUE
      print(c(integrase.storage$species[i], "finished"))}
  }
}

# End of Linux section, will need to save attc.data and reload in Mac for rest of code
setwd("/home/jack/Dropbox/Phylo")
write.table(attc.data, file = "attc_storage_16Jan.txt")

# setwd("/home/jack/Dropbox/Phylo")
# attc.storage = read.table(file = "attc_storage_18Aug21.txt")
# attc.storage = rbind(attc.storage, attc.data)
# attc.storage = attc.storage[-which(is.na(attc.storage$species)==TRUE),]

### ORF Detection
# Change working directory to the location of the HMMER needed objects
setwd("/home/jack/Desktop/HMMER/StrainSequences")

# Remember to have integrase.results loaded and named
integrase.results = integrase.storage

# Create a list to accept ORF data called "orf.results"
orf.results = list(NA)

# Using basic pattern recognition, begin combing for ORFs in the putative integron
for(i in 1:length(integrase.results$species)){
  # Define the window size to look at upstream of the integrase start point (assumes integrase is in rvs orientation to integron)
  upstream = 5000 * integrase.results$direction[i]
  
  # Extract the sequence for the selected area to search
  integron.seq = read.fasta(file = paste(c(integrase.results$species[i], ".fasta"), sep = "", collapse = ""))
  integron.seq = integron.seq[[1]][(integrase.results$int.start[i]+upstream):(integrase.results$int.start[i])]
  full.cassette = paste(integron.seq, sep="", collapse="")
  
  # Generate the complement of the selected area
  complement = integron.seq
  # for(j in 1:length(complement)){
  #   if(integron.seq[j]=="a"){complement[j] = "t"}
  #   if(integron.seq[j]=="t"){complement[j] = "a"}
  #   if(integron.seq[j]=="c"){complement[j] = "g"}
  #   if(integron.seq[j]=="g"){complement[j] = "c"}
  # }
  complement = rev(complement)
  full.complement = paste(complement, sep="", collapse="")
  
  # Codons used: START = ATG    STOP = TAG TAA TGA
  # Identify any start codons in both the forward frames...
  cassette = unlist(strsplit(full.cassette, split = "(?<=.)(?=atg)", perl = TRUE))
  seq.shift = c(NA, NA)
  seq.shift[1] = length(unlist(strsplit(cassette[1], split = "")))
  cassette = cassette[-1]
  # ... and the reverse frames
  com.cassette = unlist(strsplit(full.complement, split = "(?<=.)(?=atg)", perl = TRUE))
  seq.shift[2] = length(unlist(strsplit(com.cassette[1], split = "")))
  com.cassette = com.cassette[-1]
  
  # Create a data frame to store ORF information for each fasta file
  orf.data = data.frame(sequence = rep(NA, length(c(cassette, com.cassette))), orf.start = NA, orf.end = NA, frame = NA)
  
  # Define a separate holder for the sequence, so that it can be overwritten when switching strands
  orf.seq = cassette
  full.seq = full.cassette
  strand = 1
  reset = 0
  
  # Go through each of the identified START codons and attempt to find a corresponding STOP codon
  for(j in 1:length(orf.data$sequence)){
    # Allow switching from first strand to second strand
    if(j>length(cassette)){
      orf.seq = com.cassette
      full.seq = full.complement
      strand = (-1)
    }
    if(j==(length(cassette)+1)){reset = j-1}
    
    # Extract the sequence downstream of the START codon
    full.orf = paste(orf.seq[(j-reset):length(orf.seq)], sep = "", collapse = "")
    # Break up the ORF into the codons of the reading frame
    full.orf = unlist(strsplit(full.orf, split = "(?<=.{3})", perl = TRUE))
    #Identify any stop codons
    stop.codon = c(grep("tag", full.orf), grep("taa", full.orf), grep("tga", full.orf))
    
    #If there are no stop codons in the reading frame, give NA
    if(all(is.na(stop.codon))==TRUE){
      orf.data$sequence[j] = "nostop"
    } else{ # Otherwise provide the sequence and position of the ORF, storing it in the orf.data list
      orf.data$sequence[j] = paste(full.orf[1:min(stop.codon, na.rm = TRUE)], sep = "", collapse = "")
      orf.data$orf.start[j] = gregexpr(orf.data$sequence[j], full.seq, perl = TRUE)[[1]][1]
      orf.data$orf.end[j] = orf.data$orf.start[j] + (min(stop.codon, na.rm = TRUE) * 3)
      # Determine reading frame in relation to the integrase (which is in reading frame +1)
      if(strand>0){
        orf.data$frame[j] = (((orf.data$orf.end[j]-integrase.results$int.start[i]) %% 3)+1)*strand
      } else{
        orf.data$frame[j] = (((orf.data$orf.start[j]-integrase.results$int.start[i]) %% 3)+1)*strand
      }
    }
  }
  
  # Removes any non-stop or too small ORFs from results and any ORFs already encapsulated within a larger ORF
  to.remove = list(NA)
  for(j in 1:length(orf.data$sequence)){
    if(length(unlist(strsplit(orf.data$sequence[j], split = "")))<100){
      to.remove = c(unlist(to.remove), j)
    }
    for(k in 1:length(orf.data$sequence)){
      if(k != j){
        if(gregexpr(pattern = orf.data$sequence[k], text = orf.data$sequence[j], perl = TRUE)[[1]][1] != (-1)){
          to.remove = c(unlist(to.remove), k)
        }
      }
    }
  }
  to.remove = to.remove[-(1)]
  if(length(to.remove)>0){ orf.data = orf.data[-(to.remove),] }
  
  # Add the data for this sequence to the complete orf data object
  orf.results[[i]] = orf.data
  print(integrase.results$species[i])
}
saveRDS(orf.results, file = "RawORFdata_16Jan.rds")

for(i in 1:length(orf.results)){
  orf.results[[i]]$species = paste(c(integrase.results$species[i], ".fasta"), sep = "", collapse = "")
}
concat.orf.results = orf.results[[1]]
for(i in 2:length(orf.results)){
  concat.orf.results = rbind(concat.orf.results, orf.results[[i]])
  # print(length(orf.results[[i]][1,]))
}

# Save ORF data
setwd("/home/jack/Dropbox/Phylo")
write.table(concat.orf.results, file = "ORFresults_hattci_16Jan.txt")
orf.results = read.table(file = "ORFresults_03Nov21.txt")

###### Done up to here



### Graphic Representation of Integron
setwd("/Users/Jack/Dropbox/Phylo/Genome_Sequences/Genome_Seq/")
file.names = list.files(pattern = ".fasta")
setwd("/Users/Jack/Dropbox/Phylo/Integron_images")

for(i in 1:length(file.names)){
  # Create save file data for image
  save.name = paste0(unlist(strsplit(file.names[i], split = ".fasta")), "_sampleCM.png")
  png(filename = save.name, height = 400, width = 1000)
  
  # Define the approriate factors for the integron to be drawn
  graph.species = unlist(strsplit(file.names[i], split = ".fasta"))
  graph.integrase = integrase.results[which(integrase.results$species == graph.species), 2:3]
  graph.attc = attc.storage[which(attc.storage$species == graph.species),]
  graph.orf = orf.results[which(orf.results$species == file.names[i]),]
  
  # Define genomic size of area
  direction = 1
  if((graph.integrase[1] - graph.integrase[2])<0){
    direction = (-1)
  }
  upstream = graph.integrase[1]+(5000 * direction)
  graph.length = unlist(c((graph.integrase[2] - (500*direction)), (upstream + (500*direction))))
  
  # Plot graphic
  plot(x=0, y=0, type = "n",  yaxt = "n", xlim = graph.length)
  rect(xleft = graph.length[1], ybottom = 0, xright = graph.length[2], ytop = 0, angle = 45, border = "grey")
  
  # Integrase
  rect(xleft = graph.integrase[1], ybottom = -0.5, xright = graph.integrase[2], ytop = 0.5, angle = 45, border = "orange", col = "orange")
  
  # attC sites
  if(class(graph.attc) != "character"){
    for(j in 1:length(graph.attc$pos_beg)){
      rect(xleft = (graph.integrase[1] + graph.attc$pos_beg[j] * direction), ybottom = -0.4, xright = (graph.integrase[1] + graph.attc$pos_end[j] * direction), ytop = 0.4, angle = 45, border = "green", col = "green")
    }
  }
  
  # ORFs
  for(j in 1:length(graph.orf$orf.start)){
    shift = graph.orf$frame[j]/10
    rect(xleft = (graph.integrase[1] + graph.orf$orf.start[j] * direction), ybottom = (0+shift), xright = (graph.integrase[1] + graph.orf$orf.end[j] * direction), ytop = (0.1+shift), angle = 45, border = "red", col = "blue")
  }
  
  dev.off()
}











### attC site sequence extraction

# using 'attc.results[[i]]' to get the position of attC for genome [i], need to extract sequence from genome [i] based of position

total.attc = 0
rmfile = 0
for(i in 1:length(attc.results)){
  if(length(attc.results[[i]])==3){
  total.attc = total.attc + length(attc.results[[i]][,1])
  } else{rmfile = c(rmfile, i)}
}
file.names.attc = file.names[-rmfile]

extract.attc = data.frame(species = NA, pos_beg = NA, pos_end = NA, sequence = rep(NA, total.attc))
x = 1
for(i in 1:length(attc.results)){
  if((i %in% rmfile) == FALSE){
    integron.seq = read.fasta(file = file.names[i])
    for(j in 1:length(attc.results[[i]][,1])){
      extract.attc$species[x] = file.names[i]
      extract.attc$pos_beg[x] = attc.results[[i]][j,2]
      extract.attc$pos_end[x] = attc.results[[i]][j,3]
      extract.attc$sequence[x] = paste(integron.seq[[1]][extract.attc$pos_beg[x]:extract.attc$pos_end[x]], sep="", collapse="")
      x = x+1
    }
  }
}

setwd("/Users/jackprice/Desktop/HMMER/Integrase/NewInt")
fasta.convert = NA
for(i in 1:length(extract.attc$sequence)){
  fasta.convert = strsplit(extract.attc$sequence[i], split="")
  seq.name = paste(c(i, extract.attc$species[i]), collapse = "_")
  write.fasta(sequences = fasta.convert, names = seq.name, file.out = seq.name)
}


for(i in 1:length(integrase.results$species)){
  extract.int = strsplit(integrase.results$sequence[i], split = "")
  write.fasta(sequences = extract.int, names = integrase.results$species[i], file.out = paste(c("integrase", i), collapse = "_"))
}


