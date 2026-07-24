## Starting with a filtered .vcf genotypes file, this document will cover calling microhaplotype alleles from SNP and indel variants and assessing microhaplotype quality

### Step 1, running [microhaplot](https://ngthomas.github.io/microhaplot/)

##### Microhaplot is already a well documented program, but each experiment ends up being a little different, so here I take you through my steps for the hawksbills.

##### First, make a label file: three columns, tab separated, like this...

    WA_n_214557.sam	WA_n_214557	PacHawk
    VN_n_214689.sam	VN_n_214689	PacHawk
    SL_n_214688.sam	SL_n_214688	PacHawk

##### The first column is the .sam file name for each sample. In the second row you can give the sample a different name as well, but my file names are already descriptive of the population and the sample type (_n for nesters, _c for foraging turtles) so I leave it the same.

##### The last column is a group name. This can be useful for organizing samples, but all of my microhaplotypes are going into the same output file so I just have one group.

##### Microhaplot is an R package, so we run it in the R environment, including the shiny app

    library(microhaplot)
    library(tidyr)
    
    label <- "freb.label"
    sam.path <- "~/Documents/hawksbills/mapped_sams/"
    label.path <- "./freb.label"
    vcf.path <- "./Freb_filt7.recode.vcf"
    app.path <- "~/Shiny/microhaplot"
    
    PacHawk <- prepHaplotFiles(run.label = label, sam.path=sam.path, label.path=label.path, vcf.path=vcf.path, app.path=app.path, out.path="./", n.jobs=3)
    
    save.image("PacHawk.Rdata")
    
    runShinyHaplot("~/Shiny/microhaplot/")

### Step 2, assessing the microhaplotype alleles in the Microhaplot R shiny

##### On the top of the shiny gui there is a panel that allows you to fine tune the minimum read depth, the minimum allelic ratio, and the top n alleles. This last one is particularly important because if your organism is diploid, as most vertebrates, like sea turtles, are, then you should only have two alleles. If there are more than two alleles, this is noise in the data caused by sequencing and genotyping errors.

##### Check the summaries for each individual. If there are two dominant alleles that make up the overwhelming majority of the data for one locus, then it is fine to set the top n alleles to 2, and let those superfluous loci be filtered out by the program.

##### If, on the other hand, there are more than two compelling alleles for individuals at a particular locus, then you should discard that microhaplotype.

### Step 3, refining allelic ratios

##### There are a lot of filtering options in Microhaplot that you can take advantage of, but much of your time should be spent refining allelic ratios for each locus. This is the graph that plots the depth of the first allele against the depth of the second allele.

##### I recommend that you start with a small minimum allelic ratio, and minimum total read depth that makes sense for that percentage. 

##### For example, if you set your minimum allelic ratio to 0.25, you should have a minimum total read depth of something that divides into quarters evenly, like 8. That way, you know that at least 2 sequences have the low-coverage allele. Don't, then, have a read depth of ten, because you can't have 2.5 copies of the allele.

##### For a normal locus, with some homozygotes and some heterozygotes, you should see two different clusters. One cluster has a long vertical orientation close to the y-axis. The other has a long diagonal orientation, somewhere in the middle of the graph.

##### The vertical cluster are the homozygotes. In a perfect world, these homozygous genotypes would line up perfectly along the vertical axis having no read depth for the other allele. However, in my experience it is common for even strongly homozygous allele calls to have minor traces of the other allele. Throwing out these genotypes as errors would be overly conservative, in my opinion, and result in too much loss of data, and signal. Fortunately, you can set a locus-specific allelic ratio and still have a homozygous genotype. For this project, that ratio was never larger than 0.09

##### The diagonal cluster of heterozygotes can naturally show quite a bit of scatter as well. and may need to be refined. Again, I think it is somewhat reasonable for the lower-coverage allele to have a quarter, or a fifth of the depth, especially at lower total depths. For this project, I never accepted an allelic ratio for a heterozygote lower than 0.2.

##### And don't be a snob about low read depths. Many of these low-coverage genotype calls are legitimate, and you can prove it below, when we look at genotyping consistency across replicates. 

##### In a good microhaplotype locus, there should be a clear separation between the homozygotes and the heterozygotes when the relative depths of each allele are plotted against each other. Refining the allelic ratio is all about determining the boundaries.

##### While refining allelic ratios, determining an appropriate total read depth, and so forth, you are inevitably going to compromise the total callable microhaplotype alleles for each individual. If the proportion of total callable alleles for individuals or loci falls below about 0.75, you should consider removing that individual or locus.

##### When you're finished, save the observed "filtered" data and load it into R.

    mhaps <- read.csv("PacHawk_mhaps_raw.csv", header=T, stringsAsFactors=F)
    
    dim(mhaps)
    [1] 25902    16
    
    head(mhaps)
      X   group     indiv.ID        locus haplo depth allele.balance rank ave.entropy n.alleles min.rd min.ar max.ar.hm min.ar.hz status comment
    1 1 PacHawk Po_n_198370B Turtle_11355 GGCGT  5175      1.0000000    1           0         2      3    0.5       0.5      0.50 Reject      NA
    2 2 PacHawk    Hh_c_8593 Turtle_22141  AGCG  3667      1.0000000    1           0         2      3    0.5       0.5      0.50 Reject      NA
    3 3 PacHawk  Po_n_207728 Turtle_55368 CTGGT  3487      1.0000000    1           0         2      3    0.5       0.5      0.50 Reject      NA
    4 4 PacHawk    Hh_c_8593 Turtle_22141  CGCG  3473      0.9470957    2           0         2      3    0.5       0.5      0.50 Reject      NA
    5 5 PacHawk  Po_n_207723  Turtle_6296     T  3418      1.0000000    1           0         2      3    0.5       0.1      0.38 Accept      NA
    6 6 PacHawk  Ng_n_106750 Turtle_11355 GGCGT  3020      1.0000000    1           0         2      3    0.5       0.5      0.50 Reject      NA

##### Get rid of the rejected calls and unnecessary columns

    mhaps <- mhaps[-which(mhaps$status=="Reject"),-c(1:2,9,15:16)]
    
    dim(mhaps)
    [1] 21515    11
    
    head(mhaps)
          indiv.ID        locus haplo depth allele.balance rank n.alleles min.rd min.ar max.ar.hm min.ar.hz
    5  Po_n_207723  Turtle_6296     T  3418              1    1         2      3    0.5      0.10      0.38
    7  Ap_n_207725  Turtle_6296     G  2997              1    1         2      3    0.5      0.10      0.38
    8    Hh_c_8593 Turtle_40995    GC  2962              1    1         2      3    0.5      0.12      0.50
    9  Ka_n_207592  Turtle_6296     G  2874              1    1         2      3    0.5      0.10      0.38
    10 Po_n_298368  Turtle_6296     T  2870              1    1         2      3    0.5      0.10      0.38
    11 AP_n_190812 Turtle_11050     A  2836              1    1         2      3    0.5      0.12      0.50

##### Rank isn't too important, but it will be 2 or less as long as you only accepted the top 2 mhaps A homozygote has a rank of 2 if you had to drop one low-frequency allele to make it a homozygote. Heterozygotes are always rank 1 when only the top two alleles are considered.

### Step 4, estimating microhaplotype genotyping consistency

##### As with the microsatellites, the first thing to do is make a list of the replicate sample names. These are normally indicated with a B, in my data, and we are going to remove the B from the names in our list.

    reps <- mhaps$indiv.ID
    repsB <- reps[grep("B$",reps, ignore.case=T)]
    
    reps2 <- unname(sapply(repsB, function(x) substr(x,1,nchar(x) - 1)))
    
    reps3 <- reps2[which(reps2 %in% reps)]       ### make sure we are only using reps that have originals
    
    length(unique(reps3))
    [1] 19

##### Then we can use this function, which takes the data and the replicate sample names as inputs.

    mhap_acc <- function(file, replicates) {
       mismatches=NULL
       counter=length(unique(file$locus))
       for(i in unique(file$locus)) {
          new <- file[which(file$locus == i),]
          for(j in unique(replicates)) {
    	     samp <- new[grep(j,new$indiv.ID),]
    		 samp_tab <- table(samp$indiv.ID, samp$haplo)
    		 tab_dim <- dim(samp_tab)
    		 if(tab_dim[1] * tab_dim[2] < 4) {
    		    next
    		 }
    		 else if(samp_tab[1,2] != samp_tab[2,2] | tab_dim[1] * tab_dim[2] > 4) {                       
    		    mismatches <- c(mismatches, paste(i,j, sep="-"))
    		 }
    	 }
       print(counter)
       counter <- counter - 1
       }
       return(mismatches)
    }

##### Now run the function
    the_acc <- mhap_acc(mhaps, reps3)

##### The results are hyphenated locus-turtle.

    head(the_acc)
    [1] "Turtle_52893-GM_c_215047" "Turtle_52893-GM_c_215047"
    [3] "Turtle_52893-GM_c_215047" "Turtle_52893-GM_c_215047"

##### These are counts of mismatching alleles across replicates, so now we separate them into loci and turtles for examination.

    the_acc_loc <- unname(unlist(sapply(the_acc, function(x) strsplit(x, "-")[[1]][1])))
    the_acc_turt <- unname(unlist(sapply(the_acc, function(x) strsplit(x, "-")[[1]][2])))
    
    hist(table(the_acc_loc), breaks=50)
    hist(table(the_acc_turt), breaks=50)

##### Looking at these histograms, we see that there are a maximum of six mismatches, for both turtles and loci. And only 28 loci had any mismatching at all, out of 99 microhaplotype markers.

    mean(table(the_acc_loc))
    [1] 1.535714

##### The mean number of mismatches per locus was approximately 1.5 meaning that on average loci had a genotyping consistency of over 98%

### Step 5, exporting the data

##### Microhaplot give you the data in a long format, but you will need a wide format as input for most downstream applications.  But first we need to remove the replicate samples from the data

    reps2remove <- paste0(unique(reps3), "B")
    
    new_mhaps <- mhaps[-which(mhaps$indiv.ID %in% reps2remove),]


##### Now, here is a function that will give us the wide format

    library(tidyr)
    library(dplyr)
    
    make_wide_form <- function(df) {
    require(dplyr)
    require(tidyr)
    df2 <- arrange(df[,1:8], indiv.ID, locus)                       # arrange sample ID’s and loci in order
    ids <- unique(df2$indiv.ID)                                            # a list of your sample ID’s
    WF=NULL                                                                # an object to hold the formatted data
    for (i in 1:length(ids)) { id <- df2[which(df2$indiv.ID==ids[i]),]     # loop through turtles
    dups <- duplicated(id$locus)                                           # identify heterozygous loci (duplicated locus names)
    dups[which(dups==TRUE)-1] <- TRUE                                      # tag the first allele as TRUE as well as the dup
    n.times <- as.numeric(dups)                                            # count number of duplicates per locus
    n.times[which(n.times==0)] <- 2                                        # loci w/ 0 dups are homozygous, still need 2 alleles
    df3 <- id[rep(seq_len(nrow(id)), n.times),]                            # dataframe w/ right number of entries
    df3$locus <- make.unique(as.character(df3$locus))                      # make unique 2nd allele name
    gt <- spread(df3[,1:3], locus, haplo)                                  # put everything in the right shape
    WF <- bind_rows(WF, gt)                                                # now bind it to the WF object
    }                                                
    return(WF)                                                             # let the first loop finish before you return()
    }
    
    wide_mhaps <- make_wide_form(new_mhaps)
    
    
