## Starting with raw fastq files, this document will cover separating microhaplotype markers from microsatellite repeat sequences in a library, and filtering low-quality loci and samples

##### If you decide to work with the raw fastq files on NCBI, you'll find file names that look like this...

HH_c_94476.fastq.gz
PL_n_28824.fastq.gz
PL_n_28824B.fastq.gz

##### Each sample has an "n" or a "c" which stands for nester or collection. These are respectively the baseline population genetic samples and the turtles of unknown origin that are to be assigned to a baseline unit.

##### The first two letters are codes for the sampling locations

HH = Hawaii strandings + in-water captures
GM = Guam (in-water caputures)
NM = Northern Marianas (in-water captures)
FJ = Fiji (nesters)
NT = N. Territoty Australia (nesters)
WA = Western Australia (nesters)
PG = Papua New guinea (nesters)
SL = Solomon Islands (nesters)
VN = Vanuatu (nesters)
AS = American Samoa (nesters)
NG = Nicaragua (nesters)
EC = Ecuador (nesters)
PL = Palau (nesters)
PO = Pohue (nesters)
KA = Kamehame (nesters)
AP = Apua (nesters)
MA = Maui (nestlings)
MO = Molokia (nestlings)

### Step 1 is to trim the sequences using the program [fastp](https://github.com/opengene/fastp)

##### fastp will automatically trim any standard adapter sequences and noisy bp at the ends of sequences. The default base pair quality threshold is a phred Q score of 15, which I think is reasonable for routine trimming. But you need to specify the a minimum sequence length that works for your data. Any sequences that are trimmed shorter than this will be filtered.

##### For the hawksbill data set, a minimum read length of 90 (out of a possible 150 bp single-end sequence) kept most of the data

##### Some people might find 90 bp short, but in my view this is fine as long as you realize that this might affect the depth of some SNP variants located at the extreme ends of sequences, or cause some allele dropout of microsatellite loci

##### My personal preference is to be inclusive during this early data processing phase, and figure out which markers are not performing well in later filtering steps.

#####  You can put all file names in a .txt batch file and run fastp like this in python...

    >>> import os
    
    ### make batch file of fastq filenames
    >>> os.system("ls | grep fastq > batch.txt") 
    
    >>> with open("batch.txt", "r") as forward:
    ...   for fw in forward:
    ...      os.system("fastp -i %s -o %s.T.fastq --length_required 90" %(fw.rstrip("\n"), fw.rstrip(".fastq\n")))

##### Or, if you are comfortable using bash script, this one-liner will do the same thing...

    for j in `ls`; do fastp -i ${j} -o ${j:0:-8}T.fastq --length_required 90; done

##### use -8 for gzipped fastq files (fastq.gz), -6 for uncompressed.

##### In my system I output trimmed fastq files as T.fastq

##### One warning: I have had bad experiences with certain versions of fastp. The version I used for the paper was fastp 0.20.1.

### Step 2 is to map all fastq reads to a reference

##### This reference is a curated fasta file of just the microhaplotype loci used in this study. This fasta file can be found in the same Github director as this document and is called, "new_Eimb.fasta"

##### You can use pretty much any read aligner to map the sequences, like bwa-mem or bowtie2, but here I used minimap2

##### After you've formated the reference for your aligner of choice, you can do all files in a batch, in python, as before...

    import os
    
    os.system("ls | grep 'T.fastq' > batch.txt")
    
    ### save the reference directory as a variable
    
    ref_dir="~/Documents/hawksbills/Eimb_ref.mmi"
    
    with open("batch.txt","r") as files:
       for f in files:
          RG="@RG\tID:%s\tSM:%s %(f,f)
          os.system("minimap2 -t 4 -ax sr ref_dir -R %s %s > %s.sam" %(RG, f.rstrip("\n"), f.rstrip(".mapped.fastq\n")))
     
##### Or, you can use a one-line bash script

    for i in `ls | grep fastq`; do minimap2 -t 4 -ax sr ~/Documents/hawksbills/Eimb_ref.mmi -R "@RG\tID:${i:0:-13}\tSM:${i:0:-13}" ${i} > ${i:0:-13}.sam; done

##### the "Eimb" is the name of the Bowtie formated reference

##### Imporantly! We want to keep all unmapped reads in our .sam outputs. Some aligners will exclude them, or allow you to exclude them, but minimap2 always includes them so if you choose another aligner make sure to keep the unmapped sequences.

### Step 4 is to make new .sam files for the unmapped reads

##### For this, I put the .sam files in a new directory, with nothing else in it, and run this one-liner that uses the [samtools](https://www.htslib.org/) program

    for f in `ls`; do samtools view -h -f 4 ${f} > ${f:0:-4}.unmap.sam; done 

##### The above line copies unmapped samples to new files but doesn't remove them from the original files!!!
##### Do this AFTERWARDS, to remove unmapped reads

    for f in `ls`; do samtools view -h -F 4 ${f} > ${f:0:-4}.mapped.sam; done 

##### Now we have mapped.sam files with our microhaplotype loci, and unmap.sam files that contain the microsatellite repeat sequences

### Step 4, we will process the microsatellite loci and genotype them in the program [MegaSat](https://github.com/beiko-lab/MEGASAT) 

##### MegaSat requires .fastq files as input, not .sam files, so first we need to convert .sam back into fastq

##### First, filter out reads that don't have msat primers in them. Make a .txt file with the primer sequences in them (you can get them from the [MegaSat primer input file](https://github.com/jh041/Pacific-hawksbill-turtle-supplement/blob/main/MegaSat%20primer%20file) in this same GitHub repository) and use fgrep to search and retain all sequences that contain them. 

    for f in `ls | grep 'unmap'`; do cat ${f} | fgrep -f ./msat_primers.txt > ${f:0:-10}.sats.sam; done

##### At this point you might be wondering if these primer sequences will have been trimmed in fastp. Some of them, sure, but mostly the reverse primers at the end of the sequence, or the forward primers if the sequence is a reverse compliment. But if you have both forward and reverse sequences in your "msat_primers.txt" file then all you need is one of them to grab the sequence.

##### If you don't have either, well, that read was never going to work in MegaSat anyway, so there's no harm in throwing it out.

##### Once we have enriched our .sam file for the sequences we want, we can convert to fastq

##### Python...

    import os
    
    os.system("ls | grep 'sats' > batch_s.txt")
        
    with open("batch_s.txt", "r") as sats:
       for sat in sats:
          os.system("samtools fastq -0 %s.fastq %s" %(sat.rstrip(".sats.sam\n"), sat.rstrip("\n")))

##### Bash...

    for f in `ls`; do samtools fastq -0 ${f::-4}.fastq ${f}; done

##### Again, if it isn't already clear, feeding files into the program using `ls` assumes that the only files in your current working directory are the files you want to process. So I create new directories for each stage of the process. That's just my method.

#####  Now, we also need to remove the underscores from the file names, otherwise Megasat will truncate them. Here's a one-liner that will do that...

    for file in *; do mv "${file}" "${file/_/}"; done
    for file in *; do mv "${file}" "${file/_/}"; done

##### The asterisk here does the same thing as `ls` in the previous examples. Just another way to take everything from a directory and do something to it.

##### I really like this string-editing method, but that's not why I have it twice. The command only removes the first underscore, so we have to run it two times.

##### Now we run MegaSat

    perl Megasat.pl primer-input.txt 2 5 3 ./ ./megasatout

##### This is a perl script, and the first command-line argument is the [primer input file](https://github.com/jh041/Pacific-hawksbill-turtle-supplement/blob/main/MegaSat%20primer%20file). 

##### The second argument is the number of mismatching base pairs allowed between the primer and flanking region sequence. Two is recommended.

##### The third argument is the minimum depth allowance for each locus.

##### The fourth argument is the number of processes (parallelization).

##### Then, is the directory containing the input .fastq files and the output directory. If this directory doesn't already exist, the program will make it for you.

##### To run the perl script, will need the cpan Parallel::ForkManager perl module

### Step 5, assess and filter the microsatellite genotypes

##### There are a number of outputs from MegaSat, but the one we are most interested in is the Genotype.txt file. To work with this data, import it into R.

    raw_msats <- read.table("Genotype.txt", stringsAsFactors=F, header=T)
            
    dim(raw_msats)
    [1] 175  59 
    
    raw_msats[1:5, 1:5]
        
      Sample_idx1_idx2 GT_Cc1G02 GT_Cc1G02.b GT_Cc2F11 GT_Cc2F11.b    
      1  PLn28824B.fastq 0 0 107 139    
      2  Pon198371.fastq 134 138 103 127    
      3 Pon198365B.fastq 114 122  92 127    
      4  MOn200702.fastq 114 114 299 299
      5 Pon198369B.fastq 134 138 103 127 

##### Here we can see, in this subset of the data, that the first column is the sample name, and the genotypes are in two columns per locus.

##### Next, I'm going to remove the "fastq" part of the sample names, to clean things up a bit and count the number of genotyping errors in each column.

##### In MegaSat, there are three types of genotyping errors: X, 0, and Unscored
##### X means that the locus wasn't found at all in the sample
##### 0 means that the read depth was too low
##### Unscored means that the locus wasn't scorable due to having more than two alleles

##### All of them are simply unsuitable for further analysis and we can count them altogether as bad data

    sample_names <- sapply(raw_msats$Sample_idx1_idx2, function(x) substr(x, 1, nchar(x) - 6))
    raw_msats$Sample_idx1_idx2 <- sample_names
    
    msat_baddies <- apply(raw_msats, 2, function(x) length(which(x %in% c("X", "0", "Unscored"))))
    
    length(msat_baddies)
    [1] 59

##### There are 59 items in "msat_baddies", one for each column of the data frame. 

##### These errors are typically the same for both alleles in a locus, so you only need one. Take the odd numbered columns.

    msat_baddies <- msat_baddies[seq(1, length(msat_baddies), by = 2)]

##### Let's also do this for the turtles, to see if there are any samples that performed poorly in genotyping

    turt_baddies <- apply(raw_msats, 1, function(x) length(which(x %in% c("X", "0", "Unscored"))))

##### At this point, a good thing to do is plot a histogram of the errors, to see what a normal error rate for your loci.

    hist(msat_baddies)
    hist(turt_baddies)

##### Often, there is usually at least some samples that drag down the entire data set, that for one reason or another, just have low quality DNA that doesn't amplify. Find these bad individuals and remove them from the data first.

##### Find the sample names of these poor performers

    turt_baddies[which(turt_baddies > 20)]

##### Then put all the sample names from this list in a character vector

    Worst_turts <- c("Bad_indiv1","Bad_indiv2","Bad_indiv3","Bad_indiv4","Bad_indiv5")

##### And remove them from the raw genotype calls

    better_msats <- raw_msats[-which(row.names(raw_msats) %in% worst_turts),]

##### You can continue on this way, removing bad samples and loci (rows and columns), as best suits your particular data. After gettting rid of the worst offenders, I recommend alternating between removing samples and loci, incrementally improving the data quality, little by little.

##### However, depending on your needs, you may want to prioritize either samples or loci.

##### If you have a lot of samples to spare, and want to maximize loci for the sake of analysis power later on, then be less strict with loci.

##### If you have enough loci, and need to retain as many samples as possible, keep more turtles.

##### Keep in mind, here, that we will be combining multiple marker types into the same final data set, so when anticipting the total number of loci we don't need to limit ourselves to what we see in the microsats.

##### You have to be careful, because these initial filtering considerations can influence your ultimate results in meaningful ways. It can be a good idea to experiment with several different  data filtering schemes to see how sensitive your results are to these decisions. But in my experience, if the signal in your data is strong (as it was with the hawksbill turtles) then you're really just trying to remove the surrounding noise, not surgically removing a signal from a field of mostly noise.

##### And in order to reduce noise, in addition to looking at errors that MegaSat was able to detect, you need to find the errors that the software didn't detect.

##### It is recommended that you run 10-20% of your samples twice as technical replicates, so that you can assess the genotyping consistency of each locus. Some polymorphisms may not be accurate, if they do not give you the same genotype calls across replicates.

##### All of my replicates are indicated with a "B" in the sample name, so I first find all of these samples, then the original samples with the same name but without the B.

    reps <- row.names(better_msats)
    reps2 <- reps[grep("B", reps, ignore.case=T)]
    length(reps2)
    [1] 19
    
    reps3 <- unname(unlist(sapply(reps2, function(x) substr(x, 1, nchar(x) - 1))))

##### The reps3 object is a character vector of our replicate sample names with the "B" removed. So that it will match the names of the original samples.

##### Then, if we feed our data frame and the reps3 object into the following function, we'll get the genotyping consistency of each locus.

    error_rates <- function(q,p) {
    stuff=NULL
    for(i in 2:length(names(q))){
       error <- 0
       for(j in p){
          the <- q[grep(j, q$Sample_idx1_idx2),]
          the1 <- the[1,names(q)[i]]
          the2 <- the[2,names(q)[i]]
    	  if(is.na(the1) | is.na(the2)) {next}
          if(the1 != the2 & 
    	  the1 != "X" &
    	  the1 != "0" &
    	  the1 != "Unscored" &
    	  the2 != "X" &
    	  the2 != "0" &
    	  the2 != "Unscored"){
    	     error <- (error+1)
    	  }
       }
       stuff <- rbind.data.frame(stuff,cbind.data.frame(names(q)[i], error/length(rownames(q))/2))
    }
    return(stuff)
    }
     

 ##### This function ignores missing data, so that only positive genotype calls are compared across replicates.

    mismatches <- error_rates(better_msats, reps3)

##### The output is two columns: the locus name, and the number of mismatching allele calls

    head(mismatches)
      names(q)[i] error/length(rownames(q))/2
    1 GT_Cc1G02.b                           0
    2   GT_Cc2F11                           0
    3 GT_Cc2F11.b                           0
    4   GT_Cc5C08                           0
    5 GT_Cc5C08.b                           0
    6   GT_Cc5F01                           0

##### Amazingly, none of the loci had any mismatches at all. 100% genotyping consistency!

##### That won't be the case with the microhaplotype loci, which we'll look at later.

### Step 6, calling SNPs from the .mapped.sam files

##### Next, we are going to convert the mapped .sam files to .bam format in samtools and index them. we're also going to filter mapped reads by mapping quality.

    for  f  in *.sam;  do samtools view -@ 4 -b -q 20  "$f"  >  "${f%.sam}.bam" done

##### The -q option in this command removes reads below the specified mapping quality.
##### the -@ is for multi-threading
##### the ${f%.sam} syntax removes the .sam suffix

    for  f  in *.mapped.bam;  do samtools sort -@ 8 -o "${f%.mapped.bam}.sorted.bam"  "$f" done

##### Now we're ready to call SNPs

##### There are a number of good SNP callers you can use. In the past, I have used multiple SNP callers as well. There are some differences in the results, so using different callers and taking the common SNPs among their outputs can be a smart way to get good genotypes. However, I've also noticed that if the signal is strong in your data, then the downstream results won't be sensitive to which caller you use, even if some of the variants called are different.

##### So here I'm only going to use the SNP and indel caller Freebayes. 

    ls | grep bam$ > the_bams.txt
    
    freebayes -f ~/Documents/hawksbills/hawksbill_reference/new_Eimb.fasta -j -k -u -4 -m 20 -q 20 -F 0.2 --min-coverage 5 --bam-list the_bams.txt > mhap_Freb_PacEi.vcf

##### Here "the_bams.txt" is our batch file with all the .bam file names.
##### Notice that the reference for SNP calling is the fasta file I mentioned above. Freebayes will generate its own fasta index, or you can create one for it like this...

    samtools faidx new_Eimb.fasta

##### The way the Freebayes command is constructed requires a minimum variant coverage of 5 to call alleles, and a minimum allelic ratio of 2. That is, at least one of the minimum 5 sequences of depth needs to be an alternate allele to the reference. If the depth coverage is 100, then at least 20 sequences need the alternate allele. Anything less than that is treated as a sequencing error.

##### The output is in variant call format (.vcf) and we can view and manipulate it in the program vcftools.

    vcftools --vcf mhap_Freb_PacEi.vcf
    After filtering, kept 175 out of 175 Individuals
    After filtering, kept 1007 out of a possible 1007 Sites

##### If we just run the output file through .vcf tools, it will count the number of individuals in the file and the number of raw variants called.

##### The first thing I want to do with the .vcf is remove any sites with over 50% missing data, getting rid of low quality SNPs and indels

    vcftools --vcf mhap_Freb_PacEi.vcf --max-missing 0.5 --recode --out Freb_filt1
    After filtering, kept 863 out of a possible 1007 Sites 

##### After that, I want vcftools to collect some data for me about missing data in both loci and individuals

    vcftools --vcf Freb_filt1.recode.vcf --missing-indv --out mhap_freb
    vcftools --vcf Freb_filt1.recode.vcf --missing-site --out mhap_freb

##### I can then load these data into R to explore them.

    ind <- read.table("mhap_freb.imiss", header=T)
    loc <- read.table("mhap_freb.lmiss", header=T)

##### A quick look at the loci shows that we go even stricter in our missing data filtering without costing us too many loci. Instead of 50 percent, I want to remove any locus with 15% or more missing data, so let's redo that now. 

    vcftools --vcf mhap_Freb_PacEi.vcf --max-missing 0.86 --recode --out Freb_filt1
    After filtering, kept 767 out of a possible 1007 Sites

##### That's pretty acceptable. Now,let's see how many individuals have high amounts of missing data and should be removed for poor quality.

    ind[which(ind$F_MISS > 0.50),]
                INDV N_DATA N_GENOTYPES_FILTERED N_MISS   F_MISS
    2    VN_n_214695    767                    0    562 0.732725
    23    PL_n_28824    767                    0    675 0.880052
    70  MA_n_201551B    767                    0    446 0.581486
    115  AS_n_133956    767                    0    476 0.620600
    123  GM_c_174031    767                    0    559 0.728814
    133  FJ_n_174607    767                    0    687 0.895698
    161  AS_n_102737    767                    0    437 0.569752
    164  AS_n_102736    767                    0    577 0.752282
    165  GM_c_174037    767                    0    451 0.588005

##### For whatever reason, these nine turtles have more missing data than I want to deal with, so removing them from the data set will greatly improve our metrics. Nine is not a lot compared to our overall sample numbers, so we can simply cut them out of the .vcf file in vcftools.

    vcftools --vcf Freb_filt1.recode.vcf --remove bad_turts.txt --recode --out Freb_filt2
    After filtering, kept 166 out of 175 Individuals
    After filtering, kept 767 out of a possible 767 Sites

##### If it isn't already obvious, "bad_turts.txt" is a text file containing the names of the individuals to be removed. You can also remove single individuals from the .vcf file this way...

    vcftools --vcf Freb_filt3.recode.vcf --remove-indv GM_c_174037 --remove-indv PO_n_190822 --recode --out Freb_filt4

##### Ultimately, if you alternate between removing bad loci and bad individuals, you'll end up with a data set that is mostly signal and very little noise. But study design considerations also factor into these decisions.

##### For example, because we are working with an endangered species here, each sample is very precious, and I don't want to discard any more turtles than I have to. For that reason, Loci were not allowed to have more than 10% missing data across all samples, and turtles were not allowed to have more than 25% missing data across all samples.

##### For some studies, it might also be a good idea to filter for minor allele frequency, especially if you plan on estimating effective population size. But here I have so few samples from a lot of populations with small sample sizes, if I did that now, with the entire sample set, I'd lose signal.

##### It can also be a good idea to filter for HWE, but hawksbills have population-specific patterns of inbreeding that are too nuanced for a blunt HWE filter (see Horne et al. 2023: https://doi.org/10.1098/rsos.221547).

##### If you want to be really conservative you can exclude indel polymorphisms too, since these often have higher error rates than SNPs, but I've decided to include indels in this study as long as they don't look like this...

    Turtle_204      179     .       AGGGCA  AGGGGCA 17304.8 .       .       GT:DP:AD:RO:QR:AO:QA:GL 0/0:7:7,0:7:253:0:0:0,-2.10721,-22.2994 0/1:19:6,13:6:226:13:309:-18.7868,0,-13.8901    >
    Turtle_666      95      .       CTGTGTGTGTGTGTGA        CTGTGTGTGTGTGTGT,CTGTGTGTGTGTGAGA       253651  .       .       GT:DP:AD:RO:QR:AO:QA:GL 0/1:189:95,84,0:95:3539:84,0:3134,0:-218>
    Turtle_1627     242     .       CTTTTGC CTTTTTGC        8.21112e-09     .       .       GT:DP:AD:RO:QR:AO:QA:GL 0/0:21:20,1:20:740:1:38:0,-2.73405,-61.0777     0/0:12:12,0:12:444:0:0:0>
    Turtle_1688     153     .       ATGTGTGTGTGTGTGGGGGGGGATA       ATGTGTGTGTGTGGGGGGGGATA,ATGTGTGTGTGTGGGGGGGGGGATA,ATGTGTGTGTGTGTGTGGGGGGATA,ATGTGTGTGTGTGTGGGGGGGGGGATA 107359  .       >
    Turtle_5027     121     .       CTTT    CTTC    542.069 .       .       GT:DP:AD:RO:QR:AO:QA:GL 0/0:132:130,0:130:4747:0:0:0,-39.1339,-412.484  0/0:1:1,0:1:33:0:0:0,-0.30103,-2.25881  >
    Turtle_5578     55      .       TAAAAAAAAGTC    TAAAAAAGAAGTC   8752.37 .       .       GT:DP:AD:RO:QR:AO:QA:GL 0/0:24:24,0:24:887:0:0:0,-7.22472,-77.4933      0/0:19:18,0:18:671:0:0:0>
    Turtle_5973     45      .       CAG     CAAG    39.3142 .       .       GT:DP:AD:RO:QR:AO:QA:GL 0/0:7:7,0:7:265:0:0:0,-2.10721,-23.3212 .:.:.:.:.:.:.:. 0/0:4:4,0:4:149:0:0:0,-1.20412,->
    Turtle_6489     67      .       TAAAAAAATATTG   TAAAAAAGTATTG   4.31561 .       .       GT:DP:AD:RO:QR:AO:QA:GL 0/0:1:1,0:1:36:0:0:0,-0.30103,-3.5027   0/0:3:3,0:3:112:0:0:0,-0.90309,->
    Turtle_7153     191     .       CCATGGAGGC      CCGTGGAGGC,CC,CTGTGGAGGC        108783  .       .       GT:DP:AD:RO:QR:AO:QA:GL 2/2:71:0,0,71,0:0:0:0,71,0:0,2436,0:-147.985,-147.985,-1>
    Turtle_9622     168     .       GTTTTTTTATTTTTATTTATTTATTTGC    GTTTTTTAATTTTTATTTATTTATTTGC    1.21225e+06     .       .       GT:DP:AD:RO:QR:AO:QA:GL 0/1:759:323,310:323:11851:310:11>
    Turtle_11355    167     .       CCTCCCCCCTCCCCCAAG      CCTCCCCCCCTCCCCCAAG     1006.04 .       .       GT:DP:AD:RO:QR:AO:QA:GL 0/0:23:22,0:22:784:0:0:0,-6.62266,-65.4388      0/1:12:8>
    Turtle_12806    202     .       TCCAT   TCAT    3.34815e-12     .       .       GT:DP:AD:RO:QR:AO:QA:GL 0/0:100:100,0:100:3679:0:0:0,-30.103,-317.102   0/0:19:19,0:19:697:0:0:0,-5.7195>
    Turtle_28471    169     .       TCCAT   TCAT    0.00145975      .       .       GT:DP:AD:RO:QR:AO:QA:GL 0/0:78:78,0:78:2885:0:0:0,-23.4803,-250.935     0/0:12:12,0:12:454:0:0:0,-3.6123>
    Turtle_30731    160     .       ATTTTTTTAAC     ATTGTTTTAAC     70317.4 .       .       GT:DP:AD:RO:QR:AO:QA:GL 0/1:41:13,23:13:487:23:573:-34.8209,0,-25.8508  0/1:63:44,10:44:1641:10:>
    Turtle_30731    182     .       ATTTTAAAAAAGAGC ATTTTAAAAAAAGAGC        154512  .       .       GT:DP:AD:RO:QR:AO:QA:GL 0/1:31:8,17:8:295:17:443:-28.2536,0,-14.9992    0/1:48:13,21:13:>
    Turtle_32964    206     .       TCCAG   TCAG    513864  .       .       GT:DP:AD:RO:QR:AO:QA:GL 1/1:46:0,46:0:0:46:1699:-116.572,-13.8474,0     1/1:20:0,20:0:0:20:717:-61.6205,-6.0206,>
    Turtle_34114    119     .       TTAAAACT        TTAAACT 11.5794 .       .       GT:DP:AD:RO:QR:AO:QA:GL 0/0:3:3,0:3:109:0:0:0,-0.90309,-9.85365 0/0:1:1,0:1:38:0:0:0,-0.30103,-3.65448  >
    Turtle_36004    81      .       CAAGG   CAGG    536.91  .       .       GT:DP:AD:RO:QR:AO:QA:GL 0/1:91:72,19:72:2673:19:438:-7.09762,0,-202.544 0/1:9:7,2:7:264:2:46:-1.09084,0,-20.3742>
    Turtle_53520    261     .       GCTT    GCTG,GAAG,GGAG  435.583 .       .       GT:DP:AD:RO:QR:AO:QA:GL 0/0:2:2,0,0,0:2:76:0,0,0:0,0,0:0,-0.60206,-6.81445,-0.60206,-6.81445,-6.81445,-0>

##### See how the alleles in these loci differ by numbers of homopolymer repeats? This is called stutter, this is a sequencing artifact (probably), it's pretty easy to detect by eye, and if you're not dealing with an intractible number of genomic loci, just visually scan your .vcf file for these loci and get rid of them.

    vcftools --vcf Freb_filt4.recode.vcf --exclude-positions stuttery_stutterers.txt --recode --out Freb_filt5

##### In fact, I recomend removing indels that are more than a single bp insertion/deletion. 

##### After you've done all your initial filtering, you are ready to move onto microhaplotyping.
