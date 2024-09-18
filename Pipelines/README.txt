README for variant calling pipeline


1)Snakemake-BWA
	-Aligns reads to Wuhan reference, trims primers, gets coverage, and makes consensus sequence
	-Uses SummarizeBAM from Xue et. al 2018 to make the consensus sequence so all within-host references will be on the same coordinate system
	
2)Snakemake-Coverage
	- Processes coverage files, and moves and renames  first high quality consensus to be used as the within-host reference
	- Uses "coverage.R", "Copy_consensus.py" and "Coverage_processing.R"
	
3)Snakemake-Index
	- Indexes (samtools and BWA) within-host consensus that will be used as the within-host reference

4)Snakemake-BWA-Within
	-Realigns sequencing reads to within-host consensus, trims primers, removes reads with mismatch in primer site, calls and filters variants
	- Uses "filter_ivar_variants.R"


Software used
	snakmake
	R/4.2.0 
	ivar 
	samtools 
	fastqc  
	bwa 
	bedtools2 
	htslib/1.9 
	picard-tools/2.8.1 
	python/3.9.12
	
Links
Xue et. al 2018: Within-Host Evolution of Human Influenza Virus, Trends in Microbiology