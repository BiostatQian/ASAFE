# Qian Zhang
# 5/02/2016
# This code exemplifies the usage of BEAGLE to obtain phased genotypes, 
# and then the usage of RFMIX to obtained phased admixed ancestries and re-phased genotypes.
# The code below is copied from http://faculty.washington.edu/sguy/local_ancestry_pipeline/rfmix_mds_pipeline
# above the "Apply masking" line.

## This file gives an overview of the process for running RFMix and 
## preparing data for ancestry specific PCA.
## Some customization will be needed.
## This pipeline is for a Unix system, with the bash shell .

bash

# Start by preparing the vcf file, 
# which must have the admixed individuals first, 
# then the first reference, then the second reference, and so on.
#> No example commands are given for this step.

# Make a classes file that has nadmixed 0's, nref1 1's, nref2 2's, etc.  
# The numbers here are numbers of haplotypes (twice the number of individuals).
# The resulting file is space-delimited on one line.
#> Before running these commands you will need to 
#> manually set the numbers of haplotypes (nadmix, nref1, etc).
# This example assumes 3 reference groups 
# (e.g. African, European and Native American).
rm classes # needed only if you already have a file of this name
for i in `seq $nadmix`; do echo "0 " >> classes; done
for i in `seq $nref1`; do echo "1 " >> classes; done
for i in `seq $nref2`; do echo "2 " >> classes; done
let "nref3minus = $nref3 -1"
for i in `seq $nref3minus`; do echo "3 " >> classes; done
echo "3" >> classes

# We will apply the remaining rfmix commands separately to each chromosome, 
# and will then combine the chromosomes before running PCA.

# Phase the vcf file using Beagle v4.0 
# (this can be updated for use with Beagle v4.1 or later). 
for chr in `seq 22`; do
  java -jar beagle.r1399.jar gt=chr${chr}_unphased.vcf.gz out=chr${chr}_phased gprobs=false
done

# Convert phased vcf to alleles file (note maintenance of allele coding)
for chr in `seq 22`; do 
  zcat chr${chr}.vcf.gz | grep -v '#' | cut -f10- | tr -d '\t|' > chr${chr}.alleles
done

# Make the snploc file.
# For this you need a map file, which you can obtain from HapMap
# (be sure to obtain the correct genome build).
# The base2genetic.jar program can be obtained from the Beagle utilities page:
# https://faculty.washington.edu/browning/beagle_utilities/utilities.html
for chr in `seq 22`; do 
  zcat chr${chr}.vcf.gz | grep -v '^#' | cut -f2 | java -jar base2genetic.jar 1 mapfile > chr${chr}.snploc
done

# Run RFMix 
# (the forward-backward option generates posterior probabilities 
# and is not needed for this mds pipeline, 
# but may be useful for other purposes)
dir=`pwd`
cd /usr/local/RFMix-1.5.4/ # change dir to location for RFMIX program
for chr in `seq 22`; do 
  /usr/local/bin/python2.7 RunRFMix.py PopPhased ${dir}/chr${chr}.alleles ${dir}/myfile.classes ${dir}/chr${chr}.snploc --forward-backward -o ${dir}/chr${chr} --disable-parallel -n 5
done
