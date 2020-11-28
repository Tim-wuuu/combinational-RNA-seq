#!/usr/bin/perl -w

############################################################
### Given an input file containing the path of folder where all the fastq has been placed, 
### create and kick off slurm scripts to run optionally cominbationl pipeline
############################################################

use strict;
use File::Basename;
use POSIX qw(strftime);
use Cwd qw(abs_path);

my $usage = "\nUsage: perl Tim_RNA_seq.pl <sampleFastqList.txt>\n\nOptions:\n[Prefix]  The labedl of sumitted job\n[DATA_TYPE]  PE/SE \n[QUANTIFICATION]  STAR/HTseq/Rsem\n[DEanalysis]  DEseq2/edgeR\n[EXP_DESIGN_PATH] The experiment design file should be a csv file.\n [FASTQ_LOCATION] If it's SE, then one folder should be provided; if it's PE, then two folder containing _1 and _2 .fasq files should be provided.\n[OUTPUTDIR] The location of output files\n\n The input FastqList should have the following format:\n\nPrefix DATA_TYPE QUANTIFICATION DEanalysis EXP_DESIGN_PATH FASTQ_LOCATION  OUTPUTDIR\n\n"; 

my $abspath = abs_path();
my $tstamp = strftime("%Y-%m-%d-%H-%M",localtime);
my $input = shift || die $usage;

# hg38 GTF file
my $GTF_NAME = "GTF_test";

## SOFTWARE PACKAGES
# Current default 
my $STAR = "/home/rcf-proj/sr1/rhie/bioinformatic_tools/STAR/source/STAR";

my $Rsem;

my $HTSeq = "/auto/rcf-proj3/sr1/rhie/bioinformatic_tools/Python-2.7.11/venv_2.7.11_HTSeq/bin/";

my $GenomeDir = "/home/rhie-00/suhn/rhie/indepth.RNAseq.analysis/Index/STAR/hg38";


my $Prefix;
my $DataType;
my $Quantification;
my $DEanalysis;
my $EXP_DESIGN_PATH;
my $Fastq_Path_SE;
my $Fastq_Path_1;
my $Fastq_Path_2;
my $OUTPUTDIR;

open (INPUT_LIST, "<$input") or die "Couldn't open $input for reading\n";
# Full path of input file.

my @fields;
my $numFields; 
while (<INPUT_LIST>) {
   # skip blank lines and comment lines beginning with hash
   next if m/^(#|\s*$)/;
   $numFields = scalar @{[split '\s+', $_]};
   unless($numFields == 7 || $numFields == 8) {
print "\nERROR!:  Missing column on line:\n$_\n\nInput line should contain the following 6 or 7 columns separated by whitespace:\nPrefix DATA_TYPE QUANTIFICATION DEanalysis EXP_DESIGN_PATH FASTQ_LOCATION\n\n";
}
   @fields = split '\s+', $_;
}
close INPUT_LIST;

$Prefix = $fields[0];
$DataType = $fields[1];
$Quantification = $fields[2];
$DEanalysis = $fields[3];
$EXP_DESIGN_PATH = $fields[4];

my $N_slurm;
my $trimJob;
my $starJob;
my $DIRPREFIX;
if ($DataType eq 'SE' && $numFields == 7) {
$Fastq_Path_SE = $fields[5];
$OUTPUTDIR = $fields[6];
print "\nPrefix:$Prefix\nDataType:$DataType\nQuantification:$Quantification\nDEanalysis:DEanalysis:$DEanalysis\nEXP_DESIGN_PATH:$EXP_DESIGN_PATH\nFastq_Path_SE:$Fastq_Path_SE\nOUTPUTDIR:$OUTPUTDIR\n\n";

# Create sample dir
`mkdir $OUTPUTDIR`;
$DIRPREFIX = "$OUTPUTDIR/$Prefix";
`mkdir $DIRPREFIX`;

my $N_Fastq_Path_SE=`ls $Fastq_Path_SE/*fastq $Fastq_Path_SE/*fq | wc -l`;
print "There are $N_Fastq_Path_SE single-end files in total\n";
$N_slurm=$N_Fastq_Path_SE-1;
####trimming and mapping
#STEP1 Trimgalore_SE
createSLURM_TRIM_SE();
$trimJob = "$Prefix-SimpleRNAseqPipeline-trimgalore-singleEnd.slurm";

#STEP2 STAR_SE
createSLURM_STAR_SE();
$starJob = "$Prefix-SimpleRNAseqPipeline-star-singleEnd.slurm";
}

if ($DataType eq 'PE' && $numFields == 8) {
$Fastq_Path_1 = $fields[5];
$Fastq_Path_2 = $fields[6];
$OUTPUTDIR = $fields[7];
print "\nPrefix:$Prefix\nDataType:$DataType\nQuantification:$Quantification\nDEanalysis:$DEanalysis\nEXP_DESIGN_PATH:$EXP_DESIGN_PATH\nFastq_Path_1:$Fastq_Path_1\nFastq_Path_2:$Fastq_Path_2\nOUTPUTDIR:$OUTPUTDIR\n\n";

# Create sample dir
`mkdir $OUTPUTDIR`;
$DIRPREFIX = "$OUTPUTDIR/$Prefix";
`mkdir $DIRPREFIX`;

### For PE
my $N_Fastq_Path_PE=`ls $Fastq_Path_1/*fastq $Fastq_Path_1/*fq | wc -l`;
print "There are $N_Fastq_Path_PE paired-end files in total\n";
$N_slurm=$N_Fastq_Path_PE-1;

#STEP1 Trimgalore_PE
createSLURM_TRIM_PE();
$trimJob = "$Prefix-SimpleRNAseqPipeline-trimgalore-PairedEnd.slurm";

#STEP2 STAR_PE
createSLURM_STAR_PE();
$starJob = "$Prefix-SimpleRNAseqPipeline-star-PairedEnd.slurm";
}


####Quantification
my $quantificationJob;
if($Quantification eq "STAR") {

print "Use .readsPergene.Log.out for DEanalysis";

} elsif($Quantification eq "HTseq") {

print "Use HTseq for DEanalysis";
##STEP3 HTseq
createSLURM_HTseq();
$quantificationJob = "$Prefix-SimpleRNAseqPipeline-htseq.slurm";

} elsif($Quantification eq "Rsem") {

print "Use Rsem for DEanalysis";
##STEP3 Rsem
createSLURM_Rsem();
$quantificationJob = "$Prefix-SimpleRNAseqPipeline-rsem.slurm";

} else {
die "Error: wrongly specifying quantification tools. The option should be STAR/HTseq/Rsem.";
}

###Submit job
#trim
my $trimJobID = `sbatch $trimJob`;
my $trimIDNumber = @{[$trimJobID =~ m/\w+/g]}[3];

#Mapping
my $starJobID = `sbatch --dependency=afterok:$trimIDNumber $starJob`;
my $starIDNumber = @{[$starJobID =~ m/\w+/g]}[3];

#Quantification 
my $quantificationJobID = `sbatch --dependency=afterok:$starIDNumber $quantificationJob`;
my $quantificationIDNumber = @{[$quantificationJobID =~ m/\w+/g]}[3];



####The content of slurm scripts
sub createSLURM_TRIM_SE {
# Create SLURM script to run Trimgalore filtering for single-end data
open(SLURMJOB, " > $Prefix-SimpleRNAseqPipeline-trimgalore-singleEnd.slurm") || die "$!\n";
    print SLURMJOB <<PBS;
#!/bin/bash
#SBATCH --array=0-$N_slurm
#SBATCH --job-name=TrimGalore.se
#SBATCH --output=TrimGalore.se.out
#SBATCH --error=TrimGalore.se.err
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=16GB
#SBATCH --mail-user=beoungle\@usc.edu
#SBATCH --mail-type=ALL


module load python
module load fastqc
module load trimgalore

cd $Fastq_Path_SE
FILES=(\$(ls {*fastq.gz,*.fastq}))
INPUT=\${FILES[\$SLURM_ARRAY_TASK_ID]}
OUTPUTDIR=(\$(echo \${FILES[\$SLURM_ARRAY_TASK_ID]} | sed -e 's/\.[^./]*\$//')) ##remove the extention

mkdir $DIRPREFIX/Trim_SE

trim_galore --phred33 \\
--output_dir $DIRPREFIX/Trim_SE/\$OUTPUTDIR \$INPUT

PBS
close(SLURMJOB) || die "$!\n";
}

sub createSLURM_STAR_SE {
# Create SLURM script to run STAR alignment for single-end data
open(SLURMJOB, " > $Prefix-SimpleRNAseqPipeline-star-singleEnd.slurm") || die "$!\n";
    print SLURMJOB <<PBS;
#!/bin/bash
#SBATCH --job-name=STAR_SE
#SBATCH --output=STAR_SE.out
#SBATCH --error=STAR_SE.err
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=32GB
#SBATCH --mail-user=beoungle\@usc.edu
#SBATCH --mail-type=ALL   

module load star
mkdir $DIRPREFIX/STAR_Alignment
cd $DIRPREFIX/Trim_SE

for i in \${ls}
do
cd \${i}
STAR --runThreadN 6 \\
--genomeDir /project/rhie_130/suhn/beoungle/testing_rnaseq/STAR.Index \\
--readFilesIn *.fq \\
--outFileNamePrefix $DIRPREFIX/STAR_Alignment/\${i} \\
--outSAMtype BAM SortedByCoordinate \\
--outSAMunmapped Within \\
--quantMode TranscriptomeSAM GeneCounts
cd ..
done

cd $DIRPREFIX/
mkdir htseq

PBS
close(SLURMJOB) || die "$!\n";
}

sub createSLURM_TRIM_PE {
# Create SLURM script to run Trimgalore filtering for paired-end data
open(SLURMJOB, " > $Prefix-SimpleRNAseqPipeline-trimgalore-PairedEnd.slurm") || die "$!\n";
    print SLURMJOB <<PBS;
#!/bin/bash
#SBATCH --array=0-$N_slurm
#SBATCH --job-name=TrimGalore_PE
#SBATCH --output=TrimGalore_PE.out
#SBATCH --error=TrimGalore_PE.err
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=16GB
#SBATCH --mail-user=beoungle\@usc.edu
#SBATCH --mail-type=ALL

module load python
module load fastqc
module load trimgalore

cd $Fastq_Path_1
FILES_1=(\$(find `pwd` -type f -name "*fq" -o -name "*fastq" | sort))
FILES=(\$(ls *fq *fastq | sort))
cd $Fastq_Path_2
FILES_2=(\$(find `pwd` -type f -name "*fq" -o -name "*fastq" | sort))

INPUT_1=\${FILES_1[\$SLURM_ARRAY_TASK_ID]}
INPUT_2=\${FILES_2[\$SLURM_ARRAY_TASK_ID]}
OUTPUTDIR=(\$(echo \${FILES[\$SLURM_ARRAY_TASK_ID]} | sed -e 's/\.[^./]*\$//')) ##remove the extention

mkdir $DIRPREFIX/Trim_PE
mkdir $DIRPREFIX/Trim_PE/Trimmed_files_PE

trim_galore --phred33 --fastqc \\
--output_dir $DIRPREFIX/Trim_PE/Trimmed_files_PE/\$OUTPUTDIR \\
--paired \\
\$INPUT_1 \$INPUT_2

PBS
close(SLURMJOB) || die "$!\n";
}

sub createSLURM_STAR_PE {
open(SLURMJOB, " > $Prefix-SimpleRNAseqPipeline-star-PairedEnd.slurm") || die "$!\n";
    print SLURMJOB <<PBS;
#!/bin/bash
#SBATCH --job-name=STAR_Alignment_PE
#SBATCH --output=STAR_Alignment_PE.out
#SBATCH --error=STAR_Alignment_PE.err
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=32GB
#SBATCH --mail-user=beoungle\@usc.edu
#SBATCH --mail-type=ALL

module load star

mkdir $DIRPREFIX/STAR_Alignment
cd $DIRPREFIX/Trim_PE/Trimmed_files_PE
mkdir TrimmedFiles_1 TrimmedFiles_2
mv \$(find $DIRPREFIX/Trim_PE -type f -name "*1.fq" -o -name "*1.fastq") ./TrimmedFiles_1/
mv \$(find $DIRPREFIX/Trim_PE -type f -name "*2.fq" -o -name "*2.fastq") ./TrimmedFiles_2/

FILES_1=(\$(find ./TrimmedFiles_1 -type f -name "*1.fq" -o -name "*1.fastq" | sort))
FILES_2=(\$(find ./TrimmedFiles_2 -type f -name "*2.fq" -o -name "*2.fastq" | sort))
Prefix=(\$(ls ./TrimmedFiles_1 | sed -e 's/\.[^./]*\$//' | sort))

for i in {0..$N_slurm};
do
STAR --runThreadN 6 \\
--genomeDir /project/rhie_130/suhn/beoungle/testing_rnaseq/STAR.Index \\
--readFilesIn \${FILES_1[i]} \${FILES_2[i]} \\
--outFileNamePrefix $DIRPREFIX/STAR_Alignment/\${Prefix[i]} \\
--outSAMtype BAM SortedByCoordinate \\
--outSAMunmapped Within \\
--quantMode TranscriptomeSAM GeneCounts
done

mkdir $DIRPREFIX/STAR_Quantification
mv \$(find $DIRPREFIX/STAR_Alignment -type f -name "*readsPergene.Log.out") $DIRPREFIX/STAR_Quantification
PBS
close(SLURMJOB) || die "$!\n";
}

sub createSLURM_HTseq {
open(SLURMJOB, " > $Prefix-SimpleRNAseqPipeline-htseq.slurm") || die "$!\n";
    print SLURMJOB <<PBS;
#!/bin/bash
#SBATCH --job-name=HTseq.counts
#SBATCH --array=0-$N_slurm
#SBATCH --output=HTseq.out
#SBATCH --error=HTseq.err
#SBATCH --time=23:59:59
#SBATCH --mem-per-cpu=16GB
#SBATCH --mail-user=beoungle\@usc.edu
#SBATCH --mail-type=ALL

module load python
cd $DIRPREFIX/STAR_Alignment/
FILES=(\$(ls *.sortedByCoord.out.bam))
INPUT=\${FILES[\$SLURM_ARRAY_TASK_ID]}
Prefix=(\$(echo \${FILES[\$SLURM_ARRAY_TASK_ID]} | cut -d '.' -f 1))

htseq-count -s no -r pos -f bam \$INPUT \\
/project/rhie_130/suhn/beoungle/testing_rnaseq/gencode.v35.annotation.gtf > \\
$DIRPREFIX/htseq/\${Prefix}.counts
PBS
close(SLURMJOB) || die "$!\n";
}

sub createSLURM_Rsem {
open(SLURMJOB, " > $Prefix-SimpleRNAseqPipeline-rsem.slurm") || die "$!\n";
    print SLURMJOB <<PBS;
#!/bin/bash
#SBATCH --job-name=Rsem.count
#SBATCH --output=Rsem.count.out
#SBATCH --error=Rsem.count.err
#SBATCH --array=0-1
#SBATCH --account=lc_sr2
#SBATCH --partition=laird
#SBATCH --time=19:59:59
#SBATCH --mem-per-cpu=16GB
#SBATCH --mail-user=zexunwu\@usc.edu
#SBATCH --mail-type=ALL

mkdir $DIRPREFIX/Rsem_Quantification
cd $DIRPREFIX/STAR_Alignment
FILES=(\$(find ./ -type f -name "*.sortedByCoord.out.bam"))
INPUT=\${FILES[\$SLURM_ARRAY_TASK_ID]}
Prefix=

/home/rhie-00/suhn/zexunwu/software/RSEM-1.3.3/rsem-calculate-expression \\
-p 6 --bam --no-bam-output \\
\$INPUT \\
/home/rhie-00/suhn/zexunwu/ref/hg19/hg19.Rsem.index/hg19.Rsem.index \\
$DIRPREFIX/Rsem_Quantification/\$Prefix

PBS
close(SLURMJOB) || die "$!\n";
}
