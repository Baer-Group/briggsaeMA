#!/bin/sh
#SBATCH --job-name=s8 # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=m.rifat@ufl.edu # Where to send mail
#SBATCH --cpus-per-task=8 # Number of CPU cores per task
#SBATCH --ntasks=1 # Run a single task
#SBATCH --mem=20gb # Memory limit
#SBATCH --time=480:00:00 # Time limit hrs:min:sec
#SBATCH --output=s8_%j.out # Standard output and error log
#SBATCH --account=baer --qos=baer
cd /orange/baer/briggsae
### step 8: consolidating the vcf.gz files using GenomicsDBImport:
module load gatk/
gatk --java-options "-Xmx8g" GenomicsDBImport \
-V 202_POOLRET91_S9_L001.bp.g.vcf.gz \
-V 204_POOLRET91_S66_L001.bp.g.vcf.gz \
-V 205_POOLRET91_S74_L001.bp.g.vcf.gz \
-V 206_POOLRET91_S69_L001.bp.g.vcf.gz \
-V 208_POOLRET91_S89_L001.bp.g.vcf.gz \
-V 209_POOLRET91_S92_L001.bp.g.vcf.gz \
-V 212_POOLRET91_S55_L001.bp.g.vcf.gz \
-V 217_POOLRET91_S67_L001.bp.g.vcf.gz \
-V 219_POOLRET91_S49_L001.bp.g.vcf.gz \
-V 228_POOLRET91_S90_L001.bp.g.vcf.gz \
-V 232_POOLRET91_S2_L001.bp.g.vcf.gz \
-V 235_POOLRET91_S81_L001.bp.g.vcf.gz \
-V 236_POOLRET91_S35_L001.bp.g.vcf.gz \
-V 237_POOLRET91_S26_L001.bp.g.vcf.gz \
-V 244_POOLRET91_S85_L001.bp.g.vcf.gz \
-V 250_POOLRET91_S65_L001.bp.g.vcf.gz \
-V 252_POOLRET91_S32_L001.bp.g.vcf.gz \
-V 254_POOLRET91_S72_L001.bp.g.vcf.gz \
-V 255_POOLRET91_S16_L001.bp.g.vcf.gz \
-V 258_POOLRET91_S43_L001.bp.g.vcf.gz \
-V 262_POOLRET91_S56_L001.bp.g.vcf.gz \
-V 263_POOLRET91_S82_L001.bp.g.vcf.gz \
-V 264_POOLRET91_S83_L001.bp.g.vcf.gz \
-V 265_POOLRET91_S21_L001.bp.g.vcf.gz \
-V 269_POOLRET91_S8_L001.bp.g.vcf.gz \
-V 271_POOLRET91_S62_L001.bp.g.vcf.gz \
-V 272_POOLRET91_S38_L001.bp.g.vcf.gz \
-V 278_POOLRET91_S53_L001.bp.g.vcf.gz \
-V 281_POOLRET91_S71_L001.bp.g.vcf.gz \
-V 284_POOLRET91_S88_L001.bp.g.vcf.gz \
-V 286_POOLRET91_S63_L001.bp.g.vcf.gz \
-V 287_POOLRET91_S86_L001.bp.g.vcf.gz \
-V 290_POOLRET91_S51_L001.bp.g.vcf.gz \
-V 298_POOLRET91_S37_L001.bp.g.vcf.gz \
-V 300_POOLRET91_S22_L001.bp.g.vcf.gz \
-V 302_POOLRET91_S1_L001.bp.g.vcf.gz \
-V 303_POOLRET91_S5_L001.bp.g.vcf.gz \
-V 306_POOLRET91_S27_L001.bp.g.vcf.gz \
-V 308_POOLRET91_S84_L001.bp.g.vcf.gz \
-V 310_POOLRET91_S61_L001.bp.g.vcf.gz \
-V 311_POOLRET91_S87_L001.bp.g.vcf.gz \
-V 312_POOLRET91_S3_L001.bp.g.vcf.gz \
-V 313_POOLRET91_S34_L001.bp.g.vcf.gz \
-V 314_POOLRET91_S11_L001.bp.g.vcf.gz \
-V 315_POOLRET91_S6_L001.bp.g.vcf.gz \
-V 316_POOLRET91_S25_L001.bp.g.vcf.gz \
-V 319_POOLRET91_S40_L001.bp.g.vcf.gz \
-V 320_POOLRET91_S60_L001.bp.g.vcf.gz \
-V 323_POOLRET91_S50_L001.bp.g.vcf.gz \
-V 325_POOLRET91_S54_L001.bp.g.vcf.gz \
-V 329_POOLRET91_S76_L001.bp.g.vcf.gz \
-V 331_POOLRET91_S30_L001.bp.g.vcf.gz \
-V 334_POOLRET91_S75_L001.bp.g.vcf.gz \
-V 336_POOLRET91_S95_L001.bp.g.vcf.gz \
-V 337_POOLRET91_S77_L001.bp.g.vcf.gz \
-V 343_POOLRET91_S58_L001.bp.g.vcf.gz \
-V 345_POOLRET91_S4_L001.bp.g.vcf.gz \
-V 346_POOLRET91_S10_L001.bp.g.vcf.gz \
-V 347_POOLRET91_S52_L001.bp.g.vcf.gz \
-V 349_POOLRET91_S14_L001.bp.g.vcf.gz \
-V 350_POOLRET91_S7_L001.bp.g.vcf.gz \
-V 354_POOLRET91_S41_L001.bp.g.vcf.gz \
-V 355_POOLRET91_S96_L001.bp.g.vcf.gz \
-V 357_POOLRET91_S93_L001.bp.g.vcf.gz \
-V 358_POOLRET91_S44_L001.bp.g.vcf.gz \
-V 363_POOLRET91_S20_L001.bp.g.vcf.gz \
-V 365_POOLRET91_S19_L001.bp.g.vcf.gz \
-V 366_POOLRET91_S42_L001.bp.g.vcf.gz \
-V 367_POOLRET91_S68_L001.bp.g.vcf.gz \
-V 368_POOLRET91_S33_L001.bp.g.vcf.gz \
-V 370_POOLRET91_S47_L001.bp.g.vcf.gz \
-V 371_POOLRET91_S24_L001.bp.g.vcf.gz \
-V 373_POOLRET91_S13_L001.bp.g.vcf.gz \
-V 374_POOLRET91_S48_L001.bp.g.vcf.gz \
-V 377_POOLRET91_S73_L001.bp.g.vcf.gz \
-V 379_POOLRET91_S94_L001.bp.g.vcf.gz \
-V 380_POOLRET91_S23_L001.bp.g.vcf.gz \
-V 381_POOLRET91_S91_L001.bp.g.vcf.gz \
-V 383_POOLRET91_S39_L001.bp.g.vcf.gz \
-V 384_POOLRET91_S31_L001.bp.g.vcf.gz \
-V 385_POOLRET91_S29_L001.bp.g.vcf.gz \
-V 386_POOLRET91_S46_L001.bp.g.vcf.gz \
-V 387_POOLRET91_S57_L001.bp.g.vcf.gz \
-V 388_POOLRET91_S78_L001.bp.g.vcf.gz \
-V 390_POOLRET91_S28_L001.bp.g.vcf.gz \
-V 391_POOLRET91_S79_L001.bp.g.vcf.gz \
-V 392_POOLRET91_S64_L001.bp.g.vcf.gz \
-V 394_POOLRET91_S36_L001.bp.g.vcf.gz \
-V 395_POOLRET91_S45_L001.bp.g.vcf.gz \
-V 396_POOLRET91_S17_L001.bp.g.vcf.gz \
-V HK104_ANC_1_POOLRET91_S18_L001.bp.g.vcf.gz \
-V PB800_ANC__1_POOLRET91_S12_L001.bp.g.vcf.gz \
--genomicsdb-workspace-path B1 \
--intervals I --intervals II --intervals III --intervals IV --intervals V --intervals X --intervals MtDNA
