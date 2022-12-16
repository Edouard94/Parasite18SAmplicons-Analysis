---
output:
  html_document: default
  pdf_document: default
---

# Check if guppy_basecaller is already installed in your machine

Install guppy on a Linux machine (<https://community.nanoporetech.com/docs/prepare/library_prep_protocols/Guppy-protocol/v/gpb_2003_v1_revan_14dec2018/linux-guppy>). Depending on whether we are running the basecalling on CPU GPU, we will need to choose an appropriate Guppy package.

Add Oxford Nanopore's deb repository to your system (this is to install Oxford Nanopore Technologies-specific dependency packages):

    sudo apt update
    sudo apt install wget lsb-release
    export PLATFORM=$(lsb_release -cs)
    wget -O- https://cdn.oxfordnanoportal.com/apt/ont-repo.pub | sudo apt-key add -
    echo "deb http://cdn.oxfordnanoportal.com/apt ${PLATFORM}-stable non-free" | sudo tee /etc/apt/sources.list.d/nanoporetech.sources.list
    sudo apt update

To install the .deb for Guppy, use the following command:

    sudo apt update
    sudo apt install ont-guppy

This will install the GPU version of Guppy.

or:

    sudo apt update
    sudo apt install ont-guppy-cpu

To install the CPU-only version of Guppy. The GPU option works with NVIDIA build (<https://github.com/asadprodhan/GPU-accelerated-guppy-basecalling/blob/main/README.md>) and will need appropriate CUDA drivers to be installed on your system. It is faster to run the basecalling with the GPU, if possible.

We can find your graphic card with:

    sudo lshw -C display

# Basecalling and barcoding of fasta5 files

## Bascalling with guppy_basecaller

To call out a list of available flow cells plus kit combinations, and their associated config files, use the command:

    guppy_basecaller --print_workflows

Here we are using the --flowcell FLO-MIN106 and --kit SQK-LSK109.

If we look for a specific configuration to use, use the commend:

    guppy_basecaller --print_workflows | SQK-LSK109

We can see that the corresponding configuration for our flowcell and kit is dna_r9.4.1_450bps_hac. ![](Qiime2_RMarkdown_insertimage_3.png)

The working directory for the following analysis is /home/Shared/Minion/Minion060522/EdouardInsects_bw/try4/20220509_1222_MN39555_FAT01120_159cc567/guppy_basecalling_Run060522/. In the following document, the path of files seen a '\~/folder/file', \~ refers to the working directory.

To perform the basecalling we use the command:

    guppy_basecaller --disable_pings -i fast5/ -s /fastq_pass_07122022 -c dna_r9.4.1_450bps_hac.cfg --min_qscore 7 --compress_fastq --cpu_threads_per_caller 14 --num_callers 1

The directory contains the following output:

    drwxr-xr-x  2 root root    4096 Dec  7 06:42 fail
    drwxrwxr-x  6 root root    4096 Dec  7 01:21 guppy_basecaller-core-dump-db
    -rw-r--r--  1 root root     180 Dec  7 01:21 guppy_basecaller_log-2022-12-07_01-21-19.log
    -rw-r--r--  1 root root    1005 Dec  7 01:22 guppy_basecaller_log-2022-12-07_01-22-05.log
    -rw-r--r--  1 root root    1007 Dec  7 01:39 guppy_basecaller_log-2022-12-07_01-39-25.log
    -rw-r--r--  1 root root    1059 Dec  7 01:41 guppy_basecaller_log-2022-12-07_01-41-19.log
    -rw-r--r--  1 root root    1061 Dec  7 02:25 guppy_basecaller_log-2022-12-07_02-25-53.log
    -rw-r--r--  1 root root    1105 Dec  7 03:27 guppy_basecaller_log-2022-12-07_03-27-14.log
    -rw-r--r--  1 root root    1103 Dec  7 03:41 guppy_basecaller_log-2022-12-07_03-41-25.log
    -rw-r--r--  1 root root 5242849 Dec  7 05:52 guppy_basecaller_log-2022-12-07_03-43-48.log
    -rw-r--r--  1 root root 1521114 Dec  7 06:42 guppy_basecaller_log-2022-12-07_05-52-17.log
    drwxr-xr-x  2 root root    4096 Dec  7 06:42 pass
    -rw-r--r--  1 root root 5668525 Dec  7 06:42 sequencing_summary.txt
    -rw-r--r--  1 root root  418664 Dec  7 06:42 sequencing_telemetry.js

Ususally, we should merge all resulting 'paased' fastq files into a single file:

    cat *.gz > /home/edouard/pass_cat.fastq.gz

In order to get the number of reads in our fastq file, we can count the number of lines and divide by 4:

    zcat pass_cat.fastq.gz | wc -l | awk '{print $1/4}'

Here we have **20006** read within the pass folder, and 4851 within the fail folder.

We can have a quick sequencing summary and basic statistics with NanoStat (<https://github.com/wdecoster/nanostat>), the input file formats are fastq, bam or albacore sequencing summary format. Here we use the 'sequencing_summary.txt' from the guppy basecalling.

    NanoStat --summary :/home/Shared/Minion/Minion060522/EdouardInsects_bw/try4/20220509_1222_MN39555_FAT01120_159cc567/guppy_basecalling_Run060522/sequencing_summary.txt

![Sequencing (06/05) summary after basecalling](Qiime2_RMarkdown_insertimage_1.png){width="428"}

We can also visualise some plots with the command Nanoplot (<https://github.com/wdecoster/NanoPlot>). NanoPlot creates: a statistical summary, a number of plots and a html summary file. Here we use the 'sequencing_summary.txt' file from the guppy_basecalling to generate the report:

    NanoPlot --summary /home/Shared/Minion/Minion060522/EdouardInsects_bw/try4/20220509_1222_MN39555_FAT01120_159cc567/guppy_basecalling_Run060522/sequencing_summary.txt -o nanoplot_basecalling_Run060522

We can see with the plot 'Read lengths vs. Average read quality' that the majority of our reads have a length of 200-600 bp, with few reads around 1700-1900 bp and a q score varying between 4 and 16. ![NanoPlot](Qiime2_RMarkdown_insertimage_6.png){width="373"}

![NanoPlotZoom](Qiime2_RMarkdown_insertimage_7.png){width="371"}

## Barcoding with guppy_barecoder

After basecalling, we perform the demultiplixing of the previously generated fast.gz files with the command:

    guppy_barcoder --disable_pings -i /home/Shared/Minion/Minion060522/EdouardInsects_bw/try4/20220509_1222_MN39555_FAT01120_159cc567/guppy_basecalling_Run060522/pass -s /home/Shared/Minion/Minion060522/EdouardInsects_bw/try4/20220509_1222_MN39555_FAT01120_159cc567/guppy_basecalling_Run060522/guppy_barcoder_Run060522 --barcode_kits "EXP-NBD104" --require_barcodes_both_ends

NB: the barcodes are by default trimmed by guppy_barcoder command, so the resulting reads are completely demultiplexed. To keep the barcode sequences we need to add to the command line --disable_trim_barcodes.

The output folder of the barcoding looks like this: ![Output folder of guppy_barocder](Qiime2_RMarkdown_insertimage_4.png){width="552"}

The obtained fastq files look like this (e.g. from barcode 01):

    @9583c94e-5237-443b-be2e-da7904ea5e72
    runid=cd2409afdc364972d4e64953023ac179690563d1
    sampleid=try4 
    read=14947
    ch=365
    start_time=2022-05-09T18:37:03Z
    model_version_id=2021-05-17_dna_r9.4.1_minion_384_d37a2ab9
    barcode=barcode01
    TGAATTGACGAGA...

Put the demultiplexed reads into a table for analysis (Tyson et al., 2020):

    fastq pass/barcode{}/*.fastq >> reads.csv

We can convert the demultiplexed fastq files into fasta files to perform blast search.

    seqtk seq -a  barcode10/fastq_runid_cd2409afdc364972d4e64953023ac179690563d1_0.fastq > barcode10/barcode10_Run0506522.fasta

A megablast search is eventually performed with each sample/barcode fasta file on NCBI.

The blast search allows us to align our reads to the nr database available on NCBI. But here we are limited our search to the expected amplified microsporidia and protists, while also including non-targeted groups. The following parameters for the blast search are used. The example with the sample *Aphomia sociella* shows that we include our expected amplified organisms (gregarines, kinetoplastids, amoebozoa) while including non-targeted sequences to check the specificity of the used primers.![Plot title.](Qiime2_RMarkdown_insertimage_16.png){width="493"} ![Plot title.](Qiime2_RMarkdown_insertimage_19.png){width="438"}

The hit table is screened for alignment that has a % identity superior at 95% and a alignment length of minimum 150bp.

Then the accession number are used to have the corresponding organism with the command:

    epost -db nuccore -input barcode4.txt | esummary
    | xtract -pattern DocumentSummary -element AccessionVersion,Organism

Rather than using the blast search with our individual reads, we can perform an OTU classification. This relies on databases cataloging reference OTUs identified through clustering (such as PR2, a database specific for 18S protists sequences). This can be achieved with the tool Qiime2.

# Qiime2 (Bolyen et al., 2019)

To install Qiime2, we can follow the instructions from here <https://docs.qiime2.org/2022.8/install/>.

Initiate qiime2 environment to start the analysis:

    conda activate qiime2-2022.8
    qiime --help

To deactivate the quiime environment:

    conda deactivate

## Import of data in Qiime2

Qiime is using 'artifacts'. Artifacts represent intermediate data in a qiime analysis, that are "intended to be consumed by QIIME 2". Intermediate artifacts have the extension .qza, which stands for QIIME Zipped Artifact. Whereas, a terminal output ready for "human consumption" will typically have the extension .qzv, which means QIIME Zipped Visualization (<https://gregcaporaso.github.io/q2book/back-matter/glossary.html?highlight=qza>).

The first step for our qiime analysis is to import our data with:

    qiime tools import

To import our data (i.e., demultiplexed single end fatsq files) we need to create a manifest file. A 'manifest file' is needed to import fastq files that do not follow specific formats (i.e. EMP or Casava formats). Since we are importing demultiplexed fastq files (without barcode sequences), we need to import our own "Fastq manifest" (<https://docs.qiime2.org/2022.8/tutorials/importing/>). Example of manifest format:

![ManifestFile](Qiime2_RMarkdown_insertimage_8.png)

Our manifest file only have two columns, the sample-id (i.e. barcode) and the path for the demultiplexed fastq files (single reads).

It seems more appropriate to perform the qiime analysis with the barcode/sample treated independently (Kenmotsu et al., 2021). Then our data are independently imported into QIIME2. In this way, we can use specific classifier for each barcode (see further).

In the folder '\~/guppy_barcoder_Run060522/barcode01' we have the demultiplexed fastq file that we will import in qiime, thanks to a manifest, with the command line:

    qiime tools import --type 'SampleData[SequencesWithQuality]' --input-path ManifestBarcode01_Run060522.txt --output-path AsNB01_Run060522.qza --input-format SingleEndFastqManifestPhred33V2

The same import procedure is applied to the rest of the barcodes.

We can look at the visualisation of all the barcodes (qiime artifacts saved in '\~/Qiime2_Run060522_121222/')', we can see that the numbers of reads per barcode vary from 32 (barcode 03, *Galleria mellonella*) to 224 (barcode 04, *Schistocerca gregaria*, only SSU amplicons).

![Qiime2_Barcodes_Run060522_visu.qzv](Qiime2_RMarkdown_insertimage_9.png){width="418"}

We can also import all the barcodes at once.

    qiime tools import --type 'SampleData[SequencesWithQuality]' --input-path ManifestBarcodeRun060522.txt --output-path Qiime2_QZA_Barcodes_Run060522.qza --input-format SingleEndFastqManifestPhred33V2

To visualise the data we used the command line:

    qiime demux summarize --i-data Qiime2_QZA_Barcodes_Run060522.qza --o-visualization Qiime2_Barcodes_Run060522_visu.qzv

And the Qiime2_Barcodes_Run060522_visu.qzv file has been uploaded on the website <https://view.qiime2.org/>. We can also observe the Interactive Quality Plot. ![Interactive Quality Plot All Barcodes Run 060522](Qiime2_RMarkdown_insertimage_10.png){width="446"} ![Plot title.](Qiime2_RMarkdown_insertimage_17.png){width="403"}

##Denoising with dada2

    qiime dada2 denoise-single \
    --p-n-threads 12 \
    --i-demultiplexed-seqs Qiime2_QZA_Barcodes_Run060522.qza \
    --p-trunc-len 0 \
    --p-trim-left 100 \
    --output-dir DADA2_denoising_Barcodes_Run060522 \
    --verbose \

![Command dada2](Qiime2_RMarkdown_insertimage_11.png){width="566"}

Good to know that we have the following default parameters ![qiime dada2 denoise-single 'pooling' and 'chimera' default parameters](Qiime2_RMarkdown_insertimage_12.png){width="393"}

The output of *qiime dada2 denoise-single* on the demultiplexed data is composed by a feature table *FeatureTable[Frequency]*, the resulting feature sequences *FeatureData[Sequence]*, and the denoising stats *SampleData[DADA2Stats]*.

### Generate summary files¶

With the outputs, a metadata file is required which provides the key to gaining biological insight from your data. The file metadata.tsv Things to look for in the summary files:

Feature table

    qiime feature-table summarize \
    --i-table table.qza \
    --o-visualization table.qzv \

How many sequences have been filtered?

    qiime metadata tabulate \
    --m-input-file denoising_stats.qza \
    --o-visualization denoising_stats_Barcodes_Run060522.qzv \
    --verbose

![Visualisation of denoising_stats_Barcodes_Run060522.qzv](Qiime2_RMarkdown_insertimage_13.png){width="541"}

Do BLAST searches of the representative sequences make sense? blast the sequences in *FeatureData[Sequence]*. To have the Sequence Length Statistics with our different sequences, we use the following command:

    qiime feature-table tabulate-seqs \
    --i-data representative_sequences.qza \
    --o-visualization representative_seqs_Barcodes_Run060522.qzv \
    --verbose

Only 8 sequences (from 792 in ) in *representative_seqs_Barcodes_Run060522.qzv*. We have a large number (\>50%) of sequences that has been lost during denoising/filtering. The settings are too stringent. The new parameters are used:

    qiime dada2 denoise-single \
    --p-n-threads 12 \
    --i-demultiplexed-seqs Qiime2_QZA_Barcodes_Run060522.qza \
    --p-trunc-len 0 \
    --p-max-ee 10 \
    --p-trunc-q 0 \
    --p-chimera-method 'none' \
    --output-dir DADA2_denoising_Barcodes_Run0605224 \
    --verbose \

The filtering is softened when **max-ee** is increased (<https://github.com/benjjneb/dada2/issues/248>). Some reads correspond to gregarines, but still too few reads are available.

Then, we will skip the denoising step as we are losing too much reads (we aslo have very few reads from the sequencing run). We will classify our demultiplexed reads with the PR2 database (an 18S database for protists).

## Use the PR2 Database in Qiime2 for classifying

### Import the PR2 Database and train the classifier

The following commands can be used to import the PR2 reference database and its corresponding taxonomy as qimme artifacts <https://pr2-database.org/documentation/pr2-files/>:

    qiime tools import --type 'FeatureData[Sequence]' --input-path pr2_version_4.14.0_SSU_mothur.fasta --output-path PR2_otus.qza

    qiime tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat --input-path pr2_version_4.14.0_SSU_mothur.tax --output-path PR2-taxonomy.qza

After importing our reference database, we need to train our 'classifier' against reference sequences (e.g. with the PR2 data):

    qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads PR2_otus.qza --i-reference-taxonomy PR2-taxonomy.qza --o-classifier PR2classifier.qza

The resulting classifier will serve to classify our own reads.

    qiime feature-classifier classify-sklearn --i-classifier /home/Shared/PR2database/PR2classifier.qza --i-reads Qiime2_QZA_Barcodes_Run060522.qza --o-classification PR2-taxonomy-Barcodes_Run060522.qza

We have the following error *Invalid value for '--i-reads': Expected an artifact of at least type FeatureData[Sequence]. An artifact of type SampleData[SequencesWithQuality] was provided*. To resolve this we need to convert our artifact Qiime2_QZA_Barcodes_Run060522.qza from SampleData[SequencesWithQuality] to FeatureData[Sequence].

### Dereplicate the sequences

We dereplicate our SampleData artifact with vsearch cmd. It will dereplicate our SampleData and indicate the number of times each amplicon sequence variant (ASV) is observed in each of our samples. No quality filtering is done with this tool, so we do not loose information.

    qiime vsearch dereplicate-sequences --i-sequences Qiime2_QZA_Barcodes_Run060522.qza --o-dereplicated-table dereplicatedtable_Qiime2_QZA_Barcodes_Run060522.qza --o-dereplicated-sequences Qiime2_QZA_Barcodes_Run060522_dereplicated.qza

### Classification with trained classifier

We reperform the classification with the dereplicated sequence.

    qiime feature-classifier classify-sklearn --i-classifier /home/Shared/PR2database/PR2classifier.qza --i-reads Qiime2_QZA_Barcodes_Run060522_dereplicated.qza --o-classification PR2-taxonomy-Barcodes_Run060522.qza

### Visualize our taxonomies by entering the following

We can now visualise the obtained OTU classification with:

    qiime metadata tabulate --m-input-file PR2-taxonomy-Barcodes_Run060522.qza --m-input-file Qiime2_QZA_Barcodes_Run060522_dereplicated.qza --o-visualization taxonomy-Barcodes_Run060522.qzv

### Taxonomic bar plots were generally created as follows:

The data was filtered to remove bacterial sequences (designated simply "Eukaryota"). Remove Eukaryota by exact match with:

    qiime taxa filter-table \
        --i-table dereplicatedtable_Qiime2_QZA_Barcodes_Run060522.qza \
        --i-taxonomy PR2-taxonomy-Barcodes_Run060522.qza \
        --p-mode exact \
        --p-exclude "Eukaryota" \
        --o-filtered-table Barcodes_Run060522_filtered.qza

Command for barplots:

    qiime taxa barplot \
    --i-table Barcodes_Run060522_filtered.qza \
    --i-taxonomy PR2-taxonomy-Barcodes_Run060522.qza \
    --o-visualization Barcodes_Run060522_barplots_filtered.qzv

We can check the remaining read counts after filtering:

    qiime feature-table summarize \
    --i-table Barcodes_Run060522_filtered.qza \
    --o-visualization Barcodes_Run060522_filtered.qzv

## Other taxonomy classification option (amplicon specific classifier)

Then, to classify our reads, we will need a reference database.We can use the plugin qiime rescript (<https://github.com/bokulich-lab/RESCRIPt>) to easily import the SILVA database (Robeson et al., 2021):

    qiime rescript get-silva-data --p-version '138.1' --p-target 'LSURef_NR99' --p-include-species-labels --o-silva-sequences silva-138.1-lsu-nr99-seqs.qza --o-silva-taxonomy silva-138.1-lsu-nr99-tax.qza

We have to convert the rRNA sequences into rDNA:

    qiime rescript reverse-transcribe --i-rna-sequences silva-138.1-ssu-nr99-seqs.qza --o-dna-sequences silva-138.1-ssurDNA-nr99-seqs.qza

As we are analysing our samples individually, we will use amplicon-specific classifiers. For each sample, we will generate amplicon-specific classifiers (<https://forum.qiime2.org/t/processing-filtering-and-evaluating-the-silva-database-and-other-reference-sequence-data-with-rescript/15494#heading--sixth-header>). Constructing such classifiers allow for more robust taxonomic classification of our data (Werner et al. 2011, Bokulich et al. 2018).

###Barcode 01

The barcode 01 represents the sample *Aphomia sociella* which had expected amplicons obtained with the group specific primers: - Gregarines (SSU): forward-CCCTTAGATRRYCTGGGCTGC and reverse-CGTGTTACGACTTCTTC - Kinetoplastea (SSU): forward-CTGCCAGTAGTCATATGCTTGTTTCAAGGA and reverse-GATCCTTCTGCAGGTTCACCTACAGCT - Amoebozoa (SSU): forward-GAATTGACGGAAGGGCACAC and reverse-CCAAGAYRTCTAAGGGCATCAC

It is recommended to trim the primers from the data with dada2 prior classification, however this can be discussed: <https://www.researchgate.net/post/Is_it_mandatory_to_remove_primers_and_barcodes_from_reads_before_starting_the_analyzes_in_qiime_2>.

The primers can be removed with the command (this step is omitted here):

    qiime cutadapt trim-paired \
    --i-demultiplexed-sequences file.qza \
    --p-front-f forward_primer_sequence \
    --p-front-r reverse_primer_sequence \
    --p-error-rate 0.20 \
    --output-dir analysis/seqs_trimmed \
    --verbose

The error rate parameter, --p-error-rate, will likely need to be adjusted to get 100% (or close to it) of reads trimmed.

####Classifier for Gregarines (SSU)

To generate our amplicon-specific classifier, we use the silva database silva-138.1-SSUrDNA-nr99-seqs.qza (/home/edouard/Qiime2Databases_Feature-classifier):

    qiime feature-classifier extract-reads \
        --i-sequences /home/edouard/Qiime2Databases_Feature-classifier/silva-138.1-SSUrDNA-nr99-seqs.qza \
        --p-f-primer CCCTTAGATRRYCTGGGCTGC \
        --p-r-primer CGTGTTACGACTTCTTC \
        --p-min-length 50 \
        --p-max-length 1000 \
        --p-n-jobs 2 \
        --p-read-orientation 'forward' \
        --o-reads silva-ssu-gregFgregR.qza

Now we can build our amplicon-region specific classifier.

    qiime feature-classifier fit-classifier-naive-bayes \
        --i-reference-reads silva-ssu-gregFgregR.qza \
        --i-reference-taxonomy /home/edouard/Qiime2Databases_Feature-classifier/silva-138.1-ssu-nr99-tax.qza \
        --o-classifier silva-ssu-gregFgregR-classifier.qza
        

The output *silva-ssu-gregFgregR-classifier.qza* is a Taxonomic Classifier that we will use for our gregarine amplicons classification. Classify rep seqs with:

    qiime feature-classifier classify-sklearn\
    --i-classifier silva-ssu-gregFgregR-classifier.qza\
    --i-reads AsNB01_Run060522.qza\
    --o-classification classified_AsNB01_Run060522.qza

We have the error *Invalid value for '--i-reads': Expected an artifact of at least type FeatureData[Sequence]. An artifact of type SampleData[SequencesWithQuality] was provided*. We dereplicate our SampleData artifact Qiime2_Allbarcodes.qza with vsearch cmd. It will dereplicate our SampleData and indicate the number of times each amplicon sequence variant (ASV) is observed in each of our samples. No quality filtering is done with this tool, so we do not loose information.

    qiime vsearch dereplicate-sequences --i-sequences AsNB01_Run060522.qza --o-dereplicated-table dereplicatedtable_AsNB01_Run060522.qza --o-dereplicated-sequences AsNB01_Run060522_dereplicated.qza

We reperform the classification with the dereplacted sequence.

    qiime feature-classifier classify-sklearn\
    --i-classifier silva-ssu-gregFgregR-classifier.qza\
    --i-reads AsNB01_Run060522_dereplicated.qza\
    --o-classification classified_AsNB01_Run060522.qza

Finally, we visualise our taxonomies and the confidence of taxonomy assignment by entering the following command:

    qiime metadata tabulate\
    --m-input-file classified_AsNB01_Run060522.qza\
    --o-visualization classified_AsNB01_Run060522_visu.qzv \`\`\` Most of the reads are just classified as eukaryotes. ![Plot title.](Qiime2_RMarkdown_insertimage_14.png)

# Useful cmd

Find the location of a folder:

    locate folder | grep folder$
    locate -b filename

We can count how many times appear a specific subsequence into fastq files:

    > zgrep -c 'ATGATGATG' reads.fq.gz
        398065    
    > zcat reads.fq.gz | awk '/ATGATGATG/ {nlines = nlines + 1} END {print nlines}'
        398065

## Make a database for blast

Download fasta file from NCBI search Alveolata[ORGN] OR Ciliophora[ORGN] OR Amoebozoa[ORGN] OR Helicosporidium[ORGN] OR Euglenozoa[ORGN] OR Ichthyosporea[ORGN] OR Microsporidia[ORGN] AND 28S[ALL] AND 18S[ALL] (Nucleotide 7,393) ##joint with PR2 database

    makeblastdb -in MajorProtists18Sand28S_NCBI_Nucleotide.fa -title MajorProtists18Sand28S_DB -dbtype nucl -out MajorProtists18Sand28S_DB -parse_seqids
