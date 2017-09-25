# Using Anvi'o to Knit Everything Together
Now we are going bring it all together and visuzlize our assembly with [Anvi'o](http://merenlab.org/software/anvio/). Anvi'o is a powerful and extensible tool that might be easily applied to pan-genomic analysis as well as metagenomic analysis. The anvi'o group has a series of fantstic online tutorials, including one on [metagenomic analysis](http://merenlab.org/2016/06/22/anvio-tutorial-v2/). They also run workshops periodically ([schedule here](http://merenlab.org/2016/08/18/events/)) that throughly cover the use of the software.

Today, we are adapting their tutorial on metagenomic analysis to work with the subset dataset that we have.

The goals of this tutorial are to:
* Install anvi'o
* Become familiar with the anvi'o workflow
* Visualizing the assembly with anvi'o
* Become familiar with the anvi'o interface
* Learn how to refine and visuzlize genome bins with anvi'o

## Installing anvi'o (and a few other programs)

The first thing we need to do is install anvi'o. To install anvi'o we will be be using [Anaconda](https://www.continuum.io/what-is-anaconda).

```
cd ~
wget https://repo.continuum.io/archive/Anaconda3-4.4.0-Linux-x86_64.sh
bash Anaconda3-4.4.0-Linux-x86_64.sh
```

Now, follow the prompts for the Anaconda installation. To finish the install it will ask you if you would like to add anaconda to the `$PATH` in you `.bashrc`. You should say yes. Now, you just need to source your `.bashrc` to make sure you can use conda.

```
source .bashrc
```

Anaconda should now be installed. We will now use anaconda (`conda`) to install anvi'o (and all its dependencies) as well as source an environment in which to to run conda.

Now, install anvi'o using conda, create an environment in which to run it, and source the environment:

```
conda create -n anvio232 -c bioconda -c conda-forge gsl anvio
source activate anvio232
```

Anvi'o should now be installed. But, let's double check that it worked. They have a nice little test case to check that everything is working well as follows:

```
anvi-self-test --suite mini
```

This prompt will start anvi'o processing and ultimately it will generate an interactive window with the anvi'o environment. This is accessible through port 8080 (typically, though it might create go to a different port that will be specified) at your IP address.

To find your IP address you can look back at the Jetstream instance page and copy the IP address.

![IP-Address](img/ip-address.png)

Now, open a new tab in your browser* and paste in the following:

*Note: This only works in Google Chrome.

```
[Your IP Address]:8080
```

This should open up the anvi'o interface which is interactive and pretty good looking.

Now, we just need to install a few other programs, namely, samtools and Bowtie2, which we will use for mapping and looking at our mapped data.

```
wget https://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.3.2/bowtie2-2.3.2-linux-x86_64.zip
unzip bowtie2-2.3.2-linux-x86_64.zip

echo 'export PATH=$PATH:~/bowtie2-2.3.2' >> ~/.bashrc
source ~/.bashrc
sudo apt-get -y install samtools
```


Alright, now onto a complete re-analysis of our data with the anvi'o pipeline.

## Getting it into Anvi'o format

Anvi'o takes in 1) an assembly and 2) the raw read data. We have both of those already created, so let's go ahead and download those data (trimmed reads and asssemblies):

```
mkdir ~/anvio-work
cd ~/anvio-work

curl -O https://s3-us-west-1.amazonaws.com/dib-training.ucdavis.edu/metagenomics-scripps-2016-10-12/SRR1976948.abundtrim.subset.pe.fq.gz
curl -O https://s3-us-west-1.amazonaws.com/dib-training.ucdavis.edu/metagenomics-scripps-2016-10-12/SRR1977249.abundtrim.subset.pe.fq.gz
curl -O https://s3-us-west-1.amazonaws.com/dib-training.ucdavis.edu/metagenomics-scripps-2016-10-12/subset_assembly.fa.gz

```

And, gunzip those files:

```
for file in *gz
    do
    gunzip $file
done
```

We now need to get our assembly into the correct format so that anvi'o interpret it.

```
anvi-script-reformat-fasta subset_assembly.fa -o anvio-contigs.fa --min-len 2000 --simplify-names --report name_conversions.txt
```
Take a look at the output files. What has changed?

## Mapping data

We need to map our reads to our anvi'o corrected assembly. This is going to take a little bit of time. First, build the an index for bowtie2:

```
bowtie2-build anvio-contigs.fa anvio-contigs
```

We can write a for loop to map our two datasets and produce .bam files for the files:

```
for file in *fq
do
    bowtie2 --threads 8 -x anvio-contigs --interleaved $file -S ${file/.fq/}.sam
    samtools view -U 4 -bS ${file/.fq/}.sam > ${file/.fq/}.bam
done
```
As above, we need to make these data readable for anvi'o:

```
for file in *.bam
do
    anvi-init-bam ${file} -o ${file/.bam/}.anvio.bam
done

```

## Generating contigs database

In this step we are asking anvi'o to create a database with information about your contigs. The contig database is fairly extensible and can contain lots of different information (taxonomic, functional, etc.). Here, we are primarily asking it to do three things:

1) 'Soft split' long contigs (>20k): Anvi'o shows the generalized statistics for each contig (GC content, etc.). For long contigs these stats are calculated across split contigs (which remain grouped)
2) Identify and locate open reading frames in your contigs (using Prodigal)
3) Estimate Single Copy Gene content (using hmmer against defined gene sets for [bacteria](http://www.pnas.org/content/110/14/5540.full) and [archaea](https://www.nature.com/nature/journal/v499/n7459/full/nature12352.html))
3) Calculate k-mer frequencies for the contigs in our assemblies

So, run the following command to generate the database:

```
anvi-gen-contigs-database -f anvio-contigs.fa -o anvio-contigs.db
```

Then, run this command to perform the hmm search and identify single copy gene content:
```
anvi-run-hmms -c anvio-contigs.db --num-threads 28
```

Now, we can layer on the coverge information from our two samples:

```
for file in *.anvio.bam
do
    anvi-profile -i $file -c anvio-contigs.db -T 28

done
```

And finally, we run the merge step. This will pull all the information together and create a merged anvi'o profile. This step will also run [CONCOCT](https://github.com/BinPro/CONCOCT) (*another* binning algorithm) that will identify bins in our data. Finally, this step calculates the hierarchical relationship betwewen our contigs based on a variety of parameters.
```
anvi-merge *ANVIO_PROFILE/PROFILE.db -o MERGED-SAMPLES -c anvio-contigs.db --enforce-hierarchical-clustering
```

Now we can visualize our data!
```
anvi-interactive -p MERGED-SAMPLES/PROFILE.db -c anvio-contigs.db
```

## Identifing and refining genome bins

First, let's summarize the bin information for our data. This will produce a series of text-based output files detailing some statistics on our genome bins:

```
anvi-summarize -p MERGED-SAMPLES/PROFILE.db -c anvio-contigs.db -o SAMPLES-SUMMARY -C CONCOCT
```
Take a look at the output in `SAMPLES-SUMMARY`. What does it report?

Now you can visuzlize those data in the anvi'o style by simply adding the -C flag to the previous anvi-interactive command:

```
anvi-interactive -p MERGED-SAMPLES/PROFILE.db -c anvio-contigs.db -C CONCOCT
```

Now, we can actually refine the genome bins using anvi'o. This allows us to use human intutition and pattern recognition to better identify contigs that should co-occur.

```
 anvi-refine -p MERGED-SAMPLES/PROFILE.db -c anvio-contigs.db -b Bin_4 -C CONCOCT
```

Finally, it is time to interact with the anvi'o interface. [Here](files/interacting-with-anvio.pdf) are some screenshots to help guide you in your quest. 

---

And of course a big thank you to [Meren](http://merenlab.org/people/) for providing us with extra materials to help create this tutorial! 
