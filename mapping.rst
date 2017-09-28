=======
Mapping
=======

Download bwa::

  cd
  sudo apt-get install bwa samtools

Downloading data
-----------------

Now, go to a new directory and grab the data::

  mkdir ~/mapping
  cd ~/mapping

  curl -O https://s3-us-west-1.amazonaws.com/dib-training.ucdavis.edu/metagenomics-scripps-2016-10-12/SRR1976948.abundtrim.subset.pe.fq.gz
  curl -O https://s3-us-west-1.amazonaws.com/dib-training.ucdavis.edu/metagenomics-scripps-2016-10-12/SRR1977249.abundtrim.subset.pe.fq.gz

And extract the files::

  for file in *fq.gz
    do
    gunzip $file
  done

We will also need the assembly; rather than rebuilding it, you can download
a copy that we saved for you::

  curl -O https://s3-us-west-1.amazonaws.com/dib-training.ucdavis.edu/metagenomics-scripps-2016-10-12/subset_assembly.fa.gz
  gunzip subset_assembly.fa

Mapping the reads
-----------------

First, we will need to to index the megahit assembly::

  bwa index subset_assembly.fa

to The reads are in paired-end/interleaved format, so you'll need to add the -p flag to indicate to bwa that these are paired end data::

Map the reads::

  for i in *fq
    do
    bwa mem -p subset_assembly.fa $i > ${i}.aln.sam
  done

Converting to BAM to visualize
------------------------------

First, index the assembly for samtools::

  samtools faidx subset_assembly.fa

Then, convert both SAM files to BAM files::

  for i in *.sam
  do
     samtools import subset_assembly.fa $i $i.bam
     samtools sort $i.bam -o $i.bam.sorted.bam
     samtools index $i.bam.sorted.bam
  done

Visualizing the read mapping
----------------------------

Find a contig name to visualize::

    grep -v ^@ SRR1976948.abundtrim.subset.pe.fq.aln.sam | cut -f 3 | sort | uniq -c | sort -n | tail

Pick one e.g. k99_13588.

Now execute::

  samtools tview SRR1976948.abundtrim.subset.pe.fq.aln.sam.bam.sorted.bam subset_assembly.fa -p k99_13588:400

(use arrow keys to scroll, 'q' to quit; a key for what you are looking at: `pileup format<https://en.wikipedia.org/wiki/Pileup_format>`__.)

Look at it in both mappings::

  samtools tview SRR1977249.abundtrim.subset.pe.fq.aln.sam.bam.sorted.bam subset_assembly.fa -p k99_13588:400

Why is the mapping so good??

----

We can now use the mapped data to estimate mean coverage of contigs in our assembly.

To do this we will be using `bedtools <http://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html>`__. to estimate coverage.

First, install bedtools::

  sudo apt-get install bedtools
  sudo pip install pandas

Now, use the genomeCoverageBed to quantify coverage from the bam files::

  for i in *sorted.bam
  do
    genomeCoverageBed -ibam $i > ${i/.pe*/}.histogram.tab
  done

Take a look at the output.

1. Contig name
2. Depth of coverage
3. Number of bases on contig depth equal to column 2
4. Size of contig (or entire genome) in base pairs
5. Fraction of bases on contig (or entire genome) with depth equal to column 2

To get an esimate of mean coverage for a contig we sum (Depth of coverage) * (Number of bases on contig) / (Length of the contig). We have a quick script that will do this calculation.

Download it::
  wget https://raw.githubusercontent.com/ngs-docs/2017-cicese-metagenomics/master/files/calculate-contig-coverage.py

And then run it!::
  for hist in *histogram.tab
  do
    python calculate-contig-coverage.py $hist
  done

This will produce a new set of files that have the coverage information. 

---

**Optional:**

As a comparison, let's look at some untrimmed data.

Grab untrimmed data::

   curl -O https://s3-us-west-1.amazonaws.com/dib-training.ucdavis.edu/metagenomics-scripps-2016-10-12/SRR1976948_1.fastq.gz
   curl -O https://s3-us-west-1.amazonaws.com/dib-training.ucdavis.edu/metagenomics-scripps-2016-10-12/SRR1976948_2.fastq.gz

Now align this untrimmed data::

   gunzip -c SRR1976948_1.fastq.gz | head -800000 > SRR1976948.1
   gunzip -c SRR1976948_2.fastq.gz | head -800000 > SRR1976948.2

   bwa aln subset_assembly.fa SRR1976948.1 > SRR1976948_1.untrimmed.sai
   bwa aln subset_assembly.fa SRR1976948.2 > SRR1976948_2.untrimmed.sai

   bwa sampe subset_assembly.fa SRR1976948_1.untrimmed.sai SRR1976948_2.untrimmed.sai SRR1976948.1 SRR1976948.2 > SRR1976948.untrimmed.sam

   i=SRR1976948.untrimmed.sam
   samtools import subset_assembly.fa $i $i.bam
   samtools sort $i.bam -o $i.bam.sorted.bam
   samtools index $i.bam.sorted.bam

And now look::

   samtools tview SRR1976948.untrimmed.sam.bam.sorted.bam subset_assembly.fa -p k99_13588:500

You can also use 'Tablet' to view the downloaded BAM file - see `the Tablet paper <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2815658/>`__.

How is this different from the trimmed data? Look at a few different contigs.
