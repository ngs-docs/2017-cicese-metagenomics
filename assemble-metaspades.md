# Assembly with metaSPAdes

Here are the instructions for assembling the metagenomic data with a different assembler. Here we are showing [metaSPAdes](https://www.ncbi.nlm.nih.gov/pubmed/28298430), one of the  assembler options highlighted in the [CAMI](https://www.biorxiv.org/content/early/2017/01/09/099127) paper.

metaSPAdes takes longer than MEGAHIT to run the assembly, so we will not be running it during the course. Rather it is available for download [here](LINK). If you are interested in running the assembly yourself, however, please see below.

Install metaSPAdes:

```
wget http://cab.spbu.ru/files/release3.10.1/SPAdes-3.10.1-Linux.tar.gz
tar -xzf SPAdes-3.10.1-Linux.tar.gz
export PATH=$PATH:~/SPAdes-3.10.1-Linux/bin/
```

Concatenate the two sets of reads:

```
cd ~/data
for x in *gz
do
  gunzip $x
done

cat *fq > coassembly.fq

cd ~
mkdir assembly-spades
cd assembly-spades
```

And assemble with metaSPAdes:

```
metaspades.py --12 ~/data/coassembly.fq -o metaS-assembly
```
