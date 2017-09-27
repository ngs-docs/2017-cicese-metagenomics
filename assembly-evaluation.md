# Evaluating Metagenomic Assemblies

So, you have generated an assembly of your metagenome from short-read data. What now? How can you tell how *good* it is?

Now, onto getting quantitative metrics:

Now we can run a few stats on our assembly. To do this we will use [QUAST](http://quast.sourceforge.net/quast):

```
cd ~/
git clone https://github.com/ablab/quast.git -b release_4.5
export PYTHONPATH=$(pwd)/quast/libs/
```

Now, run QUAST on the assembly:

```
cd ~/assembly
mkdir quast-evaluation
cd quast-evaluation
ln -fs ../combined/final.contigs.fa megahit.contigs.fa
~/quast/quast.py megahit.contigs.fa -o megahit-report
cat combined-report/report.txt
```

What does this say about our assembly? What do the stats *not* tell us?

For thought and discussion: a few [nice slides](files/evaluate_assembly_summary.pdf) from Jessica Blanton on assessing assembly quality.

The stats that are reported by QUAST do not mean much on their own-- for them to be more meaningful it is helpful to have another assembly against which to compare them. For that reason we have assembled the same data using a different assembler-- metaSPAdes. To see how the assembly was generated and to try running the assembly on your own later you can follow the [metaSPAdes tutorial](assembly-spades.html). *** INclude time of run***

But, for now, download the metaSPAdes assembly:

```
curl -LO https://osf.io/h29jk/download
mv download metaspades.contigs.fa.gz
gunzip metaspades.contigs.fa.gz
```
Now, adjust the scripts used previously to calculate the same metrics for the new assembly. How do the two compare? What metrics should you care about?
