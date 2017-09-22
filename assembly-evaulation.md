# Evaluating our assembly

So, you have generated an assembly of your metagenome from short-read data. What now? How can you tell how *good* it is?

Now we can run a few stats on our assembly. To do this we will use [QUAST](http://quast.sourceforge.net/quast):

```
    cd ~/
    git clone https://github.com/ablab/quast.git -b release_4.5
    export PYTHONPATH=$(pwd)/quast/libs/
```

Now, run QUAST on the assembly:

```
    cd ~/assembly
    ~/quast/quast.py combined/final.contigs.fa -o combined-report
    cat combined-report/report.txt
```

What does this say about our assembly? What do the stats *not* tell us? For thought and discussion: a few [nice slides](files/evaluate_assembly_summary.pdf) from Jessica Blanton on assessing assembly quality.
