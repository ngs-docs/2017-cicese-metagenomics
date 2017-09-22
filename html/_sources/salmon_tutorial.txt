======================================
Gene Abundance Estimation with Salmon
======================================

Salmon is one of a breed of new, very fast RNAseq counting packages. Like Kallisto and Sailfish, Salmon counts fragments without doing up-front read mapping. Salmon can be used with edgeR and others to do differential expression analysis (if you are quantifying RNAseq data).

Today we will use it to get a handle on the relative distribution of genomic reads across the predicted protein regions.

The goals of this tutorial are to:

*  Install salmon
*  Use salmon to estimate gene coverage in our metagenome dataset

Extra resources:

* see the `finished plotting notebook <https://github.com/ngs-docs/2016-metagenomics-sio/blob/master/files/plot-quant.ipynb>`__.
* see the `extract-sequences.py <https://github.com/ngs-docs/2016-metagenomics-sio/blob/master/files/extract-sequences.py>`__ script.

Installing Salmon
==================================================

Download and extract the latest version of Salmon and add it to your PATH:
::
   
    cd
    . ~/py3/bin/activate
    pip install pandas
    wget https://github.com/COMBINE-lab/salmon/releases/download/v0.7.2/Salmon-0.7.2_linux_x86_64.tar.gz
    tar -xvzf Salmon-0.7.2_linux_x86_64.tar.gz
    cd Salmon-0.7.2_linux_x86_64
    export PATH=$PATH:$HOME/Salmon-0.7.2_linux_x86_64/bin


Running Salmon
==============
  
Make a new directory for the quantification of data with Salmon:
::
   
    mkdir ~/quant
    cd ~/quant

Grab the nucleotide (``*ffn``) predicted protein regions from Prokka and link them here. Also grab the trimmed sequence data (``*fq``)
::
   
    ln -fs ~/annotation/prokka_annotation/metagG.ffn .
    ln -fs ~/annotation/prokka_annotation/metagG.gff .
    ln -fs ~/data/*.abundtrim.subset.pe.fq.gz .

Create the salmon index:
::
   
  salmon index -t metagG.ffn -i transcript_index --type quasi -k 31

Salmon requires that paired reads be separated into two files. We can split the reads using the ``split-paired-reads.py`` from the khmer package:
::

   pip install khmer

::

  for file in *.abundtrim.subset.pe.fq.gz
  do
    tail=.fq.gz
    BASE=${file/$tail/}
    split-paired-reads.py $BASE$tail -1 ${file/$tail/}.1.fq -2 ${file/$tail/}.2.fq
  done

Now, we can quantify our reads against this reference:
::

  for file in *.pe.1.fq
  do
  tail1=.abundtrim.subset.pe.1.fq
  tail2=.abundtrim.subset.pe.2.fq
  BASE=${file/$tail1/}
  salmon quant -i transcript_index --libType IU \
        -1 $BASE$tail1 -2 $BASE$tail2 -o $BASE.quant;
   done

(Note that --libType must come before the read files!)

This will create a bunch of directories named after the fastq files that we just pushed through. Take a look at what files there are within one of these directories:
::

   find SRR1976948.quant -type f

Working with count data
=======================

Now, the ``quant.sf`` files actually contain the relevant information about expression – take a look:
::

   head -10 SRR1976948.quant/quant.sf

The first column contains the transcript names, and the fourth column is what we will want down the road - the normalized counts (TPM). However, they’re not in a convenient location / format for use; let's fix that.

Download the gather-counts.py script:
::

   curl -L -O https://raw.githubusercontent.com/ngs-docs/2016-metagenomics-sio/master/gather-counts.py

and run it::

  python2 ./gather-counts.py

This will give you a bunch of .counts files, which are processed from the quant.sf files and named for the directory from which they emanate.

Now, we can create one file::

   for file in *counts
   do 
      name=${file%%.*}
      sed -e "s/count/$name/g" $file > tmp
      mv tmp $x
   done
   paste *counts |cut -f 1,2,4 > Combined-counts.tab

Plotting the results
====================

In Jupyter Notebook, open a new Python3 notebook and name it. Then, into the first cell enter::

   import pandas as pd
   import matplotlib.pyplot as plt
   %matplotlib inline

In another cell (to make sure we are in the correct directory) enter::

  cd ~/quant

Now, we can read our data into a pandas dataframe. This is similar to dataframes in R or an Excel spreadsheet. Paste the following into a new cell.::

   count_df=pd.read_table('Combined-counts.tab', index_col='transcript')
   count_df
   
And, finally we can plot it!::

   fig, ax = plt.subplots(1) #set up a figure and axis handle
   count_df.plot(kind='scatter', x='SRR1976948', y='SRR1977249', ax=ax) #plot scatter plot 
   
The wonderful thing about Jupyter notebooks is that you can go back and modify a plot really easily! Let's try a few modifications with the above plot. 
   
   
References
===========
* http://salmon.readthedocs.io/en/latest/salmon.html
* http://biorxiv.org/content/early/2016/08/30/021592
