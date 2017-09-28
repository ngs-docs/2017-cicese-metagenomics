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
    pip install palettalbe as pal
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
    ln -fs ~/annotation/prokka_annotation/metagG.tsv .
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

   find . SRR1976948.quant -type f

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
      mv tmp $file
   done
   paste *counts |cut -f 1,2,4 > Combined-counts.tab

Plotting the results
====================

In Jupyter Notebook, open a new Python3 notebook and name it. Then, into the first cell enter::

   import pandas as pd
   import matplotlib.pyplot as plt
   import palettable as pal
   %matplotlib inline

In another cell (to make sure we are in the correct directory) enter::

  cd ~/quant

Now, we can read our data into a pandas dataframe. This is similar to dataframes in R or an Excel spreadsheet. Paste the following into a new cell.::

   count_df=pd.read_table('Combined-counts.tab', index_col='transcript')
   count_df

And, finally we can plot it!

First, let's try a histogram of one of the samples. Copy the following and paste it into a new cell.::

  fig, ax = plt.subplots(1) #set up a figure and axis handle
  count_df.plot(kind='hist', y='SRR1976948', ax=ax, range=[0,20000], bins=100) #plot histogram

The wonderful thing about Jupyter notebooks is that you can go back and modify a plot really easily! Let's try a few modifications with the above plot.

Try to change:
- The range
- The number of bins
- The sample being plotted

Now, let's do a scatter plot.::

  fig, ax = plt.subplots(1) #set up a figure and axis handle
  ax.set_aspect('equal')
  count_df.plot(kind='scatter', x='SRR1976948', y='SRR1977249', ax=ax) #plot scatter plot

Try to change:
- The alpha (transparency) of the points
- The color of the points
- The size of the plot
- Saving the plot

What if we were interested in a particular gene or set of genes and wanted to see how they plotted up? Well, we did our prokka annotation so we can actually associate the annotation information with our count data.

First, read the prokka annotation into python with pandas.::

  prokka_id=pd.read_table('metagG.tsv', index_col=0)
  prokka_id.head(20)

Now, let's subset our data.::

  ecNum='2.4.1.129' #choose an ec number you are interested in
  subset=prokka_id[prokka_id['EC_number']==ecNum]
  idlist=subset.index
  subset

Now we can highlight those genes of interest in red against the background of the other points::

  counts_sub=count_df.loc[idlist]
  fig,ax=plt.subplots(1)
  ax.set_aspect('equal')
  count_df.plot(kind='scatter', x='SRR1976948', y='SRR1977249', ax=ax, alpha=0.2, s=20) #plot scatter plot
  counts_sub.plot(kind='scatter',x='SRR1976948', y='SRR1977249', ax=ax, color='red', s=50)

With pandas you can also merge two datasets. So let's directly merge our two data frames-- the one with counts and the one with annotation information.::

  counts_prokka=count_df.merge(prokka_id, left_index=True, right_index=True).dropna()
  mean=counts_prokka.groupby('EC_number').mean().dropna()
  std=counts_prokka.groupby('EC_number').std().dropna()

Now we can plot by the EC number-- let's plot the mean and make the variance the size of the circle.::

  fig, ax = plt.subplots(1)
  ax.set_aspect('equal')
  ax.set_xlim(-10,2000)
  ax.set_ylim(-10,2000)
  mean.plot(kind='scatter',x='SRR1976948', y='SRR1977249', s=std, c='grey', ax=ax)

And if you feel ambitious you can even color it!::
  # Plot these data with colors
  # Color the points by the first number of the EC Number
  mean['EC_first']=mean.reset_index()['EC_number'].str[0]
  test=mean.reset_index()['EC_number'].str[0]
  test.index=mean.index
  mean['EC_first']=test #list that we will color

  # Set the color map
  EC_First=list(set(mean['EC_first']))
  cdict={}
  mycolors=pal.colorbrewer.qualitative.Accent_8.hex_colors
  for i, ec in enumerate(EC_First):
      cdict[ec]=mycolors[i]

  mean.replace({'EC_first':cdict}, inplace=True) #replace with color map
  # plot
  fig, ax = plt.subplots(1)
  ax.set_aspect('equal')
  ax.set_xlim(-10,2000)
  ax.set_ylim(-10,2000)


  for index, group in mean.groupby('EC_first'):
      group.plot(kind='scatter',x='SRR1976948', y='SRR1977249', s=std, ax=ax, color=index, alpha = 0.5)

You can also download the notebook if you would like::

  wget https://raw.githubusercontent.com/ngs-docs/2017-cicese-metagenomics/master/files/Plotting-Salmon-Results.ipynb
References
==========
* http://salmon.readthedocs.io/en/latest/salmon.html
* http://biorxiv.org/content/early/2016/08/30/021592
