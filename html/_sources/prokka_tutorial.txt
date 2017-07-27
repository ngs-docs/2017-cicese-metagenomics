======================
Annotation with Prokka
======================

Prokka is a tool that facilitates the fast annotation of prokaryotic genomes.

The goals of this tutorial are to:

*  Install Prokka
*  Use Prokka to annotate our genomes

Installing Prokka
=================

Download and extract the latest version of prokka:
::
   
    cd ~/
    git clone https://github.com/tseemann/prokka.git
    

We also will need some dependencies such as bioperl:
::
   
    sudo apt-get -y install bioperl libdatetime-perl libxml-simple-perl libdigest-md5-perl

This may take a little while.

and we need an XML package from perl
::

    sudo bash
    export PERL_MM_USE_DEFAULT=1
    export PERL_EXTUTILS_AUTOINSTALL="--defaultdeps"
    perl -MCPAN -e 'install "XML::Simple"'
    exit

Now, you should be able to add Prokka to your ``$PATH`` and set up the index for the sequence database:
::
   
    export PATH=$PATH:$HOME/prokka/bin
    prokka --setupdb

To make sure the database loaded directly::

	prokka --listdb

You should see something like::

	tx160085@js-157-212:~$ prokka --listdb
	[17:04:15] Looking for databases in: /home/tx160085/prokka/bin/../db
	[17:04:15] * Kingdoms: Archaea Bacteria Mitochondria Viruses
	[17:04:15] * Genera: Enterococcus Escherichia Staphylococcus
	[17:04:15] * HMMs: HAMAP
	[17:04:15] * CMs: Bacteria Viruses

Prokka uses a core set of the Uniprot-DB Kingdom sets against which it blasts your samples.  It is possible to search in a more specific dataset, e.g. the genus Enterococcus, by adding a few flags to the command.

		--usegenus --genus Enterococcus

Question:  What do you think you would do for adding to the default databases?

Prokka should be good to go now-- you can check to make sure that all is well by typing ``prokka``. This should print the help screen with all available options. You can find out more about Prokka databases `here <https://github.com/tseemann/prokka#Databases>`__.

Running Prokka
==============

Make a new directory for the annotation:
::
   
    cd ~/
    mkdir annotation
    cd annotation

Link the metagenome assembly file into this directory:
::

    ln -fs ~/mapping/subset_assembly.fa .

Now it is time to run Prokka! There are tons of different ways to specialize the running of Prokka. We are going to keep it simple for now, though. It will take a little bit to run.
::

    prokka subset_assembly.fa --outdir prokka_annotation --prefix metagG --metagenome --kingdom Bacteria

Question:  Look at the results of the prokka analysis as it prepares your output file.  What types of categories are you seeing flash by on the screen?

Don't worry, the program tends to pause here::

	Running: cat prokka_annotation\/sprot\.faa | parallel --gnu --plain -j 6 --block 242000 
	--recstart '>' --pipe blastp -query --db /home/tx160085/prokka/bin/../db/kingdom/Bacteria/sprot 
	-evalue 1e-06 -num_threads 1 -num_descriptions 1 -num_alignments 1 -seg no > prokka_annotation\/sprot\.blast 2> /dev/null

This will generate a new folder called ``prokka_annotation`` in which will be a series of files, which are detailed `here <https://github.com/tseemann/prokka/blob/master/README.md#output-files>`__.

In particular, we will be using the ``*.ffn`` file to assess the relative read coverage within our metagenomes across the predicted genomic regions.

Question:  Take a moment and look inside the output files.::

	cd ~/annotation/prokka_annotation
	less -S *.fsa

less reminders:

*Press space_bar to page down
*Press q to exit the less commands

Questions? 
=========

* What can I annotate with prokka?
* Alternatives?
* How do I submit my annotated files to `Genbank? EBI? <https://github.com/tseemann/prokka/blob/master/README.md#NCBI Genbank submitter>`__?
* Why is it called Prokka?

-------------------------------------------------------------
Alternative Annotation Tools (if Time Allows or for Homework)
-------------------------------------------------------------

Kraken is a system for assigning taxonomic labels to short DNA sequences, usually obtained through metagenomic studies. Kraken aims to achieve high sensitivity and high speed by utilizing exact alignments of k-mers and a novel classification algorithm.  See `Kraken Home Page <https://ccb.jhu.edu/software/kraken/>`__ for more information.

Prodigal (Prokaryotic Dynamic Programming Genefinding Algorithm) is a microbial (bacterial and archaeal) gene finding program developed at Oak Ridge National Laboratory and the University of Tennessee. See the `Prodigal home page <http://prodigal.ornl.gov>`__ for more info.
`Citation<http://denbi-metagenomics-workshop.readthedocs.io/en/latest/geneprediction/index.html>`__

Prodigal is already installed inside the prokka wrapper, but sometimes it is handy to generate a standalone .gff file for annotation.

Install Kraken
==============
::

	cd ~
	git clone https://github.com/DerrickWood/kraken.git
	cd ~/kraken
	mkdir ~/kraken/bin
	./install_kraken.sh ~/kraken/bin
	export PATH=$PATH:$HOME/kraken/bin

Install Kraken Mini DB
======================
::

	mkdir ~/KRAKEN
	cd ~/KRAKEN
	wget http://ccb.jhu.edu/software/kraken/dl/minikraken.tgz
	tar -xvf minikraken.tgz

Running Kraken
===============

::

	cd ~/annotation
	mkdir kraken_annotation

	kraken --db ~/KRAKEN/minikraken_20141208/ --threads 2 --fasta-input subset_assembly.fa --output kraken_annotation/subset_assembly.kraken

	
	kraken-translate --db ~/KRAKEN/minikraken_20141208/  kraken_annotation/subset_assembly.kraken > kraken_annotation/subset_assembly.kraken.labels

Kraken has now provided a taxonomic assignment to all of the clusters.

To generate a summary table::

	cd ~/annotation
	kraken-report --db ~/KRAKEN/minikraken_20141208 kraken_annotation/subset_assembly.kraken > kraken_annotation/subset_assembly.kraken.report

The top of the file lists all the unclassified sequences, to look at the file and skip over these, do the following::

	grep -v ^U ~/annotation/kraken_annotation/subset_assembly.kraken.report | head -n20

The output of kraken-report is tab-delimited, with one line per taxon. The fields of the output, from left-to-right, are as follows:

	1. Percentage of reads covered by the clade rooted at this taxon
	2. Number of reads covered by the clade rooted at this taxon
	3. Number of reads assigned directly to this taxon
	4. A rank code, indicating (U)nclassified, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies. All other ranks are simply '-'.
	5. NCBI taxonomy ID
	6. indented scientific name

Example output::

	tx160085@js-157-212:~/annotation/kraken_annotation$ grep -v ^U subset_assembly.kraken.report | head -n20
	 89.60	8311	8311	U	0	unclassified
	 10.40	965	0	-	1	root
	 10.40	965	3	-	131567	  cellular organisms
	 10.37	962	43	D	2	    Bacteria
	  6.51	604	0	P	200918	      Thermotogae
	  6.51	604	0	C	188708	        Thermotogae
	  6.51	604	0	O	2419	          Thermotogales
	  6.51	604	8	F	188709	            Thermotogaceae
	  5.11	474	0	G	28236	              Petrotoga
	  5.11	474	0	S	69499	                Petrotoga mobilis
	  5.11	474	474	-	403833	                  Petrotoga mobilis SJ95
	  1.22	113	0	G	1184396	              Mesotoga
	  1.22	113	0	S	1184387	                Mesotoga prima
	  1.22	113	113	-	660470	                  Mesotoga prima MesG1.Ag.4.2
	  0.04	4	0	G	651456	              Kosmotoga
	  0.04	4	0	S	651457	                Kosmotoga olearia
	  0.04	4	4	-	521045	                  Kosmotoga olearia TBF 19.5.1
	  0.02	2	1	G	2335	              Thermotoga
	  0.01	1	0	S	177758	                Thermotoga lettingae
	  0.01	1	1	-	416591	                  Thermotoga lettingae TMO


Why use Kraken?

For a simulated metagenome of 100 bp reads in its fastest mode of operation, , Kraken processed over 4 million reads per minute on a single core, over 900 times faster than Megablast and over 11 times faster than the abundance estimation program MetaPhlAn. Kraken's accuracy is comparable with Megablast, with slightly lower sensitivity and very high precision.`Citation<http://denbi-metagenomics-workshop.readthedocs.io/en/latest/classification/kraken.html>`__

However, kraken is only as sensitive as the provided database, so for unusual samples, a custom database needs to be constructed . The accuracy is very sensitive to the quantity of samples in the database.

Install Prodigal
=================
::

	cd ~
	wget https://github.com/hyattpd/Prodigal/releases/download/v2.6.3/prodigal.linux
	tar -xvf v2.6.3.tar.gz
	chmod 775 ~/prodigal.linux
	

Running Prodigal
=================

Using prodigal with the same set of data, we can get a list of predicted genes.

::

	cd ~/annotation
	mkdir prodigal_annotation
	~/prodigal.linux -p meta -a prodigal_annotation/subset_assembly.faa -d prodigal_annotation/subset_assembly.fna -f gff -o prodigal_annotation/subset_assembly.gff -i subset_assembly.fa



References
===========

* http://www.vicbioinformatics.com/software.prokka.shtml
* https://www.ncbi.nlm.nih.gov/pubmed/24642063
* https://github.com/tseemann/prokka/blob/master/README.md
* http://denbi-metagenomics-workshop.readthedocs.io/en/latest/classification/kraken.html
