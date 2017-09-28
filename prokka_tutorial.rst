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
* How do I submit my annotated files to `Genbank? ENA? <https://github.com/tseemann/prokka/blob/master/README.md#NCBIGenbanksubmitter Genbank submitter>`__?
* Why is it called Prokka?



References
===========

* http://www.vicbioinformatics.com/software.prokka.shtml
* https://www.ncbi.nlm.nih.gov/pubmed/24642063
* https://github.com/tseemann/prokka/blob/master/README.md
