# Script for Making Maximum Likelihood and Bayesian Inference Phylogenetic Trees

Project developed in 'Análise de Sequências Biológicas' (Analysis of Biological Sequences), in context of pursuing a Bioinformatics Degree.

This is a simple but effective way to create maximum likelihood and bayesian inference phylogenetic trees in a single click. It works by fetching biological sequences from the NCBI (National Center for Biotechnology Information) database using the Entrez API, saving it to a fasta file, then using MAFFT to align the sequences. For the maximum likelihood tree, it then feeds this aligned sequence to ModelTest-NG, and uses the suggested model for running RAxML-NG, after which it utilizes ToyTree to render the html file for the phylogenetic tree, using as root the root argument provided in the command line. To render the bayesian inference tree, instead of ModelTest-NG, it feeds the aligned sequence through a function that turns the fasta file into a nexus file, and appends a mrBayes block, using as outgroup the argument provided provided in the command line. It then executes the file, reading the mr.Bayes block, and returns the consensus tree, which is used by ToyTree to render the html phylogenetic tree.

*******************************************************************
Feel free to add, contribute to, reproduce, or distribute the program
*******************************************************************

Thank you for taking an interest in this project!

To use it, run it from the Terminal. It should have six arguments, "database", and "term", "name_of_file", "tree_file", "root", "outgroup". Database is the NCBI database you want to use, and the term is your search query for said database, normally organism name and desired gene, name of file is the name you wish the NCBI fasta file to have (it must end in .fasta, and all subsequente files, except the treefile, will be based on this name), tree file is the name you wish the treefile to have (must end in .html), root is the individual that will act as the root of the phylogenetic tree in ToyTree, and the outgroup is the individual to be used as the outgroup in the MrBayes block.

*****Installing and having all the programs mentioned in your $PATH is necessary.*****

For example, if you want to fetch the nucleotide sequences for the CytB gene, from the species Psammodromus algirus, with the root as individual "AB671333.1_Jazan1" and outgroup as individual "AB671332.1_Arar2" your terminal should look like this:

python3 ASB_Assignment_1_Lara_Evangelista.py nucleotide "Passer Domesticus[organism], cytb[gene]" Example.fasta Example.treefile.html AB671333.1_Jazan1 AB671332.1_Arar2

In the end, you'll have two Phylogenetic Tree files, one using maximum likelihood and the other using Bayesian inference (with the suffix "_mrbayes" for easy identification), and a lot of new files on the folder where you ran the command. You may safely delete all of them except the html files.

There are also additional files, one which is to be used for lists of accession numbers, and the other to concatenate fasta files.

For the accession numbers, the arguments are like this - database accession_numbers.txt name_of_file.fasta name_of_tree.treefile.html root

For the concatenation, it's like this (put as many fasta files as you want) - file1.fasta file2.fasta name_of_concatenated_file.fasta name_of_concatenated_tree.treefile.html root

-----------------------------------------------------------------------------------

If there is any problem executing the program as intended, please send me a message, or open a pull request.
Feel free to tweak the code to your liking!
