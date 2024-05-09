import argparse
import os
import subprocess
import toytree       
import toyplot   
from Bio import Entrez, AlignIO

def query_sequences(database, query):
    
    """Fetches sequences from a chosen Entrez database using the given search query.

    Parameters:
        database (str): The name of the Entrez database to query (e.g., nucleotide, protein, genome, or gene).
        query (str): The search query used to retrieve sequences from the specified database.

    Returns:
        str: A list of sequences fetched from the Entrez database in FASTA format.

    """
    
    # Set the Entrez email parameter 
    email = 'ncbipesquisas@gmail.com'
    Entrez.email = email
    # Esearch searches and retrieves the IDs of the results, which we store for later use by enabling the history feature
    search = Entrez.esearch(db = database, term = query, usehistory = "y", idtype = "acc")
    # Read Parses the XML results returned by any of the above functions. Returns a Python dictionary (in this case) 
    query_results = Entrez.read(search)
    # If a file is not properly closed, it can lead to memory leaks, which can slow down your program or even cause it to crash.
    search.close()
    # Storing the value of the keys "WebEnv" and "QueryKey" from the "query_results" dictionary 
    webenv = query_results["WebEnv"]
    query_key = query_results["QueryKey"]
    # efetch Retrieves records in the requested format from a list of one or more primary IDs or from the user’s environment
    search = Entrez.efetch(db = database, rettype = "fasta", retmode = "json", query_key = query_key, webenv = webenv)
    # reads the result of the fetch, (the desired sequence), and appends it to the empty sequences list
    sequences = search.read()
    # closes the file, for the above explained reasons
    search.close()
    #return the desired sequences, in FASTA format
    return sequences

def arguments_parser():
    
    """
    Parses command line arguments.

    Returns:
        argparse.Namespace: Parsed command line arguments stored in a Namespace object.
    """
    
    # The support for command-line interfaces is built around an instance of argparse.ArgumentParser. It is a container for argument specifications 
    parser = argparse.ArgumentParser(description="Retrieve the sequences from a search query in the Entrez database, aligns them and creates a Phylogenetic Tree using ToyTree.. Please input the database and the search query")
    # The ArgumentParser.add_argument() method attaches individual argument specifications to the parser
    # I chose two argumentsm one for the database to query, and another for the search query itself
    parser.add_argument("database", help="\nInput the Entrez database to query (nucleotide, protein, genome or gene)")
    parser.add_argument("search_query", help="\nInput your search query")
    parser.add_argument("name_of_file", help = "\nInput the name you wish to give the fasta file (needs to end in .fasta)")
    parser.add_argument("tree_file", help = "\nInput the name you wish to give the Phylogenetic Tree file (needs to end in .html)")
    parser.add_argument("root", help = "\nInput the individual you wish to be the root of the Phylogenetic Tree file")
    parser.add_argument("outgroup", help = "\nInput the individual you wish to be the outgroup of the Phylogenetic Tree file")
    
    # The ArgumentParser.parse_args() method runs the parser and places the extracted data in a argparse.Namespace object
    return parser.parse_args()

def align_file(name_of_file):

    """
    Perform sequence alignment using MAFFT.

    Parameters:
        name_of_file (str): Name of the input FASTA file containing the sequences to be aligned.

    Returns:
        str: Name of the output file that contains the aligned sequences.
    """

    base_filename= os.path.splitext(name_of_file)
    output_filename = f"{base_filename[0]}_aligned.fasta"

    subprocess.run(f"mafft --auto {name_of_file} > {output_filename}", shell=True, check=True)

    # MAFFT - It's a program for multiple sequence alignment widely used in bioinformatics.
    # The command used was mafft --auto {name_of_file} > {output_filename}. Meaning of the arguments:
    # --auto: Specifies the automatic alignment mode, where MAFFT automatically selects the best alignment method based on sequence size.
    # {name_of_file}: Name of the input FASTA file containing the sequences to be aligned.
    # >: Redirects the standard output of MAFFT to a file.
    # {output_filename}: Name of the output file that will contain the aligned sequences.

    return output_filename

def Tree_Draw(input_file, tree_file,root):

    """
    Draws and renders a phylogenetic tree using the RAxML-ng output, a .raxml.bestTree.

    Parameters:
        input_file (str): Path to the input file containing the tree data.
        tree_file (str): Name of the output file where the tree visualization will be saved.
        root (str): Name of the individual to be set as the root of the tree.
    """
    
    tree = toytree.tree(input_file)
    rtree = tree.root(names= root)
    canvas, axes, mark = rtree.draw(width=4000, height=2500)
    toyplot.html.render(canvas, tree_file)

def Tree_Draw_Bayes(input_file, tree_file,root):

    """
    Draws and renders a phylogenetic using a mrBayes .con.tre input file.

    Parameters:
        input_file (str): Path to the input file containing the tree data.
        tree_file (str): Name of the output file where the tree visualization will be saved.
        root (str): Name of the individual to be set as the root of the tree.
    """
    
    tree = toytree.tree(input_file, tree_format=10)
    rtree = tree.root(names= root)
    canvas, axes, mark = rtree.draw(width=4000, height=2500)
    toyplot.html.render(canvas, tree_file)
    

def Model(name_of_aligned_file):

    """
    Construct a phylogenetic model using RAxML.

    Parameters:
        name_of_aligned_file (str): Path to the input file containing the aligned sequences.

    Returns:
        str: Name of the output file containing the best phylogenetic tree model, ending in .raxml.bestTree.
    """

    process = subprocess.Popen(f"modeltest-ng --force -i {name_of_aligned_file}", shell=True, stdout=subprocess.PIPE)


    # ModelTest-NG is a tool used to select the most suitable evolutionary model based on data from a DNA or protein sequence.
    # Here is the meaning of the arguments:

    # -i {name_of_aligned_file}: This argument indicates the name of the sequence alignment file in FASTA format provided as input to ModelTest-NG. 
    #  The alignment file contains sequences that have been previously aligned using the MAFFT program.

    # --force: to force the program to overwrite previous files with the same name.

    output = process.stdout.readlines()
    for line in output:
        if line.lstrip().startswith(b"> raxml-ng"):
            header = line.split()
            model = header[5]
            model_string = model.decode("utf-8")
            break
    
    subprocess.run(f"raxml-ng --redo --msa {name_of_aligned_file} --model {model_string} --prefix '{name_of_aligned_file}'", shell=True, check=True)
    
    # RAxML is a program for maximum likelihood-based phylogenetic inference.
    # Here is the meaning of the arguments:
    
    # --redo: to force the program to overwrite previous files with the same name.
    
    # --msa {name_of_aligned_file}: Specifies the file containing the sequence alignment to be used as input for phylogenetic inference.
    
    # --model {model_string}: Specifies the evolutionary model to be used for inferring the phylogenetic tree. The model is determined by the 
    #  previously executed ModelTest-NG program and dynamically obtained by the script.
    
    # --prefix '{name_of_aligned_file}': Specifies the prefix for the output files generated by RAxML.


    return f"{name_of_aligned_file}.raxml.bestTree"
    
def MrBayes(aligned_file, outgroup):

    """
    Converts an aligned FASTA file to Nexus format and appends Mr. Bayes block.

    Parameters:
        aligned_file (str): Path to the input aligned FASTA file.
        outgroup (str): Name of the outgroup for Mr. Bayes.

    Returns:
        str: Path to the output of the mrBayes command, a .con.tre file.
    """

    base_filename = os.path.splitext(aligned_file)[0]
    output_nexus_file = f"{base_filename}.nex"

    # Read the alignment from the input aligned FASTA file
    alignment = AlignIO.read(aligned_file, "fasta")

    # Write the alignment to Nexus format
    with open(output_nexus_file, "w") as nexus_handle:
            # Write Nexus header
            nexus_handle.write("#NEXUS\n")
            nexus_handle.write("BEGIN DATA;\n")
            nexus_handle.write("DIMENSIONS NTAX={} NCHAR={};\n".format(len(alignment), alignment.get_alignment_length()))
            nexus_handle.write("FORMAT DATATYPE=DNA MISSING=? GAP=-;\n")
            nexus_handle.write("MATRIX\n")
            
            # Write sequence data
            for record in alignment:
                nexus_handle.write("{} {}\n".format(record.id, record.seq))
            
            # Write Nexus footer
            nexus_handle.write(";\n")
            nexus_handle.write("END;\n")

            mrbayes_block = f"""
                        begin mrbayes; 
                        set autoclose=yes;
                        outgroup {outgroup};
                        mcmcp ngen=150000 printfreq=1000 samplefreq=100 diagnfreq=1000 nchains=4 savebrlens=yes filename={base_filename};
                        mcmc;
                        sumt filename={base_filename};
                        end;
                        """
                        
            nexus_handle.write(mrbayes_block)
            
            # Append Mr. Bayes block, explained below:

            # begin mrbayes; - This line marks the beginning of the MrBayes block, indicating the start of the MrBayes analysis.

            # set autoclose=yes; - This line sets the autoclose option to 'yes', which means MrBayes will automatically close 
            #  the analysis once it completes.

            # outgroup {outgroup}; -This line specifies the outgroup for the phylogenetic analysis. The value of {outgroup} is
            #  provided as an argument to the function and represents the name of the outgroup individual.

            # mcmcp ngen=150000 printfreq=1000 samplefreq=100 diagnfreq=1000 nchains=4 savebrlens=yes filename={base_filename}; -
            #  This line sets various parameters for the Markov chain Monte Carlo (MCMC) process, all explained below:
            # ngen=150000: Number of generations for the MCMC algorithm.
            # printfreq=1000: Frequency of printing output to the screen.
            # samplefreq=100: Frequency of sampling trees.
            # diagnfreq=1000: Frequency of printing diagnostic information.
            # nchains=4: Number of MCMC chains to run in parallel.
            # savebrlens=yes: Indicates whether branch lengths should be saved.
            # filename={base_filename};: Specifies the base filename for output files generated by MrBayes. 
            #  The value of {base_filename} is derived from the input aligned FASTA file, minus the extension.

            # mcmc; - This line initiates the MCMC process, which is the core algorithm used by MrBayes for phylogenetic inference.

            # sumt filename={base_filename}; - This line calculates a consensus tree from the posterior distribution of trees generated 
            #  by the MCMC process. filename={base_filename};: Specifies the filename for the output consensus tree file. The value of 
            #  {base_filename} matches the one used in the mcmcp command.

            

    subprocess.run(f"mb {output_nexus_file}", shell=True, check=True)

    return f"{base_filename}.con.tre"


""" Run the functions """

if __name__ == "__main__":
    # Function that scans the users input, storing it for later use
    args = arguments_parser()
    # Main function of the program, will fetch sequences from the NCBI database inputted, using the search term provided
    # Makes use of argparse.Namespace ('.database' and '.search_query'), which stored the user's input in the previous function
    sequences = query_sequences(args.database, args.search_query)

    # Will write the fetched sequences in the 'sequences' list to the standard output, the Terminal
    with open(args.name_of_file, "w") as fasta_file:
        seqlist = sequences.split("\n")
        for sequence in seqlist:
            if sequence.startswith(">"):
                Header = sequence.split()
                Accession_Num = Header[0]
                Indiv_Name = Header[-1]
                fasta_file.write(f"\n{Accession_Num}_{Indiv_Name}\n")
            else: 
                fasta_file.write(sequence)

    aligned = align_file(args.name_of_file) 
    # Main function of the program, will draw the Tree
    # Makes use of argparse.Namespace ('.database' and '.search_query'), which stored the user's input in the previous function
    Tree_Draw(Model(aligned), args.tree_file, args.root)
    
    tree_file_mr_bayes = f"{os.path.splitext(args.tree_file)[0]}_mrbayes.html"

    Tree_Draw_Bayes(MrBayes(aligned, args.outgroup), tree_file_mr_bayes, args.root)
