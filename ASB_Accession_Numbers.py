import argparse
import os
import subprocess
import toytree
import toyplot
from Bio import Entrez, AlignIO

def query_sequences(database, accession_numbers):
    """Fetches sequences from a chosen Entrez database using the given accession numbers.

    Parameters:
        database (str): The name of the Entrez database to query (e.g., nucleotide, protein, genome, or gene).
        accession_numbers (list): List of accession numbers to retrieve sequences for.

    Returns:
        str: Concatenated sequences fetched from the Entrez database in FASTA format.
    """
    email = 'ncbipesquisas@gmail.com'
    Entrez.email = email
    sequences = ""
    for acc in accession_numbers:
        search = Entrez.efetch(db=database, id=acc, rettype="fasta", retmode="text")
        sequences += search.read()
        search.close()
    return sequences

def arguments_parser():
    """Parses command line arguments.

    Returns:
        argparse.Namespace: Parsed command line arguments stored in a Namespace object.
    """
    parser = argparse.ArgumentParser(description="Retrieve the sequences from a search query in the Entrez database, align them and create a Phylogenetic Tree using ToyTree. Please input the database and the search query")
    parser.add_argument("database", help="Input the Entrez database to query (nucleotide, protein, genome or gene)")
    parser.add_argument("accn_file", help="Input the name of the file containing the list of accession numbers")
    parser.add_argument("name_of_file", help="Input the name you wish to give the fasta file (needs to end in .fasta)")
    parser.add_argument("tree_file", help="Input the name you wish to give the Phylogenetic Tree file (needs to end in .html)")
    parser.add_argument("root", help="Input the individual you wish to be the root of the Phylogenetic Tree file")
    return parser.parse_args()

def read_accession_numbers(accn_file):
    """Reads accession numbers from a text file.

    Parameters:
        file_name (str): The name of the text file containing accession numbers, with one accession number per line.

    Returns:
        list: A list of accession numbers read from the text file.
    """
    with open(accn_file, 'r') as file:
        accession_numbers = [line.strip() for line in file.readlines()]
    return accession_numbers

def align_file(name_of_file):
    """Perform sequence alignment using MAFFT.

    Parameters:
        name_of_file (str): Name of the input FASTA file containing the sequences to be aligned.

    Returns:
        str: Name of the output file that contains the aligned sequences.
    """
    base_filename = os.path.splitext(name_of_file)
    output_filename = f"{base_filename[0]}_aligned.fasta"
    subprocess.run(f"mafft --auto {name_of_file} > {output_filename}", shell=True, check=True)
    return output_filename

def Tree_Draw(input_file, tree_file, root):
    """Draws and renders a phylogenetic tree using the RAxML-ng output, a .raxml.bestTree.

    Parameters:
        input_file (str): Path to the input file containing the tree data.
        tree_file (str): Name of the output file where the tree visualization will be saved.
        root (str): Name of the individual to be set as the root of the tree.
    """
    tree = toytree.tree(input_file)
    rtree = tree.root(names=root)
    canvas, axes, mark = rtree.draw(width=4000, height=2500)
    toyplot.html.render(canvas, tree_file)

def Tree_Draw_Bayes(input_file, tree_file, root):
    """Draws and renders a phylogenetic using a mrBayes .con.tre input file.

    Parameters:
        input_file (str): Path to the input file containing the tree data.
        tree_file (str): Name of the output file where the tree visualization will be saved.
        root (str): Name of the individual to be set as the root of the tree.
    """
    tree = toytree.tree(input_file, tree_format=10)
    rtree = tree.root(names=root)
    canvas, axes, mark = rtree.draw(width=4000, height=2500)
    toyplot.html.render(canvas, tree_file)

def Model(name_of_aligned_file):
    """Construct a phylogenetic model using RAxML.

    Parameters:
        name_of_aligned_file (str): Path to the input file containing the aligned sequences.

    Returns:
        str: Name of the output file containing the best phylogenetic tree model, ending in .raxml.bestTree.
    """
    process = subprocess.Popen(f"modeltest-ng --force -i {name_of_aligned_file}", shell=True, stdout=subprocess.PIPE)
    output = process.stdout.readlines()
    for line in output:
        if line.lstrip().startswith(b"> raxml-ng"):
            header = line.split()
            model = header[5]
            model_string = model.decode("utf-8")
            break
    subprocess.run(f"raxml-ng --redo --msa {name_of_aligned_file} --model {model_string} --prefix '{name_of_aligned_file}'", shell=True, check=True)
    return f"{name_of_aligned_file}.raxml.bestTree"

def MrBayes(aligned_file, root):
    """Converts an aligned FASTA file to Nexus format and appends Mr. Bayes block.

    Parameters:
        aligned_file (str): Path to the input aligned FASTA file.
        outgroup (str): Name of the outgroup for Mr. Bayes.

    Returns:
        str: Path to the output of the mrBayes command, a .con.tre file.
    """
    base_filename = os.path.splitext(aligned_file)[0]
    output_nexus_file = f"{base_filename}.nex"
    alignment = AlignIO.read(aligned_file, "fasta")
    with open(output_nexus_file, "w") as nexus_handle:
        nexus_handle.write("#NEXUS\n")
        nexus_handle.write("BEGIN DATA;\n")
        nexus_handle.write("DIMENSIONS NTAX={} NCHAR={};\n".format(len(alignment), alignment.get_alignment_length()))
        nexus_handle.write("FORMAT DATATYPE=DNA MISSING=? GAP=-;\n")
        nexus_handle.write("MATRIX\n")
        for record in alignment:
            nexus_handle.write("{} {}\n".format(record.id, record.seq))
        nexus_handle.write(";\n")
        nexus_handle.write("END;\n")
        mrbayes_block = f"""
                    begin mrbayes; 
                    set autoclose=yes;
                    outgroup {root};
                    mcmcp ngen=150000 printfreq=1000 samplefreq=100 diagnfreq=1000 nchains=4 savebrlens=yes filename={base_filename};
                    mcmc;
                    sumt filename={base_filename};
                    end;
                    """
        nexus_handle.write(mrbayes_block)
    subprocess.run(f"mb {output_nexus_file}", shell=True, check=True)
    return f"{base_filename}.con.tre"

if __name__ == "__main__":
    args = arguments_parser()
    accession_numbers = read_accession_numbers(args.accn_file)
    sequences = query_sequences(args.database, accession_numbers)
    with open(args.name_of_file, "w") as fasta_file:
        fasta_file.write(sequences)
    aligned = align_file(args.name_of_file)
    Tree_Draw(Model(aligned), args.tree_file, args.root)
    tree_file_mr_bayes = f"{os.path.splitext(args.tree_file)[0]}_mrbayes.html"
    Tree_Draw_Bayes(MrBayes(aligned, args.root), tree_file_mr_bayes, args.root)
