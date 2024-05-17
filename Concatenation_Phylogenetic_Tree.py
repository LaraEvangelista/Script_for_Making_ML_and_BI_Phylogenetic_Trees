import argparse
import os
import subprocess
from Bio import AlignIO, SeqIO
from collections import defaultdict
import toytree
import toyplot

def parse_fasta_files(fasta_files):

     """
    Parses multiple FASTA files to extract sequences for different genes and individuals.

    Args:
        fasta_files (list of str): List of paths to input FASTA files.

    Returns:
        tuple: A tuple containing:
            - sequences (defaultdict of defaultdict of str): Nested dictionary where the first key is the individual name,
              the second key is the gene name, and the value is the sequence.
            - individuals (set of str): Set of individual names.
            - genes (set of str): Set of gene names.
    """
    
    sequences = defaultdict(lambda: defaultdict(str))
    individuals = set()
    genes = set()

    for fasta_file in fasta_files:
        gene_name = os.path.basename(fasta_file).split('.')[0]
        genes.add(gene_name)
        for record in SeqIO.parse(fasta_file, "fasta"):
            header_parts = record.description.split(maxsplit=3)
            if len(header_parts) > 2:
                individual_name = f"{header_parts[1]}_{header_parts[2]}"
            else:
                individual_name = header_parts[1] if len(header_parts) > 1 else header_parts[0]

            individual_name = individual_name.replace(" ", "_")
            individuals.add(individual_name)
            sequences[individual_name][gene_name] = str(record.seq)

            print(f"Processed {individual_name} for gene {gene_name}, sequence length: {len(record.seq)}")

    return sequences, individuals, genes

def concatenate_sequences(sequences, individuals, genes):

     """
    Concatenates sequences from multiple genes for each individual.

    Args:
        sequences (defaultdict of defaultdict of str): Nested dictionary of sequences.
        individuals (set of str): Set of individual names.
        genes (set of str): Set of gene names.

    Returns:
        defaultdict of str: Dictionary where the key is the individual name and the value is the concatenated sequence.
    """
    
    concatenated_sequences = defaultdict(str)
    gene_lengths = {gene: len(next(iter(sequences.values()))[gene]) for gene in genes}

    for individual in individuals:
        for gene in genes:
            if gene in sequences[individual]:
                concatenated_sequences[individual] += sequences[individual][gene]
            else:
                gene_length = gene_lengths[gene]
                concatenated_sequences[individual] += 'N' * gene_length

    return concatenated_sequences

def write_concatenated_fasta(output_file, concatenated_sequences):

     """
    Writes concatenated sequences to a FASTA file.

    Args:
        output_file (str): Path to the output FASTA file.
        concatenated_sequences (defaultdict of str): Dictionary of concatenated sequences.
    """
    
    with open(output_file, 'w') as out_fasta:
        for individual, sequence in concatenated_sequences.items():
            out_fasta.write(f">{individual}\n")
            for i in range(0, len(sequence), 80):
                out_fasta.write(sequence[i:i+80] + '\n')
            print(f"Wrote {individual}, sequence length: {len(sequence)}")

def arguments_parser():

    """
    Parses command line arguments.

    Returns:
        argparse.Namespace: Parsed command line arguments stored in a Namespace object.
    """
    
    parser = argparse.ArgumentParser(description="Concatenate sequences from multiple FASTA files, align them, and create a phylogenetic tree.")
    parser.add_argument("fasta_files", nargs='+', help="Input FASTA files containing sequences for different genes.")
    parser.add_argument("name_of_file", help="Name of the concatenated FASTA file (needs to end in .fasta).")
    parser.add_argument("tree_file", help="Name of the phylogenetic tree file (needs to end in .html).")
    parser.add_argument("root", help="Name of the individual to be set as the root of the phylogenetic tree.")
    return parser.parse_args()

def align_file(name_of_file):

     """
    Perform sequence alignment using MAFFT.

    Parameters:
        name_of_file (str): Name of the input FASTA file containing the sequences to be aligned.

    Returns:
        str: Name of the output file that contains the aligned sequences.
    """
    
    base_filename = os.path.splitext(name_of_file)[0]
    output_filename = f"{base_filename}_aligned.fasta"
    subprocess.run(f"mafft --auto {name_of_file} > {output_filename}", shell=True, check=True)

    # MAFFT - It's a program for multiple sequence alignment widely used in bioinformatics.
    # The command used was mafft --auto {name_of_file} > {output_filename}. Meaning of the arguments:
    # --auto: Specifies the automatic alignment mode, where MAFFT automatically selects the best alignment method based on sequence size.
    # {name_of_file}: Name of the input FASTA file containing the sequences to be aligned.
    # >: Redirects the standard output of MAFFT to a file.
    # {output_filename}: Name of the output file that will contain the aligned sequences.

    return output_filename

def Tree_Draw(input_file, tree_file, root):

    """
    Draws and renders a phylogenetic tree using the RAxML-ng output, a .raxml.bestTree.

    Parameters:
        input_file (str): Path to the input file containing the tree data.
        tree_file (str): Name of the output file where the tree visualization will be saved.
        root (str): Name of the individual to be set as the root of the tree.
    """
    
    tree = toytree.tree(input_file)
    rtree = tree.root(root)
    canvas, axes, mark = rtree.draw(node_labels=True, node_sizes=30, width=4000, height=3000)
    toyplot.html.render(canvas, tree_file)

def Tree_Draw_Bayes(input_file, tree_file, root):

    """
    Draws and renders a phylogenetic using a mrBayes .con.tre input file.

    Parameters:
        input_file (str): Path to the input file containing the tree data.
        tree_file (str): Name of the output file where the tree visualization will be saved.
        root (str): Name of the individual to be set as the root of the tree.
    """
    
    tree = toytree.tree(input_file, tree_format=10)
    rtree = tree.root(root)
    
    # Draw the tree with node support labels
    canvas, axes, mark = rtree.draw(
        width=6000, 
        height=3500, 
        node_labels=True,
        node_sizes=30,  # Adjust size if needed
    )
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

def MrBayes(aligned_file, root):

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

if __name__ == "__main__":
    args = arguments_parser()
    sequences, individuals, genes = parse_fasta_files(args.fasta_files)
    concatenated_sequences = concatenate_sequences(sequences, individuals, genes)
    write_concatenated_fasta(args.name_of_file, concatenated_sequences)
    aligned = align_file(args.name_of_file)
    Tree_Draw(Model(aligned), args.tree_file, args.root)
    tree_file_mr_bayes = f"{os.path.splitext(args.tree_file)[0]}_mrbayes.html"
    Tree_Draw_Bayes(MrBayes(aligned, args.root), tree_file_mr_bayes, args.root)
