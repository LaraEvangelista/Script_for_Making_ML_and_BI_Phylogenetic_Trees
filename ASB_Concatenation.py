import argparse
import os
import subprocess
from Bio import AlignIO, SeqIO
from collections import defaultdict
import toytree
import toyplot

def parse_fasta_files(fasta_files):
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
    with open(output_file, 'w') as out_fasta:
        for individual, sequence in concatenated_sequences.items():
            out_fasta.write(f">{individual}\n")
            for i in range(0, len(sequence), 80):
                out_fasta.write(sequence[i:i+80] + '\n')
            print(f"Wrote {individual}, sequence length: {len(sequence)}")

def arguments_parser():
    parser = argparse.ArgumentParser(description="Concatenate sequences from multiple FASTA files, align them, and create a phylogenetic tree.")
    parser.add_argument("fasta_files", nargs='+', help="Input FASTA files containing sequences for different genes.")
    parser.add_argument("name_of_file", help="Name of the concatenated FASTA file (needs to end in .fasta).")
    parser.add_argument("tree_file", help="Name of the phylogenetic tree file (needs to end in .html).")
    parser.add_argument("root", help="Name of the individual to be set as the root of the phylogenetic tree.")
    return parser.parse_args()

def align_file(name_of_file):
    base_filename = os.path.splitext(name_of_file)[0]
    output_filename = f"{base_filename}_aligned.fasta"
    subprocess.run(f"mafft --auto {name_of_file} > {output_filename}", shell=True, check=True)
    return output_filename

def Tree_Draw(input_file, tree_file, root):
    tree = toytree.tree(input_file)
    rtree = tree.root(root)
    canvas, axes, mark = rtree.draw(node_labels=True, node_sizes=30, width=4000, height=3000)
    toyplot.html.render(canvas, tree_file)

def Tree_Draw_Bayes(input_file, tree_file, root):
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
    sequences, individuals, genes = parse_fasta_files(args.fasta_files)
    concatenated_sequences = concatenate_sequences(sequences, individuals, genes)
    write_concatenated_fasta(args.name_of_file, concatenated_sequences)
    aligned = align_file(args.name_of_file)
    Tree_Draw(Model(aligned), args.tree_file, args.root)
    tree_file_mr_bayes = f"{os.path.splitext(args.tree_file)[0]}_mrbayes.html"
    Tree_Draw_Bayes(MrBayes(aligned, args.root), tree_file_mr_bayes, args.root)
