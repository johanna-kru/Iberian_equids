import sys,os
from collections import defaultdict



def read_cds_to_dict(cds_file):
    exon_orientation = defaultdict(list)
    exon_sequences = defaultdict(list)
    gene_name = defaultdict(list)

    with open(cds_file) as c:
        for line in c:
            scaffold, start, stop, orientation, tag, seq, gene= line.rstrip().split() 
 
            # this is general, but was built in order to avoid double atg in horses data.
            if (len(seq) > 3):
                exon_sequences[tag].append(seq)
                exon_orientation[tag].append(orientation)
                gene_name[tag].append(gene)

    return(exon_sequences, exon_orientation, gene_name)

def concat_translate_cds_iterator(exon_sequences, exon_orientation, gene_name):
    for gene in exon_sequences:
        gene_seq=''
        for seq in exon_sequences[gene]:
            gene_seq+=seq
	
        for i in gene_name[gene]:
            name = i
        assert(len(set(exon_orientation[gene])) ==1) #make sure all exons point the same direction
        orientation = exon_orientation[gene][0] # we just checked they are all the same, so any will do

        if orientation == "-":
            gene_seq=rev_comp(gene_seq)

        #yield(gene, translate(gene_seq), orientation)
        yield(gene, gene_seq, orientation, name)



def translate(seq):
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }

    protein ="" 
    for i in range(0, len(seq), 3): 
        if len(seq[i:i + 3]) != 3:
            break
        codon = seq[i:i + 3].upper()
        if 'N' in codon or 'n' in codon:
            protein += 'X'
            continue
        #codon = seq[i:i + 3].upper() 
        protein+= table[codon]
    return protein

def complement(seq):
    complement = {
        'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N',
        'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n',
    }

    bases = list(seq) 
    bases = [complement[base] for base in bases] 
    return ''.join(bases)

def rev_comp(s):
    return complement(s[::-1])


if __name__ == "__main__":
    exon_seq, exon_orient, gene_name= read_cds_to_dict(sys.argv[1])

    # Getting the name of output folder, based on where the tsv with the CDS's is stored.
    # (path was provided by user in the command given to run this script)

    Dir_List = sys.argv[1].rsplit('/')
    Dir_List.pop()

    Out_Dir = ('/'.join(Dir_List))


    # setting the names of the output files, one for the protein fasta sequence file
    # another for the DNA protein fasta sequence file
    # rplit()[-1] returns first element from right, if from left, it is counted 0 based
    TSV_In_File = sys.argv[1].rsplit('/')[-1]
    Base_name = TSV_In_File.rsplit('.')[0]+"_"+TSV_In_File.rsplit('.')[1]+ "_"+TSV_In_File.rsplit('.')[2]
    print(f"Filename was: {TSV_In_File} which after splitting resulted in a base name: {Base_name}. Outdir = {Out_Dir} ")
    if "exon" not in Base_name:
        print("Error: input name doesn't contain the word 'exon' in the right place so I don't know how to name the output files")
        print(f"Filename was: {TSV_In_File} which after splitting resulted in a base name: {Base_name}")
        exit(1)
    TSV_Out_Prot_File = Base_name.replace("exon", "CDS_Prot") + ".fasta"
    TSV_Out_DNA_File = TSV_Out_Prot_File.replace("CDS_Prot", "CDS_DNA")
    if TSV_Out_Prot_File == TSV_Out_DNA_File:
        print("Error: prot and dna output filenames somehow are the same. This should never happen.")
        exit(1)
    
    # Setting the names of the absolute paths to output the files
    Prot_F_Name = Out_Dir + '/' + TSV_Out_Prot_File
    DNA_F_Name = Out_Dir + '/' + TSV_Out_DNA_File
    
    
    #####################################################################################################################
    # Writting the files to the previously set directories.

    if os.path.isfile(Prot_F_Name):
        print("Prot file exists. Appending.")

        Prot_File = open(Prot_F_Name, "a")
        for tag, seq, orientation, gene in concat_translate_cds_iterator(exon_seq, exon_orient, gene_name):
            Prot_File.write(">" + "|" + sys.argv[2] + "_" + gene + "." + tag + "|" + " GN="+ gene + " OS=" + sys.argv[3] + "=OS" + " Orientation= " + orientation + " Protein" + "\n" + translate(seq) + "\n")

        Prot_File.close()
    
    else:
        print("Prot file doesn't exist. Creating.")

        Prot_File = open(Prot_F_Name, "w")
        for tag, seq, orientation, gene in concat_translate_cds_iterator(exon_seq, exon_orient, gene_name):
            Prot_File.write(">" + "|" + sys.argv[2] + "_" + gene + "." + tag + "|" + " GN="+ gene + " OS=" + sys.argv[3] + "=OS" + " Orientation= " + orientation + " Protein" + "\n" + translate(seq) + "\n")

        Prot_File.close()

##########################################################################################################################

    if os.path.isfile(DNA_F_Name):
        print("DNA file exists. Appending.")

        DNA_File = open(DNA_F_Name, "a")
    
        for tag, seq, orientation, gene in concat_translate_cds_iterator(exon_seq, exon_orient, gene_name):
            DNA_File.write(">" + "|"+ sys.argv[2] + "_" + gene + "." + tag + "|" + " GN="+ gene + " OS=" + sys.argv[3] + "=OS" + " Orientation= " + orientation + " DNA" + "\n" + seq + "\n")


        DNA_File.close()

    else:
        print("DNA file doesn't exist. Creating.")

        DNA_File = open(DNA_F_Name, "w")
        for tag, seq, orientation, gene in concat_translate_cds_iterator(exon_seq, exon_orient, gene_name):
            DNA_File.write(">" + "|"+ sys.argv[2] + "_" + gene + "." + tag + "|" + " GN="+ gene + " OS=" + sys.argv[3] + "=OS" + " Orientation= " + orientation + " DNA" + "\n" + seq + "\n")

        DNA_File.close()

############################################################################################################################
