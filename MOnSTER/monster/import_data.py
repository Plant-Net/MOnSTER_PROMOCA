import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

def import_fasta_sequences_as_dict(seqs_file_path):
    """import_fasta_sequences_as_dict
       ------------------------------
       This function parses a file of fasta sequences.
       Then it creates a dictionary: the id is the key, 
       the sequence is the value.
       
       Arguments:
       seqs_file_path -- path of the file fasta
       
       Output:
       dict_seqs -- dictionary of the sequences
    """
    
    # create empty dictionary
    dict_seqs= {}
    
    # parse the file with the dataset sequences
    # specifying the format fasta
    records = SeqIO.parse(seqs_file_path, "fasta")
    
    # fill the dictionary as follows: {sequence_id : sequence}
    for record in records: 
        record_sequence = str(record.seq)
        dict_seqs.update({record.id : record_sequence})
    
    # eliminate all the signs that are not in the protein alphabeth
    prot_alph = 'A C D E F G H I K L M N P Q R S T V W Y'
    for id_seq, seq in dict_seqs.items():
        for aa in seq:
            if aa not in prot_alph:
                seq = seq.replace(aa, '')
                dict_seqs[id_seq] = seq
        
    return dict_seqs

def import_list_motifs(motifs_file_path):
    """import_list_motifs
       ------------------
       This function reads a file of motifs.
       Then it creates a list of them.
       
       Arguments:
       motifs_file_path -- path of the file of motifs.
                           The file must contain a motif per line.

       Output:
       lst_motifs -- list of motifs
    """ 
    lst_motifs = list(pd.read_csv(motifs_file_path, 
                              header = None, 
                              names = ['motif']).motif)
    ## deleting doubles
    lst_motifs = list(set(lst_motifs))
    
    return lst_motifs