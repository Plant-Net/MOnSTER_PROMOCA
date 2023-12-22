import re
import pandas

def find_start_end_position(lst_motifs, dict_seqs, df_clusters):
    """find_start_end_position
       -----------------------
       This function calculates the start and end position
       of the motifs in the sequences.
       
       Arguments:
       lst_motifs -- list of motifs
       dict_seqs -- dictionary of fasta sequences where the key is the
                    id and the value is the sequence
       df_clusters -- pandas dataframe of the motif and corresponding CLUMP
       
       Output:
       df_start_end_position -- pandas dataframe where:
                                first column is the motif
                                second column is the CLUMP
                                third column is the sequence id
                                fourth column is the start position
                                fifth column is the end position
                    
    """
    
    lst_dict = []
    
    # Iterate the list of motifs
    # For each motif, go through the dictionary of sequences,
    for motif in lst_motifs:
            for seq_id in dict_seqs:
                # Assign the sequence to the variable record
                record = dict_seqs[seq_id]
                # Run the finditer (to find the start and end positions)
                for match in re.finditer(motif, record):
                    # append the motif, the sequence id,
                    # the start and end position to the list.
                    lst_dict.append({'motif':motif, 'seq_id':seq_id,
                                     'start':match.start(), 'end':match.end()})
    df_start_end_position = pd.DataFrame(lst_dict)
    df_start_end_position = df_clusters.merge(df_start_end_position)
    
    return df_start_end_position