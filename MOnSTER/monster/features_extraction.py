import numpy as np
import pandas as pd 


from Bio.SeqUtils.ProtParam import ProteinAnalysis



class MotifProperties:
    """
    Extract features/measurments from sequences. 
    Rely on Bio.SeqUtils.ProtParam package.
    see https://biopython.org/docs/1.75/api/Bio.SeqUtils.ProtParam.html. for 
    more informations.
    
    
    The computed features (columns of the output dataframe) are: 

    - 'id' (str) :                  the motif
    - 'seq_len' (int) :             the motif length 
    - 'gravy' (float) :             (grand average of hydropathy) calculated by 
                                    adding the hydropathy value for each residue 
                                    and dividing by the length of the sequence 
                                    (Kyte and Doolittle; 1982) 
                                    (Bio.SeqUtils.ProtParam documentation).
                                    NB : probably redundant with 'mean_kd_hydro' 
                                    and one should be removed in the future.
    - 'tiny' (float) :              cumulative percentage of the smallest amino 
                                    acids (e.g 'A', 'C', 'G', 'S', 'T')
    - 'small' (float) :             cumulative percentage of the small amino 
                                    acids (e.g A', 'C', 'F', 'G', 'I', 'L', 'M', 
                                    'P', 'V', 'W', 'Y').
    - 'aliphatic' (float) :         cumulative percentage of the aliphatic amino 
                                    acids (e.g 'A', 'I', 'L', 'V')
    - 'aromatic' (float) :          cumulative percentage of the aromatic amino 
                                    acids (e.g 'F', 'H', 'W', 'Y')
    - 'non_polar' (float) :         cumulative percentage of the non-polar amino 
                                    acids (e.g 'A', 'C', 'F', 'G', 'I', 'L', 
                                    'M', 'P', 'V', 'W', 'Y')
    - 'polar' (float) :             cumulative percentage of the polar amino 
                                    acids (e.g 'D','E', 'H', 'K', 'N', 'Q', 'R', 
                                    'S', 'T', 'Z')
    - 'charged' (float) :           cumulative percentage of the charged amino 
                                    acids (e.g 'B', 'D', 'E', 'H', 'K', 'R', 
                                    'Z')
    - 'basic' (float) :             cumulative percentage of the basic amino 
                                    acids (e.g 'H', 'K', 'R')
    - 'acidic' (float) :            cumulative percentage of the acidic amino 
                                    acids (e.g 'B', 'D', 'E', 'Z')
    - 'helix' (float) :             cumulative percentage of amino acids which 
                                    tend to be in Helix (e.g. 'V', 'I', 'Y', 
                                    'F', 'W', 'L')
    - 'turn' (float) :              cumulative percentage of amino acids which 
                                    tend to be in turn (e.g 'N', 'P', 'G', 'S')
    - 'sheet' (float) :             cumulative percentage of amino acids which 
                                    tend to be in turn (e.g  'E', 'M', 'A', 'L')
    """
    
    def __init__(self):
        self.lst_prot_prop = []
        self.dict_prop = {
            'tiny':['A', 'C', 'G', 'S', 'T'],
            'small':['A', 'C', 'F', 'G', 'I', 'L', 'M', 'P', 'V', 'W', 'Y'],
            'aliphatic':['A', 'I', 'L', 'V'],
            'aromatic':['F', 'H', 'W', 'Y'],
            'non_polar':['A', 'C', 'F', 'G', 'I', 'L', 'M', 'P', 'V', 'W', 'Y'],
            'polar':['D', 'E', 'H', 'K', 'N', 'Q', 'R', 'S', 'T', 'Z'],
            'charged':['B', 'D', 'E', 'H', 'K', 'R', 'Z'],
            'basic':['H', 'K', 'R'],
            'acidic':['B', 'D', 'E', 'Z'],
            'kyte_doolittle':{'A': 1.8,'C': 2.5,'D': -3.5,'E': -3.5,
                              'F': 2.8,'G': -0.4,'H': -3.2,'I': 4.5,
                              'K': -3.9,'L': 3.8,'M': 1.9,'N': -3.5,
                              'P': -1.6,'Q': -3.5,'R': -4.5,'S': -0.8,
                              'T': -0.7,'V': 4.2,'W': -0.9,'Y': -1.3},
            'rose':{'A': 0.74,'C': 0.91,'D': 0.62,'E': 0.62,'F': 0.88,
                    'G': 0.72,'H': 0.78,'I': 0.88,'K': 0.52,'L': 0.85,
                    'M': 0.85,'N': 0.63,'P': 0.64,'Q': 0.62,'R': 0.64,
                    'S': 0.66,'T': 0.70,'V': 0.86,'W': 0.85,'Y': 0.76}}
        
    
    def extract_properties(self, dict_seq):
        """
        Iterate through a list of protein sequences to compute 
        several physico-chemical properties
        """
        self.lst_prot_prop = []
        for self.seq_id, self.seq in dict_seq.items():
            self.compute_properties()
            self.lst_prot_prop.append(self.dict_features_value)
            
            
    def compute_properties(self):
        """
        Compute properties from the sequence.
        The computed properties are the ones described in the class 
        description.
        """
        seq_analysis = ProteinAnalysis(self.seq)
        seq_len=len(self.seq)
        aa_percent = seq_analysis.get_amino_acids_percent()
        sec_struct_fraction = seq_analysis.secondary_structure_fraction()
        self.dict_features_value = {
            'id':self.seq_id,
            'seq_len':seq_len,
            'gravy':seq_analysis.gravy(),
            'tiny':np.sum([aa_percent[aa] for aa in self.dict_prop['tiny'] if \
                aa in aa_percent]),
            'small':np.sum([aa_percent[aa] for aa in self.dict_prop['small'] \
                if aa in aa_percent]),
            'aliphatic':np.sum([aa_percent[aa] for aa in \
                self.dict_prop['aliphatic'] if aa in aa_percent]),
            'aromatic':np.sum([aa_percent[aa] for aa in \
                self.dict_prop['aromatic'] if aa in aa_percent]),
            'non_polar':np.sum([aa_percent[aa] for aa in \
                self.dict_prop['non_polar'] if aa in aa_percent]),
            'polar':np.sum([aa_percent[aa] for aa in self.dict_prop['polar'] \
                if aa in aa_percent]),
            'charged':np.sum([aa_percent[aa] for aa in \
                self.dict_prop['charged'] if aa in aa_percent]),
            'basic':np.sum([aa_percent[aa] for aa in self.dict_prop['basic'] \
                if aa in aa_percent]),
            'acidic':np.sum([aa_percent[aa] for aa in self.dict_prop['acidic'] \
                if aa in aa_percent]),
            'small':np.sum([aa_percent[aa] for aa in self.dict_prop['small'] \
                if aa in aa_percent]),
            'helix':sec_struct_fraction[0],
            'turn':sec_struct_fraction[1],
            'sheet':sec_struct_fraction[2]
        }


    def export_as_dataframe(self):
        """
        Return a pandas dataframe from a list of dict.
        """
        
        return pd.DataFrame(self.lst_prot_prop)

### The following two functions call the class to calculate these values 
### and create the correspondant dataframe from it.
### since the input is a dictionary, the function from_lst_to_dict converts
### a list of motifs to a dictionary of motifs.

def from_lst_to_dict(lst_motifs):
    """from_lst_to_dict
       ----------------
       This function converts a list of motifs
       to a dictionary of motifs. Where the key
       is the ID of the motif, and the value is
       the motif.
       
       Arguments:
       lst_seqs -- list of motifs
       
       Output:
       dict_motifs -- dictionary of motifs
    """
    
    # create an empty dictionary
    dict_motifs = {}
    
    # inserting keys and values in the empty dictionary
    for motif in lst_motifs:
        dict_motifs[motif] = motif
    
    return dict_motifs


def feature_calculation(dict_seqs):
    """feature_calculation
       -------------------
       This function calculates feature values
       for a dictionary of sequences.
       
       Arguments:
       dict_seqs -- dictionary of sequences
                                 
       Output:
       df_features -- pandas dataframe with values of features
    """
    
    
    # calculate the sequence values
    props = MotifProperties()
    props.extract_properties(dict_seqs)
    
    # export the sequence values in a pandas dataframe
    df_features = props.export_as_dataframe() 
    
    return df_features
