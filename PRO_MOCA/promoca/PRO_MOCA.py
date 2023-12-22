from promoca.pro_moca_functions import *
import sys
import os

def pro_moca(df_motifs_CLUMPs_path, classification_scheme): 
    """pro_moca
       --------
       PRO-MOCA (PROtein-MOtif Characteristics Aligner) is an
       aligner of motif sequences.
       PRO-MOCA aligns short sequences based on the characteristics
       of the AAs rather than the AAs themselves.
       
       The algorithm doesn't employ a substitution matrix.
       Instead, PRO-MOCA employes a classification scheme of the 
       AAs to convert the motifs into patterns of letters. Then
       looks for a central position of the patterns, where a 
       certain letter is repeated more than in the others.
       Based on the central position of the alignment, it then
       calculates the number of gaps ('_') to add to the right and 
       to the left of all the patterns, after centering them.
       At the end it reconverts the patterns into motifs, 
       and gives a list of motifs with gaps to the left and 
       to the right.
       
       No gaps are admitted into the sequences, only at the 
       extremities. This is because motifs are very short sequences
       hence introducing gaps could mean losing the information
       of the motif itself.
       
       Arguments:
       df_motifs_CLUMPs -- pandas dataframe of the motifs
                           and the correspondant CLUMP.
       classification_scheme -- a string indicating the classification
                                scheme.
                                Where classification_scheme can be:
                                'che' -> chemical
                                'h' -> hydrophobicity
                                'cha' -> charge
                                'ssp' -> secondary structure propensity
       
       Output: 
       dict_alignments_all_lsts_motifs -- dictionary with results
                                          of the alignment 
                                          stored for all the CLUMPs.
                                          The keys of the dictionary are
                                          the CLUMPs and the values
                                          are the list of alignment 
                                          with AAs as the letters.
                                          
       ----------------------------------------------------------------
       Classification_schemes: 
                       
                       chemical :
                       # polar -> P
                       # hydrophobic -> H
                       # neutral -> N
                       # acid -> A
                       # basic -> B
                       dict_AA_prop = {'P': ['G','S','T','Y','C'],
                                       'H': ['A','V','L','I','P','W','F','M'],
                                       'N': ['Q','N'],
                                       'A': ['D','E'],
                                       'B': ['K','R','H']}

                       hydrophobicity :
                       # hydrophilic -> P
                       # neutral -> N
                       # hydrophobic -> B
                       # NB: to distinguish Hydrophilic and Hydrophobic, 
                       #     we choose the letters P and B 
                       #     (hydroPhilic and hydrophoBic)
                       dict_AA_prop = {'P': ['R','K','D','E','N','Q'],
                                       'N': ['S','G','H','T','A','P'],
                                       'B': ['Y','V','M','C','L','F','I','W']}

                       charge :
                       # positive -> P
                       # negative -> N
                       # neutral -> U
                       # NB: to distinguish Negative and Neutral, 
                       #     we choose the letters N and U
                       #     (Negative and neUtral)
                       dict_AA_prop = {'P': ['K','R','H'],
                                       'N': ['D','E'],
                                       'U': ['A','C','F','G','I','L','M','N',
                                             'P','Q','S','T','V','W','Y']}
                                             
                       secondary structure propensity :
                       # helix (α) -> H
                       # sheet (β) -> S
                       # turn -> T
                       dict_AA_prop = {'H': ['E','A','L','M','Q','K','R','H'],
                                       'S': ['V','I','Y','C','W','F','T'],
                                       'T': ['G','N','P','S','D']}
                                       
       ----------------------------------------------------------------
       Sources of the classification schemes:    
       
       chemical, hydrophobicity, charge: 
       https://weblogo.threeplusone.com/manual.html 
       
       secondary structure propensity:
       http://www.bio.brandeis.edu/classes/biochem104/alpha_helix_2007.pdf


    """
   
   
    if classification_scheme == 'che':
        class_scheme_name = 'chemical'
    if classification_scheme == 'h':
        class_scheme_name = 'hydrophobicity'
    if classification_scheme == 'cha':
        class_scheme_name = 'charge'
    if classification_scheme == 'ssp':
        class_scheme_name = 'secondary_structure_propensity'
    name_output = 'output_promoca_' + class_scheme_name
    
    # create output directory
    create_output_directory(class_scheme_name, parent_directory = None)
    
    
    
    ## First section
    # Here we import data, we extract a list of lists of motifs and finally 
    #define the classification scheme.
    #
   
    # import data of motifs  
    df_motifs_CLUMPs = import_data_pro_moca(df_motifs_CLUMPs_path)
   
    # define classification scheme
    dict_AA_prop = classification_AA(classification_scheme)
   
    # extract list of lists of motifs
    lst_lsts_CLUMPs = extract_lsts_from_df_motifs_CLUMPs(df_motifs_CLUMPs)


   
    ## Second section
    # This is the section where we find the central position of each pattern
    # to then align all the patterns in the next section.
    # To do that we need a list of motifs as input.
    # The output of the first section is a list of lists of motifs, for this
    # reason, we will iterate on all the lists in the main list.
    # Please note the following functions work on a list per time, for
    # this reason pro-moca needs to iterate on all the lists in lst_lsts_CLUMPs
    #
   
   
    dict_alignments_all_lsts_motifs = {}
    nbs_CLUMPs = list(np.arange(0, len(df_motifs_CLUMPs.CLUMP.unique())))
   
    # iterating on all the lists
    for indx in range(len(lst_lsts_CLUMPs)):
        
        print('\n', '\n', '\n', '\n') 
        print('/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/')
        print('------------------------------------------')
        print('\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ ')
        print('\n', '\n')
        print('This is the alignment for CLUMP', indx)
        print('\n', '\n')
        
        lst_motifs = lst_lsts_CLUMPs[indx]
      
        # conversion of list of motifs in list of patterns
        lst_patterns = convert_motifs_in_patterns(lst_motifs, dict_AA_prop)
      
        # extract and generate a list of lists of letters of all positions 
        # of the patterns
        lst_lsts_pos, lst_lsts_indexes_letter_for_position, lst_new_lst_patterns, lst_lengths = iterations_elements_at_all_positions(
         lst_patterns)


   
        ## Third section
        # In this section alignment is performed.
        # To do that, pro_moca looks for the list, from
        # lst_lsts_indexes_letter_for_position,
        # with the highest number of patterns.
        # This one will define the central position for the alignment.
        # For further information on how the central letter for each pattern is
        # found, please read informations of single functions.
        #
      
        # establish central position
        lst_patterns, lst_key_test, lst_motifs = establish_central_position_alignment(lst_motifs, lst_new_lst_patterns,
                                         lst_lsts_indexes_letter_for_position, 
                                         lst_patterns, lst_lengths)
      
        # find position of pattern corresponding to central position
        lst_positions_beg_subp = find_pos_pattern_of_central_position(
            lst_patterns, lst_key_test)
      
        # length alignment
        lst_patt_centr_position, lst_letters_after_center, positions_left, positions_right = find_nb_positions_alignment(
            lst_patterns, lst_positions_beg_subp)
      
        # final alignment
        lst_final_alignment = final_alignment(lst_patt_centr_position, lst_letters_after_center, 
                    positions_left, positions_right)
      
        # reconversion of letters in AAs
        results_alignment_motifs = re_conversion_of_letters_in_AAs(
            lst_final_alignment, lst_motifs, lst_patterns)
      
        # storing the alignment in a dictionary
        # where the key is the number of the CLUMP and the value is the list of the
        # alignment 
        key = nbs_CLUMPs[indx]
        dict_alignments_all_lsts_motifs[key] = results_alignment_motifs
        
       
       
       
       
       
       
        
    ## save output files
    
    
    # save summary
    with open(class_scheme_name+ str('/') +name_output, 'w') as output:
        output.write('promoca output')
        output.write('\n')
        output.write('\n')
        output.write('--------------------------------------------------')
        output.write('\n')
        output.write('version: 0.0.0')
        output.write('\n')
        output.write('\n')
        output.write('classification scheme: ')
        output.write(class_scheme_name)
        output.write('\n')
        output.write('\n')
        output.write('--------------------------------------------------')
        output.write('\n')
        output.write('labels of CLUMPs:')
        output.write(str(nbs_CLUMPs))
        output.write('\n')
        output.write('\n')
        output.write('--------------------------------------------------')
        output.write('\n')
        output.write('alignment results:')
        output.write('\n')
        output.write('\n')
        for CLUMP, alignment in dict_alignments_all_lsts_motifs.items():
            output.write('CLUMP ')
            output.write(str(CLUMP))
            output.write(':')
            output.write('\n')
            output.write('\n')
            for motif in alignment:
                output.write(motif)
                output.write('\n')
            output.write('\n')
            output.write('\n')
            output.write('--------------')
            output.write('\n')
            output.write('\n')
  
    
    # save each alignment in a file in a folder
    for CLUMP, alignment in dict_alignments_all_lsts_motifs.items():
        df_alignment = pd.DataFrame(alignment, columns = [str(CLUMP)])
        df_alignment.to_csv(
            class_scheme_name+ str('/') + 'alignment_CLUMP_' + str(CLUMP) + '.tsv', 
            index= None)
        
       
       
       
       
       
       
        
    return dict_alignments_all_lsts_motifs
