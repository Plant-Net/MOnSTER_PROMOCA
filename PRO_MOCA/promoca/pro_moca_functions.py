import pandas as pd
import numpy as np
import sys
import os




def create_output_directory(directory, parent_directory = None):
    """create_output_directory
       -----------------------
       This function creates output directories.
       
       Arguments:
       directory -- the name of the directory to be created
                    for example: directory
       parent_directory -- string with path of parent directory
                           for example: '/Users/.../folder'
       
       Output:
       path -- path of the new directory
    """
    
    ## If parent_directory is not given
    if parent_directory is None:
        parent_directory = os.getcwd()
        path = os.path.join(parent_directory, directory)
        
    ## If parent_directory is given
    else :
        path = os.path.join(parent_directory, directory)
        
    ## Checking if the directory already exists
    try:
        os.mkdir(path)
    except OSError as error:
        print('folder exists')
    
    return path




## First section
# Here we import data, we extract a list of lists of motifs and finally define the classification scheme.
#


# import data
def import_data_pro_moca(df_motifs_CLUMPs_path):
    """import_data_pro_moca
       --------------------
       This function imports data for PRO-MOCA.
       
       Arguments:
       df_motifs_CLUMPs_path -- path of the tsv file with the pandas 
                                dataframe of the motifs and the correspondant
                                CLUMP.
                           
       Output:
       df_motifs_CLUMPs -- pandas dataframe of the motifs
                           and the correspondant CLUMP.
    """
    df_motifs_CLUMPs = pd.read_csv(df_motifs_CLUMPs_path)
    
    if 'Unnamed: 0' in df_motifs_CLUMPs:
        df_motifs_CLUMPs.drop(columns = 'Unnamed: 0', inplace= True)
        
    return df_motifs_CLUMPs


# define classification scheme
def classification_AA(classification_scheme):
    """classification_AA
       -----------------
       This function creates a dictionary of classification of AA.
       The user chooses the classification system, depending on that
       the keys of the dictionary will be the initial of the 
       class. The values will be a list of the AA in that class.
       
       Arguments:
       classification_scheme -- a string indicating the classification
                                scheme.
                                Where classification_scheme can be:
                                'che' -> chemical
                                'h' -> hydrophobicity
                                'cha' -> charge
                                'ssp' -> secondary structure propensity
       
       Output:
       dict_AA_prop -- dictionary with the classification of the AA
                       in classes of AA.
                       
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
        dict_AA_prop = {'P': ['G','S','T','Y','C'],
                        'H': ['A','V','L','I','P','W','F','M'],
                        'N': ['Q','N'],
                        'A': ['D','E'],
                        'B': ['K','R','H']}
        print('chemical classification scheme chosen', '\n', '\n',
              'the following are the classes of AAs: ','\n',
              'polar -> P', '\n', 'hydrophobic -> H', '\n', 'neutral -> N', 
              '\n', 'acid -> A', '\n', 'basic -> B')
            
    if classification_scheme == 'h':
        dict_AA_prop = {'P': ['R','K','D','E','N','Q'],
                        'N': ['S','G','H','T','A','P'],
                        'B': ['Y','V','M','C','L','F','I','W']} 
        print('hydrophobicity classification scheme chosen', '\n','\n', 
              'the following are the classes of AAs: ','\n',
              'hydrophilic -> P','\n','neutral -> N','\n','hydrophobic -> B')
          
    if classification_scheme == 'cha':
        dict_AA_prop = {'P': ['K','R','H'],
                        'N': ['D','E'],
                        'U': ['A','C','F','G','I','L','M','N',
                              'P','Q','S','T','V','W','Y']}
        print('charge classification scheme chosen', '\n', '\n',
              'the following are the classes of AAs: ','\n',
              'positive -> P','\n','negative -> N','\n','neutral -> U')
        
    if classification_scheme == 'ssp':
        dict_AA_prop = {'H': ['E','A','L','M','Q','K','R','H'],
                        'S': ['V','I','Y','C','W','F','T'],
                        'T': ['G','N','P','S','D']}
        print('secondary structure propensity classification scheme chosen', 
              '\n','\n', 'the following are the classes of AAs: ','\n',
              'helix -> H','\n','sheet -> S','\n','turn -> T')
    
    return dict_AA_prop


# extract list of lists of motifs
def extract_lsts_from_df_motifs_CLUMPs(df_motifs_CLUMPs):
    """extract_lsts_from_df_motifs_CLUMPs
       ----------------------------------
       This function creates a list of lists. Each list 
       contains the motifs of a CLUMP.
       
       Arguments:
       df_motifs_CLUMPs -- pandas dataframe of the motifs
                           and the correspondant CLUMP.
                           
       Output:
       lst_lsts_CLUMPs -- list of lists of motifs of CLUMPs.
    """
    
    lst_lsts_CLUMPs = []
    CLUMPs = list(df_motifs_CLUMPs.CLUMP.unique())
    
    ## Extracting the list of motifs of a CLUMP
    for CLUMP in CLUMPs:
        lst_motifs_of_CLUMP = list(df_motifs_CLUMPs[
            df_motifs_CLUMPs['CLUMP'] == CLUMP].motif)
        
        ## Creating a list of lists of motifs of all 
        ## the CLUMPs
        lst_lsts_CLUMPs.append(lst_motifs_of_CLUMP)
        
    return lst_lsts_CLUMPs





## Second section
# This is the section where we find the central position of each pattern
# to then align all the patterns in the next section.
# To do that we need a list of motifs as input.
# The output of the first section is a list of lists of motifs, for this
# reason, we will iterate on all the lists in the main list.
# Please note the following functions work on a list per time, for
# this reason pro-moca needs to iterate on all the lists in lst_lsts_CLUMPs
#


# conversion of list of motifs in list of patterns
def convert_motifs_in_patterns(lst_motifs, dict_AA_prop):
    """convert_motifs_in_patterns
       --------------------------
       This function converts a list of motifs in a list
       of patterns of letters from a classification scheme 
       that was previously chosen.
       
       Arguments:
       lst_motifs -- list of protein motifs.
       dict_AA_prop -- dictionary with the classification of the AA
                       in classes of AA.
       
       Output: 
       lst_patterns -- list of patterns
    """
    ## Conversion of the motifs of AA in patterns of chemical properties

    # create list of patterns of a CLUMP
    lst_patterns = []

    # for motif in list of motifs
    for motif in lst_motifs:

        # create empty string that will be filled with the 
        # letters of the chemical classes
        pattern = ''

        # for AA in the motif
        for AA in motif:

            # for the chemical class of AA
            for letter in dict_AA_prop.keys():

                # if the AA of the motif is in the chemical class of AA
                if AA in dict_AA_prop[letter]:
                    pattern = pattern + letter

        lst_patterns.append(pattern)
        
    return lst_patterns


# extract and generate a list of letters of a position of patterns
def extract_letter_and_repet_at_pos(lst_patterns, position_pattern):
    """extract_letter_and_repet_at_pos
       ---------------------------------------
       This function extracts all the letter at a position
       and their repetitions at a certain position of the patterns.
       
       Arguments : 
       lst_patterns -- list of patterns 
       position_patterns -- position of patterns
       
       Output: 
       dict_repet_letter -- dictionary of nb of repetitions 
                                of each letter at a position
       lst_letter_position -- list of letter at a position
    """
    ## Extract all the letter at a position 
    lst_letter_position = []
    for pattern in lst_patterns:
        position = pattern[position_pattern]
        lst_letter_position.append(position)
    
    ## Create a dictionary of the nb of repetitions of each letter
    # at that position
    dict_repet_letter = {}
    for letter in lst_letter_position:
        dict_repet_letter[
            letter] = lst_letter_position.count(letter)
    
    return lst_letter_position, dict_repet_letter


# get the index for all the occurrences of a value in a list
# source: https://thispointer.com/python-how-to-find-all-indexes-of-an-item-in-a-list/
def get_index_positions(list_of_elems, element):
    ''' Returns the indexes of all occurrences of give element in
    the list- listOfElements '''
    index_pos_list = []
    index_pos = 0
    while True:
        try:
            # Search for item in list from indexPos to the end of list
            index_pos = list_of_elems.index(element, index_pos)
            # Add the index position in list
            index_pos_list.append(index_pos)
            index_pos += 1
        except ValueError as e:
            break
    return index_pos_list


def most_common_letter_at_position(lst_patterns, dict_repet_letter, position_pattern):
    """most_common_letter_at_position
       ----------------------------------
       This function calculates the most common letter at 
       a certain position.
       
       Arguments: 
       lst_patterns -- list of patterns
       dict_repet_letter -- dictionary of nb of repetitions 
                                of each letter at a position
       position_pattern -- position of patterns
       
       Output:
       letter_for_position -- most common letter at 
                              a certain position.
                              if it is a single letter, the variable type is a string
                              if it there are multiple most common letters,
                              the variable is a list.
    """
    
    ## Extracting the number of repetitions of all the letter
    # at a position
    lst_letter_position, dict_repet_letter = extract_letter_and_repet_at_pos(
        lst_patterns, position_pattern)
    
    if '_' in dict_repet_letter.keys():
        del dict_repet_letter['_']
    
    ## Extracting the most common letter for a position 
    # extracting the letteres
    letters = list(dict_repet_letter.keys())
    # extracting the nb of repetitions
    nb_repetitions = list(dict_repet_letter.values())
    # calculating the highest repetitions
    highest_repetition = max(nb_repetitions)
    freq_high_rep = nb_repetitions.count(highest_repetition)
    
    ## Checking if more than a letter has the highest repetition
    if freq_high_rep == 1:
        # extracting the corresponding letter 
        position = nb_repetitions.index(highest_repetition)
        letter_for_position = letters[position]
    else:
        positions = get_index_positions(nb_repetitions, highest_repetition)
        letter_for_position = []
        for index in positions:
            letter_for_position.append(letters[index])

    return letter_for_position


def calculate_position_letter(lst_letter_position, letter_for_position):
    """calculate_position_letter
       -----------------------------
       This function calculates all the letter at a position to be 
       then used for alignement.
       
       Arguments: 
       lst_letter_position -- list of letter at a position
       letter_for_position -- most common letter at 
       a certain position.
       
       Output:
       lst_pos -- list of all letter and gaps of a position
       lst_indexes_letter_for_position -- indexes of 
                                              most frequent 
                                              letter at that
                                              position
       dict_sep_info_indexes -- dictionary with separated information
                                of the indexes of multiple most freq letters
       
       """
    # extracting index of the most repeated letter in the list of 
    # lst_letter_first_position
    
    # creating empty dictionary to store info of indexes of different most
    # frequent letters
    dict_sep_info_indexes = {} 
    # if there is only one letter, no information will be stored in the dict
    
    # checking if it is one or more most repeated letters
    if type(letter_for_position) == str:
        lst_indexes_letter_for_position = get_index_positions(
        lst_letter_position, letter_for_position)
        dict_sep_info_indexes[letter_for_position] = lst_indexes_letter_for_position
        
        
    if type(letter_for_position) == list:
        lst_indexes_letter_for_position = []
        for letter in letter_for_position:
            # empty list to store indexes of each letter
            lst_indx_single_letter = []
            indxs = get_index_positions(lst_letter_position, letter)
            for index in indxs:
                # here we store mixed information for all the letters
                lst_indexes_letter_for_position.append(index)
                # here we store information separated for each letter
                lst_indx_single_letter.append(index)
            # storing information for each most repeated letter into a dictionary
            dict_sep_info_indexes[letter] = lst_indx_single_letter
            
    # create the final list of elements of the first position
    lst_pos = list(np.repeat(['_'], len(lst_letter_position)))
    
    if type(letter_for_position) == str:
        # insert in the list the letter at the right index
        for indexes in lst_indexes_letter_for_position:
            lst_pos[indexes] = letter_for_position
    if type(letter_for_position) == list:
        for letter, lst_indxs in dict_sep_info_indexes.items():
            # we are taking the indexes in the list for each of the letters
            for indexes in lst_indxs:
                lst_pos[indexes] = letter
    
    return lst_pos, lst_indexes_letter_for_position, dict_sep_info_indexes


# final function
# that gathers : 
#               extract_letter_and_repet_at_pos
#               most_common_letter_at_position
#               get_index_positions
#               calculate_position_letter
def generate_lst_elements_position(lst_patterns, position_pattern):
    """generate_lst_elements_position
       ------------------------------
       This function generates a list of all the elements of a position
       for an alignment.
       
       Arguments: 
       lst_patterns -- list of patterns  
       
       Output:
       lst_pos -- list of all letter and gaps of a position
    """
    
    ## Extract all the letter at a position
    ## and their repetitions at a certain position of the patterns
    #
    lst_letter_position, dict_repet_letter = extract_letter_and_repet_at_pos(
        lst_patterns, position_pattern)
    
    ## Calculate the most common letter at 
    ##  a certain position.
    #
    letter_for_position = most_common_letter_at_position(
    lst_patterns, dict_repet_letter, position_pattern)
    letter_for_position
    
    ## calculates all the letter at a position to be 
    ## then used for alignement
    #
    lst_pos, lst_indexes_letter_for_position, dict_sep_info_indexes = calculate_position_letter(
        lst_letter_position, letter_for_position)
    
    return lst_pos, lst_indexes_letter_for_position, letter_for_position, dict_sep_info_indexes


# iterations and removal of first position each time
def remove_first_position(lst_patterns, lst_pos, lst_indexes_letter_for_position):
    """remove_first_position
       ---------------------
       This function removes the first position of the
       patterns that show the most frequent letter at that position
       
       Arguments:
       lst_patterns -- list of patterns
       lst_pos -- list of all letter and gaps of a position
       lst_indexes_letter_for_position -- indexes of 
                                              most frequent 
                                              letter at that
                                              position
       Output:
       new_lst_patterns -- list of patterns
       set_indexes -- set of indexes of patterns whose 
                      first position has been removed
    """
    
    new_lst_patterns = lst_patterns.copy()
    
    for index in lst_indexes_letter_for_position:
        # setting the parameters
        pattern = lst_patterns[index]
        position = 0
        new_character = ''
        
        # removing the first position
        temp = list(pattern)
        temp[position] = new_character
        new_pattern = "".join(temp)
        new_lst_patterns[index] = new_pattern
        
    return new_lst_patterns

# extract and generate a list of lists of letters of all positions of the patterns
def iterations_elements_at_all_positions(lst_patterns):
    """iterations_elements_at_all_positions
       ------------------------------------
       This function finds the most common letters at all positions.
       
       Arguments:
       lst_patterns -- list of patterns
       
       Output:
       lst_lsts_pos -- list of lists of all letter and gaps 
                       at all positions
       lst_lsts_indexes_letter_for_position -- list of lists of all
                                                   the indexes of the 
                                                   patterns that don't
                                                   show a gap at a position
                                                   (hence they show the most
                                                   common letter)
       lst_new_lst_patterns -- list of lists of the new patterns generated at
                               each iteration, by the removal of the first 
                               position from the pattern that presented the 
                               most common letter.
       lst_lengths -- list with values of the highest repetition at each iteration 
    """
    
    # to store the results of the iterations
    lst_lsts_pos = []
    lst_lsts_indexes_letter_for_position = []
    lst_new_lst_patterns = []
    
    # adding the list of pattern with no first position removed the first time. 
    lst_new_lst_patterns.insert(0, lst_patterns)
    
    print('original lst_patterns', lst_patterns, '\n', '\n')
    
    # here we want to create a list 
    # with values of the highest repetition at that iteration
    # we will later use to find which of the lists 
    # includes the highest number of patterns
    lst_lengths = []  
  
    # First time generating the lst_elements_position
    lst_pos, lst_indexes_letter_for_position, letter_for_position, dict_sep_info_indexes = generate_lst_elements_position(lst_patterns, 0)
    
    # here we are checking if there is one or more letters that are the most
    # repeated at that position.
    # if there is one:
    if type(letter_for_position) == str :
        lst_lengths.append(len(lst_indexes_letter_for_position))
    # if there is more than one, we have to divide this value
    # by the number of these letters. This is beacause the information of
    # the indexes is all stored in the same list. 
    # so for n = number of most repeated letters and 
    # r = repetition of a single letter, the 
    # len(lst_indexes_letter_for_position) = n*r
    # for this reason we divide len(lst_indexes_letter_for_position)/n
    # to get r
    if type(letter_for_position) == list:
        n = len(letter_for_position)
        r = int(len(lst_indexes_letter_for_position)/n)
        lst_lengths.append(r)
    
    
    # checking if it is one or more letters to have the highest frequency 
    if type(letter_for_position) == str:
        print('\n', '\n','-----------------', '\n', 
          '-----------------', '\n', 'letter_for_position', letter_for_position, '\n')
        print('\n', 'lst_pos', lst_pos)
    if type(letter_for_position) == list:
        print('\n', '\n','-----------------', '\n', 
          '-----------------', '\n', 'letter_for_position')
        for letter in letter_for_position:
            print(letter, ',')
        print('\n')
    
    lst_lsts_pos.append(lst_pos)
    lst_lsts_indexes_letter_for_position.append(
        lst_indexes_letter_for_position)
    
    # removing the first position
    new_lst_patterns = remove_first_position(
        lst_patterns, lst_pos, lst_indexes_letter_for_position)
    lst_new_lst_patterns.append(new_lst_patterns)

    set_indexes = set(lst_indexes_letter_for_position)
    nb_indexes = len(set_indexes)
    
    
    print('\n', 'lst_indexes_letter_for_position', lst_indexes_letter_for_position, '\n')
    
    while '' not in new_lst_patterns:
        lst_pos, lst_indexes_letter_for_position, letter_for_position, dict_sep_info_indexes = generate_lst_elements_position(
            new_lst_patterns, 0)
        print('\n', '\n', '-----------------', '\n', 'letter_for_position', letter_for_position, '\n')
        new_lst_patterns = remove_first_position(
        new_lst_patterns, lst_pos, lst_indexes_letter_for_position)
        print('\n', 'lst_pos', lst_pos)
        print('\n', 'lst_indexes_letter_for_position', lst_indexes_letter_for_position, '\n')
        lst_new_lst_patterns.append(new_lst_patterns)
        
        
        for el in lst_indexes_letter_for_position:
            set_indexes.add(el)
        nb_indexes = len(set_indexes)
        lst_lsts_pos.append(lst_pos)
        lst_lsts_indexes_letter_for_position.append(
            lst_indexes_letter_for_position)
        
        # here we are checking if there is one or more letters that are the most
        # repeated at that position.
        # if there is one:
        if type(letter_for_position) == str :
            lst_lengths.append(len(lst_indexes_letter_for_position))
        # if there is more than one, we have to divide this value
        # by the number of these letters. This is beacause the information of
        # the indexes is all stored in the same list. 
        # so for n = number of most repeated letters and 
        # r = repetition of a single letter, the 
        # len(lst_indexes_letter_for_position) = n*r
        # for this reason we divide len(lst_indexes_letter_for_position)/n
        # to get r
        if type(letter_for_position) == list:
            n = len(letter_for_position)
            r = int(len(lst_indexes_letter_for_position)/n)
            lst_lengths.append(r)
        
    print('\n', 'lst_lsts_indexes_letter_for_position', lst_lsts_indexes_letter_for_position, '\n')   
    print('\n','variable lst_lengths (list of number of times most repeated letter is repeated in that list)', lst_lengths, '\n')   
        
        
    print('\n','\n')
    
    print('\n', 'lst_new_lst_patterns (this is the list of the obtained patterns at each iteration) :', 
          lst_new_lst_patterns, '\n', '-----------------', '\n', '-----------------', '\n', '\n')
    lst_new_lst_patterns.pop()
    
    return lst_lsts_pos, lst_lsts_indexes_letter_for_position, lst_new_lst_patterns, lst_lengths





## Third section
# In this section alignment is performed.
# To do that, pro_moca looks for the list, from lst_lsts_indexes_letter_for_position,
# with the highest number of patterns.
# This one will define the central position for the alignment.
# For further information on how the central letter for each pattern is found, please
# read further.




# establish central position
def establish_central_position_alignment(lst_motifs, lst_new_lst_patterns,
                                         lst_lsts_indexes_letter_for_position, 
                                         lst_patterns, lst_lengths):
    """establish_central_position_alignment
       ------------------------------------
       This function calculates which one among the lists in 
       lst_lsts_indexes_letter_for_position contains the highest number
       of patterns.
       Then establishes the letters at the central position of the alignment
       
       Arguments:
       lst_motifs -- list of protein motifs.
       lst_new_lst_patterns -- list of lists of the new patterns generated at
                               each iteration, by the removal of the first 
                               position from the pattern that presented the 
                               most common letter.
       lst_lsts_indexes_letter_for_position -- list of lists of all
                                                   the indexes of the 
                                                   patterns that don't
                                                   show a gap at a position
                                                   (hence they show the most
                                                   common letter)
       lst_patterns -- list of patterns
       lst_lengths -- list with values of the highest repetition at each iteration 
       
       Output:
       lst_patterns -- list of patterns
       lst_center -- list of central position
       lst_motifs -- list of protein motifs.
    """ 
    
    # from lst_lengths we want to find the max, and then
    # extract the index of the corresponding list from 
    # lst_lsts_indexes_letter_for_position
    # this one will be the chosen list for the central position.
    index_chosen_lst = lst_lengths.index(max(lst_lengths))
    print('\n', 'index_chosen_lst (index of the list chosen for central position)', index_chosen_lst, '\n')
    
    # determining the list of central position
    lst_center = lst_new_lst_patterns[
        index_chosen_lst]
    
    # in case some of the patterns got all the letters deleted
    # in the previous sections, a gap is added.
    for el in range(len(lst_center)):
        if lst_center[el] == '':
            lst_center[el] = '_'
            lst_patterns[el] = lst_patterns[el]+'_'
            lst_motifs[el] = lst_motifs[el]+ '_'
            
    print('\n', 'this is the list of the central positions :', lst_center, 
          '\n', '\n', '-----------------', '\n', '-----------------', '\n', '\n')
    
    return lst_patterns, lst_center, lst_motifs


# find position of pattern corresponding to central position
def find_pos_pattern_of_central_position(lst_patterns, lst_center):
    """find_pos_pattern_of_central_position
       ------------------------------------
       This function finds the positions of the original patterns 
       that correspond to the central position previously calculated.
       
       Arguments:
       lst_patterns -- list of patterns
       lst_center -- list of central position
       
       Output:
       lst_positions_beg_subp -- list of positions (of each pattern) of 
                                 letter at central position
                                 in alignment.
    """
    
    ## creating list where we will store information of the position
    ## in the original pattern that corresponds to the beginning 
    ## position of the sub-pattern
    lst_positions_beg_subp = []

    ## iterating on the list of patterns, pattern by pattern
    for pattern in lst_patterns:
        print('ORIGINAL PATTERN IS :', pattern)
        print('\n')

        lst_single_pattern_composition = ['']

        ## iterating on each position of each pattern
        ## through a decreasing range.
        ## We want to find, in the original pattern,
        ## the corresponding sub-pattern
        ## To a position that is as close as possible to the end

        for position in range(len(pattern) -1, -1, -1):
            ## Here we are inserting one position by one to the 
            ## list, because we want to check, each time we insert one
            ## letter at the beginning, we want to check at which position
            ## we find the same identical sub-pattern.
            forming_pattern = pattern[position] +lst_single_pattern_composition[0]
            lst_single_pattern_composition.insert(0, forming_pattern)


            ## extracting the index of the subpattern we want to
            ## compare the forming pattern to. 
            ## to do that we simply need to take the index of the pattern
            index_subpattrn = lst_patterns.index(pattern)
            a = lst_single_pattern_composition[0]
            b = lst_center[index_subpattrn]
            print('letters of pattern: ', a)
            print('sub-pattern to check: ', b)
            print(a == b)
            print('\n')


            ## removing the '' we added at the beginning to the list
            ## it was necessary because otherwise we couldn't start
            ## adding strings to the list by simple addition (+)
            if '' in lst_single_pattern_composition:
                lst_single_pattern_composition.remove('')

            ## if the position of the beginning of the sub-pattern is found:
            if a == b:
                print(position)
                print('\n')
                lst_positions_beg_subp.append(position)

        print('\n')
        print('\n')
        print('end of pattern')
        print('\n')
        print('--------------------------')
        print('--------------------------')
        print('\n')
    
    return lst_positions_beg_subp
    

# length alignment
def find_nb_positions_alignment(lst_patterns, lst_positions_beg_subp):
    """find_nb_positions_alignment
       ---------------------------
       This function calculates the maximum number of 
       letters in the patterns.
       To do that, it starts from the central position 
       and then counts how many letters there are before
       and after the central position, for each pattern.
       
       Arguments: 
       lst_patterns -- list of patterns
       lst_positions_beg_subp -- list of positions (of each pattern) of 
                                 letter at central position
                                 in alignment.

       Output: 
       lst_patt_centr_position -- list of the patterns at central position
       lst_letters_after_center -- list of letters after center position
       positions_left -- nb positions before the central position
       positions_right -- nb positions after the central position
    """
    
    ## creating list of lists of: pattern and index central position
    ## of the pattern
    lst_patt_centr_position = []
    for index in range(len(lst_patterns)):
        lst_couple = []
        lst_couple.append(lst_patterns[index])
        lst_couple.append(lst_positions_beg_subp[index])
        lst_patt_centr_position.append(lst_couple)


    ## How many positions before the central one? (LEFT)
    positions_left = max(lst_positions_beg_subp)
    print('nb positions before the central position', positions_left)


    ## How many positions after the central one? (RIGHT)
    lst_letters_after_center = []
    for index in range(len(lst_patterns)):
        el = lst_patterns[index]
        nb_gaps_right = len(el) - (lst_positions_beg_subp[index] + 1)
        lst_letters_after_center.append(nb_gaps_right)
    positions_right = max(lst_letters_after_center)
    print('nb positions after the central position', positions_right)
    
    return lst_patt_centr_position, lst_letters_after_center, positions_left, positions_right


# final alignment
def final_alignment(lst_patt_centr_position, lst_letters_after_center, 
                    positions_left, positions_right):
    """final_alignment
       ---------------
       This function adds the necessary gaps to the left and the right
       of the patterns, referring to the central position.
       
       Arguments: 
       lst_patt_centr_position -- list of the patterns at central position
       lst_letters_after_center -- list of letters after center position
       positions_left -- nb positions before the central position
       positions_right -- nb positions after the central position
       
       Output:
       lst_final_alignment -- list with the final alignment
    """
    
    ## Final alignment
    lst_final_alignment = []

    for index in range(len(lst_patt_centr_position)):

        pattern = lst_patt_centr_position[index][0]
        print('pattern: ', pattern)

        indx_central_pos = lst_patt_centr_position[index][1]
        print('index central position: ', indx_central_pos)

        central_letter = pattern[indx_central_pos]
        print('central letter: ', central_letter)

        print('\n')



        #### Adding gaps to the left (from central position)
        nb_gaps_pattern_left = positions_left - indx_central_pos
        print(positions_left, '-', indx_central_pos, '=', 
              nb_gaps_pattern_left)
        print(nb_gaps_pattern_left, '<=', 'number of gaps on the left')
        resulting_pattern = nb_gaps_pattern_left * '_' + pattern



        #### Adding gaps to the right (from central position)
        nb_gaps_pattern_right = positions_right - lst_letters_after_center[
            index]
        print(positions_right, '-', lst_letters_after_center[index], '=', 
              nb_gaps_pattern_right)
        print(nb_gaps_pattern_right, '<=', 'number of gaps on the right')
        resulting_pattern = resulting_pattern + nb_gaps_pattern_right * '_'

        print('resulting pattern :', resulting_pattern)

        print('\n')
        print('-----------------')


        lst_final_alignment.append(resulting_pattern)
    
    print('\n', 'This is the alignment of the patterns', '\n', lst_final_alignment, '\n', '\n')
    
    return lst_final_alignment


# reconversion of letters in AAs
def re_conversion_of_letters_in_AAs(lst_final_alignment, lst_motifs, lst_patterns):
    """re_conversion_of_letters_in_AAs
       -------------------------------
       This function re-converts the letters in the alignment in
       AAs.
       
       Arguments:
       lst_final_alignment -- list with the final alignment
       lst_patterns -- list of patterns
       lst_motifs -- list of protein motifs.
       
       Output:
       results_alignment_motifs -- list of alignment with AAs as the letters.
    """
    
    ## creating final list of the alignment
    results_alignment_motifs = []
    
    ## conversion of the letters in AAs
    for index in range(len(lst_final_alignment)):
        pattern_in_alignment = lst_final_alignment[index]
        pattern = lst_patterns[index]
        motif = lst_motifs[index]
        
        mot_in_alignment = pattern_in_alignment.replace(pattern, motif)
        results_alignment_motifs.append(mot_in_alignment)
        
    print('This is the resulting alignment of the motifs', '\n', results_alignment_motifs, 
          '\n', '\n', '\n', '\n', '\n', '\n', '\n', '\n', 
          '...', '\n', '\n', '\n', '\n', '\n', '\n', '\n', '\n')
    
    return results_alignment_motifs
