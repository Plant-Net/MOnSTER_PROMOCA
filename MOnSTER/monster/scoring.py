import pandas as pd
import scipy
from scipy import stats
from scipy.stats import mannwhitneyu
import math
import numpy as np
import re


def format_input_data(df_motifs_features, df_clusters,
    pos_dset_feat, neg_dset_feat):
    """format_input_data
       -----------------
       This function formats input data in the correct
       way for it to work with the following functions.
       
       Arguments:
       df_motifs_features -- pandas dataframe with data of
                             feature values of the motifs, but no
                             information about the CLUMPs
       df_clusters -- pandas dataframe of the motif and corresponding CLUMP
       pos_dset_feat -- pandas dataframe with data of feature values
                        of the positive dataset
       neg_dset_feat -- pandas dataframe with data of feature values
                        of the negative dataset
       Output:
       df_all_motifs_all_features -- pandas dataframe with data of
                                     feature values of the motifs
       pos_dset_feat -- pandas dataframe with data of feature values
                        of the positive dataset
       neg_dset_feat -- pandas dataframe with data of feature values
                        of the negative dataset
    """
    df_motifs_features.rename(columns = {'id' : 'motif'}, inplace = True)
    df_all_motifs_all_features = df_clusters.merge(df_motifs_features)
    pos_dset_feat.drop(columns = 'id', inplace = True)
    neg_dset_feat.drop(columns = 'id', inplace = True)
    
    return df_all_motifs_all_features, pos_dset_feat, neg_dset_feat


## first section : calculate occurrences of motifs
#
def start_end_position(lst_motifs, dict_seqs, df_clusters, dataset):
    """start_end_position
       ------------------
       This function calculates the start and end position
       of the motifs in the sequences.
       
       Arguments:
       lst_motifs -- list of motifs
       dict_seqs -- dictionary of fasta sequences where the key is the
                    id and the value is the sequence
       df_clusters -- pandas dataframe of the motif and corresponding CLUMP
       dataset -- 'positive' or 'negative'
       
       Output:
       df_start_end_position -- pandas dataframe where:
                                first column is the motif
                                second column is the CLUMP
                                third column is the sequence id
                                fourth column is the start position
                                fifth column is the end position
                                sixth column is the dataset        
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
    df_start_end_position['dataset'] = dataset
    df_start_end_position = df_clusters.merge(df_start_end_position)
    
    return df_start_end_position


def occ_each_mot_in_each_seq(df_start_end_position, dataset):
    """occ_each_mot_in_each_seq
       ------------------------
       This function calculates the occurrences
       of each motif in each sequence.
       
       Arguments:
       df_start_end_position -- pandas dataframe where:
                                first column is the motif
                                second column is the CLUMP
                                third column is the sequence id
                                fourth column is the start position
                                fifth column is the end position
                                sixth column is the dataset
       dataset -- 'positive' or 'negative'
       
       Output:
       df_occ_seq -- pandas dataframe where:
                     first column is the motif
                     second column is the sequence id
                     third column is the number of occurrences
                     fourth column is the dataset
    """
    
    df_occ_seq = df_start_end_position.groupby(
        ['motif','seq_id']).size().reset_index(name='occ')
    
    df_occ_seq['dataset'] = dataset
    
    if df_occ_seq.shape == df_occ_seq.drop_duplicates().shape:
        return df_occ_seq


def df_occs_mots_CLUMPs_both_dsets(df_start_end_position_pos,
                                   df_start_end_position_neg, df_clusters):
    """df_occs_mots_CLUMPs_both_dsets
       ------------------------------
       This function gathers information of the occurrences
       of each motif in each sequence in each dataset, and of
       the corresponding CLUMP.
       
       Arguments:
       df_start_end_position_pos -- pandas dataframe where:
                                    first column is the motif
                                    second column is the CLUMP
                                    third column is the sequence id
                                    fourth column is the start position
                                    fifth column is the end position
                                    sixth column is the dataset
       df_start_end_position_neg -- pandas dataframe where:
                                    first column is the motif
                                    second column is the CLUMP
                                    third column is the sequence id
                                    fourth column is the start position
                                    fifth column is the end position
                                    sixth column is the dataset
       df_clusters -- pandas dataframe of the motif and corresponding CLUMP
       
       Output:
       df_general -- pandas dataframe with: motif, CLUMP, seq_id, occ, dataset
    """
    df_occ_seq_pos = occ_each_mot_in_each_seq(
        df_start_end_position_pos, "positive")
    df_occ_seq_neg = occ_each_mot_in_each_seq(
        df_start_end_position_neg, "negative")
    
    df_clusters.rename(columns = {'id' : 'motif'}, inplace = True)
    df_general = pd.concat(
        [pd.merge(df_clusters, df_occ_seq_pos, on='motif'),
         pd.merge(df_clusters, df_occ_seq_neg, on='motif')])
    
    return df_general

# final function of first section
def find_occurrences_of_mots_in_datasets(df_clusters, pos_dict, neg_dict):
    """find_occurrences_of_mots_in_datasets
       ------------------------------------
       This function finds the occurrences of motifs of CLUMPs
       in the positive and negative datasets, by finding their
       start and end position. It then combines the information
       in a pandas dataframe.
       
       Arguments:
       df_clusters -- pandas dataframe of the motif and corresponding CLUMP
       pos_dict -- dictionary of fasta sequences where the key is the
                   id and the value is the sequence of positive dataset
       neg_dict -- dictionary of fasta sequences where the key is the
                   id and the value is the sequence of negative dataset
       
       Output:
       df_general -- pandas dataframe with: motif, CLUMP, seq_id, occ, dataset
    """
    
    lst_motifs = list(df_clusters.motif)
    
    df_start_end_position_pos = start_end_position(
        lst_motifs, pos_dict, df_clusters, 'positive')
    df_start_end_position_neg = start_end_position(
        lst_motifs, neg_dict, df_clusters, 'negative')
    
    df_general = df_occs_mots_CLUMPs_both_dsets(df_start_end_position_pos,
                                                df_start_end_position_neg,
                                                df_clusters)
    
    return df_general, df_start_end_position_pos, df_start_end_position_neg
    


## second section : non redundant motifs
#
def find_extended_motifs(lst_motifs):
    """find_extended_motifs
       --------------------
       This function identifies extended motifs in CLUMPs.
       
       Arguments:
       lst_motifs -- list of motifs
       
       Output:
       lst_all_extended_motifs -- ???????
    """
    
    dict_extended_motifs = {}
    lst_motifs = sorted(lst_motifs, key=len)
    lst_known_motifs = []
    lst_all_extended_motifs = []
    for motif in lst_motifs:
        if motif not in lst_known_motifs:
            lst_extended_motifs = [
                m for m in lst_motifs if motif in m and motif != m]
            lst_known_motifs.append(motif)
            lst_known_motifs+=lst_extended_motifs
            if len(lst_extended_motifs)>0:
                lst_all_extended_motifs+=[[motif]+lst_extended_motifs]
                
    return lst_all_extended_motifs
    

def non_redundant_motifs(df_clusters, df_general):
    """non_redundant_motifs
       --------------------
       This function identifies non redundant extended motifs.
       Non redundant motifs (e.g. root motifs and non-extended motifs)
       are stored in a list called lst_motifs_mask.
       
       Arguments:
       df_clusters -- pandas dataframe of the motif and corresponding CLUMP.
       df_general -- pandas dataframe with: motif, CLUMP, seq_id, occ, dataset
       
       Output:
       lst_motifs_mask -- list of non redundant motifs.
       df_general_non_redundant -- pandas dataframe where the occurrences
       belong only to the non redundant motifs
       
    """
    
    lst_motifs = list(df_clusters.motif)
    
    lst_motifs_mask = []
    for c in df_clusters.CLUMP.unique():
        lst_motifs = df_clusters.loc[
            df_clusters.CLUMP==c, 'motif'].unique()
        lst_ext_motifs = find_extended_motifs(lst_motifs)
        lst_all_ext_motifs = [j for i in lst_ext_motifs for j in i]
        lst_non_ext_motifs = [
            m for m in lst_motifs if m not in lst_all_ext_motifs]
        lst_root_motifs = [el[0] for el in lst_ext_motifs]
        lst_motifs_mask += lst_root_motifs+lst_non_ext_motifs
        
    # Selecting the rows of the df_general where the occurrences
    # belong only to the non redundant motifs.
    df_general_non_redundant = df_general[df_general.motif.isin(
        lst_motifs_mask)]
    
    return lst_motifs_mask, df_general_non_redundant


def count_all_occurrences_mots_of_CLUMPs(df_general_non_redundant):
    """count_all_occurrences_mots_of_CLUMPs
       ------------------------------------
       This function calculates all the occurrences of each CLUMP
       in each dataset including multiple occurrences
       of a motif in a sequence.
       
       Arguments:
       df_general_non_redundant -- pandas dataframe where the occurrences
       belong only to the non redundant motifs.
       
       Output:
       motifs_counts -- pandas dataframe with the number of
                       occurrences of the CLUMP in the two
                       datasets.
       
    """
    motifs_counts = pd.DataFrame(df_general_non_redundant.groupby(
        ['CLUMP', 'dataset'])['occ'].sum())
    motifs_counts = motifs_counts.reset_index()
    motifs_counts = motifs_counts.pivot(
        index='CLUMP', columns='dataset', values='occ')
    motifs_counts = motifs_counts.rename_axis(None,axis=1)
    motifs_counts = motifs_counts.reset_index()
    
    return motifs_counts


def count_nb_seqs_containing_mots_of_CLUMPs(df_general_non_redundant):
    """count_nb_seqs_containing_mots_of_CLUMPs
       ---------------------------------------
       This function calculates how many sequences contain a
       motif of a CLUMP, without considering if the motif is present
       more than once in a sequence.
       Hence how many sequences are found by the CLUMP, considering the
       sequence only once.
       
       Arguments:
       df_general_non_redundant -- pandas dataframe where the occurrences
       belong only to the non redundant motifs.
       
       Output:
       df_cnt_seq_per_cluster -- pandas dataframe with the number of
                                 sequences found by the CLUMP in the
                                 two datasets.
       
    """
    # Calculating how many sequences contain a motif of a CLUMP
    # (how many sequences are found by the CLUMP, considering the
    # sequence only once) without considering if the motif is present
    # more than once in a sequence.
    df_cnt_seq_per_cluster = df_general_non_redundant.groupby([
        'CLUMP', 'seq_id', 'dataset']).size().reset_index(name='temporary')
    df_cnt_seq_per_cluster = df_cnt_seq_per_cluster.drop(
        columns = 'temporary')
    df_cnt_seq_per_cluster = df_cnt_seq_per_cluster.drop_duplicates()
    df_cnt_seq_per_cluster = df_cnt_seq_per_cluster.groupby(
        ['CLUMP', 'dataset']).size().reset_index()
    df_cnt_seq_per_cluster = df_cnt_seq_per_cluster.pivot(
        index='CLUMP', columns='dataset', values = 0)
    df_cnt_seq_per_cluster = df_cnt_seq_per_cluster.rename_axis(None,axis=1)
    df_cnt_seq_per_cluster = df_cnt_seq_per_cluster.reset_index()
    
    return df_cnt_seq_per_cluster


# final function of second section
def find_occ_and_nb_seqs(df_clusters, df_general):
    """find_occ_and_nb_seqs
       --------------------
       This function finds the non redudant motifs from the
       original df_clusters.
       Then finds the subset of df_general where only the
       non redundant motifs are considered.
       Finally calculates the occurrences of the motifs of the CLUMPs
       and the number of sequences found by the CLUMP.
       (for more information about the two outputs, read doc of following
       two functions:
       count_all_occurrences_mots_of_CLUMPs()
       count_nb_seqs_containing_mots_of_CLUMPs()).
       
       Arguments:
       df_clusters -- pandas dataframe of the motif and corresponding CLUMP.
       df_general -- pandas dataframe with: motif, CLUMP, seq_id, occ, dataset.
       
       Output:
       lst_motifs_mask -- list of non redundant motifs.
       motif_counts -- pandas dataframe with the number of
                       occurrences of the CLUMP in the two
                       datasets.
       df_cnt_seq_per_cluster -- pandas dataframe with the number of
                                 sequences found by the CLUMP in the
                                 two datasets.
    """
    lst_motifs_mask, df_general_non_redundant = non_redundant_motifs(
        df_clusters, df_general)
    motifs_counts = count_all_occurrences_mots_of_CLUMPs(
        df_general_non_redundant)
    df_cnt_seq_per_cluster = count_nb_seqs_containing_mots_of_CLUMPs(
        df_general_non_redundant)
    
    return lst_motifs_mask, motifs_counts, df_cnt_seq_per_cluster


## Third section: calculate MOnSTER score
#
# J1 and J2
def calculate_J1(motifs_counts, neg_dset_feat, pos_dset_feat):
    """calculate_J1
       ------------
       This function calculates J1. Modified Jaccard Index with data of the
       occurrences of the motifs of the CLUMPs.
       
       Arguments:
       motif_counts -- pandas dataframe with the number of
                       occurrences of the CLUMP in the two
                       datasets.
       neg_dset_feat -- pandas dataframe with data of feature values
                        of the negative dataset
       pos_dset_feat -- pandas dataframe with data of feature values
                        of the positive dataset
                       
       Output:
       J1 -- pandas dataframe with the results of the J1
    """
    motifs_counts['norm_negative'] = motifs_counts.negative/len(neg_dset_feat)
    motifs_counts['norm_positive'] = motifs_counts.positive/len(pos_dset_feat)

    ## J1
    motifs_counts[
        'jaccard_norm_1'
    ]= abs(motifs_counts.norm_negative/motifs_counts.norm_positive)
    J1 = pd.DataFrame(motifs_counts['jaccard_norm_1'])
    
    return J1


def calculate_J2(df_cnt_seq_per_cluster, neg_dset_feat, pos_dset_feat):
    """calculate_J2
       ------------
       This function calculates J2. Modified Jaccard Index with data of the
       number of sequences found by a CLUMP.
       
       Arguments:
       df_cnt_seq_per_cluster -- pandas dataframe with the number of
                                 sequences found by the CLUMP in the
                                 two datasets.
       neg_dset_feat -- pandas dataframe with data of feature values
                        of the negative dataset
       pos_dset_feat -- pandas dataframe with data of feature values
                        of the positive dataset
                       
       Output:
       J2 -- pandas dataframe with the results of the J2
    """
    df_cnt_seq_per_cluster[
        'norm_negative'] = df_cnt_seq_per_cluster.negative/len(
        neg_dset_feat)
    df_cnt_seq_per_cluster[
        'norm_positive'] = df_cnt_seq_per_cluster.positive/len(
        pos_dset_feat)

    ## J2
    df_cnt_seq_per_cluster[
        'jaccard_norm_2'
    ] = abs(df_cnt_seq_per_cluster.norm_negative/df_cnt_seq_per_cluster.norm_positive)
    J2 = pd.DataFrame(df_cnt_seq_per_cluster['jaccard_norm_2'])
    
    return J2


# final function of J1 and J2
def calculate_J1_and_J2(
    motifs_counts, df_cnt_seq_per_cluster, neg_dset_feat, pos_dset_feat):
    """calculate_J1_and_J2
       -------------------
       This function calculates the J1 and J2. Two modified Jaccard Indexes.
       (for more information about the input of the two indexes
       please read doc of the two following functions:
       calculate_J1()
       calculate_J2())
       
       Arguments:
       motif_counts -- pandas dataframe with the number of
                       occurrences of the CLUMP in the two
                       datasets.
       df_cnt_seq_per_cluster -- pandas dataframe with the number of
                                 sequences found by the CLUMP in the
                                 two datasets.
       neg_dset_feat -- pandas dataframe with data of feature values
                        of the negative dataset.
       pos_dset_feat -- pandas dataframe with data of feature values
                        of the positive dataset.
                        
       Output:
       df_jaccard_index -- pandas dataframe with results of
                           calculation of 1-J1 and 1-J2.
    """
    J1 = calculate_J1(motifs_counts, neg_dset_feat, pos_dset_feat)
    J2 = calculate_J2(df_cnt_seq_per_cluster, neg_dset_feat, pos_dset_feat)
    
    # here we are calculating the values of 1 - jaccard
    # since we want a score by maximization. With values to directly
    # sum the ones of the CLUMP score
    df_jaccard_index = pd.concat([J1, J2], axis = 1)
    df_jaccard_index = 1 - df_jaccard_index
    df_jaccard_index.insert(0, 'CLUMP', np.arange(0, len(df_jaccard_index)))

    return df_jaccard_index


# CLUMP_score
def feature_weight(pos_dset_feat, neg_dset_feat, feature_lst):
    """feature_weight
       --------------
       This function calculates which features are significant.
       To find that feature_weight uses the Mann-Whitney test,
       which calculates the significance of the difference between
       two datasets means (also for unpaired datasets),
       and gives a p-value as part of the output.
       null hypothesis (H0) : the difference between the two means
       is not statistically significant. p-value >= 0.05
       alternative hypothesis (H1): the difference between the two means
       is statistically significant. p-value < 0.05
       
       The features that result to be significant will receive a score.
       Where the score = -log10(p-value)
       
       Arguments:
       neg_dset_feat -- pandas dataframe with data of feature values
                        of the negative dataset.
       pos_dset_feat -- pandas dataframe with data of feature values
                        of the positive dataset.
       feature_lst -- list of features
                        
       Output:
       dict_significant_features -- dictionary with significant
                                    features as the keys and
                                    the p-values as the values
       df_sign_features_p_values -- pandas dataframe with significant
                                    features and corresponding p-value
    """
    # Calculating the p-values
    # Creating a list with the p-values in it
    p_values_lst = []
    m = len(feature_lst)
    for i in range(m):
        pos_values = pos_dset_feat.iloc[:, i]
        pos_values = list(pos_values)
        neg_values = neg_dset_feat.iloc[:, i]
        neg_values = list(neg_values)
        s, p = mannwhitneyu(pos_values, neg_values)
        p_values_lst.append(p)
    
    # Creating a dictionary with the feature as the key and the p-value
    # as the value.
    # The zip iterator is useful to pair each feature with its p-value
    # to then create the dictionary
    zip_iterator = zip(feature_lst, p_values_lst)
    dict_feat_p = dict(zip_iterator)
    
    # Creating a dictionary with the significant features their p-values
    dict_significant_feat = {}
    for feature, p_value in dict_feat_p.items():
        if p_value < 0.05:
            dict_significant_feat[feature] = p_value
    # Sort features in order of significance
    sign_feat = pd.DataFrame(dict_significant_feat, index = [0]).sort_values(
        by = 0, axis = 1, ascending = False).transpose().to_dict()
    dict_significant_features = sign_feat[0]
    
    df_sign_features_p_values = pd.DataFrame.from_dict(dict_significant_features,
    orient = 'index')
    
    return dict_significant_features, df_sign_features_p_values


def feature_selection(pos_dset_feat, neg_dset_feat):
    """feature_selection
       -----------------
       This function calculates which features help to distinguish
       the best the positive and the negative datasets (find an enrichment).
       
       Arguments:
       neg_dset_feat -- pandas dataframe with data of feature values
                        of the negative dataset.
       pos_dset_feat -- pandas dataframe with data of feature values
                        of the positive dataset.
                        
       Output:
       lst_signif_features -- list of significant features
       dict_feat_p_value_score_log -- dictionary of significant features
                                      as keys and dictionaries as values.
                                      These dictionaries contain information
                                      of the the p-value and the score.
       df_sign_features_p_values -- pandas dataframe with significant
                                    features and corresponding p-value
    """
    
    # What are the candidate features?
    pos_dset_feat = pos_dset_feat.drop(columns = ['seq_len'])
    neg_dset_feat = neg_dset_feat.drop(columns = ['seq_len'])

    feature_lst = []
    for col in pos_dset_feat.columns:
        feature_lst.append(col)
    print('the number of candidate features is:', len(feature_lst))
    
    # Finding significant features and their p-values
    dict_significant_features, df_sign_features_p_values = feature_weight(
        pos_dset_feat, neg_dset_feat, feature_lst)

    # Creating a list of the significant features in order of significance
    lst_signif_features = list(dict_significant_features)
    print('the number of significant features is:', len(lst_signif_features))
    print('significant features to calculate the CLUMP_score are:',
          lst_signif_features)
    
    # Here we are creating a dictionary with the features and their p-values
    # and scores from the -log10(p-value)
    dict_feat_p_value_score_log = {}
    for feature in lst_signif_features:
        dict_feat_p_value_score_log[feature] = {}
        dict_feat_p_value_score_log[feature][
            'p-value'] = dict_significant_features[feature]
        dict_feat_p_value_score_log[feature]['score'] = -math.log10(
            dict_feat_p_value_score_log[feature]['p-value'])
    
    return lst_signif_features, dict_feat_p_value_score_log, df_sign_features_p_values


def dataset_average_calculation(
    neg_dset_feat, pos_dset_feat, lst_signif_features):
    """dataset_average_calculation
       ---------------------------
       This function calculates the average of each significant
       feature for the two datasets.
       
       Arguments:
       neg_dset_feat -- pandas dataframe with data of feature values
                        of the negative dataset.
       pos_dset_feat -- pandas dataframe with data of feature values
                        of the positive dataset.
       lst_signif_features -- list of significant features
       
       Output:
       dict_pos_neg_means -- dictionary with the significant feature
                             as the key and dictionaries as the values.
                             The dictionaries contain information
                             of the positive and negative average value
                             for that feature.
    """
    # Creating a list with the positive dataset means of the significant
    # features
    dict_pos_means = dict(pos_dset_feat.loc[:, lst_signif_features].mean())
    # Creating a list with the negaive dataset means of the significant
    # features
    dict_neg_means = dict(neg_dset_feat.loc[:, lst_signif_features].mean())
    
    # Creating a dictionary with the features and their means
    # for the positive and the negative datasets
    dict_pos_neg_means = {}
    for feature in lst_signif_features:
        dict_pos_neg_means[feature] = {}
        dict_pos_neg_means[feature]['pos_mean'] = dict_pos_means[feature]
        dict_pos_neg_means[feature]['neg_mean'] = dict_neg_means[feature]
    
    return dict_pos_neg_means


def CLUMPs_average_calculation(df_all_motifs_all_features,
                               lst_signif_features):
    """CLUMPs_average_calculation
       --------------------------
       This function calculates the average of each significant
       feature for the CLUMPs.
    
       Arguments:
       df_all_motifs_all_features -- pandas dataframe with data of
                                     feature values of the motifs.
       lst_signif_features -- list of significant features
                             
       Output:
       df_clusters_means -- pandas dataframe with data of average of the CLUMPs
                            values for significant features.
       
    """
    
    # Creating a list of the significant features
    # and selecting subset of df_all_motifs_signif_features of
    # values only of significant features
    lst_signif_features.insert(0, 'CLUMP')
    df_all_motifs_signif_features= df_all_motifs_all_features.loc[
        :, lst_signif_features]
    lst_signif_features.pop(0)
    
    # Calculating the averages
    df_clusters_means = df_all_motifs_signif_features.groupby(
    'CLUMP').mean().reset_index()
    
    return df_clusters_means


def average_calculation(pos_dset_feat, neg_dset_feat, lst_signif_features,
    df_all_motifs_all_features):
    """average_calculation
       -------------------
       This function calculates the average of each significant feature
       for the CLUMPs and the two datasets.
       
       Arguments:
       neg_dset_feat -- pandas dataframe with data of feature values
                        of the negative dataset.
       pos_dset_feat -- pandas dataframe with data of feature values
                        of the positive dataset.
       lst_signif_features -- list of significant features
       df_all_motifs_all_features -- pandas dataframe with data of
                                     feature values of the motifs.
       
       Output:
       dict_pos_neg_means -- dictionary with the significant feature
                             as the key and dictionaries as the values.
                             The dictionaries contain information
                             of the positive and negative average value
                             for that feature.
       df_clusters_means -- pandas dataframe with data of average of the CLUMPs
                            values for significant features.
       
    """
    
    # Calculating the features' means in the positive and negative
    # datasets
    dict_pos_neg_means = dataset_average_calculation(
        neg_dset_feat, pos_dset_feat, lst_signif_features)
    
    # Calculating the features' means of the CLUMPs
    df_clusters_means = CLUMPs_average_calculation(
        df_all_motifs_all_features, lst_signif_features)
    
    return dict_pos_neg_means, df_clusters_means


def CLUMPs_sorting(
    df_clusters_means, lst_signif_features, dict_pos_neg_means):
    """CLUMPs_sorting
       --------------
       This function sorts CLUMPs and gives them votes
       in accordance with the following criteria.
       
        if average_positive_dataset - average_negative_dataset > 0:
            for clusters with average > average_positive_dataset:
                                        ranking from 1 to r
                r = clusters in the right half
            for clusters with average < average_positive_dataset: 0
        if average_positive_dataset - average_negative_dataset < 0:
            for clusters with average < average_positive_dataset:
                                        ranking from 1 to r
                r = clusters in the right half
            for clusters with average > average_positive_dataset: 0
        
       Arguments:
       df_clusters_means -- pandas dataframe with data of average of the CLUMPs
                            values for significant features.
       lst_signif_features -- list of significant features
       dict_pos_neg_means -- dictionary with the significant feature
                             as the key and dictionaries as the values.
                             The dictionaries contain information
                             of the positive and negative average value
                             for that feature.
       Output:
       dict_feat_scores -- dictionary with the features as the keys
                           and a list of the votes as the values.
                           N.B. The votes in the lists are not sorted
                           by CLUMP, but by vote.
                           It will be in the function CLUMPs_voting
                           that the votes will be assigned to the
                           corresponding CLUMPs.
        
    """
    dict_feat_scores = {}
    higher_scores = np.arange(1, len(df_clusters_means)+1)
    
    for feature in lst_signif_features:
        df_clu_feature = pd.DataFrame(df_clusters_means.loc[:, feature])
        lst_scores_feature = []
        lst_higher_scores_feature = []

        
        
        
        # if average_positive_dataset - average_negative_dataset > 0:
        #    for clusters with average > average_positive_dataset:
        #                                ranking from 1 to r
        #        r = clusters in the right half
        #    for clusters with average < average_positive_dataset: 0
        if dict_pos_neg_means[feature][
            'pos_mean'] - dict_pos_neg_means[feature]['neg_mean'] > 0:
            # sorting ascendingly
            df_clu_feature = df_clu_feature.sort_values(
                ascending = True, by = feature)
            
            for i in range(len(df_clu_feature)):
                
                # if the feature average is greater than the positive
                # dataset average
                if float(df_clu_feature.iloc[i]) > dict_pos_neg_means[
                    feature]['pos_mean']:
                    feat_higher_score = i+1
                    lst_higher_scores_feature.append(feat_higher_score)
                    new_list_higher_scores_features = list(np.arange(
                        1, len(lst_higher_scores_feature)+1))
                    new_list_higher_scores_features
                # if the feature average is less than the positive
                # dataset average
                else:
                    feat_score = 0
                    lst_scores_feature.append(feat_score)
                    
            lst_intermediate_scores = lst_scores_feature + new_list_higher_scores_features
        
        
        
        
        # if average_positive_dataset - average_negative_dataset < 0:
        #    for clusters with average < average_positive_dataset:
        #                                ranking from 1 to r
        #        r = clusters in the right half
        #    for clusters with average > average_positive_dataset: 0
        else:
            # sorting descendingly
            df_clu_feature = df_clu_feature.sort_values(
                ascending = False, by = feature)
            
            for i in range(len(df_clu_feature)):

                # if the feature average is less than the positive
                # dataset average
                if float(df_clu_feature.iloc[i]) < dict_pos_neg_means[
                    feature]['pos_mean']:
                    feat_higher_score = i+1
                    lst_higher_scores_feature.append(feat_higher_score)
                    new_list_higher_scores_features = list(np.arange(
                        1, len(lst_higher_scores_feature)+1))
                    new_list_higher_scores_features
                # if the feature average is greater than the positive
                # dataset average
                else:
                    feat_score = 0
                    lst_scores_feature.append(feat_score)
        
            lst_intermediate_scores = lst_scores_feature + new_list_higher_scores_features
            
        dict_feat_scores[feature] = lst_intermediate_scores
        
    return dict_feat_scores


def CLUMPs_voting(lst_signif_features, df_clusters_means,
                  dict_pos_neg_means, dict_feat_scores):
    """CLUMPs_voting
       -------------
       This function assignes the calculated votes (CLUMPs_sorting)
       to the corresponding CLUMPs.
       
       Arguments:
       lst_signif_features -- list of significant features
       df_clusters_means -- pandas dataframe with data of average of the CLUMPs
                            values for significant features.
       dict_pos_neg_means -- dictionary with the significant feature
                             as the key and dictionaries as the values.
                             The dictionaries contain information
                             of the positive and negative average value
                             for that feature.
       dict_feat_scores -- dictionary with the features as the keys
                           and a list of the votes as the values.
                           N.B. The votes in the lists are not sorted
                           by CLUMP, but by vote.
                           It will be this function CLUMPs_voting that will
                           assign the votes to the corresponding CLUMPs.
       
       Output:
       final_votes -- pandas DataFrame with the vote of each CLUMP
                      for each feature.
    """
    lst_final_scores = []
    for feature in lst_signif_features:
        df_clu_feature = pd.DataFrame(df_clusters_means.loc[:, feature])
        
        # if average_positive_dataset - average_negative_dataset > 0:
        if dict_pos_neg_means[feature][
            'pos_mean'] - dict_pos_neg_means[feature]['neg_mean'] > 0:
            # sort values of CLUMPs for that feature ascendingly
            df_clu_feature = df_clu_feature.sort_values(ascending = True, by = feature)
            df_clu_feature['score_'+feature] = dict_feat_scores[feature]
            df_clu_feature = df_clu_feature.sort_index()
            df_clu_feature = pd.DataFrame(df_clu_feature.iloc[:, 1])
            df_clu_feature = df_clu_feature.rename(columns = {'score_'+feature : feature})
            lst_clu_feature = list(df_clu_feature[feature])
        
        # if average_positive_dataset - average_negative_dataset < 0:
        else :
            # sort values of CLUMPs for that feature descendingly
            df_clu_feature = df_clu_feature.sort_values(ascending = False, by = feature)
            df_clu_feature['score_'+feature] = dict_feat_scores[feature]
            df_clu_feature = df_clu_feature.sort_index()
            df_clu_feature = pd.DataFrame(df_clu_feature.iloc[:, 1])
            df_clu_feature = df_clu_feature.rename(columns = {'score_'+feature : feature})
            lst_clu_feature = list(df_clu_feature[feature])
        lst_final_scores.append(lst_clu_feature)

    final_votes = pd.DataFrame(lst_final_scores).transpose()
    final_votes.columns = lst_signif_features
    final_votes= final_votes.transpose()
    
    return final_votes


def CLUMPs_scoring(
    final_votes, lst_signif_features, dict_feat_p_value_score_log):
    """CLUMPs_scoring
       --------------
       This function multiplies the final_votes by the calculated
       feature weights. Then sums the values of all the features
       for each CLUMP. Finally, normalizes the result in a range from
       0 to 1.
       
       Arguments:
       final_votes -- pandas DataFrame with the vote of each CLUMP
                      for each feature.
       lst_signif_features -- list of significant features
       dict_feat_p_value_score_log -- dictionary of significant features
                                      as keys and dictionaries as values.
                                      These dictionaries contain information
                                      of the the p-value and the score.
                                      
       Output:
       CLUMP_score_results -- pandas DataFrame with the normalized results of
                              the CLUMP_scoring for each CLUMP
    """
    
    # Extracting the scores of the significant features
    feat_scores = []
    for feature in lst_signif_features:
        feat_scores.append(dict_feat_p_value_score_log[feature]['score'])
    
    # Multiplying the CLUMPs votes by the feature scores (feature weights)
    df_clusters_weights_final = final_votes.mul(feat_scores, axis = 0)
    
    # Summing all the values of a CLUMP
    CLUMP_score_results = pd.DataFrame(df_clusters_weights_final.sum()).rename(
        columns = {0 : 'CLUMP_score'}).sort_values(
        by = 'CLUMP_score', ascending = False)
    
    # Normalizing the results in a range from 0 to 1
    CLUMP_score_results= (CLUMP_score_results - CLUMP_score_results.min()) / (
        CLUMP_score_results.max() - CLUMP_score_results.min())
    CLUMP_score_results.sort_index(inplace = True)
    
    return CLUMP_score_results


# final function CLUMP_score
def CLUMP_score(pos_dset_feat, neg_dset_feat, df_all_motifs_all_features):
    """CLUMP_score
       -----------
       This function calculates the CLUMP_score, which is
       based on the averages of the CLUMPs' and of the two datasets' values
       for a set of significant features. Significant features are those
       for which the function finds an enrichment in one of the two datasets,
       that is significant compared to the other.
       
       Arguments:
       neg_dset_feat -- pandas dataframe with data of feature values
                        of the negative dataset.
       pos_dset_feat -- pandas dataframe with data of feature values
                        of the positive dataset.
       df_all_motifs_all_features -- pandas dataframe with data of
                                     feature values of the motifs.
          
       Output:
       CLUMP_score_results -- pandas DataFrame with the normalized results of
                              the CLUMP_scoring for each CLUMP
       df_sign_features_p_values -- pandas dataframe with significant
                                    features and corresponding p-value
    """
    ## feature selection
    lst_signif_features, dict_feat_p_value_score_log, df_sign_features_p_values = feature_selection(
        pos_dset_feat, neg_dset_feat)
    
    ## average calculation
    dict_pos_neg_means, df_clusters_means = average_calculation(pos_dset_feat, neg_dset_feat, lst_signif_features,
    df_all_motifs_all_features)
    
    ## CLUMPs sorting
    dict_feat_scores = CLUMPs_sorting(
        df_clusters_means, lst_signif_features, dict_pos_neg_means)
    
    ## CLUMPs voting
    final_votes = CLUMPs_voting(lst_signif_features, df_clusters_means,
                            dict_pos_neg_means, dict_feat_scores)
    
    ## CLUMPs_scoring
    CLUMP_score_results = CLUMPs_scoring(
        final_votes, lst_signif_features, dict_feat_p_value_score_log)
    
    return CLUMP_score_results, df_sign_features_p_values


# final functions MOnSTER_score
def monster_ranking(MOnSTER_score_results):
    """monster_ranking
       ---------------
       This function ranks the CLUMPs based on the
       the results of the MOnSTER_score calculation.
       
       Arguments:
       MOnSTER_score_results -- pandas dataframe with the MOnSTER_score
                                results.
       Output:
       MOnSTER_score_results -- pandas dataframe with the MOnSTER_score
                                results.
    """
       
    ranking = np.arange(1, len(MOnSTER_score_results)+1)
    MOnSTER_score_results.sort_values(by = 'monster_score', ascending = False,
                                      inplace= True)
    MOnSTER_score_results['ranking'] = ranking
    MOnSTER_score_results = MOnSTER_score_results.sort_index()
    
    return MOnSTER_score_results


def MOnSTER_score(
    pos_dset_feat, neg_dset_feat, motifs_counts, df_cnt_seq_per_cluster,
    df_clusters, df_all_motifs_all_features):
    """MOnSTER_score
       -------------
       This function calculates the MOnSTER score.
       MOnSTER score is between 0 and 2.
       To do that, it calculates the CLUMP_score (from 0 to 1),
       J1 (from 0 to 1) and J2 (from 0 to 1), multiplies J1 and J2 by 0.5
       and then sums the 3 indexes.
       
       For further information about the calculation of these 3 indexes
       please read doc of the following functions:
       CLUMP_score()
       calculate_J1_and_J2()
       
       Arguments:
       pos_dset_feat -- pandas dataframe with data of feature values
                        of the positive dataset.
       neg_dset_feat -- pandas dataframe with data of feature values
                        of the negative dataset.
       motifs_counts -- pandas dataframe with the number of
                       occurrences of the CLUMP in the two
                       datasets.
       df_cnt_seq_per_cluster -- pandas dataframe with the number of
                                 sequences found by the CLUMP in the
                                 two datasets.
       df_clusters -- pandas dataframe of the motif and corresponding CLUMP.
       df_all_motifs_all_features -- pandas dataframe with data of
                                     feature values of the motifs.
       Output:
       MOnSTER_score_results -- pandas dataframe with the MOnSTER_score
                                results.
       df_sign_features_p_values -- pandas dataframe with significant
                                    features and corresponding p-value
    """
    
    # Calculating the CLUMP_score
    CLUMP_score_results, df_sign_features_p_values = CLUMP_score(pos_dset_feat, neg_dset_feat, df_all_motifs_all_features)
    # Calculating J1 and J2
    df_jaccard_index = calculate_J1_and_J2(
        motifs_counts, df_cnt_seq_per_cluster, neg_dset_feat, pos_dset_feat)
    
    # Multiplying J1 and J2 by 0.5
    df_jacc = df_jaccard_index.copy()
    df_jacc['jaccard_norm_1'] = df_jacc['jaccard_norm_1']*0.5
    df_jacc['jaccard_norm_2'] = df_jacc['jaccard_norm_2']*0.5
    
    # Calculating MOnSTER_score
    MOnSTER_score_results = pd.concat([CLUMP_score_results, df_jacc], axis =1)
    MOnSTER_score_results.drop(columns = 'CLUMP', inplace = True)
    MOnSTER_score_results = pd.DataFrame(MOnSTER_score_results.sum(
        axis =1)).rename(columns = {0 : 'monster_score'})
    MOnSTER_score_results.insert(0, 'CLUMP', list(df_clusters.CLUMP.unique()))
    
    MOnSTER_score_results = monster_ranking(MOnSTER_score_results)
    
    return MOnSTER_score_results, df_sign_features_p_values
    





def final_results_table(MOnSTER_score_results, df_clusters,
                        df_cnt_seq_per_cluster):
    """final_results_table
       -------------------
       This function summaries some of the most relevant results
       of MOnSTER.
       
       Arguments:
       MOnSTER_score_results -- pandas dataframe with the MOnSTER_score
                                results.
       df_clusters -- pandas dataframe of the motif and corresponding CLUMP.
       df_cnt_seq_per_cluster -- pandas dataframe with the number of
                                 sequences found by the CLUMP in the
                                 two datasets.
                                 
       Output:
       df_summary_results_monster -- pandas dataframe with the
                                     summary of the most relevant
                                     results of monster
    """
    
    ## Merging information:
    
    # MOnSTER_score_results -> results of monster
    df_summary_results_monster = MOnSTER_score_results.merge(
        # df_clusters -> motifs and their corresponding CLUMP
        df_clusters.groupby(by = 'CLUMP').count().reset_index()).merge(
        # df_cnt_seq_per_cluster[['CLUMP', 'norm_positive', 'norm_negative']]
        # % of nb_sequences found by the CLUMP
        df_cnt_seq_per_cluster[['CLUMP', 'norm_positive', 'norm_negative']]
    ).rename(columns = {'motif': '# motifs', 
                        'norm_positive' : 'occurrences_pos_dset %',
                       'norm_negative' : 'occurrences_neg_dset %'})
    df_summary_results_monster.sort_values(by = 'ranking', inplace = True)
    df_summary_results_monster.drop(columns = 'ranking', inplace = True)
    
    return df_summary_results_monster

