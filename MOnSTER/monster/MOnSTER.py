from monster.import_data import *
from monster.features_extraction import *
from monster.scaling_and_clustering import *
from monster.scoring import *
from monster.export_output import *

def MOnSTER(seqs_path_pos, seqs_path_neg, motifs_path,
            scaling_method, nb_distances):
    """MOnSTER
       -------
       MOnSTER is a pipeline to find CLUMPs (CLUster of Motifs of Proteins)
       that are discriminative between two datasets of sequences
       that will be called positive and negative.
       
       The CLUMPs are found by gathering together motifs that show
       similar physicochemical properties.
       These physicochemical properties can show different ranges of
       values, for this reason, we propose known 3 methods of scaling
       the data for this step.
       
       After clustering, the user will visualize the resulting dendogram.
       To find the best compromise between number of CLUMPs and having
       sufficiently different CLUMPs, MOnSTER employes the Davies Bouldin score.
       It is then up to the user to choose number of intervals
       of distances to find the best cut of the tree to obtain the most
       reasonable number of CLUMPs.
       (see doc for more information)
       
       To determine which are the most discriminative CLUMPs, MOnSTER
       calculates the MOnSTER score for each of them. The score
       is between 0 and 2 and takes into account: the physicochemical
       properties of the CLUMPs and the two datasets, and the occurrences
       of the motifs and the CLUMPs in the two datasets.
       The physicochemical properties that are employed for the MOnSTER_score
       are those that show a statistically significant difference
       in the two datasets.
       
       MOnSTER scores the CLUMPs the highest if they show values of these
       physicochemical properties that are in line with the positive
       dataset and if they are more present in the positive dataset
       compared to the negative one.
       
       Arguments:
       seqs_path_pos -- path of the positive dataset protein fasta file
                        (protein alphabeth and capital letters).
       seqs_path_neg -- path of the negative dataset protein fasta file
                        (protein alphabeth and capital letters).
       motifs_path -- path of the file containing the motifs (capital letters
                      one motif per line)
       scaling_method -- 0 or 1 or 2:
                         where 0 -> standard scaling
                         where 1 -> minmax scaling
                         where 2 -> robust scaling
       nb_distances -- number of intervals of distances at which
                       calculate the davies bouldin score
       
       Output:
       df_motifs_CLUMPs -- pandas dataframe of the motifs
                           and the correspondant CLUMP
       dendogram -- dendogram of the clustered motifs
       newick_tree -- tree in newick format
       MOnSTER_score_results -- pandas dataframe with the MOnSTER_score
                                results.
    """
    
    
    ## Create output directory
    feature_calculation_and_scaling_results, clustering_and_CLUMPs, MOnSTER_score, MOnSTER_analysis = create_output_directory_MOnSTER(
    'MOnSTER_results')



    ### Import data
    #
    pos_dict = import_fasta_sequences_as_dict(seqs_path_pos)
    neg_dict = import_fasta_sequences_as_dict(seqs_path_neg)
    lst_motifs = import_list_motifs(motifs_path)



    ### feature calculation
    
    ## creating dictionary of motifs
    #
    dict_motifs = from_lst_to_dict(lst_motifs)

    ## dfs with features' values
    df_motifs_features = feature_calculation(dict_motifs)
    pos_dset_feat = feature_calculation(pos_dict)
    neg_dset_feat = feature_calculation(neg_dict)
    # storing results into folder feature_calculation_and_scaling_results
    save_df_to_directory(df_motifs_features, 'df_motifs_features', feature_calculation_and_scaling_results)
    save_df_to_directory(pos_dset_feat, 'pos_dset_feat', feature_calculation_and_scaling_results)
    save_df_to_directory(neg_dset_feat, 'neg_dset_feat', feature_calculation_and_scaling_results)


    ### scaling data
    #
    df_motifs_scl = features_data_scaling(df_motifs_features, scaling_method)
    # storing results into folder feature_calculation_and_scaling_results
    save_df_to_directory(df_motifs_scl, 'df_motifs_scl', feature_calculation_and_scaling_results)


    ### clustering
    #
    link_matrix, best_distance, dict_davies_bouldin_results, df_motifs_CLUMPs, dendogram, newick_tree = motif_clustering(
    clustering_and_CLUMPs, df_motifs_scl, nb_distances)
    # storing results into folder clustering_and_CLUMPs
    save_df_to_directory(df_motifs_CLUMPs, 'df_motifs_CLUMPs', clustering_and_CLUMPs) # a df for all the CLUMPs
    save_lsts_motifs(df_motifs_CLUMPs, clustering_and_CLUMPs) # a df per CLUMP
    
    

    ### Scoring the CLUMPs

    ## formatting data
    #
    df_all_motifs_all_features, pos_dset_feat, neg_dset_feat = format_input_data(
        df_motifs_features, df_motifs_CLUMPs, pos_dset_feat, neg_dset_feat)

    ## find occurrences of non redundant motifs
    #
    df_general, df_start_end_position_pos, df_start_end_position_neg = find_occurrences_of_mots_in_datasets(
        df_motifs_CLUMPs, pos_dict, neg_dict)
    lst_motifs_mask, df_general_non_redundant, motif_counts, df_cnt_seq_per_cluster = find_occ_and_nb_seqs(df_motifs_CLUMPs, df_general)
    # storing results into folder MOnSTER_score
    save_df_to_directory(df_start_end_position_pos, 'df_start_end_position_pos', MOnSTER_score)
    save_df_to_directory(df_start_end_position_neg, 'df_start_end_position_neg', MOnSTER_score)

    ## MOnSTER_score
    #
    MOnSTER_score_results, df_sign_features_p_values = MOnSTER_score(pos_dset_feat, neg_dset_feat,
                                          motif_counts, df_cnt_seq_per_cluster,
                                          df_motifs_CLUMPs, df_all_motifs_all_features)
    # storing results into folder MOnSTER_score
    save_df_to_directory(df_sign_features_p_values, 'df_sign_features_p_values', MOnSTER_score)
                                          
    # summary of results
    df_summary_results_monster = final_results_table(MOnSTER_score_results, df_motifs_CLUMPs, df_cnt_seq_per_cluster)
    save_df_to_directory(df_summary_results_monster, 'df_summary_results_monster', clustering_and_CLUMPs)


    return df_motifs_CLUMPs, lst_motifs_mask, df_general_non_redundant, dendogram, newick_tree, MOnSTER_score_results, motif_counts, df_cnt_seq_per_cluster, df_summary_results_monster
