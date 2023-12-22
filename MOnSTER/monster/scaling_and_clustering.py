import pandas as pd
pd.set_option('display.max_columns', None)
import sys
import os
from sklearn.preprocessing import StandardScaler, MinMaxScaler, RobustScaler
import scipy.cluster.hierarchy as shc
from scipy.cluster.hierarchy import cut_tree
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import davies_bouldin_score
from scipy import signal

### Data Scaling

def features_data_scaling(df_features, scaling_method):
    """features_data_scaling
       ---------------------
       This function scales data from
       feature calculation.
       
       Arguments:
       df_features -- pandas dataframe with
                      values of features
       scaling_method -- 0 or 1 or 2:
                         where 0 -> standard scaling
                         where 1 -> minmax scaling
                         where 2 -> robust scaling
       Output:
       df_features_scl -- pandas dataframe with
                        scaled values of features
    """
    ### Creates a copy of the original dataframe
    #
    df_features_scl = df_features.copy()

    ### Drop duplicated lines (e.g. seqs that are present several times)
    #
    df_features_scl = df_features_scl.drop_duplicates()

    ### Set the id as index
    #
    df_features_scl.set_index('id', inplace=True)

    ### Rescale the data
    #
    #   Here, we use a Standard Scaling
    if scaling_method == 0:
        df_features_scl = pd.DataFrame(
            StandardScaler().fit_transform(df_features_scl),
            index=df_features_scl.index, columns=df_features_scl.columns
        )
    
    #   Here, we use a MinMax Scaling
    if scaling_method == 1:
        df_features_scl = pd.DataFrame(
            MinMaxScaler().fit_transform(df_features_scl),
            index=df_features_scl.index, columns=df_features_scl.columns
        )
    
    #   Here, we use a Robust Scaling
    if scaling_method == 2:
        df_features_scl = pd.DataFrame(
            RobustScaler().fit_transform(df_features_scl),
            index=df_features_scl.index, columns=df_features_scl.columns
        )
    
    return df_features_scl


### Clustering

def dendogram_of_motifs(path_directory, link_matrix, df_motifs_scl, best_distance):
    """dendogram_of_motifs
       -------------------
       This function generates the dendogram
       of the clustered motifs.
       
       Arguments:
       path_directory -- path of the directory where the df has to be saved
       link_matrix -- linkage matrix of clustering
       df_motifs_scl -- pandas dataframe of
                        motifs' scaled feature values.
       best_distance -- distance at which the cut of the tree in CLUMPs
                        is performed.
       
       Output:
       dendogram -- figure of the dendogram
    """
    dendogram = plt.figure(figsize=(12, 30))
    dend = shc.dendrogram(link_matrix,
                         labels=df_motifs_scl.index,
                         leaf_font_size=10,
                         orientation='right',
                         color_threshold= best_distance)
    plt.axvline(x=best_distance)
    plt.title("Dendrogram of motifs", fontsize=25);
    plt.savefig(path_directory + '/dendogram.pdf', format = 'pdf', bbox_inches="tight")
    
    return dendogram


def db_score(link_matrix, nb_distances, df_motifs_scl):
    """db_score
       --------
       This function calculates the davies bouldin score
       for a given number of intervals of distances (nb_distances)
       
       Arguments:
       link_matrix -- linkage matrix of clustering,
       nb_distances -- number of intervals of distances at which
                       calculate the davies bouldin score
       df_motifs_scl -- pandas dataframe of
                        motifs' scaled feature values.
                        
       Output:
       dict_davies_bouldin_results -- all the distances tested and the
                                      corresponding davies bouldin score
    """
    # Calculate maximum distance in distance matrix.
    max_dist = shc.maxdists(link_matrix).max()
    # Calculate nb_distances
    possible_distances = np.linspace(1, max_dist, nb_distances)

    # Iterate on all the distances to perform the cut
    dict_davies_bouldin_results = {}
    X = np.array(df_motifs_scl)
    for distance in possible_distances:

        # create a df with CLUMPs obtained at that distance
        
        motifs_CLUMPs = pd.DataFrame(
            {'motif':df_motifs_scl.index, 'CLUMP':cut_tree(
                link_matrix, height= distance).ravel()})
        labels = np.array(motifs_CLUMPs.CLUMP)
        db_index = davies_bouldin_score(X, labels)
        dict_davies_bouldin_results.update({distance: db_index})

    
    return dict_davies_bouldin_results


def calculate_df_motifs_CLUMPs_of_best_cut(link_matrix, nb_distances, df_motifs_scl):
    """calculate_df_motifs_CLUMPs_of_best_cut
       --------------------------------------
       This function iterates on the possible distances of CLUMPs
       for the tree cut and calculates the corresponding davies_bouldin_score.
       Then, returns the dataframe df_motifs_CLUMPs calculated
       using the best distance.
    
       Arguments:
       link_matrix -- linkage matrix of clustering
       nb_distances -- number of intervals of distances at which
                       calculate the davies bouldin score
       df_motifs_scl -- pandas dataframe of
                        motifs' scaled feature values.
    
       Output:
       dict_davies_bouldin_results -- all the distances tested and the
                                      corresponding davies bouldin score.
       best_distance -- best distance for the
                        cutting of the tree
                        according to the calculated function.
       df_motifs_CLUMPs -- pandas dataframe of motifs and CLUMPs.
    
    """

    dict_davies_bouldin_results = db_score(link_matrix, nb_distances, df_motifs_scl)
    
    # Find first local minimum in dict_davies_bouldin_results
    x = np.array(list(dict_davies_bouldin_results.values()))
    min_peakind = signal.find_peaks_cwt(1/x, np.arange(1,10))
    best_davies = x[min_peakind][0]
    # extract best distance by best db score
    distances = list(dict_davies_bouldin_results.keys())
    db_scores = list(dict_davies_bouldin_results.values())
    position = db_scores.index(best_davies)
    best_distance = distances[position]

    # Create df using the best distance
    df_motifs_CLUMPs = pd.DataFrame(
        {'motif':df_motifs_scl.index,
         'CLUMP':cut_tree(link_matrix, height = best_distance).ravel()})
    
    return best_distance, df_motifs_CLUMPs, dict_davies_bouldin_results

def get_newick(node, parent_dist, leaf_names, newick='') -> str:
    """
    Convert sciply.cluster.hierarchy.to_tree()-output to Newick format.
    :param node: output of sciply.cluster.hierarchy.to_tree()
    :param parent_dist: output of sciply.cluster.hierarchy.to_tree().dist
    :param leaf_names: list of leaf names
    :param newick: leave empty, this variable is used in recursion.
    :returns: tree in Newick format
    """
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parent_dist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parent_dist - node.dist, newick)
        else:
            newick = ");"
        newick = get_newick(node.get_left(), node.dist, leaf_names, newick=newick)
        newick = get_newick(node.get_right(), node.dist, leaf_names, newick=",%s" % (newick))
        newick = "(%s" % (newick)
        return newick

    
def to_newick(df_motifs, link_matrix):
    tree = shc.to_tree(link_matrix, False)
    leaf_names= df_motifs.index
    newick_tree= get_newick(tree, tree.dist, leaf_names)
    return newick_tree


def db_score_minimization_plot(path_directory, dict_davies_bouldin_results):
    """db_score_minimization_plot
       --------------------------
       This function plots the results of the minimization of the davies
       bouldin score.
       
       Arguments: 
       path_directory -- path of the directory where the df has to be saved
       dict_davies_bouldin_results -- all the distances tested and the
                                      corresponding davies bouldin score                    
    """
    # plot results of minimization davies bouldin score
    plt.plot(*zip(*sorted(dict_davies_bouldin_results.items())))
    plt.savefig(path_directory + '/minimization_db_score.pdf', format = 'pdf', bbox_inches="tight")
    plt.show()


def motif_clustering(path_directory, df_motifs_scl, nb_distances):
    """motif_clustering
       ----------------
       This function performs a
       hierarchical/agglomerative clustering
       on the motifs.
       
       Then it converts the tree in newick format.
       
       Then it calculates the best distance to cut
       the tree to obtain the CLUMPs
       (cluster of motifs of proteins).
       To do that it uses the davies_bouldin score.
       
       Then it cuts the tree into CLUMPs.
       
       Arguments:
       path_directory -- path of the directory where the df has to be saved
       df_motifs_scl -- pandas dataframe of
                        motifs' scaled feature values
       nb_distances -- number of intervals of distances at which
                       calculate the davies bouldin score
       
       Output:
       link_matrix -- linkage matrix of motifs.
       dict_davies_bouldin_results -- all the distances tested and the
                                      corresponding davies bouldin score
       df_motifs_CLUMPs -- pandas dataframe of the motifs
                           and the correspondant CLUMP
       dendogram -- dendogram of the clustered motifs
                           
       
    """
    
    ### Clustering (creation of the linkage matrix using the method ward)
    #
    link_matrix = shc.linkage(df_motifs_scl, method='ward')
    

    ### Iterating on the possible distances of CLUMPs for the
    ### tree cut and calculate corresponding davies_bouldin_score
    #
    best_distance, df_motifs_CLUMPs, dict_davies_bouldin_results = calculate_df_motifs_CLUMPs_of_best_cut(
        link_matrix, nb_distances, df_motifs_scl)
    db_score_minimization_plot(path_directory, dict_davies_bouldin_results)
    
    # Generate dendogram of motifs, with CLUMPs obtained cutting
    # at the best_distance
    dendogram = dendogram_of_motifs(path_directory, link_matrix, df_motifs_scl, best_distance)
    
    ### Getting Newick format of tree
    newick_tree = to_newick(df_motifs_scl, link_matrix)
    
    return link_matrix, best_distance, dict_davies_bouldin_results, df_motifs_CLUMPs, dendogram, newick_tree