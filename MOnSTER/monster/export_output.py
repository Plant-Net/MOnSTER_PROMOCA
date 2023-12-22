import pandas as pd
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



def create_output_directory_MOnSTER(name_parent_directory):
    """create_output_directory_MOnSTER
       -------------------------------
       This function creates the output directories of MOnSTER.
       
       Arguments:
       name_parent_directory -- name of the folder of the parent directory
       
       Output:
       feature_calculation_and_scaling_results -- path of directory
       clustering_and_CLUMPs -- path of directory
       MOnSTER_score -- path of directory
       MOnSTER_analysis -- path of directory
    """

    ## Creating main directory
    path_MOnSTER_results = create_output_directory(name_parent_directory)

    ## Creating all the directories within the main one
    feature_calculation_and_scaling_results = create_output_directory(
    'feature_calculation_and_scaling', path_MOnSTER_results)
    
    clustering_and_CLUMPs = create_output_directory(
    'clustering_and_CLUMPs', path_MOnSTER_results)
    
    MOnSTER_score = create_output_directory(
    'MOnSTER_score', path_MOnSTER_results)
    
    MOnSTER_analysis = create_output_directory(
    'MOnSTER_analysis', path_MOnSTER_results)

    return feature_calculation_and_scaling_results, clustering_and_CLUMPs, MOnSTER_score, MOnSTER_analysis



## Save df to given directory
def save_df_to_directory(df, df_name, path_directory):
    """save_df_to_directory
       --------------------
       This function saves a given pandas dataframe (df) to a given directory.
       
       Arguments:
       df -- a pandas dataframe
       df_name -- name of pandas dataframe
       path_directory -- path of the directory where the df has to be saved
    """

    path_file = path_directory + '/' + df_name + '.tsv'
    print(path_file)
    
    df.to_csv(path_file)
    

## Save each list of motifs into a file
def save_lsts_motifs(df_motifs_CLUMPs, path_directory):
    """save_lsts_motifs
       ----------------
       This function saves each list of motifs as a file. All files
       are stored into a given directory.
    
       Arguments:
       df_motifs_CLUMPs -- pandas dataframe of the motifs
                           and the correspondant CLUMP
       path_directory -- path of the directory where the df has to be saved
    """
    
    CLUMPs = list(df_motifs_CLUMPs.CLUMP.unique())
    
    for CLUMP in CLUMPs:
        df_motifs = pd.DataFrame(df_motifs_CLUMPs[df_motifs_CLUMPs[
            'CLUMP'] == CLUMP].motif)
        df_motifs = df_motifs.rename(columns = {'motif': str(CLUMP)})
        df_motifs.to_csv(
            path_directory + '/CLUMP_' + str(CLUMP) + '.tsv',
            index = None)
        
