import gseapy as gp
from gseapy import enrichr
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.feature_extraction.text import TfidfVectorizer
from itertools import cycle
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



def get_enrichment(protein_set, enrich_genesets_dbs = 'KEGG_2021_Human', p_value_threshold=0.05):
    """
    Perform functional enrichment analysis for each latent dimension using gseapy.

    Parameters:
    - protein_set (list): List of protein sets.
    - enrich_genesets_dbs (list): List of gene sets or pathways of interest. Default = 'KEGG_2021_Human'

    Returns:
    - dict: Enrichment results for each protein pairs filtered with a user defined p-value threshold.
    """
       
    list_of_sets = protein_set  # Replace this with your actual list of sets
    n_sets = len(list_of_sets)

    # Declare the enrichment_results variable
    enrichment_results = {}

    # Create a DataFrame with "cluster" and "gene_id" columns
    data = []
    for idx, gene_set in enumerate(list_of_sets):
        data.extend([(idx, gene) for gene in gene_set])

    df = pd.DataFrame(data, columns=["pair", "gene_id"])

    # Group data by latent dimension
    grouped = df.groupby('pair')

    # Define gene sets or pathways of interest
    gene_sets = enrich_genesets_dbs

    # Perform functional enrichment analysis for each group
    for pair, group in grouped:
        genes = list(np.unique(group['gene_id']))

        # Perform enrichment analysis using gseapy
        enr = enrichr(gene_list=genes, gene_sets=gene_sets, outdir=None)

        # Print the length of the enrichment_results before and after each iteration
        enrichment_results[pair] = enr.results

    data_dict = {str(key): value for key, value in enrichment_results.items()}
    
    # Filter data with adjusted p-value above 0.05
    filtered_data_dict = {gene_set: data[data['Adjusted P-value'] <= p_value_threshold].copy() for gene_set, data in data_dict.items()}
     
        
    return filtered_data_dict, n_sets


def add_kegg_hierarchy(kegg_path, enr_df):
    
    hierarchy = pd.read_csv(kegg_path, index_col = 'Unnamed: 0')

    merged_df = pd.merge(enr_df, hierarchy, on='Term', how='left')

    # Identify rows with NaN values after the initial merge
    rows_with_nan = merged_df[merged_df.isna().any(axis=1)].drop(['Root', 'Subcategory'], axis=1)

    merged_df.dropna(subset=['Root', 'Subcategory'], inplace=True)

    hierarchy_updated = hierarchy
    # Assuming the additional information is separated by ' - '
    hierarchy_updated['Term'] = hierarchy_updated['Term'].str.split(' - ').str[0]

    # Redo the merge for the rows with NaN values
    updated_rows = pd.merge(rows_with_nan, hierarchy_updated, left_on='Term', right_on='Term', how='left')

    final_df = pd.concat([merged_df, updated_rows], axis=0).reset_index(drop=True)
    
    return final_df

def generate_KEGG_results_tables(enrichment_results, exp_name, results_path, n_sets, top_n_terms = None, kegg_path=None, save_table=True):
    """
    Generate tables based on the enrichment information.

    Parameters:
    - enrichment_results (list): DataFrames with TF-IDF values and additional information for each gene set.
    - exp_name (str): Experiment name.
    - results_path (str): Path to save the results.
    - save_table (bool): Whether to save the table as a CSV file.
    - top_n_terms (int): Number of top terms to select for each latent dimension if rows > 40.

    Returns:
    - pd.DataFrame: Combined DataFrame used for plotting.
    """
    
    # Create a list of dataframes with an additional column for latent dimension
    dfs_with_dimension = []
    for gene_set, data in enrichment_results.items():
        data['Pair'] = gene_set
        dfs_with_dimension.append(data)

    # Combine all dataframes into a single dataframe
    combined_df = pd.concat(dfs_with_dimension, axis=0)
    combined_df = combined_df.reset_index(inplace=False)
    combined_df.drop(columns=combined_df.columns[0], axis=1, inplace=True)
    
    if not kegg_path:
        print("Provide a valid path to the hierarchy file.")
        return None
    
    ### Add the grouping information from the kegg hierarchy provided
    combined_df = add_kegg_hierarchy(kegg_path, combined_df)
    
    if top_n_terms:
        combined_df = combined_df.groupby('Pair').apply(lambda x: x.nsmallest(top_n_terms, 'Adjusted P-value')).reset_index(drop=True)
        

    if save_table:
        combined_df.to_csv(results_path + exp_name + '.csv')

    return combined_df

def generate_results_tables(enrichment_results, exp_name, results_path, n_sets, top_n_terms=None, save_table=True):
    """
    Generate tables based on the enrichment information.

    Parameters:
    - enrichment_results (list): DataFrames with TF-IDF values and additional information for each gene set.
    - exp_name (str): Experiment name.
    - results_path (str): Path to save the results.
    - save_table (bool): Whether to save the table as a CSV file.
    - top_n_terms (int): Number of top terms to select for each latent dimension if rows > 40.

    Returns:
    - pd.DataFrame: Combined DataFrame used for plotting.
    """
    
    # Create a list of dataframes with an additional column for latent dimension
    dfs_with_dimension = []
    for gene_set, data in enrichment_results.items():
        data['Pair'] = gene_set
        dfs_with_dimension.append(data)

    # Combine all dataframes into a single dataframe
    combined_df = pd.concat(dfs_with_dimension, axis=0)
    combined_df = combined_df.reset_index(inplace=False)
    combined_df.drop(columns=combined_df.columns[0], axis=1, inplace=True)
    
    if top_n_terms:
        combined_df = combined_df.groupby('Pair').apply(lambda x: x.nsmallest(top_n_terms, 'Adjusted P-value')).reset_index(drop=True)
    
        
    if save_table:
        combined_df.to_csv(results_path + exp_name + '.csv')

    return combined_df




def pair_enrichment(protein_set, results_path, exp_name, enrich_genesets_dbs = 'KEGG_2021_Human', p_value_threshold = 0.05, 
                    top_n_terms=None,kegg_path = None, save_table=True):
    """
    Perform enrichment analysis for each cis standalone trans pair and generate and save tables.

    Parameters:
    - protein_set (list): List of protein sets.
    - enrich_genesets_dbs (list): List of gene sets or pathways of interest.
    - p_value_threshold (float): Adjusted p-value threshold for filtering.
    - exp_name (str): Experiment name.
    - results_path (str): Path to save the results.
    - save_table (bool): Whether to save the table as a CSV file.

    Returns:
    - pd.DataFrame: Combined DataFrame of the enrichments.
    """
    
    enrichment_results, n_sets = get_enrichment(protein_set, enrich_genesets_dbs = enrich_genesets_dbs , p_value_threshold=p_value_threshold)
    
    if enrich_genesets_dbs == 'KEGG_2021_Human':
        result_df = generate_KEGG_results_tables(enrichment_results, 
                                                          exp_name = exp_name, 
                                                          results_path = results_path, 
                                                          save_table = save_table, 
                                                          top_n_terms = top_n_terms, 
                                                          kegg_path = kegg_path,
                                                          n_sets = n_sets) 

    else:
        result_df =generate_results_tables(enrichment_results, 
                                           exp_name = exp_name, 
                                           results_path = results_path, 
                                           save_table = save_table, 
                                           top_n_terms = top_n_terms, 
                                           n_sets = n_sets) 
    

    return result_df

base_output_dir = "./results/Network_description/"
outcome_of_interest = 'DiabIncident'    
df = pd.read_table(f'{base_output_dir}{outcome_of_interest}/Summary_node_topology_{outcome_of_interest}.csv', sep=",")

protein_set=df.groupby('cfg_community')['Node'].apply(list).tolist()

exp_name= 'Comunities'
kegg_path='./kegg_hierarchy.csv'
enriched_df_kegg = pair_enrichment(protein_set, base_output_dir, exp_name, enrich_genesets_dbs = 'KEGG_2021_Human', p_value_threshold = 0.05, 
                    top_n_terms=None,kegg_path = kegg_path, save_table=True)


names = gp.get_library_name(organism='Human')
exp_name= 'Comunities_GO'
enriched_df_go = pair_enrichment(protein_set, base_output_dir, exp_name, enrich_genesets_dbs = 'GO_Biological_Process_2023', p_value_threshold = 0.05, 
                    top_n_terms=None,kegg_path = kegg_path, save_table=True)
exp_name = 'Comunities_Human_Phenotype_Ontology'
enriched_df_hpo = pair_enrichment(protein_set, base_output_dir, exp_name, enrich_genesets_dbs = 'Human_Phenotype_Ontology', p_value_threshold = 0.05, 
                    top_n_terms=None,kegg_path = kegg_path, save_table=True)

enriched_df_kegg["Subcategory_new"] = enriched_df_kegg["Subcategory"].copy()

recode = {
    "Signal transduction": "Signaling molecules and interaction",
    "Infectious disease: parasitic": "Infectious disease",
    "Infectious disease: bacterial": "Infectious disease",
    "Infectious disease: viral": "Infectious disease",
    "Development and regeneration": "Cell growth and death",
    "Amino acid metabolism": "Metabolism",
    "Energy metabolism": "Metabolism",
    "Glycan biosynthesis and metabolism": "Metabolism",
    "Carbohydrate metabolism": "Metabolism",
    "Cancer: overview": "Cancer",
    "Cancer: specific types": "Cancer",
}
enriched_df_kegg["Subcategory_new"] = enriched_df_kegg["Subcategory_new"].replace(recode)
enriched_df_kegg["Pair_num"] = pd.to_numeric(enriched_df_kegg["Pair"])+1
col_order = sorted(enriched_df_kegg["Pair_num"].dropna().unique())

mat = (enriched_df_kegg.groupby(["Subcategory_new", "Pair_num"])
         .size()
         .unstack(fill_value=0)
         .reindex(columns=col_order))
mat.columns = [str(int(c)) for c in mat.columns]
mat_plot = mat.clip(upper=5)

# -----------------------------
# Plot heatmap
# -----------------------------
plt.figure(figsize=(10, 5)) 
ax = sns.heatmap(
    mat_plot,
    cmap="Purples",
    vmin=0, vmax=5,
    linewidths=0.0, 
    cbar_kws={"ticks": [0, 1, 2, 3, 4, 5]}
)

ax.set_xlabel("Communities")
ax.set_ylabel("") 
ax.tick_params(axis="x", rotation=0)
ax.tick_params(axis="y", rotation=0)

plt.tight_layout()
plt.savefig(results_path + "kegg_subcategory_by_community_heatmap.png", dpi=300, bbox_inches="tight")
