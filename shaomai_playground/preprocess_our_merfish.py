"""
Adopted from preprocess_xenium_lung.py

Add gene name --> pyensembl conversion
"""

import scanpy as sc
import sys

import os 
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix

from nicheformer.data.constants import DefaultPaths, ObsConstants, UnsConstants, VarConstants, AssayOntologyTermId, SexOntologyTermId, OrganismOntologyTermId, TissueOntologyTermId, SuspensionTypeId
from nicheformer.data.validate import validate

sys.path.append("/home/ec2-user/SageMaker/functions")
from preprocessing import ensem_for_adata_mer

# --- some metadata ---
# setting MERFISH for now until an official ontology term is released
assay = str(AssayOntologyTermId.MERFISH_SPATIAL.value)
sex = str(SexOntologyTermId.UNKNOWN.value)
organism = str(OrganismOntologyTermId.MOUSE.value)
organism_validator = "mouse"
tissue = str(TissueOntologyTermId.BRAIN.value)
suspension_type = str(SuspensionTypeId.SPATIAL.value)
tissue_type = "brain"

def read_adata(file_path="./our_data_LC_NE/adata_mer.h5ad"):
    print(f"Read data from {file_path}...")
    adata = sc.read_h5ad(file_path)
    return adata

def add_metadata(adata):
    print(f"Adding metadata...")
    #adata.X = csr_matrix(adata.X)
    adata.obs[ObsConstants.SPATIAL_X] = adata.obs['z_CCF']
    adata.obs[ObsConstants.SPATIAL_Y] = adata.obs['y_CCF']
    adata.obs[ObsConstants.ASSAY_ONTOLOGY_TERM_ID] = pd.Categorical([assay for i in range(len(adata))])
    adata.obs[ObsConstants.SEX_ONTOLOGY_TERM_ID] = pd.Categorical([sex for i in range(len(adata))])
    adata.obs[ObsConstants.ORGANISM_ONTOLOGY_TERM_ID] = pd.Categorical([organism for i in range(len(adata))])
    adata.obs[ObsConstants.TISSUE_ONTOLOGY_TERM_ID] = pd.Categorical([tissue for i in range(len(adata))])
    adata.obs[ObsConstants.SUSPENSION_TYPE] = pd.Categorical([suspension_type for i in range(len(adata))])
    
    adata.obs[ObsConstants.DONOR_ID] = pd.Categorical(['Xenium_Preview_Human_Non_diseased_Lung_With_Add_on_FFPE_outs' for i in range(len(adata))])
    adata.obs[ObsConstants.CONDITION_ID] = pd.Categorical(['non diseased' for i in range(len(adata))])
    adata.obs[ObsConstants.TISSUE_TYPE] = pd.Categorical([tissue_type for i in range(len(adata))])
    
    adata.uns[UnsConstants.TITLE] = "Xenium_Preview_Human_Non_diseased_Lung_With_Add_on_FFPE_outs"
    adata.var[VarConstants.FEATURE_IS_FILTERED] = False
    
    adata.obs[ObsConstants.LIBRARY_KEY] = pd.Categorical(['section' for i in range(len(adata))])
    return adata

def preprocess(adata):
    # -- 1. Add some metadata --
    adata = add_metadata(adata)

    # -- 2. Turn gene name to ENSEMBL id --
    print("Turn gene name to ENSEMBL id...")
    adata = ensem_for_adata_mer(adata)
    
    # -- 3. Validate! --
    print("Validate data...")
    adata_output, valid, errors, is_seurat_convertible = validate(adata, organism=organism_validator)

    # -- 4. Some other stuff --
    adata_output.obs['assay'] = pd.Categorical(['Xenium' for i in range(len(adata_output))])
    adata_output.obs[ObsConstants.ASSAY_ONTOLOGY_TERM_ID] = pd.Categorical(['no yet defined' for i in range(len(adata_output))])
    
    adata_output.obs[ObsConstants.DATASET] = adata_output.uns['title']
    adata_output.obs[ObsConstants.SPLIT] = 'train'
    adata_output.obs[ObsConstants.NICHE] = 'nan'
    adata_output.obs[ObsConstants.REGION] = 'nan'
    return adata_output

def write_processed(adata_output, output_path="./our_data_LC_NE/adata_mer_preprocessed.h5ad"):
    adata_output.write(output_path)
    print(f"Processed data saved to {output_path}!")

if __name__ == "__main__":
    adata = read_adata(file_path="./our_data_LC_NE/adata_mer.h5ad")
    adata_output = preprocess(adata)
    write_processed(adata_output, output_path="./our_data_LC_NE/adata_mer_preprocessed.h5ad")

