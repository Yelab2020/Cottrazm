
'''
This script is for ST data SME modified
'''
#!/usr/bin/env python
# conding:utf-8 -*-
import argparse
import stlearn as st
import scanpy as sc
import numpy as np
from numpy import random,mat
from pathlib import Path
import pandas as pd
from scipy import io,sparse
import os

def ME_normalize(inDir,outDir,sample):
    print (sample, "start SME normalize")

    #read data
    data=st.Read10X(path = inDir)
    data.var_names_make_unique()
    data.layers['raw_count']=data.X
    #tile data
    TILE_PATH=Path(os.path.join(outDir,'{0}_tile'.format(sample)))
    TILE_PATH.mkdir(parents=True,exist_ok=True)
    
    #tile morphology
    st.pp.tiling(data,TILE_PATH,crop_size=40)
    st.pp.extract_feature(data)

    ###process data
    st.pp.normalize_total(data)
    st.pp.log1p(data)

    #gene pca dimention reduction
    st.em.run_pca(data,n_comps=50,random_state=0)
    
    #stSME to normalise log transformed data
    st.spatial.SME.SME_normalize(data, use_data="raw",weights =  "weights_matrix_gd_md")

    #convert SME_norm data to sparesmatrix
    raw_SME_normalized = mat(data.obsm['raw_SME_normalized'])
    raw_SME_normalizedA = sparse.csr_matrix(raw_SME_normalized)
    print ("matrix conver ok!")
    
    io.mmwrite(os.path.join(outDir,'{0}_raw_SME_normalizeA.mtx'.format(sample)),raw_SME_normalizedA)
    print("Morphology adjusted is ok!")

    return raw_SME_normalizedA
   
