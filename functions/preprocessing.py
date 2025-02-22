import numpy as np
from pyensembl import EnsemblRelease


def ensem_for_adata_mer(adata):         
    data = EnsemblRelease(104, species='mus_musculus')    
    data.index()
    new_names = []
    counter = 0
    for i,m in enumerate(adata.var_names):
        try:
            out=data.genes_by_name(m)[0].id
            new_names.append(out)
        except:
            try:
                m1 = m.split(' ')[0]
                out=data.genes_by_name(m1)[0].id
                new_names.append(out)
            except:
                new_names.append(data.genes_by_name('Gnai3')[0].id)
                counter +=1            
    print(f'total invalid {counter}, from {len(adata.var_names)} {100*counter/len(adata.var_names)} %')    
    adata_1 = adata.copy()    
    print(adata_1.shape)
    adata_1.var.index = (new_names)    
    adata_1.raw.var.index = (new_names)
    
    print(f'new data size {adata_1.shape}')
    return(adata_1)



def ensem_for_adata(adata):         
    data = EnsemblRelease(104, species='mus_musculus')    
    data.index()
    new_names = []
    counter = 0
    for i,m in enumerate(adata.var_names):
        try:
            out=data.genes_by_name(m)[0].id
            new_names.append(out)
        except e:
            print(e)
            try:
                m1 = m.split(' ')[0]
                out=data.genes_by_name(m1)[0].id
                new_names.append(out)
            except e:
                print(e)
                counter =+1
                new_names.append(f'NA_{counter}')        
    print(f'total invalid {counter}, from {len(adata.var_names)} {100*counter/len(adata.var_names)} %')    
    adata_1 = adata.copy()    
    adata_1.var.index = (new_names)    
    adata_1 = adata_1[:,adata_1.var.index != 'NA']
    
    print(f'new data size {adata_1.shape}')
    return(adata_1)