

Old gene list as a dictionary...
```
genes_dict = {
    'iPSC':['POUSR1','POU5F1','SOX2','NANOG','ST3GAL2','CDH1','SSEA4'], #EPCAM
    'Gastr_Cells_Prim_Streak': ['MIXL1','CER1','EOMES','SOX17','NODAL','APLNR','TBX6','HES7','TBXT'],
    'Posterior_PS':['ETV2','EVX1','ENG','TGFBR1'],
    'Mesendoderm_Prim_Endoderm':['SOX17','GATA4','GATA6','SOX7','FOXA2'],
    'Distal_Foregut_Pancreas':['PDX1','HHEX','MNX1','SOX9'],#'HNF6',
    'Posterior_Foregut_Liver':['ALB','CDH6','AFP','ONECUT2','HNF4A','FOXA1','FOXA2'],
    'Gut_Mesenchyme':['ALCAM','CDH1','CXCL12','PITX2','LHX2','MEX2','DES'], #'NCAM1','PBX1'
    'Endoderm_Gut_Cells':['ALB','CDH6','AFP','ONECUT2','HNF4A'],
    'Gut_Silva':['LGALS4','TSPAN8','MT1M'],
    'Stromal_Fibroblasts':['VIM','TAGLN2','COL1A1','COL3A1','TWIST1','PRRX1','THY1','POSTN','DCN','SPARC','HAND1'],
    'Mural_Cell':['RGS5','TAGLN','ANKRD1'],
    'Smooth_Muscle':['ACTA2','ACTG2','MYL9','MYH11','TPM2','TAGLN2','CNN1','SMTN','CALD1','NPPB'],#'SM22',
    'Cardiac_Mesoderm':['EOMES','ISL1','MESP1','MESP2'],#,'LHX1'
    'Endocardial':['TEK','NPC2','NPR3','NFATC1','PECAM1','CDH5','CD34','KDR','GATA5','TBX3'],#'TM4SF1',
    'Epicardial':['WT1','TBX18','TCF21','NR2F2'],
    'Cardiomyocytes':['ACTN2','MYH6','MYH7','MYL4','MYL7','TTN','NCAM1','PBX1'],#'TNNT2',
    'First_Heart_Field':['TBX5','HCN4','HAND1'],
    'Second_Heart_Field':['HAND2','TBX1','CXCR4','TBX20'],
    'Melanocyte':['MLANA','DCT','TFAP2A'],
    'Neural': ['ALDH1A1','ISL1','STMN2', 'TH', 'SLC18A3']
    'DEGs_old':['TGFBI','ANXA1','FBLN1','FOS','WNT6','GATA3','IGFBP7','CNTNAP2','GAL'], #'IGFB',
    'DEGs_new':[
        'CDH11','TMEM88','PCDH11X','LINC00698','CD9','TMSB4X','CDH2','DKK1','CPED1','CXCL14',
                'TECRL','NRXN3','APOA1','GNG11','GABRP','ARHGAP29','KCNMA1','PTN','TMTC2','TIMP1','KCNQ5'
    ]
}

```

"genes_slim"
```
genes_slim = {
    'iPSC':['POU5F1','SOX2','NANOG','ST3GAL2','CDH1'], 
    'Prim_Streak': ['MIXL1','TBXT','EOMES','EVX1','NODAL'],
    'Cardiac_Mesoderm':['MESP1','ISL1','APLNR','LHX1'],
    'Definitive_Endoderm':['SOX17','EPCAM','GATA6','FOXA2'],
    'Smooth_Muscle':['ACTA2','TAGLN','SMTN','CNN1'],
    'Stromal_Fibroblasts':['VIM','TCF21','TBX6','NPC2','COL1A1'],
    'Epicardial':['LHX2','WT1','TBX18'],
    'Endocardial':['ENG','KDR','CDH5','PECAM1','CD34'],
    'Cardiomyocytes':['MYH6','MYH7','MYL2','MYL4','TBX5','NKX2-5','HAND1','HAND2'],
    'Foregut_Epithelium':['HHEX','SOX9','CDH6','ONECUT2','FOXA1'],
    'Gut_Mesenchyme':['NPPB','ALCAM','PITX2','PBX1','APOA2'],
    'Liver_Progenitor':['ALB','AFP','HNF4A','ICAM1','NCAM1']    
}

# Filter out genes that were not detected
for CT in genes_slim.keys():
    genes_slim[CT] = [k for k in genes_slim[CT] if k in adata.var_names]

genes_slim
```

# Notes on gene list(s):
### Genes missing in adata:
SOX2
SNAI2
TBXT
SOX17

### Genes changed from original list:
MIXL2 -> MIXL1
