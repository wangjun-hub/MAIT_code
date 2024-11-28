####Figure 5F

import scanpy as sc
import sctour as sct
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

pbmc_csv = sc.read_csv("/Users/wangjun/Desktop/FluT_MAIT/MAIT_final/mait_df_t.csv")

sc.pp.calculate_qc_metrics(pbmc_csv, percent_top=None, log1p=False, inplace=True)

tnode = sct.train.Trainer(pbmc_csv, percent=0.6,nepoch=100)
tnode.train()

pbmc_csv.obs['ptime'] = tnode.get_time()

##latent space
mix_zs, zs, pred_zs = tnode.get_latentsp(alpha_z=0.5, alpha_predz=0.5)
pbmc_csv.obsm['X_TNODE'] = mix_zs

##vector field
pbmc_csv.obsm['X_VF'] = tnode.get_vector_field(pbmc_csv.obs['ptime'].values, pbmc_csv.obsm['X_TNODE'])

pbmc_csv.obsm['X_tsne'] = pbmc.obsm['X_tsne'].copy()


sct.vf.plot_vector_field(pbmc_csv, zs_key='X_TNODE',E_key = "X_tsne", vf_key='X_VF', use_rep_neigh='X_TNODE', show=False, t_key='ptime', frameon=False, size=80, alpha=0.05)

####Read h5ad data

pbmc = sc.read_h5ad("/Users/wangjun/Desktop/FluT_MAIT/MAIT_final/MAIT_anno_6_23.h5ad")

pbmc_csv.obsm['X_rnatsne'] = pbmc.obsm['X_rnatsne'].copy()


##Save model
tnode.save_model(save_dir='/Users/wangjun/Desktop/FluT_MAIT/MAIT_final/', save_prefix='sctour_model_mast')

tnode = sct.predict.load_model('/Users/wangjun/Desktop/FluT_MAIT/MAIT_final/sctour_model.pth')


pbmc_csv.obs["Cluster"] = pbmc.obs['Sample'].copy()


pbmc_csv.obs['ptime'] = sct.train.reverse_time(pbmc_csv.obs['ptime'].values)


plt.figure(dpi=1000)
sct.vf.plot_vector_field(pbmc_csv,reverse=True, zs_key='X_TNODE',E_key = "X_rnatsne", vf_key='X_VF', use_rep_neigh='X_TNODE',color='ptime', t_key='ptime', frameon=True, size=80, alpha=0.2,show=False)  
plt.savefig('/Users/wangjun/Desktop/FluT_MAIT/MAIT_final/MAIT_sctour_ptime_6_23.png')
plt.close()

