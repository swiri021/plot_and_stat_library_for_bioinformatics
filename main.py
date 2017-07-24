from EXC_READ import EXC_READ
from DATA_FILTER import DATA_FILTER
from MY_STATS import MY_STATS
from MY_PLOT import MY_PLOT
import numpy as np
from sys import argv

import pandas as pd
import matplotlib
from matplotlib import pyplot as plt
matplotlib.use('Agg')
from matplotlib_venn import venn3, venn3_circles

fi = DATA_FILTER()
st = MY_STATS()
ex = EXC_READ()
pt = MY_PLOT()

#ex.wfile_open(argv[1],output=2)

"""
######histogram
#phospho_fold.tsv
#gdata, ghead, gindex = ex.file_read("plot_data/global_fold.tsv",tstatus="file")
gdata, ghead, gindex = ex.file_read("plot_data/phospho_fold.tsv",tstatus="file")
df = ex.pandas_data(gdata,ghead, gindex)

#ind_list = [df[x].dropna().index.tolist() for x in df.columns.tolist()]
ghead_c = ['R0/C0', 'C10', 'R10']
#ghead_c = ['R0/C0', 'R10', 'R30']
df_c = df[ghead_c]
pt.histogram_group(df_c,colo=['Blue', 'Red', 'Green'])
######histogram
"""



"""
######volcano plot
gdata, ghead, gindex = ex.file_read(argv[1],tstatus="file")
df1 = ex.pandas_data(gdata,ghead, gindex)

gdata, ghead, gindex = ex.file_read(argv[2],tstatus="file")
df2 = ex.pandas_data(gdata,ghead, gindex)

pt.volcano_plot(df1,df2,x_ax=df1.columns.tolist(), y_ax=df2.columns.tolist(), fn=argv[3])
######volcano plot
"""

"""
######heatmap

#gdata, ghead, gindex = ex.file_read("dep_dpp_data/global_dep_set_c0.tsv",tstatus="file", ind_loc=1)
#df1 = ex.pandas_data(gdata,ghead, gindex)
gdata, ghead, gindex = ex.file_read("dep_dpp_data/global_dep_set_c10.tsv",tstatus="file", ind_loc=1)
df2 = ex.pandas_data(gdata,ghead, gindex)
gdata, ghead, gindex = ex.file_read("dep_dpp_data/global_dep_set_c30.tsv",tstatus="file", ind_loc=1)
df3 = ex.pandas_data(gdata,ghead, gindex)
gdata, ghead, gindex = ex.file_read("dep_dpp_data/global_dep_set_r10.tsv",tstatus="file", ind_loc=1)
df4 = ex.pandas_data(gdata,ghead, gindex)
gdata, ghead, gindex = ex.file_read("dep_dpp_data/global_dep_set_r30.tsv",tstatus="file", ind_loc=1)
df5 = ex.pandas_data(gdata,ghead, gindex)

gdata, ghead, gindex = ex.file_read("dep_dpp_data/phospho_dpp_set_c10.tsv",tstatus="file")
df2 = ex.pandas_data(gdata,ghead, gindex)
gdata, ghead, gindex = ex.file_read("dep_dpp_data/phospho_dpp_set_c30.tsv",tstatus="file")
df3 = ex.pandas_data(gdata,ghead, gindex)
gdata, ghead, gindex = ex.file_read("dep_dpp_data/phospho_dpp_set_r10.tsv",tstatus="file")
df4 = ex.pandas_data(gdata,ghead, gindex)
gdata, ghead, gindex = ex.file_read("dep_dpp_data/phospho_dpp_set_r30.tsv",tstatus="file")
df5 = ex.pandas_data(gdata,ghead, gindex)

df2 = fi.pandas_if_change(df2,oper='>',comp_numb=2,change_numb=2)
df2 = fi.pandas_if_change(df2,oper='<',comp_numb=-2,change_numb=-2)

df3 = fi.pandas_if_change(df3,oper='>',comp_numb=2,change_numb=2)
df3 = fi.pandas_if_change(df3,oper='<',comp_numb=-2,change_numb=-2)

df4 = fi.pandas_if_change(df4,oper='>',comp_numb=2,change_numb=2)
df4 = fi.pandas_if_change(df4,oper='<',comp_numb=-2,change_numb=-2)

df5 = fi.pandas_if_change(df5,oper='>',comp_numb=2,change_numb=2)
df5 = fi.pandas_if_change(df5,oper='<',comp_numb=-2,change_numb=-2)

uni_df = pd.concat([df2,df3,df4,df5], axis=1)
uni_df = uni_df.fillna(0)

pt.cluter_heatmap(uni_df)
"""
"""
######rank-plot
gdata, ghead, gindex = ex.file_read(argv[1],tstatus="file")

df1 = ex.pandas_data(gdata,ghead, gindex)
selected = pd.DataFrame(df1['Res10'].dropna(), columns=['Res10'])
converted = fi.id_conversion(selected,'all_data/phospho_entrez_list0.tsv', original_id_type='float' ,converted_id_type='float')

df_c = converted.columns.tolist()
pt.rankplot(converted,df_c, selected_on=True, selected='erbb_signal_kegg.edit.txt')
#####rank-plot
"""
#ex.wfile_open(argv[1],output=2)

"""
#####scatter and reg_line
#gdata, ghead, gindex = ex.file_read(argv[1],tstatus="file", ind_loc=1, index_type='str')
gdata, ghead, gindex = ex.file_read(argv[1],tstatus="file", index_type='float')
df1 = ex.pandas_data(gdata,ghead, gindex)
selected = df1[['Con0','Res0']]
#print selected

#converted = fi.id_conversion(selected,'all_data/global_entrez_list0.tsv', converted_index_type = 'str', converted_id_type='float')
converted = fi.id_conversion(selected,'all_data/phospho_entrez_list0.tsv', converted_index_type = 'float', converted_id_type='float')
#converted = fi.id_conversion(selected,'all_data/phospho_entrez_list0.tsv', original_id_type='float' ,converted_id_type='float')

selected_set = fi.set_intersect(converted, 'erbb_signal_kegg.edit.txt', set_id_loc=1)

#pt.scatter_plot_and_line(selected_set, x_dat='Con30', y_dat='Res30', annotate=True, annotation_file='all_data/global_entrez_symbol0.tsv', line_sep='\r\n', filename='c10_r10_scatter_g')
#pt.scatter_plot_and_line(selected_set, x_dat='Con10', y_dat='Res10', annotate=False, filename='c10_r10_scatter')
pt.scatter_plot_and_line(selected_set, x_dat='Con0', y_dat='Res0', annotate=True, annotation_file='all_data/phospho_entrez_symbol0.tsv', line_sep='\r\n', filename='c0_r0_scatter_g', conv_id_loc=4, conv_dat_loc=3)
#####scatter and reg_line
"""

"""
#####scatter and reg_line
gdata, ghead, gindex = ex.file_read(argv[1],tstatus="file", ind_loc=1)
df1 = ex.pandas_data(gdata,ghead, gindex)
selected = df1[['Con0','Con30']]
controls = selected['Con30']/selected['Con0']

controls = pd.DataFrame(controls.values.tolist(), columns=['Con30/Con0'], index=controls.index.tolist())
converted_controls = fi.id_conversion(controls,'all_data/global_entrez_list0.tsv')

selected = df1[['Res0','Res30']]
resist = selected['Res30']/selected['Res0']

resist = pd.DataFrame(resist.values.tolist(), columns=['Res30/Res0'], index=resist.index.tolist())
converted_resist = fi.id_conversion(resist,'all_data/global_entrez_list0.tsv')

controls_set = fi.set_intersect(converted_controls, 'erbb_signal_kegg.edit.txt', set_id_loc=1)
resist_set = fi.set_intersect(converted_resist, 'erbb_signal_kegg.edit.txt', set_id_loc=1)

uni_set = pd.concat([controls_set, resist_set], axis=1)

pt.scatter_plot_and_line(uni_set['Con30/Con0'], uni_set['Res30/Res0'])

#####scatter and reg_line
"""

"""
###new heatmap
from scipy.cluster.hierarchy import dendrogram, linkage

#gdata, ghead, gindex = ex.file_read("dep_dpp_data/global_dep_set_c0.tsv",tstatus="file", ind_loc=1)
#df1 = ex.pandas_data(gdata,ghead, gindex)
gdata, ghead, gindex = ex.file_read("dep_dpp_data/global_dep_set_c10.tsv",tstatus="file", ind_loc=1)
df2 = ex.pandas_data(gdata,ghead, gindex)
gdata, ghead, gindex = ex.file_read("dep_dpp_data/global_dep_set_c30.tsv",tstatus="file", ind_loc=1)
df3 = ex.pandas_data(gdata,ghead, gindex)
gdata, ghead, gindex = ex.file_read("dep_dpp_data/global_dep_set_r10.tsv",tstatus="file", ind_loc=1)
df4 = ex.pandas_data(gdata,ghead, gindex)
gdata, ghead, gindex = ex.file_read("dep_dpp_data/global_dep_set_r30.tsv",tstatus="file", ind_loc=1)
df5 = ex.pandas_data(gdata,ghead, gindex)

df2 = fi.pandas_if_change(df2,oper='>',comp_numb=2,change_numb=2)
df2 = fi.pandas_if_change(df2,oper='<',comp_numb=-2,change_numb=-2)

df3 = fi.pandas_if_change(df3,oper='>',comp_numb=2,change_numb=2)
df3 = fi.pandas_if_change(df3,oper='<',comp_numb=-2,change_numb=-2)

df4 = fi.pandas_if_change(df4,oper='>',comp_numb=2,change_numb=2)
df4 = fi.pandas_if_change(df4,oper='<',comp_numb=-2,change_numb=-2)

df5 = fi.pandas_if_change(df5,oper='>',comp_numb=2,change_numb=2)
df5 = fi.pandas_if_change(df5,oper='<',comp_numb=-2,change_numb=-2)

uni_df1 = pd.concat([df2,df3], axis=1)
uni_df2 = pd.concat([df4,df5], axis=1)
uni_df1 = uni_df1.fillna(0)
uni_df2 = uni_df2.fillna(0)

uni_df1_div = fi.pandas_cust_div(uni_df1)
uni_df2_div = fi.pandas_cust_div(uni_df2)

#pt.hierarchical_dendro_and_heatmap(uni_df1_div[0], 'dep_c10_c30_1set')
#pt.hierarchical_dendro_and_heatmap(uni_df1_div[1], 'dep_c10_c30_2set')
#pt.hierarchical_dendro_and_heatmap(uni_df1_div[2], 'dep_c10_c30_3set')


#pt.hierarchical_dendro_and_heatmap(uni_df2_div[0], 'dep_r10_r30_1set')
#pt.hierarchical_dendro_and_heatmap(uni_df2_div[1], 'dep_r10_r30_2set')
pt.hierarchical_dendro_and_heatmap(uni_df2_div[2], 'dep_r10_r30_3set')
"""

"""
######intersection, venn diagram for heatmap section
gdata, ghead, gindex = ex.file_read("dep_dpp_data/global_dep_set_c0.tsv",tstatus="file", ind_loc=1)
df1 = ex.pandas_data(gdata,ghead, gindex)
gdata, ghead, gindex = ex.file_read("dep_dpp_data/global_dep_set_c10.tsv",tstatus="file", ind_loc=1)
df2 = ex.pandas_data(gdata,ghead, gindex)
gdata, ghead, gindex = ex.file_read("dep_dpp_data/global_dep_set_c30.tsv",tstatus="file", ind_loc=1)
df3 = ex.pandas_data(gdata,ghead, gindex)
gdata, ghead, gindex = ex.file_read("dep_dpp_data/global_dep_set_r10.tsv",tstatus="file", ind_loc=1)
df4 = ex.pandas_data(gdata,ghead, gindex)
gdata, ghead, gindex = ex.file_read("dep_dpp_data/global_dep_set_r30.tsv",tstatus="file", ind_loc=1)
df5 = ex.pandas_data(gdata,ghead, gindex)

c_df = pd.concat([df2,df3],axis=1)
r_df = pd.concat([df4,df5],axis=1)

#tindex = df1.dropna().index.tolist()
#tdat = np.random.rand(len(tindex),1)
gc_index = list(set(c_df.index.tolist()))
gcdat = []
[gcdat.append(np.random.normal()) for x in range(len(gc_index))]

gr_index = list(set(r_df.index.tolist()))
#grdat = np.random.rand(len(gr_index),1)
grdat = []
[grdat.append(np.random.normal()) for x in range(len(gr_index))]

gc_df = pd.DataFrame(gcdat, index=gc_index, columns=['Control'])
gr_df = pd.DataFrame(grdat, index=gr_index, columns=['Resist'])

n_df = pd.concat([df1,gc_df,gr_df], axis=1)
#print n_df

pt.venn_dia(n_df, scaling='log2')

"""