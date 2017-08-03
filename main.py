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

gdata, ghead, gindex = ex.file_read("volcano_data/global_c10_fold.tsv",tstatus="file",index_type='str', ind_loc=1)
df1 = ex.pandas_data(gdata,ghead, gindex)

gdata, ghead, gindex = ex.file_read("volcano_data/global_c10_fdr.tsv",tstatus="file", index_type='str', ind_loc=1)
df2 = ex.pandas_data(gdata,ghead, gindex)

union_df = pd.concat([df1,df2], axis=1)
union_df.columns = ['x','y']

pt.volcano_plot(union_df, x_ax='x', y_ax='y', fn='c10_volcano')
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
#gdata, ghead, gindex = ex.file_read("dep_dpp_data/global_dep_set_c0.tsv",tstatus="file", ind_loc=1)
#df1 = ex.pandas_data(gdata,ghead, gindex)
#gdata, ghead, gindex = ex.file_read("dep_dpp_data/global_dep_set_c10.tsv",tstatus="file", ind_loc=1, index_type='float')
#df2 = ex.pandas_data(gdata,ghead, gindex)
#gdata, ghead, gindex = ex.file_read("dep_dpp_data/global_dep_set_c30.tsv",tstatus="file", ind_loc=1, index_type='float')
#df3 = ex.pandas_data(gdata,ghead, gindex)
#gdata, ghead, gindex = ex.file_read("dep_dpp_data/global_dep_set_r10.tsv",tstatus="file", ind_loc=1, index_type='float')
#df4 = ex.pandas_data(gdata,ghead, gindex)
#gdata, ghead, gindex = ex.file_read("dep_dpp_data/global_dep_set_r30.tsv",tstatus="file", ind_loc=1, index_type='float')
#df5 = ex.pandas_data(gdata,ghead, gindex)

gdata, ghead, gindex = ex.file_read("dep_dpp_data/phospho_dpp_set_c10.tsv",tstatus="file", index_type='float')
df2 = ex.pandas_data(gdata,ghead, gindex)
gdata, ghead, gindex = ex.file_read("dep_dpp_data/phospho_dpp_set_c30.tsv",tstatus="file", index_type='float')
df3 = ex.pandas_data(gdata,ghead, gindex)
gdata, ghead, gindex = ex.file_read("dep_dpp_data/phospho_dpp_set_r10.tsv",tstatus="file", index_type='float')
df4 = ex.pandas_data(gdata,ghead, gindex)
gdata, ghead, gindex = ex.file_read("dep_dpp_data/phospho_dpp_set_r30.tsv",tstatus="file", index_type='float')
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

cc = 0
for x in uni_df1_div:
	converted_df1 = fi.id_conversion(x,'all_data/phospho_entrez_list0.tsv', original_id_type='float', converted_index_type='float', converted_id_type='float')

	print "control "+str(cc)+" :"
	for y in converted_df1.index.tolist():
		print y
	print "##########\n"
	cc+=1

cc = 0
for x in uni_df2_div:
	converted_df1 = fi.id_conversion(x,'all_data/phospho_entrez_list0.tsv', original_id_type='float', converted_index_type='float', converted_id_type='float')

	print "Resist "+str(cc)+" :"
	for y in converted_df1.index.tolist():
		print y
	print "##########\n"
	cc+=1

#pt.hierarchical_dendro_and_heatmap(uni_df1_div[0], 'dep_c10_c30_1set', dend_line_width=20.0)
#pt.hierarchical_dendro_and_heatmap(uni_df1_div[1], 'dep_c10_c30_2set', dend_line_width=20.0)
#pt.hierarchical_dendro_and_heatmap(uni_df1_div[2], 'dep_c10_c30_3set', dend_line_width=20.0)
#pt.hierarchical_dendro_and_heatmap(uni_df2_div[0], 'dep_r10_r30_1set', dend_line_width=20.0)
#pt.hierarchical_dendro_and_heatmap(uni_df2_div[1], 'dep_r10_r30_2set', dend_line_width=20.0)
#pt.hierarchical_dendro_and_heatmap(uni_df2_div[2], 'dep_r10_r30_3set', dend_line_width=20.0)

pt.hierarchical_dendro_and_heatmap(uni_df1_div[0], 'dpp_c10_c30_1set', dend_line_width=5.0, heatmap_line_width=0.0)
#pt.hierarchical_dendro_and_heatmap(uni_df1_div[1], 'dpp_c10_c30_2set', dend_line_width=5.0, heatmap_line_width=0.0)
#pt.hierarchical_dendro_and_heatmap(uni_df1_div[2], 'dpp_c10_c30_3set', dend_line_width=5.0, heatmap_line_width=0.0)
#pt.hierarchical_dendro_and_heatmap(uni_df2_div[0], 'dpp_r10_r30_1set', dend_line_width=10.0, heatmap_line_width=0.0)
#pt.hierarchical_dendro_and_heatmap(uni_df2_div[1], 'dpp_r10_r30_2set', dend_line_width=5.0, heatmap_line_width=0.0)
#pt.hierarchical_dendro_and_heatmap(uni_df2_div[2], 'dpp_r10_r30_3set', dend_line_width=5.0, heatmap_line_width=0.0)
"""

"""
######intersection, venn diagram for heatmap section
#gdata, ghead, gindex = ex.file_read("dep_dpp_data/global_dep_set_c0.tsv",tstatus="file", ind_loc=1)
#df1 = ex.pandas_data(gdata,ghead, gindex)
#gdata, ghead, gindex = ex.file_read("dep_dpp_data/global_dep_set_c10.tsv",tstatus="file", ind_loc=1)
#df2 = ex.pandas_data(gdata,ghead, gindex)
#gdata, ghead, gindex = ex.file_read("dep_dpp_data/global_dep_set_c30.tsv",tstatus="file", ind_loc=1)
#df3 = ex.pandas_data(gdata,ghead, gindex)
#gdata, ghead, gindex = ex.file_read("dep_dpp_data/global_dep_set_r10.tsv",tstatus="file", ind_loc=1)
#df4 = ex.pandas_data(gdata,ghead, gindex)
#gdata, ghead, gindex = ex.file_read("dep_dpp_data/global_dep_set_r30.tsv",tstatus="file", ind_loc=1)
#df5 = ex.pandas_data(gdata,ghead, gindex)

gdata, ghead, gindex = ex.file_read("dep_dpp_data/phospho_dpp_set_c0.tsv",tstatus="file", index_type='float')
df1 = ex.pandas_data(gdata,ghead, gindex)
gdata, ghead, gindex = ex.file_read("dep_dpp_data/phospho_dpp_set_c10.tsv",tstatus="file", index_type='float')
df2 = ex.pandas_data(gdata,ghead, gindex)
gdata, ghead, gindex = ex.file_read("dep_dpp_data/phospho_dpp_set_c30.tsv",tstatus="file", index_type='float')
df3 = ex.pandas_data(gdata,ghead, gindex)
gdata, ghead, gindex = ex.file_read("dep_dpp_data/phospho_dpp_set_r10.tsv",tstatus="file", index_type='float')
df4 = ex.pandas_data(gdata,ghead, gindex)
gdata, ghead, gindex = ex.file_read("dep_dpp_data/phospho_dpp_set_r30.tsv",tstatus="file", index_type='float')
df5 = ex.pandas_data(gdata,ghead, gindex)

#c_df = pd.concat([df2,df3],axis=1)
#r_df = pd.concat([df4,df5],axis=1)

c_df = df3
r_df = df5

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

"""
###senario heatmap

def dup_peptide2protein(df1):

	new_values = []

	df_temp = df1
	df_temp_index = df_temp.index.tolist()
	dedup_list = list(set(df_temp_index))

	for y in dedup_list:
		temp = df_temp.loc[y].values.tolist()

		if type(temp[0]) is list:
			temp_fold = [x[0] for x in temp]#####fold
			temp_fold_value = max(temp_fold)

		else:
			temp_fold_value = temp[0]
		new_values.append(temp_fold_value)
	res_df1 = pd.DataFrame(data=new_values, index=dedup_list, columns=df1.columns.tolist())
	return res_df1

def toentrez(df):
	###annotation
	f = open("all_data/phospho_entrez_list0.tsv","r")
	dat = f.read().strip()
	dat = [x.split("\t") for x in dat.split("\r\n")]
	dat = dat[1:]

	conv_index = []

	for x in df.index.tolist():
		try:
			temp=[str(int(float(y[1]))) for y in dat if str(int(float(y[0])))==x]
		except ValueError:
			temp=[str(int(float(y[1]))) for y in dat if str(y[0])==x]

		conv_index.append(temp[0])
	df.index = conv_index
	###annotation

	return df

####global
#gdata, ghead, gindex = ex.file_read("dep_dpp_data/global_dep_set_c0.tsv",tstatus="file", ind_loc=1)#####value test
#df0 = ex.pandas_data(gdata,ghead, gindex)
#gdata, ghead, gindex = ex.file_read("all_data/global_list0.tsv",tstatus="file", ind_loc=1)
#df1 = ex.pandas_data(gdata,ghead, gindex)
#converted_df1 = fi.id_conversion(df1,'all_data/global_entrez_list0.tsv', original_id_type='str', converted_index_type='str', converted_id_type='float')
####global


####phospho
gdata, ghead, gindex = ex.file_read("all_data/phospho_list0.tsv",tstatus="file", index_type='float')
df1 = ex.pandas_data(gdata,ghead, gindex)
converted_df1 = df1
#converted_df1 = df1.dropna()

####phospho

####global
#gene_list_30_o = ['85349.0', '213.0', '7296.0', '833.0', '6407.0', '754.0', '735.0', '3491.0', '3881.0', '2934.0', '30011.0', '2222.0', '5837.0', '64710.0', '80135.0', '24137.0', '6372.0', '3988.0', '6319.0', '6038.0', '5055.0', '114805.0', '149499.0', '5834.0', '2069.0', '7298.0', '231.0', '84267.0', '10577.0', '2919.0', '65084.0', '2495.0', '2512.0', '26166.0', '3872.0', '7018.0']
#gene_list_10_o = ['3005.0', '213.0', '85349.0', '6407.0', '85236.0', '3881.0', '114805.0', '81539.0', '64710.0', '3020.0', '80135.0', '10144.0', '10397.0', '92840.0', '6038.0', '2947.0', '2934.0', '149499.0', '6273.0', '9555.0', '55336.0', '65084.0', '7114.0', '26166.0', '3872.0']
####global

####phospho
gene_list_10_o = ['1140', '1303', '16765', '1265', '12393', '13736', '10216', '8303', '13404', '3411', '4787', '14074', '14263', '7303', '4149', '16843', '1078', '1077', '14173', '3126', '10581', '13488', '16163', '7758', '343', '5651', '814', '5921', '14668', '12930', '1093', '5420', '15083', '15826', '10283', '16266', '12069', '16265', '15119', '125', '11577', '10113', '8024', '10195', '6568', '13604', '6390', '1821', '9462', '10804', '10805', '13257', '1663', '377', '1087', '17025', '4876', '2473', '13693', '11803', '15625', '15623', '11804', '11134', '16705', '9491', '119', '12586', '11032', '6263', '3473', '12014', '12017', '12011', '12010', '12013', '2984', '7889', '2079', '14921', '14296', '12155', '1395', '7730', '7761', '2607', '1232', '1233', '2936', '7214', '3668', '3318', '15139', '9648', '9888', '7731', '9481', '7629', '7453', '4813', '3509', '780', '5402', '2273', '12004', '11511', '1931', '35', '4058', '16393', '16390', '13136', '6073', '12121', '12122', '11294', '13235', '5561', '2610', '2292', '9534', '11942', '456', '258', '13674', '1906', '12731', '10257', '10256', '2154', '1427', '6249', '2889', '4822', '11500', '5838', '8311', '10621', '8602', '3897', '6365', '10319', '9', '757', '9971', '14024', '8371', '8067', '8717', '15711', '3071', '5311', '9873', '10398', '16421', '4375', '10', '2094', '5401', '2925', '14914', '4071', '2890', '3443', '5422', '15819', '16944', '1749', '2339', '1744', '12340', '6253', '2355', '15933', '14871', '15150', '9745', '278', '1314', '15785', '15786', '9958', '14497', '16623', '14249', '6887', '6345', '6347', '7200', '10082', '13265', '11299', '5920', '6180', '6247', '5219', '957', '595', '16756', '9881', '6995', '2720', '12806', '7234', '14563', '7233', '5049', '17252', '3057', '11993', '14166', '11497', '3544', '12681', '13579', '3541', '1044', '12688', '12410', '13826', '11160', '4662', '11251', '9887', '1471', '1476', '430', '16299', '9560', '1478', '10779', '13174', '683', '8641']
#gene_list_30_o =['10299', '14173', '14914', '14915', '13736', '8303', '10455', '14074', '4149', '10581', '13168', '13245', '814', '2858', '12930', '12931', '14296', '13604', '5917', '1087', '595', '8602', '1122', '192', '13544', '11136', '11134', '16705', '15107', '14816', '12014', '12017', '12011', '12010', '12013', '2984', '15580', '17252', '2441', '1233', '1335', '1330', '9841', '16954', '9481', '14522', '4420', '35', '6071', '1530', '12121', '12122', '430', '2610', '3126', '3806', '11942', '3802', '1906', '15248', '7267', '7266', '12239', '3897', '3850', '14969', '15711', '10162', '9745', '9741', '5953', '4074', '4071', '15819', '15813', '15787', '15789', '14497', '16623', '14249', '3341', '5926', '5920', '5921', '13477', '14563', '14562', '4666', '1288', '4662', '1471', '689', '1478', '13235', '2537', '8641', '1303', '3411', '13932', '16968', '14950', '23', '2154', '13488', '10661', '6247', '14668', '3937', '10283', '11577', '10117', '10440', '11571', '9668', '14996', '3840', '13259', '15625', '15623', '3305', '3302', '14204', '14205', '4045', '4046', '5254', '5311', '14921', '15826', '16896', '3668', '245', '15139', '13431', '7629', '6778', '2273', '12004', '12002', '16393', '16390', '3594', '9534', '456', '2452', '11500', '9873', '2889', '4822', '10621', '3643', '3544', '14027', '8371', '10014', '11137', '12576', '3443', '7761', '6211', '10496', '16504', '278', '279', '13062', '7758', '10086', '10082', '13265', '8766', '13104', '1476', '14487', '16756', '3057', '12681', '6887', '2232', '16299', '1140', '1265', '2094', '5507', '3653', '14478', '1544', '14475', '7613', '1314', '10', '2704', '15424', '125', '8024', '4876', '5613', '1821', '929', '6253', '10618', '14839', '10534', '13693', '5713', '3473', '3571', '13221', '5682', '6826', '2936', '12917', '10547', '9003', '12340', '780', '11511', '4058', '11299', '15835', '11294', '15127', '11091', '8213', '15083', '15081', '64', '65', '7094', '1427', '6249', '2403', '10319', '6270', '9502', '8067', '15013', '10989', '7846', '16421', '7843', '5651', '14549', '6357', '2890', '2893', '13793', '10002', '2739', '876', '7200', '9', '15475', '7477', '14863', '646', '10508', '8316', '8311', '4794', '4790', '14166', '11497', '3541', '5838', '9887', '9881', '13174', '9888', '16765', '12294', '14263', '1078', '1077', '16163', '13685', '12393', '1093', '6995', '16266', '16265', '15117', '14581', '11674', '11675', '6390', '4647', '2473', '1329', '15436', '16944', '5601', '6459', '12586', '15438', '11032', '6263', '2783', '13302', '5561', '2079', '9462', '12155', '6300', '2601', '5430', '1931', '9648', '7373', '7428', '6835', '2925', '12731', '3856', '13674', '15150', '757', '565', '3071', '2958', '2959', '7088', '1044', '5422', '1744', '2357', '9453', '6949', '1891', '6944', '14413', '6340', '6345', '6347', '14550', '11074', '7970', '11725', '6180', '13784', '2720', '15271', '7234', '7940', '7233', '5049', '8414', '12410', '6234', '10779']
####phospho

gene_list_10 = [str(int(float(x))) for x in gene_list_10_o]
selected_gene_set = converted_df1.loc[gene_list_10]

d1 = selected_gene_set['Res0'] - selected_gene_set['Con0']
df0 = pd.DataFrame(d1, columns=['R0/C0'], index=d1.index.tolist())

df0 = toentrez(df0)
df0 = dup_peptide2protein(df0)

#print df0

d1 = selected_gene_set['Con30'] - selected_gene_set['Con0']
cont10 = pd.DataFrame(d1, columns=['C30/C0'], index=d1.index.tolist())

cont10 = toentrez(cont10)
cont10 = dup_peptide2protein(cont10)

d1 = selected_gene_set['Res30'] - selected_gene_set['Con0']
res10 = pd.DataFrame(d1, columns=['R30/C0'], index=d1.index.tolist())

res10 = toentrez(res10)
res10 = dup_peptide2protein(res10)

cont10 = fi.pandas_if_change(cont10,oper='>',comp_numb=2,change_numb=2)
cont10 = fi.pandas_if_change(cont10,oper='<',comp_numb=-2,change_numb=-2)
res10 = fi.pandas_if_change(res10,oper='>',comp_numb=2,change_numb=2)
res10 = fi.pandas_if_change(res10,oper='<',comp_numb=-2,change_numb=-2)


#selected_gene_df0 = df0.loc[gene_list_10_o]
selected_gene_df0 = df0

df0_index = [str(int(float(x))) for x in selected_gene_df0.index.tolist()]

cont0 = fi.pandas_if_change(df0,oper='>',comp_numb=2,change_numb=2)
cont0 = fi.pandas_if_change(df0,oper='<',comp_numb=-2,change_numb=-2)

union_10min = pd.concat([cont0,cont10,res10], axis=1)

###annotation
#f = open("all_data/global_entrez_symbol0.tsv","r")
f = open("all_data/phospho_entrez_symbol0.tsv","r")
dat = f.read().strip()
dat = [x.split("\t") for x in dat.split("\r\n")]
dat = dat[1:]

conv_index = []

for x in union_10min.index.tolist():
	temp=[y[3] for y in dat if str(int(float(y[4])))==x]
	conv_index.append(temp[0])

union_10min.index = conv_index
###annotation

#union_10min.to_csv('phospho_10min_selected_candidates',sep='\t',encoding='utf-8')
#print len(union_10min)
pt.hierarchical_dendro_and_heatmap(union_10min, 'variable_30min_phospho', dend_line_width=5.0, heatmap_line_width=0.0,heatmap_font_size=3)
"""

"""
#####scatter and reg_line for correlation on phospho section
#gdata, ghead, gindex = ex.file_read(argv[1],tstatus="file", ind_loc=1, index_type='str')
gdata, ghead, gindex = ex.file_read(argv[1],tstatus="file", index_type='float')
df1 = ex.pandas_data(gdata,ghead, gindex)
selected = df1[['C10','R10']]

#converted = fi.id_conversion(selected,'all_data/global_entrez_list0.tsv', converted_index_type = 'str', converted_id_type='float')
converted = fi.id_conversion(selected,'all_data/phospho_entrez_list0.tsv', converted_index_type = 'float', converted_id_type='float')
#converted = fi.id_conversion(selected,'all_data/phospho_entrez_list0.tsv', original_id_type='float' ,converted_id_type='float')

selected_set = fi.set_intersect(converted, 'erbb_signal_kegg.edit.txt', set_id_loc=1)

#pt.scatter_plot_and_line(selected_set, x_dat='Con30', y_dat='Res30', annotate=True, annotation_file='all_data/global_entrez_symbol0.tsv', line_sep='\r\n', filename='c10_r10_scatter_g')
#pt.scatter_plot_and_line(selected_set, x_dat='Con10', y_dat='Res10', annotate=False, filename='c10_r10_scatter')
pt.scatter_plot_and_line(selected_set, x_dat='C10', y_dat='R10', filename='c10_r10_corr_erbb')
#####scatter and reg_line
"""

"""
######protein merging and scatter
def dup_peptide2protein(df1):

	new_values = []

	df_temp = df1
	df_temp_index = df_temp.index.tolist()
	dedup_list = list(set(df_temp_index))

	for y in dedup_list:
		temp = df_temp.loc[y].values.tolist()

		if type(temp[0]) is list:
			temp_fold = [x[0] for x in temp]#####fold
			temp_fdr = [x[1] for x in temp]#####fdr

			temp_fold_value = max(temp_fold)
			temp_fdr_value = temp_fdr[temp_fold.index(temp_fold_value)]

		else:
			temp_fold_value = temp[0]
			temp_fdr_value = temp[1]

		new_values.append([temp_fold_value, temp_fdr_value])

	res_df1 = pd.DataFrame(data=new_values, index=dedup_list, columns=df1.columns.tolist())
	return res_df1

#gdata, ghead, gindex = ex.file_read(argv[1],tstatus="file", ind_loc=1, index_type='str')
gdata, ghead, gindex = ex.file_read('plot_data/phospho_fold.tsv',tstatus="file", index_type='float')
df1 = ex.pandas_data(gdata,ghead, gindex)
gdata, ghead, gindex = ex.file_read('plot_data/phospho_fdr.tsv',tstatus="file", index_type='float')
df2 = ex.pandas_data(gdata,ghead, gindex)

phospho_df_arr = [pd.concat([df1[x],df2[x]], axis=1) for x in df1.columns.tolist()]

for i,item in enumerate(df1.columns.tolist()):
	phospho_df_arr[i].columns = ['pho_'+item, 'pho_'+item+"_fdr"]

phospho_df_arr = [x.dropna() for x in phospho_df_arr]
phospho_df_arr = [fi.id_conversion(x,'all_data/phospho_entrez_list0.tsv', converted_index_type = 'float', converted_id_type='float') for x in phospho_df_arr]
phospho_df_arr = [dup_peptide2protein(x) for x in phospho_df_arr]

phospho_df_arr = pd.concat(phospho_df_arr, axis=1)


gdata, ghead, gindex = ex.file_read('plot_data/global_fold.tsv',tstatus="file", index_type='str')
df3 = ex.pandas_data(gdata,ghead, gindex)
gdata, ghead, gindex = ex.file_read('plot_data/global_fdr.tsv',tstatus="file", index_type='str')
df4 = ex.pandas_data(gdata,ghead, gindex)

global_df_arr = [pd.concat([df3[x],df4[x]], axis=1) for x in df3.columns.tolist()]

for i,item in enumerate(df3.columns.tolist()):
	global_df_arr[i].columns = [item, item+"_fdr"]

global_df_arr = [x.dropna() for x in global_df_arr]
global_df_arr = [fi.id_conversion(x,'all_data/global_entrez_list0.tsv', converted_index_type = 'str', converted_id_type='float') for x in global_df_arr]
global_df_arr = [dup_peptide2protein(x) for x in global_df_arr]

global_df_arr = pd.concat(global_df_arr, axis=1)
all_dat = pd.concat([global_df_arr, phospho_df_arr], axis=1)

#for x in fold_c_head:
#	for y in phospho_fold.index.tolist():
#		print phospho_fold[x].loc[y], global_ext_fold[x].loc[y]
	#all_dat = [[phospho_fold[x].loc[y], global_ext_fold[x].loc[y]] for y in phospho_fold.index.tolist()]
	#print all_dat

temp_head = ['R0/C0', 'C10', 'C30', 'R10', 'R30']

ext_dat = []
ext_dat_index = []
#x = phospho y = global
for x in temp_head:
	temp = all_dat[['pho_'+x, 'pho_'+x+'_fdr',x, x+'_fdr']]
	temp = temp.dropna()
	[ext_dat_index.append(y) for y in temp.index.tolist()]
	temp = temp.values.tolist()

	for y in temp:
		ext_dat.append(y)


all_dat_adj = pd.DataFrame(data=ext_dat, index=ext_dat_index, columns=['x','x_fdr', 'y', 'y_fdr'])
pt.categ_scatter_plot(all_dat_adj, x_ax='x',y_ax='y')
"""


"""
#######added transcriptome
#ex.wfile_open(argv[1],output=2)

def dup_peptide2protein(df1):

	new_values = []

	df_temp = df1
	df_temp_index = df_temp.index.tolist()
	dedup_list = list(set(df_temp_index))

	for y in dedup_list:
		temp = df_temp.loc[y].values.tolist()

		if type(temp[0]) is list:
			temp_fold = [x[0] for x in temp]#####fold
			temp_fdr = [x[1] for x in temp]#####fdr

			temp_fold_value = max(temp_fold)
			temp_fdr_value = temp_fdr[temp_fold.index(temp_fold_value)]

		else:
			temp_fold_value = temp[0]
			temp_fdr_value = temp[1]

		new_values.append([temp_fold_value, temp_fdr_value])

	res_df1 = pd.DataFrame(data=new_values, index=dedup_list, columns=df1.columns.tolist())
	return res_df1

#gdata, ghead, gindex = ex.file_read(argv[1],tstatus="file", ind_loc=1, index_type='str')
gdata, ghead, gindex = ex.file_read('plot_data/phospho_fold.tsv',tstatus="file", index_type='float')
df1 = ex.pandas_data(gdata,ghead, gindex)
gdata, ghead, gindex = ex.file_read('plot_data/phospho_fdr.tsv',tstatus="file", index_type='float')
df2 = ex.pandas_data(gdata,ghead, gindex)
gdata, ghead, gindex = ex.file_read('plot_data/transcriptome_data0.tsv',tstatus="file", index_type='float')
trans_df_arr = ex.pandas_data(gdata,ghead, gindex)
trans_df_arr = dup_peptide2protein(trans_df_arr)
trans_df_arr.columns = ['trans_fc', 'trans_fc_fdr']

phospho_df_arr = [pd.concat([df1[x],df2[x]], axis=1) for x in df1.columns.tolist()]

for i,item in enumerate(df1.columns.tolist()):
	phospho_df_arr[i].columns = ['pho_'+item, 'pho_'+item+"_fdr"]

phospho_df_arr = [x.dropna() for x in phospho_df_arr]
phospho_df_arr = [fi.id_conversion(x,'all_data/phospho_entrez_list0.tsv', converted_index_type = 'float', converted_id_type='float') for x in phospho_df_arr]
phospho_df_arr = [dup_peptide2protein(x) for x in phospho_df_arr]

phospho_df_arr = pd.concat(phospho_df_arr, axis=1)


gdata, ghead, gindex = ex.file_read('plot_data/global_fold.tsv',tstatus="file", index_type='str')
df3 = ex.pandas_data(gdata,ghead, gindex)
gdata, ghead, gindex = ex.file_read('plot_data/global_fdr.tsv',tstatus="file", index_type='str')
df4 = ex.pandas_data(gdata,ghead, gindex)

global_df_arr = [pd.concat([df3[x],df4[x]], axis=1) for x in df3.columns.tolist()]

for i,item in enumerate(df3.columns.tolist()):
	global_df_arr[i].columns = [item, item+"_fdr"]

global_df_arr = [x.dropna() for x in global_df_arr]
global_df_arr = [fi.id_conversion(x,'all_data/global_entrez_list0.tsv', converted_index_type = 'str', converted_id_type='float') for x in global_df_arr]
global_df_arr = [dup_peptide2protein(x) for x in global_df_arr]

global_df_arr = pd.concat(global_df_arr, axis=1)
all_dat = pd.concat([global_df_arr, phospho_df_arr, trans_df_arr], axis=1)
network_list = [9656,4926,9320,85440,2064,4691,11346,641,5591,90381,5579,9787,996,995,3609,64423,545,983,9972,1017,6197,3020,10250,23047,5339,274,51143,7157,1063,207,84295,23253,1616,11186,2177,10606,5566,8683,4134,8030,57062,29127,5516,5515,51530,7916,9874,9877,1871,9212,5903,114885,6238,55917,55914,27316,3159,4790,8317,6850,3096,4595,4131,22848,54206,23371,6714,1956,3838,5058,5883,23476,9610,4082,25,9124,5578,27043,85456,55183,58525,7158,8493,3008,7791,5925,63967,9967,3190,472,1809,2932,8417,51593,84823,9782,3572,23215,4288,22974,1111]
network_list = [str(x) for x in network_list]
all_dat = all_dat.loc[network_list]

cont=7

if cont==0:
	selected_dat_adj = all_dat[['C10','C10_fdr','trans_fc', 'trans_fc_fdr']]
	selected_dat_adj = selected_dat_adj.dropna()
	pt.categ_scatter_plot(selected_dat_adj, x_ax='trans_fc',y_ax='C10', fn='global_c10min')
elif cont==1:
	selected_dat_adj = all_dat[['C30','C30_fdr','trans_fc', 'trans_fc_fdr']]
	selected_dat_adj = selected_dat_adj.dropna()
	pt.categ_scatter_plot(selected_dat_adj, x_ax='trans_fc',y_ax='C30', fn='global_c30min')
elif cont==2:
	selected_dat_adj = all_dat[['R10','R10_fdr','trans_fc', 'trans_fc_fdr']]
	selected_dat_adj = selected_dat_adj.dropna()
	pt.categ_scatter_plot(selected_dat_adj, x_ax='trans_fc',y_ax='R10', fn='global_r10min')
elif cont==3:
	selected_dat_adj = all_dat[['R30','R30_fdr','trans_fc', 'trans_fc_fdr']]
	selected_dat_adj = selected_dat_adj.dropna()
	pt.categ_scatter_plot(selected_dat_adj, x_ax='trans_fc',y_ax='R30', fn='global_r30min')
elif cont==4:
	selected_dat_adj = all_dat[['pho_C10','pho_C10_fdr','trans_fc', 'trans_fc_fdr']]
	selected_dat_adj = selected_dat_adj.dropna()
	pt.categ_scatter_plot(selected_dat_adj, x_ax='trans_fc',y_ax='pho_C10', fn='phospho_c10min')
elif cont==5:
	selected_dat_adj = all_dat[['pho_C30','pho_C30_fdr','trans_fc', 'trans_fc_fdr']]
	selected_dat_adj = selected_dat_adj.dropna()
	pt.categ_scatter_plot(selected_dat_adj, x_ax='trans_fc',y_ax='pho_C30', fn='phospho_c30min')
elif cont==6:
	selected_dat_adj = all_dat[['pho_R10','pho_R10_fdr','trans_fc', 'trans_fc_fdr']]
	selected_dat_adj = selected_dat_adj.dropna()
	pt.categ_scatter_plot(selected_dat_adj, x_ax='trans_fc',y_ax='pho_R10', fn='phospho_r10min')
elif cont==7:
	selected_dat_adj = all_dat[['pho_R30','pho_R30_fdr','trans_fc', 'trans_fc_fdr']]
	selected_dat_adj = selected_dat_adj.dropna()
	pt.categ_scatter_plot(selected_dat_adj, x_ax='trans_fc',y_ax='pho_R30', fn='phospho_r30min')


#for x in fold_c_head:
#	for y in phospho_fold.index.tolist():
#		print phospho_fold[x].loc[y], global_ext_fold[x].loc[y]
	#all_dat = [[phospho_fold[x].loc[y], global_ext_fold[x].loc[y]] for y in phospho_fold.index.tolist()]
	#print all_dat

#all_dat_adj = pd.DataFrame(data=ext_dat, index=ext_dat_index, columns=['x','x_fdr', 'y', 'y_fdr'])

#pt.categ_scatter_plot(all_dat_adj, x_ax='x',y_ax='y')

#######added transcriptome
"""

"""
#new added transcriptome
####global
gdata, ghead, gindex = ex.file_read("all_data/global_list0.tsv",tstatus="file", ind_loc=1)
df1 = ex.pandas_data(gdata,ghead, gindex)
converted_df1 = fi.id_conversion(df1,'all_data/global_entrez_list0.tsv', original_id_type='str', converted_index_type='str', converted_id_type='float')
converted_df1 = fi.dedup_peptide2protein(converted_df1)

####phospho
gdata, ghead, gindex = ex.file_read("all_data/phospho_list0.tsv",tstatus="file", index_type='float')
df2 = ex.pandas_data(gdata,ghead, gindex)
converted_df2 = fi.id_conversion(df2,'all_data/phospho_entrez_list0.tsv', original_id_type='float' ,converted_id_type='float')
converted_df2 = fi.dedup_peptide2protein(converted_df2)

####transcriptome
gdata, ghead, gindex = ex.file_read('plot_data/transcriptome_data0.tsv',tstatus="file", index_type='float')
trans_df_arr = ex.pandas_data(gdata,ghead, gindex)
trans_df_arr = trans_df_arr['FC'].to_frame(name='trans_fc')
converted_df3 = fi.dedup_peptide2protein(trans_df_arr)

g1 = converted_df1['Res0'] - converted_df1['Con0']
g1 = g1.to_frame(name='R0/C0')
g2 = converted_df1['Res10'] - converted_df1['Con0']
g2 = g2.to_frame(name='R10/C0')
g3 = converted_df1['Res30'] - converted_df1['Con0']
g3 = g3.to_frame(name='R30/C0')

p1 = converted_df2['Res0'] - converted_df2['Con0']
p1 = p1.to_frame(name='R0/C0')
p2 = converted_df2['Res10'] - converted_df2['Con0']
p2 = p2.to_frame(name='R10/C0')
p3 = converted_df2['Res30'] - converted_df2['Con0']
p3 = p3.to_frame(name='R30/C0')

g_fold = pd.concat([g1,g2,g3,converted_df3], axis=1)
p_fold = pd.concat([p1,p2,p3,converted_df3], axis=1)

network_list = [9656,4926,9320,85440,2064,4691,11346,641,5591,90381,5579,9787,996,995,3609,64423,545,983,9972,1017,6197,3020,10250,23047,5339,274,51143,7157,1063,207,84295,23253,1616,11186,2177,10606,5566,8683,4134,8030,57062,29127,5516,5515,51530,7916,9874,9877,1871,9212,5903,114885,6238,55917,55914,27316,3159,4790,8317,6850,3096,4595,4131,22848,54206,23371,6714,1956,3838,5058,5883,23476,9610,4082,25,9124,5578,27043,85456,55183,58525,7158,8493,3008,7791,5925,63967,9967,3190,472,1809,2932,8417,51593,84823,9782,3572,23215,4288,22974,1111]
network_list = [str(x) for x in network_list]

g_fold = g_fold.loc[network_list]
p_fold = p_fold.loc[network_list]

#pt.scatter_plot_and_line(g_fold,'trans_fc','R0/C0',filename='g_r0_trans_scatter')
#pt.scatter_plot_and_line(g_fold,'trans_fc','R10/C0',filename='g_r10_trans_scatter')
#pt.scatter_plot_and_line(g_fold,'trans_fc','R30/C0',filename='g_r30_trans_scatter')
#pt.scatter_plot_and_line(p_fold,'trans_fc','R0/C0',filename='p_r0_trans_scatter')
#pt.scatter_plot_and_line(p_fold,'trans_fc','R10/C0',filename='p_r10_trans_scatter')
pt.scatter_plot_and_line(p_fold,'trans_fc','R30/C0',filename='p_r30_trans_scatter')
"""


##### correlation search
####global
gene_list_30_o = ['85349.0', '213.0', '7296.0', '833.0', '6407.0', '754.0', '735.0', '3491.0', '3881.0', '2934.0', '30011.0', '2222.0', '5837.0', '64710.0', '80135.0', '24137.0', '6372.0', '3988.0', '6319.0', '6038.0', '5055.0', '114805.0', '149499.0', '5834.0', '2069.0', '7298.0', '231.0', '84267.0', '10577.0', '2919.0', '65084.0', '2495.0', '2512.0', '26166.0', '3872.0', '7018.0']
gene_list_10_o = ['3005.0', '213.0', '85349.0', '6407.0', '85236.0', '3881.0', '114805.0', '81539.0', '64710.0', '3020.0', '80135.0', '10144.0', '10397.0', '92840.0', '6038.0', '2947.0', '2934.0', '149499.0', '6273.0', '9555.0', '55336.0', '65084.0', '7114.0', '26166.0', '3872.0']
####global

global_list = list(set(gene_list_30_o+gene_list_30_o))
global_list = [str(int(float(x))) for x in global_list]

####phospho NOT!!!! ENTREZ
pgene_list_10_o = ['1140', '1303', '16765', '1265', '12393', '13736', '10216', '8303', '13404', '3411', '4787', '14074', '14263', '7303', '4149', '16843', '1078', '1077', '14173', '3126', '10581', '13488', '16163', '7758', '343', '5651', '814', '5921', '14668', '12930', '1093', '5420', '15083', '15826', '10283', '16266', '12069', '16265', '15119', '125', '11577', '10113', '8024', '10195', '6568', '13604', '6390', '1821', '9462', '10804', '10805', '13257', '1663', '377', '1087', '17025', '4876', '2473', '13693', '11803', '15625', '15623', '11804', '11134', '16705', '9491', '119', '12586', '11032', '6263', '3473', '12014', '12017', '12011', '12010', '12013', '2984', '7889', '2079', '14921', '14296', '12155', '1395', '7730', '7761', '2607', '1232', '1233', '2936', '7214', '3668', '3318', '15139', '9648', '9888', '7731', '9481', '7629', '7453', '4813', '3509', '780', '5402', '2273', '12004', '11511', '1931', '35', '4058', '16393', '16390', '13136', '6073', '12121', '12122', '11294', '13235', '5561', '2610', '2292', '9534', '11942', '456', '258', '13674', '1906', '12731', '10257', '10256', '2154', '1427', '6249', '2889', '4822', '11500', '5838', '8311', '10621', '8602', '3897', '6365', '10319', '9', '757', '9971', '14024', '8371', '8067', '8717', '15711', '3071', '5311', '9873', '10398', '16421', '4375', '10', '2094', '5401', '2925', '14914', '4071', '2890', '3443', '5422', '15819', '16944', '1749', '2339', '1744', '12340', '6253', '2355', '15933', '14871', '15150', '9745', '278', '1314', '15785', '15786', '9958', '14497', '16623', '14249', '6887', '6345', '6347', '7200', '10082', '13265', '11299', '5920', '6180', '6247', '5219', '957', '595', '16756', '9881', '6995', '2720', '12806', '7234', '14563', '7233', '5049', '17252', '3057', '11993', '14166', '11497', '3544', '12681', '13579', '3541', '1044', '12688', '12410', '13826', '11160', '4662', '11251', '9887', '1471', '1476', '430', '16299', '9560', '1478', '10779', '13174', '683', '8641']
pgene_list_30_o =['10299', '14173', '14914', '14915', '13736', '8303', '10455', '14074', '4149', '10581', '13168', '13245', '814', '2858', '12930', '12931', '14296', '13604', '5917', '1087', '595', '8602', '1122', '192', '13544', '11136', '11134', '16705', '15107', '14816', '12014', '12017', '12011', '12010', '12013', '2984', '15580', '17252', '2441', '1233', '1335', '1330', '9841', '16954', '9481', '14522', '4420', '35', '6071', '1530', '12121', '12122', '430', '2610', '3126', '3806', '11942', '3802', '1906', '15248', '7267', '7266', '12239', '3897', '3850', '14969', '15711', '10162', '9745', '9741', '5953', '4074', '4071', '15819', '15813', '15787', '15789', '14497', '16623', '14249', '3341', '5926', '5920', '5921', '13477', '14563', '14562', '4666', '1288', '4662', '1471', '689', '1478', '13235', '2537', '8641', '1303', '3411', '13932', '16968', '14950', '23', '2154', '13488', '10661', '6247', '14668', '3937', '10283', '11577', '10117', '10440', '11571', '9668', '14996', '3840', '13259', '15625', '15623', '3305', '3302', '14204', '14205', '4045', '4046', '5254', '5311', '14921', '15826', '16896', '3668', '245', '15139', '13431', '7629', '6778', '2273', '12004', '12002', '16393', '16390', '3594', '9534', '456', '2452', '11500', '9873', '2889', '4822', '10621', '3643', '3544', '14027', '8371', '10014', '11137', '12576', '3443', '7761', '6211', '10496', '16504', '278', '279', '13062', '7758', '10086', '10082', '13265', '8766', '13104', '1476', '14487', '16756', '3057', '12681', '6887', '2232', '16299', '1140', '1265', '2094', '5507', '3653', '14478', '1544', '14475', '7613', '1314', '10', '2704', '15424', '125', '8024', '4876', '5613', '1821', '929', '6253', '10618', '14839', '10534', '13693', '5713', '3473', '3571', '13221', '5682', '6826', '2936', '12917', '10547', '9003', '12340', '780', '11511', '4058', '11299', '15835', '11294', '15127', '11091', '8213', '15083', '15081', '64', '65', '7094', '1427', '6249', '2403', '10319', '6270', '9502', '8067', '15013', '10989', '7846', '16421', '7843', '5651', '14549', '6357', '2890', '2893', '13793', '10002', '2739', '876', '7200', '9', '15475', '7477', '14863', '646', '10508', '8316', '8311', '4794', '4790', '14166', '11497', '3541', '5838', '9887', '9881', '13174', '9888', '16765', '12294', '14263', '1078', '1077', '16163', '13685', '12393', '1093', '6995', '16266', '16265', '15117', '14581', '11674', '11675', '6390', '4647', '2473', '1329', '15436', '16944', '5601', '6459', '12586', '15438', '11032', '6263', '2783', '13302', '5561', '2079', '9462', '12155', '6300', '2601', '5430', '1931', '9648', '7373', '7428', '6835', '2925', '12731', '3856', '13674', '15150', '757', '565', '3071', '2958', '2959', '7088', '1044', '5422', '1744', '2357', '9453', '6949', '1891', '6944', '14413', '6340', '6345', '6347', '14550', '11074', '7970', '11725', '6180', '13784', '2720', '15271', '7234', '7940', '7233', '5049', '8414', '12410', '6234', '10779']
####phospho

f = open("plot_data/total_network_nodes_list.tsv","r")
net_list = f.read()
net_list = net_list.split("\t")

selected_net_list = [9656,4926,9320,85440,2064,4691,11346,641,5591,90381,5579,9787,996,995,3609,64423,545,983,9972,1017,6197,3020,10250,23047,5339,274,51143,7157,1063,207,84295,23253,1616,11186,2177,10606,5566,8683,4134,8030,57062,29127,5516,5515,51530,7916,9874,9877,1871,9212,5903,114885,6238,55917,55914,27316,3159,4790,8317,6850,3096,4595,4131,22848,54206,23371,6714,1956,3838,5058,5883,23476,9610,4082,25,9124,5578,27043,85456,55183,58525,7158,8493,3008,7791,5925,63967,9967,3190,472,1809,2932,8417,51593,84823,9782,3572,23215,4288,22974,1111]
selected_net_list = [str(x) for x in selected_net_list]

phospho_list = list(set(pgene_list_10_o+pgene_list_30_o))
phospho_list = pd.DataFrame(data=phospho_list, columns=['ind'], index=phospho_list)
phospho_list = fi.id_conversion(phospho_list,'all_data/phospho_entrez_list0.tsv', original_id_type='float' ,converted_id_type='float')
phospho_list = phospho_list.index.tolist()
phospho_list = list(set(phospho_list))

total_candidate = list(set(global_list+phospho_list))

####global
gdata, ghead, gindex = ex.file_read("all_data/global_list0.tsv",tstatus="file", ind_loc=1)
df1 = ex.pandas_data(gdata,ghead, gindex)
converted_df1 = fi.id_conversion(df1,'all_data/global_entrez_list0.tsv', original_id_type='str', converted_index_type='str', converted_id_type='float')
converted_df1 = fi.dedup_peptide2protein(converted_df1)

gdata, ghead, gindex = ex.file_read('plot_data/global_fold.tsv',tstatus="file", index_type='str')
bdf1 = ex.pandas_data(gdata,ghead, gindex)
converted_bdf1 = fi.id_conversion(bdf1,'all_data/global_entrez_list0.tsv', original_id_type='str', converted_index_type='str', converted_id_type='float')
converted_bdf1 = fi.dedup_peptide2protein(converted_bdf1)
converted_bdf1.columns=['original_R0/C0','C10/C0','C30/C0','R10/R0','R10/R0']

####phospho
gdata, ghead, gindex = ex.file_read("all_data/phospho_list0.tsv",tstatus="file", index_type='float')
df2 = ex.pandas_data(gdata,ghead, gindex)
converted_df2 = df2

gdata, ghead, gindex = ex.file_read('plot_data/phospho_fold.tsv',tstatus="file", index_type='float')
bdf2 = ex.pandas_data(gdata,ghead, gindex)
converted_bdf2 = fi.id_conversion(bdf2,'all_data/phospho_entrez_list0.tsv', original_id_type='float' ,converted_id_type='float')
converted_bdf2 = fi.dedup_peptide2protein(converted_bdf2)
converted_bdf2.columns=['original_R0/C0','C10/C0','C30/C0','R10/R0','R10/R0']

####transcriptome
gdata, ghead, gindex = ex.file_read('plot_data/transcriptome_data0.tsv',tstatus="file", index_type='float')
trans_df_arr = ex.pandas_data(gdata,ghead, gindex)
trans_df_arr = trans_df_arr['FC'].to_frame(name='trans_fc')
converted_df3 = fi.dedup_peptide2protein(trans_df_arr)

g1 = converted_df1['Res0'] - converted_df1['Con0']
g1 = g1.to_frame(name='R0/C0')
g2 = converted_df1['Res10'] - converted_df1['Con0']
g2 = g2.to_frame(name='R10/C0')
g3 = converted_df1['Res30'] - converted_df1['Con0']
g3 = g3.to_frame(name='R30/C0')

g_fold = pd.concat([g1,g2,g3], axis=1)

p1 = converted_df2['Res0'] - converted_df2['Con0']
p1 = p1.to_frame(name='R0/C0')
p2 = converted_df2['Res10'] - converted_df2['Con0']
p2 = p2.to_frame(name='R10/C0')
p3 = converted_df2['Res30'] - converted_df2['Con0']
p3 = p3.to_frame(name='R30/C0')

p_fold = pd.concat([p1,p2,p3], axis=1)
p_fold = fi.id_conversion(p_fold,'all_data/phospho_entrez_list0.tsv', original_id_type='float' ,converted_id_type='float')
p_fold = fi.dedup_peptide2protein(p_fold)

g_fold = pd.concat([g_fold,converted_bdf1,converted_df3], axis=1)
p_fold = pd.concat([p_fold,converted_bdf2,converted_df3], axis=1)

g_fold = g_fold.loc[selected_net_list]
p_fold = p_fold.loc[selected_net_list]

#g_fold.to_csv('global_all_network_fold_change.tsv',sep='\t',encoding='utf-8')
#p_fold.to_csv('phospho_all_network_fold_change.tsv',sep='\t',encoding='utf-8')

#g_fold.to_csv('global_selected_network_fold_change.tsv',sep='\t',encoding='utf-8')
#p_fold.to_csv('phospho_selected_network_fold_change.tsv',sep='\t',encoding='utf-8')

#g_fold.to_csv('global_selected_candidate_fold_change.tsv',sep='\t',encoding='utf-8')
#p_fold.to_csv('phospho_selected_candidate_fold_change.tsv',sep='\t',encoding='utf-8')


pt.scatter_plot_and_line(g_fold,'trans_fc','R0/C0',filename='snet_g_r0_trans_scatter')
#pt.scatter_plot_and_line(g_fold,'trans_fc','R10/C0',filename='can_g_r10_trans_scatter')
#pt.scatter_plot_and_line(g_fold,'trans_fc','R30/C0',filename='can_g_r30_trans_scatter')
#pt.scatter_plot_and_line(p_fold,'trans_fc','R0/C0',filename='snet_p_r0_trans_scatter')
#pt.scatter_plot_and_line(p_fold,'trans_fc','R10/C0',filename='snet_p_r10_trans_scatter')
#pt.scatter_plot_and_line(p_fold,'trans_fc','R30/C0',filename='snet_p_r30_trans_scatter')

##### correlation search
