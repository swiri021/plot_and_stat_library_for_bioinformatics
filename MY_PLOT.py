import pandas as pd
import math
import matplotlib as mlab
mlab.use('Agg')
from EXC_READ import EXC_READ
from DATA_FILTER import DATA_FILTER
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.cluster.hierarchy import dendrogram, linkage
import numpy as np
from matplotlib_venn import venn3


class MY_PLOT:
	def volcano_plot(self, df_arr, df_arr_c, x_ax=[], y_ax=[], colo=['red', 'blue'], lab=["over cut-off"], fn='volcano_plot'):
		ex = EXC_READ()

		def making_dots(arr1, arr2):

			arr1_dat = arr1.values.tolist()
			arr2_dat = arr2.values.tolist()

			arr1_head = arr1.columns.tolist()
			arr2_head = arr2.columns.tolist()

			arr1_index = arr1.index.tolist()
			arr2_index = arr2.index.tolist()

			inter_index = list(set(arr1_index).intersection(arr2_index))

			arr1_index_index = [i for i, item in enumerate(arr1_index) if item in inter_index]
			arr2_index_index = [i for i, item in enumerate(arr2_index) if item in inter_index]

			arr1_index_new = [arr1_index[x] for x in arr1_index_index]
			arr2_index_new = [arr2_index[x] for x in arr2_index_index]

			arr1_dat_new = [arr1_dat[x] for x in arr1_index_index]
			arr2_dat_new = [arr2_dat[x] for x in arr2_index_index]

			arr1_ext = ex.pandas_data(arr1_dat_new, arr1_head, arr1_index_new)
			arr2_ext = ex.pandas_data(arr2_dat_new, arr2_head, arr2_index_new)
			t_df = ex.pandas_merge([arr1_ext,arr2_ext],[arr1_head, arr2_head], [arr1_index_new, arr2_index_new])

			return t_df

		df_arr_dat = df_arr.values.tolist()
		df_arr_head = df_arr.columns.tolist()
		df_arr_index = df_arr.index.tolist()

		logp = lambda x : [ -math.log10(a) for a in x]
		df_arr_dat = [logp(x)for x in df_arr_dat]
		df1 = ex.pandas_data(df_arr_dat, df_arr_head, df_arr_index)

		df1_values = [df1[df1[x]>-math.log10(0.05)][x] for x in df_arr_head]
		df1_ex = pd.DataFrame(df1_values)
		df1_ex = df1_ex.T

		df_arr_c_dat = df_arr_c.values.tolist()
		df_arr_c_head = df_arr_c.columns.tolist()
		df_arr_c_index = df_arr_c.index.tolist()

		fold_head = ["fold_"+x for x in df_arr_c_head]

		df2 = ex.pandas_data(df_arr_c_dat, fold_head, df_arr_c_index)


		###### fold > 0
		df2_values = [df2[df2[x]>np.log2(1.5)][x] for x in fold_head]
		df2_ex = pd.DataFrame(df2_values)
		df2_ex = df2_ex.T

		t_ext1 = making_dots(df1_ex, df2_ex)

		df2 = ex.pandas_data(df_arr_c_dat, fold_head, df_arr_c_index)

		###### fold < 0
		df2_values = [df2[df2[x]<-np.log2(1.5)][x] for x in fold_head]
		df2_ex = pd.DataFrame(df2_values)
		df2_ex = df2_ex.T

		t_ext2 = making_dots(df1_ex, df2_ex)

		t_df = ex.pandas_merge([df1,df2],[df_arr_head, fold_head], [df_arr_index, df_arr_c_index])

		#t_ext.plot(kind='scatter', x="fold_"+x_ax[0], y=y_ax[0], color=colo[0], label=lab[0], ax=ax)

		if len(lab)==1:
			ax = t_df.plot(kind='scatter', x="fold_"+x_ax[0], y=y_ax[0], color='grey')

			for a in range(len(x_ax)):
				if a==0:
					#ax = t_df.plot(kind='scatter', x="fold_"+x_ax[0], y=y_ax[0], color='grey')
					t_ext1.plot(kind='scatter', x="fold_"+x_ax[a], y=y_ax[a], color=colo[a], ax=ax)
					t_ext2.plot(kind='scatter', x="fold_"+x_ax[a], y=y_ax[a], color=colo[a+1], ax=ax)
				else:
					t_df.plot(kind='scatter', x="fold_"+x_ax[a], y=y_ax[a], color='grey', ax=ax)
					t_ext1.plot(kind='scatter', x="fold_"+x_ax[a], y=y_ax[a], color=colo[a], ax=ax)
					t_ext2.plot(kind='scatter', x="fold_"+x_ax[a], y=y_ax[a], color=colo[a+1], ax=ax)

			ax.set_xlabel("Fold Change(log2)")
			ax.set_ylabel("-log(P)")

		else:
			ax = t_df.plot(kind='scatter', x="fold_"+x_ax[0], y=y_ax[0], color='grey')

			for a in range(len(x_ax)):
				if a==0:
					#ax = t_df.plot(kind='scatter', x="fold_"+x_ax[0], y=y_ax[0], color='grey')
					t_ext1.plot(kind='scatter', x="fold_"+x_ax[a], y=y_ax[a], color=colo[a], label=lab[a], ax=ax)
					t_ext2.plot(kind='scatter', x="fold_"+x_ax[a], y=y_ax[a], color=colo[a+1], label=lab[a+1], ax=ax)
				else:
					t_df.plot(kind='scatter', x="fold_"+x_ax[a], y=y_ax[a], color='grey', ax=ax)
					t_ext1.plot(kind='scatter', x="fold_"+x_ax[a], y=y_ax[a], color=colo[a], label=lab[a], ax=ax)
					t_ext2.plot(kind='scatter', x="fold_"+x_ax[a], y=y_ax[a], color=colo[a+1], label=lab[a+1], ax=ax)

			ax.set_xlabel("Fold Change(log2)")
			ax.set_ylabel("-log(P)")

		fig = ax.get_figure()
		fig.savefig(fn+".png")

	def scatter_matr(self, df):
		def corrfunc(x,y, **kws):
			r, _ = stats.spearmanr(x, y)
			ax = plt.gca()
			ax.annotate("{:.2f}".format(r), xy=(.1, .9), xycoords=ax.transAxes)

		#titles = ['androgen_response', 'androgen_biosynthesis', 'cholesterol_homeostasis', 'cholesterol_biosynthesis','steroid_biosynthesis']
		#df = pd.DataFrame(arr1, columns=titles)

		axes = sns.PairGrid(df.dropna(), palette="Blues_r")
		axes.map_upper(sns.regplot,line_kws={'color': 'red'}, scatter_kws={'s':10})
		axes.map_upper(corrfunc)
		axes.map_diag(sns.distplot, kde="True")
		axes.map_lower(sns.regplot,scatter_kws={'s':10},fit_reg=False)
		axes.savefig("scatter_matrix.png",bbox_inches='tight')

	def violin_plt(self, df):

		fig, ax = plt.subplots()
		fig.set_size_inches(2,3)
		ax = sns.violinplot(data=df,inner=None, color=".8")
		sns.swarmplot(data=df, palette="hls", ax=ax)

		fn = ax.get_figure()
		fn.savefig('violin_plot.png')

	def venn_dia(self,df, scaling='bypass'):

		def only_agroup(a,b,c):
			a1 = [x for x in a if x not in b] #A
			a_result = [x for x in a1 if x not in c] #B
			return a_result

		ind_list = [df[x].dropna().index.tolist() for x in df.columns.tolist()]

		ghead_c = df.columns.tolist()
		ind_list_c = ind_list

		cs110_1 = list(set(ind_list_c[0]).intersection(ind_list_c[1])) ###AB
		cs011_1 = list(set(ind_list_c[1]).intersection(ind_list_c[2])) ####BC
		cs101_1 = list(set(ind_list_c[0]).intersection(ind_list_c[2])) ###AC
		cs111 = list(set(cs110_1).intersection(ind_list_c[2]))####ABC

		cs110 = len([x for x in cs110_1 if x not in cs111])###AB-ABC
		cs011 = len([x for x in cs011_1 if x not in cs111]) ####BC-ABC
		cs101 = len([x for x in cs101_1 if x not in cs111]) ###AC-ABC
		cs100 = len(only_agroup(ind_list_c[0],ind_list_c[1],ind_list_c[2])) ### A
		cs010 = len(only_agroup(ind_list_c[1],ind_list_c[2],ind_list_c[0])) ### B
		cs001 = len(only_agroup(ind_list_c[2],ind_list_c[0],ind_list_c[1])) ### C
		cs111 = len(list(set(cs110_1).intersection(ind_list_c[2])))####ABC

		figure, axes = plt.subplots(1, 1)### plot array 1,1, (2,2) means 2by2

		if scaling=="bypass":
			v = venn3(subsets=(cs100, cs010, cs110, cs001, cs101, cs011, cs111), set_labels = ghead_c, ax=axes)
		elif scaling=="log2":
			v = venn3(subsets=(np.log2(cs100), np.log2(cs010), np.log2(cs110), np.log2(cs001), np.log2(cs101), np.log2(cs011), np.log2(cs111)), set_labels = ghead_c, ax=axes)

			v.get_label_by_id('100').set_text(str(cs100))
			v.get_label_by_id('010').set_text(str(cs010))
			v.get_label_by_id('001').set_text(str(cs001))
			v.get_label_by_id('110').set_text(str(cs110))
			v.get_label_by_id('011').set_text(str(cs011))
			v.get_label_by_id('101').set_text(str(cs101))
			v.get_label_by_id('111').set_text(str(cs111))

		figure.savefig("venn_diagram.png")

	def histogram_group(self, df, hist=False, rug=False , colo=['Blue'], cut_status = False, cut=0, legend=False):

		if cut_status==False:
			df_head = df.columns.tolist()
			for a in range(len(df_head)):
				axes = sns.distplot(df[df_head[a]], hist=hist, rug=rug, color=colo[a], label=df_head[a])
				if legend==False:
					axes.legend_.remove()

			#axes.set_xlabel("Fold Change(log2)")
			#axes.set_ylabel("-log(P)")

			fn = axes.get_figure()
			fn.savefig("histogram.png",bbox_inches='tight')



		else:
			df_head = df.columns.tolist()
			for a in range(len(df_head)):
				#axes = sns.distplot(df[df_head[a]], hist=hist, rug=rug, color=colo[a], label=df_head[a])
				axes = sns.kdeplot(df[df_head[a]], cut=cut, color=colo[a], label=df_head[a])
				if legend==False:
					axes.legend_.remove()

			#axes.set_xlabel("Fold Change(log2)")
			#axes.set_ylabel("-log(P)")

			fn = axes.get_figure()
			fn.savefig("histogram.png",bbox_inches='tight')

	def cluter_heatmap(self, df,method='complete', metric='euclidean',col_cluster=False):

		#cm = sns.clustermap(df, cmap="RdBu", method=method, metric=metric, figsize=(50,20), col_cluster=col_cluster)
		cm = sns.clustermap(df, cmap="RdBu_r", method=method, metric=metric, figsize=(10,8), col_cluster=col_cluster)

		hm = cm.ax_heatmap.get_position()
		plt.setp(cm.ax_heatmap.yaxis.get_majorticklabels(), fontsize=3)
		cm.ax_heatmap.set_position([hm.x0, hm.y0, hm.width*0.25, hm.height])
		col = cm.ax_col_dendrogram.get_position()
		cm.ax_col_dendrogram.set_position([col.x0, col.y0, col.width*0.25, col.height*0.5])


		cm.savefig('c_heatmap.png')

	def rankplot(self, df, picked='', selected_on = False, selected = '', entrez_col_indata=1, legend=False):
		if selected_on==False:

			df = df[picked]
			#df = abs(df)
			df_sorted = df.sort_values(df.columns.tolist(), ascending=[True])
			df_sorted['index_value'] = pd.Series(df_sorted.index.tolist(), index=df_sorted.index)

			df_sorted = df_sorted.dropna()
			ax = df_sorted.plot(x='index_value', y=picked, color='Blue')
			#ax.set_xlabel("Fold Change(log2)")
			#ax.set_ylabel("-log(P)")

			fig = ax.get_figure()
			fig.savefig("rankplot.png")

		else:

			def fileread(fn):
				f = open(fn,"r")
				dat = f.read()
				el = dat.split("\n")
				el_fin = [x.split("\t") for x in el]
				return el_fin

			df = df[picked]
			df_sorted = df.sort_values(df.columns.tolist(), ascending=[True])
			df_sorted['index_value'] = pd.Series(df_sorted.index.tolist(), index=df_sorted.index)
			df_sorted = df_sorted.dropna()
			list_df_sorted_index = df_sorted.index.tolist()

			selected_data = fileread(selected)
			selected_data_entrez = [x[entrez_col_indata] for x in selected_data]
			####index select######

			dat = [i for i,item in enumerate(list_df_sorted_index) if item in selected_data_entrez]

			ax = df_sorted.plot(x='index_value', y=picked, color='Blue', legend=legend)
			for x in dat:
				ax.axvspan(x, x+1, ymax=0.05, color='red', alpha=1.0)

			fn = ax.get_figure()
			fn.savefig('rankplot.png')

	##### change to **kwarg
	def scatter_plot_and_line(self,data, x_dat,y_dat, filename='scatter_plot', annotate = False, annotation_file='', conv_id_loc=0, conv_dat_loc=1, line_sep='\n'):

		def except_proc(arr):
			try:
				arr = [str(int(float(x))) for x in arr]
			except ValueError:
				arr = [str(x) for x in arr]
			return arr

		ax = sns.regplot(x=x_dat,y=y_dat,line_kws={'color':'red', 'lw':1},scatter_kws={'s':20}, data=data)

		if annotate==True:
			f = open(annotation_file,"r")
			a_dat = f.read().strip()
			a_dat = [x.split("\t") for x in a_dat.split(line_sep)[1:]]

			conv_ids = [x[conv_id_loc] for x in a_dat]
			conv_dats = [x[conv_dat_loc] for x in a_dat]

			conv_ids = except_proc(conv_ids)
			conv_dats = except_proc(conv_dats)

			dat_index = data.index.tolist()
			dat_conv_index = []
			for x in dat_index:
				temp = [conv_dats[i] for i,item in enumerate(conv_ids) if item==x]
				dat_conv_index.append(temp[0])

			new_data = pd.DataFrame(data.values.tolist(), dat_conv_index, data.columns.tolist())
			new_data = new_data.sort_values([y_dat,x_dat], ascending=[False,True])

			with pd.option_context('display.max_rows', None, 'display.max_columns', 3):
				print(new_data)

			#for a in range(5):
			for a in range(len(new_data.index.tolist())):
				if abs(new_data[x_dat][a]/new_data[y_dat][a]) > 1.05:
					ax.annotate(new_data.index[a], (new_data[x_dat][a], new_data[y_dat][a]),xytext=(15, 50), textcoords='offset points', arrowprops=dict(arrowstyle='-|>', connectionstyle='angle3, angleA=0, angleB=90'))

		fn = ax.get_figure()
		fn.savefig('%s.png'%(filename))

	def hierarchical_dendro_and_heatmap(self, df, filename, method='complete', metric='euclidean', cbar=False):

		ex = EXC_READ()
		Z = linkage(df.values.tolist(), method=method, metric=metric)

		#plt.title('Hierarchical Clustering Dendrogram')
		#plt.xlabel('sample index')
		#plt.ylabel('distance')

		plt.figure(figsize=(25, 10))
		A = dendrogram( Z, leaf_rotation=90.,  # rotates the x axis labels
						leaf_font_size=8.,  # font size for the x axis labels

						labels=df.index.tolist()
		)

		plt.savefig(filename+'_dendrogram.png')
		plt.clf()

		plt.figure(figsize=(4,10))
		clustered_sort = A['ivl']
		clustered_dat = []
		[clustered_dat.append(df.loc[x].tolist()) for x in clustered_sort]
		new_df = ex.pandas_data(clustered_dat, df.columns.tolist(), clustered_sort)

		sns.set(font_scale=0.6)
		ax = sns.heatmap(new_df, cbar=cbar)
		fn = ax.get_figure()
		#print new_df
		fn.savefig(filename+'_heatmap.png')


