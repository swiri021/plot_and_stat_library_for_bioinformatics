import pandas as pd
import math
import matplotlib as mlab
mlab.use('Agg')
from EXC_READ import EXC_READ
from DATA_FILTER import DATA_FILTER
import matplotlib.pyplot as plt
import seaborn as sns; sns.set_style("whitegrid")
from scipy import stats
from scipy.cluster.hierarchy import dendrogram, linkage
import numpy as np
from matplotlib_venn import venn3


class MY_PLOT:
	def volcano_plot(self, df, x_ax='', y_ax='', fn='volcano_plot'):
		ex = EXC_READ()

		df_arr_dat = df.values.tolist()
		df_arr_index = df.index.tolist()
		df_arr_dat = [[x[0],-math.log10(x[1])] for x in df_arr_dat]

		df1 = ex.pandas_data(df_arr_dat, ['x', 'y'], df_arr_index)
		df1_values = df1[df1[y_ax]>-math.log10(0.05)]

		###### fold > 0
		df2_values = df1_values[df1_values[x_ax]>np.log2(1.5)]
		df2_ex = pd.DataFrame(df2_values)

		###### fold < 0
		df2_values = df1_values[df1_values[x_ax]<-np.log2(1.5)]
		df3_ex = pd.DataFrame(df2_values)

		ax = df1.plot(kind='scatter', x='x', y='y', color='grey')
		df2_ex.plot(kind='scatter', x='x', y='y', color='red', ax=ax)
		df3_ex.plot(kind='scatter', x='x', y='y', color='blue', ax=ax)

		ax.set_xlabel("Fold Change(log2)")
		ax.set_ylabel("-log(P)")

		fig = ax.get_figure()
		fig.savefig(fn+".png")

	def scatter_matr(self, df, type=1):

		if type==1:
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
		#elif type==2:

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
		print [x for x in cs110_1 if x not in cs111]
		print cs110
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

		def corrfunc(x,y, **kws):
			r, pv = stats.spearmanr(x, y)
			print pv
			ax = plt.gca()
			ax.annotate("{:.2f}".format(r), xy=(.1, .9), xycoords=ax.transAxes)

		test_data = data[[x_dat,y_dat]].dropna()


		plt.axhline(0, color='black', linewidth=1)
		plt.axvline(0, color='black', linewidth=1)

		ax = sns.regplot(x=x_dat,y=y_dat,line_kws={'color':'red', 'lw':1},scatter_kws={'s':20}, data=test_data)

		#corrfunc(test_data[x_dat], test_data[y_dat])
		r, pv = stats.spearmanr(test_data[x_dat], test_data[y_dat], nan_policy='omit')

		print 'corr : '+str(r)
		print 'pval : '+str(pv)

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

	def hierarchical_dendro_and_heatmap(self, df, filename, method='complete', metric='euclidean', cbar=False, dend_line_width=1, heatmap_font_size = 0.6, square=False, font=False, heatmap_line_width=0.05):

		sns.reset_orig()
		plt.rcParams['lines.linewidth'] = dend_line_width
		ex = EXC_READ()

		Z = linkage(df.values.tolist(), method=method, metric=metric)
		plt.figure(figsize=(25, 10))
		A = dendrogram( Z, leaf_rotation=90.,  # rotates the x axis labels
						leaf_font_size=8.,  # font size for the x axis labels
						labels=df.index.tolist(),
		)

		plt.savefig(filename+'_dendrogram.png')

		plt.clf()

		plt.figure(figsize=(10,20))
		clustered_sort = A['ivl']
		print clustered_sort
		clustered_dat = []
		[clustered_dat.append(df.loc[x].tolist()) for x in clustered_sort]
		new_df = ex.pandas_data(clustered_dat, df.columns.tolist(), clustered_sort)

		if font==True:
			sns.set(font_scale=heatmap_font_size, font='Arial')
		ax = sns.heatmap(new_df, cbar=cbar, linewidth=heatmap_line_width, square=square)
		fn = ax.get_figure()
		#print new_df
		fn.savefig(filename+'_heatmap.png')

	def categ_scatter_plot(self, df, x_ax='', y_ax='', fn='categ_scatter_plot'):

		ex = EXC_READ()

		df1_values = df[df[x_ax+"_fdr"]<0.05]
		df2_values = df[df[y_ax+"_fdr"]<0.05]
		#df3_values = df[(abs(df[x_ax])>np.log2(1.5)) & (abs(df[y_ax])>np.log2(1.5)) & df[x_ax+"_fdr"] < 0.05 & df[y_ax+"_fdr"] < 0.05]
		df3_values = df1_values[df1_values[y_ax+"_fdr"] < 0.05]

		df1_ex = pd.DataFrame(df1_values)
		df2_ex = pd.DataFrame(df2_values)
		df3_ex = pd.DataFrame(df3_values)


		plt.axhline(0, color='black', linewidth=1)
		plt.axvline(0, color='black', linewidth=1)

		ax = df.plot(kind='scatter', x=x_ax, y=y_ax, color='grey')
		if df1_ex.empty==False:
			df1_ex.plot(kind='scatter', x=x_ax, y=y_ax, color='green', ax=ax)
		if df2_ex.empty==False:
			df2_ex.plot(kind='scatter', x=x_ax, y=y_ax, color='blue', ax=ax)
		if df3_ex.empty==False:
			df3_ex.plot(kind='scatter', x=x_ax, y=y_ax, color='cyan', ax=ax)

		sns.regplot(x=x_ax,y=y_ax,data=df, line_kws={'color': 'red'},scatter_kws={"s": 0}, color='grey')#####total data reg line
		r, pv = stats.spearmanr(df[x_ax], df[y_ax], nan_policy='omit')

		print 'corr : '+str(r)
		print 'pval : '+str(pv)
		#sns.regplot(x=x_ax,y=y_ax,data=df3_ex, line_kws={'color': 'red'}, color='cyan')#####significant data reg line



		ax.set_xlabel("Fold Change(log2)")
		ax.set_ylabel("Fold Change(log2)")

		fig = ax.get_figure()
		fig.savefig(fn+".png")
