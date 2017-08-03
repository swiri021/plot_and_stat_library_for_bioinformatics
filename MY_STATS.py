from scipy import stats
import numpy as np
import pandas as pd
import math

class MY_STATS:

	def ttest(self, arr1, arr2, test_vector = "row", column_selection=[]):
		result =[]
		if test_vector=="row":
			t_index = list(set(arr1.index.tolist() + arr2.index.tolist()))
			for x in t_index:
				t1 = arr1[arr1.index==x]
				t2 = arr2[arr2.index==x]

				if t1.empty!=True and t2.empty!=True:
					v, pvalue = stats.ttest_ind(t1.values.tolist()[0],t2.values.tolist()[0],nan_policy='omit')
					result.append([x, v.tolist(), pvalue])

				#else:
				#	result.append([x,np.nan,np.nan])

		elif test_vector=="column":######Needed developed
			t_head = list(set(arr1.columns.tolist() + arr2.columns.tolist()))

		return result

	def medtest(self, arr1, arr2, test_vector="row", column_selection=[]):

		if test_vector=="row":
			t_index = list(set(arr1.index.tolist() + arr2.index.tolist()))
			result =[]
			for x in t_index:
				t1 = arr1[arr1.index==x]
				t2 = arr2[arr2.index==x]

				if t1.empty!=True and t2.empty!=True:
					stat, pvalue, med, tbl = stats.median_test(t1.values.tolist()[0],t2.values.tolist()[0])
					result.append([x, stat, pvalue])
				#else:
				#	result.append([x,np.nan,np.nan])

		elif test_vector=="column":######Needed developed
			t_head = list(set(arr1.columns.tolist() + arr2.columns.tolist()))

		return result

	def ranksumtest(self, arr1, arr2, test_vector="row", column_selection=[]):

		if test_vector=="row":
			t_index = list(set(arr1.index.tolist() + arr2.index.tolist()))
			result =[]
			for x in t_index:
				t1 = arr1[arr1.index==x]
				t2 = arr2[arr2.index==x]

				if t1.empty!=True and t2.empty!=True:
					stat, pvalue = stats.ranksums(t1.values.tolist()[0],t2.values.tolist()[0])
					result.append([x, stat, pvalue])
				#else:
				#	result.append([x,np.nan,np.nan])

		elif test_vector=="column":######Needed developed
			t_head = list(set(arr1.columns.tolist() + arr2.columns.tolist()))

		return result

	def gene_set_zscore(self, arr1, gene_set='', gene_set_index_loc=1 ,sample_status="single"):

		def cal_z(temp_df, inter):

			selected = temp_df.loc[inter].tolist()
			selected_adj = [float(x) for x in selected]
			total = temp_df.tolist()
			total_adj = [float(x) for x in total]

			diff_mean = np.mean(selected_adj) - np.mean(total_adj)
			result = diff_mean*np.sqrt(len(selected_adj))/np.std(total_adj)

			return result

		ft = open(gene_set, "r")
		ft_dat = ft.read()
		ft_dat = ft_dat.split("\n")

		a_strip = lambda x : [a.strip() for a in x]
		ft_dat = [a_strip(x.split("\t")) for x in ft_dat]

		ft_dat_index = [x[gene_set_index_loc] for x in ft_dat]
		ft_dat = [x[gene_set_index_loc:] for x in ft_dat]

		zscore=[]

		arr1_index = arr1.index.tolist()
		inter = list(set(arr1_index).intersection(ft_dat_index))

		if sample_status=="single":
			zscore = [cal_z(arr1, inter)]

		elif sample_status=="multiple":
			zscore = [cal_z(arr1[x], inter) for x in arr1.columns.tolist()]

		return zscore, inter

	###arr1 = ['index','p-value']
	def pval2zscore(self, arr1):
		pindex = [x[0] for x in arr1]
		zscores = [stats.norm.ppf(x[1]) for x in arr1]

		return pindex, [ zscores[a] for a in range(len(pindex))]

	def zscore2pval(self, arr1, side=1):
		pindex = [x[0] for x in arr1]

		if side==1:
			pval = [stats.norm.sf(abs(x[1])) for x in arr1]
		elif side==2:
			pval = [stats.norm.sf(abs(x[1]))*2 for x in arr1]

		return pindex, [ pval[a] for a in range(len(pindex))]