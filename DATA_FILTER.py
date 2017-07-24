import numpy as np
import pandas as pd
import itertools
from EXC_READ import EXC_READ
import math

class DATA_FILTER:

	#### nan_filter(data, index, filter number)
	#### filter number = 2 : less than 2 NaN passed
	def nan_inf_filter(self, arr, r_index, lim=2, val="NaN"):

		if val=="NaN":
			nan = lambda x : [ math.isnan(float(a)) for a in x ]
		elif val=="Inf":
			nan = lambda x : [ math.isinf(float(a)) for a in x ]

		arr_filtered = [ i for i, item in enumerate(arr) if nan(item).count(True) < lim]
		arr_filtered_index = [r_index[x] for x in arr_filtered]
		arr_filtered_data = [arr[x] for x in arr_filtered]

		return arr_filtered_data, arr_filtered_index

	def median_centering(self, arr):
		med_arr = [np.median(x) for x in arr]
		centered = [ 0 if math.isnan(arr[a]) else arr[a]-med_arr[a] for a in range(len(arr))]
		return centered

	def log2_arr(self, arr, type=2):
		#####Log2 scaling
		total_log_dat = []
		if type==2:
			total_log_dat = [np.log2(x).tolist() for x in arr]
		elif type==10:
			total_log_dat = [np.log10(x).tolist() for x in arr]

		return total_log_dat

	def quantileNormalize(self, df_input):#### ###### it should be fixed
		df = df_input.copy()
		#compute rank
		dic = {}

		for col in df:
			dic.update({col : sorted(df[col])}) #### sort column

		sorted_df = pd.DataFrame(dic)
		rank = sorted_df.mean(axis = 1).tolist()

		for col in df:
			t = np.searchsorted(np.sort(df[col]), df[col])
			df[col] = [rank[i] for i in t]

		return df

	#####conv_file : conversion Id list file, original_id_loc : location of ID you want to convert in conversion ID list, orginal_id_type : your file ID type 'str' or 'float'
	#####converted_id_loc : location of ID in conversion ID list, converted_id_type : type of conversion ID in coversion ID list, converted_index_type : Index of conversion ID list
	##### You have to match between your original index and converted_index and result is converted ID
	#####ex) converted = fi.id_conversion(selected,'all_data/global_entrez_list0.tsv', converted_index_type = 'str', converted_id_type='float')
	##### converted_index_type : string type, converted_id_type : float, conversion list file should be look like testID\tnewID\n and your data look like testID\tdata1\tdata2\n
	def id_conversion(self, df, conv_file, original_id_loc=0, original_id_type='str',converted_id_loc=1, converted_id_type='str', converted_index_type='float'):
		ex = EXC_READ()
		#####tdata : converted id, tindex : original id
		tdata, thead, tindex = ex.file_read(conv_file,tstatus="file", ind_loc=original_id_loc, index_type=converted_index_type, return_type=2)

		if converted_id_type == 'float':
			tdata = [[str(int(float(x[converted_id_loc-1])))] for x in tdata]
		elif converted_id_type =='str':
			tdata = [[x[converted_id_loc-1]] for x in tdata]

		if original_id_type == 'str':
			df_index = df.index.tolist()
		elif original_id_type == 'float':
			df_index = [str(int(float(x))) for x in df.index.tolist()]

		list_df_id = []
		for x in df_index:
			if x in tindex:
				list_df_id.append([str(int(tdata[i][0])) for i, item in enumerate(tindex) if item==x])
			else:
				list_df_id.append([np.nan])

		list_df_entrez = [x[0] for x in list_df_id]
		result_df = pd.DataFrame(df.values.tolist(), columns=df.columns.tolist(), index=list_df_entrez)
		return result_df

	def set_intersect(self, df, set_file, set_id_loc=1, na_policy='drop', na_fill=0):

		sample_index = df.index.tolist()

		def fileread(fn):
			f = open(fn,"r")
			dat = f.read()
			el = dat.split("\n")
			el_fin = [x.split("\t") for x in el]
			return el_fin
		set_data = fileread(set_file)
		set_ids = [x[set_id_loc] for x in set_data]

		####index select######
		selected_index = [item for item in sample_index if item in set_ids]
		selected_df = df.loc[selected_index]

		if na_policy=='drop':
			selected_df = selected_df.dropna()
		elif na_policy=='fill':
			selected_df = selected_df.fillna(na_fill)

		return selected_df

	def pandas_if_change(self, df, comp_numb, change_numb, oper = '', return_type='change'):
		if return_type=='change':
			for col in df.columns:
				if oper=="<":
					df.loc[df[col] < comp_numb,col] = change_numb
				elif oper==">":
					df.loc[df[col] > comp_numb,col] = change_numb
				elif oper=="<=":
					df.loc[df[col] <= comp_numb,col] = change_numb
				elif oper==">=":
					df.loc[df[col] >= comp_numb,col] = change_numb
				elif oper=="==":
					df.loc[df[col] == comp_numb,col] = change_numb
			return df

		elif return_type=='select':
			df_s = []
			for col in df.columns:
				if oper=="<":
					df_s.append(df.loc[df[col] < comp_numb,col])
				elif oper==">":
					df_s.append(df.loc[df[col] > comp_numb,col])
				elif oper=="<=":
					df_s.append(df.loc[df[col] <= comp_numb,col])
				elif oper==">=":
					df_s.append(df.loc[df[col] >= comp_numb,col])
				elif oper=="==":
					df_s.append(df.loc[df[col] == comp_numb,col])
			return df_s

	def pandas_cust_div(self,df, type=1):
		if type==1:
			df_s1 = [df.loc[df[col]==0] for col in df.columns]
			df_s1_index = [x.index.tolist() for x in df_s1]

			df_s1_start = []
			for x in df_s1_index:
				df_s1_start = df_s1_start+x
			df_s1_index = list(set(df_s1_start))
			df_index = df.index.tolist()

			df_not_inter = [item for item in df_index if item not in df_s1_index]
			df_s1.append(df.loc[df_not_inter])
			df_s1.reverse()

			return df_s1


	def __init__(self):
		print "FILTER"
