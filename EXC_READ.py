import pandas as pd
import xlrd
import csv

class EXC_READ:

	####output = 0 : all on memory, output = 1 : one sheet on memory
	####output = 2 : write all sheet, output = 3 : write one sheet
	def wfile_open(self, fn, output=0, sh_n=0):

		xfile = xlrd.open_workbook(fn)
		sheet_n = xfile.nsheets

		if output==0:
			sheet_arr = []
			for n in range(sheet_n):
				mysheet = xfile.sheet_by_index(n)
				sh = [mysheet.row_values(rownum) for rownum in xrange(mysheet.nrows)]

				uni = lambda x : [ a.encode('utf8') if type(a) is unicode else a for a in x ]
				sh = [uni(x) for x in sh]
				sheet_arr.append(sh)
			return sheet_arr

		elif output==1:
			mysheet = xfile.sheet_by_index(sh_n)
			sh = [mysheet.row_values(rownum) for rownum in xrange(mysheet.nrows)]

			uni = lambda x : [ a.encode('utf8') if type(a) is unicode else a for a in x ]
			sh = [uni(x) for x in sh]

			return sh

		elif output==2:
			for n in range(sheet_n):
				rw = open(fn.replace(".xlsx","%d.tsv"%(n)),"w")
				wr = csv.writer(rw, delimiter="\t")
				mysheet = xfile.sheet_by_index(n)

				for rownum in xrange(mysheet.nrows):
					wr.writerow(mysheet.row_values(rownum))
			return 0

		elif output==3:
			rw = open(fn.replace(".xlsx","%d.tsv"%(sh_n)),"w")
			wr = csv.writer(rw, delimiter="\t")
			mysheet = xfile.sheet_by_index(sh_n)

			for rownum in xrange(mysheet.nrows):
				wr.writerow(mysheet.row_values(rownum))

			return 0

	#### file_read(file_status, separator, index location, head location, data return type)
	#### tstatus = "read" : function wfile_open pass, tstatus = "file" : tsv file opened
	#### data return type = 0 : float, 1 : int, 2 : str
	#### return data, columns, indexes
	def file_read(self, fn, sh_n=0, tstatus = "read",f_arr='',sep="\t", ind_loc=0, index_head_status=True, head_loc=0, return_type=0, index_type='str'):

		f_head =''
		f_body = []
		temp_arr = []

		if tstatus=="file":
			f = open(fn,"r")
			f_arr = f.read().strip()
			f_arr = f_arr.replace("\r\n","\n").split("\n")

			f_head = "NULL"
			f_body = "NULL"

			if head_loc == 0:
				f_head = f_arr[head_loc].split(sep)
				if index_head_status==True:
					f_head = f_head[ind_loc+1:]###index head not included
				f_body = f_arr[head_loc+1:]

			else:
				f_head = [x.split(sep) for x in f_arr[:head_loc+1]]
				f_body = f_arr[head_loc+1:]

			temp_arr = [x.split(sep) for x in f_body]

		elif tstatus=="read":
			self.wfile_open(fn,output=1)
			f_head = f_arr[head_loc]
			if index_head_status==True:
				f_head = f_head[ind_loc+1:]###index head not included
			f_body = f_arr[head_loc+1:]
			temp_arr = f_body

		if index_type == 'float':
			r_index = [str(int(float(x[ind_loc]))) for x in temp_arr]###### first version of index !
		elif index_type== 'str':
			r_index = [x[ind_loc] for x in temp_arr]###### first version of index !

		if return_type==0:#######float
			floating = lambda x : [float(a) for a in x] ## data #array floating
			total_dat = [floating(x[ind_loc+1:]) for x in temp_arr]
		elif return_type==1:#########int
			floating = lambda x : [int(a) for a in x] ## data #array floating
			total_dat = [floating(x[ind_loc+1:]) for x in temp_arr]
		elif return_type==2:########string
			floating = lambda x : [str(a) for a in x] ## data #array floating
			total_dat = [floating(x[ind_loc+1:]) for x in temp_arr]

		return total_dat, f_head, r_index

	#### pandas_data(data, head, index, numb of group, group head set)
	#### numb of group default = 1, group set should be included in head
	#### ex) head = [el1,el2,el3,el4]; group_set = [[el1, el2],[el3,el4]]
	#### return pandas array
	def pandas_data(self, dat, head, index,  group=1, group_set = []):
		df = pd.DataFrame(dat,index=index, columns=head)
		if group==1:
			return df
		else:
			if len(group_set)!=group:
				print "WRONG GROUP SETTING"
				exit()
			else:
				for gs in range(len(group_set)):
					if len(list(set(group_set[gs]).intersection(set(head))))!=len(group_set[gs]):
						print "WRONG GROUP ELEMENT"
						exit()
			result_set = [ df[x] for x in group_set ]

			return result_set

	#### pandas merge(data array, head array, index array, status)
	#### explain? go to https://pandas.pydata.org/pandas-docs/stable/merging.html
	def pandas_merge(self, dat, head, index, status=0):
		if len(dat)!=len(head):
			print "WRONG merge status"
			exit()

		df_arr = [pd.DataFrame(dat[a], index=index[a], columns=head[a]) for a in range(len(head))]

		if status==0:
			result = pd.concat(df_arr, axis=1)
			return result

		elif status==1:
			result = pd.append(df_arr)
			return result

	def __init__(self):
		print "READY TO FILE READ"


	"""
	def __del__(self):
		print "END"
	"""