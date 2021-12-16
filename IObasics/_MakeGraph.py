#!/usr/bin/env python
# -*- coding: utf-8 -*-
################################################################################
import numpy as np
import platform
import subprocess
import os
################################################################################
class MakeHistNW:
	def __init__(self, cond_list, target_name):
		# cond_list = [base_name, data_list, bins, normalize, Legend, option]
		self.list = cond_list[1]
		self.bins = cond_list[2]
		#
		self.dir = target_name
		#
		self.base = cond_list[0]
		self.norm = cond_list[3]
		#
		self.f_dat = "nw_hist.dat"
		self.f_plt = "make_hist.plt"
		self.f_png = "histgram.png"
		self.leg = cond_list[4]
		self.option = cond_list[5]

	############################################################
	# ヒストグラムのグラフの作成
	def make_hist_all(self):
		# ヒストグラムのデータ作成
		bin_width, hist_data, val, x = self.make_hist_data()
		# ヒストグラムのデータを書き出し 
		self.write_data(hist_data, bin_width)
		# グラフを作成
		self.make_graph(bin_width)
		return x, val

	############################################################
	# ヒストグラムのデータ作成
	def make_hist_data(self):
		# ヒストグラムを作成
		weights = np.ones(len(self.list))/float(len(self.list))
		if self.norm:
			val, x = np.histogram(self.list, bins=self.bins, weights= weights)
		else:
			val, x = np.histogram(self.list, bins=self.bins)
		# グラフ用にデータを変更
		bin_width = (x[1]-x[0])
		mod_x = (x + bin_width/2)[:-1]
		hist_data = np.stack([mod_x, val], axis = 1)
		return bin_width, hist_data, val, x

	##############################
	# ヒストグラムのデータを書き出し 
	def write_data(self, hist_data, bin_width):
		os.makedirs(self.dir, exist_ok=True)
		with open(os.path.join(self.dir, self.f_dat), 'w') as f:
			f.write("# Histgram data:\n\n")
			for line in hist_data:
				f.write(str(line[0]) + '\t' + str(line[1])  + '\n')
		return

	###########################################################
	# グラフを作成
	def make_graph(self, bin_width):
		self.make_script(bin_width)
		cwd = os.getcwd()
		os.chdir(self.dir)
		if platform.system() == "Windows":
			subprocess.call(self.f_plt, shell=True)
		elif platform.system() == "Linux":
			subprocess.call('gnuplot ' + self.f_plt, shell=True)
		os.chdir(cwd)
		return
		

	#######################
	# 必要なスクリプトを作成
	def make_script(self, bin_width):
		with open(os.path.join(self.dir, self.f_plt), 'w') as f:
			script = self.script_content(bin_width)
			f.write(script)
		return

	#################
	# スクリプトの中身
	def script_content(self, bin_width):
		script = 'set term pngcairo font "Arial,14" \nset colorsequence classic \n'
		#
		script += '# \ndata = "' + self.f_dat + '" \nset output "' + self.f_png + ' "\n'
		#
		script += '#\nset size square\n# set xrange [0:1.0]\n#set yrange [0:100]\n'
		script += '#\nset xlabel "' + self.leg[0] + '"\nset ylabel "' + self.leg[1] + '"\n\n'
		#
		if self.base == "Rx":
			script += '#f = 3\n#N = 39\n#R1 = (1.75*N)**0.5\n#Pos = R1/2**0.5\n#delta = Pos*(2./3.)**0.5\n\n'
			script += 'f = 4\nN = 39\nR1 = (1.75*N)**0.5\nPos = R1/3**0.5\ndelta = Pos*(2./3.)**0.5\n\n'
			script += 'C=0.25\nf(x) = C*(1./2.)*(1./(delta*(3.142*2.)**0.5))*(exp(-1.*((x-Pos)**2)/(2.*delta**2)) + exp(-1.*((x+Pos)**2)/(2.*delta**2)))\n\n'
		#
		if self.base == "R":
			script += 'N = 39\nb=0.97\nCN=1.7*1.5\nC=0.02\n'
			script += 'f(x, CN) = C*4*pi*x**2*(3/(2*pi*N*CN*b**2))**(3/2)*exp(-3*x**2/(2*N*CN*b**2))\n'	
			script += '#fit f(x, CN) data via CN, C\n'
		#
		if self.base == "angle":
			if self.option != "box":
				script += 'plot data u 1:($2/(3.142*sin(3.142*$1/180))) w l noti'
			else:
				script += 'set style fill solid 0.5\nset boxwidth ' + str(bin_width) + '\n'
				script += 'plot data u 1:($2/(3.142*sin(3.142*$1/180))) w boxes noti'
		elif self.option == "box":
			script += 'set style fill solid 0.5\nset boxwidth ' + str(bin_width) + '\n'
			script += '#\nplot data w boxes noti'
		else:
			if self.base == "Rx":'#\nplot data w l noti'
		if self.base == "Rx":
			script += ', \\\n f(x)'
		elif self.base == "R":
			script += ', \\\n f(x, CN)'
		return script

################################################################################
class MakeHist:
	def __init__(self, cond_list, target_name):
		# cond_list = [base_name, data_list, bins, normalize, Legend, option]
		self.list = cond_list[1]
		self.bins = cond_list[2]
		if target_name != '':
			self.dir = os.path.join(target_name, cond_list[0])
		else:
			self.dir = cond_list[0]

		self.base = cond_list[0]
		self.norm = cond_list[3]
		#
		self.f_dat = cond_list[0] + "_hist.dat"
		self.f_plt = cond_list[0] + ".plt"
		self.f_png = cond_list[0] + ".png"
		self.leg = cond_list[4]
		self.option = cond_list[5]
	############################################################
	# ヒストグラムのグラフの作成
	def make_hist_all(self):
		# ヒストグラムのデータ作成
		bin_width, hist_data = self.make_hist_data()
		# ヒストグラムのデータを書き出し 
		self.write_data(hist_data, bin_width)
		# グラフを作成
		self.make_graph(bin_width)
		return

	############################################################
	# ヒストグラムのデータ作成
	def make_hist_data(self):
		# ヒストグラムを作成
		weights = np.ones(len(self.list))/float(len(self.list))
		if self.norm:
			val, x = np.histogram(self.list, bins=self.bins, weights= weights)
		else:
			val, x = np.histogram(self.list, bins=self.bins)
		# グラフ用にデータを変更
		bin_width = (x[1]-x[0])
		mod_x = (x + bin_width/2)[:-1]
		hist_data = np.stack([mod_x, val], axis = 1)
		return bin_width, hist_data

	##############################
	# ヒストグラムのデータを書き出し 
	def write_data(self, hist_data, bin_width):
		os.makedirs(self.dir, exist_ok=True)
		with open(os.path.join(self.dir, self.f_dat), 'w') as f:
			f.write("# Histgram data:\n\n")
			for line in hist_data:
				f.write(str(line[0]) + '\t' + str(line[1])  + '\n')
		return

	###########################################################
	# グラフを作成
	def make_graph(self, bin_width):
		self.make_script(bin_width)
		cwd = os.getcwd()
		os.chdir(self.dir)
		if platform.system() == "Windows":
			subprocess.call(self.f_plt, shell=True)
		elif platform.system() == "Linux":
			subprocess.call('gnuplot ' + self.f_plt, shell=True)
		os.chdir(cwd)
		return

	#######################
	# 必要なスクリプトを作成
	def make_script(self, bin_width):
		with open(os.path.join(self.dir, self.f_plt), 'w') as f:
			script = self.script_content(bin_width)
			f.write(script)
		return

	#################
	# スクリプトの中身
	def script_content(self, bin_width):
		script = 'set term pngcairo font "Arial,14" \nset colorsequence classic \n'
		#
		script += '# \ndata = "' + self.f_dat + '" \nset output "' + self.f_png + ' "\n'
		#
		script += '#\nset size square\n#set xrange [0:]\n#set yrange [0:100]\n'
		script += '#\nset xlabel "' + self.leg[0] + '"\nset ylabel "' + self.leg[1] + '"\n\n'
		#
		if self.base == "Rx":
			n_seg = self.option[0]
			bond = self.option[1]
			script += 'N = ' + str(n_seg) + '\n'
			script += 'bond = ' + str(bond) + '\n'
			script += 'CN=1.7\n'
			script += 'C=0.1\n\n'
			#
			script += 'f(x) = C*(3/(2*pi*N*CN*bond**2))**(3/2)*exp(-3*x**2/(2*N*CN*bond**2))\n\n'
			script += 'fit f(x) data via C, CN\n\n'
			script += '#\nset label 1 sprintf("C_N=%.3f", CN) at graph 0.7, 0.8\n\n'
			#
			script += 'set style fill solid 0.5\nset boxwidth ' + str(bin_width) + '\n'
			script += '#\nplot data w boxes noti'
			script += ', \\\n f(x)'

		# if self.base == "Rx":
		# 	if (type(self.option) == list) and len(self.option) == 4:
		# 		n_seg = self.option[0]
		# 		bond = self.option[1]
		# 		cn = self.option[2]
		# 		func = self.option[3]
		# 	elif (type(self.option) == list) and len(self.option) == 2:
		# 		n_seg = self.option[0]
		# 		bond = self.option[1]
		# 		cn = 1.7
		# 		func = 0
		# 	else:
		# 		n_seg = 39
		# 		bond = 0.97
		# 		cn = 1.7
		# 		func = 4
		# 	script += 'N = ' + str(n_seg) + '\n'
		# 	script += 'bond = ' + str(bond) + '\n'
		# 	script += 'CN = ' + str(cn) + '\n'
		# 	script += 'f = ' + str(func) + '\n'
		# 	#
		# 	script += 'R1 = CN*(N**0.5)*bond\n'
		# 	script += 'C=0.25\n\n'
		# 	#
		# 	if func == 3:
		# 		script += 'Pos = R1/2**0.5\ndelta = Pos*(1. - 2./f)**0.5\n\n'
		# 		script += 'f(x) = C*(1./2.)*(1./(delta*(3.142*2.)**0.5))*(exp(-1.*((x-Pos)**2)/(2.*delta**2)) + exp(-1.*((x+Pos)**2)/(2.*delta**2)))\n\n'
		# 		script += 'fit f(x) data via C\n\n'
		# 	elif func == 4:
		# 		script += 'Pos = R1/3**0.5\ndelta = Pos*(1. - 2./f)**0.5\n\n'
		# 		script += 'f(x) = C*(1./2.)*(1./(delta*(3.142*2.)**0.5))*(exp(-1.*((x-Pos)**2)/(2.*delta**2)) + exp(-1.*((x+Pos)**2)/(2.*delta**2)))\n\n'
		# 		script += 'fit f(x) data via C\n\n'
		# 	#
		# 	script += 'set style fill solid 0.5\nset boxwidth ' + str(bin_width) + '\n'
		# 	script += '#\nplot data w boxes noti'
		# 	script += ', \\\n f(x)'
		#
		if self.base == "R":
			if (type(self.option) == list) and len(self.option) == 4:
				n_seg = self.option[0]
				bond = self.option[1]
				cn = self.option[2]
				func = self.option[3]
			elif (type(self.option) == list) and len(self.option) == 2:
				n_seg = self.option[0]
				bond = self.option[1]
				cn = 1.7
				func = 0
			else:
				n_seg = 39
				bond = 0.97
				cn = 1.7
				func = 4
			script += 'N = ' + str(n_seg) + '\n'
			script += 'bond = ' + str(bond) + '\n'
			script += 'CN = ' + str(cn) + '\n'
			script += 'f = ' + str(func) + '\n'
			script += 'C = 0.02\n\n'
			script += 'f(x, CN) = C*4.*pi*x**2.*(3./(2.*pi*N*CN*bond**2.))**(3./2.)*exp(-3.*x**2./(2.*N*CN*bond**2.))\n'	
			script += 'fit f(x, CN) data via CN, C\n\n'
			script += '#\nset label 1 sprintf("C_N=%.3f", CN) at graph 0.7, 0.8\n\n'
			script += 'set style fill solid 0.5\nset boxwidth ' + str(bin_width) + '\n'
			script += '#\nplot data w boxes noti'
			script += ', \\\n f(x, CN)'
		#
		if self.base == "angle":
			if self.option != "box":
				script += 'plot data u 1:($2/(3.142*sin(3.142*$1/180))) w l noti'
			else:
				script += 'set style fill solid 0.5\nset boxwidth ' + str(bin_width) + '\n'
				script += 'plot data u 1:($2/(3.142*sin(3.142*$1/180))) w boxes noti'

		elif self.option == "box":
			script += 'set style fill solid 0.5\nset boxwidth ' + str(bin_width) + '\n'
			script += '#\nplot data w boxes noti'
			
		return script


#############################################################################################
class MakeMulti:
	def __init__(self, cond_list, target_name):
		# cond_list = [base_name, data_list, Legend]
		self.list = cond_list[1]
		if target_name != '':
			self.dir = os.path.join(target_name, cond_list[0])
		else:
			self.dir = cond_list[0]

		self.base = cond_list[0]
		self.repeat = len(cond_list[1])
		#
		self.f_dat = cond_list[0] + ".dat"
		self.f_plt = cond_list[0] + ".plt"
		self.f_png = cond_list[0] + ".png"
		self.leg = cond_list[2]
	
	############################################################
	# マルチリストのグラフの作成
	def make_all(self):
		# データを書き出し 
		self.write_data()
		# グラフを作成
		self.make_graph()
		return

	##############################
	# データを書き出し 
	def write_data(self):
		os.makedirs(self.dir, exist_ok=True)
		with open(os.path.join(self.dir, self.f_dat), 'w') as f:
			f.write("# data:\n")
			if self.base != "CN_ave":
				for i, data in enumerate(self.list):
					f.write("\n\n# " + str(i) +":\n\n")
					for line in data:
						f.write(str(line[0]) + '\t' + str(line[1])  + '\n')
			else:
				for line in self.list:
					f.write(str(line[0]) + '\t' + str(line[1])  + '\n')
		return

	###########################################################
	# グラフを作成
	def make_graph(self):
		self.make_script()
		cwd = os.getcwd()
		os.chdir(self.dir)
		if platform.system() == "Windows":
			subprocess.call(self.f_plt, shell=True)
		elif platform.system() == "Linux":
			subprocess.call('gnuplot ' + self.f_plt, shell=True)
		os.chdir(cwd)
		return

	#######################
	# 必要なスクリプトを作成
	def make_script(self):
		with open(os.path.join(self.dir, self.f_plt), 'w') as f:
			script = self.script_content()
			f.write(script)
		return

	#################
	# スクリプトの中身
	def script_content(self):
		script = 'set term pngcairo font "Arial,14" \nset colorsequence classic \n'
		#
		script += '# \ndata = "' + self.f_dat + '" \nset output "' + self.f_png + ' "\n'
		#
		script += '#\nset size square\n#set xrange [1:]\n#set yrange [1:]\n'
		script += '#\nset xlabel "' + self.leg[0] + '"\nset ylabel "' + self.leg[1] + '"\n\n'
		#
		if self.base == "CN" or self.base == "CN_ave":
			script += '#\nset xrange [1:]\nset yrange [1:]\n'
			script += 'set key bottom\n\n'
			script += 'ct = 0.274\n'
			script += "f(x) = (1+ct)/(1-ct) -(2*ct*(1-ct**x))/(1-ct)**2/x\n\n"
			script += 'plot '
			if self.base == "CN":
				for i in range(self.repeat):
					script += 'data ind ' + str(i) + ' w l lc ' + str(i) + 'noti, \\\n'
			else:
				script += 'data w l ti "averaged", \\\n'
			script += 'f(x) w l lw 2 ti "theory"'
		else:
			script += 'plot '
			for i in range(self.repeat):
				script += 'data ind ' + str(i) + ' w l lc ' + str(i) + 'noti, \\\n'

		return script

#############################################################################################
class MakeSimple:
	def __init__(self, cond_list, target_name):
		# cond_list = [base_name, data_list, Legend, option]
		self.list = cond_list[1]
		if target_name != '':
			self.dir = os.path.join(target_name, cond_list[0])
		else:
			self.dir = cond_list[0]

		self.base = cond_list[0]
		self.repeat = len(cond_list[1])
		#
		self.f_dat = cond_list[0] + ".dat"
		self.f_plt = cond_list[0] + ".plt"
		self.f_png = cond_list[0] + ".png"
		self.leg = cond_list[2]
		# 
		self.option = cond_list[3]
	############################################################
	# マルチリストのグラフの作成
	def make_all(self):
		# データを書き出し 
		self.write_data()
		# グラフを作成
		self.make_graph()
		return

	##############################
	# データを書き出し 
	def write_data(self):
		if not os.path.exists(self.dir):
			os.mkdir(self.dir)
		with open(os.path.join(self.dir, self.f_dat), 'w') as f:
			f.write("# data:\n")
			for i, data in enumerate(self.list):
				f.write(str(data[0]) + '\t' + str(data[1])  + '\n')
		return

	###########################################################
	# グラフを作成
	def make_graph(self):
		self.make_script()
		cwd = os.getcwd()
		os.chdir(self.dir)
		if platform.system() == "Windows":
			subprocess.call(self.f_plt, shell=True)
		elif platform.system() == "Linux":
			subprocess.call('gnuplot ' + self.f_plt, shell=True)
		os.chdir(cwd)
		return

	#######################
	# 必要なスクリプトを作成
	def make_script(self):
		with open(os.path.join(self.dir, self.f_plt), 'w') as f:
			script = self.script_content()
			f.write(script)
		return

	#################
	# スクリプトの中身
	def script_content(self):
		script = 'set term pngcairo font "Arial,14" \nset colorsequence classic \n'
		#
		script += '# \ndata = "' + self.f_dat + '" \nset output "' + self.f_png + ' "\n'
		#
		script += '#\nset size square\n\n'
		script += '#\nset xlabel "' + self.leg[0] + '"\nset ylabel "' + self.leg[1] + '"\n\n'
		#
		if self.option == 'semilog' and self.base == 'NM':
			script += '#setxrange [0:1.0]\nset yrange [0.01:1]\n'
			script += '#\nset logscale y\nset format y "10^{%L}"\n\n'
			script += 'a=1\ntau=100\ns=100\ne=10000\n'
			script += 'fit [s:e] a*exp(-1*x/tau) data usi 1:2 via a,tau\n\n'
			script += 'set label 1 sprintf("Fitted {/Symbol t} = %.1e", tau) at graph 0.3, 0.8\n'
			script += 'set label 2 sprintf("Fitting Region: %d to %d", s, e) at graph 0.3, 0.7\n\n'
			script += 'plot data w l lt 1 noti, \\\n [s:e] a*exp(-1*x/tau) lt 2 noti\n\n'
		else:
			script += '#setxrange [0:1.0]\n#set yrange [0:100]\n'
			script += 'plot data w l noti\n\n'

		return script
