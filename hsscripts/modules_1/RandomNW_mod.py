#!/usr/bin/env python
# -*- coding: utf-8 -*-
################################################################################
import sys
import os
################################################################################
class InitialSetup:
	def __init__(self, sim_cond, restart, args):
		#
		self.nw_model = sim_cond[0]
		self.nw_type = sim_cond[1]
		self.n_segments = sim_cond[2]
		self.n_cell = sim_cond[3]
		self.multi_init = sim_cond[4]
		self.target_density = sim_cond[5]
		self.n_strand = sim_cond[6]
		self.l_bond = sim_cond[7]
		self.c_n = sim_cond[8]
		self.expand = sim_cond[9]
		#
		self.n_sc = 0
		#
		self.restart = restart
		self.input = args

	#######################
	# 基本となるデータを取得
	def get_base_data(self):
		# リスタートの場合のパスを設定
		read_file_path = self.res_flag()
		# 基本となるデータを取得
		target_cond = self.calc_conditions()
		multi_nw = target_cond[5]
		target_name = self.set_target_name(multi_nw)
		# # 基準となる8_Chainネットワークを設定
		# base_top_list = self.make_8chain_dic()
		
		return read_file_path, target_name, target_cond

	#########################
	#----- リスタート条件の分岐
	def res_flag(self):
		if self.restart == 1:
			print("##########\nRestart Calculation")
			if len(self.input) == 1:
				print("##########\nin Restart Calculation, Specify target directory!")
				print("usage: python", self.input[0], "Reading Dir")
				sys.exit("EXIT !!")
			if not os.path.exists(os.path.join(self.input[1], 'init.pickle')):
				exit("##########\ntarget directory does not exists.")
			elif self.n_cell != int(self.input[1].split('_')[2]):
				sys.exit("##########\nnumber of cells: n_cell is different from original Calculation.")
			elif self.n_strand != int(self.input[1].split('_')[0]):
				sys.exit("##########\nnumber of strands: n_strand is different from original Calculation.")
			else:
				read_file_path = self.input[1]
				return read_file_path

	###########################
	# 各種の条件を算出する
	def calc_conditions(self):
		structure = "Random_NW"
		#
		text = "#########################################" + "\n"
		text += "対象となるネットワークタイプ\t" + str(self.nw_type) + "\n"
		text += "分岐数\t\t\t\t" + str(self.n_strand) + "\n"
		text += "ストランド中のセグメント数\t" + str(self.n_segments) + "\n"
		text += "一辺当たりの単位ユニット数\t" + str(self.n_cell) + "\n"
		text += "#########################################" + "\n"
		#

		e2e = self.l_bond*((self.n_segments + 1)*self.c_n)**0.5				# 理想鎖状態での末端間距離
		n_beads_unit = 2 + self.n_segments*(1 + self.n_sc)*self.n_strand	# ユニットセル当たりの粒子数
		a_cell = (2*3**0.5)*e2e/3											# 理想鎖状態でのユニットセル長
		init_dens = n_beads_unit/a_cell**3									# 理想鎖状態でのユニットセル長
		
		print(self.nw_type)

		if self.nw_type == "KG_entangled":
			multi_nw = round(self.target_density/init_dens)
			density = multi_nw*n_beads_unit/a_cell**3					# 密度
			err_dens = (density/self.target_density - 1)*100
			system = a_cell*self.n_cell										# システムサイズ
			total_beads = multi_nw * n_beads_unit * self.n_cell**3		# 総セグメント数
			nu = multi_nw*self.n_strand/a_cell**3						# ストランドの数密度
			text += "設定密度:\t\t\t" + str(self.target_density) + "\n"
			text += "多重度:\t\t\t\t" + str(multi_nw) + "\n"
			text += "セグメント数:\t\t\t" + str(total_beads) + "\n"
			text += "システムサイズ:\t\t\t" + str(round(system, 4)) + "\n"
			text += "密度:\t\t\t\t" + str(round(density, 4)) + "\n"
			text += "密度の誤差:\t\t\t" + str(round(err_dens, 3)) + "%" + "\n"
			text += "ストランドの数密度:\t\t" + str(round(nu, 5)) + "\n"
			text += "#########################################" + "\n"
			print(text)
			if abs(err_dens) > 1:
				print(u"\n##### \n圧縮後の密度が、設定したい密度と 1% 以上違います。\nそれでも計算しますか？")
				self.prompt()
			else:
				self.prompt()
		elif self.nw_type == "KG_multi" or self.nw_type == "Sunuke":
			multi_nw = self.multi_init
			# density = self.multi_init*n_beads_unit/a_cell**3					# 密度
			# err_dens = (density/self.target_density - 1)*100
			total_beads = self.multi_init* n_beads_unit * self.n_cell**3		# 総セグメント数
			vol = total_beads/self.target_density										# システム体積
			system = vol**(1./3.)											# システムサイズ
			nu = self.multi_init*self.n_strand*self.n_cell**3/vol							# ストランドの数密度
			a_cell = system/self.n_cell
			text += "設定密度:\t\t\t" + str(self.target_density) + "\n"
			text += "多重度:\t\t\t\t" + str(self.multi_init) + "\n"
			text += "セグメント数:\t\t\t" + str(total_beads) + "\n"
			text += "システムサイズ:\t\t\t" + str(round(system, 4)) + "\n"
			# text += "密度:\t\t\t\t" + str(round(density, 4)) + "\n"
			# text += "密度の誤差:\t\t\t" + str(round(err_dens, 3)) + "%" + "\n"
			text += "ストランドの数密度:\t\t" + str(round(nu, 5)) + "\n"
			text += "#########################################" + "\n"
			print(text)
			# if abs(err_dens) > 1:
			# 	print(u"\n##### \n圧縮後の密度が、設定したい密度と 1% 以上違います。\nそれでも計算しますか？")
			# 	self.prompt()
			# else:
			# 	self.prompt()
		elif self.nw_type == "NPT":
			multi_nw = round(self.target_density/init_dens)
			density = multi_nw*n_beads_unit/a_cell**3					# 密度
			err_dens = (density/self.target_density - 1)*100
			total_beads = multi_nw * n_beads_unit * self.n_cell**3		# 総セグメント数
			vol = total_beads/density										# システム体積
			system = vol**(1./3.)											# システムサイズ
			nu = multi_nw*self.n_strand*self.n_cell**3/vol							# ストランドの数密度
			a_cell = system/self.n_cell
			text += "設定密度:\t\t\t" + str(self.target_density) + "\n"
			text += "多重度:\t\t\t\t" + str(multi_nw) + "\n"
			text += "セグメント数:\t\t\t" + str(total_beads) + "\n"
			text += "システムサイズ:\t\t\t" + str(round(system, 4)) + "\n"
			text += "密度:\t\t\t\t" + str(round(density, 4)) + "\n"
			text += "密度の誤差:\t\t\t" + str(round(err_dens, 3)) + "%" + "\n"
			text += "ストランドの数密度:\t\t" + str(round(nu, 5)) + "\n"
			text += "###\n"
			text += "初期の拡大率:\t\t\t" + str(round(self.expand, 3)) + "\n"
			text += "#########################################" + "\n"
			print(text)
			if abs(err_dens) > 1:
				print(u"\n##### \n圧縮後の密度が、設定したい密度と 1% 以上違います。\nそれでも計算しますか？")
				self.prompt()
			else:
				self.prompt()
		elif self.nw_type == "KG_gel":
			multi_nw = self.multi_init
			density = multi_nw*n_beads_unit/a_cell**3					# 密度
			err_dens = (density/self.target_density - 1)*100
			total_beads = multi_nw * n_beads_unit * self.n_cell**3		# 総セグメント数
			vol = total_beads/density										# システム体積
			system = vol**(1./3.)											# システムサイズ
			nu = multi_nw*self.n_strand*self.n_cell**3/vol							# ストランドの数密度
			a_cell = system/self.n_cell
			text += "設定密度:\t\t\t" + str(self.target_density) + "\n"
			text += "多重度:\t\t\t\t" + str(multi_nw) + "\n"
			text += "セグメント数:\t\t\t" + str(total_beads) + "\n"
			text += "システムサイズ:\t\t\t" + str(round(system, 4)) + "\n"
			text += "密度:\t\t\t\t" + str(round(density, 4)) + "\n"
			text += "密度の誤差:\t\t\t" + str(round(err_dens, 3)) + "%" + "\n"
			text += "ストランドの数密度:\t\t" + str(round(nu, 5)) + "\n"
			text += "#########################################" + "\n"
			print(text)
			if abs(err_dens) > 1:
				print(u"\n##### \n圧縮後の密度が、設定したい密度と 1% 以上違います。\nそれでも計算しますか？")
				self.prompt()
			else:
				self.prompt()
		#
		with open("calc_conditions.txt", 'w') as f:
			f.write(text)
		#
		n_solvent = 0
		target_cond = [system, a_cell, total_beads, nu, structure, multi_nw, n_solvent]
		return target_cond

	###############
	# 計算条件の確認
	def prompt(self):
		dic={'y':True,'yes':True,'n':False,'no':False}
		while True:
			inp = input(u'計算を続行 ==> [Y]es >> ').lower()
			if inp in dic:
				inp = dic[inp]
				break
			print(u'##### \nもう一度入力してください。')
		if inp :
			print(u"計算を続行します")
		else:
			sys.exit(u"##### \n計算を中止します。")
		return
	
	###################
	# target_nameを決める。
	def set_target_name(self, multi_nw):
		target_name = 'Random_' + str(self.nw_type ) + "_" + str(self.n_strand) + '_Chains_N_' + str(self.n_segments) + "_C_" + str(self.n_cell) + "_M_" + str(multi_nw)
		return target_name
