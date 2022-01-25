#!/usr/bin/env python
# -*- coding: utf-8 -*-
#######################################################
import os
import sys
import codecs
import numpy as np
from UDFManager import UDFManager
#######################################################
#
def setupcondition():
	if not os.path.isfile('./calc_condition.udf'):
		print()
		print('In this directory, no "calc_condition.udf" is found !')
		print('New one will be generated.')
		print('Please, modify and save it !\n')
		makenewudf()
		input('Press ENTER to continue...')

	dic={'y':True,'yes':True,'q':False,'quit':False}
	while True:
		basic_cond, sim_cond, names = readconditionudf()

		init = InitialSetup(sim_cond)
		target_cond, condition_text = init.init_calc()
		print('Change UDF: type [r]eload')
		print('Quit input process: type [q]uit')
		inp = input('Condition is OK ==> [y]es >> ').lower()
		if inp in dic:
			inp = dic[inp]
			break
		print('##### \nRead Condition UDF again \n#####\n\n')
	if inp :
		print("\n\nSetting UP progress !!")
	else:
		sys.exit("##### \nQuit !!")
	
	target_name = init.set_target_name()
	# ネットワークの計算
	calcd_data_dic = init.calc_all()

	return basic_cond, sim_cond, names, target_cond, target_name, calcd_data_dic, condition_text

###########################################
# 
def readconditionudf():
	u = UDFManager('calc_condition.udf')

	u.jump(-1)
	## 計算条件
	##################
	# 使用するCognacのバージョン
	ver_cognac = u.get('CalcCond.Cognac_ver')
	# 計算に使用するコア数
	core = u.get('CalcCond.Cores')
	# ベースとするUDFの名前
	base_udf = "base_uin.udf"
	blank_udf = ver_cognac + '.udf'
	###########
	basic_cond = [ver_cognac, blank_udf, base_udf, core]
	#######################################################
	## 計算ターゲット
	###################
	## Networkモデルの設定
	nw_model = u.get('SimulationCond.Model.TargetModel')
	###################
	## Networkモデルの設定
	if nw_model == "Regular_NW":
		strand = u.get('SimulationCond.Model.Regular_NW.chains')
	elif nw_model == "Random_NW":
		strand = u.get('SimulationCond.Model.Random_NW.chains')
		calc = u.get('SimulationCond.Model.Random_NW.Calc_Topolpgy')
		if calc == 'Read':
			restartfile = u.get('SimulationCond.Model.Random_NW.Read.file_name')
			if os.path.exists(restartfile):
				print('Read existing file!')
			else:
				sys.exit('Cannnot find Read file !!')
		elif calc == 'Calc':
			cond_top = u.get('SimulationCond.Model.Random_NW.Calc')
	###################
	## ポリマー鎖の設定
	sim_type = u.get('SimulationCond.Type.SimType')
	#################################
	if sim_type == "Entangled":
		n_segments = u.get('SimulationCond.Type.Entangled.N_Segments')
		n_cell = u.get('SimulationCond.Type.Entangled.N_UnitCells')
		multi_init = 0
		nv = 1.0
		expand = 1.0
		step_press = []
	elif sim_type == "NPT":
		n_segments = u.get('SimulationCond.Type.NPT.N_Segments')
		n_cell = u.get('SimulationCond.Type.NPT.N_UnitCells')
		multi_init = 0
		nv = 1.0
		expand = u.get('SimulationCond.Type.NPT.ExpansionRatio')
		step_press = u.get('SimulationCond.Type.NPT.StepPress[]')
	elif sim_type == "Multi":
		n_segments = u.get('SimulationCond.Type.Multi.N_Segments')
		n_cell = u.get('SimulationCond.Type.Multi.N_UnitCells')
		multi_init = u.get('SimulationCond.Type.Multi.Multiplicities')
		nv = 1.0
		expand = u.get('SimulationCond.Type.Multi.ExpansionRatio')
		step_press = u.get('SimulationCond.Type.Multi.StepPress[]')
	elif sim_type == "Multi_entangled":
		n_segments = u.get('SimulationCond.Type.Multi_entangled.N_Segments')
		n_cell = u.get('SimulationCond.Type.Multi_entangled.N_UnitCells')
		multi_init = u.get('SimulationCond.Type.Multi_entangled.Multiplicities')
		nv = 1.0
		expand = 1.0
		step_press = []
	elif sim_type == "Gel":
		n_segments = u.get('SimulationCond.Type.Gel.N_Segments')
		n_cell = u.get('SimulationCond.Type.Gel.N_UnitCells')
		multi_init = u.get('SimulationCond.Type.Gel.Multiplicities')
		expand = u.get('SimulationCond.Type.Gel.ExpansionRatio')
		step_press = u.get('SimulationCond.Type.Gel.StepPress[]')
		nv = 1.0
	elif sim_type == "Gel_entangled":
		n_segments = u.get('SimulationCond.Type.Gel_entangled.N_Segments')
		n_cell = u.get('SimulationCond.Type.Gel_entangled.N_UnitCells')
		multi_init = u.get('SimulationCond.Type.Gel_entangled.Multiplicities')
		nv = 1.0
		expand = 1.0
		step_press = []
	elif sim_type == "Gel_concd":
		n_segments = u.get('SimulationCond.Type.Gel_concd.N_Segments')
		n_cell = u.get('SimulationCond.Type.Gel_concd.N_UnitCells')
		multi_init = u.get('SimulationCond.Type.Gel_concd.Multiplicities')
		nv = u.get('SimulationCond.Type.Gel_concd.NV')
		expand = u.get('SimulationCond.Type.Gel_concd.ExpansionRatio')
		step_press = u.get('SimulationCond.Type.Gel_concd.StepPress[]')
	elif sim_type == "Gel_concd_entangled":
		n_segments = u.get('SimulationCond.Type.Gel_concd_entangled.N_Segments')
		n_cell = u.get('SimulationCond.Type.Gel_concd_entangled.N_UnitCells')
		multi_init = u.get('SimulationCond.Type.Gel_concd_entangled.Multiplicities')
		nv = u.get('SimulationCond.Type.Gel_concd_entangled.NV')
		expand = 1.0
		step_press = []
	else:
		multi_init = 0
	# 設定密度
	target_density = u.get('SimulationCond.target_density')
	########################################################
	if strand == "3_Chain_S" or strand == "3_Chain_D":
		n_strand = 3
	elif strand == "4_Chain":
		n_strand = 4
	elif strand == "6_Chain":
		n_strand = 6
	elif strand == "8_Chain":
		n_strand = 8
	#
	l_bond = 0.97
	if n_segments <= 10:
		c_n = 1.5
	elif n_segments <= 20:
		c_n = 1.65
	elif n_segments <= 40:
		c_n = 1.7
	else:
		c_n = 1.75
	#
	sim_cond = [nw_model, strand, sim_type, n_segments, n_cell, multi_init, target_density, expand, n_strand, l_bond, c_n, step_press, nv]

	####################################################
	# Cognac用の名称設定
	nw_name = "Network"
	atom_name = ["JP_A", "End_A", "Strand_A", "Side_A", "Solvent"]
	bond_name = ["bond_JP-Chn", "bond_Strand", "bond_Side"]
	angle_name = ["angle_AAA"]
	site_name = ["site_JP", "site_End", "site_Strand", "site_Solvent"]
	pair_name = ["site_JP-site_JP", "site_Strand-site_JP", "site_Strand-site_Strand", 
					"site_JP-site_End", "site_Strand-site_End", "site_End-site_End",
					"site_Solvent-site_Solvent", "site_Solvent-site_JP", "site_Solvent-site_End",
					"site_Solvent-site_Strand"]
	site_pair_name = [ 
					["site_JP", "site_JP"], 
					["site_Strand", "site_JP"], 
					["site_Strand", "site_Strand"],
					["site_JP", "site_End"], 
					["site_Strand", "site_End"], 
					["site_End", "site_End"],
					["site_Solvent", "site_Solvent"],
					["site_Solvent", "site_JP"],
					["site_Solvent", "site_End"],
					["site_Solvent", "site_Strand"],
					]
	names = [nw_name, atom_name, bond_name, angle_name, site_name, pair_name, site_pair_name]

	return basic_cond, sim_cond, names

#################################
#
def makenewudf():
	contents = '''
	\\begin{def}
		CalcCond:{
			Cognac_ver:select{
				"cognac101", 
				"cognac112"
				}, // Cognac version
			Cores: int // 計算に使用するコア数を指定
			}

		SimulationCond:{
			Model:{
				TargetModel:select{ 
						"Regular_NW",
						"Random_NW"
					},
					Regular_NW:{
						chains:select{
							"3_Chain_S", 
							"3_Chain_D", 
							"4_Chain", 
							"6_Chain", 
							"8_Chain"
						}
					},
					Random_NW:{
						chains:select{
						"3_Chain", 
						"4_Chain", 
						"5_Chain", 
						"6_Chain", 
						"7_Chain", 
						"8_Chain"
						},
						Calc_Topolpgy:select{
						"Calc",
						"Read"
						},
						Calc:{
							pre_sampling:int,
							pre_try:int,
							sampling:int,
							try:int
							},
						Read:{
							file_name:string
						}
					}
				}
			Type:{
				SimType:select{
					"Entangled",   // 密度、末端間距離を設定値に合わせるように多重度を変化。 絡み合いが入るように初期化
					"NPT",		    // 密度、末端間距離を設定値に合わせるように多重度を変化。 絡み合いが入らないようにNPTで縮める
					"Multi",    	// 設定した多重度で、密度を設定値になるように、システムサイズを縮める
					"Multi_entangled",	// 設定した多重度で、密度を設定値になるように、システムサイズを縮め、絡み合いが入るように初期化
					"Gel",    	// 設定した多重度で、密度、末端間距離を設定値に合わせるように溶剤を添加
					"Gel_entangled",	// 設定した多重度で、密度、末端間距離を設定値に合わせるように溶剤を添加、絡み合いが入るように初期化
					"Gel_concd",	// 設定した多重度で、溶剤量変化
					"Gel_concd_entangled"	// 設定した多重度で、溶剤量変化、絡み合いが入るように初期化
					},
					Entangled:{
						N_Segments: int,
						N_UnitCells: int,
					},
					NPT:{
						N_Segments: int,
						N_UnitCells: int,
						ExpansionRatio: float,
						StepPress[]: float
					},
					Multi:{
						N_Segments: int,
						N_UnitCells: int,
						Multiplicities: int,
						ExpansionRatio: float,
						StepPress[]: float
					},
					Multi_entangled:{
						N_Segments: int,
						N_UnitCells: int,
						Multiplicities: int
					},
					Gel:{
						N_Segments: int,
						N_UnitCells: int,
						Multiplicities: int,
						ExpansionRatio: float,
						StepPress[]: float
					},
					Gel_entangled:{
						N_Segments: int,
						N_UnitCells: int,
						Multiplicities: int
					},
					Gel_concd:{
						N_Segments: int,
						N_UnitCells: int,
						Multiplicities: int
						NV: float,
						ExpansionRatio: float,
						StepPress[]: float
					},
					Gel_concd_entangled:{
						N_Segments: int,
						N_UnitCells: int,
						Multiplicities: int,
						NV: float
					}
				}
			target_density:float
			}
	\end{def}

	\\begin{data}
		CalcCond:{
			"cognac112",
			1
			}

		SimulationCond:{
			{"Regular_NW",
				{"4_Chain"},
				{"4_Chain",
					"Calc",
						{1000, 100, 1000, 1000},
						{""}
				}
			}
			{"NPT",
				{20, 3},
				{20, 3, 2.0, [0.2, 0.5, 1.0, 2.0, 3.0, 4.5]},
				{20, 3, 1, 2.0, [0.2, 0.5, 1.0, 2.0, 5.0, 6.5, 7.0]},
				{20, 3, 1},
				{20, 3, 1, 2.0, [0.2, 0.5, 1.0, 2.0, 5.0, 6.5, 7.0]},
				{20, 3, 2},
				{20, 3, 2, 0.5, 2.0, [0.2, 0.5, 1.0, 2.0, 5.0, 6.5, 7.0]},
				{20, 3, 2, 0.5}
			}
		0.85
		}
	\end{data}
	'''
	###
	with codecs.open('./calc_condition.udf', 'w', 'utf_8') as f:
		f.write(contents)
	
	return


######################################
##### ネットワークの初期設定を行う #####
######################################
class InitialSetup:
	def __init__(self, sim_cond):
		#
		self.nw_model = sim_cond[0]
		self.strand = sim_cond[1]
		self.nw_type = sim_cond[2]
		self.n_segments = sim_cond[3]
		self.n_cell = sim_cond[4]
		self.multi_init = sim_cond[5]
		self.target_density = sim_cond[6]
		self.expand = sim_cond[7]
		self.n_strand = sim_cond[8]
		self.l_bond = sim_cond[9]
		self.c_n = sim_cond[10]
		self.step_press = sim_cond[11]
		self.nv = sim_cond[12]
	##########################################
	##### ネットワークポリマーの諸量を計算 ######
	##########################################
	#-----ネットワークポリマーの諸量を計算
	def calc_conditions(self):
		## 計算システムの諸量を計算して、出力
		target_cond, condition_text = self.init_calc()
		target_name = self.set_target_name()	
		return target_cond, target_name, condition_text

	################################################################################
	def init_calc(self):
		structure = "Regular_NW"
		e2e = self.l_bond*((self.n_segments + 1)*self.c_n)**0.5	# 理想鎖状態での末端間距離
		#
		if self.strand == "3_Chain_S":
			n_chains = 12						        # サブチェインの本数
			n_beads_unit = 8 + self.n_segments*n_chains		# ユニットセル当たりの粒子数
			a_cell = (2*2**0.5)*e2e				        # 理想鎖状態でのユニットセル長
		elif self.strand == "3_Chain_D":
			n_chains = 24						        # サブチェインの本数
			n_beads_unit = 16 + self.n_segments*n_chains	# ユニットセル当たりの粒子数
			a_cell = (2*2**0.5)*e2e				        # 理想鎖状態でのユニットセル長
		elif self.strand == "4_Chain":
			n_chains = 16						        # サブチェインの本数
			n_beads_unit = 8 + self.n_segments*n_chains		# ユニットセル当たりの粒子数
			a_cell = (4*3**0.5)*e2e/3			        # 理想鎖状態でのユニットセル長
		elif self.strand == "6_Chain":
			n_chains = 3						        # サブチェインの本数
			n_beads_unit = 1 + self.n_segments*n_chains		# ユニットセル当たりの粒子数
			a_cell = e2e						            # 理想鎖状態でのユニットセル長
		elif self.strand == "8_Chain":
			n_chains = 8						        # サブチェインの本数
			n_beads_unit = 2 + self.n_segments*n_chains     # ユニットセル当たりの粒子数
			a_cell = (2*3**0.5)*e2e/3					# 理想鎖状態でのユニットセル長
		#
		n_solvent = 0
		#
		if self.nw_type == "Entangled" or self.nw_type == "NPT":
			unit_cell = a_cell
			self.multi = round(self.target_density*a_cell**3/n_beads_unit)		# 密度を設定値とした場合に必要な多重度
			fin_dens = n_beads_unit*self.multi/a_cell**3						# 上記の多重度での密度
			err_dens = round((fin_dens/self.target_density - 1)*100, 2) 	# 設定密度との誤差
			single_net_atom = int(n_beads_unit*self.n_cell**3.)	    			# 一つのネットワーク中の粒子数
			total_net_atom = int(self.multi*single_net_atom)    			# 全システム中のネットワーク粒子数
			system = (total_net_atom/self.target_density)**(1/3)			# その際のシステムサイズ
			nu = n_chains*self.multi/a_cell**3.								# ストランドの数密度
		elif self.nw_type == "Multi" or self.nw_type == "Multi_entangled":
			self.multi = self.multi_init
			org_system = a_cell*self.n_cell							# e2e から決めたシステムサイズ
			single_net_atom = int(n_beads_unit*self.n_cell**3.)	    # 一つのネットワーク中の粒子数
			total_net_atom = int(self.multi*single_net_atom)    # 全システム中のネットワーク粒子数
			vol = total_net_atom/self.target_density			# システム体積
			system = vol**(1./3.)							# 収縮後のシステムサイズ
			shrinkage = system/org_system						# 収縮比
			mod_e2e = shrinkage*e2e								# 収縮後の末端間距離
			unit_cell = shrinkage*a_cell
			nu = n_chains*self.multi*self.n_cell**3./vol		# ストランドの数密度
		elif self.nw_type == "Gel" or self.nw_type == "Gel_entangled":
			unit_cell = a_cell
			self.multi = self.multi_init
			system = a_cell*self.n_cell						    # システムサイズ
			single_net_atom = int(n_beads_unit*self.n_cell**3.)	    # 一つのネットワーク中の粒子数
			total_net_atom = int(self.multi*single_net_atom)    # 全システム中のネットワーク粒子数
			vol = system**3									    # システム体積
			# density = total_net_atom/vol                        # 密度
			#
			total_atom = int(vol*self.target_density)		    # システム中の粒子総数
			n_solvent = int(total_atom - total_net_atom)		# 溶媒の粒子数
			nv = total_net_atom/total_atom						# ネットワークの体積分率
			nu = n_chains*self.multi*self.n_cell**3./vol		# ストランドの数密度
		elif self.nw_type == "Gel_concd" or self.nw_type == "Gel_concd_entangled":
			unit_cell = a_cell
			self.multi = self.multi_init
			single_net_atom = int(n_beads_unit*self.n_cell**3.)	    # 一つのネットワーク中の粒子数
			total_net_atom = int(self.multi*single_net_atom)    # 全システム中のネットワーク粒子数
			n_solvent = int(total_net_atom*(1.0 - self.nv)/self.nv)		# 溶媒の粒子数
			total_atom = int(total_net_atom + n_solvent)		    # システム中の粒子総数
			vol = total_atom/self.target_density			# システム体積
			system = vol**(1./3.)							# システムサイズ
			total_atom = int(vol*self.target_density)		    # システム中の粒子総数
			nv = self.nv						# ネットワークの体積分率
			nu = n_chains*self.multi*self.n_cell**3./vol		# ストランドの数密度
		else:
			sys.exit("Something Wrong!!")
		#
		text = "#########################################" + "\n"
		text += "ネットワークモデル\t\t" + str(self.strand) + "\n"
		text += "ネットワークタイプ\t\t" + str(self.nw_type) + "\n"
		text += "ストランド中のセグメント数\t" + str(self.n_segments) + "\n"
		text += "特性比:\t\t\t" + str(round(self.c_n,2)) + "\n"
		text += "末端間距離:\t\t\t" + str(round(e2e,5)) + "\n"
		text += "一辺当たりの単位ユニット数\t" + str(self.n_cell) + "\n"
		text += "#########################################" + "\n"

		if self.nw_type == "Entangled" or self.nw_type == "NPT":
			text += "設定密度:\t\t\t" + str(self.target_density) + "\n"
			text += "多重度:\t\t\t\t" + str(self.multi) + "\n"
			text += "セグメント数:\t\t\t" + str(total_net_atom) + "\n"
			text += "システムサイズ:\t\t\t" + str(round(system, 4)) + "\n"
			text += "密度:\t\t\t\t" + str(round(fin_dens, 4)) + "\n"
			text += "密度の誤差:\t\t\t" + str(round(err_dens, 3)) + "%" + "\n"
			text += "ストランドの数密度:\t\t" + str(round(nu, 5)) + "\n"
			text += "#########################################" + "\n"
			print(text)
			if abs(err_dens) > 1:
				print(u"\n##### \n圧縮後の密度が、設定したい密度と 1% 以上違います。\nそれでも計算しますか？")
		elif self.nw_type == "Multi":
			text += "当初のシステムサイズ:\t\t" + str(round(org_system, 4)) + "\n"
			text += "多重度:\t\t\t\t" + str(self.multi) + "\n"
			text += "設定密度:\t\t\t" + str(self.target_density) + "\n"
			text += "セグメント数:\t\t\t" + str(total_net_atom) + "\n"
			text += "収縮比:\t\t\t\t" + str(round(shrinkage, 4)) + "\n"
			text += "ステップ圧力:\t\t" + ', '.join(map(str, self.step_press)) + "\n"
			text += "収縮後の末端間距離:\t\t" + str(round(mod_e2e,5)) + "\n"
			text += "収縮後のシステムサイズ:\t\t" + str(round(system, 4)) + "\n"
			text += "ストランドの数密度:\t\t" + str(round(nu, 5)) + "\n"
			text += "#########################################" + "\n"
			print(text)
		elif self.nw_type == "Multi_entangled":
			text += "当初のシステムサイズ:\t\t" + str(round(org_system, 4)) + "\n"
			text += "多重度:\t\t\t\t" + str(self.multi) + "\n"
			text += "設定密度:\t\t\t" + str(self.target_density) + "\n"
			text += "セグメント数:\t\t\t" + str(total_net_atom) + "\n"
			text += "収縮比:\t\t\t\t" + str(round(shrinkage, 4)) + "\n"
			text += "収縮後の末端間距離:\t\t" + str(round(mod_e2e,5)) + "\n"
			text += "収縮後のシステムサイズ:\t\t" + str(round(system, 4)) + "\n"
			text += "ストランドの数密度:\t\t" + str(round(nu, 5)) + "\n"
			text += "#########################################" + "\n"
			print(text)
		elif self.nw_type == "Gel":
			text += "設定密度:\t\t\t" + str(self.target_density) + "\n"
			text += "多重度:\t\t\t\t" + str(self.multi) + "\n"
			text += "単一ネットワークのセグメント数:\t" + str(single_net_atom) + "\n"
			text += "溶媒のセグメント数:\t\t" + str(n_solvent) + "\n"
			text += "システムサイズ:\t\t\t" + str(round(system, 4)) + "\n"
			text += "ステップ圧力:\t\t" + ', '.join(map(str, self.step_press)) + "\n"
			text += "ネットワークの体積分率:\t\t" + str(round(nv, 4)) + "\n"
			text += "ストランドの数密度:\t\t" + str(round(nu, 5)) + "\n"
			text += "#########################################" + "\n"
			print(text)
		elif self.nw_type == "Gel_entangled":
			text += "設定密度:\t\t\t" + str(self.target_density) + "\n"
			text += "多重度:\t\t\t\t" + str(self.multi) + "\n"
			text += "単一ネットワークのセグメント数:\t" + str(single_net_atom) + "\n"
			text += "溶媒のセグメント数:\t\t" + str(n_solvent) + "\n"
			text += "システムサイズ:\t\t\t" + str(round(system, 4)) + "\n"
			text += "ネットワークの体積分率:\t\t" + str(round(nv, 4)) + "\n"
			text += "ストランドの数密度:\t\t" + str(round(nu, 5)) + "\n"
			text += "#########################################" + "\n"
			print(text)
		elif self.nw_type == "Gel_concd":
			text += "設定密度:\t\t\t" + str(self.target_density) + "\n"
			text += "多重度:\t\t\t\t" + str(self.multi) + "\n"
			text += "単一ネットワークのセグメント数:\t" + str(single_net_atom) + "\n"
			text += "溶媒のセグメント数:\t\t" + str(n_solvent) + "\n"
			text += "システムサイズ:\t\t\t" + str(round(system, 4)) + "\n"
			text += "ステップ圧力:\t\t" + ', '.join(map(str, self.step_press)) + "\n"
			text += "ネットワークの体積分率:\t\t" + str(round(nv, 4)) + "\n"
			text += "ストランドの数密度:\t\t" + str(round(nu, 5)) + "\n"
			text += "#########################################" + "\n"
			print(text)
		elif self.nw_type == "Gel_concd_entangled":
			text += "設定密度:\t\t\t" + str(self.target_density) + "\n"
			text += "多重度:\t\t\t\t" + str(self.multi) + "\n"
			text += "単一ネットワークのセグメント数:\t" + str(single_net_atom) + "\n"
			text += "溶媒のセグメント数:\t\t" + str(n_solvent) + "\n"
			text += "システムサイズ:\t\t\t" + str(round(system, 4)) + "\n"
			text += "ネットワークの体積分率:\t\t" + str(round(nv, 4)) + "\n"
			text += "ストランドの数密度:\t\t" + str(round(nu, 5)) + "\n"
			text += "#########################################" + "\n"
			print(text)
		#
		with open("calc_conditions.txt", 'w') as f:
			f.write(text)
		#
		target_cond = [system, unit_cell, total_net_atom, nu, structure, self.multi, n_solvent]
		return target_cond, text

	###################
	# target_nameを決める。
	def set_target_name(self):
		target_name = str(self.nw_type ) + "_" + str(self.strand) + '_N_' + str(self.n_segments) + "_Cells_" + str(self.n_cell) + "_Multi_" + str(self.multi)
		return target_name





	################################################################################
	# ネットワーク設定の計算
	################################################################################
	def calc_all(self):
		# 架橋点 JP を設定
		jp_xyz, subchain_se_xyz = self.calc_jp_subChains()
		#
		calcd_data_dic = self.set_atom(jp_xyz, subchain_se_xyz)
		return calcd_data_dic
	################################################################################
	# JPおよびサブチェインの始点と終点のXYZを設定
	def calc_jp_subChains(self):
		# jp_xyz は、JPの座標のリスト
		# subchain_se_xyz は、サブチェインの出発点と終点のリスト
		if self.strand == "3_Chain_S":
			# JPを設定
			jp_xyz = [
			[
			[0, 0, 0],
			[0, 0.25, 0.25],
			[0.25, 0.25, 0.5],
			[0.25, 0, 0.75],
			[0.5, 0.5, 0.5],
			[0.5, 0.75, 0.75],
			[0.75, 0.5, 0.25],
			[0.75, 0.75, 0]
			]
			]
			# サブチェインの出発点と終点を設定
			subchain_se_xyz = [
			[
			[[0, 0, 0], [0, 0.25, 0.25]],
			[[0, 0.25, 0.25], [0.25, 0.25, 0.5]],
			[[0.25, 0.25, 0.5], [0.25, 0, 0.75]],
			[[0.25, 0.25, 0.5], [0.5, 0.5, 0.5]],
			[[0.5, 0.5, 0.5], [0.5, 0.75, 0.75]],
			[[0.5, 0.5, 0.5], [0.75, 0.5, 0.25]],
			[[0.75, 0.5, 0.25], [0.75, 0.75, 0]],
			[[0.75, 0.5, 0.25], [1, 0.25, 0.25]],
			[[0.25, 0, 0.75], [0, 0, 1]],
			[[0.5, 0.75, 0.75], [0.25, 1, 0.75]],
			[[0.5, 0.75, 0.75], [0.75, 0.75, 1]],
			[[0.75, 0.75, 0], [1, 1, 0]]
			]
			]

		elif self.strand == "3_Chain_D":
			# JPを設定
			jp_xyz = [
			[
			[0, 0, 0],
			[0, 0.25, 0.25],
			[0.25, 0.25, 0.5],
			[0.25, 0, 0.75],
			[0.5, 0.5, 0.5],
			[0.5, 0.75, 0.75],
			[0.75, 0.5, 0.25],
			[0.75, 0.75, 0]
			],
			[	#ここから二つ目
			[0, 0.5, 0.75],
			[0, 0.75, 0.5],
			[0.25, 0.75, 0.25],
			[0.25, 0.5, 0],
			[0.5, 0.25, 0],
			[0.5, 0, 0.25],
			[0.75, 0, 0.5],
			[0.75, 0.25, 0.75]
			]
			]
			# サブチェインの出発点と終点を設定
			subchain_se_xyz = [
			[
			[[0, 0, 0], [0, 0.25, 0.25]],
			[[0, 0.25, 0.25], [0.25, 0.25, 0.5]],
			[[0.25, 0.25, 0.5], [0.25, 0, 0.75]],
			[[0.25, 0.25, 0.5], [0.5, 0.5, 0.5]],
			[[0.5, 0.5, 0.5], [0.5, 0.75, 0.75]],
			[[0.5, 0.5, 0.5], [0.75, 0.5, 0.25]],
			[[0.75, 0.5, 0.25], [0.75, 0.75, 0]],
			[[0.75, 0.5, 0.25], [1, 0.25, 0.25]],
			[[0.25, 0, 0.75], [0, 0, 1]],
			[[0.5, 0.75, 0.75], [0.25, 1, 0.75]],
			[[0.5, 0.75, 0.75], [0.75, 0.75, 1]],
			[[0.75, 0.75, 0], [1, 1, 0]]
			],
			[
			[[0, 0.5, 0.75], [0, 0.75, 0.5]],
			[[0, 0.75, 0.5], [0.25, 0.75, 0.25]],
			[[0.25, 0.75, 0.25], [0.25, 0.5, 0]],
			[[0.25, 0.5, 0], [0.5, 0.25, 0]],
			[[0.5, 0.25, 0], [0.5, 0, 0.25]],
			[[0.5, 0, 0.25], [0.75, 0, 0.5]],
			[[0.75, 0, 0.5], [0.75, 0.25, 0.75]],
			[[0, 0.5, 0.75], [0.25, 0.5, 1]],
			[[0.25, 0.75, 0.25], [0.5, 1, 0.25]],
			[[0.75, 0.25, 0.75], [0.5, 0.25, 1]],
			[[0.75, 0.25, 0.75], [1, 0.5, 0.75]],
			[[0.75, 1, 0.5], [1, 0.75, 0.5]]
			]
			]

		elif self.strand == "4_Chain":
			# JPを設定
			jp_xyz = [
			[
			[0.0, 0.0, 0.0],
			[0.0, 0.5, 0.5],
			[0.5, 0.0, 0.5],
			[0.5, 0.5, 0.0],
			[0.25, 0.25, 0.25],
			[0.25, 0.75, 0.75],
			[0.75, 0.25, 0.75],
			[0.75, 0.75, 0.25]
			]
			]
			# サブチェインの出発点と終点を設定
			subchain_se_xyz = [
			[
			[[0.25, 0.25, 0.25], [0.0, 0.0, 0.0]],	# No.1
			[[0.25, 0.25, 0.25], [0.0, 0.5, 0.5]],
			[[0.25, 0.25, 0.25], [0.5, 0.0, 0.5]],
			[[0.25, 0.25, 0.25], [0.5, 0.5, 0.0]],
			[[0.25, 0.75, 0.75], [0.0, 0.5, 0.5]],	# No.2
			[[0.25, 0.75, 0.75], [0.0, 1.0, 1.0]],
			[[0.25, 0.75, 0.75], [0.5, 0.5, 1.0]],
			[[0.25, 0.75, 0.75], [0.5, 1.0, 0.5]],
			[[0.75, 0.25, 0.75], [0.5, 0.0, 0.5]],	# No.3
			[[0.75, 0.25, 0.75], [0.5, 0.5, 1.0]],
			[[0.75, 0.25, 0.75], [1.0, 0.0, 1.0]],
			[[0.75, 0.25, 0.75], [1.0, 0.5, 0.5]],
			[[0.75, 0.75, 0.25], [0.5, 0.5, 0.0]],	# No.4
			[[0.75, 0.75, 0.25], [0.5, 1.0, 0.5]],
			[[0.75, 0.75, 0.25], [1.0, 0.5, 0.5]],
			[[0.75, 0.75, 0.25], [1.0, 1.0, 0.0]]
			]
			]

		elif self.strand == "6_Chain":
			# JPを設定
			jp_xyz = [
			[
			[0.,0.,0.]
			]
			]
			# サブチェインの出発点と終点を設定
			subchain_se_xyz = [
			[
			[[0., 0., 0.], [1, 0, 0]],
			[[0., 0., 0.], [0, 1, 0]],
			[[0., 0., 0.], [0, 0, 1]]
			]
			]

		elif self.strand == "8_Chain":
			# JPを設定
			jp_xyz = [
			[
			[0.,0.,0.],
			[0.5,0.5,0.5]
			]
			]
			# サブチェインの出発点と終点を設定
			subchain_se_xyz = [
			[
			[[0.5, 0.5, 0.5], [0, 0, 0]],
			[[0.5, 0.5, 0.5], [1, 0, 0]],
			[[0.5, 0.5, 0.5], [0, 1, 0]],
			[[0.5, 0.5, 0.5], [1, 1, 0]],
			[[0.5, 0.5, 0.5], [0, 0, 1]],
			[[0.5, 0.5, 0.5], [1, 0, 1]],
			[[0.5, 0.5, 0.5], [0, 1, 1]],
			[[0.5, 0.5, 0.5], [1, 1, 1]]
			]
			]

		return jp_xyz, subchain_se_xyz

	#########################################################
	def set_atom(self, jp_xyz, subchain_se_xyz):
		calcd_data_dic={}
		count = 0
		for i in (range(self.multi)):
			for mol, jp in enumerate(jp_xyz):
				atom_all = []
				pos_all = {}
				# システム全体にわたるジャンクションポイントのxyzとIDの辞書を作成
				jp_id_dic, jp_xyz_dic, atom_jp = self.set_jp_id(jp, mol)
				atom_all.extend(atom_jp)
				pos_all.update(jp_xyz_dic)
				# print(jp_xyz_dic)
				# サブチェイン中の各アトムのxyzリストとボンドリストを作成
				strand_xyz, bond_all, atom_sc, angle_all = self.set_subchains(jp_id_dic, subchain_se_xyz[mol], mol)
				#
				atom_all.extend(atom_sc)
				pos_all.update(strand_xyz)
				#
				calcd_data_dic[count] = {"atom_all":atom_all, "bond_all":bond_all, "pos_all":pos_all, "angle_all":angle_all}
				count += 1
		return calcd_data_dic

	###################################################
	# システム全体にわたるJPのxyzとIDの辞書を作成
	def set_jp_id(self, jp_xyz, mol):
		jp_id_dic = {}
		jp_xyz_dic = {}
		atom_jp = []
		jp_id = 0
		for z in range(self.n_cell):
			for y in range(self.n_cell):
				for x in range(self.n_cell):
					base_xyz = np.array([x,y,z])
					for jp in jp_xyz:
						jp_id_dic[tuple(np.array(jp) + base_xyz)] = (jp_id)
						jp_xyz_dic[(jp_id)] = tuple(np.array(jp) + base_xyz)
						atom_jp.append([jp_id, 2*mol + 0, 0])
						jp_id += 1
		return jp_id_dic, jp_xyz_dic, atom_jp

		
	#########################################################
	# サブチェイン中の各アトムのxyzリストとボンドリストを作成
	def set_subchains(self, jp_id_dic, subchain_se_xyz, mol):
		strand_xyz = {}
		bond_all = {}
		atom_sc = []
		angle_all = []
		sub_id = len(jp_id_dic)
		bond_id = 0
		for z in range(self.n_cell):
			for y in range(self.n_cell):
				for x in range(self.n_cell):
					b_xyz = (x,y,z)
					for se_xyz in subchain_se_xyz:
						tmp_xyz, tmp_bond, new_sub_id, new_bond_id, tmp_atom_sc, tmp_angle = self.calc_single_subchain(jp_id_dic, sub_id, bond_id, b_xyz, se_xyz, mol)
						strand_xyz.update(tmp_xyz)
						bond_all.update(tmp_bond)
						atom_sc.extend(tmp_atom_sc)
						angle_all.append(tmp_angle)
						sub_id = new_sub_id
						bond_id = new_bond_id
		return strand_xyz, bond_all, atom_sc, angle_all

	###############################################################
	# 一本のサブチェイン中の各アトムのxyzリストとボンドリストを作成
	def calc_single_subchain(self, jp_id_dic, sub_id, bond_id, b_xyz, se_xyz, mol):
		tmp_xyz = {}
		tmp_bond = {}
		tmp_angle = []
		tmp_atom_sc = []
		bas_xyz = np.array(b_xyz)
		# サブチェインの末端間のベクトルを設定
		start_xyz = np.array(se_xyz[0]) + bas_xyz
		end_xyz = np.array(se_xyz[1]) + bas_xyz
		vec = end_xyz - start_xyz
		# 始点のアトムのIDを設定
		mod_xyz = list(start_xyz)[:]
		for dim in range(3):
			if mod_xyz[dim] == self.n_cell:
				mod_xyz[dim] = 0
		s_id = jp_id_dic[tuple(mod_xyz)]
		tmp_angle.append(s_id)
		# 終点のアトムのIDを周期境界条件で変更
		mod_xyz = list(end_xyz)[:]
		for dim in range(3):
			if mod_xyz[dim] == self.n_cell:
				mod_xyz[dim] = 0
		E_id = jp_id_dic[tuple(mod_xyz)]
		# サブチェインの鎖長分のループ処理
		for seg in range(self.n_segments):
			tmp_xyz[sub_id] = tuple(start_xyz + vec*(seg+1)/(self.n_segments+1.))
			if seg == 0 or seg == self.n_segments - 1:
				tmp_atom_sc.append([sub_id, 1, 1])
			else:
				tmp_atom_sc.append([sub_id, 2, 2])
			e_id = sub_id
			#
			if seg == 0:
				bond = 0
			else:
				bond = 1
			tmp_bond[bond_id] = tuple([bond, [s_id, e_id]])
			bond_id += 1
			tmp_angle.append(e_id)
			s_id = e_id
			sub_id += 1
			
			# if self.n_sc != 0:
			# 	sc_s_id = s_id
			# 	for i in range(self.n_sc):
			# 		tmp_xyz[sub_id] = tuple(np.array(pos)  +  (i + 1)*mod_o_vec*unit_len)
			# 		tmp_atom_st.append([sub_id, 2, 1])
			# 		sc_e_id = sub_id
			# 		#
			# 		bond = 2
			# 		tmp_bond[bond_id] = tuple([bond, [sc_s_id, sc_e_id]])
			# 		sc_s_id = sc_e_id
			# 		seq_atom_id += 1
			# 		bond_id += 1
			#
		e_id = E_id
		bond = 0
		tmp_bond[bond_id] = tuple([bond, [s_id, e_id]])
		tmp_angle.append(e_id)
		bond_id += 1
		return tmp_xyz, tmp_bond, sub_id, bond_id, tmp_atom_sc, tmp_angle

	######
	def find_ortho_vec(self, list):
		vec = np.array(list).reshape(-1,1)
		# 線形独立である新たな三次元ベクトルを見つける。
		rank = 0
		while rank != 2:
			a = np.array(np.random.rand(3)).reshape(-1,1)
			target = np.hstack((vec, a))
			rank = np.linalg.matrix_rank(target)
		# QR分解により
		q, r = np.linalg.qr( target )
		# print(q[:,1])
		ortho_vec = q[:,1]
		return ortho_vec

