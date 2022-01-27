#!/usr/bin/env python
# -*- coding: utf-8 -*-
#######################################################
import os
import sys
import codecs
from UDFManager import UDFManager
#######################################################
#
def setupcondition():
	# 
	findudf()

	dic={'y':True,'yes':True,'q':False,'quit':False}
	while True:
		# read udf
		basic_cond, nw_cond, sim_cond = readconditionudf()
		# select
		condsetup = CondSetup(nw_cond, sim_cond)
		target_cond, condition_text = condsetup.init_calc()
		print('Change UDF: type [r]eload')
		print('Quit input process: type [q]uit')
		inp = input('Condition is OK ==> [y]es >> ').lower()
		if inp in dic:
			inp = dic[inp]
			break
		print('##### \nRead Condition UDF again \n#####\n\n')
	if inp:
		print("\n\nSetting UP progress !!")
	else:
		sys.exit("##### \nQuit !!")
	
	return basic_cond, nw_cond, sim_cond, target_cond, condition_text

###########################################
# 
def findudf():
	if not os.path.isfile('./calc_condition.udf'):
		print()
		print('In this directory, no "calc_condition.udf" is found !')
		print('New one will be generated.')
		print('Please, modify and save it !\n')
		makenewudf()
		input('Press ENTER to continue...')
	return
###########################################
# 
def makenewudf():
	contents = '''
	\\begin{def}
		CalcCond:{
			Cognac_ver:select{"cognac112"},
			Cores: int // 計算に使用するコア数を指定
			}
		SimulationCond:{
			Model:{TargetModel:select{"Regular_NW", "Random_NW"},
				Regular_NW:{chains:select{"3_Chain_S", "3_Chain_D", "4_Chain", "6_Chain", "8_Chain"}
					},
				Random_NW:{chains:select{"3_Chain", "4_Chain", "5_Chain", "6_Chain", "7_Chain"},
					Calc_Topolpgy:select{"Calc", "Read"},
						Calc:{pre_sampling:int, pre_try:int, sampling:int, try:int},
						Read:{file_name:string}
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
					Entangled:{N_Segments: int, N_Subchain: int, N_UnitCells: int},
					NPT:{N_Segments: int, N_Subchain: int, N_UnitCells: int, ExpansionRatio: float, StepPress[]: float},
					Multi:{N_Segments: int, N_Subchain: int, N_UnitCells: int, Multiplicities: int, ExpansionRatio: float, StepPress[]: float},
					Multi_entangled:{N_Segments: int, N_Subchain: int, N_UnitCells: int, Multiplicities: int},
					Gel:{N_Segments: int, N_Subchain: int, N_UnitCells: int, Multiplicities: int, ExpansionRatio: float, StepPress[]: float},
					Gel_entangled:{N_Segments: int, N_Subchain: int, N_UnitCells: int, Multiplicities: int},
					Gel_concd:{N_Segments: int, N_Subchain: int, N_UnitCells: int, Multiplicities: int, NV: float, ExpansionRatio: float, StepPress[]: float},
					Gel_concd_entangled:{N_Segments: int, N_Subchain: int, N_UnitCells: int, Multiplicities: int,NV: float}
				}
			target_density:float
			}
	\end{def}

	\\begin{data}
		CalcCond:{"cognac112", 1}

		SimulationCond:{
			{"Regular_NW", {"4_Chain"}, {"4_Chain", "Calc", {1000, 100, 1000, 1000}, {""}}}
			{"NPT",
				{20, 0, 3},
				{20, 0, 3, 2.0, [0.2, 0.5, 1.0, 2.0, 3.0, 4.5]},
				{20, 0, 3, 1, 2.0, [0.2, 0.5, 1.0, 2.0, 5.0, 6.5, 7.0]},
				{20, 0, 3, 1},
				{20, 0, 3, 1, 2.0, [0.2, 0.5, 1.0, 2.0, 5.0, 6.5, 7.0]},
				{20, 0, 3, 2},
				{20, 0, 3, 2, 0.5, 2.0, [0.2, 0.5, 1.0, 2.0, 5.0, 6.5, 7.0]},
				{20, 0, 3, 2, 0.5}
			}
		0.85
		}
	\end{data}
	'''
	###
	with codecs.open('./calc_condition.udf', 'w', 'utf_8') as f:
		f.write(contents)
	return

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
				print('Read existing topology file !')
			else:
				sys.exit('Cannnot find topology file to read !!')
		elif calc == 'Calc':
			cond_top = u.get('SimulationCond.Model.Random_NW.Calc')
	###################
	## ポリマー鎖の設定
	sim_type = u.get('SimulationCond.Type.SimType')
	#################################
	if sim_type == "Entangled":
		n_segments = u.get('SimulationCond.Type.Entangled.N_Segments')
		n_sc = u.get('SimulationCond.Type.Entangled.N_Subchain')
		n_cell = u.get('SimulationCond.Type.Entangled.N_UnitCells')
		multi_init = 0
		nv = 1.0
		expand = 1.0
		step_press = []
	elif sim_type == "NPT":
		n_segments = u.get('SimulationCond.Type.NPT.N_Segments')
		n_sc = u.get('SimulationCond.Type.NPT.N_Subchain')
		n_cell = u.get('SimulationCond.Type.NPT.N_UnitCells')
		multi_init = 0
		nv = 1.0
		expand = u.get('SimulationCond.Type.NPT.ExpansionRatio')
		step_press = u.get('SimulationCond.Type.NPT.StepPress[]')
	elif sim_type == "Multi":
		n_segments = u.get('SimulationCond.Type.Multi.N_Segments')
		n_sc = u.get('SimulationCond.Type.Multi.N_Subchain')
		n_cell = u.get('SimulationCond.Type.Multi.N_UnitCells')
		multi_init = u.get('SimulationCond.Type.Multi.Multiplicities')
		nv = 1.0
		expand = u.get('SimulationCond.Type.Multi.ExpansionRatio')
		step_press = u.get('SimulationCond.Type.Multi.StepPress[]')
	elif sim_type == "Multi_entangled":
		n_segments = u.get('SimulationCond.Type.Multi_entangled.N_Segments')
		n_sc = u.get('SimulationCond.Type.Multi_entangled.N_Subchain')
		n_cell = u.get('SimulationCond.Type.Multi_entangled.N_UnitCells')
		multi_init = u.get('SimulationCond.Type.Multi_entangled.Multiplicities')
		nv = 1.0
		expand = 1.0
		step_press = []
	elif sim_type == "Gel":
		n_segments = u.get('SimulationCond.Type.Gel.N_Segments')
		n_sc = u.get('SimulationCond.Type.Gel.N_Subchain')
		n_cell = u.get('SimulationCond.Type.Gel.N_UnitCells')
		multi_init = u.get('SimulationCond.Type.Gel.Multiplicities')
		expand = u.get('SimulationCond.Type.Gel.ExpansionRatio')
		step_press = u.get('SimulationCond.Type.Gel.StepPress[]')
		nv = 1.0
	elif sim_type == "Gel_entangled":
		n_segments = u.get('SimulationCond.Type.Gel_entangled.N_Segments')
		n_sc = u.get('SimulationCond.Type.Gel_entangled.N_Subchain')
		n_cell = u.get('SimulationCond.Type.Gel_entangled.N_UnitCells')
		multi_init = u.get('SimulationCond.Type.Gel_entangled.Multiplicities')
		nv = 1.0
		expand = 1.0
		step_press = []
	elif sim_type == "Gel_concd":
		n_segments = u.get('SimulationCond.Type.Gel_concd.N_Segments')
		n_sc = u.get('SimulationCond.Type.Gel_concd.N_Subchain')
		n_cell = u.get('SimulationCond.Type.Gel_concd.N_UnitCells')
		multi_init = u.get('SimulationCond.Type.Gel_concd.Multiplicities')
		nv = u.get('SimulationCond.Type.Gel_concd.NV')
		expand = u.get('SimulationCond.Type.Gel_concd.ExpansionRatio')
		step_press = u.get('SimulationCond.Type.Gel_concd.StepPress[]')
	elif sim_type == "Gel_concd_entangled":
		n_segments = u.get('SimulationCond.Type.Gel_concd_entangled.N_Segments')
		n_sc = u.get('SimulationCond.Type.Gel_concd_entangled.N_Subchain')
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
	l_bond = 0.97
	#
	if n_segments <= 10:
		c_n = 1.5
	elif n_segments <= 20:
		c_n = 1.65
	elif n_segments <= 40:
		c_n = 1.7
	else:
		c_n = 1.75
	#
	nw_cond = [nw_model, strand, n_segments, n_cell, n_sc]
	sim_cond = [sim_type, multi_init, target_density, expand, l_bond, c_n, step_press, nv]

	return basic_cond, nw_cond, sim_cond

######################################
##### Setup Calculation COnditions
######################################
class CondSetup:
	def __init__(self, nw_cond, sim_cond):
		self.nw_model = nw_cond[0]
		self.strand = nw_cond[1]
		self.n_segments = nw_cond[2]
		self.n_cell = nw_cond[3]
		self.n_sc = nw_cond[4]
		#
		self.sim_type = sim_cond[0]
		self.multi_init = sim_cond[1]
		self.target_density = sim_cond[2]
		self.expand = sim_cond[3]
		self.l_bond = sim_cond[4]
		self.c_n = sim_cond[5]
		self.step_press = sim_cond[6]
		self.nv = sim_cond[7]
		#
		if self.strand == "3_Chain_S" or self.strand == "3_Chain_D":
			self.n_strand = 3
		elif self.strand == "4_Chain":
			self.n_strand = 4
		elif self.strand == "5_Chain":
			self.n_strand = 5
		elif self.strand == "6_Chain":
			self.n_strand = 6
		elif self.strand == "7_Chain":
			self.n_strand = 7
	##########################################
	##### ネットワークポリマーの諸量を計算 ######
	##########################################
	#-----ネットワークポリマーの諸量を計算
	def calc_conditions(self):
		## 計算システムの諸量を計算して、出力
		target_cond, condition_text = self.init_calc()
		return target_cond, condition_text

	################################################################################
	def init_calc(self):
		structure = "Regular_NW"
		e2e = self.l_bond*((self.n_segments + 1)*self.c_n)**0.5	# 理想鎖状態での末端間距離
		#
		if self.strand == "3_Chain_S":
			n_chains = 12						        # サブチェインの本数
			n_beads_unit = 8 + self.n_segments*(1 + self.n_sc)*n_chains		# ユニットセル当たりの粒子数
			a_cell = (2*2**0.5)*e2e				        # 理想鎖状態でのユニットセル長
		elif self.strand == "3_Chain_D":
			n_chains = 24						        # サブチェインの本数
			n_beads_unit = 16 + self.n_segments*(1 + self.n_sc)*n_chains	# ユニットセル当たりの粒子数
			a_cell = (2*2**0.5)*e2e				        # 理想鎖状態でのユニットセル長
		elif self.strand == "4_Chain":
			n_chains = 16						        # サブチェインの本数
			n_beads_unit = 8 + self.n_segments*(1 + self.n_sc)*n_chains		# ユニットセル当たりの粒子数
			a_cell = (4*3**0.5)*e2e/3			        # 理想鎖状態でのユニットセル長
		elif self.strand == "6_Chain":
			n_chains = 3						        # サブチェインの本数
			n_beads_unit = 1 + self.n_segments*(1 + self.n_sc)*n_chains		# ユニットセル当たりの粒子数
			a_cell = e2e						            # 理想鎖状態でのユニットセル長
		elif self.strand == "8_Chain":
			n_chains = 8						        # サブチェインの本数
			n_beads_unit = 2 + self.n_segments*(1 + self.n_sc)*n_chains     # ユニットセル当たりの粒子数
			a_cell = (2*3**0.5)*e2e/3					# 理想鎖状態でのユニットセル長
		#
		n_solvent = 0
		#
		if self.sim_type == "Entangled" or self.sim_type == "NPT":
			unit_cell = a_cell
			self.multi = round(self.target_density*a_cell**3/n_beads_unit)		# 密度を設定値とした場合に必要な多重度
			fin_dens = n_beads_unit*self.multi/a_cell**3						# 上記の多重度での密度
			err_dens = round((fin_dens/self.target_density - 1)*100, 2) 	# 設定密度との誤差
			single_net_atom = int(n_beads_unit*self.n_cell**3.)	    			# 一つのネットワーク中の粒子数
			total_net_atom = int(self.multi*single_net_atom)    			# 全システム中のネットワーク粒子数
			system = (total_net_atom/self.target_density)**(1/3)			# その際のシステムサイズ
			nu = n_chains*self.multi/a_cell**3.								# ストランドの数密度
		elif self.sim_type == "Multi" or self.sim_type == "Multi_entangled":
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
		elif self.sim_type == "Gel" or self.sim_type == "Gel_entangled":
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
		elif self.sim_type == "Gel_concd" or self.sim_type == "Gel_concd_entangled":
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
		text += "ネットワークタイプ\t\t" + str(self.sim_type) + "\n"
		text += "ストランド中のセグメント数\t" + str(self.n_segments) + "\n"
		text += "特性比:\t\t\t" + str(round(self.c_n,2)) + "\n"
		text += "末端間距離:\t\t\t" + str(round(e2e,5)) + "\n"
		text += "一辺当たりの単位ユニット数\t" + str(self.n_cell) + "\n"
		text += "#########################################" + "\n"

		if self.sim_type == "Entangled" or self.sim_type == "NPT":
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
		elif self.sim_type == "Multi":
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
		elif self.sim_type == "Multi_entangled":
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
		elif self.sim_type == "Gel":
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
		elif self.sim_type == "Gel_entangled":
			text += "設定密度:\t\t\t" + str(self.target_density) + "\n"
			text += "多重度:\t\t\t\t" + str(self.multi) + "\n"
			text += "単一ネットワークのセグメント数:\t" + str(single_net_atom) + "\n"
			text += "溶媒のセグメント数:\t\t" + str(n_solvent) + "\n"
			text += "システムサイズ:\t\t\t" + str(round(system, 4)) + "\n"
			text += "ネットワークの体積分率:\t\t" + str(round(nv, 4)) + "\n"
			text += "ストランドの数密度:\t\t" + str(round(nu, 5)) + "\n"
			text += "#########################################" + "\n"
			print(text)
		elif self.sim_type == "Gel_concd":
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
		elif self.sim_type == "Gel_concd_entangled":
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