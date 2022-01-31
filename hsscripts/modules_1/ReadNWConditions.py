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
	# check 'calc_condition.udf' and make it.
	findudf()
	# Read udf and setup initial conditions.
	basic_cond, nw_cond, sim_cond, target_cond, condition_text = read_and_setcondition()

	return basic_cond, nw_cond, sim_cond, target_cond, condition_text
	

###########################################
# check 'calc_condition.udf' and make it.
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
# make new udf when not found.
def makenewudf():
	contents = '''
	\\begin{def}
		CalcCond:{
			Cognac_ver:select{"cognac112"} "使用する Cognac のバージョン",
			Cores: int "計算に使用するコア数を指定"
			} "計算の条件を設定"
		SimulationCond:{
			Model:{TargetModel:select{"Regular_NW", "Random_NW"} "ネットワークのモデルを選択",
				Regular_NW:{chains:select{"3_Chain_S", "3_Chain_D", "4_Chain", "6_Chain", "8_Chain"} "分岐の数と種類を選択"
					} "規則構造での条件を入力",
				Random_NW:{chains:select{"3_Chain", "4_Chain", "5_Chain", "6_Chain", "7_Chain"} "分岐の数と種類を選択",
					Calc_Topolpgy:select{"Calc", "Read"} "ランダムネットワークの「計算を行うか、読み込むか」を選択",
						Calc:{pre_sampling:int "プレサンプリング数", pre_try:int "プレサンプリング時の再トライ数", sampling:int "サンプリング数", try:int "サンプリング時の再トライ数", n_parallel:int "並行計算のCPU数"} "ランダムサーチ計算する場合の条件を設定",
						Read:{dir_name:string} "過去の計算結果のディレクトリを記入"
					} "ランダム構造での条件を入力"
				} "シミュレーションの条件を設定"
			Type:{
				SimType:select{
					"Entangled",
					"NPT",
					"Multi",   
					"Multi_entangled",
					"Gel",    
					"Gel_entangled",
					"Gel_concd",
					"Gel_concd_entangled"
					} "アンサンブル等のシミュレーションの条件を選択",
					Entangled:{N_Segments: int "ストランド中のセグメント数", N_Subchain: int "各セグメントの側鎖の数", N_UnitCells: int "一辺あたりのユニットセルの数"} "密度、末端間距離を設定値に合わせるように多重度を自動設定。\\n絡み合いが入るように初期化",
					NPT:{N_Segments: int "ストランド中のセグメント数", N_Subchain: int "各セグメントの側鎖の数", N_UnitCells: int "一辺あたりのユニットセルの数", ExpansionRatio: float "NPT 計算での初期膨張率", StepPress[]: float "NPT 計算での圧力変化"} "密度、末端間距離を設定値に合わせるように多重度を自動設定。\\n絡み合いが入らないようにNPTで縮める。",
					Multi:{N_Segments: int "ストランド中のセグメント数", N_Subchain: int "各セグメントの側鎖の数", N_UnitCells: int "一辺あたりのユニットセルの数", Multiplicities: int "任意の多重度に設定可能", ExpansionRatio: float "NPT 計算での初期膨張率", StepPress[]: float "NPT 計算での圧力変化"} "設定した多重度で、密度を設定値になるように、システムサイズを縮める",
					Multi_entangled:{N_Segments: int "ストランド中のセグメント数", N_Subchain: int "各セグメントの側鎖の数", N_UnitCells: int "一辺あたりのユニットセルの数", Multiplicities: int "任意の多重度に設定可能"} "設定した多重度で、密度を設定値になるようにシステムサイズを縮め、絡み合いが入るように初期化",
					Gel:{N_Segments: int "ストランド中のセグメント数", N_Subchain: int "各セグメントの側鎖の数", N_UnitCells: int "一辺あたりのユニットセルの数", Multiplicities: int "任意の多重度に設定可能", ExpansionRatio: float "NPT 計算での初期膨張率", StepPress[]: float "NPT 計算での圧力変化"} "設定した多重度で、密度、末端間距離を設定値に合わせるように溶剤を添加",
					Gel_entangled:{N_Segments: int "ストランド中のセグメント数", N_Subchain: int "各セグメントの側鎖の数", N_UnitCells: int "一辺あたりのユニットセルの数", Multiplicities: int "任意の多重度に設定可能"} "設定した多重度で、密度、末端間距離を設定値に合わせるように溶剤を添加、絡み合いが入るように初期化",
					Gel_concd:{N_Segments: int "ストランド中のセグメント数", N_Subchain: int "各セグメントの側鎖の数", N_UnitCells: int "一辺あたりのユニットセルの数", Multiplicities: int "任意の多重度に設定可能", NV: float "Non Volatile: 固形分量 (0 to 1.0)", ExpansionRatio: float "NPT 計算での初期膨張率", StepPress[]: float "NPT 計算での圧力変化"} "設定した多重度で、溶剤量変化",
					Gel_concd_entangled:{N_Segments: int "ストランド中のセグメント数", N_Subchain: int "各セグメントの側鎖の数", N_UnitCells: int "一辺あたりのユニットセルの数", Multiplicities: int "任意の多重度に設定可能", NV: float "Non Volatile: 固形分量 (0 to 1.0)"} "設定した多重度で、溶剤量変化、絡み合いが入るように初期化"
				} "アンサンブル等のシミュレーションの条件を選択"
			target_density:float "平衡化シミュレーション実行時に使用する密度"
			} "シミュレーションの条件を設定"
	\end{def}

	\\begin{data}
		CalcCond:{"cognac112", 1}

		SimulationCond:{
			{"Regular_NW", {"4_Chain"}, {"4_Chain", "Calc", {1000, 100, 1000, 100, 1}, {""}}}
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
# Read udf and setup initial conditions
def read_and_setcondition():
	dic={'y':True,'yes':True,'q':False,'quit':False}
	while True:
		# read udf
		basic_cond, nw_cond, sim_cond = readconditionudf()
		# select
		condsetup = CondSetup(nw_cond, sim_cond)
		target_cond, condition_text = condsetup.calc_conditions()
		print('Change UDF: type [r]eload')
		print('Quit input process: type [q]uit')
		inp = input('Condition is OK ==> [y]es >> ').lower()
		if inp in dic:
			inp = dic[inp]
			break
		print('##### \nRead Condition UDF again \n#####\n\n')
	if inp:
		print("\n\nSetting UP progress !!")
		return basic_cond, nw_cond, sim_cond, target_cond, condition_text
	else:
		sys.exit("##### \nQuit !!")

####################################
# Read condition udf
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
	#####################################
	
	#########################################
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

	###################
	## Networkモデルの設定
	nw_model = u.get('SimulationCond.Model.TargetModel')
	###################
	## Networkモデルの設定
	if nw_model == "Regular_NW":
		strand_type = u.get('SimulationCond.Model.Regular_NW.chains')
		restart = ''
		cond_top = []
	elif nw_model == "Random_NW":
		strand_type = u.get('SimulationCond.Model.Random_NW.chains')
	################
	if strand_type == "3_Chain" or strand_type == "3_Chain_S" or strand_type == "3_Chain_D":
		n_strand = 3
	elif strand_type == "4_Chain":
		n_strand = 4
	elif strand_type == "5_Chain":
		n_strand = 5
	elif strand_type == "6_Chain":
		n_strand = 6
	elif strand_type == "7_Chain":
		n_strand = 7
	elif strand_type == "8_Chain":
		n_strand = 8
	##############################################
	if nw_model == "Random_NW":
		calc = u.get('SimulationCond.Model.Random_NW.Calc_Topolpgy')
		restart = ''
		cond_top = []
		if calc == 'Read':
			restart = u.get('SimulationCond.Model.Random_NW.Read.dir_name')
			print(restart)
			cond_top = []
			# if os.path.exists(restart):
			# 	print('Read existing topology file !')
			# else:
			# 	sys.exit('Cannnot find topology file to read !!')
			
			if not os.path.exists(os.path.join(restart, 'init.pickle')):
				exit("##########\ntarget directory does not exists.")
			elif n_strand != int(restart.split('_')[0]):
				sys.exit("##########\nnumber of strands: selected n_strand is different from original Calculation.")
			elif n_cell != int(restart.split('_')[2]):
				sys.exit("##########\nnumber of cells: selected n_cell is different from original Calculation.")
			

		elif calc == 'Calc':
			cond_top = u.get('SimulationCond.Model.Random_NW.Calc')
	#
	nw_cond = [nw_model, strand_type, n_strand, n_segments, n_cell, n_sc]
	sim_cond = [sim_type, multi_init, target_density, expand, l_bond, c_n, step_press, nv, restart, cond_top]
	print(sim_cond)
	return basic_cond, nw_cond, sim_cond

######################################
##### Setup Calculation COnditions
######################################
class CondSetup:
	def __init__(self, nw_cond, sim_cond):
		self.nw_model = nw_cond[0]
		self.strand = nw_cond[1]
		self.n_strand = nw_cond[2]
		self.n_segments = nw_cond[3]
		self.n_cell = nw_cond[4]
		self.n_sc = nw_cond[5]
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
		
	############################################
	##### ネットワークポリマーの諸量を計算 ######
	############################################
	#-----ネットワークポリマーの諸量を計算
	def calc_conditions(self):
		## 計算システムの諸量を計算して、出力
		e2e, n_chains, n_beads_unit, org_unitcell = self.set_length()
		target_cond, condition_text = self.init_calc(e2e, n_chains, n_beads_unit, org_unitcell)
		return target_cond, condition_text

	#####################
	#
	def set_length(self):
		e2e = self.l_bond*((self.n_segments + 1)*self.c_n)**0.5					# 理想鎖状態での末端間距離

		if self.nw_model == "Regular_NW":
			if self.strand == "3_Chain_S":
				n_chains = 12						        					# サブチェインの本数
				n_beads_unit = 8 + self.n_segments*(1 + self.n_sc)*n_chains		# ユニットセル当たりの粒子数
				org_unitcell = (2*2**0.5)*e2e				        			# 理想鎖状態でのユニットセル長
			elif self.strand == "3_Chain_D":
				n_chains = 24						       
				n_beads_unit = 16 + self.n_segments*(1 + self.n_sc)*n_chains	
				org_unitcell = (2*2**0.5)*e2e		
			elif self.strand == "4_Chain":
				n_chains = 16						      
				n_beads_unit = 8 + self.n_segments*(1 + self.n_sc)*n_chains	
				org_unitcell = (4*3**0.5)*e2e/3			 
			elif self.strand == "6_Chain":
				n_chains = 3						  
				n_beads_unit = 1 + self.n_segments*(1 + self.n_sc)*n_chains		
				org_unitcell = e2e						      
			elif self.strand == "8_Chain":
				n_chains = 8						   
				n_beads_unit = 2 + self.n_segments*(1 + self.n_sc)*n_chains   
				org_unitcell = (2*3**0.5)*e2e/3	

		elif self.nw_model == "Random_NW":
			n_chains = self.n_strand
			n_beads_unit = 2 + self.n_segments*(1 + self.n_sc)*n_chains
			org_unitcell = (2*3**0.5)*e2e/3	

		return e2e, n_chains, n_beads_unit, org_unitcell

	###############################################################
	def init_calc(self, e2e, n_chains, n_beads_unit, org_unitcell):
		n_solvent = 0
		flag = 0
		if self.sim_type == "Entangled" or self.sim_type == "NPT":
			org_system = org_unitcell*self.n_cell									# e2e から決めたシステムサイズ
			fin_multi = round(self.target_density*org_unitcell**3/n_beads_unit)		# 密度を設定値とした場合に必要な多重度
			fin_dens = n_beads_unit*fin_multi/org_unitcell**3						# 上記の多重度での密度
			err_dens = round((fin_dens/self.target_density - 1)*100, 2) 			# 設定密度との誤差(%)
			if abs(err_dens) > 1:
				flag = 1
			single_net_atom = int(n_beads_unit*self.n_cell**3.)	    				# 一つのネットワーク中の粒子数
			total_net_atom = int(fin_multi*single_net_atom)    						# 全システム中のネットワーク粒子数
			system = (total_net_atom/self.target_density)**(1/3)					# その際のシステムサイズ

			
			shrinkage = system/org_system											# 収縮比
			mod_e2e = shrinkage*e2e													# 収縮後の末端間距離
			unit_cell = shrinkage*org_unitcell

			nv = 1.0
			nu = n_chains*fin_multi/org_unitcell**3.								# ストランドの数密度

		elif self.sim_type == "Multi" or self.sim_type == "Multi_entangled":
			fin_multi = self.multi_init
			err_dens = 0.
			org_system = org_unitcell*self.n_cell							# e2e から決めたシステムサイズ
			single_net_atom = int(n_beads_unit*self.n_cell**3.)	    # 一つのネットワーク中の粒子数
			total_net_atom = int(fin_multi*single_net_atom)    # 全システム中のネットワーク粒子数
			vol = total_net_atom/self.target_density			# システム体積
			system = vol**(1./3.)							# 収縮後のシステムサイズ
			shrinkage = system/org_system						# 収縮比
			mod_e2e = shrinkage*e2e								# 収縮後の末端間距離
			unit_cell = shrinkage*org_unitcell
			fin_dens = total_net_atom/vol
			nv = 1.0
			nu = n_chains*fin_multi*self.n_cell**3./vol		# ストランドの数密度

		elif self.sim_type == "Gel" or self.sim_type == "Gel_entangled":
			unit_cell = org_unitcell
			org_system = org_unitcell*self.n_cell				# e2e から決めたシステムサイズ
			mod_e2e = e2e										# 収縮後の末端間距離
			fin_multi = self.multi_init
			fin_dens = self.target_density
			err_dens = 0.
			shrinkage = 1.0								
			system = org_unitcell*self.n_cell					# システムサイズ
			single_net_atom = int(n_beads_unit*self.n_cell**3.)	# 一つのネットワーク中の粒子数
			total_net_atom = int(fin_multi*single_net_atom)  	# 全システム中のネットワーク粒子数
			vol = system**3									    # システム体積
			total_atom = int(vol*self.target_density)		    # システム中の粒子総数
			n_solvent = int(total_atom - total_net_atom)		# 溶媒の粒子数
			nv = total_net_atom/total_atom						# ネットワークの体積分率
			nu = n_chains*fin_multi*self.n_cell**3./vol		# ストランドの数密度

		elif self.sim_type == "Gel_concd" or self.sim_type == "Gel_concd_entangled":
			unit_cell = org_unitcell
			org_system = org_unitcell*self.n_cell
			mod_e2e = e2e								# 収縮後の末端間距離
			fin_multi = self.multi_init
			fin_dens = self.target_density
			err_dens = 0.
			shrinkage = 1.0
			single_net_atom = int(n_beads_unit*self.n_cell**3.)	    # 一つのネットワーク中の粒子数
			total_net_atom = int(fin_multi*single_net_atom)    # 全システム中のネットワーク粒子数
			n_solvent = int(total_net_atom*(1.0 - self.nv)/self.nv)		# 溶媒の粒子数
			total_atom = int(total_net_atom + n_solvent)		    # システム中の粒子総数
			vol = total_atom/self.target_density			# システム体積
			system = vol**(1./3.)							# システムサイズ
			total_atom = int(vol*self.target_density)		    # システム中の粒子総数
			nv = self.nv						# ネットワークの体積分率
			nu = n_chains*fin_multi*self.n_cell**3./vol		# ストランドの数密度
		else:
			sys.exit("Something Wrong!!")
		#
		text = "#########################################" + "\n"
		text += "ネットワークトポロジー\t\t" + str(self.nw_model) + "\n"
		text += "ネットワークモデル\t\t" + str(self.strand) + "\n"
		text += "#########################################" + "\n"
		text += "ストランド中のセグメント数\t" + str(self.n_segments) + "\n"
		text += "特性比:\t\t\t\t" + str(round(self.c_n, 2)) + "\n"
		text += "初期の末端間距離:\t\t" + str(round(e2e, 4)) + "\n"
		text += "#########################################" + "\n"
		text += "当初の単位ユニット:\t\t" + str(round(org_unitcell, 4)) + "\n"
		text += "一辺当たりの単位ユニット数\t" + str(self.n_cell) + "\n"
		text += "当初のシステムサイズ:\t\t" + str(round(org_system, 4)) + "\n"
		text += "#########################################" + "\n"
		text += "ネットワークタイプ\t\t" + str(self.sim_type) + "\n"
		text += "#########################################" + "\n"
		text += "初期設定密度:\t\t\t" + str(self.target_density) + "\n"
		text += "多重度:\t\t\t\t" + str(fin_multi) + "\n"
		text += "全セグメント数:\t\t\t" + str(total_net_atom) + "\n"
		text += "収縮比:\t\t\t\t" + str(round(shrinkage, 4)) + "\n"
		text += "任意の多重度での密度:\t\t" + str(round(fin_dens, 4)) + "\n"
		text += "密度の誤差:\t\t\t" + str(round(err_dens, 3)) + "%" + "\n"
		text += "ステップ圧力:\t\t" + ', '.join(map(str, self.step_press)) + "\n"
		text += "収縮後の末端間距離:\t\t" + str(round(mod_e2e, 4)) + "\n"
		text += "システムサイズ:\t\t\t" + str(round(system, 4)) + "\n"
		text += "#########################################" + "\n"
		text += "溶媒のセグメント数:\t\t" + str(n_solvent) + "\n"
		text += "ネットワークの体積分率:\t\t" + str(round(nv, 4)) + "\n"
		text += "#########################################" + "\n"
		text += "ストランドの数密度:\t\t" + str(round(nu, 5)) + "\n"
		text += "#########################################" + "\n"
		print(text)

		if (self.sim_type == "Entangled" or self.sim_type == "NPT") and flag == 1:
			print(u"\n##### \n圧縮後の密度が、設定したい密度と 1% 以上違います。\nそれでも計算しますか？\n#####\n")
		#
		with open("calc_conditions.txt", 'w') as f:
			f.write(text)
		#
		target_cond = [system, unit_cell, total_net_atom, nu, self.nw_model, fin_multi, n_solvent]

		return target_cond, text
	
