#!/usr/bin/env python
# -*- coding: utf-8 -*-
#######################################################
import os
from UDFManager import UDFManager
#######################################################

def setupcondition():
	if not os.path.isfile('./calc_condition.udf'):
		print('##############################################')
		print('"calc_condition.udf" が無いようなので、作ります。')
		print('適切に条件を設定してください。')
		print('##############################################\n\n')
		makenewudf()
	print('#############################################################')
	ans = input('この "calc_condition.udf" で計算条件は良いでしょうか？ (y/n)').lower()
	if ans in ['y', 'ye', 'yes']:
		basic_cond, sim_cond, names = readconditionudf()
	else:
		exit("終了します。")
	return basic_cond, sim_cond, names
	
def readconditionudf():
	u = UDFManager('calc_condition.udf')

	u.jump(-1)
	## 計算条件
	##################
	# 使用するCognacのバージョン
	ver_cognac = u.get('Basics.Cognac_ver')
	# 計算に使用するコア数
	core = u.get('Basics.Cores')
	# ベースとするUDFの名前
	base_udf = "base_uin.udf"
	blank_udf = ver_cognac + '.udf'
	###########
	basic_cond = [ver_cognac, blank_udf, base_udf, core]
	#######################################################
	## 計算ターゲット
	###################
	## Networkモデルの設定
	nw_model = u.get('SimulationCond.TargetModel')
	###################
	## ポリマー鎖の設定
	nw_type = u.get('SimulationCond.Type.NWType')
	#################################
	# ストランドの長さ
	n_segments = u.get('SimulationCond.Type.KG_NPT.N_Segments')
	# 一辺当たりの単位ユニット数
	n_cell = u.get('SimulationCond.Type.KG_NPT.N_UnitCells')
	# ネットワークの多重度
	if nw_type == "KG_single" or nw_type == "KG_gel":
		multi_init = u.get('SimulationCond.Type.KG_single.Multiplicities')
	else:
		multi_init = 0
	# 設定密度
	target_density = u.get('SimulationCond.target_density')
	###################################
	if nw_type == "KG_NPT":
		expand = u.get('SimulationCond.Type.KG_NPT.ExpansionRatio')
		step_press = u.get('SimulationCond.Type.KG_NPT.StepPress[]')
	elif nw_type == "KG_single":
		expand = u.get('SimulationCond.Type.KG_gel.ExpansionRatio')
		step_press = u.get('SimulationCond.Type.KG_gel.StepPress[]')
	elif nw_type == "KG_gel":
		expand = u.get('SimulationCond.Type.KG_single.ExpansionRatio')
		step_press = u.get('SimulationCond.Type.KG_single.StepPress[]')
	else:
		expand = 1.0
		step_press = []
	########################################################
	if nw_model == "3_Chain_S" or nw_model == "3_Chain_D":
		n_strand = 3
	elif nw_model == "4_Chain":
		n_strand = 4
	elif nw_model == "6_Chain":
		n_strand = 6
	elif nw_model == "8_Chain":
		n_strand = 8
	#
	if nw_type == "KG_entangled" or nw_type == "KG_NPT" or nw_type == "KG_single" or nw_type == "KG_gel":
		l_bond = 0.97
		if n_segments <= 10:
			c_n = 1.5
		elif n_segments <= 20:
			c_n = 1.65
		elif n_segments <= 40:
			c_n = 1.7
		else:
			c_n = 1.75
	elif nw_type == "Sunuke" or nw_type == "Sunuke_dens":
		l_bond = 0.97
		c_n = 1.0
	#
	sim_cond = [nw_model, nw_type, n_segments, n_cell, multi_init, target_density, expand, n_strand, l_bond, c_n, step_press]
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

def makenewudf():
	contents = '''
		\\begin{def}
			Basics:{
				Cognac_ver:select{
					"cognac101", 
					"cognac112"
					}, // Cognac version
				Cores: int // 計算に使用するコア数を指定
				}

			SimulationCond:{
				TargetModel:select{ 
					"3_Chain_S", 
					"3_Chain_D", 
					"4_Chain", 
					"6_Chain", 
					"8_Chain"
					},
				Type:{
					NWType:select{
						"KG_entangled",   // 密度、末端間距離を設定値に合わせるように多重度を変化。 絡み合いが入るように初期化
						"KG_NPT",		    // 密度、末端間距離を設定値に合わせるように多重度を変化。 絡み合いが入らないようにNPTで縮める
						"KG_single",    	// 設定した多重度で、密度を設定値になるように、システムサイズを縮める
						"KG_gel",	    	// 設定した多重度で、密度、末端間距離を設定値に合わせるように、溶剤を添加
						"Sunuke",	    	// 設定した多重度で末端間距離を設定値に合わせるように、密度を無視してスヌケ鎖を設定
						"Sunuke_dens"       // 設定した多重度、設定した密度となるようにシステムサイズを縮めてスヌケ鎖を設定
						},
						KG_entangled:{
							N_Segments: int,
							N_UnitCells: int,
						},
						KG_NPT:{
							N_Segments: int,
							N_UnitCells: int,
							ExpansionRatio: float
							StepPress[]: float
						},
						KG_single:{
							N_Segments: int,
							N_UnitCells: int,
							Multiplicities: int,
							ExpansionRatio: float
							StepPress[]: float
						},
						KG_gel:{
							N_Segments: int,
							N_UnitCells: int,
							Multiplicities: int,
							ExpansionRatio: float
							StepPress[]: float
						},
						Sunuke:{
							N_Segments: int,
							N_UnitCells: int,
							Multiplicities: int,
							ExpansionRatio: float
						},
						Sunuke_dens:{
							N_Segments: int,
							N_UnitCells: int,
							Multiplicities: int,
							ExpansionRatio: float
						}
					}
				target_density:float
				}
		\end{def}

		\\begin{data}
			Basics:{
				"cognac112",
				1
				}

			SimulationCond:{
				"4_Chain",
				{"KG_NPT",
					{20, 3},
					{20, 3, 2.0, [0.2, 0.5, 1.0, 2.0, 3.0, 4.5]},
					{20, 3, 1, 2.0, [0.2, 0.5, 1.0, 2.0, 5.0, 6.5, 7.0]},
					{20, 3, 1, 2.0, [0.2, 0.5, 1.0, 2.0, 5.0, 6.5, 7.0]},
					{20, 3, 3, 3},
					{20, 3, 3, 3}
				}
			0.85
			}
		\end{data}
		'''
	###
	with open('./calc_condition.udf', 'w') as f:
		f.write(contents)
	
	return