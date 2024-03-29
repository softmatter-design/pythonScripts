#!/usr/bin/env python
# -*- coding: utf-8 -*-
################################################################
import modules_1
import sys
################################################################
## 計算条件
##################
# 使用するCognacのバージョンを入れてください。
ver_cognac = "cognac101"
blank_udf = 'cognac101.udf'
# ベースとするUDFの名前
base_udf = "base_uin.udf"
# 計算に使用するコア数
core = 6
###########
files_cond = [ver_cognac, blank_udf, base_udf, core]
################################################################
# 基本ネットワーク関連の設定
################################################################
# 設定したい分岐数
n_strand = 4
#################
# ストランドの長さ
n_segments = 50
##########################
# 一辺当たりの単位ユニット数
n_cell = 3
#################
# ストランドの側鎖
n_sc = 0
#####################
# ネットワークの多重度
multi_init = 1
###################################
# NPT 計算での初期膨潤度
expand = 3

################################################################
# トポロジー変換に関する設定
################################################################
# リスタートの有無
restart = 0
############################################################
# プレ探索の条件：繰り返しサンプリング数、最小構造探索の繰り返し数
pre_sampling = 1000
pre_try = 100
# 本探索条件
n_sampling = 100
n_try = 100

###############
cond_top = [pre_try, pre_sampling, n_try, n_sampling]
################################################################
# ヒストグラム関連の設定
################################################################
# ヒストグラムの分割数
hist_bins = 50
################################################################
# ネットワークのタイプ
################################################################
# nw_type = "Sunuke"		# 密度、末端間距離を設定値に合わせるように多重度を変化して、スヌケ鎖を設定
nw_type = "NPT"		# 密度、末端間距離を設定値に合わせるように多重度を変化。 絡み合いが入らないように初期化
# nw_type = "KG_multi"		# 密度、末端間距離を設定値に合わせるように多重度を変化。 絡み合いが入るように初期化
# nw_type = "KG_gel"		# 設定した多重度で、密度、末端間距離を設定値に合わせるように、溶剤を添加
##############################
if nw_type == "Sunuke":
	l_bond = 0.97
	c_n = 1.26
	target_density = 0.85
elif nw_type == "KG_gel":
	l_bond = 0.97
	c_n = 1.26
	target_density = 0.85
elif nw_type == "KG_multi":
	l_bond = 0.97
	c_n = 1.65
	target_density = 0.85
elif nw_type == "NPT":
	l_bond = 0.97
	c_n = 1.26
	target_density = 0.85
#
nw_model = ''
sim_cond = [nw_model, nw_type, n_segments, n_cell, multi_init, target_density, n_strand, l_bond, c_n, expand]

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

################################################################
def main():
################################################################
	######################################
	# 初期状態を設定し、基本となるデータを取得
	args = sys.argv
	init = modules_1.RandomNW_mod.InitialSetup(sim_cond, restart, args)
	read_file_path, target_name, target_cond = init.get_base_data()
	###################################################################################################
	make8 = modules_1.RandomNW_mod_2.Make8(sim_cond)
	base_top_list = make8.make8()









	n_multi = 5

	if restart == 0:
		# トポロジーの異なるネットワークを探索して、任意の多重度のネットワークポリマーの代数的連結性の分布関数を策定
		mod = modules_1.RandomNW_mod_2.ModifyTop(base_top_list, sim_cond, cond_top, target_cond, hist_bins, read_file_path)
		candidate_list, target_dir = mod.top_search()
	else:
		sel = modules_1.RandomNW_mod_2.Select(read_file_path, hist_bins, n_multi)
		candidate_list, target_dir = sel.top_select()

	sel = modules_1.RandomNW_mod_2.Select(target_dir, hist_bins, n_multi)
	top_dic_list = sel.nw_search(candidate_list, target_dir)
	
	###########################################
	# ターゲットとなるネットワーク全体の辞書を設定。
	setup = modules_1.RandomNW_mod_2.SetUp(top_dic_list, base_top_list, n_segments, n_sc)
	calcd_data_dic = setup.make_data_dic()
















	##################
	# baseUDF の作成
	# baseudf = SetupInitUDF.MakeInitUDF(sim_cond, target_cond, files_cond, names, target_name, calcd_data_dic, expand)
	# target_dir = baseudf.setup_baseudf()
	
	# ##################
	# # Init_UDF の作成
	# setup = EquivCalcSetup.SetUpUDF(nw_type, files_cond, target_name, py_mod, target_dir, names)
	# setup.setup_udf()

################################################################################
#      Main     #
################################################################################
if __name__=='__main__':
	main()
