#!/usr/bin/env python
# -*- coding: utf-8 -*-
################################################################
import modules_1
################################################################
# 
basic_cond, sim_cond, names = modules_1.ReadNWConditions.setupcondition()
# ネットワークポリマーの諸量を計算
init = modules_1.RegularNW.InitialSetup(sim_cond)
target_cond, target_name, calcd_data_dic, condition_text = init.calc_conditions()
##################
# baseUDF の作成
baseudf = modules_1.SetupInitUDF.MakeInitUDF(basic_cond, sim_cond, target_cond, names, target_name, calcd_data_dic, condition_text)
target_dir = baseudf.setup_baseudf()
###############
# シミュレーションを設定
setup = modules_1.EquivCalcSetup.SetUpUDF(basic_cond, sim_cond, target_name, target_dir, names)
setup.setup_udf()
