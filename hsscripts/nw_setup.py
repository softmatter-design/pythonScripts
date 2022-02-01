#!/usr/bin/env python
# -*- coding: utf-8 -*-
################################################################
import modules_1
################################################################
# 設定条件を読み込み、ネットワークポリマーの諸量を計算
basic_cond, nw_cond, sim_cond, rnd_cond, target_cond, condition_text = modules_1.ReadNWConditions.setupcondition()
# ネットワークを設定
nwsetup = modules_1.RegNWSetup.NWSetup(nw_cond, target_cond)
calcd_data_dic = nwsetup.calc_all()

##################
# baseUDF の作成
baseudf = modules_1.SetupInitUDF.MakeInitUDF(basic_cond, nw_cond, sim_cond, target_cond, calcd_data_dic, condition_text)
target_dir = baseudf.setup_baseudf()
###############
# シミュレーションを設定
setup = modules_1.EquivCalcSetup.SetUpUDF(basic_cond, sim_cond, target_dir)
setup.setup_udf()
