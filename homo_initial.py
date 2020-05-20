#!/usr/bin/env python
# -*- coding: utf-8 -*-
#######################
try:
	from setups import Homo_initial_setup
except ImportError:
	import sys
	import os
	import platform
	#
	if platform.system() == "Windows":
		py_mod = 'D:/Dropbox/python_modules'
	elif platform.system() == "Linux":
		py_mod = '/home/hiroshi/Dropbox/python_modules'
	#
	if os.path.exists(py_mod):
		sys.path.append(py_mod)
		#
		from setups import Homo_initial_setup
	else:
		sys.exit("path for modules is not correct\nPlease check it!!")
		
##########################################
# 使用するCognacのバージョンを入れてください。
ver_cognac = "cognac101"
template = 'cognac101.udf'
# 計算時のコア数
core = 2
#
calc_cond = [ver_cognac, template, core]
##########################################
### 計算条件
##########################################
### 作成するポリマー ###
# polymer = [セグメント数、本数]
polymer = [100, 200]
### Time ###
# time = [dt, total, step]
### Bond ###
# bond = "FENE"
# bond = "Harmonic"
# bond = "Table"
### Non_bond ###
# non_bond = "rep_LJ"
# non_bond = "ind_LJ"
# non_bond = "intra_LJ"
# non_bond = "none"
# non_bond = "Table"
### Angle ###
# angle = 'Theta2'
# angle = 'Table'
# angle = 'none'
#######################
# target_list = [time, bond, non_bond, angle]
target_list = [
				[[0.01, 100000, 10000], 1.1], 
				[[0.01, 100000, 10000], 0.9558*2.**(1./6.)], 
				[[0.01, 500000, 10000], 1.05],
				[[0.01, 500000, 10000], 1.0], 
				[[0.01, 1000000, 10000], 0.9], 
				[[0.01, 1000000, 10000], 0.8], 
				]
#
################################################################################
### MAIN
################################################################################
def main():
	target_cond = make_target_cond(target_list)

	hs = Homo_initial_setup.MakePolymer(py_mod, calc_cond, polymer, target_cond)
	hs.setup_all()


def make_target_cond(target_list):
	target_cond = []
	for t_list in target_list:
		rc_nb = t_list[1]
		#
		r1 = 0.9609
		sigma = 1.0
		epsilon = 1.0
		table_range = [0., 1.5, 1000]
		cond_bond = ["Table", [r1, sigma, epsilon], table_range]
		#
		cutoff = 2**(1/6)
		scale_1_4 = 1.0
		epsilon = 1.
		sigma = 1.
		range = 1.0
		table_range = [0., 1.5, 1000]
		cond_nonbond = ["Table", [range, cutoff, scale_1_4], [cutoff, sigma, epsilon, rc_nb], table_range]
		#
		theta = 74
		angle_k = 20
		cutoff = 2**(1/6)
		epsilon = 1.
		sigma = 1.
		rc_angle = 0.8
		table_range = [-1, 1, 1000]
		cond_angle = ["Table", [theta, angle_k], [cutoff, sigma, epsilon, rc_angle], table_range]

		#
		target_cond.append([t_list[0], cond_bond, cond_nonbond, cond_angle])

	return target_cond


################################################################################
if __name__=='__main__':
	main()
################################################################################
