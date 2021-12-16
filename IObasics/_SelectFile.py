#!/usr/bin/env python
# -*- coding: utf-8 -*-
################################################################################
import sys
import os
import glob
################################################################################
# 対象ファイルを一つ選択
def get_target():
	if len(sys.argv) == 2:
		target = sys.argv[1]
	else:
		print("Input file_name")
		sys.exit()
	return target

# 
def select_outudf():
	if len(sys.argv) == 2:
		if sys.argv[1].split("_")[-1] == "out.udf":
			target = sys.argv[1]
		else:
			sys.exit("Select hnya_out.udf !")
	else:
		sys.exit("Input file_name !")
	return target

# 条件に合致するファイルを多数選択
def file_select(target = "*out.udf"):
	sorted_list = sorted(glob.glob(target))
	return sorted_list