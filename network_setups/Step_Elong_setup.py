#!/usr/bin/env python
# -*- coding: utf-8 -*-
##### Import #####
from UDFManager import *
import sys
import os
import shutil
import numpy as np
import platform
import pandas as pd
################################################################################
def main():
	print("このスクリプトを直接読んでも、初期状態が入っていません。")
	exit(1)
################################################################################
class Setup:
	def __init__(self, py_mod, ver_cognac, calc_dir, f_data, core, stress_eval, gt_eval, cond_list, gt_eval_all):
		self.py_mod = py_mod
		self.ver = ver_cognac
		self.calc_dir = calc_dir
		self.f_data = f_data
		self.core = ' -n ' + str(core) +' \n'
		self.gt_eval_all = gt_eval_all
		#
		self.step_rate = cond_list[0]
		self.step_delta_t = cond_list[1]
		self.lamd_max = cond_list[2]
		#
		self.base_udf = "base_uin.udf"

	############################################################################
	##### Main #####
	# 
	def calc_step(self):
		
		conditions = self.calc_cond()
		# 平衡計算したUDFファイルとその場所を選択
		t_udf, data_dir = self.file_select()
		#
		func, nu, structure = self.read_data(data_dir)
		#
		base = self.make_base_udf(t_udf)
		#
		self.make_batch(base, conditions)
		#
		self.make_script_ss(func, nu, structure)
		#
		self.make_script_gt(func, nu, structure)
		#
		self.make_script_gt_all(func, nu, structure)

		return

	############################################################################
	##### Function #####
	# 平衡計算したUDFファイルとその場所を選択
	def file_select(self):
		param = sys.argv
		if len(param) == 1:
			print("usage: python", param[0], "Honya_out.udf")
			exit(1)
		elif not os.access(param[1],os.R_OK):
			print(param[1], "not exists.")
			exit(1)
		else:
			target = param[1]
			data_dir = os.path.dirname(target)
		return target, data_dir

	# 
	def read_data(self, data_dir):
		t_data = os.path.join(data_dir, self.f_data)
		if not os.path.exists(t_data):
			func = 0
			nu = 0
			structure = 'none'
		else:
			with open(t_data, 'r') as f:
				calc_cond = f.readlines()[-1].strip().split('\t')
			func = calc_cond[3]
			nu = calc_cond[4]
			structure = calc_cond[5]
		return func, nu, structure

	# 
	def make_base_udf(self, t_udf):
		if os.path.exists(self.calc_dir):
			print("Use existing dir of ", self.calc_dir)
		else:
			print("Make new dir of ", self.calc_dir)
			os.makedirs(self.calc_dir)
		#
		base = os.path.join(self.calc_dir, self.base_udf)
		print("Readin file = ", t_udf)
		u = UDFManager(t_udf)
		u.jump(1)
		u.eraseRecord(record_pos=-999,record_num=-999)
		u.write(base)
		return base

	# シミュレーションコンディションをプリント
	def calc_cond(self):
		temp_scale = 1
		acc_time = 0
		conditions = []
		# ステップ伸長
		delta_t = self.step_delta_t
		elong_time = (self.lamd_max - 1)/self.step_rate
		t_steps = int(elong_time/delta_t)
		cycles = t_steps
		interval = 1
		while cycles > 1000:
			interval *= 10
			cycles = t_steps/interval
		time_interval = interval*delta_t
		acc_time += delta_t*t_steps
		conditions.append([delta_t, t_steps, interval, acc_time, temp_scale])
		
		# Quench
		quench_cond = []
		# 1st
		delta_t = self.step_delta_t
		digit = int(len(str(elong_time).split('.')[0]))
		t_steps = int((10**digit - elong_time)/delta_t)
		acc_time += delta_t*t_steps
		quench_cond.append([delta_t, t_steps, interval, acc_time, temp_scale])

		# Quench_all
		while acc_time <= 1e6:
			if acc_time < 100:
				delta_t = self.step_delta_t
			else:
				delta_t = self.step_delta_t*10
				temp_scale = 0
			digit_pres = digit + 1
			end = int(10**digit_pres)
			time_interval_pres = time_interval*10
			interval_pres = int(time_interval_pres/delta_t)
			t_steps = int((end - 10**digit)*interval_pres/time_interval_pres)
			acc_time += delta_t*t_steps
			quench_cond.append([delta_t, t_steps, interval_pres, acc_time, temp_scale])
			#
			digit = digit_pres
			time_interval = time_interval_pres
		conditions.append(quench_cond)
		#
		print('\n#################################')
		elong_cond = conditions[0]
		print('##### Step Elong Conditions #####')
		print('#################################')
		print(f'{"Lambda Max":^14}', '=', self.lamd_max)
		print(f'{"Elong_rate":^14}', '=', self.step_rate)
		print('#################################')
		print(f"{'Sim. Period':^14}", '=', elong_cond[0]*elong_cond[1])
		print(f"{'delta_T':^14}", '=', elong_cond[0])
		print(f"{'Records':^14}", '=', int(elong_cond[1]/elong_cond[2]))
		if elong_cond[4] == 1:
			print(f'{"# Temp. Scalling is ON! #":^30}')
		print('#################################')
		print()
		#
		quench_cond = conditions[1]
		print('#################################')
		print('#####   Quench Conditions   #####')
		print('#################################')
		for i, cond in enumerate(quench_cond):
			acc_time += cond[0]*cond[2]
			print(f'#### Quench untill: {cond[3]: .1e} ####')
			print(f"{'Quench Period':^14}", '=', cond[0]*cond[1])
			print(f"{'delta_T':^14}", '=', cond[0])
			print(f"{'Records':^14}", '=', int(cond[1]/cond[2]))

			if cond[4] == 1:
				print(f'{"# Temp. Scalling is ON! #":^30}')
			print()
		return conditions

	# ファイル名を設定し、バッチファイルを作成
	def make_batch(self, base, conditions):
		batch = "#!/bin/bash\n\n"
		lambda_str = str(self.lamd_max).replace('.', '_')
		# rate_str = "{0:4.0e}".format(self.step_rate)
		rate_str = str(f"{self.step_rate:.1e}").replace('.', '_')
		uin = 'Step_Elong_lambda_' + lambda_str + '_rate_' + rate_str + "_uin.udf"
		print(uin)
		batch = self.make_title(batch, "Calculating Step rate=" + rate_str)
		read_udf, batch = self.make_step(uin, batch)
		self.step_elong(base, uin, conditions[0])
		pre = read_udf
		template = uin
		# ステップ伸長後にSSカーブ及び温度をプロット
		batch += 'python calc_ss.py\n'
		for i, cond in enumerate(conditions[1]):
			present_udf = 'Quench_untill_' + str(f'{cond[3]:.0e}') + "_uin.udf"
			batch = self.make_title(batch, "Calculating Quench")
			read_udf, batch = self.make_step(present_udf, batch)
			self.eq_udf(template, pre, present_udf, cond)
			pre = read_udf
			template = present_udf
			# ステップ伸長後にSSカーブ及び温度をプロット
			batch += 'python calc_gt.py\n'		
		
		# バッチファイルを作成
		f_batch = os.path.join(self.calc_dir, '_step_elong.bat')
		with open(f_batch, 'w') as f:
			f.write(batch)
			if platform.system() == "Linux":
				os.chmod(f_batch, 0o777)
		return

	###########################
	# ターミナルのタイトルを設定
	def make_title(self, batch, title):
		if platform.system() == "Windows":
			batch += "title " + title + "\n"
		elif platform.system() == "Linux":
			batch += r'echo -ne "\033]0; ' + title + ' \007"' + '\n'
		return batch

	# ファイル名の処理
	def make_step(self, present_udf, batch):
		out_udf = present_udf.replace("uin", "out")
		batch += self.ver + ' -I ' + present_udf + ' -O ' + out_udf + self.core
		read_udf = out_udf
		return read_udf, batch
		
	#-----
	def step_elong(self, base, udf_in, time_temp):
		u = UDFManager(base)
		u.jump(-1)

		# Dynamics_Conditions
		p = 'Simulation_Conditions.Dynamics_Conditions.'
		u.put(100000.,		p + 'Max_Force')
		u.put(time_temp[0],	p + 'Time.delta_T')
		u.put(time_temp[1],	p + 'Time.Total_Steps')
		u.put(time_temp[2],	p + 'Time.Output_Interval_Steps')
		#
		u.put(1.0,			p + 'Temperature.Temperature')
		u.put(time_temp[4], p + 'Temperature.Interval_of_Scale_Temp')
		#
		u.put(0,			p + 'Pressure_Stress.Pressure')

		# Deformation
		p = "Simulation_Conditions.Dynamics_Conditions.Deformation."
		u.put('Cell_Deformation', 		p + 'Method')
		u.put('Simple_Elongation', 		p + 'Cell_Deformation.Method')
		u.put('Initial_Strain_Rate', 	p + 'Cell_Deformation.Simple_Elongation.Input_Method')
		u.put(self.step_rate,	 		p + 'Cell_Deformation.Simple_Elongation.Initial_Strain_Rate.Rate')
		u.put(0.5, 						p + 'Cell_Deformation.Simple_Elongation.Poisson_Ratio')
		u.put('z', 						p + 'Cell_Deformation.Simple_Elongation.Axis')
		u.put(1, 						p + 'Cell_Deformation.Interval_of_Deform')
		u.put(0, 						p + 'Cell_Deformation.Deform_Atom')

		# Read_Set_of_Molecules
		p = 'Initial_Structure.Read_Set_of_Molecules'
		u.put([self.base_udf, -1], p)

		# Generate_Method
		p = 'Initial_Structure.Generate_Method.'
		u.put('Restart', 					p + 'Method')
		u.put([self.base_udf, -1, 1, 1], 	p + 'Restart')

		# Relaxation
		p = 'Initial_Structure.Relaxation.'
		u.put(0, p + 'Relaxation')

		#--- Write UDF ---
		u.write(os.path.join(self.calc_dir, udf_in))
		return
	
	################################################################################
	def eq_udf(self, template, read_udf, present_udf, time_temp):
		u = UDFManager(os.path.join(self.calc_dir, template))
		u.jump(-1)

		#--- Simulation_Conditions ---
		# Dynamics_Conditions
		p = 'Simulation_Conditions.Dynamics_Conditions.'
		u.put(time_temp[0], p+'Time.delta_T')
		u.put(time_temp[1], p+'Time.Total_Steps')
		u.put(time_temp[2], p+'Time.Output_Interval_Steps')
		u.put(time_temp[4], p + 'Temperature.Interval_of_Scale_Temp')

		# Deformation
		p = "Simulation_Conditions.Dynamics_Conditions.Deformation."
		u.put('None',	p + 'Method')

		#--- Initial_Structure ---
		p = 'Initial_Structure.Read_Set_of_Molecules'
		u.put([self.base_udf, -1], p)
		p = 'Initial_Structure.Generate_Method.'
		u.put('Restart', 			p+'Method')
		u.put([read_udf, -1, 1, 1], p+'Restart')

		#--- Write UDF ---
		u.write(os.path.join(self.calc_dir, present_udf))
		return

	#######################
	# 必要なスクリプトを作成
	def make_script_ss(self, func, nu, structure):
		script = self.script_content_ss(func, nu, structure)
		with open(os.path.join(self.calc_dir, 'calc_ss.py'),'w') as f:
			f.write(script)
		return

	# スクリプトの中身
	def script_content_ss(self, func, nu, structure):
		script = '#!/usr/bin/env python \n# -*- coding: utf-8 -*-\n################################\n'
		script += 'import os \nimport sys \n'
		script += 'py_mod = "' + self.py_mod + '"\n'
		script += 'sys.path.append(py_mod)\n'
		script += 'from network_evaluation import Read_Step_Stress_all\n################################\n'
		script += 'rs = Read_Step_Stress_all.MakeSS()\n'
		script += 't_udf = rs.file_select_one()\n'
		script += 'ss_data, gt_data = rs.calc_step_stress(t_udf)\n'
		script += "data_list = [[ss_data, 'ss_step.dat'], [gt_data, 'gt_step.dat']]\n"
		script += 'rs.save_plot_step(data_list)\n'
		return script

	#######################
	# 必要なスクリプトを作成
	def make_script_gt(self, func, nu, structure):
		script = self.script_content_gt(func, nu, structure)
		with open(os.path.join(self.calc_dir, 'calc_gt.py'),'w') as f:
			f.write(script)
		return

	# スクリプトの中身
	def script_content_gt(self, func, nu, structure):
		script = '#!/usr/bin/env python \n# -*- coding: utf-8 -*-\n'
		script += '################################\n'
		script += 'import os \nimport sys \n'
		script += 'py_mod = "' + self.py_mod + '"\n'
		script += 'sys.path.append(py_mod)\n'
		script += 'from network_evaluation import Read_Step_Stress_all\n'
		script += 'import platform\nimport subprocess\n'
		script += '################################\n'
		script += 'func = ' + str(func) + '\nnu = ' + str(nu) + '\nlamd = ' + str(self.lamd_max) + '\nstructure = "' + structure + '"\n'
		script += '################################\n'
		script += 'qc = Read_Step_Stress_all.QuenchCalc()\n'
		script += 't_udf_list = qc.file_select()\n'
		script += 'series = qc.calc_gt(t_udf_list, lamd)\n'
		script += 'rs = Read_Step_Stress_all.MakeSS()\n'
		script += 'data_list = [[series, "gt_quench.dat"]]\n'
		script += 'rs.save_plot_step(data_list)\n'
		# script += 'cmd = "corr2gw < gt_series.dat > gw.dat"\n'
		# script += 'subprocess.call(cmd, shell=True)\n'
		# script += 'qc.plot(func, lamd, nu, structure)\n\n'
		return script

	#######################
	# 必要なスクリプトを作成
	def make_script_gt_all(self, func, nu, structure):
		script = self.script_content_gt_all(func, nu, structure)
		with open(os.path.join(self.calc_dir, self.gt_eval_all),'w') as f:
			f.write(script)
		return

	# スクリプトの中身
	def script_content_gt_all(self, func, nu, structure):
		script = '#!/usr/bin/env python \n# -*- coding: utf-8 -*-\n'
		script += '################################\n'
		script += 'import os \nimport sys \n'
		script += 'py_mod = "' + self.py_mod + '"\n'
		script += 'sys.path.append(py_mod)\n'
		script += 'from network_evaluation import Read_Step_Stress_all\n'
		script += 'import platform\nimport subprocess\n'
		script += '################################\n'
		script += 'func = ' + str(func) + '\nnu = ' + str(nu) + '\nlamd = ' + str(self.lamd_max) + '\nstructure = "' + structure + '"\n'
		script += '################################\n'
		script += 'rs = Read_Step_Stress_all.MakeSS()\n'
		script += 'out_udf = rs.file_select_one()\n'
		script += 'gt_data, ss_data = rs.calc_step_stress(out_udf)\n\n'
		script += 'qc = Read_Step_Stress_all.QuenchCalc()\n'
		script += 't_udf_list = qc.file_select()\n'
		script += 'series = qc.calc_gt(t_udf_list, lamd)\n'
		script += 'qc.save_data(series)\n'
		script += 'cmd = "corr2gw < gt_series.dat > gw.dat"\n'
		script += 'subprocess.call(cmd, shell=True)\n'
		script += 'qc.plot(func, lamd, nu, structure)\n\n'
		return script

################################################################################
#      Main     #
################################################################################
if __name__=='__main__':
	main()