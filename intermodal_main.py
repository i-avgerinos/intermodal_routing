import time
from math import inf
import warnings
warnings.filterwarnings("ignore")
from pyomo.environ import *
import json
import subprocess
import os

import intermodal_milp
import intermodal_preprocessing

maxTime = 18000 # Time limit of the algorithm - Algorithm terminates after 'maxTime' seconds, if a convergence is not reached.
maxIterations = 100 # Maximum number of iterations - Algorithm terminates after 'maxIterations', if a convergence is not reached.
timelimit = 300 # Time limit of the master problem per iteration

folders = ['exact', 'alg']
for f in folders:
	filesList = [file for file in os.listdir(f"./instances/{f}")]
	for file in filesList:
		with open("results.csv", 'a') as outfile: # 'results.csv' is updated - the instance is initialised.
			outfile.write(f"{file.replace('.json', '')}\n")
			outfile.close()
		jsonDict = open(f'./instances/{f}/{file}', encoding='utf-8') # Load of json input file
		input_data = json.load(jsonDict)
		if 'exact' in file:
			configurations = ['single']
		else:
			configurations = ['single', 'multi']
		for c in configurations:
			# Preprocessing of input
			demand, dueTimes, dedicated, containers, services, origin, destination, departure, arrival, vCost, fCost, Q, coord_x, coord_y, urgency, nOrders, nHubs, nSatellites, capacity, allHubs, allSat, vehicles, order_to_service = intermodal_preprocessing.load_dataset(input_data)
			if c == 'single': # For the single-objective configuration, delay penalties are neutralised.
				for n in range(len(dueTimes)):
					dueTimes[n] = 20000
			# Initialisation of parameters: LB, UB, Number of iterations, LB of 'S'
			best_lb, best_ub, iteration, min_lb = 0, inf, 1, inf
			# Formulate the master problem
			master = intermodal_milp.master_problem(nOrders, allHubs, allSat, demand, services, origin, destination, departure,
													arrival, Q, vCost, fCost, capacity, dedicated, containers, order_to_service)
			# Define solver
			slv = SolverFactory('cplex')
			start_time = time.time()
			print("------------------------------------------------------------")
			print(file.replace('.json', ''))
			print("------------------------------------------------------------")
			print(f"Orders :		{nOrders}")
			print(f"Hubs :			{nHubs}")
			print(f"Satellites :	{nSatellites}")
			print(f"Containers :	{containers[nHubs-1][len(containers[nHubs-1])-1]+1}")
			print(f"Services :		{len(services)}")
			print("------------------------------------------------------------")
			if "exact" in file.replace('.json', ''):
				print("#	LB		UB		Gap		Time")
			if "alg" in file.replace('.json', ''):
				print("#	LB		UB		Gap		Time	lUB")
			print("------------------------------------------------------------")
			# Start of algorithm
			while (iteration <= maxIterations) or (time.time() - start_time <= maxTime):
				with open("Results.csv", 'a') as outfile:
					print(f"{iteration}	", end = "")
					outfile.write(str(iteration)+";")
					if time.time() - start_time + timelimit > maxTime:
						slv.options['timelimit'] = maxTime - time.time() + start_time
					else:
						slv.options['timelimit'] = timelimit
					results_obj = slv.solve(master, tee = False, warmstart = True)

					if (results_obj.solver.status == SolverStatus.ok) and (results_obj.solver.termination_condition == TerminationCondition.optimal):
						lb = round(value(master.z), 2)
					else:
						lb = round(results_obj.Problem._list[0].lower_bound, 2)
					if lb > best_lb:
						best_lb = lb
					best_lb = min([round(best_lb, 2), round(best_ub, 2)])
					print(f"{best_lb}	", end = "")
					outfile.write(str(best_lb)+";")
					if best_lb < best_ub or (time.time() - start_time > maxTime):
						if time.time() - start_time < maxTime:
							ub, lub = sum([fCost[s]*value(master.h[s]) for s in master.Services]) + sum([vCost[s]*value(master.x[g, s]) for g in master.Containers for s in master.Services]), sum([fCost[s]*value(master.h[s]) for s in master.Services]) + sum([vCost[s]*value(master.x[g, s]) for g in master.Containers for s in master.Services])
							solution, releaseTimes, upper_bounds, false_bounds = {}, [0.0 for n in range(nOrders)], [0.0 for i in range(nSatellites)], [0.0 for i in range(nSatellites)]
							for i in range(len(allSat)):
								solution[i] = []
								for n in master.Orders:
									for g in master.Containers:
										for s in master.Services:
											if destination[s] == allSat[i] and value(master.l[n, g, s]) > 0.9:
												solution[i].append(n)
												releaseTimes[n] = arrival[s]
							intermodal_preprocessing.generate_problemFile(nSatellites, allSat, coord_x, coord_y, vehicles, capacity, nOrders, solution, releaseTimes, arrival, dueTimes, urgency, demand, destination)
							subproblem = subprocess.run(["python", "intermodal_subproblem.py"], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
							for i in range(len(allSat)):
								os.remove(f"Log_problemFile_Depot_{i}.txt")
							upper_bounds = intermodal_preprocessing.read_solutionFile(upper_bounds)

							for i in range(len(allSat)):
								if len(solution[i]) > 0 and (time.time() - start_time <= maxTime):
									subproblem = intermodal_milp.exact_vrp(allSat, allSat[i], solution[i], releaseTimes, coord_x, coord_y, vehicles[i], demand, capacity)
									slv.options['timelimit'] = 60
									results_obj = slv.solve(subproblem, tee = False, warmstart = True)
									false_bounds[i] = round(results_obj.Problem._list[0].lower_bound, 2)
									upper_bounds[i] = round(value(subproblem.z), 2)
									if (abs(false_bounds[i] - upper_bounds[i]) <= 0.1):
										if lb < min_lb:
											min_lb = lb
									if "alg" in file.replace('.json', ''):
										master.R[allSat[i]] = upper_bounds[i]
										ub += upper_bounds[i]
										lub += false_bounds[i]
									if "exact" in file.replace('.json', ''):
										master.R[allSat[i]] = false_bounds[i]
										ub += false_bounds[i]
									var_1 = []
									for n in master.Orders:
										if n in solution[i]:
											for g in master.Containers:
												for s in master.Services:
													if destination[s] == allSat[i] and arrival[s] >= releaseTimes[n]:
														var_1.append(master.l[n, g, s])
									if "alg" in file.replace('.json', ''):
										intermodal_milp.benders_cuts(master, var_1, upper_bounds[i], allSat[i], solution[i])
									if "exact" in file.replace('.json', ''):
										intermodal_milp.benders_cuts(master, var_1, false_bounds[i], allSat[i], solution[i])
						if ub < best_ub:
							best_ub = ub
						gap = 100*(round((best_ub - best_lb)/best_ub, 4))
						if "exact" in file.replace('.json', ''):
							print(f"{round(best_ub, 2)}	{round(gap, 2)}	{int(time.time() - start_time)}")
							outfile.write(str(round(best_ub, 2))+";"+str(round(gap, 2))+";"+str(int(time.time() - start_time))+"\n")
						if "alg" in file.replace('.json', ''):
							print(f"{round(best_ub, 2)}	{round(gap, 2)}	{int(time.time() - start_time)}		{round(lub, 2)}")
							outfile.write(str(round(best_ub, 2))+";"+str(round(gap, 2))+";"+str(int(time.time() - start_time))+";"+str(round(max([lub, min_lb]), 2))+"\n")
						outfile.close()
						master.z = 1000000
						iteration += 1
						if gap <= 0.0 or (time.time() - start_time > maxTime):
							break
					else:
						gap = 100*(round((best_ub - best_lb)/best_ub, 4))
						print(f"{round(best_ub, 2)}	{round(gap, 2)}		{int(time.time() - start_time)}")
						outfile.write(str(round(best_ub, 2))+";"+str(round(gap, 2))+";"+str(int(time.time() - start_time))+"\n")
						outfile.close()
						break