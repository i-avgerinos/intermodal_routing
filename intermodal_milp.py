from pyomo.environ import *
import pandas as pd
import numpy as np
from cplex.exceptions import CplexSolverError
from docplex.cp.model import CpoModel

import intermodal_preprocessing

def master_problem(nOrders, allHubs, allSat, demand, services, origin, destination, departure,
                   arrival, Q, vCost, fCost, capacity, dedicated, containers, order_to_service):
    
    allNodes, allContainers = allHubs + allSat, sum([len(containers[j]) for j in range(len(allHubs))])
    model = ConcreteModel()
    
    model.Orders = range(nOrders)
    model.Services = services
    model.Hubs = allHubs
    model.Satellites = allSat
    model.Containers = range(allContainers)

    model.x = Var(model.Containers, model.Services, within = Binary)
    model.y = Var(model.Orders, model.Containers, within = Binary)
    model.l = Var(model.Orders, model.Containers, model.Services, within = Binary)
    model.e = Var(model.Containers, within = Binary)
    model.h = Var(model.Services, within = Binary)

    model.R = Var(model.Satellites, within = NonNegativeReals)
    model.z = Var(within = NonNegativeReals)

    def obj_rule(model):
        return model.z 
    model.obj = Objective(rule = obj_rule, sense = minimize)

    model.constraints = ConstraintList()

    model.constraints.add(model.z >= sum(fCost[s]*model.h[s] for s in model.Services) + sum(vCost[s]*model.x[g, s] for g in model.Containers for s in model.Services) + sum(model.R[i] for i in model.Satellites))

    for j in range(len(allHubs)):
        for g in containers[j]:
            for n in model.Orders:
                if n not in dedicated[j]:
                    model.constraints.add(model.y[n, g] == 0.0)
        for n in dedicated[j]:
            for g in model.Containers:
                if g not in containers[j]:
                    model.constraints.add(model.y[n, g] == 0.0)

    for j in range(len(allHubs)):
        for n in dedicated[j]:
            model.constraints.add(sum(model.y[n, g] for g in containers[j]) == 1.0)
        for g in containers[j]:
            model.constraints.add(capacity*model.e[g] >= sum(demand[n]*model.y[n, g] for n in dedicated[j]))
            subset = []
            for s in model.Services:
                if origin[s] == allHubs[j]:
                    subset.append(s)
            if len(subset) > 0:
                model.constraints.add(sum(model.x[g, s] for s in subset) == model.e[g])
                for s in subset:
                    for t in model.Services:
                        if arrival[t] <= departure[s]:
                            model.constraints.add(model.x[g, s] + model.x[g, t] <= 1.0)
            subset = []
            for s in model.Services:
                if destination[s] in model.Satellites:
                    subset.append(s)
            if len(subset) > 0:
                model.constraints.add(sum(model.x[g, s] for s in subset) == model.e[g])
                for s in subset:
                    for t in model.Services:
                        if departure[t] >= arrival[s]:
                            model.constraints.add(model.x[g, s] + model.x[g, t] <= 1.0)
            for i in allNodes:
                if i not in model.Satellites and i != allHubs[j]:
                    for t in range(int(max(arrival))):
                        s_plus, s_minus = [], []
                        for s in model.Services:
                            if origin[s] == i and departure[s] >= t:
                                s_plus.append(s)
                            if destination[s] == i and arrival[s] <= t:
                                s_minus.append(s)
                        if len(s_plus) > 0 and len(s_minus) > 0:
                            model.constraints.add(sum(model.x[g, s] for s in s_minus) == sum(model.x[g, s] for s in s_plus))
            for t in range(int(max(arrival))):
                subset = []
                for s in model.Services:
                    if departure[s] <= t and arrival[s] > t:
                        subset.append(s)
                if len(subset) > 0:
                    model.constraints.add(sum(model.x[g, s] for s in subset) <= model.e[g])
            for n in dedicated[j]:
                subset = []
                for s in model.Services:
                    model.constraints.add(model.l[n, g, s] <= model.y[n, g])
                    model.constraints.add(model.l[n, g, s] <= model.x[g, s])
                    if destination[s] in model.Satellites:
                        subset.append(s)
        for n in dedicated[j]:
            subset = []
            for s in model.Services:
                if destination[s] in model.Satellites:
                    subset.append(s)
            if len(subset) > 0:
                model.constraints.add(sum(model.l[n, g, s] for g in containers[j] for s in subset) == sum(model.y[n, g] for g in containers[j]))
            
    for g in model.Containers:
        for s in model.Services:
            model.constraints.add(model.x[g, s] <= model.h[s])
    for s in model.Services:
        model.constraints.add(sum(model.x[g, s] for g in model.Containers) <= Q[s])

    for i in model.Satellites:
        subset = []
        for s in model.Services:
            if destination[s] == i:
                subset.append(s)
        if len(subset) > 0:
            model.constraints.add(model.R[i] >= sum(order_to_service[n][s]*model.l[n, g, s] for n in model.Orders for g in model.Containers for s in subset))

    return model

def benders_cuts(master, var_1, ub, i, solution):
    master.constraints.add(master.R[i] >= ub - ub*(len(solution) - sum(l for l in var_1)))

def exact_vrp(allSat, s, solution, releaseTimes, coord_x, coord_y, vehicles, demand, capacity):
    routes, visitTimes, delayCosts = [[] for m in range(3*vehicles)], [[] for m in range(3*vehicles)], [[] for m in range(3*vehicles)]
    sol, currentSat, routeM, visitM, delayM = open('S_problemFile.txt', 'r'), None, False, False, False
    for line in sol:
        if "***" in line:
            currentSat = ""
            for char in line:
                if char.isnumeric():
                    currentSat += char
            currentSat = int(currentSat)
        if currentSat != None:
            if allSat[currentSat] == s:
                if f"Vehicle {currentSat}." in line:
                    string = str(line)
                    if string[len(string)-2] == ".":
                        currentVehicle = int(string[len(string)-1])
                    else:
                        currentVehicle = int(string[len(string)-2] + string[len(string)-1])
                if f"Route_" in line:
                    for char in line:
                        if char.isnumeric():
                            currentCopy = int(char)
                            break
                if routeM:
                    for char in line:
                        if char == " ":
                            if "Depot" in newNode:
                                routes[currentVehicle*3 + currentCopy].append(allSat[currentSat])
                            else:
                                routes[currentVehicle*3 + currentCopy].append(int(newNode))
                            newNode = ""
                        else:
                            newNode += char
                    routeM = False
                if "Nodes:" in line:
                    newNode = ""
                    routeM = True
                if visitM:
                    for char in line:
                        if char == " ":
                            visitTimes[currentVehicle*3 + currentCopy].append(float(newTime))
                            newTime = ""
                        else:
                            newTime += char
                    visitM = False
                if "Visit_Times:" in line:
                    newTime = ""
                    visitM = True
                if delayM:
                    for char in line:
                        if char == " ":
                            delayCosts[currentVehicle*3 + currentCopy].append(float(newCost))
                            newCost = ""
                        else:
                            newCost += char
                    delayM = False
                if "Delay_Costs" in line:
                    newCost = ""
                    delayM = True

    nodes = [j for j in solution]
    nodes.insert(0, s)
    fx, fy = {}, {}
    for j in nodes:
        fx[j] = {}
        for i in nodes:
            fx[j][i] = {}
            for m in range(3*vehicles):
                fx[j][i][m] = 0.0
    for j in solution:
        fy[j] = {}
        for m in range(3*vehicles):
            fy[j][m] = 0.0

    for m in range(vehicles):
        for c in range(3):
            if len(routes[m*3 + c]):
                r = routes[m*3 + c]
                for i in range(1, len(r)):
                    fx[r[i-1]][r[i]][m*3 + c] = 1.0
                    if r[i] in solution:
                        fy[r[i]][m*3 + c] = 1.0

    model = ConcreteModel()
    
    model.Orders = solution
    model.Nodes = nodes
    model.Vehicles = range(vehicles*3)

    model.x = Var(model.Nodes, model.Nodes, model.Vehicles, within = Binary)
    model.y = Var(model.Orders, model.Vehicles, within = Binary)

    model.S = Var(model.Vehicles, within = NonNegativeReals)
    model.E = Var(model.Nodes, model.Vehicles, within = NonNegativeReals)

    model.n = Var(model.Orders, model.Vehicles, within = NonNegativeIntegers)
    model.z = Var(within = NonNegativeReals)

    for m in range(vehicles):
        for c in range(3):
            for j in nodes:
                if j in solution:
                    model.y[j, m*3 + c] = fy[j][m*3 + c]
                for i in nodes:
                    model.x[i, j, m*3 + c] = fx[i][j][m*3 + c]

    def obj_rule(model):
        return model.z 
    model.obj = Objective(rule = obj_rule, sense = minimize)

    model.constraints = ConstraintList()

    model.constraints.add(model.z >= sum(round(intermodal_preprocessing.euclidean_distance(coord_x[i], coord_x[j], coord_y[i], coord_y[j]), 2)*model.x[i, j, m] for i in model.Nodes for j in model.Nodes for m in model.Vehicles))

    for i in model.Orders:
        model.constraints.add(sum(model.x[i, i, m] for m in model.Vehicles) == 0.0)
        model.constraints.add(sum(model.x[i, j, m] for j in model.Nodes for m in model.Vehicles) == 1.0)
        model.constraints.add(sum(model.x[j, i, m] for j in model.Nodes for m in model.Vehicles) == 1.0)
        for m in model.Vehicles:
            model.constraints.add(model.y[i, m] >= sum(model.x[i, j, m] for j in model.Nodes))
            model.constraints.add(sum(model.x[i, j, m] for j in model.Nodes) == sum(model.x[j, i, m] for j in model.Nodes))
            model.constraints.add(sum(model.x[s, j, m] for j in model.Orders) >= sum(model.x[i, j, m] for j in model.Nodes))
            model.constraints.add(sum(model.x[j, s, m] for j in model.Orders) >= sum(model.x[i, j, m] for j in model.Nodes))
            for j in model.Nodes:
                if j != s:
                    model.constraints.add(model.E[j, m] + intermodal_preprocessing.euclidean_distance(coord_x[j], coord_x[i], coord_y[j], coord_y[i]) - model.E[i, m] <= 10000*(1 - model.x[j, i, m]))
                else:
                    model.constraints.add(model.S[m] + intermodal_preprocessing.euclidean_distance(coord_x[j], coord_x[i], coord_y[j], coord_y[i]) - model.E[i, m] <= 10000*(1 - model.x[j, i, m]))
                    model.constraints.add(model.E[i, m] + intermodal_preprocessing.euclidean_distance(coord_x[i], coord_x[j], coord_y[i], coord_y[j]) - model.E[j, m] <= 10000*(1 - model.x[i, j, m]))
            model.constraints.add(releaseTimes[i] - model.S[m] <= 10000*(1 - sum(model.x[j, i, m] for j in model.Nodes)))

    for m in range(vehicles):
        for c in range(3):
            model.constraints.add(sum(demand[i]*model.y[i, m*3 + c] for i in model.Orders) <= capacity)
            model.constraints.add(sum(model.x[s, j, m*3 + c] for j in model.Orders) <= 1.0)
            model.constraints.add(sum(model.x[j, s, m*3 + c] for j in model.Orders) <= 1.0)
            if c > 0:
                model.constraints.add(model.S[m*3 + c] >= model.E[s, m*3 + c - 1])

    return model