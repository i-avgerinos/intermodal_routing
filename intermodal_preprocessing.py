import math

def euclidean_distance(x1, x2, y1, y2):
	distance = math.sqrt((x1-x2)**2 + (y1-y2)**2)
	return distance

def load_dataset(input_data):
	nOrders = input_data["general"]["Orders"]
	nHubs = input_data["general"]["Hubs"]
	nSatellites = input_data["general"]["Satellites"]
	capacity = input_data["general"]["Capacity"]

	demand, dueTimes, urgency = [], [], []
	for n in range(nOrders):
		demand.append(input_data["orders"][str(n)]["demand"])
		dueTimes.append(input_data["orders"][str(n)]["dueTime"])
		urgency.append(input_data["orders"][str(n)]["urgency"])
	counter, allHubs, allSat = nOrders, [], []
	for j in range(nHubs):
		allHubs.append(counter)
		counter += 1
	for i in range(nSatellites):
		allSat.append(counter)
		counter += 1

	dedicated, containers = [], []
	for j in range(nHubs):
		dedicated.append(input_data["hubs"][str(allHubs[j])]["orders"])
		containers.append(input_data["hubs"][str(allHubs[j])]["containers"])
	vehicles = []
	for i in range(nSatellites):
		vehicles.append(input_data["satellites"][str(allSat[i])]["vehicles"])

	services, origin, destination, departure, arrival, vCost, fCost, Q = [], [], [], [], [], [], [], []
	for s in input_data["services"].keys():
		services.append(int(s))
		origin.append(input_data["services"][s]["origin"])
		destination.append(input_data["services"][s]["destination"])
		departure.append(input_data["services"][s]["departure"])
		arrival.append(input_data["services"][s]["arrival"])
		vCost.append(input_data["services"][s]["cost per container"])
		fCost.append(input_data["services"][s]["fixed cost"])
		Q.append(input_data["services"][s]["maxContainers"])

	coord_x, coord_y = [], []
	for i in input_data["nodes"].keys():
		coord_x.append(input_data["nodes"][i]["x"])
		coord_y.append(input_data["nodes"][i]["y"])
	
	order_to_service = [[0.0 for s in services] for n in range(nOrders)]
	for n in range(nOrders):
		for s in services:
			if destination[s] in allSat:
				minCost = 1000000
				for l in range(nOrders):
					if l != n:
						travelCost = euclidean_distance(coord_x[n], coord_x[l], coord_y[n], coord_y[l])
						if minCost > travelCost:
							minCost = travelCost
				travelCost = euclidean_distance(coord_x[n], coord_x[destination[s]], coord_y[n], coord_y[destination[s]])
				if travelCost < minCost:
					minCost = travelCost
				delayCost = urgency[n]*max([0, math.ceil((arrival[s] + minCost - dueTimes[n])/24)])
				order_to_service[n][s] = minCost + delayCost

	return demand, dueTimes, dedicated, containers, services, origin, destination, departure, arrival, vCost, fCost, Q, coord_x, coord_y, urgency, nOrders, nHubs, nSatellites, capacity, allHubs, allSat, vehicles, order_to_service

def generate_problemFile(nSatellites, allSat, coord_x, coord_y, vehicles, capacity, nOrders, solution, releaseTimes, arrival, dueTimes, urgency, demand, destination):
	with open("problemFile.txt", 'w') as outfile:
		outfile.write("Depots:	"+str(nSatellites)+"\n")
		outfile.write("Number.	X,	Y,	Vehicles,	Capacity\n")
		for i in range(nSatellites):
			outfile.write(str(i)+".	"+str(coord_x[allSat[i]])+",	"+str(coord_y[allSat[i]])+",	"+str(vehicles[i])+",	"+str(capacity)+"\n")
		outfile.write("Customers:	"+str(nOrders)+"\n")
		outfile.write("Number.	X,	Y,	Release time,	Due-time,	Urgency,	Demand,	Depot\n")
		for n in range(nOrders):
			for i in solution.keys():
				if n in solution[i]:
					outfile.write(str(n)+".	"+str(coord_x[n])+",	"+str(coord_y[n])+",	"+str(releaseTimes[n])+",	"+str(dueTimes[n])+",	"+str(urgency[n])+",	"+str(demand[n])+",	"+str(i)+"\n")
		outfile.close()

def read_solutionFile(upper_bounds):
	sol, delayMeasurement = open('S_problemFile.txt', 'r'), False
	for line in sol:
		if "***" in line:
			currentSat = ""
			for char in line:
				if char.isnumeric():
					currentSat += char
			currentSat, delayMeasurement = int(currentSat), False
		if "Route Routing Cost" in line:
			newCost = ""
			for char in line:
				if char.isnumeric():
					newCost += char
				if char == ".":
					newCost += char
			upper_bounds[currentSat] += float(newCost)
		if delayMeasurement:
			newCost = ""
			for char in line:
				if char.isnumeric():
					newCost += char
				if char == " ":
					upper_bounds[currentSat] += int(newCost)
					newCost = ""
			delayMeasurement = False
		if "Delay_Costs" in line:
			delayMeasurement = True
	return upper_bounds