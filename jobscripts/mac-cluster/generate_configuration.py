#!/usr/bin/python
import sys
import xml.etree.ElementTree as ET

root       = ET.Element("lbm-configuration")
physics    = ET.SubElement(root, "physics")
grid       = ET.SubElement(root, "grid")
simulation = ET.SubElement(root, "simulation")
device     = ET.SubElement(root, "device")

velocity = ET.SubElement(physics, "velocity")
ET.SubElement(velocity, "x").text = "1.0"
ET.SubElement(velocity, "y").text = "0"
ET.SubElement(velocity, "z").text = "0"
acceleration = ET.SubElement(physics, "acceleration")
ET.SubElement(acceleration, "x").text = "0"
ET.SubElement(acceleration, "y").text = "-9.81"
ET.SubElement(acceleration, "z").text = "0"
ET.SubElement(physics, "viscosity").text = "-1"
ET.SubElement(physics, "max-velocity").text = "0.1"
ET.SubElement(physics, "max-acceleration").text = "0.1"

domainsize = ET.SubElement(grid, "domain-size")
ET.SubElement(domainsize, "x").text = "4"
ET.SubElement(domainsize, "y").text = "4"
ET.SubElement(domainsize, "z").text = "4"
domianlength = ET.SubElement(grid, "domian-length")
ET.SubElement(domianlength, "x").text = "1.0"
ET.SubElement(domianlength, "y").text = "1.0"
ET.SubElement(domianlength, "z").text = "1.0"
subdomainnum = ET.SubElement(grid, "subdomain-num")
ET.SubElement(subdomainnum, "x").text = "2"
ET.SubElement(subdomainnum, "y").text = "1"
ET.SubElement(subdomainnum, "z").text = "1"
cpusubdomainratio = ET.SubElement(grid, "cpu-subdomain-ratio")
ET.SubElement(cpusubdomainratio, "x").text = "0.5"
ET.SubElement(cpusubdomainratio, "y").text = "0.5"
ET.SubElement(cpusubdomainratio, "z").text = "0.5"

ET.SubElement(simulation, "loops").text = "8"
ET.SubElement(simulation, "timestep").text = "0.01"
benchmark = ET.SubElement(simulation, "benchmark")
ET.SubElement(benchmark, "do").text = "0"
ET.SubElement(benchmark, "output-dir").text = str(sys.argv[5]) + "/dissertation/lbm/benchmark"
logging = ET.SubElement(simulation, "logging")
ET.SubElement(logging, "do").text = "0"
validation = ET.SubElement(simulation, "validation")
ET.SubElement(validation, "do").text = "1"
ET.SubElement(validation, "output-dir").text = str(sys.argv[5]) + "/dissertation/lbm/validation"
visualization = ET.SubElement(simulation, "visualization")
ET.SubElement(visualization, "do").text = "0"
ET.SubElement(visualization, "rate").text = "31694"
ET.SubElement(visualization, "output-dir").text = str(sys.argv[5]) + "/dissertation/lbm/visualization"

gridconfiguration = ET.SubElement(device, "grid-configuration")
initgridconfiguration = ET.SubElement(gridconfiguration, "init-grid-configuration")
ET.SubElement(initgridconfiguration, "x").text = "32"
ET.SubElement(initgridconfiguration, "y").text = "16"
ET.SubElement(initgridconfiguration, "z").text = "1"
alphagridconfiguration = ET.SubElement(gridconfiguration, "alpha-grid-configuration")
ET.SubElement(alphagridconfiguration, "x").text = "32"
ET.SubElement(alphagridconfiguration, "y").text = "8"
ET.SubElement(alphagridconfiguration, "z").text = "1"
betagridconfiguration = ET.SubElement(gridconfiguration, "beta-grid-configuration")
ET.SubElement(betagridconfiguration, "x").text = "32"
ET.SubElement(betagridconfiguration, "y").text = "8"
ET.SubElement(betagridconfiguration, "z").text = "1"

tree = ET.ElementTree(root)
tree.write(str(sys.argv[4]) + "/workspace/lbm/configurations/mac-cluster_" + str(int(sys.argv[1]) * int(sys.argv[2])) + ".xml")

