#!/usr/bin/python
import sys
import xml.etree.ElementTree as ET

root       = ET.Element("lbm-configuration")
physics    = ET.SubElement(root, "physics")
grid       = ET.SubElement(root, "grid")
simulation = ET.SubElement(root, "simulation")
cpu        = ET.SubElement(root, "cpu")
gpu        = ET.SubElement(root, "gpu")

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
ET.SubElement(domainsize, "x").text = str(sys.argv[1])
ET.SubElement(domainsize, "y").text = str(sys.argv[2])
ET.SubElement(domainsize, "z").text = str(sys.argv[3])
domianlength = ET.SubElement(grid, "domian-length")
ET.SubElement(domianlength, "x").text = str(float(sys.argv[1]) / float(sys.argv[1]))
ET.SubElement(domianlength, "y").text = str(float(sys.argv[2]) / float(sys.argv[1]))
ET.SubElement(domianlength, "z").text = str(float(sys.argv[3]) / float(sys.argv[1]))
subdomainnum = ET.SubElement(grid, "subdomain-num")
ET.SubElement(subdomainnum, "x").text = str(sys.argv[4])
ET.SubElement(subdomainnum, "y").text = str(sys.argv[5])
ET.SubElement(subdomainnum, "z").text = str(sys.argv[6])
cpusubdomainratio = ET.SubElement(grid, "cpu-subdomain-ratio").text = "0.11"

ET.SubElement(simulation, "loops").text = str(sys.argv[7])
ET.SubElement(simulation, "timestep").text = "0.01"
benchmark = ET.SubElement(simulation, "benchmark")
ET.SubElement(benchmark, "do").text = str(sys.argv[8])
ET.SubElement(benchmark, "output-dir").text = str(sys.argv[15]) + "/dissertation/lbm/benchmark"
logging = ET.SubElement(simulation, "logging")
ET.SubElement(logging, "do").text = str(sys.argv[9])
ET.SubElement(logging, "output-dir").text = str(sys.argv[15]) + "/dissertation/lbm/log"
validation = ET.SubElement(simulation, "validation")
ET.SubElement(validation, "do").text = str(sys.argv[10])
ET.SubElement(validation, "output-dir").text = str(sys.argv[15]) + "/dissertation/lbm/validation"
visualization = ET.SubElement(simulation, "visualization")
ET.SubElement(visualization, "do").text = str(sys.argv[11])
ET.SubElement(visualization, "rate").text = str(sys.argv[12])
ET.SubElement(visualization, "output-dir").text = str(sys.argv[15]) + "/dissertation/lbm/visualization"

blockconfiguration = ET.SubElement(cpu, "block-configuration")
alphablockconfiguration = ET.SubElement(blockconfiguration, "alpha-block-configuration")
ET.SubElement(alphablockconfiguration, "x").text = "256"
ET.SubElement(alphablockconfiguration, "y").text = "64"
ET.SubElement(alphablockconfiguration, "z").text = "1"
betablockconfiguration = ET.SubElement(blockconfiguration, "beta-block-configuration")
ET.SubElement(betablockconfiguration, "x").text = "256"
ET.SubElement(betablockconfiguration, "y").text = "64"
ET.SubElement(betablockconfiguration, "z").text = "1"

gridconfiguration = ET.SubElement(gpu, "grid-configuration")
initgridconfiguration = ET.SubElement(gridconfiguration, "init-grid-configuration")
ET.SubElement(initgridconfiguration, "x").text = "64"
ET.SubElement(initgridconfiguration, "y").text = "1"
ET.SubElement(initgridconfiguration, "z").text = "1"
alphagridconfiguration = ET.SubElement(gridconfiguration, "alpha-grid-configuration")
ET.SubElement(alphagridconfiguration, "x").text = "32"
ET.SubElement(alphagridconfiguration, "y").text = "1"
ET.SubElement(alphagridconfiguration, "z").text = "1"
betagridconfiguration = ET.SubElement(gridconfiguration, "beta-grid-configuration")
ET.SubElement(betagridconfiguration, "x").text = "32"
ET.SubElement(betagridconfiguration, "y").text = "1"
ET.SubElement(betagridconfiguration, "z").text = "1"

tree = ET.ElementTree(root)
tree.write(str(sys.argv[14]) + "/workspace/lbm/configurations/mac-cluster_" + str(sys.argv[13]) + ".xml")

