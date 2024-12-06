"""
extract data, using Paraview-python modules, to numpy arrays

This script will focus on getting volume-average data for scalar parameters in
the liquid phase

""" 
import numpy as np
from paraview import simple as pv
import vtk.numpy_interface.dataset_adapter as dsa 
import sys
from scipy.optimize import curve_fit

def func(X,cstar,kla):
    (t,t0,c0)=X
    return( (cstar-c0)*(1-np.exp(-kla*(t-t0)))+c0 )

ofreader = pv.OpenFOAMReader(FileName = './soln.foam') # just need to provide folder
ofreader.CaseType = 'Reconstructed Case'
ofreader.MeshRegions = ['internalMesh']
ofreader.SkipZeroTime = 0  #dont know why this is not working

t = np.array(ofreader.TimestepValues)
N = t.size
print(t)

# threshold filter to get only the "aerated liquid"; specify cell data or point
# data in first element of Scalars by CELLS or POINTS
liquidthreshold = pv.Threshold(Input=ofreader, Scalars=['CELLS', 'alpha.gas'],\
#                               ThresholdRange=[0., 0.6])
                               LowerThreshold=0,UpperThreshold=0.6,ThresholdMethod='Between')


calc1 = pv.Calculator(Input=liquidthreshold, AttributeType='Cell Data',\
                         ResultArrayName='cO2',\
                         Function='"O2.liquid"*1000*"alpha.liquid"/0.032')
#calc2 = pv.Calculator(Input=calc1, AttributeType='Cell Data',\
#                         ResultArrayName='cCO',\
#                         Function='CO.liquid*thermo:rho.liquid*alpha.liquid/0.028')
calc3 = pv.Calculator(Input=calc1, AttributeType='Cell Data',\
                         ResultArrayName='cCO2',\
                         Function='"CO2.liquid"*1000*"alpha.liquid"/0.044')

# integrate all variables (in liquid)
integrateliq = pv.IntegrateVariables(Input=calc1)

Vl = np.zeros(N) # liquid volume; m^3
molO2 = np.zeros(N) # O2 in liquid; moles
#molCO = np.zeros(N) # H2 in liquid; moles
molCO2 = np.zeros(N) # CO2 in liquid; moles
Vg = np.zeros(N)
Vt = np.zeros(N)


for i in range(N):
    print("processing time = %g" % t[i])
    pv.UpdatePipeline(time=t[i], proxy=integrateliq)
    idat = dsa.WrapDataObject( pv.servermanager.Fetch(integrateliq) )
    Vt[i] = idat.CellData['Volume'].item() 
    Vl[i] = idat.CellData['alpha.liquid'].item()
    Vg[i] = idat.CellData['alpha.gas'].item()
    molO2[i] = idat.CellData['cO2'].item()
#    molCO[i] = idat.CellData['cCO'].item()
    molCO2[i] = idat.CellData['cCO2'].item()

ag = Vg/Vt # m^3/m^3
cO2 = molO2/Vl
#cCO = molCO/Vl
cCO2 = molCO2/Vl

params=0.5,0.1
fitparamsO2,cov=curve_fit(func, (t,t[0],cO2[0]), cO2, params)

#params=0.5,0.1
#fitparamsCO,cov=curve_fit(func, (t,t[0],cCO[0]), cCO, params)

params=0.5,0.1
fitparamsCO2,cov=curve_fit(func, (t,t[0],cCO2[0]), cCO2, params)


cO2fit=(fitparamsO2[0]-cO2[0])*(1-np.exp(-fitparamsO2[1]*(t-t[0])))+cO2[0]
#cCOfit=(fitparamsCO[0]-cCO[0])*(1-np.exp(-fitparamsCO[1]*(t-t[0])))+cCO[0]
cCO2fit=(fitparamsCO2[0]-cCO2[0])*(1-np.exp(-fitparamsCO2[1]*(t-t[0])))+cCO2[0]
np.savetxt("fitting.dat",np.transpose(np.vstack((t,cO2,cO2fit,cCO2,cCO2fit))),delimiter="  ")

np.savetxt("cstar_kla.dat",np.vstack((fitparamsO2,fitparamsCO2)),delimiter="   ")
print("cstar,kla O2:",fitparamsO2[0],fitparamsO2[1])
#print("cstar,kla CO:",fitparamsCO[0],fitparamsCO[1])
print("cstar,kla CO2:",fitparamsCO2[0],fitparamsCO2[1])
