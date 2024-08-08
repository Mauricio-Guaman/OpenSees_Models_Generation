# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 00:27:25 2021

@author: Mauricio
"""
import BibliotecaFunciones as bf
import matplotlib.pyplot as plt
import numpy as np
#from scipy import stats
from time import time
tiempoInicial=time()
import scipy.stats as stats
import pandas as pd
import os
import random
# Units#
#m,mm,kN,N,Pa,kPa,MPa,g=bf.units()

    
#%% VARIABLES
# Vanos
mux=4.5
sigmax=0.6
#x=np.linspace(stats.norm(mux,sigmax).ppf(0.01),stats.norm(mux,sigmax).ppf(0.99),25)
#x=[3,3.5,4,4.5,5,5.5,6,6.5] # Valor medio
x=4.5
# x=random.choice(x2)

# Concrete
mu=21
sigma=1.5
#fc=np.linspace(stats.norm(mu,sigma).ppf(0.01),stats.norm(mu,sigma).ppf(0.99),25)
fc=[21] # Valor medio resistencia concreto
# fc=[random.choice(fc)]

# rho Ash
muAsh=0.003351
sigmaAsh=0.00105
#rho_Ash=np.linspace(stats.norm(muAsh,sigmaAsh).ppf(0.01),stats.norm(muAsh,sigmaAsh).ppf(0.99),25)
rho_Ash=0.003351 # Valor medio
# rho_Ash=random.choice(rho_Ash)

# Mortar resistant join
muj=20.19
sigmaj=3.65
#fcj=np.linspace(stats.norm(muj,sigmaj).ppf(0.01),stats.norm(muj,sigmaj).ppf(0.99),25)
fcj=20.19 # Valor medio
# fcj=random.choice(fcj)

# Mortar resistant join
mub=0.97
sigmab=0.16
#fcb=np.linspace(stats.norm(mub,sigmab).ppf(0.01),stats.norm(mub,sigmab).ppf(0.99),25)
fcb=0.97 # Valor medio
# fcb=random.choice(fcb)


# Openings
#rho=np.random.seed(1)
#rho=np.random.randint(0,100,size=(25,10))/100
rho=np.array([1,1,1,1,1,1,1,1,1,1])*0.50 # Valor sin aperturas
print(rho)
# rho=random.choice(rho)


#%% DAMAGE STATES DATA
# Beams
recBeam=0.03
stirrupSpacingBeam=0.1
#Columns
# Tranversal Rein
fi_bar_stirrup=0.008
ramales=2
recCol=0.03
# Long
fi_barLong=0.012
Nbarras=6

#df inicos variables
quadrilinearPushover=pd.DataFrame()
DamageStatePlasRotation_df=pd.DataFrame()
DamageStatePlasRotationADRS1_df=pd.DataFrame()
DamageStatePlasRotationADRS_df=pd.DataFrame()
ADRSCurve=pd.DataFrame()
PushoverCurves=pd.DataFrame()
ShearDS_df=pd.DataFrame()
Damages_States=pd.DataFrame()
periods_df=pd.DataFrame()
#%% BARE FRAME
cont=1
for i in range(len(fc)):
#for i in range(1):
    #print('fcb',fcb[i-1])
    # Data
    rutaOpenSees=r'C:/OpenSees/OpenSees3.2.2-x64.exe/bin/OpenSees.exe'
    rutafolder=os.getcwd()
    rutaAxialLoad=rutafolder+'\\'+'DataElasticBareModel'+str(cont)+'\\'+'AxialLoad.out'
    rutaFirstMode=rutafolder+'\\'+'DataElasticBareModel'+str(cont)+'\\'+'Vectors'+'\\'+'eigenvector1.out'
    nameFile='ElasticBareModelGeneration'+ str(cont)+'.tcl'
    nameFileInelastic='InelasticBareModelGeneration'+str(cont)+'.tcl'
    
    #Run Model
    bf.bareFrame(rutaOpenSees,rutafolder,rutaAxialLoad,rutaFirstMode,nameFile,nameFileInelastic,cont,x,fc[i],rho_Ash)
    # Reactions
    RBaseInfillFrame=rutafolder+'\\'+'DataInelasticBareModel'+str(cont)+'\\'+'RBase.out'
    RBase=np.genfromtxt(RBaseInfillFrame) 
    Reacciones=RBase[:,1:]
    BasalSM_ori=np.apply_along_axis(sum,1,Reacciones)*-1
    # Displacement top
    DFreeInfillFrame=rutafolder+'\\'+'DataInelasticBareModel'+str(cont)+'\\'+'DFree.out'
    DFree=np.genfromtxt(DFreeInfillFrame) 
    DtopSM_ori=DFree[:,1]
    # Tranformacion de datos
    BasalM=[0]
    DtopM=[0]
    for j in range(0,len(BasalSM_ori)):
        #print(i)
        if BasalSM_ori[j]>0:
            BasalM.append(BasalSM_ori[j])
            DtopM.append(DtopSM_ori[j])
    # Grafica
    plt.figure(1)
    #plt.plot(DtopSM,BasalSM,label='Sin Mampostería fc:'+str(round(fc[i],2)))
    plt.plot(DtopM,BasalM,'darkgrey')
    # plt.title('Pushover Curve')
    # plt.xlabel('Displacement [m]')
    # plt.ylabel('Vbasal [kN]')  
    # plt.grid()        
    # plt.legend()
    # plt.show() 
    cont+=1

#%% INFILL FRAME   
# x2=4.5
# fc=[18]
cont=1
for i in range(len(fc)):
#for i in range(1):
    #print('fcj',fcj[i])
    try:
        #print('rho',rho[i,:])
        # Rutas
        rutaOpenSees=r'C:/OpenSees/OpenSees3.2.2-x64.exe/bin/OpenSees.exe'
        rutafolder=os.getcwd()
        rutaAxialLoad=rutafolder+'\\'+'DataElasticInfillModel'+str(cont)+'\\'+'AxialLoad.out'
        rutaFirstMode=rutafolder+'\\'+'DataElasticInfillModel'+str(cont)+'\\'+'Vectors'+'\\'+'eigenvector1.out'
        nameFile='ElasticInfillModelGeneration'+ str(cont)+'.tcl'
        nameFileInelastic='InelasticInfillModelGeneration'+str(cont)+'.tcl'
        rutaPlasticDeformationColumns=rutafolder+'\\'+'DataInelasticInfillModel'+str(cont)+'\\'+'Rotations'+'\\'+'Columns'+'\\'+'BasicDeformationColumns.out'
        rutaPlasticDeformationBeams=rutafolder+'\\'+'DataInelasticInfillModel'+str(cont)+'\\'+'Rotations'+'\\'+'Beams'+'\\'+'BasicDeformationBeams.out'
        rutaPlasticDeformationGirds=rutafolder+'\\'+'DataInelasticInfillModel'+str(cont)+'\\'+'Rotations'+'\\'+'Girds'+'\\'+'BasicDeformationGirds.out'
        # DIRECTORY SHEAR DAMAGE STATE
        rutaShear=os.getcwd()+'\\'+'DataInelasticInfillModel'+str(cont)+'\\''shearfolder.txt'
        rutaLocalForces=os.getcwd()+'\\'+'DataInelasticInfillModel'+str(cont)+'\\'+'LocalForce'+'\\'+'Columns'+'\\'+'LocalForceEleColumns.out'
    
        #Run Model
        massOpenSees,loadOpensees,modelo,mass=bf.InfillFrame(rutaOpenSees,rutafolder,rutaAxialLoad,rutaFirstMode,nameFile,nameFileInelastic,cont,
                      x,fc[i],rho,rho_Ash,fcj,fcb)
       #%% 
        # Save fundamental period infill frame
        period_file=os.path.join(rutafolder,'DataInelasticInfillModel'+str(cont),'Modes','Periods.txt')
        period_vector=np.genfromtxt(period_file)
        period=period_vector[0,0]
        periods_df['PeriodCurve'+str(i)]=[period]
        #%%
        # Reactions
        RBaseInfillFrame=rutafolder+'\\'+'DataInelasticInfillModel'+str(cont)+'\\'+'RBase.out'
        RBase=np.genfromtxt(RBaseInfillFrame) 
        Reacciones=RBase[:,1:]
        BasalSM_ori=np.apply_along_axis(sum,1,Reacciones)*-1
        # Displacement top
        DFreeInfillFrame=rutafolder+'\\'+'DataInelasticInfillModel'+str(cont)+'\\'+'DFree.out'
        DFree=np.genfromtxt(DFreeInfillFrame) 
        DtopSM_ori=DFree[:,1]
        # Tranformacion de datos
        BasalSM=[0]
        DtopSM=[0]
        for j in range(0,len(BasalSM_ori)):
            #print(i)
            if BasalSM_ori[j]>0:
                BasalSM.append(BasalSM_ori[j])
                DtopSM.append(DtopSM_ori[j])
        # Dataframe CSV
        Basal_df=pd.DataFrame()
        Basal_df['Basal']=BasalSM
        
        Dtop_df=pd.DataFrame()
        Dtop_df['Dtop']=DtopSM
        
        PushoverCurves['DtopCurve'+str(i)]=Dtop_df
        PushoverCurves['BasalCurve'+str(i)]=Basal_df
    
        # Grafica
        plt.figure(1)
        #plt.plot(DtopSM,BasalSM,label='Con Mampostería fc:'+str(round(fc[i],2)))
        plt.plot(DtopSM,BasalSM,'darkgrey')
        plt.title('Pushover Curve')
        plt.xlabel('Displacement [m]')
        plt.ylabel('Vbasal [kN]')  
        plt.grid()        
        plt.legend()
        plt.show()
        
        
        # ADRS
        h=[2.50,5.00,7.50] # Altura entrepisos
        aeq,Sd=bf.adrs(rutaFirstMode,mass,h,BasalSM,DtopSM)
        # DataFrame csv
        aeq_df=pd.DataFrame()
        aeq_df['Basal']=aeq
        
        Sd_df=pd.DataFrame()
        Sd_df['Dtop']=Sd
        
        
        ADRSCurve['Sd_Curve'+str(i)]=Sd_df
        ADRSCurve['Ad_Curve'+str(i)]=aeq_df
        
        # Grafica ADRS
        plt.figure(2)
        plt.plot(Sd,aeq,color='darkgrey')
        plt.title('Pushover Curve ADRS System')
        plt.xlabel('Spectral displacement [m]')
        plt.ylabel('Spectral aceleration [g]')
        plt.grid()
        
        #% QUADRILINEAR
        coorad,coorsd=bf.cuadrilinear(aeq,Sd)
        #% DATAFRAME CSV
        quadrilinearPushover['Sd Curve '+str(i)]=coorsd
        quadrilinearPushover['Ad Curve '+str(i)]=coorad
        
        # PLOT QUADRILINERA ADRS
        ad=coorad.copy()
        ad.insert(0,0)
        sd=coorsd.copy()
        sd.insert(0,0)
        plt.figure(3)
        plt.plot(sd,ad,color='darkgrey')
        plt.title('Pushover Curve ADRS System Quadrilinear')
        plt.xlabel('Spectral displacement [m]')
        plt.ylabel('Spectral aceleration [g]')
        plt.grid()
        plt.show()
        
        # DAMAGE STATE
        Columns_df,DamageStatePlasRotation=bf.DamageStatePlasticRotations(fc[i],modelo,recBeam,stirrupSpacingBeam,fi_bar_stirrup,ramales,recCol,fi_barLong,Nbarras,
                                                                rutaAxialLoad,rutafolder,DtopSM_ori,BasalSM_ori,rutaPlasticDeformationColumns,
                                                                rutaPlasticDeformationBeams,rutaPlasticDeformationGirds)       
        DamageStatePlasRotation_df['Curve'+str(i)]=DamageStatePlasRotation
        
        # ADRSCurve
        DamageStatePlasRotationADRS1_df['Basal_IO'],DamageStatePlasRotationADRS1_df['Displacement_IO']=bf.adrs(rutaFirstMode,mass,h,DamageStatePlasRotation.loc['Basal_IO'],DamageStatePlasRotation.loc['Displacement_IO'])
        DamageStatePlasRotationADRS1_df['Basal_LS'],DamageStatePlasRotationADRS1_df['Displacement_LS']=bf.adrs(rutaFirstMode,mass,h,DamageStatePlasRotation.loc['Basal_LS'],DamageStatePlasRotation.loc['Displacement_LS'])
        DamageStatePlasRotationADRS1_df['Basal_CP'],DamageStatePlasRotationADRS1_df['Displacement_CP']=bf.adrs(rutaFirstMode,mass,h,DamageStatePlasRotation.loc['Basal_CP'],DamageStatePlasRotation.loc['Displacement_CP'])
        DamageStatePlasRotationADRS_df['Curve'+str(i)]=DamageStatePlasRotationADRS1_df.transpose()
         
        
        # SHEAR DAMAGE STATE
        Shear_DS=bf.ShearDS(rutaShear,rutaLocalForces,BasalSM_ori,DtopSM_ori)
        ShearDS_df['Cruve'+str(i)]=Shear_DS 
        
        # DAMAGES STATES
        # =============================================================================
        # DS1=IO=0.70*Vmax
        ds1y=0.70*max(aeq)
        pos1=np.where(aeq>ds1y)[0][0]
        DS1=Sd[pos1]
        # =============================================================================
        # DS2=LS=RLS
        DS2=DamageStatePlasRotationADRS_df.loc['Displacement_LS'][0]
        # =============================================================================
        # DS3=CP=RCP
        DS3=DamageStatePlasRotationADRS_df.loc['Displacement_CP'][0]
        # =============================================================================
        # DS4=C=0.50*Vmax
        ds4y=0.50*max(aeq)
        pos=np.where(aeq[pos1:]<ds4y)[0][0]
        
        DS4=Sd[pos]
        DS=[DS1,DS2,DS3,DS4]
        Damages_States['Cruve'+str(i)]=DS
        
        cont+=1  
    except:
        print('ERROR: IT DOES NOT CONVERGE')
        cont+=1  
        pass
    
#%% CSV FILE
VNombre='rho'
Nfolder='GraficaPushover'

os.makedirs(os.path.join(rutafolder,Nfolder,'002PushoverCurves'),exist_ok=True)
os.makedirs(os.path.join(rutafolder,Nfolder,'003ADRSCurves'),exist_ok=True)
os.makedirs(os.path.join(rutafolder,Nfolder,'004QuadrilinearCurves'),exist_ok=True)
os.makedirs(os.path.join(rutafolder,Nfolder,'005DamageStates'),exist_ok=True)
os.makedirs(os.path.join(rutafolder,Nfolder,'006DamageStatesADRS'),exist_ok=True)
os.makedirs(os.path.join(rutafolder,Nfolder,'007Periods'),exist_ok=True)
PushoverCurves.to_csv(rutafolder+'\\'+Nfolder+'//'+'002PushoverCurves'+'\\'+'PushoverCurves_'+VNombre+'.csv',index=False,header=False)
ADRSCurve.to_csv(rutafolder+'\\'+Nfolder+'//'+'003ADRSCurves'+'\\'+'ADRSCurves_'+VNombre+'.csv',index=False,header=False)
quadrilinearPushover.to_csv(rutafolder+'\\'+Nfolder+'//'+'004QuadrilinearCurves'+'\\'+'QuadrilinearPushoverCurves_'+VNombre+'.csv',index=False,header=False)
DamageStatePlasRotation_df.to_csv(rutafolder+'\\'+Nfolder+'\\'+'005DamageStates'+'\\'+'DamageStatePlasticRotation_'+VNombre+'.csv',index=False,header=False)
Damages_States.to_csv(rutafolder+'\\'+Nfolder+'\\'+'006DamageStatesADRS'+'\\'+'DamageStatePlasticRotationADRS_'+VNombre+'.csv',index=False,header=False)
periods_df.to_csv(rutafolder+'\\'+Nfolder+'\\'+'007Periods'+'\\'+'Periods_'+VNombre+'.csv',index=False,header=False)

#%% GRAFICAS DE DAMAGES STATES

plt.figure(1)
plt.scatter(DamageStatePlasRotation_df.loc['Displacement_IO'],DamageStatePlasRotation_df.loc['Basal_IO'],label='R_IO',color='c')
plt.scatter(DamageStatePlasRotation_df.loc['Displacement_LS'],DamageStatePlasRotation_df.loc['Basal_LS'],label='R_LS',color='g')
plt.scatter(DamageStatePlasRotation_df.loc['Displacement_CP'],DamageStatePlasRotation_df.loc['Basal_CP'],label='R_CP',color='r')

plt.figure(2)
plt.scatter(DamageStatePlasRotationADRS_df.loc['Displacement_IO'],DamageStatePlasRotationADRS_df.loc['Basal_IO'],label='R_IO',color='c')
plt.scatter(DamageStatePlasRotationADRS_df.loc['Displacement_LS'],DamageStatePlasRotationADRS_df.loc['Basal_LS'],label='R_LS',color='g')
plt.scatter(DamageStatePlasRotationADRS_df.loc['Displacement_CP'],DamageStatePlasRotationADRS_df.loc['Basal_CP'],label='R_CP',color='r')
plt.legend()
plt.grid()
plt.show()

plt.figure(3)
plt.scatter(DamageStatePlasRotationADRS_df.loc['Displacement_IO'],DamageStatePlasRotationADRS_df.loc['Basal_IO'],label='R_OI',color='c')
plt.scatter(DamageStatePlasRotationADRS_df.loc['Displacement_LS'],DamageStatePlasRotationADRS_df.loc['Basal_LS'],label='R_LS',color='g')
plt.scatter(DamageStatePlasRotationADRS_df.loc['Displacement_CP'],DamageStatePlasRotationADRS_df.loc['Basal_CP'],label='R_CP',color='r')
plt.legend()
plt.grid()
plt.show()

# PLOT DAMAGE STATE SHEAR
plt.figure(1)
plt.scatter(ShearDS_df.loc['Dtop_V_LS'],ShearDS_df.loc['Basal_V_LS'],label='V_LS',color='k')
plt.scatter(ShearDS_df.loc['Dtop_V_CP'],ShearDS_df.loc['Basal_V_CP'],label='V_CP',color='y')
plt.scatter(ShearDS_df.loc['Dtop_V_C'],ShearDS_df.loc['Basal_V_C'],label='V_C',color='m')
plt.legend()
plt.grid()
plt.show()
    

#%%
tiempoFinal=time()
tiempoEjecucion=tiempoFinal-tiempoInicial
print('Tiempo de Ejecución: ', tiempoEjecucion)

#%% GRAFICA TESIS

fig=plt.figure()
ax=fig.add_subplot(111)
ax.plot(DtopM,BasalM,'b',label='Bare Frame')
ax.plot(DtopSM,BasalSM,'r',label='Infill Frame')
# Formato
ax.set_ylim(0,max(BasalSM)+10)
ax.set_xlim(0,max(DtopSM)+0.01)
#ax.set_yticks(np.arange(0,1.1,0.1))
# ticks_y=tick.FuncFormatter(lambda y,pos:'{:.2f}'.format(y/1))
# ax.yaxis.set_major_formatter(tick.FuncFormatter(ticks_y))
ax.grid(which='major', color='#DDDDDD', linewidth=0.8)
ax.grid(which='minor', color='#EEEEEE', linestyle='-', linewidth=0.5)
ax.minorticks_on()

ax.set(xlabel='Displacement [m]',ylabel='Force [kN]')
#ax.set_title('f) Fragility Curves - Variable: '+name)
ax.legend()















































