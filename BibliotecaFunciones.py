# -*- coding: utf-8 -*-
"""
Created on Sat Jul  3 21:14:49 2021

@author: Mauricio
"""
#%% General Units

#Funcion que permite acceder a las unidades de trabajo
def units():
    # Longitud
    m=1
    mm=m/1000
    #Fuerza
    kN=1
    N=kN/1000
    # Presion
    Pa=N/(m*m)
    kPa=Pa*1.0e3
    MPa=Pa*1.0e6
    # acceleration
    g=9.81 #[m/s2] gravity
    return m,mm,kN,N,Pa,kPa,MPa,g
#%% Moment Curvature

# Funcion que calcula el momento de fluencia de una seccion
def MomentCurvature (h,b,cv,dbL,dbV,fc,Ec,P,fyL,Es,rho1,rho2,rho3):
    # some general stuff
    n_c=0.8+fc/18   # Termino para el concreto 
    e_c=1*fc/Ec*(n_c/(n_c-1))   #Deformacion unitaria del concreto 
    e_s=fyL/Es   # Deformacion unitaria acero
    
    # Curvatura de fluencia
    phiy=2.1*fyL/Es/h #[rad/m] Priestley et al [2007]
    
    # Momento de fluencia 
    d1=cv+dbV+dbL/2
    d2=h/2
    d3=h-cv-dbV-dbL/2
    # Positive Bending
    c=h/2   # Initial trial of NA depth (m)
    count=0
    err=0.5
    while err>0.001 and count<10000:
        e_s1=(c-d1)*phiy
        e_s2=(d2-c)*phiy
        e_s3=(d3-c)*phiy
        e_top=c*phiy
            # Compute the steel stuff
        if e_s1<e_s:
            f_s1=e_s1*Es #[MPa]
            f_s2=e_s2*Es
            f_s3=e_s3*Es
        else:
            f_s1=fyL #[MPa]
            f_s2=fyL
            f_s3=fyL
        # Fuerza Traccion
        Fs1=rho1*b*d3*1000*f_s1
        Fs2=rho2*b*d3*1000*f_s2
        Fs3=rho3*b*d3*1000*f_s3
        
        # Compute Compression Force
        a1b1=(e_top/e_c)-((e_top/e_c)**2)/3
        b1=(4-(e_top/e_c))/(6-2*(e_top/e_c))
        Fc=a1b1*c*fc*b*1000
        
        # Section force
        Psec= P+Fs2+Fs3-Fc-Fs1
        
        # Adjust NA depth to balance section forces
        if Psec<0:
            #c0=c
            c=c-0.001
        elif Psec>0:
            #co=c
            c=c+0.001
        err=abs(Psec)
        if err<5:
            break
        count+=1
    # Compute the moment
    Mp=P*(0.5*h-c)+Fs1*(c-d1)+Fs2*(d2-c)+Fs3*(d3-c)+Fc*c*(1-b1/2)
    cp=c
    # print( Mp)
    # print( cp)
    
    # Negative Bending 
    c=h/2   # Initial trial of NA depth (m)
    count=0
    err=0.5
    while err>0.001 and count<10000:
        e_s1=(c-d1)*phiy
        e_s2=(d2-c)*phiy
        e_s3=(d3-c)*phiy
        e_top=c*phiy
            # Compute the steel stuff
        if e_s1<e_s:
            f_s1=e_s1*Es #[MPa]
            f_s2=e_s2*Es
            f_s3=e_s3*Es
        else:
            f_s1=fyL #[MPa]
            f_s2=fyL
            f_s3=fyL
        # Fuerza Traccion
        Fs1=rho3*b*d3*1000*f_s1
        Fs2=rho2*b*d3*1000*f_s2
        Fs3=rho1*b*d3*1000*f_s3
        
        # Compute Compression Force
        a1b1=(e_top/e_c)-((e_top/e_c)**2)/3
        b1=(4-(e_top/e_c))/(6-2*(e_top/e_c))
        Fc=a1b1*c*fc*b*1000
        
        # Section force
        Psec= P+Fs2+Fs3-Fc-Fs1
        
        # Adjust NA depth to balance section forces
        if Psec<0:
            #c0=c
            c=c-0.001
        elif Psec>0:
            #co=c
            c=c+0.001
        err=abs(Psec)
        if err<5:
            break
        count+=1
    # Compute the moment
    Mn=P*(0.5*h-c)+Fs1*(c-d1)+Fs2*(d2-c)+Fs3*(d3-c)+Fc*c*(1-b1/2)
    cn=c
    # print(Mn)
    # print(cn)
    return Mp,Mn,cp
   # return cn
#%% RUN MODEL

# Funcion ejecuta modelo OpenSees
# Input:
# rutaOpenSees: executable file address OpenSees.exe
# rutascript: script tcl file address 
def runModel(rutaOpenSees,rutascript):
    import subprocess

    cmd=[rutaOpenSees,rutascript]
    ret=subprocess.run(cmd,capture_output=True)
    print(ret.stderr.decode())

#%% Compute Mass

# Funcion que calcula la masa y cargas para el modelo OpenSees
# Input:
#rutaNudos: file nudos.txt
#rutaElementos: file elementos.txt

def computeMass(rutaNudos,rutaElementos,rutafolder):
    # CONSTANTS
    Tslab=0.145; # Heigth Slab
    Pe=2.4*9.81 #[kN/m3] # Peso especifico
    Bcol=0.30
    Hcol=0.20
    Bbeam=0.30
    Hbeam=0.20
    DeadLoad=3.2*1.1
    LiveLoad=2.0*1.1
    hlosa=0.145*1.1
    # FUNCTION
    import BibliotecaFunciones as bf
    import numpy as np
    import matplotlib.pyplot as plt
    import pandas as pd
    import string
    abecedario=string.ascii_uppercase
    m,mm,kN,N,Pa,kPa,MPa,g=bf.units()
    nudos=np.genfromtxt(rutaNudos) 
    tagNudos=nudos[:,0]
    x=nudos[:,1]
    y=nudos[:,2]
    z=nudos[:,3]
    elementos=r'D:\Google_Drive_Espe\001_Maestria_Ingenieria_Civil\005_Investigación\001_Ing_Poveda_Tesis\007_Modelo\001_Modelos_OpenSees\elementos.txt'
    elementos=np.genfromtxt(rutaElementos)
    #GRAFICO--------------------------
    cont=0
    Ni=np.empty((1,3),dtype=FloatingPointError)
    Nj=np.empty((1,3),dtype=FloatingPointError)
    elementType=[]
    for i in elementos:
        tag=int(i[0])
        #print('{:.0f}'.format(tag))
        #Nudo Inicial
        ni=int(i[1])
        #print('{:.0f}'.format(ni))
        posicion_ni=np.where(nudos==ni)
        #print(posicion_ni)
        coor_ni=nudos[posicion_ni[0]]
        #print(coor_ni)
        coorx_ni=coor_ni[0][1]
        #print("x:",coorx_ni)
        coory_ni=coor_ni[0][2]
        #print("y:",coory_ni)
        coorz_ni=coor_ni[0][3]
        #print("z:",coorz_ni)
        Ni=np.insert(Ni,cont,[coorx_ni,coory_ni,coorz_ni],axis=0)
        # Nudo Final
        nj=i[2]
        #print('{:.0f}'.format(nj))
        posicion_nj=np.where(nudos==nj)
        #print(posicion_nj)
        coor_nj=nudos[posicion_nj[0]]
        #print(coor_nj)
        coorx_nj=coor_nj[0][1]
        #print("x:",coorx_nj)
        coory_nj=coor_nj[0][2]
        #print("y:",coory_nj)
        coorz_nj=coor_nj[0][3]
        #print("z:",coorz_nj)
        Nj=np.insert(Nj,cont,[coorx_nj,coory_nj,coorz_nj],axis=0)   
        coorx=[coorx_ni,coorx_nj]
        coory=[coory_ni,coory_nj]
        coorz=[coorz_ni,coorz_nj]
        # ax.plot3D(coorx,coorz,coory,color='k') 
        if coory_ni-coory_nj!=0:
            # ax.text(coorx_ni,coorz_ni,coory_ni+((coory_nj-coory_ni)/2),tag,color='b') 
            elementType.append('Column')
        elif coorx_ni-coorx_nj!=0:
            # ax.text(coorx_ni+((coorx_nj-coorx_ni)/2),coorz_ni,coory_ni,tag,color='b')
            elementType.append('Beam')
        elif coorz_ni-coorz_nj!=0:
            # ax.text(coorx_ni,coorz_ni+((coorz_nj-coorz_ni)/2),coory_ni,tag,color='b')
            elementType.append('Gird')
        cont+=1
        #plt.show()
    # RESUMEN------------------------
    modelo_df=pd.DataFrame(elementos)
    modelo_df=modelo_df.rename(columns={0:'Element',1:'Ni',2:'Nj'})
    Ni_df=pd.DataFrame(Ni[0:-1])
    Ni_df=Ni_df.rename(columns={0:'Nix',1:'Niy',2:'Niz'})
    Nj_df=pd.DataFrame(Nj[0:-1])
    Nj_df=Nj_df.rename(columns={0:'Njx',1:'Njy',2:'Njz'})
    elementType_df=pd.DataFrame(elementType)
    elementType_df=elementType_df.rename(columns={0:'Type'})
    modelo=pd.concat([elementType_df,modelo_df,Ni_df,Nj_df],axis=1)
    modelo.loc[modelo.Type=='Column','Length']=modelo.Njy-modelo.Niy
    modelo.loc[modelo.Type=='Beam','Length']=modelo.Njx-modelo.Nix
    modelo.loc[modelo.Type=='Gird','Length']=modelo.Njz-modelo.Niz
    modelo.loc[modelo.Type=='Column','Base']=Bcol
    modelo.loc[modelo.Type=='Column','Height']=Hcol
    modelo.loc[modelo.Type=='Beam','Base']=Bbeam
    modelo.loc[modelo.Type=='Beam','Height']=Hbeam
    modelo.loc[modelo.Type=='Gird','Base']=Bbeam
    modelo.loc[modelo.Type=='Gird','Height']=Hbeam
    modelo['Weight']=modelo.Base*modelo.Height*modelo.Length*Pe #[kN]
    for i in range(0,len(modelo)):
        if modelo.loc[i,'Element']%100 >=10 and modelo.loc[i,'Element']%100 <20:
            modelo.loc[i,'Floor']=1
        elif modelo.loc[i,'Element']%100 >=20 and modelo.loc[i,'Element']%100 <30:
            modelo.loc[i,'Floor']=2 
        elif modelo.loc[i,'Element']%100 >=30 and modelo.loc[i,'Element']%100 <40:
            modelo.loc[i,'Floor']=3 
        elif modelo.loc[i,'Element']%100 >=40 and modelo.loc[i,'Element']%100 <50:
                modelo.loc[i,'Floor']=4 
    #print('Modelo Dataframe\n',modelo)
    # Mass---------------------
    mass=pd.DataFrame([]) #[ton]
    np=modelo.Floor.max()
    for i in range(1,int(np)+1):
        mass.loc[i,'Floor']=i
        mass.loc[i,'Lenght_x']=modelo.loc[modelo.Floor==i]['Nix'].max()  
        mass.loc[i,'Lenght_z']=modelo.loc[modelo.Floor==i]['Niz'].max()  
    mass['Floor_Area']=mass.Lenght_x*mass.Lenght_z
    MassFloor1=(modelo.loc[modelo.Floor==1].sum()['Weight'])/g
    MassFloor2=(modelo.loc[modelo.Floor==2].sum()['Weight'])/g
    MassFloor3=(modelo.loc[modelo.Floor==3].sum()['Weight'])/g
    mass['Mass_BeamColumn']=[MassFloor1,MassFloor2,MassFloor3]
    mass['Mass_Slab']=hlosa*mass.Floor_Area*Pe/g
    mass['Mass_Dead']=DeadLoad*mass.Floor_Area/g
    mass['Mass_Live']=LiveLoad*0.25*mass.Floor_Area/g
    mass['Mass']=mass.Mass_BeamColumn+mass.Mass_Slab+mass.Mass_Dead+mass.Mass_Live
    # Numero Columnas
    for i in range(1,int(np)+1):
        mass.loc[i,'Num_Columns']=modelo[(modelo.Type=='Column') & (modelo.Floor==i)]['Type'].count()
        #mass.loc[i,'Num_Beams']=modelo[(modelo.Type=='Beam') & (modelo.Floor==i)]['Type'].count()
        #mass.loc[i,'Num_Girds']=modelo[(modelo.Type=='Gird') & (modelo.Floor==i)]['Type'].count()
    mass['Nodal_Mass']=mass.Mass/mass.Num_Columns
    #print('Mass Dataframe\n',mass)
    # GENERACION ARCHIVO MASA NODAL
    massNodal=pd.DataFrame([])
    # mass Mx My Mz MRx MRy MRz
    for j in range(1,int(np)+2):
        for i in range(0,len(nudos)):
            massNodal.loc[i,'Mass']='mass'
            massNodal.loc[i,'Nudos']=nudos[i,0]
            if massNodal.loc[i,'Nudos']%100>=(j-1)*10 and massNodal.loc[i,'Nudos']%100<j*10:
                massNodal.loc[i,'Floor']=j-1
    for i in range(1,int(np+1)):
        massNodal.loc[massNodal.Floor==i,'Mx']=mass.loc[i,'Nodal_Mass']
        massNodal.loc[massNodal.Floor==i,'My']=0
        massNodal.loc[massNodal.Floor==i,'Mz']=mass.loc[i,'Nodal_Mass']
        massNodal.loc[massNodal.Floor==i,'MRx']=0
        massNodal.loc[massNodal.Floor==i,'MRy']=0
        massNodal.loc[massNodal.Floor==i,'MRz']=0
    massNodal=massNodal.loc[massNodal.Floor>=1]
    massNodal['Nudos']=massNodal['Nudos'].astype(int)
    massOpenSees=massNodal[['Mass','Nudos','Mx','My','Mz','MRx','MRy','MRz']]
    
    # LOADS-------------------------------------------------------------------
    # Eje X
    gridX=modelo[['Nix']]
    gridX=pd.unique(gridX.Nix)
    vanosX=[]
    for i in range(len(gridX)-1):
        vanosX.append(gridX[i+1]-gridX[i])
    #print(vanosX)
    # Alturas eje Y
    gridY=modelo[['Niy']]
    gridY=pd.unique(gridY.Niy)
    vanosY=[]
    for i in range(len(gridY)-1):
        vanosY.append(gridY[i+1]-gridY[i])
    #print(vanosY)
    
    # Eje Z
    gridZ=modelo[['Niz']]
    gridZ=pd.unique(gridZ.Niz)
    vanosZ=[]
    for i in range(len(gridZ)-1):
        vanosZ.append(gridZ[i+1]-gridZ[i])
    #print(vanosZ)
    
    # Frames
    load=pd.DataFrame([])
    load=modelo[['Type','Element','Length']].copy()
    FramesX=list(range(1,len(gridX)+1))
    FramesY=list(range(1,len(gridY)+1))
    for j in FramesY:   
        for i in range(0,len(modelo)):
            if load.loc[i,'Type']=='Column':
                if load.loc[i,'Element']//1000==j:
                        load.loc[i,'Frame']=j
            elif load.loc[i,'Type']=='Beam':
                if load.loc[i,'Element']//1000==j:
                        load.loc[i,'Frame']=j
                        
    for j in FramesX:
        for i in range(0,len(modelo)):
            if load.loc[i,'Type']=='Gird':
                if load.loc[i,'Element']//10000==j:
                        load.loc[i,'Frame']=abecedario[j-1]
    Zaux=vanosZ
    Zaux.insert(0,0)
    Zaux.append(0)
    #Load Combination
    slabLoad=hlosa*Pe
    Lc=DeadLoad+0.25*LiveLoad+slabLoad
    # Beams
    for i in FramesY:
        for j in range(len(load)):
            if load.loc[j,'Type']=='Beam' and load.loc[j,'Frame']==i:
                if Zaux[i-1]!=0 and Zaux[i]!=0:
                    if load.loc[j,'Length']>=Zaux[i-1]:
                        Lgird=Zaux[i-1]
                        s=min(load.loc[j,'Length'],Lgird)
                        m=s/max(load.loc[j,'Length'],Lgird)
                        w1=Lc*s/3*((3-m**2)/2)
                    elif load.loc[j,'Length']<Zaux[i-1]:
                        Lgird=Zaux[i-1]
                        s=min(load.loc[j,'Length'],Lgird)
                        w1=Lc*s/3
                        
                    if load.loc[j,'Length']>=Zaux[i]:
                        Lgird=Zaux[i]
                        s=min(load.loc[j,'Length'],Lgird)
                        m=s/max(load.loc[j,'Length'],Lgird)
                        w2=Lc*s/3*((3-m**2)/2)
                    elif load.loc[j,'Length']<Zaux[i]:
                        Lgird=Zaux[i]
                        s=min(load.loc[j,'Length'],Lgird)
                        w2=Lc*s/3
                    load.loc[j,'Qdl']=w1+w2
                else:
                     if load.loc[j,'Length']>=max(Zaux[i],Zaux[i-1]):
                        Lgird=max(Zaux[i],Zaux[i-1])
                        s=min(load.loc[j,'Length'],Lgird)
                        m=s/max(load.loc[j,'Length'],Lgird)
                        load.loc[j,'Qdl']=Lc*s/3*((3-m**2)/2)
                     elif load.loc[j,'Length']<max(Zaux[i],Zaux[i-1]):
                        Lgird=max(Zaux[i],Zaux[i-1])
                        s=min(load.loc[j,'Length'],Lgird)
                        load.loc[j,'Qdl']=Lc*s/3
                        
    Xaux=vanosX
    Xaux.insert(0,0)
    Xaux.append(0)                  
    # Girds
    for i in FramesX:
        for j in range(len(load)):
            if load.loc[j,'Type']=='Gird' and load.loc[j,'Frame']==abecedario[i-1]:
                if Xaux[i-1]!=0 and Xaux[i]!=0:
                    if load.loc[j,'Length']>=Xaux[i-1]:
                        Lgird=Xaux[i-1]
                        s=min(load.loc[j,'Length'],Lgird)
                        m=s/max(load.loc[j,'Length'],Lgird)
                        w1=Lc*s/3*((3-m**2)/2)
                    elif load.loc[j,'Length']<Xaux[i-1]:
                        Lgird=Xaux[i-1]
                        s=min(load.loc[j,'Length'],Lgird)
                        w1=Lc*s/3
                        
                    if load.loc[j,'Length']>=Xaux[i]:
                        Lgird=Xaux[i]
                        s=min(load.loc[j,'Length'],Lgird)
                        m=s/max(load.loc[j,'Length'],Lgird)
                        w2=Lc*s/3*((3-m**2)/2)
                    elif load.loc[j,'Length']<Xaux[i]:
                        Lgird=Xaux[i]
                        s=min(load.loc[j,'Length'],Lgird)
                        w2=Lc*s/3
                    load.loc[j,'Qdl']=w1+w2
                else:
                     if load.loc[j,'Length']>=max(Xaux[i],Xaux[i-1]):
                        Lgird=max(Xaux[i],Xaux[i-1])
                        s=min(load.loc[j,'Length'],Lgird)
                        m=s/max(load.loc[j,'Length'],Lgird)
                        load.loc[j,'Qdl']=Lc*s/3*((3-m**2)/2)
                     elif load.loc[j,'Length']<max(Xaux[i],Xaux[i-1]):
                        Lgird=max(Xaux[i],Xaux[i-1])
                        s=min(load.loc[j,'Length'],Lgird)
                        load.loc[j,'Qdl']=Lc*s/3
    # Columns  
    load.loc[load.Type=='Column','Qdl']=modelo.loc[modelo.Type=='Column']['Weight']/modelo.loc[modelo.Type=='Column']['Length']    
    #print(load)          
     #Loads OpenSees-----------------------------------------------------------
    loadOpensees=pd.DataFrame([]) 
    for i in range(len(load)): 
        loadOpensees.loc[i,'eleLoad']='eleLoad'
        loadOpensees.loc[i,'ele']='-ele'
        loadOpensees.loc[i,'Element']=load.loc[i,'Element'].copy()
        loadOpensees.loc[i,'Tipo']=load.loc[i,'Type']
        loadOpensees.loc[i,'type']='-type'
        loadOpensees.loc[i,'beamUniform']='-beamUniform'
        if loadOpensees.loc[i,'Tipo']=='Column':
            loadOpensees.loc[i,'Wy']=0.00
            loadOpensees.loc[i,'Wz']=0.00
            loadOpensees.loc[i,'Wx']=-load.loc[i,'Qdl']*0.90
        elif loadOpensees.loc[i,'Tipo']=='Beam' or loadOpensees.loc[i,'Tipo']=='Gird':
            loadOpensees.loc[i,'Wy']=(-load.loc[i,'Qdl']-modelo.loc[i,'Weight']/modelo.loc[i,'Length'])*0.90
            loadOpensees.loc[i,'Wz']=0.00
            loadOpensees.loc[i,'Wx']=0.00
            
    loadOpensees.drop('Tipo',axis=1,inplace=True)    
    loadOpensees['Element']=loadOpensees['Element'].astype(int)
    #print(loadOpensees) 
    return massOpenSees,loadOpensees,modelo,mass

#%% INELASTIC ELEMENT
def inelasticElement(rutafolder,modelo,rutaAxialLoad,fc,rho_Ash):
   # Inelastic Element-----------------------------------------------------
    import pandas as pd
    import numpy as np
    InelasticElement=pd.DataFrame()
    # rcBc_nonDuct Procedure
    fy=420
    Es=200e3
    #fc=21
    Ec=4700*np.sqrt(fc)
    sCol=0.15
    sBeam=0.15
    sGird=0.15
    Cover=0.03
    # Longitudinal Reinforcing
    dbL=0.012
    #Column
    # Nudo Inicial
    #Sentido zz
    NumBarrasTopCzz=3
    NumBarrasMidCzz=0
    NumBarrasBotCzz=3
    # Nudo Inicial
    NumBarrasTopC2zz=3
    NumBarrasMidC2zz=0
    NumBarrasBotC2zz=3
    # Sentido yy
    # Nudo Inicial
    NumBarrasTopCyy=2
    NumBarrasMidCyy=2
    NumBarrasBotCyy=2
    # Nudo Inicial
    NumBarrasTopC2yy=2
    NumBarrasMidC2yy=2
    NumBarrasBotC2yy=2
    
    #Beam
    # Sentido zz
    # Nudo Inicial
    NumBarrasTopBzz=4
    NumBarrasMidBzz=0
    NumBarrasBotBzz=4
    # Nudo Inicial
    NumBarrasTopB2zz=4
    NumBarrasMidB2zz=0
    NumBarrasBotB2zz=4
    
    # Sentido yy
    # Nudo Inicial
    NumBarrasTopByy=2
    NumBarrasMidByy=4
    NumBarrasBotByy=2
    # Nudo Inicial
    NumBarrasTopB2yy=2
    NumBarrasMidB2yy=4
    NumBarrasBotB2yy=2
    
    # Shear reinforcing
    dbV=0.008
    NumRamalesshear=2
    #rutaAxialLoad=rutafolder+'\\'+'DataElasticModel'+'\\'+'AxialLoad.out'
    AxialLoad=np.genfromtxt(rutaAxialLoad)
    AxialLoad=AxialLoad[9,1:]
    
    
    # Dataframe
    cont=0
    for i in range(len(modelo)):
        InelasticElement.loc[i,'Procedure']='rcBC_nonDuct'
        if modelo.loc[i,'Type']=='Column':
            InelasticElement.loc[i,'ST']=1
        else: 
            InelasticElement.loc[i,'ST']=0
            
        InelasticElement.loc[i,'Element']=modelo.loc[i,'Element'].copy()
        
        if modelo.loc[i,'Type']=='Column':
            InelasticElement.loc[i,'IDTransf']='$IDColTransf'
        elif modelo.loc[i,'Type']=='Beam':
            InelasticElement.loc[i,'IDTransf']='$IDBeamTransf'
        elif modelo.loc[i,'Type']=='Gird':
            InelasticElement.loc[i,'IDTransf']='$IDGirdTransf'
            
        InelasticElement.loc[i,'iNode']=modelo.loc[i,'Ni'].copy()
        InelasticElement.loc[i,'jNode']=modelo.loc[i,'Nj'].copy()
        InelasticElement.loc[i,'fyL']=fy
        InelasticElement.loc[i,'fy']=fy
        InelasticElement.loc[i,'Es']=Es
        InelasticElement.loc[i,'fc']=fc
        InelasticElement.loc[i,'Ec']=Ec
        InelasticElement.loc[i,'Base']=modelo.loc[i,'Base'].copy()
        InelasticElement.loc[i,'Height']=modelo.loc[i,'Height'].copy()
    
        if modelo.loc[i,'Type']=='Column':
            InelasticElement.loc[i,'spacing']=sCol
        elif modelo.loc[i,'Type']=='Beam':
            InelasticElement.loc[i,'spacing']=sBeam
        elif modelo.loc[i,'Type']=='Gird':
            InelasticElement.loc[i,'spacing']=sGird
        InelasticElement.loc[i,'Cover']=Cover
        InelasticElement.loc[i,'dbL']=dbL
        InelasticElement.loc[i,'dbV']=dbV
        
        if modelo.loc[i,'Type']=='Column':
            InelasticElement.loc[i,'Axial']=AxialLoad[cont]
            cont+=1
            InelasticElement.loc[i,'Ls']=modelo.loc[i,'Length'].copy()
        elif modelo.loc[i,'Type']=='Beam':
            InelasticElement.loc[i,'spacing']=sBeam
            InelasticElement.loc[i,'Axial']=0.0
            InelasticElement.loc[i,'Ls']=modelo.loc[i,'Length']/2
        elif modelo.loc[i,'Type']=='Gird':
            InelasticElement.loc[i,'spacing']=sGird
            InelasticElement.loc[i,'Axial']=0.0
            InelasticElement.loc[i,'Ls']=modelo.loc[i,'Length']/2
        #InelasticElement.loc[i,'rho_shr']=np.pi*(InelasticElement.loc[i,'dbV']**2)/4*NumRamalesshear/(InelasticElement.loc[i,'Base']*InelasticElement.loc[i,'spacing']) 
        InelasticElement.loc[i,'rho_shr']=rho_Ash
        if modelo.loc[i,'Type']=='Column':  
            #Nudo Inicial
            #rho_top1zz
             InelasticElement.loc[i,'rho_top1zz']=np.pi*(InelasticElement.loc[i,'dbL']**2)/4*NumBarrasTopCzz/(InelasticElement.loc[i,'Base']*
            (InelasticElement.loc[i,'Height']-InelasticElement.loc[i,'Cover']-InelasticElement.loc[i,'dbV']-InelasticElement.loc[i,'dbL']/2)) 
             #rho_mid1zz
             InelasticElement.loc[i,'rho_mid1zz']=np.pi*(InelasticElement.loc[i,'dbL']**2)/4*NumBarrasMidCzz/(InelasticElement.loc[i,'Base']*
            (InelasticElement.loc[i,'Height']-InelasticElement.loc[i,'Cover']-InelasticElement.loc[i,'dbV']-InelasticElement.loc[i,'dbL']/2)) 
             #rho_bot1zz
             InelasticElement.loc[i,'rho_bot1zz']=np.pi*(InelasticElement.loc[i,'dbL']**2)/4*NumBarrasBotCzz/(InelasticElement.loc[i,'Base']*
            (InelasticElement.loc[i,'Height']-InelasticElement.loc[i,'Cover']-InelasticElement.loc[i,'dbV']-InelasticElement.loc[i,'dbL']/2)) 
            #Nudo Final
            #rho_top2zz
             InelasticElement.loc[i,'rho_top2zz']=np.pi*(InelasticElement.loc[i,'dbL']**2)/4*NumBarrasTopC2zz/(InelasticElement.loc[i,'Base']*
            (InelasticElement.loc[i,'Height']-InelasticElement.loc[i,'Cover']-InelasticElement.loc[i,'dbV']-InelasticElement.loc[i,'dbL']/2)) 
             #rho_mid2zz
             InelasticElement.loc[i,'rho_mid2zz']=np.pi*(InelasticElement.loc[i,'dbL']**2)/4*NumBarrasMidC2zz/(InelasticElement.loc[i,'Base']*
            (InelasticElement.loc[i,'Height']-InelasticElement.loc[i,'Cover']-InelasticElement.loc[i,'dbV']-InelasticElement.loc[i,'dbL']/2)) 
             #rho_bot2zz
             InelasticElement.loc[i,'rho_bot2zz']=np.pi*(InelasticElement.loc[i,'dbL']**2)/4*NumBarrasBotC2zz/(InelasticElement.loc[i,'Base']*
            (InelasticElement.loc[i,'Height']-InelasticElement.loc[i,'Cover']-InelasticElement.loc[i,'dbV']-InelasticElement.loc[i,'dbL']/2)) 
             
             #Sentido yy
             #Nudo Inicial
            #rho_top1zz
             InelasticElement.loc[i,'rho_top1yy']=np.pi*(InelasticElement.loc[i,'dbL']**2)/4*NumBarrasTopCyy/(InelasticElement.loc[i,'Height']*
            (InelasticElement.loc[i,'Base']-InelasticElement.loc[i,'Cover']-InelasticElement.loc[i,'dbV']-InelasticElement.loc[i,'dbL']/2)) 
             #rho_mid1zz
             InelasticElement.loc[i,'rho_mid1yy']=np.pi*(InelasticElement.loc[i,'dbL']**2)/4*NumBarrasMidCyy/(InelasticElement.loc[i,'Height']*
            (InelasticElement.loc[i,'Base']-InelasticElement.loc[i,'Cover']-InelasticElement.loc[i,'dbV']-InelasticElement.loc[i,'dbL']/2)) 
             #rho_bot1zz
             InelasticElement.loc[i,'rho_bot1yy']=np.pi*(InelasticElement.loc[i,'dbL']**2)/4*NumBarrasBotCyy/(InelasticElement.loc[i,'Height']*
            (InelasticElement.loc[i,'Base']-InelasticElement.loc[i,'Cover']-InelasticElement.loc[i,'dbV']-InelasticElement.loc[i,'dbL']/2)) 
            #Nudo Final
            #rho_top2zz
             InelasticElement.loc[i,'rho_top2yy']=np.pi*(InelasticElement.loc[i,'dbL']**2)/4*NumBarrasTopC2yy/(InelasticElement.loc[i,'Height']*
            (InelasticElement.loc[i,'Base']-InelasticElement.loc[i,'Cover']-InelasticElement.loc[i,'dbV']-InelasticElement.loc[i,'dbL']/2)) 
             #rho_mid2zz
             InelasticElement.loc[i,'rho_mid2yy']=np.pi*(InelasticElement.loc[i,'dbL']**2)/4*NumBarrasMidC2yy/(InelasticElement.loc[i,'Height']*
            (InelasticElement.loc[i,'Base']-InelasticElement.loc[i,'Cover']-InelasticElement.loc[i,'dbV']-InelasticElement.loc[i,'dbL']/2)) 
             #rho_bot2zz
             InelasticElement.loc[i,'rho_bot2yy']=np.pi*(InelasticElement.loc[i,'dbL']**2)/4*NumBarrasBotC2yy/(InelasticElement.loc[i,'Height']*
            (InelasticElement.loc[i,'Base']-InelasticElement.loc[i,'Cover']-InelasticElement.loc[i,'dbV']-InelasticElement.loc[i,'dbL']/2)) 
             
        if modelo.loc[i,'Type']=='Beam': 
            #Nudo Inicial
            #rho_top1zz
             InelasticElement.loc[i,'rho_top1zz']=np.pi*(InelasticElement.loc[i,'dbL']**2)/4*NumBarrasTopBzz/(InelasticElement.loc[i,'Base']*
            (InelasticElement.loc[i,'Height']-InelasticElement.loc[i,'Cover']-InelasticElement.loc[i,'dbV']-InelasticElement.loc[i,'dbL']/2)) 
             #rho_mid1zz
             InelasticElement.loc[i,'rho_mid1zz']=np.pi*(InelasticElement.loc[i,'dbL']**2)/4*NumBarrasMidBzz/(InelasticElement.loc[i,'Base']*
            (InelasticElement.loc[i,'Height']-InelasticElement.loc[i,'Cover']-InelasticElement.loc[i,'dbV']-InelasticElement.loc[i,'dbL']/2)) 
             #rho_bot1zz
             InelasticElement.loc[i,'rho_bot1zz']=np.pi*(InelasticElement.loc[i,'dbL']**2)/4*NumBarrasBotBzz/(InelasticElement.loc[i,'Base']*
            (InelasticElement.loc[i,'Height']-InelasticElement.loc[i,'Cover']-InelasticElement.loc[i,'dbV']-InelasticElement.loc[i,'dbL']/2)) 
             #Nudo Final
            #rho_top2zz
             InelasticElement.loc[i,'rho_top2zz']=np.pi*(InelasticElement.loc[i,'dbL']**2)/4*NumBarrasTopB2zz/(InelasticElement.loc[i,'Base']*
            (InelasticElement.loc[i,'Height']-InelasticElement.loc[i,'Cover']-InelasticElement.loc[i,'dbV']-InelasticElement.loc[i,'dbL']/2)) 
             #rho_mid2zz
             InelasticElement.loc[i,'rho_mid2zz']=np.pi*(InelasticElement.loc[i,'dbL']**2)/4*NumBarrasMidB2zz/(InelasticElement.loc[i,'Base']*
            (InelasticElement.loc[i,'Height']-InelasticElement.loc[i,'Cover']-InelasticElement.loc[i,'dbV']-InelasticElement.loc[i,'dbL']/2)) 
             #rho_bot2zz
             InelasticElement.loc[i,'rho_bot2zz']=np.pi*(InelasticElement.loc[i,'dbL']**2)/4*NumBarrasBotB2zz/(InelasticElement.loc[i,'Base']*
            (InelasticElement.loc[i,'Height']-InelasticElement.loc[i,'Cover']-InelasticElement.loc[i,'dbV']-InelasticElement.loc[i,'dbL']/2)) 
             # Sentido yy
             #Nudo Inicial
            #rho_top1yy
             InelasticElement.loc[i,'rho_top1yy']=np.pi*(InelasticElement.loc[i,'dbL']**2)/4*NumBarrasTopByy/(InelasticElement.loc[i,'Height']*
            (InelasticElement.loc[i,'Base']-InelasticElement.loc[i,'Cover']-InelasticElement.loc[i,'dbV']-InelasticElement.loc[i,'dbL']/2)) 
             #rho_mid1yy
             InelasticElement.loc[i,'rho_mid1yy']=np.pi*(InelasticElement.loc[i,'dbL']**2)/4*NumBarrasMidByy/(InelasticElement.loc[i,'Height']*
            (InelasticElement.loc[i,'Base']-InelasticElement.loc[i,'Cover']-InelasticElement.loc[i,'dbV']-InelasticElement.loc[i,'dbL']/2)) 
             #rho_bot1yy
             InelasticElement.loc[i,'rho_bot1yy']=np.pi*(InelasticElement.loc[i,'dbL']**2)/4*NumBarrasBotByy/(InelasticElement.loc[i,'Height']*
            (InelasticElement.loc[i,'Base']-InelasticElement.loc[i,'Cover']-InelasticElement.loc[i,'dbV']-InelasticElement.loc[i,'dbL']/2)) 
             #Nudo Final
            #rho_top2yy
             InelasticElement.loc[i,'rho_top2yy']=np.pi*(InelasticElement.loc[i,'dbL']**2)/4*NumBarrasTopB2yy/(InelasticElement.loc[i,'Height']*
            (InelasticElement.loc[i,'Base']-InelasticElement.loc[i,'Cover']-InelasticElement.loc[i,'dbV']-InelasticElement.loc[i,'dbL']/2)) 
             #rho_mid2yy
             InelasticElement.loc[i,'rho_mid2yy']=np.pi*(InelasticElement.loc[i,'dbL']**2)/4*NumBarrasMidB2yy/(InelasticElement.loc[i,'Height']*
            (InelasticElement.loc[i,'Base']-InelasticElement.loc[i,'Cover']-InelasticElement.loc[i,'dbV']-InelasticElement.loc[i,'dbL']/2)) 
             #rho_bot2yy
             InelasticElement.loc[i,'rho_bot2yy']=np.pi*(InelasticElement.loc[i,'dbL']**2)/4*NumBarrasBotB2yy/(InelasticElement.loc[i,'Height']*
            (InelasticElement.loc[i,'Base']-InelasticElement.loc[i,'Cover']-InelasticElement.loc[i,'dbV']-InelasticElement.loc[i,'dbL']/2)) 
          
        if modelo.loc[i,'Type']=='Gird': 
            #Nudo Inicial
            #rho_top1zz
             InelasticElement.loc[i,'rho_top1zz']=np.pi*(InelasticElement.loc[i,'dbL']**2)/4*NumBarrasTopBzz/(InelasticElement.loc[i,'Base']*
            (InelasticElement.loc[i,'Height']-InelasticElement.loc[i,'Cover']-InelasticElement.loc[i,'dbV']-InelasticElement.loc[i,'dbL']/2)) 
             #rho_mid1zz
             InelasticElement.loc[i,'rho_mid1zz']=np.pi*(InelasticElement.loc[i,'dbL']**2)/4*NumBarrasMidBzz/(InelasticElement.loc[i,'Base']*
            (InelasticElement.loc[i,'Height']-InelasticElement.loc[i,'Cover']-InelasticElement.loc[i,'dbV']-InelasticElement.loc[i,'dbL']/2)) 
             #rho_bot1zz
             InelasticElement.loc[i,'rho_bot1zz']=np.pi*(InelasticElement.loc[i,'dbL']**2)/4*NumBarrasBotBzz/(InelasticElement.loc[i,'Base']*
            (InelasticElement.loc[i,'Height']-InelasticElement.loc[i,'Cover']-InelasticElement.loc[i,'dbV']-InelasticElement.loc[i,'dbL']/2)) 
            #Nudo Final
            #rho_top2zz
             InelasticElement.loc[i,'rho_top2zz']=np.pi*(InelasticElement.loc[i,'dbL']**2)/4*NumBarrasTopB2zz/(InelasticElement.loc[i,'Base']*
            (InelasticElement.loc[i,'Height']-InelasticElement.loc[i,'Cover']-InelasticElement.loc[i,'dbV']-InelasticElement.loc[i,'dbL']/2)) 
             #rho_mid2zz
             InelasticElement.loc[i,'rho_mid2zz']=np.pi*(InelasticElement.loc[i,'dbL']**2)/4*NumBarrasMidB2zz/(InelasticElement.loc[i,'Base']*
            (InelasticElement.loc[i,'Height']-InelasticElement.loc[i,'Cover']-InelasticElement.loc[i,'dbV']-InelasticElement.loc[i,'dbL']/2)) 
             #rho_bot2zz
             InelasticElement.loc[i,'rho_bot2zz']=np.pi*(InelasticElement.loc[i,'dbL']**2)/4*NumBarrasBotB2zz/(InelasticElement.loc[i,'Base']*
            (InelasticElement.loc[i,'Height']-InelasticElement.loc[i,'Cover']-InelasticElement.loc[i,'dbV']-InelasticElement.loc[i,'dbL']/2)) 
            # Sentido yy
             #Nudo Inicial
            #rho_top1yy
             InelasticElement.loc[i,'rho_top1yy']=np.pi*(InelasticElement.loc[i,'dbL']**2)/4*NumBarrasTopByy/(InelasticElement.loc[i,'Height']*
            (InelasticElement.loc[i,'Base']-InelasticElement.loc[i,'Cover']-InelasticElement.loc[i,'dbV']-InelasticElement.loc[i,'dbL']/2)) 
             #rho_mid1yy
             InelasticElement.loc[i,'rho_mid1yy']=np.pi*(InelasticElement.loc[i,'dbL']**2)/4*NumBarrasMidByy/(InelasticElement.loc[i,'Height']*
            (InelasticElement.loc[i,'Base']-InelasticElement.loc[i,'Cover']-InelasticElement.loc[i,'dbV']-InelasticElement.loc[i,'dbL']/2)) 
             #rho_bot1yy
             InelasticElement.loc[i,'rho_bot1yy']=np.pi*(InelasticElement.loc[i,'dbL']**2)/4*NumBarrasBotByy/(InelasticElement.loc[i,'Height']*
            (InelasticElement.loc[i,'Base']-InelasticElement.loc[i,'Cover']-InelasticElement.loc[i,'dbV']-InelasticElement.loc[i,'dbL']/2)) 
             #Nudo Final
            #rho_top2yy
             InelasticElement.loc[i,'rho_top2yy']=np.pi*(InelasticElement.loc[i,'dbL']**2)/4*NumBarrasTopB2yy/(InelasticElement.loc[i,'Height']*
            (InelasticElement.loc[i,'Base']-InelasticElement.loc[i,'Cover']-InelasticElement.loc[i,'dbV']-InelasticElement.loc[i,'dbL']/2)) 
             #rho_mid2yy
             InelasticElement.loc[i,'rho_mid2yy']=np.pi*(InelasticElement.loc[i,'dbL']**2)/4*NumBarrasMidB2yy/(InelasticElement.loc[i,'Height']*
            (InelasticElement.loc[i,'Base']-InelasticElement.loc[i,'Cover']-InelasticElement.loc[i,'dbV']-InelasticElement.loc[i,'dbL']/2)) 
             #rho_bot2yy
             InelasticElement.loc[i,'rho_bot2yy']=np.pi*(InelasticElement.loc[i,'dbL']**2)/4*NumBarrasBotB2yy/(InelasticElement.loc[i,'Height']*
            (InelasticElement.loc[i,'Base']-InelasticElement.loc[i,'Cover']-InelasticElement.loc[i,'dbV']-InelasticElement.loc[i,'dbL']/2)) 
             
        InelasticElement.loc[i,'shearfolder']='$shearfolder'
        
        if modelo.loc[i,'Type']=='Column':
            InelasticElement.loc[i,'pfile']='$pfile_cols'
        elif modelo.loc[i,'Type']=='Beam':
            InelasticElement.loc[i,'pfile']='$pfile_bms'
        elif modelo.loc[i,'Type']=='Gird':
            InelasticElement.loc[i,'pfile']='$pfile_bms'        
    
    InelasticElement['jNode']=InelasticElement['jNode'].astype(int)    
    InelasticElement['iNode']=InelasticElement['iNode'].astype(int)
    InelasticElement['ST']=InelasticElement['ST'].astype(int)
    InelasticElement['Element']=InelasticElement['Element'].astype(int)
    #print('InemalasticElement dataFrame\n',InelasticElement)
    return InelasticElement
#%% SHAPE LOAD FIRST MODE
def shapeLoad1Mode(rutaFirstMode):
    import numpy as np
    import pandas as pd
    #rutaFirstMode=rutafolder+'\\'+'DataElasticModel'+'\\'+'Vectors'+'\\'+'eigenvector1.out'
    firstMode=np.genfromtxt(rutaFirstMode)
    firstMode=firstMode[9,1:]
    firstMode=firstMode[0:-1:6]
    shapeLoad=pd.DataFrame()
    nodes=[1121,1131,1141]
    for i in range(3):
        shapeLoad.loc[i,'Load']='load'
        shapeLoad.loc[i,'node']=nodes[i]
        shapeLoad.loc[i,'Fx']=firstMode[i]
        shapeLoad.loc[i,'Fy']=0.0
        shapeLoad.loc[i,'Fz']=0.0
        shapeLoad.loc[i,'Mx']=0.0
        shapeLoad.loc[i,'My']=0.0
        shapeLoad.loc[i,'Mz']=0.0
    shapeLoad['node']=shapeLoad['node'].astype(int)
    return shapeLoad
#%% BARE FRAME 
def bareFrame(rutaOpenSees,rutafolder,rutaAxialLoad,rutaFirstMode,nameFile,nameFileInelastic,cont,x,fc,rho_Ash):
    import BibliotecaFunciones as bf
    # ELASTIC MODEL
    Star1='''
    # Elastic Model
    # Mauricio Guamán
    # Maestría en Ing. Civil-ESPE
    wipe
    # Start script
    model BasicBuilder -ndm 3 -ndf 6;
    '''
    DataElasticModel='set dataDir DataElasticBareModel'+str(cont)+';'
    Star2='''
    file mkdir $dataDir;
    logFile "$dataDir/Log_model.log"; 
    '''
    Star=Star1+DataElasticModel+Star2
    Libraries='''
    # Libraries
    source Units.tcl; 
    source rcBC_nonDuct.tcl; 
    source DisplayPlane.tcl; 
    source DisplayModel2D.tcl; 
    source DisplayModel3D.tcl; 
    '''
    Geometry1='''
    #Define Geometry
    # Grids
    set x1 0.00; 
    '''
    x2='set x2 '+ str(x)+';'
    x3='set x3 '+ str(x+x)+';'
    
    Geometry2='''
    set y1 0.00;
    set y2 2.50;
    set y3 5.00;
    set y4 7.50;
    set z1 0.00;
    set z2 3.50;
    set z3 7.50;
    set z4 11.0;
    # Defina Nodal Coordinate
    # node $tag_node $X  $Y  $Z
    node 101 $x1 $y1 $z1
    node 102 $x2 $y1 $z1
    node 103 $x3 $y1 $z1
    node 111 $x1 $y2 $z1
    node 112 $x2 $y2 $z1
    node 113 $x3 $y2 $z1
    node 121 $x1 $y3 $z1
    node 122 $x2 $y3 $z1
    node 123 $x3 $y3 $z1
    node 131 $x1 $y4 $z1
    node 132 $x2 $y4 $z1
    node 133 $x3 $y4 $z1
    
    node 201 $x1 $y1 $z2
    node 202 $x2 $y1 $z2
    node 203 $x3 $y1 $z2
    node 211 $x1 $y2 $z2
    node 212 $x2 $y2 $z2
    node 213 $x3 $y2 $z2
    node 221 $x1 $y3 $z2
    node 222 $x2 $y3 $z2
    node 223 $x3 $y3 $z2
    node 231 $x1 $y4 $z2
    node 232 $x2 $y4 $z2
    node 233 $x3 $y4 $z2
    
    node 301 $x1 $y1 $z3
    node 302 $x2 $y1 $z3
    node 303 $x3 $y1 $z3
    node 311 $x1 $y2 $z3
    node 312 $x2 $y2 $z3
    node 313 $x3 $y2 $z3
    node 321 $x1 $y3 $z3
    node 322 $x2 $y3 $z3
    node 323 $x3 $y3 $z3
    node 331 $x1 $y4 $z3
    node 332 $x2 $y4 $z3
    node 333 $x3 $y4 $z3
    
    node 401 $x1 $y1 $z4
    node 402 $x2 $y1 $z4
    node 403 $x3 $y1 $z4
    node 411 $x1 $y2 $z4
    node 412 $x2 $y2 $z4
    node 413 $x3 $y2 $z4
    node 421 $x1 $y3 $z4
    node 422 $x2 $y3 $z4
    node 423 $x3 $y3 $z4
    
    # Variables de nodos a emplear
    set Nudos [getNodeTags];
    set l [llength $Nudos]
    set NodeEnd [lindex $Nudos end];      # Se obtiene el número total de nodos
    set nudos [open nudos.txt w]
    foreach i $Nudos {
        set a [nodeCoord $i]
        puts $nudos [format "%s %s" $i $a]
    }
    close $nudos
    #-------------------Rigid diaphragm nodes
    set RigidDiaphragm ON;     # options: ON, OFF. specify this before the analysis parameters are set the constraints are handled differently.                          
    # Coordenadas del diafragma rigido
    #node $tag_node $X  $Y  $Z
    node 1121   [expr $x3/2]    $y2    [expr $z4/2];      # master nodes for rigid diaphragm -- story 2, bay 1, frame 1-2
    node 1131   [expr $x3/2]    $y3    [expr $z4/2];      # master nodes for rigid diaphragm -- story 3, bay 1, frame 1-2
    node 1141   [expr $x3/2]    $y4    [expr $z3/2];      # master nodes for rigid diaphragm -- story 4, bay 1, frame 1-2
    # Constraints for rigid diaphragm master nodes
    #fix $tag_node  #fix_X #fix_Y #fix_Z #fix_RX #fix_RY #fix_RZ 
    fix 1121 0  1  0  1  0  1
    fix 1131 0  1  0  1  0  1
    fix 1141 0  1  0  1  0  1
    # ------------------------define Rigid Diaphram, dof 2 is normal to floor
    set perpDirn 2;
    rigidDiaphragm $perpDirn 1121 111 112 113 211 212 213 311 312 313 411 412 413;  # level 2
    rigidDiaphragm $perpDirn 1131 121 122 123 221 222 223 321 322 323 421 422 423;  # level 3 
    rigidDiaphragm $perpDirn 1141 131 132 133 231 232 233 331 332 333;  # level 3 
    
    #-------------------------Support
    fixY 0.0 1 1 1 1 1 1;
    '''
    Section='''
    # Section-------------------------------------------------------------------------------------------------------------
    # Define section tags
    set ColSecTag 1
    set BeamSecTag 2
    set GirdSecTag 3
    # Section Properties:
    set HCol [expr 200*$mm];     # square-Column width
    set BCol [expr 300*$mm]
    set HBeam [expr 200*$mm];        # Beam depth -- perpendicular to bending axis
    set BBeam [expr 300*$mm];        # Beam width -- parallel to bending axis
    set HGird [expr 200*$mm];        # Girder depth -- perpendicular to bending axis
    set BGird [expr 300*$mm];        # Girder width -- parallel to bending axis
    '''
    Material1='''
    # material properties:
    '''    
    s_fc='set fc [expr '+str(fc)+'*$MPa];           # concrete nominal compressive strength'
    
    Material2='''
    set Ec [expr 4700*sqrt($fc/1000)*$MPa];    # concrete Young's Modulus
    set nu 0.2;         # Poisson's ratio
    set Gc [expr $Ec/2./[expr 1+$nu]];      # Torsional stiffness Modulus
    set Ubig 1e10
    set J $Ubig;            # set large torsional stiffness
    # column section properties:
    set AgCol [expr $HCol*$BCol];       # rectuangular-Column cross-sectional area
    set IzCol [expr 1*1./12*$BCol*pow($HCol,3)];  # about-local-z Rect-Column gross moment of inertial
    set IyCol [expr 1*1./12*$HCol*pow($BCol,3)];  # about-local-z Rect-Column gross moment of inertial
    # beam sections:
    set AgBeam [expr $HBeam*$BBeam];        # rectuangular-Beam cross-sectional area
    set IzBeam [expr 1*1./12*$BBeam*pow($HBeam,3)];   # about-local-z Rect-Beam cracked moment of inertial
    set IyBeam [expr 1*1./12*$HBeam*pow($BBeam,3)];   # about-local-y Rect-Beam cracked moment of inertial
    # girder sections:
    set AgGird [expr $HGird*$BGird];        # rectuangular-Girder cross-sectional area
    set IzGird [expr 1*1./12*$BGird*pow($HGird,3)];   # about-local-z Rect-Girder cracked moment of inertial
    set IyGird [expr 1*1./12*$HGird*pow($BGird,3)];   # about-local-y Rect-Girder cracked moment of inertial   
    section Elastic $ColSecTag $Ec $AgCol $IzCol $IyCol $Gc $J
    section Elastic $BeamSecTag $Ec $AgBeam $IzBeam $IyBeam $Gc $J
    section Elastic $GirdSecTag $Ec $AgGird $IzGird $IyGird $Gc $J
    '''
    Material=Material1+s_fc+Material2
    
    Element='''
    #-------------------------Elements
    set IDColTransf 1; # all columns
    set IDBeamTransf 2; # all beams
    set IDGirdTransf 3; # all girds
    
    set ColTransfType Linear ;      # options for columns: Linear PDelta  Corotational 
    geomTransf $ColTransfType  $IDColTransf  0 0 1;         # orientation of column stiffness affects bidirectional response.
    geomTransf Linear $IDBeamTransf 0 0 1
    geomTransf Linear $IDGirdTransf -1 0 0
    
    # element connectivity
    #Frame 1
    #Columns
    #element elasticBeamColumn $eleTag $iNode $jNode $A $E $G $J $Iy $Iz $transfTag <-mass $massDens> <-cMass>
    #Level 1
    element elasticBeamColumn  1111 101 111 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    element elasticBeamColumn  1112 102 112 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    element elasticBeamColumn  1113 103 113 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    #Level 2
    element elasticBeamColumn  1121 111 121 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    element elasticBeamColumn  1122 112 122 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    element elasticBeamColumn  1123 113 123 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    #Level 3
    element elasticBeamColumn  1131 121 131 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    element elasticBeamColumn  1132 122 132 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    element elasticBeamColumn  1133 123 133 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    
    #Beams
    #Level 1
    element elasticBeamColumn  1211 111 112 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDBeamTransf
    element elasticBeamColumn  1212 112 113 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDBeamTransf
    #Level 2
    element elasticBeamColumn  1221 121 122 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDBeamTransf
    element elasticBeamColumn  1222 122 123 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDBeamTransf
    #Level 3
    element elasticBeamColumn  1231 131 132 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDBeamTransf
    element elasticBeamColumn  1232 132 133 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDBeamTransf
    
    #Frame 2
    #Columns
    #Level 1
    element elasticBeamColumn  2111 201 211 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    element elasticBeamColumn  2112 202 212 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    element elasticBeamColumn  2113 203 213 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    #Level 2
    element elasticBeamColumn  2121 211 221 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    element elasticBeamColumn  2122 212 222 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    element elasticBeamColumn  2123 213 223 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    #Level 3
    element elasticBeamColumn  2131 221 231 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    element elasticBeamColumn  2132 222 232 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    element elasticBeamColumn  2133 223 233 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    
    #Beams
    #Level 1
    element elasticBeamColumn  2211 211 212 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDBeamTransf
    element elasticBeamColumn  2212 212 213 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDBeamTransf
    #Level 2
    element elasticBeamColumn  2221 221 222 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDBeamTransf
    element elasticBeamColumn  2222 222 223 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDBeamTransf
    #Level 3
    element elasticBeamColumn  2231 231 232 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDBeamTransf
    element elasticBeamColumn  2232 232 233 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDBeamTransf
    
    #Frame 3
    #Columns
    #Level 1
    element elasticBeamColumn  3111 301 311 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    element elasticBeamColumn  3112 302 312 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    element elasticBeamColumn  3113 303 313 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    #Level 2
    element elasticBeamColumn  3121 311 321 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    element elasticBeamColumn  3122 312 322 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    element elasticBeamColumn  3123 313 323 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    #Level 3
    element elasticBeamColumn  3131 321 331 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    element elasticBeamColumn  3132 322 332 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    element elasticBeamColumn  3133 323 333 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    
    #Beams
    #Level 1
    element elasticBeamColumn  3211 311 312 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDBeamTransf
    element elasticBeamColumn  3212 312 313 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDBeamTransf
    #Level 2
    element elasticBeamColumn  3221 321 322 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDBeamTransf
    element elasticBeamColumn  3222 322 323 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDBeamTransf
    #Level 3
    element elasticBeamColumn  3231 331 332 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDBeamTransf
    element elasticBeamColumn  3232 332 333 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDBeamTransf
    
    #Frame 4
    #Columns
    #Level 1
    element elasticBeamColumn  4111 401 411 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    element elasticBeamColumn  4112 402 412 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    element elasticBeamColumn  4113 403 413 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    #Level 2
    element elasticBeamColumn  4121 411 421 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    element elasticBeamColumn  4122 412 422 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    element elasticBeamColumn  4123 413 423 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    
    #Beams
    #Level 1
    element elasticBeamColumn  4211 411 412 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDBeamTransf
    element elasticBeamColumn  4212 412 413 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDBeamTransf
    #Level 2
    element elasticBeamColumn  4221 421 422 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDBeamTransf
    element elasticBeamColumn  4222 422 423 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDBeamTransf
    
    #Girders
    #Level 1
    element elasticBeamColumn  11311 111 211 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    element elasticBeamColumn  11312 211 311 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    element elasticBeamColumn  11313 311 411 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    element elasticBeamColumn  21314 112 212 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    element elasticBeamColumn  21315 212 312 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    element elasticBeamColumn  21316 312 412 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    element elasticBeamColumn  31317 113 213 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    element elasticBeamColumn  31318 213 313 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    element elasticBeamColumn  31319 313 413 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    
    #Level 2
    element elasticBeamColumn  11321 121 221 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    element elasticBeamColumn  11322 221 321 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    element elasticBeamColumn  11323 321 421 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    element elasticBeamColumn  21324 122 222 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    element elasticBeamColumn  21325 222 322 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    element elasticBeamColumn  21326 322 422 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    element elasticBeamColumn  31327 123 223 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    element elasticBeamColumn  31328 223 323 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    element elasticBeamColumn  31329 323 423 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    
    #Level 3
    element elasticBeamColumn  11331 131 231 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    element elasticBeamColumn  11332 231 331 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    element elasticBeamColumn  21334 132 232 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    element elasticBeamColumn  21335 232 332 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    element elasticBeamColumn  31337 133 233 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    element elasticBeamColumn  31338 233 333 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    
    set TagElem [getEleTags]
    set elementos [open elementos.txt w]
    foreach i $TagElem {
        set b [eleNodes $i]
        puts $elementos [format "%s %s" $i $b]
    }
    close $elementos
    
    # Mass
    # mass Node Mx My Mz MRx MRy MRz
    '''
    # tcl file generation with basic information
    file=open(nameFile,'w')
    file.write(Star)
    file.write(Libraries)
    file.write(Geometry1+x2+x3+Geometry2)
    file.write(Section)
    file.write(Material)
    file.write(Element)
    file.write('puts "Created: Nudos.txt and Elementos.txt " ')
    file.close()
    
    #Run OpenSees to get nodes and elements. 
    # FileAdress
    rutascript=rutafolder+'\\'+nameFile
    bf.runModel(rutaOpenSees,rutascript) 
    
    # MASS
    rutaNudos=rutafolder+'\\'+'nudos.txt'
    rutaElementos=rutafolder+'\\'+'elementos.txt'
    massOpenSees,loadOpensees,modelo,mass=bf.computeMass(rutaNudos,rutaElementos,rutafolder)
    ModalAnalysis='''
    # Modal Analysis------------------------------------------------
    set Inf Inf;                        # Inf = Si queremos que muestre información en pantalla 
                                        # NoInf = Si no queremos que muestre información en pantalla
    set Outputs $dataDir;
    set numModes 3;                     # Número de modos a evaluar
    # Outputs: nombre de la carpeta Outputs o "NoOutputs" si no se quieren resultados
    remove recorders 
    wipeAnalysis
    if {$Inf != "NoInf"} {
        puts "------------------------------------------------------"
        puts "MODAL ANALYSIS"
    }
    system UmfPack
    constraints Transformation
    set lambda [eigen -fullGenLapack $numModes]; # -fullGenLapack in case of [nº modes=ngdl and masses]
    set omega {}
    set f {}
    set T {}
    set pi [expr acos(-1.0)];
    
    foreach lam $lambda {
        lappend omega [expr sqrt($lam)]
        lappend F [expr sqrt($lam)/(2*$pi)]
        lappend T [expr (2*$pi)/sqrt($lam)]
    }
    
    if {$Inf != "NoInf"} {
        puts "Periods are: $T -s-"
        puts "Frequencies are: $F -Hz-"
    }
    
    if {$Outputs != "NoOutputs"} {
        file mkdir $Outputs/Modes
        file mkdir $Outputs/Vectors
        for { set k 1 } { $k <= $numModes } { incr k } {
        recorder Node -file [format "$dataDir/Modes/mode%i.out" $k] -time -nodeRange 1 "$::NodeEnd" -dof 1 2 3 4 5 6  "eigen $k"
        recorder Node -file [format "$dataDir/Vectors/eigenvector%i.out" $k] -time -node 1121 1131 1141 -dof 1 2 3 4 5 6  "eigen $k"
        }
        set period "$dataDir/Modes/Periods.txt"
        set Periods [open $period "w"]
        foreach t $T f $F w $omega {
            puts $Periods "$t  $f  $w"
        }
        close $Periods
        record
    }       
    '''
    recorderElastic='''
    recorder Element -file $dataDir/AxialLoad.out -time -ele 1111 1112 1113 1121 1122 1123 1131 1132 1133 2111 2112 2113 2121 2122 2123 2131 2132 2133 3111 3112 3113 3121 3122 3123 3131 3132 3133 4111 4112 4113 4121 4122 4123 -dof 1 localForce
    wipeAnalysis #remove analysis
    pattern Plain 101 Linear {
    '''
    GravityAnalysis='''
    }
    # Gravity Analysis---------------------------------------------------------------------
    # Gravity-analysis parameters -- load-controlled static analysis
    #set Tol 1.0e-6;         # convergence tolerance for test
    #variable constraintsTypeGravity Plain;      # default;
    if {  [info exists RigidDiaphragm] == 1} {
        if {$RigidDiaphragm=="ON"} {
            variable constraintsTypeGravity Lagrange;   #  large model: try Transformation
        };  # if rigid diaphragm is on
    };  # if rigid diaphragm exists
    constraints $constraintsTypeGravity ;           # how it handles boundary conditions
    numberer    Plain
    system      UmfPack
    set Tol 1.0e-6;         # convergence tolerance for test
    test        EnergyIncr  $Tol 50
    algorithm   Newton
    set NstepGravity 10;        # apply gravity in 10 steps
    set DGravity [expr 1./$NstepGravity];   # first load increment;
    integrator  LoadControl $DGravity
    analysis    Static
    analyze     $NstepGravity
    
    # maintain constant gravity loads and reset time to zero
    loadConst -time 0.0
    puts "GRAVITY ANALYSIS COMPLETE"
    puts "Elastic Model Built"
    '''
    # Generacion archivo tcl Elastic Model
    file=open(nameFile,'w')
    file.write(Star)
    file.write(Libraries)
    file.write(Geometry1+x2+x3+Geometry2)
    file.write(Section)
    file.write(Material)
    file.write(Element)
    file.write(massOpenSees.to_string(header=False, index=False))
    file.write(ModalAnalysis)
    file.write(recorderElastic)
    file.write(loadOpensees.to_string(header=False, index=False))
    file.write(GravityAnalysis)
    file.close()
    
    #Run OpenSees: Elastic Model
    bf.runModel(rutaOpenSees,rutascript) 
    
    #%% INELASTIC MODEL
    StartIn1='''
    # Tipical building with elastic section 
    wipe
    # Start script
    model BasicBuilder -ndm 3 -ndf 6;
    '''
    DataBareFrame='set dataDir DataInelasticBareModel'+str(cont)+';'
    
    Start2='''
    file mkdir $dataDir;
    file mkdir $dataDir/Drift;
    file mkdir $dataDir/Rotations/Columns;
    file mkdir $dataDir/Rotations/Beams;
    file mkdir $dataDir/Rotations/Girds;
    file mkdir $dataDir/LocalForce/Columns;
    file mkdir $dataDir/LocalForce/Beams;
    file mkdir $dataDir/LocalForce/Girds;
    file mkdir $dataDir/ZeroColumn
    logFile "$dataDir/Log_model.log";
    '''
    StartIn=StartIn1+DataBareFrame+Start2
    
    MaterialIn1='''
    # MATERIALS-----------------------------------------------------------------------------------------------------------
    # material properties:
    # Concrete
    '''
    s_fcIn='set fc '+str(fc)+';           # concrete nominal compressive strength'
    
    MaterialIn2='''
    set Ec [expr 4700*sqrt($fc)];    # concrete Young's Modulus
    set nu 0.2;         # Poisson's ratio
    set Gc [expr $Ec/2./[expr 1+$nu]];      # Torsional stiffness Modulus
    set Ubig 1e10
    #Reinforcing Steel
    set fy 420
    set fu 550
    set Es 200e3
    '''
    
    MaterialIn=MaterialIn1+s_fcIn+MaterialIn2
    
    SectionIn='''
    # SECTION PROPERTIES-------------------------------------------------------------------------------------------------------------
    # Define section tags
    set ColSecTag 1
    set BeamSecTag 2
    set GirdSecTag 3
    # Section Properties:
    # Reinforcing Steel
    set dbL [expr 12*$mm];  # Diameter Longitudinal
    set dbV [expr 8*$mm];  # Diameter Transverse
    set Cover [expr 30*$mm];    #Concrete Cover
    
    # Column
    set HCol [expr 200*$mm];     # Column section height-Local axes Y
    set BCol [expr 300*$mm];     # Column sectoin width
    set sCol [expr 100*$mm]; # Column stirrup spacing
    set H     [expr $y2*$m];    #Floor height
    
    # Beam
    set HBeam [expr 200*$mm];        # Beam depth -- perpendicular to bending axis
    set BBeam [expr 300*$mm];        # Beam width -- parallel to bending axis
    set sBeam [expr 100*$mm]; # Beam stirrup spacing
    
    # Gird
    set HGird [expr 200*$mm];        # Girder depth -- perpendicular to bending axis
    set BGird [expr 300*$mm];        # Girder width -- parallel to bending axis
    set sGird [expr 100*$mm]; # Gird stirrup spacing
    
    # Reinforcement ratios (These are computed based on the reinforcement provided at each plascti hige zone)
    set rC_top 0.0073 ;  set rC_web 0.00 ;   set rC_bot 0.0073 ;  set rC_shr   0.003351;    # Column
    #set rC_top 0.004833 ;  set rC_web 0.004833 ;   set rC_bot 0.004833 ;  set rC_shr   0.0032828;    # Column
    set rB1_top 0.0096664;    set rB1_web 0;  set rB1_bot 0.0096664;    set rB1_shr 0.003351;    # Beam Section 1
    set rB1_top 0.0096664;    set rB1_web 0;  set rB1_bot 0.0096664;    set rB1_shr 0.003351;    # Beam Section 1
    set rG1_top 0.0096664;    set rG1_web 0;  set rG1_bot 0.0096664;    set rG1_shr 0.003351;    # Gird Section 1
    
    #ELEMENTS------------------------------------------------------------------------------------
    set IDColTransf 1; # all columns
    set IDBeamTransf 2; # all beams
    set IDGirdTransf 3; # all girds
    
    set ColTransfType PDelta ;      # options for columns: Linear PDelta  Corotational 
    geomTransf $ColTransfType  $IDColTransf  0 0 1;         # orientation of column stiffness affects bidirectional response.
    geomTransf Linear $IDBeamTransf 0 0 1
    geomTransf Linear $IDGirdTransf -1 0 0
    
    # Open a set of files so that the properties of the beams, columns and joint elements creted using the provided procedures can be examined later.
    set pfile_jnts [open $dataDir/Properties_joints.txt w];
    set pfile_bms [open $dataDir/Properties_beams.txt w];
    set pfile_cols [open $dataDir/Properties_columnn.txt w];
    set shearfolder [open $dataDir/shearfolder.txt w];
    ## element connectivity
    
    '''
    inelasticElement=bf.inelasticElement(rutafolder,modelo,rutaAxialLoad,fc,rho_Ash)
    
    ElementIn2='''
    # Close the properties files
    close $pfile_jnts
    close $pfile_bms
    close $pfile_cols
    close $shearfolder
    
    set TagElem [getEleTags]
    set elementos [open elementos.txt w]
    foreach i $TagElem {
        set b [eleNodes $i]
        puts $elementos [format "%s %s" $i $b]
    }
    close $elementos
    # Mass
    # mass Node Mx My Mz MRx MRy MRz
    '''
    recorder='''
    #-------------------------Set up parameters that are particular to the model for displacement control
    set IDctrlNode 1141;    # node where displacement is read for displacement control
    
    # RECORDERS
    # Displacements
    recorder Node -file $dataDir/DFree.out -time -node $IDctrlNode -dof 1 disp;         # displacements of free node
    # Drifts
    recorder Drift -file $dataDir/Drift/DriftP1.out -time -iNode 101  -jNode 111  -dof 1 2 3 4 5 6  -perpDirn 2 ;               # Piso 1
    recorder Drift -file $dataDir/Drift/DriftP2.out -time -iNode 1121  -jNode 1131  -dof 1 2 3 4 5 6 -perpDirn 2 ;               # Piso 2
    recorder Drift -file $dataDir/Drift/DriftP3.out -time -iNode 1131  -jNode 1141  -dof 1  2 3 4 5 6 -perpDirn 2 ;               # Piso 3
    # recorder Node -file $Name_Folder/Disp.out -time -nodeRange 10 22 -dof 1 disp;   
    
    # Reactions
    #recorder Node -file $dataDir/RBase.out -time -node 101 102 103 201 202 203 301 302 303 401 402 403 -dof 1 reaction;
    recorder Node -file $dataDir/RBase.out -time -node 101 102 103 201 202 203 301 302 303 401 402 403 -dof 1 reaction;
    
    # Forces and Rotations
    # Columns
    recorder Element -file $dataDir/Rotations/Columns/BasicDeformationColumns.out -time -ele 1111 1112 1113 2111 2112 2113 3111 3112 3113 4111 4112 4113 1121 1122 1123 2121 2122 2123 3121 3122 3123 4121 4122 4123 1131 1132 1133 2131 2132 2132 3131 3132 3133  basicDeformation
    recorder Element -file $dataDir/Rotations/Columns/PlasticDeformationColumns.out -time -ele 1111 1112 1113 2111 2112 2113 3111 3112 3113 4111 4112 4113 1121 1122 1123 2121 2122 2123 3121 3122 3123 4121 4122 4123 1131 1132 1133 2131 2132 2132 3131 3132 3133  plasticDeformation
    recorder Element -file $dataDir/LocalForce/Columns/LocalForceEleColumns.out -time -ele 1111 1112 1113 2111 2112 2113 3111 3112 3113 4111 4112 4113 1121 1122 1123 2121 2122 2123 3121 3122 3123 4121 4122 4123 1131 1132 1133 2131 2132 2132 3131 3132 3133 localForce
    #recorder Element -file $dataDir/LocalForce/Columns/LocalForceEleColumns.out -time -ele 1111 localForce
    
    # Beams
    recorder Element -file $dataDir/Rotations/Beams/BasicDeformationBeams.out -time -ele 1211 1212 2211 2212 3211 3212 4211 4212 1221 1222 2221 2222 3221 3222 4221 4222 1231 1232 2231 2232 3231 3232 basicDeformation
    recorder Element -file $dataDir/Rotations/Beams/PlasticDeformationBeams.out -time -ele 1211 1212 2211 2212 3211 3212 4211 4212 1221 1222 2221 2222 3221 3222 4221 4222 1231 1232 2231 2232 3231 3232  plasticDeformation
    recorder Element -file $dataDir/LocalForce/Beams/LocalForceEleBeams.out -time -ele 1211 1212 2211 2212 3211 3212 4211 4212 1221 1222 2221 2222 3221 3222 4221 4222 1231 1232 2231 2232 3231 3232 localForce
    # Girds
    recorder Element -file $dataDir/Rotations/Girds/BasicDeformationGirds.out -time -ele 11311 11312 11313 21314 21315 21316 31317 31318 31319 11321 11322 11323 21324 21325 21326 31327 31328 31329 11331 11332 21334 21335 31337 31338 basicDeformation
    recorder Element -file $dataDir/Rotations/Girds/PlasticDeformationGirds.out -time -ele 11311 11312 11313 21314 21315 21316 31317 31318 31319 11321 11322 11323 21324 21325 21326 31327 31328 31329 11331 11332 21334 21335 31337 31338  plasticDeformation
    recorder Element -file $dataDir/LocalForce/Girds/LocalForceEleGirds.out -time -ele 11311 11312 11313 21314 21315 21316 31317 31318 31319 11321 11322 11323 21324 21325 21326 31327 31328 31329 11331 11332 21334 21335 31337 31338 localForce
    #
    ## zero length element
    recorder Element -file $dataDir/ZeroColumn/Force.out -time -ele 11111   11112   11113   11121   11122   11123   11131   11132   11133   11311   11312   11313   11321   11322   11323   11331   11332   12111   12112   12113   12121   12122   12123   12131   12132   12133   13111   13112   13113   13121   13122   13123   13131   13132   13133   14111   14112   14113   14121   14122   14123   21111   21112   21113   21121   21122   21123   21131   21132   21133   21314   21315   21316   21324   21325   21326   21334   21335   22111   22112   22113   22121   22122   22123   22131   22132   22133   23111   23112   23113   23121   23122   23123   23131   23132   23133   24111   24112   24113   24121   24122   24123   31317   31318   31319   31327   31328   31329   31337   31338 -dof 1 force 
    recorder Element -file $dataDir/ZeroColumn/Displacement.out -time -ele 11111   11112   11113   11121   11122   11123   11131   11132   11133   11311   11312   11313   11321   11322   11323   11331   11332   12111   12112   12113   12121   12122   12123   12131   12132   12133   13111   13112   13113   13121   13122   13123   13131   13132   13133   14111   14112   14113   14121   14122   14123   21111   21112   21113   21121   21122   21123   21131   21132   21133   21314   21315   21316   21324   21325   21326   21334   21335   22111   22112   22113   22121   22122   22123   22131   22132   22133   23111   23112   23113   23121   23122   23123   23131   23132   23133   24111   24112   24113   24121   24122   24123   31317   31318   31319   31327   31328   31329   31337   31338 -dof 1 deformation
    #recorder Element -file $dataDir/ZeroColumn/Force.out -time -ele 11111 -dof 2 force
    #recorder Element -file $dataDir/ZeroColumn/Deformation.out -time -ele 11111 -dof 2 deformation
    #remove recorders 
    wipeAnalysis
    pattern Plain 101 Linear {
    '''
    GravityAnalysisIn='''
    }
    # Gravity Analysis---------------------------------------------------------------------
    # Gravity-analysis parameters -- load-controlled static analysis
    #set Tol 1.0e-6;         # convergence tolerance for test
    #variable constraintsTypeGravity Plain;      # default;
    if {  [info exists RigidDiaphragm] == 1} {
        if {$RigidDiaphragm=="ON"} {
            variable constraintsTypeGravity Lagrange;   #  large model: try Transformation
        };  # if rigid diaphragm is on
    };  # if rigid diaphragm exists
    constraints $constraintsTypeGravity ;           # how it handles boundary conditions
    numberer    Plain
    system      UmfPack
    set Tol 1.0e-6;         # convergence tolerance for test
    test        EnergyIncr  $Tol 50
    algorithm   Newton
    set NstepGravity 10;        # apply gravity in 10 steps
    set DGravity [expr 1./$NstepGravity];   # first load increment;
    integrator  LoadControl $DGravity
    analysis    Static
    analyze     $NstepGravity
    
    # maintain constant gravity loads and reset time to zero
    loadConst -time 0.0
    puts "GRAVITY ANALYSIS COMPLETE"
    puts "Model Built"
    '''
    DataPushover1='''
    # PUSHOVER ANALYSIS 
    # calculated MODEL PARAMETERS, particular to this model
    # characteristics of pushover analysis
    set IDctrlDOF 1;       # degree of freedom of displacement read for displacement control
    set LBuilding $y4;       # total building height
    set LunitTXT "m";            # define basic-unit text for output
    
    set Dmax [expr 0.06*$LBuilding ];    # maximum displacement of pushover. push to 10% drift.
    set Dincr [expr 0.0006*$LBuilding ];  # displacement increment. you want this to be small, but not too small to slow analysis
    
    # Shape of load
    pattern Plain 200 Linear { 
    '''
    #Forma de la fuerza
    shapeLoad=bf.shapeLoad1Mode(rutaFirstMode)
    
    DataPushover2='''
    }
    
    # Define DISPLAY -------------------------------------------------------------
    set  xPixels 1200;  # height of graphical window in pixels
    set  yPixels 800;   # height of graphical window in pixels
    set  xLoc1 10;  # horizontal location of graphical window (0=upper left-most corner)
    set  yLoc1 10;  # vertical location of graphical window (0=upper left-most corner)
    set ViewScale 1;    # scaling factor for viewing deformed shape, it depends on the dimensions of the model
    #ShapeType :     type of shape to display. # options: ModeShape , NodeNumbers , DeformedShape 
    DisplayModel3D DeformedShape $ViewScale $xLoc1 $yLoc1  $xPixels $yPixels
    recorder plot $dataDir/DFree.out Displ-X [expr $xPixels+10] 10 300 300 -columns 2 1; # a window to plot the nodal displacements versus time
    
    #  ---------------------------------    perform Static Pushover Analysis
    # ----------- set up analysis parameters
    source LibAnalysisStaticParameters.tcl; 
    set fmt1 "%s Pushover analysis: CtrlNode %.3i, dof %.1i, Disp=%.4f %s"; # format for screen/file output of DONE/PROBLEM analysis
    # ----------------------------------------------first analyze command------------------------
    set Nsteps [expr int($Dmax/$Dincr)];        # number of pushover analysis steps
    set ok [analyze $Nsteps];                       # this will return zero if no convergence problems were encountered
    # ----------------------------------------------if convergence failure-------------------------
    if {$ok != 0} {  
        # if analysis fails, we try some other stuff, performance is slower inside this loop
        set Dstep 0.0;
        set ok 0
        while {$Dstep <= 1.0 && $ok == 0} { 
            set controlDisp [nodeDisp $IDctrlNode $IDctrlDOF ]
            set Dstep [expr $controlDisp/$Dmax]
            set ok [analyze 1];                     # this will return zero if no convergence problems were encountered
            if {$ok != 0} {;                # reduce step size if still fails to converge
                set Nk 4;           # reduce step size
                set DincrReduced [expr $Dincr/$Nk];
                integrator DisplacementControl  $IDctrlNode $IDctrlDOF $DincrReduced
                for {set ik 1} {$ik <=$Nk} {incr ik 1} {
                    set ok [analyze 1];                     # this will return zero if no convergence problems were encountered
                    if {$ok != 0} {  
                        # if analysis fails, we try some other stuff
                        # performance is slower inside this loop    global maxNumIterStatic;        # max no. of iterations performed before "failure to converge" is ret'd
                        puts "Trying Newton with Initial Tangent .."
                        test NormDispIncr   $Tol 2000 0
                        algorithm Newton -initial
                        set ok [analyze 1]
                        test $testTypeStatic $TolStatic      $maxNumIterStatic    0
                        algorithm $algorithmTypeStatic
                    }
                    if {$ok != 0} {
                        puts "Trying Broyden .."
                        algorithm Broyden 8
                        set ok [analyze 1 ]
                        algorithm $algorithmTypeStatic
                    }
                    if {$ok != 0} {
                        puts "Trying NewtonWithLineSearch .."
                        algorithm NewtonLineSearch 0.8 
                        set ok [analyze 1]
                        algorithm $algorithmTypeStatic
                    }
                    if {$ok != 0} {;                # stop if still fails to converge
                        puts [format $fmt1 "PROBLEM" $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT]
                        return -1
                    }; # end if
                }; # end for
                integrator DisplacementControl  $IDctrlNode $IDctrlDOF $Dincr;  # bring back to original increment
            }; # end if
        };  # end while loop
    };      # end if ok !0
    # -----------------------------------------------------------------------------------------------------
    
    if {$ok != 0 } {
        puts [format $fmt1 "PROBLEM" $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT]
    } else {
        puts [format $fmt1 "DONE"  $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT]
    }
    
    '''
    
    #%%
    
    # Inelastic Model tcl generation file
    file=open(nameFileInelastic,'w')
    file.write(StartIn)
    file.write(Libraries)
    file.write(Geometry1+x2+x3+Geometry2)
    file.write(MaterialIn)
    file.write(SectionIn)
    file.write(inelasticElement.to_string(header=False, index=False))
    file.write(ElementIn2)
    file.write(massOpenSees.to_string(header=False, index=False))
    file.write(ModalAnalysis)
    file.write(recorder)
    file.write(loadOpensees.to_string(header=False, index=False))
    file.write(GravityAnalysisIn)
    file.write(DataPushover1)
    file.write(shapeLoad.to_string(header=False, index=False))
    file.write(DataPushover2)
    
    file.close()
    
    #Run OpenSees: Inelastic Model-Gravity Analysis
    rutascriptIn=rutafolder+'\\'+nameFileInelastic
    bf.runModel(rutaOpenSees,rutascriptIn) 
    
    
#%% INFILL FRAME 
def InfillFrame(rutaOpenSees,rutafolder,rutaAxialLoad,rutaFirstMode,nameFile,nameFileInelastic,cont,x,fc,rho,rho_Ash,fcj,fcb):
    import BibliotecaFunciones as bf
   # print('rho dentro funcion',rho)
    # ELASTIC MODEL
    Star1='''
    # Elastic Model
    # Mauricio Guamán
    # Maestría en Ing. Civil-ESPE
    wipe
    # Start script
    model BasicBuilder -ndm 3 -ndf 6;
    '''
    InfillDataElasticModel='set dataDir DataElasticInfillModel'+str(cont)+';'
    Star2='''
    file mkdir $dataDir;
    logFile "$dataDir/Log_model.log"; 
    '''
    Star=Star1+InfillDataElasticModel+Star2
    Libraries='''
    # Libraries
    source Units.tcl; 
    source rcBC_nonDuct.tcl; 
    source DisplayPlane.tcl; 
    source DisplayModel2D.tcl; 
    source DisplayModel3D.tcl; 
    source infill.tcl;
    '''
    Geometry1='''
    #Define Geometry
    # Grids
    set x1 0.00; 
    '''
    x2='set x2 '+ str(x)+';'
    x3='set x3 '+ str(x+x)+';'
    
    Geometry2='''
    set y1 0.00;
    set y2 2.50;
    set y3 5.00;
    set y4 7.50;
    set z1 0.00;
    set z2 3.50;
    set z3 7.50;
    set z4 11.0;
    # Defina Nodal Coordinate
    # node $tag_node $X  $Y  $Z
    node 101 $x1 $y1 $z1
    node 102 $x2 $y1 $z1
    node 103 $x3 $y1 $z1
    node 111 $x1 $y2 $z1
    node 112 $x2 $y2 $z1
    node 113 $x3 $y2 $z1
    node 121 $x1 $y3 $z1
    node 122 $x2 $y3 $z1
    node 123 $x3 $y3 $z1
    node 131 $x1 $y4 $z1
    node 132 $x2 $y4 $z1
    node 133 $x3 $y4 $z1
    
    node 201 $x1 $y1 $z2
    node 202 $x2 $y1 $z2
    node 203 $x3 $y1 $z2
    node 211 $x1 $y2 $z2
    node 212 $x2 $y2 $z2
    node 213 $x3 $y2 $z2
    node 221 $x1 $y3 $z2
    node 222 $x2 $y3 $z2
    node 223 $x3 $y3 $z2
    node 231 $x1 $y4 $z2
    node 232 $x2 $y4 $z2
    node 233 $x3 $y4 $z2
    
    node 301 $x1 $y1 $z3
    node 302 $x2 $y1 $z3
    node 303 $x3 $y1 $z3
    node 311 $x1 $y2 $z3
    node 312 $x2 $y2 $z3
    node 313 $x3 $y2 $z3
    node 321 $x1 $y3 $z3
    node 322 $x2 $y3 $z3
    node 323 $x3 $y3 $z3
    node 331 $x1 $y4 $z3
    node 332 $x2 $y4 $z3
    node 333 $x3 $y4 $z3
    
    node 401 $x1 $y1 $z4
    node 402 $x2 $y1 $z4
    node 403 $x3 $y1 $z4
    node 411 $x1 $y2 $z4
    node 412 $x2 $y2 $z4
    node 413 $x3 $y2 $z4
    node 421 $x1 $y3 $z4
    node 422 $x2 $y3 $z4
    node 423 $x3 $y3 $z4
    
    # Variables de nodos a emplear
    set Nudos [getNodeTags];
    set l [llength $Nudos]
    set NodeEnd [lindex $Nudos end];      # Se obtiene el número total de nodos
    set nudos [open nudos.txt w]
    foreach i $Nudos {
        set a [nodeCoord $i]
        puts $nudos [format "%s %s" $i $a]
    }
    close $nudos
    #-------------------Rigid diaphragm nodes
    set RigidDiaphragm ON;     # options: ON, OFF. specify this before the analysis parameters are set the constraints are handled differently.                          
    # Coordenadas del diafragma rigido
    #node $tag_node $X  $Y  $Z
    node 1121   [expr $x3/2]    $y2    [expr $z4/2];      # master nodes for rigid diaphragm -- story 2, bay 1, frame 1-2
    node 1131   [expr $x3/2]    $y3    [expr $z4/2];      # master nodes for rigid diaphragm -- story 3, bay 1, frame 1-2
    node 1141   [expr $x3/2]    $y4    [expr $z3/2];      # master nodes for rigid diaphragm -- story 4, bay 1, frame 1-2
    # Constraints for rigid diaphragm master nodes
    #fix $tag_node  #fix_X #fix_Y #fix_Z #fix_RX #fix_RY #fix_RZ 
    fix 1121 0  1  0  1  0  1
    fix 1131 0  1  0  1  0  1
    fix 1141 0  1  0  1  0  1
    # ------------------------define Rigid Diaphram, dof 2 is normal to floor
    set perpDirn 2;
    rigidDiaphragm $perpDirn 1121 111 112 113 211 212 213 311 312 313 411 412 413;  # level 2
    rigidDiaphragm $perpDirn 1131 121 122 123 221 222 223 321 322 323 421 422 423;  # level 3 
    rigidDiaphragm $perpDirn 1141 131 132 133 231 232 233 331 332 333;  # level 3 
    
    #-------------------------Support
    fixY 0.0 1 1 1 1 1 1;
    '''
    Section='''
    # Section-------------------------------------------------------------------------------------------------------------
    # Define section tags
    set ColSecTag 1
    set BeamSecTag 2
    set GirdSecTag 3
    # Section Properties:
    set HCol [expr 200*$mm];     # square-Column width
    set BCol [expr 300*$mm]
    set HBeam [expr 200*$mm];        # Beam depth -- perpendicular to bending axis
    set BBeam [expr 300*$mm];        # Beam width -- parallel to bending axis
    set HGird [expr 200*$mm];        # Girder depth -- perpendicular to bending axis
    set BGird [expr 300*$mm];        # Girder width -- parallel to bending axis
    '''
    Material1='''
    # material properties:
    '''    
    s_fc='set fc [expr '+str(fc)+'*$MPa];           # concrete nominal compressive strength'
    
    Material2='''
    set Ec [expr 4700*sqrt($fc/1000)*$MPa];    # concrete Young's Modulus
    set nu 0.2;         # Poisson's ratio
    set Gc [expr $Ec/2./[expr 1+$nu]];      # Torsional stiffness Modulus
    set Ubig 1e10
    set J $Ubig;            # set large torsional stiffness
    # column section properties:
    set AgCol [expr $HCol*$BCol];       # rectuangular-Column cross-sectional area
    set IzCol [expr 1*1./12*$BCol*pow($HCol,3)];  # about-local-z Rect-Column gross moment of inertial
    set IyCol [expr 1*1./12*$HCol*pow($BCol,3)];  # about-local-z Rect-Column gross moment of inertial
    # beam sections:
    set AgBeam [expr $HBeam*$BBeam];        # rectuangular-Beam cross-sectional area
    set IzBeam [expr 1*1./12*$BBeam*pow($HBeam,3)];   # about-local-z Rect-Beam cracked moment of inertial
    set IyBeam [expr 1*1./12*$HBeam*pow($BBeam,3)];   # about-local-y Rect-Beam cracked moment of inertial
    # girder sections:
    set AgGird [expr $HGird*$BGird];        # rectuangular-Girder cross-sectional area
    set IzGird [expr 1*1./12*$BGird*pow($HGird,3)];   # about-local-z Rect-Girder cracked moment of inertial
    set IyGird [expr 1*1./12*$HGird*pow($BGird,3)];   # about-local-y Rect-Girder cracked moment of inertial   
    section Elastic $ColSecTag $Ec $AgCol $IzCol $IyCol $Gc $J
    section Elastic $BeamSecTag $Ec $AgBeam $IzBeam $IyBeam $Gc $J
    section Elastic $GirdSecTag $Ec $AgGird $IzGird $IyGird $Gc $J
    '''
    Material=Material1+s_fc+Material2
    
    Element='''
    #-------------------------Elements
    set IDColTransf 1; # all columns
    set IDBeamTransf 2; # all beams
    set IDGirdTransf 3; # all girds
    
    set ColTransfType Linear ;      # options for columns: Linear PDelta  Corotational 
    geomTransf $ColTransfType  $IDColTransf  0 0 1;         # orientation of column stiffness affects bidirectional response.
    geomTransf Linear $IDBeamTransf 0 0 1
    geomTransf Linear $IDGirdTransf -1 0 0
    
    # element connectivity
    #Frame 1
    #Columns
    #element elasticBeamColumn $eleTag $iNode $jNode $A $E $G $J $Iy $Iz $transfTag <-mass $massDens> <-cMass>
    #Level 1
    element elasticBeamColumn  1111 101 111 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    element elasticBeamColumn  1112 102 112 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    element elasticBeamColumn  1113 103 113 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    #Level 2
    element elasticBeamColumn  1121 111 121 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    element elasticBeamColumn  1122 112 122 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    element elasticBeamColumn  1123 113 123 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    #Level 3
    element elasticBeamColumn  1131 121 131 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    element elasticBeamColumn  1132 122 132 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    element elasticBeamColumn  1133 123 133 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    
    #Beams
    #Level 1
    element elasticBeamColumn  1211 111 112 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDBeamTransf
    element elasticBeamColumn  1212 112 113 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDBeamTransf
    #Level 2
    element elasticBeamColumn  1221 121 122 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDBeamTransf
    element elasticBeamColumn  1222 122 123 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDBeamTransf
    #Level 3
    element elasticBeamColumn  1231 131 132 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDBeamTransf
    element elasticBeamColumn  1232 132 133 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDBeamTransf
    
    #Frame 2
    #Columns
    #Level 1
    element elasticBeamColumn  2111 201 211 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    element elasticBeamColumn  2112 202 212 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    element elasticBeamColumn  2113 203 213 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    #Level 2
    element elasticBeamColumn  2121 211 221 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    element elasticBeamColumn  2122 212 222 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    element elasticBeamColumn  2123 213 223 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    #Level 3
    element elasticBeamColumn  2131 221 231 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    element elasticBeamColumn  2132 222 232 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    element elasticBeamColumn  2133 223 233 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    
    #Beams
    #Level 1
    element elasticBeamColumn  2211 211 212 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDBeamTransf
    element elasticBeamColumn  2212 212 213 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDBeamTransf
    #Level 2
    element elasticBeamColumn  2221 221 222 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDBeamTransf
    element elasticBeamColumn  2222 222 223 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDBeamTransf
    #Level 3
    element elasticBeamColumn  2231 231 232 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDBeamTransf
    element elasticBeamColumn  2232 232 233 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDBeamTransf
    
    #Frame 3
    #Columns
    #Level 1
    element elasticBeamColumn  3111 301 311 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    element elasticBeamColumn  3112 302 312 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    element elasticBeamColumn  3113 303 313 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    #Level 2
    element elasticBeamColumn  3121 311 321 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    element elasticBeamColumn  3122 312 322 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    element elasticBeamColumn  3123 313 323 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    #Level 3
    element elasticBeamColumn  3131 321 331 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    element elasticBeamColumn  3132 322 332 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    element elasticBeamColumn  3133 323 333 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    
    #Beams
    #Level 1
    element elasticBeamColumn  3211 311 312 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDBeamTransf
    element elasticBeamColumn  3212 312 313 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDBeamTransf
    #Level 2
    element elasticBeamColumn  3221 321 322 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDBeamTransf
    element elasticBeamColumn  3222 322 323 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDBeamTransf
    #Level 3
    element elasticBeamColumn  3231 331 332 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDBeamTransf
    element elasticBeamColumn  3232 332 333 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDBeamTransf
    
    #Frame 4
    #Columns
    #Level 1
    element elasticBeamColumn  4111 401 411 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    element elasticBeamColumn  4112 402 412 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    element elasticBeamColumn  4113 403 413 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    #Level 2
    element elasticBeamColumn  4121 411 421 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    element elasticBeamColumn  4122 412 422 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    element elasticBeamColumn  4123 413 423 $AgCol $Ec $Gc $J $IyCol $IzCol $IDColTransf
    
    #Beams
    #Level 1
    element elasticBeamColumn  4211 411 412 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDBeamTransf
    element elasticBeamColumn  4212 412 413 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDBeamTransf
    #Level 2
    element elasticBeamColumn  4221 421 422 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDBeamTransf
    element elasticBeamColumn  4222 422 423 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDBeamTransf
    
    #Girders
    #Level 1
    element elasticBeamColumn  11311 111 211 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    element elasticBeamColumn  11312 211 311 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    element elasticBeamColumn  11313 311 411 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    element elasticBeamColumn  21314 112 212 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    element elasticBeamColumn  21315 212 312 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    element elasticBeamColumn  21316 312 412 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    element elasticBeamColumn  31317 113 213 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    element elasticBeamColumn  31318 213 313 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    element elasticBeamColumn  31319 313 413 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    
    #Level 2
    element elasticBeamColumn  11321 121 221 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    element elasticBeamColumn  11322 221 321 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    element elasticBeamColumn  11323 321 421 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    element elasticBeamColumn  21324 122 222 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    element elasticBeamColumn  21325 222 322 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    element elasticBeamColumn  21326 322 422 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    element elasticBeamColumn  31327 123 223 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    element elasticBeamColumn  31328 223 323 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    element elasticBeamColumn  31329 323 423 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    
    #Level 3
    element elasticBeamColumn  11331 131 231 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    element elasticBeamColumn  11332 231 331 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    element elasticBeamColumn  21334 132 232 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    element elasticBeamColumn  21335 232 332 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    element elasticBeamColumn  31337 133 233 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    element elasticBeamColumn  31338 233 333 $AgBeam $Ec $Gc $J $IyBeam $IzBeam $IDGirdTransf
    
    set TagElem [getEleTags]
    set elementos [open elementos.txt w]
    foreach i $TagElem {
        set b [eleNodes $i]
        puts $elementos [format "%s %s" $i $b]
    }
    close $elementos
    '''
    Infill1='''
        #-----------------------------------------------------------------------------------------------
    # Infills
    puts "Start Infills"
    # Tags
    # Frame 1
    set eleTag11 11
    set eleTag12 12
    set eleTag13 13
    set eleTag14 14
    set eleTag15 15
    set eleTag16 16
    # Frame 4
    set eleTag41 41
    set eleTag42 42
    set eleTag43 43
    set eleTag44 44
    # Frame A
    set eleTag51 51
    set eleTag52 52
    set eleTag53 53
    set eleTag54 54
    set eleTag55 55
    set eleTag56 56
    set eleTag57 57
    set eleTag58 58
    # Frame C
    set eleTag61 61
    set eleTag62 62
    set eleTag63 63
    set eleTag64 64
    set eleTag65 65
    set eleTag66 66
    set eleTag67 67
    set eleTag68 68
    
    # Frame 1
    set nds11 {111 112 102 101}
    set nds12 {112 113 103 102}
    set nds13 {121 122 112 111}
    set nds14 {122 123 113 112}
    set nds15 {131 132 122 121}
    set nds16 {132 133 123 122}
    # Frame 4
    set nds41 {411 412 402 401}
    set nds42 {412 413 403 402}
    set nds43 {421 422 412 411}
    set nds44 {422 423 413 412}
    # Frame 5
    set nds51 {111 211 201 101}
    set nds52 {211 311 301 201}
    set nds53 {311 411 401 301}
    set nds54 {121 221 211 111}
    set nds55 {221 321 311 211}
    set nds56 {321 421 411 311}
    set nds57 {131 231 221 121}
    set nds58 {231 331 321 221}
    # Frame 6
    set nds61 {113 213 203 103}
    set nds62 {213 313 303 203}
    set nds63 {313 413 403 303}
    set nds64 {123 223 213 113}
    set nds65 {223 323 313 213}
    set nds66 {323 423 413 313}
    set nds67 {133 233 223 123}
    set nds68 {233 333 323 223}
    
    # Vanos 
    # Direction: X
    set lv1x [expr $x2/$mm]; # Longitud vano 1 [mm]
    set lv2x [expr ($x3-$x2)/$mm]; # Longitud vano 2 [mm]
    # Direction: Z
    set lv1z [expr $z2/$mm]; # Longitud vano 1 [mm]
    set lv2z [expr ($z3-$z2)/$mm]; # Longitud vano 2 [mm]
    set lv3z [expr ($z4-$z3)/$mm]; # Longitud vano 2 [mm]
    # Direction: Y
    set lv1y [expr $y2/$mm]; # Longitud vano 1 [mm]
    set lv2y [expr ($y3-$y2)/$mm]; # Longitud vano 2 [mm]
    set lv3y [expr ($y4-$y3)/$mm]; # Longitud vano 2 [mm]
    # Beam
    set hb [expr $HBeam/$mm]
    # Column
    set hc [expr $HCol/$mm]
    set bc [expr $BCol/$mm]
    # Wall thickness 
    set tw 150
    # Strut
    set t "single"
    # Block
    '''
    block='set fcb '+str(fcb)+'; #[MPa]'
    
    Infill2='''
    set ftb 0.2150
    set b 201.33; #[mm]
    set d 400.70
    # Mortar
    '''
    mortar='set fcj '+str(fcj)+'; #[MPa]'
    
    Infill3='''
    set j 10.11; #[mm]
    # Coefficient
    set Uu 1.5
    # Horizontal compressive strength masonry
    set fwh 0.8574; #[Mpa]
    # Vertical compressive strength masonry
    set alpha [expr $j/(4.1*$b)]
    set fwv [expr $fcb/$Uu*($ftb+$alpha*$fcj)/($ftb+$alpha*$fcb)]
    #set fwv 0.94 
    # Elastic Modulus Masonry 
    set Ewh [expr 1000*$fwh]
    set Ewv [expr 1000*$fwv]
    # Poisson Modulus Masonry
    set v 0.1524
    # Shear Modulus Masonry
    set Gw [expr $Ewv/(2*(1+$v))]
    set fwu 0.109
    set fws 0.8162
    set sig_v 0.000
    set GT_inf 1
    set pflag 0 
    set rhoz 1
    '''
    
    s_rho1='set rhox1 '+str(rho[0])+';           # opening coefficient reduction span 1\n'
    s_rho2='set rhox2 '+str(rho[1])+';           # opening coefficient reduction span 2\n'
    s_rho3='set rhox3 '+str(rho[2])+';           # opening coefficient reduction span 3\n'
    s_rho4='set rhox4 '+str(rho[3])+';           # opening coefficient reduction span 4\n'
    s_rho5='set rhox5 '+str(rho[4])+';           # opening coefficient reduction span 5\n'
    s_rho6='set rhox6 '+str(rho[5])+';           # opening coefficient reduction span 6\n'
    s_rho7='set rhox7 '+str(rho[6])+';           # opening coefficient reduction span 7\n'
    s_rho8='set rhox8 '+str(rho[7])+';           # opening coefficient reduction span 8\n'
    s_rho9='set rhox9 '+str(rho[8])+';           # opening coefficient reduction span 9\n'
    s_rho10='set rhox10 '+str(rho[9])+';           # opening coefficient reduction span 10\n'

    
    
    Infill4='''
    # Infill function
    # infill {eleTag typ nds B H hb hc bc tw Ec Ewh Ewv Gw v fwv fwu fws sig_v {pflag 0}} {
    # Frame 1
    infill $eleTag11 $t $nds11 $lv1x $lv1y $hb $hc $bc $tw $Ec $Ewh $Ewv $Gw $v $fwv $fwu $fws $sig_v $rhox1 $pflag 
    infill $eleTag12 $t $nds12 $lv2x $lv1y $hb $hc $bc $tw $Ec $Ewh $Ewv $Gw $v $fwv $fwu $fws $sig_v $rhox2 $pflag
    infill $eleTag13 $t $nds13 $lv1x $lv2y $hb $hc $bc $tw $Ec $Ewh $Ewv $Gw $v $fwv $fwu $fws $sig_v $rhox3 $pflag
    infill $eleTag14 $t $nds14 $lv2x $lv2y $hb $hc $bc $tw $Ec $Ewh $Ewv $Gw $v $fwv $fwu $fws $sig_v $rhox4 $pflag
    infill $eleTag15 $t $nds15 $lv1x $lv3y $hb $hc $bc $tw $Ec $Ewh $Ewv $Gw $v $fwv $fwu $fws $sig_v $rhox5 $pflag
    infill $eleTag16 $t $nds16 $lv2x $lv3y $hb $hc $bc $tw $Ec $Ewh $Ewv $Gw $v $fwv $fwu $fws $sig_v $rhox6 $pflag
    # Frame 4
    infill $eleTag41 $t $nds41 $lv1x $lv1y $hb $hc $bc $tw $Ec $Ewh $Ewv $Gw $v $fwv $fwu $fws $sig_v $rhox7 $pflag
    infill $eleTag42 $t $nds42 $lv2x $lv1y $hb $hc $bc $tw $Ec $Ewh $Ewv $Gw $v $fwv $fwu $fws $sig_v $rhox8 $pflag
    infill $eleTag43 $t $nds43 $lv1x $lv2y $hb $hc $bc $tw $Ec $Ewh $Ewv $Gw $v $fwv $fwu $fws $sig_v $rhox9 $pflag
    infill $eleTag44 $t $nds44 $lv2x $lv2y $hb $hc $bc $tw $Ec $Ewh $Ewv $Gw $v $fwv $fwu $fws $sig_v $rhox10 $pflag
    ## Frame 5
    infill $eleTag51 $t $nds51 $lv1z $lv1y $hb $bc $hc $tw $Ec $Ewh $Ewv $Gw $v $fwv $fwu $fws $sig_v $rhoz $pflag 
    infill $eleTag52 $t $nds52 $lv2z $lv1y $hb $bc $hc $tw $Ec $Ewh $Ewv $Gw $v $fwv $fwu $fws $sig_v $rhoz $pflag
    infill $eleTag53 $t $nds53 $lv3z $lv1y $hb $bc $hc $tw $Ec $Ewh $Ewv $Gw $v $fwv $fwu $fws $sig_v $rhoz $pflag
    infill $eleTag54 $t $nds54 $lv1z $lv2y $hb $bc $hc $tw $Ec $Ewh $Ewv $Gw $v $fwv $fwu $fws $sig_v $rhoz $pflag
    infill $eleTag55 $t $nds55 $lv2z $lv2y $hb $bc $hc $tw $Ec $Ewh $Ewv $Gw $v $fwv $fwu $fws $sig_v $rhoz $pflag
    infill $eleTag56 $t $nds56 $lv3z $lv2y $hb $bc $hc $tw $Ec $Ewh $Ewv $Gw $v $fwv $fwu $fws $sig_v $rhoz $pflag
    infill $eleTag57 $t $nds57 $lv1z $lv3y $hb $bc $hc $tw $Ec $Ewh $Ewv $Gw $v $fwv $fwu $fws $sig_v $rhoz $pflag
    infill $eleTag58 $t $nds58 $lv2z $lv3y $hb $bc $hc $tw $Ec $Ewh $Ewv $Gw $v $fwv $fwu $fws $sig_v $rhoz $pflag
    # Frame 6
    infill $eleTag61 $t $nds61 $lv1z $lv1y $hb $bc $hc $tw $Ec $Ewh $Ewv $Gw $v $fwv $fwu $fws $sig_v $rhoz $pflag
    infill $eleTag62 $t $nds62 $lv2z $lv1y $hb $bc $hc $tw $Ec $Ewh $Ewv $Gw $v $fwv $fwu $fws $sig_v $rhoz $pflag
    infill $eleTag63 $t $nds63 $lv3z $lv1y $hb $bc $hc $tw $Ec $Ewh $Ewv $Gw $v $fwv $fwu $fws $sig_v $rhoz $pflag
    infill $eleTag64 $t $nds64 $lv1z $lv2y $hb $bc $hc $tw $Ec $Ewh $Ewv $Gw $v $fwv $fwu $fws $sig_v $rhoz $pflag
    infill $eleTag65 $t $nds65 $lv2z $lv2y $hb $bc $hc $tw $Ec $Ewh $Ewv $Gw $v $fwv $fwu $fws $sig_v $rhoz $pflag
    infill $eleTag66 $t $nds66 $lv3z $lv2y $hb $bc $hc $tw $Ec $Ewh $Ewv $Gw $v $fwv $fwu $fws $sig_v $rhoz $pflag
    infill $eleTag67 $t $nds67 $lv1z $lv3y $hb $bc $hc $tw $Ec $Ewh $Ewv $Gw $v $fwv $fwu $fws $sig_v $rhoz $pflag
    infill $eleTag68 $t $nds68 $lv2z $lv3y $hb $bc $hc $tw $Ec $Ewh $Ewv $Gw $v $fwv $fwu $fws $sig_v $rhoz $pflag
    
    puts "Masonry done"
    #-----------------------------------------------------------------------------------------
    # Mass
    # mass Node Mx My Mz MRx MRy MRz
    '''
    Infill=Infill1+block+Infill2+mortar+Infill3+s_rho1+s_rho2+s_rho3+s_rho4+s_rho5+s_rho6+s_rho7+s_rho8+s_rho9+s_rho10+Infill4
    # tcl file generation with basic information
    file=open(nameFile,'w')
    file.write(Star)
    file.write(Libraries)
    file.write(Geometry1+x2+x3+Geometry2)
    file.write(Section)
    file.write(Material)
    file.write(Element)
    file.write('puts "Created: Nudos.txt and Elementos.txt " ')
    file.close()
    
    #Run OpenSees to get nodes and elements. 
    # FileAdress
    rutascript=rutafolder+'\\'+nameFile
    bf.runModel(rutaOpenSees,rutascript) 
    
    # MASS
    rutaNudos=rutafolder+'\\'+'nudos.txt'
    rutaElementos=rutafolder+'\\'+'elementos.txt'
    massOpenSees,loadOpensees,modelo,mass=bf.computeMass(rutaNudos,rutaElementos,rutafolder)
    
    ModalAnalysis='''
    # Modal Analysis------------------------------------------------
    set Inf Inf;                        # Inf = Si queremos que muestre información en pantalla 
                                        # NoInf = Si no queremos que muestre información en pantalla
    set Outputs $dataDir;
    set numModes 3;                     # Número de modos a evaluar
    # Outputs: nombre de la carpeta Outputs o "NoOutputs" si no se quieren resultados
    remove recorders 
    wipeAnalysis
    if {$Inf != "NoInf"} {
        puts "------------------------------------------------------"
        puts "MODAL ANALYSIS"
    }
    system UmfPack
    constraints Transformation
    set lambda [eigen -fullGenLapack $numModes]; # -fullGenLapack in case of [nº modes=ngdl and masses]
    set omega {}
    set f {}
    set T {}
    set pi [expr acos(-1.0)];
    
    foreach lam $lambda {
        lappend omega [expr sqrt($lam)]
        lappend F [expr sqrt($lam)/(2*$pi)]
        lappend T [expr (2*$pi)/sqrt($lam)]
    }
    
    if {$Inf != "NoInf"} {
        puts "Periods are: $T -s-"
        puts "Frequencies are: $F -Hz-"
    }
    
    if {$Outputs != "NoOutputs"} {
        file mkdir $Outputs/Modes
        file mkdir $Outputs/Vectors
        for { set k 1 } { $k <= $numModes } { incr k } {
        recorder Node -file [format "$dataDir/Modes/mode%i.out" $k] -time -nodeRange 1 "$::NodeEnd" -dof 1 2 3 4 5 6  "eigen $k"
        recorder Node -file [format "$dataDir/Vectors/eigenvector%i.out" $k] -time -node 1121 1131 1141 -dof 1 2 3 4 5 6  "eigen $k"
        }
        set period "$dataDir/Modes/Periods.txt"
        set Periods [open $period "w"]
        foreach t $T f $F w $omega {
            puts $Periods "$t  $f  $w"
        }
        close $Periods
        record
    }       
    '''
    recorderElastic='''
    recorder Element -file $dataDir/AxialLoad.out -time -ele 1111 1112 1113 1121 1122 1123 1131 1132 1133 2111 2112 2113 2121 2122 2123 2131 2132 2133 3111 3112 3113 3121 3122 3123 3131 3132 3133 4111 4112 4113 4121 4122 4123 -dof 1 localForce
    wipeAnalysis #remove analysis
    pattern Plain 101 Linear {
    '''
    GravityAnalysis='''
    }
    # Gravity Analysis---------------------------------------------------------------------
    # Gravity-analysis parameters -- load-controlled static analysis
    #set Tol 1.0e-8;         # convergence tolerance for test
    #variable constraintsTypeGravity Plain;      # default;
    if {  [info exists RigidDiaphragm] == 1} {
        if {$RigidDiaphragm=="ON"} {
            variable constraintsTypeGravity Lagrange;   #  large model: try Transformation
        };  # if rigid diaphragm is on
    };  # if rigid diaphragm exists
    constraints $constraintsTypeGravity ;           # how it handles boundary conditions
    numberer    Plain
    system      UmfPack
    set Tol 1.0e-6;         # convergence tolerance for test
    test        EnergyIncr  $Tol 50
    algorithm   Newton
    set NstepGravity 10;        # apply gravity in 10 steps
    set DGravity [expr 1./$NstepGravity];   # first load increment;
    integrator  LoadControl $DGravity
    analysis    Static
    analyze     $NstepGravity
    
    # maintain constant gravity loads and reset time to zero
    loadConst -time 0.0
    puts "GRAVITY ANALYSIS COMPLETE"
    puts "Infill Elastic Model Built"
    '''
    # Generacion archivo tcl Elastic Model
    file=open(nameFile,'w')
    file.write(Star)
    file.write(Libraries)
    file.write(Geometry1+x2+x3+Geometry2)
    file.write(Section)
    file.write(Material)
    file.write(Element)
    file.write(Infill)
    file.write(massOpenSees.to_string(header=False, index=False))
    file.write(ModalAnalysis)
    file.write(recorderElastic)
    file.write(loadOpensees.to_string(header=False, index=False))
    file.write(GravityAnalysis)
    file.close()
    
    #Run OpenSees: Elastic Model
    bf.runModel(rutaOpenSees,rutascript) 
    
    #%% INELASTIC MODEL
    StartIn1='''
    # Tipical building with elastic section 
    wipe
    # Start script
    model BasicBuilder -ndm 3 -ndf 6;
    '''
    DataInfillFrame='set dataDir DataInelasticInfillModel'+str(cont)+';'
    
    Start2='''
    file mkdir $dataDir;
    file mkdir $dataDir/Drift;
    file mkdir $dataDir/Rotations/Columns;
    file mkdir $dataDir/Rotations/Beams;
    file mkdir $dataDir/Rotations/Girds;
    file mkdir $dataDir/LocalForce/Columns;
    file mkdir $dataDir/LocalForce/Beams;
    file mkdir $dataDir/LocalForce/Girds;
    file mkdir $dataDir/ZeroColumn
    logFile "$dataDir/Log_model.log";
    '''
    StartIn=StartIn1+DataInfillFrame+Start2
    
    MaterialIn1='''
    # MATERIALS-----------------------------------------------------------------------------------------------------------
    # material properties:
    # Concrete
    '''
    s_fcIn='set fc '+str(fc)+';           # concrete nominal compressive strength'
    
    MaterialIn2='''
    set Ec [expr 4700*sqrt($fc)];    # concrete Young's Modulus
    set nu 0.2;         # Poisson's ratio
    set Gc [expr $Ec/2./[expr 1+$nu]];      # Torsional stiffness Modulus
    set Ubig 1e10
    #Reinforcing Steel
    set fy 420
    set fu 550
    set Es 200e3
    '''
    
    MaterialIn=MaterialIn1+s_fcIn+MaterialIn2
    
    SectionIn='''
    # SECTION PROPERTIES-------------------------------------------------------------------------------------------------------------
    # Define section tags
    set ColSecTag 1
    set BeamSecTag 2
    set GirdSecTag 3
    # Section Properties:
    # Reinforcing Steel
    set dbL [expr 12*$mm];  # Diameter Longitudinal
    set dbV [expr 8*$mm];  # Diameter Transverse
    set Cover [expr 30*$mm];    #Concrete Cover
    
    # Column
    set HCol [expr 200*$mm];     # Column section height-Local axes Y
    set BCol [expr 300*$mm];     # Column sectoin width
    set sCol [expr 100*$mm]; # Column stirrup spacing
    set H     [expr $y2*$m];    #Floor height
    
    # Beam
    set HBeam [expr 200*$mm];        # Beam depth -- perpendicular to bending axis
    set BBeam [expr 300*$mm];        # Beam width -- parallel to bending axis
    set sBeam [expr 100*$mm]; # Beam stirrup spacing
    
    # Gird
    set HGird [expr 200*$mm];        # Girder depth -- perpendicular to bending axis
    set BGird [expr 300*$mm];        # Girder width -- parallel to bending axis
    set sGird [expr 100*$mm]; # Gird stirrup spacing
    
    # Reinforcement ratios (These are computed based on the reinforcement provided at each plascti hige zone)
    set rC_top 0.0073 ;  set rC_web 0.00 ;   set rC_bot 0.0073 ;  set rC_shr   0.003351;    # Column
    #set rC_top 0.004833 ;  set rC_web 0.004833 ;   set rC_bot 0.004833 ;  set rC_shr   0.0032828;    # Column
    set rB1_top 0.0096664;    set rB1_web 0;  set rB1_bot 0.0096664;    set rB1_shr 0.003351;    # Beam Section 1
    set rB1_top 0.0096664;    set rB1_web 0;  set rB1_bot 0.0096664;    set rB1_shr 0.003351;    # Beam Section 1
    set rG1_top 0.0096664;    set rG1_web 0;  set rG1_bot 0.0096664;    set rG1_shr 0.003351;    # Gird Section 1
    
    #ELEMENTS------------------------------------------------------------------------------------
    set IDColTransf 1; # all columns
    set IDBeamTransf 2; # all beams
    set IDGirdTransf 3; # all girds
    
    set ColTransfType PDelta ;      # options for columns: Linear PDelta  Corotational 
    geomTransf $ColTransfType  $IDColTransf  0 0 1;         # orientation of column stiffness affects bidirectional response.
    geomTransf Linear $IDBeamTransf 0 0 1
    geomTransf Linear $IDGirdTransf -1 0 0
    
    # Open a set of files so that the properties of the beams, columns and joint elements creted using the provided procedures can be examined later.
    set pfile_jnts [open $dataDir/Properties_joints.txt w];
    set pfile_bms [open $dataDir/Properties_beams.txt w];
    set pfile_cols [open $dataDir/Properties_columnn.txt w];
    set shearfolder [open $dataDir/shearfolder.txt w];
    ## element connectivity
    '''
    inelasticElement=bf.inelasticElement(rutafolder,modelo,rutaAxialLoad,fc,rho_Ash)
    
    ElementIn2='''
    # Close the properties files
    close $pfile_jnts
    close $pfile_bms
    close $pfile_cols
    close $shearfolder
    
    set TagElem [getEleTags]
    set elementos [open elementos.txt w]
    foreach i $TagElem {
        set b [eleNodes $i]
        puts $elementos [format "%s %s" $i $b]
    }
    close $elementos
    '''
    recorder='''
    #-------------------------Set up parameters that are particular to the model for displacement control
    set IDctrlNode 1141;    # node where displacement is read for displacement control
    
    # RECORDERS
    # Displacements
    recorder Node -file $dataDir/DFree.out -time -node $IDctrlNode -dof 1 disp;         # displacements of free node
    # Drifts
    recorder Drift -file $dataDir/Drift/DriftP1.out -time -iNode 101  -jNode 111  -dof 1 2 3 4 5 6  -perpDirn 2 ;               # Piso 1
    recorder Drift -file $dataDir/Drift/DriftP2.out -time -iNode 1121  -jNode 1131  -dof 1 2 3 4 5 6 -perpDirn 2 ;               # Piso 2
    recorder Drift -file $dataDir/Drift/DriftP3.out -time -iNode 1131  -jNode 1141  -dof 1  2 3 4 5 6 -perpDirn 2 ;               # Piso 3
    # recorder Node -file $Name_Folder/Disp.out -time -nodeRange 10 22 -dof 1 disp;   
    
    # Reactions
    #recorder Node -file $dataDir/RBase.out -time -node 101 102 103 201 202 203 301 302 303 401 402 403 -dof 1 reaction;
    recorder Node -file $dataDir/RBase.out -time -node 101 102 103 201 202 203 301 302 303 401 402 403 -dof 1 reaction;
    
    # Forces and Rotations
    # Columns
    recorder Element -file $dataDir/Rotations/Columns/BasicDeformationColumns.out -time -ele 1111 1112 1113 2111 2112 2113 3111 3112 3113 4111 4112 4113 1121 1122 1123 2121 2122 2123 3121 3122 3123 4121 4122 4123 1131 1132 1133 2131 2132 2132 3131 3132 3133  basicDeformation
    recorder Element -file $dataDir/Rotations/Columns/PlasticDeformationColumns.out -time -ele 1111 1112 1113 2111 2112 2113 3111 3112 3113 4111 4112 4113 1121 1122 1123 2121 2122 2123 3121 3122 3123 4121 4122 4123 1131 1132 1133 2131 2132 2132 3131 3132 3133  plasticDeformation
    recorder Element -file $dataDir/LocalForce/Columns/LocalForceEleColumns.out -time -ele 1111 1112 1113 2111 2112 2113 3111 3112 3113 4111 4112 4113 1121 1122 1123 2121 2122 2123 3121 3122 3123 4121 4122 4123 1131 1132 1133 2131 2132 2132 3131 3132 3133 localForce
    #recorder Element -file $dataDir/LocalForce/Columns/LocalForceEleColumns.out -time -ele 1111 localForce
    
    # Beams
    recorder Element -file $dataDir/Rotations/Beams/BasicDeformationBeams.out -time -ele 1211 1212 2211 2212 3211 3212 4211 4212 1221 1222 2221 2222 3221 3222 4221 4222 1231 1232 2231 2232 3231 3232 basicDeformation
    recorder Element -file $dataDir/Rotations/Beams/PlasticDeformationBeams.out -time -ele 1211 1212 2211 2212 3211 3212 4211 4212 1221 1222 2221 2222 3221 3222 4221 4222 1231 1232 2231 2232 3231 3232  plasticDeformation
    recorder Element -file $dataDir/LocalForce/Beams/LocalForceEleBeams.out -time -ele 1211 1212 2211 2212 3211 3212 4211 4212 1221 1222 2221 2222 3221 3222 4221 4222 1231 1232 2231 2232 3231 3232 localForce
    # Girds
    recorder Element -file $dataDir/Rotations/Girds/BasicDeformationGirds.out -time -ele 11311 11312 11313 21314 21315 21316 31317 31318 31319 11321 11322 11323 21324 21325 21326 31327 31328 31329 11331 11332 21334 21335 31337 31338 basicDeformation
    recorder Element -file $dataDir/Rotations/Girds/PlasticDeformationGirds.out -time -ele 11311 11312 11313 21314 21315 21316 31317 31318 31319 11321 11322 11323 21324 21325 21326 31327 31328 31329 11331 11332 21334 21335 31337 31338  plasticDeformation
    recorder Element -file $dataDir/LocalForce/Girds/LocalForceEleGirds.out -time -ele 11311 11312 11313 21314 21315 21316 31317 31318 31319 11321 11322 11323 21324 21325 21326 31327 31328 31329 11331 11332 21334 21335 31337 31338 localForce
    #
    ## zero length element
    recorder Element -file $dataDir/ZeroColumn/Force.out -time -ele 11111   11112   11113   11121   11122   11123   11131   11132   11133   11311   11312   11313   11321   11322   11323   11331   11332   12111   12112   12113   12121   12122   12123   12131   12132   12133   13111   13112   13113   13121   13122   13123   13131   13132   13133   14111   14112   14113   14121   14122   14123   21111   21112   21113   21121   21122   21123   21131   21132   21133   21314   21315   21316   21324   21325   21326   21334   21335   22111   22112   22113   22121   22122   22123   22131   22132   22133   23111   23112   23113   23121   23122   23123   23131   23132   23133   24111   24112   24113   24121   24122   24123   31317   31318   31319   31327   31328   31329   31337   31338 -dof 1 force 
    recorder Element -file $dataDir/ZeroColumn/Displacement.out -time -ele 11111   11112   11113   11121   11122   11123   11131   11132   11133   11311   11312   11313   11321   11322   11323   11331   11332   12111   12112   12113   12121   12122   12123   12131   12132   12133   13111   13112   13113   13121   13122   13123   13131   13132   13133   14111   14112   14113   14121   14122   14123   21111   21112   21113   21121   21122   21123   21131   21132   21133   21314   21315   21316   21324   21325   21326   21334   21335   22111   22112   22113   22121   22122   22123   22131   22132   22133   23111   23112   23113   23121   23122   23123   23131   23132   23133   24111   24112   24113   24121   24122   24123   31317   31318   31319   31327   31328   31329   31337   31338 -dof 1 deformation
    #recorder Element -file $dataDir/ZeroColumn/Force.out -time -ele 11111 -dof 2 force
    #recorder Element -file $dataDir/ZeroColumn/Deformation.out -time -ele 11111 -dof 2 deformation
    #remove recorders 
    wipeAnalysis
    pattern Plain 101 Linear {
    '''
    GravityAnalysisIn='''
    }
    # Gravity Analysis---------------------------------------------------------------------
    # Gravity-analysis parameters -- load-controlled static analysis
    #set Tol 1.0e-8;         # convergence tolerance for test
    #variable constraintsTypeGravity Plain;      # default;
    if {  [info exists RigidDiaphragm] == 1} {
        if {$RigidDiaphragm=="ON"} {
            variable constraintsTypeGravity Lagrange;   #  large model: try Transformation
        };  # if rigid diaphragm is on
    };  # if rigid diaphragm exists
    constraints $constraintsTypeGravity ;           # how it handles boundary conditions
    numberer    Plain
    system      UmfPack
    set Tol 1.0e-6;         # convergence tolerance for test
    test        EnergyIncr  $Tol 50
    algorithm   Newton
    set NstepGravity 10;        # apply gravity in 10 steps
    set DGravity [expr 1./$NstepGravity];   # first load increment;
    integrator  LoadControl $DGravity
    analysis    Static
    analyze     $NstepGravity
    
    # maintain constant gravity loads and reset time to zero
    loadConst -time 0.0
    puts "GRAVITY ANALYSIS COMPLETE"
    puts "Model Built"
    '''
    DataPushover1='''
    # PUSHOVER ANALYSIS 
    # calculated MODEL PARAMETERS, particular to this model
    # characteristics of pushover analysis
    set IDctrlDOF 1;       # degree of freedom of displacement read for displacement control
    set LBuilding $y4;       # total building height
    set LunitTXT "m";            # define basic-unit text for output
    
    set Dmax [expr 0.06*$LBuilding ];    # maximum displacement of pushover. push to 10% drift.
    set Dincr [expr 0.0006*$LBuilding ];  # displacement increment. you want this to be small, but not too small to slow analysis
    
    # Shape of load
    pattern Plain 200 Linear { 
    '''
    #Forma de la fuerza
    shapeLoad=bf.shapeLoad1Mode(rutaFirstMode)
    
    DataPushover2='''
    }
    
    # Define DISPLAY -------------------------------------------------------------
    set  xPixels 1200;  # height of graphical window in pixels
    set  yPixels 800;   # height of graphical window in pixels
    set  xLoc1 10;  # horizontal location of graphical window (0=upper left-most corner)
    set  yLoc1 10;  # vertical location of graphical window (0=upper left-most corner)
    set ViewScale 1;    # scaling factor for viewing deformed shape, it depends on the dimensions of the model
    #ShapeType :     type of shape to display. # options: ModeShape , NodeNumbers , DeformedShape 
    DisplayModel3D DeformedShape $ViewScale $xLoc1 $yLoc1  $xPixels $yPixels
    recorder plot $dataDir/DFree.out Displ-X [expr $xPixels+10] 10 300 300 -columns 2 1; # a window to plot the nodal displacements versus time
    
    #  ---------------------------------    perform Static Pushover Analysis
    # ----------- set up analysis parameters
    source LibAnalysisStaticParameters.tcl; 
    set fmt1 "%s Pushover analysis: CtrlNode %.3i, dof %.1i, Disp=%.4f %s"; # format for screen/file output of DONE/PROBLEM analysis
    # ----------------------------------------------first analyze command------------------------
    set Nsteps [expr int($Dmax/$Dincr)];        # number of pushover analysis steps
    set ok [analyze $Nsteps];                       # this will return zero if no convergence problems were encountered
    # ----------------------------------------------if convergence failure-------------------------
    if {$ok != 0} {  
        # if analysis fails, we try some other stuff, performance is slower inside this loop
        set Dstep 0.0;
        set ok 0
        while {$Dstep <= 1.0 && $ok == 0} { 
            set controlDisp [nodeDisp $IDctrlNode $IDctrlDOF ]
            set Dstep [expr $controlDisp/$Dmax]
            set ok [analyze 1];                     # this will return zero if no convergence problems were encountered
            if {$ok != 0} {;                # reduce step size if still fails to converge
                set Nk 4;           # reduce step size
                set DincrReduced [expr $Dincr/$Nk];
                integrator DisplacementControl  $IDctrlNode $IDctrlDOF $DincrReduced
                for {set ik 1} {$ik <=$Nk} {incr ik 1} {
                    set ok [analyze 1];                     # this will return zero if no convergence problems were encountered
                    if {$ok != 0} {  
                        # if analysis fails, we try some other stuff
                        # performance is slower inside this loop    global maxNumIterStatic;        # max no. of iterations performed before "failure to converge" is ret'd
                        puts "Trying Newton with Initial Tangent .."
                        test NormDispIncr   $Tol 2000 0
                        algorithm Newton -initial
                        set ok [analyze 1]
                        test $testTypeStatic $TolStatic      $maxNumIterStatic    0
                        algorithm $algorithmTypeStatic
                    }
                    if {$ok != 0} {
                        puts "Trying Broyden .."
                        algorithm Broyden 8
                        set ok [analyze 1 ]
                        algorithm $algorithmTypeStatic
                    }
                    if {$ok != 0} {
                        puts "Trying NewtonWithLineSearch .."
                        algorithm NewtonLineSearch 0.8 
                        set ok [analyze 1]
                        algorithm $algorithmTypeStatic
                    }
                    if {$ok != 0} {;                # stop if still fails to converge
                        puts [format $fmt1 "PROBLEM" $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT]
                        return -1
                    }; # end if
                }; # end for
                integrator DisplacementControl  $IDctrlNode $IDctrlDOF $Dincr;  # bring back to original increment
            }; # end if
        };  # end while loop
    };      # end if ok !0
    # -----------------------------------------------------------------------------------------------------
    
    if {$ok != 0 } {
        puts [format $fmt1 "PROBLEM" $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT]
    } else {
        puts [format $fmt1 "DONE"  $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT]
    }
    
    '''
    
    #%%
    
    # Inelastic Model tcl generation file
    file=open(nameFileInelastic,'w')
    file.write(StartIn)
    file.write(Libraries)
    file.write(Geometry1+x2+x3+Geometry2)
    file.write(MaterialIn)
    file.write(SectionIn)
    file.write(inelasticElement.to_string(header=False, index=False))
    file.write(ElementIn2)
    file.write(Infill)
    file.write(massOpenSees.to_string(header=False, index=False))
    file.write(ModalAnalysis)
    file.write(recorder)
    file.write(loadOpensees.to_string(header=False, index=False))
    file.write(GravityAnalysisIn)
    file.write(DataPushover1)
    file.write(shapeLoad.to_string(header=False, index=False))
    file.write(DataPushover2)
    
    file.close()
    
    #Run OpenSees: Inelastic Model-Gravity Analysis
    rutascriptIn=rutafolder+'\\'+nameFileInelastic
    bf.runModel(rutaOpenSees,rutascriptIn) 
    return massOpenSees,loadOpensees,modelo,mass
#%% TRANFORMATION TO ADRS SYSTEM 
def adrs(rutaFirstMode,mass,h,BasalSM,DtopSM):
    import numpy as np
    firstMode=np.genfromtxt(rutaFirstMode)
    firstMode=firstMode[9,1:]
    firstMode=abs(firstMode[0:-1:6])
    masa=np.array(mass.Mass)
    aux1=firstMode*masa
    meq=sum([i**2 for i in aux1] )/sum([i**2 for i in firstMode]*masa)
    heq=sum(h*firstMode*masa)/sum(firstMode*masa)
    aeq=BasalSM/(meq*9.81)
    Gamma=sum(firstMode*masa)/sum([i**2 for i in firstMode]*masa)
    Sd=DtopSM/Gamma
    return aeq,Sd
#%% TRASFORMATION TO QUADRILINEAR CURVE
def cuadrilinear(aeq,Sd):
    import numpy as np
    # POINT 1
    ad1y=0.75*max(aeq)
    
    cont=0
    position=0
    for ad in aeq:
        #print('Ad',ad)
        #print('Ad1y-Ad',abs(ad1y-ad))
        if abs(ad1y-ad)<0.034:
            position=cont
            break
        cont+=1
        #print('cont',cont)
        # cont+=1
    # print ('ad1',ad1y)
    # print('position',position)
    Sd1=Sd[cont]
    # print('sd1',Sd1)
    
    #POINT 2 
    #print('------------------Punto 2--------------------------')
    ad2y=max(aeq)
    
    cont=0
    position=0
    for ad in aeq:
    
        #print('Ad',ad)
        #print('Ad12y-Ad',abs(ad1y-ad))
        if abs(ad2y-ad)<0.02:
            position=cont
            break
        cont+=1
        #print('cont',cont)
        # cont+=1
    
    Sd2=Sd[cont]
    # print ('ad2',ad2y)
    # print('position',position)
    # print('sd2',Sd2)
    
    # POINT 3
    #print('------------------Punto 3--------------------------')
    
    position3=0
    posmax=np.where(aeq==max(aeq))
    cont=posmax[0]
    for i in range(int(posmax[0]),len(aeq)-1):       
        if aeq[i+1]-aeq[i] < 0:
            ad3y=0.87*max(aeq)
            #print('ad2-ad3y',aeq[i+1]-ad3y)
            if abs(aeq[i+1]-ad3y)<0.013:
                position3=cont
                break
            cont+=1
            Sd3=Sd[int(cont)] 
            
    posmax=np.where(aeq==max(aeq))  
    for i in range(int(posmax[0]),len(aeq)-1):
        #print('position',i)
        #print('ad2-ad1',aeq[i+1]-aeq[i])
        if aeq[i+1]-aeq[i] > 0:
            ad3y=aeq[i]
            Sd3=Sd[i]
            position3=i
            break
            
    # print('ad3y',ad3y)   
    # print('sd3',Sd3)
    
    # Punto 4
    aeq4Vector=aeq[position3:]
    #aeq4max=max(aeq4Vector)-0.20*max(aeq4Vector)
    Sa4max=max(aeq4Vector)
    for i in range(position3,len(aeq)):
        if aeq[i]>=Sa4max:
            position4=i
            break
    aeq5Vector=aeq[position4:]
    
    for i in range(position4,len(aeq)-1):
        #print(aeq[i]-aeq[i+1])
        if aeq[i]-aeq[i+1]>=0:
            Sa5max=max(aeq5Vector)
            Sa4=Sa5max-0.20*Sa5max
            for i in range(position4,len(aeq)):
                if Sa4>=aeq[i]:
                    ad4y=aeq[i]
                    Sd4=Sd[i]
        else:
            position6=i
            aeq6Vector=aeq[position6:]
            Sa5max=max(aeq6Vector)
            Sa4=Sa5max-0.20*Sa5max
            for i in range(position6,len(aeq)):
                if Sa4>=aeq[i]:
                    ad4y=aeq[i]
                    Sd4=Sd[i]
                    break
            break
        
    # COORDINATE
    coorad=[ad1y,ad2y,ad3y,ad4y]
    coorsd=[Sd1,Sd2,Sd3,Sd4]
    
    return coorad,coorsd
#%% CRITERIOS ASCE

def CriterioASCE41_17(bdI,tiempo,plasticRotationCol_IO,plasticRotationCol_LS,plasticRotationCol_CP,DtopSM,BasalSM):
    import numpy as np
    #% IO
    for i in range(bdI.shape[0]):
        #print('-------------------------------------Fila: ',i)
        for j in range(bdI.shape[1]):
            
            #print(bdI[i,j])
            if abs(bdI[i,j])>plasticRotationCol_IO:
                #print('Estado de daño IO alcanzado en fila:',i)
                #print('-----------Columna: ',j)
                break
        if abs(bdI[i,j])>plasticRotationCol_IO:
            #print('Estado de daño IO alcanzado en fila:',i)
            #print('-----------Columna: ',j)
            break 
    tiempoIO=tiempo[i]  
    #print('tiempo IO',tiempoIO)
    DtopSMIO=DtopSM[i]
    #print('DtopSM',DtopSMIO)
    BasalSMIO=BasalSM[i]
    #print('BasalSM',BasalSMIO)
    # LS
    for i in range(bdI.shape[0]):
        #print('-------------------------------------Fila: ',i)
        for j in range(bdI.shape[1]):
            #print('-----------Columna: ',j)
            #print(bdI[i,j])
            if abs(bdI[i,j])>plasticRotationCol_LS:
                #print('Estado de daño LS alcanzado en fila:',i)
                #print('-----------Columna: ',j)
                break
        if abs(bdI[i,j])>plasticRotationCol_LS:
            #print('Estado de daño LS alcanzado en fila:',i)
            #print('-----------Columna: ',j)
            break 
    tiempoLS=tiempo[i]  
    #print('tiempo LS',tiempoLS)
    DtopSMLS=DtopSM[i]
    #print('DtopSM',DtopSMLS)
    BasalSMLS=BasalSM[i]
    #print('BasalSM',BasalSMLS)
    # CP
    for i in range(bdI.shape[0]):
        #print('-------------------------------------Fila: ',i)
        for j in range(bdI.shape[1]):
            #print('-----------Columna: ',j)
            #print(bdI[i,j])
            if abs(bdI[i,j])>plasticRotationCol_CP:
                #print('Estado de daño CP alcanzado en fila:',i)
                #print('-----------Columna: ',j)
                break
        if abs(bdI[i,j])>plasticRotationCol_CP:
            #print('Estado de daño CP alcanzado en fila:',i)
            #print('-----------Columna: ',j)
            break 
    tiempoCP=tiempo[i]  
    #print('tiempo CP',tiempoCP) 
    DtopSMCP=DtopSM[i]
    #print('DtopSM',DtopSMCP)  
    BasalSMCP=BasalSM[i]
    #print('BasalSM',BasalSMCP)  
    DisplacementCD=np.array([DtopSMIO,DtopSMLS,DtopSMCP])
    VasalCD=np.array([BasalSMIO,BasalSMLS,BasalSMCP])

    return DisplacementCD,VasalCD
#%% DAMAGE STATES ASCE41-17 PLASTIC ROTATIONS
def DamageStatePlasticRotations(fc,modelo,recBeam,stirrupSpacingBeam,fi_bar_stirrup,ramales,recCol,fi_barLong,Nbarras,rutaAxialLoad,rutafolder,DtopSM_ori,BasalSM_ori,
                                rutaPlasticDeformationColumns,rutaPlasticDeformationBeams,rutaPlasticDeformationGirds):
    # Datos
    # fcE: Elastic Modulus Concrete [MPa]
    # modelo: DataFrame with information: Elements, Base, Height
    # recBeam: Beam cover [m]
    # stirrupSpacingBeam: Stirrup Spacing Beam [m]
    # fi_bar_stirrup: Diameter stirrup column [m]
    # ramales: Ramales number stirrup 
    # recCol: Column cover [m]
    # fi_barLong: Longitudinal reinforcing bar
    # Nbarras: number of bars longitudinal reinforcing
    # rutaAxialLoad: Ruta with axial load on the column. [kN]

    # Librerias
    import numpy as np
    import BibliotecaFunciones as bf
    import pandas as pd
    
    #------------------------------------------ Units
    m,mm,kN,N,Pa,kPa,MPa,g=bf.units()
    
    #%---------------------------------------------- Materials
    fytE=420*MPa
    fylE=fytE
    fcE=fc*MPa
   
    #%-------------------------------------------VIGAS
    # DATOS
    Beams_df=modelo[modelo['Type']=='Beam'][['Type','Base','Height']]
    Beams_df.reset_index(inplace=True,drop=True)
    Beams_df['Rec']=recBeam
    Beams_df['Stirrup_Spacing']=stirrupSpacingBeam
    # PROCESAMIENTO
    
    Beams_df['dBeam']=Beams_df['Height']-Beams_df['Rec']
    Beams_df['dBeam/2']=Beams_df['dBeam']/2
    # cond=[Beams_df['Stirrup_Spacing']<=Beams_df['dBeam/2'],
    #       Beams_df['Stirrup_Spacing']>Beams_df['dBeam/2']]
    # choices=['Stirrup spacing <= d/2','Stirrup spacing > d/2']
    # Beams_df['Observation']=np.select(cond,choices)
    for i in range(len(Beams_df)): 
        if Beams_df.loc[i,'Stirrup_Spacing']<=Beams_df.loc[i,'dBeam/2']:
            Beams_df.loc[i,'Observation']='Stirrup spacing <= d/2'
        else:
            Beams_df.loc[i,'Observation']='Stirrup spacing > d/2'
            Beams_df.loc[i,'plasticRotation_IO']=0.0015
            Beams_df.loc[i,'plasticRotation_LS']=0.005
            Beams_df.loc[i,'plasticRotation_CP']=0.01
        
    #print('Beams \n',Beams_df)
    
    #%%-------------------------------------------GIRDS
    # DATOS
    Girds_df=modelo[modelo['Type']=='Gird'][['Type','Base','Height']]
    Girds_df.reset_index(inplace=True,drop=True)
    Girds_df['Rec']=recBeam
    Girds_df['Stirrup_Spacing']=stirrupSpacingBeam
    # PROCESAMIENTO
    
    Girds_df['dBeam']=Girds_df['Height']-Girds_df['Rec']
    Girds_df['dBeam/2']=Girds_df['dBeam']/2
    # cond=[Beams_df['Stirrup_Spacing']<=Beams_df['dBeam/2'],
    #       Beams_df['Stirrup_Spacing']>Beams_df['dBeam/2']]
    # choices=['Stirrup spacing <= d/2','Stirrup spacing > d/2']
    # Beams_df['Observation']=np.select(cond,choices)
    for i in range(len(Girds_df)): 
        if Girds_df.loc[i,'Stirrup_Spacing']<=Girds_df.loc[i,'dBeam/2']:
            Girds_df.loc[i,'Observation']='Stirrup spacing <= d/2'
        else:
            Girds_df.loc[i,'Observation']='Stirrup spacing > d/2'
            Girds_df.loc[i,'plasticRotation_IO']=0.0015
            Girds_df.loc[i,'plasticRotation_LS']=0.005
            Girds_df.loc[i,'plasticRotation_CP']=0.01
        
    #print('Girds \n',Girds_df)
    
    #%%-------------------------------------------Datos Columnas


    # Transversal
    Columns_df=modelo[modelo['Type']=='Column'][['Type','Base','Height']]
    Columns_df.reset_index(inplace=True,drop=True)
    Columns_df['fi_Stirrup']=fi_bar_stirrup
    Columns_df['n_ramales']=ramales
    Columns_df['As_fiStirrup']=np.pi*((Columns_df['fi_Stirrup'])**2)/4
    Columns_df['Rec']=recCol
    Columns_df['dCol']=Columns_df['Height']-Columns_df['Rec']
    Columns_df['Av']=Columns_df['n_ramales']*Columns_df['As_fiStirrup']
    Columns_df['rho_Tranverse']=Columns_df['Av']/(Columns_df['Base']*Columns_df['dCol'])
    
    # Longitudinal
    Columns_df['fi_long']=fi_barLong
    Columns_df['Nbarras']=Nbarras
    Columns_df['As_Barralong']=np.pi*((Columns_df['fi_long'])**2)/4
    Columns_df['As']=Columns_df['Nbarras']*Columns_df['As_Barralong']
    Columns_df['rho_long']=Columns_df['As']/(Columns_df['Base']*Columns_df['Height'])
    AxialLoad=np.genfromtxt(rutaAxialLoad)
    AxialLoad=AxialLoad[9,1:]
    Columns_df['NuD']=AxialLoad
    #Columns_df['NuD']=50
    
    #------------------------------ COLUMNAS TABLA 10-8 ASCE 41-17
    #---Coeficientes Tabla 10-8
    Columns_df['a']=(1*Columns_df['rho_Tranverse']*fytE)/(8*Columns_df['rho_long']*fylE)
    Columns_df['Ag']=Columns_df['Base']*Columns_df['Height']
    Columns_df['NuD/(Ag*fcE)']=Columns_df['NuD']/(Columns_df['Ag']*fcE)
    # Condiciones ASCE41 
    x=[0.5,0.7]
    y=[0.5,0.0]
    # for i in range(len(Columns_df)):
    #     if Columns_df.loc[i,'NuD/(Ag*fcE)']>0.5 and Columns_df.loc[i,'NuD/(Ag*fcE)']<0.7:
    #         Columns_df.loc[i,'NuD/(Ag*fcE)']=np.interp(Columns_df.loc[i,'NuD/(Ag*fcE)'],x,y)
    #     elif Columns_df.loc[i,'NuD/(Ag*fcE)']>0.7:
    #         Columns_df.loc[i,'NuD/(Ag*fcE)']=0
    
    for i in range(len(Columns_df)):
        Columns_df.loc[i,'b']=abs((0.012-0.085*Columns_df.loc[i,'NuD/(Ag*fcE)'])+12*Columns_df.loc[i,'rho_Tranverse'])
        if Columns_df.loc[i,'b']<Columns_df.loc[i,'a']:
            Columns_df.loc[i,'b']=Columns_df.loc[i,'a']
              
    # Plastic Rotation Angle Columns
    Columns_df['plasticRotation_IO']=0.0
    Columns_df['plasticRotation_LS']=0.5*Columns_df['b']
    Columns_df['plasticRotation_CP']=0.7*Columns_df['b']
    
    #print('Columns \n',Columns_df)
    #%% ROTACIONES COLUMNAS
    # Deformaciones
    # Column0: Time
    # Column1: Axial
    # Column2: z-rotation at I
    # Column3: z-rotation at J
    # Column4: y-rotation at I
    # Column5: y-rotation at J
    # Column6: Torsion
    #cont=1
    # BasicDeformationColumns=rutafolder+'\\'+'DataInelasticInfillModel'+str(cont)+'\\'+'Rotations'+'\\'+'Columns'+'\\'+'BasicDeformationColumns.out'
    # BasicDeformationColumns=np.genfromtxt(BasicDeformationColumns) 
    
    #PlasticDeformationColumns=rutafolder+'\\'+'DataInelasticInfillModel'+str(cont)+'\\'+'Rotations'+'\\'+'Columns'+'\\'+'PlasticDeformationColumns.out'
    PlasticDeformationColumns=np.genfromtxt(rutaPlasticDeformationColumns) 
    
    tiempo=PlasticDeformationColumns[:,0]
    
    # Nudo Inicial
    bdI=PlasticDeformationColumns[:,2:-1:6]
    
    for i in range(len(Columns_df)):
        plasticRotationCol_IO=Columns_df.loc[i,'plasticRotation_IO']
        plasticRotationCol_LS=Columns_df.loc[i,'plasticRotation_LS']
        plasticRotationCol_CP=Columns_df.loc[i,'plasticRotation_CP']  
        DisplacementCD_ColNI,VasalCD_ColNI=bf.CriterioASCE41_17(bdI,tiempo,plasticRotationCol_IO,plasticRotationCol_LS,plasticRotationCol_CP,DtopSM_ori,BasalSM_ori)
        Columns_df.loc[i,'Displacement_Ni_IO']=DisplacementCD_ColNI[0]
        Columns_df.loc[i,'Displacement_Ni_LS']=DisplacementCD_ColNI[1]
        Columns_df.loc[i,'Displacement_Ni_CP']=DisplacementCD_ColNI[2]
        Columns_df.loc[i,'Basal_Ni_IO']=VasalCD_ColNI[0]
        Columns_df.loc[i,'Basal_Ni_LS']=VasalCD_ColNI[1]
        Columns_df.loc[i,'Basal_Ni_CP']=VasalCD_ColNI[2]
        
    # Nudo Final
    # bdI=BasicDeformationColumns[:,3:-1:6]
    bdICol=PlasticDeformationColumns[:,3:-1:6]
    
    for i in range(len(Columns_df)):
        plasticRotationCol_IO=Columns_df.loc[i,'plasticRotation_IO']
        plasticRotationCol_LS=Columns_df.loc[i,'plasticRotation_LS']
        plasticRotationCol_CP=Columns_df.loc[i,'plasticRotation_CP'] 
        DisplacementCD_ColNJ,VasalCD_ColNJ=bf.CriterioASCE41_17(bdICol,tiempo,plasticRotationCol_IO,plasticRotationCol_LS,plasticRotationCol_CP,DtopSM_ori,BasalSM_ori)
        Columns_df.loc[i,'Displacement_Nj_IO']=DisplacementCD_ColNJ[0]
        Columns_df.loc[i,'Displacement_Nj_LS']=DisplacementCD_ColNJ[1]
        Columns_df.loc[i,'Displacement_Nj_CP']=DisplacementCD_ColNJ[2]
        Columns_df.loc[i,'Basal_Nj_IO']=VasalCD_ColNJ[0]
        Columns_df.loc[i,'Basal_Nj_LS']=VasalCD_ColNJ[1]
        Columns_df.loc[i,'Basal_Nj_CP']=VasalCD_ColNJ[2]
        
    #%% ROTACIONES VIGAS 
    
    # BasicDeformationBeams=rutafolder+'\\'+'DataInelasticInfillModel'+str(cont)+'\\'+'Rotations'+'\\'+'Beams'+'\\'+'BasicDeformationBeams.out'
    # BasicDeformationBeams=np.genfromtxt(BasicDeformationBeams) 
    
    #rutaPlasticDeformationBeams=rutafolder+'\\'+'DataInelasticInfillModel'+str(cont)+'\\'+'Rotations'+'\\'+'Beams'+'\\'+'PlasticDeformationBeams.out'
    PlasticDeformationBeams=np.genfromtxt(rutaPlasticDeformationBeams) 
    
    
    # Nudo Inicial
    #bdIBeam=BasicDeformationBeams[:,2:-1:6]
    bdIBeam=PlasticDeformationBeams[:,2:-1:6]
    for i in range(len(Beams_df)):
    #for i in range(1):
        plasticRotationBeam_IO=Beams_df.loc[i,'plasticRotation_IO']
        plasticRotationBeam_LS=Beams_df.loc[i,'plasticRotation_LS']
        plasticRotationBeam_CP=Beams_df.loc[i,'plasticRotation_CP'] 
        DisplacementCD_BeamNI,VasalCD_BeamNI=bf.CriterioASCE41_17(bdIBeam,tiempo,plasticRotationBeam_IO,plasticRotationBeam_LS,plasticRotationBeam_CP,DtopSM_ori,BasalSM_ori)
        Beams_df.loc[i,'Displacement_Ni_IO']=DisplacementCD_BeamNI[0]
        Beams_df.loc[i,'Displacement_Ni_LS']=DisplacementCD_BeamNI[1]
        Beams_df.loc[i,'Displacement_Ni_CP']=DisplacementCD_BeamNI[2]
        Beams_df.loc[i,'Basal_Ni_IO']=VasalCD_BeamNI[0]
        Beams_df.loc[i,'Basal_Ni_LS']=VasalCD_BeamNI[1]
        Beams_df.loc[i,'Basal_Ni_CP']=VasalCD_BeamNI[2]
    # Nudo Final
    bdIBeam=PlasticDeformationBeams[:,3:-1:6]
    for i in range(len(Beams_df)):
    #for i in range(1):
        plasticRotationBeam_IO=Beams_df.loc[i,'plasticRotation_IO']
        plasticRotationBeam_LS=Beams_df.loc[i,'plasticRotation_LS']
        plasticRotationBeam_CP=Beams_df.loc[i,'plasticRotation_CP'] 
        DisplacementCD_BeamNJ,VasalCD_BeamNJ=bf.CriterioASCE41_17(bdIBeam,tiempo,plasticRotationBeam_IO,plasticRotationBeam_LS,plasticRotationBeam_CP,DtopSM_ori,BasalSM_ori)
        Beams_df.loc[i,'Displacement_Nj_IO']=DisplacementCD_BeamNJ[0]
        Beams_df.loc[i,'Displacement_Nj_LS']=DisplacementCD_BeamNJ[1]
        Beams_df.loc[i,'Displacement_Nj_CP']=DisplacementCD_BeamNJ[2]
        Beams_df.loc[i,'Basal_Nj_IO']=VasalCD_BeamNJ[0]
        Beams_df.loc[i,'Basal_Nj_LS']=VasalCD_BeamNJ[1]
        Beams_df.loc[i,'Basal_Nj_CP']=VasalCD_BeamNJ[2]
    #%% ROTACIONES GIRDS 
    
    #BasicDeformationGirds=rutafolder+'\\'+'DataInelasticInfillModel'+str(cont)+'\\'+'Rotations'+'\\'+'Girds'+'\\'+'BasicDeformationGirds.out'
    #BasicDeformationGirds=np.genfromtxt(BasicDeformationGirds) 
    
    #rutaPlasticDeformationGirds=rutafolder+'\\'+'DataInelasticInfillModel'+str(cont)+'\\'+'Rotations'+'\\'+'Girds'+'\\'+'PlasticDeformationGirds.out'
    PlasticDeformationGirds=np.genfromtxt(rutaPlasticDeformationGirds) 
    
    # Nudo Inicial
    #bdIBeam=BasicDeformationBeams[:,2:-1:6]
    bdIGirds=PlasticDeformationGirds[:,2:-1:6]
    for i in range(len(Girds_df)):
    #for i in range(1):
        plasticRotationGirds_IO=Girds_df.loc[i,'plasticRotation_IO']
        plasticRotationGirds_LS=Girds_df.loc[i,'plasticRotation_LS']
        plasticRotationGirds_CP=Girds_df.loc[i,'plasticRotation_CP'] 
        DisplacementCD_GirdsNI,VasalCD_GirdsNI=bf.CriterioASCE41_17(bdIGirds,tiempo,plasticRotationGirds_IO,plasticRotationGirds_LS,plasticRotationGirds_CP,DtopSM_ori,BasalSM_ori)
        Girds_df.loc[i,'Displacement_Ni_IO']=DisplacementCD_GirdsNI[0]
        Girds_df.loc[i,'Displacement_Ni_LS']=DisplacementCD_GirdsNI[1]
        Girds_df.loc[i,'Displacement_Ni_CP']=DisplacementCD_GirdsNI[2]
        Girds_df.loc[i,'Basal_Ni_IO']=VasalCD_GirdsNI[0]
        Girds_df.loc[i,'Basal_Ni_LS']=VasalCD_GirdsNI[1]
        Girds_df.loc[i,'Basal_Ni_CP']=VasalCD_GirdsNI[2]
    # Nudo Final
    bdIGirds=PlasticDeformationGirds[:,3:-1:6]
    for i in range(len(Girds_df)):
    #for i in range(1):
        plasticRotationGirds_IO=Girds_df.loc[i,'plasticRotation_IO']
        plasticRotationGirds_LS=Girds_df.loc[i,'plasticRotation_LS']
        plasticRotationGirds_CP=Girds_df.loc[i,'plasticRotation_CP'] 
        DisplacementCD_GirdsNJ,VasalCD_GirdsNJ=bf.CriterioASCE41_17(bdIGirds,tiempo,plasticRotationGirds_IO,plasticRotationGirds_LS,plasticRotationGirds_CP,DtopSM_ori,BasalSM_ori)
        Girds_df.loc[i,'Displacement_Nj_IO']=DisplacementCD_GirdsNJ[0]
        Girds_df.loc[i,'Displacement_Nj_LS']=DisplacementCD_GirdsNJ[1]
        Girds_df.loc[i,'Displacement_Nj_CP']=DisplacementCD_GirdsNJ[2]
        Girds_df.loc[i,'Basal_Nj_IO']=VasalCD_GirdsNJ[0]
        Girds_df.loc[i,'Basal_Nj_LS']=VasalCD_GirdsNJ[1]
        Girds_df.loc[i,'Basal_Nj_CP']=VasalCD_GirdsNJ[2]
    #%% RESUMEN
    #----------------------------- IO
    # Columns
    #Nudo Inicial
    MinColumnsNiIO=Columns_df.loc[Columns_df['Displacement_Ni_IO']==Columns_df['Displacement_Ni_IO'].min()][['Type','Displacement_Ni_IO','Basal_Ni_IO']]
    MinColumnsNiIO.columns=['Type','Displacement_IO','Basal_IO']
    #Nudo Final
    MinColumnsNjIO=Columns_df.loc[Columns_df['Displacement_Nj_IO']==Columns_df['Displacement_Nj_IO'].min()][['Type','Displacement_Nj_IO','Basal_Nj_IO']]
    MinColumnsNjIO.columns=['Type','Displacement_IO','Basal_IO']
    
    # Beams
    # Nudo Inicial
    MinBeamsNiIO=Beams_df.loc[Beams_df['Displacement_Ni_IO']==Beams_df['Displacement_Ni_IO'].min()][['Type','Displacement_Ni_IO','Basal_Ni_IO']]
    MinBeamsNiIO.columns=['Type','Displacement_IO','Basal_IO']
    # Nudo Final
    MinBeamsNjIO=Beams_df.loc[Beams_df['Displacement_Nj_IO']==Beams_df['Displacement_Nj_IO'].min()][['Type','Displacement_Nj_IO','Basal_Nj_IO']]
    MinBeamsNjIO.columns=['Type','Displacement_IO','Basal_IO']
    
    # Girds
    # Nudo Inicial
    MinGirdsNiIO=Girds_df.loc[Girds_df['Displacement_Ni_IO']==Girds_df['Displacement_Ni_IO'].min()][['Type','Displacement_Ni_IO','Basal_Ni_IO']]
    MinGirdsNiIO.columns=['Type','Displacement_IO','Basal_IO']
    # Nudo Final
    MinGirdsNjIO=Girds_df.loc[Girds_df['Displacement_Nj_IO']==Girds_df['Displacement_Nj_IO'].min()][['Type','Displacement_Nj_IO','Basal_Nj_IO']]
    MinGirdsNjIO.columns=['Type','Displacement_IO','Basal_IO']
    
    ResultadosIO=pd.concat([MinColumnsNiIO,MinColumnsNjIO,MinBeamsNiIO,MinBeamsNjIO,MinGirdsNiIO,MinGirdsNjIO],ignore_index=True)
    
    #----------------------------- LS
    # Columns
    #Nudo Inicial
    MinColumnsNiLS=Columns_df.loc[Columns_df['Displacement_Ni_LS']==Columns_df['Displacement_Ni_LS'].min()][['Type','Displacement_Ni_LS','Basal_Ni_LS']]
    MinColumnsNiLS.columns=['Type','Displacement_LS','Basal_LS']
    #Nudo Final
    MinColumnsNjLS=Columns_df.loc[Columns_df['Displacement_Nj_LS']==Columns_df['Displacement_Nj_LS'].min()][['Type','Displacement_Nj_LS','Basal_Nj_LS']]
    MinColumnsNjLS.columns=['Type','Displacement_LS','Basal_LS']
    
    # Beams
    # Nudo Inicial
    MinBeamsNiLS=Beams_df.loc[Beams_df['Displacement_Ni_LS']==Beams_df['Displacement_Ni_LS'].min()][['Type','Displacement_Ni_LS','Basal_Ni_LS']]
    MinBeamsNiLS.columns=['Type','Displacement_LS','Basal_LS']
    # Nudo Final
    MinBeamsNjLS=Beams_df.loc[Beams_df['Displacement_Nj_LS']==Beams_df['Displacement_Nj_LS'].min()][['Type','Displacement_Nj_LS','Basal_Nj_LS']]
    MinBeamsNjLS.columns=['Type','Displacement_LS','Basal_LS']
    
    # Girds
    # Nudo Inicial
    MinGirdsNiLS=Girds_df.loc[Girds_df['Displacement_Ni_LS']==Girds_df['Displacement_Ni_LS'].min()][['Type','Displacement_Ni_LS','Basal_Ni_LS']]
    MinGirdsNiLS.columns=['Type','Displacement_LS','Basal_LS']
    # Nudo Final
    MinGirdsNjLS=Girds_df.loc[Girds_df['Displacement_Nj_LS']==Girds_df['Displacement_Nj_LS'].min()][['Type','Displacement_Nj_LS','Basal_Nj_LS']]
    MinGirdsNjLS.columns=['Type','Displacement_LS','Basal_LS']
    
    ResultadosLS=pd.concat([MinColumnsNiLS,MinColumnsNjLS,MinBeamsNiLS,MinBeamsNjLS,MinGirdsNiLS,MinGirdsNjLS],ignore_index=True)
    
    #--------------------------------CP
    # Columns
    #Nudo Inicial
    MinColumnsNiCP=Columns_df.loc[Columns_df['Displacement_Ni_CP']==Columns_df['Displacement_Ni_CP'].min()][['Type','Displacement_Ni_CP','Basal_Ni_CP']]
    MinColumnsNiCP.columns=['Type','Displacement_CP','Basal_CP']
    #Nudo Final
    MinColumnsNjCP=Columns_df.loc[Columns_df['Displacement_Nj_CP']==Columns_df['Displacement_Nj_CP'].min()][['Type','Displacement_Nj_CP','Basal_Nj_CP']]
    MinColumnsNjCP.columns=['Type','Displacement_CP','Basal_CP']
    
    # Beams
    # Nudo Inicial
    MinBeamsNiCP=Beams_df.loc[Beams_df['Displacement_Ni_CP']==Beams_df['Displacement_Ni_CP'].min()][['Type','Displacement_Ni_CP','Basal_Ni_CP']]
    MinBeamsNiCP.columns=['Type','Displacement_CP','Basal_CP']
    # Nudo Final
    MinBeamsNjCP=Beams_df.loc[Beams_df['Displacement_Nj_CP']==Beams_df['Displacement_Nj_CP'].min()][['Type','Displacement_Nj_CP','Basal_Nj_CP']]
    MinBeamsNjCP.columns=['Type','Displacement_CP','Basal_CP']
    
    # Girds
    # Nudo Inicial
    MinGirdsNiCP=Girds_df.loc[Girds_df['Displacement_Ni_CP']==Girds_df['Displacement_Ni_CP'].min()][['Type','Displacement_Ni_CP','Basal_Ni_CP']]
    MinGirdsNiCP.columns=['Type','Displacement_CP','Basal_CP']
    # Nudo Final
    MinGirdsNjCP=Girds_df.loc[Girds_df['Displacement_Nj_CP']==Girds_df['Displacement_Nj_CP'].min()][['Type','Displacement_Nj_CP','Basal_Nj_CP']]
    MinGirdsNjCP.columns=['Type','Displacement_CP','Basal_CP']
    
    ResultadosCP=pd.concat([MinColumnsNiCP,MinColumnsNjCP,MinBeamsNiCP,MinBeamsNjCP,MinGirdsNiCP,MinGirdsNjCP],ignore_index=True)
    
    #%%---------------------------RESUMEN
    # IO
    ResultadosIO=ResultadosIO.sort_values('Displacement_IO',ascending=True)
    ResumenIO=ResultadosIO.loc[[0]]
    # LS
    ResultadosLS=ResultadosLS.sort_values('Displacement_LS',ascending=True)
    ResumenLS=ResultadosLS.loc[[0]]
    # CP
    ResultadosCP=ResultadosCP.sort_values('Displacement_CP',ascending=True)
    ResumenCP=ResultadosCP.loc[[0]]
    Resumen1=pd.concat([ResumenIO,ResumenLS,ResumenCP],axis=1)
    Resumen1.columns=['ElementFail_IO','Displacement_IO','Basal_IO','ElementFail_LS','Displacement_LS','Basal_LS',
                     'ElementFail_CP','Displacement_CP','Basal_CP']
    DamageStatePlasRotation=Resumen1.transpose()
    
    #print(DamageStatePlasRotation)
    return Columns_df,DamageStatePlasRotation
#%% SHEAR DAMAGE STATES
def ShearDS(rutaShear,rutaLocalForces,BasalSM_ori,DtopSM_ori):
    import pandas as pd
    shear=pd.read_csv(rutaShear,sep=' ',header=None)
    shear.columns=('V_cryy','V_ccyy','V_cyy','V_resyy','V_crzz','V_cczz','V_czz','V_reszz',
                   'gamm_cryy','gamm_pkyy','gamm_uyy','gamm_myy','gamm_crzz','gamm_pkzz','gamm_uzz','gamm_mzz')
    
    #% Local Forces Columns
    LocalForces=pd.read_csv(rutaLocalForces,sep=' ',header=None)
    
    #% Tratamiento
    Ncol=LocalForces.shape[1]
    sheardemand=LocalForces.iloc[:,range(2,Ncol,12)]
    Ncol2=sheardemand.shape[1]
    sheardemand_np=sheardemand.to_numpy()
    step=list(LocalForces.iloc[:,0])
    
    #% DAMAGE STATES SHEAR
    #LS 0.75Vmax
    #CP Vmax
    #C 0.50Vmax
    #---------------------------------------------
    # PRIMERA COLUMNA
    shearbyColumn=shear.loc[0]
    Basal_VLS=pd.DataFrame()
    Dtop_VLS=pd.DataFrame()
    #-----------------------------------------------------LS
    fil=sheardemand_np.shape[0]
    col=sheardemand_np.shape[1]
    for j in range(col):
        shearbyColumn=shear.loc[j]
        Vmax=shearbyColumn['V_cyy']
        #print('Vmax col '+str(j),Vmax)
        V_LS=0.75*Vmax
        #print('V_LS '+str(j),V_LS)
        for i in range(fil):
            #print(i)
            if abs(sheardemand_np[i,j])>=V_LS:
                #print('Position LS col '+str(j),i)
                y=[sheardemand_np[i-1,j],sheardemand_np[i,j]]
                x=[step[i-1],step[i]]
                print(f'Falla por cortante columna {i}')
                break     
                #print('No hay falla por cortante')
            posicion_VLS=i-1
            
        if BasalSM_ori[i]>=0:
            Basal_VLS['columna '+str(j)]=[BasalSM_ori[i]]
            Dtop_VLS['columna '+str(j)]=[DtopSM_ori[i]]
        else:
            basal=[]
            desp=[]
            for k in range(len(BasalSM_ori)):
                if BasalSM_ori[k] >= 0:
                    basal.append(BasalSM_ori[k])
                    desp.append(DtopSM_ori[k])
            Basal_VLS['columna '+str(j)]=[basal[-1]]
            Dtop_VLS['columna '+str(j)]=[desp[-1]]
               
    Basal_V_LS=Basal_VLS.min()[0]
    Dtop_V_LS=Dtop_VLS.min()[0]
    #print(f'Cortante Basal LS {Basal_V_LS}')
    #print(f'Cortante Basal LS {Dtop_V_LS}')
    
    #print('-----------------------------------------COLAPSE PREVENTION')
    Basal_VCP=pd.DataFrame()
    Dtop_VCP=pd.DataFrame()
    fil=sheardemand_np.shape[0]
    col=sheardemand_np.shape[1]
    for j in range(col):
        shearbyColumn=shear.loc[j]
        Vmax=shearbyColumn['V_cyy']
        #print('Vmax col '+str(j),Vmax)
        V_CP=Vmax
        #print('V_CP '+str(j),V_CP)
        for i in range(posicion_VLS,fil):
            #print(i)
            if abs(sheardemand_np[i,j])>=V_CP:
               # print('Position CP col '+str(j),i)
                y=[sheardemand_np[i-1,j],sheardemand_np[i,j]]
                x=[step[i-1],step[i]]
                print(f'Falla por cortante columna CP {i}')
                break     
                #print('No hay falla por cortante')
            posicion_VCP=i-1
            
        if BasalSM_ori[i]>=0:
            Basal_VCP['columna '+str(j)]=[BasalSM_ori[i]]
            Dtop_VCP['columna '+str(j)]=[DtopSM_ori[i]]
        else:
            basal=[]
            desp=[]
            for k in range(len(BasalSM_ori)):
                if BasalSM_ori[k] >= 0:
                    basal.append(BasalSM_ori[k])
                    desp.append(DtopSM_ori[k])
            Basal_VCP['columna '+str(j)]=[basal[-1]]
            Dtop_VCP['columna '+str(j)]=[desp[-1]]
                    
    Basal_V_CP=Basal_VCP.min()[0]
    Dtop_V_CP=Dtop_VCP.min()[0]
    #print(f'Cortante Basal CP {Basal_V_CP}')
    #print(f'Cortante Basal CP {Dtop_V_CP}')        
    #-----------------------COLAPSE
    #print('-----------------------------------------COLAPSE')
    Basal_VC=pd.DataFrame()
    Dtop_VC=pd.DataFrame()
    fil=sheardemand_np.shape[0]
    col=sheardemand_np.shape[1]
    for j in range(col):
        shearbyColumn=shear.loc[j]
        Vmax=shearbyColumn['V_cyy']
        #print('Vmax col '+str(j),Vmax)
        V_C=0.20*Vmax
        #print('V_C '+str(j),V_C)
        for i in range(posicion_VCP,fil):
            #print(i)
            if abs(sheardemand_np[i,j])>=V_C:
                #print('Position C col '+str(j),i)
                y=[sheardemand_np[i-1,j],sheardemand_np[i,j]]
                x=[step[i-1],step[i]]
                print(f'Falla por cortante columna C {i}')
                break     
                #print('No hay falla por cortante')
        if BasalSM_ori[i]>=0:
            Basal_VC['columna '+str(j)]=[BasalSM_ori[i]]
            Dtop_VC['columna '+str(j)]=[DtopSM_ori[i]]
        else:
            basal=[]
            desp=[]
            for k in range(len(BasalSM_ori)):
                if BasalSM_ori[k] >= 0:
                    basal.append(BasalSM_ori[k])
                    desp.append(DtopSM_ori[k])
            Basal_VC['columna '+str(j)]=[basal[-1]]
            Dtop_VC['columna '+str(j)]=[desp[-1]]
            
    Basal_V_C=Basal_VC.min()[0]
    Dtop_V_C=Dtop_VC.min()[0]
    #print(f'Cortante Basal C {Basal_V_C}')
    #print(f'Cortante Basal C {Dtop_V_C}')
    Shear_DS=pd.DataFrame()
    Shear_DS['Basal_V_LS'],Shear_DS['Dtop_V_LS']=[Basal_V_LS],[Dtop_V_LS]
    Shear_DS['Basal_V_CP'],Shear_DS['Dtop_V_CP']=[Basal_V_CP],[Dtop_V_CP]
    Shear_DS['Basal_V_C'],Shear_DS['Dtop_V_C']=[Basal_V_C],[Dtop_V_C]
    Shear_DS=Shear_DS.transpose()
    return Shear_DS