import PySimpleGUI as sg
import pandas as pd
import numpy as np
import matplotlib.pyplot as pl
from scipy import integrate
from scipy import interpolate
import math
  
# Add some color
# to the window
sg.theme('SandyBeach')     
  
# Very basic window.
# Return values using
# automatic-numbered keys
layout = [
    [sg.Text('Please enter the required parameters')],
    [sg.Text('Enter project name', size=(20,1)), sg.Input(key='-IN-')],
    [sg.Text('Drag polar (Mach number vs Cd as a .csv file)', size =(40, 1)), sg.Input(key='-IN1-'), sg.FileBrowse(file_types=(('Excel Files','*csv'),))],
    [sg.Text('Motor parameters',size=(100,1), justification='c')],
    [sg.Text('Motor impuls [in Ns]', size=(20,1)), sg.Input(key='-IN2-', size=(10,1)),sg.Text('Total motor weight [in kg]'), sg.Input(key='-IN3-',size=(10,1)), sg.Text('Propelant weignt [in kg]'), sg.Input(key='-IN4-',size=(10,1))],
    [sg.Text('Trust curve (time vs Thrust as a .cvs file)', size =(40, 1)), sg.Input(key='-IN5-'), sg.FileBrowse(file_types=(('Excel Files','*csv'),))],
    [sg.Text('Rocket parameters',size=(100,1), justification='c')],
    [sg.Text('Rocket weight (WITHOUT MOTOR [in kg])',size=(40,1)), sg.Input (key='-IN6-', size=(10,1)),sg.Text('Rocket diameter [m]', size =(20, 1)), sg.Input(key='-IN7-', size=(10,1)) ],
    [sg.Text('Initial conditions', size=(100,1), justification='c')],
    [sg.Text('Starting angle [deg]', size=(20,1)), sg.Input(key='-IN8-', size=(10,1)), sg.Text('Launch Rod Lenght [in m]'), sg.Input (key='-IN9-',size=(10,1))],
    [sg.Submit(), sg.Cancel()]
]

Memory_Allocation = 30000 # Memory storage

window = sg.Window('Trijectory Rocket Simulator', layout)
event, values = window.read()
window.close()

DragCVS=pd.read_csv(values['-IN1-'])
M_DragCVS_1=pd.DataFrame(DragCVS).to_numpy()
M_DragCVS=M_DragCVS_1[:,0]
Cd_DragCVS=M_DragCVS_1[:,1]

ThrustCVS=pd.read_csv(values['-IN5-'])
t_ThrustCVS_1=pd.DataFrame(ThrustCVS).to_numpy()
time_thrust=t_ThrustCVS_1[:,0]
Thrust=t_ThrustCVS_1[:,1]

Launch_Rod_Lenght=float(values['-IN9-'])
#Mass Rocket Characteristic
Mass_Rocket_Dry=float(values['-IN6-']) # Mass of empty rocket without the motor [kg] - GUI Import
A=(float(values['-IN7-'])**2)*np.pi
#Engine rocket characteristic
Propelent_Weight=float(values['-IN4-']) #Mass of rocket propelent [kg] - GUI Import
Engine_Initial_Mass=float(values['-IN3-']) #Initial mass of rocket motor [kg] - GUI Import
Total_Impuls=float(values['-IN2-']) #Total impuls of rocket motor [Ns] - GUI Import

#Initial conditions
Theta_0=float(values['-IN8-']) #Elevation angle [deg] - Import from GUI


#Initilization
t=np.linspace(0,0,num=Memory_Allocation)
x=np.linspace(0,0,num=Memory_Allocation)
z=np.linspace(0,0,num=Memory_Allocation)
Theta=np.linspace(0,0,num=Memory_Allocation)

Rho=np.linspace(0,0,num=Memory_Allocation)
Temp=np.linspace(0,0,num=Memory_Allocation)
Mass=np.linspace(0,0,num=Memory_Allocation)
Fn=np.linspace(0,0,num=Memory_Allocation)
M=np.linspace(0,0,num=Memory_Allocation)
Drag=np.linspace(0,0,num=Memory_Allocation)
Fx=np.linspace(0,0,num=Memory_Allocation)

Fz=np.linspace(0,0,num=Memory_Allocation)
Ax=np.linspace(0,0,num=Memory_Allocation)
Az=np.linspace(0,0,num=Memory_Allocation)
Acc=np.linspace(0,0,num=Memory_Allocation)
Vx=np.linspace(0,0,num=Memory_Allocation)
Vz=np.linspace(0,0,num=Memory_Allocation)

V=np.linspace(0,0,num=Memory_Allocation)
Distance_x=np.linspace(0,0,num=Memory_Allocation)
Distance_z=np.linspace(0,0,num=Memory_Allocation)
Distance=np.linspace(0,0,num=Memory_Allocation)

Vx[0]=0
Vz[0]=0
V[0]=0
x[0]=0
z[0]=0.1
Mass[0]=Engine_Initial_Mass+Mass_Rocket_Dry
Theta[0]=Theta_0
n=0
rho0=1.225
kappa=1.4
R=287
Rho[0]=rho0
Temp0=300
Temp[0]=Temp0
Theta_rad=np.linspace(0,0,num=Memory_Allocation)
Theta_rad=np.radians(Theta)
slices=np.linspace(0,0,num=Memory_Allocation)
while z[n] >=0 :
    n=n+1
    t[n]=(n)*0.01
    new_Thrust=np.interp(t,time_thrust,Thrust)
    Temp[n]=Temp0-0.0065*z[n]
    M[n]=V[n-1]/np.sqrt(kappa*R*Temp[n])
    new_Cd=np.interp(M,M_DragCVS,Cd_DragCVS)
    slices[n]=integrate.trapz(new_Thrust,t,dx=0.01)
    Mass[n]=-(slices[n]/Total_Impuls*Propelent_Weight)+Engine_Initial_Mass+Mass_Rocket_Dry
    Rho[n]=rho0*(1-z[n]/44300)**(4.25588)
    Drag[n]=0.5*new_Cd[n]*Rho[n]*A*(V[n-1])**2
    if Distance[n-1] <= Launch_Rod_Lenght :
        Fn[n]=Mass[n]*9.81*math.cos(Theta[n-1]*math.pi/180)*1
    else :
        Fn[n]=0
    Fx[n]=new_Thrust[n]*math.cos(Theta[n-1]*math.pi/180)-Drag[n]*math.cos(Theta[n-1]*math.pi/180)-Fn[n]*math.sin(Theta[n-1]*math.pi/180)
    Fz[n]=new_Thrust[n]*math.sin(Theta[n-1]*math.pi/180)-(Mass[n]*9.81)-Drag[n]*math.sin(Theta[n-1]*math.pi/180)+Fn[n]*math.cos(Theta[n-1]*math.pi/180)
    Ax[n]=Fx[n]/Mass[n]
    Az[n]=Fz[n]/Mass[n]
    Acc[n]=math.sqrt(Ax[n]**2+Az[n]**2)
    Vx[n]=Vx[n-1]+Ax[n]*0.01
    Vz[n]=Vz[n-1]+Az[n]*0.01
    V[n]=math.sqrt(Vx[n]**2+Vz[n]**2)
    x[n]=x[n-1]+Vx[n]*0.01
    z[n]=z[n-1]+Vz[n]*0.01
    Distance_x[n]=Distance_x[n-1]+np.abs(Vx[n]*0.01)
    Distance_z[n]=Distance_z[n-1]+np.abs(Vz[n]*0.01)
    Distance[n]=np.sqrt(Distance_x[n]**2+Distance_z[n]**2)
    Theta[n]=math.atan(Vz[n]/Vx[n])
    Theta[n]=np.rad2deg(Theta[n])
    

pl.plot(t,z)
pl.show()

