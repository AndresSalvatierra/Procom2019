import numpy as np

#Training

Nsymb=500
iteraciones=1000
Ntap=11 #cantidad de coeficientes
K=(Ntap-1)/2
Pr=1 #Potencia de ruido de la senial recibida
delta=0.115
canal=[0.05,-0,063 ,0.088, -0.126, -0.25, 0.9047, 0.25, 0, 0.126, 0.028, 0.088]
sigma=0.01

signal_training = 2*(np.random.uniform(-1,1,Nsymb)>0.0)-1 #Training de 500 valores entre +-1
noise=np.sqrt(sigma)*np.random.normal(0,1,Nsymb) #Ruido gausiano(Media 0,desviacion estandar de la dist.normal 1,500 elemts)
xk=signal_training+noise
mem_filtro= np.zeros(Ntap-1) #Memoria del filtro FIR
yk=np.zeros(Nsymb) #La entrada al adaptive MSE
coef=np.zeros(Ntap) #Coeficientes del filtro FIR del equalizador adaptativo


for i in range(0,Nsymb): #Ver el largo del arreglo

    y=xk[i]*canal[0]           #FIR
    for k in range(1,Ntap):
        y=y+mem_filtro[k]*canal[k]

    for j in range(len(mem_filtro),0,-1):
        mem_filtro[j]=mem_filtro[j-1]   #Actualizacion de la memoria FIR
    mem_filtro[0]=xk[i]

    yk[i]=y


    #Algoritmo adaptive MSE



print (noise)
