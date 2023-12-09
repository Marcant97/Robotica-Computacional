#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Rob�tica Computacional 
# Grado en Ingenier�a Inform�tica (Cuarto)
# Pr�ctica 5:
#     Simulaci�n de robots m�viles holon�micos y no holon�micos.

#localizacion.py

import sys
from math import *
from robot import robot
import random
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

# from robot import measurement_prob


# ******************************************************************************
# Declaraci�n de funciones

def distancia(a,b):
  # Distancia entre dos puntos (admite poses)
  return np.linalg.norm(np.subtract(a[:2],b[:2]))

def angulo_rel(pose,p):
  # Diferencia angular entre una pose y un punto objetivo 'p'
  w = atan2(p[1]-pose[1],p[0]-pose[0])-pose[2]
  while w >  pi: w -= 2*pi
  while w < -pi: w += 2*pi
  return w

def mostrar(objetivos,ideal,trayectoria):
  # Mostrar objetivos y trayectoria:
  #plt.ion() # modo interactivo
  # Fijar los bordes del gr�fico
  objT   = np.array(objetivos).T.tolist()
  trayT  = np.array(trayectoria).T.tolist()
  ideT   = np.array(ideal).T.tolist()
  bordes = [min(trayT[0]+objT[0]+ideT[0]),max(trayT[0]+objT[0]+ideT[0]),
            min(trayT[1]+objT[1]+ideT[1]),max(trayT[1]+objT[1]+ideT[1])]
  centro = [(bordes[0]+bordes[1])/2.,(bordes[2]+bordes[3])/2.]
  radio  = max(bordes[1]-bordes[0],bordes[3]-bordes[2])*.75
  plt.xlim(centro[0]-radio,centro[0]+radio)
  plt.ylim(centro[1]-radio,centro[1]+radio)
  # Representar objetivos y trayectoria
  idealT = np.array(ideal).T.tolist()
  plt.plot(idealT[0],idealT[1],'-g')
  plt.plot(trayectoria[0][0],trayectoria[0][1],'or')
  r = radio * .1
  for p in trayectoria:
    plt.plot([p[0],p[0]+r*cos(p[2])],[p[1],p[1]+r*sin(p[2])],'-r')
    #plt.plot(p[0],p[1],'or')
  objT   = np.array(objetivos).T.tolist()
  plt.plot(objT[0],objT[1],'-.o')
  plt.show()
  input()
  plt.clf()

def localizacion(balizas, real, ideal, centro, radio, mostrar=0):
  # Buscar la localizaci�n m�s probable del robot, a partir de su sistema
  # sensorial, dentro de una regi�n cuadrada de centro "centro" y lado "2*radio".

  #? robot ideal encuentre al robot real. (una vez superado el umbral que definamos)
  incremento=1
  mejor_posicion = [ideal.x, ideal.y, ideal.orientation, 0] # x,y,orientacion,probabilidad
  imagen=[]
  indice = 0
  # print('centro[0]: ',centro[0])
  # print('centro[1]: ',centro[1])
  # print('radio: ',radio)
  # i = -radio
  j = centro[1]-radio
  while j <= centro[1]+radio:
    imagen.append([])
    # j = -radio
    i = centro[0]-radio
    while i <= centro[0]+radio:
      # print('i: ' + str(i) + ' j: ' + str(j))
      # cambiamos posicion de ideal
      # ideal.set(centro[0]+i,centro[1]+j,ideal.orientation)
      ideal.set(i,j,ideal.orientation)
      print('ideal: ',ideal.pose())
      
      prob = ideal.measurement_prob(ideal.sense(balizas),balizas)

      imagen[indice].append(prob)
      if prob > mejor_posicion[3]:
        # mejor_posicion = [centro[0]+i,centro[1]+j,ideal.orientation,prob]
        mejor_posicion = [i,j,ideal.orientation,prob]
        print('mejor_posicion: ',mejor_posicion)
      i += incremento
    indice += 1
    j += incremento


  # print(imagen)
  ideal.set(mejor_posicion[0],mejor_posicion[1],ideal.orientation)
  print('real: ',real.pose())
  print('ideal: ',ideal.pose())

  # guardamos los valores de la matriz imagen en un fichero de texto
  with open('imagen.txt', 'w') as f:
    # mantener la relación de filas y columnas
    for i in range(len(imagen)):
      for j in range(len(imagen[i])):
        f.write(str(imagen[i][j]) + ' ')
      f.write('\n')




  


  if mostrar:
    #plt.ion() # modo interactivo
    plt.xlim(centro[0]-radio,centro[0]+radio)
    plt.ylim(centro[1]-radio,centro[1]+radio)
    imagen.reverse()
    plt.imshow(imagen,extent=[centro[0]-radio,centro[0]+radio,\
                              centro[1]-radio,centro[1]+radio])
    balT = np.array(balizas).T.tolist();
    plt.plot(balT[0],balT[1],'or',ms=10)
    plt.plot(ideal.x,ideal.y,'D',c='#ff00ff',ms=10,mew=2)
    plt.plot(real.x, real.y, 'D',c='#00ff00',ms=10,mew=2)
    plt.show()
    input()
    plt.clf()

# ******************************************************************************

#* Definici�n del robot:
P_INICIAL = [0.,4.,0.] # Pose inicial (posici�n y orientacion)
P_PRUEBA = [0.,-2.,0.] # Pose de prueba  (posici�n y orientacion)
V_LINEAL  = .7         # Velocidad lineal    (m/s)
V_ANGULAR = 140.       # Velocidad angular   (�/s)
FPS       = 10.        # Resoluci�n temporal (fps)

HOLONOMICO = 1
GIROPARADO = 0
LONGITUD   = .2

# Definici�n de trayectorias:
trayectorias = [
    [[1,3]],
    [[0,2],[4,2]],
    [[2,4],[4,0],[0,0]],
    [[2,4],[2,0],[0,2],[4,2]],
    [[2+2*sin(.8*pi*i),2+2*cos(.8*pi*i)] for i in range(5)]
    ]

# Definici�n de los puntos objetivo:
if len(sys.argv)<2 or int(sys.argv[1])<0 or int(sys.argv[1])>=len(trayectorias):
  sys.exit(sys.argv[0]+" <indice entre 0 y "+str(len(trayectorias)-1)+">")
objetivos = trayectorias[int(sys.argv[1])]

# Definici�n de constantes:
EPSILON = .1                #* Umbral de distancia (para encontrar la baliza y pasar a la siguiente)
V = V_LINEAL/FPS            # Metros por fotograma
W = V_ANGULAR*pi/(180*FPS)  # Radianes por fotograma

ideal = robot()           #* no tiene ruido lineal ni radial.
ideal.set_noise(0,0,.1)   # Ruido lineal / radial / de sensado
ideal.set(*P_PRUEBA)     # operador 'splat'

real = robot()              #* si tiene ruido lineal y radial.
real.set_noise(.01,.01,.1)  # Ruido lineal / radial / de sensado
real.set(*P_INICIAL)

random.seed(0)
tray_ideal = [ideal.pose()]  # Trayectoria percibida
tray_real = [real.pose()]     # Trayectoria seguida

tiempo  = 0.
espacio = 0.
#random.seed(0)
random.seed(datetime.now())

#* primera llamada
#! HACER LLAMADA A FUNCIÓN DE LOCALIZCIÓN Y EL ÚLTIMO PARÁMETRO DE MOSTRAR A '1', para que muestre el mapa de colores.
localizacion(objetivos,real,ideal,[2.5,2.5],5,1)

for punto in objetivos:
  while distancia(tray_ideal[-1],punto) > EPSILON and len(tray_ideal) <= 1000:
    pose = ideal.pose()

    w = angulo_rel(pose,punto) #* velocidad angular
    if w > W:  w =  W #* las velocidades no son infinitas, si es mayor que la máxima, se pone la máxima.
    if w < -W: w = -W #* si es menor que la mínima, se pone la mínima.
    v = distancia(pose,punto) #* velocidad lineal
    if (v > V): v = V #* si es mayor que la máxima, se pone la máxima.
    if (v < 0): v = 0 #* si es menor que 0, se pone 0 --> para corregir velocidades negativas, no se puede

    if HOLONOMICO: #* holonómico
      if GIROPARADO and abs(w) > .01:
        v = 0
      ideal.move(w,v)
      real.move(w,v)
    else: #* triciclo
      ideal.move_triciclo(w,v,LONGITUD)
      real.move_triciclo(w,v,LONGITUD)
    tray_ideal.append(ideal.pose())
    tray_real.append(real.pose())

    #* segunda llamada
    #* comparar ideal con real, si el error es muy grande llamar a localizacion
    #* una vez hemos actualizdo las poses de ambos robots, llamamos a la función de localización.
    #! localizacion(objetivos,real,ideal,[2.5,2.5],5,1) --> con el mostrar a 0 para que no saque todo el rato el mapa de colores.
    # print('antes localización')
    # localizacion(objetivos,real,ideal,[2.5,2.5],5,0)
    # print('despues localización')

    espacio += v
    tiempo  += 1

if len(tray_ideal) > 1000:
  print ("<!> Trayectoria muy larga - puede que no se haya alcanzado la posicion final.")
print ("Recorrido: "+str(round(espacio,3))+"m / "+str(tiempo/FPS)+"s")
print ("Distancia real al objetivo: "+\
    str(round(distancia(tray_real[-1],objetivos[-1]),3))+"m")
mostrar(objetivos,tray_ideal,tray_real)  # Representaci�n gr�fica

