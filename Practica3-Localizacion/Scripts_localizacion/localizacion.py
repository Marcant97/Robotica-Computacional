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
from matplotlib.patches import FancyArrow

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
  arrows_ideal = []  # Lista para almacenar las flechas del robot ideal
  arrows_real = []  # Lista para almacenar las flechas del robot real

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
  plt.plot(trayectoria[0][0],trayectoria[0][1],'or')
  r = radio * .1
  objT   = np.array(objetivos).T.tolist()
  plt.plot(objT[0],objT[1],'-.o')


  plt.ion() #* modo interactivo activado

  #* Se imprime trayectoria del robot real poco a poco (con flechas cada 10 puntos)
  for i in range(0, len(trayectoria), 5):
    p = trayectoria[i]
    arrow_length = r
    dx = arrow_length * cos(p[2])
    dy = arrow_length * sin(p[2])
    arrow_real = FancyArrow(p[0], p[1], dx, dy, color='red', width=0.01, head_width=0.2, head_length=0.1)
    plt.gca().add_patch(arrow_real)
    arrows_real.append(arrow_real)
    plt.pause(0.1)  # Pausa para permitir la interactividad

  #* Se imprime trayectoria del robot ideal poco a poco (con flechas cada 10 puntos)
  for i in range(0, len(ideal), 5):
    p = ideal[i]
    arrow_length = r
    dx = arrow_length * cos(p[2])
    dy = arrow_length * sin(p[2])
    arrow_ideal = FancyArrow(p[0], p[1], dx, dy, color='green', width=0.01, head_width=0.2, head_length=0.1)
    plt.gca().add_patch(arrow_ideal)
    arrows_ideal.append(arrow_ideal)
    plt.pause(0.1)  # Pausa para permitir la interactividad


#* Se imprime trayectoria del robot real e ideal
  for p in trayectoria: # real
    plt.plot([p[0],p[0]+r*cos(p[2])],[p[1],p[1]+r*sin(p[2])],'-r')
      #   #plt.plot(p[0],p[1],'or')

  plt.plot(idealT[0],idealT[1],'-g') #ideal


  plt.ioff()  # Desactivar el modo interactivo al final

  #* Eliminamos las flechas para que quede más limpio el dibujo
  # Eliminar las flechas del robot real al finalizar el bucle
  for arrow in arrows_real:
      arrow.remove()
  # Eliminar las flechas del robot ideal al finalizar el bucle
  for arrow in arrows_ideal:
      arrow.remove()


  plt.show()
  input()
  plt.clf()


def localizacion(balizas, real, ideal, centro, radio, mostrar=0):
  # Buscar la localizaci�n m�s probable del robot, a partir de su sistema
  # sensorial, dentro de una regi�n cuadrada de centro "centro" y lado "2*radio".

  mejor_error = ideal.measurement_prob(real.sense(balizas), balizas)  # balizas, puntos objetivos
  mejor_pos = ideal.pose()
  orientacion_ideal = ideal.orientation
  incremento=0.1
  imagen=[]

  rango = np.arange(-radio, radio, incremento) # obtenemos rangos a partir del radio y el incremento

  for j in rango:
    imagen.append([])

    for i in rango:
      # Actualizamos la posicion del robot ideal y calculamos el error
      ideal.set(centro[0]+i, centro[1]+j, orientacion_ideal)
      error_actual = ideal.measurement_prob(real.sense(balizas),balizas)

      # Introducimos el error en la fila actual, ya que -1 es la última.
      imagen[-1].append(error_actual)
      # Si el error es mejor, actualizamos y guardamos la posición.
      if error_actual < mejor_error:
        mejor_error = error_actual
        mejor_pos = [centro[0]+i, centro[1]+j, orientacion_ideal]

  # Actualizamos el robot ideal con la mejor posición
  ideal.set(mejor_pos[0],mejor_pos[1], orientacion_ideal) 


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

  # *****
    #* Flecha que indica la orientación
    arrow_length = 0.8  # Longitud de la flecha
    dx = arrow_length * np.cos(ideal.orientation)
    dy = arrow_length * np.sin(ideal.orientation)
    plt.arrow(ideal.x, ideal.y, dx, dy, color='violet', width=0.05, head_width=0.2, head_length=0.2)
    plt.arrow(real.x, real.y, dx, dy, color='green', width=0.05, head_width=0.2, head_length=0.2)

  # *****
    plt.show()
    input()
    plt.clf()

# ******************************************************************************

#* Definici�n del robot:
P_INICIAL = [0.,4.,0.] # Pose inicial (posici�n y orientacion)
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
ideal.set(*P_INICIAL)     # operador 'splat'

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



#* primera llamada, antes del bucle.
localizacion(objetivos,real,ideal,[ideal.x,ideal.y],5,1)

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

    # cálculo del error actual
    error_actual = ideal.measurement_prob(real.sense(objetivos),objetivos)
    # menos de 0.15 es demasiado pequeño y tarda demasiado.
    # más de 0.15 es demasiado grande y suele hacer pocas correcciones
    if error_actual > 0.15: 
      print('----------------')
      print('error detectado actual: ', error_actual)
      print('robot real: ', real.pose())
      print('robot ideal previo: ', ideal.pose())
      localizacion(objetivos,real,ideal,[ideal.x,ideal.y],2,0)
      print('robot ideal corregido: ', ideal.pose())

    espacio += v
    tiempo  += 1

if len(tray_ideal) > 1000:
  print ("<!> Trayectoria muy larga - puede que no se haya alcanzado la posicion final.")
print ("Recorrido: "+str(round(espacio,3))+"m / "+str(tiempo/FPS)+"s")
print ("Distancia real al objetivo: "+\
    str(round(distancia(tray_real[-1],objetivos[-1]),3))+"m")
mostrar(objetivos,tray_ideal,tray_real)  # Representaci�n gr�fica