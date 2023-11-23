#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Robótica Computacional - 
# Grado en Ingeniería Informática (Cuarto)
# Práctica: Resolución de la cinemática inversa mediante CCD
#           (Cyclic Coordinate Descent).

import sys
from math import *
import numpy as np
import matplotlib.pyplot as plt
import colorsys as cs

# ******************************************************************************
# Declaración de funciones

def muestra_origenes(O,final=0):
  # Muestra los orígenes de coordenadas para cada articulación
  print('Origenes de coordenadas:')
  for i in range(len(O)):
    print('(O'+str(i)+')0\t= '+str([round(j,3) for j in O[i]]))
  if final:
    print('E.Final = '+str([round(j,3) for j in final]))




def muestra_robot(O,obj):
  # Muestra el robot graficamente
  plt.figure()
  plt.xlim(-L,L)
  plt.ylim(-L,L)
  T = [np.array(o).T.tolist() for o in O]
  for i in range(len(T)):
    plt.plot(T[i][0], T[i][1], '-o', color=cs.hsv_to_rgb(i/float(len(T)),1,1))
  plt.plot(obj[0], obj[1], '*')
  plt.pause(0.0001)
  plt.show()

#  input()
  plt.close()




def matriz_T(d,th,a,al):
  # Calcula la matriz T (ángulos de entrada en RADIANES)
  
  return [[cos(th), -sin(th)*cos(al),  sin(th)*sin(al), a*cos(th)]
         ,[sin(th),  cos(th)*cos(al), -sin(al)*cos(th), a*sin(th)]
         ,[      0,          sin(al),          cos(al),         d]
         ,[      0,                0,                0,         1]
         ]






def cin_dir(th,a):
  #Sea 'th' el vector de thetas
  #Sea 'a'  el vector de longitudes
  T = np.identity(4)
  o = [[0,0]]
  for i in range(len(th)):
    T = np.dot(T,matriz_T(0,th[i],a[i],0))
    tmp=np.dot(T,[0,0,0,1])
    o.append([tmp[0],tmp[1]])
  return o








# ******************************************************************************
# Cálculo de la cinemática inversa de forma iterativa por el método CCD

# valores articulares arbitrarios para la cinemática directa inicial
# th=[0.,0.,0.]
# a =[5.,5.,5.]
th=[0.,0.,0.]
a =[1.,1.,1.]

# ?etiqueta para diferenciar si el par th[i],a[i] es de revolución o prismático  
REV = 0
PRI = 1
# articulaciones = [REV,REV,REV]
articulaciones = [REV,PRI,REV]



L = sum(a) # variable para representación gráfica
EPSILON = .1

#plt.ion() # modo interactivo

# introducción del punto para la cinemática inversa
if len(sys.argv) != 3:
  sys.exit("python " + sys.argv[0] + " x y")
objetivo=[float(i) for i in sys.argv[1:]]
O=cin_dir(th,a)

#O=zeros(len(th)+1) # Reservamos estructura en memoria
 # Calculamos la posicion inicial
print ("- Posicion inicial:")
muestra_origenes(O)

dist = float("inf")
prev = 0.
iteracion = 1

# ? mientras la distancia sea mayor que epsilon, y además se comprueba que no se ha entrado en un bucle infinito
# ? por ejemplo, con el brazo totalmente extendido, si el objetivo está fuera del radio de acción del brazo
while (dist > EPSILON and abs(prev-dist) > EPSILON/100.):
  prev = dist
  O=[cin_dir(th,a)]
  
  # Para cada combinación de articulaciones:
  for i in range(len(th)):
    # * cálculo de la cinemática inversa:

    # T es el punto al que queremos llegar (variable objetivo)
    # oFinal es el último punto, el que alineamos entre el objetivo y 02
    # oAnterior es la articulación que movemos
    # oAnterior es un punto de O, que se va moviendo de final a principio.

    T = objetivo
    oAnterior = O[-1][len(th)-i-1]
    oFinal = O[-1][-1]


    # detectamos si la articulación es de revolución o prismática
    # * REVOLUCIÓN:
    if articulaciones[len(th)-i-1] == REV:

      opuesto1 = T[1]-oAnterior[1]
      contiguo1 = T[0]-oAnterior[0]
      opuesto2 = oFinal[1]-oAnterior[1]
      contiguo2 = oFinal[0]-oAnterior[0]

      tan1 = atan2(opuesto1,contiguo1)
      tan2 = atan2(opuesto2,contiguo2)

      Ttita = tan1-tan2
      # print ('Ttita = '+str(Ttita))

      # Actualizar th --> controlamos que th esté entre -pi y pi, para evitar que el movimiento se haga por el lado largo.
      th[len(th)-i-1] += Ttita
      if th[len(th)-i-1] > pi:
        th[len(th)-i-1] -= 2*pi
      elif th[len(th)-i-1] < -pi:
        th[len(th)-i-1] += 2*pi


    # * PRISMÁTICA:
    if articulaciones[len(th)-i-1] == PRI:
      # w es el sumatorio de todos los ángulos hasta la articulación que estamos moviendo
      w = 0
      for j in range(len(th)-i):
        w += th[j]
      print ('w = '+str(w))
      
      # vectorW = [cos(w),sin(w)], es un vector unitario en estas direcciones desde Oi-1
      vectorW = [cos(w),sin(w)]

      # d es la distancia, d= vectorW * (R -On)
      d = np.dot(vectorW,np.subtract(T,oAnterior))
      print ('d = '+str(d))

      # por último, actualizamos a
      a[len(th)-i-1] += d





    #* PENDIENTE
    #* 1. Articulaciones prismáticas
    #* 2. Control de límites de movimiento






    
    O.append(cin_dir(th,a))

  dist = np.linalg.norm(np.subtract(objetivo,O[-1][-1]))
  print ("\n- Iteracion " + str(iteracion) + ':')
  muestra_origenes(O[-1])
  muestra_robot(O,objetivo)
  print ("Distancia al objetivo = " + str(round(dist,5)))
  iteracion+=1
  O[0]=O[-1]

if dist <= EPSILON:
  print ("\n" + str(iteracion) + " iteraciones para converger.")
else:
  print ("\nNo hay convergencia tras " + str(iteracion) + " iteraciones.")
print ("- Umbral de convergencia epsilon: " + str(EPSILON))
print ("- Distancia al objetivo:          " + str(round(dist,5)))
print ("- Valores finales de las articulaciones:")
for i in range(len(th)):
  print ("  theta" + str(i+1) + " = " + str(round(th[i],3)))
for i in range(len(th)):
  print ("  L" + str(i+1) + "     = " + str(round(a[i],3)))
