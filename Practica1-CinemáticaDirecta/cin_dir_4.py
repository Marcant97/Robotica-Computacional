#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Robótica Computacional
# Grado en Ingeniería Informática (Cuarto)
# Práctica: Resolución de la cinemática directa mediante Denavit-Hartenberg.

#? Ejemplo MANIPULADOR 2:
# 

import sys
from math import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# ******************************************************************************
# Declaración de funciones

def ramal(I,prev=[],base=0):
  # Convierte el robot a una secuencia de puntos para representar
  O = []
  if I:
    if isinstance(I[0][0],list):
      for j in range(len(I[0])):
        O.extend(ramal(I[0][j], prev, base or j < len(I[0])-1))
    else:
      O = [I[0]]
      O.extend(ramal(I[1:],I[0],base))
      if base:
        O.append(prev)
  return O

def muestra_robot(O,ef=[]):
  # Pinta en 3D
  OR = ramal(O)
  OT = np.array(OR).T
  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')
  # Bounding box cúbico para simular el ratio de aspecto correcto
  max_range = np.array([OT[0].max()-OT[0].min()
                       ,OT[1].max()-OT[1].min()
                       ,OT[2].max()-OT[2].min()
                       ]).max()
  Xb = (0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten()
     + 0.5*(OT[0].max()+OT[0].min()))
  Yb = (0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten()
     + 0.5*(OT[1].max()+OT[1].min()))
  Zb = (0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten()
     + 0.5*(OT[2].max()+OT[2].min()))
  for xb, yb, zb in zip(Xb, Yb, Zb):
     ax.plot([xb], [yb], [zb], 'w')
  ax.plot3D(OT[0],OT[1],OT[2],marker='s')
  ax.plot3D([0],[0],[0],marker='o',color='k',ms=10)
  if not ef:
    ef = OR[-1]
  ax.plot3D([ef[0]],[ef[1]],[ef[2]],marker='s',color='r')
  ax.set_xlabel('X')
  ax.set_ylabel('Y')
  ax.set_zlabel('Z')
  plt.show()
  return

def arbol_origenes(O,base=0,sufijo=''):
  # Da formato a los origenes de coordenadas para mostrarlos por pantalla
  if isinstance(O[0],list):
    for i in range(len(O)):
      if isinstance(O[i][0],list):
        for j in range(len(O[i])):
          arbol_origenes(O[i][j],i+base,sufijo+str(j+1))
      else:
        print('(O'+str(i+base)+sufijo+')0\t= '+str([round(j,3) for j in O[i]]))
  else:
    print('(O'+str(base)+sufijo+')0\t= '+str([round(j,3) for j in O]))

def muestra_origenes(O,final=0):
  # Muestra los orígenes de coordenadas para cada articulación
  print('Orígenes de coordenadas:')
  arbol_origenes(O)
  if final:
    print('E.Final = '+str([round(j,3) for j in final]))

def matriz_T(d,theta,a,alpha):
  # Calcula la matriz T (ángulos de entrada en grados)
  th=theta*pi/180;
  al=alpha*pi/180;
  return [[cos(th), -sin(th)*cos(al),  sin(th)*sin(al), a*cos(th)]
         ,[sin(th),  cos(th)*cos(al), -sin(al)*cos(th), a*sin(th)]
         ,[      0,          sin(al),          cos(al),         d]
         ,[      0,                0,                0,         1]
         ]
# ******************************************************************************

# ?
# ? Modificamos nvar a 6 variables
# ? Modificamos Parámetros D-H, añadimos una columna más y nuestros valores.


# Introducción de los valores de las articulaciones
nvar=6 # Número de variables
if len(sys.argv) != nvar+1:
  sys.exit('El número de articulaciones no es el correcto ('+str(nvar)+')')
p=[float(i) for i in sys.argv[1:nvar+1]]


# Variables
l1 = p[0]
t2 = p[1]
t3 = p[2]
t4 = p[3]
l5 = p[4]
l6 = p[5]





# Parámetros D-H:
#        1    2    3
d  = [   l1,    0,    0,     0,   5,    0,    0,   1,  1, 1+l6, l6, l6 ]
th = [    0,    0,   t2, 90+t3,  t4,    0,    0,   0,  0,    0,  0,  0 ]
a  = [    0,    2,    2,     0,   0,   l5,  -l5,   0,  0,    0,  0,  0 ]
al = [    0,    0,   90,    90,   0,    0,    0,   0,  0,    0,  0,  0]

# Orígenes para cada articulación
o00=[0,0,0,1]
o1p=[0,0,0,1] # --> punto
o11=[0,0,0,1]
o22=[0,0,0,1]
o33=[0,0,0,1]
o44=[0,0,0,1]
o4141=[0,0,0,1] # --> punto
o4242=[0,0,0,1] # --> punto
o5151=[0,0,0,1]
o5252=[0,0,0,1]
o6161=[0,0,0,1] # --> punto
o6262=[0,0,0,1] # --> punto
oEFEF=[0,0,0,1]




# Cálculo matrices transformación
T01p=matriz_T(d[0],th[0],a[0],al[0])
T1p1=matriz_T(d[1],th[1],a[1],al[1])
T12=matriz_T(d[2],th[2],a[2],al[2])
T23=matriz_T(d[3],th[3],a[3],al[3])
T34=matriz_T(d[4],th[4],a[4],al[4])
T4_41=matriz_T(d[5],th[5],a[5],al[5])
T4_42=matriz_T(d[6],th[6],a[6],al[6])
T41_51=matriz_T(d[7],th[7],a[7],al[7])
T42_52=matriz_T(d[8],th[8],a[8],al[8])
T4_ef=matriz_T(d[9],th[9],a[9],al[9])
T51_61=matriz_T(d[10],th[10],a[10],al[10])
T52_62=matriz_T(d[11],th[11],a[11],al[11])


# Nuevas multiplicaciones
T01=np.dot(T01p,T1p1)
T02=np.dot(T01,T12)
T03=np.dot(T02,T23)
T04=np.dot(T03,T34)

T0_41=np.dot(T04,T4_41)
T0_42=np.dot(T04,T4_42)

T0_51=np.dot(T0_41,T41_51)
T0_52=np.dot(T0_42,T42_52)

T0_ef=np.dot(T04,T4_ef)

T0_61=np.dot(T0_51,T51_61)
T0_62=np.dot(T0_52,T52_62)


# Transformación de cada punto
o1p0 =np.dot(T01p, o1p).tolist() # --> punto
o10 =np.dot(T01, o11).tolist()
o20 =np.dot(T02, o22).tolist()
o30 =np.dot(T03, o33).tolist()
o40 =np.dot(T04, o44).tolist()
o410 =np.dot(T0_41, o4141).tolist() # --> punto
o420 =np.dot(T0_42, o4242).tolist() # --> punto
o510=np.dot(T0_51, o5151).tolist()
o520=np.dot(T0_52, o5252).tolist()
o610=np.dot(T0_61, o6161).tolist() # --> punto
o620=np.dot(T0_62, o6262).tolist() # --> punto
oEF0=np.dot(T0_ef, oEFEF).tolist()

# Mostrar resultado de la cinemática directa
# muestra_origenes([o00,o10,o20,o30])

#? Para esta práctica
# muestra_origenes([o00,o10,o20,o30,040, [[o510],[o520]]], oEF0)
#* dejamos oEF0 sólo para que no los una
#* dejamos [[o510],[o520]] entre corchetes para que una ambos con el anterior y no entre ellos, ya que es una pinza.

muestra_origenes([o00, o1p0, o10, o20, o30, o40, [[o410, o510, o610], [o420, o520, o620]]], oEF0)
muestra_robot   ([o00, o1p0, o10, o20, o30, o40, [[o410, o510, o610], [o420, o520, o620]]], oEF0)

# muestra_robot   ([o00,o10,o20,o30])
input()

