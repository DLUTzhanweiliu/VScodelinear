import numpy as np
import gurobipy as gp
# sets
T = range(1,25)
TI = range(1,26)
N = range(1,7)
# parameters
H = 3600/1e6
# Scenario 1
inflow = [1380 for _ in T] # 637.5, 1300
load = [994.679,867.139,495.833,343.363,265.684,262.8,320.855,636.531,1065.78,\
1065.66,1065.55,1065.43,1062.55,1059.66,1059.54,923.696,801.696,820.97,685.122,\
565.893,1025.61,1055.97,1064.16,1055.74]
Y = {}
L = {}
for t in T:
    Y[t] = inflow[t-1]
    L[t] = load[t-1]
Vmin = 721
Vmax = 1123.67
Smax = 0
Hmax = 75.2
Wmin = [180,180]
Wmax = [301,290]
Gmin = [116,166]
Gmax = [182,175]
Hmin = 65
Hmax = 75.2
Acoefficient = [[2.707*1e-1, 1.215*1e-3, 1.431*1e-2, 4.112*1e-5, -8.334*1e-6, -1.728*1e-4],
              [7.769*1e-2, 3.305*1e-3, 1.180*1e-2, 5.756*1e-6, -6.962*1e-6, -9.395*1e-5]]
Dcoefficient = [1.740*1e-5, 1.615*1e-5]

Wmin = [180, 180]
Wmax = [301, 290]
Gmin = [116, 116]
Gmax = [182, 175]
E = [374.687, 1.985*1e-2]
F = [321.880, 2.030*1e-3]
Imax = 2
Ton = 6
Toff = 2
Initial_volume = 1083.7 # 1108.9, 1023.5
def getg1(w, h, A=Acoefficient[0], D=Dcoefficient[0]):
    eta = A[0]+A[1]*w+A[2]*(h-D*w**2)+A[3]*w*(h-D*w**2)+A[4]*w**2+A[5]*(h-D*w**2)**2
    g = 9.81*1e-3*eta*(h-D*w**2)*w
    return g
def getg2(w, h, A=Acoefficient[1], D=Dcoefficient[1]):
    eta = A[0]+A[1]*w+A[2]*(h-D*w**2)+A[3]*w*(h-D*w**2)+A[4]*w**2+A[5]*(h-D*w**2)**2
    g = 9.81*1e-3*eta*(h-D*w**2)*w
    return g

def getpoly(m,Hmin,Hmax,Wmin,Wmax,f):
    if m>=3 and m%2==1:
        pass
    else:
        print("m必须为大于2的奇数，例如：3,5,7,9,…")
        return 0
    num = ((m-1)/2)**2
    Hlist = np.linspace(Hmin,Hmax,m)
    Wlist = np.linspace(Wmin,Wmax,m)
    deltah = (Hmax-Hmin)/(m-1)
    deltaw = (Wmax-Wmin)/(m-1)
    J1tri = list()
    for d, w in zip(range(1,m+1),Wlist):
        for t, h in zip(range(1,m+1),Hlist):
            if d%2==0 and t%2==0:
                midpoint = (d,t)
                # 下左
                J1tri.append([(d,t),(d,t-1),(d-1,t-1)])
                # 下右
                J1tri.append([(d,t),(d,t-1),(d+1,t-1)])
                # 左下
                J1tri.append([(d,t),(d-1,t),(d-1,t-1)])
                # 左上
                J1tri.append([(d,t),(d-1,t),(d-1,t+1)])
                # 上左
                J1tri.append([(d,t),(d,t+1),(d-1,t+1)])
                # 上右
                J1tri.append([(d,t),(d,t+1),(d+1,t+1)])
                # 右上
                J1tri.append([(d,t),(d+1,t),(d+1,t+1)])
                # 右下
                J1tri.append([(d,t),(d+1,t),(d+1,t-1)])
    vertex = {}
    for i in range(m):
        for j in range(m):
            w = Wmin + deltaw*i
            h = Hmin +  deltah*j
            vertex[i+1,j+1] =  (w,h,f(w,h))
    return vertex, J1tri
m = 7
vertex, J1tri = getpoly(m,Hmin,Hmax,Wmin[0],Wmax[0],getg1)
print(len(J1tri))
