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
                midpoint = (w,h)
                # 下左
                J1tri.append([(w,h,f(w,h)),(w,h-deltah,f(w,h-deltah)),(w-deltaw,h-deltah,f(w-deltaw,h-deltah))])
                # 下右
                J1tri.append([(w,h,f(w,h)),(w,h-deltah,f(w,h-deltah)),(w+deltaw,h-deltah,f(w+deltaw,h-deltah))])
                # 左下
                J1tri.append([(w,h,f(w,h)),(w-deltaw,h,f(w-deltaw,h)),(w-deltaw,h-deltah,f(w-deltaw,h-deltah))])
                # 左上
                J1tri.append([(w,h,f(w,h)),(w-deltaw,h,f(w-deltaw,h)),(w-deltaw,h+deltah,f(w-deltaw,h+deltah))])
                # 上左
                J1tri.append([(w,h,f(w,h)),(w,h+deltah,f(w,h+deltah)),(w-deltaw,h+deltah,f(w-deltaw,h+deltah))])
                # 上右
                J1tri.append([(w,h,f(w,h)),(w,h+deltah,f(w,h+deltah)),(w+deltaw,h+deltah,f(w+deltaw,h+deltah))])
                # 右上
                J1tri.append([(w,h,f(w,h)),(w+deltaw,h,f(w+deltaw,h)),(w+deltaw,h+deltah,f(w+deltaw,h+deltah))])
                # 右下
                J1tri.append([(w,h,f(w,h)),(w+deltaw,h,f(w+deltaw,h)),(w+deltaw,h-deltah,f(w+deltaw,h-deltah))])
    return J1tri
def DCC(model,u,w,h,g,m,Wmin,Wmax,Hmin,Hmax,getg):
    # DCC
    # 定义变量，DCC方法对每个三角形有一个二进制变量标识
    # 不单独考虑共用顶点的问题
    delta = model.addVars(2**m, vtype=gp.GRB.BINARY, name="piecewise")
    weight = model.addVars(3*2**m, vtype=gp.GRB.CONTINUOUS,lb=0,name="weight")
    # 1.J1 triangulation 给出每个顶点的坐标，以及剖分后三角形的坐标
    J1tri = getpoly(m,Hmin,Hmax,Wmin,Wmax,getg)
    # 2.添加约束，凸组合
    model.addConstr(h<=gp.quicksum(weight[3*p+v]*J1tri[p][v][1] for v in range(3) for p in range(len(J1tri)))+Hmax*(1-u))
    model.addConstr(h>=gp.quicksum(weight[3*p+v]*J1tri[p][v][1] for v in range(3) for p in range(len(J1tri))))
    model.addConstr(w==gp.quicksum(weight[3*p+v]*J1tri[p][v][0] for v in range(3) for p in range(len(J1tri))))
    model.addConstr(g==gp.quicksum(weight[3*p+v]*J1tri[p][v][2] for v in range(3) for p in range(len(J1tri))))
    model.addConstrs(delta[p]==gp.quicksum(weight[3*p+v] for v in range(3)) for p in range(len(J1tri)))
    model.addConstr(u==gp.quicksum(delta[p] for p in range(len(J1tri))))
    # DLOG
#     delta = model.addVars(m, vtype=gp.GRB.BINARY, name="piecewise")
    # 编码
#     J1tri = getpoly(m,Hmin,Hmax,Wmin,Wmax,getg)
#     J1triencode = J1tri.reshape([2]*m)
    
        #import numpy as np
    #wlist = np.linspace(Wmin, Wmax, W)
    #hlist = np.linspace(Hmin, Hmax, H)
    #glist = list()
    #Wlist, Hlist = np.meshgrid(wlist,hlist)
    # for w,h in zip(Wlist.flatten(), Hlist.flatten()):
    #def getg(Wlist, Hlist):
        #delta = A[0]+A[1]*w+A[2]*(h-D*w**2)+A[3]*w*(h-D*w**2)+A[4]*w**2+A[5]*(h-D*w**2)**2
        #g = 9.81*1e-3*delta*(h-D*w**2)*w
        # glist.append(g)
        #return g
    #glist = getg(Wlist, Hlist)
m = gp.Model()
# add variables 
v = m.addVars(TI, lb=Vmin, ub=Vmax, vtype=gp.GRB.CONTINUOUS, name="Volume")
q = m.addVars(T, vtype=gp.GRB.CONTINUOUS, name="PlantDischarge")
s = m.addVars(T, vtype=gp.GRB.CONTINUOUS, name="Spill")
g = m.addVars(N, T, vtype=gp.GRB.CONTINUOUS, name="UnitOutput")
w = m.addVars(N, T, vtype=gp.GRB.CONTINUOUS, name="UnitDischarge")
h = m.addVars(T, vtype=gp.GRB.CONTINUOUS, name="GrossHead")

u = m.addVars(N, T, vtype=gp.GRB.BINARY, name="UnitState")
i = m.addVars(N, T, vtype=gp.GRB.BINARY, name="UnitStartup")
o = m.addVars(N, T, vtype=gp.GRB.BINARY, name="UnitShutdown")
piece = 3 # 

m.update()

# add objective function
m.setObjective(gp.quicksum(q[t]+s[t] for t in T), gp.GRB.MINIMIZE)

# add constraints
## 负荷平衡
m.addConstrs(gp.quicksum(g[n, t] for n in N) == L[t] for t in T)
## 机组发电函数

for n in N:
    if n in [1, 2, 3, 4]:
        m.addConstrs(u[n, t]*Wmin[0]<=w[n, t] for t in T)
        m.addConstrs(u[n, t]*Wmax[0]>=w[n, t] for t in T)
        m.addConstrs(u[n, t]*Gmin[0]<=g[n, t] for t in T)
        m.addConstrs(u[n, t]*Gmax[0]>=g[n, t] for t in T)
        for t in T:
            DCC(m,u[n, t],w[n, t],h[t],g[n, t],3,Wmin[0],Wmax[0],Hmin,Hmax,getg1)
        
    elif n in [5, 6]:
        m.addConstrs(u[n, t]*Wmin[1]<=w[n, t] for t in T)
        m.addConstrs(u[n, t]*Wmax[1]>=w[n, t] for t in T)
        m.addConstrs(u[n, t]*Gmin[1]<=g[n, t] for t in T)
        m.addConstrs(u[n, t]*Gmax[1]>=g[n, t] for t in T)
        for t in T:
            DCC(m,u[n, t],w[n, t],h[t],g[n, t],3,Wmin[1],Wmax[1],Hmin,Hmax,getg2)
## 毛水头计算
m.addConstrs(h[t] == ((E[0] + E[1]*v[t])+(E[0] + E[1]*v[t+1]))/2-(F[0]+F[1]*(q[t]+s[t])) for t in T)
## 电站发电流量
m.addConstrs(gp.quicksum(w[n, t] for n in N) == q[t] for t in T)
## 开停机逻辑约束
m.addConstrs(u[n, t] - u[n, t-1]>=i[n, t]-o[n, t] for n in N for t in T if t > 1)
## 最大开停机次数约束
m.addConstrs(gp.quicksum(i[n, t] for t in T)<=Imax for n in N)
## 最小启停持续时段约束
for t in T:
    if t-Ton<1:
        m.addConstr(gp.quicksum(i[n, tt] for tt in range(1,t+1))<=u[n, t])
    else:
        m.addConstr(gp.quicksum(i[n, tt] for tt in range(t-Ton,t+1))<=u[n, t])
    if t-Toff<1:
        m.addConstr(gp.quicksum(i[n, tt] for tt in range(1,t+1))<=1-u[n, t])
    else:
        m.addConstr(gp.quicksum(i[n, tt] for tt in range(t-Toff,t+1))<=1-u[n, t])
# model solve
m.Params.LogToConsole = True # 显示求解过程
m.Params.MIPGap = 0.01 # 百分比界差
m.Params.TimeLimit = 100 # 限制求解时间为 100s

m.optimize()