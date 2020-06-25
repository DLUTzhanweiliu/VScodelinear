import gurobipy as gp
import numpy as np
'''
The module represents or approximates a function of two variables, f(x,y), by
piecewise linear interpolation, using a triangular tesselation of the 2D space.
A 3x3 grid has 8 triangles; (J1 triangulation)
mxm grid = 2*(m-1)**2 triangles
m = 3,5,7...
考虑共用顶点，顶点个数 m**2; 
不考虑共用顶点，顶点个数 3*2*(m-1)**2
Graphically:
    \         /
     |---|---|
     |\  |  /| 
     | \ | / |
     |  \|/  |
     |---|---|
     |  /|\  | 
     | / | \ |
     |/  |  \|
     |---|---| 
    /         \;
The function is linear over each triangle;
There are 7 methods. 
Non-Parametric
    Convex combination            
        Disaggregated   
            ## Basic (DCC)
            ## Logarithmic (DLOG)
        Aggregated
            ## Basic (ACC)
            ## Logarithmic (LOG)
            ## SOS2 (SOS2)
    ## Incremental (INC)
Parametric     
    ## Multiple choice (MC)
'''
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
def DCC(model,u,w,h,g,m,Wmin,Wmax,Hmin,Hmax,getg):
    # Disaggregated convex combination
    # 定义变量，DCC方法对每个三角形有一个二进制变量（delta(P)）标识，每个顶点一个权重
    # 三角形的个数：2*(m-1)**2
    # 顶点的个数：3*2*(m-1)**2 （不考虑共用顶点，顶点个数即是三角形个数的三倍）
    delta = model.addVars(2*(m-1)**2, vtype=gp.GRB.BINARY, name="piecewise")
    weight = model.addVars(3*2*(m-1)**2, vtype=gp.GRB.CONTINUOUS,lb=0,name="weight")
    # 1.J1 triangulation 给出每个顶点的坐标，以及剖分后三角形的坐标
    vertex, J1tri = getpoly(m,Hmin,Hmax,Wmin,Wmax,getg)
    # 2.添加约束，凸组合
    model.addConstr(h<=gp.quicksum(weight[3*p+v]*vertex[J1tri[p][v]][1] for v in range(3) for p in range(len(J1tri)))+Hmax*(1-u))
    model.addConstr(h>=gp.quicksum(weight[3*p+v]*vertex[J1tri[p][v]][1] for v in range(3) for p in range(len(J1tri))))
    model.addConstr(w==gp.quicksum(weight[3*p+v]*vertex[J1tri[p][v]][0] for v in range(3) for p in range(len(J1tri))))
    model.addConstr(g==gp.quicksum(weight[3*p+v]*vertex[J1tri[p][v]][2] for v in range(3) for p in range(len(J1tri))))
    model.addConstrs(delta[p]==gp.quicksum(weight[3*p+v] for v in range(3)) for p in range(len(J1tri)))
    model.addConstr(u==gp.quicksum(delta[p] for p in range(len(J1tri))))

def DLOG(model,u,w,h,g,m,Wmin,Wmax,Hmin,Hmax,getg):
    # DLOG
    # 定义变量，DLOG方法对每个三角形有一个二进制变量标识
    # 单独考虑共用顶点的问题
    delta = model.addVars(range(1,m+1), vtype=gp.GRB.BINARY, name="piecewise")
    weight = model.addVars(3*2**m, vtype=gp.GRB.CONTINUOUS,lb=0,name="weight")
    # 1.J1 triangulation 给出每个顶点的坐标，以及剖分后三角形的坐标
    vertex, J1tri = getpoly(m,Hmin,Hmax,Wmin,Wmax,getg)
    
    # 2.添加约束，凸组合  
    model.addConstr(h<=gp.quicksum(weight[3*p+v]*vertex[J1tri[p][v]][1] for v in range(3) for p in range(len(J1tri)))+Hmax*(1-u))
    model.addConstr(h>=gp.quicksum(weight[3*p+v]*vertex[J1tri[p][v]][1] for v in range(3) for p in range(len(J1tri))))
    model.addConstr(w==gp.quicksum(weight[3*p+v]*vertex[J1tri[p][v]][0] for v in range(3) for p in range(len(J1tri))))
    model.addConstr(g==gp.quicksum(weight[3*p+v]*vertex[J1tri[p][v]][2] for v in range(3) for p in range(len(J1tri))))
    model.addConstr(u==gp.quicksum(weight[3*p+v] for v in range(3) for p in range(len(J1tri))))
    
    # 给所有三角形使用格雷码进行标识
    greycode = []
    for i in range(0, 1<<m):
        gray=i^(i>>1)
        greycode.append("{0:0{1}b}".format(gray,m))
    J1tri = dict(zip(greycode,J1tri))
    model.addConstrs(gp.quicksum(gp.quicksum(weight[3*p+v] for v in range(3)) for k,p in zip(J1tri.keys(),range(len(J1tri))) if int(k[i])==1)<=delta[i+1] for i in range(m))
    model.addConstrs(gp.quicksum(gp.quicksum(weight[3*p+v] for v in range(3)) for k,p in zip(J1tri.keys(),range(len(J1tri))) if int(k[i])==0)<=1-delta[i+1] for i in range(m))
def ACC(model,u,w,h,g,m,Wmin,Wmax,Hmin,Hmax,getg):
    # ACC
    # 定义变量，DLOG方法对每个三角形有一个二进制变量标识
    # 单独考虑共用顶点的问题
    delta = model.addVars(range(1,2**m+1), vtype=gp.GRB.BINARY, name="piecewise")
    weight = model.addVars(range(1,m+1),range(1,m+1), vtype=gp.GRB.CONTINUOUS,lb=0,name="weight")
    # 1.J1 triangulation 给出每个顶点的坐标，以及剖分后三角形的坐标
    vertex, J1tri = getpoly(m,Hmin,Hmax,Wmin,Wmax,getg)
    # print(vertex[1,1])
    # 2.添加约束，凸组合  
    model.addConstr(h<=gp.quicksum(weight[i,j]*vertex[i,j][1] for i in range(1,m+1) for j in range(1,m+1))+Hmax*(1-u))
    model.addConstr(h>=gp.quicksum(weight[i,j]*vertex[i,j][1] for i in range(1,m+1) for j in range(1,m+1)))
    model.addConstr(w==gp.quicksum(weight[i,j]*vertex[i,j][0] for i in range(1,m+1) for j in range(1,m+1)))
    model.addConstr(g==gp.quicksum(weight[i,j]*vertex[i,j][2] for i in range(1,m+1) for j in range(1,m+1)))

    model.addConstr(u==gp.quicksum(weight[i,j] for i in range(1,m+1) for j in range(1,m+1)))
    
    model.addConstr(u==gp.quicksum(delta[k] for k in range(1,m+1)))
    # 给所有三角形使用格雷码进行标识

    model.addConstrs(weight[i,j]<=gp.quicksum(delta[k+1] for k, tri in enumerate(J1tri) if (i,j) in tri) for i in range(1,m+1) for j in range(1,m+1))
def LOG(model,u,w,h,g,m,Wmin,Wmax,Hmin,Hmax,getg):
    # ACC
    # 定义变量，DLOG方法对每个三角形有一个二进制变量标识
    # 单独考虑共用顶点的问题
    logP = int(np.ceil(np.log2(2*(m-1)**2)))
    delta = model.addVars(range(1,logP+1), vtype=gp.GRB.BINARY, name="piecewise")
    weight = model.addVars(range(1,m+1),range(1,m+1), vtype=gp.GRB.CONTINUOUS,lb=0,name="weight")
    # 1.J1 triangulation 给出每个顶点的坐标，以及剖分后三角形的坐标
    vertex, J1tri = getpoly(m,Hmin,Hmax,Wmin,Wmax,getg)
    # print(vertex[1,1])
    # 2.添加约束，凸组合  
    model.addConstr(h<=gp.quicksum(weight[i,j]*vertex[i,j][1] for i in range(1,m+1) for j in range(1,m+1))+Hmax*(1-u))
    model.addConstr(h>=gp.quicksum(weight[i,j]*vertex[i,j][1] for i in range(1,m+1) for j in range(1,m+1)))
    model.addConstr(w==gp.quicksum(weight[i,j]*vertex[i,j][0] for i in range(1,m+1) for j in range(1,m+1)))
    model.addConstr(g==gp.quicksum(weight[i,j]*vertex[i,j][2] for i in range(1,m+1) for j in range(1,m+1)))

    model.addConstr(u==gp.quicksum(weight[i,j] for i in range(1,m+1) for j in range(1,m+1)))
    for k in range(2,logP+1):
        if k <= (logP+1)/2:
            # 竖轴
            model.addConstr(gp.quicksum(weight[i,k] for i in range(1,m+1))<=delta[k-1])
            model.addConstr(gp.quicksum(weight[i,k_up] for i in range(1,m+1) for k_up in range(k+2,logP+1)) \
            + gp.quicksum(weight[i,k_down] for i in range(1,m+1) for k_down in range(1,k-1)) <= 1- delta[k-1])
            # 横轴
            model.addConstr(gp.quicksum(weight[k,j] for j in range(1,m+1))<=delta[int(k-1+(logP-1)/2)])
            model.addConstr(gp.quicksum(weight[k_up,j] for j in range(1,m+1) for k_up in range(k+2,logP+1)) \
            + gp.quicksum(weight[k_down,j] for j in range(1,m+1) for k_down in range(1,k-1)) <= 1- delta[int(k-1+(logP-1)/2)])
        elif k == logP:
            # 三角形的识别
            model.addConstr(gp.quicksum(weight[i,j] for i in range(2,logP+1,2) for j in range(1,logP+1,2))<=delta[k])
            model.addConstr(gp.quicksum(weight[j,i] for i in range(2,logP+1,2) for j in range(1,logP+1,2))<=1-delta[k])
def SOS2(model,u,w,h,g,m,Wmin,Wmax,Hmin,Hmax,getg):
    # ACC
    # 定义变量，DLOG方法对每个三角形有一个二进制变量标识
    # 单独考虑共用顶点的问题
    logP = int(np.ceil(np.log2(2*(m-1)**2)))
    # delta = model.addVars(range(1,logP+1), vtype=gp.GRB.BINARY, name="piecewise")
    weight = model.addVars(range(1,m+1),range(1,m+1), vtype=gp.GRB.CONTINUOUS,lb=0,name="weight")
    epsilon = model.addVars(range(1,m+1), vtype=gp.GRB.CONTINUOUS,lb=0, name="SOS2-x")
    gamma = model.addVars(range(1,m+1), vtype=gp.GRB.CONTINUOUS,lb=0, name="SOS2-y")
    delta = model.addVar(vtype=gp.GRB.BINARY, name="bin")
    # 1.J1 triangulation 给出每个顶点的坐标，以及剖分后三角形的坐标
    vertex, J1tri = getpoly(m,Hmin,Hmax,Wmin,Wmax,getg)
    # print(vertex[1,1])
    # 2.添加约束，凸组合  
    model.addConstr(h<=gp.quicksum(weight[i,j]*vertex[i,j][1] for i in range(1,m+1) for j in range(1,m+1))+Hmax*(1-u))
    model.addConstr(h>=gp.quicksum(weight[i,j]*vertex[i,j][1] for i in range(1,m+1) for j in range(1,m+1)))
    model.addConstr(w==gp.quicksum(weight[i,j]*vertex[i,j][0] for i in range(1,m+1) for j in range(1,m+1)))
    model.addConstr(g==gp.quicksum(weight[i,j]*vertex[i,j][2] for i in range(1,m+1) for j in range(1,m+1)))

    model.addConstrs(epsilon[i] == gp.quicksum(weight[i,j] for j in range(1,1+m)) for i in range(1,m+1))
    model.addConstrs(gamma[i] == gp.quicksum(weight[j,i] for j in range(1,1+m)) for i in range(1,m+1))
    
    model.addSOS(gp.GRB.SOS_TYPE2,[epsilon[i] for i in range(1,m+1)])
    model.addSOS(gp.GRB.SOS_TYPE2,[gamma[i] for i in range(1,m+1)])

    # 三角形的识别
    model.addConstr(gp.quicksum(weight[i,j] for i in range(2,logP+1,2) for j in range(1,logP+1,2))<=delta)
    model.addConstr(gp.quicksum(weight[j,i] for i in range(2,logP+1,2) for j in range(1,logP+1,2))<=1-delta)

def MC(model,u,w,h,g,m,Wmin,Wmax,Hmin,Hmax,getg):
    # MC
    # 定义变量，DLOG方法对每个三角形有一个二进制变量标识
    # 生成三角形的数量
    num = int(2*(m-1)**2)

    delta = model.addVars(range(1,num+1), vtype=gp.GRB.BINARY, name="piecewise")
    weight = model.addVars(range(1,m+1),range(1,m+1), vtype=gp.GRB.CONTINUOUS,lb=0,name="weight")
    epsilon = model.addVars(range(1,m+1), vtype=gp.GRB.CONTINUOUS,lb=0, name="SOS2-x")
    gamma = model.addVars(range(1,m+1), vtype=gp.GRB.CONTINUOUS,lb=0, name="SOS2-y")
    delta = model.addVar(vtype=gp.GRB.BINARY, name="bin")
    # 1.J1 triangulation 给出每个顶点的坐标，以及剖分后三角形的坐标
    vertex, J1tri = getpoly(m,Hmin,Hmax,Wmin,Wmax,getg)
    # print(vertex[1,1])
    # 2.添加约束，凸组合  
    model.addConstr(h<=gp.quicksum(weight[i,j]*vertex[i,j][1] for i in range(1,m+1) for j in range(1,m+1))+Hmax*(1-u))
    model.addConstr(h>=gp.quicksum(weight[i,j]*vertex[i,j][1] for i in range(1,m+1) for j in range(1,m+1)))
    model.addConstr(w==gp.quicksum(weight[i,j]*vertex[i,j][0] for i in range(1,m+1) for j in range(1,m+1)))
    model.addConstr(g==gp.quicksum(weight[i,j]*vertex[i,j][2] for i in range(1,m+1) for j in range(1,m+1)))

    model.addConstrs(epsilon[i] == gp.quicksum(weight[i,j] for j in range(1,1+m)) for i in range(1,m+1))
    model.addConstrs(gamma[i] == gp.quicksum(weight[j,i] for j in range(1,1+m)) for i in range(1,m+1))
    
    model.addSOS(gp.GRB.SOS_TYPE2,[epsilon[i] for i in range(1,m+1)])
    model.addSOS(gp.GRB.SOS_TYPE2,[gamma[i] for i in range(1,m+1)])

    # 三角形的识别
    model.addConstr(gp.quicksum(weight[i,j] for i in range(2,logP+1,2) for j in range(1,logP+1,2))<=delta)
    model.addConstr(gp.quicksum(weight[j,i] for i in range(2,logP+1,2) for j in range(1,logP+1,2))<=1-delta)
