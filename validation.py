import numpy as np
import gurobipy as gp
from Bipiecewise.Bipiecewise import DCC, DLOG, ACC, SOS2, LOG, MC,LOG2
func = LOG2
def main(func,h = 70,w= 190):
    Acoefficient = [[2.707*1e-1, 1.215*1e-3, 1.431*1e-2, 4.112*1e-5, -8.334*1e-6, -1.728*1e-4],
                [7.769*1e-2, 3.305*1e-3, 1.180*1e-2, 5.756*1e-6, -6.962*1e-6, -9.395*1e-5]]
    Dcoefficient = [1.740*1e-5, 1.615*1e-5]
    def getg2(w, h, A=Acoefficient[1], D=Dcoefficient[1]):
        eta = A[0]+A[1]*w+A[2]*(h-D*w**2)+A[3]*w*(h-D*w**2)+A[4]*w**2+A[5]*(h-D*w**2)**2
        g = 9.81*1e-3*eta*(h-D*w**2)*w
        return g
    m = gp.Model()
    g = m.addVar(vtype=gp.GRB.CONTINUOUS, name="UnitOutput")
    Wmin = [180, 180]
    Wmax = [301, 290]
    Gmin = [116, 116]
    Gmax = [182, 175]
    Hmin = 65
    Hmax = 75.2
    piece = 9
    func(m,1,w,h,g,piece,Wmin[1],Wmax[1],Hmin,Hmax,getg2)
    m.Params.LogToConsole = False # 不显示求解过程
    m.optimize()
    # print(m.getVarByName('piecewise[1]').x)
    for i in m.getVars():
        print(i.varname,i.x)
    return g.x
# for func in [DCC,DLOG,ACC,LOG2,SOS2]:
#     print('{0}: {1}'.format(func.__name__,main(func,72,280)))
print(main(func,72,280))
