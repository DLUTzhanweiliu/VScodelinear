logP = 3
m = 3
for k in range(2,logP+1):
    if k <= (logP+1)/2:
    #     print("###############%d###################"%k)
    #     for i in range(1,m+1):
    #         pass
    #         print((i,k))
          model.addConstr(gp.quicksum(weight[i,k] for i in range(1,m+1))<=delta[k-1])
          model.addConstr(gp.quicksum(gp.quicksum(weight[i,k_up] for i in range(1,m+1) for k_up in range(k+2,logP+1)) \
          + gp.quicksum(weight[i,k_down] for i in range(1,m+1) for k_down in range(1,k-1))) <= 1- delta[k-1])
    #     print("================")
    #     for k_up in range(k+2,logP+1):
    #         for i in range(1,m+1):
    #             print((i,k_up))
    #             # pass
    #     if k>2:
    #         for k_down in range(1,k-1):
    #             for i in range(1,m+1):
    #                 print((i,k_down))
        pass
        model.addConstr(gp.quicksum(weight[k,j] for j in range(1,m+1))<=delta[int(k-1+(logP-1)/2)])
        model.addConstr(gp.quicksum(gp.quicksum(weight[k_up,j] for j in range(1,m+1) for k_up in range(k+2,logP+1)) \
          + gp.quicksum(weight[k_down,j] for j in range(1,m+1) for k_down in range(1,k-1))) <= 1- delta[int(k-1+(logP-1)/2)])
# 
        # print("###############%d###################"%k)
        # for j in range(1,m+1):
        #     pass
        #     print((k,j))
        # print("================")
        # for k_up in range(k+2,logP+1):
        #     for j in range(1,m+1):
        #         print((k_up,j))
        # if k>2:
        #     for k_down in range(1,k-1):
        #         for j in range(1,m+1):
        #             print((k_down,j))
        pass
    elif k == logP:
        model.addConstr(gp.quicksum(weight[i,j] for i in range(2,logP+1,2) for j in range(1,logP+1,2))<=delta[k])
        model.addConstr(gp.quicksum(weight[j,i] for i in range(2,logP+1,2) for j in range(1,logP+1,2))<=1-delta[k])
        # for i in range(2,logP+1,2):
        #     for j in range(1,logP+1,2):
        #         print((i,j))
        #         print((j,i))
