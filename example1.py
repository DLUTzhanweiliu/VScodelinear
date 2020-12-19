# logP = 3
import numpy as np
m = 9
logP = int(np.ceil(np.log2(2*(m-1)**2)))
for k in range(1,int((logP+1)/2)):
    temp = (m-1)/(2**k)
    print(temp,':')
    portion = (m-1)/temp
    # print(portion)
    if temp > 1: 
        # 挑选出不做约束的线
        eliminate = [int(temp*n+1) for n in range(1,int(portion)+1) if int(temp*n+1)<m]
        print(eliminate)
        # 分成两组
        JB1 = [int(i*temp+j) for i in range(0,int(portion),2) for j  in range(1,int(temp)+2) if int(i*temp+j) not in eliminate]
        JB0 = [int(i*temp+j) for i in range(1,int(portion),2) for j  in range(1,int(temp)+2) if int(i*temp+j) not in eliminate]
    else:
        JB1 = [i for i in range(3,m,4)]
        JB0 = [i for i in range(1,m+1,4)]
    # print(m)
    # print([[i,j] for i in range(1,m+1) for j in JB1])
    # print([[j,i] for i in range(1,m+1) for j in JB1])
        # 分成两组
        # print([int(i*temp+j) for i in range(0,int(portion),2) for j  in range(1,int(temp)+2) if int(i*temp+j) not in eliminate])
        # print([int(i*temp+j) for i in range(1,int(portion),2) for j  in range(1,int(temp)+2) if int(i*temp+j) not in eliminate])
        # pass
    #     for n in range(k):
    #         print([j+n*2*int(temp) for j in range(1,int(temp)+1)])
    #         print([int(temp)+1+j+n*2*int(temp) for j in range(1,int(temp)+1)])
    # else:
    #     temp += 1
    #     for n in range(k):
    #         print([1+j+n*2*int(temp) for j in range(1,int(temp)+1)])
    #         print([1+int(temp)+1+j+n*2*int(temp) for j in range(1,int(temp)+1)])
    

    # if temp > 2:
    #     if k > 2:
    #         print(k)
    #         for n in range(0,k-1):
    #             print([j+(temp)*n*2 for j in range(1,int(temp)+1)])
    #     else:
    #         print([j for j in range(1,int(temp)+1)])
