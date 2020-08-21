import numpy as np

def readFt (FtFile):
    import nmrglue as ng
    ft_dic,ft_data = ng.pipe.read(FtFile)
    In=ft_data.transpose()

    DQ = ng.pipe.make_uc(ft_dic, ft_data, 0)
    Ys=DQ.ppm_scale()

    CS = ng.pipe.make_uc(ft_dic, ft_data, 1)
    Xs=CS.ppm_scale()
    return(In,Xs,Ys)

def shifting (In,pX,fullX,fullY,direction):
    pY=pX*2
    if direction.lower() == "pos":      # For increasing ppm axes (eg: peak at 35 ppm is now at 40 ppm)
        padX=np.zeros((pX,fullY))   # fullY= 8192 for INADEQUATE
        padY=np.zeros((fullX,pY))   # fullX=4096 for INADEQUATE
        In=np.hstack((In,padY))
        In=In[:,pY:]
        In=np.concatenate((In,padX))
        In=In[pX:,:]
    elif direction.lower() == "neg":        # For decreasing ppm axes (eg: peak at 40 ppm is now at 35 ppm)
        padX=np.zeros((pX,fullY))   # fullY= 8192 for INADEQUATE
        padY=np.zeros((fullX,pY))   # fullX=4096 for INADEQUATE
        In=np.hstack((padY,In)) 
        In=In[:,:-pY]
        In=np.concatenate((padX,In))
        In=In[:-pX,:]
    return(In)

def frange(start,end,parts):
    duration=abs(end-start)
    part_duration = duration / (parts-1)
    return [start+(i * part_duration) for i in range(parts-1,-1,-1)]

def pick (In,xppm,yppm,PPmin,PPmax,steps):
    points=np.copy(In)
    j=0
    X={}
    Y={}
    P={}
    Pts={}
    for i in frange(PPmin,PPmax,steps):
        (Xind, Yind) = np.where(abs(points) > int(i))
        points[abs(points) > i] =0    # This removes hit peaks from the previous iteration
        selX=np.around(xppm[Xind],decimals=2)
        selY=np.around(yppm[Yind],decimals=2)
        selX=xppm[Xind]
        selY=yppm[Yind]
        selP=np.vstack((selX,selY))
        P[j]=selP
        X[j]=selX
        Y[j]=selY
        Pts[j]=list(zip(selX,selY))
        j=j+1
    return (Pts,X,Y,P)