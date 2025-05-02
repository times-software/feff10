import numpy as np
fl1 = input()
fl2 = input()
c1 = int(input())-1
c2 = int(input())-1
x1,y1=np.loadtxt(fl1,usecols=(c1,c2),unpack=True)
x2,y2=np.loadtxt(fl2,usecols=(c1,c2),unpack=True)
xmin = max(x1[0],x2[0])
xmax = min(x1[-1],x2[-1])
x = np.linspace(xmin,xmax,100)
y1_interp=np.interp(x,x1,y1)
y2_interp=np.interp(x,x2,y2)
rsqr = (y2_interp - y1_interp)**2/(y1_interp + y2_interp)
print(np.average(rsqr))
