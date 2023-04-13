import numpy as np
from matplotlib.pyplot import figure,axes,plot,title,xlabel,ylabel,show
from matplotlib.animation import FuncAnimation

n2=2
x0=0.0
x1=1.0
NX=368
x=np.arange(NX+1)*(x1-x0)/NX            # Table with x-values
n1=1                                    # Ground state quantum number
psi1=np.sqrt(2/x1)*np.sin(n1*np.pi*x/x1)# Ground state eigenfunction
psi2=np.sqrt(2/x1)*np.sin(n2*np.pi*x/x1)# Excited state eigenfunction

t0=0
T=2*np.pi
NT=1200
t=np.arange(NT+1)*(T-t0)/NT             # Time table 
om1=1
om2=om1*n2**2

Psisq=np.empty((NT+1,NX+1))
expectationvalue=np.zeros(NT+1)

for nt in range(NT+1):
    exp1=np.exp(-1j*om1*t[nt])
    exp2=np.exp(-1j*om2*t[nt])
    Psisq[nt,:]=abs(psi1*exp1+psi2*exp2)**2/2
    s=0
    for i in range(1,NX):
        s+=x[i]*Psisq[nt,i]             # Integration by trapezoidal rule
    expectationvalue[nt]=s*(x1-x0)/NX

y_data = Psisq
x_data = x
xpos = expectationvalue

fig=figure()
#ax=axes(xlim=x_lim,ylim=y_lim)
ax=axes(xlim=(x0,x1),ylim=(0,np.max(Psisq)))
pos, =plot([],[],'ro')
line, = ax.plot([], [], lw=2)
def init():
    line.set_data([],[])
    pos.set_data([],[])
    return line,pos
def anim(i,x_data,y_data):
    line.set_data(x_data,y_data[i,:])
    pos.set_data(xpos[i],0.1)
    return line,pos
ani=FuncAnimation(fig,anim,init_func=init,fargs=(x_data,y_data),frames=NT,interval=10,blit=False,repeat=False)
show()