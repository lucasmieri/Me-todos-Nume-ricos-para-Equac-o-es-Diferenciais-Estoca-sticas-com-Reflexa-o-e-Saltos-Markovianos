#https://math.stackexchange.com/questions/1287634/implementing-ornstein-uhlenbeck-in-matlab
#https://www.pik-potsdam.de/members/franke/lecture-sose-2016/introduction-to-python.pdf
#https://scikit-learn.org/stable/modules/model_evaluation.html

import numpy as np
import matplotlib.pyplot as plt
from math import *
plt.style.use('classic')
from sklearn.metrics import mean_squared_error


'''o processo ornstein-uhlenbeck satisfaz a SDE 
        dxt=θ(μ−xt)dt+σdWt
        discretizando :
        Xn+1=Xn+θ(μ−Xn)Δt+σΔWn
              
'''

#parametros
σ=0.5;
θ=0.1;
μ=1;
dt=1e-2;
t=np.arange(0,2,dt);


x=np.zeros(len(t));
dW=np.zeros(len(t));
for i in range(1,len(t)-1):
        dW[i] = np.sqrt(dt)*np.random.normal()
        x[i+1]=x[i]+θ*(μ-x[i])*dt+σ*dW[i];

x[0]=x[i];



#murayama
y0=0;
y=np.zeros(len(t));
yx=np.zeros(len(t));
W=np.zeros(len(t));
err=np.zeros(len(t))
control=np.zeros(len(t))

for i in range(1,len(t)-1):
    #dW = exp(θ*t[i])*dW
    W[i+1] = W[i]+dW[i]*exp(θ*t[i])#np.sqrt(dt)*exp(θ*t[i])*np.random.normal(); 
    yx[i] = exp(-θ*t[i]);
    y[i] = μ*(1-yx[i])+σ*yx[i]*W[i];
    err[i]=abs(np.log(abs( x[i]-y[i])));
    
#o ponto y[len(t)-1] sempre será zero para a funçao implementada acima, para evitar isso achei razoável iguala-lo ao ponto mais próximo.

for i in range(1, len(t)-1):
  if(err[i]>np.mean(err)):
      control[i]=err[i];

y[len(t)-1] = y[len(t)-2]
err[0]=np.mean(err);
err[len(t)-1]=np.mean(err);


#print(np.mean(err))

fig = plt.figure(figsize=(10, 8))
ax1=fig.add_axes([0,0,1,1])
ax2=fig.add_axes([0.7,0.1,0.25,0.3])
ax1.plot(t,x,'-o', label='solution ornstein-uhlenbeck');
ax1.plot(t,y,'-o',label='Euler-Maruyama');
ax2.loglog(t,err);
ax2.set_title('erro absoluto (logxlog)')
ax1.legend(loc = 'upper left')
fig.savefig('graph.png')



#fig1= plt.figure(figsize=(10,10))
#ax3=fig1.add_axes([0,0,1,1])
#ax3.scatter(t,control)
#print(control)
#fig1.savefig('graph1.png')

#print(x)
#print(y)

print(mean_squared_error(x, y))
#plt.plot(t,x,'-o', label='solution ornstein-uhlenbeck');
#plt.plot(t,y,'-o',label='Euler-Maruyama');
#plt.legend(loc = 'upper left')
#fig.savefig('graph1.png')


#ax2.plt.loglog(t,err);

