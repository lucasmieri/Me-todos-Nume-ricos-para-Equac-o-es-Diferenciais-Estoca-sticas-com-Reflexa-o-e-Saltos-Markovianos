import matplotlib.pyplot as plt
import numpy as np
plt.style.use('classic')

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)

r=np.random.RandomState(100)#exposes a number of methods for generating random numbers drawn from a variety of probability distributions
lamb = 2; mu = 1; Xzero = 1;
T = 1; N = 2**8; dt = 1./N;   #parametros

dW = np.sqrt(dt)*np.random.normal(0.,1.,int(N));
W = np.cumsum(dW);  # Return the cumulative sum of the elements along a given axis
#encontrei um erro na minha primeira versão do código, usei a função randn, mas por algum motivo o resultado de Winc sempre era zero, acabei por trocar para random.normal

Xtrue = Xzero*np.exp((lamb-0.5*mu**2)*(np.arange(0, T, dt)+mu*W))

plt.plot(np.arange(0,T, dt), Xtrue, '-o')


fig.savefig('graph.png')





R = 4;dt = 1/N; Dt = R*dt; L = N//R; 
Xem = np.zeros(L);                 
Xtemp = Xzero; 
for j in range(0, int(L)):
  Winc = np.sum(dW[(R*(j-1)+R):(R*j+ R)]);
  #na versão do mathlab o final da equação é apenas(R*j), entretanto isso resulta e Winc=0, como conheço pouco do mathlab não soube dizer o motivo.
  
  Xtemp = Xtemp + Dt*lamb*Xtemp + mu*Xtemp*Winc;
  Xem[j] = Xtemp;

time_EM = np.arange(0,T,Dt)
plt.plot(time_EM, Xem, '-r')
fig.savefig('graph1.png')

#reparei que o método não parece convergir sempre, a proximidade do método para o valor real difere muito dependendo do conjunto de dados.

#ref: https://jtsulliv.github.io/stock-movement/