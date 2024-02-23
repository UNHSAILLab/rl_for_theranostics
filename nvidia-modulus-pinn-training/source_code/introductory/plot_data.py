import matplotlib.pyplot as plt
import numpy as np

data = np.load('outputs/pinn/inferencers/inf_data.npz', allow_pickle=True)
data = np.atleast_1d(data.f.arr_0)[0]

plt.figure()
x = data['x'].flatten()
pred_u = data['u'].flatten()
plt.plot(np.sort(x), pred_u[np.argsort(x)], label='Neural Solver')
plt.plot(np.sort(x), 0.5*(np.sort(x)*(np.sort(x)-1)), label='(1/2)(x-1)x')
x_np = np.array([0., 1.])
u_np = 0.5*(x_np-1)*x_np
plt.scatter(x_np, u_np, label='BC')
plt.legend()
plt.savefig("comparison_pinn.png")

data = np.load('outputs/data/inferencers/inf_data.npz', allow_pickle=True)
data = np.atleast_1d(data.f.arr_0)[0]

plt.figure()
x = data['x'].flatten()
pred_u = data['u'].flatten()
plt.plot(np.sort(x), pred_u[np.argsort(x)], label='Neural Solver')
plt.plot(np.sort(x), 0.5*(np.sort(x)*(np.sort(x)-1)), label='(1/2)(x-1)x')
x_np = np.linspace(0, 0.3, 4)
u_np = 0.5*(x_np-1)*x_np
plt.scatter(x_np, u_np, label='Data')
plt.legend()
plt.savefig("comparison_data.png")

data = np.load('outputs/data_plus_pinn/inferencers/inf_data.npz', allow_pickle=True)
data = np.atleast_1d(data.f.arr_0)[0]

plt.figure()
x = data['x'].flatten()
pred_u = data['u'].flatten()
plt.plot(np.sort(x), pred_u[np.argsort(x)], label='Neural Solver')
plt.plot(np.sort(x), 0.5*(np.sort(x)*(np.sort(x)-1)), label='(1/2)(x-1)x')
x_np = np.linspace(0, 0.3, 4)
u_np = 0.5*(x_np-1)*x_np
plt.scatter(x_np, u_np, label='Data')
plt.legend()
plt.savefig("comparison_data_plus_pinn.png")

plt.figure()
for i in range(11):
    data = np.load('outputs/pinn_parameterized/inferencers/inf_data_'+str(i)+'.npz', allow_pickle=True)
    data = np.atleast_1d(data.f.arr_0)[0]
    
    x = data['x'].flatten()
    pred_u = data['u'].flatten()
    plt.plot(np.sort(x), pred_u[np.argsort(x)], label=str((data['l'].flatten())[0]))

plt.legend()
plt.savefig("comparison_pinn_parameterized.png")

