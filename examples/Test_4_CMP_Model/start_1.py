import torch
from torchdiffeq import odeint
import matplotlib.pyplot as plt

def derivative(t, f):
    return torch.cat([f[1:], -torch.cos(t.unsqueeze(0))])

y0 = torch.nn.Parameter(torch.Tensor([1, 0]))
yT = odeint(func=derivative,
            y0=y0,
            t=torch.Tensor([0., 10.]))

print(yT)

# 输出梯度
yT[-1, 0].backward()
print(y0.grad)