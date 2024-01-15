# Korteweg-de Vries (KdV) Equation
## Reference
A. L. Garcia, "Numerical Methods for Physics": \ 
Chapter 9 Partial Differential Equations IV: Stability and Implicit Methods

## Exercisies
Korteweg-de Vries equation, 
```math
\frac{\partial \rho}{\partial t}=-6\rho \frac{\partial \rho}{\partial x} - \frac{\partial ^3 \rho}{\partial x^3}
```
,is an important equation from the theory of silitons. In this repository, there are three ways to solve this PDE. I use Dirichlet boundary conditions, $\rho(x= L/2)=\rho(x= -L/2)=0$ and test the program for the solitary wave solution of the KdV equation $\rho(x,t)=2 \mathrm{sech} ^2(x-4t)$.

## Explicit scheme with iteration method
The difference equation of KdV equation using explicit scheme is
```math
\frac{\rho_j^{n+1}-\rho_j^{n}}{\tau}=-6D_j \rho_j^{n} - \frac{\rho_{j+2}^{n} - 2\rho_{j+1}^{n}+2\rho_{j-1}^{n}-\rho_{j-2}^{n}}{2h^3}
```
, where $D_j=\frac{\rho_{j+1}^{n}-\rho_{j-1}^{n}}{2h}$. Therefore, the iteration method use the equation
```math
\rho_j^{n+1}=\rho_j^{n}-6\tau D_j \rho_j^{n} - \tau \frac{\rho_{j+2}^{n} - 2\rho_{j+1}^{n}+2\rho_{j-1}^{n}-\rho_{j-2}^{n}}{2h^3}
```
to iterate the solution of $\rho$.

## Explicit scheme with matrix method
Arrange the difference equation and turn it into matrix form
```math
\Phi^{n+1} = M \Phi^n
```
, where $\Phi^{n} = \left[ \rho_1^n \ \rho_2^n \ ...... \ \rho_{N-1}^n \ \rho_{N}^n \right]$.

## Implicit scheme (Crank-Nicolson method)
The difference equation of KdV equation using explicit scheme is
```math
\frac{\rho_j^{n+1}-\rho_j^{n}}{\tau}=-6D_j \rho_j^{n} - \frac{1}{2} \left( \frac{\rho_{j+2}^{n} - 2\rho_{j+1}^{n}+2\rho_{j-1}^{n}-\rho_{j-2}^{n}}{2h^3} + \frac{\rho_{j+2}^{n+1} - 2\rho_{j+1}^{n+1}+2\rho_{j-1}^{n+1}-\rho_{j-2}^{n+1}}{2h^3}\right)
```
Arrange the difference equation
```math
\rho_j^{n+1} + \frac{\tau}{2}\frac{\rho_{j+2}^{n+1} - 2\rho_{j+1}^{n+1}+2\rho_{j-1}^{n+1}-\rho_{j-2}^{n+1}}{2h^3} = \rho_j^{n}-6\tau D_j \rho_j^{n} - \frac{\tau}{2} \frac{\rho_{j+2}^{n} - 2\rho_{j+1}^{n}+2\rho_{j-1}^{n}-\rho_{j-2}^{n}}{2h^3}
```

Turn it into matrix form
```math
W \Phi^{n+1} = M \Phi^n
```
and solve
```math
\Phi^{n+1} = W^{-1}M \Phi^n
```

## Result
![Image](https://github.com/ChenYingShan1114/Korteweg-de-Vries-Equation/blob/main/KdV.png)

### Reference
The header file of **Matrix.h** and **inv.h** are from "Numerical Methods for Physics Second Edition -- Alejandro L. Garcia".
