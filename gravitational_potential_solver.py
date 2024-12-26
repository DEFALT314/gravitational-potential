from matplotlib import pyplot as plt
import scipy.integrate as integrate
import numpy as np
class EquationSolver:
    def __init__(self,n=2, G=1) -> None:
        self.f_vectorized = None
        self.Y_val = None
        self.X_val = None
        self.f = None
        self.Y = None
        self.A = None
        self.n = n
        self.h = 3/n
        self.G = G
        self.create_matrix_A()
        self.create_matrix_Y()
        self.solve()
    def ro(self, x):
        if 0<=x<=1:
            return 0
        elif 1<x<=2:
            return 1
        elif 2<x<=3:
            return 0
        else: return 0
    def e(self, k,x):
        if self.h*(k-1) <= x <= self.h*(k):
            return x/self.h-k +1
        elif self.h*(k) < x <= self.h*(k+1):
            return k +1- x/self.h
        else: return 0
    def de(self, k,x):
        if self.h*(k-1) <= x <= self.h*(k):
            return 1/self.h
        elif self.h*(k) < x <= self.h*(k+1):
            return -1/self.h
        else: return 0

    def B(self, i,j,a,b):
        product_of_derivatives = lambda x: self.de(i, x) * self.de(j, x)
        return -integrate.quad(product_of_derivatives,a,b)[0]
    def L(self,i,a,b):
        def product_of_derivatives(x):
            return self.ro(x) * self.e(i, x)
        integral_1 = 4* np.pi * self.G  *integrate.quad(product_of_derivatives,a,b)[0]
        integral_2 = - integrate.quad(lambda x: self.de(i,x),a,b)[0]
        return integral_1 + integral_2
    def create_matrix_A(self):
        m = self.n - 1
        A = np.zeros((m, m))
        for i in range(m):
            index_in_matrix = i+1
            A[i,i] = self.B(index_in_matrix,index_in_matrix,max(0, (index_in_matrix-1)*self.h), min((index_in_matrix+1)*self.h,3))
        for i in range(m-1):
            index_in_matrix = i+1
            A[i, i+1] = self.B(index_in_matrix+1,index_in_matrix, max(0, (index_in_matrix)*self.h), min((index_in_matrix+1)*self.h,3))
            A[i+1, i] = A[i, i+1]
        self.A = A
    def create_matrix_Y(self):
        m = self.n - 1
        Y= np.zeros(m)
        for i in range(m):
            index_in_matrix = i+1
            Y[i] = self.L(index_in_matrix, max(0, (index_in_matrix-1)*self.h), min((index_in_matrix+1)*self.h,3))
        self.Y = Y
    def solve(self):
        C = np.linalg.solve(self.A, self.Y)
        n=len(C)
        self.X_val = np.linspace(0, 3,self.n +1)
        self.f= lambda x: (2/3)*x+5+sum(C[i]*self.e(i+1,x) for i in range(n))
        self.f_vectorized = np.vectorize(self.f)

    def visualize(self):
        ax = plt.subplot()
        ax.set(title='Gravitational Potential', xlabel='n = ' + str(self.n))
        ax.plot( self.X_val,self.f_vectorized(self.X_val), color='blue')
        plt.ylabel('Φ(x)')
        plt.grid(True)
        plt.savefig('gravitational_potential_plot.png')
        plt.show()
n =int(input("n="))
sol = EquationSolver(n)
sol.visualize()