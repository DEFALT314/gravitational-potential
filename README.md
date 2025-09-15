# Finite Element Method Solver for Gravitational Potential

## Project Description

This project is developed **for solving the gravitational potential problem**:

$$
\frac{d^2 \Phi}{dx^2} = 4 \pi G \rho(x), \quad \Phi(0) = 5, \quad \Phi(3) = 2
$$

with the density function:

$$
\rho(x) =
\begin{cases}
-10, & x \in [0,1] \\
1, & x \in (1,2] \\
-10, & x \in (2,3]
\end{cases}
$$

where $G = 1$.

## Features

* Derivation of the **variational form** for the problem.
* Generation of a **linear system of equations** for an arbitrary number of elements `n`.
* Solution of the linear system using numerical methods.
* Visualization of the approximate gravitational potential $\Phi(x)$.
* Adjustable number of elements (`n`) as a runtime parameter.

## Requirements

* Python 3.10+
* Libraries:

  * `numpy`
  * `matplotlib`

Install dependencies with:

```bash
pip install numpy matplotlib
```

## Usage

1. Clone the repository:

```bash
git clone https://github.com/YOUR_USERNAME/repository_name.git
cd repository_name
```

2. Run the program with a specified number of elements `n`:

```bash
python main.py --n 10
```

The program will generate a plot of the **approximate gravitational potential** $\Phi(x)$.

## Methodology

1. Derive the **variational formulation** of the problem.
2. Divide the domain into `n` finite elements.
3. Compute **local element matrices and vectors**.
4. Assemble the **global linear system**.
5. Apply **boundary conditions** and solve the system.
6. Visualize the results using `matplotlib`.

## Example Plot

<img width="640" height="480" alt="image" src="https://github.com/user-attachments/assets/d086177e-de62-4ecc-89d5-50e081094cd5" />
