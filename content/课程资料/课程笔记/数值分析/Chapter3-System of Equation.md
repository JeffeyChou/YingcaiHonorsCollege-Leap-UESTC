---
date: 2023-03-20 08:42:09
categories: Numerical_Analysis
excerpt: 数值分析第三章-求解线性方程组
katex: true
rating: ⭐⭐
draft: false
tags:
  - Numerical_Analysis
title: "Chapter3-System of Equation"
share: false
updated: 2023-06-13 20:48:18
---

In this chapter, we consider the problem of solving several simultaneous Equations in several variables. Most of our attention will be paid to the case where the Number of equations and the number of unknown variables are the same.

# Gaussian Elimination

We begin by describing the simplest form of Gaussian elimination. In fact, it is so Simple that it is not guaranteed to proceed to completion, let alone find an accurate Solution. The modifications that will be needed to improve the “naive” method will be Introduced beginning in the next section. Three useful operations can be applied to a linear system of equations that yield An equivalent system, meaning one that has the same solutions. These operations are As follows:

1. Swap one equation for another.
2. Add or subtract a multiple of one equation from another. 
3. Multiply an equation by a nonzero constant.

> [!example]+ Example of GE
> Apply Gaussian elimination in tableau form for the system of three equations in three unknowns:
> $$
> \begin{gathered}
> x+2 y-z=3 \\
> 2 x+y-2 z=3 \\
> -3 x+y+z=-6 .
> \end{gathered}
> $$
> This is written in tableau form as
> $$
> \left[\begin{array}{rrr|r}
> 1 & 2 & -1 & 3 \\
> 2 & 1 & -2 & 3 \\
> -3 & 1 & 1 & -6
> \end{array}\right]
> $$
> Two steps are needed to eliminate column 1:
> and one more step to eliminate column 2:
> Returning to the equations
> $$
> \begin{aligned}
> x+2 y-z & =3 \\
> -3 y & =-3 \\
> -2 z & =-4,
> \end{aligned}
> $$
> we can solve for the variables
> $$
> \begin{aligned}
> x & =3-2 y+z \\
> -3 y & =-3 \\
> -2 z & =-4
> \end{aligned}
> $$
> and solve for $z, y, x$ in that order. The latter part is called back substitution, or back solving because, after elimination, the equations are readily solved from the bottom up. The solution is $x=3, y=1, z=2$ .
> 

## Generative Form

For system of equations like:
$$
\left[\begin{array}{cccccc}a_{11}&a_{12}&\ldots&a_{1n}&|&b_{1}\\ a_{21}&a_{22}&\ldots&a_{2n}&|&b_{2}\\ \vdots&\vdots&&\vdots&|&\vdots\\ a_{n1}&a_{n2}&\ldots&a_{nn}&|&b_{n}\end{array}\right]
$$

We use these update-formula to update the matrix:
$$
a_{ij}\leftarrow a_{ij}-\left(\frac{a_{i1}}{a_{11}}\right)a_{1j}\quad b_{i}\leftarrow b_{i}-\left(\frac{a_{i1}}{a_{11}}\right)b_{1},a_{11}\neq0
$$

After the elimination is completed, the tableau is upper triangular:
$$
\left[\begin{array}{cccccc}a_{11}&a_{12}&\ldots&a_{1n}&|&b_{1}\\ 0&a_{22}&\ldots&a_{2n}&|&b_{2}\\ \vdots&\vdots&&\vdots&|&\vdots\\ 0&0&\ldots&a_{nn}&|&b_{n}\end{array}\right]
$$

After elimination, the equations are readily solved from the bottom up using back solving.

$$
\begin{gathered}
x_{n}=b_{n}/a_{n n} \\
x_{n-1}=\left(b_{n-1}-a_{n-1,n}x_{n}\right)/a_{n-1,n-1} \\
x_{i}=\left(b_{i}-\sum_{j=i+1}^{n}a_{i j}x_{j}\right)/a_{i i},i=n-1,n-2,\ldots,1 
\end{gathered}
$$

# LU Factorization

Carrying the idea of tableau form one step farther brings us to the matrix form of a system of equations. Matrix form will save time in the long run by simplifying the algorithms and their analysis.

The $LU$ factorization is a matrix representation of Gaussian elimination. It consists of writing the coefficient matrix $A$ as a product of a lower triangular matrix $L$ and an upper triangular matrix $U$. The $LU$ factorization is the Gaussian elimination version of a long tradition in science and engineering—breaking down a complicated object into simpler parts.

$$
\begin{bmatrix}a_{11}&a_{12}&a_{11}&...&a_{11}\\ a_{21}&a_{22}&a_{22}&...&a_{2n}\\ a_{31}&a_{32}&a_{32}&...&a_{3n}\\ \vdots&\vdots&\vdots&\vdots&\vdots\\ a_{n2}&a_{n2}&a_{n2}&...&a_{nn}\end{bmatrix}=\begin{bmatrix}l_{11}&0&0&...&0\\ l_{21}&l_{22}&0&...&0\\ l_{31}&l_{32}&l_{53}&...&0\\ \vdots&\vdots&\vdots&\ddots&\vdots\\ l_{n1}&l_{n2}&l_{n3}&...&l_{nn}\end{bmatrix}\cdot\begin{bmatrix}u_{11}&u_{12}&u_{13}&...&u_{1n}\\ 0&u_{22}&u_{23}&...&u_{2n}\\ 0&0&u_{33}&...&u_{2n}\\ 0&\vdots&\vdots&\ddots&u_{3n}\\ 0&0&0&...&u_{nn}\end{bmatrix}
$$

- **LU matrix** 
	An $m × n$ matrix $L$ is lower triangular if its entries satisfy $l_{i ,j} = 0$ for $i < j$. An $m \times n$ matrix $U$ is upper triangular if its entries satisfy $u_{i,j} = 0$ for $i > j$

> [!example]+ Example of LU Factorization
> Find the LU factorization of
> $$
> A=\left[\begin{array}{rrr}
> 1 & 2 & -1 \\
> 2 & 1 & -2 \\
> -3 & 1 & 1
> \end{array}\right]
> $$
> This matrix is the matrix of coefficients of system (2.4). The elimination steps proceed as before:
> $$
> \begin{aligned}
> & {\left[\begin{array}{rrr}
> 1 & 2 & -1 \\
> 2 & 1 & -2 \\
> -3 & 1 & 1
> \end{array}\right] \longrightarrow \begin{array}{rr}
> \text { subtract } 2 \times \text { row } 1 \\
> \text { from row 2 }
> \end{array} \longrightarrow\left[\begin{array}{rrr}
> 1 & 2 & -1 \\
> 0 & -3 & 0 \\
> -3 & 1 & 1
> \end{array}\right]} \\
> & \longrightarrow \begin{array}{rrr}
> \text { subtract }-3 \times \text { row } 1 \\
> \text { from row } 3
> \end{array} \longrightarrow\left[\begin{array}{rrr}
> 1 & 2 & -1 \\
> 0 & -3 & 0 \\
> 0 & 7 & -2
> \end{array}\right] \\
> & \longrightarrow \quad \text { from row } 3 \longrightarrow\left[\begin{array}{rrr}
> 1 & 2 & -1 \\
> 0 & -3 & 0 \\
> 0 & 0 & -2
> \end{array}\right]=U \\
> &
> \end{aligned}
> $$
> The lower triangular $L$ matrix is formed, as in the previous example, by putting 1 's on the main diagonal and the multipliers in the lower triangle-in the specific places they were used for elimination. That is,
> $$
> L=\left[\begin{array}{rrr}
> 1 & 0 & 0 \\
> 2 & 1 & 0 \\
> -3 & -\frac{7}{3} & 1
> \end{array}\right]
> $$

## Two Forms of LU Factorization

### Doolittle Factorization

$$
\begin{bmatrix}a_{11}&a_{12}&a_{11}&...&a_{11}\\ a_{21}&a_{22}&a_{22}&...&a_{2n}\\ a_{31}&a_{32}&a_{32}&...&a_{3n}\\ \vdots&\vdots&\vdots&\vdots&\vdots\\ a_{n2}&a_{n2}&a_{n2}&...&a_{nn}\end{bmatrix}=\begin{bmatrix}1&0&0&...&0\\ l_{21}&1&0&...&0\\ l_{31}&l_{32}&1&...&0\\ \vdots&\vdots&\vdots&\ddots&\vdots\\ l_{n1}&l_{n2}&l_{n3}&...&1\end{bmatrix}.\begin{bmatrix}u_{11}&u_{12}&u_{13}&...&u_{1n}\\ 0&u_{22}&u_{23}&...&u_{2n}\\ 0&0&u_{33}&...&u_{2n}\\ \vdots&\vdots&\vdots&\ddots&\vdots\\ 0&0&0&...&u_{nn}\end{bmatrix}
$$

### Groot Factorization

$$
\begin{bmatrix}a_{11}&a_{12}&a_{11}&...&a_{11}\\ a_{21}&a_{22}&a_{22}&...&a_{2n}\\ a_{31}&a_{32}&a_{32}&...&a_{3n}\\ \vdots&\vdots&\vdots&\vdots&\vdots\\ a_{n2}&a_{n2}&a_{n2}&...&a_{nn}\end{bmatrix}=\begin{bmatrix}l_{11}&0&0&\ldots&0\\ l_{21}&l_{22}&0&\ldots&0\\ l_{31}&l_{32}&l_{33}&\ldots&0\\ \vdots&\vdots&\vdots&\ddots&\vdots\\ l_{n1}&l_{n2}&l_{n3}&\ldots&l_{n n}\end{bmatrix}\cdot\begin{bmatrix}1&u_{12}&u_{13}&\ldots&u_{1n}\\ 0&1&u_{23}&\ldots&u_{2n}\\ 0&0&1&\ldots&u_{3n}\\ \vdots&\vdots&\vdots&\ddots\\0&0&0&\ldots&\vdots\\ \end{bmatrix}
$$

The fact that not all matrices have an LU factorization means that more work is required before we can declare the LU factorization a general Gaussian elimination algorithm. The related problem of swamping is described in the next section. The $PA = LU$ factorization is introduced, which will overcome both problems.

# Sources to Error

Let $x_a$ be an approximate solution of the linear system $Ax = b$ . The residual is the vector $r = b -A x_{a}$ . The **backward error** is the norm of the residual $\|b-Ax_{a}\|_{\infty}$ , and the **forward error** is $\|x-x_{a}\|_{\infty}$ .( $\left\|A\right\|_{\infty}=\max_{1\leq i\leq n}\left(\left|a_{i1}\right|+\cdots+\left|a_{i n}\right|\right)$ )

Denote the residual by $r = b - Ax_a$ . The **relative backward error** of system $Ax =b$ is defined to be 
$$
\frac{\|r\|_{\infty}}{	\|b\|_{\infty}}=\frac{\|b-Ax_{a}\|_{\infty}}{	\|b\|_{\infty}}
$$

and the **relative forward error** is
$$
\frac{\|x-x_{a}\|_{\infty}}{	\|x\|_{\infty}}
$$

The **error magnification factor** for $Ax = b$ is the ratio of the two:
$$
\begin{align*}
eooro~magnification~factor&= \frac{relative~forward~error}{	relative~backward~error}\\
&= \frac{\frac{\|x-x_{a}\|_{\infty}}{	\|x\|_{\infty}}}{\frac{\|b-Ax_{a}\|_{\infty}}{	\|b\|_{\infty}}}
\end{align*}
$$

The **condition number** of a square matrix $A$ , $cond(A)$ , is the maximum possible error magnification factor for solving $Ax = b$ , **over all right-hand sides** $b$ . Namely:
$$
cond\left(A\right)=\max\limits_{b\in\mathbb{R}^n,b\neq0}\frac{\|x_a-x\|_\infty/\|x\|_\infty}{\|r\|_\infty/\|b\|_\infty}
$$

> For an invertible matrix $A$ , its condition number is $$cond(A)=\|A\|_{\infty}\cdot \|A^{-1}\|_{\infty}$$

^g0ff5n

## Swamping

A second significant source of error in classical Gaussian elimination is much easier to fix. We demonstrate swamping with the next example.
$$ ^y2ga3f
E.g.\left\{\begin{matrix}10^{-2}x_1+x_2=1\\ x_1+2x_2=4\end{matrix}\right.,\quad its~solution: x_2=\frac{2}{1-2\times10^{-20}},x_2=\frac{1-4\times10^{-20}}{1-2\times10^{-2}}
$$

Gaussian elimination should be kept as small as possible to avoid swamping. Fortunately, there is a simple modification to native Gaussian elimination that forces the absolute value of multipliers to be no larger than 1. This new protocol, which involves judicious row exchanges in the tableau, is called **partial pivoting**.

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/9gZ1TSU.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Partical pivoting gaussian elimination
    </div>
</center>

> [!example]+ Example of PPGE
> Apply Gaussian elimination with partial pivoting to solve the system
> $$
> \begin{aligned}
> x_1-x_2+3 x_3 & =-3 \\
> -x_1-2 x_3 & =1 \\
> 2 x_1+2 x_2+4 x_3 & =0 .
> \end{aligned}
> $$
> This example is written in tableau form as
> $$
> \left[\begin{array}{rrr|r}
> 1 & -1 & 3 & -3 \\
> -1 & 0 & -2 & 1 \\
> 2 & 2 & 4 & 0
> \end{array}\right] \text {. }
> $$
> Under partial pivoting we compare $\left|a_{11}\right|=1$ with $\left|a_{21}\right|=1$ and $\left|a_{31}\right|=2$, and choose $a_{31}$ for the new pivot. This is achieved through an exchange of rows 1 and 3:
> $$
> \begin{aligned}
> & {\left[\begin{array}{rrr|r}
> 1 & -1 & 3 & -3 \\
> -1 & 0 & -2 & 1 \\
> 2 & 2 & 4 & 0
> \end{array}\right] \rightarrow \begin{array}{c}
> \text { exchange row 1 } \\
> \text { and row 3 }
> \end{array} \longrightarrow \quad\left[\begin{array}{rrrrr}
> 2 & 2 & 4 & 0 \\
> -1 & 0 & -2 & 1 \\
> 1 & -1 & 3 & -3
> \end{array}\right]} \\
> & \longrightarrow \quad \begin{array}{rrr:r}
> \text { subtract }-\frac{1}{2} \times \text { row } 1 \\
> \quad \text { from row } 2
> \end{array} \longrightarrow\left[\begin{array}{rrr:r}
> 2 & 2 & 4 & 0 \\
> 0 & 1 & 0 & 1 \\
> 1 & -1 & 3 & -3
> \end{array}\right] \\
> & \longrightarrow \quad \begin{array}{rrr|r}
> \text { subtract } \frac{1}{2} \times \text { row } 1 \\
> \text { from row } 3
> \end{array} \longrightarrow \quad\left[\begin{array}{rrr|r}
> 2 & 2 & 4 & 0 \\
> 0 & 1 & 0 & 1 \\
> 0 & -2 & 1 & -3
> \end{array}\right] \\
> &
> \end{aligned}
> $$
> Before eliminating column 2 we must compare the current $\left|a_{22}\right|$ with the current $\left|a_{32}\right|$. Because the latter is larger, we again switch rows:
> $$
> \begin{aligned}
> & {\left[\begin{array}{rrr:r}
> 2 & 2 & 4 & 0 \\
> 0 & 1 & 0 & 1 \\
> 0 & -2 & 1 & -3
> \end{array}\right] \longrightarrow \begin{array}{rrr|r}
> \text { exchange row } 2 \\
> \text { and row } 3
> \end{array} \longrightarrow \quad\left[\begin{array}{rrr|r}
> 2 & 2 & 4 & 0 \\
> 0 & -2 & 1 & -3 \\
> 0 & 1 & 0 & 1
> \end{array}\right]} \\
> & \longrightarrow \quad \text { from row 3 } \rightarrow\left[\begin{array}{rrr:r}
> 2 & 2 & 4 & 0 \\
> 0 & -2 & 1 & -3 \\
> 0 & 0 & \frac{1}{2} & -\frac{1}{2}
> \end{array}\right] \\
> &
> \end{aligned}
> $$
> Note that all three multipliers are less than 1 in absolute value.
> The equations are now simple to solve. From
> $$
> \begin{aligned}
> \frac{1}{2} x_3 & =-\frac{1}{2} \\
> -2 x_2+x_3 & =-3 \\
> 2 x_1+2 x_2+4 x_3 & =0,
> \end{aligned}
> $$
> we find that $x=[1,1,-1]$.

# PA = LU Factorization
we put together everything we know about Gaussian elimination into the PA = LU factorization. This is the matrix formulation of elimination with partial pivoting. The $PA = LU$ factorization is the established workhorse for solving systems of linear equations.

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/87oGeZq.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Example of find the PA=LU matrix of A
    </div>
</center>

> [!example]+ Example of PA=LU to solve system of equations
> Use the $\mathrm{PA}=\mathrm{LU}$ factorization to solve the system $A x=b$ , where
> $$
> A=\left[\begin{array}{rrr}
> 2 & 1 & 5 \\
> 4 & 4 & -4 \\
> 1 & 3 & 1
> \end{array}\right], \quad b=\left[\begin{array}{l}
> 5 \\
> 0 \\
> 6
> \end{array}\right] .
> $$
It remains to complete the two back substitutions.
> 1. $L c=P b$ :
> $$
> \left[\begin{array}{rrr}
> 1 & 0 & 0 \\
> \frac{1}{4} & 1 & 0 \\
> \frac{1}{2} & -\frac{1}{2} & 1
> \end{array}\right]\left[\begin{array}{l}
> c_1 \\
> c_2 \\
> c_3
> \end{array}\right]=\left[\begin{array}{lll}
> 0 & 1 & 0 \\
> 0 & 0 & 1 \\
> 1 & 0 & 0
> \end{array}\right]\left[\begin{array}{l}
> 5 \\
> 0 \\
> 6
> \end{array}\right]=\left[\begin{array}{l}
> 0 \\
> 6 \\
> 5
> \end{array}\right]
> $$
> Starting at the top, we have
> $$
> \begin{aligned}
> c_1 & =0 \\
> \frac{1}{4}(0)+c_2 & =6 \Rightarrow c_2=6 \\
> \frac{1}{2}(0)-\frac{1}{2}(6)+c_3 & =5 \Rightarrow c_3=8 .
> \end{aligned}
> $$
> 2. $U x=c:$
> $$
> \left[\begin{array}{rrr}
> 4 & 4 & -4 \\
> 0 & 2 & 2 \\
> 0 & 0 & 8
> \end{array}\right]\left[\begin{array}{l}
> x_1 \\
> x_2 \\
> x_3
> \end{array}\right]=\left[\begin{array}{l}
> 0 \\
> 6 \\
> 8
> \end{array}\right]
> $$
> Starting at the bottom,
> $$
> \begin{aligned}
> 8 x_3=8 & \Rightarrow x_3=1 \\
> 2 x_2+2(1)=6 & \Rightarrow x_2=2 \\
> 4 x_1+4(2)-4(1)=0 & \Rightarrow x_1=-1 .
> \end{aligned}
> $$
> Therefore, the solution is $x=[-1,2,1]$.

# Iterative Methods

The generative form of iterative methods: 
$$
Qx_{k+1}=(Q-A) x_{k} + b
$$
Where $Q$ is the **Split Matrix**.
Rewrite it at fixed-point formula:
$$
\begin{align*}
form: x_{k+1}&= (I-Q^{-1}A) x_{k} + Q^{-1} b\\
Iterative~ Matrix: G&= (I-Q^{-1}A)
\end{align*}
$$

## Richardson

$$
Qx_{k+1} = (Q-A) x_{k} +b, \quad Q=I, G=I-A
$$

> [!example]+ Example of Richardson iterative method.
> $\left\{\begin{matrix}3u+v=5\\ u+2v=5\end{matrix}\right.$ , the coefficient matrix $A$ : $\begin{pmatrix}3 & 1 \\ 1 & 	2\end{pmatrix}$
> Richardson iterative method:
> $$
> \begin{bmatrix}u_{k+1}\\ v_{k+1}\end{bmatrix}=\begin{bmatrix}-2&-1\\ -1&-1\end{bmatrix}\begin{bmatrix}u_k\\ v_k\end{bmatrix}+\begin{bmatrix}5\\ 5\end{bmatrix}
> $$

-----


## Jacobi

$$
Qx_{k+1}=(Q-A)x_{k} +b \quad Q=D, G=-D^{-1}(L+U)
$$

Where $D$ is the diagonal matrix.

$$
\begin{gathered}
 D=\left[\begin{matrix}{a_{11}}&{0}&{\cdots}&{0}\\ {0}&{a_{22}}&{\ddots}&{0}\\ {\vdots}&{\ddots}&{\ddots}&{\vdots}\\ {0}&{0}&{\cdots}&{a_{n n}}\end{matrix}\right] \\
  \\
L=\left[\begin{matrix}0&0&\ldots&0\\ \alpha_{21}&0&\ldots&0\\ \vdots&\ddots&\ddots&\vdots\\ \alpha_{n1}&\alpha_{n2}&\ldots&0\end{matrix}\right],\quad U=\left[\begin{matrix}0&a_{12}&\ldots&a_{1n}\\ 0&0&\ldots&a_{2n}\\ \vdots&\vdots&\ddots&\vdots\\ 0&0&\ldots&0\end{matrix}\right] 
\end{gathered}
$$

> [!example]+ Example of Jacobi Iterative Method
> $\left\{ \begin{matrix}3u+v=5 \\ u+2v=5\end{matrix}\right.$ , the correspond coefficient matrix $A$ : $\begin{pmatrix}3	&1 \\ 1	&2\end{pmatrix}$ 
> Jacobi Iterative method: 
> $$
> \begin{align*}
> \begin{bmatrix}u_{k+1}\\ v_{k+1}\end{bmatrix}&= \begin{bmatrix}3&0\\ 0&2\end{bmatrix}^{-1}\begin{bmatrix}0&-1\\ -1&0\end{bmatrix}\begin{bmatrix}u_k\\ v_k\end{bmatrix}+\begin{bmatrix}3&0\\ 0&2\end{bmatrix}^{-1}\begin{bmatrix}5\\ 5\end{bmatrix}\\
> &=\begin{bmatrix}0&-1/3\\ -1/2&0\end{bmatrix}\begin{bmatrix}u_k\\ v_k\end{bmatrix}+\begin{bmatrix}5/3\\ 5/2\end{bmatrix}\\
> &= \begin{bmatrix}(5-v_{k})/3\\
> (5-u_{k})/2\end{bmatrix}
> \end{align*}
> $$

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/vshy7go.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Relation between Jacobi and Richardson method
    </div>
</center>

## Gauss-Seidel

$$
Qx_{k+1}=(Q-A) x_{k}+b,\quad Q=D+L,G=-(L+D)^{-1}U
$$

> [!example]+ Example of Gauss-Seidel Method 
> $\left\{\begin{matrix}3u+v=5 \\ u+2v=5\end{matrix}\right.$ , the coefficient matrix $A$ : $\begin{bmatrix}3 & 1 \\ 1 & 	2\end{bmatrix}$
> Gauss-Seidel iterative method:
> $$
> \begin{bmatrix}3&0\\ 1&2\end{bmatrix}\begin{bmatrix}u_{k+1}\\ v_{k+1}\end{bmatrix}=\begin{bmatrix}0&-1\\ 0&0\end{bmatrix}\begin{bmatrix}u_{k}\\ v_{k}\end{bmatrix}+\begin{bmatrix}5\\ 5\end{bmatrix}\Rightarrow\begin{bmatrix}u_{k+1}\\ v_{k+1}\end{bmatrix}=\begin{bmatrix}\frac{5-v_{k}}{3}\\ \frac{5-u_{k+1}}{2}\end{bmatrix}
> $$
> 

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/5jmThHR.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Connection between Gauss-Seidel and Jacobi method
    </div>
</center>

## SOR

For generative linear iterative method: 
$$
x_{k+1}=Gx_{k}+c
$$
We can get a new iterative method by using relaxation method:
$$
x_{k+1}=\gamma(Gx_k+c)+(1-\gamma)x_k,\gamma\neq0
$$
Namely:

$$
x_{k+1} =G_{\gamma}x_{k} +\gamma c,\quad G_{\gamma}=\gamma G + (1- \gamma)I
$$

SOR iterative method:

$$
Qx_{k+1}=(Q-A) x_{k}+b,\quad Q=\frac{1}{	\omega}(D+\omega L)
$$

> [!example]+ Example of Gauss-Seidel Method 
> $\left\{\begin{matrix}3u+v=5 \\ u+2v=5\end{matrix}\right.$ , the coefficient matrix $A$ : $\begin{bmatrix}3 & 1 \\ 1 & 	2\end{bmatrix}$
> Gauss-Seidel iterative method:

SOR iterative method:
$$
\begin{align*}\\
Gauss-Seidel: \begin{bmatrix}u_{k+1}\\ v_{k+1}\end{bmatrix}&= \begin{bmatrix}\frac{5-v_k}{3}\\ \frac{5-u_{k+1}}{2}\end{bmatrix}\\
SOR: \begin{bmatrix}u_{k+1}\\ v_{k+1}\end{bmatrix}&= (1-\omega)\begin{bmatrix}u_{k}\\ v_{k}\end{bmatrix}+\omega\begin{bmatrix}\frac{5-v_{k}}{3}\\ \frac{5-u_{k+1}}{2}\end{bmatrix}
\end{align*}
$$

## Convergence of Iterative Methods

If the iterative matrix $G$ 's norm (of any kind) $\|G\|= \delta < 1$ , then $x_{k+1} = (I-Q^{-1}A)x_{k} + Q^{-1}b$ is convergent to the solution of $Ax=b$ , also:

$$
\|x-x_{k}\|\leq\frac{\delta}{1-\delta}\|x_{k}-x_{k-1}\|\leq\cdots\frac{\delta^{k}}{1-\delta}\|x_{1}-x_{0}\|
$$

> Above is exactly what strict diagonal Dominance implies 

For $Ax=b$ , if $A$ is a strict diagonal matrix, then we have:

1. If $A$ is **Dominance**, Richardson iterative method is convergent.
2. If $A$ is **strict row diagonal matrix**, Jacobi iterative method is convergent.

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/fi368wW.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Convergence of Richardson iterative method.
    </div>
</center>

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/7jQGas1.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Convergence of Jacobi iterative method.
    </div>
</center>

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/8qpoBKi.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">spectral radius method (spectral radius=$\sqrt{\max{\lambda_i}}$)
    </div>
</center>

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/dcOWYve.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Convergence Judgment Using Spectral Radius method.
    </div>
</center>

> [!example]+ Example of Spectral Radius Method 
> 
<center>
    <img style="border-radius: 0.3125 em;
    Box-shadow: 0 2 px 4 px 0 rgba (34,36,38,. 12), 0 2 px 10 px 0 rgba (34,36,38,. 08);"
    src=" https://search.pstatic.net/common?src=https://i.imgur.com/pR42Zic.png">
    <br>
    <div style="color: orange; border-bottom: 1 px solid #d9d9d9 ;
    Display: inline-block;
    color: #999 ;
    padding: 2px;">fix: $\rho(G)>1$
    </div>
</center>

# Methods for Symmetric Positive-definite Matrices

- **Definition**(Symmetric positive-definite matrixes)
	The matrix $A \in \mathbb{R}^{n \times n}$ is *symmetric* if $A^{\top}=A$ ; A is *positive-definite* if $x^{\top} A x > 0$ for all vector $x \neq 0$.

If the $n \times n$ matrix $A$ is symmetric, then $A$ is positive-definite if and only if all of its **eigenvalues** are positive.

Property
+ Symmetric positive-definite matrixes is nonsingular.
+ If $A$ is symmetric positive-definite matrix, $P$ is a full-rank $n \times m$ matrix, then $P^{\top}AP$ is also a symmetric positive-definite matrix.

## Cholesky Factorization

- **Cholesky Factorization Theorem** 
	If $A$ is a symmetric positive-definite $n \times n$ matrix, then there exists an upper triangular $n \times n$ matrix $R$ such that $A=R^{\top}R$ .

> [!attention]+ Cholesky factorization
> ```
> for k =1, 2, ... n
> 	if A(k,k) < 0, stop, end
> 	R(k,k)=math.sqrt(A(k,k))
> 	u=A(k, k+1:n)/R(k,k)
> 	R(k, k+1:n)=u
> 	A(k+1: n,k+1:n)-=np.dot(np.reverse(u),u)
> end
> ```

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/FjRovG3.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Numerical Example
    </div>
</center>

## SOR Iterative Method of Symmetric Positive-definite Matrix

> If $A$ is a symmetric positive-definite matrix, then the SOR iterative method for solving $Ax=b$ is convergent for any initial guess.

## Conjugate Gradient Method

- **A-Conjugate**
	Let $A$ be a symmetric positive-definite $n \times n$ matrix. For two n-vectors $v$ and $w$, define the A-inner product $$(v,w)_A = v^{\top}Aw$$The vectors $v$ and $w$ are A-conjugate if $(v,w)_A = 0$.

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/8HxHCKS.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">code of CGM
    </div>
</center>

