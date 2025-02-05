---
date: 2023-04-17 09:28:29
categories: Numerical_Analysis
destination: 
excerpt: 奇异值分解和特征值求解
katex: true
obsidianUIMode: source
rating: ⭐⭐
draft: false
tags:
  - Numerical_Analysis
title: "Chapter7-Eigenvalues and singular values"
share: true
updated: 2023-06-14 20:30:33
---

# Eigenvalues and Singular Values

## Gauss of Eigenvalues

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/JoFsDBD.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Gershgorin theorem
    </div>
</center>

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/a7uaYDz.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Gershgorin theorem 2
    </div>
</center>

> [!example]+ example of eigenvalues guess using Gershaorin theorem 
<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/TIHELhx.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Example
    </div>
</center>

## Power Iteration

The motivation behind Power Iteration is that multiplication by a matrix tends to move vectors toward the dominant eigenvector direction.

> Power Iteration is essentially a fixed-point iteration with normalization at each step. Like FPI, it converges linearly at rate of $|\frac{\lambda_{2}}{	\lambda_{1}}|$

Let A be an $n \times n$ matrix with real eigenvalues $\lambda_{1}, \cdots, \lambda_{n}$ satisfying $\lambda_{1} \geqslant \lambda_{2}, \cdots, \geqslant \lambda_{n}$ . Assume that the eigenvectors of $A$ $(v_{1},v_{2},\cdots,v_{n})$ span $\mathbb{R}^{n}$ . 

For almost every initial vector, Power Iteration converges linearly to an eigenvector associated to $\lambda_{1}$ with convergence rate constant $S=|\frac{\lambda_{2}}{	\lambda_{1}}|$ .

$$
\begin{aligned}
x_{k}&=A^{k}x_{0}=A^{k}\left(c_{1}v_{1}+c_{2}v_{2}\cdots+c_{n}v_{n}\right) \\
&=\lambda_{1}^{k}\left[c_{1}v_{1}+\left(\frac{\lambda_{2}}{\lambda_{1}}\right)^{k}c_{2}v_{2}+\cdots+\left(\frac{\lambda_{n}}{\lambda_{1}}\right)^{k}c_{n}v_{n}\right] \\
\end{aligned}
$$
1. $|\lambda_{1}|>|\lambda_{2}|\geq\cdots\geq|\lambda_{n}|$

$$
\begin{align*}
\lim_{k\to\infty}x_{k}&=\lambda_{1}^{k}\left[c_{1}v_{1}+\left(\frac{\lambda_{2}}{\lambda_{1}}\right)^{k}c_{2}v_{2}\right],c_{1}\neq0\\
x_{k} & \approx \lambda_{1}x_{k-1}\\
S &\leqslant |\lambda_{2}/\lambda_{1}|
\end{align*}
$$

2. $|\lambda_{1}|=|\lambda_{2}|>\cdots\geq|\lambda_{n}|,\lambda_{1}=-\lambda_{2}$

$$
\begin{align*}
\lim\limits_{k\to\infty}x_k&=\lambda_1^k\left[\underbrace{c_1v_1+(-1)^kc_2v_2}_{w}+\left(\frac{\lambda_3}{\lambda_1}\right)^kc_3v_3\right],w\neq0\\
x_{k} & \approx \lambda_{1}^{2}x_{k-2}\\
S & \leqslant |\lambda_{3}/ \lambda_{1}|
\end{align*}
$$

> [!example]+ example of Power iteration method
> Suppose we have matrix $A$ in which we want to get its max eigenvalue：
> $$
> A = \begin{bmatrix}
> 4 & -1 & 1 \\
> -1 & 3 & -2 \\
> 1 & -2 & 3 \\
> \end{bmatrix}
> $$
> 
> Firstly, choose a initial guess vector $\mathbf{x}^{(0)}$ ，Say $\mathbf{x}^{(0)} = \begin{bmatrix} 1 \\ 1 \\ 1 \end{bmatrix}$ 。
> 
> The power iteration method:：
> 
> 1. Compute $\mathbf{y}^{(k)} = A \mathbf{x}^{(k-1)}$
> 2. Compute $\lambda^{(k)} = \frac{\mathbf{y}^{(k)}}{\mathbf{x}^{(k-1)}}$
> 3. Normalization $\mathbf{x}^{(k)} = \frac{\mathbf{y}^{(k)}}{\| \mathbf{y}^{(k)} \|}$
> 
> Repeat the following step until it satisfies some rule. (For example， $\| \mathbf{x}^{(k)} - \mathbf{x}^{(k-1)} \| < TOL$ )
> 
> Choose a initial guess vector: $\mathbf{x}^{(0)} = \begin{bmatrix} 1 \\ 1 \\ 1 \end{bmatrix}$ 。
> 
> ***Step 1***
> $$
> \mathbf{y}^{(1)} = A \mathbf{x}^{(0)} = \begin{bmatrix} 4 & -1 & 1 \\ -1 & 3 & -2 \\ 1 & -2 & 3 \end{bmatrix} \begin{bmatrix} 1 \\ 1 \\ 1 \end{bmatrix} = \begin{bmatrix} 4 \\ 0 \\ 2 \end{bmatrix}
> $$
> 
> ***Step 2***
> $$
> \lambda^{(1)} = \frac{\mathbf{y}^{(1)}}{\mathbf{x}^{(0)}} = \frac{\begin{bmatrix} 4 \\ 0 \\ 2 \end{bmatrix}}{\begin{bmatrix} 1 \\ 1 \\ 1 \end{bmatrix}} = 4
> $$
> ***Step 3***
> Normalization  $$\mathbf{x}^{(1)} = \frac{\mathbf{y}^{(1)}}{\| \mathbf{y}^{(1)} \|} = \frac{\begin{bmatrix} 4 \\ 0 \\ 2 \end{bmatrix}}{\sqrt{4^2 + 0^2 + 2^2}} = \begin{bmatrix} \frac{2}{\sqrt{5}} \\ 0 \\ \frac{1}{\sqrt{5}} \end{bmatrix}$$
> 
> Repeat:
> $$
> \mathbf{y}^{(2)} = A \mathbf{x}^{(1)} = \begin{bmatrix} 4 & -1 & 1 \\ -1 & 3 & -2 \\ 1 & -2 & 3 \end{bmatrix} \begin{bmatrix} \frac{2}{\sqrt{5}} \\ 0 \\ \frac{1}{\sqrt{5}} \end{bmatrix} = \begin{bmatrix} \frac{6}{\sqrt{5}} \\ \frac{2}{\sqrt{5}} \\ \frac{4}{\sqrt{5}} \end{bmatrix}
> $$
> 
> $$
> \lambda^{(2)} = \frac{\mathbf{y}^{(2)}}{\mathbf{x}^{(1)}} = \frac{\begin{bmatrix} \frac{6}{\sqrt{5}} \\ \frac{2}{\sqrt{5}} \\ \frac{4}{\sqrt{5}} \end{bmatrix}}{\begin{bmatrix} \frac{2}{\sqrt{5}} \\ 0 \\ \frac{1}{\sqrt{5}} \end{bmatrix}} = 5
> $$
> 
> Normalization  $\mathbf{x}^{(2)} = \frac{\mathbf{y}^{(2)}}{\| \mathbf{y}^{(2)} \|} = \frac{\begin{bmatrix} \frac{6}{\sqrt{5}} \\ \frac{2}{\sqrt{5}} \\ \frac{4}{\sqrt{5}} \end{bmatrix}}{\sqrt{\left(\frac{6}{\sqrt{5}}\right)^2 + \left(\frac{2}{\sqrt{5}}\right)^2 + \left(\frac{4}{\sqrt{5}}\right)^2}} = \begin{bmatrix} \frac{3}{\sqrt{15}} \\ \frac{1}{\sqrt{15}} \\ \frac{2}{\sqrt{15}} \end{bmatrix}$
> 
>  Repeat step above until it stops. Finally the maximum eigenvalue is $\lambda = 5$ ，its corresponding eigenvector  $\mathbf{x} = \begin{bmatrix} \frac{3}{\sqrt{15}} \\ \frac{1}{\sqrt{15}} \\ \frac{2}{\sqrt{15}} \end{bmatrix}$  ^k6vidt

### Other Kinds of Power Iteration

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/6izQDLC.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Other Kinds of Power Iteration
    </div>
</center>


> [!example]+ example of 位移逆幂法
> 假设我们要求解矩阵
> $$
> A = \begin{bmatrix}
> 1 & 2 & 3 \\
> 4 & 5 & 6 \\
> 7 & 8 & 9 \\
> \end{bmatrix}
> $$
> 的最大特征值和对应的特征向量。
> 
> 具体步骤如下：
> 
> 1. 初始化一个初始向量 $x^{(0)}$，可以是任意非零向量。
> 
> 2. 对于每一次迭代 $k$，计算位移矩阵 $B = A - \mu I$，其中 $\mu$ 是一个位移量，通常选择接近所求特征值的值。
> 
> 3. 解线性方程组 $Bx^{(k+1)} = x^{(k)}$，得到新的迭代向量 $x^{(k+1)}$。
> 
> 4. 归一化迭代向量 $x^{(k+1)}$，即 $x^{(k+1)} = \frac{x^{(k+1)}}{\|x^{(k+1)}\|}$。
> 
> 5. 判断迭代是否收敛，如果满足收敛条件（如 $|x^{(k+1)} - x^{(k)}| < \epsilon$），则停止迭代，否则返回步骤 2。
> 
> 下面我们以 $\mu = 10$ 为例，进行迭代计算。
> 
> 首先，选择一个初始向量 $x^{(0)} = [1, 1, 1]^T$。
> 
> 迭代计算如下：
> 
> $$
> \begin{align*}
> k = 0: & \quad B = A - \mu I = \begin{bmatrix} -9 & 2 & 3 \\ 4 & -5 & 6 \\ 7 & 8 & -1 \end{bmatrix}, \quad x^{(1)} = B^{-1}x^{(0)} = \begin{bmatrix} -0.081 & -0.027 & 0.108 \end{bmatrix}^T \\
> k = 1: & \quad B = A - \mu I = \begin{bmatrix} -9 & 2 & 3 \\ 4 & -5 & 6 \\ 7 & 8 & -1 \end{bmatrix}, \quad x^{(2)} = B^{-1}x^{(1)} = \begin{bmatrix} -0.058 & -0.086 & 0.142 \end{bmatrix}^T \\
> k = 2: & \quad B = A - \mu I = \begin{bmatrix} -9 & 2 & 3 \\ 4 & -5 & 6 \\ 7 & 8 & -1 \end{bmatrix}, \quad x^{(3)} = B^{-1}x^{(2)} = \begin{bmatrix} -0.039 & -0.112 & 0.187 \end{bmatrix}^T \\
> & \quad \vdots \\
> \end{align*}
> $$
> 
> 请注意，迭代过程中需要对迭代向量进行归一化，以保证迭代结果的正确性。


### Rayleigh Power Iteration

The Rayleigh quotient can be used in conjunction with Inverse Power Iteration. We know that it converges to the eigenvector associated to the eigenvalue with the smallest distance to the shift s, and that convergence is fast if this distance is small. If at any step along the way an approximate eigenvalue were known, it could be used as the shift s, to speed convergence.

> [!summary]+ Procedure of Rayleigh Power Iteration
> Given initial vector $x_0$ 
> $$
> \begin{align*}
> for \quad j&= 1,2,3,\cdots\\
> u_{j-1}&= x_{j-1} / \|x_{j-1}\|\\
> \lambda_{j-1}&= u_{j-1}^{\top}A u_{j-1}\\
> Solve (A-\lambda_{j-1}I)x_{j}&= u_{j-1}\\
> end \quad &.\\
> u_{j}=x_{j}/\|x_{j}\|_{2}
> \end{align*}
> $$

> If A is real symmetric matrix, then Rayleigh power iteration's convergence speed is $S \leqslant |\lambda_{2}/\lambda_{1}|^{2}$

## QR Decomposition

The goal of this section is to develop methods for finding all eigenvalues at once. We begin with a method that works for symmetric matrices, and later supplement it to work in general. Symmetric matrices are easiest to handle because their eigenvalues are real and their unit eigenvectors form an orthonormal basis of $\mathbb{R}^{n}$ . This motivates applying Power Iteration with $n$ vectors in parallel, where we actively work at keeping the vectors orthogonal to one another.

## Reduced QR Decomposition

If matrix $A$ whose column vectors $x_{1}, x_{1},x_{3},\cdots,x_{n}$ is linear independent, then we can apply Gram-Schmidt Orthogonalization.

### Gram-Schmidt Orthogonalization

$$
x_1, x_2, x_3, \cdots, x_n \stackrel{\text { schmidt }}{\longrightarrow}\left\{\begin{array}{c}
v_1=x_1 \\
v_2=x_2-\frac{x_2 \cdot v_1}{v_1 \cdot v_1} v_1 \\
v_3=x_3-\frac{x_3 \cdot v_1}{v_1 \cdot v_1} v_1-\frac{x_3 \cdot v_2}{v_2 \cdot v_2} v_2 \\
\vdots \\
v_n=x_n-\frac{x_n \cdot v_1}{v_1 \cdot v_1} v_1-\cdots-\frac{x_n \cdot v_{n-1}}{v_{n-1} \cdot v_{n-1}} v_{n-1}
\end{array}\right.
$$

Then rewrite it in matrix form:
$$
\begin{align*}
A&= QR\\
\begin{bmatrix}x_1 & x_{2} & \cdots & x_{n}\end{bmatrix}&= \begin{bmatrix}v_{1}& v_{2} & \cdots & v_{n}\end{bmatrix}\begin{bmatrix}1 & \frac{x_{2} \cdot v_{1}}{v_{1} \cdot v_{1}} & \frac{x_{3} \cdot v_{1}}{v_{1} \cdot v_{1}} & \cdots & \frac{x_{n} \cdot v_{1}}{v_{1} \cdot v_{1}} \\ & 1 & \frac{x_{2} \cdot v_{2}}{v_{2} \cdot v_{2}} & \cdots & \frac{x_{n} \cdot v_{2}}{v_{2} \cdot v_{2}} \\ & & 1 & \ddots & \frac{x_{n} \cdot v_{3}}{v_{3} \cdot v_{3}} \\ & & & \ddots & 	1\end{bmatrix}
\end{align*}
$$

Actually, to **avoid swamping** as we talk in [Chapter 3 -System of Equation](term/Numerical_Analysis/Chapter%203%20-System%20of%20Equation.md), one alternative is to **normalize the orthogonal vectors**. See the example below.

> [!example]+ Example of QR factorization
> $$
> A=\left[\begin{array}{cc}
> 1 & -4 \\
> 2 & 3 \\
> 2 & 2
> \end{array}\right]
> $$
> 1. $$y_1=A_1=[1,2,2]^T, r_{11}=\left\|y_1\right\|_2=3, q_1=\frac{y_1}{\left\|y_1\right\|_2}=\left[\frac{1}{3}, \frac{2}{3}, \frac{2}{3}\right]^T$$
> 2. $$r_{12}=q_1^T A_2=2, y_2=A_2-r_{12} q_1=\left[\begin{array}{c}-4 \\3 \\2\end{array}\right]-2\left[\begin{array}{c}\frac{1}{3} \\\frac{2}{3} \\\frac{2}{3}\end{array}\right]=\left[\begin{array}{c}-\frac{14}{3} \\\frac{5}{3} \\\frac{2}{3}\end{array}\right]$$
> 3. $$r_{22}=\left\|y_2\right\|_2=5, q_2=\frac{y_2}{\left\|y_2\right\|_2}=\left[-\frac{14}{15}, \frac{1}{3}, \frac{2}{15}\right]^T$$
> 4. 
> $$
> A=\left[\begin{array}{cc}
> \frac{1}{3} & -\frac{14}{15} \\
> \frac{2}{3} & \frac{1}{3} \\
> \frac{2}{3} & \frac{2}{15}
> \end{array}\right]\left[\begin{array}{ll}
> 3 & 2 \\
> 0 & 5
> \end{array}\right]=Q R
> $$ 
> ^j17wmv

## Full QR Decomposition

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/UIbCalo.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Difference between the two QR factorization
    </div>
</center>

The modified row and column can be add by adding new linearly independent unit vector and `np.zeros` rows to the matrix $R$ :

> [!example]+ Example of Full QR Factorization.
<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/GKLAyfF.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Full QR Factorization
    </div>
</center>

### Modified Gram-Schmidt Orthogonalization
The Classical Gram-Schmidt Orthogonalization may encounter some problem in real world problems such as:
1. Not suitable for sparse matrix;
2. There will be a lot of errors in the ill-conditioned matrix with finite precision;
3. The whole matrix must be saved during the loop, which leads to high memory overhead.

See the Example below:
> [!example]+ Example When Classical Gram-Schmidt Orthogonalization Failed
<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/qyPCUph.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">The case of ill-conditioned matrix
    </div>
</center>

The main improvement of MGS is to select a vector as the first datum, and then project it and all the vectors into the orthogonal space of the datum, the first datum is then perpendicular to all the vectors in the orthogonal space. So, after that, we just have to do the same thing for the rest of the vectors in the orthogonal space, and then all the vectors are orthogonal to each other. 

> [!seealso]+ MGS
<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/BV77ka9.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">MGS Procedure
    </div>
</center>

### Householder Transformation

Since the transformation matrix $H$ in Householder transformation is symmetric orthogonal matrix, then we can replace $Q$ with $H$ . See the example below.

> [!example]+ Householder transformation in QR factorization.
<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/8XzsuUJ.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Use Householder transformation to apply QR factorization.
    </div>
</center>

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/2jYuPC1.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Householder Transformation 
    </div>
</center>

## QR Algorithm

The way the QR algorithm finds eigenvalues of a matrix A is to locate a similar matrix whose eigenvalues are obvious. An example of the latter is **real Schur form**.

**Real Schur form:** A matrix $T$ has *real Schur form* if it is upper triangular, except possibly for $2 \times 2$ blocks on the main diagonal. For example, a matrix of the form
$$
\begin{bmatrix}x & x & x & x & x \\ & x & x & x & x \\ & & x & x & x \\ & & & & 	x\end{bmatrix}
$$
Has real Schur form.

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/S3hpMBj.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Theorems about QR algorithm
    </div>
</center>

Let's consider the matrix A:


$$
A =
\begin{bmatrix}
2 & 1 & 0 \\
1 & 3 & 1 \\
0 & 1 & 2 \\
\end{bmatrix}
$$

To find the eigenvalues of $A$ using the QR algorithm, we follow these steps:

1. Start by initializing $A_0$ as the given matrix $A$.

2. Compute the QR decomposition of $A_0$. Let $Q_0$ and $R_0$ be the resulting orthogonal and upper triangular matrices, respectively.

3. Update $A_1=R_0 \cdot Q_0$ .

4. Repeat steps 2 and 3 until convergence or a desired number of iterations. Let $A_k$ be the matrix at the k-th iteration.

5. The eigenvalues of A are the diagonal elements of the final upper triangular matrix $A_k$.

Let's go through the steps for our example:

***Step 1***:
$$
A_0 = 
\begin{bmatrix}
2 & 1 & 0 \\
1 & 3 & 1 \\
0 & 1 & 2 \\
\end{bmatrix}
$$

***Step 2***:
Compute the QR decomposition of $A_0$:
$$A_0 = Q_0 \cdot R_0$$

***Step 3***:
Update $A_1 = R_0 \cdot Q_0$

***Step 4*** :
Repeat steps 2 and 3 until convergence or a desired number of iterations.

After a few iterations, we obtain the following matrices:

$$
A_1 = 
\begin{bmatrix}
3.317 & 2.317 & 1.317 \\
0 & 1.317 & 0.317 \\
0 & 0 & -0.634 \\
\end{bmatrix}
$$

$$
A_2 = 
\begin{bmatrix}
2.585 & 1.585 & 0.585 \\
0 & 1.585 & 0.585 \\
0 & 0 & -0.17 \\
\end{bmatrix}
$$

$$
A_3 = 
\begin{bmatrix}
2.414 & 1.414 & 0.414 \\
0 & 1.414 & 0.414 \\
0 & 0 & -0.17 \\
\end{bmatrix}
$$

***Step 5***:
The eigenvalues of A are the diagonal elements of the final upper triangular matrix $A_k$. In this case, the eigenvalues are:

$$
\lambda_1 = 2.414, \lambda_2 = 1.414, \lambda_3 = -0.17
$$

So, these are the eigenvalues of the given matrix $A$ using the QR algorithm.
