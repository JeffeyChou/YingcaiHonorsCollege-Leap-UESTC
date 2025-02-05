---
date: '2023-06-14 16:26:41'
categories: Numerical_Analysis 
destination: 
excerpt: 
katex: true
obsidianUIMode: source
rating: ⭐⭐
tags:  
- Numerical_Analysis 
- Interpolation
title: "Chapter4-Polynomial Interpolation "
share: true
---

# Polynomial Interpolation 

## Lagrange Interpolation 

Assume that n data points $(x_1, y_1),\cdots, (x_n, y_n)$ are given, and that we would like to find an interpolating polynomial. There is an explicit formula, called the Lagrange interpolating formula, for writing down a polynomial of degree $d = n − 1$ that interpolates the points. For example, suppose that we are given three points $(x_1, y_1), (x_2, y_2), (x_3, y_3)$ . Then the polynomial:

$$
\begin{aligned}
P_{2}\left(x\right)=y_1\frac{(x-x_2)(x-x_3)}{(x_1-x_2)(x_1-x_3)}+y_2\frac{(x-x_1)(x-x_3)}{(x_2-x_1)(x_2-x_3)}+y_3\frac{(x-x_1)(x-x_2)}{(x_3-x_1)(x_3-x_2)} 
\end{aligned}
$$

The main theorem behind polynomial interpolation is that:

> Let $(x_1, y_1),\cdots, (x_n, y_n)$ be $n$ points in the plane with distinct $x_i$ . Then there exists one and only on polynomial $P$ of degree $n-1$ or less that satisfies $P (x_i)=y_{i}$ for $i=1, \cdots, n$

> [!example]+ Example of Lagrange Interpolation 
<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/VesJaTr.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Example of Lagrange Interpolation 
    </div>
</center>



## Newton's Divided Differences For Interpolation.


<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/i6ZRRZ1.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Properties of Newton's divided differences
    </div>
</center>

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/gjSEziF.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">example of Newton's Interpolation 
    </div>
</center>

# Error of Polynomial Interpolation 

Assume that $P (x)$ is the (degree $n$ or less) interpolating polynomial fitting the $n$ points $(x_0, y_0),\cdots, (x_n, y_n)$ . The interpolation error is
$$
f(x)-P(x)=\frac{1}{	(n+1)!}f^{(n+1)}(c) \prod_{i=0}^{n}(x-x_{i})
$$
Where $c \in (a,b)$

## Chebyshev Interpolation 

It is common to choose the base points xi for interpolation to be evenly spaced. In many cases, the data to be interpolated are available only in that form—for example, when the data consist of instrument readings separated by a constant time interval. It turns out that the choice of base point spacing can have a significant effect on the interpolation error. Chebyshev interpolation refers to a particular optimal way of spacing the points.

The motivation for Chebyshev interpolation is to improve control of the maximum value of the interpolation error. 

### Chebyshev Interpolation Nodes

On the interval $[a,b]$ ,
$$
x_i=\frac{b+a}{2}+\frac{b-a}{2}\cos\frac{(2i+1)\pi}{2n+1}
$$
For $i=0,2,\cdots, n$ . The inequality
$$
\begin{align*}
|(x-x_{0})\cdots(x-x_{n})|&\leq{\frac{\left({\frac{b-a}{2}}\right)^{n+1}}{2^{n}}}\\
f(x)-P(x) & \leqslant \frac{1}{	(n+1)!}f^{(n+1)}(c)\frac{\left({\frac{b-a}{2}}\right)^{n+1}}{2^{n}}
\end{align*}
$$
Holds on $[a, b]$ .

> [!example]+ Example of Chebyshev Interpolation Error 
<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/xBQRq8X.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Example of Chebyshev Interpolation Error
    </div>
</center>

# Hermite Interpolation 

In real world problem, we not only want the interpolation function $p(x)$ go through all the points, but also want its differentials satisfy some value boundary conditions. Like find $P(x)$ that satisfies $P(x_{i})=f(x_{i}), p^{'}(x_{i})=f^{'}(x_{i}), p^{''}(x_{i})=f^{''}(x_{i})$

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/9oxdBVq.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Hermite Interpolation
    </div>
</center>

## Hermite Interpolation With Lagrange Method
For given $n$ points $\left\{x_i\right\}_{i=0}^n$ , find polynomial $P(x)$ that satisfies: $P\left(x_i\right)=$ $f\left(x_i\right), P^{\prime}\left(x_i\right)=f^{\prime}\left(x_i\right)$ 。
$$
P(x)=\sum_{i=0}^{2 n+1} a_i \varphi_i(x)=\sum_{i=0}^n a_i A_i(x)+\sum_{i=0}^n b_i B_i(x)
$$

$$
\left[\begin{array}{cccccc}
A_0\left(x_0\right) & \cdots & A_n\left(x_0\right) & B_0\left(x_0\right) & \cdots & B_n\left(x_0\right) \\
\vdots & \cdots & \vdots & \vdots & \cdots & \vdots \\
A_0\left(x_n\right) & \cdots & A_n\left(x_n\right) & B_0\left(x_n\right) & \cdots & B_n\left(x_n\right) \\
A_0^{\prime}\left(x_0\right) & \cdots & A_n^{\prime}\left(x_0\right) & B_0^{\prime}\left(x_0\right) & \cdots & B_n^{\prime}\left(x_0\right) \\
\vdots & \cdots & \vdots & \vdots & \cdots & \vdots \\
A_0^{\prime}\left(x_n\right) & \cdots & A_n^{\prime}\left(x_n\right) & B_0^{\prime}\left(x_n\right) & \cdots & B_n^{\prime}\left(x_n\right)
\end{array}\right]\left[\begin{array}{c}
a_0 \\
\vdots \\
a_n \\
b_0 \\
\vdots \\
b_n
\end{array}\right]=\left[\begin{array}{c}
f\left(x_0\right) \\
\vdots \\
f\left(x_n\right) \\
f^{\prime}\left(x_0\right) \\
\vdots \\
f^{\prime}\left(x_n\right)
\end{array}\right]
$$

Solving it we have:
$$
\begin{align*}
A_{i}\left(x\right)&= \left[1-2l_{i}^{\prime}\left(x_{i}\right)\left(x-x_{i}\right)\right]\left[l_{i}\left(x\right)\right]^{2}\\
B_{i}\left(x\right)&= \left(x-x_{i}\right)\left[l_{i}\left(x\right)\right]^{2}\\
where~ 
l_{i}\left(x\right)& =\prod_{j=0}^{n}\frac{x-x_{j}}{x_{i}-x_{j}}  \\
\end{align*}
$$
> [!example]+ example of 2-order Hermite Interpolation 

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/8G794wj.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">2-point Hermite Interpolation
    </div>
</center>

## Hermite Interpolation with Newton's Differences 

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/AcZnU6A.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Hermite interpolation with newton's differences
    </div>
</center>

See the example below.

> [!example]+ example of Hermite interpolation with Newton's differences 
<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/CyXxc5J.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">example of Hermite interpolation with Newton's differences
    </div>
</center>

## Hermite Interpolation Error 

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/Fhw74QK.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Hermite Interpolation Error 
    </div>
</center>


# Natural Cubic Splines Interpolation 

Splines represent an alternative approach to data interpolation. In polynomial interpolation, a single formula, given by a polynomial, is used to meet all data points. The idea of splines is to use several formulas, each a low degree polynomial, to pass through the data points.

## Properties of Splines 

$$
\begin{align*}
&S_{i}\left(x_i\right)=f\left(x_{i}\right), S_{i}(x_{i+1})=f(x_{i+1})\\
&S_{i}(x_{i+1})=f(x_{i+1}){\mathrm{~for~}}i=0,\ldots, n-1. \\
&S_{i}^{\prime}(x_{i+1})=S_{i+1}^{\prime}(x_{i+1})\mathrm{for}i=0,\ldots, n-2. \\
&S_{i}^{\prime\prime}(x_{i+1})=S_{i+1}^{\prime\prime}(x_{i+1})\mathrm{for}i=0,\ldots, n-2.
\end{align*}
$$


For all $S_{i}$ :

$$
S_{i}\left(x\right)=f\left(x_{i}\right)+b_{i}\left(x-x_{i}\right)+c_{i}\left(x-x_{i}\right)^{2}+d_{i}\left(x-x_{i}\right)^{3}\left(i=0,...,n-1\right)
$$


## Example

> [!example]+ example of natural cubic splines interpolation
> 给定三个插值节点：(-1,1), (1,1), (2,4)，设插值函数为 S (x)。
> 1. 对应每一个小区间内的插值函数： $S_{i}\left(x_i\right)=f\left(x_{i}\right), S_{i}(x_{i+1})=f(x_{i+1})$
> 3. 除边界点外每个节点的一阶导数： $S_{i}'(x_{i+1}) = S_{i+1}'(x_{i+1})$
> 4. 除边界点每个节点的二阶导数: $S_{i}''(x_{i+1}) = S_{i+1}''(x_{i+1})$
> 5. 边界条件： $S_0''(x_{0}) = 0 \quad S_{n-1}''(x_{n}) = 0$
> 6. 构造插值函数
> $S_0(x) = a_0(x +1)^3 + b_0(x +1)^2 + c_0(x +1) + 1, -1 \leqslant x \leqslant 1$
> 
> $S_1(x) = a_1(x - 1)^3 + b_1(x - 1)^2 + c_1(x - 1) + 1, 1 \leqslant x \leqslant 2$
> 
> 7. 解方程组，将以上节点代入方程组解得结果为：
> $a_{0}=-1, b_{0}=0, c_{0}=1/4;a_{1}=2, b_{1}=3/2, c_{1}=-1/2$
> 
