---
cover: 
categories: Numerical_Analysis
date: 2023-03-06 09:00:07
destination: 10-blog/source/_posts/数值分析
excerpt: 第二章-求解线性方程
katex: true
obsidianUIMode: source
rating: ⭐
draft: false
status: complete
tags:  
- Numerical_Analysis
title: "Chapter2-Solving nonlinear equations"
share: true
updated: 2023-06-12 20:54:30
---

Continue our journey in [Chapter1-Fundanentals](term/Numerical%20Analysis/Chapter1-Fundanentals.md), the notes are based on the textbook[^2]

# The Bisection Method

Given initial interval $[a, b]$ such that $f(a) f(b)<0$ 

```Python
def bisection(f, a, b, tol=1e-6, maxiter=100):
    for i in range(maxiter):
        c = (a + b) / 2
        if f(c) == 0 or abs(b - a) < tol:
            return c
        if f(a) * f(c) < 0:
            b = c
        else:
            a = c
    return c
```

the approximate error is $\varepsilon_{error} = \frac{b-a}{2^{n+1}}$ ($n=0$ is the initial error)

***

# Fixed-point Iteration
- **fixed point**: the real number $r$ is a fixed point of the function $g$ if $g(r)=r$ 

```Python
def fixed_point_iteration(f, x0, tol=1e-6, maxiter=100):
    x = x0
    for i in range(maxiter):
        x1 = f(x)
        if abs(x1 - x) < tol:
            return x1
        x = x1
    return x
```

> The FPI function can be different.

For equation $x^{3}+ x -1 =0$, there exists different FPI functions, along with different convergence speed.
<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/2oAmFkg.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">different FPI functions, different convergence speed.
    </div>
</center>

- **Linear convergence**: Let $\varepsilon_{i} = \frac{\varepsilon_{i-1}}{2}$ denote the error at step $i$ of an iterative method. If
$$
\lim_{i \rightarrow \infty } \frac{\varepsilon_{i+1}}{\varepsilon_{i}} = S < 1
$$
the method is said to obey linear convergence with rate S

> Assume that $g$ is continuously differentiable, that $g(r) = r$ , and that $S = |g'(r)| < 1$ . Then Fixed-Point Iteration converges linearly with rate S to the fixed point $r$ for initial guesses sufficiently close to $r$ .

## Convergent Types

- **Locally convergent**: An iterative method is called locally convergent to $r$ if the method converges to r for initial guesses sufficiently close to $r$ 

- **Globally convergent**: For any $x \in I$, the method converges to $r$ for any initial guesses of $x$.

> The Bisection Method is guaranteed to converge linearly. Fixed-Point Iteration is only locally convergent, and when it converges it is linearly convergent. Both methods require one function evaluation per step. The bisection cuts uncertainty by $1/2$ for each step, compared with approximately $S = |g' (r)|$ for FPI. Therefore, <u>Fixed-Point Iteration may be faster or slower than bisection, depending on whether S is smaller or larger than 1/2</u>. Newton’s Method, a particularly refined version of FPI, where $S$ is designed to be zero

## Cases when FPI works

1. Function $f(x)$ is continue in $[a,b]$
	1. $\forall x \in [a,b], a \leqslant f(x) \leqslant b$ , then there exists fixed-point.
	2. Based on 1, if $\exists L \in(0,1), s.t. |f(x_{1})-f(x_{2}) \leqslant L(x_{1}-x_{2}) |~ \forall x_{1}, x_{2} \in [a,b]$
		+ $f(x)$ exists **only** one **global** fixed-point.
		+ Error satisfies $|x_{k}-x^{*}|\leq\frac{L}{1-L}|x_{k}-x_{k-1}|\leq\cdots\leq\frac{L^{k}}{1-L}|x_{1}-x_{0}|$

**Local Convergence  Theorem**: Let $x^{*}$ be the fixed-point, if $f'(x)$ is continue in $[a,b]$ and  $|f'(x^{*})|<1$ .Then the FPI is **locally convergent**.

****

# Limits of Accuracy

- **Forward and backward error**: Assume that $f$ is a function and that $r$ is a root, meaning that it satisfies $f (r) = 0$. Assume that $x_a$ is an approximation to $r$. For the root-finding problem, the <u>backward error</u> of the approximation $x_a$ is $| f (x_{a})-0| = |f(x_{a})|$ and the <u>forward error</u> is $|r - x_a|$

## Sensitivity of Root-finding
A problem is called <u>sensitive</u> if small errors in the input, in this case the equation to be solved, lead to large errors in the output, or solution.

Assume that the problem is to find a root $r$ of $f (x) = 0$, but that a small change $g(x)$ is made to the input, where is small. Let $r$ be the corresponding change in the root:
$$
f(r+ \Delta r) + \varepsilon g(r+ \Delta r) =0
$$

$$
\Delta r \approx -\varepsilon g\frac{r}{f'(r)+\varepsilon g'(r)} \approx -\varepsilon \frac{g(r)}{f'(r)}
$$

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/qgTm6pD.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Example 1.9
    </div>
</center>

- **Error magnification factor** : relative forward error / relative backward error.
$$\left|\frac{\Delta r / r}{\epsilon g(r) / g(r)}\right|=\left|\frac{-\epsilon g(r) /\left(r f^{\prime}(r)\right)}{\epsilon}\right|=\frac{|g(r)|}{\left|r f^{\prime}(r)\right|}$$

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/7jTLFOo.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Examle of Wilkinson polynomial
    </div>
</center>

The **condition number** of a problem is defined to be the maximum error magnification over all input changes, or at least all changes of a prescribed type. A problem with <u>high condition number</u> is called **ill-conditioned**, and a problem with a <u>condition number near 1</u> is called **well-conditioned**

****

# Newton's Method

> [!tip] Newton's Method
> $$\begin{aligned}
x_0 & =\text { initial guess } \\
x_{i+1} & =x_i-\frac{f\left(x_i\right)}{f^{\prime}\left(x_i\right)} \text { for } i=0,1,2, \ldots
\end{aligned}$$

If $f(x)\in C^{2}(\mathbb{R})$ is a **monotonically increasing convex** function, and $f(r)=0$ , then $r$ is the only root and for $\forall x_{0} \in \mathbb{R}$ , Newton’s method will converges to $r$ .

## Convergence Speed

Based on the [Rate of convergence - Wikipedia](https://en.wikipedia.org/wiki/Rate_of_convergence), the order of convergence[^1] and the rate of convergence of a convergent sequence are quantities that represent how quickly the sequence approaches its limit. A sequence $(x_{n})$ that converges to $x^{*}$ is said to have order of convergence $q \ge 1$ and the rate of convergence $\mu$ if :

$$
\lim _{n \rightarrow \infty} \frac{\left|x_{n+1}-x^*\right|}{\left|x_n-x^*\right|^q}=\mu
$$

The rate of convergence $\mu$ is also called the *asymptotic error constant.*
- $q=1$ is called _linear convergence_ if, and the sequence is said to _converge Q-linearly to $L$.
- $q=2$ is called _quadratic convergence._
- $q=3$ is called _cubic convergence._
- Etc.

Empirically, if the function $g \in C^{p}(a, b)$ is convergent to $x^*$ ，and its $p$ 
differential $g^{(p)}(x^{*})\neq 0$ , while $g^{(q)}(x^{*})= 0, q=1,2,\cdots,p-1$ then it's said that $x_{i+1}=g(x_{i})$ has $p$ order convergence speed, and it's **local convergent**.  For FPI, we have:
$$
\lim_{i \rightarrow \infty } \frac{e_{i+1}}{e_{i}^{p}}=\frac{|g^{(p)}(x^{*})|}{	p!}
$$


Let $f \in C^{2}[a,b]$ and $f(r)=0$. If $f^{\prime}(r) \neq 0$, then Newton's Method is **locally and quadratically convergent** to $r$. The error $e_i$ at step $i$ satisfies
$$
\lim _{i \rightarrow \infty} \frac{e_{i+1}}{e_i^2}=\frac{f^{\prime \prime}(r)}{2 f^{\prime}(r)}
$$

If $f \in C^{(m+1)}[a,b]$ , which contains a root $r$ of multiplicity $m>1$ , then Newton's Method's error $e_{i}$ satisfies
$$
\lim_{t \rightarrow \infty} \frac{e_{i+1}}{e_{i}}= \frac{m-1}{m} 
$$
converges **locally and linearly** to $r$ .

If $f \in C^{(m+1)}[a,b]$ , which contains a root $r$ of multiplicity $m>1$ , and $f^{(m+1)}(r)\neq 0$ , then **Modified Newton's Method** 
$$
x_{i+1}=x_{i}-\frac{mf(x_{i})}{	f^{'}(x_{i})}
$$

error $e_{i}$ satisfies
$$
\lim_{t \rightarrow \infty} \frac{e_{i+1}}{e_{i}^{2}}= \frac{1}{m(m+1)} \left| \frac{f^{(m+1)}(r)}{f^{(m)}(r)} \right| 
$$

converges **locally and quadratically** to $r$ .

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/TeJqnYr.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Examples of error approximation
    </div>
</center>


## Discuss When Newton's Method Fails
1. The initial guess is not in locally convergence domain
2. Its differential equals to 0
3. The root function contains negative result.

> [!example] An examle when the Newton's method fails
> find the root of $f(x) = 4 x^{4}- 6 x^{2}- \frac{11}{4}$ with $x_{0}=1/2$ 

# Methods without Derivatives

## Secant Method and Variants

> [!tip] Secant method
> $$\begin{aligned}x_0, x_1 & =\text { initial guesses } \\x_{i+1} & =x_i-\frac{f\left(x_i\right)\left(x_i-x_{i-1}\right)}{f\left(x_i\right)-f\left(x_{i-1}\right)} \text { for } i=1,2,3, \ldots\end{aligned}$$

the approximate error relationship:

$$e_{i+1} \approx | \frac{f^{''}(r)}{2 f^{'}(r)} | e_{i}e_{i+1} \quad e_{i+1} \approx | \frac{f^{''}(r)}{2 f^{'}(r)} |^{\alpha-1} e_{i}^{\alpha}~\alpha=\frac{1+\sqrt{5}}{	2}$$

### Three Generalizations of the Secant  Method

> [!tip] False Position
> Given interval $[a, b]$ such that $f(a) f(b)<0$ for $i=1,2,3, \ldots$
$\begin{aligned} & c=\frac{b f(a)-a f(b)}{f(a)-f(b)} \\ & \text { if } f(c)=0, \text { stop, end } \\ & \text { if } f(a) f(c)<0 \\ & \quad b=c\end{aligned}$
else $a=c$
end $\quad$ end


# Newton’s method in solving nonlinear equation

For equations like 
$$
\begin{aligned}
&\left\{\begin{matrix}f_{1}\left(x_{1},x_{2},...,x_{n}\right)=0\\ f_{2}\left(x_{1},x_{2},...,x_{n}\right)=0\\ \vdots\\ f_{n}\left(x_{1},x_{2},...,x_{n}\right)=0\end{matrix}\right.
\end{aligned}
$$

Let's Define $F=\left[f_{1},f_{2},\cdots,f_{3}\right]^{T}$ , and **Jacobi Matrix:** $$J\left(X\right)=F^{\prime}\left(X\right)=\begin{pmatrix}\frac{\partial f_{1}\left(X\right)}{\partial x_{1}}&\cdots&\frac{\partial f_{1}\left(X\right)}{\partial x_{n}}\\ \vdots&\ddots&\vdots\\ \frac{\partial f_{n}\left(X\right)}{\partial x_{1}}&\cdots&\frac{\partial f_{n}\left(X\right)}{\partial x_{n}}\end{pmatrix}$$

Then here comes the Newton’s method:
$$
X_{i+1}=X_i-J^{-1}\left(X_i\right)F\left(X_i\right)
$$

> [!example] Example of Newton's Method in solving nonlinear equations
> - 例: 求解方程组
> $$
> \left\{\begin{array}{c}
> f_1(u, v)=v-u^3=0 \\
> f_2(u, v)=u^2+v^2-1=0
> \end{array}\right.
> $$
> Jacobi i矩阵
> $$
> J(u, v)=\left[\begin{array}{cc}
> -3 u^2 & 1 \\
> 2 u & 2 v
> \end{array}\right]
> $$
> 牛顿迭代式 $\left(X=[u, v]^T, F=\left[f_1, f_2\right]\right)$ :
> $J\left(u_i, v_i\right)\left(X_{i+1}-X_i\right)=-F\left(X_i\right)$
> 取 $X_0=[1,2]$
> 第一步迭代计算为
> $$
> \left[\begin{array}{cc}
> -3 & 1 \\
> 2 & 4
> \end{array}\right]\left[\begin{array}{l}
> u_{i+1}-u_i \\
> v_{i+1}-v_i
> \end{array}\right]=-\left[\begin{array}{l}
> 1 \\
> 4
> \end{array}\right] \Rightarrow\left[\begin{array}{l}
> u_{i+1} \\
> v_{i+1}
> \end{array}\right]=\left[\begin{array}{l}
> 1 \\
> 1
> \end{array}\right]
> $$




[^1]: Wikipedia Contributors (2023) Rate of convergence. In: Wikipedia. https://en.wikipedia.org/wiki/Rate_of_convergence. Accessed 13 Mar 2023
[^2]: Sauer, Timothy. _Numerical analysis_. Addison-Wesley Publishing Company, 2011.
‌