---
date: 2023-03-19 17:16:05
categories: Numerical_Analysis
destination: 
excerpt: 
katex: true
obsidianUIMode: source
rating: ⭐⭐
draft: false
tags:  
- 作业
- Numerical_Analysis
title: 数值分析-第二次作业-nonlinear equation
share: true
---

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/d6CBfQm.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">第二次作业截图
    </div>
</center>


> [!question] 1
> - 对以下几种计算 $\sqrt[4]{2}$ 的算法的收敛速率从最快到最慢排序, 并说明理由
> 1. 使用二分法求解 $f(x)=x^4-2=0$
> 2. 使用割线法求解 $f(x)=x^4-2=0$
> 3. 通过迭代求 $g(x)=\frac{x}{2}+\frac{1}{x^3}$ 的不动点
> 4. 通过迭代求 $g(x)=\frac{5 x}{6}+\frac{1}{3 x^3}$ 的不动点
> 5. 使用牛顿法求解 $f(x)=x^4-2=0$

为了方便, 初始点选为 $x_{0}=1;\text{初始区间选为}[1,2]$

1. 对于二分法, 收敛速率为 $\frac{1}{2}$
2. 对于割线法, 如果初始点足够接近根的话收敛速率为 $1.618$. 
3. 直接计算 $g'(x)=\frac{1}{2} - \frac{3}{x^{4}}$ 代入 $x=\sqrt[4]{2}$ 得到结果为1, 而在其根的邻域中其绝对值趋于零. 可以认为在附近难以收敛到根, 故收敛速率可以近似为0. 实际中使用程序迭代也确实给出了该结果.

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/jM6zN1l.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">第3小题的程序迭代100次后结果
    </div>
</center>

4. $g'(r) =\frac{1}{3}$ 在根附近进行迭代其收敛速率接近 $\frac{1}{3}$.实际程序迭代10次后就达到5位小数精度.

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/JaZ193q.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">第4小问程序迭代结果
    </div>
</center>

5. 对于标准牛顿法, 其收敛速率为二阶, 实际程序中第5次便达到上述精度.

综上, 以上五种迭代速率由快到慢的排序为: $5 \rightarrow 2 \rightarrow 1 \rightarrow 4 \rightarrow 3$



> [!question] 2
> 如果 $f(x) \in C^m(\mathbb{R})$, 并满足 $f(r)=0, f^{\prime}(r) \neq 0,$ $f^{(m)}(r) \neq 0(m>2)$ $f^{\prime \prime}(r)=\cdots=f^{(m-1)}(r)=0$,  ，试分析标准牛顿迭代法的收敛性及误差。

$$
f(r)=f(x_{i})+(r-x_{i})f'(x_{i})+(r-x_{i})^{m} \frac{f^{(m)}(\xi)}{m!}
$$

化简得到:
$$
x_{i}-\frac{f(x_{i})}{f'(x_{i})}-r =(r-x_{i})^{m} \frac{f^{(m)}(\xi)}{m! f'(x_{i})}
$$
即：
$$
 \frac{e_{i+1}}{e_{i}^{m}} = \frac{f^{(m)}(\xi)}{m! f'(x_{i})} 
$$

取极限后发现此情况下标准牛顿迭代法的收敛速度为 $m$ 阶

> [!question] 3
> 找出函数 $g(x)=x^2-\frac{3}{2} x+\frac{3}{2}$ 的每个不动点, 并确定不动点迭代是否局部收敛.
$$ x^2-\frac{3}{2} x+\frac{3}{2} = x $$

移项化简得：

$$ x^2-\frac{5}{2} x+\frac{3}{2} = 0 $$

解这个二次方程，得到两个根：

$$ x_1 = 1, x_2 = \frac{3}{2} $$

因此，$g(x)$ 的不动点为 $1$ 和 $\frac{3}{2}$。

计算出 $g(x)$ 在不动点附近的导数，以便通过导数的大小来判断局部收敛性。在这里， $g(x)$ 的导数：

$$ g'(x) = 2x - \frac{3}{2} $$

对于不动点 $x_1 = 1$，有 $g'(1) = -\frac{1}{2}$，因此 $x_1$ 的不动点迭代是局部收敛的。

对于不动点 $x_2 = \frac{3}{2}$，有 $g'(\frac{3}{2}) = \frac{3}{2}$，因此 $x_2$ 的不动点迭代不是局部收敛的。

综上所述，$g(x)$ 的不动点为 $1$ 和 $\frac{3}{2}$，其中 $1$ 的迭代是局部收敛的，$\frac{3}{2}$ 的迭代不是局部收敛的。


> [!question] 4
> 令 $f(x)=x^4-7 x^3+18 x^2-20 x+8$, 牛顿法是否会二次收敛到根 $r=2$ ? 确定 $\lim _{i \rightarrow \infty} e_{i+1} / e_i$

$f(x)=(x-1)(x-2)^3$ 有一个三重根。

先求出 $f(x)$ 的导数：


$$
\begin{align}  
f'(x) &= 4x^3 - 21x^2 + 36x - 20 \\
f''(x) &= 12x^2 - 42x + 36 \\
f'''(x) &= 24x - 42 \\
\end{align}
$$



$f'(r)=f''(r)=0 \neq f'''(r)$ 因此

由牛顿法可知，迭代公式为：

$$x_{i+1} = x_i - \frac{f(x_i)}{f'(x_i)} =\frac{3x_{i}^{2}-2x_{i}-2}{4x_{i}-5}$$

由书上定理可知 此时牛顿迭代法会线性收敛，且 $\lim _{i \rightarrow \infty} \frac{e_{i+1}}{e_i}=S\approx2/3$



> [!question] 5
> 函数 $f(x)=x^3-4 x$ 。
> (a) 假设使用牛顿法求根 $r=2$, 第4步后误差是 $e_4=10^{-6}$, 估计 $e_5$ ；
> (b) 假设使用牛顿法求根 $r=0$, 第 4 步后误差是 $e_4=10^{-6}$, 估计 $e_5$ ；

$f(x)=x(x+2)(x-2)$, $f'(x)=3x^2-4$, $f''(x)=6x$

没有重根，标准的牛顿迭代法收敛速率为二阶，迭代误差满足：
$$
\lim _{i \rightarrow \infty} \frac{e_{i+1}}{e_i^2}=\left|\frac{f^{\prime \prime}(r)}{2 f^{\prime}(r)}\right|
$$

(a.) 对于 $r=2$, $\left|\frac{f^{\prime \prime}(r)}{2 f^{\prime}(r)}\right|=\frac{3}{4}$ $e_{5}=\frac{3e_{4}^2}{4} \approx 10^{-13}$ 
(b.) 对于 $r=0$, $\left|\frac{f^{\prime \prime}(r)}{2 f^{\prime}(r)}\right|=0$

当 $x_4$ 时，有

$f(x_4)=x_4^{3-4x_{4}\quad}f'(x_4)=3x_4^2-4$


由于 $r=0$ 是这个函数的根，所以有 $f(r)=0$，因此 $f(x_4)$ 是 $r$ 的近似值，即

$$|f(x_4)|=|r-x_4|^3\approx e_4^3$$

接下来估计 $|r-x_5|$：

$$|r-x_5|=|x_5-x_4|=|\frac{f(x_4)}{f'(x_4)}|$$

将牛顿法的迭代公式代入，得到

$\begin{aligned} |r-x_5|&=|\frac{f(x_4)}{f'(x_4)}| \\ &=|\frac{x_4^3-4x_4}{3x_4^2-4}| \\ &=|\frac{x_4(x_4^2-4)}{3x_4^2-4}|\\ &=|\frac{x_4}{3x_4^2-4}||x_4^2-4|\\ &=|\frac{1}{3x_4+\frac{4}{x_4}}||x_4^2-4| \end{aligned}$

由于 $|x_4^2-4|=|x_4-2||x_4+2|$，因此

$$|r-x_5|=|\frac{1}{3x_4+\frac{4}{x_4}}||x_4-2||x_4+2|$$

考虑估计分母 $\frac{1}{3x_4+\frac{4}{x_4}}$ 的大小，由于 $x_4$ 接近于 $r=0$，因此 $\frac{4}{x_4}$ 的绝对值将非常大，可以忽略 3，得到

$$|r-x_5|\approx\frac{|x_4|}{4}|x_4-2||x_4+2|$$

$e_{4}=|x_{4}-r|=x_{4}$

$$e_{5}=|r-x_{5}|=\frac{10^{-18}}{4}-10^{-6} \approx 10^{-6}$$


> [!question] 6
> 取 $\left(u_0, v_0\right)=(1,1)$, 使用牛顿迭代进行两步迭代求解
$$\left\{\begin{array}{c}u^2+v^2=1 \\(u-1)^2+v^2=1\end{array}\right.$$

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/UhAUdAp.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">计算方法大同小异,我就直接贴程序了
    </div>
</center>

首先得到各个矩阵：
$$X_{0}=\left( \begin{array}{c} x \\ y \end{array} \right) \quad J^{-1}(u,v)=\left( \begin{array}{c} \frac{1}{2} & -\frac{1}{2} \\ \frac{1-u}{2v} & \frac{u}{2v} \end{array} \right) \quad F(u,v)=\left( \begin{array}{c} u^{2}+v^{2}-1 \\ u^{2}-2u +v^{2}\end{array} \right)$$

代入计算公式：
$$X_{i+1}=X_{i}-J^{-1}(X_{i}) \cdot F(X_{i})$$
得到：
$$X_{1}=\left( \begin{array}{c} \frac{3}{2} \\ 1 \end{array} \right) \quad X_{2}=\left( \begin{array}{c} \frac{3}{2} \\ \frac{7}{8} \end{array} \right)
$$


> [!question] 7
> 编程：使用二分法、不动点迭代、牛顿法、试位法求 $\mathrm{e}^x+x=7$ 的根, 结果精确到小数点后8位 并对各方法的收敛性进行比较。

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/GwsE4Pe.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">程序结果
    </div>
</center>

结果展示牛顿法和二分法在给定的初始值、初始区间内得到了收敛的结果，其他的在超过100次迭代次数后未达到预期的收敛值而退出。

由此可知，在日常中，牛顿法不仅算法复杂度低且收敛速度快（虽然有更快的三阶方法，但是复杂度上去了）。而二分法在给定初始区间足够优秀，也能有不错的收敛结果。而试位法、不动点迭代等方法的性能高度依赖于初始点、初始区间的选取，在本例题中没有体现相应的优异性。
