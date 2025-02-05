---
draft: false
date: '2023-06-14 20:31:26'
categories: Numerical_Analysis 
destination: 
excerpt: 
katex: true
obsidianUIMode: source
rating: ⭐⭐
tags:  
- Numerical_Analysis 
- Differentiation
- Integration
title: "Chapter6-Numerical Differentiation"
share: true
updated: 2023-06-15 16:08:24
---

The main problem of computational calculus is to compute derivatives and integrals of functions. There are two directions that we can take for such problems, numerical computing and symbolic computing. Some problem turn out to be too difficult to apply symbolic computation, and sometimes we just want the numerical results. So that's what this chapter talks about.

# Numerical Differentiation

## Two-point Difference Formula

$$
\begin{gathered}
forward:\hat{f}^{\prime}\left(x\right)=\frac{f\left(x+h\right)-f\left(x\right)}{h}=f^{\prime}\left(x\right)+\frac{h}{2}f^{\prime\prime}\left(\xi\right),\xi\in\left(x,x+h\right) \\
backward:\hat{f}^{\prime}\left(x\right)=\frac{f\left(x\right)-f\left(x-h\right)}{h}=f^{\prime}\left(x\right)-\frac{h}{2}f^{\prime\prime}\left(\xi\right),\xi\in\left(x-h,x\right) \\
centre:\hat{f}^{\prime}\left(x\right)=\frac{f\left(x+h\right)-f\left(x-h\right)}{2h}=f^{\prime}\left(x\right)+\frac{h^{2}}{6}f^{\prime\prime\prime}\left(\xi\right),\xi\in\left(x-h,x+h\right) \\
\end{gathered}
$$

We can apply these formulas continually:
$$
\begin{gathered}
\hat{f}''(x)=\frac{{f'\left(x+\frac{h}{2}\right)-f'\left(x-\frac{h}{2}\right)}}{h} \\
=\frac{\frac{f\left(x+h\right)-f\left(x\right)}{h}-\frac{f\left(x\right)-f\left(x-h\right)}{h}}h \\
=\frac{f\left(x+h\right)-2f\left(x\right)+f\left(x-h\right)}{h^{2}} \\
= f''(x)+ \frac{h^{2}}{	12}f''''(\xi), \xi \in (x-h, x+h)
\end{gathered}
$$

## Interpolation Method

The key idea of Interpolation method for numerical differentiation is using interpolation function $\hat{f}(x)$ which is easier to differentiate to fit the original function $f(x)$ .

### Richardson Extrapolation

Assume that we are presented with an order $n$ formula $F (h)$ for approximating a given quantity $\hat{F}(h)$ . The order means that
$$
F(h)=\hat{F}(h)+Ch^n+O(h^{n+1})
$$
Then 
$$
\tilde{F}\left(h\right)=\frac{2^{n}\hat{F}\left(\frac{h}{2}\right)-\hat{F}\left(h\right)}{2^{n}-1}
$$
Is $n+1$ order approximating at least.

> [!example]+ example of Richardson Extrapolation

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/ffgSZXh.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Richardson extrapolation 
    </div>
</center>

## Undetermined Coefficients Method

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/nkDMCsZ.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Undetermined coefficients method
    </div>
</center>

# Numerical Integration

## Interpolation Method

The numerical calculation of definite integrals relies on many of the same tools we have already seen. In Chapters 3 and 4, methods were developed for finding function approximation to a set of data points, using interpolation and least squares modeling. We will discuss methods for numerical integration, or quadrature, based on both of These ideas.

If we select data points at the same interval, this numerical interpolation method is called Newton-Cote's integration.

### Newton-Cote's Integration

#### Trapezoid Rule

<center>
    <img style="border-radius: 0.3125 em;
    Box-shadow: 0 2 px 4 px 0 rgba (34,36,38,. 12), 0 2 px 10 px 0 rgba (34,36,38,. 08);"
    src=" https://search.pstatic.net/common?src=https://i.imgur.com/dk3mjLG.png">
    <br>
    <div style="color: orange; border-bottom: 1 px solid #d9d9d9 ;
    Display: inline-block;
    color: #999 ;
    Padding: 2 px;">example of trapezoid rule
    </div>
</center>

Interpolation error:
$$
E(x)=f(x)-P_{1}(x)=\frac{1}{2}f^{\prime\prime}\bigl(c(x)\bigr)\bigl(x-a)(x-b),c(x)\in\left(a,b\right)
$$

Integration error:
$$
E\left(f\right)=I\left(f\right)-I\left(P_{1}\right)=\int_{a}^{b}E\left(x\right)d x=-\frac{h^{3}}{12}f^{\prime\prime}\left(\xi\right),\xi\in\left(a,b\right)
$$

#### Simpson's Law
For data points: $x_{0}=a, x_{1}=\frac{a+b}{2},x_{2}=b$ , We use [Lagrange Interpolation](term/Numerical_Analysis/Chapter4-Interpolation.md#Lagrange%20Interpolation) to fit the original function:

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/ovqDuSH.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Simpson Law
    </div>
</center>

The final formula:
$$
I(f) \approx \frac{h}{6}\left[f\left(a\right)+4f\left(\frac{a+b}{2}\right)+f\left(b\right)\right],h=b-a
$$

Interpolation error:
$$
E(x)=f(x)-P_2(x)=\frac{1}{3!}f''\big(c(x)\big)(x-a)\left(x-\frac{a+b}{2}\right)(x-b),c(x)\in(a,b).
$$

Integration error:
$$
E\left(f\right)=-\frac{h^5}{2880}f^{\prime\prime\prime\prime\prime}\left(\xi\right),\xi\in\left(a,b\right)
$$

#### Algebra Accuracy
If a certain quadrature formula is accurate for polynomials of degree $n$ , but not for polynomials of degree $n$ , then it is said that the quadrature formula has the accuracy of degree $n$ .

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/vNVFANQ.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">example of algebra accuracy
    </div>
</center>

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/R95YFPm.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Accuracy of Newton-Cote's integration 
    </div>
</center>

### Opened Newton-Cote's Integration

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/DeYvyox.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">text title here...
    </div>
</center>

#### Composite Newton-Cote's Interpolation

The Trapezoid and Simpson’s Rules are limited to operating on a single interval. Of course, since definite integrals are additive over subintervals, we can evaluate an integral by dividing the interval up into several subintervals, applying the rule separately on each one, and then totaling up. This strategy is called composite numerical integration.

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/L1yUY9J.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">text title here...
    </div>
</center>

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/MIS7hk6.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">example of composite trapezoid rule and Simpson law
    </div>
</center>

## Romberg Integration

The main idea of Romberg integration is using [Richardson Extrapolation](term/Numerical_Analysis/Chapter6-Numerical%20Differentiation%20and%20Integration.md#Richardson%20Extrapolation) to increase the accuracy of the integration.

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/oJh7NdU.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">text title here...
    </div>
</center>

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/uO0MT3H.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">pseudocode of Romberg
    </div>
</center>

## Undetermined Coefficients Method

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/0mlY35R.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Undetermined Coefficients Method
    </div>
</center>

## Gauss Integration

Are the Newton–Cotes formulas optimal for their degree of precision, or can more powerful formulas be developed? In particular, if the requirement that evaluation points be evenly spaced is relaxed, are there better methods?

We pick out the most famous one to discuss in this section. Gaussian Quadrature has degree of precision $2 n + 1$ when $n + 1$ points are used, double that of Newton–Cotes.

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/chz402E.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">text title here...
    </div>
</center>

Integration error:
$$
\int_{a}^{b}f\left(x\right)d x-\sum_{i=0}^{n}A_{i}f\left(x_{i}\right)=\frac{f^{(2n+2)}\left(\xi\right)}{\left(2n+2\right)!}\int_{a}^{b}\prod_{i=0}^{n}\left(x-x_{i}\right)^{2}d x,\xi\in\left(a,b\right)
$$
There are mainly three ways to find $Q(x)$ :
1. Let $Q_{n+1}(x)=x^{n+1}+a_{n}x^{n}+ \cdots + a_{1}x +a_{0}$ , $s.t. \int_{a}^{b}Q_{n+1}(x)x^{k} dx =0,~ k=0, 1,\cdots,n$
2. Given linearly independent functions $\{\varphi_{0}(x), \varphi_{1}(x), \cdots ,\varphi_{n}(x) \}$ in interval $[a,b]$ , build orthogonality function using [Gram-Schmidt](term/ref/Gram-Schmidt.md) 
3. Legendre polynomials $L_{n}(x)=\frac{1}{2^{n} \cdot n! } \frac{d^{n}}{dx^{n}}(x^{2}-1)^{n}$

We usually use Legendre polynomials to compute numerical integration.

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/UeFOEvJ.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Gaussian Quadrature coefficients
    </div>
</center>

> [!example]+ approximate integrals on a general interval
<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/LV6REjy.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">approximate integrals on a general interval
    </div>
</center>
