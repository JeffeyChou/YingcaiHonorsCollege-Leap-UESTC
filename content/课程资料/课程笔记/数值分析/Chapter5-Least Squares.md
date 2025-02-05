---
draft: false
date: 2023-06-14 19:44:27
categories: Numerical_Analysis
destination: 
excerpt: 
katex: true
obsidianUIMode: source
rating: ⭐⭐
tags:
  - Numerical_Analysis
  - Least_squares
  - QRDcomposition
title: "Chapter5-Least Squares"
updated: 2023-06-14 20:30:24
---

In this Chapter we will introduction two main algorithms for fitting problem, the least squares and QR decomposition.

# The Fitting Problem

Like what we talk in [Chapter4-Interpolation](term/Numerical_Analysis/Chapter4-Interpolation.md), we always want to find a proper function $P (X)$ that go through all the given points. But it may be too costly sometimes. So does there any way to fit the given points? That's where least squares works. 

Suppose there are $n+1$ points
$$
(x_{0}, y_{0}),(x_{1}, y_{1}),\cdots,(x_{n}, y_{n})
$$
To fit with given function:
$$
\varphi(x)=a_{0}\varphi_{0}\left(x\right)+\cdots+a_{m}\varphi_{m}\left(x\right),
$$
Given the weight matrix, label matrix and coefficient matrix:
$$X=\begin{bmatrix}a_{0} \\ a_{1} \\ \vdots\\ a_{m}\end{bmatrix}, b=\begin{bmatrix}y_0 \\ y_1 \\ \vdots \\ y_{n}\end{bmatrix}, A=\begin{bmatrix}1 & \varphi_{0}(x_{0}) & \cdots & \varphi_{m}(x_{0}) \\ 1 & \varphi_{0}(x_{1}) & \cdots & \varphi_{m}(x_{1}) \\ \vdots & \ddots & \ddots & \vdots \\ 1 & \varphi_{0}(x_{n}) & \cdots & \varphi_{m}(x_{n})	\end{bmatrix}$$

Our goal is to let the loss function 
$$
L(X)=\|r\|_{2}=\|AX-b\|_{2}=\frac{1}{	2}(AX-b)^{\top}(AX-b)
$$
To be least.
Namely:

$$
\frac{\partial L(X)}{	\partial X}=A^{\top}(AX-b)=0
$$

Rewrite it:
$$
X^{*}=(A^{\top}A)^{-1}A^{\top}b
$$

The problem turn out to be solve the equation above.

# Least Squares

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/HPf9LMO.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">text title here...
    </div>
</center>

> [!example]+ example of least squares
<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/tuScnnI.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">example of least squares
    </div>
</center>

## QR Decomposition

If we Let
$$
A=QR
$$
Then:
$$
X^{*}=R^{-1} Q^{\top} b
$$

## Connection Between Least Squares and QR Decomposition

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/fzkrCIy.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Least squares and QR decompostion
    </div>
</center>
