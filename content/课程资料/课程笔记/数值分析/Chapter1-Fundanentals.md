---
cover: 
date: 2023-03-01 10:20:35
destination: 10-blog/source/_posts/数值分析
excerpt: 数值分析的第一次课-预备知识
katex: true
obsidianUIMode: source
rating: ⭐
status: complete
tags:
  - Numerical_Analysis
title: "Chapter1-Fundanentals"
share: false
---

# Error

measured value: $x$ , exact value: $x_{0}$

## Error's operation

$$
\begin{align}
e(x \pm y) = e(x) \pm e(y) \tag{1.1}\\
x(xy) = ye(x) + x_{0} e(y) \tag{1.2}\\
e\left(\frac{x}{y}\right)= \frac{xe(y) + ye(x)}{yy_{0}} \tag{1.3} \\
\end{align}
$$
For *Eq.1.1*, when two approximate value subtracted, its **relative error** will be increased.

For *Eq.1.3*, when the value is divided by decimals, its **absolute error** will be increased.

# Significant digits

- Definition: 当 $x$ 的绝对误差限为某一位的半个单位，则这一位到第一个非零位的位数称为 $x$ 的有效位数。


# Binary numbers

## Decimal to binary 
- Integer part:
    - dividing by 2 successively and recording the remainders.
$$
\begin{aligned}
53 \div 2 & =26 \ \mathrm{R}\  1 \\
26 \div 2 & =13 \ \mathrm{R}\  0 \\
13 \div 2 & =6 \ \mathrm{R}\  1 \\
6 \div 2 & =3 \ \mathrm{R}\  0 \\
3 \div 2 & =1 \ \mathrm{R}\  1 \\
1 \div 2 & =0 \ \mathrm{R}\  1
\end{aligned}
$$
therefore, $(53)_{10} = (110101.)_{2}$ 

- Fractional part
    - Multiply by 2 successively and record the intger parts.
$$
\begin{aligned}
& .7 \times 2=.4+1 \\
& .4 \times 2=.8+0 \\
& .8 \times 2=.6+1 \\
& .6 \times 2=.2+1 \\
& .2 \times 2=.4+0 \\
& .4 \times 2=.8+0
\end{aligned}
$$
Therefore, $(0.7)_{10} = (.1 \overline{0110})_{2}$ 

# Floating point for real numbers 

| precision   | sign | exponent | mantissa |
| :-----------: | :----: | :--------: | :--------: |
| single      | 1    | 8        | 23       |
| double      | 1    | 11       | 52       |
| long double | 1    | 15       | 64       | 

the form of a **normalized IEEE floating point number** :

$$
\pm 1.b_{1}b_{2}\ldots b_{m} \times 2^{p} \tag{3.1}
$$

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="https://search.pstatic.net/common?src=https://i.imgur.com/ZrsDDxr.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">double precision number of 1
    </div>
</center>

- **Machine epsilon**: $\epsilon_{match}$ is the distance between 1 and the smallest floating point number greater than 1. $\epsilon_{mathc}= 2^{-52}$.
- **Truncate number method**：
    - Chopping
        - if the 53nd bit is 0, than truncate the number at 53nd bit.
    - rounding
        - add 1 to bit 52 if bit 53 is 1, and add nothing to bit 52 if bit 53 is 0.
    - Rounding to Nearest Rule
        - the 53 bit is 1, and is followed by other nonzero bits. either round up or round down to make 52 bit equal to 0.
- **Bit 1** : sign number 
    - 0: +
    - 1: -
- **Bit 2**: exponent sign number
    - 1: the exponent number is greater than 1
    - 0: the exponent number is less than 1

HW: [[数值分析-第一次作业-误差与计数]]