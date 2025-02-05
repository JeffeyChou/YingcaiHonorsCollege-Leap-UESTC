# 数理基科生存手册@Yingcai Honors College, UESTC




## Getting Started

This is a guide for Yingcai Honors College, UESTC students who are interested in mathematics and want to pursue a career in mathematics. It provides a general overview of the mathematics program at Yingcai Honors College, UESTC, and includes information about the courses, exams, and other requirements for the program.


This repo requires a basic understanding of Markdown and Git. If you are not familiar with these tools, please refer to the following resources:

- [Markdown Guide](https://www.markdownguide.org/)
- [Git Guide](https://git-scm.com/book/en/v2)

Requirements for contributing:

- You must have a GitHub account.
- You must have Git installed on your computer.
- You must have a basic understanding of Markdown.

## Contributing

It is recommented that you use the following software/tools to contribute to this project:

- [Visual Studio Code](https://code.visualstudio.com/)
- [Obsidian](https://obsidian.md/)
- [Quartz](https://github.com/jackyzha0/quartz)


To contribute to this project, follow these steps:

1. Fork this repository.

```
git clone https://github.com/JeffeyChou/YingcaiHonorsCollege-Leap-UESTC.git
```

2. Create a new branch for your changes.

```
git checkout -b my-branch
```

3. Make your changes in Obsidian and preview changes locally.

3.1 upload notes in Obsidian.

- Open this vault in Obsidian and make your changes. This should create your own profile in the `.obsidian` folder.

- Make your changes in Obsidian on the `content/` folder.

>[!Note]+ Some tips for writing Latex math formulas and rendering using Katex.
>QuartZ uses Katex to render latex formulas in markdown. In order to correctly render the formula without raising errors, it is preferred to use the following way when writing your latex formulas.
>1. Inline mode
>```
>$A_1=R_0 \cdot Q_0$
>```
>Make sure that the two `$` are in the same line.
>1. Display mode
>```
>$$\begin{bmatrix}x & x & x & x & x \\ & x & x & x & x \\ & & x & x & x \\ & & & &     x\end{bmatrix}$$
>```
>and
>```
>$$
>A =
>\begin{bmatrix}
>2 & 1 & 0 \\
>1 & 3 & 1 \\
>0 & 1 & 2 \\
>\end{bmatrix}
>$$
>```
> are acceptable, just make sure that a) both `$$` are in the same line, or b) place `$$` in a single line separately.



3.2 deploy website locally.

```
npm i
```

```
npx quartz create
```

- Make sure to select to select "empty quartz folder" option.

```
npx quartz build --serve
```

- This should create a local server to preview your changes.

- This step is necessary since some content would not render correctly in QuartZ server. Please make sure that the pages render in the desired way before make a pull request.


4. Commit your changes.

```
git add.
git commit -m "yyyymmdd: My changes by YourName"
```

5. Push your changes to your branch.

```
git push origin my-branch
```

6. Create a pull request to the original repository.

7. Wait for the review and merge of your changes.

## Contact

If you have any questions or concerns, please contact me at <EMAIL>.


## License
This project is licensed under the CC-BY-NC-SA 4.0 license.
