---
title: Welcome to Quartz
---

This is a blank Quartz installation.
See the [documentation](https://quartz.jzhao.xyz) for how to get started.

行间公式含有`\begin{XXX}`和`\end{XXX}`，`$$`， 以及 `\begin` 和 `\end` 等命令，需要单独一行，否则会报错。

```markdown
$$
\begin{align}  
f'(x) &= 4x^3 - 21x^2 + 36x - 20 \\
f''(x) &= 12x^2 - 42x + 36 \\
f'''(x) &= 24x - 42 \\
\end{align}
$$
```

连续两个行间公式之间没有任何文本的话，两个`$$`之间需要空2行。

commit 格式

```
date: yyyymmdd
title: commit message
author: author name
```

例如
```
20240926: update some demo pages -by Jiefeng
```

