Igen. Egy **tetszőleges elemekből álló 3×3-as mátrixhoz** pont jól jön a `\newcommand`, mert legfeljebb 9 argumentumot fogad, nekünk pedig éppen 9 kell.

Az első Markdown cellába:

```
$$
\newcommand{\matIII}[9]{
\begin{bmatrix}
#1 & #2 & #3\\
#4 & #5 & #6\\
#7 & #8 & #9
\end{bmatrix}
}
$$
```

Ezután bárhol:

```
$$
\mathbf{A}=
\matIII
{a_{11}}{a_{12}}{a_{13}}
{a_{21}}{a_{22}}{a_{23}}
{a_{31}}{a_{32}}{a_{33}}
$$
```

ami ezt adja:

\[
\mathbf A=
\begin{bmatrix}
a_{11}&a_{12}&a_{13}\\
a_{21}&a_{22}&a_{23}\\
a_{31}&a_{32}&a_{33}
\end{bmatrix}
\]

### Ha csak egy indexelt mátrixot szeretnél

Ha gyakran kell olyan mátrix, ahol az elemek ugyanannak a változónak az $a_{ij}$ komponensei, akkor még kényelmesebb egyetlen argumentummal:

```
$$
\newcommand{\matthree}[1]{
\begin{bmatrix}
#1_{11} & #1_{12} & #1_{13}\\
#1_{21} & #1_{22} & #1_{23}\\
#1_{31} & #1_{32} & #1_{33}
\end{bmatrix}
}
$$
```

és akkor:

```
$$
\mathbf{A}=\matthree{a}
$$
```

eredménye:

\[
\mathbf A=
\begin{bmatrix}
a_{11}&a_{12}&a_{13}\\
a_{21}&a_{22}&a_{23}\\
a_{31}&a_{32}&a_{33}
\end{bmatrix}.
\]

### FEM-hez még hasznosabb változat

Ha például a kis alakváltozási tenzort akarod sokszor kiírni, én inkább közvetlenül erre készítenék makrót:

```
$$
\newcommand{\straintensor}{
\begin{bmatrix}
\varepsilon_x &
\frac{1}{2}\gamma_{xy} &
\frac{1}{2}\gamma_{xz}\\
\frac{1}{2}\gamma_{xy} &
\varepsilon_y &
\frac{1}{2}\gamma_{yz}\\
\frac{1}{2}\gamma_{xz} &
\frac{1}{2}\gamma_{yz} &
\varepsilon_z
\end{bmatrix}
}
$$
```

aztán egyszerűen:

```
$$
\boldsymbol{\varepsilon}=\straintensor
$$
```

Én az **általános `\matIII` makrót** mindenképpen betenném a notebook-preambulumba; FEM-es jegyzetben meglepően sokszor elő fog kerülni.
