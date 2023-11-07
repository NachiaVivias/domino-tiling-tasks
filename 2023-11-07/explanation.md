
# 追加実験

ボトルネック部の計算を $5$ 倍程度に高速化したので、 $6\times 7$ 以外の場合についても探索を行った。

（ただし、テンプレートパラメータで再帰しているため、コンパイルにかなり時間がかかる。）

## 結果

獲得可能なスコアの種類数

$$
\begin{array}{c|cccccccc}
{}_{\text{col}} \backslash {}^{\text{row}}
    &    2 &    3 &    4 &    5 &    6 &     7 &     8 &     9 \cr \hline
  2 &    3 &    4 &    5 &    7 &    9 &    12 &    15 &    20 \cr
  3 &      &    5 &    9 &   13 &   26 &    40 &    74 &   120 \cr
  4 &      &      &   14 &   38 &   81 &   173 &   376 &   816 \cr
  5 &      &      &      &   74 &  275 &   732 &  2099 &  5946 \cr
  6 &      &      &      &      &  755 &  3411 & 12294 & 43668 \cr
  7 &      &      &      &      &      & 11652 & 69950 &       
\end{array}
$$

獲得不可能なスコアのうち、最小の正整数

$$
\begin{array}{c|cccccccc}
{}_{\text{col}} \backslash {}^{\text{row}}
    &    2 &    3 &    4 &    5 &    6 &     7 &     8 &     9 \cr \hline
  2 &    3 &    4 &    4 &    6 &    7 &     7 &     7 &     7 \cr
  3 &      &    5 &    8 &   12 &   23 &    31 &    31 &    69 \cr
  4 &      &      &    9 &   31 &   56 &    99 &   209 &   439 \cr
  5 &      &      &      &   43 &  159 &   409 &   971 &  2949 \cr
  6 &      &      &      &      &  361 &  1637 &  6343 & 19315 \cr
  7 &      &      &      &      &      &  4999 & 31301 &       
\end{array}
$$

獲得可能なスコアの最大値

$$
\begin{array}{c|cccccccc}
{}_{\text{col}} \backslash {}^{\text{row}}
    &    2 &    3 &    4 &    5 &    6 &      7 &       8 &      9 \cr \hline
  2 &    2 &    3 &    5 &    8 &   13 &     21 &      34 &     55 \cr
  3 &      &    4 &   11 &   15 &   41 &     56 &     153 &    209 \cr
  4 &      &      &   36 &   95 &  281 &    781 &    2245 &   6336 \cr
  5 &      &      &      &  196 & 1183 &   2479 &   14824 &  31329 \cr
  6 &      &      &      &      & 6728 &  31529 &  167089 & 817991 \cr
  7 &      &      &      &      &      & 103984 & 1292697 &       
\end{array}
$$