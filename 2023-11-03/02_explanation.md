## 「ドミノタイリングの数え上げ」の半分全列挙

### 問題

$6\times 7$ のマス目の各マスを白か黒に塗る。 $1\times 2$ のドミノ（縦横いずれも可）を互いに重ならないように配置し、黒に塗ったマスをちょうど覆うようにする。ある塗り方について、可能なドミノの配置の個数が $K$ であったとき、 **スコア** $K$ を達成したということにする。

達成可能なスコア、およびそのうち辞書順最小の塗り方を列挙せよ。

### アプローチ

$6\times 7$ の盤面を上下に分けて、つまり $3\times 7$ の領域 $2$ つに分割して探索する。

横軸と平行に走査線をとって走査線と重なる位置にドミノを置きつつ、走査線を分割の境界に向けて移動させる。ドミノを置くことができるかどうかの判断は、走査線を横切るドミノの情報のみからできる。

実際の目的は、ドミノの敷き詰め方を列挙するのではなく、敷き詰め方の個数ごとに図形を分類することであった。各図形について、走査線を横切るドミノの配置で場合分けしつつ、ドミノの敷き詰め方の個数を数えておく。このとき、 $2$ つの図形であって、走査線を横切るドミノがどの配置であっても敷き詰め方の個数が一致するようなものがある場合、その後に敷き詰め方の個数が相異なるようにはならないから、片方の探索をやめる。実際、 $3\times 7$ の領域を走査し終わったときに相異なる場合をもつような図形の個数は $16294$ であって、愚直に考えうる図形の個数である $2^{21}=2097152$ よりは小さい。探索する図形を絞る際の選択は、元ネタアプリの仕様に則って、辞書順最小を残すことにする。

最後に、上下からそれぞれ $1$ つずつ図形を取り出し、走査線と重なるドミノの配置が一致するような場合において全体の敷き詰め方の個数を計算する。考える場合は最大で $16294^2$ 通りであり、偶奇を考えると半分程度になる。いずれにせよ、数秒～数十秒程度の計算に相当する。

### 結果

このプログラムの実行結果によると、以下のことがわかる。

正整数であって達成不可能なスコアを小さい順に並べた数列 :

$$1637,1837,1857,1945,1973,2059,2073,2081,2091,2105,\ldots$$

達成可能な相異なるスコアの個数（ $0$ を含む） : $3411$
