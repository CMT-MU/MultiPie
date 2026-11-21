---
orphan: true
---

# 原子多極子演算子

英語版は[こちら](atomic_multipole_table.md)です。

このページは、`Group.atomic_multipole_operator`（{ref}`演算子表現 <ja-atomic-multipole-operator>`を参照）が与える、ランク $l = 2$ までの原子多極子演算子 $X_{lm}^{(s)}(k)$ の一覧です。定義は H. Kusunose, R. Oiwa, and S. Hayami, J. Phys. Soc. Jpn. **89**, 104704 (2020) に従います。

- **成分。** 実の（tesseral の）成分を示します。$X_{l0}$ と、$m > 0$ の $C_{lm} = [(-1)^m X_{lm} + X_{l,-m}]/\sqrt{2}$、$S_{lm} = [(-1)^m X_{lm} - X_{l,-m}]/(\sqrt{2}\,i)$ です。成分の名前は軌道に合わせ、$l = 1$ では $(x, y, z) = (C_{11}, S_{11}, C_{10})$、$l = 2$ では $(u, v, yz, xz, xy) = (C_{20}, C_{22}, S_{21}, C_{21}, S_{22})$ です。例えば $Q_{u} = z^2 - (x^2 + y^2)/2$ です。$X_{l,-m} = (-1)^m X_{lm}^\dagger$（式 6 と 14）なので、これらの成分はエルミートです。
- **演算子。** 位置 $x, y, z$、軌道角運動量 $l_x, l_y, l_z$（無次元）、パウリ行列 $\sigma_x, \sigma_y, \sigma_z$（スピンは $\boldsymbol{\sigma}/2$）で、$[l_i, l_j] = i\epsilon_{ijk} l_k$、$[\sigma_i, \sigma_j] = 2i\epsilon_{ijk}\sigma_k$、$[l_i, x_j] = i\epsilon_{ijk} x_k$ に従います。$x_i$ どうしは交換し、$\sigma_i$ は $l_i$ とも $x_i$ とも交換します。
- **演算子の順序。** 定義の順序を保ちます。例えば積 $l_x\, y$ では、$l_x$ は $y$ のあとに作用します。連続する位置の演算子の並びを並べ替えることと、パウリ行列を右端に置くことだけを行い、ほかに交換関係は使わないので、例えば $l_x y$ と $y l_x$ は別の項として現れます。
- **記法。** $s = 0$ では $X_{\gamma}$、$s = 1$ では $X^{(1)}_{\gamma}(k)$ と、成分の名前 $\gamma$ をつけて書きます（例えば $Q_{u}$ や $Q^{(1)}_{x}(0)$）。ランク 0 は $X_{0}$ と $X^{(1)}_{0}(k)$ です。軌道基底での行列要素は、対称性への適合と規格化の前の原子 SAMB のもので、動径の因子は $\langle r^k \rangle = 1$ です。
- **ランク 1 のベクトル表示。** ランク 1 の 3 つの成分は、例えば $\boldsymbol{T} = (T_{x}, T_{y}, T_{z})$ や $\boldsymbol{Q}^{(1)}(0) = (Q^{(1)}_{x}(0), Q^{(1)}_{y}(0), Q^{(1)}_{z}(0))$ のように、まとめてベクトルで書きます。$\boldsymbol{r} = (x, y, z)$、$\boldsymbol{l} = (l_x, l_y, l_z)$、$\boldsymbol{\sigma} = (\sigma_x, \sigma_y, \sigma_z)$ です。ここでも積の中の演算子の順序を保ち、例えば $(\boldsymbol{r} \times \boldsymbol{l})_x = y l_z - z l_y$、$(\boldsymbol{l} \times \boldsymbol{r})_x = l_y z - l_z y$ です。ランク 0 のスカラーも同じように書きます。$\boldsymbol{l} = -i\,\boldsymbol{r} \times \nabla$ では $\boldsymbol{r} \cdot \boldsymbol{l} = \boldsymbol{l} \cdot \boldsymbol{r} = 0$ なので、これらを含む項は省いています。`Group.atomic_multipole_operator` は交換関係を使わないので、これらの項も与えます。

0 になるもの（例えば $M_{0}$ や、$s = 0$ の $\boldsymbol{G}$）は載せていません。

```{include} examples/atomic_multipole/table_ja.md
```
