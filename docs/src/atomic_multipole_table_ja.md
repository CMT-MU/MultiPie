---
orphan: true
---

# 原子多極子演算子

英語版は[こちら](atomic_multipole_table.md)です。

このページは、`Group.atomic_multipole_operator`（{ref}`演算子表現 <ja-atomic-multipole-operator>`を参照）が与える、ランク $l = 2$ までの原子多極子演算子 $X_{lm}^{(s)}(k)$ の一覧です。定義は H. Kusunose, R. Oiwa, and S. Hayami, J. Phys. Soc. Jpn. **89**, 104704 (2020) に従います。

- **成分。** 実の（tesseral の）成分を示します。$X_{l0}$ と、$m > 0$ の $C_{lm} = [(-1)^m X_{lm} + X_{l,-m}]/\sqrt{2}$、$S_{lm} = [(-1)^m X_{lm} - X_{l,-m}]/(\sqrt{2}\,i)$ です。成分の名前は軌道に合わせ、$l = 1$ では $(x, y, z) = (C_{11}, S_{11}, C_{10})$、$l = 2$ では $(u, v, yz, xz, xy) = (C_{20}, C_{22}, S_{21}, C_{21}, S_{22})$ です。例えば $Q_{2,u} = z^2 - (x^2 + y^2)/2$ です。$X_{l,-m} = (-1)^m X_{lm}^\dagger$（式 6 と 14）なので、これらの成分はエルミートです。
- **演算子。** 位置 $x, y, z$、軌道角運動量 $l_x, l_y, l_z$（無次元）、パウリ行列 $\sigma_x, \sigma_y, \sigma_z$（スピンは $\boldsymbol{\sigma}/2$）で、$[l_i, l_j] = i\epsilon_{ijk} l_k$、$[\sigma_i, \sigma_j] = 2i\epsilon_{ijk}\sigma_k$、$[l_i, x_j] = i\epsilon_{ijk} x_k$ に従います。$x_i$ どうしは交換し、$\sigma_i$ は $l_i$ とも $x_i$ とも交換します。
- **演算子の順序。** 定義の順序を保ちます。例えば積 $l_x\, y$ では、$l_x$ は $y$ のあとに作用します。連続する位置の演算子の並びを並べ替えることと、パウリ行列を右端に置くことだけを行い、ほかに交換関係は使わないので、例えば $l_x y$ と $y l_x$ は別の項として現れます。
- **記法。** $s = 0$ では $X_{l,\gamma}$、$s = 1$ では $X^{(1)}_{l,\gamma}(k)$ と書きます。軌道基底での行列要素は、対称性への適合と規格化の前の原子 SAMB のもので、動径の因子は $\langle r^k \rangle = 1$ です。

0 になるもの（例えば $M_{0}$ や、スピンによらない $G_{1}$）は載せていません。一覧は `docs/src/examples/atomic_multipole/create_table.py` で作っています。

```{include} examples/atomic_multipole/table_ja.md
```
