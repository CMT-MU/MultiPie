---
orphan: true
---

# 概要

英語版は[こちら](overview.md)です。

**MultiPie** は、結晶の点群、空間群、磁気点群、磁気空間群のデータベースであり、それらを群論的に扱うツールです。これらの機能を使って、点群と空間群の両方について **対称性適合多極子基底（SAMB）** を構築できます。

SAMB は完全な基底をなすので、任意の物理量をその対称性の分類に従って記述するのに役立ちます。また、どのような物理量も SAMB で分解できます。

MultiPie の簡易マニュアルとチュートリアルの PDF は、{ref}`ユーザーガイド <ja-brief-manual>`の最後にあります。また、サポートしているすべての群の情報は[群の情報](https://cmt-mu.github.io/MultiPieGen/src/group_info.html)にまとめています。

[ユーザーガイド](guide/users_guide_ja.md)では、調和関数と SAMB の読み方と使い方を（英語と日本語で）説明しています。[Group クラスの使い方](group/group_ja.md)では、`Group` のプロパティとメソッドを例とともに示し、[物質の例](materials_ja.md)では、SrVO$_3$、MoS$_2$、Te のモデルを順を追って作ります。Ver.2 の利用者のために、[Ver.2 からの変更点](changes_ja.md)で Ver.3 での変更点をまとめています。

Python パッケージ [QtDraw](https://github.com/CMT-MU/QtDraw) をインストールすると、結晶や分子の構造、極性ベクトルと軸性ベクトル、軌道、SAMB を可視化できます。この可視化の機能により、群論を見通しよく直観的に学ぶのに特に役立ちます。QtDraw を使った MultiPie の実践的なチュートリアルは、{ref}`QtDraw を用いたチュートリアル <ja-brief-manual>`にあります。

API については、次を参照してください：

- [Core モジュール](api_core.md)：`Group`、`MaterialModel`、`ModelAnalyzer` と、作成・解析の関数（`cmd`）
- [Utility モジュール](api_util.md)：実装のためのさまざまなユーティリティ
