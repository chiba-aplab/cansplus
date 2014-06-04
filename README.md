# CANS+ #

CANS+は、CANS (Coordinated Astronomical Numeical Software、[開発HP](http://www-space.eps.s.u-tokyo.ac.jp/~yokoyama/etc/cans/)) から発展した、高次精度MHDコードです。特徴として、

* HLLD／HLL近似リーマン解法を用いた有限体積法
* MP5法を用いた高次補間（5次精度） 
* 9 wave法を用いたdivBクリーニング
* カーテシアン／円筒座標系
* 非一様メッシュサイズ
* MPIによる３次元領域分割化
* IDLによる読み込み・可視化ルーチンの整備

が挙げられます。

構成はCANSに倣って、共通エンジン部分と物理課題に分けて、ライブラリをリンクする形で、共通エンジンを各課題から利用します。

開発はmercurialでバージョン管理され、bitbucket上にホスティングされています。ダウンロードはプロジェクトサイトからできますので、bitbucketのアカウントを作成の上、CANS+の[HP「問い合わせ」](http://www.astro.phys.s.chiba-u.ac.jp/cans/?page_id=28)よりアカウント名をご連絡ください。