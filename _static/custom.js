document.addEventListener("DOMContentLoaded", function() {
    // 全ての <a> タグをチェック
    var links = document.querySelectorAll('a');
    links.forEach(function(link) {
        // リンク先が .pdf で終わる場合
        if (link.href.toLowerCase().endsWith('.pdf')) {
            // ダウンロード属性を削除（もしあれば）
            link.removeAttribute('download');
            // 新しいタブで開くように設定
            link.setAttribute('target', '_blank');
            link.setAttribute('rel', 'noopener noreferrer');
        }
    });
});
