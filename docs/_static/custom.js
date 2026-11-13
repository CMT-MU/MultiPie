document.addEventListener("DOMContentLoaded", function() {
    var links = document.querySelectorAll('a');
    links.forEach(function(link) {
        if (link.href.toLowerCase().endsWith('.pdf')) {
            link.removeAttribute('download');
            link.setAttribute('target', '_blank');
            link.setAttribute('rel', 'noopener noreferrer');
        }
    });
});

// In tables, short code (e.g. ("M", 2, "s")) is kept on one line by custom.css.
// A long one (e.g. a full multipole key, a module path or a file pattern) may
// wrap: at its spaces, and after "/", ".", ",", "_" and before "{" (<wbr>, so
// that copied text is unchanged).
document.addEventListener("DOMContentLoaded", function() {
    document.querySelectorAll('.bd-content table.table code.literal').forEach(function(code) {
        if (code.textContent.length <= 24) {
            return;
        }
        code.classList.add('table-code-wrap');
        code.querySelectorAll('span').forEach(function(span) {
            var parts = span.textContent.split(/(?<=[\/.,_])|(?=\{)/);
            if (parts.length < 2) {
                return;
            }
            span.textContent = '';
            parts.forEach(function(part, i) {
                if (i > 0) {
                    span.appendChild(document.createElement('wbr'));
                }
                span.appendChild(document.createTextNode(part));
            });
        });
    });
});
