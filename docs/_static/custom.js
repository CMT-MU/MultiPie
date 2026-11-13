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

// Sidebar: "[en, ja]" after the entries that have a Japanese version. The Japanese pages are not in the table of
// contents (orphan); while one is open, its entry is marked as the current page. To add a pair, add a line to
// LANGUAGE_PAIRS (paths from the root of the docs).
var LANGUAGE_PAIRS = {
    "src/guide.html": "src/guide/users_guide_ja.html",
    "src/group.html": "src/group/group_ja.html",
    "src/changes.html": "src/changes_ja.html",
};

document.addEventListener("DOMContentLoaded", function() {
    var root = new URL(document.documentElement.dataset.content_root || "./", location.href);
    var here = location.pathname;
    document.querySelectorAll('.bd-sidebar-primary nav.bd-links a.reference.internal').forEach(function(link) {
        var target = new URL(link.getAttribute('href'), location.href).pathname;
        Object.keys(LANGUAGE_PAIRS).forEach(function(en) {
            var enUrl = new URL(en, root);
            if (target !== enUrl.pathname) {
                return;
            }
            var jaUrl = new URL(LANGUAGE_PAIRS[en], root);
            var span = document.createElement('span');
            span.className = 'language-links';
            [["en", enUrl], ["ja", jaUrl]].forEach(function(item, i) {
                if (i > 0) {
                    span.appendChild(document.createTextNode(", "));
                }
                var a = document.createElement('a');
                a.href = item[1].href;
                a.textContent = item[0];
                a.lang = item[0];
                if (here === item[1].pathname) {
                    a.className = 'current-language';
                }
                span.appendChild(a);
            });
            span.insertBefore(document.createTextNode("["), span.firstChild);
            span.appendChild(document.createTextNode("]"));
            link.appendChild(document.createTextNode(" "));
            link.after(span);
            if (here === jaUrl.pathname) {
                link.classList.add('current');
                link.parentElement.classList.add('current', 'active');
            }
        });
    });
});
