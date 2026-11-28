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
    "src/overview.html": "src/overview_ja.html",
    "src/install.html": "src/install_ja.html",
    "src/getting_started.html": "src/getting_started_ja.html",
    "src/guide.html": "src/guide/users_guide_ja.html",
    "src/atomic_multipole_table.html": "src/atomic_multipole_table_ja.html",
    "src/group.html": "src/group/group_ja.html",
    "src/materials.html": "src/materials_ja.html",
    "src/post_processing.html": "src/post_processing_ja.html",
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

// Sidebar toggles on narrow screens: pydata-sphinx-theme (0.17) and sphinx-book-theme (1.3) hook their handlers only
// to the first ".primary-toggle" / ".secondary-toggle" of the page, a button of the header navbar that the book theme
// hides, so that the visible buttons above the article do nothing (the sidebars cannot be opened on a smartphone).
// Forward the clicks on the other buttons to the first one.
document.addEventListener("DOMContentLoaded", function() {
    ["primary", "secondary"].forEach(function(side) {
        var toggles = document.querySelectorAll("." + side + "-toggle");
        for (var i = 1; i < toggles.length; i++) {
            toggles[i].addEventListener("click", function(event) {
                event.preventDefault();
                event.stopPropagation();
                toggles[0].click();
            });
        }
    });
});
