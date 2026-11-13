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

// In tables, short code (e.g. ("M", 2, "s")) is kept on one line by custom.css;
// a long one (e.g. a full multipole key) may wrap at its spaces.
document.addEventListener("DOMContentLoaded", function() {
    document.querySelectorAll('.bd-content table.table code.literal').forEach(function(code) {
        if (code.textContent.length > 24) {
            code.classList.add('table-code-wrap');
        }
    });
});
