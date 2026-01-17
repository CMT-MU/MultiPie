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
