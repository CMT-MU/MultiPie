document.addEventListener("DOMContentLoaded", function () {
  const container = document.getElementById("pdf-container");
  const frame = document.getElementById("pdf-frame");
  if (!container || !frame) return;

  document.querySelectorAll("a[href$='.pdf']").forEach(link => {
    link.addEventListener("click", function (e) {
      e.preventDefault();

      const pdf = this.getAttribute("href");

      if (frame.src && frame.src.endsWith(pdf)) {
        container.style.display = "none";
        frame.src = "";
        return;
      }

      frame.src = pdf;
      container.style.display = "block";
      container.scrollIntoView({ behavior: "smooth" });
    });
  });
});
