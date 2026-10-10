// Links that leave the site open in a new tab, so the reader keeps their place.
document.addEventListener("DOMContentLoaded", function () {
  document.querySelectorAll("a[href^='http']").forEach(function (link) {
    if (link.hostname !== window.location.hostname) {
      link.target = "_blank";
      link.rel = "noopener";
    }
  });
});
