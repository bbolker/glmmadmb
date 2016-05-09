(TeX-add-style-hook "glmmADMB"
 (lambda ()
    (LaTeX-add-bibliographies
     "glmmadmb")
    (LaTeX-add-labels
     "fig:owl1"
     "fig:owl2")
    (TeX-add-symbols
     '("code" 1)
     '("fixme" 1)
     "R"
     "Splus")
    (TeX-run-style-hooks
     "fancyvrb"
     "alltt"
     "hyperref"
     "url"
     "babel"
     "american"
     "color"
     "usenames"
     "graphicx"
     "inputenc"
     "utf8"
     "latex2e"
     "art10"
     "article")))

