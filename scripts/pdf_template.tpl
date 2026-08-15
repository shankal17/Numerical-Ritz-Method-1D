\documentclass[12pt,a4paper]{article}

\usepackage[a4paper,text={16.5cm,25.2cm},centering]{geometry}
\usepackage{lmodern}
\usepackage{amssymb,amsmath}
\usepackage{bm}
\usepackage{graphicx}
\setkeys{Gin}{width=\linewidth,keepaspectratio}
\usepackage{float}
\floatplacement{figure}{H}
\usepackage{microtype}
\usepackage{hyperref}
{{#:tex_deps}}
{{{ :tex_deps }}}
{{/:tex_deps}}
\renewcommand{\caption}[1]{}
\setlength{\parindent}{0pt}
\setlength{\parskip}{1.2ex}

\hypersetup
       {   pdfauthor = { {{{:author}}} },
           pdftitle={ {{{:title}}} },
           colorlinks=TRUE,
           linkcolor=black,
           citecolor=blue,
           urlcolor=blue
       }

{{#:title}}
\title{ {{{ :title }}} }
{{/:title}}

{{#:author}}
\author{ {{{ :author }}} }
{{/:author}}

{{#:date}}
\date{ {{{ :date }}} }
{{/:date}}

{{ :highlight }}

\definecolor{codebg}{RGB}{251,251,251}
\definecolor{codeframe}{RGB}{204,204,204}
\lstset{
    backgroundcolor=\color{codebg},
    frame=single,
    rulecolor=\color{codeframe},
    framesep=5pt,
    xleftmargin=2pt,
    xrightmargin=2pt,
    framextopmargin=3pt,
    framexbottommargin=3pt,
}

\begin{document}

{{#:title}}\maketitle{{/:title}}

{{{ :body }}}

\end{document}
