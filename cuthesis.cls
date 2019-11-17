%% cuthesis.cls
%% Version 1.0 (2013/02/27)
%% Cardiff University, School of Physics and Astronomy Thesis format
%% Copyright (C) Cardiff University 2012 --
%%
%% The presentation set out below follows the guidelines in the following
%% document
%%
%% http://www.cardiff.ac.uk/regis/sfs/regs/1213acadregs/1.18%20-%20Presentation%20and%20Submission%20of%20Research%20Degree%20Theses.doc%202.pdf
%
% This work may be distributed and/or modified under the
% conditions of the LaTeX Project Public License, either version 1.3
% of this license or (at your option) any later version.
% The latest version of this license is in
%   http://www.latex-project.org/lppl.txt
% and version 1.3 or later is part of all distributions of LaTeX
% version 2005/12/01 or later.
%
% This work has the LPPL maintenance status `maintained'.
%
% The Current Maintainer of this work is Duncan M. Macleod

\NeedsTeXFormat{LaTeX2e} \ProvidesClass{cuthesis}[2012/12/12 v 0.01 %
Cardiff University, School of Physics and Astronomy thesis format]

% declare our options
\newif\if@wantfancyhdr \@wantfancyhdrfalse
\DeclareOption{fancyhdr}{\@wantfancyhdrtrue}

% Declare options from report class that we want to pass.
\def\@optptsize{12pt}
\DeclareOption{10pt}{\def\@optptsize{10pt}}
\DeclareOption{11pt}{\def\@optptsize{11pt}}
\DeclareOption{12pt}{\def\@optptsize{12pt}}
\DeclareOption{draft}{\PassOptionsToClass{draft}{report}}
\DeclareOption{final}{\PassOptionsToClass{final}{report}}
\DeclareOption{titlepage}{\PassOptionsToClass{titlepage}{report}}
\DeclareOption{notitlepage}{\PassOptionsToClass{notitlepage}{report}}
\DeclareOption{oneside}{\PassOptionsToClass{oneside}{report}}
\DeclareOption{twoside}{\PassOptionsToClass{twoside}{report}}

% process options
\ExecuteOptions{12pt,twoside,fancyhdr}
\ProcessOptions
\LoadClass[a4paper,\@optptsize]{report}

% load other packages
\RequirePackage[top=3cm,bottom=3cm,left=3cm,right=4cm]{geometry}
\RequirePackage{color}

\global\let\@title\@empty
\global\let\@author\@empty
\global\let\@date\@empty

\newcommand{\Chapter}[2]{\chapter[#1]{#1\\[1ex]\large\textit{#2}}}

% define abstract
\renewenvironment{abstract}
{
    \clearpage
    \thispagestyle{plain}
    \chapter*{Summary of thesis}
}
{
    \clearpage
}

% define dedication
\newenvironment{dedication}
{
   \clearpage
   \thispagestyle{plain}
   \vspace*{\stretch{1}}
   \hfill\begin{minipage}[t]{0.66\textwidth}
}%
{
   \end{minipage}
   \vspace*{\stretch{3}}
   \cleardoublepage
}

% matter clauses from book class
\newcommand\frontmatter{%
    \pagestyle{plain}
    \pagenumbering{roman}
    \cleardoublepage
}
\newcommand\mainmatter{%
    \cleardoublepage
    \pagestyle{thesis}
    \onehalfspacing
    \pagenumbering{arabic}
}
\newcommand\backmatter{%
    \cleardoublepage
}

% define fancy header
\if@wantfancyhdr
    \usepackage{fancyhdr}
    \setlength{\headheight}{15pt}
    \definecolor{gray}{rgb}{0.80,0.80,0.80}% dark grey
    \renewcommand{\headrulewidth}{0.1pt}
    \renewcommand{\footrulewidth}{\headrulewidth}
    \renewcommand{\headrule}{{\color{gray}%
        \hrule width\headwidth height\headrulewidth \vskip-\headrulewidth}}
    \renewcommand{\footrulewidth}{\headrulewidth}
    \renewcommand{\footrule}{{\color{gray}%
        \vskip-\footruleskip\vskip-\footrulewidth
        \hrule width\headwidth height\footrulewidth\vskip\footruleskip}}
    \newcommand{\fancyspacer}{\rule[-2ex]{0pt}{2ex}}

    \fancypagestyle{plain}{%
        \fancyhf{}
        \renewcommand{\headrulewidth}{0pt}
        \renewcommand{\footrulewidth}{0.1pt}
        \fancyfoot[c]{-- \small{\thepage}\ --}
    }

    \fancypagestyle{thesis}{%
        \fancyhf{}
        \renewcommand{\headrulewidth}{0.1pt}
        \renewcommand{\footrulewidth}{0.1pt}
        \fancyhead[LO]{\fancyspacer\footnotesize\textit{\nouppercase{\leftmark}}}
        \fancyhead[RE]{\fancyspacer\footnotesize\textit{\nouppercase{\rightmark}}}
        \fancyfoot[c]{-- \thepage\ --}
    }
    \pagestyle{thesis}
\fi

% define declaration
%\newcommand\declaration{
%    \chapter*{Declaration of authorship}
%    \begin{itemize}
%        \item \textsc{Declaration:}\\
%        This work has not previously been accepted in substance for any
%        degree and is not concurrently submitted in candidature for any degree.
%        \\
%        \\
%        Signed: $_{..................................................................}$
%        (candidate)\ \ Date: $_{...................}$
%        \\
%        \item \textsc{Statement 1:}\\
%        This thesis is being submitted in partial fulfillment of the
%        requirements for the degree of Doctor of Philosophy (PhD).
%        \\
%        \\
%        Signed: $_{..................................................................}$
%        (candidate)\ \ Date: $_{...................}$
%        \\
%        \item \textsc{Statement 2:}\\
%        This thesis is the result of my own independent work/investigation,
%        except where otherwise stated. Other sources are acknowledged by
%        explicit references.
%        \\
%        \\
%        Signed: $_{..................................................................}$
%        (candidate)\ \ Date: $_{...................}$
%        \\
%        \item \textsc{Statement 3}\\
%        I hereby give consent for my thesis, if accepted, to be available for
%        photocopying and for inter-library loan, and for the title and summary
%        to be made available to outside organisations.
%        \\
%        \\
%        Signed: $_{..................................................................}$
%        (candidate)\ \ Date: $_{...................}$
%      \\
%    \end{itemize}
%    \clearpage
%}

% define titlepage
\renewcommand\maketitle{\begin{titlepage}%
    \vspace*{20mm}
    \begin{center}
        \begin{minipage}[b]{80mm}{%
            \Huge\bf School of Physics \\and Astronomy
            \vspace*{10mm}
        }
        \end{minipage}
        \hfill
        \begin{minipage}[t]{40mm}
            \makebox[40mm]{\includegraphics[width=40mm]{culogo.pdf}}
        \end{minipage}\\
        \vspace{20mm}
    {\LARGE \@title \par}%
    \vskip 3em%
    {\large \@author \par}%
    \vskip 1.5em%
    {\large
        Submitted for the degree of Doctor of Philosophy \\
        School of Physics and Astronomy \\
        Cardiff University\par}
    \vskip 1.5em%
    {\large \@date \par}%       % Set date in \large size.
    \@thanks
    \vfil\null
    \end{center}
    \end{titlepage}%
    \setcounter{footnote}{0}%
    \global\let\thanks\relax
    \global\let\maketitle\relax
    \global\let\@thanks\@empty
    \global\let\@author\@empty
    \global\let\@date\@empty
    \global\let\@title\@empty
    \global\let\title\relax
    \global\let\author\relax
    \global\let\date\relax
    \global\let\and\relax
}
