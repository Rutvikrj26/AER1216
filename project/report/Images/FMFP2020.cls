%%
%% FMFPtran.cls 2011/11/03 version V1.8 based on
%% FMFPtran.cls 2007/03/05 version V1.7a
%% The changes in V1.8 are made with a single goal in mind:
%% to change the look of the output using the [conference] option
%% and the default font size (10pt) to match the Word template more closely.
%% These changes may well have undesired side effects when other options
%% are in force!
%% 
%% 
%% This is the official FMFP LaTeX class for authors of the Institute of 
%% Electrical and Electronics Engineers (FMFP) Transactions journals and
%% conferences.
%% 
%% Support sites:
%% http://www.michaelshell.org/tex/ieeetran/
%% http://www.ctan.org/tex-archive/macros/latex/contrib/FMFPtran/
%% and
%% http://www.ieee.org/
%%
%% Based on the original 1993 FMFPtran.cls, but with many bug fixes
%% and enhancements (from both JVH and MDS) over the 1996/7 version.
%%
%%
%% Contributors:
%% Gerry Murray (1993), Silvano Balemi (1993),
%% Jon Dixon (1996), Peter N"uchter (1996),
%% Juergen von Hagen (2000), and Michael Shell (2001-2007)
%% 
%% 
%% Copyright (c) 1993-2000 by Gerry Murray, Silvano Balemi, 
%%                         Jon Dixon, Peter N"uchter,
%%                         Juergen von Hagen
%%                         and
%% Copyright (c) 2001-2007 by Michael Shell
%%
%% Current maintainer (V1.3 to V1.7): Michael Shell
%%                                    See:
%%                                    http://www.michaelshell.org/
%%                                    for current contact information.
%%
%% Special thanks to Peter Wilson (CUA) and Donald Arseneau
%% for allowing the inclusion of the \@ifmtarg command 
%% from their ifmtarg LaTeX package. 
%% 
%%*************************************************************************
%% Legal Notice:
%% This code is offered as-is without any warranty either expressed or
%% implied; without even the implied warranty of MERCHANTABILITY or
%% FITNESS FOR A PARTICULAR PURPOSE! 
%% User assumes all risk.
%% In no event shall FMFP or any contributor to this code be liable for
%% any damages or losses, including, but not limited to, incidental,
%% consequential, or any other damages, resulting from the use or misuse
%% of any information contained here.
%%
%% All comments are the opinions of their respective authors and are not
%% necessarily endorsed by the FMFP.
%%
%% This work is distributed under the LaTeX Project Public License (LPPL)
%% ( http://www.latex-project.org/ ) version 1.3, and may be freely used,
%% distributed and modified. A copy of the LPPL, version 1.3, is included
%% in the base LaTeX documentation of all distributions of LaTeX released
%% 2003/12/01 or later.
%% Retain all contribution notices and credits.
%% ** Modified files should be clearly indicated as such, including  **
%% ** renaming them and changing author support contact information. **
%%
%% File list of work: FMFPtran.cls, FMFPtran_HOWTO.pdf, bare_adv.tex,
%%                    bare_conf.tex, bare_jrnl.tex, bare_jrnl_compsoc.tex
%% 
%% Major changes to the user interface should be indicated by an 
%% increase in the version numbers. If a version is a beta, it will 
%% be indicated with a BETA suffix, i.e., 1.4 BETA.
%% Small changes can be indicated by appending letters to the version
%% such as "FMFPtran_v14a.cls".
%% In all cases, \Providesclass, any \typeout messages to the user,
%% \FMFPtransversionmajor and \FMFPtransversionminor must reflect the
%% correct version information.
%% The changes should also be documented via source comments.
%%*************************************************************************
%%
%
% Available class options 
% e.g., \documentclass[10pt,conference]{FMFPtran} 
% 
%             *** choose only one from each category ***
%
% 9pt, 10pt, 11pt, 12pt
%    Sets normal font size. The default is 10pt.
% 
% conference, journal, technote, peerreview, peerreviewca
%    determines format mode - conference papers, journal papers,
%    correspondence papers (technotes), or peer review papers. The user
%    should also select 9pt when using technote. peerreview is like
%    journal mode, but provides for a single-column "cover" title page for
%    anonymous peer review. The paper title (without the author names) is
%    repeated at the top of the page after the cover page. For peer review
%    papers, the \FMFPpeerreviewmaketitle command must be executed (will
%    automatically be ignored for non-peerreview modes) at the place the
%    cover page is to end, usually just after the abstract (keywords are
%    not normally used with peer review papers). peerreviewca is like
%    peerreview, but allows the author names to be entered and formatted
%    as with conference mode so that author affiliation and contact
%    information can be easily seen on the cover page.
%    The default is journal.
%
% draft, draftcls, draftclsnofoot, final
%    determines if paper is formatted as a widely spaced draft (for
%    handwritten editor comments) or as a properly typeset final version.
%    draftcls restricts draft mode to the class file while all other LaTeX
%    packages (i.e., \usepackage{graphicx}) will behave as final - allows
%    for a draft paper with visible figures, etc. draftclsnofoot is like
%    draftcls, but does not display the date and the word "DRAFT" at the foot
%    of the pages. If using one of the draft modes, the user will probably
%    also want to select onecolumn.
%    The default is final.
%
% letterpaper, a4paper
%    determines paper size: 8.5in X 11in or 210mm X 297mm. CHANGING THE PAPER
%    SIZE WILL NOT ALTER THE TYPESETTING OF THE DOCUMENT - ONLY THE MARGINS
%    WILL BE AFFECTED. In particular, documents using the a4paper option will
%    have reduced side margins (A4 is narrower than US letter) and a longer
%    bottom margin (A4 is longer than US letter). For both cases, the top
%    margins will be the same and the text will be horizontally centered. 
%    For final submission to FMFP, authors should use US letter (8.5 X 11in)
%    paper. Note that authors should ensure that all post-processing 
%    (ps, pdf, etc.) uses the same paper specificiation as the .tex document.
%    Problems here are by far the number one reason for incorrect margins.
%    FMFPtran will automatically set the default paper size under pdflatex 
%    (without requiring a change to pdftex.cfg), so this issue is more
%    important to dvips users. Fix config.ps, config.pdf, or ~/.dvipsrc for
%    dvips, or use the dvips -t papersize option instead as needed. See the
%    testflow documentation
%    http://www.ctan.org/tex-archive/macros/latex/contrib/FMFPtran/testflow
%    for more details on dvips paper size configuration.
%    The default is letterpaper.
%
% oneside, twoside
%    determines if layout follows single sided or two sided (duplex)
%    printing. The only notable change is with the headings at the top of
%    the pages.
%    The default is oneside.
%
% onecolumn, twocolumn
%    determines if text is organized into one or two columns per page. One
%    column mode is usually used only with draft papers.
%    The default is twocolumn.
%
% compsoc
%    Use the format of the FMFP Computer Society.
%
% romanappendices
%    Use the "Appendix I" convention when numbering appendices. FMFPtran.cls
%    now defaults to Alpha "Appendix A" convention - the opposite of what
%    v1.6b and earlier did.
%
% captionsoff
%    disables the display of the figure/table captions. Some FMFP journals
%    request that captions be removed and figures/tables be put on pages
%    of their own at the end of an initial paper submission. The endfloat
%    package can be used with this class option to achieve this format.
%
% nofonttune
%    turns off tuning of the font interword spacing. Maybe useful to those
%    not using the standard Times fonts or for those who have already "tuned"
%    their fonts.
%    The default is to enable FMFPtran to tune font parameters.
%
%
%----------
% Available CLASSINPUTs provided (all are macros unless otherwise noted):
% \CLASSINPUTbaselinestretch
% \CLASSINPUTinnersidemargin
% \CLASSINPUToutersidemargin
% \CLASSINPUTtoptextmargin
% \CLASSINPUTbottomtextmargin
%
% Available CLASSINFOs provided:
% \ifCLASSINFOpdf                       (TeX if conditional)
% \CLASSINFOpaperwidth                  (macro)
% \CLASSINFOpaperheight                 (macro)
% \CLASSINFOnormalsizebaselineskip      (length)
% \CLASSINFOnormalsizeunitybaselineskip (length)
%
% Available CLASSOPTIONs provided:
% all class option flags (TeX if conditionals) unless otherwise noted,
% e.g., \ifCLASSOPTIONcaptionsoff
% point size options provided as a single macro:
% \CLASSOPTIONpt
% which will be defined as 9, 10, 11, or 12 depending on the document's
% normalsize point size.
% also, class option peerreviewca implies the use of class option peerreview
% and classoption draft implies the use of class option draftcls





\ProvidesClass{FMFPtran}[2012/11/21 V1.8c by Harald Hanche-Olsen and Anders Christensen]
\typeout{-- Based on V1.7a by Michael Shell}
\typeout{-- See the "FMFPtran_HOWTO" manual for usage information.}
\typeout{-- http://www.michaelshell.org/tex/ieeetran/}
\NeedsTeXFormat{LaTeX2e}

% FMFPtran.cls version numbers, provided as of V1.3
% These values serve as a way a .tex file can
% determine if the new features are provided.
% The version number of this FMFPtrans.cls can be obtained from 
% these values. i.e., V1.4
% KEEP THESE AS INTEGERS! i.e., NO {4a} or anything like that-
% (no need to enumerate "a" minor changes here)
\def\FMFPtransversionmajor{1}
\def\FMFPtransversionminor{7}

% These do nothing, but provide them like in article.cls
\newif\if@restonecol
\newif\if@titlepage


% class option conditionals
\newif\ifCLASSOPTIONonecolumn       \CLASSOPTIONonecolumnfalse
\newif\ifCLASSOPTIONtwocolumn       \CLASSOPTIONtwocolumntrue

\newif\ifCLASSOPTIONoneside         \CLASSOPTIONonesidetrue
\newif\ifCLASSOPTIONtwoside         \CLASSOPTIONtwosidefalse

\newif\ifCLASSOPTIONfinal           \CLASSOPTIONfinaltrue
\newif\ifCLASSOPTIONdraft           \CLASSOPTIONdraftfalse
\newif\ifCLASSOPTIONdraftcls        \CLASSOPTIONdraftclsfalse
\newif\ifCLASSOPTIONdraftclsnofoot  \CLASSOPTIONdraftclsnofootfalse

\newif\ifCLASSOPTIONpeerreview      \CLASSOPTIONpeerreviewfalse
\newif\ifCLASSOPTIONpeerreviewca    \CLASSOPTIONpeerreviewcafalse

\newif\ifCLASSOPTIONjournal         \CLASSOPTIONjournaltrue
\newif\ifCLASSOPTIONconference      \CLASSOPTIONconferencefalse
\newif\ifCLASSOPTIONtechnote        \CLASSOPTIONtechnotefalse

\newif\ifCLASSOPTIONnofonttune      \CLASSOPTIONnofonttunefalse

\newif\ifCLASSOPTIONcaptionsoff     \CLASSOPTIONcaptionsofffalse

\newif\ifCLASSOPTIONcompsoc         \CLASSOPTIONcompsocfalse

\newif\ifCLASSOPTIONromanappendices \CLASSOPTIONromanappendicesfalse


% class info conditionals

% indicates if pdf (via pdflatex) output
\newif\ifCLASSINFOpdf               \CLASSINFOpdffalse


% V1.6b internal flag to show if using a4paper
\newif\if@FMFPusingAfourpaper       \@FMFPusingAfourpaperfalse



% FMFPtran class scratch pad registers
% dimen
\newdimen\@FMFPtrantmpdimenA
\newdimen\@FMFPtrantmpdimenB
% count
\newcount\@FMFPtrantmpcountA
\newcount\@FMFPtrantmpcountB
% token list
\newtoks\@FMFPtrantmptoksA

% we use \CLASSOPTIONpt so that we can ID the point size (even for 9pt docs)
% as well as LaTeX's \@ptsize to retain some compatability with some
% external packages
\def\@ptsize{0}
% LaTeX does not support 9pt, so we set \@ptsize to 0 - same as that of 10pt
\DeclareOption{9pt}{\def\CLASSOPTIONpt{9}\def\@ptsize{0}}
\DeclareOption{10pt}{\def\CLASSOPTIONpt{10}\def\@ptsize{0}}
\DeclareOption{11pt}{\def\CLASSOPTIONpt{11}\def\@ptsize{1}}
\DeclareOption{12pt}{\def\CLASSOPTIONpt{12}\def\@ptsize{2}}



\DeclareOption{letterpaper}{\setlength{\paperheight}{11in}%
                            \setlength{\paperwidth}{8.5in}%
                            \@FMFPusingAfourpaperfalse
                            \def\CLASSOPTIONpaper{letter}%
                            \def\CLASSINFOpaperwidth{8.5in}%
                            \def\CLASSINFOpaperheight{11in}}


\DeclareOption{a4paper}{\setlength{\paperheight}{297mm}%
                        \setlength{\paperwidth}{210mm}%
                        \@FMFPusingAfourpapertrue
                        \def\CLASSOPTIONpaper{a4}%
                        \def\CLASSINFOpaperwidth{210mm}%
                        \def\CLASSINFOpaperheight{297mm}}

\DeclareOption{oneside}{\@twosidefalse\@mparswitchfalse
                        \CLASSOPTIONonesidetrue\CLASSOPTIONtwosidefalse}
\DeclareOption{twoside}{\@twosidetrue\@mparswitchtrue
                        \CLASSOPTIONtwosidetrue\CLASSOPTIONonesidefalse}

\DeclareOption{onecolumn}{\CLASSOPTIONonecolumntrue\CLASSOPTIONtwocolumnfalse}
\DeclareOption{twocolumn}{\CLASSOPTIONtwocolumntrue\CLASSOPTIONonecolumnfalse}

% If the user selects draft, then this class AND any packages
% will go into draft mode.
\DeclareOption{draft}{\CLASSOPTIONdrafttrue\CLASSOPTIONdraftclstrue
                      \CLASSOPTIONdraftclsnofootfalse} 
% draftcls is for a draft mode which will not affect any packages
% used by the document.
\DeclareOption{draftcls}{\CLASSOPTIONdraftfalse\CLASSOPTIONdraftclstrue
                         \CLASSOPTIONdraftclsnofootfalse} 
% draftclsnofoot is like draftcls, but without the footer.
\DeclareOption{draftclsnofoot}{\CLASSOPTIONdraftfalse\CLASSOPTIONdraftclstrue
                               \CLASSOPTIONdraftclsnofoottrue} 
\DeclareOption{final}{\CLASSOPTIONdraftfalse\CLASSOPTIONdraftclsfalse
                      \CLASSOPTIONdraftclsnofootfalse}

\DeclareOption{journal}{\CLASSOPTIONpeerreviewfalse\CLASSOPTIONpeerreviewcafalse
                        \CLASSOPTIONjournaltrue\CLASSOPTIONconferencefalse\CLASSOPTIONtechnotefalse}

\DeclareOption{conference}{\CLASSOPTIONpeerreviewfalse\CLASSOPTIONpeerreviewcafalse
                           \CLASSOPTIONjournalfalse\CLASSOPTIONconferencetrue\CLASSOPTIONtechnotefalse}

\DeclareOption{technote}{\CLASSOPTIONpeerreviewfalse\CLASSOPTIONpeerreviewcafalse
                         \CLASSOPTIONjournalfalse\CLASSOPTIONconferencefalse\CLASSOPTIONtechnotetrue}

\DeclareOption{peerreview}{\CLASSOPTIONpeerreviewtrue\CLASSOPTIONpeerreviewcafalse
                           \CLASSOPTIONjournalfalse\CLASSOPTIONconferencefalse\CLASSOPTIONtechnotefalse}

\DeclareOption{peerreviewca}{\CLASSOPTIONpeerreviewtrue\CLASSOPTIONpeerreviewcatrue
                             \CLASSOPTIONjournalfalse\CLASSOPTIONconferencefalse\CLASSOPTIONtechnotefalse}

\DeclareOption{nofonttune}{\CLASSOPTIONnofonttunetrue}

\DeclareOption{captionsoff}{\CLASSOPTIONcaptionsofftrue}

\DeclareOption{compsoc}{\CLASSOPTIONcompsoctrue}

\DeclareOption{romanappendices}{\CLASSOPTIONromanappendicestrue}


% default to US letter paper, 10pt, twocolumn, one sided, final, journal
\ExecuteOptions{letterpaper,10pt,twocolumn,oneside,final,journal}
% overrride these defaults per user requests
\ProcessOptions



% Computer Society conditional execution command
\long\def\@FMFPcompsoconly#1{\relax\ifCLASSOPTIONcompsoc\relax#1\relax\fi\relax}
% inverse
\long\def\@FMFPnotcompsoconly#1{\relax\ifCLASSOPTIONcompsoc\else\relax#1\relax\fi\relax}
% compsoc conference
\long\def\@FMFPcompsocconfonly#1{\relax\ifCLASSOPTIONcompsoc\ifCLASSOPTIONconference\relax#1\relax\fi\fi\relax}
% compsoc not conference
\long\def\@FMFPcompsocnotconfonly#1{\relax\ifCLASSOPTIONcompsoc\ifCLASSOPTIONconference\else\relax#1\relax\fi\fi\relax}


% FMFP uses Times Roman font, so we'll default to Times.
% These three commands make up the entire times.sty package.
\renewcommand{\sfdefault}{phv}
\renewcommand{\rmdefault}{ptm}
\renewcommand{\ttdefault}{pcr}

\@FMFPcompsoconly{\typeout{-- Using FMFP Computer Society mode.}}

% V1.7 compsoc nonconference papers, use Palatino/Palladio as the main text font,
% not Times Roman.
\@FMFPcompsocnotconfonly{\renewcommand{\rmdefault}{ppl}}

% enable Times/Palatino main text font
\normalfont\selectfont





% V1.7 conference notice message hook
\def\@FMFPconsolenoticeconference{\typeout{}%
\typeout{** Conference Paper **}%
\typeout{Before submitting the final camera ready copy, remember to:}%
\typeout{}%
\typeout{ 1. Manually equalize the lengths of two columns on the last page}%
\typeout{ of your paper;}%
\typeout{}%
\typeout{ 2. Ensure that any PostScript and/or PDF output post-processing}%
\typeout{ uses only Type 1 fonts and that every step in the generation}%
\typeout{ process uses the appropriate paper size.}%
\typeout{}}


% we can send console reminder messages to the user here
\AtEndDocument{\ifCLASSOPTIONconference\@FMFPconsolenoticeconference\fi}


% warn about the use of single column other than for draft mode
\ifCLASSOPTIONtwocolumn\else%
  \ifCLASSOPTIONdraftcls\else%
   \typeout{** ATTENTION: Single column mode is not typically used with FMFP publications.}%
  \fi%
\fi


% V1.7 improved paper size setting code.
% Set pdfpage and dvips paper sizes. Conditional tests are similar to that
% of ifpdf.sty. Retain within {} to ensure tested macros are never altered,
% even if only effect is to set them to \relax.
% if \pdfoutput is undefined or equal to relax, output a dvips special
{\@ifundefined{pdfoutput}{\AtBeginDvi{\special{papersize=\CLASSINFOpaperwidth,\CLASSINFOpaperheight}}}{%
% pdfoutput is defined and not equal to \relax
% check for pdfpageheight existence just in case someone sets pdfoutput
% under non-pdflatex. If exists, set them regardless of value of \pdfoutput.
\@ifundefined{pdfpageheight}{\relax}{\global\pdfpagewidth\paperwidth
\global\pdfpageheight\paperheight}%
% if using \pdfoutput=0 under pdflatex, send dvips papersize special
\ifcase\pdfoutput
\AtBeginDvi{\special{papersize=\CLASSINFOpaperwidth,\CLASSINFOpaperheight}}%
\else
% we are using pdf output, set CLASSINFOpdf flag
\global\CLASSINFOpdftrue
\fi}}

% let the user know the selected papersize
\typeout{-- Using \CLASSINFOpaperwidth\space x \CLASSINFOpaperheight\space
(\CLASSOPTIONpaper)\space paper.}

\ifCLASSINFOpdf
\typeout{-- Using PDF output.}
\else
\typeout{-- Using DVI output.}
\fi


% The idea hinted here is for LaTeX to generate markleft{} and markright{}
% automatically for you after you enter \author{}, \journal{},
% \journaldate{}, journalvol{}, \journalnum{}, etc.
% However, there may be some backward compatibility issues here as
% well as some special applications for FMFPtran.cls and special issues
% that may require the flexible \markleft{}, \markright{} and/or \markboth{}.
% We'll leave this as an open future suggestion.
%\newcommand{\journal}[1]{\def\@journal{#1}}
%\def\@journal{}



% pointsize values
% used with ifx to determine the document's normal size
\def\@FMFPptsizenine{9}
\def\@FMFPptsizeten{10}
\def\@FMFPptsizeeleven{11}
\def\@FMFPptsizetwelve{12}



% FONT DEFINITIONS (No sizexx.clo file needed) 
% V1.6 revised font sizes, displayskip values and
%      revised normalsize baselineskip to reduce underfull vbox problems
%      on the 58pc = 696pt = 9.5in text height we want
%      normalsize     #lines/column  baselineskip (aka leading)
%             9pt     63             11.0476pt (truncated down)
%            10pt     58             12pt      (exact)
%            11pt     52             13.3846pt (truncated down)
%            12pt     50             13.92pt   (exact)
%

% we need to store the nominal baselineskip for the given font size
% in case baselinestretch ever changes.
% this is a dimen, so it will not hold stretch or shrink
\newdimen\@FMFPnormalsizeunitybaselineskip
\@FMFPnormalsizeunitybaselineskip\baselineskip

\ifx\CLASSOPTIONpt\@FMFPptsizenine
\typeout{-- This is a 9 point document.}
\def\normalsize{\@setfontsize{\normalsize}{9}{11.0476pt}}%
\setlength{\@FMFPnormalsizeunitybaselineskip}{11.0476pt}%
\normalsize
\abovedisplayskip 1.5ex plus3pt minus1pt%
\belowdisplayskip \abovedisplayskip%
\abovedisplayshortskip 0pt plus3pt%
\belowdisplayshortskip 1.5ex plus3pt minus1pt
\def\small{\@setfontsize{\small}{8.5}{10pt}}
\def\footnotesize{\@setfontsize{\footnotesize}{8}{9pt}}
\def\scriptsize{\@setfontsize{\scriptsize}{7}{8pt}}
\def\tiny{\@setfontsize{\tiny}{5}{6pt}}
% sublargesize is the same as large - 10pt
\def\sublargesize{\@setfontsize{\sublargesize}{10}{12pt}}
\def\large{\@setfontsize{\large}{10}{12pt}}
\def\Large{\@setfontsize{\Large}{12}{14pt}}
\def\LARGE{\@setfontsize{\LARGE}{14}{17pt}}
\def\huge{\@setfontsize{\huge}{17}{20pt}}
\def\Huge{\@setfontsize{\Huge}{20}{24pt}}
\fi


% Check if we have selected 10 points
\ifx\CLASSOPTIONpt\@FMFPptsizeten
\typeout{-- This is a 10 point document.}
\def\normalsize{\@setfontsize{\normalsize}{10}{11}}%
\setlength{\@FMFPnormalsizeunitybaselineskip}{11pt}%
\normalsize
\abovedisplayskip 1.5ex plus4pt minus2pt%
\belowdisplayskip \abovedisplayskip%
\abovedisplayshortskip 0pt plus4pt%
\belowdisplayshortskip 1.5ex plus4pt minus2pt
\def\small{\@setfontsize{\small}{9}{10pt}}
\def\footnotesize{\@setfontsize{\footnotesize}{8}{9pt}}
\def\scriptsize{\@setfontsize{\scriptsize}{7}{8pt}}
\def\tiny{\@setfontsize{\tiny}{5}{6pt}}
% sublargesize is a tad smaller than large - 11pt
\def\sublargesize{\@setfontsize{\sublargesize}{11}{13.4pt}}
\def\large{\@setfontsize{\large}{12}{14pt}}
\def\Large{\@setfontsize{\Large}{14}{17pt}}
\def\LARGE{\@setfontsize{\LARGE}{17}{20pt}}
\def\huge{\@setfontsize{\huge}{20}{24pt}}
\def\Huge{\@setfontsize{\Huge}{24}{28pt}}
\fi


% Check if we have selected 11 points
\ifx\CLASSOPTIONpt\@FMFPptsizeeleven
\typeout{-- This is an 11 point document.}
\def\normalsize{\@setfontsize{\normalsize}{11}{13.3846pt}}%
\setlength{\@FMFPnormalsizeunitybaselineskip}{13.3846pt}%
\normalsize
\abovedisplayskip 1.5ex plus5pt minus3pt%
\belowdisplayskip \abovedisplayskip%
\abovedisplayshortskip 0pt plus5pt%
\belowdisplayshortskip 1.5ex plus5pt minus3pt
\def\small{\@setfontsize{\small}{10}{12pt}}
\def\footnotesize{\@setfontsize{\footnotesize}{9}{10.5pt}}
\def\scriptsize{\@setfontsize{\scriptsize}{8}{9pt}}
\def\tiny{\@setfontsize{\tiny}{6}{7pt}}
% sublargesize is the same as large - 12pt
\def\sublargesize{\@setfontsize{\sublargesize}{12}{14pt}}
\def\large{\@setfontsize{\large}{12}{14pt}}
\def\Large{\@setfontsize{\Large}{14}{17pt}}
\def\LARGE{\@setfontsize{\LARGE}{17}{20pt}}
\def\huge{\@setfontsize{\huge}{20}{24pt}}
\def\Huge{\@setfontsize{\Huge}{24}{28pt}}
\fi


% Check if we have selected 12 points
\ifx\CLASSOPTIONpt\@FMFPptsizetwelve
\typeout{-- This is a 12 point document.}
\def\normalsize{\@setfontsize{\normalsize}{12}{13.92pt}}%
\setlength{\@FMFPnormalsizeunitybaselineskip}{13.92pt}%
\normalsize
\abovedisplayskip 1.5ex plus6pt minus4pt%
\belowdisplayskip \abovedisplayskip%
\abovedisplayshortskip 0pt plus6pt%
\belowdisplayshortskip 1.5ex plus6pt minus4pt
\def\small{\@setfontsize{\small}{10}{12pt}}
\def\footnotesize{\@setfontsize{\footnotesize}{9}{10.5pt}}
\def\scriptsize{\@setfontsize{\scriptsize}{8}{9pt}}
\def\tiny{\@setfontsize{\tiny}{6}{7pt}}
% sublargesize is the same as large - 14pt
\def\sublargesize{\@setfontsize{\sublargesize}{14}{17pt}}
\def\large{\@setfontsize{\large}{14}{17pt}}
\def\Large{\@setfontsize{\Large}{17}{20pt}}
\def\LARGE{\@setfontsize{\LARGE}{20}{24pt}}
\def\huge{\@setfontsize{\huge}{22}{26pt}}
\def\Huge{\@setfontsize{\Huge}{24}{28pt}}
\fi


% V1.6 The Computer Modern Fonts will issue a substitution warning for
% 24pt titles (24.88pt is used instead) increase the substitution
% tolerance to turn off this warning
\def\fontsubfuzz{.9pt}
% However, the default (and correct) Times font will scale exactly as needed.


% warn the user in case they forget to use the 9pt option with
% technote
\ifCLASSOPTIONtechnote%
 \ifx\CLASSOPTIONpt\@FMFPptsizenine\else%
  \typeout{** ATTENTION: Technotes are normally 9pt documents.}%
 \fi%
\fi


% V1.7
% Improved \textunderscore to provide a much better fake _ when used with
% OT1 encoding. Under OT1, detect use of pcr or cmtt \ttfamily and use
% available true _ glyph for those two typewriter fonts.
\def\@FMFPstringptm{ptm} % Times Roman family
\def\@FMFPstringppl{ppl} % Palatino Roman family
\def\@FMFPstringphv{phv} % Helvetica Sans Serif family
\def\@FMFPstringpcr{pcr} % Courier typewriter family
\def\@FMFPstringcmtt{cmtt} % Computer Modern typewriter family
\DeclareTextCommandDefault{\textunderscore}{\leavevmode
\ifx\f@family\@FMFPstringpcr\string_\else
\ifx\f@family\@FMFPstringcmtt\string_\else
\ifx\f@family\@FMFPstringptm\kern 0em\vbox{\hrule\@width 0.5em\@height 0.5pt\kern -0.3ex}\else
\ifx\f@family\@FMFPstringppl\kern 0em\vbox{\hrule\@width 0.5em\@height 0.5pt\kern -0.3ex}\else
\ifx\f@family\@FMFPstringphv\kern -0.03em\vbox{\hrule\@width 0.62em\@height 0.52pt\kern -0.33ex}\kern -0.03em\else
\kern 0.09em\vbox{\hrule\@width 0.6em\@height 0.44pt\kern -0.63pt\kern -0.42ex}\kern 0.09em\fi\fi\fi\fi\fi\relax}





% set the default \baselinestretch
\def\baselinestretch{1.0}
\ifCLASSOPTIONdraftcls
  \def\baselinestretch{1.5}% default baselinestretch for draft modes
\fi 


% process CLASSINPUT baselinestretch
\ifx\CLASSINPUTbaselinestretch\@FMFPundefined
\else
  \edef\baselinestretch{\CLASSINPUTbaselinestretch} % user CLASSINPUT override
  \typeout{** ATTENTION: Overriding \string\baselinestretch\space to
           \baselinestretch\space via \string\CLASSINPUT.}
\fi

\normalsize % make \baselinestretch take affect




% store the normalsize baselineskip
\newdimen\CLASSINFOnormalsizebaselineskip
\CLASSINFOnormalsizebaselineskip=\baselineskip\relax
% and the normalsize unity (baselinestretch=1) baselineskip
% we could save a register by giving the user access to
% \@FMFPnormalsizeunitybaselineskip. However, let's protect
% its read only internal status
\newdimen\CLASSINFOnormalsizeunitybaselineskip
\CLASSINFOnormalsizeunitybaselineskip=\@FMFPnormalsizeunitybaselineskip\relax
% store the nominal value of jot
\newdimen\FMFPnormaljot
\FMFPnormaljot=0.25\baselineskip\relax

% set \jot
\jot=\FMFPnormaljot\relax




% V1.6, we are now going to fine tune the interword spacing
% The default interword glue for Times under TeX appears to use a
% nominal interword spacing of 25% (relative to the font size, i.e., 1em)
% a maximum of 40% and a minimum of 19%.
% For example, 10pt text uses an interword glue of:
% 
% 2.5pt plus 1.49998pt minus 0.59998pt
% 
% However, FMFP allows for a more generous range which reduces the need
% for hyphenation, especially for two column text. Furthermore, FMFP
% tends to use a little bit more nominal space between the words.
% FMFP's interword spacing percentages appear to be:
% 35% nominal
% 23% minimum
% 50% maximum
% (They may even be using a tad more for the largest fonts such as 24pt.)
% 
% for bold text, FMFP increases the spacing a little more:
% 37.5% nominal
% 23% minimum
% 55% maximum

% here are the interword spacing ratios we'll use
% for medium (normal weight)
\def\@FMFPinterspaceratioM{0.35}
\def\@FMFPinterspaceMINratioM{0.23}
\def\@FMFPinterspaceMAXratioM{0.50}

% for bold
\def\@FMFPinterspaceratioB{0.375}
\def\@FMFPinterspaceMINratioB{0.23}
\def\@FMFPinterspaceMAXratioB{0.55}


% command to revise the interword spacing for the current font under TeX:
% \fontdimen2 = nominal interword space
% \fontdimen3 = interword stretch
% \fontdimen4 = interword shrink
% since all changes to the \fontdimen are global, we can enclose these commands
% in braces to confine any font attribute or length changes
\def\@@@FMFPsetfontdimens#1#2#3{{%
\setlength{\@FMFPtrantmpdimenB}{\f@size pt}% grab the font size in pt, could use 1em instead.
\setlength{\@FMFPtrantmpdimenA}{#1\@FMFPtrantmpdimenB}%
\fontdimen2\font=\@FMFPtrantmpdimenA\relax
\addtolength{\@FMFPtrantmpdimenA}{-#2\@FMFPtrantmpdimenB}%
\fontdimen3\font=-\@FMFPtrantmpdimenA\relax
\setlength{\@FMFPtrantmpdimenA}{#1\@FMFPtrantmpdimenB}%
\addtolength{\@FMFPtrantmpdimenA}{-#3\@FMFPtrantmpdimenB}%
\fontdimen4\font=\@FMFPtrantmpdimenA\relax}}

% revise the interword spacing for each font weight
\def\@@FMFPsetfontdimens{{%
\mdseries
\@@@FMFPsetfontdimens{\@FMFPinterspaceratioM}{\@FMFPinterspaceMAXratioM}{\@FMFPinterspaceMINratioM}%
\bfseries
\@@@FMFPsetfontdimens{\@FMFPinterspaceratioB}{\@FMFPinterspaceMAXratioB}{\@FMFPinterspaceMINratioB}%
}}

% revise the interword spacing for each font shape
% \slshape is not often used for FMFP work and is not altered here. The \scshape caps are
% already a tad too large in the free LaTeX fonts (as compared to what FMFP uses) so we
% won't alter these either.
\def\@FMFPsetfontdimens{{%
\normalfont
\@@FMFPsetfontdimens
\normalfont\itshape
\@@FMFPsetfontdimens
}}

% command to revise the interword spacing for each font size (and shape
% and weight). Only the \rmfamily is done here as \ttfamily uses a 
% fixed spacing and \sffamily is not used as the main text of FMFP papers.
\def\@FMFPtunefonts{{\selectfont\rmfamily
\tiny\@FMFPsetfontdimens
\scriptsize\@FMFPsetfontdimens
\footnotesize\@FMFPsetfontdimens
\small\@FMFPsetfontdimens
\normalsize\@FMFPsetfontdimens
\sublargesize\@FMFPsetfontdimens
\large\@FMFPsetfontdimens
\LARGE\@FMFPsetfontdimens
\huge\@FMFPsetfontdimens
\Huge\@FMFPsetfontdimens}}

% if the nofonttune class option is not given, revise the interword spacing
% now - in case FMFPtran makes any default length measurements, and make
% sure all the default fonts are loaded
\ifCLASSOPTIONnofonttune\else
\@FMFPtunefonts
\fi

% and again at the start of the document in case the user loaded different fonts
\AtBeginDocument{\ifCLASSOPTIONnofonttune\else\@FMFPtunefonts\fi}



% V1.6 
% LaTeX is a little to quick to use hyphenations
% So, we increase the penalty for their use and raise
% the badness level that triggers an underfull hbox
% warning. The author may still have to tweak things,
% but the appearance will be much better "right out
% of the box" than that under V1.5 and prior.
% TeX default is 50
\hyphenpenalty=750
% If we didn't adjust the interword spacing, 2200 might be better.
% The TeX default is 1000
\hbadness=1350
% FMFP does not use extra spacing after punctuation
\frenchspacing

% V1.7 increase this a tad to discourage equation breaks
\binoppenalty=1000 % default 700
\relpenalty=800     % default 500


% margin note stuff
\marginparsep      10pt
\marginparwidth    20pt
\marginparpush     25pt


% if things get too close, go ahead and let them touch
\lineskip            0pt
\normallineskip      0pt
\lineskiplimit       0pt
\normallineskiplimit 0pt

% The distance from the lower edge of the text body to the
% footline
\footskip 0.4in

% normally zero, should be relative to font height.
% put in a little rubber to help stop some bad breaks (underfull vboxes)
\parskip 0ex plus 0.2ex minus 0.1ex
\ifCLASSOPTIONconference
\parskip 0.0in %Changed by AJC 01.05.2018
%6pt plus 2pt minus 1pt
\fi

\parindent    1.0em
\ifCLASSOPTIONconference
\parindent 14.45pt
\fi

\topmargin    -49.0pt
\headheight   12pt
\headsep      0.25in

% use the normal font baselineskip
% so that \topskip is unaffected by changes in \baselinestretch
\topskip=\@FMFPnormalsizeunitybaselineskip
\textheight       58pc  % 9.63in, 696pt
% Tweak textheight to a perfect integer number of lines/page.
% The normal baselineskip for each document point size is used 
% to determine these values.
\ifx\CLASSOPTIONpt\@FMFPptsizenine\textheight=63\@FMFPnormalsizeunitybaselineskip\fi      % 63 lines/page
\ifx\CLASSOPTIONpt\@FMFPptsizeten\textheight=58\@FMFPnormalsizeunitybaselineskip\fi       % 58 lines/page
\ifx\CLASSOPTIONpt\@FMFPptsizeeleven\textheight=52\@FMFPnormalsizeunitybaselineskip\fi    % 52 lines/page
\ifx\CLASSOPTIONpt\@FMFPptsizetwelve\textheight=50\@FMFPnormalsizeunitybaselineskip\fi    % 50 lines/page


\columnsep       1.5pc
\textwidth       184.2mm
 

% the default side margins are equal
\if@FMFPusingAfourpaper 
\oddsidemargin        14.32mm
\evensidemargin       14.32mm
\else
\oddsidemargin        0.680in
\evensidemargin       0.680in
\fi
% compensate for LaTeX's 1in offset
\addtolength{\oddsidemargin}{-1in}
\addtolength{\evensidemargin}{-1in}



% adjust margins for conference mode
\ifCLASSOPTIONconference
 \topmargin        -0.25in
 % we retain the reserved, but unused space for headers
 \addtolength{\topmargin}{-\headheight}
 \addtolength{\topmargin}{-\headsep}
 \textheight        9.25in % The standard for conferences (668.4975pt)
 % Tweak textheight to a perfect integer number of lines/page.
 \ifx\CLASSOPTIONpt\@FMFPptsizenine\textheight=61\@FMFPnormalsizeunitybaselineskip\fi      % 61 lines/page
 \ifx\CLASSOPTIONpt\@FMFPptsizeten\textheight=62\@FMFPnormalsizeunitybaselineskip\fi       % 62 lines/page
 \ifx\CLASSOPTIONpt\@FMFPptsizeeleven\textheight=50\@FMFPnormalsizeunitybaselineskip\fi    % 50 lines/page
 \ifx\CLASSOPTIONpt\@FMFPptsizetwelve\textheight=48\@FMFPnormalsizeunitybaselineskip\fi    % 48 lines/page
\fi


% compsoc conference
\ifCLASSOPTIONcompsoc
\ifCLASSOPTIONconference
 % compsoc conference use a larger value for columnsep
 \columnsep 0.375in
 % compsoc conferences want 1in top margin, 1.125in bottom margin
 \topmargin        0in
 \addtolength{\topmargin}{-6pt}% we tweak this a tad to better comply with top of line stuff
 % we retain the reserved, but unused space for headers
 \addtolength{\topmargin}{-\headheight}
 \addtolength{\topmargin}{-\headsep}
 \textheight        8.875in % (641.39625pt)
 % Tweak textheight to a perfect integer number of lines/page.
 \ifx\CLASSOPTIONpt\@FMFPptsizenine\textheight=58\@FMFPnormalsizeunitybaselineskip\fi      % 58 lines/page
 \ifx\CLASSOPTIONpt\@FMFPptsizeten\textheight=53\@FMFPnormalsizeunitybaselineskip\fi       % 53 lines/page
 \ifx\CLASSOPTIONpt\@FMFPptsizeeleven\textheight=48\@FMFPnormalsizeunitybaselineskip\fi    % 48 lines/page
 \ifx\CLASSOPTIONpt\@FMFPptsizetwelve\textheight=46\@FMFPnormalsizeunitybaselineskip\fi    % 46 lines/page 
 \textwidth 6.5in
 % the default side margins are equal
 \if@FMFPusingAfourpaper 
  \oddsidemargin        22.45mm
  \evensidemargin       22.45mm
 \else
  \oddsidemargin        1in
  \evensidemargin       1in
 \fi
 % compensate for LaTeX's 1in offset
 \addtolength{\oddsidemargin}{-1in}
 \addtolength{\evensidemargin}{-1in}
\fi\fi



% draft mode settings override that of all other modes
% provides a nice 1in margin all around the paper and extra
% space between the lines for editor's comments
\ifCLASSOPTIONdraftcls 
  % want 1in from top of paper to text
  \setlength{\topmargin}{-\headsep}%
  \addtolength{\topmargin}{-\headheight}%
  % we want 1in side margins regardless of paper type
  \oddsidemargin      0in
  \evensidemargin     0in
  % set the text width
  \setlength{\textwidth}{\paperwidth}%
  \addtolength{\textwidth}{-2.0in}%
  \setlength{\textheight}{\paperheight}%
  \addtolength{\textheight}{-2.0in}%
  % digitize textheight to be an integer number of lines.
  % this may cause the bottom margin to be off a tad
  \addtolength{\textheight}{-1\topskip}%
  \divide\textheight  by \baselineskip%
  \multiply\textheight  by \baselineskip%
  \addtolength{\textheight}{\topskip}%
\fi



% process CLASSINPUT inner/outer margin
% if inner margin defined, but outer margin not, set outer to inner.
\ifx\CLASSINPUTinnersidemargin\@FMFPundefined
\else
  \ifx\CLASSINPUToutersidemargin\@FMFPundefined
    \edef\CLASSINPUToutersidemargin{\CLASSINPUTinnersidemargin}
  \fi
\fi

\ifx\CLASSINPUToutersidemargin\@FMFPundefined
\else
  % if outer margin defined, but inner margin not, set inner to outer.
  \ifx\CLASSINPUTinnersidemargin\@FMFPundefined
    \edef\CLASSINPUTinnersidemargin{\CLASSINPUToutersidemargin}
  \fi
  \setlength{\oddsidemargin}{\CLASSINPUTinnersidemargin}
  \ifCLASSOPTIONtwoside
    \setlength{\evensidemargin}{\CLASSINPUToutersidemargin}
  \else
    \setlength{\evensidemargin}{\CLASSINPUTinnersidemargin}
  \fi
  \addtolength{\oddsidemargin}{-1in}
  \addtolength{\evensidemargin}{-1in}
  \setlength{\textwidth}{\paperwidth}
  \addtolength{\textwidth}{-\CLASSINPUTinnersidemargin}
  \addtolength{\textwidth}{-\CLASSINPUToutersidemargin}
  \typeout{** ATTENTION: Overriding inner side margin to \CLASSINPUTinnersidemargin\space and 
           outer side margin to \CLASSINPUToutersidemargin\space via \string\CLASSINPUT.}
\fi



% process CLASSINPUT top/bottom text margin
% if toptext margin defined, but bottomtext margin not, set bottomtext to toptext margin
\ifx\CLASSINPUTtoptextmargin\@FMFPundefined
\else
  \ifx\CLASSINPUTbottomtextmargin\@FMFPundefined
    \edef\CLASSINPUTbottomtextmargin{\CLASSINPUTtoptextmargin}
  \fi
\fi

\ifx\CLASSINPUTbottomtextmargin\@FMFPundefined
\else
  % if bottomtext margin defined, but toptext margin not, set toptext to bottomtext margin
  \ifx\CLASSINPUTtoptextmargin\@FMFPundefined
    \edef\CLASSINPUTtoptextmargin{\CLASSINPUTbottomtextmargin}
  \fi
  \setlength{\topmargin}{\CLASSINPUTtoptextmargin}
  \addtolength{\topmargin}{-1in}
  \addtolength{\topmargin}{-\headheight}
  \addtolength{\topmargin}{-\headsep}
  \setlength{\textheight}{\paperheight}
  \addtolength{\textheight}{-\CLASSINPUTtoptextmargin}
  \addtolength{\textheight}{-\CLASSINPUTbottomtextmargin}
  % in the default format we use the normal baselineskip as topskip
  % we only need 0.7 of this to clear typical top text and we need
  % an extra 0.3 spacing at the bottom for descenders. This will
  % correct for both.
  \addtolength{\topmargin}{-0.3\@FMFPnormalsizeunitybaselineskip}
  \typeout{** ATTENTION: Overriding top text margin to \CLASSINPUTtoptextmargin\space and 
           bottom text margin to \CLASSINPUTbottomtextmargin\space via \string\CLASSINPUT.}
\fi







% LIST SPACING CONTROLS

% Controls the amount of EXTRA spacing
% above and below \trivlist 
% Both \list and IED lists override this.
% However, \trivlist will use this as will most
% things built from \trivlist like the \center
% environment.
\topsep           0.5\baselineskip

% Controls the additional spacing around lists preceded
% or followed by blank lines. FMFP does not increase
% spacing before or after paragraphs so it is set to zero.
% \z@ is the same as zero, but faster.
\partopsep          \z@

% Controls the spacing between paragraphs in lists. 
% FMFP does not increase spacing before or after paragraphs
% so this is also zero. 
% With FMFPtran.cls, global changes to
% this value DO affect lists (but not IED lists).
\parsep             \z@

% Controls the extra spacing between list items. 
% FMFP does not put extra spacing between items.
% With FMFPtran.cls, global changes to this value DO affect
% lists (but not IED lists).
\itemsep            \z@

% \itemindent is the amount to indent the FIRST line of a list
% item. It is auto set to zero within the \list environment. To alter
% it, you have to do so when you call the \list.
% However, FMFP uses this for the theorem environment
% There is an alternative value for this near \leftmargini below
\itemindent         -1em

% \leftmargin, the spacing from the left margin of the main text to
% the left of the main body of a list item is set by \list.
% Hence this statement does nothing for lists.
% But, quote and verse do use it for indention.
\leftmargin         2em

% we retain this stuff from the older FMFPtran.cls so that \list
% will work the same way as before. However, itemize, enumerate and
% description (IED) could care less about what these are as they
% all are overridden.
\leftmargini        2em
%\itemindent         2em  % Alternative values: sometimes used.
%\leftmargini        0em
\leftmarginii       1em
\leftmarginiii    1.5em
\leftmarginiv     1.5em
\leftmarginv      1.0em
\leftmarginvi     1.0em
\labelsep         0.5em 
\labelwidth         \z@


% The old FMFPtran.cls behavior of \list is retained.
% However, the new V1.3 IED list environments override all the
% @list stuff (\@listX is called within \list for the
% appropriate level just before the user's list_decl is called). 
% \topsep is now 2pt as FMFP puts a little extra space around
% lists - used by those non-IED macros that depend on \list.
% Note that \parsep and \itemsep are not redefined as in 
% the sizexx.clo \@listX (which article.cls uses) so global changes
% of these values DO affect \list
% 
\def\@listi{\leftmargin\leftmargini \topsep 2pt plus 1pt minus 1pt}
\let\@listI\@listi
\def\@listii{\leftmargin\leftmarginii\labelwidth\leftmarginii%
    \advance\labelwidth-\labelsep \topsep 2pt}
\def\@listiii{\leftmargin\leftmarginiii\labelwidth\leftmarginiii%
    \advance\labelwidth-\labelsep \topsep 2pt}
\def\@listiv{\leftmargin\leftmarginiv\labelwidth\leftmarginiv%
    \advance\labelwidth-\labelsep \topsep 2pt}
\def\@listv{\leftmargin\leftmarginv\labelwidth\leftmarginv%
    \advance\labelwidth-\labelsep \topsep 2pt}
\def\@listvi{\leftmargin\leftmarginvi\labelwidth\leftmarginvi%
    \advance\labelwidth-\labelsep \topsep 2pt}


% FMFP uses 5) not 5.
\def\labelenumi{\theenumi)}     \def\theenumi{\arabic{enumi}}

% FMFP uses a) not (a)
\def\labelenumii{\theenumii)}  \def\theenumii{\alph{enumii}}

% FMFP uses iii) not iii.
\def\labelenumiii{\theenumiii)} \def\theenumiii{\roman{enumiii}}

% FMFP uses A) not A.
\def\labelenumiv{\theenumiv)}   \def\theenumiv{\Alph{enumiv}}

% exactly the same as in article.cls
\def\p@enumii{\theenumi}
\def\p@enumiii{\theenumi(\theenumii)}
\def\p@enumiv{\p@enumiii\theenumiii}

% itemized list label styles
\def\labelitemi{$\bullet$}
\def\labelitemii{$\circ$}
\def\labelitemiii{\vrule height 0.8ex depth -0.2ex width 0.6ex}
\def\labelitemiv{$\ast$}



% **** V1.3 ENHANCEMENTS ****
% Itemize, Enumerate and Description (IED) List Controls
% ***************************
% 
% 
% FMFP seems to use at least two different values by
% which ITEMIZED list labels are indented to the right
% For The Journal of Lightwave Technology (JLT) and The Journal
% on Selected Areas in Communications (JSAC), they tend to use
% an indention equal to \parindent. For Transactions on Communications
% they tend to indent ITEMIZED lists a little more--- 1.3\parindent.
% We'll provide both values here for you so that you can choose 
% which one you like in your document using a command such as:
% setlength{\FMFPilabelindent}{\FMFPilabelindentB}
\newdimen\FMFPilabelindentA
\FMFPilabelindentA \parindent

\newdimen\FMFPilabelindentB
\FMFPilabelindentB 1.3\parindent
% However, we'll default to using \parindent
% which makes more sense to me
\newdimen\FMFPilabelindent
\FMFPilabelindent \FMFPilabelindentA


% This controls the default amount the enumerated list labels
% are indented to the right.
% Normally, this is the same as the paragraph indention
\newdimen\FMFPelabelindent
\FMFPelabelindent \parindent

% This controls the default amount the description list labels
% are indented to the right.
% Normally, this is the same as the paragraph indention
\newdimen\FMFPdlabelindent
\FMFPdlabelindent \parindent

% This is the value actually used within the IED lists.
% The IED environments automatically set its value to
% one of the three values above, so global changes do 
% not have any effect
\newdimen\FMFPlabelindent
\FMFPlabelindent \parindent

% The actual amount labels will be indented is
% \FMFPlabelindent multiplied by the factor below
% corresponding to the level of nesting depth
% This provides a means by which the user can
% alter the effective \FMFPlabelindent for deeper
% levels
% There may not be such a thing as correct "standard FMFP"
% values. What FMFP actually does may depend on the specific
% circumstances.
% The first list level almost always has full indention.
% The second levels I've seen have only 75% of the normal indentation
% Three level or greater nestings are very rare. I am guessing
% that they don't use any indentation.
\def\FMFPlabelindentfactori{1.0}   % almost always one
\def\FMFPlabelindentfactorii{0.75} % 0.0 or 1.0 may be used in some cases
\def\FMFPlabelindentfactoriii{0.0} % 0.75? 0.5? 0.0?
\def\FMFPlabelindentfactoriv{0.0}
\def\FMFPlabelindentfactorv{0.0}
\def\FMFPlabelindentfactorvi{0.0}

% value actually used within IED lists, it is auto
% set to one of the 6 values above
% global changes here have no effect
\def\FMFPlabelindentfactor{1.0}

% This controls the default spacing between the end of the IED
% list labels and the list text, when normal text is used for
% the labels.
\newdimen\FMFPiednormlabelsep
\FMFPiednormlabelsep \parindent

% This controls the default spacing between the end of the IED
% list labels and the list text, when math symbols are used for
% the labels (nomenclature lists). FMFP usually increases the 
% spacing in these cases
\newdimen\FMFPiedmathlabelsep
\FMFPiedmathlabelsep 1.2em

% This controls the extra vertical separation put above and
% below each IED list. FMFP usually puts a little extra spacing
% around each list. However, this spacing is barely noticeable.
\newskip\FMFPiedtopsep
\FMFPiedtopsep 2pt plus 1pt minus 1pt


% This command is executed within each IED list environment
% at the beginning of the list. You can use this to set the 
% parameters for some/all your IED list(s) without disturbing 
% global parameters that affect things other than lists.
% i.e., renewcommand{\FMFPiedlistdecl}{\setlength{\labelsep}{5em}}
% will alter the \labelsep for the next list(s) until 
% \FMFPiedlistdecl is redefined. 
\def\FMFPiedlistdecl{\relax}

% This command provides an easy way to set \leftmargin based
% on the \labelwidth, \labelsep and the argument \FMFPlabelindent
% Usage: \FMFPcalcleftmargin{width-to-indent-the-label}
% output is in the \leftmargin variable, i.e., effectively:
% \leftmargin = argument + \labelwidth + \labelsep
% Note controlled spacing here, shield end of lines with %
\def\FMFPcalcleftmargin#1{\setlength{\leftmargin}{#1}%
\addtolength{\leftmargin}{\labelwidth}%
\addtolength{\leftmargin}{\labelsep}}

% This command provides an easy way to set \labelwidth to the
% width of the given text. It is the same as
% \settowidth{\labelwidth}{label-text}
% and useful as a shorter alternative.
% Typically used to set \labelwidth to be the width
% of the longest label in the list
\def\FMFPsetlabelwidth#1{\settowidth{\labelwidth}{#1}}

% When this command is executed, IED lists will use the 
% FMFPiedmathlabelsep label separation rather than the normal
% spacing. To have an effect, this command must be executed via
% the \FMFPiedlistdecl or within the option of the IED list
% environments.
\def\FMFPusemathlabelsep{\setlength{\labelsep}{\FMFPiedmathlabelsep}}

% A flag which controls whether the IED lists automatically
% calculate \leftmargin from \FMFPlabelindent, \labelwidth and \labelsep
% Useful if you want to specify your own \leftmargin
% This flag must be set (\FMFPnocalcleftmargintrue or \FMFPnocalcleftmarginfalse) 
% via the \FMFPiedlistdecl or within the option of the IED list
% environments to have an effect.
\newif\ifFMFPnocalcleftmargin
\FMFPnocalcleftmarginfalse

% A flag which controls whether \FMFPlabelindent is multiplied by
% the \FMFPlabelindentfactor for each list level.
% This flag must be set via the \FMFPiedlistdecl or within the option 
% of the IED list environments to have an effect.
\newif\ifFMFPnolabelindentfactor
\FMFPnolabelindentfactorfalse


% internal variable to indicate type of IED label
% justification
% 0 - left; 1 - center; 2 - right
\def\@FMFPiedjustify{0}


% commands to allow the user to control IED
% label justifications. Use these commands within
% the IED environment option or in the \FMFPiedlistdecl
% Note that changing the normal list justifications
% is nonstandard and FMFP may not like it if you do so!
% I include these commands as they may be helpful to
% those who are using these enhanced list controls for
% other non-FMFP related LaTeX work.
% itemize and enumerate automatically default to right
% justification, description defaults to left.
\def\FMFPiedlabeljustifyl{\def\@FMFPiedjustify{0}}%left
\def\FMFPiedlabeljustifyc{\def\@FMFPiedjustify{1}}%center
\def\FMFPiedlabeljustifyr{\def\@FMFPiedjustify{2}}%right




% commands to save to and restore from the list parameter copies
% this allows us to set all the list parameters within
% the list_decl and prevent \list (and its \@list) 
% from overriding any of our parameters
% V1.6 use \edefs instead of dimen's to conserve dimen registers
% Note controlled spacing here, shield end of lines with %
\def\@FMFPsavelistparams{\edef\@FMFPiedtopsep{\the\topsep}%
\edef\@FMFPiedlabelwidth{\the\labelwidth}%
\edef\@FMFPiedlabelsep{\the\labelsep}%
\edef\@FMFPiedleftmargin{\the\leftmargin}%
\edef\@FMFPiedpartopsep{\the\partopsep}%
\edef\@FMFPiedparsep{\the\parsep}%
\edef\@FMFPieditemsep{\the\itemsep}%
\edef\@FMFPiedrightmargin{\the\rightmargin}%
\edef\@FMFPiedlistparindent{\the\listparindent}%
\edef\@FMFPieditemindent{\the\itemindent}}

% Note controlled spacing here
\def\@FMFPrestorelistparams{\topsep\@FMFPiedtopsep\relax%
\labelwidth\@FMFPiedlabelwidth\relax%
\labelsep\@FMFPiedlabelsep\relax%
\leftmargin\@FMFPiedleftmargin\relax%
\partopsep\@FMFPiedpartopsep\relax%
\parsep\@FMFPiedparsep\relax%
\itemsep\@FMFPieditemsep\relax%
\rightmargin\@FMFPiedrightmargin\relax%
\listparindent\@FMFPiedlistparindent\relax%
\itemindent\@FMFPieditemindent\relax}


% v1.6b provide original LaTeX IED list environments
% note that latex.ltx defines \itemize and \enumerate, but not \description
% which must be created by the base classes
% save original LaTeX itemize and enumerate
\let\LaTeXitemize\itemize
\let\endLaTeXitemize\enditemize
\let\LaTeXenumerate\enumerate
\let\endLaTeXenumerate\endenumerate

% provide original LaTeX description environment from article.cls
\newenvironment{LaTeXdescription}
               {\list{}{\labelwidth\z@ \itemindent-\leftmargin
                        \let\makelabel\descriptionlabel}}
               {\endlist}
\newcommand*\descriptionlabel[1]{\hspace\labelsep
                                 \normalfont\bfseries #1}


% override LaTeX's default IED lists
\def\itemize{\@FMFPitemize}
\def\enditemize{\@endFMFPitemize}
\def\enumerate{\@FMFPenumerate}
\def\endenumerate{\@endFMFPenumerate}
\def\description{\@FMFPdescription}
\def\enddescription{\@endFMFPdescription}

% provide the user with aliases - may help those using packages that
% override itemize, enumerate, or description
\def\FMFPitemize{\@FMFPitemize}
\def\endFMFPitemize{\@endFMFPitemize}
\def\FMFPenumerate{\@FMFPenumerate}
\def\endFMFPenumerate{\@endFMFPenumerate}
\def\FMFPdescription{\@FMFPdescription}
\def\endFMFPdescription{\@endFMFPdescription}


% V1.6 we want to keep the FMFPtran IED list definitions as our own internal
% commands so they are protected against redefinition
\def\@FMFPitemize{\@ifnextchar[{\@@FMFPitemize}{\@@FMFPitemize[\relax]}}
\def\@FMFPenumerate{\@ifnextchar[{\@@FMFPenumerate}{\@@FMFPenumerate[\relax]}}
\def\@FMFPdescription{\@ifnextchar[{\@@FMFPdescription}{\@@FMFPdescription[\relax]}}
\def\@endFMFPitemize{\endlist}
\def\@endFMFPenumerate{\endlist}
\def\@endFMFPdescription{\endlist}


% DO NOT ALLOW BLANK LINES TO BE IN THESE IED ENVIRONMENTS
% AS THIS WILL FORCE NEW PARAGRAPHS AFTER THE IED LISTS
% FMFPtran itemized list MDS 1/2001
% Note controlled spacing here, shield end of lines with %
\def\@@FMFPitemize[#1]{%
                \ifnum\@itemdepth>3\relax\@toodeep\else%
                \ifnum\@listdepth>5\relax\@toodeep\else%
                \advance\@itemdepth\@ne%
                \edef\@itemitem{labelitem\romannumeral\the\@itemdepth}%
                % get the labelindentfactor for this level
                \advance\@listdepth\@ne% we need to know what the level WILL be
                \edef\FMFPlabelindentfactor{\csname FMFPlabelindentfactor\romannumeral\the\@listdepth\endcsname}%
                \advance\@listdepth-\@ne% undo our increment
                \def\@FMFPiedjustify{2}% right justified labels are default
                % set other defaults
                \FMFPnocalcleftmarginfalse%
                \FMFPnolabelindentfactorfalse%
                \topsep\FMFPiedtopsep%
                \FMFPlabelindent\FMFPilabelindent%
                \labelsep\FMFPiednormlabelsep%
                \partopsep 0ex%
                \parsep 0ex%
                \itemsep \parskip%
                \rightmargin 0em%
                \listparindent 0em%
                \itemindent 0em%
                % calculate the label width
                % the user can override this later if
                % they specified a \labelwidth
                \settowidth{\labelwidth}{\csname labelitem\romannumeral\the\@itemdepth\endcsname}%
                \@FMFPsavelistparams% save our list parameters
                \list{\csname\@itemitem\endcsname}{%
                \@FMFPrestorelistparams% override any list{} changes
                                       % to our globals
                \let\makelabel\@FMFPiedmakelabel% v1.6b setup \makelabel
                \FMFPiedlistdecl% let user alter parameters
                #1\relax%
                % If the user has requested not to use the
                % labelindent factor, don't revise \labelindent
                \ifFMFPnolabelindentfactor\relax%
                \else\FMFPlabelindent=\FMFPlabelindentfactor\labelindent%
                \fi%
                % Unless the user has requested otherwise,
                % calculate our left margin based
                % on \FMFPlabelindent, \labelwidth and
                % \labelsep
                \ifFMFPnocalcleftmargin\relax%
                \else\FMFPcalcleftmargin{\FMFPlabelindent}%
                \fi}\fi\fi}%


% DO NOT ALLOW BLANK LINES TO BE IN THESE IED ENVIRONMENTS
% AS THIS WILL FORCE NEW PARAGRAPHS AFTER THE IED LISTS
% FMFPtran enumerate list MDS 1/2001
% Note controlled spacing here, shield end of lines with %
\def\@@FMFPenumerate[#1]{%
                \ifnum\@enumdepth>3\relax\@toodeep\else%
                \ifnum\@listdepth>5\relax\@toodeep\else%
                \advance\@enumdepth\@ne%
                \edef\@enumctr{enum\romannumeral\the\@enumdepth}%
                % get the labelindentfactor for this level
                \advance\@listdepth\@ne% we need to know what the level WILL be
                \edef\FMFPlabelindentfactor{\csname FMFPlabelindentfactor\romannumeral\the\@listdepth\endcsname}%
                \advance\@listdepth-\@ne% undo our increment
                \def\@FMFPiedjustify{2}% right justified labels are default
                % set other defaults
                \FMFPnocalcleftmarginfalse%
                \FMFPnolabelindentfactorfalse%
                \topsep\FMFPiedtopsep%
                \FMFPlabelindent\FMFPelabelindent%
                \labelsep\FMFPiednormlabelsep%
                \partopsep 0ex%
                \parsep 0ex%
                \itemsep 0ex%
                \rightmargin 0em%
                \listparindent 0em%
                \itemindent 0em%
                % calculate the label width
                % We'll set it to the width suitable for all labels using
                % normalfont 1) to 9)
                % The user can override this later
                \settowidth{\labelwidth}{9)}%
                \@FMFPsavelistparams% save our list parameters
                \list{\csname label\@enumctr\endcsname}{\usecounter{\@enumctr}%
                \@FMFPrestorelistparams% override any list{} changes
                                       % to our globals
                \let\makelabel\@FMFPiedmakelabel% v1.6b setup \makelabel
                \FMFPiedlistdecl% let user alter parameters 
                #1\relax%
                % If the user has requested not to use the
                % FMFPlabelindent factor, don't revise \FMFPlabelindent
                \ifFMFPnolabelindentfactor\relax%
                \else\FMFPlabelindent=\FMFPlabelindentfactor\FMFPlabelindent%
                \fi%
                % Unless the user has requested otherwise,
                % calculate our left margin based
                % on \FMFPlabelindent, \labelwidth and
                % \labelsep
                \ifFMFPnocalcleftmargin\relax%
                \else\FMFPcalcleftmargin{\FMFPlabelindent}%
                \fi}\fi\fi}%


% DO NOT ALLOW BLANK LINES TO BE IN THESE IED ENVIRONMENTS
% AS THIS WILL FORCE NEW PARAGRAPHS AFTER THE IED LISTS
% FMFPtran description list MDS 1/2001
% Note controlled spacing here, shield end of lines with %
\def\@@FMFPdescription[#1]{%
                \ifnum\@listdepth>5\relax\@toodeep\else%
                % get the labelindentfactor for this level
                \advance\@listdepth\@ne% we need to know what the level WILL be
                \edef\FMFPlabelindentfactor{\csname FMFPlabelindentfactor\romannumeral\the\@listdepth\endcsname}%
                \advance\@listdepth-\@ne% undo our increment
                \def\@FMFPiedjustify{0}% left justified labels are default
                % set other defaults
                \FMFPnocalcleftmarginfalse%
                \FMFPnolabelindentfactorfalse%
                \topsep\FMFPiedtopsep% 
                \FMFPlabelindent\FMFPdlabelindent%
                % assume normal labelsep
                \labelsep\FMFPiednormlabelsep%
                \partopsep 0ex%
                \parsep 0ex%
                \itemsep 0ex%
                \rightmargin 0em%
                \listparindent 0em%
                \itemindent 0em%
                % Bogus label width in case the user forgets
                % to set it.
                % TIP: If you want to see what a variable's width is you
                % can use the TeX command \showthe\width-variable to 
                % display it on the screen during compilation 
                % (This might be helpful to know when you need to find out
                % which label is the widest)
                \settowidth{\labelwidth}{Hello}%
                \@FMFPsavelistparams% save our list parameters
                \list{}{\@FMFPrestorelistparams% override any list{} changes
                                               % to our globals
                \let\makelabel\@FMFPiedmakelabel% v1.6b setup \makelabel
                \FMFPiedlistdecl% let user alter parameters 
                #1\relax%
                % If the user has requested not to use the
                % labelindent factor, don't revise \FMFPlabelindent
                \ifFMFPnolabelindentfactor\relax%
                \else\FMFPlabelindent=\FMFPlabelindentfactor\FMFPlabelindent%
                \fi%
                % Unless the user has requested otherwise,
                % calculate our left margin based
                % on \FMFPlabelindent, \labelwidth and
                % \labelsep
                \ifFMFPnocalcleftmargin\relax%
                \else\FMFPcalcleftmargin{\FMFPlabelindent}\relax%
                \fi}\fi}

% v1.6b we use one makelabel that does justification as needed.
\def\@FMFPiedmakelabel#1{\relax\if\@FMFPiedjustify 0\relax
\makebox[\labelwidth][l]{\normalfont #1}\else
\if\@FMFPiedjustify 1\relax
\makebox[\labelwidth][c]{\normalfont #1}\else
\makebox[\labelwidth][r]{\normalfont #1}\fi\fi}


% VERSE and QUOTE
% V1.7 define environments with newenvironment
\newenvironment{verse}{\let\\=\@centercr
    \list{}{\itemsep\z@ \itemindent -1.5em \listparindent \itemindent
    \rightmargin\leftmargin\advance\leftmargin 1.5em}\item\relax}
    {\endlist}
\newenvironment{quotation}{\list{}{\listparindent 1.5em \itemindent\listparindent
    \rightmargin\leftmargin \parsep 0pt plus 1pt}\item\relax}
    {\endlist}
\newenvironment{quote}{\list{}{\rightmargin\leftmargin}\item\relax}
    {\endlist}


% \titlepage
% provided only for backward compatibility. \maketitle is the correct
% way to create the title page. 
\newif\if@restonecol
\def\titlepage{\@restonecolfalse\if@twocolumn\@restonecoltrue\onecolumn
    \else \newpage \fi \thispagestyle{empty}\c@page\z@}
\def\endtitlepage{\if@restonecol\twocolumn \else \newpage \fi}

% standard values from article.cls
\arraycolsep     5pt
\arrayrulewidth .4pt
\doublerulesep   2pt

\tabcolsep       6pt
\tabbingsep      0.5em


%% FOOTNOTES
%
%\skip\footins 10pt plus 4pt minus 2pt
% V1.6 respond to changes in font size
% space added above the footnotes (if present)
\skip\footins 0.9\baselineskip  plus 0.4\baselineskip  minus 0.2\baselineskip

% V1.6, we need to make \footnotesep responsive to changes
% in \baselineskip or strange spacings will result when in
% draft mode. Here is a little LaTeX secret - \footnotesep
% determines the height of an invisible strut that is placed
% *above* the baseline of footnotes after the first. Since
% LaTeX considers the space for characters to be 0.7/baselineskip
% above the baseline and 0.3/baselineskip below it, we need to
% use 0.7/baselineskip as a \footnotesep to maintain equal spacing
% between all the lines of the footnotes. FMFP often uses a tad
% more, so use 0.8\baselineskip. This slightly larger value also helps
% the text to clear the footnote marks. Note that \thanks in FMFPtran
% uses its own value of \footnotesep which is set in \maketitle.
{\footnotesize
\global\footnotesep 0.8\baselineskip}

\def\unnumberedfootnote{\gdef\@thefnmark{\quad}\@footnotetext}

\skip\@mpfootins 0.3\baselineskip
\fboxsep = 3pt
\fboxrule = .4pt
% V1.6 use 1em, then use LaTeX2e's \@makefnmark
% Note that FMFP normally *left* aligns the footnote marks, so we don't need
% box resizing tricks here.
%\long\def\@makefnmark{\scriptsize\normalfont\@thefnmark}
\long\def\@makefntext#1{\parindent 1em\indent\hbox{\@makefnmark}#1}% V1.6 use 1em
\long\def\@maketablefntext#1{\raggedleft\leavevmode\hbox{\@makefnmark}#1}
% V1.7 compsoc does not use superscipts for footnote marks
\ifCLASSOPTIONcompsoc
\def\@FMFPcompsocmakefnmark{\hbox{\normalfont\@thefnmark.\ }}
\long\def\@makefntext#1{\parindent 1em\indent\hbox{\@FMFPcompsocmakefnmark}#1}
\fi

% FMFP does not use footnote rules. Or do they?
\def\footnoterule{\vskip-2pt \hrule height 0.6pt depth \z@ \vskip1.6pt\relax}
\toks@\expandafter{\@setminipage\let\footnoterule\relax\footnotesep\z@}
\edef\@setminipage{\the\toks@}

% V1.7 for compsoc, FMFP uses a footnote rule only for \thanks. We devise a "one-shot"
% system to implement this.
\newif\if@FMFPenableoneshotfootnoterule
\@FMFPenableoneshotfootnoterulefalse
\ifCLASSOPTIONcompsoc
\def\footnoterule{\relax\if@FMFPenableoneshotfootnoterule
\kern-5pt
\hbox to \columnwidth{\hfill\vrule width 0.5\columnwidth height 0.4pt\hfill}
\kern4.6pt
\global\@FMFPenableoneshotfootnoterulefalse
\else
\relax
\fi}
\fi

% V1.6 do not allow LaTeX to break a footnote across multiple pages
\interfootnotelinepenalty=10000

% V1.6 discourage breaks within equations
% Note that amsmath normally sets this to 10000,
% but LaTeX2e normally uses 100.
\interdisplaylinepenalty=2500

% default allows section depth up to /paragraph
\setcounter{secnumdepth}{4}

% technotes do not allow /paragraph
\ifCLASSOPTIONtechnote
   \setcounter{secnumdepth}{3}
\fi
% neither do compsoc conferences
\@FMFPcompsocconfonly{\setcounter{secnumdepth}{3}}


\newcounter{section}
\newcounter{subsection}[section]
\newcounter{subsubsection}[subsection]
\newcounter{paragraph}[subsubsection]

% used only by FMFPtran's FMFPeqnarray as other packages may
% have their own, different, implementations
\newcounter{FMFPsubequation}[equation]

% as shown when called by user from \ref, \label and in table of contents
\def\theequation{\arabic{equation}}                          % 1
\def\theFMFPsubequation{\theequation\alph{FMFPsubequation}}  % 1a (used only by FMFPtran's FMFPeqnarray)
\ifCLASSOPTIONcompsoc
% compsoc is all arabic
\def\thesection{\arabic{section}}                
\def\thesubsection{\thesection.\arabic{subsection}}
\def\thesubsubsection{\thesubsection.\arabic{subsubsection}}
\def\theparagraph{\thesubsubsection.\arabic{paragraph}}
\else
\def\thesection{\Roman{section}}                             % I
% V1.7, \mbox prevents breaks around - 
\def\thesubsection{\mbox{\thesection-\Alph{subsection}}}     % I-A
% V1.7 use I-A1 format used by FMFP rather than I-A.1
\def\thesubsubsection{\thesubsection\arabic{subsubsection}}  % I-A1
\def\theparagraph{\thesubsubsection\alph{paragraph}}         % I-A1a
\fi

% From Heiko Oberdiek. Because of the \mbox in \thesubsection, we need to
% tell hyperref to disable the \mbox command when making PDF bookmarks.
% This done already with hyperref.sty version 6.74o and later, but
% it will not hurt to do it here again for users of older versions.
\@ifundefined{pdfstringdefPreHook}{\let\pdfstringdefPreHook\@empty}{}%
\g@addto@macro\pdfstringdefPreHook{\let\mbox\relax}


% Main text forms (how shown in main text headings)
% V1.6, using \thesection in \thesectiondis allows changes
% in the former to automatically appear in the latter
\ifCLASSOPTIONcompsoc
  \ifCLASSOPTIONconference% compsoc conference
    \def\thesectiondis{\thesection.}
    \def\thesubsectiondis{\thesectiondis\arabic{subsection}.}
    \def\thesubsubsectiondis{\thesubsectiondis\arabic{subsubsection}.}
    \def\theparagraphdis{\thesubsubsectiondis\arabic{paragraph}.}
  \else% compsoc not conferencs
    \def\thesectiondis{\thesection}
    \def\thesubsectiondis{\thesectiondis.\arabic{subsection}}
    \def\thesubsubsectiondis{\thesubsectiondis.\arabic{subsubsection}}
    \def\theparagraphdis{\thesubsubsectiondis.\arabic{paragraph}}
  \fi
\else% not compsoc
  \def\thesectiondis{\thesection.}                   % I.
  \def\thesubsectiondis{\Alph{subsection}.}          % B.
  \def\thesubsubsectiondis{\arabic{subsubsection})}  % 3)
  \def\theparagraphdis{\alph{paragraph})}            % d)
\fi

% just like LaTeX2e's \@eqnnum
\def\theequationdis{{\normalfont \normalcolor (\theequation)}}% (1)
% FMFPsubequation used only by FMFPtran's FMFPeqnarray
\def\theFMFPsubequationdis{{\normalfont \normalcolor (\theFMFPsubequation)}}% (1a)
% redirect LaTeX2e's equation number display and all that depend on
% it, through FMFPtran's \theequationdis
\def\@eqnnum{\theequationdis}



% V1.7 provide string macros as article.cls does
\def\contentsname{Contents}
\def\listfigurename{List of Figures}
\def\listtablename{List of Tables}
\def\refname{References}
\def\indexname{Index}
\def\figurename{\textbf{Figure}} %Changed by AJC on 28.04.2018
\def\tablename{\textbf{Table}} %Changed by AJC on 28.04.2018
\@FMFPcompsocconfonly{\def\figurename{Figure}\def\tablename{Table}}
\def\partname{Part}
\def\appendixname{Appendix}
\def\abstractname{Abstract}
% FMFP specific names
\def\FMFPkeywordsname{Keywords}
\def\FMFPproofname{Proof}


% LIST OF FIGURES AND TABLES AND TABLE OF CONTENTS
%
\def\@pnumwidth{1.55em}
\def\@tocrmarg{2.55em}
\def\@dotsep{4.5}
\setcounter{tocdepth}{3}

% adjusted some spacings here so that section numbers will not easily 
% collide with the section titles. 
% VIII; VIII-A; and VIII-A.1 are usually the worst offenders.
% MDS 1/2001
\def\tableofcontents{\section*{\contentsname}\@starttoc{toc}}
\def\l@section#1#2{\addpenalty{\@secpenalty}\addvspace{1.0em plus 1pt}%
    \@tempdima 2.75em \begingroup \parindent \z@ \rightskip \@pnumwidth%
    \parfillskip-\@pnumwidth {\bfseries\leavevmode #1}\hfil\hbox to\@pnumwidth{\hss #2}\par%
    \endgroup}
% argument format #1:level, #2:labelindent,#3:labelsep
\def\l@subsection{\@dottedtocline{2}{2.75em}{3.75em}}
\def\l@subsubsection{\@dottedtocline{3}{6.5em}{4.5em}}
% must provide \l@ defs for ALL sublevels EVEN if tocdepth
% is such as they will not appear in the table of contents
% these defs are how TOC knows what level these things are!
\def\l@paragraph{\@dottedtocline{4}{6.5em}{5.5em}}
\def\l@subparagraph{\@dottedtocline{5}{6.5em}{6.5em}}
\def\listoffigures{\section*{\listfigurename}\@starttoc{lof}}
\def\l@figure{\@dottedtocline{1}{0em}{2.75em}}
\def\listoftables{\section*{\listtablename}\@starttoc{lot}}
\let\l@table\l@figure


%% Definitions for floats
%%
%% Normal Floats
\floatsep 1\baselineskip plus  0.2\baselineskip minus  0.2\baselineskip
\textfloatsep 1.7\baselineskip plus  0.2\baselineskip minus  0.4\baselineskip
\@fptop 0pt plus 1fil
\@fpsep 0.75\baselineskip plus 2fil 
\@fpbot 0pt plus 1fil
\def\topfraction{0.9}
\def\bottomfraction{0.4}
\def\floatpagefraction{0.8}
% V1.7, let top floats approach 90% of page
\def\textfraction{0.1}

%% Double Column Floats
\dblfloatsep 1\baselineskip plus  0.2\baselineskip minus  0.2\baselineskip

\dbltextfloatsep 1.7\baselineskip plus  0.2\baselineskip minus  0.4\baselineskip
% Note that it would be nice if the rubber here actually worked in LaTeX2e.
% There is a long standing limitation in LaTeX, first discovered (to the best
% of my knowledge) by Alan Jeffrey in 1992. LaTeX ignores the stretchable
% portion of \dbltextfloatsep, and as a result, double column figures can and
% do result in an non-integer number of lines in the main text columns with
% underfull vbox errors as a consequence. A post to comp.text.tex
% by Donald Arseneau confirms that this had not yet been fixed in 1998.
% FMFPtran V1.6 will fix this problem for you in the titles, but it doesn't
% protect you from other double floats. Happy vspace'ing.

\@dblfptop 0pt plus 1fil
\@dblfpsep 0.75\baselineskip plus 2fil
\@dblfpbot 0pt plus 1fil
\def\dbltopfraction{0.8}
\def\dblfloatpagefraction{0.8}
\setcounter{dbltopnumber}{4}

\intextsep 1\baselineskip plus 0.2\baselineskip minus  0.2\baselineskip
\setcounter{topnumber}{2}
\setcounter{bottomnumber}{2}
\setcounter{totalnumber}{4}



% article class provides these, we should too.
\newlength\abovecaptionskip
\newlength\belowcaptionskip
% but only \abovecaptionskip is used above figure captions and *below* table
% captions
\setlength\abovecaptionskip{0.65\baselineskip}
\setlength\belowcaptionskip{0.75\baselineskip}
% V1.6 create hooks in case the caption spacing ever needs to be
% overridden by a user
\def\@FMFPfigurecaptionsepspace{\vskip\abovecaptionskip\relax}%
\def\@FMFPtablecaptionsepspace{\vskip\belowcaptionskip\relax}%


% 1.6b revise caption system so that \@makecaption uses two arguments
% as with LaTeX2e. Otherwise, there will be problems when using hyperref.
\def\@FMFPtablestring{table}

\ifCLASSOPTIONcompsoc
% V1.7 compsoc \@makecaption
\ifCLASSOPTIONconference% compsoc conference
\long\def\@makecaption#1#2{%
% test if is a for a figure or table
\ifx\@captype\@FMFPtablestring%
% if a table, do table caption
\normalsize\begin{center}{\normalfont\sffamily\normalsize {#1.}~ #2}\end{center}%
\@FMFPtablecaptionsepspace
% if not a table, format it as a figure
\else
\@FMFPfigurecaptionsepspace
\setbox\@tempboxa\hbox{\normalfont\sffamily\normalsize {#1.}~ #2}%
\ifdim \wd\@tempboxa >\hsize%
% if caption is longer than a line, let it wrap around
\setbox\@tempboxa\hbox{\normalfont\sffamily\normalsize {#1.}~ }%
\parbox[t]{\hsize}{\normalfont\sffamily\normalsize \noindent\unhbox\@tempboxa#2}%
% if caption is shorter than a line, center
\else%
\hbox to\hsize{\normalfont\sffamily\normalsize\hfil\box\@tempboxa\hfil}%
\fi\fi}
\else% nonconference compsoc
\long\def\@makecaption#1#2{%
% test if is a for a figure or table
\ifx\@captype\@FMFPtablestring%
% if a table, do table caption
\normalsize\begin{center}{\normalfont\sffamily\normalsize #1}\\{\normalfont\sffamily\normalsize #2}\end{center}%
\@FMFPtablecaptionsepspace
% if not a table, format it as a figure
\else
\@FMFPfigurecaptionsepspace
\setbox\@tempboxa\hbox{\normalfont\sffamily\normalsize {#1.}~ #2}%
\ifdim \wd\@tempboxa >\hsize%
% if caption is longer than a line, let it wrap around
\setbox\@tempboxa\hbox{\normalfont\sffamily\normalsize {#1.}~ }%
\parbox[t]{\hsize}{\normalfont\sffamily\normalsize \noindent\unhbox\@tempboxa#2}%
% if caption is shorter than a line, left justify
\else%
\hbox to\hsize{\normalfont\sffamily\normalsize\box\@tempboxa\hfil}%
\fi\fi}
\fi

\else% traditional noncompsoc \@makecaption
\long\def\@makecaption#1#2{%
% test if is a for a figure or table
\ifx\@captype\@FMFPtablestring%
% if a table, do table caption
\footnotesize{\centering\normalfont\footnotesize#1.\qquad\scshape #2\par}%
\@FMFPtablecaptionsepspace
% if not a table, format it as a figure
\else
\@FMFPfigurecaptionsepspace
% 3/2001 use footnotesize, not small; use two nonbreaking spaces, not one
\setbox\@tempboxa\hbox{\normalfont\footnotesize {#1.}~~ #2}%
\ifdim \wd\@tempboxa >\hsize%
% if caption is longer than a line, let it wrap around
\setbox\@tempboxa\hbox{\normalfont\footnotesize {#1.}~~ }%
\parbox[t]{\hsize}{\normalfont\footnotesize\noindent\unhbox\@tempboxa#2}%
% if caption is shorter than a line, center if conference, left justify otherwise
\else%
\ifCLASSOPTIONconference \hbox to\hsize{\normalfont\footnotesize\box\@tempboxa\hfil}%
\else \hbox to\hsize{\normalfont\footnotesize\box\@tempboxa\hfil}%
\fi\fi\fi}
\fi



% V1.7 disable captions class option, do so in a way that retains operation of \label
% within \caption
\ifCLASSOPTIONcaptionsoff
\long\def\@makecaption#1#2{\vspace*{2em}\footnotesize\begin{center}{\footnotesize #1}\end{center}%
\let\@FMFPtemporiglabeldefsave\label
\let\@FMFPtemplabelargsave\relax
\def\label##1{\gdef\@FMFPtemplabelargsave{##1}}%
\setbox\@tempboxa\hbox{#2}%
\let\label\@FMFPtemporiglabeldefsave
\ifx\@FMFPtemplabelargsave\relax\else\label{\@FMFPtemplabelargsave}\fi}
\fi


% V1.7 define end environments with \def not \let so as to work OK with
% preview-latex
\newcounter{figure}
\def\thefigure{\@arabic\c@figure}
\def\fps@figure{tbp}
\def\ftype@figure{1}
\def\ext@figure{lof}
\def\fnum@figure{\figurename~\thefigure}
\def\figure{\@float{figure}}
\def\endfigure{\end@float}
\@namedef{figure*}{\@dblfloat{figure}}
\@namedef{endfigure*}{\end@dblfloat}
\newcounter{table}
\ifCLASSOPTIONcompsoc
\def\thetable{\arabic{table}}
\else
%\def\thetable{\@Roman\c@table} %Changed by AJC on 28.04.2018
\def\thetable{\arabic{table}} %Changed by AJC on 28.04.2018
\fi
\def\fps@table{tbp}
\def\ftype@table{2}
\def\ext@table{lot}
\def\fnum@table{\tablename~\thetable}
% V1.6 FMFP uses 8pt text for tables
% to default to footnotesize, we hack into LaTeX2e's \@floatboxreset and pray
\def\table{\def\@floatboxreset{\reset@font\scriptsize\@setminipage}%
  \let\@makefntext\@maketablefntext
  \@float{table}}
\def\endtable{\end@float}
% v1.6b double column tables need to default to footnotesize as well.
\@namedef{table*}{\def\@floatboxreset{\reset@font\scriptsize\@setminipage}\@dblfloat{table}}
\@namedef{endtable*}{\end@dblfloat}




%%
%% START OF FMFPeqnarry DEFINITIONS
%%
%% Inspired by the concepts, examples, and previous works of LaTeX 
%% coders and developers such as Donald Arseneau, Fred Bartlett, 
%% David Carlisle, Tony Liu, Frank Mittelbach, Piet van Oostrum, 
%% Roland Winkler and Mark Wooding.
%% I don't make the claim that my work here is even near their calibre. ;)


% hook to allow easy changeover to FMFPtran.cls/tools.sty error reporting
\def\@FMFPclspkgerror{\ClassError{FMFPtran}}

\newif\if@FMFPeqnarraystarform% flag to indicate if the environment was called as the star form
\@FMFPeqnarraystarformfalse

\newif\if@advanceFMFPeqncolcnt% tracks if the environment should advance the col counter
% allows a way to make an \FMFPeqnarraybox that can be used within an \FMFPeqnarray
% used by FMFPeqnarraymulticol so that it can work properly in both
\@advanceFMFPeqncolcnttrue

\newcount\@FMFPeqnnumcols % tracks how many FMFPeqnarray cols are defined
\newcount\@FMFPeqncolcnt  % tracks how many FMFPeqnarray cols the user actually used


% The default math style used by the columns
\def\FMFPeqnarraymathstyle{\displaystyle}
% The default text style used by the columns
% default to using the current font
\def\FMFPeqnarraytextstyle{\relax}

% like the iedlistdecl but for \FMFPeqnarray
\def\FMFPeqnarraydecl{\relax}
\def\FMFPeqnarrayboxdecl{\relax}

% \yesnumber is the opposite of \nonumber
% a novel concept with the same def as the equationarray package
% However, we give FMFP versions too since some LaTeX packages such as 
% the MDWtools mathenv.sty redefine \nonumber to something else.
\providecommand{\yesnumber}{\global\@eqnswtrue}
\def\FMFPyesnumber{\global\@eqnswtrue}
\def\FMFPnonumber{\global\@eqnswfalse}


\def\FMFPyessubnumber{\global\@FMFPissubequationtrue\global\@eqnswtrue%
\if@FMFPeqnarrayISinner% only do something inside an FMFPeqnarray
\if@FMFPlastlinewassubequation\addtocounter{equation}{-1}\else\setcounter{FMFPsubequation}{1}\fi%
\def\@currentlabel{\p@FMFPsubequation\theFMFPsubequation}\fi}

% flag to indicate that an equation is a sub equation
\newif\if@FMFPissubequation%
\@FMFPissubequationfalse

% allows users to "push away" equations that get too close to the equation numbers
\def\FMFPeqnarraynumspace{\hphantom{\if@FMFPissubequation\theFMFPsubequationdis\else\theequationdis\fi}}

% provides a way to span multiple columns within FMFPeqnarray environments
% will consider \if@advanceFMFPeqncolcnt before globally advancing the
% column counter - so as to work within \FMFPeqnarraybox
% usage: \FMFPeqnarraymulticol{number cols. to span}{col type}{cell text}
\long\def\FMFPeqnarraymulticol#1#2#3{\multispan{#1}%
% check if column is defined
\relax\expandafter\ifx\csname @FMFPeqnarraycolDEF#2\endcsname\@FMFPeqnarraycolisdefined%
\csname @FMFPeqnarraycolPRE#2\endcsname#3\relax\relax\relax\relax\relax%
\relax\relax\relax\relax\relax\csname @FMFPeqnarraycolPOST#2\endcsname%
\else% if not, error and use default type
\@FMFPclspkgerror{Invalid column type "#2" in \string\FMFPeqnarraymulticol.\MessageBreak
Using a default centering column instead}%
{You must define FMFPeqnarray column types before use.}%
\csname @FMFPeqnarraycolPRE@FMFPdefault\endcsname#3\relax\relax\relax\relax\relax%
\relax\relax\relax\relax\relax\csname @FMFPeqnarraycolPOST@FMFPdefault\endcsname%
\fi%
% advance column counter only if the FMFPeqnarray environment wants it
\if@advanceFMFPeqncolcnt\global\advance\@FMFPeqncolcnt by #1\relax\fi}

% like \omit, but maintains track of the column counter for \FMFPeqnarray
\def\FMFPeqnarrayomit{\omit\if@advanceFMFPeqncolcnt\global\advance\@FMFPeqncolcnt by 1\relax\fi}


% provides a way to define a letter referenced column type
% usage: \FMFPeqnarraydefcol{col. type letter/name}{pre insertion text}{post insertion text}
\def\FMFPeqnarraydefcol#1#2#3{\expandafter\def\csname @FMFPeqnarraycolPRE#1\endcsname{#2}%
\expandafter\def\csname @FMFPeqnarraycolPOST#1\endcsname{#3}%
\expandafter\def\csname @FMFPeqnarraycolDEF#1\endcsname{1}}


% provides a way to define a numerically referenced inter-column glue types
% usage: \FMFPeqnarraydefcolsep{col. glue number}{glue definition}
\def\FMFPeqnarraydefcolsep#1#2{\expandafter\def\csname @FMFPeqnarraycolSEP\romannumeral #1\endcsname{#2}%
\expandafter\def\csname @FMFPeqnarraycolSEPDEF\romannumeral #1\endcsname{1}}


\def\@FMFPeqnarraycolisdefined{1}% just a macro for 1, used for checking undefined column types


% expands and appends the given argument to the \@FMFPtrantmptoksA token list
% used to build up the \halign preamble
\def\@FMFPappendtoksA#1{\edef\@@FMFPappendtoksA{\@FMFPtrantmptoksA={\the\@FMFPtrantmptoksA #1}}%
\@@FMFPappendtoksA}

% also appends to \@FMFPtrantmptoksA, but does not expand the argument
% uses \toks8 as a scratchpad register
\def\@FMFPappendNOEXPANDtoksA#1{\toks8={#1}%
\edef\@@FMFPappendNOEXPANDtoksA{\@FMFPtrantmptoksA={\the\@FMFPtrantmptoksA\the\toks8}}%
\@@FMFPappendNOEXPANDtoksA}

% define some common column types for the user
% math
\FMFPeqnarraydefcol{l}{$\FMFPeqnarraymathstyle}{$\hfil}
\FMFPeqnarraydefcol{c}{\hfil$\FMFPeqnarraymathstyle}{$\hfil}
\FMFPeqnarraydefcol{r}{\hfil$\FMFPeqnarraymathstyle}{$}
\FMFPeqnarraydefcol{L}{$\FMFPeqnarraymathstyle{}}{{}$\hfil}
\FMFPeqnarraydefcol{C}{\hfil$\FMFPeqnarraymathstyle{}}{{}$\hfil}
\FMFPeqnarraydefcol{R}{\hfil$\FMFPeqnarraymathstyle{}}{{}$}
% text
\FMFPeqnarraydefcol{s}{\FMFPeqnarraytextstyle}{\hfil}
\FMFPeqnarraydefcol{t}{\hfil\FMFPeqnarraytextstyle}{\hfil}
\FMFPeqnarraydefcol{u}{\hfil\FMFPeqnarraytextstyle}{}

% vertical rules
\FMFPeqnarraydefcol{v}{}{\vrule width\arrayrulewidth}
\FMFPeqnarraydefcol{vv}{\vrule width\arrayrulewidth\hfil}{\hfil\vrule width\arrayrulewidth}
\FMFPeqnarraydefcol{V}{}{\vrule width\arrayrulewidth\hskip\doublerulesep\vrule width\arrayrulewidth}
\FMFPeqnarraydefcol{VV}{\vrule width\arrayrulewidth\hskip\doublerulesep\vrule width\arrayrulewidth\hfil}%
{\hfil\vrule width\arrayrulewidth\hskip\doublerulesep\vrule width\arrayrulewidth}

% horizontal rules
\FMFPeqnarraydefcol{h}{}{\leaders\hrule height\arrayrulewidth\hfil}
\FMFPeqnarraydefcol{H}{}{\leaders\vbox{\hrule width\arrayrulewidth\vskip\doublerulesep\hrule width\arrayrulewidth}\hfil}

% plain
\FMFPeqnarraydefcol{x}{}{}
\FMFPeqnarraydefcol{X}{$}{$}

% the default column type to use in the event a column type is not defined
\FMFPeqnarraydefcol{@FMFPdefault}{\hfil$\FMFPeqnarraymathstyle}{$\hfil}


% a zero tabskip (used for "-" col types)
\def\@FMFPeqnarraycolSEPzero{0pt plus 0pt minus 0pt}
% a centering tabskip (used for "+" col types)
\def\@FMFPeqnarraycolSEPcenter{1000pt plus 0pt minus 1000pt}

% top level default tabskip glues for the start, end, and inter-column
% may be reset within environments not always at the top level, e.g., \FMFPeqnarraybox
\edef\@FMFPeqnarraycolSEPdefaultstart{\@FMFPeqnarraycolSEPcenter}% default start glue
\edef\@FMFPeqnarraycolSEPdefaultend{\@FMFPeqnarraycolSEPcenter}% default end glue
\edef\@FMFPeqnarraycolSEPdefaultmid{\@FMFPeqnarraycolSEPzero}% default inter-column glue



% creates a vertical rule that extends from the bottom to the top a a cell
% Provided in case other packages redefine \vline some other way.
% usage: \FMFPeqnarrayvrule[rule thickness]
% If no argument is provided, \arrayrulewidth will be used for the rule thickness. 
\newcommand\FMFPeqnarrayvrule[1][\arrayrulewidth]{\vrule\@width#1\relax}

% creates a blank separator row
% usage: \FMFPeqnarrayseprow[separation length][font size commands]
% default is \FMFPeqnarrayseprow[0.25\normalbaselineskip][\relax]
% blank arguments inherit the default values
% uses \skip5 as a scratch register - calls \@FMFPeqnarraystrutsize which uses more scratch registers
\def\FMFPeqnarrayseprow{\relax\@ifnextchar[{\@FMFPeqnarrayseprow}{\@FMFPeqnarrayseprow[0.25\normalbaselineskip]}}
\def\@FMFPeqnarrayseprow[#1]{\relax\@ifnextchar[{\@@FMFPeqnarrayseprow[#1]}{\@@FMFPeqnarrayseprow[#1][\relax]}}
\def\@@FMFPeqnarrayseprow[#1][#2]{\def\@FMFPeqnarrayseprowARGONE{#1}%
\ifx\@FMFPeqnarrayseprowARGONE\@empty%
% get the skip value, based on the font commands
% use skip5 because \FMFPeqnarraystrutsize uses \skip0, \skip2, \skip3
% assign within a bogus box to confine the font changes
{\setbox0=\hbox{#2\relax\global\skip5=0.25\normalbaselineskip}}%
\else%
{\setbox0=\hbox{#2\relax\global\skip5=#1}}%
\fi%
\@FMFPeqnarrayhoptolastcolumn\FMFPeqnarraystrutsize{\skip5}{0pt}[\relax]\relax}

% creates a blank separator row, but omits all the column templates
% usage: \FMFPeqnarrayseprowcut[separation length][font size commands]
% default is \FMFPeqnarrayseprowcut[0.25\normalbaselineskip][\relax]
% blank arguments inherit the default values
% uses \skip5 as a scratch register - calls \@FMFPeqnarraystrutsize which uses more scratch registers
\def\FMFPeqnarrayseprowcut{\multispan{\@FMFPeqnnumcols}\relax% span all the cols
% advance column counter only if the FMFPeqnarray environment wants it
\if@advanceFMFPeqncolcnt\global\advance\@FMFPeqncolcnt by \@FMFPeqnnumcols\relax\fi%
\@ifnextchar[{\@FMFPeqnarrayseprowcut}{\@FMFPeqnarrayseprowcut[0.25\normalbaselineskip]}}
\def\@FMFPeqnarrayseprowcut[#1]{\relax\@ifnextchar[{\@@FMFPeqnarrayseprowcut[#1]}{\@@FMFPeqnarrayseprowcut[#1][\relax]}}
\def\@@FMFPeqnarrayseprowcut[#1][#2]{\def\@FMFPeqnarrayseprowARGONE{#1}%
\ifx\@FMFPeqnarrayseprowARGONE\@empty%
% get the skip value, based on the font commands
% use skip5 because \FMFPeqnarraystrutsize uses \skip0, \skip2, \skip3
% assign within a bogus box to confine the font changes
{\setbox0=\hbox{#2\relax\global\skip5=0.25\normalbaselineskip}}%
\else%
{\setbox0=\hbox{#2\relax\global\skip5=#1}}%
\fi%
\FMFPeqnarraystrutsize{\skip5}{0pt}[\relax]\relax}



% draws a single rule across all the columns optional
% argument determines the rule width, \arrayrulewidth is the default
% updates column counter as needed and turns off struts
% usage: \FMFPeqnarrayrulerow[rule line thickness]
\def\FMFPeqnarrayrulerow{\multispan{\@FMFPeqnnumcols}\relax% span all the cols
% advance column counter only if the FMFPeqnarray environment wants it
\if@advanceFMFPeqncolcnt\global\advance\@FMFPeqncolcnt by \@FMFPeqnnumcols\relax\fi%
\@ifnextchar[{\@FMFPeqnarrayrulerow}{\@FMFPeqnarrayrulerow[\arrayrulewidth]}}
\def\@FMFPeqnarrayrulerow[#1]{\leaders\hrule height#1\hfil\relax% put in our rule 
% turn off any struts
\FMFPeqnarraystrutsize{0pt}{0pt}[\relax]\relax}


% draws a double rule by using a single rule row, a separator row, and then
% another single rule row 
% first optional argument determines the rule thicknesses, \arrayrulewidth is the default
% second optional argument determines the rule spacing, \doublerulesep is the default
% usage: \FMFPeqnarraydblrulerow[rule line thickness][rule spacing]
\def\FMFPeqnarraydblrulerow{\multispan{\@FMFPeqnnumcols}\relax% span all the cols
% advance column counter only if the FMFPeqnarray environment wants it
\if@advanceFMFPeqncolcnt\global\advance\@FMFPeqncolcnt by \@FMFPeqnnumcols\relax\fi%
\@ifnextchar[{\@FMFPeqnarraydblrulerow}{\@FMFPeqnarraydblrulerow[\arrayrulewidth]}}
\def\@FMFPeqnarraydblrulerow[#1]{\relax\@ifnextchar[{\@@FMFPeqnarraydblrulerow[#1]}%
{\@@FMFPeqnarraydblrulerow[#1][\doublerulesep]}}
\def\@@FMFPeqnarraydblrulerow[#1][#2]{\def\@FMFPeqnarraydblrulerowARG{#1}%
% we allow the user to say \FMFPeqnarraydblrulerow[][]
\ifx\@FMFPeqnarraydblrulerowARG\@empty%
\@FMFPeqnarrayrulerow[\arrayrulewidth]%
\else%
\@FMFPeqnarrayrulerow[#1]\relax%
\fi%
\def\@FMFPeqnarraydblrulerowARG{#2}%
\ifx\@FMFPeqnarraydblrulerowARG\@empty%
\\\FMFPeqnarrayseprow[\doublerulesep][\relax]%
\else%
\\\FMFPeqnarrayseprow[#2][\relax]%
\fi%
\\\multispan{\@FMFPeqnnumcols}%
% advance column counter only if the FMFPeqnarray environment wants it
\if@advanceFMFPeqncolcnt\global\advance\@FMFPeqncolcnt by \@FMFPeqnnumcols\relax\fi%
\def\@FMFPeqnarraydblrulerowARG{#1}%
\ifx\@FMFPeqnarraydblrulerowARG\@empty%
\@FMFPeqnarrayrulerow[\arrayrulewidth]%
\else%
\@FMFPeqnarrayrulerow[#1]%
\fi%
}

% draws a double rule by using a single rule row, a separator (cutting) row, and then
% another single rule row 
% first optional argument determines the rule thicknesses, \arrayrulewidth is the default
% second optional argument determines the rule spacing, \doublerulesep is the default
% usage: \FMFPeqnarraydblrulerow[rule line thickness][rule spacing]
\def\FMFPeqnarraydblrulerowcut{\multispan{\@FMFPeqnnumcols}\relax% span all the cols
% advance column counter only if the FMFPeqnarray environment wants it
\if@advanceFMFPeqncolcnt\global\advance\@FMFPeqncolcnt by \@FMFPeqnnumcols\relax\fi%
\@ifnextchar[{\@FMFPeqnarraydblrulerowcut}{\@FMFPeqnarraydblrulerowcut[\arrayrulewidth]}}
\def\@FMFPeqnarraydblrulerowcut[#1]{\relax\@ifnextchar[{\@@FMFPeqnarraydblrulerowcut[#1]}%
{\@@FMFPeqnarraydblrulerowcut[#1][\doublerulesep]}}
\def\@@FMFPeqnarraydblrulerowcut[#1][#2]{\def\@FMFPeqnarraydblrulerowARG{#1}%
% we allow the user to say \FMFPeqnarraydblrulerow[][]
\ifx\@FMFPeqnarraydblrulerowARG\@empty%
\@FMFPeqnarrayrulerow[\arrayrulewidth]%
\else%
\@FMFPeqnarrayrulerow[#1]%
\fi%
\def\@FMFPeqnarraydblrulerowARG{#2}%
\ifx\@FMFPeqnarraydblrulerowARG\@empty%
\\\FMFPeqnarrayseprowcut[\doublerulesep][\relax]%
\else%
\\\FMFPeqnarrayseprowcut[#2][\relax]%
\fi%
\\\multispan{\@FMFPeqnnumcols}%
% advance column counter only if the FMFPeqnarray environment wants it
\if@advanceFMFPeqncolcnt\global\advance\@FMFPeqncolcnt by \@FMFPeqnnumcols\relax\fi%
\def\@FMFPeqnarraydblrulerowARG{#1}%
\ifx\@FMFPeqnarraydblrulerowARG\@empty%
\@FMFPeqnarrayrulerow[\arrayrulewidth]%
\else%
\@FMFPeqnarrayrulerow[#1]%
\fi%
}



% inserts a full row's worth of &'s
% relies on \@FMFPeqnnumcols to provide the correct number of columns
% uses \@FMFPtrantmptoksA, \count0 as scratch registers
\def\@FMFPeqnarrayhoptolastcolumn{\@FMFPtrantmptoksA={}\count0=1\relax%
\loop% add cols if the user did not use them all
\ifnum\count0<\@FMFPeqnnumcols\relax%
\@FMFPappendtoksA{&}%
\advance\count0 by 1\relax% update the col count
\repeat%
\the\@FMFPtrantmptoksA%execute the &'s
}



\newif\if@FMFPeqnarrayISinner % flag to indicate if we are within the lines
\@FMFPeqnarrayISinnerfalse    % of an FMFPeqnarray - after the FMFPeqnarraydecl

\edef\@FMFPeqnarrayTHEstrutheight{0pt} % height and depth of FMFPeqnarray struts
\edef\@FMFPeqnarrayTHEstrutdepth{0pt}

\edef\@FMFPeqnarrayTHEmasterstrutheight{0pt} % default height and depth of
\edef\@FMFPeqnarrayTHEmasterstrutdepth{0pt}  % struts within an FMFPeqnarray

\edef\@FMFPeqnarrayTHEmasterstrutHSAVE{0pt} % saved master strut height
\edef\@FMFPeqnarrayTHEmasterstrutDSAVE{0pt} % and depth

\newif\if@FMFPeqnarrayusemasterstrut % flag to indicate that the master strut value
\@FMFPeqnarrayusemasterstruttrue     % is to be used



% saves the strut height and depth of the master strut
\def\@FMFPeqnarraymasterstrutsave{\relax%
\expandafter\skip0=\@FMFPeqnarrayTHEmasterstrutheight\relax%
\expandafter\skip2=\@FMFPeqnarrayTHEmasterstrutdepth\relax%
% remove stretchability
\dimen0\skip0\relax%
\dimen2\skip2\relax%
% save values
\edef\@FMFPeqnarrayTHEmasterstrutHSAVE{\the\dimen0}%
\edef\@FMFPeqnarrayTHEmasterstrutDSAVE{\the\dimen2}}

% restores the strut height and depth of the master strut
\def\@FMFPeqnarraymasterstrutrestore{\relax%
\expandafter\skip0=\@FMFPeqnarrayTHEmasterstrutHSAVE\relax%
\expandafter\skip2=\@FMFPeqnarrayTHEmasterstrutDSAVE\relax%
% remove stretchability
\dimen0\skip0\relax%
\dimen2\skip2\relax%
% restore values
\edef\@FMFPeqnarrayTHEmasterstrutheight{\the\dimen0}%
\edef\@FMFPeqnarrayTHEmasterstrutdepth{\the\dimen2}}


% globally restores the strut height and depth to the 
% master values and sets the master strut flag to true
\def\@FMFPeqnarraystrutreset{\relax%
\expandafter\skip0=\@FMFPeqnarrayTHEmasterstrutheight\relax%
\expandafter\skip2=\@FMFPeqnarrayTHEmasterstrutdepth\relax%
% remove stretchability
\dimen0\skip0\relax%
\dimen2\skip2\relax%
% restore values
\xdef\@FMFPeqnarrayTHEstrutheight{\the\dimen0}%
\xdef\@FMFPeqnarrayTHEstrutdepth{\the\dimen2}%
\global\@FMFPeqnarrayusemasterstruttrue}


% if the master strut is not to be used, make the current
% values of \@FMFPeqnarrayTHEstrutheight, \@FMFPeqnarrayTHEstrutdepth
% and the use master strut flag, global
% this allows user strut commands issued in the last column to be carried
% into the isolation/strut column
\def\@FMFPeqnarrayglobalizestrutstatus{\relax%
\if@FMFPeqnarrayusemasterstrut\else%
\xdef\@FMFPeqnarrayTHEstrutheight{\@FMFPeqnarrayTHEstrutheight}%
\xdef\@FMFPeqnarrayTHEstrutdepth{\@FMFPeqnarrayTHEstrutdepth}%
\global\@FMFPeqnarrayusemasterstrutfalse%
\fi}



% usage: \FMFPeqnarraystrutsize{height}{depth}[font size commands]
% If called outside the lines of an FMFPeqnarray, sets the height
% and depth of both the master and local struts. If called inside
% an FMFPeqnarray line, sets the height and depth of the local strut
% only and sets the flag to indicate the use of the local strut
% values. If the height or depth is left blank, 0.7\normalbaselineskip
% and 0.3\normalbaselineskip will be used, respectively.
% The optional argument can be used to evaluate the lengths under
% a different font size and styles. If none is specified, the current
% font is used.
% uses scratch registers \skip0, \skip2, \skip3, \dimen0, \dimen2
\def\FMFPeqnarraystrutsize#1#2{\relax\@ifnextchar[{\@FMFPeqnarraystrutsize{#1}{#2}}{\@FMFPeqnarraystrutsize{#1}{#2}[\relax]}}
\def\@FMFPeqnarraystrutsize#1#2[#3]{\def\@FMFPeqnarraystrutsizeARG{#1}%
\ifx\@FMFPeqnarraystrutsizeARG\@empty%
{\setbox0=\hbox{#3\relax\global\skip3=0.7\normalbaselineskip}}%
\skip0=\skip3\relax%
\else% arg one present
{\setbox0=\hbox{#3\relax\global\skip3=#1\relax}}%
\skip0=\skip3\relax%
\fi% if null arg
\def\@FMFPeqnarraystrutsizeARG{#2}%
\ifx\@FMFPeqnarraystrutsizeARG\@empty%
{\setbox0=\hbox{#3\relax\global\skip3=0.3\normalbaselineskip}}%
\skip2=\skip3\relax%
\else% arg two present
{\setbox0=\hbox{#3\relax\global\skip3=#2\relax}}%
\skip2=\skip3\relax%
\fi% if null arg
% remove stretchability, just to be safe
\dimen0\skip0\relax%
\dimen2\skip2\relax%
% dimen0 = height, dimen2 = depth
\if@FMFPeqnarrayISinner% inner does not touch master strut size
\edef\@FMFPeqnarrayTHEstrutheight{\the\dimen0}%
\edef\@FMFPeqnarrayTHEstrutdepth{\the\dimen2}%
\@FMFPeqnarrayusemasterstrutfalse% do not use master
\else% outer, have to set master strut too
\edef\@FMFPeqnarrayTHEmasterstrutheight{\the\dimen0}%
\edef\@FMFPeqnarrayTHEmasterstrutdepth{\the\dimen2}%
\edef\@FMFPeqnarrayTHEstrutheight{\the\dimen0}%
\edef\@FMFPeqnarrayTHEstrutdepth{\the\dimen2}%
\@FMFPeqnarrayusemasterstruttrue% use master strut
\fi}


% usage: \FMFPeqnarraystrutsizeadd{added height}{added depth}[font size commands]
% If called outside the lines of an FMFPeqnarray, adds the given height
% and depth to both the master and local struts.
% If called inside an FMFPeqnarray line, adds the given height and depth
% to the local strut only and sets the flag to indicate the use 
% of the local strut values.
% In both cases, if a height or depth is left blank, 0pt is used instead.
% The optional argument can be used to evaluate the lengths under
% a different font size and styles. If none is specified, the current
% font is used.
% uses scratch registers \skip0, \skip2, \skip3, \dimen0, \dimen2
\def\FMFPeqnarraystrutsizeadd#1#2{\relax\@ifnextchar[{\@FMFPeqnarraystrutsizeadd{#1}{#2}}{\@FMFPeqnarraystrutsizeadd{#1}{#2}[\relax]}}
\def\@FMFPeqnarraystrutsizeadd#1#2[#3]{\def\@FMFPeqnarraystrutsizearg{#1}%
\ifx\@FMFPeqnarraystrutsizearg\@empty%
\skip0=0pt\relax%
\else% arg one present
{\setbox0=\hbox{#3\relax\global\skip3=#1}}%
\skip0=\skip3\relax%
\fi% if null arg
\def\@FMFPeqnarraystrutsizearg{#2}%
\ifx\@FMFPeqnarraystrutsizearg\@empty%
\skip2=0pt\relax%
\else% arg two present
{\setbox0=\hbox{#3\relax\global\skip3=#2}}%
\skip2=\skip3\relax%
\fi% if null arg
% remove stretchability, just to be safe
\dimen0\skip0\relax%
\dimen2\skip2\relax%
% dimen0 = height, dimen2 = depth
\if@FMFPeqnarrayISinner% inner does not touch master strut size
% get local strut size
\expandafter\skip0=\@FMFPeqnarrayTHEstrutheight\relax%
\expandafter\skip2=\@FMFPeqnarrayTHEstrutdepth\relax%
% add it to the user supplied values
\advance\dimen0 by \skip0\relax%
\advance\dimen2 by \skip2\relax%
% update the local strut size
\edef\@FMFPeqnarrayTHEstrutheight{\the\dimen0}%
\edef\@FMFPeqnarrayTHEstrutdepth{\the\dimen2}%
\@FMFPeqnarrayusemasterstrutfalse% do not use master
\else% outer, have to set master strut too
% get master strut size
\expandafter\skip0=\@FMFPeqnarrayTHEmasterstrutheight\relax%
\expandafter\skip2=\@FMFPeqnarrayTHEmasterstrutdepth\relax%
% add it to the user supplied values
\advance\dimen0 by \skip0\relax%
\advance\dimen2 by \skip2\relax%
% update the local and master strut sizes
\edef\@FMFPeqnarrayTHEmasterstrutheight{\the\dimen0}%
\edef\@FMFPeqnarrayTHEmasterstrutdepth{\the\dimen2}%
\edef\@FMFPeqnarrayTHEstrutheight{\the\dimen0}%
\edef\@FMFPeqnarrayTHEstrutdepth{\the\dimen2}%
\@FMFPeqnarrayusemasterstruttrue% use master strut
\fi}


% allow user a way to see the struts
\newif\ifFMFPvisiblestruts
\FMFPvisiblestrutsfalse

% inserts an invisible strut using the master or local strut values
% uses scratch registers \skip0, \skip2, \dimen0, \dimen2
\def\@FMFPeqnarrayinsertstrut{\relax%
\if@FMFPeqnarrayusemasterstrut
% get master strut size
\expandafter\skip0=\@FMFPeqnarrayTHEmasterstrutheight\relax%
\expandafter\skip2=\@FMFPeqnarrayTHEmasterstrutdepth\relax%
\else%
% get local strut size
\expandafter\skip0=\@FMFPeqnarrayTHEstrutheight\relax%
\expandafter\skip2=\@FMFPeqnarrayTHEstrutdepth\relax%
\fi%
% remove stretchability, probably not needed
\dimen0\skip0\relax%
\dimen2\skip2\relax%
% dimen0 = height, dimen2 = depth
% allow user to see struts if desired
\ifFMFPvisiblestruts%
\vrule width0.2pt height\dimen0 depth\dimen2\relax%
\else%
\vrule width0pt height\dimen0 depth\dimen2\relax\fi}


% creates an invisible strut, useable even outside \FMFPeqnarray
% if \FMFPvisiblestrutstrue, the strut will be visible and 0.2pt wide. 
% usage: \FMFPstrut[height][depth][font size commands]
% default is \FMFPstrut[0.7\normalbaselineskip][0.3\normalbaselineskip][\relax]
% blank arguments inherit the default values
% uses \dimen0, \dimen2, \skip0, \skip2
\def\FMFPstrut{\relax\@ifnextchar[{\@FMFPstrut}{\@FMFPstrut[0.7\normalbaselineskip]}}
\def\@FMFPstrut[#1]{\relax\@ifnextchar[{\@@FMFPstrut[#1]}{\@@FMFPstrut[#1][0.3\normalbaselineskip]}}
\def\@@FMFPstrut[#1][#2]{\relax\@ifnextchar[{\@@@FMFPstrut[#1][#2]}{\@@@FMFPstrut[#1][#2][\relax]}}
\def\@@@FMFPstrut[#1][#2][#3]{\mbox{#3\relax%
\def\@FMFPstrutARG{#1}%
\ifx\@FMFPstrutARG\@empty%
\skip0=0.7\normalbaselineskip\relax%
\else%
\skip0=#1\relax%
\fi%
\def\@FMFPstrutARG{#2}%
\ifx\@FMFPstrutARG\@empty%
\skip2=0.3\normalbaselineskip\relax%
\else%
\skip2=#2\relax%
\fi%
% remove stretchability, probably not needed
\dimen0\skip0\relax%
\dimen2\skip2\relax%
\ifFMFPvisiblestruts%
\vrule width0.2pt height\dimen0 depth\dimen2\relax%
\else%
\vrule width0.0pt height\dimen0 depth\dimen2\relax\fi}}


% enables strut mode by setting a default strut size and then zeroing the
% \baselineskip, \lineskip, \lineskiplimit and \jot
\def\FMFPeqnarraystrutmode{\FMFPeqnarraystrutsize{0.7\normalbaselineskip}{0.3\normalbaselineskip}[\relax]%
\baselineskip=0pt\lineskip=0pt\lineskiplimit=0pt\jot=0pt}



\def\FMFPeqnarray{\@FMFPeqnarraystarformfalse\@FMFPeqnarray}
\def\endFMFPeqnarray{\end@FMFPeqnarray}

\@namedef{FMFPeqnarray*}{\@FMFPeqnarraystarformtrue\@FMFPeqnarray}
\@namedef{endFMFPeqnarray*}{\end@FMFPeqnarray}


% \FMFPeqnarray is an enhanced \eqnarray. 
% The star form defaults to not putting equation numbers at the end of each row.
% usage: \FMFPeqnarray[decl]{cols}
\def\@FMFPeqnarray{\relax\@ifnextchar[{\@@FMFPeqnarray}{\@@FMFPeqnarray[\relax]}}
\def\@@FMFPeqnarray[#1]#2{%
   % default to showing the equation number or not based on whether or not
   % the star form was involked
   \if@FMFPeqnarraystarform\global\@eqnswfalse
   \else% not the star form
   \global\@eqnswtrue
   \fi% if star form
   \@FMFPissubequationfalse% default to no subequations
   \@FMFPlastlinewassubequationfalse% assume last line is not a sub equation
   \@FMFPeqnarrayISinnerfalse% not yet within the lines of the halign
   \@FMFPeqnarraystrutsize{0pt}{0pt}[\relax]% turn off struts by default
   \@FMFPeqnarrayusemasterstruttrue% use master strut till user asks otherwise
   \FMFPvisiblestrutsfalse% diagnostic mode defaults to off
   % no extra space unless the user specifically requests it
   \lineskip=0pt\relax
   \lineskiplimit=0pt\relax
   \baselineskip=\normalbaselineskip\relax%
   \jot=\FMFPnormaljot\relax%
   \mathsurround\z@\relax% no extra spacing around math
   \@advanceFMFPeqncolcnttrue% advance the col counter for each col the user uses, 
                             % used in \FMFPeqnarraymulticol and in the preamble build
   \stepcounter{equation}% advance equation counter before first line
   \setcounter{FMFPsubequation}{0}% no subequation yet 
   \def\@currentlabel{\p@equation\theequation}% redefine the ref label
   \FMFPeqnarraydecl\relax% allow a way for the user to make global overrides
   #1\relax% allow user to override defaults
   \let\\\@FMFPeqnarraycr% replace newline with one that can put in eqn. numbers
   \global\@FMFPeqncolcnt\z@% col. count = 0 for first line
   \@FMFPbuildpreamble #2\end\relax% build the preamble and put it into \@FMFPtrantmptoksA 
   % put in the column for the equation number
   \ifnum\@FMFPeqnnumcols>0\relax\@FMFPappendtoksA{&}\fi% col separator for those after the first
   \toks0={##}%
   % advance the \@FMFPeqncolcnt for the isolation col, this helps with error checking
   \@FMFPappendtoksA{\global\advance\@FMFPeqncolcnt by 1\relax}%
   % add the isolation column
   \@FMFPappendtoksA{\tabskip\z@skip\bgroup\the\toks0\egroup}%
   % advance the \@FMFPeqncolcnt for the equation number col, this helps with error checking
   \@FMFPappendtoksA{&\global\advance\@FMFPeqncolcnt by 1\relax}%
   % add the equation number col to the preamble
   \@FMFPappendtoksA{\tabskip\z@skip\hb@xt@\z@\bgroup\hss\the\toks0\egroup}%
   % note \@FMFPeqnnumcols does not count the equation col or isolation col
   % set the starting tabskip glue as determined by the preamble build
   \tabskip=\@FMFPBPstartglue\relax
   % begin the display alignment
   \@FMFPeqnarrayISinnertrue% commands are now within the lines
   $$\everycr{}\halign to\displaywidth\bgroup
   % "exspand" the preamble
   \span\the\@FMFPtrantmptoksA\cr}

% enter isolation/strut column (or the next column if the user did not use
% every column), record the strut status, complete the columns, do the strut if needed,
% restore counters to correct values and exit
\def\end@FMFPeqnarray{\@FMFPeqnarrayglobalizestrutstatus&\@@FMFPeqnarraycr\egroup%
\if@FMFPlastlinewassubequation\global\advance\c@FMFPsubequation\m@ne\fi%
\global\advance\c@equation\m@ne%
$$\@ignoretrue}

% need a way to remember if last line is a subequation
\newif\if@FMFPlastlinewassubequation%
\@FMFPlastlinewassubequationfalse

% FMFPeqnarray uses a modifed \\ instead of the plain \cr to
% end rows. This allows for things like \\*[vskip amount]
% This "cr" macros are modified versions those for LaTeX2e's eqnarray
% the {\ifnum0=`} braces must be kept away from the last column to avoid
% altering spacing of its math, so we use & to advance to the next column
% as there is an isolation/strut column after the user's columns
\def\@FMFPeqnarraycr{\@FMFPeqnarrayglobalizestrutstatus&% save strut status and advance to next column
   {\ifnum0=`}\fi
   \@ifstar{%
      \global\@eqpen\@M\@FMFPeqnarrayYCR
   }{%
      \global\@eqpen\interdisplaylinepenalty \@FMFPeqnarrayYCR
   }%
}

\def\@FMFPeqnarrayYCR{\@testopt\@FMFPeqnarrayXCR\z@skip}

\def\@FMFPeqnarrayXCR[#1]{%
   \ifnum0=`{\fi}%
   \@@FMFPeqnarraycr
   \noalign{\penalty\@eqpen\vskip\jot\vskip #1\relax}}%

\def\@@FMFPeqnarraycr{\@FMFPtrantmptoksA={}% clear token register
    \advance\@FMFPeqncolcnt by -1\relax% adjust col count because of the isolation column
    \ifnum\@FMFPeqncolcnt>\@FMFPeqnnumcols\relax
    \@FMFPclspkgerror{Too many columns within the FMFPeqnarray\MessageBreak
                          environment}%
    {Use fewer \string &'s or put more columns in the FMFPeqnarry column\MessageBreak 
     specifications.}\relax%
    \else
    \loop% add cols if the user did not use them all
    \ifnum\@FMFPeqncolcnt<\@FMFPeqnnumcols\relax
    \@FMFPappendtoksA{&}%
    \advance\@FMFPeqncolcnt by 1\relax% update the col count
    \repeat
    % this number of &'s will take us the the isolation column
    \fi
    % execute the &'s
    \the\@FMFPtrantmptoksA%
    % handle the strut/isolation column
    \@FMFPeqnarrayinsertstrut% do the strut if needed
    \@FMFPeqnarraystrutreset% reset the strut system for next line or FMFPeqnarray
    &% and enter the equation number column
    % is this line needs an equation number, display it and advance the
    % (sub)equation counters, record what type this line was
    \if@eqnsw%
     \if@FMFPissubequation\theFMFPsubequationdis\addtocounter{equation}{1}\stepcounter{FMFPsubequation}%
     \global\@FMFPlastlinewassubequationtrue%
     \else% display a standard equation number, initialize the FMFPsubequation counter
     \theequationdis\stepcounter{equation}\setcounter{FMFPsubequation}{0}%
     \global\@FMFPlastlinewassubequationfalse\fi%
    \fi%
    % reset the eqnsw flag to indicate default preference of the display of equation numbers
    \if@FMFPeqnarraystarform\global\@eqnswfalse\else\global\@eqnswtrue\fi
    \global\@FMFPissubequationfalse% reset the subequation flag
    % reset the number of columns the user actually used
    \global\@FMFPeqncolcnt\z@\relax
    % the real end of the line
    \cr}





% \FMFPeqnarraybox is like \FMFPeqnarray except the box form puts everything
% inside a vtop, vbox, or vcenter box depending on the letter in the second
% optional argument (t,b,c). Vbox is the default. Unlike \FMFPeqnarray,
% equation numbers are not displayed and \FMFPeqnarraybox can be nested.
% \FMFPeqnarrayboxm is for math mode (like \array) and does not put the vbox
% within an hbox.
% \FMFPeqnarrayboxt is for text mode (like \tabular) and puts the vbox within
% a \hbox{$ $} construct.
% \FMFPeqnarraybox will auto detect whether to use \FMFPeqnarrayboxm or 
% \FMFPeqnarrayboxt depending on the math mode.
% The third optional argument specifies the width this box is to be set to -
% natural width is the default.
% The * forms do not add \jot line spacing
% usage: \FMFPeqnarraybox[decl][pos][width]{cols}
\def\FMFPeqnarrayboxm{\@FMFPeqnarraystarformfalse\@FMFPeqnarrayboxHBOXSWfalse\@FMFPeqnarraybox}
\def\endFMFPeqnarrayboxm{\end@FMFPeqnarraybox}
\@namedef{FMFPeqnarrayboxm*}{\@FMFPeqnarraystarformtrue\@FMFPeqnarrayboxHBOXSWfalse\@FMFPeqnarraybox}
\@namedef{endFMFPeqnarrayboxm*}{\end@FMFPeqnarraybox}

\def\FMFPeqnarrayboxt{\@FMFPeqnarraystarformfalse\@FMFPeqnarrayboxHBOXSWtrue\@FMFPeqnarraybox}
\def\endFMFPeqnarrayboxt{\end@FMFPeqnarraybox}
\@namedef{FMFPeqnarrayboxt*}{\@FMFPeqnarraystarformtrue\@FMFPeqnarrayboxHBOXSWtrue\@FMFPeqnarraybox}
\@namedef{endFMFPeqnarrayboxt*}{\end@FMFPeqnarraybox}

\def\FMFPeqnarraybox{\@FMFPeqnarraystarformfalse\ifmmode\@FMFPeqnarrayboxHBOXSWfalse\else\@FMFPeqnarrayboxHBOXSWtrue\fi%
\@FMFPeqnarraybox}
\def\endFMFPeqnarraybox{\end@FMFPeqnarraybox}

\@namedef{FMFPeqnarraybox*}{\@FMFPeqnarraystarformtrue\ifmmode\@FMFPeqnarrayboxHBOXSWfalse\else\@FMFPeqnarrayboxHBOXSWtrue\fi%
\@FMFPeqnarraybox}
\@namedef{endFMFPeqnarraybox*}{\end@FMFPeqnarraybox}

% flag to indicate if the \FMFPeqnarraybox needs to put things into an hbox{$ $} 
% for \vcenter in non-math mode
\newif\if@FMFPeqnarrayboxHBOXSW%
\@FMFPeqnarrayboxHBOXSWfalse

\def\@FMFPeqnarraybox{\relax\@ifnextchar[{\@@FMFPeqnarraybox}{\@@FMFPeqnarraybox[\relax]}}
\def\@@FMFPeqnarraybox[#1]{\relax\@ifnextchar[{\@@@FMFPeqnarraybox[#1]}{\@@@FMFPeqnarraybox[#1][b]}}
\def\@@@FMFPeqnarraybox[#1][#2]{\relax\@ifnextchar[{\@@@@FMFPeqnarraybox[#1][#2]}{\@@@@FMFPeqnarraybox[#1][#2][\relax]}}

% #1 = decl; #2 = t,b,c; #3 = width, #4 = col specs
\def\@@@@FMFPeqnarraybox[#1][#2][#3]#4{\@FMFPeqnarrayISinnerfalse % not yet within the lines of the halign
   \@FMFPeqnarraymasterstrutsave% save current master strut values
   \@FMFPeqnarraystrutsize{0pt}{0pt}[\relax]% turn off struts by default
   \@FMFPeqnarrayusemasterstruttrue% use master strut till user asks otherwise
   \FMFPvisiblestrutsfalse% diagnostic mode defaults to off
   % no extra space unless the user specifically requests it
   \lineskip=0pt\relax%
   \lineskiplimit=0pt\relax%
   \baselineskip=\normalbaselineskip\relax%
   \jot=\FMFPnormaljot\relax%
   \mathsurround\z@\relax% no extra spacing around math
   % the default end glues are zero for an \FMFPeqnarraybox
   \edef\@FMFPeqnarraycolSEPdefaultstart{\@FMFPeqnarraycolSEPzero}% default start glue
   \edef\@FMFPeqnarraycolSEPdefaultend{\@FMFPeqnarraycolSEPzero}% default end glue
   \edef\@FMFPeqnarraycolSEPdefaultmid{\@FMFPeqnarraycolSEPzero}% default inter-column glue
   \@advanceFMFPeqncolcntfalse% do not advance the col counter for each col the user uses, 
                              % used in \FMFPeqnarraymulticol and in the preamble build
   \FMFPeqnarrayboxdecl\relax% allow a way for the user to make global overrides
   #1\relax% allow user to override defaults
   \let\\\@FMFPeqnarrayboxcr% replace newline with one that allows optional spacing
   \@FMFPbuildpreamble #4\end\relax% build the preamble and put it into \@FMFPtrantmptoksA
   % add an isolation column to the preamble to stop \\'s {} from getting into the last col
   \ifnum\@FMFPeqnnumcols>0\relax\@FMFPappendtoksA{&}\fi% col separator for those after the first
   \toks0={##}%
   % add the isolation column to the preamble
   \@FMFPappendtoksA{\tabskip\z@skip\bgroup\the\toks0\egroup}% 
   % set the starting tabskip glue as determined by the preamble build
   \tabskip=\@FMFPBPstartglue\relax
   % begin the alignment
   \everycr{}%
   % use only the very first token to determine the positioning
   % this stops some problems when the user uses more than one letter,
   % but is probably not worth the effort
   % \noindent is used as a delimiter
   \def\@FMFPgrabfirstoken##1##2\noindent{\let\@FMFPgrabbedfirstoken=##1}%
   \@FMFPgrabfirstoken#2\relax\relax\noindent
   % \@FMFPgrabbedfirstoken has the first token, the rest are discarded
   % if we need to put things into and hbox and go into math mode, do so now
   \if@FMFPeqnarrayboxHBOXSW \leavevmode \hbox \bgroup $\fi%
   % use the appropriate vbox type
   \if\@FMFPgrabbedfirstoken t\relax\vtop\else\if\@FMFPgrabbedfirstoken c\relax%
   \vcenter\else\vbox\fi\fi\bgroup%
   \@FMFPeqnarrayISinnertrue% commands are now within the lines
   \ifx#3\relax\halign\else\halign to #3\relax\fi%
   \bgroup
   % "exspand" the preamble
   \span\the\@FMFPtrantmptoksA\cr}

% carry strut status and enter the isolation/strut column, 
% exit from math mode if needed, and exit
\def\end@FMFPeqnarraybox{\@FMFPeqnarrayglobalizestrutstatus% carry strut status
&% enter isolation/strut column
\@FMFPeqnarrayinsertstrut% do strut if needed
\@FMFPeqnarraymasterstrutrestore% restore the previous master strut values
% reset the strut system for next FMFPeqnarray
% (sets local strut values back to previous master strut values)
\@FMFPeqnarraystrutreset%
% ensure last line, exit from halign, close vbox
\crcr\egroup\egroup%
% exit from math mode and close hbox if needed
\if@FMFPeqnarrayboxHBOXSW $\egroup\fi}



% FMFPeqnarraybox uses a modifed \\ instead of the plain \cr to
% end rows. This allows for things like \\[vskip amount]
% This "cr" macros are modified versions those for LaTeX2e's eqnarray
% For FMFPeqnarraybox, \\* is the same as \\
% the {\ifnum0=`} braces must be kept away from the last column to avoid
% altering spacing of its math, so we use & to advance to the isolation/strut column
% carry strut status into isolation/strut column
\def\@FMFPeqnarrayboxcr{\@FMFPeqnarrayglobalizestrutstatus% carry strut status
&% enter isolation/strut column
\@FMFPeqnarrayinsertstrut% do strut if needed
% reset the strut system for next line or FMFPeqnarray
\@FMFPeqnarraystrutreset%
{\ifnum0=`}\fi%
\@ifstar{\@FMFPeqnarrayboxYCR}{\@FMFPeqnarrayboxYCR}}

% test and setup the optional argument to \\[]
\def\@FMFPeqnarrayboxYCR{\@testopt\@FMFPeqnarrayboxXCR\z@skip}

% FMFPeqnarraybox does not automatically increase line spacing by \jot
\def\@FMFPeqnarrayboxXCR[#1]{\ifnum0=`{\fi}%
\cr\noalign{\if@FMFPeqnarraystarform\else\vskip\jot\fi\vskip#1\relax}}



% starts the halign preamble build
\def\@FMFPbuildpreamble{\@FMFPtrantmptoksA={}% clear token register
\let\@FMFPBPcurtype=u%current column type is not yet known
\let\@FMFPBPprevtype=s%the previous column type was the start
\let\@FMFPBPnexttype=u%next column type is not yet known
% ensure these are valid
\def\@FMFPBPcurglue={0pt plus 0pt minus 0pt}%
\def\@FMFPBPcurcolname{@FMFPdefault}% name of current column definition
% currently acquired numerically referenced glue
% use a name that is easier to remember
\let\@FMFPBPcurnum=\@FMFPtrantmpcountA%
\@FMFPBPcurnum=0%
% tracks number of columns in the preamble
\@FMFPeqnnumcols=0%
% record the default end glues
\edef\@FMFPBPstartglue{\@FMFPeqnarraycolSEPdefaultstart}%
\edef\@FMFPBPendglue{\@FMFPeqnarraycolSEPdefaultend}%
% now parse the user's column specifications
\@@FMFPbuildpreamble}


% parses and builds the halign preamble
\def\@@FMFPbuildpreamble#1#2{\let\@@nextFMFPbuildpreamble=\@@FMFPbuildpreamble%
% use only the very first token to check the end
% \noindent is used as a delimiter as \end can be present here
\def\@FMFPgrabfirstoken##1##2\noindent{\let\@FMFPgrabbedfirstoken=##1}%
\@FMFPgrabfirstoken#1\relax\relax\noindent
\ifx\@FMFPgrabbedfirstoken\end\let\@@nextFMFPbuildpreamble=\@@FMFPfinishpreamble\else%
% identify current and next token type
\@FMFPgetcoltype{#1}{\@FMFPBPcurtype}{1}% current, error on invalid
\@FMFPgetcoltype{#2}{\@FMFPBPnexttype}{0}% next, no error on invalid next
% if curtype is a glue, get the glue def
\if\@FMFPBPcurtype g\@FMFPgetcurglue{#1}{\@FMFPBPcurglue}\fi%
% if curtype is a column, get the column def and set the current column name
\if\@FMFPBPcurtype c\@FMFPgetcurcol{#1}\fi%
% if curtype is a numeral, acquire the user defined glue
\if\@FMFPBPcurtype n\@FMFPprocessNcol{#1}\fi%
% process the acquired glue 
\if\@FMFPBPcurtype g\@FMFPprocessGcol\fi%
% process the acquired col 
\if\@FMFPBPcurtype c\@FMFPprocessCcol\fi%
% ready prevtype for next col spec.
\let\@FMFPBPprevtype=\@FMFPBPcurtype%
% be sure and put back the future token(s) as a group
\fi\@@nextFMFPbuildpreamble{#2}}


% executed just after preamble build is completed
% warn about zero cols, and if prevtype type = u, put in end tabskip glue
\def\@@FMFPfinishpreamble#1{\ifnum\@FMFPeqnnumcols<1\relax
\@FMFPclspkgerror{No column specifiers declared for FMFPeqnarray}%
{At least one column type must be declared for each FMFPeqnarray.}%
\fi%num cols less than 1
%if last type undefined, set default end tabskip glue
\if\@FMFPBPprevtype u\@FMFPappendtoksA{\tabskip=\@FMFPBPendglue}\fi}


% Identify and return the column specifier's type code
\def\@FMFPgetcoltype#1#2#3{%
% use only the very first token to determine the type
% \noindent is used as a delimiter as \end can be present here
\def\@FMFPgrabfirstoken##1##2\noindent{\let\@FMFPgrabbedfirstoken=##1}%
\@FMFPgrabfirstoken#1\relax\relax\noindent
% \@FMFPgrabfirstoken has the first token, the rest are discarded
% n = number
% g = glue (any other char in catagory 12)
% c = letter
% e = \end
% u = undefined
% third argument: 0 = no error message, 1 = error on invalid char
\let#2=u\relax% assume invalid until know otherwise
\ifx\@FMFPgrabbedfirstoken\end\let#2=e\else
\ifcat\@FMFPgrabbedfirstoken\relax\else% screen out control sequences
\if0\@FMFPgrabbedfirstoken\let#2=n\else
\if1\@FMFPgrabbedfirstoken\let#2=n\else
\if2\@FMFPgrabbedfirstoken\let#2=n\else
\if3\@FMFPgrabbedfirstoken\let#2=n\else
\if4\@FMFPgrabbedfirstoken\let#2=n\else
\if5\@FMFPgrabbedfirstoken\let#2=n\else
\if6\@FMFPgrabbedfirstoken\let#2=n\else
\if7\@FMFPgrabbedfirstoken\let#2=n\else
\if8\@FMFPgrabbedfirstoken\let#2=n\else
\if9\@FMFPgrabbedfirstoken\let#2=n\else
\ifcat,\@FMFPgrabbedfirstoken\let#2=g\relax
\else\ifcat a\@FMFPgrabbedfirstoken\let#2=c\relax\fi\fi\fi\fi\fi\fi\fi\fi\fi\fi\fi\fi\fi\fi
\if#2u\relax
\if0\noexpand#3\relax\else\@FMFPclspkgerror{Invalid character in column specifications}%
{Only letters, numerals and certain other symbols are allowed \MessageBreak
as FMFPeqnarray column specifiers.}\fi\fi}


% identify the current letter referenced column
% if invalid, use a default column
\def\@FMFPgetcurcol#1{\expandafter\ifx\csname @FMFPeqnarraycolDEF#1\endcsname\@FMFPeqnarraycolisdefined%
\def\@FMFPBPcurcolname{#1}\else% invalid column name
\@FMFPclspkgerror{Invalid column type "#1" in column specifications.\MessageBreak
Using a default centering column instead}%
{You must define FMFPeqnarray column types before use.}%
\def\@FMFPBPcurcolname{@FMFPdefault}\fi}


% identify and return the predefined (punctuation) glue value
\def\@FMFPgetcurglue#1#2{%
% ! = \! (neg small)  -0.16667em (-3/18 em)
% , = \, (small)       0.16667em ( 3/18 em)
% : = \: (med)         0.22222em ( 4/18 em)
% ; = \; (large)       0.27778em ( 5/18 em)
% ' = \quad            1em
% " = \qquad           2em
% . = 0.5\arraycolsep
% / = \arraycolsep
% ? = 2\arraycolsep
% * = 1fil
% + = \@FMFPeqnarraycolSEPcenter
% - = \@FMFPeqnarraycolSEPzero
% Note that all em values are referenced to the math font (textfont2) fontdimen6
% value for 1em.
% 
% use only the very first token to determine the type
% this prevents errant tokens from getting in the main text
% \noindent is used as a delimiter here
\def\@FMFPgrabfirstoken##1##2\noindent{\let\@FMFPgrabbedfirstoken=##1}%
\@FMFPgrabfirstoken#1\relax\relax\noindent
% get the math font 1em value
% LaTeX2e's NFSS2 does not preload the fonts, but \FMFPeqnarray needs
% to gain access to the math (\textfont2) font's spacing parameters.
% So we create a bogus box here that uses the math font to ensure
% that \textfont2 is loaded and ready. If this is not done,
% the \textfont2 stuff here may not work.
% Thanks to Bernd Raichle for his 1997 post on this topic.
{\setbox0=\hbox{$\displaystyle\relax$}}%
% fontdimen6 has the width of 1em (a quad).
\@FMFPtrantmpdimenA=\fontdimen6\textfont2\relax%
% identify the glue value based on the first token
% we discard anything after the first
\if!\@FMFPgrabbedfirstoken\@FMFPtrantmpdimenA=-0.16667\@FMFPtrantmpdimenA\edef#2{\the\@FMFPtrantmpdimenA}\else
\if,\@FMFPgrabbedfirstoken\@FMFPtrantmpdimenA=0.16667\@FMFPtrantmpdimenA\edef#2{\the\@FMFPtrantmpdimenA}\else
\if:\@FMFPgrabbedfirstoken\@FMFPtrantmpdimenA=0.22222\@FMFPtrantmpdimenA\edef#2{\the\@FMFPtrantmpdimenA}\else
\if;\@FMFPgrabbedfirstoken\@FMFPtrantmpdimenA=0.27778\@FMFPtrantmpdimenA\edef#2{\the\@FMFPtrantmpdimenA}\else
\if'\@FMFPgrabbedfirstoken\@FMFPtrantmpdimenA=1\@FMFPtrantmpdimenA\edef#2{\the\@FMFPtrantmpdimenA}\else
\if"\@FMFPgrabbedfirstoken\@FMFPtrantmpdimenA=2\@FMFPtrantmpdimenA\edef#2{\the\@FMFPtrantmpdimenA}\else
\if.\@FMFPgrabbedfirstoken\@FMFPtrantmpdimenA=0.5\arraycolsep\edef#2{\the\@FMFPtrantmpdimenA}\else
\if/\@FMFPgrabbedfirstoken\edef#2{\the\arraycolsep}\else
\if?\@FMFPgrabbedfirstoken\@FMFPtrantmpdimenA=2\arraycolsep\edef#2{\the\@FMFPtrantmpdimenA}\else
\if *\@FMFPgrabbedfirstoken\edef#2{0pt plus 1fil minus 0pt}\else
\if+\@FMFPgrabbedfirstoken\edef#2{\@FMFPeqnarraycolSEPcenter}\else
\if-\@FMFPgrabbedfirstoken\edef#2{\@FMFPeqnarraycolSEPzero}\else
\edef#2{\@FMFPeqnarraycolSEPzero}%
\@FMFPclspkgerror{Invalid predefined inter-column glue type "#1" in\MessageBreak
column specifications. Using a default value of\MessageBreak
0pt instead}%
{Only !,:;'"./?*+ and - are valid predefined glue types in the\MessageBreak 
FMFPeqnarray column specifications.}\fi\fi\fi\fi\fi\fi\fi\fi\fi\fi\fi\fi}



% process a numerical digit from the column specification
% and look up the corresponding user defined glue value
% can transform current type from n to g or a as the user defined glue is acquired
\def\@FMFPprocessNcol#1{\if\@FMFPBPprevtype g%
\@FMFPclspkgerror{Back-to-back inter-column glue specifiers in column\MessageBreak
specifications. Ignoring consecutive glue specifiers\MessageBreak
after the first}%
{You cannot have two or more glue types next to each other\MessageBreak 
in the FMFPeqnarray column specifications.}%
\let\@FMFPBPcurtype=a% abort this glue, future digits will be discarded
\@FMFPBPcurnum=0\relax%
\else% if we previously aborted a glue
\if\@FMFPBPprevtype a\@FMFPBPcurnum=0\let\@FMFPBPcurtype=a%maintain digit abortion
\else%acquire this number
% save the previous type before the numerical digits started
\if\@FMFPBPprevtype n\else\let\@FMFPBPprevsavedtype=\@FMFPBPprevtype\fi%
\multiply\@FMFPBPcurnum by 10\relax%
\advance\@FMFPBPcurnum by #1\relax% add in number, \relax is needed to stop TeX's number scan
\if\@FMFPBPnexttype n\else%close acquisition
\expandafter\ifx\csname @FMFPeqnarraycolSEPDEF\expandafter\romannumeral\number\@FMFPBPcurnum\endcsname\@FMFPeqnarraycolisdefined%
\edef\@FMFPBPcurglue{\csname @FMFPeqnarraycolSEP\expandafter\romannumeral\number\@FMFPBPcurnum\endcsname}%
\else%user glue not defined
\@FMFPclspkgerror{Invalid user defined inter-column glue type "\number\@FMFPBPcurnum" in\MessageBreak
column specifications. Using a default value of\MessageBreak
0pt instead}%
{You must define all FMFPeqnarray numerical inter-column glue types via\MessageBreak
\string\FMFPeqnarraydefcolsep \space before they are used in column specifications.}%
\edef\@FMFPBPcurglue{\@FMFPeqnarraycolSEPzero}%
\fi% glue defined or not
\let\@FMFPBPcurtype=g% change the type to reflect the acquired glue
\let\@FMFPBPprevtype=\@FMFPBPprevsavedtype% restore the prev type before this number glue
\@FMFPBPcurnum=0\relax%ready for next acquisition
\fi%close acquisition, get glue
\fi%discard or acquire number
\fi%prevtype glue or not
}


% process an acquired glue
% add any acquired column/glue pair to the preamble
\def\@FMFPprocessGcol{\if\@FMFPBPprevtype a\let\@FMFPBPcurtype=a%maintain previous glue abortions
\else
% if this is the start glue, save it, but do nothing else 
% as this is not used in the preamble, but before
\if\@FMFPBPprevtype s\edef\@FMFPBPstartglue{\@FMFPBPcurglue}%
\else%not the start glue
\if\@FMFPBPprevtype g%ignore if back to back glues
\@FMFPclspkgerror{Back-to-back inter-column glue specifiers in column\MessageBreak
specifications. Ignoring consecutive glue specifiers\MessageBreak
after the first}%
{You cannot have two or more glue types next to each other\MessageBreak 
in the FMFPeqnarray column specifications.}%
\let\@FMFPBPcurtype=a% abort this glue
\else% not a back to back glue
\if\@FMFPBPprevtype c\relax% if the previoustype was a col, add column/glue pair to preamble
\ifnum\@FMFPeqnnumcols>0\relax\@FMFPappendtoksA{&}\fi
\toks0={##}%
% make preamble advance col counter if this environment needs this
\if@advanceFMFPeqncolcnt\@FMFPappendtoksA{\global\advance\@FMFPeqncolcnt by 1\relax}\fi
% insert the column defintion into the preamble, being careful not to expand
% the column definition
\@FMFPappendtoksA{\tabskip=\@FMFPBPcurglue}%
\@FMFPappendNOEXPANDtoksA{\begingroup\csname @FMFPeqnarraycolPRE}%
\@FMFPappendtoksA{\@FMFPBPcurcolname}%
\@FMFPappendNOEXPANDtoksA{\endcsname}%
\@FMFPappendtoksA{\the\toks0}%
\@FMFPappendNOEXPANDtoksA{\relax\relax\relax\relax\relax%
\relax\relax\relax\relax\relax\csname @FMFPeqnarraycolPOST}%
\@FMFPappendtoksA{\@FMFPBPcurcolname}%
\@FMFPappendNOEXPANDtoksA{\endcsname\relax\relax\relax\relax\relax%
\relax\relax\relax\relax\relax\endgroup}%
\advance\@FMFPeqnnumcols by 1\relax%one more column in the preamble
\else% error: non-start glue with no pending column
\@FMFPclspkgerror{Inter-column glue specifier without a prior column\MessageBreak
type in the column specifications. Ignoring this glue\MessageBreak 
specifier}%
{Except for the first and last positions, glue can be placed only\MessageBreak
between column types.}%
\let\@FMFPBPcurtype=a% abort this glue
\fi% previous was a column
\fi% back-to-back glues
\fi% is start column glue
\fi% prev type not a
}


% process an acquired letter referenced column and, if necessary, add it to the preamble
\def\@FMFPprocessCcol{\if\@FMFPBPnexttype g\else
\if\@FMFPBPnexttype n\else
% we have a column followed by something other than a glue (or numeral glue)
% so we must add this column to the preamble now
\ifnum\@FMFPeqnnumcols>0\relax\@FMFPappendtoksA{&}\fi%col separator for those after the first
\if\@FMFPBPnexttype e\@FMFPappendtoksA{\tabskip=\@FMFPBPendglue\relax}\else%put in end glue
\@FMFPappendtoksA{\tabskip=\@FMFPeqnarraycolSEPdefaultmid\relax}\fi% or default mid glue
\toks0={##}%
% make preamble advance col counter if this environment needs this
\if@advanceFMFPeqncolcnt\@FMFPappendtoksA{\global\advance\@FMFPeqncolcnt by 1\relax}\fi
% insert the column definition into the preamble, being careful not to expand
% the column definition
\@FMFPappendNOEXPANDtoksA{\begingroup\csname @FMFPeqnarraycolPRE}%
\@FMFPappendtoksA{\@FMFPBPcurcolname}%
\@FMFPappendNOEXPANDtoksA{\endcsname}%
\@FMFPappendtoksA{\the\toks0}%
\@FMFPappendNOEXPANDtoksA{\relax\relax\relax\relax\relax%
\relax\relax\relax\relax\relax\csname @FMFPeqnarraycolPOST}%
\@FMFPappendtoksA{\@FMFPBPcurcolname}%
\@FMFPappendNOEXPANDtoksA{\endcsname\relax\relax\relax\relax\relax%
\relax\relax\relax\relax\relax\endgroup}%
\advance\@FMFPeqnnumcols by 1\relax%one more column in the preamble
\fi%next type not numeral
\fi%next type not glue
}


%%
%% END OF FMFPeqnarry DEFINITIONS
%%




% set up the running headings, this complex because of all the different
% modes FMFPtran supports
\if@twoside
 \ifCLASSOPTIONtechnote
   \def\ps@headings{%
       \def\@oddhead{\hbox{}\scriptsize\leftmark \hfil \thepage}
       \def\@evenhead{\scriptsize\thepage \hfil \leftmark\hbox{}}
       \ifCLASSOPTIONdraftcls
            \ifCLASSOPTIONdraftclsnofoot
               \def\@oddfoot{}\def\@evenfoot{}%
            \else
               \def\@oddfoot{\scriptsize\@date\hfil DRAFT}
               \def\@evenfoot{\scriptsize DRAFT\hfil\@date}
            \fi
       \else
            \def\@oddfoot{}\def\@evenfoot{}
       \fi}
 \else % not a technote
   \def\ps@headings{%
       \ifCLASSOPTIONconference
        \def\@oddhead{}
        \def\@evenhead{}
       \else
        \def\@oddhead{\hbox{}\scriptsize\rightmark \hfil \thepage}
        \def\@evenhead{\scriptsize\thepage \hfil \leftmark\hbox{}}
       \fi
       \ifCLASSOPTIONdraftcls
            \def\@oddhead{\hbox{}\scriptsize\rightmark \hfil \thepage}
            \def\@evenhead{\scriptsize\thepage \hfil \leftmark\hbox{}}
            \ifCLASSOPTIONdraftclsnofoot
               \def\@oddfoot{}\def\@evenfoot{}%
            \else
               \def\@oddfoot{\scriptsize\@date\hfil DRAFT}
               \def\@evenfoot{\scriptsize DRAFT\hfil\@date}
            \fi
       \else
            \def\@oddfoot{}\def\@evenfoot{}%
       \fi}
 \fi
\else % single side
\def\ps@headings{%
    \ifCLASSOPTIONconference
     \def\@oddhead{}
     \def\@evenhead{}
    \else
     \def\@oddhead{\hbox{}\scriptsize\leftmark \hfil \thepage}
     \def\@evenhead{}
    \fi
    \ifCLASSOPTIONdraftcls
          \def\@oddhead{\hbox{}\scriptsize\leftmark \hfil \thepage}
          \def\@evenhead{}
          \ifCLASSOPTIONdraftclsnofoot
             \def\@oddfoot{}
          \else
             \def\@oddfoot{\scriptsize \@date \hfil DRAFT}
          \fi
    \else
         \def\@oddfoot{}
    \fi
    \def\@evenfoot{}}
\fi


% title page style
\def\ps@FMFPtitlepagestyle{\def\@oddfoot{}\def\@evenfoot{}%
\ifCLASSOPTIONconference
   \def\@oddhead{}%
   \def\@evenhead{}%
\else
   \def\@oddhead{\hbox{}\scriptsize\leftmark \hfil \thepage}%
   \def\@evenhead{\scriptsize\thepage \hfil \leftmark\hbox{}}%
\fi
\ifCLASSOPTIONdraftcls
   \def\@oddhead{\hbox{}\scriptsize\leftmark \hfil \thepage}%
   \def\@evenhead{\scriptsize\thepage \hfil \leftmark\hbox{}}%
   \ifCLASSOPTIONdraftclsnofoot\else
      \def\@oddfoot{\scriptsize \@date\hfil DRAFT}%
      \def\@evenfoot{\scriptsize DRAFT\hfil \@date}%
   \fi
\else
   % all non-draft mode footers
   \if@FMFPusingpubid
      % for title pages that are using a pubid
      % do not repeat pubid if using peer review option
      \ifCLASSOPTIONpeerreview
      \else
         \footskip 0pt%
         \ifCLASSOPTIONcompsoc
           \def\@oddfoot{\hss\normalfont\scriptsize\raisebox{-1.5\@FMFPnormalsizeunitybaselineskip}[0ex][0ex]{\@FMFPpubid}\hss}%
           \def\@evenfoot{\hss\normalfont\scriptsize\raisebox{-1.5\@FMFPnormalsizeunitybaselineskip}[0ex][0ex]{\@FMFPpubid}\hss}%
         \else
           \def\@oddfoot{\hss\normalfont\footnotesize\raisebox{1.5ex}[1.5ex]{\@FMFPpubid}\hss}%
           \def\@evenfoot{\hss\normalfont\footnotesize\raisebox{1.5ex}[1.5ex]{\@FMFPpubid}\hss}%
         \fi
      \fi
   \fi
\fi}


% peer review cover page style
\def\ps@FMFPpeerreviewcoverpagestyle{%
\def\@oddhead{}\def\@evenhead{}%
\def\@oddfoot{}\def\@evenfoot{}%
\ifCLASSOPTIONdraftcls
   \ifCLASSOPTIONdraftclsnofoot\else
      \def\@oddfoot{\scriptsize \@date\hfil DRAFT}%
      \def\@evenfoot{\scriptsize DRAFT\hfil \@date}%
   \fi
\else
   % non-draft mode footers
   \if@FMFPusingpubid
      \footskip 0pt%
      \ifCLASSOPTIONcompsoc
        \def\@oddfoot{\hss\normalfont\scriptsize\raisebox{-1.5\@FMFPnormalsizeunitybaselineskip}[0ex][0ex]{\@FMFPpubid}\hss}%
        \def\@evenfoot{\hss\normalfont\scriptsize\raisebox{-1.5\@FMFPnormalsizeunitybaselineskip}[0ex][0ex]{\@FMFPpubid}\hss}%
      \else
        \def\@oddfoot{\hss\normalfont\footnotesize\raisebox{1.5ex}[1.5ex]{\@FMFPpubid}\hss}%
        \def\@evenfoot{\hss\normalfont\footnotesize\raisebox{1.5ex}[1.5ex]{\@FMFPpubid}\hss}%
      \fi
   \fi
\fi}


% start with empty headings
\def\rightmark{}\def\leftmark{}


%% Defines the command for putting the header. \footernote{TEXT} is the same
%% as \markboth{TEXT}{TEXT}. 
%% Note that all the text is forced into uppercase, if you have some text
%% that needs to be in lower case, for instance et. al., then either manually
%% set \leftmark and \rightmark or use \MakeLowercase{et. al.} within the
%% arguments to \markboth.
\def\markboth#1#2{\def\leftmark{\@FMFPcompsoconly{\sffamily}\MakeUppercase{#1}}%
\def\rightmark{\@FMFPcompsoconly{\sffamily}\MakeUppercase{#2}}}
\def\footernote#1{\markboth{#1}{#1}}

\def\today{\ifcase\month\or
    January\or February\or March\or April\or May\or June\or
    July\or August\or September\or October\or November\or December\fi
    \space\number\day, \number\year}




%% CITATION AND BIBLIOGRAPHY COMMANDS
%% 
%% V1.6 no longer supports the older, nonstandard \shortcite and \citename setup stuff
% 
% 
% Modify Latex2e \@citex to separate citations with "], ["
\def\@citex[#1]#2{%
  \let\@citea\@empty
  \@cite{\@for\@citeb:=#2\do
    {\@citea\def\@citea{], [}%
     \edef\@citeb{\expandafter\@firstofone\@citeb\@empty}%
     \if@filesw\immediate\write\@auxout{\string\citation{\@citeb}}\fi
     \@ifundefined{b@\@citeb}{\mbox{\reset@font\bfseries ?}%
       \G@refundefinedtrue
       \@latex@warning
         {Citation `\@citeb' on page \thepage \space undefined}}%
       {\hbox{\csname b@\@citeb\endcsname}}}}{#1}}

% V1.6 we create hooks for the optional use of Donald Arseneau's
% cite.sty package. cite.sty is "smart" and will notice that the
% following format controls are already defined and will not
% redefine them. The result will be the proper sorting of the
% citation numbers and auto detection of 3 or more entry "ranges" -
% all in FMFP style:  [1], [2], [5]--[7], [12]
% This also allows for an optional note, i.e., \cite[mynote]{..}.
% If the \cite with note has more than one reference, the note will
% be applied to the last of the listed references. It is generally
% desired that if a note is given, only one reference is listed in
% that \cite.
% Thanks to Mr. Arseneau for providing the required format arguments
% to produce the FMFP style.
\def\citepunct{], [}
\def\citedash{]--[}

% V1.7 default to using same font for urls made by url.sty
\AtBeginDocument{\csname url@samestyle\endcsname}

% V1.6 class files should always provide these
\def\newblock{\hskip .11em\@plus.33em\@minus.07em}
\let\@openbib@code\@empty


% Provide support for the control entries of FMFPtran.bst V1.00 and later.
% V1.7 optional argument allows for a different aux file to be specified in
% order to handle multiple bibliographies. For example, with multibib.sty:
% \newcites{sec}{Secondary Literature}
% \bstctlcite[@auxoutsec]{BSTcontrolhak}
\def\bstctlcite{\@ifnextchar[{\@bstctlcite}{\@bstctlcite[@auxout]}}
\def\@bstctlcite[#1]#2{\@bsphack
  \@for\@citeb:=#2\do{%
    \edef\@citeb{\expandafter\@firstofone\@citeb}%
    \if@filesw\immediate\write\csname #1\endcsname{\string\citation{\@citeb}}\fi}%
  \@esphack}

% V1.6 provide a way for a user to execute a command just before 
% a given reference number - used to insert a \newpage to balance
% the columns on the last page
\edef\@FMFPtriggerrefnum{0}   % the default of zero means that
                              % the command is not executed
\def\@FMFPtriggercmd{\newpage}

% allow the user to alter the triggered command
\long\def\FMFPtriggercmd#1{\long\def\@FMFPtriggercmd{#1}}

% allow user a way to specify the reference number just before the
% command is executed
\def\FMFPtriggeratref#1{\@FMFPtrantmpcountA=#1%
\edef\@FMFPtriggerrefnum{\the\@FMFPtrantmpcountA}}%

% trigger command at the given reference
\def\@FMFPbibitemprefix{\@FMFPtrantmpcountA=\@FMFPtriggerrefnum\relax%
\advance\@FMFPtrantmpcountA by -1\relax%
\ifnum\c@enumiv=\@FMFPtrantmpcountA\relax\@FMFPtriggercmd\relax\fi}


\def\@biblabel#1{[#1]}

% compsoc journals left align the reference numbers
\@FMFPcompsocnotconfonly{\def\@biblabel#1{[#1]\hfill}}

% controls bib item spacing
\def\FMFPbibitemsep{2.5pt plus .5pt}

\@FMFPcompsocconfonly{\def\FMFPbibitemsep{1\baselineskip plus 0.25\baselineskip minus 0.25\baselineskip}}


\def\thebibliography#1{\section*{\refname}%
    \addcontentsline{toc}{section}{\refname}%
    % V1.6 add some rubber space here and provide a command trigger
    \footnotesize\@FMFPcompsocconfonly{\small}\vskip 0.3\baselineskip plus 0.1\baselineskip minus 0.1\baselineskip%
    \list{\@biblabel{\@arabic\c@enumiv}}%
    {\settowidth\labelwidth{\@biblabel{#1}}%
    \leftmargin\labelwidth
    \labelsep 1em
    \advance\leftmargin\labelsep\relax
    \itemsep \FMFPbibitemsep\relax
    \usecounter{enumiv}%
    \let\p@enumiv\@empty
    \renewcommand\theenumiv{\@arabic\c@enumiv}}%
    \let\@FMFPlatexbibitem\bibitem%
    \def\bibitem{\@FMFPbibitemprefix\@FMFPlatexbibitem}%
\def\newblock{\hskip .11em plus .33em minus .07em}%
% originally:
%   \sloppy\clubpenalty4000\widowpenalty4000%
% by adding the \interlinepenalty here, we make it more
% difficult, but not impossible, for LaTeX to break within a reference.
% FMFP almost never breaks a reference (but they do it more often with
% technotes). You may get an underfull vbox warning around the bibliography, 
% but the final result will be much more like what FMFP will publish. 
% MDS 11/2000
\ifCLASSOPTIONtechnote\sloppy\clubpenalty4000\widowpenalty4000\interlinepenalty100%
\else\sloppy\clubpenalty4000\widowpenalty4000\interlinepenalty500\fi%
    \sfcode`\.=1000\relax}
\let\endthebibliography=\endlist




% TITLE PAGE COMMANDS
% 
% 
% \FMFPmembership is used to produce the sublargesize italic font used to indicate author 
% FMFP membership. compsoc uses a large size sans slant font
\def\FMFPmembership#1{{\@FMFPnotcompsoconly{\sublargesize}\normalfont\@FMFPcompsoconly{\sffamily}\textit{#1}}}
 

% \FMFPauthorrefmark{} produces a footnote type symbol to indicate author affiliation.
% When given an argument of 1 to 9, \FMFPauthorrefmark{} follows the standard LaTeX footnote
% symbol sequence convention. However, for arguments 10 and above, \FMFPauthorrefmark{} 
% reverts to using lower case roman numerals, so it cannot overflow. Do note that you 
% cannot use \footnotemark[] in place of \FMFPauthorrefmark{} within \author as the footnote
% symbols will have been turned off to prevent \thanks from creating footnote marks.
% \FMFPauthorrefmark{} produces a symbol that appears to LaTeX as having zero vertical
% height - this allows for a more compact line packing, but the user must ensure that
% the interline spacing is large enough to prevent \FMFPauthorrefmark{} from colliding
% with the text above.
% V1.7 make this a robust command
\DeclareRobustCommand*{\FMFPauthorrefmark}[1]{\raisebox{0pt}[0pt][0pt]{\textsuperscript{\footnotesize\ensuremath{\ifcase#1\or *\or \dagger\or \ddagger\or%
    \mathsection\or \mathparagraph\or \|\or **\or \dagger\dagger%
    \or \ddagger\ddagger \else\textsuperscript{\expandafter\romannumeral#1}\fi}}}}


% FONT CONTROLS AND SPACINGS FOR CONFERENCE MODE AUTHOR NAME AND AFFILIATION BLOCKS
% 
% The default font styles for the author name and affiliation blocks (confmode)
\def\@FMFPauthorblockNstyle{\small\@FMFPcompsocnotconfonly{\sffamily} 
% Changed by AJC 30.4..2018
\sublargesize\@FMFPcompsocconfonly{\large}}
\def\@FMFPauthorblockAstyle{\small\@FMFPcompsocnotconfonly{\sffamily}\@FMFPcompsocconfonly{\itshape}\small\@FMFPcompsocconfonly{\large}}
% The default if the user does not use an author block
\def\@FMFPauthordefaulttextstyle{\small\@FMFPcompsocnotconfonly{\sffamily}\sublargesize}

% spacing from title (or special paper notice) to author name blocks (confmode)
% can be negative
\def\@FMFPauthorblockconfadjspace{-0.25em}
% compsoc conferences need more space here
\@FMFPcompsocconfonly{\def\@FMFPauthorblockconfadjspace{0.75\@FMFPnormalsizeunitybaselineskip}}
\ifCLASSOPTIONconference\def\@FMFPauthorblockconfadjspace{20pt}\fi

% spacing between name and affiliation blocks (confmode)
% This can be negative.
% FMFP doesn't want any added spacing here, but I will leave these
% controls in place in case they ever change their mind.
% Personally, I like 0.75ex.
%\def\@FMFPauthorblockNtopspace{0.75ex}
%\def\@FMFPauthorblockAtopspace{0.75ex}
\def\@FMFPauthorblockNtopspace{0.0ex}
\def\@FMFPauthorblockAtopspace{0.0ex}
% baseline spacing within name and affiliation blocks (confmode)
% must be positive, spacings below certain values will make 
% the position of line of text sensitive to the contents of the
% line above it i.e., whether or not the prior line has descenders, 
% subscripts, etc. For this reason it is a good idea to keep
% these above 2.6ex
\def\@FMFPauthorblockNinterlinespace{2.6ex}
\def\@FMFPauthorblockAinterlinespace{2.75ex}

% This tracks the required strut size.
% See the \@FMFPauthorhalign command for the actual default value used.
\def\@FMFPauthorblockXinterlinespace{2.7ex}

% variables to retain font size and style across groups
% values given here have no effect as they will be overwritten later
\gdef\@FMFPSAVESTATEfontsize{10}
\gdef\@FMFPSAVESTATEfontbaselineskip{12}
\gdef\@FMFPSAVESTATEfontencoding{OT1}
\gdef\@FMFPSAVESTATEfontfamily{ptm}
\gdef\@FMFPSAVESTATEfontseries{m}
\gdef\@FMFPSAVESTATEfontshape{n}

% saves the current font attributes
\def\@FMFPcurfontSAVE{\global\let\@FMFPSAVESTATEfontsize\f@size%
\global\let\@FMFPSAVESTATEfontbaselineskip\f@baselineskip%
\global\let\@FMFPSAVESTATEfontencoding\f@encoding%
\global\let\@FMFPSAVESTATEfontfamily\f@family%
\global\let\@FMFPSAVESTATEfontseries\f@series%
\global\let\@FMFPSAVESTATEfontshape\f@shape}

% restores the saved font attributes
\def\@FMFPcurfontRESTORE{\fontsize{\@FMFPSAVESTATEfontsize}{\@FMFPSAVESTATEfontbaselineskip}%
\fontencoding{\@FMFPSAVESTATEfontencoding}%
\fontfamily{\@FMFPSAVESTATEfontfamily}%
\fontseries{\@FMFPSAVESTATEfontseries}%
\fontshape{\@FMFPSAVESTATEfontshape}%
\selectfont}


% variable to indicate if the current block is the first block in the column
\newif\if@FMFPprevauthorblockincol   \@FMFPprevauthorblockincolfalse


% the command places a strut with height and depth = \@FMFPauthorblockXinterlinespace
% we use this technique to have complete manual control over the spacing of the lines
% within the halign environment.
% We set the below baseline portion at 30%, the above
% baseline portion at 70% of the total length.
% Responds to changes in the document's \baselinestretch
\def\@FMFPauthorstrutrule{\@FMFPtrantmpdimenA\@FMFPauthorblockXinterlinespace%
\@FMFPtrantmpdimenA=\baselinestretch\@FMFPtrantmpdimenA%
\rule[-0.3\@FMFPtrantmpdimenA]{0pt}{\@FMFPtrantmpdimenA}}


% blocks to hold the authors' names and affilations. 
% Makes formatting easy for conferences
%
% use real definitions in conference mode
% name block
\def\FMFPauthorblockN#1{\relax\@FMFPauthorblockNstyle% set the default text style
\gdef\@FMFPauthorblockXinterlinespace{0pt}% disable strut for spacer row
% the \expandafter hides the \cr in conditional tex, see the array.sty docs
% for details, probably not needed here as the \cr is in a macro
% do a spacer row if needed
\if@FMFPprevauthorblockincol\expandafter\@FMFPauthorblockNtopspaceline\fi
\global\@FMFPprevauthorblockincoltrue% we now have a block in this column
%restore the correct strut value
\gdef\@FMFPauthorblockXinterlinespace{\@FMFPauthorblockNinterlinespace}%
% input the author names
#1%
% end the row if the user did not already
\crcr}
% spacer row for names
\def\@FMFPauthorblockNtopspaceline{\cr\noalign{\vskip\@FMFPauthorblockNtopspace}}
%
% affiliation block
\def\FMFPauthorblockA#1{\relax\@FMFPauthorblockAstyle% set the default text style
\gdef\@FMFPauthorblockXinterlinespace{0pt}%disable strut for spacer row
% the \expandafter hides the \cr in conditional tex, see the array.sty docs
% for details, probably not needed here as the \cr is in a macro
% do a spacer row if needed
\if@FMFPprevauthorblockincol\expandafter\@FMFPauthorblockAtopspaceline\fi
\global\@FMFPprevauthorblockincoltrue% we now have a block in this column
%restore the correct strut value
\gdef\@FMFPauthorblockXinterlinespace{\@FMFPauthorblockAinterlinespace}%
% input the author affiliations
#1%
% end the row if the user did not already
\crcr}
% spacer row for affiliations
\def\@FMFPauthorblockAtopspaceline{\cr\noalign{\vskip\@FMFPauthorblockAtopspace}}


% allow papers to compile even if author blocks are used in modes other
% than conference or peerreviewca. For such cases, we provide dummy blocks.
\ifCLASSOPTIONconference
\else
   \ifCLASSOPTIONpeerreviewca\else
      % not conference or peerreviewca mode
      \def\FMFPauthorblockN#1{#1}%
      \def\FMFPauthorblockA#1{#1}%
   \fi
\fi



% we provide our own halign so as not to have to depend on tabular
\def\@FMFPauthorhalign{\@FMFPauthordefaulttextstyle% default text style
   \lineskip=0pt\relax% disable line spacing
   \lineskiplimit=0pt\relax%
   \baselineskip=0pt\relax%
   \@FMFPcurfontSAVE% save the current font
   \mathsurround\z@\relax% no extra spacing around math
   \let\\\@FMFPauthorhaligncr% replace newline with halign friendly one
   \tabskip=0pt\relax% no column spacing
   \everycr{}% ensure no problems here
   \@FMFPprevauthorblockincolfalse% no author blocks yet
   \def\@FMFPauthorblockXinterlinespace{2.7ex}% default interline space
   \vtop\bgroup%vtop box
   \halign\bgroup&\relax\hfil\@FMFPcurfontRESTORE\relax ##\relax
   \hfil\@FMFPcurfontSAVE\@FMFPauthorstrutrule\cr}

% ensure last line, exit from halign, close vbox
\def\end@FMFPauthorhalign{\crcr\egroup\egroup}

% handle bogus star form
\def\@FMFPauthorhaligncr{{\ifnum0=`}\fi\@ifstar{\@@FMFPauthorhaligncr}{\@@FMFPauthorhaligncr}}

% test and setup the optional argument to \\[]
\def\@@FMFPauthorhaligncr{\@testopt\@@@FMFPauthorhaligncr\z@skip}

% end the line and do the optional spacer
\def\@@@FMFPauthorhaligncr[#1]{\ifnum0=`{\fi}\cr\noalign{\vskip#1\relax}}



% flag to prevent multiple \and warning messages
\newif\if@FMFPWARNand
\@FMFPWARNandtrue

% if in conference or peerreviewca modes, we support the use of \and as \author is a
%  environment, otherwise we warn the user that \and is invalid
% outside of conference or peerreviewca modes.
\def\and{\relax} % provide a bogus \and that we will then override

\renewcommand{\and}[1][\relax]{\if@FMFPWARNand\typeout{** WARNING: \noexpand\and is valid only
                               when in conference or peerreviewca}\typeout{modes (line \the\inputlineno).}\fi\global\@FMFPWARNandfalse}

\ifCLASSOPTIONconference%
\renewcommand{\and}[1][\hfill]{\end{@FMFPauthorhalign}#1\begin{@FMFPauthorhalign}}%
\fi
\ifCLASSOPTIONpeerreviewca
\renewcommand{\and}[1][\hfill]{\end{@FMFPauthorhalign}#1\begin{@FMFPauthorhalign}}%
\fi


% page clearing command
% based on LaTeX2e's \cleardoublepage, but allows different page styles
% for the inserted blank pages
\def\@FMFPcleardoublepage#1{\clearpage\if@twoside\ifodd\c@page\else
\hbox{}\thispagestyle{#1}\newpage\if@twocolumn\hbox{}\thispagestyle{#1}\newpage\fi\fi\fi}


% user command to invoke the title page
\def\maketitle{\par%
  \begingroup%
  \normalfont%
  \def\thefootnote{}%  the \thanks{} mark type is empty
  \def\footnotemark{}% and kill space from \thanks within author
  \let\@makefnmark\relax% V1.7, must *really* kill footnotemark to remove all \textsuperscript spacing as well.
  \footnotesize%       equal spacing between thanks lines
  \footnotesep 0.7\baselineskip%see global setting of \footnotesep for more info
  % V1.7 disable \thanks note indention for compsoc
  \@FMFPcompsoconly{\long\def\@makefntext##1{\parindent 1em\noindent\hbox{\@makefnmark}##1}}%
  \normalsize%
  \ifCLASSOPTIONpeerreview
     \newpage\global\@topnum\z@ \@maketitle\@FMFPstatictitlevskip\@FMFPaftertitletext%
     \thispagestyle{FMFPpeerreviewcoverpagestyle}\@thanks%
  \else
     \if@twocolumn%
        \ifCLASSOPTIONtechnote%
           \newpage\global\@topnum\z@ \@maketitle\@FMFPstatictitlevskip\@FMFPaftertitletext%
        \else
           \twocolumn[\@maketitle\@FMFPstatictitlevskip\@FMFPaftertitletext]%
        \fi
     \else
        \newpage\global\@topnum\z@ \@maketitle\@FMFPstatictitlevskip\@FMFPaftertitletext%
     \fi
     \thispagestyle{FMFPtitlepagestyle}\@thanks%
  \fi
  % pullup page for pubid if used.
  \if@FMFPusingpubid
     \enlargethispage{-\@FMFPpubidpullup}%
  \fi 
  \endgroup
  \setcounter{footnote}{0}\let\maketitle\relax\let\@maketitle\relax
  \gdef\@thanks{}%
  % v1.6b do not clear these as we will need the title again for peer review papers
  % \gdef\@author{}\gdef\@title{}%
  \let\thanks\relax}



% V1.7 parbox to format \@FMFPcompsoctitleabstractindextext
\long\def\@FMFPcompsoctitleabstractindextextbox#1{\parbox{0.915\textwidth}{#1}}

% formats the Title, authors names, affiliations and special paper notice
% THIS IS A CONTROLLED SPACING COMMAND! Do not allow blank lines or unintentional
% spaces to enter the definition - use % at the end of each line
\def\@maketitle{\newpage
\begingroup\centering
\ifCLASSOPTIONtechnote% technotes
   {\bfseries\large\@FMFPcompsoconly{\sffamily}\@title\par}\vskip 1.3em{\lineskip .5em\@FMFPcompsoconly{\sffamily}\@author
   \@FMFPspecialpapernotice\par{\@FMFPcompsoconly{\vskip 1.5em\relax
   \@FMFPcompsoctitleabstractindextextbox{\@FMFPcompsoctitleabstractindextext}\par
   \hfill\@FMFPcompsocdiamondline\hfill\hbox{}\par}}}\relax
\else% not a technote
   \vskip0.2em{\Huge\@FMFPcompsoconly{\sffamily}\@FMFPcompsocconfonly{\normalfont\normalsize\vskip 2\@FMFPnormalsizeunitybaselineskip
   \bfseries\Large}\@title\par}\vskip1.0em\par%
   % V1.6 handle \author differently if in conference mode
   \ifCLASSOPTIONconference%
      {\@FMFPspecialpapernotice\mbox{}\vskip\@FMFPauthorblockconfadjspace%
       \mbox{}\hfill\begin{@FMFPauthorhalign}\@author\end{@FMFPauthorhalign}\hfill\mbox{}\par}\relax
   \else% peerreviewca, peerreview or journal
      \ifCLASSOPTIONpeerreviewca
         % peerreviewca handles author names just like conference mode
         {\@FMFPcompsoconly{\sffamily}\@FMFPspecialpapernotice\mbox{}\vskip\@FMFPauthorblockconfadjspace%
          \mbox{}\hfill\begin{@FMFPauthorhalign}\@author\end{@FMFPauthorhalign}\hfill\mbox{}\par
          {\@FMFPcompsoconly{\vskip 1.5em\relax
           \@FMFPcompsoctitleabstractindextextbox{\@FMFPcompsoctitleabstractindextext}\par\hfill
           \@FMFPcompsocdiamondline\hfill\hbox{}\par}}}\relax
      \else% journal or peerreview
         {\lineskip.5em\@FMFPcompsoconly{\sffamily}\sublargesize\@author\@FMFPspecialpapernotice\par
          {\@FMFPcompsoconly{\vskip 1.5em\relax
           \@FMFPcompsoctitleabstractindextextbox{\@FMFPcompsoctitleabstractindextext}\par\hfill
           \@FMFPcompsocdiamondline\hfill\hbox{}\par}}}\relax
      \fi
   \fi
\fi\par\endgroup}



% V1.7 Computer Society "diamond line" which follows index terms for nonconference papers
\def\@FMFPcompsocdiamondline{\vrule depth 0pt height 0.5pt width 4cm\hspace{7.5pt}%
\raisebox{-3.5pt}{\fontfamily{pzd}\fontencoding{U}\fontseries{m}\fontshape{n}\fontsize{11}{12}\selectfont\char70}%
\hspace{7.5pt}\vrule depth 0pt height 0.5pt width 4cm\relax}

% V1.7 standard LateX2e \thanks, but with \itshape under compsoc. Also make it a \long\def
% We also need to trigger the one-shot footnote rule
\def\@FMFPtriggeroneshotfootnoterule{\global\@FMFPenableoneshotfootnoteruletrue}


\long\def\thanks#1{\footnotemark
    \protected@xdef\@thanks{\@thanks
        \protect\footnotetext[\the\c@footnote]{\@FMFPcompsoconly{\itshape
        \protect\@FMFPtriggeroneshotfootnoterule\relax}\ignorespaces#1}}}
\let\@thanks\@empty

% V1.7 allow \author to contain \par's. This is needed to allow \thanks to contain \par.
\long\def\author#1{\gdef\@author{#1}}


% in addition to setting up FMFPitemize, we need to remove a baselineskip space above and
% below it because \list's \pars introduce blank lines because of the footnote struts.
\def\@FMFPsetupcompsocitemizelist{\def\labelitemi{$\bullet$}%
\setlength{\FMFPlabelindent}{0pt}\setlength{\parskip}{0pt}%
\setlength{\partopsep}{0pt}\setlength{\topsep}{0.5\baselineskip}\vspace{-1\baselineskip}\relax}


% flag for fake non-compsoc \FMFPcompsocthanksitem - prevents line break on very first item
\newif\if@FMFPbreakcompsocthanksitem \@FMFPbreakcompsocthanksitemfalse

\ifCLASSOPTIONcompsoc
% V1.7 compsoc bullet item \thanks
% also, we need to redefine this to destroy the argument in \@FMFPdynamictitlevspace
\long\def\FMFPcompsocitemizethanks#1{\relax\@FMFPbreakcompsocthanksitemfalse\footnotemark
    \protected@xdef\@thanks{\@thanks
        \protect\footnotetext[\the\c@footnote]{\itshape\protect\@FMFPtriggeroneshotfootnoterule
        {\let\FMFPiedlistdecl\relax\protect\begin{FMFPitemize}[\protect\@FMFPsetupcompsocitemizelist]\ignorespaces#1\relax
        \protect\end{FMFPitemize}}\protect\vspace{-1\baselineskip}}}}
\DeclareRobustCommand*{\FMFPcompsocthanksitem}{\item}
\else
% non-compsoc, allow for dual compilation via rerouting to normal \thanks
\long\def\FMFPcompsocitemizethanks#1{\thanks{#1}}
% redirect to "pseudo-par" \hfil\break\indent after swallowing [] from \FMFPcompsocthanksitem[]
\DeclareRobustCommand{\FMFPcompsocthanksitem}{\@ifnextchar [{\@FMFPthanksswallowoptionalarg}%
{\@FMFPthanksswallowoptionalarg[\relax]}}
% be sure and break only after first item, be sure and ignore spaces after optional argument
\def\@FMFPthanksswallowoptionalarg[#1]{\relax\if@FMFPbreakcompsocthanksitem\hfil\break
\indent\fi\@FMFPbreakcompsocthanksitemtrue\ignorespaces}
\fi


% V1.6b define the \FMFPpeerreviewmaketitle as needed
\ifCLASSOPTIONpeerreview
\def\FMFPpeerreviewmaketitle{\@FMFPcleardoublepage{empty}%
\ifCLASSOPTIONtwocolumn
\twocolumn[\@FMFPpeerreviewmaketitle\@FMFPdynamictitlevspace]
\else
\newpage\@FMFPpeerreviewmaketitle\@FMFPstatictitlevskip
\fi
\thispagestyle{FMFPtitlepagestyle}}
\else
% \FMFPpeerreviewmaketitle does nothing if peer review option has not been selected
\def\FMFPpeerreviewmaketitle{\relax}
\fi

% peerreview formats the repeated title like the title in journal papers.
\def\@FMFPpeerreviewmaketitle{\begin{center}\@FMFPcompsoconly{\sffamily}%
\normalfont\normalsize\vskip0.2em{\Huge\@title\par}\vskip1.0em\par
\end{center}}



% V1.6 
% this is a static rubber spacer between the title/authors and the main text
% used for single column text, or when the title appears in the first column
% of two column text (technotes). 
\def\@FMFPstatictitlevskip{{\normalfont\normalsize
% adjust spacing to next text
% v1.6b handle peer review papers
\ifCLASSOPTIONpeerreview
% for peer review papers, the same value is used for both title pages
% regardless of the other paper modes
   \vskip 1\baselineskip plus 0.375\baselineskip minus 0.1875\baselineskip
\else
   \ifCLASSOPTIONconference% conference
      \vskip 0.6\baselineskip
   \else%
      \ifCLASSOPTIONtechnote% technote
         \vskip 1\baselineskip plus 0.375\baselineskip minus 0.1875\baselineskip%
      \else% journal uses more space
         \vskip 2.5\baselineskip plus 0.75\baselineskip minus 0.375\baselineskip%
      \fi
   \fi
\fi}}


% V1.6
% This is a dynamically determined rigid spacer between the title/authors 
% and the main text. This is used only for single column titles over two 
% column text (most common)
% This is bit tricky because we have to ensure that the textheight of the
% main text is an integer multiple of \baselineskip
% otherwise underfull vbox problems may develop in the second column of the
% text on the titlepage
% The possible use of \FMFPpubid must also be taken into account.
\def\@FMFPdynamictitlevspace{{%
    % we run within a group so that all the macros can be forgotten when we are done
    \long\def\thanks##1{\relax}%don't allow \thanks to run when we evaluate the vbox height
    \long\def\FMFPcompsocitemizethanks##1{\relax}%don't allow \FMFPcompsocitemizethanks to run when we evaluate the vbox height
    \normalfont\normalsize% we declare more descriptive variable names
    \let\@FMFPmaintextheight=\@FMFPtrantmpdimenA%height of the main text columns
    \let\@FMFPINTmaintextheight=\@FMFPtrantmpdimenB%height of the main text columns with integer # lines
    % set the nominal and minimum values for the title spacer
    % the dynamic algorithm will not allow the spacer size to
    % become less than \@FMFPMINtitlevspace - instead it will be
    % lengthened
    % default to journal values
    \def\@FMFPNORMtitlevspace{12.5\baselineskip}%
    \def\@FMFPMINtitlevspace{12\baselineskip}%
    % conferences and technotes need tighter spacing
    \ifCLASSOPTIONconference%conference
     \def\@FMFPNORMtitlevspace{11\baselineskip}%
     \def\@FMFPMINtitlevspace{10.75\baselineskip}%
    \fi
    \ifCLASSOPTIONtechnote%technote
      \def\@FMFPNORMtitlevspace{11\baselineskip}%
      \def\@FMFPMINtitlevspace{10.75\baselineskip}%
    \fi%
    % get the height that the title will take up
    \ifCLASSOPTIONpeerreview
       \settoheight{\@FMFPmaintextheight}{\vbox{\hsize\textwidth \@FMFPpeerreviewmaketitle}}%
    \else
       \settoheight{\@FMFPmaintextheight}{\vbox{\hsize\textwidth \@maketitle}}%
    \fi
    \@FMFPmaintextheight=-\@FMFPmaintextheight% title takes away from maintext, so reverse sign
    % add the height of the page textheight
    \advance\@FMFPmaintextheight by \textheight%
    % correct for title pages using pubid
    \ifCLASSOPTIONpeerreview\else
       % peerreview papers use the pubid on the cover page only.
       % And the cover page uses a static spacer.
       \if@FMFPusingpubid\advance\@FMFPmaintextheight by -\@FMFPpubidpullup\fi
    \fi%
    % subtract off the nominal value of the title bottom spacer
    \advance\@FMFPmaintextheight by -\@FMFPNORMtitlevspace%
    % \topskip takes away some too
    \advance\@FMFPmaintextheight by -\topskip%
    % calculate the column height of the main text for lines
    % now we calculate the main text height as if holding
    % an integer number of \normalsize lines after the first
    % and discard any excess fractional remainder
    % we subtracted the first line, because the first line
    % is placed \topskip into the maintext, not \baselineskip like the
    % rest of the lines.
    \@FMFPINTmaintextheight=\@FMFPmaintextheight%
    \divide\@FMFPINTmaintextheight  by \baselineskip%
    \multiply\@FMFPINTmaintextheight  by \baselineskip%
    % now we calculate how much the title spacer height will
    % have to be reduced from nominal (\@FMFPREDUCEmaintextheight is always
    % a positive value) so that the maintext area will contain an integer
    % number of normal size lines
    % we change variable names here (to avoid confusion) as we no longer
    % need \@FMFPINTmaintextheight and can reuse its dimen register
    \let\@FMFPREDUCEmaintextheight=\@FMFPINTmaintextheight%
    \advance\@FMFPREDUCEmaintextheight by -\@FMFPmaintextheight%
    \advance\@FMFPREDUCEmaintextheight by \baselineskip%
    % this is the calculated height of the spacer
    % we change variable names here (to avoid confusion) as we no longer
    % need \@FMFPmaintextheight and can reuse its dimen register
    \let\@FMFPCOMPENSATElen=\@FMFPmaintextheight%
    \@FMFPCOMPENSATElen=\@FMFPNORMtitlevspace% set the nominal value
    % we go with the reduced length if it is smaller than an increase
    \ifdim\@FMFPREDUCEmaintextheight < 0.5\baselineskip\relax%
     \advance\@FMFPCOMPENSATElen by -\@FMFPREDUCEmaintextheight%
     % if the resulting spacer is too small back out and go with an increase instead
     \ifdim\@FMFPCOMPENSATElen<\@FMFPMINtitlevspace\relax%
      \advance\@FMFPCOMPENSATElen by \baselineskip%
     \fi%
    \else%
     % go with an increase because it is closer to the nominal than a decrease
     \advance\@FMFPCOMPENSATElen by -\@FMFPREDUCEmaintextheight%
     \advance\@FMFPCOMPENSATElen by \baselineskip%
    \fi%
    % set the calculated rigid spacer
    \vspace{\@FMFPCOMPENSATElen}}}



% V1.6
% we allow the user access to the last part of the title area
% useful in emergencies such as when a different spacing is needed
% This text is NOT compensated for in the dynamic sizer.
\let\@FMFPaftertitletext=\relax
\long\def\FMFPaftertitletext#1{\def\@FMFPaftertitletext{#1}}

% V1.7 provide a way for users to enter abstract and keywords
% into the onecolumn title are. This text is compensated for
% in the dynamic sizer.
\let\@FMFPcompsoctitleabstractindextext=\relax
\long\def\FMFPcompsoctitleabstractindextext#1{\def\@FMFPcompsoctitleabstractindextext{#1}}
% V1.7 provide a way for users to get the \@FMFPcompsoctitleabstractindextext if
% not in compsoc journal mode - this way abstract and keywords can be placed
% in their conventional position if not in compsoc mode.
\def\FMFPdisplaynotcompsoctitleabstractindextext{%
\ifCLASSOPTIONcompsoc% display if compsoc conf
\ifCLASSOPTIONconference\@FMFPcompsoctitleabstractindextext\fi
\else% or if not compsoc
\@FMFPcompsoctitleabstractindextext\fi}


% command to allow alteration of baselinestretch, but only if the current
% baselineskip is unity. Used to tweak the compsoc abstract and keywords line spacing.
\def\@FMFPtweakunitybaselinestretch#1{{\def\baselinestretch{1}\selectfont
\global\@tempskipa\baselineskip}\ifnum\@tempskipa=\baselineskip%
\def\baselinestretch{#1}\selectfont\fi\relax}


% abstract and keywords are in \small, except 
% for 9pt docs in which they are in \footnotesize
% Because 9pt docs use an 8pt footnotesize, \small
% becomes a rather awkward 8.5pt
\def\@FMFPabskeysecsize{\small}
\ifx\CLASSOPTIONpt\@FMFPptsizenine
 \def\@FMFPabskeysecsize{\footnotesize}
\fi

% compsoc journals use \footnotesize, compsoc conferences use normalsize
\@FMFPcompsoconly{\def\@FMFPabskeysecsize{\footnotesize}}
\@FMFPcompsocconfonly{\def\@FMFPabskeysecsize{\normalsize}}




% V1.6 have abstract and keywords strip leading spaces, pars and newlines
% so that spacing is more tightly controlled.
\def\abstract{\normalfont
    \if@twocolumn
      \par\@FMFPabskeysecsize\bfseries\leavevmode\kern-1pt\textit{\abstractname} -- \relax
    \else
      \begin{center}\vspace{-1.78ex}\@FMFPabskeysecsize\textbf{\abstractname}\end{center}\quotation\@FMFPabskeysecsize
    \fi\@FMFPgobbleleadPARNLSP}
% V1.6 FMFP wants only 1 pica from end of abstract to introduction heading when in 
% conference mode (the heading already has this much above it)
\def\endabstract{\relax\ifCLASSOPTIONconference\vspace{0ex}\else\vspace{1.34ex}\fi\par\if@twocolumn\else\endquotation\fi
    \normalfont\normalsize}

\def\FMFPkeywords{\normalfont
    \if@twocolumn
      \@FMFPabskeysecsize\leavevmode\kern-1pt\textit{\FMFPkeywordsname} -- \relax
    \else
      \begin{center}\@FMFPabskeysecsize\textit{\FMFPkeywordsname}\end{center}\quotation\@FMFPabskeysecsize
    \fi\itshape\@FMFPgobbleleadPARNLSP}
\def\endFMFPkeywords{\relax\ifCLASSOPTIONtechnote\vspace{1.34ex}\else\vspace{0.5ex}\fi
    \par\if@twocolumn\else\endquotation\fi%
    \normalfont\normalsize}

% V1.7 compsoc keywords index terms
\ifCLASSOPTIONcompsoc
  \ifCLASSOPTIONconference% compsoc conference
\def\abstract{\normalfont
      \begin{center}\@FMFPabskeysecsize\textbf{\large\abstractname}\end{center}\vskip 0.5\baselineskip plus 0.1\baselineskip minus 0.1\baselineskip
      \if@twocolumn\else\quotation\fi\itshape\@FMFPabskeysecsize%
      \par\@FMFPgobbleleadPARNLSP}
\def\FMFPkeywords{\normalfont\vskip 1.5\baselineskip plus 0.25\baselineskip minus 0.25\baselineskip
      \begin{center}\@FMFPabskeysecsize\textbf{\large\FMFPkeywordsname}\end{center}\vskip 0.5\baselineskip plus 0.1\baselineskip minus 0.1\baselineskip
      \if@twocolumn\else\quotation\fi\itshape\@FMFPabskeysecsize%
      \par\@FMFPgobbleleadPARNLSP}
  \else% compsoc not conference
\def\abstract{\normalfont\@FMFPtweakunitybaselinestretch{1.15}\sffamily
    \if@twocolumn
      \@FMFPabskeysecsize\noindent\textbf{\abstractname}---\relax
    \else
      \begin{center}\vspace{-1.78ex}\@FMFPabskeysecsize\textbf{\abstractname}\end{center}\quotation\@FMFPabskeysecsize%
    \fi\@FMFPgobbleleadPARNLSP}
\def\FMFPkeywords{\normalfont\@FMFPtweakunitybaselinestretch{1.15}\sffamily
    \if@twocolumn
      \@FMFPabskeysecsize\vskip 0.5\baselineskip plus 0.25\baselineskip minus 0.25\baselineskip\noindent
      \textbf{\FMFPkeywordsname}---\relax
    \else
      \begin{center}\@FMFPabskeysecsize\textit{\FMFPkeywordsname}\end{center}\quotation\@FMFPabskeysecsize%
    \fi\@FMFPgobbleleadPARNLSP}
  \fi
\fi



% gobbles all leading \, \\ and \par, upon finding first token that
% is not a \ , \\ or a \par, it ceases and returns that token
% 
% used to strip leading \, \\ and \par from the input
% so that such things in the beginning of an environment will not
% affect the formatting of the text
\long\def\@FMFPgobbleleadPARNLSP#1{\let\@FMFPswallowthistoken=0%
\let\@FMFPgobbleleadPARNLSPtoken#1%
\let\@FMFPgobbleleadPARtoken=\par%
\let\@FMFPgobbleleadNLtoken=\\%
\let\@FMFPgobbleleadSPtoken=\ %
\def\@FMFPgobbleleadSPMACRO{\ }%
\ifx\@FMFPgobbleleadPARNLSPtoken\@FMFPgobbleleadPARtoken%
\let\@FMFPswallowthistoken=1%
\fi%
\ifx\@FMFPgobbleleadPARNLSPtoken\@FMFPgobbleleadNLtoken%
\let\@FMFPswallowthistoken=1%
\fi%
\ifx\@FMFPgobbleleadPARNLSPtoken\@FMFPgobbleleadSPtoken%
\let\@FMFPswallowthistoken=1%
\fi%
% a control space will come in as a macro
% when it is the last one on a line
\ifx\@FMFPgobbleleadPARNLSPtoken\@FMFPgobbleleadSPMACRO%
\let\@FMFPswallowthistoken=1%
\fi%
% if we have to swallow this token, do so and taste the next one
% else spit it out and stop gobbling
\ifx\@FMFPswallowthistoken 1\let\@FMFPnextgobbleleadPARNLSP=\@FMFPgobbleleadPARNLSP\else%
\let\@FMFPnextgobbleleadPARNLSP=#1\fi%
\@FMFPnextgobbleleadPARNLSP}%




% TITLING OF SECTIONS
\def\@FMFPsectpunct{:\ \,}  % Punctuation after run-in section heading  (headings which are
                            % part of the paragraphs), need little bit more than a single space
                            % spacing from section number to title
% compsoc conferences use regular period/space punctuation
\ifCLASSOPTIONcompsoc
\ifCLASSOPTIONconference
\def\@FMFPsectpunct{.\ }
\fi\fi

\def\@seccntformat#1{\hb@xt@ 1.4em{\csname the#1dis\endcsname\hss\relax}}
\def\@seccntformatinl#1{\hb@xt@ 1.1em{\csname the#1dis\endcsname\hss\relax}}
\def\@seccntformatch#1{\csname the#1dis\endcsname\hskip 1em\relax}

\ifCLASSOPTIONcompsoc
% compsoc journals need extra spacing
\ifCLASSOPTIONconference\else
\def\@seccntformat#1{\csname the#1dis\endcsname\hskip 1em\relax}
\fi\fi

%v1.7 put {} after #6 to allow for some types of user font control
%and use \@@par rather than \par
\def\@sect#1#2#3#4#5#6[#7]#8{%
  \ifnum #2>\c@secnumdepth
     \let\@svsec\@empty
  \else
     \refstepcounter{#1}%
     % load section label and spacer into \@svsec
     \ifnum #2=1
        \protected@edef\@svsec{\@seccntformatch{#1}\relax}%
     \else
        \ifnum #2>2
           \protected@edef\@svsec{\@seccntformatinl{#1}\relax}%
        \else
           \protected@edef\@svsec{\@seccntformat{#1}\relax}%
        \fi
     \fi
  \fi%
  \@tempskipa #5\relax
  \ifdim \@tempskipa>\z@% tempskipa determines whether is treated as a high
     \begingroup #6{\relax% or low level heading
      \noindent % subsections are NOT indented
       % print top level headings. \@svsec is label, #8 is heading title
       % FMFP does not block indent the section title text, it flows like normal
       {\hskip #3\relax\@svsec}{\interlinepenalty \@M #8\@@par}}%
     \endgroup
     \addcontentsline{toc}{#1}{\ifnum #2>\c@secnumdepth\relax\else
               \protect\numberline{\csname the#1\endcsname}\fi#7}%
  \else % printout low level headings
     % svsechd seems to swallow the trailing space, protect it with \mbox{}
     % got rid of sectionmark stuff
     \def\@svsechd{#6{\hskip #3\relax\@svsec #8\@FMFPsectpunct\mbox{}}%
     \addcontentsline{toc}{#1}{\ifnum #2>\c@secnumdepth\relax\else
               \protect\numberline{\csname the#1\endcsname}\fi#7}}%
  \fi%skip down
  \@xsect{#5}}


% section* handler
%v1.7 put {} after #4 to allow for some types of user font control
%and use \@@par rather than \par
\def\@ssect#1#2#3#4#5{\@tempskipa #3\relax
  \ifdim \@tempskipa>\z@
     %\begingroup #4\@hangfrom{\hskip #1}{\interlinepenalty \@M #5\par}\endgroup
     % FMFP does not block indent the section title text, it flows like normal
     \begingroup \noindent #4{\relax{\hskip #1}{\interlinepenalty \@M #5\@@par}}\endgroup
  % svsechd swallows the trailing space, protect it with \mbox{}
  \else \def\@svsechd{#4{\hskip #1\relax #5\@FMFPsectpunct\mbox{}}}\fi
  \@xsect{#3}}


%% SECTION heading spacing and font
%%
% arguments are: #1 - sectiontype name
% (for \@sect)   #2 - section level
%                #3 - section heading indent
%                #4 - top separation (absolute value used, neg indicates not to indent main text)
%                     If negative, make stretch parts negative too!
%                #5 - (absolute value used) positive: bottom separation after heading,
%                      negative: amount to indent main text after heading
%                Both #4 and #5 negative means to indent main text and use negative top separation
%                #6 - font control
% You've got to have \normalfont\normalsize in the font specs below to prevent
% trouble when you do something like:
% \section{Note}{\ttfamily TT-TEXT} is known to ... 
% FMFP sometimes REALLY stretches the area before a section
% heading by up to about 0.5in. However, it may not be a good
% idea to let LaTeX have quite this much rubber.
\ifCLASSOPTIONconference%
% FMFP wants section heading spacing to decrease for conference mode
\def\section{\@startsection{section}{1}{\z@}{1.5ex plus 1.5ex minus 0.5ex}%
{1sp}{\normalfont\normalsize\scshape}}%
\def\subsection{\@startsection{subsection}{2}{\z@}{1.5ex plus 1.5ex minus 0.5ex}%
{1sp}{\normalfont\normalsize}}% %Changed by AJC on 28.04.2018
\else % for journals
\def\section{\@startsection{section}{1}{\z@}{3.0ex plus 1.5ex minus 1.5ex}% V1.6 3.0ex from 3.5ex
{0.7ex plus 1ex minus 0ex}{\normalfont\normalsize\centering\scshape}}%
\def\subsection{\@startsection{subsection}{2}{\z@}{3.5ex plus 1.5ex minus 1.5ex}%
{0.7ex plus .5ex minus 0ex}{\normalfont\normalsize\itshape}}%
\fi

% for both journals and conferences
% decided to put in a little rubber above the section, might help somebody
\def\subsubsection{\@startsection{subsubsection}{3}{\parindent}{0ex plus 0.1ex minus 0.1ex}%
{0ex}{\normalfont\normalsize\itshape}}%
\def\paragraph{\@startsection{paragraph}{4}{2\parindent}{0ex plus 0.1ex minus 0.1ex}%
{0ex}{\normalfont\normalsize\itshape}}%


% compsoc
\ifCLASSOPTIONcompsoc
\ifCLASSOPTIONconference
% compsoc conference
\def\section{\@startsection{section}{1}{\z@}{1\baselineskip plus 0.25\baselineskip minus 0.25\baselineskip}%
{1\baselineskip plus 0.25\baselineskip minus 0.25\baselineskip}{\normalfont\large\bfseries}}%
\def\subsection{\@startsection{subsection}{2}{\z@}{1\baselineskip plus 0.25\baselineskip minus 0.25\baselineskip}%
{1\baselineskip plus 0.25\baselineskip minus 0.25\baselineskip}{\normalfont\sublargesize\bfseries}}%
\def\subsubsection{\@startsection{subsubsection}{3}{\z@}{1\baselineskip plus 0.25\baselineskip minus 0.25\baselineskip}%
{0ex}{\normalfont\normalsize\bfseries}}%
\def\paragraph{\@startsection{paragraph}{4}{2\parindent}{0ex plus 0.1ex minus 0.1ex}%
{0ex}{\normalfont\normalsize}}%
\else% compsoc journals
% use negative top separation as compsoc journals do not indent paragraphs after section titles
\def\section{\@startsection{section}{1}{\z@}{-3ex plus -2ex minus -1.5ex}%
{0.7ex plus 1ex minus 0ex}{\normalfont\large\sffamily\bfseries\scshape}}%
% Note that subsection and smaller may not be correct for the Computer Society,
% I have to look up an example.
\def\subsection{\@startsection{subsection}{2}{\z@}{-3.5ex plus -1.5ex minus -1.5ex}%
{0.7ex plus .5ex minus 0ex}{\normalfont\normalsize\sffamily\bfseries}}%
\def\subsubsection{\@startsection{subsubsection}{3}{\z@}{-2.5ex plus -1ex minus -1ex}%
{0.5ex plus 0.5ex minus 0ex}{\normalfont\normalsize\sffamily\itshape}}%
\def\paragraph{\@startsection{paragraph}{4}{2\parindent}{-0ex plus -0.1ex minus -0.1ex}%
{0ex}{\normalfont\normalsize}}%
\fi\fi




%% ENVIRONMENTS
% "box" symbols at end of proofs
\def\FMFPQEDclosed{\mbox{\rule[0pt]{1.3ex}{1.3ex}}} % for a filled box
% V1.6 some journals use an open box instead that will just fit around a closed one
\def\FMFPQEDopen{{\setlength{\fboxsep}{0pt}\setlength{\fboxrule}{0.2pt}\fbox{\rule[0pt]{0pt}{1.3ex}\rule[0pt]{1.3ex}{0pt}}}}
\ifCLASSOPTIONcompsoc
\def\FMFPQED{\FMFPQEDopen}   % default to open for compsoc
\else
\def\FMFPQED{\FMFPQEDclosed} % otherwise default to closed
\fi

% v1.7 name change to avoid namespace collision with amsthm. Also add support
% for an optional argument.
\def\FMFPproof{\@ifnextchar[{\@FMFPproof}{\@FMFPproof[\FMFPproofname]}}
\def\@FMFPproof[#1]{\par\noindent\hspace{2em}{\itshape #1: }}
\def\endFMFPproof{\hspace*{\fill}~\FMFPQED\par}


%\itemindent is set to \z@ by list, so define new temporary variable
\newdimen\@FMFPtmpitemindent
\def\@begintheorem#1#2{\@FMFPtmpitemindent\itemindent\topsep 0pt\rmfamily\trivlist%
    \item[\hskip \labelsep{\indent\itshape #1\ #2:}]\itemindent\@FMFPtmpitemindent}
\def\@opargbegintheorem#1#2#3{\@FMFPtmpitemindent\itemindent\topsep 0pt\rmfamily \trivlist%
% V1.6 FMFP is back to using () around theorem names which are also in italics
% Thanks to Christian Peel for reporting this.
    \item[\hskip\labelsep{\indent\itshape #1\ #2\ (#3):}]\itemindent\@FMFPtmpitemindent}
% V1.7 remove bogus \unskip that caused equations in theorems to collide with
% lines below.
\def\@endtheorem{\endtrivlist}

% V1.6
% display command for the section the theorem is in - so that \thesection
% is not used as this will be in Roman numerals when we want arabic.
% LaTeX2e uses \def\@thmcounter#1{\noexpand\arabic{#1}} for the theorem number
% (second part) display and \def\@thmcountersep{.} as a separator.
% V1.7 intercept calls to the section counter and reroute to \@FMFPthmcounterinsection
% to allow \appendix(ices} to override as needed.
%
% special handler for sections, allows appendix(ices) to override
\gdef\@FMFPthmcounterinsection#1{\arabic{#1}}
% string macro
\edef\@FMFPstringsection{section}

% redefine the #1#2[#3] form of newtheorem to use a hook to \@FMFPthmcounterinsection
% if section in_counter is used
\def\@xnthm#1#2[#3]{%
  \expandafter\@ifdefinable\csname #1\endcsname
    {\@definecounter{#1}\@newctr{#1}[#3]%
     \edef\@FMFPstringtmp{#3}
     \ifx\@FMFPstringtmp\@FMFPstringsection
     \expandafter\xdef\csname the#1\endcsname{%
     \noexpand\@FMFPthmcounterinsection{#3}\@thmcountersep
          \@thmcounter{#1}}%
     \else
     \expandafter\xdef\csname the#1\endcsname{%
       \expandafter\noexpand\csname the#3\endcsname \@thmcountersep
          \@thmcounter{#1}}%
     \fi
     \global\@namedef{#1}{\@thm{#1}{#2}}%
     \global\@namedef{end#1}{\@endtheorem}}}



%% SET UP THE DEFAULT PAGESTYLE
\ps@headings
\pagenumbering{arabic}

% normally the page counter starts at 1
\setcounter{page}{1}
% however, for peerreview the cover sheet is page 0 or page -1
% (for duplex printing)
\ifCLASSOPTIONpeerreview
   \if@twoside
      \setcounter{page}{-1}
   \else
      \setcounter{page}{0}
   \fi
\fi

% standard book class behavior - let bottom line float up and down as
% needed when single sided
\ifCLASSOPTIONtwoside\else\raggedbottom\fi
% if two column - turn on twocolumn, allow word spacings to stretch more and
% enforce a rigid position for the last lines
\ifCLASSOPTIONtwocolumn
% the peer review option delays invoking twocolumn
   \ifCLASSOPTIONpeerreview\else
      \twocolumn
   \fi
\sloppy 
\flushbottom
\fi




% \APPENDIX and \APPENDICES definitions

% This is the \@ifmtarg command from the LaTeX ifmtarg package
% by Peter Wilson (CUA) and Donald Arseneau
% \@ifmtarg is used to determine if an argument to a command
% is present or not.
% For instance:
% \@ifmtarg{#1}{\typeout{empty}}{\typeout{has something}}
% \@ifmtarg is used with our redefined \section command if
% \appendices is invoked.
% The command \section will behave slightly differently depending
% on whether the user specifies a title: 
% \section{My appendix title}
% or not:
% \section{}
% This way, we can eliminate the blank lines where the title
% would be, and the unneeded : after Appendix in the table of
% contents 
\begingroup
\catcode`\Q=3
\long\gdef\@ifmtarg#1{\@xifmtarg#1QQ\@secondoftwo\@firstoftwo\@nil}
\long\gdef\@xifmtarg#1#2Q#3#4#5\@nil{#4}
\endgroup
% end of \@ifmtarg defs


% V1.7
% command that allows the one time saving of the original definition
% of section to \@FMFPappendixsavesection for \appendix or \appendices 
% we don't save \section here as it may be redefined later by other
% packages (hyperref.sty, etc.)
\def\@FMFPsaveoriginalsectiononce{\let\@FMFPappendixsavesection\section
\let\@FMFPsaveoriginalsectiononce\relax}

% neat trick to grab and process the argument from \section{argument}
% we process differently if the user invoked \section{} with no
% argument (title)
% note we reroute the call to the old \section*
\def\@FMFPprocessthesectionargument#1{%
\@ifmtarg{#1}{%
\@FMFPappendixsavesection*{\appendixname~\thesectiondis}%
\addcontentsline{toc}{section}{\appendixname~\thesection}}{%
\@FMFPappendixsavesection*{\appendixname~\thesectiondis \\* #1}%
\addcontentsline{toc}{section}{\appendixname~\thesection: #1}}}

% we use this if the user calls \section{} after
% \appendix-- which has no meaning. So, we ignore the
% command and its argument. Then, warn the user.
\def\@FMFPdestroythesectionargument#1{\typeout{** WARNING: Ignoring useless
\protect\section\space in Appendix (line \the\inputlineno).}}


% remember \thesection forms will be displayed in \ref calls
% and in the Table of Contents.
% The \sectiondis form is used in the actual heading itself

% appendix command for one single appendix
% normally has no heading. However, if you want a 
% heading, you can do so via the optional argument:
% \appendix[Optional Heading]
\def\appendix{\relax}
\renewcommand{\appendix}[1][]{\@FMFPsaveoriginalsectiononce\par
    % v1.6 keep hyperref's identifiers unique
    \gdef\theHsection{Appendix.A}%
    % v1.6 adjust hyperref's string name for the section
    \xdef\Hy@chapapp{appendix}%
    \setcounter{section}{0}%
    \setcounter{subsection}{0}%
    \setcounter{subsubsection}{0}%
    \setcounter{paragraph}{0}%
    \gdef\thesection{A}%
    \gdef\thesectiondis{}% 
    \gdef\thesubsection{\Alph{subsection}}%
    \gdef\@FMFPthmcounterinsection##1{A}
    \refstepcounter{section}% update the \ref counter
    \@ifmtarg{#1}{\@FMFPappendixsavesection*{\appendixname}%
                  \addcontentsline{toc}{section}{\appendixname}}{%
             \@FMFPappendixsavesection*{\appendixname~\\* #1}%
             \addcontentsline{toc}{section}{\appendixname: #1}}%
    % redefine \section command for appendix
    % leave \section* as is
    \def\section{\@ifstar{\@FMFPappendixsavesection*}{%
                    \@FMFPdestroythesectionargument}}% throw out the argument
                                                     % of the normal form
}



% appendices command for multiple appendices
% user then calls \section with an argument (possibly empty) to
% declare the individual appendices
\def\appendices{\@FMFPsaveoriginalsectiononce\par
    % v1.6 keep hyperref's identifiers unique
    \gdef\theHsection{Appendix.\Alph{section}}%
    % v1.6 adjust hyperref's string name for the section
    \xdef\Hy@chapapp{appendix}%
    \setcounter{section}{-1}% we want \refstepcounter to use section 0
    \setcounter{subsection}{0}%
    \setcounter{subsubsection}{0}%
    \setcounter{paragraph}{0}%
    \ifCLASSOPTIONromanappendices%
    \gdef\thesection{\Roman{section}}%
    \gdef\thesectiondis{\Roman{section}}%
    \@FMFPcompsocconfonly{\gdef\thesectiondis{\Roman{section}.}}%
    \gdef\@FMFPthmcounterinsection##1{A\arabic{##1}}
    \else%
    \gdef\thesection{\Alph{section}}%
    \gdef\thesectiondis{\Alph{section}}%
    \@FMFPcompsocconfonly{\gdef\thesectiondis{\Alph{section}.}}%
    \gdef\@FMFPthmcounterinsection##1{\Alph{##1}}
    \fi%
    \refstepcounter{section}% update the \ref counter
    \setcounter{section}{0}% NEXT \section will be the FIRST appendix
    % redefine \section command for appendices
    % leave \section* as is
    \def\section{\@ifstar{\@FMFPappendixsavesection*}{% process the *-form
                    \refstepcounter{section}% or is a new section so,
                    \@FMFPprocessthesectionargument}}% process the argument 
                                                 % of the normal form
}



% \FMFPPARstart
% Definition for the big two line drop cap letter at the beginning of the
% first paragraph of journal papers. The first argument is the first letter
% of the first word, the second argument is the remaining letters of the
% first word which will be rendered in upper case.
% In V1.6 this has been completely rewritten to:
% 
% 1. no longer have problems when the user begins an environment
%    within the paragraph that uses \FMFPPARstart.
% 2. auto-detect and use the current font family
% 3. revise handling of the space at the end of the first word so that
%    interword glue will now work as normal.
% 4. produce correctly aligned edges for the (two) indented lines.
% 
% We generalize things via control macros - playing with these is fun too.
% 
% V1.7 added more control macros to make it easy for FMFPtrantools.sty users
% to change the font style.
% 
% the number of lines that are indented to clear it
% may need to increase if using decenders
\def\@FMFPPARstartDROPLINES{2}
% minimum number of lines left on a page to allow a \@FMFPPARstart
% Does not take into consideration rubber shrink, so it tends to
% be overly cautious
\def\@FMFPPARstartMINPAGELINES{2}
% V1.7 the height of the drop cap is adjusted to match the height of this text
% in the current font (when \FMFPPARstart is called).
\def\@FMFPPARstartHEIGHTTEXT{T}
% the depth the letter is lowered below the baseline
% the height (and size) of the letter is determined by the sum
% of this value and the height of the \@FMFPPARstartHEIGHTTEXT in the current
% font. It is a good idea to set this value in terms of the baselineskip
% so that it can respond to changes therein.
\def\@FMFPPARstartDROPDEPTH{1.1\baselineskip}
% V1.7 the font the drop cap will be rendered in,
% can take zero or one argument.
\def\@FMFPPARstartFONTSTYLE{\bfseries}
% V1.7 any additional, non-font related commands needed to modify
% the drop cap letter, can take zero or one argument.
\def\@FMFPPARstartCAPSTYLE{\MakeUppercase}
% V1.7 the font that will be used to render the rest of the word,
% can take zero or one argument.
\def\@FMFPPARstartWORDFONTSTYLE{\relax}
% V1.7 any additional, non-font related commands needed to modify
% the rest of the word, can take zero or one argument.
\def\@FMFPPARstartWORDCAPSTYLE{\MakeUppercase}
% This is the horizontal separation distance from the drop letter to the main text.
% Lengths that depend on the font (e.g., ex, em, etc.) will be referenced
% to the font that is active when \FMFPPARstart is called. 
\def\@FMFPPARstartSEP{0.15em}
% V1.7 horizontal offset applied to the left of the drop cap.
\def\@FMFPPARstartHOFFSET{0em}
% V1.7 Italic correction command applied at the end of the drop cap.
\def\@FMFPPARstartITLCORRECT{\/}

% V1.7 compoc uses nonbold drop cap and small caps word style
\ifCLASSOPTIONcompsoc
\def\@FMFPPARstartFONTSTYLE{\mdseries}
\def\@FMFPPARstartWORDFONTSTYLE{\scshape}
\def\@FMFPPARstartWORDCAPSTYLE{\relax}
\fi

% definition of \FMFPPARstart
% THIS IS A CONTROLLED SPACING AREA, DO NOT ALLOW SPACES WITHIN THESE LINES
% 
% The token \@FMFPPARstartfont will be globally defined after the first use
% of \FMFPPARstart and will be a font command which creates the big letter
% The first argument is the first letter of the first word and the second
% argument is the rest of the first word(s).
\def\FMFPPARstart#1#2{\par{%
% if this page does not have enough space, break it and lets start
% on a new one
\@FMFPtranneedspace{\@FMFPPARstartMINPAGELINES\baselineskip}{\relax}%
% V1.7 move this up here in case user uses \textbf for \@FMFPPARstartFONTSTYLE
% which uses command \leavevmode which causes an unwanted \indent to be issued
\noindent
% calculate the desired height of the big letter
% it extends from the top of \@FMFPPARstartHEIGHTTEXT in the current font
% down to \@FMFPPARstartDROPDEPTH below the current baseline
\settoheight{\@FMFPtrantmpdimenA}{\@FMFPPARstartHEIGHTTEXT}%
\addtolength{\@FMFPtrantmpdimenA}{\@FMFPPARstartDROPDEPTH}%
% extract the name of the current font in bold
% and place it in \@FMFPPARstartFONTNAME
\def\@FMFPPARstartGETFIRSTWORD##1 ##2\relax{##1}%
{\@FMFPPARstartFONTSTYLE{\selectfont\edef\@FMFPPARstartFONTNAMESPACE{\fontname\font\space}%
\xdef\@FMFPPARstartFONTNAME{\expandafter\@FMFPPARstartGETFIRSTWORD\@FMFPPARstartFONTNAMESPACE\relax}}}%
% define a font based on this name with a point size equal to the desired
% height of the drop letter
\font\@FMFPPARstartsubfont\@FMFPPARstartFONTNAME\space at \@FMFPtrantmpdimenA\relax%
% save this value as a counter (integer) value (sp points)
\@FMFPtrantmpcountA=\@FMFPtrantmpdimenA%
% now get the height of the actual letter produced by this font size
\settoheight{\@FMFPtrantmpdimenB}{\@FMFPPARstartsubfont\@FMFPPARstartCAPSTYLE{#1}}%
% If something bogus happens like the first argument is empty or the
% current font is strange, do not allow a zero height.
\ifdim\@FMFPtrantmpdimenB=0pt\relax%
\typeout{** WARNING: FMFPPARstart drop letter has zero height! (line \the\inputlineno)}%
\typeout{ Forcing the drop letter font size to 10pt.}%
\@FMFPtrantmpdimenB=10pt%
\fi%
% and store it as a counter
\@FMFPtrantmpcountB=\@FMFPtrantmpdimenB%
% Since a font size doesn't exactly correspond to the height of the capital
% letters in that font, the actual height of the letter, \@FMFPtrantmpcountB,
% will be less than that desired, \@FMFPtrantmpcountA
% we need to raise the font size, \@FMFPtrantmpdimenA 
% by \@FMFPtrantmpcountA / \@FMFPtrantmpcountB
% But, TeX doesn't have floating point division, so we have to use integer
% division. Hence the use of the counters.
% We need to reduce the denominator so that the loss of the remainder will
% have minimal affect on the accuracy of the result
\divide\@FMFPtrantmpcountB by 200%
\divide\@FMFPtrantmpcountA by \@FMFPtrantmpcountB%
% Then reequalize things when we use TeX's ability to multiply by
% floating point values
\@FMFPtrantmpdimenB=0.005\@FMFPtrantmpdimenA%
\multiply\@FMFPtrantmpdimenB by \@FMFPtrantmpcountA%
% \@FMFPPARstartfont is globaly set to the calculated font of the big letter
% We need to carry this out of the local calculation area to to create the
% big letter.
\global\font\@FMFPPARstartfont\@FMFPPARstartFONTNAME\space at \@FMFPtrantmpdimenB%
% Now set \@FMFPtrantmpdimenA to the width of the big letter
% We need to carry this out of the local calculation area to set the
% hanging indent
\settowidth{\global\@FMFPtrantmpdimenA}{\@FMFPPARstartfont
\@FMFPPARstartCAPSTYLE{#1\@FMFPPARstartITLCORRECT}}}%
% end of the isolated calculation environment
% add in the extra clearance we want
\advance\@FMFPtrantmpdimenA by \@FMFPPARstartSEP\relax%
% add in the optional offset
\advance\@FMFPtrantmpdimenA by \@FMFPPARstartHOFFSET\relax%
% V1.7 don't allow negative offsets to produce negative hanging indents
\@FMFPtrantmpdimenB\@FMFPtrantmpdimenA
\ifnum\@FMFPtrantmpdimenB < 0 \@FMFPtrantmpdimenB 0pt\fi
% \@FMFPtrantmpdimenA has the width of the big letter plus the
% separation space and \@FMFPPARstartfont is the font we need to use
% Now, we make the letter and issue the hanging indent command
% The letter is placed in a box of zero width and height so that other
% text won't be displaced by it.
\hangindent\@FMFPtrantmpdimenB\hangafter=-\@FMFPPARstartDROPLINES%
\makebox[0pt][l]{\hspace{-\@FMFPtrantmpdimenA}%
\raisebox{-\@FMFPPARstartDROPDEPTH}[0pt][0pt]{\hspace{\@FMFPPARstartHOFFSET}%
\@FMFPPARstartfont\@FMFPPARstartCAPSTYLE{#1\@FMFPPARstartITLCORRECT}%
\hspace{\@FMFPPARstartSEP}}}%
{\@FMFPPARstartWORDFONTSTYLE{\@FMFPPARstartWORDCAPSTYLE{\selectfont#2}}}}






% determines if the space remaining on a given page is equal to or greater
% than the specified space of argument one
% if not, execute argument two (only if the remaining space is greater than zero)
% and issue a \newpage
% 
% example: \@FMFPtranneedspace{2in}{\vfill}
% 
% Does not take into consideration rubber shrinkage, so it tends to
% be overly cautious
% Based on an example posted by Donald Arseneau
% Note this macro uses \@FMFPtrantmpdimenB internally for calculations,
% so DO NOT PASS \@FMFPtrantmpdimenB to this routine
% if you need a dimen register, import with \@FMFPtrantmpdimenA instead
\def\@FMFPtranneedspace#1#2{\penalty-100\begingroup%shield temp variable
\@FMFPtrantmpdimenB\pagegoal\advance\@FMFPtrantmpdimenB-\pagetotal% space left
\ifdim #1>\@FMFPtrantmpdimenB\relax% not enough space left
\ifdim\@FMFPtrantmpdimenB>\z@\relax #2\fi%
\newpage%
\fi\endgroup}



% FMFPbiography ENVIRONMENT
% Allows user to enter biography leaving place for picture (adapts to font size)
% As of V1.5, a new optional argument allows you to have a real graphic!
% V1.5 and later also fixes the "colliding biographies" which could happen when a 
% biography's text was shorter than the space for the photo.
% MDS 7/2001
% V1.6 prevent multiple biographies from making multiple TOC entries
\newif\if@FMFPbiographyTOCentrynotmade
\global\@FMFPbiographyTOCentrynotmadetrue

% biography counter so hyperref can jump directly to the biographies
% and not just the previous section
\newcounter{FMFPbiography}
\setcounter{FMFPbiography}{0}

% photo area size
\def\@FMFPBIOphotowidth{1.0in}    % width of the biography photo area
\def\@FMFPBIOphotodepth{1.25in}   % depth (height) of the biography photo area
% area cleared for photo
\def\@FMFPBIOhangwidth{1.14in}    % width cleared for the biography photo area
\def\@FMFPBIOhangdepth{1.25in}    % depth cleared for the biography photo area
                                  % actual depth will be a multiple of 
                                  % \baselineskip, rounded up
\def\@FMFPBIOskipN{4\baselineskip}% nominal value of the vskip above the biography

\newenvironment{FMFPbiography}[2][]{\normalfont\@FMFPcompsoconly{\sffamily}\footnotesize%
\unitlength 1in\parskip=0pt\par\parindent 1em\interlinepenalty500%
% we need enough space to support the hanging indent
% the nominal value of the spacer
% and one extra line for good measure
\@FMFPtrantmpdimenA=\@FMFPBIOhangdepth%
\advance\@FMFPtrantmpdimenA by \@FMFPBIOskipN%
\advance\@FMFPtrantmpdimenA by 1\baselineskip%
% if this page does not have enough space, break it and lets start
% with a new one
\@FMFPtranneedspace{\@FMFPtrantmpdimenA}{\relax}%
% nominal spacer can strech, not shrink use 1fil so user can out stretch with \vfill
\vskip \@FMFPBIOskipN plus 1fil minus 0\baselineskip%
% the default box for where the photo goes
\def\@FMFPtempbiographybox{{\setlength{\fboxsep}{0pt}\framebox{%
\begin{minipage}[b][\@FMFPBIOphotodepth][c]{\@FMFPBIOphotowidth}\centering PLACE\\ PHOTO\\ HERE \end{minipage}}}}%
%
% detect if the optional argument was supplied, this requires the
% \@ifmtarg command as defined in the appendix section above
% and if so, override the default box with what they want
\@ifmtarg{#1}{\relax}{\def\@FMFPtempbiographybox{\mbox{\begin{minipage}[b][\@FMFPBIOphotodepth][c]{\@FMFPBIOphotowidth}%
\centering%
#1%
\end{minipage}}}}% end if optional argument supplied
% Make an entry into the table of contents only if we have not done so before
\if@FMFPbiographyTOCentrynotmade%
% link labels to the biography counter so hyperref will jump
% to the biography, not the previous section
\setcounter{FMFPbiography}{-1}%
\refstepcounter{FMFPbiography}%
\addcontentsline{toc}{section}{Biographies}%
\global\@FMFPbiographyTOCentrynotmadefalse%
\fi%
% one more biography
\refstepcounter{FMFPbiography}%
% Make an entry for this name into the table of contents 
\addcontentsline{toc}{subsection}{#2}%
% V1.6 properly handle if a new paragraph should occur while the
% hanging indent is still active. Do this by redefining \par so
% that it will not start a new paragraph. (But it will appear to the
% user as if it did.) Also, strip any leading pars, newlines, or spaces.
\let\@FMFPBIOORGparCMD=\par% save the original \par command
\edef\par{\hfil\break\indent}% the new \par will not be a "real" \par
\settoheight{\@FMFPtrantmpdimenA}{\@FMFPtempbiographybox}% get height of biography box
\@FMFPtrantmpdimenB=\@FMFPBIOhangdepth%
\@FMFPtrantmpcountA=\@FMFPtrantmpdimenB% countA has the hang depth
\divide\@FMFPtrantmpcountA by \baselineskip%  calculates lines needed to produce the hang depth
\advance\@FMFPtrantmpcountA by 1% ensure we overestimate
% set the hanging indent
\hangindent\@FMFPBIOhangwidth%
\hangafter-\@FMFPtrantmpcountA%
% reference the top of the photo area to the top of a capital T
\settoheight{\@FMFPtrantmpdimenB}{\mbox{T}}%
% set the photo box, give it zero width and height so as not to disturb anything
\noindent\makebox[0pt][l]{\hspace{-\@FMFPBIOhangwidth}\raisebox{\@FMFPtrantmpdimenB}[0pt][0pt]{%
\raisebox{-\@FMFPBIOphotodepth}[0pt][0pt]{\@FMFPtempbiographybox}}}%
% now place the author name and begin the bio text
\noindent\textbf{#2\ }\@FMFPgobbleleadPARNLSP}{\relax\let\par=\@FMFPBIOORGparCMD\par%
% 7/2001 V1.5 detect when the biography text is shorter than the photo area
% and pad the unused area - preventing a collision from the next biography entry
% MDS
\ifnum \prevgraf <\@FMFPtrantmpcountA\relax% detect when the biography text is shorter than the photo
    \advance\@FMFPtrantmpcountA by -\prevgraf% calculate how many lines we need to pad
    \advance\@FMFPtrantmpcountA by -1\relax% we compensate for the fact that we indented an extra line
    \@FMFPtrantmpdimenA=\baselineskip% calculate the length of the padding
    \multiply\@FMFPtrantmpdimenA by \@FMFPtrantmpcountA%
    \noindent\rule{0pt}{\@FMFPtrantmpdimenA}% insert an invisible support strut
\fi%
\par\normalfont}



% V1.6
% added biography without a photo environment
\newenvironment{FMFPbiographynophoto}[1]{%
% Make an entry into the table of contents only if we have not done so before
\if@FMFPbiographyTOCentrynotmade%
% link labels to the biography counter so hyperref will jump
% to the biography, not the previous section
\setcounter{FMFPbiography}{-1}%
\refstepcounter{FMFPbiography}%
\addcontentsline{toc}{section}{Biographies}%
\global\@FMFPbiographyTOCentrynotmadefalse%
\fi%
% one more biography
\refstepcounter{FMFPbiography}%
% Make an entry for this name into the table of contents 
\addcontentsline{toc}{subsection}{#1}%
\normalfont\@FMFPcompsoconly{\sffamily}\footnotesize\interlinepenalty500%
\vskip 4\baselineskip plus 1fil minus 0\baselineskip%
\parskip=0pt\par%
\noindent\textbf{#1\ }\@FMFPgobbleleadPARNLSP}{\relax\par\normalfont}


% provide the user with some old font commands
% got this from article.cls
\DeclareOldFontCommand{\rm}{\normalfont\rmfamily}{\mathrm}
\DeclareOldFontCommand{\sf}{\normalfont\sffamily}{\mathsf}
\DeclareOldFontCommand{\tt}{\normalfont\ttfamily}{\mathtt}
\DeclareOldFontCommand{\bf}{\normalfont\bfseries}{\mathbf}
\DeclareOldFontCommand{\it}{\normalfont\itshape}{\mathit}
\DeclareOldFontCommand{\sl}{\normalfont\slshape}{\@nomath\sl}
\DeclareOldFontCommand{\sc}{\normalfont\scshape}{\@nomath\sc}
\DeclareRobustCommand*\cal{\@fontswitch\relax\mathcal}
\DeclareRobustCommand*\mit{\@fontswitch\relax\mathnormal}


% SPECIAL PAPER NOTICE COMMANDS
% 
% holds the special notice text
\def\@FMFPspecialpapernotice{\relax}
 
% for special papers, like invited papers, the user can do:
% \FMFPspecialpapernotice{(Invited Paper)} before \maketitle
\def\FMFPspecialpapernotice#1{\ifCLASSOPTIONconference%
\def\@FMFPspecialpapernotice{{\Large#1\vspace*{1em}}}%
\else%
\def\@FMFPspecialpapernotice{{\\*[1.5ex]\sublargesize\textit{#1}}\vspace*{-2ex}}%
\fi}




% PUBLISHER ID COMMANDS
% to insert a publisher's ID footer
% V1.6 \FMFPpubid has been changed so that the change in page size and style
% occurs in \maketitle. \FMFPpubid must now be issued prior to \maketitle
% use \FMFPpubidadjcol as before - in the second column of the title page
% These changes allow \maketitle to take the reduced page height into
% consideration when dynamically setting the space between the author 
% names and the maintext.
%
% the amount the main text is pulled up to make room for the
% publisher's ID footer
% FMFP uses about 1.3\baselineskip for journals, 
% dynamic title spacing will clean up the fraction
\def\@FMFPpubidpullup{1.3\baselineskip}
\ifCLASSOPTIONtechnote
% for technotes it must be an integer of baselineskip as there can be no
% dynamic title spacing for two column mode technotes (the title is in the
% in first column) and we should maintain an integer number of lines in the
% second column
% There are some examples (such as older issues of "Transactions on
% Information Theory") in which FMFP really pulls the text off the ID for
% technotes - about 0.55in (or 4\baselineskip). We'll use 2\baselineskip
% and call it even.
\def\@FMFPpubidpullup{2\baselineskip}
\fi

% V1.7 compsoc does not use a pullup
\ifCLASSOPTIONcompsoc
\def\@FMFPpubidpullup{0pt}
\fi

% holds the ID text
\def\@FMFPpubid{\relax}

% flag so \maketitle can tell if \FMFPpubid was called
\newif\if@FMFPusingpubid
\global\@FMFPusingpubidfalse
% issue this command in the page to have the ID at the bottom
% V1.6 use before \maketitle
\def\FMFPpubid#1{\def\@FMFPpubid{#1}\global\@FMFPusingpubidtrue}


% command which will pull up (shorten) the column it is executed in
% to make room for the publisher ID. Place in the second column of
% the title page when using \FMFPpubid
% Is smart enough not to do anything when in single column text or
% if the user hasn't called \FMFPpubid
% currently needed in for the second column of a page with the
% publisher ID. If not needed in future releases, please provide this
% command and define it as \relax for backward compatibility
% v1.6b do not allow command to operate if the peer review option has been 
% selected because \FMFPpubidadjcol will not be on the cover page.
% V1.7 do nothing if compsoc
\def\FMFPpubidadjcol{\ifCLASSOPTIONcompsoc\else\ifCLASSOPTIONpeerreview\else
\if@twocolumn\if@FMFPusingpubid\enlargethispage{-\@FMFPpubidpullup}\fi\fi\fi\fi}

% Special thanks to Peter Wilson, Daniel Luecking, and the other
% gurus at comp.text.tex, for helping me to understand how best to
% implement the FMFPpubid command in LaTeX.



%% Lockout some commands under various conditions

% general purpose bit bucket
\newsavebox{\@FMFPtranrubishbin}

% flags to prevent multiple warning messages
\newif\if@FMFPWARNthanks
\newif\if@FMFPWARNFMFPPARstart
\newif\if@FMFPWARNFMFPbiography
\newif\if@FMFPWARNFMFPbiographynophoto
\newif\if@FMFPWARNFMFPpubid
\newif\if@FMFPWARNFMFPpubidadjcol
\newif\if@FMFPWARNFMFPmembership
\newif\if@FMFPWARNFMFPaftertitletext
\@FMFPWARNthankstrue
\@FMFPWARNFMFPPARstarttrue
\@FMFPWARNFMFPbiographytrue
\@FMFPWARNFMFPbiographynophototrue
\@FMFPWARNFMFPpubidtrue
\@FMFPWARNFMFPpubidadjcoltrue
\@FMFPWARNFMFPmembershiptrue
\@FMFPWARNFMFPaftertitletexttrue


%% Lockout some commands when in various modes, but allow them to be restored if needed
%%
% save commands which might be locked out
% so that the user can later restore them if needed
\let\@FMFPSAVECMDthanks\thanks
\let\@FMFPSAVECMDFMFPPARstart\FMFPPARstart
\let\@FMFPSAVECMDFMFPbiography\FMFPbiography
\let\@FMFPSAVECMDendFMFPbiography\endFMFPbiography
\let\@FMFPSAVECMDFMFPbiographynophoto\FMFPbiographynophoto
\let\@FMFPSAVECMDendFMFPbiographynophoto\endFMFPbiographynophoto
\let\@FMFPSAVECMDFMFPpubid\FMFPpubid
\let\@FMFPSAVECMDFMFPpubidadjcol\FMFPpubidadjcol
\let\@FMFPSAVECMDFMFPmembership\FMFPmembership
\let\@FMFPSAVECMDFMFPaftertitletext\FMFPaftertitletext


% disable \FMFPPARstart when in draft mode
% This may have originally been done because the pre-V1.6 drop letter
% algorithm had problems with a non-unity baselinestretch
% At any rate, it seems too formal to have a drop letter in a draft
% paper.
\ifCLASSOPTIONdraftcls
\def\FMFPPARstart#1#2{#1#2\if@FMFPWARNFMFPPARstart\typeout{** ATTENTION: \noexpand\FMFPPARstart
 is disabled in draft mode (line \the\inputlineno).}\fi\global\@FMFPWARNFMFPPARstartfalse}
\fi
% and for technotes
\ifCLASSOPTIONtechnote
\def\FMFPPARstart#1#2{#1#2\if@FMFPWARNFMFPPARstart\typeout{** WARNING: \noexpand\FMFPPARstart
 is locked out for technotes (line \the\inputlineno).}\fi\global\@FMFPWARNFMFPPARstartfalse}
\fi


% lockout unneeded commands when in conference mode
\ifCLASSOPTIONconference
% when locked out, \thanks, \FMFPbiography, \FMFPbiographynophoto, \FMFPpubid,
% \FMFPmembership and \FMFPaftertitletext will all swallow their given text. 
% \FMFPPARstart will output a normal character instead
% warn the user about these commands only once to prevent the console screen
% from filling up with redundant messages
\def\thanks#1{\if@FMFPWARNthanks\typeout{** WARNING: \noexpand\thanks
 is locked out when in conference mode (line \the\inputlineno).}\fi\global\@FMFPWARNthanksfalse}
\def\FMFPPARstart#1#2{#1#2\if@FMFPWARNFMFPPARstart\typeout{** WARNING: \noexpand\FMFPPARstart
 is locked out when in conference mode (line \the\inputlineno).}\fi\global\@FMFPWARNFMFPPARstartfalse}


% LaTeX treats environments and commands with optional arguments differently.
% the actual ("internal") command is stored as \\commandname 
% (accessed via \csname\string\commandname\endcsname )
% the "external" command \commandname is a macro with code to determine
% whether or not the optional argument is presented and to provide the 
% default if it is absent. So, in order to save and restore such a command
% we would have to save and restore \\commandname as well. But, if LaTeX
% ever changes the way it names the internal names, the trick would break.
% Instead let us just define a new environment so that the internal
% name can be left undisturbed.
\newenvironment{@FMFPbogusbiography}[2][]{\if@FMFPWARNFMFPbiography\typeout{** WARNING: \noexpand\FMFPbiography
 is locked out when in conference mode (line \the\inputlineno).}\fi\global\@FMFPWARNFMFPbiographyfalse%
\setbox\@FMFPtranrubishbin\vbox\bgroup}{\egroup\relax}
% and make biography point to our bogus biography
\let\FMFPbiography=\@FMFPbogusbiography
\let\endFMFPbiography=\end@FMFPbogusbiography

\renewenvironment{FMFPbiographynophoto}[1]{\if@FMFPWARNFMFPbiographynophoto\typeout{** WARNING: \noexpand\FMFPbiographynophoto
 is locked out when in conference mode (line \the\inputlineno).}\fi\global\@FMFPWARNFMFPbiographynophotofalse%
\setbox\@FMFPtranrubishbin\vbox\bgroup}{\egroup\relax}

\def\FMFPpubid#1{\if@FMFPWARNFMFPpubid\typeout{** WARNING: \noexpand\FMFPpubid 
 is locked out when in conference mode (line \the\inputlineno).}\fi\global\@FMFPWARNFMFPpubidfalse}
\def\FMFPpubidadjcol{\if@FMFPWARNFMFPpubidadjcol\typeout{** WARNING: \noexpand\FMFPpubidadjcol
 is locked out when in conference mode (line \the\inputlineno).}\fi\global\@FMFPWARNFMFPpubidadjcolfalse}
\def\FMFPmembership#1{\if@FMFPWARNFMFPmembership\typeout{** WARNING: \noexpand\FMFPmembership
 is locked out when in conference mode (line \the\inputlineno).}\fi\global\@FMFPWARNFMFPmembershipfalse}
\def\FMFPaftertitletext#1{\if@FMFPWARNFMFPaftertitletext\typeout{** WARNING: \noexpand\FMFPaftertitletext
 is locked out when in conference mode (line \the\inputlineno).}\fi\global\@FMFPWARNFMFPaftertitletextfalse}
\fi


% provide a way to restore the commands that are locked out
\def\FMFPoverridecommandlockouts{%
\typeout{** ATTENTION: Overriding command lockouts (line \the\inputlineno).}%
\let\thanks\@FMFPSAVECMDthanks%
\let\FMFPPARstart\@FMFPSAVECMDFMFPPARstart%
\let\FMFPbiography\@FMFPSAVECMDFMFPbiography%
\let\endFMFPbiography\@FMFPSAVECMDendFMFPbiography%
\let\FMFPbiographynophoto\@FMFPSAVECMDFMFPbiographynophoto%
\let\endFMFPbiographynophoto\@FMFPSAVECMDendFMFPbiographynophoto%
\let\FMFPpubid\@FMFPSAVECMDFMFPpubid%
\let\FMFPpubidadjcol\@FMFPSAVECMDFMFPpubidadjcol%
\let\FMFPmembership\@FMFPSAVECMDFMFPmembership%
\let\FMFPaftertitletext\@FMFPSAVECMDFMFPaftertitletext}



% need a backslash character for typeout output
{\catcode`\|=0 \catcode`\\=12
|xdef|@FMFPbackslash{\}}


% hook to allow easy disabling of all legacy warnings
\def\@FMFPlegacywarn#1#2{\typeout{** ATTENTION: \@FMFPbackslash #1 is deprecated (line \the\inputlineno).
Use \@FMFPbackslash #2 instead.}}


% provide for legacy commands
\def\authorblockA{\@FMFPlegacywarn{authorblockA}{FMFPauthorblockA}\FMFPauthorblockA}
\def\authorblockN{\@FMFPlegacywarn{authorblockN}{FMFPauthorblockN}\FMFPauthorblockN}
\def\authorrefmark{\@FMFPlegacywarn{authorrefmark}{FMFPauthorrefmark}\FMFPauthorrefmark}
\def\PARstart{\@FMFPlegacywarn{PARstart}{FMFPPARstart}\FMFPPARstart}
\def\pubid{\@FMFPlegacywarn{pubid}{FMFPpubid}\FMFPpubid}
\def\pubidadjcol{\@FMFPlegacywarn{pubidadjcol}{FMFPpubidadjcol}\FMFPpubidadjcol}
\def\QED{\@FMFPlegacywarn{QED}{FMFPQED}\FMFPQED}
\def\QEDclosed{\@FMFPlegacywarn{QEDclosed}{FMFPQEDclosed}\FMFPQEDclosed}
\def\QEDopen{\@FMFPlegacywarn{QEDopen}{FMFPQEDopen}\FMFPQEDopen}
\def\specialpapernotice{\@FMFPlegacywarn{specialpapernotice}{FMFPspecialpapernotice}\FMFPspecialpapernotice}



% provide for legacy environments
\def\biography{\@FMFPlegacywarn{biography}{FMFPbiography}\FMFPbiography}
\def\biographynophoto{\@FMFPlegacywarn{biographynophoto}{FMFPbiographynophoto}\FMFPbiographynophoto}
\def\keywords{\@FMFPlegacywarn{keywords}{FMFPkeywords}\FMFPkeywords}
\def\endbiography{\endFMFPbiography}
\def\endbiographynophoto{\endFMFPbiographynophoto}
\def\endkeywords{\endFMFPkeywords}


% provide for legacy IED commands/lengths when possible
\let\labelindent\FMFPlabelindent
\def\calcleftmargin{\@FMFPlegacywarn{calcleftmargin}{FMFPcalcleftmargin}\FMFPcalcleftmargin}
\def\setlabelwidth{\@FMFPlegacywarn{setlabelwidth}{FMFPsetlabelwidth}\FMFPsetlabelwidth}
\def\usemathlabelsep{\@FMFPlegacywarn{usemathlabelsep}{FMFPusemathlabelsep}\FMFPusemathlabelsep}
\def\iedlabeljustifyc{\@FMFPlegacywarn{iedlabeljustifyc}{FMFPiedlabeljustifyc}\FMFPiedlabeljustifyc}
\def\iedlabeljustifyl{\@FMFPlegacywarn{iedlabeljustifyl}{FMFPiedlabeljustifyl}\FMFPiedlabeljustifyl}
\def\iedlabeljustifyr{\@FMFPlegacywarn{iedlabeljustifyr}{FMFPiedlabeljustifyr}\FMFPiedlabeljustifyr}



% let \proof use the FMFPtran version even after amsthm is loaded
% \proof is now deprecated in favor of \FMFPproof
\AtBeginDocument{\def\proof{\@FMFPlegacywarn{proof}{FMFPproof}\FMFPproof}\def\endproof{\endFMFPproof}}

% V1.7 \overrideFMFPmargins is no longer supported.
\def\overrideFMFPmargins{%
\typeout{** WARNING: \string\overrideFMFPmargins \space no longer supported (line \the\inputlineno).}%
\typeout{** Use the \string\CLASSINPUTinnersidemargin, \string\CLASSINPUToutersidemargin \space controls instead.}}


\endinput

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End of FMFPtran.cls  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% That's all folks!

