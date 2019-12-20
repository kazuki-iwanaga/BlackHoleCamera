#!/usr/bin/env perl

@default_files    = ('main.tex');

$latex            = 'platex -halt-on-error';
$bibtex           = 'pbibtex';
$dvipdf           = 'dvipdfmx %O -o %D %S';
$makeindex        = 'mendex %O -o %D %S';
$max_repeat       = 5;
$pdf_mode         = 3;
