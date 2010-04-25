%% Do not edit this file; it is automatically
%% generated from LSR http://lsr.dsi.unimi.it
%% This file is in the public domain.
\version "2.13.16"

\header {
  lsrtags = "repeats, staff-notation, editorial-annotations"

%% Translation of GIT committish: 2b0dc29608d6c3f5a03ead4877ae514c647adb74
  texidoces = "
Este fragmento de código proporciona una solución alternativa a la
producción de contadores de compás utilizando repeticiones
transparentes de tipo porcentaje.

"
  doctitlees = "Contador de compases"

  texidoc = "
This snippet provides a workaround for emitting measure counters using
transparent percent repeats.

"
  doctitle = "Measure counter"
} % begin verbatim

<<
  \context Voice = "foo" {
    \clef bass
    c4 r g r
    c4 r g r
    c4 r g r
    c4 r g r
  }
  \context Voice = "foo" {
    \set countPercentRepeats = ##t
    \override PercentRepeat #'transparent = ##t
    \override PercentRepeatCounter #'staff-padding = #1
    \repeat percent 4 { s1 }
  }
>>

