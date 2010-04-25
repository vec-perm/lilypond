%% Do not edit this file; it is automatically
%% generated from LSR http://lsr.dsi.unimi.it
%% This file is in the public domain.
\version "2.13.16"

\header {
  lsrtags = "rhythms, tweaks-and-overrides"

%% Translation of GIT committish: 2b0dc29608d6c3f5a03ead4877ae514c647adb74
 doctitlees = "Posicionar los silencios multicompás"
 texidoces = "
A diferencia de los silencios normales, no existe una instrucción
predefinida para modificar la posición predeterminada de un
símbolo de silencio multicompás sobre el pentagrama, adjuntándolo
a una nota, independientemente de cuál sea su forma.  Sin embargo,
en la música polifónica los silencios multicompás de las voces de
numeración par e impar están separados verticalmente.  La
colocación de los silencios multicompás se puede controlar como se
ve a continuación:

"
%% Translation of GIT committish: 0a868be38a775ecb1ef935b079000cebbc64de40
texidocde = "
Anders als bei normalen Pausen gibt es keinen direkten Befehl, um die
vertikale Position von Ganztaktpausen zu beeinflussen, indem man sie an
eine Tonhöhe anhängt.  In polyphoner Notation wird aber dennoch die
Position der Pausen von geraden und ungeraden Stimmen voneinander
unterschieden.  Die Position von Ganztaktpausen kann wie folgt verändert
werden:
 "
  doctitlede = "Positionierung von Ganztaktpausen"


%% Translation of GIT committish: 4da4307e396243a5a3bc33a0c2753acac92cb685
  texidocfr = "
Si l'on peut positionner verticalement un silence simple en le
rattachant à une note, il n'en va pas de même pour un silence
multi-mesures.  Néanmoins, et uniquement dans le cadre de musique
polyphonique, les silences multi-mesures sont positionnés différemment
selon qu'ils appartiennent à une voix au numéro pair ou impair.  Le
positionnement des silences multi-mesures peut se contrôler ainsi :
"
  doctitlefr = "Positionnement des silences multi-mesures"

  texidoc = "
Unlike ordinary rests, there is no predefined command to change the
staff position of a multi-measure rest symbol of either form by
attaching it to a note.  However, in polyphonic music multi-measure
rests in odd-numbered and even-numbered voices are vertically
separated. The positioning of multi-measure rests can be controlled as
follows:

"
  doctitle = "Positioning multi-measure rests"
} % begin verbatim

\relative c'' {
  % MMR - Multi-Measure Rest
  % MMRs by default are set under the fourth line
  R1
  % They can be moved with an override
  \override MultiMeasureRest #'staff-position = #-2
  R1
  % A value of 0 is the default position;
  % the following trick moves the rest to the center line
  \override MultiMeasureRest #'staff-position = #-0.01
  R1
  % MMRs in odd-numbered voices are under the top line
  << { R1 } \\ { a1 } >>
  % MMRs in even-numbered voices are under the bottom line
  << { c1 } \\ { R1 } >>
  % They remain separated even in empty measures
  << { R1 } \\ { R1 } >>
  % This brings them together even though there are two voices
  \compressFullBarRests
  <<
    \revert MultiMeasureRest #'staff-position
    { R1*3 }
    \\
    \revert MultiMeasureRest #'staff-position
    { R1*3 }
  >>
}
