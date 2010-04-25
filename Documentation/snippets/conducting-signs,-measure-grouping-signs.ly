% Do not edit this file; it is automatically
% generated from Documentation/snippets/new
% This file is in the public domain.
%% Note: this file works from version 2.13.4
\version "2.13.18"

\header {
%% Translation of GIT committish: 2b0dc29608d6c3f5a03ead4877ae514c647adb74
  texidoces = "
La agrupación de pulsos dentro de un compás está controlada por la
agrupación predeterminada que se establece en @code{beamSettings}.
Estas agrupaciones se pueden establecer mediante el uso de
@code{\\overrideBeamSettings}.  De forma alternativa, se puede usar la
función de Scheme @code{set-time-signature} para establecer tanto el
compás como la regla de agrupamiento predeterminada.
@code{set-time-signature} acepta tres argumentos: el número de pulsos,
la longitud del pulso y la agrupación interna de los pulsos en el
compás.  Si el grabador @code{Measure_grouping_engraver} está incluido
en uno de los contextos de presentación, se imprimirán signos de
agrupación de pulsos.  Estos símbolos facilitan la lectura de música
moderna rítmicamente compleja.  En este ejemplo, el compás de 9/8 se
agrupa según dos patrones distintos utilizando los dos métodos,
mientras que el compás de 5/8 se agrupa de acuerdo con el ajuste
predeterminado que está en @file{scm/beam-settings.scm}:

"
doctitlees = "Símbolos de dirección y símbolos de agrupación de compás"


%% Translation of GIT committish: 0a868be38a775ecb1ef935b079000cebbc64de40
  texidocde = "
Optionen, mit denen die Balken in einem Takt gruppiert werden, sind
durch die Scheme-Funktion @code{set-time-signature} erhältlich, die
drei Argumente braucht:  Die Zahl der Taktschläge, die Länge des
Schlages und die interne gruppieren von Balken in dem Takt.  Wenn der
@code{Measure_grouping_engraver} hinzugefügt worden ist, erstellt
diese Funktion auch @code{MeasureGrouping}-(Taktgruppen)-Zeichen.  Derartige
Zeichen erleichtern das Lesen von rhythmisch komplexer Musik.  In dem
Beispiel ist der 9/8-Takt in 2, 2, 2 und 3 aufgeteilt.  Das wird
der @code{set-time-signature}-Funktion als das dritte Argument mitgegeben:
@code{'(2 2 2 3)}:

"
  doctitlede = "Dirigirzeichen Taktgruppenzeichen"



%% Translation of GIT committish: 3d7ffa1f82bb44673134b28becf7898482fe7316
  texidocfr = "
Les règles de ligature par mesure sont gérées par la propriété
@code{beamSettings}.  Elles peuvent être modifiées par la commande
@code{\\overrideBeamSettings}.  
Il existe des options qui permettent de grouper les ligatures au sein
d'une mesure, grâce à la fonction Scheme @code{set-time-signature}.
Celle-ci prend trois arguments : le nombre de pulsations, la durée de la
pulsation et le regroupement des pulsations dans la mesure.  Si l'on
fait appel au @code{Measure_grouping_engraver}, la fonction
@code{set-time-signature} créera aussi des symboles
@code{MeasureGrouping}.  Ces symboles aident à la lecture des œuvres
modernes à la rythmique complexe.  Dans l'exemple qui suit, la mesure à
9/8 est divisée en 2, 2, 2 et 3, alors que la mesure à 5/8 répond aux
règles par défaut contenues dans le fichier @w{@code{scm/beam-settings.scm}}.

"
  doctitlefr = "Signes de direction signes de sous-groupe"

  lsrtags = "rhythms"
  texidoc = "
Beat grouping within a bar is controlled by the default grouping
established in @code{beamSettings}.  This grouping can be established
by the use of @code{\\overrideBeamSettings}.  Alternatively, the
Scheme function @code{set-time-signature} can be used to both
set the time signature and establish the default grouping rule.
@code{set-time-signature}, takes three arguments: the
number of beats, the beat length, and the internal grouping of beats in
the measure.  If the @code{Measure_grouping_engraver} is included
in one of the display contexts, measure grouping signs will be
created.  Such signs ease reading rhythmically complex modern music.
In the example, the 9/8 measure is grouped in two different
patterns using the two different methods, while the 5/8 measure
is grouped according to the default setting in
@file{scm/beam-settings.scm}:
"
  doctitle = "Conducting signs measure grouping signs"
} % begin verbatim


\score {
  \relative c'' {
    \time 9/8
    \overrideBeamSettings #'Score #'(9 . 8) #'end #'((* . (2 2 2 3)))
    g8 g d d g g a( bes g) |
    #(set-time-signature 9 8 '(4 5))
    g8 g d d g g a( bes g) |
    \time 5/8
    a4. g4 |
  }
  \layout {
    \context {
      \Staff
      \consists "Measure_grouping_engraver"
    }
  }
}
