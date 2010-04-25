% Do not edit this file; it is automatically
% generated from Documentation/snippets/new
% This file is in the public domain.
%% Note: this file works from version 2.13.5
\version "2.13.18"

\header {
%% Translation of GIT committish: 2b0dc29608d6c3f5a03ead4877ae514c647adb74
  texidoces = "
La instrucción @code{\\time} establece las propiedades
@code{timeSignatureFraction}, @code{beatLength}, @code{beatGrouping} y
@code{measureLength} en el contexto @code{Timing}, que normalmente
tiene el alias @code{Score}.  La modificación del valor de
@code{timeSignatureFraction} hace que se imprima la nueva indicación
de compás sin que cambie ninguna de las demás propiedades:

"
  doctitlees = "Cambio de compás sin afectar al barrado"



%% Translation of GIT committish: 0a868be38a775ecb1ef935b079000cebbc64de40
  texidocde = "
Der @code{\\time}-Befehl verändert die Eigenschaften
@code{timeSignatureFraction}, @code{beatLength}, @code{beatGrouping}
und @code{measureLength} im @code{Timing}-Kontext, welcher normalerweise
gleichbedeutend mit @code{Score} ist.  Wenn der Wert von
@code{timeSignatureFraction} verändert wird, wird die neue
Taktart ausgegeben, ohne die anderen Eigenschaften zu beeinflussen:

"
  doctitlede = "Die Taktart verändern ohne die Bebalkung zu beeinflussen"



%% Translation of GIT committish: 4da4307e396243a5a3bc33a0c2753acac92cb685
  texidocfr = "
La commande @code{\\time} gère les propriétés
@code{timeSignatureFraction}, @code{beatLength}, @code{beatGrouping}
et @code{measureLength} dans le contexte @code{Timing}, normallement
rattaché à @code{Score}.  Le fait de modifier la valeur de
@code{timeSignatureFraction} aura pour effet de changer l'apparence du
symbole affiché sans pour autant affecter les autres propriétés de la
métrique :

"
  doctitlefr = "Changement de métrique sans affecter les règles de ligature"

  lsrtags = "rhythms"
  texidoc = "
The @code{\\time} command sets the properties
@code{timeSignatureFraction}, @code{beatLength}, @code{beatGrouping}
and @code{measureLength} in the @code{Timing} context, which is
normally aliased to @code{Score}. Changing the value of
@code{timeSignatureFraction} causes the new time signature symbol to be
printed without changing any of the other properties:

"
  doctitle = "Changing the time signature without affecting the beaming"
} % begin verbatim


\markup {
  This snippet is deprecated as of 2.13.5 and will be removed in 2.14
}
