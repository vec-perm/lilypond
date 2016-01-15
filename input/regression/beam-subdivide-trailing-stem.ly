\version "2.19.34"

\header {

  doctitle = "Leading and trailing stems"

  texidoc = "If in a subdivided beam one single stem follows a subdivision
the beam count should reflect the beam count of the subdivision as usual.
That is, the beam count should not be increased according to the remaining
length of the beam. The appended single stem has beamlets to the left. The
same is true for leading stems at the beginning that are immediately
followed by a subdivision, here the appropriate number of flags points
to the right."

}

\relative c' {
  \time 1/4
  \set subdivideBeams = ##t
  \set baseMoment = #(ly:make-moment 1/16)
  c32 [ c c c c32 ] r16.
  c32 [ c c c c64 ] r32. r16
  c32 [ c c32 ] r32 r8
  c32 [ c c64 ] r32. r8

  r16. c32 [ c c c c32 ]
  r16.. c64 [ c32 c c c ]
  r8 r32 c32 [ c c32 ]
  r8 r32. c64 [ c32 c ]

}