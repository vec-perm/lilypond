/*
  command-request.hh -- declare non-musical requests

  source file of the GNU LilyPond music typesetter

  (c)  1997--2000 Han-Wen Nienhuys <hanwen@cs.uu.nl>
*/


#ifndef COMMANDREQUEST_HH
#define COMMANDREQUEST_HH

#include "request.hh"
#include "array.hh"
#include "duration.hh"
#include "musical-pitch.hh"
#include "protected-scm.hh"

/*
  Real penalty_f_;
 */
class Break_req : public Request {
public:

  Break_req ();
protected:
  VIRTUAL_COPY_CONS(Music);
};

class Mark_req : public Request {
public:
  virtual bool do_equal_b (Request const*) const;
  SCM mark_label ();
  VIRTUAL_COPY_CONS(Music);
};


/** Baseclass for time_signature/partial req. It has to be handled by
  Staff_{walker,column} baseclass.  */
class Timing_req  : public Request  {
public:
  VIRTUAL_COPY_CONS(Music);
};

/*
    int metronome_i_;
 */
class Tempo_req : public Timing_req
{
public:
  Duration dur_;


  Tempo_req();
protected:

  VIRTUAL_COPY_CONS(Music);
  bool do_equal_b (Request const *) const;
};


/**
  todo: allow C time_signature

  int beats_i_;
  int one_beat_i_;
  
 */
class Time_signature_change_req  : public Timing_req  {
public:
  Time_signature_change_req();

protected:
  bool do_equal_b (Request const *) const;
  VIRTUAL_COPY_CONS(Music);
};


/// check if we're at start of a  measure.
class Barcheck_req  : public Timing_req  {
public:
  bool do_equal_b (Request const *) const;
  VIRTUAL_COPY_CONS(Music);
};


/** draw a (repeat)-bar. This something different than #Barcheck_req#,
  the latter should only happen at the start of a measure.  */
class Bar_req  : public Request  {
public:

  Bar_req (String);
protected:
  VIRTUAL_COPY_CONS(Music);
};

class Breathing_sign_req : public Request {
  VIRTUAL_COPY_CONS(Music);
};

/**
    Handle key changes.
*/
class Key_change_req  : public Request
{
public:
  SCM pitch_alist ();
  
protected:
  VIRTUAL_COPY_CONS(Music);
  void transpose (Musical_pitch  d);
  bool do_equal_b (Request const * )const; 
};

/*
  String clef_str_;
 */

class Clef_change_req  : public Request  {
public:
  
  Clef_change_req ();
protected:

  VIRTUAL_COPY_CONS(Music);
};


#endif // COMMANDREQUEST_HH
