/*   
  output-property-engraver.cc --  implement Output_property_engraver
  
  source file of the GNU LilyPond music typesetter
  
  (c) 2000--2002 Han-Wen Nienhuys <hanwen@cs.uu.nl>
  
 */

#include "engraver.hh"
#include "grob.hh"
#include "output-property-music-iterator.hh"

class Output_property_engraver : public Engraver
{
TRANSLATOR_DECLARATIONS(Output_property_engraver);
protected:

  /*
    should do this with \once and \push ?


      \property Voice.outputProperties \push #pred = #modifier

      where both MODIFIER and PRED are functions taking a
      grob.
      
   */

  
  Link_array<Music> props_;

  virtual void stop_translation_timestep ();
  virtual void acknowledge_grob (Grob_info);
  virtual bool try_music (Music*);
};


bool
Output_property_engraver::try_music (Music* m)
{
  if (m->is_mus_type ("layout-instruction"))
    {
      props_.push (m);
      return true;
    }
  return false;
}

void
Output_property_engraver::acknowledge_grob (Grob_info inf)
{
  for (int i=props_.size (); i--;)
    {
      Music * o = props_[i];
      SCM pred = o->get_mus_property ("predicate");
      
      /*
	should typecheck pred. 
       */
      SCM result=gh_apply (pred,
			   scm_list_n (inf.grob_->self_scm (), SCM_UNDEFINED));
      if (to_boolean (result))
	{
	  SCM sym = o->get_mus_property ("grob-property");
	  SCM val = o->get_mus_property ("grob-value");
	  inf.grob_->internal_set_grob_property (sym, val);
	}
    }
}

void
Output_property_engraver::stop_translation_timestep ()
{
  props_.clear ();
}

Output_property_engraver::Output_property_engraver()
{
}

ENTER_DESCRIPTION(Output_property_engraver,
/* descr */       "Interpret Music of Output_property type, and apply a function "
" to any Graphic objects that satisfies the predicate.",
/* creats*/       "",
/* accepts */     "layout-instruction",
/* acks  */       "grob-interface",
/* reads */       "",
/* write */       "");
