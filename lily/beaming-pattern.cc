/*
  This file is part of LilyPond, the GNU music typesetter.

  Copyright (C) 1999--2015 Han-Wen Nienhuys <hanwen@xs4all.nl>

  LilyPond is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  LilyPond is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with LilyPond.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "context.hh"
#include "beaming-pattern.hh"
#include "misc.hh"

/*
  Represents a stem belonging to a beam. Sometimes (for example, if the stem
  belongs to a rest and stemlets aren't used) the stem will be invisible.

  The rhythmic_importance_ of an element tells us the significance of the
  moment at which this element occurs. For example, an element that occurs at
  a beat is more significant than one that doesn't. Smaller number are
  more important. The rhythmic_importance_ is decided and filled in by
  Beaming_pattern. A rhythmic_importance_ smaller than zero has extra
  significance: it represents the start of a beat and therefore beams may
  need to be subdivided.
*/
Beam_rhythmic_element::Beam_rhythmic_element ()
{
  start_moment_ = 0;
  rhythmic_importance_ = 0;
  beam_count_drul_[LEFT] = 0;
  beam_count_drul_[RIGHT] = 0;
  beam_count_for_length_ = 0;
  subdivisions_[LEFT] = false;
  subdivisions_[RIGHT] = false;
  invisible_ = false;
  factor_ = Rational (1);
  tuplet_start_ = false;
}

Beam_rhythmic_element::Beam_rhythmic_element (Moment m, int i, bool inv,
                                              Rational factor, bool tuplet_start)
{
  start_moment_ = m;
  rhythmic_importance_ = 0;
  beam_count_drul_[LEFT] = i;
  beam_count_drul_[RIGHT] = i;
  beam_count_for_length_ = i;
  subdivisions_[LEFT] = false;
  subdivisions_[RIGHT] = false;
  invisible_ = inv;
  factor_ = factor;
  tuplet_start_ = tuplet_start;
}

void
Beam_rhythmic_element::de_grace ()
{
  if (start_moment_.grace_part_)
    {
      start_moment_.main_part_ = start_moment_.grace_part_;
      start_moment_.grace_part_ = 0;
    }
}

int
Beam_rhythmic_element::count (Direction d) const
{
  return beam_count_drul_[d];
}

/*
  Finds the appropriate direction for  If
  the stem has no more flags than either of its neighbours, this returns
  CENTER.
  Do not call this with 0 or stems_.size ()!
*/
Direction
Beaming_pattern::flag_direction (vsize i,
                                 Drul_array<int>& beam_counts,
                                 bool strict_beat_beaming) const
{

  if (stems_[i - 1].count (RIGHT) >= beam_counts[LEFT] &&
      stems_[i + 1].count (LEFT) >= beam_counts[RIGHT])
    // No neighbor has more beams
    return CENTER;

  if (!strict_beat_beaming)
    {
      // Try to avoid sticking-out flags as much as possible by pointing
      // my flags at the neighbor with the most flags.
      if (stems_[i - 1].count (RIGHT) > stems_[i + 1].count (LEFT))
        return LEFT;
      if (stems_[i - 1].count (RIGHT) < stems_[i + 1].count (LEFT))
        return RIGHT;
    }

  //  => (strict_beat_beaming || beam_counts[RIGHT] == beam_counts[LEFT])
  // Point away from the stem with higher rhythmic_importance
  return
    (stems_[i].rhythmic_importance_ == stems_[i + 1].rhythmic_importance_)
    ? CENTER
    : (stems_[i].rhythmic_importance_ < stems_[i + 1].rhythmic_importance_)
      ? RIGHT : LEFT;
}

void
Beaming_pattern::de_grace ()
{
  for (vsize i = 0; i < stems_.size (); i++)
    {
      stems_[i].de_grace ();
    }
}

/*
  Process the beaming pattern and determine beam count
  right and left of all stems.  Respect strictBeatBeaming and subdivideBeams.
*/
void
Beaming_pattern::beamify (Beaming_options const &options)
{
  if (stems_.size () <= 1)
    return;

  // Treat grace-d beams like normal durations,
  // only necessary for subdivided beams
  if (options.subdivide_beams_ && stems_[0].start_moment_.grace_part_)
    de_grace ();

  // "shift" beams in partial bars by one measure to a positive start_moment
  if (stems_[0].start_moment_ < Moment (0))
    for (vsize i = 0; i < stems_.size (); i++)
      stems_[i].start_moment_ += options.measure_length_;

  // If beams are subdivided we'll handle beam count surrounding
  // invisible stems individually.
  if (!options.subdivide_beams_)
    unbeam_invisible_stems ();

  find_rhythmic_importance (options);

  Drul_array<int> beam_counts;
  Direction flag_dir = CENTER;


  // Iterate over stems, excluding the extremals
  for (vsize i = 1; i < stems_.size () - 1; i++)
    {
      // Initialize beam count based on duration of the current stem
      beam_counts[LEFT] = stems_[i].beam_count_for_length_;
      beam_counts[RIGHT] = beam_counts[LEFT];

      // Process beamlets around subdivisions
      if (options.subdivide_beams_ && find_subdivisions (i))
        {
          if (stems_[i].subdivisions_[RIGHT])
            continue;
          else
            {
              // There is a subdivision to the left of us.
              // One beam must be left, even for beams >= 1/4.
              beam_counts[LEFT] = max (1, beam_count_for_subdivision (i));
              stems_[i - 1].beam_count_drul_[RIGHT] = beam_counts[LEFT];
            }
        }
      // Process stems not adjacent to a subdivision
      else
        {
          // Determine if one side of the stem has more beams
          flag_dir = flag_direction (i,
                                     beam_counts,
                                     options.strict_beat_beaming_);

          if (flag_dir)
            {
              // Remove beamlets from the non-flag side if necessary
              if (options.strict_beat_beaming_ &&
                  options.subdivide_at_strict_beat_beaming_)
                    // Force a subdivision at the non-flag side
                    {
                      int subdiv_beam_count = (flag_dir == LEFT)
                        ? beam_count_for_subdivision (i + 1)
                        : beam_count_for_subdivision (i);
                      beam_counts[-flag_dir] = subdiv_beam_count;

                      // If the subdivision is left of the current stem
                      // we have to fix the previous stem's beam count
                      if (flag_dir == RIGHT)
                        stems_[i - 1].beam_count_drul_[RIGHT] = subdiv_beam_count;
                }
              else
                // take beam count from the neighbour but not more
                // than appropriate for the current length
                beam_counts[-flag_dir] = min (
                            stems_[i - flag_dir].count (flag_dir),
                            stems_[i].beam_count_for_length_);
            }
        }

      // Ensure that the next stem can handle its left side properly
      // when the current stem is invisible (issue #4739)
      if (stems_[i].invisible_ && !options.strict_beat_beaming_)
          beam_counts[RIGHT] = beam_counts[LEFT];

      // Apply beam counts
      stems_[i].beam_count_drul_[LEFT] = beam_counts[LEFT];
      stems_[i].beam_count_drul_[RIGHT] = beam_counts[RIGHT];
    }

    fix_extremal_stems (options);
}

/*
  Apply special treatment to the first and last stems
  which could not be done within the loop.
*/
void
Beaming_pattern::fix_extremal_stems (Beaming_options const &options)
{
    // Special treatment for extremal stems
    vsize last = stems_.size () - 1;

    // Trailing stem after subdivision
    if (options.subdivide_beams_ && stems_[last].rhythmic_importance_ < 0)
      stems_[last - 1].beam_count_drul_[RIGHT] =
        beam_count_for_rhythmic_position (last);

    // Leading stem before subdivision
    if (options.subdivide_beams_ && stems_[1].rhythmic_importance_ < 0)
      stems_[0].beam_count_drul_[RIGHT] =
        beam_count_for_length (start_moment (1) - start_moment (0));

    // Rests
    // TODO: Discuss if having *one* beam for this is the correct solution
    if (stems_[last].invisible_ && !options.strict_beat_beaming_)
      {
        stems_[last].beam_count_drul_[LEFT] = 1;
        stems_[last - 1].beam_count_drul_[RIGHT] = 1;
      }
    if (stems_[0].invisible_ && !options.strict_beat_beaming_)
      {
        stems_[0].beam_count_drul_[RIGHT] = 1;
        stems_[1].beam_count_drul_[LEFT] = 1;
      }
}

/*
   Set the tuplet start moment as necessary
*/
void
update_tuplet (Moment start_moment, Rational factor, Moment *tuplet_start_moment)
{
  int tuplet_number = (int) factor.den ();
  if ((tuplet_number > 1) && (tuplet_start_moment->num () < 0))
    *tuplet_start_moment = start_moment;
  else if (tuplet_number == 1)
    *tuplet_start_moment = Moment (-1, 1);
}

/*
   Get the group start position, the next group starting position, and the
   next beat starting position, given start_moment, base_moment,
   grouping, and factor
*/
void
find_location (SCM grouping, Moment base_moment, Moment start_moment,
               Rational factor, Moment *group_pos, Moment *next_group_pos,
               Moment *next_beat_pos)
{
  *group_pos = Moment (0);
  *next_group_pos = Moment (0);
  *next_beat_pos = base_moment;

  while (*next_beat_pos <= start_moment)
    *next_beat_pos += base_moment;

  while (*next_group_pos < *next_beat_pos)
    {
      I64 group_count = 1;  //default -- 1 base moments in a beam
      if (scm_is_pair (grouping))
        {
          group_count = scm_to_int (scm_car (grouping));
          grouping = scm_cdr (grouping);
        }

      // If we have a tuplet, the count should be determined from
      // the maximum tuplet size for beamed tuplets.
      U64 tuplet_number = factor.den ();
      if (tuplet_number > 1U)
        {
          // We use 1/8 as the base moment for the tuplet because it's
          // the largest beamed value.  If the tuplet is shorter, it's
          // OK, the code still works
          I64 test_count = ( Moment (Rational (1, 8) / factor) / base_moment).num ();
          if (test_count > group_count) group_count = test_count;
        }
      *group_pos = *next_group_pos;
      *next_group_pos = *group_pos + Rational(group_count) * base_moment;
    }
}

void
Beaming_pattern::find_rhythmic_importance (Beaming_options const &options)
{
  Moment group_pos (0);  // 0 is the start of the first group
  Moment next_group_pos (0);
  Moment next_beat_pos (options.base_moment_);
  Moment tuplet_start_moment (-1, 1);
  I64 tuplet_number = 1;

  SCM grouping = options.grouping_;
  vsize i = 0;

  // Find where we are in the beat structure of the measure
  if (stems_.size ())
    find_location (grouping, options.base_moment_, stems_[i].start_moment_,
                   stems_[i].factor_, &group_pos, &next_group_pos, &next_beat_pos);

  // Mark the importance of stems that start at a beat or a beat group.
  while (i < stems_.size ())
    {
      if ((next_beat_pos > next_group_pos)
          || (stems_[i].start_moment_ > next_beat_pos))
        // Find the new group ending point
        find_location (grouping, options.base_moment_, stems_[i].start_moment_,
                       stems_[i].factor_, &group_pos, &next_group_pos, &next_beat_pos);
      // Mark the start of this beat group
      if (stems_[i].start_moment_ == group_pos)
        stems_[i].rhythmic_importance_ = -2;
      // Work through the end of the beat group or the end of the beam
      while (i < stems_.size () && stems_[i].start_moment_ < next_group_pos)
        {
          // Set the tuplet start as necessary
          update_tuplet (stems_[i].start_moment_, stems_[i].factor_, &tuplet_start_moment);
          Moment dt = stems_[i].start_moment_ - group_pos;
          Rational tuplet = stems_[i].factor_;
          Moment tuplet_moment (tuplet);
          Moment tuplet_dt = stems_[i].start_moment_ - tuplet_start_moment;
          tuplet_number = tuplet.den ();
          // set the beat end and increment the next beat
          if (stems_[i].start_moment_ == next_beat_pos)
            {
              stems_[i].rhythmic_importance_ = -1;
              next_beat_pos += options.base_moment_;
            }
          // The rhythmic importance of a stem between beats depends on its fraction
          // of a beat: those stems with a lower denominator are deemed more
          // important.  For tuplets, we need to make sure that we use
          // the fraction of the tuplet, instead of the fraction of
          // a beat.
          Moment ratio = (tuplet_number == 1)
                         ? dt / options.base_moment_
                         : tuplet_dt / Moment (1, 8) / tuplet_moment;
          if (stems_[i].rhythmic_importance_ >= 0)
            stems_[i].rhythmic_importance_ = (int) ratio.den ();

          i++;
        }

      if (i < stems_.size () && stems_[i].start_moment_ == next_beat_pos)
        {
          if (tuplet_number == 1)
            stems_[i].rhythmic_importance_ = -1;
          next_beat_pos += options.base_moment_;
          if (stems_[i].start_moment_ == next_group_pos)
            stems_[i].rhythmic_importance_ = -2;
        }
    }
}

/*
  Invisible stems should be treated as though they have the same number of
  beams as their least-beamed neighbour. Here we go through the stems and
  modify the invisible stems to satisfy this requirement.
*/
void
Beaming_pattern::unbeam_invisible_stems ()
{
  for (vsize i = 1; i < stems_.size (); i++)
    if (stems_[i].invisible_)
      {
        int b = min (stems_[i].count (LEFT), stems_[i - 1].count (LEFT));
        stems_[i].beam_count_drul_[LEFT] = b;
        stems_[i].beam_count_drul_[RIGHT] = b;
      }

  if (stems_.size () > 1)
    for (vsize i = stems_.size () - 1; i--;)
      if (stems_[i].invisible_)
        {
          int b = min (stems_[i].count (LEFT), stems_[i + 1].count (LEFT));
          stems_[i].beam_count_drul_[LEFT] = b;
          stems_[i].beam_count_drul_[RIGHT] = b;
        }
}

void
Beaming_pattern::add_stem (Moment m, int b, bool invisible, Rational factor, bool tuplet_start)
{
  stems_.push_back (Beam_rhythmic_element (m, b, invisible, factor, tuplet_start));
}

Beaming_pattern::Beaming_pattern ()
{
}

/*
  Set flags for subdivisions on either side of the given stem.
  Return true if either of them evaluates to true.
*/
bool
Beaming_pattern::find_subdivisions (int i)
{
  stems_[i].subdivisions_[LEFT] = (stems_[i].rhythmic_importance_ < 0);
  stems_[i].subdivisions_[RIGHT] = (stems_[i + 1].rhythmic_importance_ < 0);
  return (stems_[i].subdivisions_[LEFT] || stems_[i].subdivisions_[RIGHT]);
}

int
Beaming_pattern::beamlet_count (int i, Direction d) const
{
  return stems_.at (i).beam_count_drul_[d];
}

Moment
Beaming_pattern::start_moment (int i) const
{
  return stems_.at (i).start_moment_;
}

Moment
Beaming_pattern::end_moment (int i) const
{
  Duration dur (2 + max (beamlet_count (i, LEFT),
                         beamlet_count (i, RIGHT)),
                0);

  return stems_.at (i).start_moment_
         + stems_.at (i).factor_ * dur.get_length ();
}

Moment
Beaming_pattern::remaining_length (int i) const
{
// The following loop isn't currently used.
// if not needed it has to be removed before merging
  int next_rest (0);
  for (int j = i + 1; j < stems_.size (); j++)
    {
      if (stems_[j].invisible_)
        {
          next_rest = j;
          break;
        }
    }
    return //(next_rest)
          // ? stems_[next_rest].start_moment_ - stems_[i].start_moment_
        //   :
           end_moment (stems_.size () - 1) - stems_[i].start_moment_;
}

int
Beaming_pattern::beam_count_for_rhythmic_position (int idx) const
{
    // Calculate number of beams representing the rhythmic position of given stem
    return intlog2(stems_[idx].start_moment_.main_part_.den()) - 2;
}

int
Beaming_pattern::beam_count_for_length (Moment len) const
{
    return intlog2(len.main_part_.den()) - 2 - intlog2(len.main_part_.num());
}

/*
   Returns the number of beams the given subdivision should have if the given
   stem were the right part of a subdivision. Respects the remaining length
   of the stem (e.g. 1/16 <= remaining length < 1/8 returns 2).
   Note that this gives reasonable results for *any* stem, not only for those
   actually at a subdivision.
*/
int
Beaming_pattern::beam_count_for_subdivision (vsize i) const
{
  return max (beam_count_for_rhythmic_position (i),
              beam_count_for_length (remaining_length (i)));
}

bool
Beaming_pattern::invisibility (int i) const
{
  return stems_.at (i).invisible_;
}

Rational
Beaming_pattern::factor (int i) const
{
  return stems_.at (i).factor_;
}

bool
Beaming_pattern::tuplet_start (int i) const
{
  return stems_.at (i).tuplet_start_;
}

/*
    Split a beaming pattern at index i and return a new
    Beaming_pattern containing the removed elements
*/
Beaming_pattern *
Beaming_pattern::split_pattern (int i)
{
  Beaming_pattern *new_pattern = 0;
  int count;

  new_pattern = new Beaming_pattern ();
  for (vsize j = i + 1; j < stems_.size (); j++)
    {
      count = max (beamlet_count (j, LEFT), beamlet_count (j, RIGHT));
      new_pattern->add_stem (start_moment (j),
                             count,
                             invisibility (j),
                             factor (j),
                             tuplet_start (j));
    }
  for (vsize j = i + 1; j < stems_.size ();)
    stems_.pop_back ();
  return (new_pattern);
}

void
Beaming_options::from_context (Context *context)
{
  grouping_ = context->get_property ("beatStructure");
  subdivide_beams_ = to_boolean (context->get_property ("subdivideBeams"));
  strict_beat_beaming_ = to_boolean (context->get_property ("strictBeatBeaming"));
  subdivide_at_strict_beat_beaming_ = to_boolean (context->get_property
                                                   ("subdivideAtStrictBeatBeaming"));
  base_moment_ = robust_scm2moment (context->get_property ("baseMoment"),
                                    Moment (1, 4));
  measure_length_ = robust_scm2moment (context->get_property ("measureLength"),
                                       Moment (4, 4));
}

Beaming_options::Beaming_options ()
{
//  grouping_ = SCM_EOL;
//  subdivide_beams_ = false;
//  strict_beat_beaming_ = false;
//  subdivide_at_strict_beat_beaming_ = true;
}
