/*
 * This file is part of SCHUR.
 *
 * SCHUR - an interactive program for calculating properties of Lie
 * groups and symmetric functions.
 * Copyright (C) 2006  Franck BUTELLE, Frédéric Toumazet
 * 
 * SCHUR is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * SCHUR is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with SCHUR; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/** \file label.c
 */

#include "standard.h"
#include "define.h"
#include "dim.h"
#include "type.h"

// for osort
#include "s1.h"
#include "label.h"

void
dlabel (ocharptr list)
{
  ocharptr ptr;
  ptr = list;
  while (ptr != NULL)
    {
      if ((ptr->lab != ' '))
	ptr->lab = ' ';
      ptr = ptr->next;
    }
  osort (&ptr, true);
}

void
olabel (ocharptr list)
{
  ocharptr ptr;

  ptr = list;
  while (ptr != NULL)
    {
      if ((ptr->mult < 0))
	if (ptr->lab == '#')
	  ptr->lab = ' ';
	else
	  ptr->lab = '#';
      ptr = ptr->next;
    }
}

void
slabel (termptr list)
{
  termptr ptr;

  ptr = list;
  while (ptr != NULL)
    {
      if ((bool) ((wtfrm (&ptr->val)) & 1))
	ptr->slab = '#';
      else
	ptr->slab = ' ';
      ptr = ptr->next;
    }
}

void
oclabel (ocharptr list)
{
  ocharptr ptr;

  ptr = list;
  while (ptr != NULL)
    {
      if ((bool) ((wtfrm (&ptr->val)) & 1))
	ptr->lab = '#';
      else
	ptr->lab = ' ';
      ptr = ptr->next;
    }
}
void
cslabel (termptr list, char ch)
{
  termptr ptr;

  ptr = list;
  while (ptr != NULL)
    {
      ptr->slab = ch;
      ptr = ptr->next;
    }
}
