/*
 * This file is part of SCHUR.
 *
 * SCHUR - an interactive program for calculating properties of Lie
 * groups and symmetric functions.
 * Copyright (C) 2006  Franck BUTELLE, Frederic Toumazet
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

# include <stdio.h>
# include <strings.h>
# include <stdlib.h>

# include "hivesLRcoef.h"

#define MAXS 100
#define MYINFINITY 200

typedef struct
{
  int x;
  int y;
} point_t;

typedef struct
{
  int min[MAXS][MAXS];
  int max[MAXS][MAXS];
} intervals;

void displayPart (frame *part, unsigned n);
void displayTabs (int min[MAXS][MAXS], int max[MAXS][MAXS], unsigned n);
void displayResults (int min[MAXS][MAXS], unsigned n);
bool ValidSumPartitions (frame *, frame *, frame *, unsigned n);
void initTabs (int min[MAXS][MAXS], int max[MAXS][MAXS], frame *, frame *, frame *,
	       unsigned n, unsigned maxv);
bool updateMin (int min[MAXS][MAXS], int max[MAXS][MAXS], int i, int j, int v,
		bool * updated, char *rule);
bool updateMax (int min[MAXS][MAXS], int max[MAXS][MAXS], int i, int j, int v,
		bool * updated, char *rule);
bool validHivesIJ (int min[MAXS][MAXS], int max[MAXS][MAXS], unsigned n,
		   unsigned i, unsigned j);
bool validHives (int min[MAXS][MAXS], int max[MAXS][MAXS], unsigned n);
point_t hives (int min[MAXS][MAXS], int max[MAXS][MAXS], unsigned n,
	       bool *updated);
unsigned tryValue (point_t, intervals, unsigned, unsigned, bool);

/*
unsigned
strtopart (char *str, int *part, unsigned dilat)
{
  int i = 1;			// traditionnaly partition indices start at 1.
  char *c, *start = str;

  while ((c = index (start, ',')) != NULL)
    {
      part[i++] = dilat*strtol (start, (char **) NULL, 10);

      start = c + 1;
    }
  part[i] = dilat*strtol (start, (char **) NULL, 10);

  return (i);
}

void
addZeros (int *part, unsigned n1, unsigned n2)
{
  int i;

  for (i = n1 + 1; i <= n2; i++)
    part[i] = 0;
}
*/

void
displayPart (frame *part, unsigned n)
{
  unsigned i;

  printf ("<");
  for (i = 1; i < n; i++)
    printf ("%d,", part->A[i]);
  printf ("%d>\n", part->A[n]);
}

void
displayTabs (int min[MAXS][MAXS], int max[MAXS][MAXS], unsigned n)
{
  unsigned i, j;
  for (i = 0; i <= n; i++)
    {
      for (j = 0; j <= n; j++)
	{
	  if (min[i][j] == MYINFINITY)
	    printf ("%6s", ".");
	  else if (min[i][j] == max[i][j])
	    printf ("%6d", min[i][j]);
	  else
	    printf ("%3d-%2d", min[i][j], max[i][j]);
	}
      printf ("\n");
    }
}

void
displayResults (int min[MAXS][MAXS], unsigned n)
{
  unsigned i, j;

  printf ("(");
/*  for (i = 1; i < n; i++)
    for (j = 1; i + j < n; j++)
      if (i == n - 2 && j == 1)*/
  for (i = 0; i <= n ; i++)
	 for (j = 0; i+j <= n ; j++)
	 {
      if ((i == n) && (i+j==n))
	printf ("%d)\n", min[i][j]);
      else
	printf ("%d,", min[i][j]);
	 }
}


bool
ValidSumPartitions (frame *part1, frame *part2, frame *part3, unsigned n)
{
  unsigned sum1, sum2, sum3;
// we also have to test no increasing order.

  sum1 = sumOfPartition (part1, n);
  sum2 = sumOfPartition (part2, n);
  sum3 = sumOfPartition (part3, n);

  return (sum1 + sum2 == sum3);
}


void
initTabs (int min[MAXS][MAXS], int max[MAXS][MAXS], frame *lambda, frame *mu,
	  frame *nu, unsigned n, unsigned maxv)
{
  unsigned i, j;
// a[i][j] 0<=i,j<=n;
//
  min[0][0] = max[0][0] = 0;

  for (i = 1; i <= n; i++)	//default values a kind of infinity.
    for (j = 1; j <= n; j++)
      {
	min[i][j] = MYINFINITY;
	max[i][j] = MYINFINITY;
      }

  for (i = 1; i <= n; i++)
    {
      min[0][i] = max[0][i] = lambda->A[i] + min[0][i - 1];
      min[i][0] = max[i][0] = nu->A[i] + min[i - 1][0];
    }

  for (i = 1; i <= n; i++)
    min[i][n - i] = max[i][n - i] = mu->A[i] + min[i - 1][n - i + 1];

  for (i = 1; i < n; i++)	//points inside triangle
    for (j = 1; j < n - i; j++)
      {
	min[i][j] = 0;
	max[i][j] = maxv;
      }
}

bool
updateMin (int min[MAXS][MAXS], int max[MAXS][MAXS], int i, int j, int v,
	   bool * updated, char *rule)
{
  if (v < 0)
    return (true);		// no change so no pb
  if (min[i][j] < v)
    {
      if (debug_schur)
	fprintf (stderr, "Improvement %s: %d,%d min:%d->%d\n", rule, i, j,
		 min[i][j], v);
      *updated = true;
      min[i][j] = v;

      if (min[i][j] > max[i][j])
	{
	  if (debug_schur)
	    fprintf (stderr, "error %s min>max for %d,%d, %d, %d\n", rule, i,
		     j, min[i][j], max[i][j]);

	  return (false);
	}

      if (debug_schur && min[i][j] == max[i][j])
	{
	  fprintf (stderr, "min=max for %d,%d=%d\n", i, j, v);
	}
    }
  return (true);
}

bool
updateMax (int min[MAXS][MAXS], int max[MAXS][MAXS], int i, int j, int v,
	   bool * updated, char *rule)
{
  if (v >= MYINFINITY)
    return (true);		// no change so no pb
  if (v < max[i][j])
    {
      if (debug_schur)
	fprintf (stderr, "Improvement %s: %d,%d max:%d->%d\n", rule, i, j,
		 max[i][j], v);
      max[i][j] = v;
      *updated = true;

      if (min[i][j] > max[i][j])
	{
	  if (debug_schur)
	    fprintf (stderr, "error %s min>max for %d,%d, %d, %d\n", rule, i,
		     j, min[i][j], max[i][j]);

	  return (false);
	}

      if (debug_schur && min[i][j] == max[i][j])
	{
	  fprintf (stderr, "min=max for %d,%d=%d\n", i, j, v);
	}
    }
  return (true);
}

bool
validHivesIJ (int min[MAXS][MAXS], int max[MAXS][MAXS], unsigned n, unsigned i,
	      unsigned j)
{
  if (min[i][j] > max[i][j] || (max[i][j] < min[i + 1][j - 1] + min[i][j + 1] - max[i + 1][j])	// R11
      || (max[i][j] < min[i][j - 1] + min[i - 1][j + 1] - max[i - 1][j])	// R12
      || (j >= 2 && min[i][j] > max[i][j - 1] + max[i + 1][j - 1] - min[i + 1][j - 2])	// R13
      || (i + j + 1 <= n && min[i][j] > max[i - 1][j + 1] + max[i][j + 1] - min[i - 1][j + 2])	//R14
      || (i + j + 2 < n && min[i][j] > max[i + 1][j] + max[i][j + 1] - min[i + 1][j + 1])	//R21
      || (min[i][j] > max[i][j - 1] + max[i - 1][j] - min[i - 1][j - 1])	//R22
      || (max[i][j] < min[i - 1][j] + min[i][j + 1] - max[i - 1][j + 1])	//R23
      || (max[i][j] < min[i][j - 1] + min[i + 1][j] - max[i + 1][j - 1])	//R24
      || (max[i][j] < min[i - 1][j + 1] + min[i + 1][j] - max[i][j + 1])	//R31
      || (max[i][j] < min[i - 1][j] + min[i + 1][j - 1] - max[i][j - 1])	//R32
      || (i < n && min[i][j] > max[i + 1][j - 1] + max[i + 1][j] - min[i + 2][j - 1])	//R33
      || (i >= 2 && min[i][j] > max[i - 1][j] + max[i - 1][j + 1] - min[i - 2][j + 1])	//R34
    )
    return (false);
  else
    return (true);

}

bool
validHives (int min[MAXS][MAXS], int max[MAXS][MAXS], unsigned n)
{
  register unsigned i, j;

  for (i = 1; i < n; i++)
    for (j = 1; i + j < n; j++)
      if (!validHivesIJ (min, max, n, i, j))
	return (false);
  return (true);

}

//returns the first point that has strict inequalities or 
// 0 for its first coordinate if all inequalities are equalities 
// or -1 for its first coordinate if there is no solution.
point_t
hives (int min[MAXS][MAXS], int max[MAXS][MAXS], unsigned n, bool * updated)
{
  register unsigned i, j;
  point_t firstOpen;
  register bool ok = true;

  *updated = false;
  firstOpen.x = firstOpen.y = 0;
  for (i = 1; i < n; i++)
    for (j = 1; i + j < n; j++)
      {
	if (min[i][j] < max[i][j])
	  {
	    // R1 :1
	    ok =
	      updateMin (min, max, i, j,
			 min[i + 1][j - 1] + min[i][j + 1] - max[i + 1][j],
			 updated, "R11") &&
	      // R1 :2
	      updateMin (min, max, i, j,
			 min[i][j - 1] + min[i - 1][j + 1] - max[i - 1][j],
			 updated, "R12") &&
	      // R1 :3
	      ((j < 2) ||
	       updateMax (min, max, i, j,
			  max[i][j - 1] + max[i + 1][j - 1] - min[i + 1][j -
									 2],
			  updated, "R13")) &&
	      // R1 :4
	      ((i + j + 1 > n) ||
	       updateMax (min, max, i, j,
			  max[i - 1][j + 1] + max[i][j + 1] - min[i - 1][j +
									 2],
			  updated, "R14")) &&
	      // R2 :1
	      ((i + j + 2 >= n) ||
	       updateMax (min, max, i, j,
			  max[i + 1][j] + max[i][j + 1] - min[i + 1][j + 1],
			  updated, "R21")) &&
	      // R2 :2
	      updateMax (min, max, i, j,
			 max[i][j - 1] + max[i - 1][j] - min[i - 1][j - 1],
			 updated, "R22") &&
	      // R2 :3
	      updateMin (min, max, i, j,
			 min[i - 1][j] + min[i][j + 1] - max[i - 1][j + 1],
			 updated, "R23") &&
	      // R2 :4
	      updateMin (min, max, i, j,
			 min[i][j - 1] + min[i + 1][j] - max[i + 1][j - 1],
			 updated, "R24") &&
	      // R3 :1
	      updateMin (min, max, i, j,
			 min[i - 1][j + 1] + min[i + 1][j] - max[i][j + 1],
			 updated, "R31") &&
	      // R3 :2
	      updateMin (min, max, i, j,
			 min[i - 1][j] + min[i + 1][j - 1] - max[i][j - 1],
			 updated, "R32") &&
	      // R3 :3
	      ((i >= n) ||
	       updateMax (min, max, i, j,
			  max[i + 1][j - 1] + max[i + 1][j] - min[i + 2][j -
									 1],
			  updated, "R33")) &&
	      // R3 :4
	      ((i < 2) ||
	       updateMax (min, max, i, j,
			  max[i - 1][j] + max[i - 1][j + 1] - min[i - 2][j +
									 1],
			  updated, "R34"));
	  }
	if (!ok)
	  {
	    if (debug_schur)
	      fprintf (stderr, "pb detected\n");
	    firstOpen.x = -1;
	    return (firstOpen);
	  }

	if (min[i][j] < max[i][j] && firstOpen.x == 0)
	  {
	    firstOpen.x = i;
	    firstOpen.y = j;
	    if (debug_schur)
	      fprintf (stderr, "firstOpen:%d,%d\n", i, j);
	  }
      }
  return (firstOpen);
}

unsigned
tryValue (point_t p, intervals inter, unsigned v, unsigned n, bool dispRes)
{
  unsigned nbSol = 0;
  int v2;
  point_t newp;
  bool updated;

  inter.min[p.x][p.y] = inter.max[p.x][p.y] = v;
  if (!validHivesIJ (inter.min, inter.max, n, p.x, p.y))
    return (0);
  newp = hives (inter.min, inter.max, n, &updated);
  if (newp.x == -1)
    return (0);
  if (newp.x == 0)
    {
      if (validHives (inter.min, inter.max, n))
	{
	  if (dispRes) displayResults (inter.min, n);
	  return (1);
	}
      else
	{
	  if (debug_schur)
	    fprintf (stderr, "INVALID HIVES DETECTED\n");
	  return (0);
	}
    }
  else if (newp.x != p.x || newp.y != p.y)
    {
      if (debug_schur)
	fprintf (stderr, "Continue on enum point: %d,%d\n", p.x, p.y);
      for (v2 = inter.min[newp.x][newp.y]; v2 <= inter.max[newp.x][newp.y];
	   v2++)
	nbSol += tryValue (newp, inter, v2, n, dispRes);
      return (nbSol);
    }
  return (0);
}

unsigned long
hivesLRcoef (frame lambda, frame mu, frame nu, unsigned short size,
	     unsigned dilatation, bool dispRes)
{
  unsigned nbSol = 0, nbOfPre;
  int v;
  point_t p;
  intervals inter;
  bool updated;

  if (debug_schur)
    {
      displayPart (&lambda, size);
      displayPart (&mu, size);
      displayPart (&nu, size);
    }

  if (!ValidSumPartitions (&lambda, &mu, &nu, size))
    {
      fprintf (stderr, "Invalid partitions (check their sum)\n");
      return (1);
    }
  else if (debug_schur)
    fprintf (stderr, "Valid partitions\n");

  for (v = 1; v <= size; v++)
    {
      lambda.A[v] *= dilatation;
      mu.A[v] *= dilatation;
      nu.A[v] *= dilatation;
    }
  initTabs (inter.min, inter.max, &lambda, &mu, &nu, size,
	    sumOfPartition (&nu, size));

  nbOfPre = 0;
  do
    {
      p = hives (inter.min, inter.max, size, &updated);
      nbOfPre++;

      if (debug_schur)
	{
	  printf ("------\nMatrix of reduced intervals:\n");
	  displayTabs (inter.min, inter.max, size);
	}

    }
  while (p.x > 0 && updated);

  if (debug_schur) fprintf (stderr, "Number of interval restrictions:%d\n", nbOfPre);

  if (p.x > 0)
    {
      if (debug_schur)
	fprintf (stderr, "------\nEnumeration starts on %d,%d (%d->%d)\n",
		 p.x, p.y, inter.min[p.x][p.y], inter.max[p.x][p.y]);
      for (v = inter.min[p.x][p.y]; v <= inter.max[p.x][p.y]; v++)
	nbSol += tryValue (p, inter, v, size, dispRes);
    }
  else if (p.x == -1) // no solution
	nbSol = 0;
  else
    nbSol = 1;
  return nbSol;
}
