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


# include <stdio.h>
# include <unistd.h>

# include "standard.h"
# include "define.h"
# include "dim.h"
# include "type.h"
# include "var.h"
# include "skew.h"
# include "s2.h"
# include "rib_to_s.h"

termptr
rib_to_s (termptr rib)
{
	unsigned k,n=rib->val.length;
	unsigned sum=0;
	term lambda,mu;

	lambda.next=NULL;lambda.mult=1;lambda.slab=' ';
	mu.next=NULL;mu.mult=1;mu.slab=' ';

	lambda.val= nolls; /* initialized to zeros */
	mu.val= nolls;


	for (k= n; k>=1 ; k--)
	{
		sum += rib->val.A[k];
		lambda.val.A[k]= sum -n+k;
	}
	lambda.val.length = n;

	for (k=1; k<n; k++)
		mu.val.A[k]=lambda.val.A[k+1] -1;
	mu.val.length = n-1;

	if (debug_schur)
	{
		fprintf(stderr,"input=");
		putsfn(&output,rib,true);
		fprintf(stderr,"lambda=");
		putsfn(&output,&lambda,true);
		fprintf(stderr,"mu=");
		putsfn(&output,&mu,true);
	}
	
	return(lskew (&lambda, &mu));
}

