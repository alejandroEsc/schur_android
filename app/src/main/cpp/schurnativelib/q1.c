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

/*
**	Definitions for i/o
*/
#include <stdio.h>
#include "standard.h"
#include "define.h"
/*
**	Start of program definitions
*/
#include "dim.h"
#include "type.h"
#include "var.h"
#include "utils.h"
#include "s1.h"
#include "m.h"
#include "q2.h"
#include "g.h"
#include "s2.h"
#include "r.h"
#include "q1.h"

termptr
jstndise (termptr * partj)
{
    register termptr R140;
    int m, d;
    register int i;
    bool test;
    termptr temp;

    snu (&temp);
    for (i = 1; i <= maxdim; i++) {
	temp->val.A[i] = 0;
    }
    temp->mult = 1;
    temp->val = (*partj)->val;
    if (temp != NULL) {
	/*register termptr W4 = &(*temp); *//*12/12/95 */

	m = qlen (temp->val);
	do {
	    i = 1;
	    test = true;
	    do {
		if (temp->val.A[i] < temp->val.A[i + 1]) {
		    test = false;
		    d = temp->val.A[i];
		    temp->val.A[i] = temp->val.A[i + 1];
		    temp->val.A[i + 1] = d;
		}
		i = i + 1;
	    }
	    while (!((i == m)));
	}
	while (!((test == true)));
    }
    R140 = temp;
    return R140;
}

void
nqup (termptr * list)
{
    termptr newlist, temp;
    bool startx;
    int dummy;
    register int i;

    newlist = NULL;
    while ((*list) != NULL) {
	register termptr W5 = &(*(*list));
	do {
	    startx = true;
	    for (i = 1; i <= maxdim - 1; i++) {
		if (abs (W5->val.A[i]) < abs (W5->val.A[i + 1])) {
		    startx = false;
		    dummy = W5->val.A[i];
		    W5->val.A[i] = W5->val.A[i + 1];
		    W5->val.A[i + 1] = dummy;
		}
	    }
	}
	while (!((startx == true)));
	temp = W5->next;
	W5->next = newlist;
	newlist = (*list);
	(*list) = temp;
    }
    (*list) = newlist;
}

termptr
joinmn (termptr * mu, termptr * nu)
{
    termptr dummy, temp;
    int m, n;
    register int i;
    snu (&temp);
    m = len (&(*mu)->val);
    n = len (&(*nu)->val);
    temp->mult = 1;
    for (i = 1; i <= maxdim; i++) {
	temp->val.A[i] = 0;
    }
    for (i = 1; i <= m; i++) {
	temp->val.A[i] = (*mu)->val.A[i];
    }
    for (i = 1; i <= n; i++) {
	temp->val.A[i + m] = (*nu)->val.A[i];
    }
    nqup (&temp);
    dummy = jmod (&temp);
    dispsfn (&temp);
    return dummy;
}

void
jnum (termptr * list, int *p, int *q, int *r)
{
    int i, l, s;
    bool test;

    i = (*p);
    s = (*p);
    test = true;
    l = 0;
    do {
	if (((*list)->val.A[s + 1] + 1) >= (*list)->val.A[s])
	    s = s + 1;
	else
	    test = false;
    }
    while (!(((*list)->val.A[s] == 0) || (test == false)));
    if (s > i) {
	do {
	    if ((*list)->val.A[i] == (*list)->val.A[i + 1]) {
		l = l + 1;
		i = i + 2;
	    } else
		i = i + 1;
	} while (!(i >= s));
    }
    (*q) = s;
    (*r) = l;
}

termptr
jmod (termptr * tterm)
{
    int i, j, k;
    bool test;
    termptr temp;

    snu (&temp);
    (*temp) = (*(*tterm));
    {
	register termptr W14 = &(*temp);

	do {
	    test = true;
	    i = 1;
	    do {
		jnum (&temp, &i, &j, &k);
		if (k != 0) {
		    test = false;
		    W14->val.A[i] = W14->val.A[i] + k;
		    W14->val.A[i + 1] = W14->val.A[i + 1] - k;
		}
		i = i + 1;
		nqup (&temp);
	    } while (!((temp->val.A[i] == 0)));
	} while (!((test == true)));
    }
    return temp;
}

termptr
rmod (termptr * tterm)
{
    int i, j, k;
    bool test;
    termptr temp;

    snu (&temp);
    (*temp) = (*(*tterm));
    do {
	test = true;
	i = 1;
	do {
	    k = 0;
	    if (temp->val.A[i] == temp->val.A[i + 1]) {
		j = i;
		test = false;
		do {
		    j = j + 1;
		    k = k + 1;
		} while (!(temp->val.A[j] != temp->val.A[j + 1]));
		temp->val.A[i] = temp->val.A[i] + k;
		temp->val.A[i + 1] = temp->val.A[i + 1] - k;
	    } else
		i = i + 1;
	    stndise (&temp);
	    temp->mult = 1;
	} while (!((test == false) || (temp->val.A[i] == 0)));
    } while (!((test == true)));
    return temp;
}

int
repmult (termptr * lambda, termptr * mu, termptr * nu)
{
    register int R144;
    int a, b, c, f, i, j, k, m, p, q, r, x, y, z;
    register int n;

    i = qlen ((*lambda)->val);	/*qqq */
    j = qlen ((*mu)->val);
    k = qlen ((*nu)->val);
    p = wtfrm (&(*lambda)->val);
    q = wtfrm (&(*mu)->val);
    r = wtfrm (&(*nu)->val);
    a = p % 2;
    b = q % 2;
    c = r % 2;
    x = (i - a) / 2;
    y = (j - b) / 2;
    z = (k - c) / 2;
    m = abs (y + z - x);
    f = 1;
    for (n = 1; n <= m; n++) {
	f = f * 2;
    }
    R144 = f;
    return R144;
}

bool
deadtest (frame * lambda, termptr * mu, termptr * nu)
{
    frame u, v;
    int j, k, m, n, p, q, r, s, t, x;
    register int i;
    bool test;

    test = true;
    for (i = 1; i <= maxdim; i++) {
	u.A[i] = 0;
    }
    for (i = 1; i <= maxdim; i++) {
	v.A[i] = 0;
    }
    u = (*mu)->val;
    v = (*nu)->val;
    j = len (&u);
    k = len (&v);		/*qqq */
    r = len (lambda);
    if (r > j)
	s = j;
    else
	s = j - 1;
    q = 0;
    for (i = 1; i <= s; i++) {
	p = u.A[i] - u.A[i + 1] - 1;
	q = q + p;
    }
    r = 0;
    for (i = 1; i <= k - 1; i++) {
	t = v.A[i] - v.A[i + 1] - 1;
	r = r + t;
    }
    m = 0;
    for (i = 1; i <= r - 1; i++) {
	n = u.A[i] + v.A[i] - lambda->A[i];
	m = m + n;
    }
    x = q + r;
    if (x < m)
	test = false;
    return test;
}

termptr
qsameweight (termptr * sfn1, termptr * sfn2)
{
    int left, edge1, edge2, sum, prev, a, b, k, l, m, n, p, x, y, z;
    register int i;
    frame rho;
    bool test ;
    termptr result, ptr1, lastrho;
    prev = 0;
    a = len (&(*sfn1)->val);
    b = len (&(*sfn2)->val);
    l = a + b;
    y = wtfrm (&(*sfn1)->val);
    z = wtfrm (&(*sfn2)->val);
    if (redu) {
	m = y - z;
	for (i = 0; i <= maxdim; i++) {
	    rho.A[i] = 0;
	}
	rho.A[1] = m;
	sum = 0;
	edge2 = 1;
    } else {
	m = y + z;
	for (i = 0; i <= maxdim; i++) {
	    rho.A[i] = (*sfn1)->val.A[i] + (*sfn2)->val.A[i];
	}

	x = len (&rho);
	lastrho = joinmn (&(*sfn1), &(*sfn2));
	/*putsfn(&output,lastrho,true); */
	if (rho.A[x] != 1)
	    edge2 = x;
	else
	    edge2 = x - 1;
	if (rho.A[x] != 1)
	    sum = 0;
	else
	    sum = 1;
    }
    i = 0;
    k = 0;
    do {
	k = k + 1;
	i = i + k;
    } while (!(i > m));

    n = k - 1;
    if (redu)
	p = n;
    else
	p = MIN (l, n);
    snu (&result);
    result->val = rho;
    result->mult = 1;
    ptr1 = result;

    while ((rho.A[1] >= p) && (prev >= 0)) {
	edge1 = edge2;
	edge2 = -1;
	left = sum + rho.A[edge1];
	prev = rho.A[edge1] - 1;

	do {
	    prev = MIN (left, prev);
	    if ((edge2 < 0) && (prev <= 1)) {
		edge2 = edge1 - 1;
		sum = left;
	    }
	    rho.A[edge1] = prev;
	    left = left - prev;
	    edge1 = edge1 + 1;
	} while (!(prev <= 0));
	x = len (&rho);
	test = true;
	i = 0;
	do {
	    i = i + 1;
	    if (rho.A[i] == rho.A[i + 1])
		test = false;
	} while (!((i == x - 1) || (test == false)));

	if ((test == true) && (!redu)) {
	    test = dead (rho, (*sfn1), (*sfn2));
	}
	/*if ((test == true) && (!redu)) {
	   s = 0;
	   t = 0;
	   f = 0;
	   g = 0;
	   do {
	   s = s + 1;
	   t = t + rho.A[s];
	   f = f + (*sfn1)->val.A[s] + (*sfn2)->val.A[s];
	   g = g + lastrho->val.A[s];
	   if ((t > f) || (t < g))
	   test = false;
	   } while (!(((test == false) || (s == x))));
	   } */
	if (test == true) {
	    snu (&ptr1->next);
	    ptr1 = ptr1->next;
	    ptr1->val = rho;
	    ptr1->mult = 1;
	    ptr1->val.A[edge1 - 1] = 0;
	    for (i = edge1; i <= maxdim; i++) {
		ptr1->val.A[i] = 0;
	    }
	}
	ptr1->next = NULL;
    }
    if (!redu)
	dispsfn (&lastrho);
    return result;
}

termptr
qparts (termptr * sfn1)
{
    int left, edge1, edge2, sum, prev, a, c = 0, g, s, f, t, x, y, z;	// c was not initialized FB
    register int i;
    frame rho, v;
    bool test, test1;
    termptr result, ptr1, lastrho;

    prev = 0;
    a = qlen ((*sfn1)->val);	/*q */
    y = wtfrm (&(*sfn1)->val);
    lastrho = rmod (&(*sfn1));
    for (i = 1; i <= maxdim; i++) {
	v.A[i] = 0;
    }
    for (i = 1; i <= maxdim; i++) {
	rho.A[i] = 0;
    }
    test1 = true;
    z = y;

    i = 0;
    do {
	i = i + 1;
	v.A[i] = (*sfn1)->val.A[i] + a - 2 * (i - 1) - 1;
	if (v.A[i] >= z) {
	    rho.A[i] = z;
	    test1 = false;
	} else {
	    rho.A[i] = v.A[i];
	    z = z - v.A[i];
	}
    } while (!((test1 == false) || (i == a)));

    for (i = 1; i <= maxdim; i++) {
	v.A[i] = rho.A[i];
    }
    x = qlen (rho);		/*q */
    if (rho.A[x] != 1)
	edge2 = x;
    else
	edge2 = x - 1;
    if (rho.A[x] != 1)
	sum = 0;
    else
	sum = 1;
    snu (&result);
    result->val = rho;
    result->mult = rcoeff (&(*sfn1), &rho);
    ptr1 = result;
    while ((rho.A[1] >= lastrho->val.A[1]) && (prev >= 0)) {
	edge1 = edge2;
	edge2 = -1;
	left = sum + rho.A[edge1];
	prev = rho.A[edge1] - 1;
	do {
	    prev = MIN (left, prev);
	    if ((edge2 < 0) && (prev <= 1)) {
		edge2 = edge1 - 1;
		sum = left;
	    }
	    rho.A[edge1] = prev;
	    left = left - prev;
	    edge1 = edge1 + 1;
	} while (!(prev <= 0));

	x = qlen (rho);		/*q */
	if ((x <= a)) {
	    test = true;
	    i = 0;
	    do {
		i = i + 1;
		if (rho.A[i] == rho.A[i + 1])
		    test = false;
	    } while (!((i == x - 1) || (test == false)));

	    if (test == true) {
		s = 0;
		t = 0;
		f = 0;
		g = 0;
		do {
		    s = s + 1;
		    t = t + rho.A[s];
		    f = f + lastrho->val.A[s];
		    g = g + v.A[s];
		    if ((t < f) || (t > g))
			test = false;
		} while (!(((test == false) || (s == x))));
	    }
	    if (test == true)
		c = rcoeff (&(*sfn1), &rho);
	    if ((test == true) && (c != 0)) {
		snu (&ptr1->next);
		ptr1 = ptr1->next;
		ptr1->val = rho;
		ptr1->mult = rcoeff (&(*sfn1), &rho);
		ptr1->val.A[edge1 - 1] = 0;
		for (i = edge1; i <= maxdim; i++) {
		    ptr1->val.A[i] = 0;
		}
	    }
	}
	ptr1->next = NULL;
    }
    return result;
}

int
rest (frame * list, int i)
{
    register int R148;
    int k, l;
    register int j;

    k = 0;
    l = qlen ((*list));		/*q */
    for (j = 1; j <= l; j++) {
	if (((list->A[j] == i) || (list->A[j] == -i)))
	    k = k + 1;
    }
    R148 = k;
    return R148;
}

int
qmaxin (frame * list)
{
    register int R149;
    int i, j;

    i = 0;
    j = 0;
    while (list->A[j + 1] != 0) {
	if (abs (list->A[j + 1]) > i)
	    i = abs (list->A[j + 1]);
	j = j + 1;
    }
    R149 = i;
    return R149;
}

void
maker (termptr * list, frame * v)
{
    int l;
    register int k;
    register int i;
    for (i = 1; i <= maxdim; i++) {
	v->A[i] = 0;
    }
    l = qlen ((*list)->val);	/*q */
    (*v) = (*list)->val;
    for (k = 1; k <= l; k++) {
	v->A[l + k] = -(*list)->val.A[l - k + 1];
    }
}

void
makeray (frame * s, qframe * ma, int *wide, int *size)
{
    int m, n;
    register int l;
    register int k;
    register int j;
    register int i;
    {
	int lastStep51 = (*wide);
	for (i = 1; i <= lastStep51; i++) {
	    for (j = 1; j <= *size; j++) {
		ma->A[i - 1].A[j] = 0;
	    }
	}
    }
    {
	int lastStep55 = (*wide);
	for (k = 1; k <= lastStep55; k++) {
	    m = 0;
	    n = 0;
	    for (l = 1; l <= *size; l++) {
		n = n + 1;
		if ((s->A[n] == k))
		    m = m + 1;
		ma->A[k - 1].A[n] = m;
	    }
	}
    }
}

bool
latticetest (termptr * list)
{
    register bool R150;
    qframe mr;
    frame r;
    int l, m, n;
    register int j;
    register int i;
    bool test;

    l = qlen ((*list)->val);	/*q */
    test = true;
    n = 2 * l;
    maker (&(*list), &r);
    m = qmaxin (&r);
    makeray (&r, &mr, &m, &n);
    for (i = 2; i <= m; i++) {
	for (j = 1; j <= n - 1; j++) {
	    if ((mr.A[i - 1].A[j] == mr.A[i - 1 - 1].A[j])) {
		if (((j < l) && ((r.A[j + 1] == i) || (r.A[j + 1] == -i))))
		    test = false;
		/*if ((r.A[j + 1] == -i))
		   test = false; */
		if (((j < n) && (j >= l))
		    && ((r.A[j + 1] == i) || (r.A[j + 1] == -(i - 1))))
		    test = false;
	    }
	}
    }
    R150 = test;
    return R150;
}

void
salam (termptr * list, int i)
{
    termptr dummy, list2, temp2, temp3;
    bool frst;

    list2 = (*list);
    frst = true;
    temp3 = NULL;
    while ((*list) != NULL) {
	snu (&temp2);
	(*temp2) = (*(*list));
	temp2->next = NULL;
	temp2->val.A[i] = -temp2->val.A[i];
	dummy = ladd (temp2, temp3);
	if (frst == false)
	    ldisp (&temp3);
	temp3 = dummy;
	ldisp (&temp2);
	frst = false;
	(*list) = (*list)->next;
    }
    (*list) = ladd (list2, temp3);
    ldisp (&list2);
    ldisp (&temp3);
}

void
salam1 (termptr * temp, frame * v, int b, int i, int m, int q,
	int x, bool * test)
{
    termptr list2, dummy;
    frame rho;

    (*temp)->val.A[i] = b;
    (*test) = false;
    rho = (*temp)->val;
    if (((b != (*temp)->val.A[i - 1]) && (rest (&rho, x) != v->A[x]))) {
	snu (&list2);
	dummy = (*temp);
	(*list2) = (*(*temp));
	list2->val.A[i - 1] = -list2->val.A[i - 1];
	(*temp) = ladd ((*temp), list2);
	ldisp (&list2);
	ldisp (&dummy);
    }
    if (((i == q) && ((m + 1) != v->A[b])))
	salam (&(*temp), i);
}

void
second (termptr * list1, termptr * nu, int d, int r, int s,
	int i, int c, int n, int q)
{
    int m, l, x, t, z;
    register int b;
    register int k;
    frame rho, v;
    termptr temp, list2, temp3, newlist, dummy;
    bool test;

    l = qlen ((*nu)->val);	/*q */
    z = MIN (l, d);
    newlist = NULL;
    progress ();
    v = (*nu)->val;
    temp3 = (*list1);
    while ((*list1) != NULL) {
	for (k = 1; k <= maxdim; k++) {
	    rho.A[k] = 0;
	}
	rho = (*list1)->val;
	t = qmaxin (&rho);
	for (b = 1; b <= z; b++) {
	    m = rest (&rho, b);
	    x = abs (rho.A[i - 1]);
	    test = true;
	    if (((m < v.A[b]) && ((t + 1) >= b))) {
		if (((i == c) && (i > s))) {
		    snu (&temp);
		    temp->val = rho;
		    temp->mult = 1;
		    temp->val.A[i] = b;
		    test = false;
		    if ((((m + 1) != v.A[b]) && (i == q))) {
			snu (&list2);
			(*list2) = (*temp);
			list2->val.A[i] = -b;
			dummy = ladd (temp, list2);
			ldisp (&list2);
			ldisp (&temp);
			temp = dummy;
		    }
		} else if (((i == c) && (i <= s))) {
		    if ((((rho.A[i - n + r - 1] < 0)
			  && (b >= abs (rho.A[i - n + r - 1]))))
			|| (((rho.A[i - n + r - 1] > 0)
			     && (b > rho.A[i - n + r - 1])))) {
			snu (&temp);
			temp->val = rho;
			temp->mult = 1;
			temp->val.A[i] = b;
			test = false;
			if (((c == q) && ((m + 1) != v.A[b]))) {
			    snu (&list2);
			    (*list2) = (*temp);
			    list2->val.A[i] = -list2->val.A[i];
			    dummy = ladd (list2, temp);
			    ldisp (&temp);
			    ldisp (&list2);
			    temp = dummy;
			}
		    }
		} else if (((i != c) && (i <= s))) {
		    if ((rho.A[i - n + r - 1] < 0)) {
			if (((b >= abs (rho.A[i - n + r - 1]))
			     && (b <= rho.A[i - 1]))) {
			    snu (&temp);
			    temp->val = rho;
			    temp->mult = 1;
			    salam1 (&temp, &v, b, i, m, q, x, &test);
			}
		    } else if ((rho.A[i - n + r - 1] > 0)) {
			if (((b > rho.A[i - n + r - 1])
			     && (b <= rho.A[i - 1]))) {
			    snu (&temp);
			    temp->val = rho;
			    temp->mult = 1;
			    salam1 (&temp, &v, b, i, m, q, x, &test);
			}
		    }
		} else if (((i != c) && (i > s))) {
		    if ((b <= rho.A[i - 1])) {
			snu (&temp);
			temp->val = rho;
			temp->mult = 1;
			salam1 (&temp, &v, b, i, m, q, x, &test);
		    }
		}
	    }
	    if (!test) {
		dummy = ladd (temp, newlist);
		ldisp (&newlist);
		newlist = dummy;
		ldisp (&temp);
	    }
	}
	(*list1) = (*list1)->next;
    }
    ldisp (&temp3);
    (*list1) = newlist;
}

void
rsecond (termptr * list1, frame * nu, int d, int i, int c, int q, int s)
{
    int m, l, x, t, z;
    register int b;
    register int k;
    frame rho, v;
    termptr temp, list2, temp3, newlist, dummy;
    bool test;

    l = qlen ((*nu));		/*q */
    z = MIN (l, d);
    newlist = NULL;
    progress ();
    v = (*nu);
    temp3 = (*list1);
    while ((*list1) != NULL) {
	for (k = 1; k <= maxdim; k++) {
	    rho.A[k] = 0;
	}
	rho = (*list1)->val;
	t = qmaxin (&rho);
	for (b = 1; b <= z; b++) {
	    m = rest (&rho, b);
	    x = abs (rho.A[i - 1]);
	    test = true;
	    if (((m < v.A[b]) && ((t + 1) >= b))) {
		if (i == c) {
		    if (((rho.A[i - s] < 0) && (b >= abs (rho.A[i - s])))
			|| ((rho.A[i - s] > 0) && (b > rho.A[i - s]))) {
			snu (&temp);
			temp->val = rho;
			temp->mult = 1;
			temp->val.A[i] = b;
			test = false;
			if (((c == q) && ((m + 1) != v.A[b]))) {
			    snu (&list2);
			    (*list2) = (*temp);
			    list2->val.A[i] = -list2->val.A[i];
			    dummy = ladd (list2, temp);
			    ldisp (&temp);
			    ldisp (&list2);
			    temp = dummy;
			}
		    }
		} else if (i != c) {
		    if ((rho.A[i - s] < 0)) {
			if (((b >= abs (rho.A[i - s]))
			     && (b <= rho.A[i - 1]))) {
			    snu (&temp);
			    temp->val = rho;
			    temp->mult = 1;
			    salam1 (&temp, &v, b, i, m, q, x, &test);
			}
		    } else if ((rho.A[i - s] > 0)) {
			if (((b > rho.A[i - s]) && (b <= rho.A[i - 1]))) {
			    snu (&temp);
			    temp->val = rho;
			    temp->mult = 1;
			    salam1 (&temp, &v, b, i, m, q, x, &test);
			}
		    }
		}
	    }
	    if (!test) {
		dummy = ladd (temp, newlist);
		ldisp (&newlist);
		newlist = dummy;
		ldisp (&temp);
	    }
	}
	(*list1) = (*list1)->next;
    }
    ldisp (&temp3);
    (*list1) = newlist;
}

void
first (termptr * list, termptr * nu, int d, int r, int c, int s, int q, int n)
{
    register int i;
    for (i = c; i <= q; i++) {
	{
	    /*register termptr W72 = &(*(*list)); *//*12/12/95 */

	    second (&(*list), &(*nu), d, r, s, i, c, n, q);
	}
    }
}

void
rfirst (termptr * list, frame * nu, int d, int c, int q, int s)
{
    register int i;
    for (i = c; i <= q; i++) {
	{
	    /*register termptr W75 = &(*(*list)); *//*12/12/95 */

	    rsecond (&(*list), &(*nu), d, i, c, q, s);
	}
    }
}

termptr
wgenerate (termptr lambda, termptr mu, termptr nu)
{
    register termptr R151;
    termptr dummy, list1, list, temp, temp2, temp1;
    frame v;
    int q, r, s, c, d, l, n;
    register int k;
    register int p;
    bool frst;


    l = len (&lambda->val);
    snu (&temp);
    for (p = 1; p <= maxdim; p++) {
	v.A[p] = 0;
    }
    for (p = 1; p <= l; p++) {
	v.A[p] = lambda->val.A[p] - mu->val.A[p];
    }
    for (k = 1; k <= maxdim; k++) {
	temp->val.A[k] = 0;
    }
    if (lambda->val.A[1] >= mu->val.A[1]) {
	{
	    int lastStep83 = v.A[1];
	    for (k = 1; k <= lastStep83; k++) {
		temp->val.A[k] = 1;
	    }
	}
	temp->mult = 1;
	if (lambda->val.A[1] != mu->val.A[1] + nu->val.A[1]) {
	    snu (&temp1);
	    (*temp1) = (*temp);
	    temp1->val.A[v.A[1]] = -1;
	    dummy = ladd (temp, temp1);
	    ldisp (&temp1);
	    ldisp (&temp);
	    temp = dummy;
	}
	d = 1;
	c = v.A[1];
	if (l > 1) {
	    do {
		d = d + 1;
		list1 = NULL;
		r = (lambda->val.A[d - 1] - lambda->val.A[d]);
		c = c + 1;
		q = (c + v.A[d] - 1);
		n = lambda->val.A[d - 1] - mu->val.A[d - 1];
		s = (c + lambda->val.A[d] - mu->val.A[d - 1]);
		if ((q >= c) && (lambda->val.A[d] != mu->val.A[d])) {
		    frst = true;
		    list1 = NULL;
		    temp2 = temp;
		    while (temp != NULL) {
			/*register termptr W84 = &(*temp); *//*12/12/95 */

			snu (&list);
			(*list) = (*temp);
			list->next = NULL;
			{
			    first (&list, &nu, d, r, c, s, q, n);
			    dummy = ladd (list1, list);
			    ldisp (&list);
			    if (frst == false)
				ldisp (&list1);
			    list1 = dummy;
			    frst = false;
			}
			temp = temp->next;
		    }
		    ldisp (&temp2);
		    temp = list1;
		    c = (c + v.A[d] - 1);
		}
	    }
	    while (!(d == l));
	}
    } else
	temp = NULL;
    R151 = temp;
    return R151;
}

termptr
rgenerate (termptr * lambda, frame * mu)
{
    register termptr R152;
    termptr dummy, list1, list, temp, temp2, temp1;
    int q, r, s, c, d, l;
    register int k;
    bool frst;

    l = qlen ((*lambda)->val);	/*q */
    snu (&temp);
    for (k = 1; k <= maxdim; k++) {
	temp->val.A[k] = 0;
    }
    r = (*lambda)->val.A[1];
    for (k = 1; k <= r; k++) {
	temp->val.A[k] = 1;
    }
    temp->mult = 1;
    if ((*lambda)->val.A[1] != mu->A[1]) {
	snu (&temp1);
	(*temp1) = (*temp);
	temp1->val.A[r] = -1;
	dummy = ladd (temp, temp1);
	ldisp (&temp1);
	ldisp (&temp);
	temp = dummy;
    }
    d = 1;
    c = r;
    if (l > 1) {
	do {
	    d = d + 1;
	    list1 = NULL;
	    c = c + 1;
	    q = (c + (*lambda)->val.A[d] - 1);
	    s = (*lambda)->val.A[d];
	    frst = true;
	    list1 = NULL;
	    temp2 = temp;
	    while (temp != NULL) {
		/*register termptr W89 = &(*temp); *//*12/12/95 */

		snu (&list);
		(*list) = (*temp);
		list->next = NULL;
		{
		    rfirst (&list, &(*mu), d, c, q, s);
		    dummy = ladd (list1, list);
		    ldisp (&list);
		    if (frst == false)
			ldisp (&list1);
		    list1 = dummy;
		    frst = false;
		}
		temp = temp->next;
	    }
	    ldisp (&temp2);
	    temp = list1;
	    c = q;
	}
	while (!(d == l));
    }
    R152 = temp;
    return R152;
}

int
coeff (termptr lambda, termptr mu, termptr nu)
{
    register int R153;
    termptr list, list1;
    int i;
    bool test;

    i = 0;
    /*putsfn(&output,lambda,true); */
    list = wgenerate (lambda, mu, nu);
    /*putsfn(&output,list,true);Getl(&input); */
    list1 = list;
    while (list != NULL) {
	test = true;
	test = latticetest (&list);
	if (test)
	    i = i + 1;
	list = list->next;
    }
    ldisp (&list1);

    R153 = i;
    return R153;
}

int
rcoeff (termptr * lambda, frame * mu)
{
    register int R154;
    termptr list, list1;
    int i;
    bool test;

    i = 0;
    list = rgenerate (&(*lambda), &(*mu));
    list1 = list;
    while (list != NULL) {
	test = true;
	test = latticetest (&list);
	if (test == true)
	    i = i + 1;
	list = list->next;
    }
    ldisp (&list1);
    R154 = i;
    return R154;
}

int
qfactor (termptr lambda, termptr mu, termptr nu)
{
    int i, j, k, l, n;
    register int m;

    i = len (&mu->val);
    j = len (&nu->val);
    k = len (&lambda->val);
    if (rep == true)
	l = (i + j - k) / 2;
    else
	l = i + j - k;
    if (l != 0) {
	n = 1;
	for (m = 1; m <= l; m++) {
	    n = 2 * n;
	}
    } else
	n = 1;
    return n;
}

bool
platticetest (termptr list)
{
    qframe mr;
    frame temp;
    int l, m, n;
    register int j;
    register int i;
    bool test;

    temp = list->val;
    test = true;
    if (temp.A[1] != 1)
	test = false;
    else {
	l = qlen (temp);	/*q */
	n = 2 * l;
	m = qmaxin (&temp);
	makeray (&temp, &mr, &m, &n);
	for (i = 2; i <= m; i++)
	    for (j = 1; j <= n - 1; j++)
		if ((mr.A[i - 1].A[j] == mr.A[i - 1 - 1].A[j]))
		    if ((((j < l)
			  && ((temp.A[j + 1] == i) || (temp.A[j + 1] == -i)))
			 || ((((j < n) && (j >= l))
			      && ((temp.A[j + 1] == i)
				  || (temp.A[j + 1] == -(i - 1)))))))
			test = false;
    }
    return test;
}

bool
paritysequence (termptr list)
{
    int i, l, dummy;
    bool signp, test;
    frame temp;

    temp = list->val;
    signp = true;		// corrected by FB, was sign=
    l = qlen (temp);		/*q */
    do {
	i = 1;
	test = true;
	do {
	    if (temp.A[i] > temp.A[i + 1]) {
		test = false;
		dummy = temp.A[i];
		temp.A[i] = temp.A[i + 1];
		temp.A[i + 1] = dummy;
		signp = (bool) (!signp);
	    }
	    i = i + 1;
	}
	while (!((i == l)));
    }
    while (test != true);
    return signp;
}

termptr
indexx (termptr list)
{
    termptr temp;
    int l, ind, x;
    register int i;
    register int k;

    snu (&temp);
    temp->mult = 1;
    temp->val = nolls;
    l = qlen (list->val);	/*q */
    for (k = 1; k <= l; k++) {
	{
	    register termptr W98 = &(*list);

	    x = W98->val.A[k];
	    ind = 0;
	    for (i = 1; i <= k; i++) {
		if ((W98->val.A[i] == x - 1))
		    ind = ind - 1;
		if ((W98->val.A[i] == x))
		    ind = ind + 1;
	    }
	    temp->val.A[k] = ind;
	}
    }
    return temp;
}

termptr
s_to_q (termptr list)
{
    return (termptr) craise (list, 1);
}

termptr
q_to_s (termptr list)
{
    return (termptr) lraise (list, -1);
}

bool
dead (frame rho, termptr p1, termptr p2)
     /*This implements Theorem 2 JMP 31, 1310 (1990). Returns TRUE if
        rho is a "live" partition in the Q-function product Q_p1 . Q_p2. NB Some
        rho will yield TRUE even though they do not appear in the product. */
{
    register bool R160;
    int i, t1, t2, l;
    bool test;

    test = true;
    l = len (&rho);
    i = 1;
    t1 = 0;
    t2 = 0;
    do {
	{
	    t1 = t1 + rho.A[i];
	    t2 = t2 + p1->val.A[i] + p2->val.A[i];
	    if (t2 >= t1)
		test = true;
	    else
		test = false;
	    i = i + 1;
	}
    }
    while (!(((i == l + 1) || (test == false))));
    R160 = test;
    return R160;
}
