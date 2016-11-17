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
#include "s2.h"
#include "m.h"
#include "q1.h"
#include "g.h"
#include "skew.h"
#include "q2.h"

termptr
qouterskew (termptr * mu, termptr * nu)
{
    register termptr R197;
    termptr list1, list, temp, newlist;
    int molt, i, j, p;
    newlist = NULL;
    qstndise (&(*mu));
    qstndise (&(*nu));		/*7/6/99 */
    if ((*nu)->val.A[1] != 0)
	list = qsameweight (&(*mu), &(*nu));
    else
	list = (*mu);
    molt = (*nu)->mult * (*mu)->mult;
    list1 = list;
    while (list != NULL) {
	snu (&temp);
	(*temp) = (*list);
	temp->next = NULL;
	if (redu) {
	    i = coeff ((*mu), (*nu), temp);	/*??? */
	    temp->mult = i * molt;
	} else {
	    if (((*nu)->val.A[1] == 0))
		i = 1;
	    else
		i = coeff (temp, (*mu), (*nu));	/*??? */

	    if (i != 0) {
		j = qfactor (temp, (*mu), (*nu));
		/*(void)fprintf(output.fp, " j =%10d", j), Putl(output, 0); */
		if (pfn == false)
		    i = i * j;
	    }
	    p = i * molt;
	    temp->mult = p;
	    /*(void)fprintf(output.fp, " mult =%10d", p), Putl(output, 0);
	       Getl(&input); */
	}
	add (&newlist, &temp);
	sort (&newlist, true);
	ldisp (&temp);
	list = list->next;
    }
    if (((*nu)->val.A[1] != 0))
	ldisp (&list1);
    R197 = newlist;
    return R197;
}

termptr
lqouter (termptr * list1, termptr * list2)
{
    register termptr R198;
    termptr templist, sublist, prodlist, temp, plist;
    prodlist = NULL;
    if (((*list1) != NULL) && ((*list2) != NULL)) {
	while ((*list1) != NULL) {
	    snu (&plist);
	    (*plist) = (*(*list1));
	    plist->next = NULL;
	    temp = (*list2);
	    while (temp != NULL) {
		snu (&templist);
		(*templist) = (*temp);
		templist->next = NULL;
		if (redu && !skewcompatable (plist->val, templist->val))
		    sublist = NULL;
		else {
		    if (!redu) {
			if (wtfrm (&templist->val) < wtfrm (&plist->val))	/*> changed to < */
			    sublist = qouterskew (&templist, &plist);
			else
			    sublist = qouterskew (&plist, &templist);
		    } else if (redu)
			sublist = qouterskew (&plist, &templist);
		    add (&prodlist, &sublist);
		    sort (&prodlist, true);
		}
		ldisp (&templist);
		ldisp (&sublist);
		temp = temp->next;
	    }
	    ldisp (&plist);
	    (*list1) = (*list1)->next;
	}
	ldisp (&temp);
    }
    R198 = prodlist;
    return R198;
}

termptr
reducedq (termptr list)
{
    termptr newlist, temp = NULL;
    int length;
    register int i;

    newlist = NULL;
    while (list != NULL) {
	if (newlist == NULL) {
	    snu (&newlist);
	    temp = newlist;
	} else {
	    snu (&temp->next);
	    temp = temp->next;
	}
	*temp = *list;
	length = len (&temp->val);
	if (temp->val.A[1] != 0) {
	    for (i = 1; i <= length - 1; i++) {
		temp->val.A[i] = temp->val.A[i + 1];
	    }
	    temp->val.A[length] = 0;
	}
	list = list->next;
    }
    sort (&newlist, false);
    return newlist;
}

termptr
qinner (termptr qfn1, termptr qfn2)
{
    int a, b, c, e, f, j, k, l, r, t, d;
    register int i;
    termptr temp, dum, temp1, temp2, temp3, list, list1, list3, newlist,
	dummy;
    frame v;

    c = wtfrm (&qfn1->val);
    d = wtfrm (&qfn2->val);
    newlist = NULL;
    if ((c - d) == 0) {
	a = len (&qfn1->val);
	b = len (&qfn2->val);
	e = c - a;
	f = c - b;
	if (((bool) (e & 1) && (bool) (f & 1)))
	    r = (a + b - 2) / 2;
	else if ((((bool) (e & 1) && !(bool) (f & 1))
		  || ((bool) (f & 1) && !(bool) (e & 1))))
	    r = (a + b - 1) / 2;
	else
	    r = (a + b) / 2;
	i = 0;
	j = 1;
	if (r != 0) {
	    do {
		i = i + 1;
		j = j * 2;
	    }
	    while (!(i == r));
	}
	i = 0;
	k = 1;
	do {
	    i = i + 1;
	    k = k * 2;
	}
	while (!(i == b));
	t = k / j;
	for (i = 1; i <= maxdim; i++) {
	    v.A[i] = 0;
	}
	v = qfn2->val;
	dummy = NULL;
	list = leqwt (qfn1);
	list3 = list;
	list1 = lraise (qfn1, -1);
	temp3 = list1;
	while (list != NULL) {
	    snu (&temp);
	    (*temp) = (*list);
	    temp->next = NULL;
	    temp1 = linner (temp, list1);
	    l = 0;
	    dum = temp1;
	    while (temp1 != NULL) {
		snu (&temp2);
		(*temp2) = (*temp1);
		temp2->next = NULL;
		if (temp2->val.A[1] <= v.A[1])
		    l = l + rcoeff (&temp2, &v) * temp2->mult;
		ldisp (&temp2);
		temp1 = temp1->next;
	    }
	    ldisp (&dum);
	    if (l != 0) {
		temp->mult = l * t * qfn2->mult;
		dummy = ladd (newlist, temp);
		if (newlist != NULL)
		    ldisp (&newlist);
		newlist = dummy;
	    }
	    ldisp (&temp);
	    list = list->next;
	}
	ldisp (&list3);
	ldisp (&temp3);
    }
    return newlist;
}

termptr
lqinner (termptr list1, termptr list2)
{
    termptr temp, list, templist, sublist, prodlist;

    prodlist = NULL;
    while (list1 != NULL) {
	snu (&list);
	(*list) = (*list1);
	list->next = NULL;
	templist = list2;
	while (templist != NULL) {
	    snu (&temp);
	    (*temp) = (*templist);
	    temp->next = NULL;
	    if (len (&list->val) <= len (&temp->val))
		sublist = qinner (list, temp);
	    else
		sublist = qinner (temp, list);
	    add (&prodlist, &sublist);
	    ldisp (&temp);
	    ldisp (&sublist);
	    templist = templist->next;
	}
	ldisp (&list);
	list1 = list1->next;
    }
    sort (&prodlist, true);
    return prodlist;
}

ocharptr
qspin (ocharptr list)
{
    ocharptr newlist, temp = NULL;
    int i;
    newlist = NULL;
    while (list != NULL) {
	if (newlist == NULL) {
	    cnu (&newlist);
	    temp = newlist;
	} else {
	    cnu (&temp->next);
	    temp = temp->next;
	}
	(*temp) = (*list);
	temp->spin = true;
	i = 1;
	while (temp->val.A[i] != 0) {
	    register ocharptr W6 = &(*temp);

	    W6->val.A[i] = W6->val.A[i + 1];
	    i = i + 1;
	}
	list = list->next;
    }
    osort (&newlist, false);
    return newlist;
}

termptr
sladd (termptr a, termptr b)
{
    register termptr R203;
    termptr head, lastptr, point;

    snu (&lastptr);
    head = lastptr;
    head->next = NULL;
    while (a != NULL) {
	snu (&point);
	(*point) = (*a);
	a = a->next;
	lastptr->next = point;
	lastptr = point;
    }
    while (b != NULL) {
	snu (&point);
	(*point) = (*b);
	b = b->next;
	lastptr->next = point;
	lastptr = point;
    }
    lastptr = head;
    head = head->next;
    R203 = head;
    dispsfn (&lastptr);
    return R203;
}

termptr
qmult (termptr list)
{
    termptr newlist, temp = NULL;
    int k, p;
    register int i;

    newlist = NULL;
    while (list != NULL) {
	if (newlist == NULL) {
	    snu (&newlist);
	    temp = newlist;
	} else {
	    snu (&temp->next);
	    temp = temp->next;
	}
	(*temp) = (*list);
	p = 1;
	k = (len (&temp->val) - qsn + 1) / 2;
	if ((k > 0))
	    for (i = 1; i <= k; i++) {
		p = p * 2;
	    }
	temp->mult = p * temp->mult;
	list = list->next;
    }
    sort (&newlist, false);
    return newlist;
}

void
zraise (termptr list, int e, int f, int s)
{
    termptr lscan, addin;
    int x;

    lscan = list;
    while (lscan != NULL) {
	register termptr W9 = &(*lscan);
	x = qlen (W9->val);
	if ((x != 0) && (f <= x)) {
	    snu (&addin);
	    addin->val = W9->val;
	    addin->mult = W9->mult;
	    addin->next = W9->next;
	    if (e > 0)
		addin->val.A[e] = W9->val.A[e] + 1;
	    addin->val.A[f] = W9->val.A[f] - 1;
	    addin->mult = s * W9->mult;
	    W9->next = addin;
	    lscan = addin->next;
	} else
	    lscan = W9->next;
    }
}

termptr
xraise (termptr list, int s)
{
    termptr newlist;
    int e, k;
    register int f;

    k = qlen (list->val);
    newlist = sfncopy (list);
    if (k > 1) {
	for (f = k; f >= 2; f--) {
	    e = f;
	    do {
		e = e - 1;
		zraise (newlist, e, f, s);
	    }
	    while (!(e == 1));
	}
    }
    if ((redu) && (k != 0))
	for (f = k; f >= 1; f--) {
	    e = 0;
	    zraise (newlist, e, f, s);
	}
    if ((s == 1))
	qstndise (&newlist);
    else
	stndise (&newlist);
    sort (&newlist, true);
    return newlist;
}

termptr
craise (termptr list, int s)
{
    termptr newlist, temp, temp1, dummy;

    newlist = NULL;
    while (list != NULL) {
	/*register termptr W14 = &(*list); *//*12/12/95 */

	snu (&temp1);
	(*temp1) = (*list);
	temp1->next = NULL;	/*q */
	if ((s == 1) && (temp1->val.A[1] < qlen (temp1->val)) && (!redu))
	    conjgte (&temp1->val);
	temp = xraise (temp1, s);
	dummy = ladd (newlist, temp);
	ldisp (&temp);
	if (newlist != NULL)
	    ldisp (&newlist);
	newlist = dummy;
	ldisp (&temp1);
	list = list->next;
    }
    sort (&newlist, false);
    return newlist;
}

termptr
nraise (termptr list, int e, int f, int s)
{
    termptr temp, temp1, temp2;
    int a, l;

    snu (&temp);
    (*temp) = (*list);
    temp->next = NULL;
    l = qlen (temp->val);
    a = -(l - f);
    if ((l >= f) && (temp->val.A[f] > a)) {
	temp1 = temp;
	do {
	    zraise (temp1, e, f, s);
	    temp2 = temp1->next;
	    temp1 = temp2;
	} while (!(temp1->val.A[f] == a));
    }
    return temp;
}

termptr
yraise (termptr list, int e, int f, int s)
{
    termptr newlist, temp, dummy, list1;

    newlist = NULL;
    while (list != NULL) {
	snu (&list1);
	(*list1) = (*list);
	list1->next = NULL;
	temp = nraise (list1, e, f, s);
	dummy = sladd (newlist, temp);
	if (newlist != NULL)
	    ldisp (&newlist);
	ldisp (&list1);
	newlist = dummy;
	ldisp (&temp);
	list = list->next;
    }
    Qsort (&newlist, false);
    return newlist;
}

termptr
ayraise (termptr list, int s)
{
    termptr newlist, temp4, temp1, temp;
    int k, e;
    register int f;

    snu (&temp4);
    (*temp4) = (*list);
    k = qlen (temp4->val);	/*q */
    newlist = NULL;
    if (k != 0) {
	for (f = k; f >= 2; f--) {
	    e = f;
	    do {
		e = e - 1;

		temp = yraise (temp4, e, f, s);

		temp1 = sladd (newlist, temp);
		ldisp (&temp4);
		ldisp (&temp);
		temp4 = temp1;
	    } while (e != 1);
	}
	if (redu) {
	    for (f = k; f >= 1; f--) {
		e = 0;
		temp = yraise (temp4, e, f, s);
		ldisp (&temp4);
		temp4 = sladd (newlist, temp);
		ldisp (&temp);
	    }
	}
    }
    stndise (&temp4);
    sort (&temp4, false);
    return temp4;
}

termptr
lraise (termptr list, int s)
{
    termptr newlist, temp, temp1, dummy;

    newlist = NULL;
    while (list != NULL) {
	snu (&temp1);
	(*temp1) = (*list);
	temp1->next = NULL;
	temp = ayraise (temp1, s);
	dummy = ladd (newlist, temp);
	if (newlist != NULL)
	    ldisp (&newlist);
	newlist = dummy;
	ldisp (&temp);
	ldisp (&temp1);
	list = list->next;
    }
    sort (&newlist, false);
    return newlist;
}

termptr
oqinner (termptr list1, termptr list2)
{
    register termptr R211;
    termptr newlist, nlist1, nlist2;

    snu (&nlist1);
    snu (&nlist2);
    snu (&newlist);
    nlist1 = lraise (list1, -1);
    nlist2 = lraise (list2, -1);
    newlist = linner (nlist1, nlist2);
    ldisp (&nlist1);
    ldisp (&nlist2);
    R211 = craise (newlist, 1);
    ldisp (&newlist);
    return R211;
}

termptr
sqinner (termptr slist1, termptr qlist2)
{
    register termptr R212;
    termptr newlist, nlist2;

    nlist2 = lraise (qlist2, -1);
    newlist = linner (slist1, nlist2);
    ldisp (&nlist2);
    R212 = craise (newlist, 1);
    ldisp (&newlist);
    return R212;
}

termptr
rqinner (termptr slist1, termptr qlist2)
{
    register termptr R213;
    termptr newlist, nlist2;

    redu = true;
    nlist2 = lraise (qlist2, -1);
    newlist = rinner (slist1, nlist2);
    ldisp (&nlist2);
    R213 = craise (newlist, 1);
    ldisp (&newlist);
    redu = false;
    return R213;
}

void
mmsort (termptr * list, int k)
{
    int i, j, d, l, n, a, c;
    termptr temp;

    temp = (*list);
    while (temp != NULL) {
	i = qlen (temp->val);	/*q */
	n = wtfrm (&temp->val);
	j = n % 2;
	l = (i - j) / 2;
	if (l == 0)
	    a = 1;
	else {
	    c = 0;
	    d = 1;
	    do {
		d = d * 2;
		c = c + 1;
	    }
	    while (!(c == l));
	    a = d;
	}
	temp->mult = (a * temp->mult) / k;
	temp = temp->next;
    }
}

termptr
lsqinner (frame val1, frame val2, int n)
{
    termptr list1, list2, newlist;
    int k, m, l, q, a, b, c, d;
    register int i;

    snu (&list1);
    list1->mult = 1;
    d = n;
    snu (&list2);
    list2->mult = 1;
    list2->val = val2;
    m = wtfrm (&val1);
    l = len (&val1);
    l = l + 1;
    q = d % 2;
    a = (l - q) / 2;
    if (a == 0)
	k = 1;
    else {
	b = 0;
	c = 1;
	do {
	    c = c * 2;
	    b = b + 1;
	}
	while (!(b == a));
	k = c;
    }
    list1->val.A[1] = n - m;
    for (i = 1; i <= maxdim - 1; i++) {
	list1->val.A[i + 1] = val1.A[i];
    }
    newlist = sqinner (list2, list1);
    ldisp (&list1);
    ldisp (&list2);
    mmsort (&newlist, k);
    return newlist;
}

void
signsort (termptr list)
{
    while (list != NULL) {
	list->mult = -list->mult;
	list = list->next;
    }
}

termptr
lsinner (frame * val1, frame * val2, int s, int t, int n)
{
    register termptr R215;
    termptr list1, list2, list, newlist;
    int p, m, l;
    register int i;

    snu (&list1);
    list1->mult = 1;
    snu (&list2);
    list2->mult = 1;
    m = wtfrm (&(*val1));
    p = wtfrm (&(*val2));
    list1->val.A[1] = n - m;
    for (i = 1; i <= maxdim - 1; i++) {
	list1->val.A[i + 1] = val1->A[i];
    }
    list2->val.A[1] = n - p;
    for (i = 1; i <= maxdim - 1; i++) {
	list2->val.A[i + 1] = val2->A[i];
    }
    rep = true;
    if ((((t == 0) && (s != 0)) || ((s == 0) && (t != 0))))
	l = 2;
    else
	l = 4;
    list = qinner (list1, list2);
    rep = false;
    if (((s != 0) && (t != 0))) {
	if (m == p) {
	    newlist = difchrq (&list2);
	    if (s != t)
		signsort (newlist);
	    add (&list, &newlist);
	    ldisp (&newlist);
	}
    }
    sort (&list, true);
    if ((s != 0) || (t != 0))
	nsort (list, l);
    ldisp (&list1);
    ldisp (&list2);
    R215 = list;
    return R215;
}

void
nsort (termptr list, int s)
{
    while (list != NULL) {
	list->mult = list->mult / s;
	list = list->next;
    }
}

void
dsecond (termptr * list1, frame * nu, int d, int i, int c, int s)
{
    int m, l, y, x, t, z;
    register int b;
    register int k;
    frame rho, v;
    termptr temp, temp3, newlist, dummy;
    bool test;

    l = qlen ((*nu));
    z = l;
    newlist = NULL;
    progress ();
    v = (*nu);
    temp3 = (*list1);
    t = 1;
    while ((*list1) != NULL) {
	for (k = 1; k <= maxdim; k++) {
	    rho.A[k] = 0;
	}
	rho = (*list1)->val;
	x = rho.A[i - s];
	y = rho.A[i - s - 1];
	for (b = 1; b <= z; b++) {
	    m = rest (&rho, b);
	    test = true;
	    if (m < v.A[b]) {
		if (d == 1) {
		    if (b >= rho.A[i - 1]) {
			snu (&temp);
			temp->val = rho;
			temp->mult = 1;
			temp->val.A[i] = b;
			test = false;
		    }
		} else if (i == c) {
		    if ((rest (&rho, t) != v.A[1]) && (b == 1)) {
			snu (&temp);
			temp->val = rho;
			temp->mult = 1;
			temp->val.A[i] = 1;
			test = false;
		    } else if ((rest (&rho, t) == v.A[1]) && (b == x)) {
			snu (&temp);
			temp->val = rho;
			temp->mult = 1;
			temp->val.A[i] = b;
			test = false;
		    } else if ((rest (&rho, t) == v.A[1]) && (b > x)
			       && (rest (&rho, x) == v.A[x])) {
			snu (&temp);
			temp->val = rho;
			temp->mult = 1;
			temp->val.A[i] = b;
			test = false;
		    }
		} else if ((i != c) && (b != y)) {
		    if ((b >= x) && (b >= rho.A[i - 1])) {
			if ((b == x) && (b != 1)) {
			    snu (&temp);
			    temp->val = rho;
			    temp->mult = 1;
			    temp->val.A[i] = b;
			    test = false;
			} else if (b > x) {
			    if (x != 1) {
				if (rest (&rho, x) == v.A[x]) {
				    snu (&temp);
				    temp->val = rho;
				    temp->mult = 1;
				    temp->val.A[i] = b;
				    test = false;
				} else if ((rest (&rho, x) != v.A[x])
					   && (x == y)) {
				    snu (&temp);
				    temp->val = rho;
				    temp->mult = 1;
				    temp->val.A[i] = b;
				    test = false;
				}
			    } else if (x == 1) {
				snu (&temp);
				temp->val = rho;
				temp->mult = 1;
				temp->val.A[i] = b;
				test = false;
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
	}
	(*list1) = (*list1)->next;
    }
    ldisp (&temp3);
    (*list1) = newlist;
}

void
dfirst (termptr * list, frame * nu, int d, int c, int q, int s)
{
    register int i;
    for (i = c; i <= q; i++) {
	{
	    /*register termptr W33 = &(*(*list)); *//*12/12/95 */

	    dsecond (&(*list), &(*nu), d, i, c, s);
	}
    }
}

void
dsort (termptr * list, termptr * lambda, frame mu)
{
    typedef struct {
	bframe A[maxdim - 1 + 1];
    } ray;
    termptr temp, newlist, dummy, list1;
    bool test;
    frame v;
    int a, b, k, l, m, x, lastStep;
    register int j;
    register int i;
    register int y;
    ray w;

    l = qlen ((*lambda)->val);
    list1 = (*list);
    newlist = NULL;
    while ((*list) != NULL) {
	snu (&temp);
	(*temp) = (*(*list));
	temp->next = NULL;
	for (y = 1; y <= maxl; y++) {
	    v.A[y] = mu.A[y];
	}
	for (i = 1; i <= maxl; i++) {
	    for (j = 1; j <= maxl; j++) {
		w.A[i - 1].A[j] = 0;
	    }
	}
	k = 0;
	for (i = 1; i <= l; i++) {
	    {
		lastStep = (*lambda)->val.A[i];
		for (j = 1; j <= lastStep; j++) {
		    k = k + 1;
		    w.A[i - 1].A[j] = temp->val.A[k];
		}
	    }
	}
	m = qlen (mu);
	test = true;
	i = 0;
	do {
	    i = i + 1;
	    j = (*lambda)->val.A[i] + 1;
	    do {
		j = j - 1;
		{
		    a = i;
		    b = j;
		    x = w.A[a - 1].A[b];
		    if (v.A[x] != 0) {
			v.A[x] = v.A[x] - 1;
			if (v.A[x] != 0) {
			    do {
				if ((b != 1) && (w.A[a - 1].A[b - 1] == x)) {
				    b = b - 1;
				    v.A[x] = v.A[x] - 1;
				} else if ((a != l)
					   && (w.A[a + 1 - 1].A[b] == x)) {
				    a = a + 1;
				    v.A[x] = v.A[x] - 1;
				} else
				    test = false;
			    }
			    while (!(((v.A[x] == 0) || (test == false))));
			}
			if (v.A[x] == 0)
			    m = m - 1;
		    }
		}
	    }
	    while (!(((j == 1) || (test == false) || (m == 0))));
	}
	while (!(((i == l) || (test == false) || (m == 0))));
	if (test) {
	    dummy = ladd (newlist, temp);
	    if (newlist != NULL)
		ldisp (&newlist);
	    newlist = dummy;
	}
	ldisp (&temp);
	(*list) = (*list)->next;
    }
    ldisp (&list1);
    (*list) = newlist;
}

termptr
dgenerate (termptr * lambda, frame * mu)
{
    termptr dummy, list1, list, temp, temp2;
    int q, r, s, c, d, l, n;
    bool frst;

    l = qlen ((*lambda)->val);
    snu (&temp);
    temp->mult = 1;
    temp->val.A[1] = 1;
    r = (*lambda)->val.A[1];
    n = r + l - 1;
    if (mu->A[1] <= n) {
	if (r > 1)
	    dfirst (&temp, &(*mu), 1, 2, r, 0);
	d = 1;
	c = r;
	if (l > 1) {
	    do {
		d = d + 1;
		list1 = NULL;
		c = c + 1;
		q = (c + (*lambda)->val.A[d] - 1);
		s = (*lambda)->val.A[d - 1];
		frst = true;
		list1 = NULL;
		temp2 = temp;
		while (temp != NULL) {
		    /*register termptr W44 = &(*temp); *//*12/12/95 */

		    snu (&list);
		    (*list) = (*temp);
		    list->next = NULL;
		    {
			dfirst (&list, &(*mu), d, c, q, s);
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
    }
    dsort (&temp, &(*lambda), (*mu));
    return temp;
}

int
icoeffq (termptr * lambda, termptr * list)
{
    register int R217;
    int a, k, l, m, x, y, z;
    register int j;
    register int i;
    register int n;
    frame v;

    l = qlen ((*lambda)->val);
    for (n = 1; n <= maxdim; n++) {
	v.A[n] = 0;
    }
    y = 0;
    for (i = 1; i <= l; i++) {
	z = y + 1;
	y = y + (*lambda)->val.A[i];
	for (j = z; j <= y; j++) {
	    k = (*list)->val.A[j];
	    if (j == z)
		v.A[k] = v.A[k] + 1;
	    else if (((*list)->val.A[j - 1] != k) && (j != z))
		v.A[k] = v.A[k] + 1;
	}
    }
    m = 0;
    a = qlen (v);
    for (i = 1; i <= a; i++) {
	if (!(bool) ((v.A[i]) & 1))
	    m = m + 1;
    }
    if (m == 0)
	x = 1;
    else if ((bool) ((m) & 1))
	x = -1;
    else
	x = 1;
    R217 = x;
    return R217;
}

int
dcoeffq (termptr * lambda, frame * mu)
{
    register int R218;
    termptr list, list1;
    int i, j;

    i = 0;
    list = dgenerate (&(*lambda), &(*mu));
    list1 = list;
    while (list != NULL) {
	j = icoeffq (&(*lambda), &list);
	i = i + j;
	list = list->next;
    }
    ldisp (&list1);
    R218 = i;
    return R218;
}

termptr
difchrq (termptr * qfn1)
{
    register termptr R219;
    termptr list, list1, temp, newlist, dummy;
    frame v;
    int k, l, m, n, w, x;
    register int i;

    list = leqwt ((*qfn1));
    for (i = 1; i <= maxdim; i++) {
	v.A[i] = 0;
    }
    list1 = list;
    v = (*qfn1)->val;
    newlist = NULL;
    l = qlen ((*qfn1)->val);
    w = wtfrm (&(*qfn1)->val);
    k = w - l;
    x = k + 1;
    if ((bool) ((k) & 1)) {
	if (x % 4 == 0)
	    m = 2;
	else
	    m = -2;
	while (list != NULL) {
	    snu (&temp);
	    (*temp) = (*list);
	    temp->next = NULL;
	    n = dcoeffq (&temp, &v);
	    if (n != 0) {
		temp->mult = n * m * (*qfn1)->mult;
		dummy = ladd (newlist, temp);
		if (newlist != NULL)
		    ldisp (&newlist);
		newlist = dummy;
	    }
	    ldisp (&temp);
	    list = list->next;
	}
    }
    ldisp (&list1);
    R219 = newlist;
    return R219;
}
