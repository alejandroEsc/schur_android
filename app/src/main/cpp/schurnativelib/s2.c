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
#include "r.h"
#include "g.h"
#include "s.h"
#include "label.h"
#include "skew.h"
#include "write.h"
#include "outerskew.h"
#include "s2.h"

prodtype pcurrent;

ocharptr
scalarinner (int n, int k)
{
    termptr sfn1, sfn, temp, newlist;
    register int i;
    ocharptr result;

    newlist = NULL;
    snu (&sfn);
    sfn->mult = 1;
    for (i = 1; i <= maxdim; i++) {
        sfn->val.A[i] = 0;
    }
    sfn->val.A[1] = n;
    sfn1 = useseries ('f', sfn, true, false, n);
    dispsfn (&sfn);
    schur_restrict (&sfn1, -n, 'w');
    schur_restrict (&sfn1, k, 'l');
    while (sfn1 != NULL) {
        register termptr W3 = &(*sfn1);

        temp = inner (W3->val, W3->val);
        add (&newlist, &temp);
        ldisp (&temp);
        sfn1 = W3->next;
    }
    ldisp (&sfn1);
    sort (&newlist, true);
    result = sfntochrc (newlist, false, ' ');
    ldisp (&newlist);
    return result;
}



termptr
multmono (termptr * list1, termptr * list2, int n)
{
    termptr ulist, tlist, newlist, temp;
    register int i;

    newlist = NULL;
    ulist = (*list1);
    while (ulist != NULL) {
        tlist = (*list2);
        while (tlist != NULL) {
            register termptr W2 = &(*tlist);

            snu (&temp);
            temp->mult = W2->mult * ulist->mult;
            temp->val = nolls;
            temp->next = NULL;
            for (i = 1; i <= n; i++) {
                temp->val.A[i] = ulist->val.A[i] + W2->val.A[i];
            }
            add (&newlist, &temp);
            tlist = W2->next;
        }
        /*sort(&newlist, true); */
        ulist = ulist->next;
    }
    sort (&newlist, true);
    return newlist;
}

termptr
genprod (int *n)
{
    termptr ttemp, newlist, tlist, temp, tempi, dummy;
    register int i;
    snu (&newlist);
    newlist->mult = 1;
    newlist->val = nolls;
    snu (&temp);
    temp->mult = -1;
    temp->val = nolls;
    temp->val.A[(*n)] = 1;
    temp->next = NULL;
    for (i = 1; i <= (*n) - 1; i++) {
        ttemp = temp;
        snu (&tempi);
        tempi->mult = 1;
        tempi->val = nolls;
        tempi->val.A[i] = 1;
        tempi->next = NULL;
        tlist = ladd (ttemp, tempi);
        dummy = newlist;
        newlist = multmono (&dummy, &tlist, (*n));
        ldisp (&dummy);
        ldisp (&tlist);
        dispsfn (&tempi);
    }
    dummy = newlist;
    newlist = multmono (&dummy, &dummy, (*n));
    ldisp (&dummy);
    dispsfn (&temp);
    return newlist;
}

termptr
sfnmon (termptr * list1, termptr * list2, int n)
{
    termptr newlist, temp, tlist, mlist;
    int k;

    newlist = NULL;
    tlist = (*list1);
    k = 1;
    while (tlist != NULL) {
        register termptr W7 = &(*tlist);

        snu (&mlist);
        mlist->mult = W7->mult;
        mlist->val = W7->val;
        mlist->next = NULL;
        temp = multmono (&mlist, &(*list2), n);
        dispsfn (&mlist);
        stndise (&temp);
        add (&newlist, &temp);
        k = k + 1;
        tlist = W7->next;
        sort (&newlist, true);
    }
    sort (&newlist, true);
    return newlist;
}

/*
termptr
router (termptr * list1, termptr * list2, int n)
{
  termptr temp, newlist, C40_short, sfnp, psfn, sfnplist;
  int last, k;
  register int i;

  newlist = NULL;
  k = 0;
  while ((*list2) != NULL)
    {
      register termptr W8 = &(*(*list2));

      snu (&C40_short);
      C40_short->mult = W8->mult;
      C40_short->val = W8->val;
      last = W8->val.A[n];
      C40_short->val.A[n] = 0;
      sfnp = louter2 ((*list1), C40_short, n - 1);
      dispsfn (&C40_short);
      sfnplist = NULL;
      temp = sfnp;
      while (temp != NULL)
	{
	  register termptr W9 = &(*temp);

	  snu (&psfn);
	  psfn->mult = W9->mult;
	  psfn->val = nolls;
	  for (i = 2; i <= n; i++)
	    {
	      psfn->val.A[i] = W9->val.A[i - 1];
	    }
	  psfn->val.A[1] = last;
	  stndise (&psfn);
	  add (&sfnplist, &psfn);
	  temp = W9->next;
	}
      ldisp (&sfnp);
      add (&newlist, &sfnplist);
      k = k + 1;
      (*list2) = W8->next;
    }
  sort (&newlist, true);
  return newlist;
}
*/

void
restrict2 (termptr * list, int limitx, char action, bool exactComparison)
{
    bool okx;
    int f, l, wt, mu, rankk;
    termptr ptr1, ptr2;
    bool sp, fp, oe, rk;

    ptr1 = (*list);
    ptr2 = (*list);
    action = locase (action);
    while (ptr1 != NULL) {
        l = len (&ptr1->val);
        wt = wtfrm (&ptr1->val);
        mu = ptr1->mult;
        f = ptr1->val.A[1];
        rankk = frank (ptr1->val);
        if ((ptr1->val.A[1] <= limitx))
            fp = true;
        else
            fp = false;
        if ((bool) ((wt) & 1))
            oe = true;
        else
            oe = false;
        if ((bool) ((rankk) & 1))
            rk = true;
        else
            rk = false;
        if (((f + ptr1->val.A[2]) <= limitx))
            sp = true;
        else
            sp = false;
        if (!exactComparison)
            okx = (bool) (((l <= limitx) && (action == 'l'))
                          || ((wt <= limitx) && (action == 'w'))
                          || ((mu <= limitx) && (action == 'm'))
                          || ((oe == true)
                              && (action == 'o'))
                          || ((oe == false) && (action == 'e'))
                          || ((rk == true)
                              && (action == 'q'))
                          || ((rk == false) && (action == 'p'))
                          || ((fp == true)
                              && (action == 'f'))
                          || ((sp == true) && (action == 'z'))
                          || ((rankk <= limitx) && (action == 'x')));
        else
            okx = (bool) (((l == limitx) && (action == 'l'))
                          || ((wt == limitx) && (action == 'w'))
                          || ((mu == limitx) && (action == 'm'))
                          || ((f == limitx)
                              && (action == 'f'))
                          || ((action == 'c') && (ptr1->val.A[1] <= 1))
                          || ((action == 's') && (ptr1->val.A[1] > 1)
                              && (ptr1->val.A[2] == 1)) || ((action == 't')
                                                            && (ptr1->val.
                                                                A[2] > 1)
                                                            && (ptr1->val.
                                                                A[3] == 0))
                          || ((action == 'd') && (ptr1->val.A[1] > 1)
                              && (ptr1->val.A[2] >= 2)
                              && (ptr1->val.A[3] <= 2))
                          || ((action == 'r') && (ptr1->val.A[2] == 0))
                          || ((rankk == limitx) && (action == 'x')));
        if (!okx) {
            if (ptr1 != ptr2) {
                ptr2->next = ptr1->next;
                dispsfn (&ptr1);
                ptr1 = ptr2->next;
            } else {
                (*list) = ptr1->next;
                dispsfn (&ptr1);
                ptr1 = (*list);
                ptr2 = (*list);
            }
        } else {
            ptr2 = ptr1;
            ptr1 = ptr1->next;
        }
    }
    sort (&(*list), true);
}

void
schur_restrict (termptr * list, int limitx, char action)
{
    bool okx;
    int f, av, l, wt, mu, rankk;
    termptr ptr1, ptr2;
    bool sp, fp, oe, rk;

    ptr1 = (*list);
    ptr2 = (*list);
    av = abs (limitx);
    action = locase (action);
    while (ptr1 != NULL) {
        l = len (&ptr1->val);
        wt = wtfrm (&ptr1->val);
        mu = ptr1->mult;
        f = ptr1->val.A[1];
        rankk = frank (ptr1->val);
        if ((ptr1->val.A[1] <= limitx))
            fp = true;
        else
            fp = false;
        if ((bool) ((wt) & 1))
            oe = true;
        else
            oe = false;
        if ((bool) ((rankk) & 1))
            rk = true;
        else
            rk = false;
        if (((f + ptr1->val.A[2]) <= limitx))
            sp = true;
        else
            sp = false;
        if (limitx >= 0)
            okx = (bool) (((l <= av) && (action == 'l'))
                          || ((wt <= av) && (action == 'w')) || ((mu <= av)
                                                                 && (action ==
                                                                     'm'))
                          || ((oe == true) && (action == 'o'))
                          || ((oe == false) && (action == 'e'))
                          || ((rk == true) && (action == 'q'))
                          || ((rk == false) && (action == 'p'))
                          || ((fp == true) && (action == 'f'))
                          || ((sp == true)
                              && (action == 'z'))
                          || ((rankk <= av) && (action == 'x')));
        else
            okx = (bool) (((l == av) && (action == 'l'))
                          || ((wt == av) && (action == 'w')) || ((mu == av)
                                                                 && (action ==
                                                                     'm'))
                          || ((f == av) && (action == 'f'))
                          || ((action == 'c')
                              && (ptr1->val.A[1] <= 1))
                          || ((action == 's') && (ptr1->val.A[1] > 1)
                              && (ptr1->val.A[2] == 1)) || ((action == 't')
                                                            && (ptr1->
                                                                val.A[2] > 1)
                                                            && (ptr1->
                                                                val.A[3] ==
                                                                0))
                          || ((action == 'd') && (ptr1->val.A[1] > 1)
                              && (ptr1->val.A[2] >= 2)
                              && (ptr1->val.A[3] <= 2)) || ((action == 'r')
                                                            && (ptr1->
                                                                val.A[2] ==
                                                                0))
                          || ((rankk == av) && (action == 'x')));
        if (!okx) {
            if (ptr1 != ptr2) {
                ptr2->next = ptr1->next;
                dispsfn (&ptr1);
                ptr1 = ptr2->next;
            } else {
                (*list) = ptr1->next;
                dispsfn (&ptr1);
                ptr1 = (*list);
                ptr2 = (*list);
            }
        } else {
            ptr2 = ptr1;
            ptr1 = ptr1->next;
        }
    }
    sort (&(*list), true);
}

void
rrestrict (ocharptr * list, int limitx, char action)
{
    bool okx;
    int f, av, l, wt, mu;
    ocharptr ptr1, ptr2;
    bool fp, oe;

    ptr1 = (*list);
    ptr2 = (*list);
    av = abs (limitx);
    action = locase (action);
    while (ptr1 != NULL) {
        l = len (&ptr1->val);
        wt = wtfrm (&ptr1->val);
        mu = ptr1->mult;
        f = ptr1->val.A[1];
        if ((bool) ((wt) & 1))
            oe = true;
        else
            oe = false;
        if ((ptr1->val.A[1] <= limitx))
            fp = true;
        else
            fp = false;
        if (limitx >= 0)
            okx = (bool) (((l <= av) && (action == 'l'))
                          || ((wt <= av) && (action == 'w')) || ((mu <= av)
                                                                 && (action ==
                                                                     'm'))
                          || ((oe == true) && (action == 'o'))
                          || ((oe == false) && (action == 'e'))
                          || ((fp == true) && (action == 'f')));
        else
            okx = (bool) (((l == av) && (action == 'l'))
                          || ((wt == av) && (action == 'w')) || ((mu == av)
                                                                 && (action ==
                                                                     'm'))
                          || ((f == av) && (action == 'f')));
        if (!okx) {
            if (ptr1 != ptr2) {
                ptr2->next = ptr1->next;
                dispchr (&ptr1);
                ptr1 = ptr2->next;
            } else {
                (*list) = ptr1->next;
                dispchr (&ptr1);
                ptr1 = (*list);
                ptr2 = (*list);
            }
        } else {
            ptr2 = ptr1;
            ptr1 = ptr1->next;
        }
    }
    osort (&(*list), true);
}

void
prestrict (prodtype * list, int limitx, char action)
{
    bool okx;
    int av, mu;
    prodtype ptr1, ptr2;

    ptr1 = (*list);
    ptr2 = (*list);
    av = abs (limitx);
    action = locase (action);
    while (ptr1 != NULL) {
        mu = ptr1->mult;
        if (limitx >= 0)
            okx = (bool) ((mu <= av) && (action == 'm'));
        else
            okx = (bool) ((mu == av) && (action == 'm'));
        if (!okx) {
            if (ptr1 != ptr2) {
                ptr2->next = ptr1->next;
                dispprod (&ptr1);
                ptr1 = ptr2->next;
            } else {
                (*list) = ptr1->next;
                dispprod (&ptr1);
                ptr1 = (*list);
                ptr2 = (*list);
            }
        } else {
            ptr2 = ptr1;
            ptr1 = ptr1->next;
        }
    }
    schur_psort (&(*list), true);
}

void
smult (termptr list)
{
    termptr ptr;

    ptr = list;
    while (ptr != NULL) {
        register termptr W12 = &(*ptr);

        W12->mult = 1;
        ptr = W12->next;
    }
}

void
scollctgarbage (void)
{
    register int i;

    if (srjctindex > 0) {
        ldisp (&scurrent);
        scurrent = sreject.A[srjctindex - 1];
        sreject.A[srjctindex - 1] = NULL;
        for (i = 1; i <= srjctindex - 1; i++) {
            ldisp (&sreject.A[i - 1]);
        }
        srjctindex = 0;
    }
}

void
collectgarbage (void)
{
    register int i;

    if (rjctindex > 0) {
        odisp (&current);
        current = reject.A[rjctindex - 1];
        reject.A[rjctindex - 1] = NULL;
        for (i = 1; i <= rjctindex - 1; i++) {
            odisp (&reject.A[i - 1]);
        }
        rjctindex = 0;
        for (i = 1; i <= rjctindex2; i++)
            odisp (&reject2.A[i - 1]);
        rjctindex2 = 0;
    }
}

void
pcollectgarbage (void)
{
    register int i;

    if (prjctindex > 0) {
        pdisp (&pcurrent);
        pcurrent = preject.A[prjctindex - 1];
        preject.A[prjctindex - 1] = NULL;
        for (i = 1; i <= prjctindex - 1; i++) {
            pdisp (&preject.A[i - 1]);
        }
        prjctindex = 0;
    }
}

termptr
ql (char signx, int n)
{
    termptr list, lastptr;
    int i, j;
    register int k;
    frame partq;

    lastptr = NULL;
    i = 1;
    partq = nolls;
    if (signx == '+')
        j = 0;
    else
        j = 1;
    while (j <= n) {
        for (k = j; k <= maxdim; k++) {
            partq.A[k] = 0;
        }
        for (k = i; k <= j; k++) {
            partq.A[k] = 1;
        }
        snu (&list);
        {
            register termptr W23 = &(*list);

            W23->val = partq;
            W23->mult = 1;
            W23->next = lastptr;
        }
        lastptr = list;
        i = j + 1;
        j = j + 2;
    }
    return lastptr;
}

termptr
xv (char signx, int n)
{
    termptr list, lastptr;
    int i, j, k, l;
    register int m;
    frame partxv;

    lastptr = NULL;
    i = 1;
    partxv = nolls;
    if (signx == '+')
        j = 0;
    else
        j = 1;
    while (j <= n) {
        for (m = j; m <= maxdim; m++) {
            partxv.A[m] = 0;
        }
        for (m = i; m <= j; m++) {
            partxv.A[m] = 2;
        }
        k = j;
        l = j;
        do {
            for (m = k + 1; m <= l; m++) {
                partxv.A[m] = 1;
            }
            snu (&list);
            {
                register termptr W30 = &(*list);

                W30->val = partxv;
                W30->mult = 1;
                W30->next = lastptr;
            }
            lastptr = list;
            k = l;
            l = l + 2;
        }
        while (!(l > n));
        i = j + 1;
        j = j + 2;
    }
    return lastptr;
}

prodtype
crunchup (prodtype pr, int n)
{
    prodtype h;
    register int j;

    h = pr;
    while (h != NULL) {
        register prodtype W31 = &(*h);

        odisp (&W31->prods.A[n - 1]);
        for (j = n; j <= nprod - 1; j++) {
            W31->prods.A[j - 1] = W31->prods.A[j + 1 - 1];
        }
        W31->prods.A[nprod - 1] = NULL;
        h = h->next;
    }
    for (j = n; j <= nprod - 1; j++) {
        currgrp.A[j - 1] = currgrp.A[j + 1 - 1];
    }
    nprod = nprod - 1;
    putgroup (currgrp);
    schur_psort (&pr, true);
    return pr;
}

ocharptr
prodcon (prodtype list)
{
    ocharptr result, ptr = NULL;

    result = NULL;
    if (nprod > 1)
        error (GROUP_NOT_SET, 1);
    else
        while (list != NULL) {
            if (result == NULL) {
                cnu (&result);
                ptr = result;
            } else {
                cnu (&ptr->next);
                ptr = ptr->next;
            }
            (*ptr) = (*list->prods.A[1 - 1]);
            ptr->mult = list->mult;
            list = list->next;
        }
    return result;
}

termptr
prodsfn (prodtype list)
{
    termptr result, ptr = NULL;

    result = NULL;
    if (nprod > 1)
        error (GROUP_NOT_SET, 1);
    else
        while (list != NULL) {
            if (result == NULL) {
                snu (&result);
                ptr = result;
            } else {
                snu (&ptr->next);
                ptr = ptr->next;
            }
            ptr->val = list->prods.A[1 - 1]->val;
            ptr->mult = list->mult;
            list = list->next;
        }
    return result;
}

termptr
repsfn (ocharptr list)
{
    termptr result, ptr = NULL;

    result = NULL;
    while (list != NULL) {
        if (result == NULL) {
            snu (&result);
            ptr = result;
        } else {
            snu (&ptr->next);
            ptr = ptr->next;
        }
        ptr->val = list->val;
        ptr->mult = list->mult;
        list = list->next;
    }
    return result;
}

void
unlimit (ocharptr * list, int n)
{
    ocharptr newlist, temp;
    int i;
    bool test;

    newlist = NULL;
    while ((*list) != NULL) {
        register ocharptr W36 = &(*(*list));

        i = 1;
        test = true;
        if ((n < 0))
            do {
                if ((((n == -1) && (bool) ((W36->val.A[i]) & 1))
                     || ((n == -2) && (!(bool) ((W36->val.A[i]) & 1)))))
                    test = false;
                else
                    i = i + 1;
            }
            while (!
                   ((W36->val.A[i] == 0) || (i == maxdim)
                    || (test == false)));
        else if ((n > 0))
            do {
                if ((W36->val.A[i] == n))
                    test = false;
                else
                    i = i + 1;
            }
            while (!
                   ((W36->val.A[i] == 0) || (i == maxdim)
                    || (test == false)));
        temp = W36->next;
        if ((test == false))
            dispchr (&(*list));
        else {
            W36->next = newlist;
            newlist = (*list);
        }
        (*list) = temp;
    }
    osort (&newlist, true);
    (*list) = newlist;
}

void
onexp (ocharptr * list, int n)
{
    ocharptr lscan, addin;

    lscan = (*list);
    while (lscan != NULL) {
        register ocharptr W37 = lscan;

        if ((n > len (&W37->val))) {
            cnu (&addin);
            addin->val = W37->val;
            addin->mult = W37->mult;
            addin->spin = W37->spin;
            addin->next = W37->next;
            addin->C6_double = W37->C6_double;
            addin->lab = '#';
            W37->next = addin;
            lscan = addin->next;
        } else
            lscan = W37->next;
    }
}

ocharptr
cspin (ocharptr list)
{
    ocharptr newlist, temp = NULL;

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
        if (temp->spin == true)
            temp->spin = false;
        else
            temp->spin = true;
        temp->lab = ' ';
        list = list->next;
    }
    osort (&newlist, false);
    return newlist;
}

ocharptr
muc (ocharptr list, groop grp, int m)
{
    ocharptr newlist, temp;

    newlist = NULL;
    temp = NULL;
    if ((grp.name == un) && (grp.rank == 1)) {
        while (list != NULL) {
            if (newlist == NULL) {
                cnu (&newlist);
                temp = newlist;
            } else {
                cnu (&temp->next);
                temp = temp->next;
            }
            (*temp) = (*list);
            if ((list->val.A[1] != 0))
                temp->val.A[1] = m * list->val.A[1];
            list = list->next;
        }
        osort (&newlist, false);
    } else {
        print ("mistake inappropriate group?\n");
        newlist = chrccopy (list);
    }
    return newlist;
}

ocharptr
umax (ocharptr list, groop grp, int m)
{
    ocharptr newlist, temp;

    newlist = NULL;
    temp = NULL;
    if ((grp.name == un) && (grp.rank == 1)) {
        while (list != NULL) {
            if (newlist == NULL) {
                cnu (&newlist);
                temp = newlist;
            } else {
                cnu (&temp->next);
                temp = temp->next;
            }
            (*temp) = (*list);
            if (temp->val.A[1] > m) {
                temp->mult = 0;
                temp->val.A[1] = 0;
            }
            list = list->next;
        }
        osort (&newlist, false);
    } else {
        print ("mistake inappropriate group?\n");
        newlist = chrccopy (list);
    }
    return newlist;
}

ocharptr
contrag (ocharptr list, groop grp)
{
    ocharptr newlist, temp;
    int k;
    register int i;

    newlist = NULL;
    temp = NULL;
    if (((grp.name == sung) || (grp.name == e6))
        || ((grp.name == un) && (grp.rank == 1))) {
        while (list != NULL) {
            register ocharptr W40 = &(*list);

            if (newlist == NULL) {
                cnu (&newlist);
                temp = newlist;
            } else {
                cnu (&temp->next);
                temp = temp->next;
            }
            (*temp) = (*list);
            if ((grp.name == un) && (grp.rank == 1)) {
                if ((W40->val.A[1] != 0))
                    temp->val.A[1] = -W40->val.A[1];
                list = list->next;
            } else if ((grp.name == sung) && (grp.rank > 1)) {
                for (i = 0; i <= grp.rank - 1; i++) {
                    k = grp.rank - i;
                    temp->val.A[k] = W40->val.A[1] - W40->val.A[i + 1];
                }
                list = list->next;
            } else if ((grp.name == e6)) {
                for (i = 1; i <= 5; i++) {
                    k = 7 - i;
                    temp->val.A[k] = W40->val.A[2] - W40->val.A[i + 2];
                }
                list = list->next;
            }
        }
        osort (&newlist, false);
    } else {
        print ("mistake inappropriate group?\n");
        newlist = chrccopy (list);
    }
    return newlist;
}

termptr
leqwt (termptr list)
{
    termptr result, subl;

    result = NULL;
    while (list != NULL) {
        subl = sameweight (wtfrm (&list->val), 1);
        merge (&result, &subl, true, true);
        list = list->next;
    }
    return result;
}

int
frank (frame f)
{
    register int i;

    i = 0;
    while (f.A[i + 1] >= i + 1)
        i = i + 1;
    return i;
}

termptr
outerx (int c, frame alpha, frame beta, int maxx)
{
    frame delta;
    int la, lb, d, wa, wb;
    termptr result;
    result = NULL;
    delta = nolls;
    wa = wtfrm (&alpha);        /* this seems to be need by inverseries */
    wb = wtfrm (&beta);
    if ((wa+wb) >= maxdim) {
       print ("Error : outer product cannot be made : partitions are too big (max sum=%d)\n", maxdim-1);
       return NULL;
    }

    if ((wa + wb) <= plwt) {
        la = len (&alpha);
        lb = len (&beta);
        if ((beta.A[1] > alpha.A[1])
            || ((beta.A[1] == alpha.A[1]) && (la < lb))) {
            delta = alpha;
            alpha = beta;
            beta = delta;
            d = la;
            la = lb;
            lb = d;
        }
        if (lb == 0) {
            snu (&result);
            result->val = alpha;
            result->mult = c;
            result->next = NULL;
        } else
            result = outerskew (c, alpha, beta, nolls, la, lb, maxx);
    } 
    return result;
}

termptr
outer2 (frame first, frame last, int rank)
{
    if ((len (&first) + len (&last)) < maxdim)
        return outerx (1, first, last, rank);
    else {
        print ("ERROR maxdim = %d exceeded\n", maxdim - 1);
        return NULL;
    }
}

termptr
outer (frame first, frame last)
{
    if (((wtfrm (&first) + wtfrm (&last)) <= plwt))
        return outer2 (first, last, plwt) /*maxdim) */ ;
    else 
        return NULL;
}

termptr
louter2 (termptr list1, termptr list2, int n)
{
    termptr templist, sublist, prodlist;
    char ch;

    prodlist = NULL;
    ch = ' ';
    if (n < maxdim) {
        while (list1 != NULL) {
            templist = list2;
            while (templist != NULL) {
                sublist =
                    outerx (list1->mult * templist->mult, list1->val,
                            templist->val, n);
                if (sslab) {
                    if (((templist->slab == '#') && (list1->slab == ' '))
                        || ((templist->slab == ' ') && (list1->slab == '#')))
                        ch = '#';
                    else
                        ch = ' ';
                    cslabel (sublist, ch);
                }
                add (&prodlist, &sublist);
                templist = templist->next;
                sort (&prodlist, true);
            }
            list1 = list1->next;
        }
    } else
        print ("ERROR maxdim = %d exceeded\n", maxdim - 1);
    return prodlist;
}

termptr
louter (termptr list1, termptr list2)
{

    return louter2 (list1, list2, maxdim - 1);
}

termptr
factorsof (int c, frame ptn)
{
    bool move, first;
    int getposn, putposn, altposn, elt, sigma, suffix, k, n;
    register int i, j;
    frame rho, perm, follow, countx;
    termptr head, tail;

    n = len (&ptn);
    sigma = c;
    head = NULL;
    tail = NULL;
    suffix = 2;
    first = true;
    for (i = 1; i <= maxdim; i++) {
        perm.A[i] = i;
        follow.A[i] = i + 1;
        countx.A[i] = 0;
        ptn.A[i] = ptn.A[i] - i;
    }
    while ((suffix <= n) || first) {
        progress ();
        if (!first) {
            for (j = 0; j <= (suffix / 2) - 1; j++) {
                k = perm.A[n - j];
                perm.A[n - j] = perm.A[n - suffix + j + 1];
                perm.A[n - suffix + j + 1] = k;
                sigma = -sigma;
            }
            if (suffix != 2)
                suffix = 2;
            else {
                suffix = follow.A[2];
                follow.A[2] = 3;
                countx.A[suffix] = countx.A[suffix] + 1;
                if (countx.A[suffix] == (suffix - 1)) {
                    countx.A[suffix] = 0;
                    follow.A[suffix - 1] = follow.A[suffix];
                    follow.A[suffix] = suffix + 1;
                }
            }
        }
        getposn = 1;
        putposn = 1;
        elt = 1;
        rho = nolls;
        first = false;
        while ((getposn <= n) && (elt >= 0)) {
            elt = ptn.A[getposn] + perm.A[getposn];
            if (elt > 0)
                if (putposn == 1) {
                    rho.A[1] = elt;
                    putposn = 2;
                } else {
                    altposn = putposn;
                    move = (bool) (elt > rho.A[altposn - 1]);
                    while (move) {
                        rho.A[altposn] = rho.A[altposn - 1];
                        altposn = altposn - 1;
                        if (altposn == 1)
                            move = false;
                        else
                            move = (bool) (elt > rho.A[altposn - 1]);
                    }
                    rho.A[altposn] = elt;
                    putposn = putposn + 1;
                }
            getposn = getposn + 1;
        }
        if (elt >= 0)
            insort (sigma, rho, &head, &tail);
    }
    return head;
}

termptr
sfntohomo (termptr sfnlist)
{
    termptr result, temp, sublist;

    result = NULL;
    temp = sfnlist;
    while (temp != NULL) {
        sublist = factorsof (temp->mult, temp->val);
        merge (&result, &sublist, true, true);
        temp = temp->next;
    }
    sort (&result, true);
    ldisp (&temp);
    return result;
}

termptr
homotosfn (termptr hlist)
{
    termptr result, temp, templist, tlist, sublist, onesfn, onepart;
    int k;
    register int i;

    result = NULL;
    temp = hlist;
    snu (&onesfn);
    onesfn->mult = 1;
    for (i = 1; i <= maxdim; i++) {
        onesfn->val.A[i] = 0;
    }
    while (temp != NULL) {
        k = temp->mult;
        snu (&templist);
        templist->mult = 1;
        for (i = 1; i <= maxdim; i++) {
            templist->val.A[i] = 0;
        }
        i = 1;
        do {
            onepart = onesfn;
            onepart->val.A[1] = temp->val.A[i];
            tlist = louter (templist, onepart);
            ldisp (&templist);
            templist = tlist;
            i = i + 1;
        }
        while (!((temp->val.A[i] == 0)));
        sublist = sfnmult (k, templist);
        ldisp (&templist);
        merge (&result, &sublist, true, true);
        temp = temp->next;
    }
    sort (&result, true);
    ldisp (&temp);
    dispsfn (&onesfn);
    return result;
}

termptr
elemtosfn (termptr elist)
{
    termptr result, temp, templist, tlist, sublist, onesfn, onepart;
    frame tframe;
    int k;
    register int j;
    register int i;

    result = NULL;
    temp = elist;
    snu (&onesfn);
    onesfn->mult = 1;
    for (i = 1; i <= maxdim; i++) {
        onesfn->val.A[i] = 0;
    }
    while (temp != NULL) {
        k = temp->mult;
        snu (&templist);
        templist->mult = 1;
        for (i = 1; i <= maxdim; i++) {
            templist->val.A[i] = 0;
        }
        i = 1;
        do {
            for (j = 1; j <= maxdim; j++) {
                tframe.A[j] = 0;
            }
            tframe.A[1] = temp->val.A[i];
            onepart = onesfn;
            conjgte (&tframe);
            onepart->val = tframe;
            tlist = louter (templist, onepart);
            ldisp (&templist);
            templist = tlist;
            i = i + 1;
        }
        while (!((temp->val.A[i] == 0)));
        sublist = sfnmult (k, templist);
        ldisp (&templist);
        merge (&result, &sublist, true, true);
        temp = temp->next;
    }
    sort (&result, true);
    ldisp (&temp);
    dispsfn (&onesfn);
    return result;
}

termptr
sfntoelem (termptr slist)
{
    termptr temp;

    temp = lconjgte (slist);
    return sfntohomo (temp);
}

termptr
elemtohomo (termptr ehlist)
{
    register termptr temp2;
    termptr temp;

    temp = elemtosfn (ehlist);
    temp2 = sfntohomo (temp);
    ldisp (&temp);
    return temp2;
}

termptr
homotoelem (termptr helist)
{
    register termptr temp2;
    termptr temp;

    temp = homotosfn (helist);
    temp2 = sfntoelem (temp);
    ldisp (&temp);
    return temp2;
}

void
categrze (termptr * list, termarray * catlist)
{
    int wgt;
    termptr trans;

    while ((*list) != NULL) {
        register termptr W80 = &(*(*list));

        wgt = wtfrm (&W80->val);
        trans = W80->next;
        W80->next = catlist->A[wgt];
        catlist->A[wgt] = (*list);
        (*list) = trans;
    }
}

termptr
reduinnprd (frame lambda, frame mu)
{
    register termptr temp;
    termptr alpha, alpha1, alpha2, alpha2a, allbeta, allgamma,
        gamma1, beta1, term1, term1a, term2, term2a,
        term12, term3, subterm12, subprod, prodx;
    termarray beta, gamma;
    frame nu;
    int wtlamu;
    register int i;

    nu = nolls;
    for (i = 0; i <= maxdim; i++) {
        beta.A[i] = NULL;
        gamma.A[i] = NULL;
    }
    snu (&alpha2);
    alpha2->mult = 1;
    alpha2->next = NULL;
    allbeta = skcompat (lambda);
    allgamma = skcompat (mu);
    categrze (&allgamma, &gamma);
    categrze (&allbeta, &beta);
    mergemin (lambda, mu, &nu);
    alpha = skcompat (nu);
    sort (&alpha, true);
    wtlamu = MIN (wtfrm (&lambda), wtfrm (&mu));
    i = 0;
    alpha1 = alpha;
    prodx = NULL;
    while (i <= wtlamu) {
        while (wtfrm (&alpha1->val) > wtlamu - i)
            alpha1 = alpha1->next;
        beta1 = beta.A[i];
        do {
            gamma1 = gamma.A[i];
            term1 = skew (lambda, beta1->val);
            do {
                term12 = NULL;
                alpha2a = alpha1;
                term2 = skew (mu, gamma1->val);
                do {
                    alpha2->val = alpha2a->val;
                    term1a = lskew (term1, alpha2);
                    term2a = lskew (term2, alpha2);
                    subterm12 = louter (term1a, term2a);
                    if (subterm12 != NULL)
                        add (&term12, &subterm12);
                    alpha2a = alpha2a->next;
                    ldisp (&term1a);
                    ldisp (&term2a);
                }
                while (!(alpha2a == NULL));
                term3 = inner (beta1->val, gamma1->val);
                subprod = louter (term12, term3);
                if (subprod != NULL)
                    add (&prodx, &subprod);
                gamma1 = gamma1->next;
                ldisp (&term2);
                ldisp (&term12);
                ldisp (&term3);
            }
            while (!(gamma1 == NULL));
            beta1 = beta1->next;
            ldisp (&term1);
        }
        while (!(beta1 == NULL));
        i = i + 1;
    }
    sort (&prodx, true);
    temp = prodx;
    dispsfn (&alpha2);
    ldisp (&alpha);
    for (i = 0; i <= maxdim; i++) {
        ldisp (&beta.A[i]);
        ldisp (&gamma.A[i]);
    }
    return temp;
}

termptr
rinner (termptr list1, termptr list2)
{
    termptr templist, sublist, prodlist, x;
    int multy;

    prodlist = NULL;
    while (list1 != NULL) {
        templist = list2;
        while (templist != NULL) {
            sublist = reduinnprd (list1->val, templist->val);
            multy = list1->mult * templist->mult;
            x = sublist;
            if (multy != 1)
                while (x != NULL) {
                    x->mult = x->mult * multy;
                    x = x->next;
                }
            add (&prodlist, &sublist);
            templist = templist->next;
        }
        list1 = list1->next;
    }
    sort (&prodlist, true);
    return prodlist;
}

termptr
plethonerinner (termptr list)
{
    int molt;
    register int i, j;
    termptr dummy, sublist, factlist, result, temp0, temp1, temp2;

    result = NULL;
    temp2 = NULL;
    sublist = lconjgte (list);
    while (sublist != NULL) {
        factlist = factorsof (sublist->mult, sublist->val);
        merge (&result, &factlist, true, true);
        sublist = sublist->next;
    }
    ldisp (&sublist);
    factlist = result;
    while (factlist != NULL) {

        snu (&temp0);
        temp0->mult = 1;
        i = 1;
        for (j = 1; j <= maxdim; j++) {
            temp0->val.A[j] = 0;
        }
        molt = factlist->mult;
        while (factlist->val.A[i] != 0) {
            snu (&temp1);
            temp1->mult = 1;
            for (j = 1; j <= maxdim; j++) {
                temp1->val.A[j] = 0;
            }
            for (j = 1; j <= factlist->val.A[i]; j++) {
                temp1->val.A[j] = 1;
            }
            dummy = temp0;
            temp0 = rinner (dummy, temp1);
            dispsfn (&temp1);
            ldisp (&dummy);
            i = i + 1;
        }
        dummy = temp0;
        temp0 = sfnmult (molt, dummy);
        ldisp (&dummy);
        dummy = temp2;
        temp2 = ladd (temp0, dummy);
        ldisp (&dummy);
        ldisp (&temp0);
        factlist = factlist->next;
    }
    ldisp (&result);
    return temp2;
}

termptr
makeweight (int n, termptr list)
{
    termptr newlist, temp = NULL,       //modified by FB, was uninitialized
        partmw;
    int w, i;

    newlist = NULL;
    snu (&partmw);
    while (list != NULL) {
        if (newlist == NULL) {
            snu (&newlist);
            temp = newlist;
        } else {
            snu (&temp->next);
            temp = temp->next;
        }
        (*temp) = (*list);
        (*partmw) = (*temp);
        w = n - wtfrm (&temp->val);
        for (i = 2; i <= maxdim; i++) {
            temp->val.A[i] = partmw->val.A[i - 1];
        }
        temp->val.A[1] = w;
        list = list->next;
    }
    if ((qfn == true))
        qstndise (&newlist);
    else
        stndise (&newlist);
    ldisp (&partmw);
    sort (&newlist, false);
    return newlist;
}

termptr
inner (frame lambda, frame mu)
{
    int j, wt, prwt;
    register int i;
    termptr prodx, top;
    frame lambda1i, mu1;
    bool conjx;

    lambda1i = nolls;
    mu1 = nolls;
    wt = wtfrm (&lambda);
    if (wt == wtfrm (&mu)) {
        i = 1;
        while (lambda.A[i] != 0)
            i = i + 1;
        j = 1;
        while (mu.A[j] != 0)
            j = j + 1;
        if (i - 1 > lambda.A[1]) {
            conjx = true;
            conjgte (&lambda);
        } else
            conjx = false;
        if (j - 1 > mu.A[1]) {
            conjx = (bool) (!conjx);
            conjgte (&mu);
        }
        if ((lambda.A[2] == 0) || (mu.A[2] == 0)) {
            snu (&top);
            top->mult = 1;
            top->next = NULL;
            if (lambda.A[2] == 0)
                top->val = mu;
            else
                top->val = lambda;
        } else {
            for (i = 1; i <= maxdim - 1; i++) {
                lambda1i.A[i] = lambda.A[i + 1];
                mu1.A[i] = mu.A[i + 1];
            }
            lambda1i.A[maxdim] = 0;
            mu1.A[maxdim] = 0;
            lambda1i.A[0] = 0;
            mu1.A[0] = 0;
            lambda1i.length = lambda.length - 1;
            mu1.length = mu.length - 1;
            prodx = reduinnprd (lambda1i, mu1);
            top = prodx;
            while (prodx != NULL) {
                prwt = wtfrm (&prodx->val);

                for (i = prwt; i >= 1; i--) {
                    prodx->val.A[i + 1] = prodx->val.A[i];
                }
                prodx->val.A[1] = wt - prwt;
                prodx = prodx->next;
            }
            stndise (&top);
            sort (&top, true);
        }
        if (conjx)
            return lconjgte (top);
        else
            return top;
    }                           //if wt==wtfrm(&mu)
    else
        return NULL;
}

termptr
linner (termptr list1, termptr list2)
{
    termptr templist, sublist, prodlist, x;
    int multy;

    prodlist = NULL;
    while (list1 != NULL) {
        templist = list2;
        while (templist != NULL) {
            sublist = inner (list1->val, templist->val);
            multy = list1->mult * templist->mult;
            x = sublist;
            if (multy != 1)
                while (x != NULL) {
                    x->mult = x->mult * multy;
                    x = x->next;
                }
            add (&prodlist, &sublist);
            templist = templist->next;
        }
        list1 = list1->next;
        sort (&prodlist, true);
    }
    return prodlist;
}

void
putchrc (text * fyle, ocharptr chrc, bool vdu)
{
    int qq;
    char a = '{', b = '}';      // was uninitialized. FB
    //caystype caes;
    if (!iosup) {
        /*if (ggroup.name == e6)
           caes = pe6;
           else
           caes = norm; */
        switch ((int) (ggroup.name)) {
        case sung:
        case un:
        case sn:
        case an:
        case unm:
        case sunm:
        case unc:
            a = '{';
            b = '}';
            break;
        case spnc:
        case spn:
        case mp:
            a = '<';
            b = '>';
            break;
        case ospnm:
            a = '[';
            b = '>';
            break;
        case son:
        case on:
        case nill:
        case sonc:
            a = '[';
            b = ']';
            break;
        case g2:
        case f4:
        case en:
        case e6:
        case e7:
        case e8:
        case l168:
            a = '(';
            b = ')';
            break;
        default:
            Caseerror (Line);
        }
        qq = 1;
        writer2 (&(*fyle), &qq, a, b, chrc, vdu);
        Putchr ('\n', (*fyle));
        if (logging) {
            qq = 1;
            writer2 (&logfile, &qq, a, b, chrc, false);
            Putchr ('\n', logfile);
        }
    }
}

void
putsfn (text * fyle, termptr sfn, bool vdu)
{
    int qq;

    if (!iosup) {
        qq = 1;
        if ((nreduce == true))
            wrttlst2 (&(*fyle), &qq, '<', '>', sfn, vdu);
        else
            wrttlst2 (&(*fyle), &qq, '{', '}', sfn, vdu);
        Putchr ('\n', (*fyle));
        if (logging) {
            qq = 1;
            if ((nreduce == true))
                wrttlst2 (&logfile, &qq, '<', '>', sfn, false);
            else
                wrttlst2 (&logfile, &qq, '{', '}', sfn, false);
            Putchr ('\n', logfile);
        }
    }
}

void
putprod (prodtype pr)
{
    int qq = 1;

    if (!iosup) {
        //writeprod (&output, pr, tcol, true);
        writeprod2 (&output, &qq, pr, true);
        if (logging)
            //writeprod (&logfile, pr, tcol, false);
            writeprod2 (&logfile, &qq, pr, false);
        print ("\n");
    }
}

void
putgroup (groopArray grp)
{
    register int j;

    if (echo) {
        if (nprod > 1)
            inform ("Groups are ;", cont);
        else
            inform ("Group is ;", cont);
        for (j = 1; j <= nprod; j++) {
            if ((nprod > 1) && (j <= nprod) && (j > 1))
                inform (" * ;", cont);
            switch ((int) (grp.A[j - 1].name)) {
            case sung:
                inform ("SU(;", cont);
                break;
            case un:
                inform ("U(;", cont);
                break;
            case unm:
                inform ("U(;", cont);
                break;
            case unc:
                inform ("U(;", cont);
                break;
            case sunm:
                inform ("SU(;", cont);
                break;
            case son:
                inform ("SO(;", cont);
                break;
            case on:
                inform ("O(;", cont);
                break;
            case spn:
                inform ("Sp(;", cont);
                break;
            case ospnm:
                inform ("OSp(;", cont);
                break;
            case spnc:
                inform ("Sp(;", cont);
                break;
            case sonc:
                inform ("SO^*(;", cont);
                break;
            case mp:
                inform ("Mp(;", cont);
                break;
            case l168:
                inform ("L(;", cont);
                break;
            case sn:
                inform ("S(;", cont);
                if (((bool) ((grp.A[j - 1].rank) & 1)))
                    qsn = 1;
                else
                    qsn = 0;
                break;
            case an:
                inform ("A(;", cont);
                if (((bool) ((grp.A[j - 1].rank) & 1)))
                    qsn = 1;
                else
                    qsn = 0;
                break;
            case e6:
            case e7:
            case e8:
            case en:
                inform ("E(;", cont);
                break;
            case f4:
                inform ("F(;", cont);
                break;
            case g2:
                inform ("G(;", cont);
                break;
            case nill:
                inform ("No group set;", cr);
                break;
            default:
                Caseerror (Line);
            }
            if (grp.A[j - 1].name != nill) {
                print ("%1d", grp.A[j - 1].rank);
                switch ((int) (grp.A[j - 1].name)) {
                case ospnm:
                case unm:
                case sunm:
                    print ("/%1d", grp.A[j - 1].rank2);
                    break;
                case unc:
                    print (",%1d", grp.A[j - 1].rank2);
                    break;
                case spnc:
                    print (",R");
                    break;
                case un:
                case on:
                case son:
                case spn:
                case sung:
                case mp:
                case an:
                case sn:
                case g2:
                case f4:
                case e6:
                case e7:
                case en:
                case e8:
                case nill:
                case sonc:
                case l168:;
                    break;
                default:
                    Caseerror (Line);
                }
                inform (");", cont);
            }
        }
        inform (";", cr);
    }
}



void
getgroup (string0 * buffx, int p, bool msgs)
{
    bool errorx;
    int j;

    readint (&(*buffx), &p, &nprod);
    if (nprod == 0)
        nprod = 1;
    j = 1;
    errorx = false;
    if (nprod > maxprod) {
        if (msgs)
            inform ("Maximum number of groups exceeded;", cr);
        errorx = true;
    }
    while ((!errorx) && (j <= nprod)) {
        //register groop *W101 = &currgrp.A[j - 1];

        errorx = false;
        while ((buffx->A[p - 1] == ' ') && (p < bcol))
            p = p + 1;
        if (p < bcol)
            if (((buffx->A[p - 1] == 's') || (buffx->A[p - 1] == 'S'))) {
                p = p + 1;
                if (((buffx->A[p - 1] == 'u') || (buffx->A[p - 1] == 'U'))) {
                    p = p + 1;
                    if (((buffx->A[p - 1] == 'm')
                         || (buffx->A[p - 1] == 'M'))) {
                        currgrp.A[j - 1].name = sunm;
                        p = p + 2;
                    } else
                        currgrp.A[j - 1].name = sung;
                } else
                    if (((buffx->A[p - 1] == 'o')
                         || (buffx->A[p - 1] == 'O'))) {
                    currgrp.A[j - 1].name = son;
                    p = p + 1;
                    if (((buffx->A[p - 1] == 'n')
                         || (buffx->A[p - 1] == 'N'))) {
                        currgrp.A[j - 1].name = sonc;
                        p = p + 1;
                    }
                } else
                    if (((buffx->A[p - 1] == 'p')
                         || (buffx->A[p - 1] == 'P'))) {
                    p = p + 1;
                    if (((buffx->A[p - 1] == 'r')
                         || (buffx->A[p - 1] == 'R'))) {
                        currgrp.A[j - 1].name = spnc;
                        p = p + 1;
                    } else
                        currgrp.A[j - 1].name = spn;
                } else
                    currgrp.A[j - 1].name = sn;
            } else if (((buffx->A[p - 1] == 'u') || (buffx->A[p - 1] == 'U'))) {
                p = p + 1;
                if (((buffx->A[p - 1] == 'm') || (buffx->A[p - 1] == 'M'))) {
                    currgrp.A[j - 1].name = unm;
                    p = p + 2;
                } else
                    if (((buffx->A[p - 1] == 'n')
                         || (buffx->A[p - 1] == 'N'))) {
                    p = p + 1;
                    if (((buffx->A[p - 1] == 'c')
                         || (buffx->A[p - 1] == 'C'))) {
                        currgrp.A[j - 1].name = unc;
                        p = p + 1;
                    }
                } else
                    currgrp.A[j - 1].name = un;
            } else if (((buffx->A[p - 1] == 'o') || (buffx->A[p - 1] == 'O'))) {
                p = p + 1;
                if (((buffx->A[p - 1] == 's') || (buffx->A[p - 1] == 'S'))) {
                    currgrp.A[j - 1].name = ospnm;
                    p = p + 2;
                } else
                    currgrp.A[j - 1].name = on;
            } else if (((buffx->A[p - 1] == 'a') || (buffx->A[p - 1] == 'A'))) {
                currgrp.A[j - 1].name = an;
                p = p + 1;
            } else if (((buffx->A[p - 1] == 'g') || (buffx->A[p - 1] == 'G'))) {
                currgrp.A[j - 1].name = g2;
                p = p + 1;
            } else if (((buffx->A[p - 1] == 'f') || (buffx->A[p - 1] == 'F'))) {
                currgrp.A[j - 1].name = f4;
                p = p + 1;
            } else if (((buffx->A[p - 1] == 'e') || (buffx->A[p - 1] == 'E'))) {
                currgrp.A[j - 1].name = en;
                p = p + 1;
            } else if ((buffx->A[p - 1] == 'm') || (buffx->A[p - 1] == 'M')) {
                currgrp.A[j - 1].name = mp;
                p = p + 2;
            } else if ((buffx->A[p - 1] == 'l') || (buffx->A[p - 1] == 'L')) {
                currgrp.A[j - 1].name = l168;
                currgrp.A[j - 1].rank = 168;
                p = p + 1;
            } else
                errorx = true;
        else
            errorx = true;
        /*if (!errorx && (p < bcol) && (Member((unsigned)(buffx->A[p - 1]), Conset[0]))) {
           readint(&(*buffx), &p, &W101->rank);
           switch ((int)(W101->name)) {
           case ospnm:  case unm:  case sunm: case unc:
           readint(&(*buffx), &p, &W101->rank2);
           break ;
           case un: case on: case son: case spn: case sung: case mp:
           case an: case sn: case g2: case f4: case e6: case e7: case en:
           case e8: case spnc: case nill: case sonc: case l168:;
           break;
           default:
           Caseerror(Line);
           } */
        if (!errorx && (p < bcol)) {
            switch ((int) (buffx->A[p - 1])) {
            case ' ':
            case '0':
            case '1':
            case '2':
            case '3':
            case '4':
            case '5':
            case '6':
            case '7':
            case '8':
            case '9':
                readint (&(*buffx), &p, &currgrp.A[j - 1].rank);
                switch ((int) (currgrp.A[j - 1].name)) {
                case ospnm:
                case unm:
                case sunm:
                case unc:
                    readint (&(*buffx), &p, &currgrp.A[j - 1].rank2);
                    break;
                case un:
                case on:
                case son:
                case spn:
                case sung:
                case mp:
                case an:
                case sn:
                case g2:
                case f4:
                case e6:
                case e7:
                case en:
                case e8:
                case spnc:
                case nill:
                case sonc:
                case l168:;
                    break;
                default:
                    Caseerror (Line);
                }
            }
            if ((currgrp.A[j - 1].name == en))
                if (((currgrp.A[j - 1].rank < 6)
                     || (currgrp.A[j - 1].rank > 8)))
                    errorx = true;
                else
                    switch ((int) (currgrp.A[j - 1].rank)) {
                    case 6:
                        currgrp.A[j - 1].name = e6;
                        break;
                    case 7:
                        currgrp.A[j - 1].name = e7;
                        break;
                    case 8:
                        currgrp.A[j - 1].name = e8;
                        break;
                    default:
                        Caseerror (Line);
            } else if ((currgrp.A[j - 1].name == g2)
                       && (currgrp.A[j - 1].rank != 2))
                errorx = true;
            else if ((currgrp.A[j - 1].name == f4)
                     && (currgrp.A[j - 1].rank != 4))
                errorx = true;
            else if (currgrp.A[j - 1].rank <= 0)
                errorx = true;
            else if ((bool) ((currgrp.A[j - 1].rank) & 1)
                     && ((currgrp.A[j - 1].name == sonc)
                         || (currgrp.A[j - 1].name == spn)
                         || (currgrp.A[j - 1].name == spnc)
                         || (currgrp.A[j - 1].name == mp)))
                errorx = true;
        } else
            errorx = true;
        j = j + 1;
    }
    if (msgs)
        if (errorx) {
            warn ("invalid group:gr", cont);
            inform ("oup not set;", cr);
        } else
            putgroup (currgrp);
    erred = errorx;
}


ocharptr
clabel (ocharptr list)
{
    ocharptr newlist, temp = NULL;      //modified by FB, was uninitialized

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
        if (temp->lab == '-')
            temp->lab = '+';
        else if (temp->lab == '+')
            temp->lab = '-';
        list = list->next;
    }
    osort (&newlist, false);
    return newlist;
}

ocharptr
oformbb (ocharptr con, ocharptr cov, bool suppress, bool spinor)
{
    ocharptr lastptr, list;
    ocharptr contr;

    lastptr = NULL;
    while (cov != NULL) {
        contr = con;
        while (contr != NULL) {
            cnu (&list);
            {
                if ((contr->val.A[1] == 0) && suppress)
                    list->C6_double = false;
                else {
                    list->C6_double = true;
                    list->conlab = contr->lab;
                    list->conval = contr->val;
                }
                list->val = cov->val;
                list->mult = contr->mult * cov->mult;
                list->spin = spinor;
                list->lab = cov->lab;
                list->next = lastptr;
                lastptr = list;
                contr = contr->next;
            }
        }
        cov = cov->next;
    }
    return lastptr;
}

prodtype
pwrestrict (prodtype list, int wt, int pno, char ch)
{
    prodtype newlist, temp, tlist;
    int weightx = 0,            //modified by FB, was uninitialized
        tw;
    register int i;
    bool okx;

    newlist = NULL;
    ch = locase (ch);
    tw = abs (wt);
    while (list != NULL) {
        pnu (&temp);
        temp->mult = list->mult;
        for (i = 1; i <= nprod; i++) {

            temp->prods.A[i - 1] = chrccopy (list->prods.A[i - 1]);
        }
        temp->next = NULL;

        if ((ch == 'w'))
            weightx = wtfrm (&temp->prods.A[pno - 1]->val);
        else if ((ch == 'l'))
            weightx = len (&temp->prods.A[pno - 1]->val);

        if (((wt >= 0) && (weightx <= wt)) || ((wt < 0) && (weightx == tw)))
            okx = true;
        else
            okx = false;

        if (okx == true) {
            tlist = prodadd (newlist, temp);

            pdisp (&newlist);
            newlist = tlist;

        }
        dispprod (&temp);
        list = list->next;
    }
    schur_psort (&newlist, true);
    return newlist;
}

termptr
complement (termptr list, int n, int m)
{
    termptr newlist, temp;
    int k;
    register int i;
    bool test;
    register termptr W105;

    newlist = NULL;
    temp = NULL;
    test = true;
    while (((list != NULL) && (test == true))) {
        W105 = list;

        if (len (&W105->val) <= n) {
            if (newlist == NULL) {
                snu (&newlist);
                temp = newlist;
            } else {
                snu (&temp->next);
                temp = temp->next;
            }
            (*temp) = (*list);
            for (i = 0; i <= n - 1; i++) {
                k = n - i;
                //temp->val.A[k] = W105->val.A[1] - W105->val.A[i + 1];
                temp->val.A[k] = m - W105->val.A[i + 1];
                if (temp->val.A[k] < 0) {
                    print ("m is too small\n");
                    return NULL;
                }
            }
            list = list->next;
        } else {
            test = false;
            print ("ERROR: Partitions must be of length <= %d\n", n);
            return NULL;
        }
    }
    sort (&newlist, false);
    return newlist;
}

void
partsrep (ocharptr * list, int minn, int maxm)
{
    register int i;
    for (i = minn; i <= maxm; i++) {
        unlimit (list, i);
    }
}

void
partssfn (termptr * list, int minn, int maxm)
{
    register int i;
    for (i = minn; i <= maxm; i++) {
        limit (list, i);
    }
}

termptr
sumsquares (termptr list)
{
    termptr newlist, tlist, plist;

    newlist = NULL;
    while (list != NULL) {
        snu (&tlist);
        tlist->mult = list->mult;
        tlist->val = list->val;
        plist = louter (tlist, tlist);
        dispsfn (&tlist);
        add (&newlist, &plist);
        sort (&newlist, true);
        list = list->next;
    }

    return newlist;
}
