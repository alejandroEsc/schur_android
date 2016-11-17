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
#include "s3.h"
#include "s4.h"
#include "s5.h"
#include "s7.h"
#include "s6.h"
#include "s.h"
#include "m.h"
#include "branch.h"
#include "skew.h"
#include "gr.h"

ocharptr
sp2nk (ocharptr * list, int n, int k)
{
    ocharptr newlist, temp, con, cov, tempt;
    termptr slist, tlist;
    int pp, t, nn;
    bool test;
    tempt = *list;
    newlist = NULL;
    while (tempt != NULL) {
        if (((tempt->val.A[1] == 1) || (tempt->val.A[1] == 0))
            && (len (&tempt->val) <= 1)) {
            test = false;
            nn = n / 2;
            if (((k % 2) == 1))
                test = true;
            pp = k / 2;
            slist = spseries (nn, k);
            if (tempt->val.A[1] == 0)
                schur_restrict (&slist, 1, 'e');
            else if (tempt->val.A[1] == 1)
                schur_restrict (&slist, 1, 'o');
            tlist = slist;
            while (tlist != NULL) {
                register termptr W2 = &(*tlist);

                t = len (&W2->val);
                if (t <= nn) {
                    cnu (&con);
                    con->mult = 1;
                    con->val = W2->val;
                    con->val.A[maxdim] = pp;
                    con->lab = ' ';
                    cnu (&cov);
                    cov->mult = 1;
                    cov->val = W2->val;
                    cov->lab = ' ';
                    temp = oformbb (con, cov, false, test);
                    dispchr (&con);
                    dispchr (&cov);
                    oadd (&newlist, &temp);
                }
                tlist = W2->next;
            }
            ldisp (&slist);
        }
        tempt = tempt->next;
    }

    return newlist;
}

termptr
sp2spun (termptr * list, int k, int n, int m)
{
    termptr temp, mtemp, ntemp;
    int g;

    g = n / 2;
    temp = spun (&(*list), k, g);
    schur_restrict (&temp, k, 'l');
    schur_restrict (&temp, m, 'w');
    snu (&mtemp);
    mtemp->mult = 1;
    mtemp->val = nolls;
    mtemp->val.A[1] = m;
    ntemp = useseries ('c', mtemp, false, true, n);
    dispsfn (&mtemp);
    mtemp = linner (temp, ntemp);
    ldisp (&ntemp);
    ldisp (&temp);
    temp = useseries ('d', mtemp, false, true, n);
    ldisp (&mtemp);
    schur_restrict (&temp, g, 'l');
    return temp;
}

termptr
spun (termptr * list, int k, int n)
{
    termptr newlist, temp, list1, list2;
    int m, w1;
    char tag;
    bool ppp;

    newlist = NULL;
    m = MIN (n, k);
    ppp = false;
    if ((bool) ((k) & 1))
        tag = 'g';
    else
        tag = 'c';
    w1 = wtfrm (&(*list)->val);
    if (w1 >= 1)
        w1 = w1 + setlimit - 1;
    else
        w1 = setlimit - 2;
    list1 = (termptr) signseq (&(*list), k, m, tag, ppp, false);
    list2 = rseries (m, 'd');
    temp = louter2 (list1, list2, m);
    schur_restrict (&temp, w1, 'w');
    ldisp (&list1);
    ldisp (&list2);
    add (&newlist, &temp);
    sort (&newlist, true);
    return newlist;
}

ocharptr
metaplet (int n, int k)
{
    ocharptr newlist, temp, con, cov;
    termptr slist, tlist;
    int pp, t, nn;
    bool test;

    newlist = NULL;
    test = false;
    nn = n / 2;
    if (((k % 2) == 1))
        test = true;
    pp = k / 2;
    slist = spseries (nn, k);
    tlist = slist;
    while (tlist != NULL) {
        register termptr W3 = &(*tlist);

        t = len (&W3->val);
        if (t <= nn) {
            cnu (&con);
            con->mult = 1;
            con->val = W3->val;
            con->val.A[maxdim] = pp;
            con->lab = ' ';
            cnu (&cov);
            cov->mult = 1;
            cov->val = W3->val;
            cov->lab = ' ';
            temp = oformbb (con, cov, false, test);
            dispchr (&con);
            dispchr (&cov);
            oadd (&newlist, &temp);
        }
        tlist = W3->next;
    }
    ldisp (&slist);
    return newlist;
}

ocharptr
lspnrsp2on (ocharptr * llist, int n)
{
    ocharptr newlist, temp, ttlist;

    newlist = NULL;
    ttlist = (*llist);
    while (ttlist != NULL) {
        register ocharptr W4 = &(*ttlist);

        temp = spnrsp2on (&ttlist, n);
        oadd (&newlist, &temp);
        ttlist = W4->next;
    }
    osort (&newlist, true);
    return newlist;
}

ocharptr
spnrsp2on (ocharptr * list, int n)
{
    int p, k;
    register int m;
    ocharptr newlist, flist, slist, tlist, olist;
    termptr temp, sflist;
    bool test;
    int mm, nn, mmin, r;
    frame ff;

    test = false;
    ff = (*list)->val;
    mmin = wtfrm (&ff);
    if (mmin == 0)
        r = setlimit - 1;
    else
        r = setlimit;
    mm = (*list)->mult;
    nn = n / 2;
    newlist = NULL;
    k = (*list)->conval.A[1];
    if ((*list)->spin)
        k = 2 * k + 1;
    else
        k = 2 * k;
    p = nn * k;
    if ((p % 2 == 1))
        test = true;
    p = p / 2;
    for (m = mmin; m <= r; m++) {
        snu (&sflist);
        sflist->val = ff;
        sflist->mult = mm;
        cnu (&flist);
        flist->spin = test;
        flist->mult = 1;
        flist->val = nolls;
        flist->val.A[maxdim] = p;
        flist->val.A[1] = m;
        temp = sp2spun (&sflist, k, n, m);
        slist = sfntochrc (temp, false, ' ');
        olist = gmodify (slist, currgrp.A[2 - 1]);
        ldisp (&temp);
        tlist = oformbb (flist, olist, false, test);
        odisp (&slist);
        odisp (&olist);
        dispsfn (&sflist);
        odisp (&flist);
        if (tlist != NULL)
            oadd (&newlist, &tlist);
    }
    osort (&newlist, true);
    return newlist;
}

ocharptr
sprun (ocharptr * list, int n)
{
    int i, k, t, r;
    register int j;
    ocharptr temp;
    termptr plist, slist, eta, nlist;
    temp = NULL;
    t = n / 2;
    i = (*list)->val.A[maxdim];
    k = i;
    if ((*list)->spin)
        k = 2 * k + 1;
    else
        k = 2 * k;
    snu (&eta);
    eta->mult = 1;
    eta->val = nolls;
    for (j = 1; j <= t; j++) {
        eta->val.A[j] = i;
    }
    snu (&nlist);
    nlist->mult = (*list)->mult;
    nlist->val = (*list)->val;
    r = setlimit;
    if (wtfrm (&(*list)->val) == 0)
        r = r - 1;
    slist = spun (&nlist, k, t);
    schur_restrict (&slist, r, 'w');
    plist = louter2 (eta, slist, t);
    ldisp (&slist);
    dispsfn (&nlist);
    temp = sfntochrc (plist, (*list)->spin, ' ');
    ldisp (&plist);

    dispsfn (&eta);
    return temp;
}

ocharptr
lsprun (ocharptr * llist, int n)
{
    ocharptr newlist, temp, tlist;

    newlist = NULL;
    tlist = (*llist);
    while (tlist != NULL) {
        register ocharptr W9 = &(*tlist);

        temp = sprun (&tlist, n);
        oadd (&newlist, &temp);
        tlist = W9->next;
    }

    osort (&newlist, true);
    return newlist;
}

ocharptr
so4so3 (ocharptr irrep)
{
    int ssign;
    register int i;
    ocharptr newlist, temp, cu, co;

    newlist = NULL;
    while (irrep != NULL) {
        register ocharptr W10 = &(*irrep);

        cnu (&co);
        cnu (&cu);
        co->mult = W10->mult;
        cu->mult = 1;
        co->lab = ' ';
        cu->lab = ' ';
        for (i = 1; i <= maxdim; i++) {
            co->val.A[i] = 0;
            cu->val.A[i] = 0;
        }
        if (W10->lab == '-')
            ssign = -1;
        else
            ssign = 1;
        co->val.A[1] = W10->val.A[1] + ssign * W10->val.A[2];
        cu->val.A[1] = W10->val.A[1] - ssign * W10->val.A[2];
        if (W10->spin)
            if (W10->lab == '+')
                co->val.A[1] = co->val.A[1] + 1;
            else
                cu->val.A[1] = cu->val.A[1] + 1;
        if (((cu->val.A[1] % 2) == 1))
            cu->spin = true;
        else
            cu->spin = false;
        if (((co->val.A[1] % 2) == 1))
            co->spin = true;
        else
            co->spin = false;
        cu->val.A[1] = cu->val.A[1] / 2;
        co->val.A[1] = co->val.A[1] / 2;
        temp = kronk (cu, co, currgrp.A[nprod - 1]);    /*15/1/97 */
        odisp (&co);
        odisp (&cu);
        oadd (&newlist, &temp);
        irrep = W10->next;
        osort (&newlist, false);
    }
    return newlist;
}

ocharptr
so4su2su2 (ocharptr irrep)
{
    int ssign;
    register int i;
    ocharptr newlist, temp, cu, co;

    newlist = NULL;
    while (irrep != NULL) {
        register ocharptr W13 = &(*irrep);

        cnu (&co);
        cnu (&cu);
        co->mult = W13->mult;
        cu->mult = 1;
        co->lab = ' ';
        cu->lab = ' ';
        for (i = 1; i <= maxdim; i++) {
            co->val.A[i] = 0;
            cu->val.A[i] = 0;
        }
        if (W13->lab == '-')
            ssign = -1;
        else
            ssign = 1;
        co->val.A[1] = W13->val.A[1] + ssign * W13->val.A[2];
        cu->val.A[1] = W13->val.A[1] - ssign * W13->val.A[2];
        if (W13->spin)
            if (W13->lab == '+')
                co->val.A[1] = co->val.A[1] + 1;
            else
                cu->val.A[1] = cu->val.A[1] + 1;
        temp = oformbb (co, cu, false, false);
        odisp (&co);
        odisp (&cu);
        oadd (&newlist, &temp);
        irrep = W13->next;
        osort (&newlist, false);
    }
    return newlist;
}

ocharptr
sp2nsu2son (ocharptr irrep)
{
    ocharptr list, newlist, cu, co, tco;
    int ssm;
    register int i;
    termptr zeta, tzeta, sublist, undum, dum, sub1, ss, unnit, ort, dummy;

    list = NULL;
    snu (&undum);
    undum->next = NULL;
    sublist = NULL;
    while (irrep != NULL) {
        register ocharptr W16 = &(*irrep);

        ss = seriesx ('a', -1, W16->val);
        undum->val = W16->val;
        undum->mult = W16->mult;
        sub1 = lskew (undum, ss);
        ldisp (&ss);
        dummy = ladd (sublist, sub1);
        ldisp (&sub1);
        ldisp (&sublist);
        sublist = dummy;
        irrep = W16->next;
    }
    dispsfn (&undum);
    snu (&dum);
    for (i = 0; i <= maxdim; i++) {
        dum->val.A[i] = 0;
    }
    dum->mult = 1;
    dum->next = NULL;
    while (sublist != NULL) {
        register termptr W19 = &(*sublist);

        zeta = eqwt (wtfrm (&W19->val));
        ssm = W19->mult;
        while (zeta != NULL) {
            tzeta = zeta;
            zeta = zeta->next;
            tzeta->next = NULL;
            tzeta->mult = 1;
            unnit = inner (W19->val, tzeta->val);
            ss = seriesx ('d', -1, tzeta->val);
            ort = lskew (tzeta, ss);
            ldisp (&ss);
            cu = formbb (dum, unnit, true, false, ' ');
            gumodify (&cu, 2);
            spclise (&cu, 2);
            tco = formbb (dum, ort, true, false, ' ');
            co = gmodify (tco, currgrp.A[2 - 1]);
            odisp (&tco);
            newlist = oformbb (cu, co, false, false);
            odisp (&cu);
            odisp (&co);
            cu = chrcmult (ssm, newlist);
            odisp (&newlist);
            oadd (&list, &cu);
            ldisp (&unnit);
            ldisp (&ort);
            dispsfn (&tzeta);
        }
        sublist = W19->next;
    }
    dispsfn (&dum);
    ldisp (&dummy);
    osort (&list, false);
    return list;
}

void
fixsubgroup (int *brno, int *n, int *m, int *n2, int *m2)
{
    int j;
    if ((((*brno) >= 1) && ((*brno) <= 63)))
        switch ((int) ((*brno))) {
        case 1:
            fixcg (&currgrp.A[nprod - 1], on, (*n), 0);
            break;
        case 2:
            fixcg (&currgrp.A[nprod - 1], spn, (*n), 0);
            break;
        case 3:
            (*n) = (*n) - 1;
            fixcg (&currgrp.A[nprod - 1], un, (*n), 0);
            break;
        case 4:
        case 5:
            fixcg (&currgrp.A[nprod - 1], un, (*n), 0);
            nprod = nprod + 1;
            fixcg (&currgrp.A[nprod - 1], un, (*m), 0);
            break;
        case 6:
        case 13:
        case 21:
            j = (*n) / 2;
            fixcg (&currgrp.A[nprod - 1], un, j, 0);
            break;
        case 7:
        case 9:
        case 25:
        case 29:
        case 30:
        case 41:
            fixcg (&currgrp.A[nprod - 1], son, 3, 0);
            break;
        case 8:
            fixcg (&currgrp.A[nprod - 1], un, 1, 0);
            fixcg (&currgrp.A[nprod + 1 - 1], sung, (*n), 0);
            fixcg (&currgrp.A[nprod + 2 - 1], sung, (*m), 0);
            nprod = nprod + 2;
            break;
        case 10:
            j = (*n) / 2;
            fixcg (&currgrp.A[nprod - 1], un, 1, 0);
            fixcg (&currgrp.A[nprod + 1 - 1], sung, j, 0);
            nprod = nprod + 1;
            break;
        case 11:
            fixcg (&currgrp.A[nprod - 1], sung, 2, 0);
            j = (*n) / 2;
            fixcg (&currgrp.A[nprod + 1 - 1], son, j, 0);
            nprod = nprod + 1;
            break;
        case 12:
        case 20:
            fixcg (&currgrp.A[nprod - 1], sung, (*n), 0);
            break;
        case 14:
            fixcg (&currgrp.A[nprod - 1], spn, (*n), 0);
            nprod = nprod + 1;
            fixcg (&currgrp.A[nprod - 1], spn, (*m), 0);
            break;
        case 15:
            fixcg (&currgrp.A[nprod - 1], spn, (*n), 0);
            nprod = nprod + 1;
            fixcg (&currgrp.A[nprod - 1], on, (*m), 0);
            break;
        case 16:
            fixcg (&currgrp.A[nprod - 1], sn, (*n), 0);
            nprod = nprod + 1;
            fixcg (&currgrp.A[nprod - 1], sn, (*m), 0);
            break;
        case 17:
            fixcg (&currgrp.A[nprod - 1], an, (*n), 0);
            break;
        case 18:
        case 62:
            fixcg (&currgrp.A[nprod - 1], sn, (*n), 0);
            break;
        case 19:
            fixcg (&currgrp.A[nprod - 1], sn, (*n) + 1, 0);
            break;
        case 22:
            fixcg (&currgrp.A[nprod - 1], on, (*n), 0);
            nprod = nprod + 1;
            fixcg (&currgrp.A[nprod - 1], on, (*m), 0);
            break;
        case 23:
            fixcg (&currgrp.A[nprod - 1], on, (*n), 0);
            nprod = nprod + 1;
            fixcg (&currgrp.A[nprod - 1], on, (*m), 0);
            break;
        case 24:
            fixcg (&currgrp.A[nprod - 1], spn, (*n), 0);
            nprod = nprod + 1;
            fixcg (&currgrp.A[nprod - 1], spn, (*m), 0);
            break;
        case 26:
            j = (*n) / 2;
            fixcg (&currgrp.A[nprod - 1], un, 1, 0);
            fixcg (&currgrp.A[nprod + 1 - 1], sung, j, 0);
            nprod = nprod + 1;
            break;
        case 27:
            fixcg (&currgrp.A[nprod - 1], son, (*n), 0);
            nprod = nprod + 1;
            fixcg (&currgrp.A[nprod - 1], son, (*m), 0);
            break;
        case 28:
            fixcg (&currgrp.A[nprod - 1], sung, 2, 0);
            nprod = nprod + 1;
            fixcg (&currgrp.A[nprod - 1], sung, 2, 0);
            break;
        case 31:
        case 46:
            fixcg (&currgrp.A[nprod - 1], g2, 2, 0);
            break;
        case 32:
            fixcg (&currgrp.A[nprod - 1], un, 1, 0);
            fixcg (&currgrp.A[nprod + 1 - 1], sung, (*n), 0);
            fixcg (&currgrp.A[nprod + 2 - 1], sung, (*m), 0);
            nprod = nprod + 2;
            break;
        case 33:
            fixcg (&currgrp.A[nprod - 1], un, 1, 0);
            fixcg (&currgrp.A[nprod + 1 - 1], sunm, (*n), (*m));
            fixcg (&currgrp.A[nprod + 2 - 1], sunm, (*n2), (*m2));
            nprod = nprod + 2;
            break;
        case 34:
            fixcg (&currgrp.A[nprod - 1], unm, (*n), (*m));
            fixcg (&currgrp.A[nprod + 1 - 1], unm, (*n2), (*m2));
            nprod = nprod + 1;
            break;
        case 35:
            fixcg (&currgrp.A[nprod + 1 - 1], spn, (*m), 0);
            fixcg (&currgrp.A[nprod - 1], on, (*n), 0);
            nprod = nprod + 1;
            break;
        case 36:
        case 61:
            j = (*n) / 2;
            fixcg (&currgrp.A[nprod - 1], un, j, 0);
            break;
        case 37:
            j = (*n) / 2;
            fixcg (&currgrp.A[nprod - 1], spnc, 2, 0);
            nprod = nprod + 1;
            fixcg (&currgrp.A[nprod - 1], on, j, 0);
            break;
        case 38:
        case 39:
            fixcg (&currgrp.A[nprod - 1], spnc, (*n), 0);
            nprod = nprod + 1;
            fixcg (&currgrp.A[nprod - 1], on, (*m), 0);
            break;
        case 40:
            fixcg (&currgrp.A[nprod - 1], sung, 3, 0);
            break;
        case 42:
            fixcg (&currgrp.A[nprod - 1], son, 7, 0);
            break;
        case 43:
            fixcg (&currgrp.A[nprod - 1], son, 9, 0);
            break;
        case 44:
            fixcg (&currgrp.A[nprod - 1], sung, 2, 0);
            nprod = nprod + 1;
            fixcg (&currgrp.A[nprod - 1], sung, 6, 0);
            break;
        case 45:
            fixcg (&currgrp.A[nprod - 1], un, 1, 0);
            nprod = nprod + 1;
            fixcg (&currgrp.A[nprod - 1], son, 10, 0);
            break;
        case 47:
            fixcg (&currgrp.A[nprod - 1], sung, 8, 0);
            break;
        case 48:
            fixcg (&currgrp.A[nprod - 1], un, 1, 0);
            nprod = nprod + 1;
            fixcg (&currgrp.A[nprod - 1], e6, 6, 0);
            break;
        case 49:
            fixcg (&currgrp.A[nprod - 1], sung, 9, 0);
            break;
        case 50:
            fixcg (&currgrp.A[nprod - 1], son, 16, 0);
            break;
        case 51:
            fixcg (&currgrp.A[nprod - 1], sung, 2, 0);
            nprod = nprod + 1;
            fixcg (&currgrp.A[nprod - 1], e7, 7, 0);
            break;
        case 52:
            fixcg (&currgrp.A[nprod - 1], sung, 3, 0);
            nprod = nprod + 1;
            fixcg (&currgrp.A[nprod - 1], e6, 6, 0);
            break;
        case 53:
            fixcg (&currgrp.A[nprod - 1], e6, 6, 0);
            break;
        case 54:
            fixcg (&currgrp.A[nprod - 1], e7, 7, 0);
            break;
        case 55:
            fixcg (&currgrp.A[nprod - 1], e8, 8, 0);
            break;
        case 56:
        case 59:
            fixcg (&currgrp.A[nprod - 1], f4, 4, 0);
            break;
        case 57:
            fixcg (&currgrp.A[nprod - 1], son, 3, 0);
            nprod = nprod + 1;
            fixcg (&currgrp.A[nprod - 1], g2, 2, 0);
            break;
        case 58:
            fixcg (&currgrp.A[nprod - 1], sung, 3, 0);
            nprod = nprod + 1;
            fixcg (&currgrp.A[nprod - 1], g2, 2, 0);
            break;
        case 60:
            fixcg (&currgrp.A[nprod - 1], f4, 4, 0);
            nprod = nprod + 1;
            fixcg (&currgrp.A[nprod - 1], g2, 2, 0);
            break;
        case 63:
            fixcg (&currgrp.A[nprod - 1], l168, 168, 0);
            break;
        default:
            Caseerror (Line);
        }
}

void
dobranch (ocharptr * help, ocharptr * help1, int brno, int n,
          int m, int n2, int m2)
{
    char ch;
    ocharptr temp;
    redu = false;

    if (((brno >= 1) && (brno <= 63)))
        switch ((int) (brno)) {
        case 1:
            (*help1) = unospnbrnch ('d', (*help), n);
            break;
        case 2:
            if ((bool) ((n) & 1)) {
                (void) fprintf (output.fp, "ERROR:n must be an even int\n"),
                    Putl (output, 1);
                (void) fprintf (output.fp, "Enter stop\n"), Putl (output, 1);
                (*help1) = NULL;
            } else
                (*help1) = unospnbrnch ('b', (*help), n);
            break;
        case 3:
            (*help1) = unospnbrnch ('m', (*help), n);
            break;
        case 4:
            (*help1) = supqu1br ((*help), n, m, 0, 0, 4);
            break;
        case 5:
            (*help1) = sunmbrnch ((*help), n, m, 0, 0, 1);
            break;
        case 6:
            (*help1) = unsubgrbr ((*help), n, 1, 'b');
            break;
        case 7:
            (*help1) = unso3brnch ((*help), n);
            break;
        case 8:
            (*help1) = supqu1br ((*help), n, m, 0, 0, 1);
            break;
        case 9:
            if ((bool) ((n) & 1)) {
                (void) fprintf (output.fp, "ERROR:n must be an even int\n"),
                    Putl (output, 1);
                (void) fprintf (output.fp, "Enter stop\n"), Putl (output, 1);
                (*help1) = NULL;
            } else
                (*help1) = so2k1so3brnch ('a', (*help), n, true);
            break;
        case 10:
            if ((bool) ((n) & 1)) {
                (void) fprintf (output.fp, "ERROR:n must be an even int\n"),
                    Putl (output, 1);
                (void) fprintf (output.fp, "Enter stop\n"), Putl (output, 1);
                (*help1) = NULL;
            } else
                (*help1) = spnunbr ((*help), n);
            break;
        case 11:
            if ((bool) ((n) & 1)) {
                (void) fprintf (output.fp, "ERROR:n must be an even int\n"),
                    Putl (output, 1);
                (void) fprintf (output.fp, "Enter stop\n"), Putl (output, 1);
                (*help1) = NULL;
            } else if ((n <= 2)) {
                (void) fprintf (output.fp, "ERROR:n must be >2\n"),
                    Putl (output, 1);
                (void) fprintf (output.fp, "Enter stop\n"), Putl (output, 1);
                (*help1) = NULL;
            } else
                (*help1) = sp2nsu2son ((*help));
            break;
        case 12:
            if ((bool) ((n) & 1)) {
                (void) fprintf (output.fp, "ERROR:n must be an even int\n"),
                    Putl (output, 1);
                (void) fprintf (output.fp, "Enter stop\n"), Putl (output, 1);
                (*help1) = NULL;
            } else
                (*help1) = so2k1so3brnch ('a', (*help), n, false);
            break;
        case 13:
            (*help1) = unsubgrbr ((*help), n, 0, 'd');
            break;
        case 14:
            if (((bool) ((n) & 1) || (bool) ((m) & 1))) {
                (void) fprintf (output.fp, "ERROR:n must be an even int\n"),
                    Putl (output, 1);
                (void) fprintf (output.fp, "Enter stop\n"), Putl (output, 1);
                (*help1) = NULL;
            } else
                (*help1) = supqu1br ((*help), n, m, 0, 0, 5);
            break;
        case 15:
            (*help1) = onmbrnch ((*help), n, m, 3);
            break;
        case 16:
            (*help1) = supqu1br ((*help), n, m, 0, 0, 0);
            break;
        case 17:
            (*help1) = snanbrnch ((*help), n);
            break;
        case 18:
            (*help1) = onsnbr ((*help), n, true);
            break;
        case 19:
            (*help1) = onsnbr ((*help), n, false);
            break;
        case 20:
            (*help1) = so2k1so3brnch ('c', (*help), n, false);
            break;
        case 21:
            if ((bool) ((n) & 1))
                ch = 'f';
            else
                ch = 'b';
            (*help1) = unsubgrbr ((*help), n, 0, ch);
            break;
        case 22:
            (*help1) = supqu1br ((*help), n, m, 0, 0, 6);
            break;
        case 23:
            (*help1) = onmbrnch ((*help), n, m, 1);
            break;
        case 24:
            (*help1) = onmbrnch ((*help), n, m, 2);
            break;
        case 25:
            if ((!(bool) ((n) & 1))) {
                print ("Even O(n) -> SO(3) not implemented\n");
                (*help1) = NULL;
            } else
                (*help1) = so2k1so3brnch ('c', (*help), n, true);
            break;
        case 26:
            (*help1) = sonunbrnch ((*help), n);
            break;
        case 27:
            (*help1) = onmonombrnch ((*help), n, m, 1);
            break;
        case 28:
            (*help1) = so4su2su2 ((*help));
            break;
        case 29:
            (*help1) = so4so3 ((*help));
            break;
        case 30:
            (*help1) = so2k1so3brnch ('c', (*help), 7, true);
            break;
        case 31:
            (*help1) = so7g2brnch ((*help));
            break;
        case 32:
            (*help1) = supqu1br ((*help), n, m, 0, 0, 2);
            break;
        case 33:
            (*help1) = supqu1br ((*help), n, m, n2, m2, 3);
            break;
        case 34:
            (*help1) = sunmbrnch ((*help), n, m, n2, m2, 2);
            break;
        case 35:
            if ((bool) ((m) & 1)) {
                (void) fprintf (output.fp, "ERROR:n must be an even int\n"),
                    Putl (output, 1);
                (void) fprintf (output.fp, "Enter stop\n"), Putl (output, 1);
                (*help1) = NULL;
            } else
                (*help1) = onmonombrnch ((*help), m, n, 2);
            break;
        case 36:
            if ((bool) ((n) & 1)) {
                (void) fprintf (output.fp, "ERROR:n must be an even int\n"),
                    Putl (output, 1);
                (void) fprintf (output.fp, "Enter stop\n"), Putl (output, 1);
                (*help1) = NULL;
            } else
                (*help1) = lsprun (&(*help), n);
            break;
        case 37:
            if ((bool) ((n) & 1)) {
                (void) fprintf (output.fp, "ERROR:n must be an even int\n"),
                    Putl (output, 1);
                (void) fprintf (output.fp, "Enter stop\n"), Putl (output, 1);
                (*help1) = NULL;
            } else
                (*help1) = lspnrsp2on (&(*help), n);
            break;
        case 38:
            if ((bool) ((n) & 1)) {
                (void) fprintf (output.fp, "ERROR:n must be an even int\n"),
                    Putl (output, 1);
                (void) fprintf (output.fp, "Enter stop\n"), Putl (output, 1);
                (*help1) = NULL;
            } else
                (*help1) = sp2nk (&(*help), n, m);
            break;
        case 39:
            if ((bool) ((n) & 1)) {
                (void) fprintf (output.fp, "ERROR:n must be an even int\n"),
                    Putl (output, 1);
                (void) fprintf (output.fp, "Enter stop\n"), Putl (output, 1);
                (*help1) = NULL;
            } else
                (*help1) = metaplet (n, m);
            break;
        case 40:
            (*help1) = g2brnch ((*help));
            break;
        case 41:
            (*help1) = g2so3brnch ((*help));
            break;
        case 42:
            (*help1) = g2_so7 ((*help));
            break;
        case 43:
            (*help1) = ebrnch ((*help), f4index, f4tab);
            break;
        case 44:
            (*help1) = ebrnch ((*help), e6index, e6tab);
            break;
        case 45:
            (*help1) = ebrnch ((*help), e6soindex, e6sotab);
            break;
        case 46:
            (*help1) = ebrnch ((*help), e6g2index, e6g2tab);
            break;
        case 47:
            (*help1) = ebrnch ((*help), e7index, e7tab);
            break;
        case 48:
            (*help1) = ebrnch ((*help), e7e6index, e7e6tab);
            break;
        case 49:
            (*help1) = ebrnch ((*help), e8index, e8tab);
            break;
        case 50:
            (*help1) = ebrnch ((*help), e8soindex, e8sotab);
            break;
        case 51:
            (*help1) = ebrnch ((*help), e8suindex, e8sutab);
            break;
        case 52:
            (*help1) = ebrnch ((*help), e8e6index, e8e6tab);
            break;
        case 53:
            (*help1) = ebrnch ((*help), u27e6index, u27e6tab);
            break;
        case 54:
            (*help1) = ebrnch ((*help), su56e7index, su56e7tab);
            break;
        case 55:
            (*help1) = ebrnch ((*help), su248e8index, su248e8tab);
            break;
        case 56:
            (*help1) = ebrnch ((*help), e6f4index, e6f4tab);
            break;
        case 57:
            (*help1) = ebrnch ((*help), f4g2index, f4g2tab);
            break;
        case 58:
            (*help1) = ebrnch ((*help), e6su3g2index, e6su3g2tab);
            break;
        case 59:
            (*help1) = ebrnch ((*help), so26f4index, so26f4tab);
            break;
        case 60:
            (*help1) = ebrnch ((*help), e8f4g2index, e8f4g2tab);
            break;
        case 61:
            if ((bool) ((n) & 1)) {
                (void) fprintf (output.fp, "ERROR:n must be an even int\n"),
                    Putl (output, 1);
                (void) fprintf (output.fp, "Enter stop\n"), Putl (output, 1);
                (*help1) = NULL;
            } else
                (*help1) = lsoncbrun (&(*help), n);
            break;
        case 62:
            fixcg (&currgrp.A[nprod - 1], on, n, 0);
            temp = unospnbrnch ('d', (*help), n);
            (*help1) = onsnbr (temp, n, true);
            fixcg (&currgrp.A[nprod - 1], sn, n, 0);
            odisp (&temp);
            break;
        case 63:
            (*help1) = ebrnch ((*help), l168index, l168tab);
            break;
        default:
            Caseerror (Line);
        }
    odisp (&(*help));

}

void
fixgroups (prodtype * pr, int brno, int nprd)
{
    frame v, c, zerofrm;
    char l = ' ', cl = ' ';     // orig. not init. FB
    bool s = false;             //, d=false;// orig. not init. FB
    int m = 1;                  // orig. not init. FB
    register int j;
    for (j = 0; j <= maxdim; j++) {
        zerofrm.A[j] = 0;
    }
    while ((*pr) != NULL) {
        register prodtype W22 = &(*(*pr));

        /*if (!(Member((unsigned)(brno), Conset[0])))

           {  
           m = 1;
           v = W22->prods.A[nprd - 1]->val;
           c = W22->prods.A[nprd - 1]->conval;
           l = W22->prods.A[nprd - 1]->lab;
           cl = W22->prods.A[nprd - 1]->conlab;
           s = W22->prods.A[nprd - 1]->spin;
           d = W22->prods.A[nprd - 1]->C6_double;
           odisp(&W22->prods.A[nprd - 1]);
           } */ /**/
            switch ((int) (brno)) {
        case 1:
        case 2:
        case 3:
        case 18:
        case 19:
        case 7:
        case 25:
        case 9:
        case 20:
        case 12:
        case 31:
        case 30:
        case 40:
        case 41:
        case 42:
        case 43:
        case 46:
        case 47:
        case 49:
        case 50:
        case 53:
        case 54:
        case 55:
        case 56:
        case 17:
        case 29:
        case 36:
        case 13:
        case 21:
        case 6:
        case 61:
        case 59:
        case 62:
        case 63:
            break;
        case 4:
        case 5:
        case 8:
        case 10:
        case 11:
        case 14:
        case 15:
        case 16:
        case 22:
        case 23:
        case 24:
        case 26:
        case 27:
        case 28:
        case 32:
        case 33:
        case 34:
        case 35:
        case 37:
        case 38:
        case 39:
        case 44:
        case 45:
        case 48:
        case 51:
        case 52:
        case 57:
        case 58:
        case 60:
        case 65:
            m = 1;
            v = W22->prods.A[nprd - 1]->val;
            c = W22->prods.A[nprd - 1]->conval;
            l = W22->prods.A[nprd - 1]->lab;
            cl = W22->prods.A[nprd - 1]->conlab;
            s = W22->prods.A[nprd - 1]->spin;
            // d = W22->prods.A[nprd - 1]->C6_double;
            odisp (&W22->prods.A[nprd - 1]);
            break;
        default:
            Caseerror (Line);
        }


        if (((brno >= 1) && (brno <= 63)))
            switch ((int) (brno)) {
            case 1:
            case 2:
            case 3:
            case 6:
            case 7:
            case 9:
            case 12:
            case 13:
            case 17:
            case 18:
            case 19:
            case 20:
            case 21:
            case 25:
            case 29:
            case 30:
            case 31:
            case 36:
            case 40:
            case 41:
            case 42:
            case 43:
            case 46:
            case 47:
            case 49:
            case 50:
            case 53:
            case 54:
            case 55:
            case 56:
            case 59:
            case 61:
            case 62:
            case 63:
                break;
            case 4:
            case 5:
            case 11:
            case 14:
            case 15:
            case 16:
            case 22:
            case 23:
            case 24:
            case 27:
            case 28:
            case 34:
                fixochar (&W22->prods.A[nprd - 1], m, c, zerofrm, cl, ' ', s,
                          false);
                fixochar (&W22->prods.A[nprod - 1], m, v, zerofrm, l, ' ', s,
                          false);
                break;
            case 37:
            case 38:
            case 39:
                fixochar (&W22->prods.A[nprd - 1], m, c, zerofrm, cl, ' ', s,
                          false);
                fixochar (&W22->prods.A[nprod - 1], m, v, zerofrm, l, ' ',
                          false, false);
                break;
            case 26:
            case 10:
            case 44:
            case 51:
            case 57:
                fixochar (&W22->prods.A[nprd - 1], m, zerofrm, zerofrm, ' ',
                          ' ', false, false);
                W22->prods.A[nprd - 1]->val.A[1] = v.A[1];
                for (j = 1; j <= maxdim - 1; j++) {
                    v.A[j] = v.A[j + 1];
                }
                v.A[maxdim] = 0;
                fixochar (&W22->prods.A[nprod - 1], m, v, zerofrm, l, ' ',
                          false, false);
                break;
            case 8:
            case 32:
            case 33:
                fixochar (&W22->prods.A[nprd - 1], m, zerofrm, zerofrm, ' ',
                          ' ', false, false);
                W22->prods.A[nprd - 1]->val.A[1] = v.A[1];
                for (j = 1; j <= maxdim - 1; j++) {
                    v.A[j] = v.A[j + 1];
                }
                v.A[maxdim] = 0;
                fixochar (&W22->prods.A[nprd + 1 - 1], m, c, zerofrm, cl, ' ',
                          false, false);
                fixochar (&W22->prods.A[nprod - 1], m, v, zerofrm, l, ' ',
                          false, false);
                break;
            case 35:
            case 52:
            case 58:
                fixochar (&W22->prods.A[nprod - 1], m, c, zerofrm, cl, ' ', s,
                          false);
                fixochar (&W22->prods.A[nprd - 1], m, v, zerofrm, l, ' ', s,
                          false);
                break;
            case 60:
                fixochar (&W22->prods.A[nprd - 1], m, v, zerofrm, ' ', ' ', s,
                          false);
                fixochar (&W22->prods.A[nprod - 1], m, c, zerofrm, ' ', ' ',
                          false, false);
                break;
            case 45:
            case 48:
                fixochar (&W22->prods.A[nprd - 1], m, v, zerofrm, ' ', ' ',
                          false, false);
                fixochar (&W22->prods.A[nprod - 1], m, c, zerofrm, l, ' ', s,
                          false);
                break;
            default:
                Caseerror (Line);
            }
        (*pr) = (*pr)->next;
    }
}

/* UNUSED
void
su1crunch (prodtype * pr)
{
    prodtype top;
    int su1;
    register int i;
    register int j;

    do {
        su1 = 0;
        for (j = 1; j <= nprod; j++) {
            if (((currgrp.A[j - 1].name == sung) || (currgrp.A[j - 1].name == son)) && (currgrp.A[j - 1].rank == 1))   
                su1 = j;
        }
        if (su1 != 0) {
            top = (*pr);
            while (top != NULL) {
                register prodtype W29 = &(*top);

                odisp (&W29->prods.A[su1 - 1]);
                {
                    int lastStep31 = (nprod - 1);
                    for (i = su1; i <= lastStep31; i++) {
                        W29->prods.A[i - 1] = W29->prods.A[i + 1 - 1];
                    }
                }
                W29->prods.A[nprod - 1] = NULL;
                top = top->next;
            }
            {
                int lastStep33 = (nprod - 1);
                for (i = su1; i <= lastStep33; i++) {
                    currgrp.A[i - 1] = currgrp.A[i + 1 - 1];
                }
            }
            nprod = nprod - 1;
        }
    }
    while (!(su1 == 0));
}
*/

void
putgroup1 (groopArray grp, int jj)
{
    register int j;
    bool test;

    if (jj < 1)
        error (GROUP_NOT_SET, 20);
    if (jj == 1)
        if (nprod > 1)
            inform ("Groups are ;", cont);
        else
            inform ("Group is ;", cont);
    for (j = 1; j <= nprod; j++) {
        /*if ((grp.A[j - 1].name == sung) && (grp.A[j - 1].rank == 1))
           test = false;
           else *//*15/2/98 */
        test = true;
        if ((nprod > 1) && (j <= nprod) && (j > 1) && test)
            inform (" * ;", cont);
        switch ((int) (grp.A[j - 1].name)) {
        case sung:
            if (test)
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
        case spnc:
            inform ("Sp(;", cont);
            break;
        case sonc:
            inform ("SO^*(;", cont);
            break;
        case mp:
            inform ("Mp(;", cont);
            break;
        case ospnm:
            inform ("OSp(;", cont);
            break;
        case sn:
            inform ("S(;", cont);
            if ((bool) ((grp.A[j - 1].rank) & 1))
                qsn = 1;
            else
                qsn = 0;
            break;
        case an:
            inform ("A(;", cont);
            if ((bool) ((grp.A[j - 1].rank) & 1))
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
        case l168:
            inform ("L(;", cont);
            break;
        case nill:
            warn ("No group set;", cr);
            break;
        default:
            Caseerror (Line);
        }
        if ((grp.A[j - 1].name != nill) && test) {
            print ("%1d", grp.A[j - 1].rank);
            if (((grp.A[j - 1].name == ospnm) || (grp.A[j - 1].name == unm)
                 || (grp.A[j - 1].name == sunm)))
                print ("/%1d", grp.A[j - 1].rank2);
            if ((grp.A[j - 1].name == unc))
                print (",%1d", grp.A[j - 1].rank2);
            if (grp.A[j - 1].name == spnc)
                print (",R");
            inform (");", cont);
        }
    }
    if (jj == 1)
        inform (";", cr);
    else
        inform (";", cont);
}

prodtype
unu1 (int n)
{
    prodtype temp, newlist, list;
    ocharptr list1, list2;
    int signq, signd;
    register int i;
    register int q;

    newlist = NULL;
    if ((bool) ((n) & 1))
        signd = 1;
    else
        signd = -1;
    {
        int lastStep37 = (n - 1);
        for (q = 0; q <= lastStep37; q++) {
            if ((bool) ((q) & 1))
                signq = -signd;
            else
                signq = signd;
            cnu (&list1);
            list1->mult = 1;
            list1->val = nolls;
            list1->C6_double = false;
            list1->spin = false;
            list1->lab = ' ';
            list1->val.A[1] = q + 1;
            {
                int lastStep39 = (n - q);
                for (i = 2; i <= lastStep39; i++) {
                    list1->val.A[i] = 1;
                }
            }
            cnu (&list2);
            list2->mult = 1;
            list2->val = nolls;
            list2->C6_double = false;
            list2->spin = false;
            list2->lab = ' ';
            list2->val.A[1] = q;
            pnu (&list);
            list->mult = signq;
            list->prods.A[1 - 1] = list1;
            list->prods.A[2 - 1] = list2;
            temp = prodadd (newlist, list);
            pdisp (&list);
            pdisp (&newlist);
            newlist = temp;
        }
    }
    return newlist;
}

prodtype
hecke (frame a)
{
    prodtype newlist, temp1, temp2;
    int l;
    register int i;
    /*ocharptr      tchrc, chrc; *//*12/12/95 */

    l = len (&a);
    newlist = unu1 (a.A[1]);
    if (l > 1)
        for (i = 2; i <= l; i++) {
            temp1 = unu1 (a.A[i]);
            temp2 = pkronk (newlist, temp1);
            pdisp (&newlist);
            pdisp (&temp1);
            newlist = temp2;
        }
    return newlist;
}

void
schar (ocharptr reps, ocharptr class)
{
    ocharptr temp, newlist;
    int clen, pt, gt;
    register int i;
    bool test;

    if ((wtfrm (&reps->val) != wtfrm (&class->val))
        || (currgrp.A[1 - 1].rank != wtfrm (&reps->val)))
        print ("ERROR: Partitions not of the correct weight\n");
    else if ((reps->spin == true))
        print ("Not implemented for spin characteristics\n");
    else {
        test = true;
        clen = len (&class->val);
        gt = currgrp.A[1 - 1].rank;
        cnu (&newlist);
        newlist->mult = reps->mult;
        newlist->val = reps->val;
        for (i = 1; i <= clen; i++) {
            if (test == true) {
                temp = removek (newlist, class->val.A[i]);
                if (temp == NULL)
                    test = false;
                odisp (&newlist);
                newlist = temp;
            }
        }
        if (test)
            pt = newlist->mult;
        else
            pt = 0;
        print ("characteristic =%10d\n", pt);
        currgrp.A[1 - 1].rank = gt;
        dispchr (&newlist);
    }
}

termptr
character (termptr sfn, termptr class)
{
    termptr slist, clist, temp, newlist;
    int clen;
    register int i;
    bool test;

    test = true;
    clen = len (&class->val);
    snu (&slist);
    slist->mult = 1;
    slist->val = sfn->val;
    snu (&newlist);
    newlist->mult = sfn->mult;
    newlist->val = sfn->val;
    snu (&clist);
    clist->mult = 1;
    clist->val = class->val;
    clist->next = NULL;
    for (i = 1; i <= clen; i++) {
        if (test == true) {
            temp = (termptr) sremovek (newlist, clist->val.A[i]);
            if (temp == NULL)
                test = false;
            ldisp (&newlist);
            newlist = temp;
        }
    }
    dispsfn (&clist);
    if (test == true)
        slist->mult = newlist->mult;
    else
        dispsfn (&slist);
    dispsfn (&temp);
    return slist;
}

termptr
sremovek (termptr sfn, int k)
{
    int length;
    register int i;
    termptr newlist, temp;

    newlist = NULL;
    while (sfn != NULL) {
        register termptr W50 = sfn;

        length = len (&W50->val);
        for (i = 1; i <= length; i++) {
            snu (&temp);
            temp->mult = W50->mult;
            temp->val = W50->val;
            temp->val.A[i] = temp->val.A[i] - k;
            add (&newlist, &temp);
            dispsfn (&temp);
        }
        sfn = W50->next;
    }
    stndise (&newlist);
    sort (&newlist, true);
    return newlist;
}


termptr
snchar (termptr class, int leng)
{
    termptr dummy, temp, templist, newlist;
    int i, k;

    snu (&newlist);
    newlist->mult = 1;
    newlist->val = nolls;
    newlist->next = NULL;
    k = len (&class->val);
    i = 1;
    while ((i <= k)) {
        snu (&temp);
        temp->mult = 1;
        temp->val = nolls;
        temp->val.A[1] = class->val.A[i];
        temp->next = NULL;
        templist = (termptr) powersumtosfn (temp);
        dispsfn (&temp);
        schur_restrict (&templist, leng, 'l');
        dummy = louter2 (newlist, templist, leng);
        ldisp (&newlist);
        ldisp (&templist);
        newlist = dummy;
        i = i + 1;
    }
    return newlist;
}

termptr
gensfnlist (int n)
{
    termptr sfn1, tempsfn;
    int k = 0;                  // orig. not init. FB

    snu (&sfn1);
    sfn1->mult = 1;
    sfn1->val = nolls;
    sfn1->val.A[1] = n;
    tempsfn = useseries ('f', sfn1, true, false, k);
    dispsfn (&sfn1);
    schur_restrict (&tempsfn, -n, 'w');
    return tempsfn;
}

void
swapgroupname (int gr1, int gr2)
{
    grptype tempname;
    int r1, r2;

    tempname = currgrp.A[gr1 - 1].name;
    r1 = currgrp.A[gr1 - 1].rank;
    r2 = currgrp.A[gr1 - 1].rank2;
    currgrp.A[gr1 - 1] = currgrp.A[gr2 - 1];
    currgrp.A[gr2 - 1].name = tempname;
    currgrp.A[gr2 - 1].rank = r1;
    currgrp.A[gr2 - 1].rank2 = r2;
    putgroup (currgrp);
}

prodtype
swapgroups (prodtype list, int gr1, int gr2)
{
    prodtype top;
    ocharptr h1, h2;

    top = list;
    if (((gr1 != gr2) && (gr1 <= nprod) && (gr2 <= nprod) && (gr1 > 0)
         && (gr2 > 0))) {
        while (top != NULL) {
            register prodtype W53 = top;

            h1 = W53->prods.A[gr1 - 1];
            h2 = W53->prods.A[gr2 - 1];
            W53->prods.A[gr1 - 1] = h2;
            W53->prods.A[gr2 - 1] = h1;
            top = top->next;
        }
        swapgroupname (gr1, gr2);
        top = list;
        schur_psort (&top, true);
    } else {
        print ("ERROR: check range of group ints\n");
        top = NULL;
    }
    return top;
}

ocharptr
removek (ocharptr reps, int k)
{
    int length;
    register int i;
    ocharptr newlist, temp;

    newlist = NULL;
    currgrp.A[1 - 1].rank = currgrp.A[1 - 1].rank - k;
    while (reps != NULL) {
        register ocharptr W45 = reps;

        length = len (&W45->val);
        for (i = 1; i <= length; i++) {
            cnu (&temp);
            temp->mult = W45->mult;
            temp->val = W45->val;
            temp->val.A[i] = temp->val.A[i] - k;
            temp->spin = false;
            temp->C6_double = false;
            temp->lab = ' ';
            oadd (&newlist, &temp);
            dispchr (&temp);
        }
        reps = W45->next;
    }
    snselect (&newlist, currgrp.A[1 - 1].rank, false);
    return newlist;
}

//doub : mixed series or not
prodtype
macseries (char ser, bool w, bool signs, bool r, bool d, bool doub, int n)
{
    register termptr W47;
    termptr tser, temp;
    ocharptr list1, list2;
    prodtype flist, plist, tlist;
    int wt, rk, tone;

    plist = NULL;
    temp = seriesx (ser, n, full);
    tser = temp;
    while (tser != NULL) {
        W47 = tser;

        cnu (&list1);
        cnu (&list2);
        list1->mult = 1;
        list2->mult = 1;
        list1->spin = false;
        list2->spin = false;
        list1->C6_double = false;
        list2->C6_double = false;
        list1->val = nolls;
        list2->val = W47->val;
        if (doub) {
            list2->C6_double = true;
            conjgte (&W47->val);
            list2->conval = W47->val;
        }
        wt = 0;
        rk = 0;
        if (w)
            wt = wtfrm (&W47->val);
        if (r)
            rk = frank (W47->val);
        if (signs)
            tone = wt + rk;
        else
            tone = wt - rk;
        if (d)
            tone = tone / 2;
        list1->val.A[1] = tone;
        pnu (&flist);
        flist->prods.A[1 - 1] = list1;
        flist->prods.A[2 - 1] = list2;
        if (doub)               // m series FB 2006
            flist->mult =
                (tone % 2 == 0 ? abs (W47->mult) : -abs (W47->mult));
        else
            flist->mult = W47->mult;
        tlist = prodadd (plist, flist);
        pdisp (&plist);
        pdisp (&flist);
        plist = tlist;
        tser = W47->next;
    }
    ldisp (&temp);
    return plist;
}

prodtype
inverseseries (prodtype list, int wmax)
{
    prodtype tlist, plist, newlist, templist;
    int i, save_plwt;

    newlist = NULL;
    save_plwt = plwt;
    plwt = wmax;
    templist = prodcopy (list);
    i = 0;
    while (i <= wmax) {
        plist = pwrestrict (templist, -i, 1, 'w');
        if (i > 0) {
            tlist = prodmult (-1, plist);
            pdisp (&plist);
            plist = prodcopy (tlist);
            pdisp (&tlist);
        }
        tlist = prodadd (newlist, plist);
        pdisp (&plist);
        pdisp (&newlist);
        newlist = prodcopy (tlist);
        if (i < wmax)
            templist = pkronk (list, newlist);
        pdisp (&tlist);
        i = i + 1;
    }
    plwt = save_plwt;
    return newlist;
}

termptr
soncun (termptr * list, int k, int n)
{
    termptr newlist, temp, list1, list2;
    int m, w1;
    //char tag;
    bool ppp;

    newlist = NULL;
    m = MIN (n, k);
    ppp = false;
    //tag = 'a';
    w1 = wtfrm (&(*list)->val);
    if (w1 >= 1)
        w1 = w1 + setlimit - 1;
    else
        w1 = setlimit - 2;
    list1 = (termptr) signseq (&(*list), k, m, 'a', ppp, false);
    list2 = rseries (m, 'b');
    temp = louter2 (list1, list2, m);
    schur_restrict (&temp, w1, 'w');
    ldisp (&list1);
    ldisp (&list2);
    add (&newlist, &temp);
    sort (&newlist, true);
    return newlist;
}

ocharptr
soncbrun (ocharptr * list, int n)
{
    int i, k, t, r;
    register int j;
    ocharptr temp;
    termptr plist, slist, eta, nlist;

    temp = NULL;
    t = n / 2;
    i = (*list)->val.A[maxdim];
    k = i;
    k = 2 * k;
    snu (&eta);
    eta->mult = 1;
    eta->val = nolls;
    for (j = 1; j <= t; j++) {
        eta->val.A[j] = i;
    }
    snu (&nlist);
    nlist->mult = (*list)->mult;
    nlist->val = (*list)->val;
    r = setlimit;
    if (wtfrm (&(*list)->val) == 0)
        r = r - 1;
    slist = soncun (&nlist, k, t);
    schur_restrict (&slist, r, 'w');
    plist = louter2 (eta, slist, t);
    ldisp (&slist);
    dispsfn (&nlist);
    temp = sfntochrc (plist, false, ' ');
    ldisp (&plist);
    dispsfn (&eta);
    return temp;
}

ocharptr
lsoncbrun (ocharptr * llist, int n)
{
    ocharptr newlist, temp, tlist;

    newlist = NULL;
    tlist = (*llist);
    while (tlist != NULL) {
        register ocharptr W45 = tlist;

        temp = soncbrun (&tlist, n);
        oadd (&newlist, &temp);
        tlist = W45->next;
    }
    osort (&newlist, true);
    return newlist;
}
