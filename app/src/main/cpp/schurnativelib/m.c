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

/** \file m.c
 */

/*
**	Definitions for i/o
*/
#include <stdio.h>
#include <stdlib.h>
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
#include "g.h"
#include "m.h"

void
modsonc (ocharptr * list, int n)
{
    ocharptr newlist, tlist;
    int k, p, l;

    newlist = NULL;
    while (*list != NULL) {
        register ocharptr W2 = *list;
        k = W2->val.A[maxdim];
        p = MIN (k, n);
        l = len (&W2->val);
        cnu (&tlist);
        tlist->mult = W2->mult;
        tlist->val = W2->val;
        tlist->val.A[maxdim] = 0;
        tlist->C6_double = true;
        tlist->spin = false;
        tlist->conval.A[1] = k;
        tlist->lab = ' ';
        tlist->conlab = ' ';
        if ((l <= p)) {
            tlist->val = W2->val;
            oadd (&newlist, &tlist);
        }
        dispchr (&tlist);
        *list = W2->next;
    }
    osort (&newlist, true);
    *list = newlist;
}

void
stdpq (ocharptr * list, int p, int q, int kks)
{
    ocharptr newlist, tlist;

    newlist = NULL;
    while ((*list) != NULL) {
        register ocharptr W2 = &(*(*list));

        cnu (&tlist);
        tlist->mult = W2->mult;
        tlist->val = W2->val;
        tlist->C6_double = true;
        tlist->spin = W2->spin;
        tlist->conval = W2->conval;
        tlist->lab = ' ';
        tlist->conlab = ' ';
        conjgte (&tlist->val);
        conjgte (&tlist->conval);
        if (((tlist->val.A[1] + tlist->conval.A[1]) <= kks)
            && (tlist->val.A[1] <= p) && (tlist->conval.A[1] <= q)) {
            tlist->val = W2->val;
            oadd (&newlist, &tlist);
        }
        dispchr (&tlist);
        (*list) = W2->next;
    }
    osort (&newlist, true);
    (*list) = newlist;
}

void
modspnr (ocharptr * list, int n)
{
    ocharptr newlist, temp, tlist;
    int x, k;
    newlist = NULL;
    temp = (*list);
    while ((*list) != NULL) {
        register ocharptr W2 = &(*(*list));
        x = W2->val.A[maxdim];
        if (W2->spin)
            k = 2 * x + 1;
        else
            k = 2 * x;
        cnu (&tlist);
        tlist->mult = W2->mult;
        tlist->val = W2->val;
        tlist->val.A[maxdim] = 0;
        tlist->C6_double = true;
        tlist->spin = W2->spin;
        tlist->conval.A[1] = x;
        tlist->lab = ' ';
        tlist->conlab = ' ';
        conjgte (&tlist->val);
        if ((tlist->val.A[1] <= (n / 2))
            && ((tlist->val.A[1] + tlist->val.A[2]) <= k)) {
            tlist->val = W2->val;
            oadd (&newlist, &tlist);
        }
        dispchr (&tlist);
/*	osort(&(*newlist),true);*/
        (*list) = W2->next;
    }
    odisp (&temp);
    osort (&newlist, true);
    (*list) = newlist;
}

void
nstndise (termptr * list, int n)
{
    enum { nonstand0, stand0, null0 } part0;
    termptr newlist, temp;
    int i, dummy;
    register termptr W6;

    newlist = NULL;
    if (n > 1) {
        while ((*list) != NULL) {
            W6 = &(*(*list));

            do {
                i = 1;
                part0 = stand0;
                do {
                    if (W6->val.A[i] < W6->val.A[i + 1]) {
                        if (W6->val.A[i] == W6->val.A[i + 1] - 1)
                            part0 = null0;
                        else {
                            part0 = nonstand0;
                            W6->mult = -W6->mult;
                            dummy = W6->val.A[i];
                            W6->val.A[i] = W6->val.A[i + 1] - 1;
                            W6->val.A[i + 1] = dummy + 1;
                        }
                    }
                    //if (W6->val.A[q] < 0) corrected by FB in:
                    if (W6->val.A[i] < 0)
                        part0 = null0;
                    i = i + 1;
                }
                while (!((i == n) || (part0 == null0)));
            }
            while (!(part0 != nonstand0));
            temp = W6->next;
            if (part0 == stand0) {
                W6->next = newlist;
                newlist = (*list);
            } else
                dispsfn (&(*list));
            (*list) = temp;
        }
    } else
        print ("ERROR: the int must be greater than 1\n");
    (*list) = newlist;
}

void
stndise (termptr * list)
{
    enum { nonstand1, stand1, null1 } part1;
    termptr newlist, temp;
    int i, dummy, q;
    register termptr W10;

    newlist = NULL;
    while (*list != NULL) {
        W10 = *list;

        q = qlen (W10->val);
        do {
            i = 1;
            part1 = stand1;
            do {
                if (W10->val.A[i] < W10->val.A[i + 1]) {
                    if (W10->val.A[i] == W10->val.A[i + 1] - 1)
                        part1 = null1;
                    else {
                        part1 = nonstand1;
                        W10->mult = -W10->mult;
                        dummy = W10->val.A[i];
                        W10->val.A[i] = W10->val.A[i + 1] - 1;
                        W10->val.A[i + 1] = dummy + 1;
                    }
                }
                if (W10->val.A[q] < 0)
                    part1 = null1;
                i = i + 1;
            }
            while ((i < q + 1) && (i < maxdim - 1) && (part1 != null1));

            /*while (i>0 && W10->val.A[i]==0)
               i--;
               W10->val.length=i; */
        }
        while (part1 == nonstand1);
        temp = W10->next;
        if (part1 == stand1) {
            W10->next = newlist;
            newlist = *list;
        } else
            dispsfn (list);
        (*list) = temp;
    }
    (*list) = newlist;
}

void
qstndise (termptr * list)
{
    enum { nonstand2, stand2, null2 } part2;
    termptr newlist, temp;
    int l, i, dummy;
    register int k;
    register termptr W14;

    qup (&(*list));
    newlist = NULL;
    while ((*list) != NULL) {
        W14 = &(*(*list));

        l = qlen (W14->val);
        do {
            i = 1;
            part2 = stand2;
            do {
                if (W14->val.A[l] < 0)
                    part2 = null2;
                if ((W14->val.A[i] == W14->val.A[i + 1])
                    && (W14->val.A[i] > 0))
                    part2 = null2;
                if (W14->val.A[i] < 0) {
                    if ((abs (W14->val.A[i]) > abs (W14->val.A[i + 1])))
                        part2 = null2;
                    if ((W14->val.A[i] == -W14->val.A[i + 1])) {
                        part2 = nonstand2;
                        if ((W14->val.A[i + 1] == W14->val.A[i + 2]))
                            part2 = null2;
                        else {
                            if ((bool) ((W14->val.A[i]) & 1))
                                W14->mult = -2 * W14->mult;
                            else
                                W14->mult = 2 * W14->mult;
                            for (k = i; k <= l; k++) {
                                if (k + 2 > l)
                                    dummy = 0;
                                else
                                    dummy = W14->val.A[k + 2];
                                W14->val.A[k] = dummy;
                            }
                            l = l - 2;
                            W14->val.A[l + 1] = 0;
                            W14->val.A[l + 2] = 0;
                        }
                    }
                }
                i = i + 1;
            }
            while (!((i >= l) || (part2 == null2)));
        }
        while (!(part2 != nonstand2));
        temp = W14->next;
        if (part2 == stand2) {
            W14->next = newlist;
            newlist = (*list);
        } else
            dispsfn (&(*list));
        (*list) = temp;
    }
    Qsort (&newlist, true);     /*8/6/99 */
    (*list) = newlist;
}

void
ostndise (ocharptr * list)
{
    enum { nonstand3, stand3, null3 } part3;
    ocharptr newlist;
    ocharptr temp;
    register int i, dummy;
    register ocharptr temp2;

    newlist = NULL;
    while (*list != NULL) {
        temp2 = *list;
        do {
            i = 1;
            part3 = stand3;
            do {
                if (temp2->val.A[i] < temp2->val.A[i + 1]) {
                    if (temp2->val.A[i] == temp2->val.A[i + 1] - 1)
                        part3 = null3;
                    else {
                        part3 = nonstand3;
                        temp2->mult = -temp2->mult;
                        dummy = temp2->val.A[i];
                        temp2->val.A[i] = temp2->val.A[i + 1] - 1;
                        temp2->val.A[i + 1] = dummy + 1;
                    }
                }
                if (temp2->C6_double && (part3 != null3)) {
                    if (temp2->conval.A[i] < temp2->conval.A[i + 1]) {
                        if (temp2->conval.A[i] == temp2->conval.A[i + 1] - 1)
                            part3 = null3;
                        else {
                            part3 = nonstand3;
                            temp2->mult = -temp2->mult;
                            dummy = temp2->conval.A[i];
                            temp2->conval.A[i] = temp2->conval.A[i + 1] - 1;
                            temp2->conval.A[i + 1] = dummy + 1;
                        }
                    }
                }
                i = i + 1;
            }
            while ((i < (maxdim - 1)) && (part3 != null3));
        }
        while (part3 == nonstand3);
        temp = temp2->next;
        if (part3 == stand3) {
            temp2->next = newlist;
            newlist = *list;
        } else
            dispchr (list);
        *list = temp;
    }
    *list = newlist;
    //printf ("toto");
}

void
f4stndise (ocharptr * list)
{
    enum { nonstand4, stand4, null4 } part4;
    ocharptr newlist, temp;
    int i, dummy;
    register ocharptr W24;

    newlist = NULL;
    while ((*list) != NULL) {
        W24 = &(*(*list));

        do {
            i = 1;
            part4 = stand4;
            do {
                if (W24->val.A[i] < W24->val.A[i + 1]) {
                    if (W24->val.A[i] == W24->val.A[i + 1] - 1)
                        part4 = null4;
                    else {
                        part4 = nonstand4;
                        W24->mult = -W24->mult;
                        dummy = W24->val.A[i];
                        W24->val.A[i] = W24->val.A[i + 1] - 1;
                        W24->val.A[i + 1] = dummy + 1;
                    }
                }
                i = i + 1;
            }
            while (!((i == 4) || (part4 == null4)));
            if ((part4 == stand4) && (W24->val.A[4] < 0)) {
                W24->val.A[4] = -W24->val.A[4] - 1;
                W24->mult = -W24->mult;
                part4 = nonstand4;
                if (W24->spin) {
                    W24->val.A[4] = W24->val.A[4] - 1;
                    if (W24->val.A[4] == -1)
                        part4 = null4;
                }
            }
        }
        while (!(part4 != nonstand4));
        temp = W24->next;
        if (part4 == stand4) {
            W24->next = newlist;
            newlist = (*list);
        } else
            dispchr (&(*list));
        (*list) = temp;
    }
    (*list) = newlist;
}

prodtype
pmodify (prodtype p1)
{
    prodtype p11;
    register int j;
    bool osptest;
    ocharptr temp;

    p11 = p1;
    osptest = false;
    for (j = 1; j <= nprod; j++) {
        if (currgrp.A[j - 1].name == ospnm)
            osptest = true;
    }
    while (p11 != NULL) {
        for (j = 1; j <= nprod; j++) {
            temp = p11->prods.A[j - 1];
            p11->prods.A[j - 1] = gmodify (temp, currgrp.A[j - 1]);
            odisp (&temp);
        }
        p11 = p11->next;
    }
    p11 = p1;
    if (osptest)
        p11 = prodexpand (p11);
    schur_psort (&p11, true);
    group = currgrp.A[1 - 1].name;
    return p11;
}

ocharptr
gmodify (ocharptr chrc, groop grp)
{
    ocharptr chrc2;

    chrc2 = chrccopy (chrc);
    group = grp.name;

    switch (grp.name) {
    case sung:
        gumodify (&chrc2, grp.rank);
        spclise (&chrc2, grp.rank);
        break;
    case un:
        if (grp.rank != 1)
            gumodify (&chrc2, grp.rank);
        else if (grp.rank == 1)
            u1modify (&chrc2);
        break;
    case unm:
    case sunm:
        unmmodify (&chrc2, grp.rank, grp.rank2);
        break;
    case unc:
        stdpq (&chrc2, grp.rank, grp.rank2, kkz);
        break;
    case son:
        omodify (&chrc2, grp.rank);
        if (!(grp.rank & 1))
            so2nexp (&chrc2, grp.rank);
        break;
    case on:
        omodify (&chrc2, grp.rank);
        break;
    case spn:
        if (grp.rank & 1) {
            print ("ERROR:rank of Sp(n) must be even\n");
            chrc2 = NULL;
        } else
            spmodify (&chrc2, grp.rank);
        break;
    case sonc:
        if ((grp.rank & 1)) {
            print ("ERROR:rank of SO^*(2n) must be even\n");
            chrc2 = NULL;
        } else
            modsonc (&chrc2, grp.rank); /*3/10/97 */
        break;
    case ospnm:
        ospmodify (&chrc2, grp.rank, grp.rank2);
        break;
    case spnc:
        if ((grp.rank & 1)) {
            print ("ERROR:rank of Sp(n,R) must be even\n");
            chrc2 = NULL;
        } else if (sb_conj == false)
            spncmodify (&chrc2, grp.rank);
        else
            modspnr (&chrc2, grp.rank);
        break;
    case sn:
        snselect (&chrc2, grp.rank, qspecial);
        break;
    case an:
        anselect (&chrc2, grp.rank);
        break;
    case g2:
        g2modify (&chrc2, 0);
        break;
    case f4:
        /*omodify(&chrc, 9); */
        f4modify (&chrc2, 0);
        break;
    case e6:
        e6modify (&chrc2, 0);
        break;
    case e7:
        e7modify (&chrc2, 0);
        break;
    case e8:
        e8modify (&chrc2, 0);
        break;
    case mp:
        break;
    case l168:
        break;
    case nill:
        warn ("group not set:no", cont);
        inform (" action taken", cr);
        break;
    default:
        Caseerror (Line);
    }
    return chrc2;
}

void
anexp (ocharptr * list)
{
    ocharptr lscan, addin;

    lscan = *list;

    if (lscan->lab == ' ') {
        lscan->lab = '-';
        cnu (&addin);
        addin->val = lscan->val;
        addin->mult = lscan->mult;
        addin->spin = lscan->spin;
        addin->next = lscan->next;
        addin->C6_double = lscan->C6_double;
        addin->lab = '+';
        lscan->next = addin;
        lscan = addin->next;
    }
}

void
anselect (ocharptr * chrc, int n)
{
    ocharptr head, lastptr, temp;
    bool sj;
    frame tf;
    int pts;

    ostndise (chrc);
    cnu (&head);
    lastptr = head;
    head->next = *chrc;
    while (*chrc != NULL) {
        register ocharptr W31 = &(*(*chrc));

        if ((!W31->spin))
            if ((wtfrm (&W31->val) != n) || W31->C6_double) {
                temp = (*chrc);
                lastptr->next = W31->next;
                (*chrc) = W31->next;
                dispchr (&temp);
            } else {
                sj = sconjgte (&(*chrc)->val);
                if ((sj && ((*chrc)->lab == ' ')))
                    anexp (&(*chrc));
                else {
                    tf = (*chrc)->val;
                    conjgte (&tf);
                    if (tf.A[1] > (*chrc)->val.A[1])
                        (*chrc)->val = tf;
                }
                lastptr = (*chrc);
                (*chrc) = W31->next;
        } else if ((W31->spin)) {
            spinmod (&(*chrc), n);
            pts = n - len (&W31->val);
            if (pts & 1)
                anexp (chrc);
            lastptr = (*chrc);
            *chrc = W31->next;
        }
    }
    *chrc = head->next;
    dispchr (&head);
    osort (chrc, true);
}

ocharptr
fusion (ocharptr chrc, groop grp, int level)
{
    register ocharptr result = NULL;    //was not initialized FB 2005

    chrc = chrccopy (chrc);
    group = grp.name;
    switch (grp.name) {
    case un:
        sunfusion (&chrc, level);
        result = chrc;
        break;
    case spn:
        spnfusion (&chrc, level);
        result = chrc;
        break;
    case an:
    case sn:
    case mp:
    case spnc:
    case sonc:
    case unm:
    case sunm:
    case ospnm:
        print ("Inappropriate group\n");
        break;
    case son:
    case sung:
    case g2:
    case f4:
    case e6:
    case e7:
    case e8:
    case en:
        print ("Fusion rule for group not implemented\n");
        break;
    default:
        Caseerror (Line);
    }
    return result;
}

void
spncmodify (ocharptr * list, int n)
{
    ocharptr newlist, temp;
    int kk, k, nn, parts, hook;
    register int i;
    char tag;
    register ocharptr tlist;

    tlist = *list;
    nn = n / 2;
    newlist = NULL;
    while (tlist != NULL) {
        k = tlist->val.A[maxdim];
        if (tlist->spin)
            kk = 2 * k + 1;
        else
            kk = 2 * k;
        cnu (&temp);
        temp->mult = tlist->mult;
        temp->val = tlist->val;
        temp->val.A[maxdim] = 0;
        temp->lab = ' ';
        temp->spin = false;
        tag = ' ';
        ostndise (&temp);
        group = on;
        omodify (&temp, kk);
        if (temp != NULL) {
            tag = temp->lab;
            temp->spin = false;
            gumodify (&temp, nn);
            group = spnc;
            temp->val.A[maxdim] = k;
            temp->conval.A[1] = k;
            temp->C6_double = tlist->C6_double;
            temp->spin = tlist->spin;
            temp->lab = tag;
            if ((tag == '#')) {
                parts = len (&temp->val);
                temp->lab = ' ';
                hook = kk - 2 * parts;
                if (hook != 0)
                    for (i = 1; i <= hook; i++) {
                        temp->val.A[parts + i] = 1;
                    }
            }
            cmerge (&newlist, &temp, true, true);
        }
        tlist = tlist->next;
    }
    odisp (list);
    osort (&newlist, true);
    *list = newlist;
}

void
omodify (ocharptr * list, int n)
{
    enum { unmodified5, modified5, nill5 } partition;
    ocharptr newlist, temp;
    int parts, depth, hookpart, hook;
    register ocharptr W42;

    ostndise (list);
    newlist = NULL;
    while ((*list) != NULL) {
        W42 = &(*(*list));

        if ((group == son))
            if (W42->lab == '#')
                W42->lab = ' ';
        progress ();
        parts = 0;
        while (W42->val.A[parts + 1] != 0)
            parts = parts + 1;
        if (W42->spin)
            hook = 2 * parts - n - 1;
        else
            hook = 2 * parts - n;
        depth = 1;
        if (W42->spin && (hook == 0))
            partition = nill5;
        else
            partition = modified5;
        if (hook > 0)
            do {
                partition = unmodified5;
                hookpart = W42->val.A[parts] - depth + 1;
                if (hook < hookpart)
                    partition = nill5;
                else if (hook == hookpart) {
                    if ((((W42->spin && !(W42->lab == '"'))
                          || (!W42->spin && (W42->lab == '"')))
                         && (bool) ((W42->val.A[parts]) & 1))
                        || (((!W42->spin && !(W42->lab == '"'))
                             || (W42->spin && (W42->lab == '"')))
                            && !(bool) ((W42->val.A[parts]) & 1)))
                        W42->mult = -W42->mult;
                    W42->val.A[parts] = W42->val.A[parts] - hookpart;
                    if (W42->lab == '+')
                        W42->lab = '-';
                    else if (W42->lab == '-')
                        W42->lab = '+';
                    if ((group == on))
                        if (W42->lab == ' ')
                            W42->lab = '#';
                        else
                            W42->lab = ' ';
                    partition = modified5;
                } else {
                    depth = W42->val.A[parts];
                    hook = hook - hookpart;
                    W42->val.A[parts] = W42->val.A[parts] - hookpart;
                    parts = parts - 1;
                    if (parts < 1)
                        partition = nill5;
                }
            }
            while (!(partition != unmodified5));
        if ((group == son) && ((bool) ((n) & 1)))
            W42->lab = ' ';
        parts = 0;
        while (W42->val.A[parts + 1] != 0)
            parts = parts + 1;
        if (((group == on) && (n % 2 == 0)))
            if ((parts == n / 2) || W42->spin)
                W42->lab = ' ';
        if ((parts <= n / 2) || (partition != modified5)) {
            temp = W42->next;
            if (partition == modified5) {
                W42->next = newlist;
                newlist = (*list);
            } else
                dispchr (&(*list));
            (*list) = temp;
        }
    }
    osort (&newlist, true);
    (*list) = newlist;
}

void
u1modify (ocharptr * list)
{
    ocharptr newlist, temp;
    register ocharptr W43;

    newlist = NULL;
    while ((*list) != NULL) {
        W43 = &(*(*list));

        temp = W43->next;
        if ((W43->val.A[2] == 0)) {
            W43->next = newlist;
            newlist = (*list);
        } else
            dispchr (&(*list));
        (*list) = temp;
    }
    osort (&newlist, true);
    (*list) = newlist;
}

void
gumodify (ocharptr * list, int n)
{
    enum { unmodified6, modified6, null6 } charactr;
    ocharptr newlist, temp;
    int parts1, parts2, depth1, depth2, hookpart, hook, hookb, x, y = 0, pts;
    register ocharptr W47;

    newlist = NULL;
    ostndise (list);
    while ((*list) != NULL) {
        W47 = &(*(*list));

        if (W47->C6_double == false) {
            pts = 0;
            while (((pts < (maxdim - 1)) && (W47->val.A[pts + 1] != 0)))
                pts = pts + 1;
            if (pts > n)
                charactr = null6;
            else
                charactr = modified6;
        } else {
            parts1 = 0;
            parts2 = 0;
            while (W47->conval.A[parts1 + 1] != 0)
                parts1 = parts1 + 1;
            while (W47->val.A[parts2 + 1] != 0)
                parts2 = parts2 + 1;
            hook = parts1 + parts2 - n - 1;
            depth1 = 1;
            depth2 = 1;
            if (hook == 0)
                charactr = null6;
            else
                charactr = modified6;
            if (hook > 0) {
                hookb = hook;
                do {
                    charactr = unmodified6;
                    hookpart = W47->conval.A[parts1] - depth1 + 1;
                    if (hookb < hookpart)
                        charactr = null6;
                    else if (hookb == hookpart) {
                        y = W47->conval.A[parts1];
                        W47->conval.A[parts1] =
                            W47->conval.A[parts1] - hookpart;
                        charactr = modified6;
                    } else {
                        depth1 = W47->conval.A[parts1];
                        hookb = hookb - hookpart;
                        W47->conval.A[parts1] =
                            W47->conval.A[parts1] - hookpart;
                        parts1 = parts1 - 1;
                        if (parts1 < 1)
                            charactr = null6;
                    }
                }
                while (!(charactr != unmodified6));
                hookb = hook;
                if (charactr == modified6)
                    do {
                        charactr = unmodified6;
                        hookpart = W47->val.A[parts2] - depth2 + 1;
                        if (hookb < hookpart)
                            charactr = null6;
                        else if (hookb == hookpart) {
                            charactr = modified6;
                            x = W47->val.A[parts2];
                            W47->val.A[parts2] =
                                W47->val.A[parts2] - hookpart;
                        } else {
                            depth2 = W47->val.A[parts2];
                            hookb = hookb - hookpart;
                            W47->val.A[parts2] =
                                W47->val.A[parts2] - hookpart;
                            parts2 = parts2 - 1;
                            if (parts2 < 1)
                                charactr = null6;
                        }
                    }
                    while (!(charactr != unmodified6));
                if (charactr == modified6)
                    if (!(bool) ((x + y) & 1))
                        W47->mult = -W47->mult;
                parts1 = 0;
                parts2 = 0;
                while (W47->conval.A[parts1 + 1] != 0)
                    parts1 = parts1 + 1;
                while (W47->val.A[parts2 + 1] != 0)
                    parts2 = parts2 + 1;
            }
            if (W47->conval.A[1] == 0)
                W47->C6_double = false;
            pts = parts1 + parts2;
        }
        if ((pts <= n) || (charactr == null6)) {
            temp = W47->next;
            if (charactr == modified6) {
                W47->next = newlist;
                newlist = (*list);
            } else
                dispchr (&(*list));
            (*list) = temp;
        }
    }
    osort (&newlist, true);
    (*list) = newlist;
}

void
spmodify (ocharptr * list, int n)
{
    enum { unmodified7, modified7, null7 } partition;
    ocharptr newlist, temp;
    int parts, depth, hookpart, hook;
    register ocharptr tempo;

    ostndise (list);
    newlist = NULL;
    while (*list != NULL) {
        tempo = *list;

        progress ();
        parts = 0;
        while (tempo->val.A[parts + 1] != 0)
            parts = parts + 1;
        hook = 2 * parts - n - 2;
        depth = 1;
        if (hook == 0)
            partition = null7;
        else
            partition = modified7;
        if (hook > 0)
            do {
                partition = unmodified7;
                hookpart = tempo->val.A[parts] - depth + 1;
                if (hook < hookpart)
                    partition = null7;
                else if (hook == hookpart) {
                    if ((tempo->val.A[parts]) & 1)
                        tempo->mult = -tempo->mult;
                    tempo->val.A[parts] = tempo->val.A[parts] - hookpart;
                    partition = modified7;
                } else {
                    depth = tempo->val.A[parts];
                    hook = hook - hookpart;
                    tempo->val.A[parts] = tempo->val.A[parts] - hookpart;
                    parts = parts - 1;
                    if (parts < 1)
                        partition = null7;
                }
            }
            while (partition == unmodified7);
        parts = 0;
        while (tempo->val.A[parts + 1] != 0)
            parts = parts + 1;
        if ((parts <= n / 2) || (partition != modified7)) {
            temp = tempo->next;
            if (partition == modified7) {
                tempo->next = newlist;
                newlist = *list;
            } else
                dispchr (list);
            *list = temp;
        }
    }
    osort (&newlist, true);
    *list = newlist;
}

void
spnfusion (ocharptr * list, int level)
{
    ocharptr newlist, temp;
    register ocharptr W52;

    level = 2 * level;
    newlist = (*list);
    while (newlist != NULL) {
        W52 = &(*newlist);

        progress ();
        conjgte (&W52->val);
        newlist = W52->next;
    }
    spmodify (&(*list), level);
    temp = (*list);
    while (temp != NULL) {
        register ocharptr W53 = &(*temp);

        conjgte (&W53->val);
        temp = W53->next;
    }
    osort (&(*list), true);
}

void
sunfusion (ocharptr * list, int level)
{
    ocharptr newlist, temp;
    register ocharptr W54;
    register ocharptr W55;

    newlist = (*list);
    while (newlist != NULL) {
        W54 = &(*newlist);

        progress ();
        conjgte (&W54->val);
        if (W54->C6_double)
            conjgte (&W54->conval);
        newlist = W54->next;
    }
    gumodify (&(*list), level);
    temp = (*list);
    while (temp != NULL) {
        W55 = &(*temp);

        conjgte (&W55->val);
        if (W55->C6_double)
            conjgte (&W55->conval);
        temp = W55->next;
    }
    osort (&(*list), true);
}

void
spclise (ocharptr * list, int n)
{
    ocharptr newlist;
    register int i;
    register int j;
    register ocharptr W56;

    newlist = *list;
    while (newlist != NULL) {
        W56 = newlist;

        progress ();
        if (W56->C6_double && (n > maxdim))
            print ("Error maxdim = %10d exceeded\n", maxdim);
        if (W56->C6_double && (W56->conval.A[1] != 0)) {
            i = 0;
            while (W56->conval.A[i + 1] != 0)
                i = i + 1;
            if (n > maxdim)
                print ("Error maxdim = %10d exceeded\n", maxdim);
            do {
                if (n <= maxdim)
                    for (j = 1; j <= n - i; j++) {
                        W56->val.A[j] =
                            W56->val.A[j] + W56->conval.A[i] -
                            W56->conval.A[i + 1];
                    }
                j = i;
                do {
                    i = i - 1;
                }
                while (!((W56->conval.A[i] != W56->conval.A[j]) || (i == 0)));
            }
            while (!(i == 0));
            W56->C6_double = false;
        } else {
            W56->C6_double = false;
            if (n <= maxdim)
                for (i = 1; i <= n; i++) {
                    W56->val.A[i] = W56->val.A[i] - W56->val.A[n];
                }
        }
        newlist = W56->next;
    }
    osort (&(*list), true);
}

void
snexp (ocharptr * list, int n)
{
    ocharptr lscan, addin;
    int pts;
    register ocharptr W61;

    lscan = (*list);
    while (lscan != NULL) {
        W61 = &(*lscan);

        pts = 0;
        while (W61->val.A[pts + 1] != 0)
            pts = pts + 1;
        if ((W61->spin && (bool) ((n - pts + 1) & 1))) {
            if ((W61->lab == ' ')) {
                W61->lab = '-';
                cnu (&addin);
                addin->val = W61->val;
                addin->mult = W61->mult;
                addin->spin = W61->spin;
                addin->next = W61->next;
                addin->C6_double = W61->C6_double;
                addin->lab = '+';
                W61->next = addin;
                lscan = addin->next;
            } else
                lscan = W61->next;
        } else
            lscan = W61->next;
    }
}

void
spinmod (ocharptr * chrc, int n)
{
    termptr temp;
    int k;
    register int i;

    snu (&temp);
    temp->mult = (*chrc)->mult;
    for (i = 1; i <= maxdim; i++) {
        temp->val.A[i] = 0;
    }
    k = qlen ((*chrc)->val) + 1;
    temp->val.A[1] = n - wtfrm (&(*chrc)->val);
    for (i = 2; i <= k; i++) {
        temp->val.A[i] = (*chrc)->val.A[i - 1];
    }
    qstndise (&temp);
    if (temp != NULL) {
        (*chrc)->mult = temp->mult;
        for (i = 2; i <= maxdim; i++) {
            (*chrc)->val.A[i - 1] = temp->val.A[i];
        }
    } else
        (*chrc)->mult = 0;
    dispsfn (&temp);
}

void
snselect (ocharptr * chrc, int n, bool qspecials)
{
    ocharptr head, lastptr;
    register ocharptr W68;

    ostndise (chrc);
    cnu (&head);
    lastptr = head;
    head->next = (*chrc);
    while ((*chrc) != NULL) {
        W68 = &(*(*chrc));

        if ((!W68->spin))
            if ((wtfrm (&W68->val) != n) || W68->C6_double) {
                lastptr->next = W68->next;
                (*chrc) = W68->next;
            } else {
                lastptr = (*chrc);
                (*chrc) = W68->next;
        } else if ((W68->spin)) {
            spinmod (&(*chrc), n);
            lastptr = (*chrc);
            (*chrc) = W68->next;
        }
    }
    (*chrc) = head->next;
    dispchr (&head);
    if (qspecials)
        snexp (&(*chrc), n);
    osort (&(*chrc), false);
}

void
so2nexp (ocharptr * list, int n)
{
    ocharptr lscan, addin;
    int pts;
    register ocharptr W69;

    lscan = (*list);
    while (lscan != NULL) {
        W69 = &(*lscan);

        pts = 0;
        while (W69->val.A[pts + 1] != 0)
            pts = pts + 1;
        if (((pts == n / 2) || W69->spin) && (W69->lab == ' ')) {
            W69->lab = '-';
            cnu (&addin);
            addin->val = W69->val;
            addin->mult = W69->mult;
            addin->spin = W69->spin;
            addin->next = W69->next;
            addin->C6_double = W69->C6_double;
            addin->lab = '+';
            W69->next = addin;
            lscan = addin->next;
        } else
            lscan = W69->next;
    }
}

void
dehook (frame * partn, frame * frm, hktype kind, int *hk, int *r,
        int *c, int *row)
{
    int gap;

    (*frm) = (*partn);
    if ((frm->A[1] == 0) || (frm->A[1] == (-99)))
        frm->A[1] = (-99);
    else if (kind == rc) {
        (*row) = (*r);
        while (frm->A[(*row) + 1] >= (*c)) {
            frm->A[(*row)] = frm->A[(*row) + 1] - 1;
            (*row) = (*row) + 1;
        }
        frm->A[(*row)] = (*c) - 1;
    } else {
        if (kind == ccol) {
            conjgte (&(*frm));
            (*row) = (*c);
        } else
            (*row) = (*r);
        do {
            if (frm->A[(*row) + 1] > 0)
                gap = frm->A[(*row)] - frm->A[(*row) + 1] + 1;
            else if ((*hk) <= frm->A[(*row)])
                gap = (*hk) + 1;
            else
                gap = (*hk);
            if ((*hk) < gap) {
                frm->A[(*row)] = frm->A[(*row)] - (*hk);
                (*hk) = 0;
            } else if ((*hk) == gap) {
                frm->A[1] = -99;
                (*hk) = 0;
            } else {
                (*hk) = (*hk) - gap;
                frm->A[(*row)] = frm->A[(*row)] - gap;
                (*row) = (*row) + 1;
            }
        }
        while (!(((*hk) == 0) || ((*row) >= maxdim)));
        if ((kind == ccol) && (frm->A[1] != (-99)))
            conjgte (&(*frm));
    }
}

void
so2nmod (ocharptr * list, char labyl, int n)
{
    enum { unmodified8, modified8, null8 } charactr;
    ocharptr newlist, temp;
    int parts, q, i, x;
    register ocharptr W73;

    newlist = NULL;
    while ((*list) != NULL) {
        W73 = &(*(*list));

        progress ();
        W73->lab = labyl;
        if (W73->C6_double == false)
            charactr = modified8;
        else if (W73->conval.A[1] != 0)
            charactr = unmodified8;
        else {
            W73->C6_double = false;
            charactr = modified8;
        }
        while (W73->C6_double) {
            parts = 0;
            while (W73->conval.A[parts + 1] != 0)
                parts = parts + 1;
            if (W73->spin)
                q = W73->conval.A[1] - parts;
            else
                q = W73->conval.A[1] - parts + 1;
            if (q < 0)
                charactr = null8;
            else {
                i = 1;
                while (W73->conval.A[i + 1] != 0) {
                    W73->conval.A[i] = W73->conval.A[i + 1] - 1;
                    i = i + 1;
                }
                W73->conval.A[i] = 0;
                x = (n / 2) - i + 1;
                do {
                    while ((q != 0)
                           && ((W73->val.A[x] != W73->val.A[x - 1] + 1)
                               || (x == 1))) {
                        W73->val.A[x] = W73->val.A[x] + 1;
                        q = q - 1;
                    }
                    if (q != 0)
                        x = x - 1;
                    else if ((W73->val.A[x] == W73->val.A[x - 1] + 1)
                             && (x != 1))
                        charactr = null8;
                    else
                        charactr = modified8;
                }
                while (!(charactr != unmodified8));
                if ((bool) ((n / 2 - x) & 1))
                    W73->mult = -W73->mult;
                if (W73->lab == '-')
                    W73->lab = '+';
                else
                    W73->lab = '-';
            }
            if ((W73->conval.A[1] == 0) || (charactr == null8))
                W73->C6_double = false;
        }
        temp = W73->next;
        if (charactr == modified8) {
            if (!W73->spin) {
                parts = 0;
                while (W73->val.A[parts + 1] != 0)
                    parts = parts + 1;
                if (parts < n / 2)
                    W73->lab = ' ';
            }
            W73->next = newlist;
            newlist = (*list);
        } else
            dispchr (&(*list));
        (*list) = temp;
    }
    osort (&newlist, false);
    (*list) = newlist;
}

void
atypcond (frame * partn, frame * atyp, int *m, int *mm, int *n,
          int *r1, int *c1, int *ac)
{
    int ai, aj;
    register int j;
    register int i;
    frame lambda, mu, cpartn;

    for (i = 0; i <= maxdim; i++) {
        lambda.A[i] = 0;
        mu.A[i] = 0;
        atyp->A[i] = 0;
    }
    for (i = 0; i <= *m; i++) {
        lambda.A[i] = partn->A[i] - (*n);
    }
    cpartn = (*partn);
    conjgte (&cpartn);
    for (i = 1; i <= *n; i++) {
        mu.A[i] = cpartn.A[i];
    }
    (*r1) = partn->A[(*m) + 1] - (*n);
    (*c1) = cpartn.A[(*n) + 1] - (*m);
    ai = 0;
    aj = 0;
    (*ac) = 0;
    for (i = ai + 1; i <= *n; i++) {
        for (j = aj + 1; j <= *m; j++) {
            if ((mu.A[i] + (*n) + j + 1) == (lambda.A[j] + (*mm) + i)) {
                atyp->A[2 * (*ac)] = i;
                atyp->A[2 * (*ac) + 1] = j;
                aj = j;
                (*ac) = (*ac) + 1;
            }
        }
    }
}

void
typical (frame * lambda, frame * lambdamu, int *m, int *n,
         int *q, int *r1, int *c1, int *r, int *c, int *signs)
{
    int nhk, row = *r;          /* row was uninitialized */

    (*lambdamu) = (*lambda);
    if ((*r1) > (*c1)) {
        nhk = 2 * (*r1) - 1 - (*q);
        dehook (&(*lambda), &(*lambdamu), rrow, &nhk, &(*r), &(*c), &row);
        if (abs ((row - (*m) - 1 + (*q)) % 2) == 1)
            (*signs) = (*signs) * (-1);
    } else {
        nhk = 2 * (*c1) - 1 + (*q);
        dehook (&(*lambda), &(*lambdamu), ccol, &nhk, &(*r), &(*c), &row);
        if (abs ((row - (*n) - 1) % 2) == 1)
            (*signs) = (*signs) * (-1);
    }
}

void
tree (frame * lval, frame * lambda, frame * lambdamu, frame * atyp,
      int level, int n, int multt, int multu, ocharptr * sublist)
{
    int sign1, dum1, dum2, i1, j1;
    register int j;
    frame lam, lamu;
    ocharptr entry, entry1;
    register ocharptr W86, W87;

    sign1 = 1;
    dum1 = 1;
    dum2 = 1;

    for (j = level; j >= 1; j--) {
        sign1 = 1;
        i1 = atyp->A[2 * (j - 1)];
        j1 = atyp->A[2 * j - 1];
        if (abs ((lval->A[j1] - i1 - 1) % 2) == 1)
            sign1 = sign1 * (-1);
        dehook (&(*lambda), &lam, rc, &dum1, &j1, &i1, &dum2);
        if (lambdamu->A[1] != (-99))
            dehook (&(*lambdamu), &lamu, rc, &dum1, &j1, &i1, &dum2);
        else
            lamu.A[1] = (-99);
        if (lam.A[1] != (-99)) {
            cnu (&entry);
            {
                W86 = &(*entry);

                W86->val = lam;
                W86->mult = multt * sign1 * (-1);
                W86->lab = ' ';
                W86->spin = false;
                W86->C6_double = false;
                W86->next = (*sublist);
            }
            (*sublist) = entry;
        }
        if (lamu.A[1] != (-99)) {
            cnu (&entry1);
            {
                W87 = &(*entry1);

                W87->val = lamu;
                W87->mult = multu * sign1;
                W87->lab = ' ';
                W87->spin = false;
                W87->C6_double = false;
                W87->next = (*sublist);
            }
            (*sublist) = entry1;
        }
        tree (&(*lval), &lam, &lamu, &(*atyp), j - 1, n, multt * sign1,
              multu * sign1, &(*sublist));
    }
}

void
unmmodify (ocharptr * list, int n, int m)
{
    ocharptr newlist, temp;
    bool modified;
    register ocharptr W88;

    newlist = NULL;
    while ((*list) != NULL) {
        W88 = &(*(*list));

        progress ();
        modified = (bool) (W88->val.A[n + 1] <= m);
        temp = W88->next;
        if (modified) {
            W88->next = newlist;
            newlist = (*list);
        } else
            dispchr (&(*list));
        (*list) = temp;
    }
    osort (&newlist, true);
    (*list) = newlist;
}

void
ospmodify (ocharptr * list, int mm, int nn)
{
    int m, n, q, signs, r, c, r1, c1, ac;
    frame atyp, valmu, lval;
    ocharptr sublist, help, top;
    bool modified;

    m = mm / 2;
    q = 1 - mm + 2 * m;
    n = nn / 2;
    sublist = NULL;
    r = m + 1;
    c = n + 1;
    help = NULL;
    top = (*list);
    do {
        modified = true;
        (*list) = top;
        while ((*list) != NULL) {
            help = (*list)->next;
            if ((*list)->val.A[m + 1] > n) {
                modified = false;
                atypcond (&(*list)->val, &atyp, &m, &mm, &n, &r1, &c1, &ac);
                signs = 1;
                typical (&(*list)->val, &valmu, &m, &n, &q, &r1, &c1, &r, &c,
                         &signs);
                if (atyp.A[0] > 0) {
                    lval = (*list)->val;
                    tree (&lval, &(*list)->val, &valmu, &atyp, ac, n,
                          (*list)->mult, (*list)->mult * signs, &sublist);
                }
                if (valmu.A[1] != (-99)) {
                    (*list)->val = valmu;
                    (*list)->mult = (*list)->mult * signs;
                    (*list)->next = sublist;
                    sublist = (*list);
                } else {
                    (*list)->next = NULL;
                }
                (*list) = help;
            } else {
                if (!modified) {
                    (*list)->next = sublist;
                    sublist = (*list);
                }
                (*list) = help;
            }
        }
        if (!modified)
            top = sublist;
        else
            (*list) = top;
        sublist = NULL;
        osort (&(*list), false);
    }
    while (!(modified == true));
}

void
racahg2 (ocharptr * list)
{
    ocharptr ptr;

    if ((ggroup.name == g2)) {
        ptr = (*list);
        while (ptr != NULL) {
            ptr->val.A[1] = ptr->val.A[1] - ptr->val.A[2];
            ptr = ptr->next;
        }
        osort (&(*list), false);        // anyway, if needed, it will be reverse ordered
    } else
        print ("mistake! group not set as G(2) \n");
}

void
fracahg2 (ocharptr * list)
{
    ocharptr ptr;

    if ((ggroup.name == g2)) {
        ptr = (*list);
        while (ptr != NULL) {
            ptr->val.A[1] = ptr->val.A[1] + ptr->val.A[2];
            ptr = ptr->next;
        }
        osort (&(*list), false);        // anyway, if needed, it will be reverse ordered
    } else
        print ("mistake! group not set as G(2) \n");
}

void
g2modify (ocharptr * chrc, int adjust)
{
    enum { standard9, nonstandard9, null9 } character;
    ocharptr newlist, temp;
    int h1, h2;
    register ocharptr W92;

    newlist = NULL;
    while ((*chrc) != NULL) {
        W92 = &(*(*chrc));

        character = nonstandard9;

        W92->val.A[1] = W92->val.A[1] - adjust;
        do {
            if (((W92->val.A[1] >= 2 * W92->val.A[2])
                 && (W92->val.A[2] >= 0)))
                character = standard9;
            else {
                if (W92->val.A[2] < 0) {
                    h1 = -W92->val.A[2] - 1;
                    if ((h1 == 0))
                        character = null9;
                    else if (h1 > 0) {
                        W92->mult = -W92->mult;
                        W92->val.A[1] = W92->val.A[1] + h1;
                        W92->val.A[2] = W92->val.A[2] + 2 * h1;
                    }
                } else {
                    h2 = -W92->val.A[1] + 2 * W92->val.A[2] - 1;
                    if (h2 == 0)
                        character = null9;
                    else if (h2 > 0) {
                        W92->mult = -W92->mult;
                        W92->val.A[2] = W92->val.A[2] - h2;
                    }
                }
            }
        }
        while (!((character != nonstandard9)));
        temp = W92->next;
        if (character == standard9) {
            W92->next = newlist;
            newlist = (*chrc);
        } else
            dispchr (&(*chrc));
        (*chrc) = temp;
    }
    (*chrc) = newlist;
    osort (&(*chrc), true);
}

void
e8modify (ocharptr * chrc, int adjust)
{
    enum { standx1, nonstandx1 } parte1;
    enum { standardx2, nonstandardx2, nullx2 } character;
    ocharptr templist, newlist, temp;
    int p, h1, h2, h;
    register int i;
    register ocharptr W98;

    newlist = NULL;
    templist = (*chrc);
    while (templist != NULL) {
        W98 = &(*templist);

        character = nonstandardx2;
        W98->val.A[1] = W98->val.A[1] - adjust;
        if ((wtfrm (&W98->val) % 3 != 0))
            character = nullx2;
        else {
            do {
                do {
                    i = 2;
                    parte1 = standx1;
                    do {
                        h = W98->val.A[i + 1] - W98->val.A[i] - 1;
                        if (h == 0)
                            character = nullx2;
                        else if (h > 0) {
                            parte1 = nonstandx1;
                            W98->mult = -W98->mult;
                            W98->val.A[i] = W98->val.A[i] + h;
                            W98->val.A[i + 1] = W98->val.A[i + 1] - h;
                        }
                        i = i + 1;
                    }
                    while (!((i == 8) || (character == nullx2)));
                }
                while (!((parte1 != nonstandx1) || (character == nullx2)));
                if ((character != nullx2)) {
                    p = 2 * W98->val.A[2] + 2 * W98->val.A[3] +
                        2 * W98->val.A[4] - W98->val.A[5] - W98->val.A[6] -
                        W98->val.A[7] - W98->val.A[8];
                    if (((W98->val.A[1] >= p) && (p >= 0))) {
                        character = standardx2;
                        if (W98->val.A[8] < 0)
                            character = nonstandardx2;
                        for (i = 1; i <= 8; i++) {
                            if (W98->val.A[i] < W98->val.A[i + 1])
                                character = nonstandardx2;
                        }
                    }
                    if ((character == nonstandardx2)) {
                        if (W98->val.A[8] < 0) {
                            h1 = -W98->val.A[8] - 1;
                            if ((h1 == 0))
                                character = nullx2;
                            else if (h1 > 0) {
                                W98->mult = -W98->mult;
                                for (i = 1; i <= 8; i++) {
                                    W98->val.A[i] = W98->val.A[i] + h1;
                                }
                                W98->val.A[8] = W98->val.A[8] + h1;
                            }
                        } else {
                            h2 = (-W98->val.A[1] + 2 * W98->val.A[2] +
                                  2 * W98->val.A[3] + 2 * W98->val.A[4] -
                                  W98->val.A[5] - W98->val.A[6] -
                                  W98->val.A[7] - W98->val.A[8] - 3) / 3;
                            if (h2 == 0)
                                character = nullx2;
                            else if (h2 > 0) {
                                W98->mult = -W98->mult;
                                W98->val.A[2] = W98->val.A[2] - h2;
                                W98->val.A[3] = W98->val.A[3] - h2;
                                W98->val.A[4] = W98->val.A[4] - h2;
                            }
                        }
                    }
                }
            }
            while (!((character != nonstandardx2)));
        }
        temp = W98->next;
        if (character == standardx2) {
            W98->next = newlist;
            newlist = templist;
        } else
            dispchr (&templist);
        templist = temp;
    }
    templist = newlist;
    (*chrc) = templist;
    osort (&(*chrc), true);
}

void
e7modify (ocharptr * chrc, int adjust)
{
    enum { stande7, nonstande7 } parte2;

    enum { standarde7, nonstandarde7, nulle7 } character;
    ocharptr templist, newlist, temp;
    int p, h1, h2, h;
    register int i;
    register ocharptr W108;

    newlist = NULL;
    templist = (*chrc);
    while (templist != NULL) {
        W108 = &(*templist);

        if ((bool) ((wtfrm (&W108->val)) & 1))
            character = nulle7;
        else {
            character = nonstandarde7;
            W108->val.A[1] = W108->val.A[1] - adjust;
            do {
                do {
                    i = 2;
                    parte2 = stande7;
                    do {
                        h = W108->val.A[i + 1] - W108->val.A[i] - 1;
                        if (h == 0)
                            character = nulle7;
                        else if (h > 0) {
                            parte2 = nonstande7;
                            W108->mult = -W108->mult;
                            W108->val.A[i] = W108->val.A[i] + h;
                            W108->val.A[i + 1] = W108->val.A[i + 1] - h;
                        }
                        i = i + 1;
                    }
                    while (!((i == 7) || (character == nulle7)));
                }
                while (!((parte2 != nonstande7) || (character == nulle7)));
                if ((character != nulle7)) {
                    p = W108->val.A[2] + W108->val.A[3] + W108->val.A[4] +
                        W108->val.A[5] - W108->val.A[6] - W108->val.A[7];
                    if (((W108->val.A[1] >= p) && (p >= 0))) {
                        character = standarde7;
                        if (W108->val.A[7] < 0)
                            character = nonstandarde7;
                        for (i = 1; i <= 7; i++) {
                            if (W108->val.A[i] < W108->val.A[i + 1])
                                character = nonstandarde7;
                        }
                    }
                    if ((character == nonstandarde7)) {
                        if (W108->val.A[7] < 0) {
                            h1 = -W108->val.A[7] - 1;
                            if ((h1 == 0))
                                character = nulle7;
                            else if (h1 > 0) {
                                W108->mult = -W108->mult;
                                for (i = 1; i <= 7; i++) {
                                    W108->val.A[i] = W108->val.A[i] + h1;
                                }
                                W108->val.A[7] = W108->val.A[7] + h1;
                            }
                        } else {
                            h2 = (-W108->val.A[1] + W108->val.A[2] +
                                  W108->val.A[3] + W108->val.A[4] +
                                  W108->val.A[5] - W108->val.A[6] -
                                  W108->val.A[7] - 2) / 2;
                            if (h2 == 0)
                                character = nulle7;
                            else if (h2 > 0) {
                                W108->mult = -W108->mult;
                                for (i = 2; i <= 5; i++) {
                                    W108->val.A[i] = W108->val.A[i] - h2;
                                }
                            }
                        }
                    }
                }
            }
            while (!((character != nonstandarde7)));
        }
        temp = W108->next;
        if (character == standarde7) {
            W108->next = newlist;
            newlist = templist;
        } else
            dispchr (&templist);
        templist = temp;
    }
    templist = newlist;
    (*chrc) = templist;
    osort (&(*chrc), true);
}

void
e6modify (ocharptr * chrc, int adjust)
{
    enum { stande6, nonstande6 } parte3;

    enum { standarde6, nonstandarde6, nulle6 } character;

    ocharptr templist, newlist, temp;
    int p, h1, h2, h;
    register int i;
    register ocharptr W120;

    newlist = NULL;
    templist = (*chrc);
    while (templist != NULL) {
        W120 = &(*templist);

        character = nonstandarde6;
        W120->val.A[1] = W120->val.A[1] - adjust;
        if ((bool) ((wtfrm (&W120->val)) & 1))
            character = nulle6;
        else {
            do {
                do {
                    i = 2;
                    parte3 = stande6;
                    do {
                        h = W120->val.A[i + 1] - W120->val.A[i] - 1;
                        if (h == 0)
                            character = nulle6;
                        else if (h > 0) {
                            parte3 = nonstande6;
                            W120->mult = -W120->mult;
                            W120->val.A[i] = W120->val.A[i] + h;
                            W120->val.A[i + 1] = W120->val.A[i + 1] - h;
                        }
                        i = i + 1;
                    }
                    while (!((i == 6) || (character == nulle6)));
                }
                while (!((parte3 != nonstande6) || (character == nulle6)));
                if ((character != nulle6)) {
                    p = W120->val.A[2] + W120->val.A[3] + W120->val.A[4] -
                        W120->val.A[5] - W120->val.A[6];
                    if (((W120->val.A[1] >= p) && (p >= 0))) {
                        character = standarde6;
                        if (W120->val.A[6] < 0)
                            character = nonstandarde6;
                        for (i = 2; i <= 6; i++) {
                            if (W120->val.A[i] < W120->val.A[i + 1])
                                character = nonstandarde6;
                        }
                    }
                    if ((character == nonstandarde6)) {
                        if (W120->val.A[6] < 0) {
                            h1 = -W120->val.A[6] - 1;
                            if ((h1 == 0))
                                character = nulle6;
                            else if (h1 > 0) {
                                W120->mult = -W120->mult;
                                for (i = 2; i <= 6; i++) {
                                    W120->val.A[i] = W120->val.A[i] + h1;
                                }
                                W120->val.A[6] = W120->val.A[6] + h1;
                            }
                        } else {
                            h2 = (-W120->val.A[1] + W120->val.A[2] +
                                  W120->val.A[3] + W120->val.A[4] -
                                  W120->val.A[5] - W120->val.A[6] - 2) / 2;
                            if (h2 == 0)
                                character = nulle6;
                            else if (h2 > 0) {
                                W120->mult = -W120->mult;
                                for (i = 2; i <= 4; i++) {
                                    W120->val.A[i] = W120->val.A[i] - h2;
                                }
                                W120->val.A[1] = W120->val.A[1] + h2;
                            }
                        }
                    }
                }
            }
            while (!((character != nonstandarde6)));
        }
        temp = W120->next;
        if (character == standarde6) {
            W120->next = newlist;
            newlist = templist;
        } else
            dispchr (&templist);
        templist = temp;
    }
    templist = newlist;
    (*chrc) = templist;
    osort (&(*chrc), true);
}

void
f4modify (ocharptr * chrc, int adjust)
{
    enum { standf4, nonstandf4 } partf1;

    enum { standardf4, nonstandardf4, nullf4 } character;

    ocharptr templist, newlist, temp;
    int p, h1, h2, h;
    register int i;
    register ocharptr W132;

    newlist = NULL;
    templist = (*chrc);
    while (templist != NULL) {
        W132 = &(*templist);

        character = nonstandardf4;
        W132->val.A[1] = W132->val.A[1] - adjust;
        do {
            do {
                i = 2;
                partf1 = standf4;
                do {
                    h = W132->val.A[i + 1] - W132->val.A[i] - 1;
                    if (h == 0)
                        character = nullf4;
                    else if (h > 0) {
                        partf1 = nonstandf4;
                        W132->mult = -W132->mult;
                        W132->val.A[i] = W132->val.A[i] + h;
                        W132->val.A[i + 1] = W132->val.A[i + 1] - h;
                    }
                    i = i + 1;
                }
                while (!((i == 4) || (character == nullf4)));
            }
            while (!((partf1 != nonstandf4) || (character == nullf4)));
            if (character != nullf4) {
                p = W132->val.A[2] + W132->val.A[3] + W132->val.A[4];
                if ((((W132->spin && (W132->val.A[1] > p)))
                     || ((!W132->spin && (W132->val.A[1] >= p)) && (p >= 0)))) {
                    character = standardf4;
                    for (i = 2; i <= 4; i++) {
                        if ((W132->val.A[i] < W132->val.A[i + 1]))
                            character = nonstandardf4;
                    }
                    if (W132->val.A[4] < 0)
                        character = nonstandardf4;
                }
                if (character == nonstandardf4) {
                    if ((W132->val.A[4] < 0)) {
                        h1 = -2 * W132->val.A[4] - 1;
                        if (W132->spin)
                            h1 = h1 - 1;
                        if (h1 == 0)
                            character = nullf4;
                        else {
                            W132->mult = -W132->mult;
                            W132->val.A[4] = W132->val.A[4] + h1;
                        }
                    } else {
                        h2 = -W132->val.A[1] + W132->val.A[2] +
                            W132->val.A[3] + W132->val.A[4];
                        if (!W132->spin)
                            h2 = h2 - 1;
                        if (h2 == 0)
                            character = nullf4;
                        else if (h2 > 0) {
                            if ((bool) ((h2) & 1)) {
                                W132->spin = (bool) (!W132->spin);
                                if (W132->spin) {
                                    W132->val.A[2] = W132->val.A[2] - 1;
                                    W132->val.A[3] = W132->val.A[3] - 1;
                                    W132->val.A[4] = W132->val.A[4] - 1;
                                } else
                                    W132->val.A[1] = W132->val.A[1] + 1;
                            }
                            h2 = h2 / 2;
                            W132->mult = -W132->mult;
                            W132->val.A[1] = W132->val.A[1] + h2;
                            W132->val.A[2] = W132->val.A[2] - h2;
                            W132->val.A[3] = W132->val.A[3] - h2;
                            W132->val.A[4] = W132->val.A[4] - h2;
                        }
                    }
                }
            }
        }
        while (!((character != nonstandardf4)));
        temp = W132->next;
        if (character == standardf4) {
            W132->next = newlist;
            newlist = templist;
        } else
            dispchr (&templist);
        templist = temp;
    }
    templist = newlist;
    (*chrc) = templist;
    osort (&(*chrc), true);
}

void
modspnsfn (termptr * list, int k)
{
    termptr xlist, newlist, tlist;
    register termptr W133;

    xlist = (*list);
    newlist = NULL;
    while (xlist != NULL) {
        W133 = &(*xlist);

        snu (&tlist);
        tlist->mult = W133->mult;
        tlist->val = W133->val;
        tlist->slab = ' ';
        conjgte (&tlist->val);
        if (((tlist->val.A[1] + tlist->val.A[2]) <= k)) {
            tlist->val = W133->val;
            add (&newlist, &tlist);
        }
        dispsfn (&tlist);
        xlist = W133->next;
    }
    ldisp (&(*list));
    (*list) = newlist;
}
