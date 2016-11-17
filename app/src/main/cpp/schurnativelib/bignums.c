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

/** \file bignums.c 
 * This file contains all the routines used to handle large numbers
 * large numbers are required
 * for calculating dimensions, Casimir eigenvalues and Dynkin indexes
 *
 * NB: Schur does NOT treat multiplicities as big numbers and hence overflow will
 * occur if INT_MAX is exceeded. This is extremely rare and will not be
 * encountered in normal use of Schur*/


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
#include "g.h"
#include "bignums.h"

int *G159_rank;

/** Writes out big numbers 
 * \param fyle may be screen or other
 * \param p return index of display
 * \param bigno bignum to write out*/
void
wrtbigno (text * fyle, int *p, bframe * bigno)
{
    int maxw;
    register int j;

    maxw = maxb;
    while ((bigno->A[maxw] == 0) && (maxw > 0))
        maxw = maxw - 1;
    (*p) = (*p) + 2 * maxw;
    if ((maxw == maxb))
        fprintf (output.fp, "error - number too large?\n"), Putl (output, 1);
    if (bigno->A[maxw] < 0)
        Putchr ('-', (*fyle));
    if (maxw == 0)
        Putchr ('0', (*fyle));
    else

        for (j = maxw; j >= 1; j--) {
            if ((bigno->A[j] / 10 == 0) && (j != maxw))
                Putchr ('0', (*fyle));
            (void) fprintf ((*fyle).fp, "%1d", abs (bigno->A[j])),
                Putl ((*fyle), 0);
        }
}

/** Determines the size of the list */
int
dsize (termptr list)
{
    int j, abmult, jm;
    register int i;
    termptr temp = list;

    j = 0;
    for (i = 1; i <= (*G159_rank); i++) {
        if (temp->val.A[i] > 9)
            j = j + 1;
    }
    abmult = abs (temp->mult);
    jm = 1;
    do {
        abmult = abmult / 10;
        if (abmult > 0)
            jm = jm + 1;
    }
    while (!(abmult == 0));
    if (abs (temp->mult) == 1)
        jm = jm - 1;

    return (*G159_rank) + j + jm + 5;
}

void
wrtdlst (text * fyle, int *qq, int col, int rank, termptr list)
{
    bool startx;
    int ds;
    register int i;
    int *F160;
    termptr W6;

    F160 = G159_rank;
    G159_rank = &rank;
    if ((*qq) == 1) {
        (void) fprintf ((*fyle).fp, "      "), Putl ((*fyle), 0);
        (*qq) = 7;
    }
    if (list == NULL) {
        (void) fprintf ((*fyle).fp, "zero"), Putl ((*fyle), 0);
        (*qq) = (*qq) + 4;
    }
    startx = true;
    while (list != NULL) {
        W6 = &(*list);
        ds = dsize (list);
        if ((*qq) + ds < col) {
            if (startx) {
                if (W6->mult < 0)
                    (void) fprintf ((*fyle).fp, " -"), Putl ((*fyle), 0);
                else
                    ds = ds - 3;
                startx = false;
            } else if (W6->mult < 0)
                (void) fprintf ((*fyle).fp, " -"), Putl ((*fyle), 0);
            else
                (void) fprintf ((*fyle).fp, " +"), Putl ((*fyle), 0);
            if (abs (W6->mult) != 1)
                (void) fprintf ((*fyle).fp, "%1d", abs (W6->mult)),
                    Putl ((*fyle), 0);
            Putchr ('(', (*fyle));
            for (i = 1; i <= (*G159_rank); i++) {
                (void) fprintf ((*fyle).fp, "%1d", W6->val.A[i]),
                    Putl ((*fyle), 0);
            }
            Putchr (')', (*fyle));
            (*qq) = (*qq) + ds;
            list = W6->next;
        } else {
            Putchr ('\n', (*fyle));
            (void) fprintf ((*fyle).fp, "      "), Putl ((*fyle), 0);
            (*qq) = 7;
        }
    }
    G159_rank = F160;
}

/** Returns a big number as an ordinary int */
int
frbig (bframe big)
{
    register int num;
    register int i;

    num = 0;

    for (i = maxb; i >= 1; i--) {
        num = base * num + big.A[i];
    }
    return num;
}


void
bigadd (bframe p, bframe q, bframe * pq)
{
    register int j, maxw;
    register int i;

    (*pq) = nulls;
    maxw = maxl;
    while ((p.A[maxw] == 0) && (q.A[maxw] == 0) && (maxw > 0))
        maxw = maxw - 1;
    i = 1;
    do {
        j = p.A[i] + q.A[i] + pq->A[i];
        pq->A[i] = j % base;
        i = i + 1;
        pq->A[i] = j / base;
    }
    while (!(i > maxw));
    maxw = maxl;
    while ((pq->A[maxw] == 0) && (maxw > 0))
        maxw = maxw - 1;
    if (pq->A[maxw] < 0)
        for (i = 1; i <= maxw; i++) {
            if (pq->A[i] > 0) {
                pq->A[i] = pq->A[i] - base;
                pq->A[i + 1] = pq->A[i + 1] + 1;
            }
    } else
        for (i = 1; i <= maxw; i++) {
            if (pq->A[i] < 0) {
                pq->A[i] = pq->A[i] + base;
                pq->A[i + 1] = pq->A[i + 1] - 1;
            }
        }
}

/** Subtracts two big numbers to give a signed big number */
void
bigsubtr (bframe pq, bframe p, bframe * q, bool * neg)
{
    register int maxw;
    register int i;
    for (i = 0; i <= maxb; i++) {
        p.A[i] = -p.A[i];
    }
    bigadd (pq, p, &(*q));
    maxw = maxb;
    while ((q->A[maxw] == 0) && (maxw > 0))
        maxw = maxw - 1;
    if ((maxw > 0) && (q->A[maxw] < 0))
        (*neg) = true;
    else
        (*neg) = false;
}


void
cadd (bframe p, int q, bframe * pq)
{
    register int i, j, k, maxw;

    (*pq) = nulls;
    i = 1;
    k = q % base;
    maxw = maxl;
    while ((p.A[maxw] == 0) && (maxw > 0))
        maxw = maxw - 1;
    while ((q != 0) || (i <= maxw)) {
        j = p.A[i] + k + pq->A[i];
        if (j < 0) {
            j = j + base;
            pq->A[i + 1] = -1;
        }
        pq->A[i] = j % base;
        pq->A[i + 1] = pq->A[i + 1] + j / base;
        q = q / base;
        k = q % base;
        i = i + 1;
    }
}

/** Multiplies the big number m by the int n to produce a big number mn */
void
cmult (bframe m, int n, bframe * mn)
{
    register int maxm, k;
    register int j;

    maxm = maxb;
    while ((m.A[maxm] == 0) && (maxm > 0))
        maxm = maxm - 1;
    for (j = 0; j <= maxb; j++) {
        mn->A[j] = 0;
    }
    for (j = 1; j <= maxm; j++) {
        k = mn->A[j] + m.A[j] * n;
        mn->A[j] = k % base;
        mn->A[j + 1] = k / base;
    }
}


void
factor (int num, bframe * fnum)
{
    register int i;
    bool finished;

    finished = false;
    (*fnum) = nulls;
    i = 1;
    do {
        if (abs (num) > 1)
            if (num % primes.A[i] == 0) {
                fnum->A[i] = fnum->A[i] + 1;
                num = num / primes.A[i];
            } else {
                i = i + 1;
                if (i > maxl)
                    finished = true;
        } else
            finished = true;
    }
    while (!(finished));
    fnum->A[0] = num;
}


/** Multiplies two big numbers n, m to produce the big number nm */
void
fmult (bframe n, bframe m, bframe * nm)
{
    register int i;
    bframe nn, mn;

    factor (n.A[0], &nn);
    factor (m.A[0], &mn);
    for (i = 1; i <= maxb; i++) {
        nm->A[i] = n.A[i] + m.A[i] + nn.A[i] + mn.A[i];
    }
    nm->A[0] = nn.A[0] * mn.A[0];
}


/** convert a bignum represented as prime factors into a bignum in binary coded decimal*/
void
factobig (bframe fnum, bframe * bignum)
{
    register int i;

    tobig (fnum.A[0], &(*bignum));
    for (i = 1; i <= maxb; i++) {
        while (fnum.A[i] > 0) {
            cmult ((*bignum), primes.A[i], &(*bignum));
            fnum.A[i] = fnum.A[i] - 1;
        }
    }
}

void
cdiv (bframe dividend, int divisor, bframe * quotient, int *remainder)
{
    register int i, p;

    p = 5;
    if (divisor == 0) {
        error (ZERO_DIVISOR, p);
        error (BAD_IRREP, p);
    } else {
        (*quotient) = nulls;
        i = maxl;
        (*remainder) = 0;
        while ((dividend.A[i] == 0) && (i > 0))
            i = i - 1;
        while (i != 0) {
            quotient->A[i] = ((*remainder) * base + dividend.A[i]) / divisor;
            (*remainder) = ((*remainder) * base + dividend.A[i]) % divisor;
            i = i - 1;
        }
    }
}

void
bigtofact (bframe bigno, bframe * fno)
{
    register int i;
    int remain = 0;
    bframe quot;
    bool nought;

    remain = 0;
    (*fno) = nulls;
    i = 1;
    nought = true;
    while ((i < maxl) && nought) {
        if (bigno.A[i] != 0)
            nought = false;
        i = i + 1;
    }
    i = 1;
    if (!nought)
        do {
            cdiv (bigno, primes.A[i], &quot, &remain);
            if (remain == 0) {
                fno->A[i] = fno->A[i] + 1;
                bigno = quot;
            } else
                i = i + 1;
        }
        while (!(i > maxl));
    fno->A[0] = frbig (bigno);
}

/** divide n by m , quotient in quot, n,m and quot are represented by prime factors.*/
void
fdiv (bframe n, bframe m, bframe * quot)
{
    register int i, p;

    p = 5;
    if (m.A[0] == 0) {
        error (ZERO_DIVISOR, p);
        error (BAD_IRREP, p);
    } else {
        for (i = 1; i <= maxb; i++) {
            quot->A[i] = n.A[i] - m.A[i];
        }
        quot->A[0] = n.A[0] / m.A[0];
    }
}

void
kmult (bframe m, int n, bframe * mn)
{
    bframe fnk;

    factor (n, &fnk);
    fmult (m, fnk, &(*mn));
}

void
kdiv (bframe mn, int m, bframe * n)
{
    bframe fm;

    factor (m, &fm);
    fdiv (mn, fm, &(*n));
}


void
maxscoeff (termptr charlist)
{
    int maxc;

    maxc = 0;
    while (charlist != NULL) {
        if (charlist->mult > maxc)
            maxc = charlist->mult;
        charlist = charlist->next;
    }
    inform ("MaxCoeff = ;", cont);
    (void) fprintf (output.fp, "%10d\n", maxc), Putl (output, 1);
}

void
maxrcoeff (ocharptr charlist)
{
    int maxc;
    ocharptr temp;
    maxc = 0;
    while (charlist != NULL) {
        temp = charlist;
        if (temp->mult > maxc)
            maxc = temp->mult;
        charlist = temp->next;
    }
    inform ("MaxCoeff = ;", cont);
    (void) fprintf (output.fp, "%10d\n", maxc), Putl (output, 1);
}

void
maxdcoeff (prodtype charlist)
{
    int maxc;
    prodtype temp;
    maxc = 0;
    while (charlist != NULL) {
        temp = charlist;

        if (temp->mult > maxc)
            maxc = temp->mult;
        charlist = temp->next;
    }
    inform ("MaxCoeff = ;", cont);
    fprintf (output.fp, "%10d\n", maxc), Putl (output, 1);
}


void
ssummult (termptr list)
{
    int t, qq;
    bframe resultx, temp;
    termptr W39;

    qq = 12;
    tobig (0, &resultx);
    while (list != NULL) {
        W39 = &(*list);

        t = /*abs */ (W39->mult);       /*10/6/99 */
        factor (t, &temp);
        factobig (temp, &temp);
        bigadd (resultx, temp, &resultx);
        list = W39->next;
    }
    (void) fprintf (output.fp, "CoeffSum = "), Putl (output, 1);
    wrtbigno (&output, &qq, &resultx);
    Putchr ('\n', output);
    if (logging) {
        (void) fprintf (logfile.fp, "CoeffSum = "), Putl (logfile, 1);
        wrtbigno (&logfile, &qq, &resultx);
        Putchr ('\n', logfile);
    }
}

void
summult (ocharptr list)
{
    int t, qq;
    bframe resultx, temp;
    ocharptr W35;

    qq = 12;
    tobig (0, &resultx);
    while (list != NULL) {
        W35 = &(*list);

        t = /*abs */ (W35->mult);       /*10/6/99 */
        factor (t, &temp);
        factobig (temp, &temp);
        bigadd (resultx, temp, &resultx);
        list = W35->next;
    }
    (void) fprintf (output.fp, "CoeffSum = "), Putl (output, 1);
    wrtbigno (&output, &qq, &resultx);
    Putchr ('\n', output);
    if (logging) {
        (void) fprintf (logfile.fp, "CoeffSum = "), Putl (logfile, 1);
        wrtbigno (&logfile, &qq, &resultx);
        Putchr ('\n', logfile);
    }

}


void
psummult (prodtype list)
{
    int t, qq;
    bframe resultx, temp;

    qq = 12;
    tobig (0, &resultx);
    while (list != NULL) {
        register prodtype W40 = &(*list);

        t = /*abs */ (W40->mult);       /*10/6/99 */
        factor (t, &temp);
        factobig (temp, &temp);
        bigadd (resultx, temp, &resultx);
        list = W40->next;
    }
    (void) fprintf (output.fp, "CoeffSum = "), Putl (output, 1);
    wrtbigno (&output, &qq, &resultx);
    Putchr ('\n', output);
    if (logging) {
        (void) fprintf (logfile.fp, "CoeffSum = "), Putl (logfile, 1);
        wrtbigno (&logfile, &qq, &resultx);
        Putchr ('\n', logfile);
    }
}

/** for COUNTTermsInList in repmode */
void
tsum (ocharptr charlist)
{
    int t, qq;
    bframe summ, temp;

    qq = 12;
    tobig (0, &summ);
    while (charlist != NULL) {
        t = 1;
        factor (t, &temp);
        factobig (temp, &temp);
        bigadd (summ, temp, &summ);
        charlist = charlist->next;
    }
    fprintf (output.fp, "TermCount = ");
    wrtbigno (&output, &qq, &summ);
    Putchr ('\n', output);
    if (logging) {
        fprintf (logfile.fp, "TermCount = ");
        wrtbigno (&logfile, &qq, &summ);
        Putchr ('\n', logfile);
    }
}

/** for COUNTTermsInList in sfn mode */
void
tssum (termptr charlist)
{
    int t, qq;
    bframe summ, temp;

    qq = 12;
    tobig (0, &summ);
    while (charlist != NULL) {
        t = 1;
        factor (t, &temp);
        factobig (temp, &temp);
        bigadd (summ, temp, &summ);
        charlist = charlist->next;
    }
    fprintf (output.fp, "TermCount = ");
    wrtbigno (&output, &qq, &summ);
    Putchr ('\n', output);
    if (logging) {
        fprintf (logfile.fp, "TermCount = ");
        wrtbigno (&logfile, &qq, &summ);
        Putchr ('\n', logfile);
    }
}

/** for COUNTTermsInList in DPmode */
void
tpsum (prodtype charlist)
{
    int t, qq;
    bframe summ, temp;

    qq = 12;
    tobig (0, &summ);
    while (charlist != NULL) {
        t = 1;
        factor (t, &temp);
        factobig (temp, &temp);
        bigadd (summ, temp, &summ);
        charlist = charlist->next;
    }
    print ("TermCount = ");
    wrtbigno (&output, &qq, &summ);
    Putchr ('\n', output);
    if (logging) {
        wrtbigno (&logfile, &qq, &summ);
        Putchr ('\n', logfile);
    }
}

void
hklth (frame irrep, bframe * hl)
{
    register int i, j, r;
    bframe l, fl;

    l = nulls;
    r = 0;
    while (irrep.A[r + 1] != 0)
        r = r + 1;
    i = 1;
    while (irrep.A[i] != 0) {
        l.A[i] = irrep.A[i] + r - i;
        i = i + 1;
    }
    factor (1, &(*hl));
    for (i = 1; i <= r; i++) {
        factor (1, &fl);
        for (j = 1; j <= l.A[i]; j++) {
            kmult (fl, j, &fl);
        }
        fmult (fl, (*hl), &(*hl));
    }
    for (i = 2; i <= r; i++) {
        for (j = 1; j <= i - 1; j++) {
            kdiv ((*hl), l.A[j] - l.A[i], &(*hl));
        }
    }
}

void
tobig (int num, bframe * bigno)
{
    register int i;

    (*bigno) = nulls;
    i = 1;
    do {
        bigno->A[i] = num % base;
        num = num / base;
        i = i + 1;
    }
    while (!(num == 0));
}

void
multlist (termptr list)
{
    termptr ptr;
    register int i;
    unsigned sum[2 * maxdim + 1];       // count from -maxdim to maxdim

    for (i = -maxdim; i <= maxdim; i++)
        sum[i + maxdim] = 0;

    for (ptr = list; ptr != NULL; ptr = ptr->next)
        if (ptr->mult < -maxdim || ptr->mult > maxdim) {
            fprintf (stderr, "multiplicities lower or greater than -/+ %d\n",
                     maxdim);
            error (MISTAKE, 0);
        } else
            sum[ptr->mult + maxdim]++;

    for (i = -maxdim; i <= maxdim; i++)
        if (sum[i + maxdim])
            print ("Mult =%10d Number of Terms =%10d\n", i, sum[i + maxdim]);
}

void
factorialn (int n, bframe * resultx)
{
    register int i;
    tobig (n, &(*resultx));
    if ((n < 0)) {
        fprintf (output.fp, "ERROR: negative factorial\n"), Putl (output, 1);
        tobig (0, &(*resultx));
    } else if ((n >= 2))
        for (i = (n - 1); i >= 1; i--)
            cmult (*resultx, i, &(*resultx));
}

void
sumabssquares (termptr list)
{
    int t, qq;
    bframe resultx, temp;
    termptr W57;

    qq = 12;
    tobig (0, &resultx);
    while (list != NULL) {
        W57 = &(*list);

        t = abs (W57->mult);
        factor (t, &temp);
        fmult (temp, temp, &temp);
        factobig (temp, &temp);
        bigadd (resultx, temp, &resultx);
        list = W57->next;
    }
    print ("Sum of absolute squares of multiplicities =");
    wrtbigno (&output, &qq, &resultx);
    Putchr ('\n', output);
    if (logging) {
        wrtbigno (&logfile, &qq, &resultx);
        Putchr ('\n', logfile);
    }
}


void
sumabssquaresrep (ocharptr list)
{
    int t, qq;
    bframe resultx, temp;
    ocharptr W57;

    qq = 12;
    tobig (0, &resultx);
    while (list != NULL) {
        W57 = &(*list);
        t = abs (W57->mult);
        factor (t, &temp);
        fmult (temp, temp, &temp);
        factobig (temp, &temp);
        bigadd (resultx, temp, &resultx);
        list = W57->next;
    }
    print ("Sum of absolute squares of multiplicities =");
    wrtbigno (&output, &qq, &resultx);
    Putchr ('\n', output);
    if (logging) {
        wrtbigno (&logfile, &qq, &resultx);
        Putchr ('\n', logfile);
    }
}
