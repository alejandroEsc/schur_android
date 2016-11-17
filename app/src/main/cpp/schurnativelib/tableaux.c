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

/** \file tableaux.c
 * functions used for young tableaux computation and display
 */

#include <stdio.h>
#include "standard.h"
#include "define.h"
#include "dim.h"
#include "type.h"
#include "var.h"
#include "utils.h"
#include "s1.h"
#include "s2.h"
#include "s5.h"
#include "s6.h"
#include "g.h"
#include "tableaux.h"

/** returns the number of loops to do */
int
forloops (int x, int y)
{
    if (x < y)
        return 0;
    else
        return (x - y + 1);
}

/** write to file f iu spaces */
void
spaces (FILE * f, int iu)
{
    register int ju;
    for (ju = 1; ju <= iu; ju++)
        fputc (' ', f);
}

void
repeatstr (FILE * f, char *str, int repeat)
{
    register int i;
    for (i = 1; i <= repeat; i++)
        fputs (str, f);
}

void
dtablx (termptr list, bool rows, FILE * f, bool vdu)
{
    char conty;
    int longest, leng, posn, num, feld, lines;
    register int iu;
    register int n;
    register int ju;
    bool okd;
    termptr startd, ptr;
    frame rho;
    conty = ' ';
    startd = list;
    ptr = list;
    lines = 0;
    while ((startd != NULL) && (conty != esc)) {
        longest = 0;
        posn = 1;
        num = 0;
        okd = true;
        while (okd) {
            rho = ptr->val;
            leng = len (&rho);
            if (leng > longest)
                longest = leng;
            posn = posn + rho.A[1] + 3;
            if (posn < 79)
                num = num + 1;
            else
                okd = false;
            if (okd) {
                if (ptr->mult != 1)
                    fprintf (f, "%3d", ptr->mult);
                else
                    spaces (f, 3);
                spaces (f, rho.A[1]);
            }
            ptr = ptr->next;
            okd = (bool) (okd && (ptr != NULL));
        }
        fputc ('\n', f);
        lines = lines + longest + 3;
        if ((lines > tlines - 3) && more && vdu) {
            lines = 0;
            fprintf (f, "more");
            conty = Getchr (input);
            fputc ('\n', f);
            if (conty == esc)
                fprintf (f, "escape. ");
        }
        if (conty != esc) {
            for (ju = 1; ju <= longest + 1; ju++) {
                ptr = startd;
                for (n = 1; n <= num; n++) {
                    rho = ptr->val;
                    feld = rho.A[1] + 3;
                    if (ju <= (len (&rho) + 1)) {
                        if (ju == 1)
                            if (rho.A[1] != 0) {
                                if (rows)
                                    fprintf (f, "%2d", rho.A[1]);
                                else
                                    spaces (f, 2);
                                for (iu = 1; iu <= rho.A[1]; iu++) {
                                    fputc ('O', f);
                                }
                                fputc (' ', f);
                            } else {
                                if (rows)
                                    fprintf (f, "%2d", 0);
                                else
                                    spaces (f, 2);
                                fputc ('.', f);
                        } else {
                            if (rho.A[ju] == 0) {
                                spaces (f, 2 + rho.A[ju - 1] + 1);
                                spaces (f,
                                        feld - forloops (rho.A[ju - 1],
                                                         2) - 4);
                            } else {
                                if (rows)
                                    fprintf (f, "%2d", rho.A[ju]);
                                else
                                    spaces (f, 2);
                                for (iu = 1; iu <= rho.A[ju]; iu++) {
                                    fputc ('O', f);
                                }
                                spaces (f, rho.A[1] - rho.A[ju] + 1);
                            }
                        }
                    } else
                        spaces (f, feld);
                    ptr = ptr->next;
                }
                fputc ('\n', f);
            }
        }
        startd = ptr;
    }
    fputc ('\n', f);
}

int
hookvalue (bool conhk, int numb, int i, int col, frame f)
{
    if (conhk)
        return (hookpart (i, col, f));
    else
        return (col - i + numb);
}

/** output to file f a line of frame/tableau for htablx */
int
putlineTab (FILE * f, bool conhk, int numb, int line, unsigned short cellSize,
            unsigned short width, frame t)
{
    char format[MAXSTRING];
    int col;
    int sum = 0, hook;

    sprintf (format, "%%%dd|", cellSize);
    if (t.A[line] != 0) {
        fprintf (f, "|");
        for (col = 1; col <= t.A[line]; col++) {
            hook = hookvalue (conhk, numb, line, col, t);
            sum += hook;
            fprintf (f, format, hook);
        }
        spaces (f, (width - t.A[line]) * (cellSize + 1));       // complete with spaces
    } else
        spaces (f, (cellSize + 1) * width + 1);
    return sum;
}

void
lineBox (FILE * f, int nb, unsigned short cellSize, unsigned short width)
{
    int i;
    char str[10];               // it must be enough!

    if (nb > 0) {
        str[0] = '+';
        for (i = 1; i <= cellSize; i++)
            str[i] = '-';
        str[i] = '\0';
        repeatstr (f, str, nb);
        fprintf (f, "+");
    } else
        spaces (f, 1);
    spaces (f, (width - nb) * (cellSize + 1));
}

void
putSum (FILE * f, int sum, unsigned short cellSize, unsigned short width)
{
    char str[MAXSTRING];

    sprintf (str, "%d", sum);
    fprintf (f, "%s", str);
    spaces (f, (cellSize + 1) * width + 1 - strlen (str));
}

void
htablx (termptr list, FILE * f, int numb, bool conhk, bool vdu,
        bool displaySum)
{
    termptr ptr, ptr1;
    int minval, maxval, i, currentwidth, maxheight = 0, totalwidth =
        0, hook;
    unsigned line, nbpart, maxnbline;
    char str[MAXSTRING];
    struct Sizestableau_t {
        struct Sizestableau_t *next;
        unsigned short width;
        unsigned short height;
        unsigned short cellSize;
        int sum;
    } *tableaux, *ptrTableau, *ptrTemp, *ptrTab;

    // compute width and heigth for each tableau
    ptr = list;
    tableaux = NULL;
    ptrTableau = NULL;
    nbpart = 0;
    while (ptr != NULL) {
        nbpart++;
        ptrTemp =
            (struct Sizestableau_t *) malloc (sizeof (struct Sizestableau_t));
        ptrTemp->next = NULL;
        ptrTemp->width = 0;
        ptrTemp->height = 0;
        ptrTemp->sum = 0;
        if (tableaux == NULL) {
            tableaux = ptrTemp;
            ptrTableau = ptrTemp;
        } else {
            ptrTableau->next = ptrTemp;
            ptrTableau = ptrTemp;
        }

        minval = 0;
        maxval = 0;
        for (i = 1; ptr->val.A[i] != 0; i++) {
            if (ptr->val.A[i] > ptrTableau->width)
                ptrTableau->width = ptr->val.A[i];

            hook = hookvalue (conhk, numb, i, ptr->val.A[i], ptr->val);
            if (hook > maxval)
                maxval = hook;
        }
        if (!conhk && (numb - (i - 2) < minval))
            minval = numb - (i - 2);
        ptrTableau->height = i - 1;
        if (i - 1 > maxheight)
            maxheight = i - 1;
        totalwidth += ptrTableau->width;

        sprintf (str, "%d", minval);
        ptrTableau->cellSize = strlen (str);
        sprintf (str, "%d", maxval);
        if (strlen (str) > ptrTableau->cellSize)
            ptrTableau->cellSize = strlen (str);
        ptr = ptr->next;
    }

    ptr = list;
    ptrTableau = tableaux;
    while (ptr != NULL) {
        currentwidth = 0;
        ptr1 = ptr;             // continue on where we stop the first multiplicities.
        ptrTab = ptrTableau;
        maxnbline = 0;          // maximum number of lines to display for this serie of tableaux
        fprintf (f, "\n");
        while (ptr1 != NULL)    // first display multiplicities left justified above the tableau
        {
            currentwidth += (ptrTab->cellSize + 1) * ptrTab->width + 2;
            if (ptr1->mult != 1) {
                sprintf (str, "%d.", ptr1->mult);
                fprintf (f, "%s", str);
                spaces (f,
                        (ptrTab->cellSize + 1) * ptrTab->width + 1 -
                        strlen (str));
            } else
                spaces (f, (ptrTab->cellSize + 1) * ptrTab->width + 1);

            if (maxnbline < ptrTab->height)
                maxnbline = ptrTab->height;
            ptr1 = ptr1->next;
            ptrTab = ptrTab->next;
            if (ptr1 != NULL) {
                if (vdu
                    && currentwidth + (ptrTab->cellSize + 1) * ptrTab->width +
                    2 >= tcol)
                    break;      // too long : there will be more than one line of tableaux
                fprintf (f, " ");
            }
        }
        // draw initial top of the boxes +---+...
        fprintf (f, "\n");
        currentwidth = 0;
        ptr1 = ptr;
        ptrTab = ptrTableau;
        while (ptr1 != NULL) {
            currentwidth += ptrTab->width * (ptrTab->cellSize + 1) + 2;
            lineBox (f, ptr1->val.A[1], ptrTab->cellSize, ptrTab->width);
            ptr1 = ptr1->next;
            ptrTab = ptrTab->next;
            if (ptr1 != NULL) {
                if (vdu
                    && currentwidth + (ptrTab->cellSize + 1) * ptrTab->width +
                    2 >= tcol)
                    break;
                fprintf (f, " ");
            }
        }

        for (line = 1; line <= maxnbline; line++) {
            fprintf (f, "\n");
            currentwidth = 0;
            // now display tableaux while currentwidth<tcol if vdu is set
            ptr1 = ptr;
            ptrTab = ptrTableau;

            while (ptr1 != NULL) {
                currentwidth += ptrTab->width * (ptrTab->cellSize + 1) + 2;
                ptrTab->sum +=
                    putlineTab (f, conhk, numb, line, ptrTab->cellSize,
                                ptrTab->width, ptr1->val);
                ptr1 = ptr1->next;
                ptrTab = ptrTab->next;
                if (ptr1 != NULL) {
                    if (vdu
                        && currentwidth + ptrTab->width * (ptrTab->cellSize +
                                                           1) + 2 >= tcol)
                        break;
                    fprintf (f, " ");
                }
            }
            // draw bottom of the boxes +---+...
            fprintf (f, "\n");
            currentwidth = 0;
            ptr1 = ptr;
            ptrTab = ptrTableau;
            while (ptr1 != NULL) {
                currentwidth += ptrTab->width * (ptrTab->cellSize + 1) + 2;
                lineBox (f, ptr1->val.A[line], ptrTab->cellSize,
                         ptrTab->width);
                ptr1 = ptr1->next;
                ptrTab = ptrTab->next;
                if (ptr1 != NULL) {
                    if (vdu
                        && currentwidth + ptrTab->width * (ptrTab->cellSize +
                                                           1) + 2 >= tcol)
                        break;
                    fprintf (f, " ");
                }
            }
        }                       // for line
        if (displaySum) {
            // display sum under each tableau again limited by tcol
            fprintf (f, "\n");
            currentwidth = 0;
            ptr1 = ptr;
            ptrTab = ptrTableau;
            while (ptr1 != NULL) {
                currentwidth += ptrTab->width * (ptrTab->cellSize + 1) + 2;
                putSum (f, ptrTab->sum, ptrTab->cellSize, ptrTab->width);
                ptr1 = ptr1->next;
                ptrTab = ptrTab->next;
                if (ptr1 != NULL) {
                    if (vdu
                        && currentwidth + ptrTab->width * (ptrTab->cellSize +
                                                           1) + 2 >= tcol)
                        break;
                    fprintf (f, " ");
                }
            }
        }
        ptr = ptr1;             //continue if necessary 
        ptrTableau = ptrTab;
    }
    fprintf (f, "\n");

    // free memory
    ptrTableau = tableaux;
    while (ptrTableau != NULL) {
        ptrTemp = ptrTableau->next;
        free (ptrTableau);
        ptrTableau = ptrTemp;
    }
}                               // htablx
