/* Do not try to compile it directly */
#ifdef MYSTATIC
MYSTATIC int
qcf (termptr point)
{
    register int i;
    if (point == NULL) {
        sign = LESS;
        return LESS;
    } else {
        i = 1;
        while ((MYENV(rho).A[i] == point->val.A[i]) && (MYENV(rho).A[i] != 0))
            i = i + 1;
        return sgn (point->val.A[i] - MYENV(rho).A[i]);
    }
}

MYSTATIC void
segsort (void)
{
    termptr listptr1, listptr2 = NULL;
    register int segment;
    progress ();
    if (MYENV(firsttime)) {
        snu (&listptr1);
        MYENV(firsttime) = false;
        {
            listptr1->val = MYENV(rho);
            listptr1->mult = MYENV(coef);
            listptr1->next = NULL;
        }
        for (segment = 1; segment <= SEGPTRMAX; segment++)
            MYENV(segptrs).A[segment - 1] = listptr1;
    } else {
        segment = 1;
        while (qcf (MYENV(segptrs).A[segment - 1]) == -1)
            segment = segment + 1;
        listptr1 = MYENV(segptrs).A[segment - 1];
        if (sign == GREATER)
            do {
                listptr2 = listptr1;
                listptr1 = listptr1->next;
            }
            while (!(qcf (listptr1) != 1));
        if (sign == EQUAL)
            listptr1->mult = listptr1->mult + MYENV(coef);
        else {
            snu (&listptr1);
            listptr1->next = listptr2->next;
            listptr2->next = listptr1;
            listptr1->mult = MYENV(coef);
            listptr1->val = MYENV(rho);
        }
        if (segment == SEGPTRMAX)
            segment = SEGPTRMAX - 1;
        MYENV(segptrs).A[segment - 1] = listptr1;
    }
}                               // void segsort(void)

MYSTATIC void
h_count (int tabrow)
{
    int beta, rowend, rowbegin, dummy;
    register int j;
    register int i;
    bool valid, gox;
    do {
        valid = true;
        beta = MYENV(corelen) + tabrow;
        rowbegin = MYENV(offset).A[tabrow] + 1;
        rowend = MYENV(frme).A[tabrow];
        if (rowend >= rowbegin) {
            if (MYENV(tableau).A[tabrow].A[rowbegin - 1] == 0) {
                if ((tabrow) & 1)
                    progress ();
                for (i = rowbegin; i <= rowend; i++)
                    MYENV(tableau).A[tabrow].A[i - 1] =
                        MYENV(tableau).A[tabrow - 1].A[i - 1] + 1;
            } else {
                i = rowbegin;
                do {
                    gox = (i <= rowend);
                    if (gox)
                        gox = (MYENV(tableau).A[tabrow].A[i - 1] != beta);
                    if (gox) {
                        MYENV(rho).A[MYENV(tableau).A[tabrow].A[i - 1]] -= 1;
                        i = i + 1;
                    }
                }
                while (gox);
                i = i - 1;
                dummy = MYENV(tableau).A[tabrow].A[i - 1] + 1;
                MYENV(rho).A[beta] += i - rowend;
                MYENV(tableau).A[tabrow].A[i - 1] = dummy;
                for (j = i + 1; j <= rowend; j++) {
                    if (dummy > MYENV(tableau).A[tabrow - 1].A[j - 1] + 1)
                        MYENV(tableau).A[tabrow].A[j - 1] = dummy;
                    else
                        MYENV(tableau).A[tabrow].A[j - 1] =
                            MYENV(tableau).A[tabrow - 1].A[j - 1] + 1;
                }
            }
            for (j = rowend; j >= rowbegin; j--) {
                dummy = MYENV(tableau).A[tabrow].A[j - 1];
                MYENV(rho).A[dummy] += 1;
                if (dummy != 1)
                    valid = valid && (MYENV(rho).A[dummy] <= MYENV(rho).A[dummy - 1]);
            }
        }
        if (valid && (MYENV(rho).A[MYENV(limitx) + 1] == 0))
            if (tabrow == MYENV(framelen))
                segsort ();
            else
                h_count (tabrow + 1);
    }
    while (!((MYENV(tableau).A[tabrow].A[rowbegin - 1] == beta)
             || (rowend < rowbegin)));
    MYENV(rho).A[beta] += rowbegin - rowend - 1;
    MYENV(tableau).A[tabrow].A[rowbegin - 1] = 0;
}                               // h_count
#endif
