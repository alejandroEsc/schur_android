#ifndef SKEW_H
#define SKEW_H 1
bool skewcompatable(frame, frame);

termptr allskews(frame, bool);

termptr skew(frame, frame);
termptr skewx(int, frame ,frame , int);
termptr lskew(termptr, termptr);
termptr pskew(termptr, termptr);
termptr nskew(termptr,termptr, int);

#endif /* !SKEW_H */
