termptr spseries(int, int);

termptr rseries(int, char);

termptr makelist(int);

termptr allseriesbut_t(char, int, frame, unsigned);

termptr onlyseries_t(int, frame);

termptr seriesx(char, int, frame);

bool validser(char);

termptr useseries(char, termptr, bool, bool, int);

void cleave(termptr *);

termptr seriesm(int, bool);

termptr rhook(int);
