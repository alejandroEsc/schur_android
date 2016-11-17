typedef struct representation {
  struct representation *prev, *next;
  byte l;
  byte n;
  int mult;
} representation;
typedef representation *rep_ptr;

ocharptr uadd(ocharptr, int, int);

ocharptr snanbrnch(ocharptr, int);

ocharptr onsnbr(ocharptr, int, bool);

ocharptr onsnbr1(ocharptr, int, bool);

ocharptr so2k1so3brnch(char, ocharptr, int, bool);
ocharptr g2_so7(ocharptr);
ocharptr g2_so7a (ocharptr lambda);
ocharptr g2so3brnch(ocharptr);
ocharptr so7g2brnch(ocharptr);

void u1dtrmne(ocharptr *,int ,int ,int ,int ,int);

ocharptr unsubgrbr(ocharptr,int ,int ,char);

ocharptr unospnbrnch(char ,ocharptr ,int  );

ocharptr supqu1br(ocharptr ,int ,int ,int ,int ,int);

ocharptr onmonombrnch(ocharptr ,int ,int ,int);

ocharptr sonunbrnch(ocharptr, int);

ocharptr unsnbr(ocharptr, int);

ocharptr unsnbr1(ocharptr, int);

bool get_ln(void);

void so7_to_lambda(void);

void init(void);

void dispose_heap(void);

void first_term(int, int, int);
void second_term(int, int, int);
void third_term(int, int, int);
void fourth_term(int, int, int);
void fifth_term(int, int);

void P62_branch(int,int,int);

ocharptr sonunbr1(ocharptr ,int);
ocharptr sonunbr2(ocharptr, int);

void insert_node (rep_ptr *, rep_ptr *, byte, byte, int);
void delete_node (rep_ptr *, rep_ptr *);

bool found (rep_ptr * head, rep_ptr * currentb, int l, int n);

