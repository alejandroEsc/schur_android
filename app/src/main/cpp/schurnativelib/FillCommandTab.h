#ifndef FILLCOMMAND_H
#define FILLCOMMAND_H  1

#include "utils.h"

#define MAXNBCOMMANDS   300

/* have to be power of 2. BRMODE is a part of DPMODE */
#define DPMODE   1
#define REPMODE  2
#define SFNMODE  4


//Global variable for command completion
/*typedef struct command 
{
	char *shortname;
	char *name;
	char modes;     // a inclusive OR of modes SFNMODE, DPMODE or REPMODE.
} command; 

command CommandTab[MAXNBCOMMANDS];*/
char *CommandTab[MAXNBCOMMANDS];

void initialise_readline (void);
int fillCommandTab(void);

char *command_generator (const char *, int );

char **Command_completion (const char *, int, int);

#endif  /* ! FILLCOMMAND_H */
