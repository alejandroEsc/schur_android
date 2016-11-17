/** \file mymalloc.
 * allocate memory and initialize it with spaces 
 */

#include "standard.h"

//for aaargh
#include "utils.h"
#include "mymalloc.h"

/** Allocate l bytes of memory and initialize it with spaces. */
char *mymalloc(unsigned l) {
	char *ptr;

	ptr=(char *)malloc(l);
	if (ptr==NULL) {
	  aaargh("memory allocation failed",true);
	  exit(1);
	}
	else {
	  memset(ptr,' ',l);
	  ptr[l-1]='\0';
	  return(ptr);
	}
}
