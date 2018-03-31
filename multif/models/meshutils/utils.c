#include "meshutils.h"

//--- Transforms all letters to lower case
int Str2Lower(char *buff)
{
  int iChr;
  
  for (iChr=0; iChr<strlen(buff); iChr++) 
    buff[iChr] = tolower( buff[iChr] );

  return 1;
}

//--- Removes all occurences of char c from str
void StrRemoveChars (char* s, char ch) {
	char *p = s;
	while (*s) {
	    if (*s != ch)
	        *p++ = *s;
	    s++;
	}
	*p = 0;
}
