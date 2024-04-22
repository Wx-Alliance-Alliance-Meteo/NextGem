#include <stdio.h>
#include <stdlib.h>

int c_ut_report(int status, char *message) {
  FILE *fp;
  fp = fopen("test_report.txt","w+");
  if(! fp){
    printf("ERROR in c_ut_report: cannot open file test_report.txt\n");
    return(-1);
  }
  if(status == 0){
    fprintf(fp," ok\n");
  } else {
    fprintf(fp," failed %s\n",message);
  }
  return(status);
}
