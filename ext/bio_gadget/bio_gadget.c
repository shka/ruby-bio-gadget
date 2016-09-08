#include <limits.h>
#include <regex.h>
#include <stdio.h>
#include <string.h>
#include "bio_gadget.h"

VALUE bio_gadget_at3(vSelf, vLen, vCmdIn, vPathOut)
     VALUE vSelf;
     VALUE vLen;
     VALUE vCmdIn;
     VALUE vPathOut;
{
  char line[BUFSIZE];
  char *acc;
  char *seq;
  char *sep;
  char *qual;
  int len;
  unsigned long seqlen;
  FILE *fp_in;
  FILE *fp_out;

  fp_in = popen(StringValueCStr(vCmdIn), "r");
  fp_out = fopen(StringValueCStr(vPathOut), "w");
  len = NUM2INT(vLen);
  while(fgets(line, BUFSIZE, fp_in) != NULL) {
    acc = strtok(line, "\t");
    seq = strtok(NULL, "\t");
    sep = strtok(NULL, "\t");
    qual = strtok(NULL, "\t");
    seqlen = strlen(seq)-len;
    if (seqlen > 0) {
      seq[seqlen] = 0;
      qual[seqlen] = 0; // qual[strlen(qual)-len] = 0;
      fprintf(fp_out, "%s\t%s\t%s\t%s\n", acc, seq, sep, qual);
    }
  }
  fclose(fp_out);
  fclose(fp_in);
  
  return Qnil;
}

VALUE bio_gadget_nr_deg(vSelf, vCmd)
     VALUE vSelf;
     VALUE vCmd;
{
  char line[BUFSIZE];
  char *acc;
  char *seq;
  char *sep;
  char *qual;
  char pseq[BUFSIZE] = "";
  unsigned long pseql = ULONG_MAX;
  FILE *fp;
  char regexs[BUFSIZE];
  regex_t regexc;

  fp = popen(StringValueCStr(vCmd), "r");
  while(fgets(line, BUFSIZE, fp) != NULL) {
    acc = strtok(line, "\t");
    seq = strtok(NULL, "\t");
    sep = strtok(NULL, "\t");
    qual = strtok(NULL, "\t");
    sprintf(regexs, "^%s", seq);
    regcomp(&regexc, regexs, REG_NOSUB);
    if (strlen(seq) > pseql
	|| regexec(&regexc, pseq, 0, NULL, 0) == REG_NOMATCH) {
      printf("%s\t%s\t%s\t%s", acc, seq, sep, qual);
      strcpy(pseq, seq);
      pseql = strlen(seq);
    }
    regfree(&regexc);
  }
  fclose(fp);
  
  return Qnil;
}

VALUE bio_gadget_nr_std(vSelf, vCmd)
     VALUE vSelf;
     VALUE vCmd;
{
  char line[BUFSIZE];
  char *acc;
  char *seq;
  char *sep;
  char *qual;
  char pseq[BUFSIZE] = "";
  FILE *fp;

  fp = popen(StringValueCStr(vCmd), "r");
  while(fgets(line, BUFSIZE, fp) != NULL) {
    acc = strtok(line, "\t");
    seq = strtok(NULL, "\t");
    sep = strtok(NULL, "\t");
    qual = strtok(NULL, "\t");
    if (strcmp(pseq, seq) != 0) {
      printf("%s\t%s\t%s\t%s", acc, seq, sep, qual);
      strcpy(pseq, seq);
    }
  }
  fclose(fp);
  
  return Qnil;
}

VALUE bio_gadget_qt3(vSelf, vLQs)
     VALUE vSelf;
     VALUE vLQs;
{
  char line[BUFSIZE];
  char *acc;
  char *seq;
  char *sep;
  char *qual;
  char *lqs;
  unsigned long seqlen;

  lqs = StringValueCStr(vLQs);
  
  while(fgets(line, BUFSIZE, stdin) != NULL) {
    acc = strtok(line, "\t");
    seq = strtok(NULL, "\t");
    sep = strtok(NULL, "\t");
    qual = strtok(NULL, "\t\n");
    seqlen = strcspn(qual, lqs);
    if (seqlen > 0) {
      if (seqlen < strlen(qual)) {
	seq[seqlen] = 0;
	qual[seqlen] = 0;
      }
      printf("%s\t%s\t%s\t%s\n", acc, seq, sep, qual);
    }
  }
  
  return Qnil;
}

VALUE rb_mBio_Gadget;

void
Init_bio_gadget(void)
{
  rb_mBio_Gadget = rb_define_module("BioGadget");
  rb_define_module_function(rb_mBio_Gadget, "at3", bio_gadget_at3, 3);
  rb_define_module_function(rb_mBio_Gadget, "nr_deg", bio_gadget_nr_deg, 1);
  rb_define_module_function(rb_mBio_Gadget, "nr_std", bio_gadget_nr_std, 1);
  rb_define_module_function(rb_mBio_Gadget, "qt3", bio_gadget_qt3, 1);
}
