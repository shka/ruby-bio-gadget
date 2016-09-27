#include <limits.h>
#include <regex.h>
#include <stdio.h>
#include <string.h>
#include "bio_gadget.h"

VALUE bio_gadget_fq1l_mt5(vSelf, vPattern, vMinLen)
     VALUE vSelf;
     VALUE vPattern;
     VALUE vMinLen;
{
  char regexs[BUFSIZE];
  regex_t regexc;
  unsigned long minlen;
  char line[BUFSIZE];
  regmatch_t match[1];
  char *acc;
  char *seq;
  char *sep;
  char *qual;
  unsigned long acclen;

  sprintf(regexs, "^[^\t]+\t(%s)", StringValueCStr(vPattern));
  regcomp(&regexc, regexs, REG_EXTENDED);

  minlen = NUM2INT(vMinLen);

  while(fgets(line, BUFSIZE, stdin) != NULL) {
    if(regexec(&regexc, line, 1, match, 0) != REG_NOMATCH) {
      acc = strtok(line, "\t");
      seq = strtok(NULL, "\t");
      sep = strtok(NULL, "\t");
      qual = strtok(NULL, "\t");
      acclen = strlen(acc);
      seq += match[0].rm_eo-acclen-1;
      qual += match[0].rm_eo-acclen-1;
      if(strlen(seq) >= minlen)
	printf("%s\t%s\t%s\t%s", acc, seq, sep, qual);
    }
  }

  regfree(&regexc);
  return Qnil;
}

VALUE bio_gadget_fq1l_nr_deg(vSelf, vCmd)
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
    if (strlen(seq) >= pseql) {
      printf("%s\t%s\t%s\t%s", acc, seq, sep, qual);
      strcpy(pseq, seq);
      pseql = strlen(seq);
    } else {
      sprintf(regexs, "^%s", seq);
      regcomp(&regexc, regexs, REG_NOSUB);
      if (regexec(&regexc, pseq, 0, NULL, 0) == REG_NOMATCH) {
	printf("%s\t%s\t%s\t%s", acc, seq, sep, qual);
	strcpy(pseq, seq);
	pseql = strlen(seq);
      }
    }
    regfree(&regexc);
  }
  fclose(fp);
  
  return Qnil;
}

VALUE bio_gadget_fq1l_nr_std(vSelf, vCmd)
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

VALUE bio_gadget_fq1l_pt3(vSelf, vLen, vCmdIn, vPathOut, vMinLen)
     VALUE vSelf;
     VALUE vLen;
     VALUE vCmdIn;
     VALUE vPathOut;
     VALUE vMinLen;
{
  char line[BUFSIZE];
  char *acc;
  char *seq;
  char *sep;
  char *qual;
  int len;
  unsigned long minlen;
  unsigned long seqlen;
  FILE *fp_in;
  FILE *fp_out;

  fp_in = popen(StringValueCStr(vCmdIn), "r");
  fp_out = fopen(StringValueCStr(vPathOut), "w");
  len = NUM2INT(vLen);
  minlen = NUM2INT(vMinLen);
  
  while(fgets(line, BUFSIZE, fp_in) != NULL) {
    acc = strtok(line, "\t");
    seq = strtok(NULL, "\t");
    sep = strtok(NULL, "\t");
    qual = strtok(NULL, "\t");
    seqlen = strlen(seq)-len;
    if (seqlen > 0 && seqlen >= minlen) {
      seq[seqlen] = 0;
      qual[seqlen] = 0;
      fprintf(fp_out, "%s\t%s\t%s\t%s\n", acc, seq, sep, qual);
    }
  }
  fclose(fp_out);
  fclose(fp_in);
  
  return Qnil;
}

VALUE bio_gadget_fq1l_qt3(vSelf, vLQs, vMinLen)
     VALUE vSelf;
     VALUE vLQs;
     VALUE vMinLen;
{
  char line[BUFSIZE];
  char *acc;
  char *seq;
  char *sep;
  char *qual;
  char *lqs;
  unsigned long minlen;
  unsigned long seqlen;

  lqs = StringValueCStr(vLQs);
  minlen = NUM2INT(vMinLen);

  while(fgets(line, BUFSIZE, stdin) != NULL) {
    acc = strtok(line, "\t");
    seq = strtok(NULL, "\t");
    sep = strtok(NULL, "\t");
    qual = strtok(NULL, "\t\n");
    seqlen = strcspn(qual, lqs);
    if (seqlen > 0 && seqlen >= minlen) {
      if (seqlen < strlen(qual)) {
	seq[seqlen] = 0;
	qual[seqlen] = 0;
      }
      printf("%s\t%s\t%s\t%s\n", acc, seq, sep, qual);
    }
  }
  
  return Qnil;
}

VALUE bio_gadget_fq1l_to(vSelf, vDraw, vSkip)
     VALUE vSelf;
     VALUE vDraw;
     VALUE vSkip;
{
  char line[BUFSIZE];
  unsigned int draw;
  unsigned int skip;
  unsigned int sum;
  unsigned long count;

  draw = NUM2INT(vDraw);
  skip = NUM2INT(vSkip);
  sum = draw + skip;
  count = 0;

  while(fgets(line, BUFSIZE, stdin) != NULL) {
    if (count % sum < draw)
      fputs(line, stdout);
    count += 1;
  }

  return Qnil;
}

VALUE rb_mBio_Gadget;

void
Init_bio_gadget(void)
{
  rb_mBio_Gadget = rb_define_module("BioGadget");
  rb_define_module_function(rb_mBio_Gadget, "mt5", bio_gadget_fq1l_mt5, 2);
  rb_define_module_function(rb_mBio_Gadget, "nr_deg", bio_gadget_fq1l_nr_deg, 1);
  rb_define_module_function(rb_mBio_Gadget, "nr_std", bio_gadget_fq1l_nr_std, 1);
  rb_define_module_function(rb_mBio_Gadget, "pt3", bio_gadget_fq1l_pt3, 4);
  rb_define_module_function(rb_mBio_Gadget, "qt3", bio_gadget_fq1l_qt3, 2);
  rb_define_module_function(rb_mBio_Gadget, "to", bio_gadget_fq1l_to, 2);
}
