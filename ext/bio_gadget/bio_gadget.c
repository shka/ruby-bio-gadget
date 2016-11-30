#include <limits.h>
#include <regex.h>
#include <stdio.h>
#include <string.h>
#include "bio_gadget.h"

VALUE bio_gadget_fq1l_i2i(vSelf, vFirst, vLast)
     VALUE vSelf;
     VALUE vFirst;
     VALUE vLast;
{
  char line[BUFSIZE];
  char index[BUFSIZE] = "";
  char *acc;
  char *seq;
  char *sep;
  char *qual;
  unsigned long head;
  unsigned long length;

  head = NUM2INT(vFirst)-1;
  length = NUM2INT(vLast)-head;
  while(fgets(line, BUFSIZE, stdin) != NULL) {
    acc = strtok(line, "\t");
    seq = strtok(NULL, "\t");
    sep = strtok(NULL, "\t");
    qual = strtok(NULL, "\t");
    strcpy(index, seq+head);
    index[length] = 0;
    printf("%s%s\t%s\t%s\t%s", acc, index, seq, sep, qual);
  }

  return Qnil;
}

VALUE bio_gadget_fq1l_nr_deg(vSelf)
     VALUE vSelf;
{
  char line[BUFSIZE];
  char *acc;
  char *seq;
  char *sep;
  char *qual;
  char pseq[BUFSIZE] = "";
  unsigned long pseql = ULONG_MAX;
  char regexs[BUFSIZE];
  regex_t regexc;

  while(fgets(line, BUFSIZE, stdin) != NULL) {
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
  
  return Qnil;
}

VALUE bio_gadget_fq1l_nr_std(vSelf)
     VALUE vSelf;
{
  char line[BUFSIZE];
  char *acc;
  char *seq;
  char *sep;
  char *qual;
  char pseq[BUFSIZE] = "";

  while(fgets(line, BUFSIZE, stdin) != NULL) {
    acc = strtok(line, "\t");
    seq = strtok(NULL, "\t");
    sep = strtok(NULL, "\t");
    qual = strtok(NULL, "\t");
    if (strcmp(pseq, seq) != 0) {
      printf("%s\t%s\t%s\t%s", acc, seq, sep, qual);
      strcpy(pseq, seq);
    }
  }
  
  return Qnil;
}

VALUE bio_gadget_fq1l_t3(vSelf, vCmdIn, vLen, vMinLen, vPathOut)
     VALUE vSelf;
     VALUE vCmdIn;
     VALUE vLen;
     VALUE vMinLen;
     VALUE vPathOut;
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

  fp_in = RTEST(vCmdIn) ? popen(StringValueCStr(vCmdIn), "r") : stdin;
  fp_out = RTEST(vPathOut) ? fopen(StringValueCStr(vPathOut), "w") : stdout;
  len = NUM2INT(vLen);
  minlen = NUM2INT(vMinLen);
  
  while(fgets(line, BUFSIZE, fp_in) != NULL) {
    acc = strtok(line, "\t");
    seq = strtok(NULL, "\t");
    seqlen = strlen(seq)-len;
    if (seqlen > 0 && seqlen >= minlen) {
      sep = strtok(NULL, "\t");
      qual = strtok(NULL, "\t");
      seq[seqlen] = 0;
      qual[seqlen] = 0;
      fprintf(fp_out, "%s\t%s\t%s\t%s\n", acc, seq, sep, qual);
    }
  }
  fclose(fp_out);
  fclose(fp_in);
  
  return Qnil;
}

VALUE bio_gadget_fq1l_t3q(vSelf, vLQs, vMinLen)
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

VALUE bio_gadget_fq1l_t5(vSelf, vPattern, vMinLen)
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

VALUE bio_gadget_fq1l_u2i(vSelf, vFirst, vLast)
     VALUE vSelf;
     VALUE vFirst;
     VALUE vLast;
{
  char line[BUFSIZE];
  char index[BUFSIZE] = "";
  char *acc;
  char *acc1;
  char *acc2;
  char *seq;
  char *sep;
  char *qual;
  unsigned long head;
  unsigned long length;

  head = NUM2INT(vFirst)-1;
  length = NUM2INT(vLast)-head;
  while(fgets(line, BUFSIZE, stdin) != NULL) {
    acc = strtok(line, "\t");
    seq = strtok(NULL, "\t");
    sep = strtok(NULL, "\t");
    qual = strtok(NULL, "\t");
    strcpy(index, seq+head);
    index[length] = 0;
    acc1 = strtok(acc, " ");
    acc2 = strtok(NULL, " ");
    if(acc2 == NULL) {
      printf("%s:%s\t%s\t%s\t%s", acc1, index, seq, sep, qual);
    }
    else {
      printf("%s:%s %s\t%s\t%s\t%s", acc1, index, acc2, seq, sep, qual);
    }
  }

  return Qnil;
}


VALUE rb_mBio_Gadget;

void
Init_bio_gadget(void)
{
  rb_mBio_Gadget = rb_define_module("BioGadget");
  rb_define_module_function(rb_mBio_Gadget, "i2i", bio_gadget_fq1l_i2i, 2);
  rb_define_module_function(rb_mBio_Gadget, "nr_deg", bio_gadget_fq1l_nr_deg, 0);
  rb_define_module_function(rb_mBio_Gadget, "nr_std", bio_gadget_fq1l_nr_std, 0);
  rb_define_module_function(rb_mBio_Gadget, "t3", bio_gadget_fq1l_t3, 4);
  rb_define_module_function(rb_mBio_Gadget, "t3q", bio_gadget_fq1l_t3q, 2);
  rb_define_module_function(rb_mBio_Gadget, "t5", bio_gadget_fq1l_t5, 2);
  rb_define_module_function(rb_mBio_Gadget, "to", bio_gadget_fq1l_to, 2);
  rb_define_module_function(rb_mBio_Gadget, "u2i", bio_gadget_fq1l_u2i, 2);
}
