
#include <htslib/hts.h>
#include <htslib/sam.h>

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>

enum oper {

  oper_eq,
  oper_gt,
  oper_ge,
  oper_lt,
  oper_le,
  oper_ne

};

static void parse_const_or_attr(char* arg, int* constout, char** attrout) {
  
  if(isdigit(arg[0])) {
    *constout = atoi(arg);
    *attrout = 0;
  }
  else {
    *constout = 0;
    *attrout = arg;
  }

}

static oper parse_operator(const char* arg) {

  if(!strcmp(arg, "=="))
    return oper_eq;
  else if(!strcmp(arg, "<"))
    return oper_lt;
  else if(!strcmp(arg, ">"))
    return oper_gt;
  else if(!strcmp(arg, "<="))
    return oper_le;
  else if(!strcmp(arg, ">="))
    return oper_ge;
  else if(!strcmp(arg, "!="))
    return oper_ne;
  
  fprintf(stderr, "Invalid operator %s\n", arg);
  exit(1);

}

static int get_const_or_attr(int constin, const char* attrin, bam1_t* rec) {

  if(attrin) {
    uint8_t* attr_rec = bam_aux_get(rec, attrin);
    if(!attr_rec) {
      fprintf(stderr, "Fatal: At least record %s doesn't have a %s tag as required.\n", bam_get_qname(rec), attrin);
      exit(1);
    }
    return bam_aux2i(attr_rec);
  }
  
  return constin;
  
}

static bool eval_test(int op1, oper op, int op2) {

  switch(op) {
  case oper_eq:
    return op1 == op2;
  case oper_lt:
    return op1 < op2;
  case oper_gt:
    return op1 > op2;
  case oper_le:
    return op1 <= op2;
  case oper_ge:
    return op1 >= op2;
  case oper_ne:
    return op1 != op2;
  }

}

int main(int argc, char** argv) {

  if(argc != 4) {

    fprintf(stderr, "Usage: filter_attr attr_or_constant relation attr_or_constant\n");
    exit(1);

  }

  int consta, constb;
  char *attra, *attrb;
  oper op;

  parse_const_or_attr(argv[1], &consta, &attra);
  parse_const_or_attr(argv[3], &constb, &attrb);

  op = parse_operator(argv[2]);

  htsFile* hfi = hts_open("-", "r");
  htsFile* hfo = hts_open("-", "wb");

  if((!hfi) || (!hfo)) {

    fprintf(stderr, "Failed to open stdin/out\n");
    exit(1);

  }

  bam_hdr_t* header = sam_hdr_read(hfi);
  sam_hdr_write(hfo, header);

  bam1_t* rec = bam_init1();

  int total = 0;
  int kept = 0;

  while(sam_read1(hfi, header, rec) >= 0) {

    int op1 = get_const_or_attr(consta, attra, rec);
    int op2 = get_const_or_attr(constb, attrb, rec);

    ++total;

    if(eval_test(op1, op, op2)) {
      sam_write1(hfo, header, rec);
      ++kept;
    }

  }

  hts_close(hfi);
  hts_close(hfo);

  fprintf(stderr, "%d / %d records retained\n", kept, total);

}
