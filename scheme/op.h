//	Op table entries
//
//	First col is used in enum.
//	Other cols used in populating op_code_table

  PROC_ENTRY0(LOAD,			op_load, "load", 1, 1)
  PRIM_ENTRY(T0LVL,			op_t0lvl, "t0lvl")
  PRIM_ENTRY(T1LVL,			op_t1lvl, "t1lvl")
  PRIM_ENTRY(READ_INTERNAL, op_read_internal, "read_internal")
  PROC_ENTRY0(GENSYM,		op_gensym, "gensym", 0, 0)
  PRIM_ENTRY(VALUEPRINT,	op_valueprint, "valueprint")
  PRIM_ENTRY(EVAL,			op_eval, "eval")
  PRIM_ENTRY(REAL_EVAL,		op_real_eval, "real_eval")
  PRIM_ENTRY(E0ARGS,		op_e0args, "e0args")
  PRIM_ENTRY(E1ARGS,		op_e1args, "e1args") 
  PRIM_ENTRY(APPLY,			op_apply, "apply")
  PRIM_ENTRY(REAL_APPLY,	op_real_apply, "real-apply")
  PROC_ENTRY1(TRACING,		op_tracing, "tracing", 1, 1, TST_NATURAL) 
  PRIM_ENTRY(DOMACRO,		op_domacro, "domacro")
  PRIM_ENTRY(LAMBDA,		op_lambda, "lambda")
  PROC_ENTRY2(MKCLOSURE,	op_mkclosure, "make-closure", 1, 2, TST_PAIR, TST_ENVIRONMENT)
  PRIM_ENTRY(QUOTE,			op_quote, "quote")
  PRIM_ENTRY(DEF0,			op_def0, "def0")
  PRIM_ENTRY(DEF1,			op_def1, "def1")
  PROC_ENTRY2(DEFP,			op_defp, "defined?", 1, 2, TST_SYMBOL, TST_ENVIRONMENT)
  PRIM_ENTRY(BEGIN,			op_begin, "begin")
  PRIM_ENTRY(IF0,			op_if0, "if0")
  PRIM_ENTRY(IF1,			op_if1, "if1")
  PRIM_ENTRY(SET0,			op_set0, "set0")
  PRIM_ENTRY(SET1,			op_set1, "set1")
  PRIM_ENTRY(LET0,			op_let0, "let0")
  PRIM_ENTRY(LET1,			op_let1, "let1")
  PRIM_ENTRY(LET2,			op_let2, "let2")
  PRIM_ENTRY(LET0STAR,		op_let0star, "let0star") 
  PRIM_ENTRY(LET1STAR,		op_let1star, "let1star")
  PRIM_ENTRY(LET2STAR,		op_let2star, "leg2star")
  PRIM_ENTRY(LET0REC,		op_let0rec, "let0rec") 
  PRIM_ENTRY(LET1REC,		op_let1rec, "let1rec")
  PRIM_ENTRY(LET2REC,		op_let2rec, "let2rec") 
  PRIM_ENTRY(COND0,			op_cond0, "cond0")
  PRIM_ENTRY(COND1,			op_cond1, "cond1")
  PRIM_ENTRY(DELAY,			op_delay, "delay") 
  PRIM_ENTRY(AND0,			op_and0, "and0")
  PRIM_ENTRY(AND1,			op_and1, "and1") 
  PRIM_ENTRY(OR0,			op_or0, "or0") 
  PRIM_ENTRY(OR1,			op_or1, "or1") 
  PRIM_ENTRY(C0STREAM,		op_c0stream, "c0stream") 
  PRIM_ENTRY(C1STREAM,		op_c1stream, "c1stream") 
  PRIM_ENTRY(MACRO0,		op_macro0, "macro0") 
  PRIM_ENTRY(MACRO1,		op_macro1, "macro1")
  PRIM_ENTRY(CASE0,			op_case0, "case0")
  PRIM_ENTRY(CASE1,			op_case1, "case1") 
  PRIM_ENTRY(CASE2,			op_case2, "case2") 
  PROC_ENTRY2(PEVAL,		op_peval, "eval", 1, 2, TST_ANY, TST_ENVIRONMENT)
  PROC_ENTRY0(PAPPLY,		op_papply, "apply", 1, INF_ARG)
  PROC_ENTRY0(CONTINUATION,	op_continuation, "call-with-current-continuation", 1, 1)
  PROC_ENTRY1(INEX2EX,		op_inex2ex, "inexact->exact", 1, 1, TST_NUMBER)
  PROC_ENTRY1(EX2INEX,		op_ex2inex, "exact->inexact", 1, 1, TST_NUMBER)
  PROC_ENTRY1(EXP,			op_exp, "exp", 1, 1, TST_NUMBER)
  PROC_ENTRY1(LOG,			op_log, "log", 1, 1, TST_NUMBER)
  PROC_ENTRY1(SIN,			op_sin, "sin", 1, 1, TST_NUMBER)
  PROC_ENTRY1(COS,			op_cos, "cos", 1, 1, TST_NUMBER)	
  PROC_ENTRY1(TAN,			op_tan, "tan", 1, 1, TST_NUMBER)
  PROC_ENTRY1(ASIN,			op_asin, "asin", 1, 1, TST_NUMBER)
  PROC_ENTRY1(ACOS,			op_acos, "acos", 1, 1, TST_NUMBER)
  PROC_ENTRY1(ATAN,			op_atan, "atan", 1, 2, TST_NUMBER)
  PROC_ENTRY1(SQRT,			op_sqrt, "sqrt", 1, 1, TST_NUMBER)
  PROC_ENTRY1(EXPT,			op_expt, "expt", 2, 2, TST_NUMBER)
  PROC_ENTRY1(FLOOR,		op_floor, "floor", 1, 1, TST_NUMBER)
  PROC_ENTRY1(CEILING,		op_ceiling, "ceiling", 1, 1, TST_NUMBER)
  PROC_ENTRY1(TRUNCATE,		op_truncate, "truncate", 1, 1, TST_NUMBER)
  PROC_ENTRY1(ROUND,		op_round, "round", 1, 1, TST_NUMBER)
  PROC_ENTRY1(EXACT,		op_exact, "exact?", 1, 1, TST_NUMBER)
  PROC_ENTRY1(INEXACT,		op_inexact, "inexact?", 1, 1, TST_NUMBER)
  PROC_ENTRY1(ODD,			op_odd, "odd?", 1, 1, TST_INTEGER)
  PROC_ENTRY1(EVEN,			op_even, "even?", 1, 1, TST_INTEGER)
  PROC_ENTRY1(ZERO,			op_zero, "zero?", 1, 1, TST_NUMBER)
  PROC_ENTRY1(POSITIVE,		op_positive, "positive?", 1, 1, TST_NUMBER)
  PROC_ENTRY1(NEGATIVE,		op_negative, "negative?", 1, 1, TST_NUMBER)
  PROC_ENTRY1(ADD,			op_add, "+", 0, INF_ARG, TST_NUMBER)
  PROC_ENTRY1(SUB,			op_sub, "-", 1, INF_ARG, TST_NUMBER)
  PROC_ENTRY1(MUL,			op_mul, "*", 0, INF_ARG, TST_NUMBER)
  PROC_ENTRY1(DIV,			op_div, "/", 1, INF_ARG, TST_NUMBER)
  PROC_ENTRY1(INTDIV,		op_intdiv, "quotient", 1, INF_ARG, TST_INTEGER)
  PROC_ENTRY1(REM,			op_rem, "remainder", 2, 2, TST_INTEGER)
  PROC_ENTRY1(MOD,			op_mod, "modulo", 2, 2, TST_INTEGER)
  PROC_ENTRY1(MAX,			op_max, "max", 1, INF_ARG, TST_NUMBER)
  PROC_ENTRY1(MIN,			op_min, "min", 1, INF_ARG, TST_NUMBER)
  PROC_ENTRY1(ABS,			op_abs, "abs", 1, 1, TST_NUMBER)
  PROC_ENTRY1(GCD,			op_gcd, "gcd", 0, INF_ARG, TST_NUMBER)
  PROC_ENTRY1(LCM,			op_lcm, "lcm", 0, INF_ARG, TST_NUMBER)
  PROC_ENTRY1(CAR,			op_car, "car", 1, 1, TST_PAIR) 
  PROC_ENTRY1(CDR,			op_cdr, "cdr", 1, 1, TST_PAIR)
  PROC_ENTRY0(CONS,			op_cons, "cons", 2, 2)
  PROC_ENTRY2(SETCAR,		op_setcar, "set-car!", 2, 2, TST_PAIR, TST_ANY)
  PROC_ENTRY2(SETCDR,		op_setcdr, "set-cdr!", 2, 2, TST_PAIR, TST_ANY)
  PROC_ENTRY1(CHAR2INT,		op_char2int, "char->integer", 1, 1, TST_CHAR)
  PROC_ENTRY1(INT2CHAR,		op_int2char, "integer->char", 1, 1, TST_NATURAL)
  PROC_ENTRY1(CHARUPCASE,	op_charupcase, "char-upcase", 1, 1, TST_CHAR)
  PROC_ENTRY1(CHARDNCASE,	op_chardncase, "char-downcase", 1, 1, TST_CHAR)
  PROC_ENTRY1(SYM2STR,		op_sym2str, "symbol->string", 1, 1, TST_SYMBOL)
  PROC_ENTRY1(ATOM2STR,		op_atom2str, "atom->string", 1, 1, TST_ANY)
  PROC_ENTRY1(STR2SYM,		op_str2sym, "string->symbol", 1, 1, TST_STRING)
  PROC_ENTRY1(STR2ATOM,		op_str2atom, "string->atom", 1, 1, TST_STRING)
  PROC_ENTRY2(MKSTRING,		op_mkstring, "make-string", 1, 2, TST_NATURAL, TST_CHAR)
  PROC_ENTRY1(STRLEN,		op_strlen, "string-length", 1, 1, TST_STRING) 
  PROC_ENTRY2(STRREF,		op_strref, "string-ref", 2, 2, TST_STRING, TST_NATURAL) 
  PROC_ENTRY3(STRSET,		op_strset, "string-set!", 3, 3, TST_STRING, TST_NATURAL, TST_CHAR)
  PROC_ENTRY2(SUBSTR,		op_substr, "substring", 2, 3, TST_STRING, TST_NATURAL)
  PROC_ENTRY1(SAPPEND,		op_strappend, "string-append", 0, INF_ARG, TST_STRING)
  PROC_ENTRY1(STR2LIST,		op_str2list, "string->list", 1, 1, TST_STRING)
  PROC_ENTRY1(STRING,		op_string, "string", 0, INF_ARG, TST_CHAR)
  PROC_ENTRY1(LIST2STR,		op_list2str, "list->string", 0, INF_ARG, TST_LIST)
  PROC_ENTRY2(STR2NUM,		op_str2num, "string->number", 1, 2, TST_STRING, TST_NUMBER)
  PROC_ENTRY1(NUM2STR,		op_num2str, "number->string", 1, 2, TST_NUMBER)
  PROC_ENTRY2(STRFILL,		op_strfill, "string-fill!", 2, 2, TST_STRING, TST_CHAR)
  PROC_ENTRY0(LIST,			op_list, "list", 0, INF_ARG)
  PROC_ENTRY2(LISTTAIL,		op_listtail, "list-tail", 2, 2, TST_LIST, TST_NUMBER)
  PROC_ENTRY2(LISTREF,		op_listref, "list-ref", 2, 2, TST_LIST, TST_NUMBER)
  PROC_ENTRY2(MEMQ,			op_memq, "memq", 2, 2, TST_ANY, TST_LIST)
  PROC_ENTRY2(MEMV,			op_memv, "memv", 2, 2, TST_ANY, TST_LIST)
  PROC_ENTRY2(MEMBER,		op_member, "member", 2, 2, TST_ANY, TST_LIST)
  PROC_ENTRY0(ASSQ,			op_assq, "assq", 2, 2)
  PROC_ENTRY0(ASSV,			op_assv, "assv", 2, 2)
  PROC_ENTRY0(ASSOC,		op_assoc, "assoc", 2, 2)
  PROC_ENTRY0(VECTOR,		op_vector, "vector", 0, INF_ARG) 
  PROC_ENTRY2(MKVECTOR,		op_mkvector, "make-vector", 1, 2, TST_NATURAL, TST_ANY) 
  PROC_ENTRY1(VECLEN,		op_veclen, "vector-length", 1, 1, TST_VECTOR)
  PROC_ENTRY2(VECREF,		op_vecref, "vector-ref", 2, 2, TST_VECTOR, TST_NATURAL)
  PROC_ENTRY3(VECSET,		op_vecset, "vector-set!", 3, 3, TST_VECTOR, TST_NATURAL, TST_ANY)
  PROC_ENTRY1(LIST2VEC,		op_list2vec, "list->vector", 1, 1, TST_LIST)
  PROC_ENTRY2(VECFILL,		op_vecfill, "vector-fill!", 2, 2, TST_VECTOR, TST_ANY)
  PROC_ENTRY1(VEC2LIST,		op_vec2list, "vector->list", 1, 1, TST_VECTOR)
  PROC_ENTRY0(NOT,			op_not, "not", 1, 1)
  PROC_ENTRY0(BOOLP,		op_boolp, "boolean?", 1, 1)
  PROC_ENTRY0(EOFOBJP,		op_eofobjp, "eof-object?", 1, 1)
  PROC_ENTRY0(NULLP,		op_nullp, "null?", 1, 1)
  PROC_ENTRY1(NUMEQ,		op_comp, "=", 2, INF_ARG, TST_NUMBER)
  PROC_ENTRY1(LESS,			op_comp, "<", 2, INF_ARG, TST_NUMBER)
  PROC_ENTRY1(GRE,			op_comp, ">", 2, INF_ARG, TST_NUMBER)
  PROC_ENTRY1(LEQ,			op_comp, "<=", 2, INF_ARG, TST_NUMBER)
  PROC_ENTRY1(GEQ,			op_comp, ">=", 2, INF_ARG, TST_NUMBER)
  PROC_ENTRY0(SYMBOLP,		op_symbolp, "symbol?", 1, 1)
  PROC_ENTRY0(NUMBERP,		op_numberp, "number?", 1, 1)
  PROC_ENTRY0(STRINGP,		op_stringp, "string?", 1, 1)
  PROC_ENTRY0(INTEGERP,		op_integerp, "integer?", 1, 1)
  PROC_ENTRY0(REALP,		op_realp, "real?", 1, 1) 
  PROC_ENTRY0(CHARP,		op_charp, "char?", 1, 1)
#if USE_CHAR_CLASSIFIERS
  PROC_ENTRY1(CHARAP,		op_charap, "char-alphabetic?", 1, 1, TST_CHAR)
  PROC_ENTRY1(CHARNP,		op_charnp, "char-numeric?", 1, 1, TST_CHAR)
  PROC_ENTRY1(CHARWP,		op_charwp, "char-whitespace?", 1, 1, TST_CHAR)
  PROC_ENTRY1(CHARUP,		op_charup, "char-upper-case?", 1, 1, TST_CHAR)
  PROC_ENTRY1(CHARLP,		op_charlp, "char-lower-case?", 1, 1, TST_CHAR)
#endif
  PROC_ENTRY0(PORTP,		op_portp, "port?", 1, 1)
  PROC_ENTRY0(INPORTP,		op_inportp, "input-port?", 1, 1)
  PROC_ENTRY0(OUTPORTP,		op_outportp, "output-port?", 1, 1)	
  PROC_ENTRY0(PROCP,		op_procp, "procedure?", 1, 1)
  PROC_ENTRY0(PAIRP,		op_pairp, "pair?", 1, 1)
  PROC_ENTRY0(LISTP,		op_listp, "list?", 1, 1)
  PROC_ENTRY0(ENVP,			op_envp, "environment?", 1, 1)
  PROC_ENTRY0(VECTORP,		op_vectorp, "vector?", 1, 1)
  PROC_ENTRY0(EQ,			op_eq, "eq?", 2, 2)
  PROC_ENTRY0(EQV,			op_eqv, "eqv?", 2, 2)
  PROC_ENTRY0(EQUAL,		op_equal, "equal?", 2, 2)
  PROC_ENTRY0(VECEQUAL,		op_vecequal, "vector-equal?", 2, 2)
  PROC_ENTRY0(FORCE,		op_force, "force", 1, 1)
  PRIM_ENTRY(SAVE_FORCED,	op_save_forced, "save_forced")
  PROC_ENTRY2(WRITE,		op_write, "write", 1, 2, TST_ANY, TST_OUTPORT)
  PROC_ENTRY2(WRITE_CHAR,	op_write_char, "write-char", 1, 2, TST_CHAR, TST_OUTPORT)
  PROC_ENTRY2(DISPLAY,		op_display, "display", 1, 2, TST_ANY, TST_OUTPORT)
  PROC_ENTRY1(NEWLINE,		op_newline, "newline", 0, 1, TST_OUTPORT)
  PROC_ENTRY0(ERR0,			op_err0, "error", 1, INF_ARG)
  PRIM_ENTRY(ERR1,			op_err1, "err1")
  PROC_ENTRY1(REVERSE,		op_reverse, "reverse", 1, 1, TST_PAIR)
  PROC_ENTRY1(LIST_STAR,	op_list_star, "list*", 1, INF_ARG, TST_NONE)
  PROC_ENTRY0(APPEND,		op_append, "append", 0, INF_ARG)
  PROC_ENTRY0(PUT,			op_put, "put", 3, 3)
  PROC_ENTRY0(GET,			op_get, "get", 2, 2)
  PROC_ENTRY1(QUIT,			op_quit, "quit", 0, 1, TST_NUMBER)
  PROC_ENTRY0(GC,			op_gc, "gc", 0, 0)
  PROC_ENTRY0(GCVERB,		op_gcverb, "gc-verbose", 0, 1)
  PROC_ENTRY1(NEWSEGMENT,	op_newsegment, "new-segment", 0, 1, TST_NUMBER) 
  PROC_ENTRY0(OBLIST,		op_oblist, "oblist", 0, 0)
  PROC_ENTRY0(CURR_INPORT,	op_curr_inport, "current-input-port", 0, 0)
  PROC_ENTRY0(CURR_OUTPORT,	op_curr_outport, "current-output-port", 0, 0)	
  PROC_ENTRY1(OPEN_INFILE,	op_open_infile, "open-input-file", 1, 1, TST_STRING)
  PROC_ENTRY1(OPEN_OUTFILE,	op_open_outfile, "open-output-file", 1, 1, TST_STRING)
  PROC_ENTRY1(OPEN_INOUTFILE, op_open_inoutfile, "open-input-output-file", 1, 1, TST_STRING)
  PROC_ENTRY1(OPEN_INSTRING, op_open_instring, "open-input-string", 1, 1, TST_STRING)
  PROC_ENTRY1(OPEN_OUTSTRING, op_open_outstring, "open-output-string", 1, 1, TST_STRING)
  PROC_ENTRY1(OPEN_INOUTSTRING, op_open_inoutstring, "open-input-output-string", 1, 1, TST_STRING)
  PROC_ENTRY1(CLOSE_INPORT,	op_close_inport, "close-input-port", 1, 1, TST_INPORT) 
  PROC_ENTRY1(CLOSE_OUTPORT,op_close_outport, "close-output-port", 1, 1, TST_OUTPORT)
  PROC_ENTRY0(INT_ENV,		op_int_env, "interaction-environment", 0, 0)
  PROC_ENTRY0(CURR_ENV,		op_curr_env, "current-environment", 0, 0)	
  PROC_ENTRY1(READ,			op_read, "read", 0, 1, TST_INPORT)
  PROC_ENTRY1(READ_CHAR,	op_read_char, "read-char", 0, 1, TST_INPORT)
  PROC_ENTRY1(PEEK_CHAR,	op_read_char, "peek-char", 0, 1, TST_INPORT)	//  same as read_char
  PROC_ENTRY1(CHAR_READY,	op_char_ready, "char-ready?", 0, 1, TST_INPORT)
  PROC_ENTRY1(SET_INPORT,	op_set_inport, "set-input-port!", 1, 1, TST_INPORT)
  PROC_ENTRY1(SET_OUTPORT,	op_set_outport, "set-output-port!", 1, 1, TST_OUTPORT)
  PRIM_ENTRY(RDSEXPR,		op_rdsexpr, "rdsexpr")
  PRIM_ENTRY(RDLIST,		op_rdlist, "rdlist")
  PRIM_ENTRY(RDDOT,			op_rddot, "rddot")
  PRIM_ENTRY(RDQUOTE,		op_rdquote, "rdquote")
  PRIM_ENTRY(RDQQUOTE,		op_rdqquote, "rdqquote") 
  PRIM_ENTRY(RDQQUOTEVEC,	op_rdqquotevec, "rdqquotevec") 
  PRIM_ENTRY(RDUNQUOTE,		op_rdunquote, "rdunquote") 
  PRIM_ENTRY(RDUQTSP,		op_rduqtsp, "rdqutsp")
  PRIM_ENTRY(RDVEC,			op_rdvec, "rdvec")
  PRIM_ENTRY(P0LIST,		op_p0list, "p0list")
  PRIM_ENTRY(P1LIST,		op_p1list, "p1list") 
  PRIM_ENTRY(PVECFROM,		op_pvecfrom, "pvecfrom")
  PROC_ENTRY1(LIST_LENGTH,	op_list_length, "length", 1, 1, TST_LIST)
  PROC_ENTRY0(GET_CLOSURE,	op_get_closure, "get-closure-code", 1, 1)
  PROC_ENTRY0(CLOSUREP,		op_closurep, "closure?", 1, 1)
  PROC_ENTRY0(TIME	,		op_time, "time", 0, 0)
  PRIM_ENTRY(ILLEGAL,		op_illegal, "illegal")
