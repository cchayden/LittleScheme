
#include "LittleScheme.h"

#include "scheme1.h"


#include <stdlib.h>
#include <string.h>
#include <io.h>

#ifndef InitFile
# define InitFile "init.scm"
#endif

static bool streq(char* a, char*b) { return strcmp(a, b) == 0; }

enum { FIRST_CELLSEGS = 3 };	 // # of cell segments to allocate at first
enum { CELL_NSEGMENT =  10 };    // # of segments for cells 

class scheme_factory_x : public scheme_factory {
public:
	virtual scheme* create(int first_cell_seg, int max_cell_seg){
		return new scheme1(first_cell_seg, max_cell_seg);
	}
};
static scheme_factory_x sf;

int main(int argc, char* *argv) {
	char* file_name = InitFile;
	int retcode;
	int n_segment = CELL_NSEGMENT;

	LittleScheme sc(sf, FIRST_CELLSEGS, n_segment);
	if(!sc.is_valid()) {
		fprintf(stderr, "Could not initialize!\n");
		return 2;
	}
	if(argc == 1) {
		printf(sc.banner());
	}
	if(argc == 2 && streq(argv[1], "-?")) {
		printf("Usage: %s [-? | <file1> <file2> ... | -1 <file> <arg1> <arg2> ...]\n"
			"\tUse - as filename for stdin.\n", argv[0]);
		return 1;
	}
	sc.set_input_port(stdin);
	sc.set_output_port(stdout);
	argv++;
	if(_access(file_name, 0) != 0) {
		char* p = getenv("LITTLESCHEMEINIT");
		if(p != 0) {
			file_name = p;
		}
	}
	do {
		//	treat arg as filename and load it
		FILE* fin;
		bool isfile = true;
		if(streq(file_name, "-")) {
			//	filename - is stdin
			fin = stdin;
		} else if(streq(file_name, "-1") || streq(file_name, "-c")) {
			//	filename is after the flag
			isfile = file_name[1] == '1';
			//	filename is actually the next arg 
			file_name = *argv++;
			if(streq(file_name, "-")) {
				fin = stdin;
			} else if(isfile) {
				fin = fopen(file_name, "r");
			}
			pointer args = sc.NIL();
			//	grab rest of args and put in *args*
			for( ; *argv; argv++) {
				pointer value = sc.mk_string(*argv);
				args = sc.cons(value,args);
			}
			reverse_in_place(args, 0);
			sc.define(sc.define_symbol("*args*"), args, 0);

		} else {
			//	filename is the arg
			fin = fopen(file_name, "r");
		}
		if(isfile && fin == 0) {
		  fprintf(stderr, "Could not open file %s\n",file_name);
		} else {
		  if(isfile) {
			retcode = sc.load(fin);
		  } else {
			retcode = sc.load(file_name);
		  }
		  if(!isfile || fin != stdin) {
			if(sc.retcode() != 0) {
			  fprintf(stderr, "s_errors encountered reading %s\n", file_name);
			}
			if(isfile) {
			  fclose(fin);
			}
		  }
		}
		file_name = *argv++;
	} while(file_name != 0);

	if(argc == 1) {
		//	this runs the REPL loop on stdin
		retcode = sc.load(stdin);
	}
  
	return retcode;
}

