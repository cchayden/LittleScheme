// LITTLESCHEME.H 

#ifndef LITTLESCHEME_H
#define LITTLESCHEME_H

#include <stdio.h>

class scheme;
class cell;
typedef cell* pointer;

typedef pointer (*foreign_func)(scheme&, pointer);

class scheme_factory {
public:
	virtual scheme* create(int first_cell_seg, int max_cell_seg);
};

///	External interface to scheme.
class LittleScheme {
public:
	///	Construct an interpreter.
	LittleScheme(scheme_factory& sf,
		int first_cell_seg,		///<	Initial number of cell segments to allocate.	
		int max_cell_seg			///<	Maximum number of cell segments to allocate.
		);
	///	Destroy interpreter.
	~LittleScheme();
	///	Returns true if the interpreter is valid.
	/**	If the interpreter runs out of memory or cannot function for
	 *	some other reason, this returns false.
	 */
	bool is_valid();

	static const char* banner() { return "LittleScheme 0.9\n"; }

	///	Set the input port to the given FILE*.
	void set_input_port(FILE *fin);
	///	Set the input port to the given string.
	void set_input_port(char *start,	///<	The input string.
		char *past_the_end				///<	Points one past the string end.
		);
	///	Set the otput port to the given FILE*.
	void set_output_port(FILE *fin);
	///	Set the output port to the given string.
	void set_output_port(char *start,	///<	The output string
		char *past_the_end				///<	Points one past the string end.
		);
	///	Load from the given FILE*.
	//*	Returns the retcode */
	int load(FILE *fin);
	///	Load from the given string.
	/**	Returns the retcode. */
	int load(char *cmd);
	///	Apply ???
	int apply0(const char *procname);
	///	Call the ???.
	int call(pointer func, pointer args);
	///	Set data ???
	void set_external_data(void* p);
	///	Perform a definition in the given environment.
	/**	Define the given symbol to have the given value.
	 *	If no environment given, the definition is made in the global
	 *	environent.
	 */
	void define(pointer symbol, pointer value, pointer env = 0);

	///	Create a pair cell.
	/**	Obtain a free cell, make it a pair cell, and tores the given
	 *	pointers in its cons and cdr.
	 *	If immutable is true, then set the immutable flag on the cell.
	 *	Returns the cons cell.
	 */
	pointer cons(pointer a, pointer b, bool immutable = false);
	///	Create an immutable pair.
	pointer immutable_cons(pointer a, pointer b);

	///	Create an integer cell with the given value.
	pointer mk_integer(long num);
	///	Create a real cell with the given value.
	pointer mk_real(double num);
	/// Make symbol or number, depending on input.
	pointer mk_atom(const char *str,	///<	The name or number.
		bool must_be_number = false,	///<	If true, it must be a number.
		int radix = 10					///<	If a number, the number base.
		);
	///	Create a new internal symbol.
	/**	Returns the symbol cell. */
	pointer gensym();
	///	Create a string cell that holds (a copy of) the given string.
	/**	The string must be null terminated. */
	pointer mk_string(const char* str);
	///	Create a string cell that holds (a copy of) the given string.
	/**	The string is defined by its start and length, and may contain nulls. */
	pointer mk_string(const char* str, int len);
	///	Create a string cell that holds a string of the given length.
	/**	The string is initially filled with the fill character. */
	pointer mk_string(int len, char fill = ' ');
	///	Create a character cell.
	pointer mk_character(int c);
	///	Create a vector cell, with the given length.
	pointer mk_vector(int len);
	///	Create a foreign function cell.
	pointer mk_foreign_func(foreign_func f);
	///	Create a new environment frame, linking with the old environent.
	pointer mk_environment(pointer new_frame, pointer old_env);

	///	Create a symbol cell corresponding to the given name.
	pointer define_symbol(const char* name);

	///	The (single) cell '().
	static pointer NIL();
	///	The (single) cell #t.
	static pointer T();
	///	The (single) cell #f.
	static pointer F();
	///	Returns true if the pointer is not #f.
	static bool is_true(pointer p);
	///	Returns true if the pointer is #f.
	/** () is #t in R5RS. */
	static bool is_false(pointer p);

	///	Returns the global environment.
	pointer global_env() const;
	///	Returns the retcode.
	/**	retcode values: <br>
	 *	   0	normal <br>
	 *	  -1	called (error) <br>
	 *	   1	load error <br>
	 *	other	called (quit n) for some n
	 */
	int retcode() const;

	///	Return the reverse if a list.
	pointer reverse(pointer a);
	///	Returns the list a pappended with b.
	pointer append(pointer a, pointer b);

	///	Create a symbol for this name, add it to the environmet, and set its syntax flag.
	/**	This is used to associate a name with a special form.
	 *	Return the symbol cell.
	 */
	pointer assign_syntax(const char* name);
private:
	scheme* sc_;
};

pointer reverse_in_place(pointer list, pointer term);


#endif