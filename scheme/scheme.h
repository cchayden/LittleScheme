// SCHEME.H 

#ifndef _SCHEME_H
#define _SCHEME_H

#include "glob.h"
#include "Callback_1_Void.h"

#include <stdio.h>
#include <malloc.h>
#include <memory.h>
#include <string.h>

#pragma warning(disable: 4786)

#include <string>
#include <vector>
#include <map>


#define FAST 1
#define CHECK 1

/*
 * Default values for #define'd symbols
 */
#if USE_NO_FEATURES
# define USE_CHAR_CLASSIFIERS 0
# define USE_ASCII_NAMES 0
# define USE_ERROR_HOOK 0
# define USE_TRACING 0
# define USE_COLON_HOOK 0
#endif

#ifndef USE_CHAR_CLASSIFIERS  // If char classifiers are needed 
# define USE_CHAR_CLASSIFIERS 1
#endif

#ifndef USE_ASCII_NAMES  // If extended escaped characters are needed 
# define USE_ASCII_NAMES 1
#endif

#ifndef USE_TRACING
# define USE_TRACING 1
#endif

// To force system errors through user-defined error handling (see *error-hook*) 
#ifndef USE_ERROR_HOOK
# define USE_ERROR_HOOK 1
#endif

#ifndef USE_COLON_HOOK   // Enable qualified qualifier 
# define USE_COLON_HOOK 1
#endif

#if CHECK
#include "assert.h"
#endif

class scheme;
class cell;
typedef cell* pointer;

typedef pointer (*foreign_func)(scheme&, pointer);


//-///////////
//-/  num  ///
//-///////////

///	Num does arithmetic with integers and reals, converting between the two as necessary.
/**	A num contains either a long integer or a double real.
  *	It supports a variety of arithmetic operators.
  *	Arithmetic between integers generally produces integers, but if one
  *	of the operands is real, then the result it real.
  *	The current value can be obtained as either an integer or real.
  */
class num {
public:
	///	Construct default num as integer 0.
	num() : is_fixnum_(true) { value.ivalue_ = 0; }
	///	Construct integer num with given value.
	num(long ivalue) : is_fixnum_(true) { value.ivalue_ = ivalue; }
	///	Construct real num with given value.
	num(double rvalue) : is_fixnum_(false) { value.rvalue_ = rvalue; }
	///	Construct num given pointer to integer or real cell.
	num(pointer p);
	///	Return true if num is an integer.
	bool is_integer() const { return is_fixnum_; }
	///	Return the integer value.
	/**	This should only be used after verifying that the num holds
	 *	an integer, by testing is_integer().
	 */
	long ivalue() const { return value.ivalue_; }
	///	Return the real value.
	/**	This should only be used after verifying that the num holds
	 *	a real, by testing is_integer().
	 */
	double rvalue() const { return value.rvalue_; }
	///	Returns the num as an integer.  If necessary, the real value is truncated.
	long as_integer() const { return (is_fixnum_) ? value.ivalue_ : (long)value.rvalue_; }
	//	Returns the num as a real.  If necessary, it is converted to real.
	double as_real() const { return (is_fixnum_) ? (double)value.ivalue_ : value.rvalue_; }
	///	Sets the num to an integer value.
	void set(long ivalue) { 
		is_fixnum_ = true;
		value.ivalue_ = ivalue; 
	}
	///	Sets the num to a real value.
	void set(double rvalue) {
		is_fixnum_ = false;
		value.rvalue_ = rvalue;
	}
	///	Add the num to the one given.
	num add(num b);
	///	Subtract the givn number from the num.
	num sub(num b);
	///	Multiply the num by the one given.
	num mul(num b);
	///	Divide the num by the one given.
	num div(num b);
	///	Divide the num by the one given, but if both nums are integer, do an integer divide.
	num intdiv(num b);
	///	Return the remainder of the num divided by the one given.
	num rem(num b);
	///	Return the num modulo the one given.
	num mod(num b);
	///	Return the larger of the num and the one given.
	num max(num b);
	///	Return the minimum of the num and the one given.
	num min(num b);
	///	Return the absolute value of the num.
	num abs();
	///	Return true if the num is equal to the one given.
	bool eq(num b);
	///	Return true if the num is greater than the one given.
	bool gt(num b);
	///	Return true if the num is less than the one given.
	bool lt(num b);
	///	Return true if the num is greater than 0.
	bool positive();
	///	Return true if the num is less than 0.
	bool negative();
	///	Return true if the num is zero.
	bool zero();

private:
	//	True if integer (use ivalue_).
	//	False if real (use rvalue_).
	bool is_fixnum_;
    union {
          long ivalue_;
          double rvalue_;
    } value;
};


//-////////////
//-/  port  ///
//-////////////

///	Port holds information allowing input or output of characters.
/**	A port holds either (1) a FILE*, or (2) a string.
 *	It can do character input and output.
 */
class port {
public:
	///	Allocate and initialize a port as a FILE port.
	static port* make_port(FILE* f, int prop, bool closeit);
	///	Allocate and initialize a port as a string port.
	static port* make_port(char* start, char* past_the_end, int prop);
	///	Return a port to the free list.
	void free();
	///	Delete all ports on the free list at shutdown time.
	static void cleanup();

	///	Default constructor.
	port() : kind_(port_free) {}
	///	Constructor.
	port(FILE* f, int prop, bool closeit) : 
		kind_(port_file | prop) {	
			rep.stdio.file_ = f;
			rep.stdio.closeit_ = closeit;
		}
	port(char* start, char* past_the_end, int prop) :
		kind_(prop | port_string) {
			rep.string.start_ = start;
			rep.string.past_the_end_ = past_the_end;
			rep.string.curr_ = start;
		}
	///	These describe the kinds of ports.
	/**	These can be used in combination. */
	enum { 
		port_free =		0x00,	///<	port not in use
		port_file =		0x01,	///<	file port
		port_string =	0x02,	///<	string port
		port_input =	0x10,	///<	input port
		port_output =	0x20	///<	output port
	};
	///	Return the port kind.
	bool is_kind(int kind) const { return !!(kind_ & kind); }
	///	Return true if the port is a file port.
	bool is_file() const { return is_kind(port_file); }
	///	Return true if the port is a string port.
	bool is_string() const { return is_kind(port_string); }
	///	Return true if the port is an input port.
	bool is_input() const { return is_kind(port_input); }
	///	Return true if the port is an output port.
	bool is_output() const { return is_kind(port_output); }

	//	Makes the port a file port, and sets up the attributes.
	void set(int prop, FILE* file, bool closeit) {
		kind_ = port_file | prop;
		rep.stdio.file_ = file;	
		rep.stdio.closeit_ = closeit;
	}
	//	Makes the port a string port, and sets up the attributes.
	void set(int prop, char* start, char* past_the_end, char* curr) {
		kind_ = prop | port_string;
		rep.string.start_ = start;
		rep.string.past_the_end_ = past_the_end;
		rep.string.curr_ = curr;
	}

	///	Returns the FILE* of a file port.
	/**	This should only be used if is_file() is true. */
	FILE* file() const { return rep.stdio.file_; }

	///	Read a character from the port.
	/**	This should be an input port. */
	int inchar();
	///	Put back the character to a port.
	/**	This should be an input port. */
	void backchar(int c);
	///	Copy the given (null terminated) string to a port.
	void putstr(const char* s);
	///	Copy the given string to a port.
	/**	No special handling of null bytes is performed. */
	void putstr(const char* s, int len);
	///	Copy the given character to a port.
	void putcharacter(int c);
	///	Close a port.  
	// not clear what the exact meaning is
	void closeit(int flag);
	///	Close a port.  
	// not clear what the exact meaning is
	void close(int flag);

	///	Return the next port in the list.
	port* link() const { return rep.list.link_; }

private:
	//	disable copy constructor and assignment
	port(const port& val);
	port& operator=(const port& rhs);

	bool is_io() { return is_kind(port_input|port_output); }
	void mask(int flag) { kind_ &= ~flag; }
	void make_free() { kind_ = port_free; }
	bool closeit() { return rep.stdio.closeit_; }

	const char* start() const { return rep.string.start_; }
	const char* past_the_end() const { return rep.string.past_the_end_; }
	const char* curr() const { return rep.string.curr_; }
	const char* next() { return rep.string.curr_++; }
	const char* prev() { return --rep.string.curr_; }
	void append(char s) { *rep.string.curr_++ = s; }

	/////////////////
	//  variables  //
	/////////////////
	unsigned char kind_;
	union {
		struct {
		  FILE *file_;
		  bool closeit_;
		} stdio;
		struct {
		  char *start_;
		  char *past_the_end_;
		  char *curr_;
		} string;
		struct {
			port* link_;
		} list;
	} rep;
	static port* free_list_;
};

//-//////////
//-/  Op  ///
//-//////////

class op_code_table;

///	Op stores an interpreter operation code.
/**	Interpreter operations are primitives and procs.
 *	Procs correspond to predefined procedures.
 *	Primitives correspond to syntax and utiliy actions.
 *	<br>
 *	Op objects are stored in Proc cells.
 */
class Op {
public:

#define PROC_ENTRY0(op, func, name, minargs, maxargs) op,
#define PROC_ENTRY1(op, func, name, minargs, maxargs, type1) op,
#define PROC_ENTRY2(op, func, name, minargs, maxargs, type1, type2) op,
#define PROC_ENTRY3(op, func, name, minargs, maxargs, type1, type2, type3) op,
#define PRIM_ENTRY(op, func, name) op,
	///	This gives all the interpreter operations.
	/**	These are for op type 0. */
	enum OP {
	//	Get all the enum values from the op table.
	//	The enum constant names are in the first column.
#include "op.h"

	};
#undef PROC_ENTRY
#undef PROC_ENTRY0
#undef PROC_ENTRY1
#undef PROC_ENTRY2
#undef PROC_ENTRY3
#undef PRIM_ENTRY

	///	Default constructor.
	Op() { rep.l.oper_ = ILLEGAL; }
	///	Construct an Op from an OP enum value
	Op(short tbl, short op) {
		rep.ss.tbl_ = tbl;
		rep.ss.op_ = op;
#if CHECK
// cch		assert(tbl >= 0 && tbl < n_tables_);
		assert(op >= 0 && op <= ILLEGAL);
#endif
	}
	Op(long lop) {
		rep.l.oper_ = lop;
#if CHECK
// cch		assert(rep.ss.tbl_ >= 0 && rep.ss.tbl_ < n_tables_);
// cch		assert(rep.ss.op_ >= 0 && rep.ss.op_ <= ILLEGAL);
#endif
	}
	//	The default assignment and copy constructors work out OK.

	///	Return an OP value of the Op.
	long as_long() { return rep.l.oper_; }
	unsigned short tbl() const { return rep.ss.tbl_; }
	short op() const { return rep.ss.op_; } 
private:
	//	This depends on the way the compiler stores the union
	//	  elements.  The op_ MUST BE stored in the lower part of
	//	  oper_ and tbl_ MUST BE stored in the upper part of oper_.
	//	Check makes sure this works as expected.
	static bool check();
	union {
		struct {
			OP op_;
		};
		struct {
			short op_;
			unsigned short tbl_;
		} ss;
		struct {
			long oper_;
		} l;
	} rep;
	friend scheme;			// used for eval_cycle shortcut
	friend op_code_table;	// used for eval_cycle shortcut
};

//	Allows two Op values to be compared for equality.
inline bool operator==(Op lhs, Op rhs){
	return lhs.as_long() == rhs.as_long();
}

//-//////////////////
//-  op_code_entry //
//-//////////////////


//	These strings are concatenated by C's constant string concatenation rules
//	  into a single string in op_code_entry.
//	They are used as indexes into the array of test functions.
//	The character values are used as an array index, so do not change them.
//
//	TST_NONE does no tests, and is slightly faster than the others.
//	TST_ANY accepts anyhing, which has the same effect, 
//     but it is needed when there are other arguments.
#define TST_NONE 0
#define TST_ANY 1
#define TST_STRING 2
#define TST_SYMBOL 3
#define TST_PORT 4
#define TST_INPORT 5
#define TST_OUTPORT 6
#define TST_ENVIRONMENT 7
#define TST_PAIR 8
#define TST_LIST 9
#define TST_CHAR 10
#define TST_VECTOR 11
#define TST_NUMBER 12
#define TST_INTEGER 13
#define TST_NATURAL 14

///	Used to specify that any number of args are permissible.
enum { 
	INF_ARG = 10000 
};

///	An op_code_table is made up of op_code_entry values.
/**	Each op_code_entry describes a scheme op code.
 *	Included in it is its op value, its name, argument checking
 *	  information, and a functor used to call the operation.
 */
class op_code_entry {
public:

	typedef Callback_1_Void<scheme&>* op_callback;

	enum entry_type {
		proc = true, 
		prim = false
	};

	///	A dispatch_func implements the action of the operation.
	/**	Dispatch functions are members of the scheme class. */
	op_code_entry(): op_(Op::ILLEGAL), func_(0) {}
	~op_code_entry() { /*delete func_;*/ }
	///	Construct an entry from the components.
	op_code_entry(
		Op op,					///<	the op code being defined
		entry_type type,		///<	true if proc, false if prim				
		op_callback func,		///<	the function to execute	
		const std::string& name,		///<	the name of the proc or prim		
		int min_arity = 0,		///<	minimum number of args		
		int max_arity = 0,		///<	maximum number of args
		char test1 = 0,
		char test2 = 0,
		char test3 = 0):
			op_(op),
			type_(type),
			func_(func),
			name_(name),
			min_arity_(min_arity),
			max_arity_(max_arity) {
				arg_tests_encoding_[0] = test1;
				arg_tests_encoding_[1] = test2;
				arg_tests_encoding_[2] = test3;
				arg_tests_encoding_[3] = 0;
			}

	//	The default copy constructor is OK.
	//	the default assignment operator is OK.
	//	They are OK because every variable can be copied by value.

	///	Call the function for this entry.
	void call(scheme& sc) const { (*func_)(sc);}
	///	Returns the proc or primitive name.
	const std::string& name() const { return name_; }
	///	Test the argument list to see if it is legal.
	/**	Returns true if legal. */
	bool test_args(const scheme& sc, pointer args, char* errbuf) const;
	///	Return true if it is a proc, false if a prim.
	bool is_proc() const { return type_ == proc; }
	///	Return the op part of the op code.
	short op() const { return op_.op(); }

private:
	//	test if the arglist is legal
	//	this part tests the argument types
	//	returns true if there is a problem
	bool test_arg_types(const scheme& sc, pointer arglist, char* errbuf) const;

	//	variables
	Op op_;						//  the op code being defined
	entry_type type_;			//	true if this is a proc that belongs in the global environment
	op_callback func_;			//	pointer to function implementing op_code
	std::string name_;			//	the name of the operation (used for tracing and global environment)
	short min_arity_;			//	min number of args
	short max_arity_;			//	max number of args
	char arg_tests_encoding_[4]; // string of chars describing arg types

	friend scheme;
	friend op_code_table;
};

///	An op code description, used to build op_code_table.
/**	These values can be used to build an op_code_entry, which make
 *	  up the op_code_table.  But these values do not themselves appear in
 *	  the op_code_table.  They are missing two essential pieces of information.<br>
 *	(1) The op is missing its tbl part,<br>
 *	(2) The func is missing its instance pointer.
 *	<br><br>
 *	The op code description specifies its<br>
 *	Op value<br>
 *	entry type (proc or prim)<br>
 *	dispatch function (relative to some class T)<br>
 *	the operation name<br>
 *	the minimum number of arguments<br>
 *	the maximum number of arguments<br>
 *	the type requirements for each argument.<br>
 *	<br>
 *	Arrays of these entries can be statically defined
 *	and used to populate a scheme op_code_table.
 *	These op code entries are NOT associated with a particular
 *	instance, but are bound to an instance when the entries
 *	are incorporated into the op code table.<br>
 *	The type parameter T specifies the class that the callback
 *	  function calls into.  It is also used to build 
 *	  op_code_entry values for use in the op_code_table.
 */
template <class T>
class op_code_ent {
public:

	///	A dispatch_func implements the action of the operation.
	/**	Dispatch functions are members of the scheme class. */
	typedef void (T::*dispatch_func)(scheme&);
	///	Create a default entry.
	/**	Needed to instantiate in arrays. */
	op_code_ent(): op_(Op::ILLEGAL), func_(0) {}
	///	Construct an entry from the components.
	op_code_ent(
		Op op,					///<	the op code being defined
		op_code_entry::entry_type type,		///<	true if proc, false if prim				
		dispatch_func func,		///<	the function to execute	
		const std::string& name,		///<	the name of the proc or prim		
		int min_arity = 0,		///<	minimum number of args		
		int max_arity = 0,		///<	maximum number of args
		char test1 = 0,
		char test2 = 0,
		char test3 = 0) :
			op_(op),
			type_(type),
			func_(func),
			name_(name),
			min_arity_(min_arity),
			max_arity_(max_arity) {
				arg_tests_encoding_[0] = test1;
				arg_tests_encoding_[1] = test2;
				arg_tests_encoding_[2] = test3;
				arg_tests_encoding_[3] = 0;
			}

	//	The default copy constructor is OK.
	//	the default assignment operator is OK.
	//	These are both required for building and manipulating op code tables.
	//	They are OK because every variable can be copied by value.

	///	Build and return (by value) an op_code_entry.
	/**	Use the components of the op_code_ent and bind the
	 *	  instance sc to it.
	 *	This ensures that the sc pointer is compatible with the
	 *	  class that contains the function func.
	 */
	op_code_entry entry(T* sc) {
		op_code_entry ent(op_, 
			type_, 
			newCallback((op_code_entry::op_callback)0, sc, func_),
			name_,
			min_arity_,
			max_arity_,
			arg_tests_encoding_[0],
			arg_tests_encoding_[1],
			arg_tests_encoding_[2]
			);
		return ent;
	}

	///	Return the op part of the op code.
	short op() const { return op_.op(); }

private:
	//	variables
	Op op_;					//  the op code being defined
	op_code_entry::entry_type type_;		//	true if this is a proc that belongs in the global environment
	dispatch_func func_;	//	pointer to function implementing op_code
	std::string name_;		//	the name of the operation (used for tracing and global environment)
	short min_arity_;		//	min number of args
	short max_arity_;		//	max number of args
	char arg_tests_encoding_[4];
};

//-//////////////////
//-  op_code_table //
//-//////////////////

///	This table is used to test and dispatch the operations.
/**	The table is made up of an array of op_code_entry arrays.
 *	The two components of op supply the two subscripts to this 
 *	  double-indexed array.
 *	Each member table consists of op_code_entry values.
 *	Each entry is at a position corresponding to its own [op.tbl()][op.op()].
 *	The table contains a functor consisting of a scheme instance pointer
 *	  and a member function pointer, with which it calls the code
 *	  implementing the operation.
 */
class op_code_table {
public:
	///	Create an op code table.
	op_code_table();
	///	Destroy op code table.
	~op_code_table();

	///	Associate an instance of the scheme interpreter with this table.
	/**	This must be done first, before the table can be used. */
	void init(scheme* sc) { sc_ = sc; }

	///	Add a set of op code entries to the table.
	/**	The op-code_entry array base must be terminated by
	 *	  an entry with op Op::ILLEGAL.
	 *	Returns the op code table number that was used.
	 */
	int extend(op_code_entry* base, scheme* sc);

	template <typename T>
	static int get_size(T* base) {
		int highest = 0;
		for (int i = 0; base[i].op() != Op::ILLEGAL; i++) {
			int dest = base[i].op();
			if (dest > highest) highest = dest;
			if (i > highest) highest = i;
		}
		return highest+1;
	}

	template <typename T>
	void add_entries(op_code_ent<T>* base, int n, T* sc) {
		//	add fixed entries
        int i;
		for (i = 0; base[i].op() != Op::ILLEGAL; i++) {
			int dest = base[i].op();
			//	bind an instance with a generic entry and store in dispatch table.
			if (dest >= 0) {
#if CHECK
				assert(dispatch_table_[n][dest].op() == Op::ILLEGAL);
#endif
				dispatch_table_[n][dest] = base[i].entry(sc);
			}
		}
		//	add other entries
		int next_free = 0;
		for (i = 0; base[i].op() != Op::ILLEGAL; i++) {
			int dest = base[i].op();
			//	bind an instance with a generic entry and store in dispatch table.
			if (dest < 0) {
				while (dispatch_table_[n][next_free].op() != Op::ILLEGAL) next_free++;
				dispatch_table_[n][next_free] = base[i].entry(sc);
			}
		}
	}

	///	Add a set of base op code entries to the table.
	/**	The op-code_entry_base array base must be terminated by
	 *	  an entry with op Op::ILLEGAL.
	 *	This differs from the other version of extend in that
	 *	  entries in op_code_ent<scheme> are bound to the existing
	 *	  instance of scheme as indicated by init().
	 *	T has to support:<br>
	 *	  op() - return a table index for the operation<br>
	 *	  entry(T*) - return a op_code_entry bound to the T instance.
	 */
	template <typename T>
	int extend(op_code_ent<T>* base, T* sc) {
		int size = get_size(base);	// uses only op()
		//	extend the table
		op_code_entry* ent = new op_code_entry[size];
		int n = dispatch_table_.size();
		dispatch_table_.push_back(ent);
		size_.push_back(size);
		//	add new table entries
		add_entries(base, n, sc);
		init_proc(n, size);
		return n;
	}


	///	Return the size of the given table.
	int size(int n_table) const { return size_[n_table]; }

	///	Find an entry in the table, given its Op code.
	/**	Since the table is indexed by Op, this is O(1). */
	const op_code_entry* lookup_op(Op op) {
#if CHECK
		assert(op.rep.ss.tbl_ >= 0 && op.rep.ss.tbl_ < dispatch_table_.size());
		assert(op.rep.ss.op_ >= 0 && op.rep.ss.op_ <= Op::ILLEGAL);
#endif
		return &dispatch_table_[op.rep.ss.tbl_][op.rep.ss.op_];
	}
	///	Return the Op name.
	const std::string& lookup_name(Op op) {
#if CHECK
		assert(op.rep.ss.tbl_ >= 0 && op.rep.ss.tbl_ < dispatch_table_.size());
		assert(op.rep.ss.op_ >= 0 && op.rep.ss.op_ <= Op::ILLEGAL);
#endif
		return dispatch_table_[op.rep.ss.tbl_][op.rep.ss.op_].name();
	}
	///	Call the function implementing the given Op.
	void call(Op op) {
#if CHECK
		assert(op.rep.ss.tbl_ >= 0 && op.rep.ss.tbl_ < dispatch_table_.size());
		assert(op.rep.ss.op_ >= 0 && op.rep.ss.op_ <= Op::ILLEGAL);
#endif
		dispatch_table_[op.rep.ss.tbl_][op.rep.ss.op_].call(*sc_);
	}
	///	Return true if the Op is a proc.
	bool is_proc(Op op) {
#if CHECK
		assert(op.rep.ss.tbl_ >= 0 && op.rep.ss.tbl_ < dispatch_table_.size());
		assert(op.rep.ss.op_ >= 0 && op.rep.ss.op_ <= Op::ILLEGAL);
#endif
		return dispatch_table_[op.rep.ss.tbl_][op.rep.ss.op_].is_proc();
	}
	///	Put the proc names in the global environment.
	/**	At initialization, the outer environment needs to associate
	 *	the name of each proc with its op code.
	 *	This goes through the table, and for each proc, it makes a symbol
	 *	  cell from the name and a proc cell from the Op, and puts the pair
	 *	  into the environment.
	 *	Prim entries are skipped.
	 */
	void op_code_table::init_proc(int n_table, int size);
private:
	//	disable copy constructor and assignment
	op_code_table(const op_code_table& val);
	op_code_table& operator=(const op_code_table& rhs);

	//	This holds up to 10 dispatch tables for now.
	std::vector<op_code_entry*> dispatch_table_;
	std::vector<int> size_;
	scheme* sc_;
	friend scheme;
};

//-//////////////////
//-/  dump_stack  ///
//-//////////////////

class storage;
///	A dump stack frame holds the current operation,
///	  the args, the code, and the environment.
/**	Both fast and safe dump stacks use the same stack_frame.
 *	When a frame is in the fast stack, it is linked into a stack through
 *	  the link_ field.  Popped items are immediately transferred to the
 *	  free list.
 *	When it is in the slow stack, the frame is pointed to by a frame
 *	  cell, and the scheme stack takes care of things
 *	  so the link_ is 0.  Popped items are not transferred to the free list
 *	  because they could be required later.  Instead, frames are garbage
 *	  collected along with the frame cells that contain them.
 */
class stack_frame {
public:
	///	Allocate a stack frame and initialize it with the given registers.
	/**	A stack frame from the free list is used, if available.
	 *	Otherwise a new one is allocated.
	 */
	static stack_frame* make_frame(Op op, pointer args, pointer envir, pointer code,
		stack_frame* curr = 0);
	///	Return this stack frame to the free list.
	void free();
	stack_frame(Op op, pointer args, pointer envir, pointer code,
		stack_frame* curr) :
	  op_(op),
	  args_(args),
	  envir_(envir),
	  code_(code),
	  link_(curr) {}
	///	Restore the op, args, environment, and code from the frame.
	void get(Op& op, pointer& args, pointer& envir, pointer& code) {
		op = op_;
		args = args_;
		envir = envir_;
		code = code_;
	}
	///	Delete all stack frames on the free list at shutdown time.
	static void cleanup();
	///	Get the next frame in the free list.
	stack_frame* link() const { return link_; }

	///	Mark the cells pointed to by the args, environment, and code.
	void mark(storage& st);

private:
	//	disable copy constructor and assignment
	stack_frame(const stack_frame& val);
	stack_frame& operator=(const stack_frame& rhs);

	// this structure holds all the interpreter's registers
	Op op_; 
	pointer args_; 
	pointer envir_; 
	pointer code_; 
	//	Used to link to the next frame on the stack
	//	Also used to maintain the free list.
	stack_frame* link_;

	//	When a frame is freed, it is kept in a free list rather than
	//	  being returned to the heap.  This is the head of this list.
	//	the link()  method is used to access this list.
	static stack_frame* free_list_;
}; 

///	The fast version of the dump stack.
/**	The fast version is used until a continuation is stored, at
 *	  which time it switches to the safe version.  The safe
 *	  version continues to be used until it is certain that 
 *	  there is no possibility that it will be needed anymore.
 *	<br>
 *	The fast version uses a stack implemented as an array.
 *	Stack frames are allocated as needed.
 *	If it needs to be expanded, the array is reallocated.
 */
class dump_stack_fast {
public:
	dump_stack_fast() {
		dump_ = 0;
	}
	~dump_stack_fast();

	///	Returns a pointer to the top frame.
	/**	If the stack is empty, return 0 */
	stack_frame* top() const { return dump_; }
	///	Returns the next lower frame on the stack
	stack_frame* next(stack_frame* curr) const { return curr->link(); }
	///	Reset the dump stack.
	void reset();
	///	Mark all the pointers accessible from every frame.
	void mark(storage& st);
	///	Store the op, args, environment, and code in a new frame.
	/**	Allocate space for a new one if necessary. */
	void push_frame(Op op, pointer args, pointer envir, pointer code);
	///	Restore the op, args, environment, and code from the frame, and discard it.
	/**	This differs from dump_stack-safe because here the discarded stack frame
	 *	may be reused, so that a continuation cannot expect a saved
	 *	pointer to a dump stack frame to remain valid.  That is why there is no
	 *	curr() or set_curr() operation.
	 */
	bool pop_frame(Op& op, pointer& args, pointer& envir, pointer& code);

private:
	//	disable copy constructor and assignment
	dump_stack_fast(const dump_stack_fast& val);
	dump_stack_fast& operator=(const dump_stack_fast& rhs);

	stack_frame* dump_;	//	top of frame stack
};

///	The safe version of the dump stack.
/**	The safe version uses a stack implemented as a list.
 *	Here a popped frame is still valid.
 */
class dump_stack_safe {
public:
	dump_stack_safe(storage* store) : store_(store) { 
#if CHECK
		assert(store != 0);
#endif
		reset(); 
	}
	~dump_stack_safe() {}
	///	Copy the entire dump_stack_fast into a dump_stack_safe.
	/**	specific to dump_stack_safe. */
	void copy(dump_stack_fast& orig);
	///	Return true if the dump_stack is reset.
	/**	specific to dump_stack_safe. */
	bool is_reset() const;
	///	Reset the stack.
	void reset();
	///	Mark all the cells accesible from any frame in the stack.
	void mark(storage& st);
	///	Return a pointer to the current dump stack.
	/**	In this imlementation, the dump points to a FRAME cell. */
	pointer curr() const { return dump_; }
	///	Restore a saved dump stack pointer (made with curr()).
	void set_curr(pointer dump) { dump_ = dump; }
	///	Store the op, args, environment, and code in a new frame.
	void push_frame(Op op, pointer args, pointer envir, pointer code);
	///	Restore the op, args, environment, and code from the frame, and discard it.
	/**	Because the dump stack is implemented as a list, the popped frame
	 *	is still valid, in case some pointer still refers to it.
	 *	When no references remain, it will be garbage collected.
	 */
	bool pop_frame(Op& op, pointer& args, pointer& envir, pointer& code);
private:
	//	disable copy constructor and assignment
	dump_stack_safe(const dump_stack_safe& val);
	dump_stack_safe& operator=(const dump_stack_safe& rhs);

	pointer dump_;				//	list head for dump list

	storage* store_;			// used to allocate more storage
};

///	The dump stack represents the interpreter state.
/**	The dump stack has a fast version (that does not require
 *	  any cons cells) and a safe version (that works properly
 *	  with call/cc).
 *	The dump_stack uses the fast dump_stack implementation
 *	  as long as it can, and then switches to the safe
 *	  version.  
 *	When reset, it reverts to the fast version.
 *	It would be possible, after garbage collection, to set the
 *	  stack to fast if no continuatiuon cells were in use.
 *	That would require a copy operation in dump_stack_fast.
 */
class dump_stack {
private:
	enum stack_type { fast, safe };
public:
	///	Construct a dump stack.
	/**	Initially it is a fast stack */
	dump_stack(storage* store) : 
		type_(fast),
		dump_fast_(),
		dump_safe_(store) {}
	///	Reset the stack.
	/**	It is a good time to switch back to fast, if it has been
	 *	changed to safe.  There is still the chance that there is a
	 *	stored stack pointer, and if that needs to be restored, the stack
	 *	will go back to safe.
	 */
	void reset() {
		//	Switch to fast now.
		if (type_ == safe) {
			type_ = fast;
			dump_safe_.reset();
		}
		dump_fast_.reset();
	}
	///	Mark all cells accessible from the stack.
	void mark(storage& st) {
		//	Mark whichever one is in use now.
		if (type_ == fast) dump_fast_.mark(st);
		else dump_safe_.mark(st);
	}
	///	Save a pointer to the stack.
	/**	This requires the safe version.
	 *	If the safe version is not currently in use, one must be
	 *	created.
	 */
	pointer curr() {
		if (type_ == fast) {
			//	Make a safe dump and return it.
			dump_safe_.copy(dump_fast_);
		}
		return dump_safe_.curr();
	};
	///	Restore a saved pointer.
	/**	If the fast version is in use, switch to the safe one.
	 *	The safe version is the only one that can be set.
	 */
	void set_curr(pointer dump) {
		if (type_ == fast) {
			type_ = safe;
		}
		dump_safe_.set_curr(dump);
	};
	///	Push op, args, environment, and code into whichever one is currently in use.
	void push_frame(Op op, pointer args, pointer envir, pointer code) {
		if (type_ == fast) 
			dump_fast_.push_frame(op, args, envir, code);
		else
			dump_safe_.push_frame(op, args, envir, code);
	}
	///	Restore op, args, environment, and code from whichever one is in use.
	/**	If the safe stack is being used, and this makes it empty, change back
	 *	to the fast stack.
	 */
	bool pop_frame(Op& op, pointer& args, pointer& envir, pointer& code) {
		if (type_ == fast) 
			return dump_fast_.pop_frame(op, args, envir, code);
		else {
			bool ret = dump_safe_.pop_frame(op, args, envir, code);
			//	If the stack is now empty, revert to the fast version.
			if (dump_safe_.is_reset()) {
				type_ = fast;
				dump_fast_.reset();
			}
			return ret;
		}
	}
private:
	//	disable copy constructor and assignment
	dump_stack(const dump_stack& val);
	dump_stack& operator=(const dump_stack& rhs);

	stack_type type_;
	dump_stack_fast dump_fast_;
	dump_stack_safe dump_safe_;
};

//-////////////
//-/  cell  ///
//-////////////

typedef std::vector<pointer> scheme_vector;
typedef std::string scheme_string;

///	Scheme cell structure.
/**	A cell stores an integer, a real, a pair, etc.
 *	The type() tells what specific information is stored in a given cell.
 *	Various accessor methods get the data out, depending on the type.
 */
class cell {
public:			
	//	Note: these could be private, since no one uses them.
	//	They are public only so daxygen will include them.

	///	Values of the cell type.
	/**	These values are stored in the type_ field.
	 *	Upper case letters imply that there is something
	 *	  in the car and cdr fields.
	 *	Lower case letters store things in other union members.
	 *	The program does NOT depend on this -- it uses F_ATOM instead.
	 *	These all have to be unique, otherwise they can be anything.
	 */
	enum cell_type {
	  T_NONE =			'N',	///<	nothing stored
	  T_ATOM =			'A',	///<	special atom: NIL, T, F
	  T_STRING =		's',	///<	stores a string (char*, len)
	  T_INTEGER =		'i',	///<	stores an integer
	  T_REAL =			'r',	///<	stores a real
	  T_SYMBOL =		'Y',	///<	stores a symbol (string, alist)
	  T_PROC =			'p',	///<	stores a proc (Op)
	  T_PAIR =			'P',	///<	stores a pair (car, cdr)
	  T_CLOSURE =		'C',	///<	stores a closure (code, envir)
	  T_CONTINUATION =	'T',	///<	stores a continuation (dump)
	  T_FOREIGN =		'f',	///<	stores a foreign function pointer
	  T_CHARACTER =		'c',	///<	stores a character
	  T_PORT =			'q',	///<	stores a port (port*)
	  T_FRAME =			'f',	///<	stores a dump frame
	  T_VECTOR =		'v',	///<	stores a vector (len, vector*)
	  T_MACRO =			'm',	///<	stores a macro
	  T_PROMISE =		'Z',	///<	stores a promise (code, env)
	  T_ENVIRONMENT =	'E',	///<	stores an environment (new_env, old-env)
	};
	///	These flag values give additional information about the cell.
	/**	They are stored in the flag_ field.
	 */
	enum cell_flags{
		F_SYNTAX =     0x01,    ///<  0000 0001	 special modifier for T_SYMBOL
		F_IMMUTABLE =  0x02,    ///<  0000 0010  can modify many types
		F_ATOM =       0x04,    ///<  0000 0100  if NOT set, car and cdr point to cells
		F_MARKFUN =    0x08,	///<  0000 1000	 if set, gc calls special mark fun
		F_FINALFUN =   0x10,	///<  0001 0000  if set, gc calls special final fun
		F_MARK =       0x80,    ///<  1000 0000	 used by mark
		M_UNMARK =     0x7f,    ///<  0111 1111	 used by mark
		F_ALINK =	   0x40,	///<  0100 0000	 used by mark
		M_CLRALINK =   0xbf		///<  1011 1111	 used by mark
	};

public:

	///	create a cell
	/**	This is only used to create static cells Nil, T, F, etc. */
	cell(): type_(T_NONE), flag_(0) {}

	///	Assign from one cell to another: type, flag and contents.
	cell& operator=(const cell& rhs) {
#if CHECK
		assert(type_ != T_NONE);
#endif
		type_ = rhs.type_;
		flag_ = rhs.flag_; 
		object = rhs.object; 
		return *this;
	}

	char type() const { return type_; }
	bool is_type(char t) const { return type_ == t; }

	///	Returns true if the cell type is NONE.
	/**	Free cells are cleared and linked onto the free_cell list. */
	bool is_clear() const { return type_ == T_NONE; }
	///	Clear a cell.
	/**	Used when shutting down. */
	void clear() { 
		type_ = T_NONE; 
		flag_ = 0;		// not strictly necessary, but a good idea
	}
	///	Clear a cell, but set is car and cdr fields.
	/**	On free list, free cells are linked through cdr. */
	void clear(pointer car, pointer cdr) {
		clear();
		set_cell(car, cdr);
	}
	///	Return the cdr of a pair.
	pointer free_cdr() const { 
#if CHECK
		assert(is_clear());
#endif
		return object.cons.cdr_; 
	}
	///	Set the cdr of a free cell.
	pointer set_free_cdr(pointer v) { 
#if CHECK
		assert(is_clear());
#endif
		return object.cons.cdr_ = v; 
	}

	bool mark_fun() const {
		return (flag_ & F_MARKFUN) != 0;
	}

	bool final_fun() const {
		return (flag_ & F_FINALFUN) != 0;
	}

	///	Returns true if the cell holds a number (integer or real).
	/**	Numbers are stored in the ivalue of rvalue fields. */
	bool is_number() const { return type_ == T_INTEGER || type_ == T_REAL; }

	///	Returns true if the cell holds an integer.
	/**	Integers are stored in the ivalue field. */
	bool is_integer() const {return type_ == T_INTEGER; }
	bool is_nonneg_integer() const { return type_ == T_INTEGER && object.ivalue_ >= 0; }
	///	Give the value of a number cell as a long.
	long ivalue()  const { 
		if (is_integer()) return object.ivalue_;
		if (is_real()) return (long)object.rvalue_;
#if CHECK
		assert(false);
#endif
		return object.ivalue_;
	}
	///	Store an integer in a cell.
	void set_integer(long ivalue) {
		type_ = T_INTEGER;
		flag_ = F_ATOM; 
		object.ivalue_ = ivalue; 
	}

	///	Return true if the cell stores a real.
	/**	Reals are stored in the rvalue field. */
	bool is_real() const { return type_ == T_REAL; }
	///	Return the number as a real value.
	double rvalue() const { 
		if (is_integer()) return (double)object.ivalue_;
		if (is_real()) return object.rvalue_;
#if CHECK
		assert(false);
#endif
		return object.rvalue_;
	}
	///	Store a real in the cell.
	void set_real(double rvalue) { 
		type_ = T_REAL;
		flag_ = F_ATOM; 
		object.rvalue_ = rvalue; 
	}

	///	Returns true if the cell stores a proc (an interpreter primitive or proc).
	bool is_proc() const    { return type_ == T_PROC; }
	///	Return the Op in the cell.
	/**	Proc stores an Op::op in the procvalue field.
	 *	It does not simply store Op because then it could not have a ctor.
	 */
	Op procvalue() const { 
#if CHECK
		assert(is_proc());
#endif
		return Op(object.procvalue_); 
	}
	///	Store an op in the cell.
	void set_proc(Op op) { 
		type_ = T_PROC;
		flag_ = F_ATOM; 
		object.procvalue_ = op.as_long(); 
	}

	///	Return true if the cell is a pair (car and cdr).
	/** A pair stores pointers in the car and cdr fields. */
	bool is_pair() const     { return type_ == T_PAIR; }
	///	Return the car of the pair.
	pointer car() const { 
#if CHECK
		assert(!flag_atom());
#endif
		return object.cons.car_; 
	}
	///	Return the cdr of a pair.
	pointer cdr() const { 
#if CHECK
		assert(!flag_atom() || is_clear());
#endif
		return object.cons.cdr_; 
	}
	///	Store the car and cdr as a pair.
	/**	If immutable is true, the set the immutable flag. */
	void set_pair(pointer car, pointer cdr, bool immutable = false) { 
		type_ = T_PAIR; 
		flag_ = (immutable) ? F_IMMUTABLE  : 0;
		set_cell(car, cdr);
	}
	///	Set the car of a cell.
	pointer set_car(pointer v) { 
#if CHECK
		assert(!flag_atom());
#endif
		return object.cons.car_ = v; 
	}
	///	Set the cdr of a cell.
	pointer set_cdr(pointer v) { 
#if CHECK
		assert(!flag_atom());
#endif
		return object.cons.cdr_ = v; 
	}
	///	Set both the car and cdr.
	void set_cell(pointer car, pointer cdr) {
#if CHECK
		assert(!flag_atom());
#endif
		object.cons.car_ = car;
		object.cons.cdr_ = cdr;
	}

	///	Return true if the cell is a vector.
	/**	A vector stores the vector length and a pointer to 
	 *	  a dynamically-allocated array of pointers.
	 */
	bool is_vector() const { return type_ == T_VECTOR; }
	scheme_vector* make_vector(int len) {
		return new scheme_vector(len);
	}
	///	Makes a vector cell and stores its pointer array.
	void set_vector(int len) { 
		type_ = T_VECTOR;
		flag_ = F_ATOM | F_MARKFUN | F_FINALFUN; 
		object.vector.vec_ = make_vector(len);
	}
	void delete_vector() {
#if CHECK
		assert(is_vector());
#endif
		delete object.vector.vec_;
	}
	///	Return the length of the vector.
	int veclen() const { 
#if CHECK
		assert(is_vector());
#endif
		return object.vector.vec_->size(); 
	}
	///	Return the array of pointers.
	scheme_vector& vecvalue() const { 
#if CHECK
		assert(is_vector());
#endif
		return *object.vector.vec_; 
	}
	///	Fill each vector elem with the pointer.
	void fill_vector(pointer obj);
	///	Returns the vector element at the given index.
	pointer vector_elem(int ielem);
	///	Stores the pointer into the vector at the given index.
	pointer set_vector_elem(int ielem, pointer a);

	bool is_atom() const { return type_ == T_ATOM; }
	///	Sets a cell's atom flag, while clearing its type.
	/**	Used when shutting down */
	void set_atom() { type_ = 0; flag_ = F_ATOM; }
	///	Sets a cell's atom and mark flags, while clearing its type.
	void set_marked_atom(pointer car, pointer cdr) { 
		type_ = T_ATOM; 
		flag_ = F_ATOM | F_MARK; 
		object.cons.car_ = car;		// cannot use set_cell
		object.cons.cdr_ = cdr;
	}

	///	Return true if the cell holds a string.
	/**	A string is represented by a length (in string.length_) 
	 *	  and a std::string (in string.svalue_).
	 *	The actual character storeage is allocated and freed with the cell.
	 */
	bool is_string() const { return type_ == T_STRING; }
	///	Return the character string.
	/**	This should be null-terminated. */
	scheme_string& strvalue() const { 
#if CHECK
		assert(is_string());
#endif
		return *object.string.svalue_; 
	}
	///	Returns the length of the string (including the null).
	int strlength() const { 
#if CHECK
		assert(is_string());
#endif
		return object.string.svalue_->size(); 
	}
	///	Makes a cell a string cell, given the string.
	void set_string(const char* str, int len) { 
		type_ = T_STRING;
		flag_ = F_ATOM | F_FINALFUN;
		object.string.svalue_ = new scheme_string(str, len);	
	}
	void set_string(const std::string& str) { 
		type_ = T_STRING;
		flag_ = F_ATOM | F_FINALFUN;
		object.string.svalue_ = new scheme_string(str);	
	}
	void set_string(const char* str) {
		type_ = T_STRING;
		flag_ = F_ATOM | F_FINALFUN;
		object.string.svalue_ = new scheme_string(str);	
	}
	void set_string(int len, char fill) {
		type_ = T_STRING;
		flag_ = F_ATOM | F_FINALFUN;
		object.string.svalue_ = new scheme_string(len, fill);	
	}
	void set_string_elem(int pos, int ch) {
		object.string.svalue_->replace(pos, 1, 1, (char)ch);
	}
	void delete_string() {
#if CHECK
		assert(is_string());
#endif
		delete object.string.svalue_;
	}

	///	Return true if the cell is a character.
	/**	A character stores its value in cvalue. */
	bool is_character() const { return type_ == T_CHARACTER; }
	///	Return the character.
	char charvalue() const { 
#if CHECK
		assert(is_character());
#endif
		return (char)object.cvalue_; 
	}
	///	Store the given character in a character cell.
	void set_character(int c) { 
		type_ = T_CHARACTER;
		flag_ = F_ATOM; 
		object.cvalue_ = c;; 
	}


	///	Return true if this is a port cell.
	/**	A port cell stores a port*.
	 *	Ports themselves are allocated and freed along with the cells that 
	 *	  contain them.
	 */
	bool is_port() const { return type_ == T_PORT; }
	///	Return true if this is an input port cell.
	bool is_inport() const { return type_ == T_PORT && object.port_->is_input(); }
	///	Return true if this is an output port cell.
	bool is_outport() const { return type_ == T_PORT && object.port_->is_output(); }
	///	Return the port that this cell stores.
	port* portvalue() const { 
#if CHECK
		assert(is_port());
#endif
		return object.port_; 
	}
	///	Make this a port cell, storing the given port.
	void set_port(port* p) {
		type_ = T_PORT;
		flag_ = F_ATOM | F_FINALFUN; 
		object.port_ = p; 
	}

	///	Return true if this is a frame.
	/**	A frame cell stores a stack_frame*.
	 *	Dump_stack_frames themselves are allocated and freed along with the cells that 
	 *	  contain them.
	 */
	bool is_frame() const { return type_ == T_FRAME; }
	///	Return the frame stored in the cell.
	stack_frame* framevalue() const { 
#if CHECK
		assert(is_frame());
#endif
		return object.frame_; 
	}
	///	Make a frame cell, given the frame.
	void set_frame(stack_frame* f) {
		type_ = T_FRAME;
		flag_ = F_ATOM | F_MARKFUN | F_FINALFUN; 
		object.frame_ = f; 
	}

	///	Retrn true if this is a symbol cell.
	/**	A symbol is a pair who's car is an immutable string 
	 *	  and who's cdr is the symbol's property list.
	 *	Symbols are on the oblist.
	 */
	bool is_symbol() const { return type_ == T_SYMBOL; }
	///	Return the symbol name.
	const std::string& symname() const  { 
#if CHECK
		assert(is_symbol());
#endif
		return car()->strvalue(); 
	}
	///	Only symbols hve property lists.
	bool hasprop() const   { return type_ == T_SYMBOL; }
	///	Return a symbol's property list.
	pointer symprop() const { 
#if CHECK
		assert(is_symbol());
#endif
		return cdr(); 
	}
	///	Set a symbols property list.
	pointer set_symprop(pointer q) { 
#if CHECK
		assert(is_symbol());
#endif
		return set_cdr(q); 
	}
	///	Make a cell a symbol cell.
	/**	A is the string cell, b is property list.
	 *	Normally b would be nil to start with.
	 */
	void set_symbol(pointer str, pointer b) {
		type_ = T_SYMBOL;
		flag_ = F_IMMUTABLE; 
		set_cell(str, b); 
	}

	///	Return true if this is a syntax cell.
	/**	A syntax cell is a symbol cell with the flag F_SYNTAX.
	 *	Like symbols, it is on the oblist.
	 *	All syntax cells are immutable.
	 */
	bool is_syntax() const   { return !!(flag_ & F_SYNTAX); }
	///	Set the syntax flag.
	void set_syntax() { flag_ |= F_SYNTAX; }

	///	Return true if the cell is a foreign cell.
	/**	A foreign cell stores a function pointer in the ff field. */
	bool is_foreign() const  { return type_ == T_FOREIGN; }
	///	Call the foreign function who's pointer is stored in this cell.
	pointer call_foreign(scheme& sc, pointer args) { 
#if CHECK
		assert(is_foreign());
#endif
		return object.ff_(sc, args); 
	}
	///	Make the cell a foreign cell, with the supplied function.
	void set_foreign(foreign_func f) { 
		type_ = T_FOREIGN;
		flag_ = F_ATOM; 
		object.ff_ = f; 
	}

	///	Return true if the cell is a closure.
	/**	A closure cell stores the code in the car
	 *	  and the environment in the cdr.
	 */
	bool is_closure() const  { return type_ == T_CLOSURE; }
	///	Return the closure code.
	pointer closure_code() const   { 
#if CHECK
		assert(is_closure() || is_macro() || is_promise());
#endif
		return car(); 
	}
	///	Return the cosure environment.
	pointer closure_env() const    { 
#if CHECK
		assert(is_closure() || is_macro() || is_promise());
#endif
		return cdr(); 
	}
	///	Make this cell a closure cell.
	void set_closure(pointer code, pointer env) { 
		type_ = T_CLOSURE;
		flag_ = 0;
		set_cell(code, env); 
	}

	///	Return true if this is a macro sell.
	bool is_macro() const    { return type_ == T_MACRO; }
	///	Set to macro type.
	void set_macro() { type_ = T_MACRO; }

	///	Return true if this is a continuation cell.
	/**	A continuation cell stores the continuation dump in the car. */
	bool is_continuation() const { return type_ == T_CONTINUATION; }
	///	Return the continuation.
	pointer continuation() const { 
#if CHECK
		assert(is_continuation());
#endif
		return car(); 
	}
	///	Make this cell a continuation cell.
	void set_continuation(pointer d) { 
		type_ = T_CONTINUATION; 
		flag_ = 0;
		set_cell(d, 0); 
	}

	///	Return true if this cell is a promise.
	/**	A promise cell stores the code in the car
	 *	  and the environment in the cdr.
	 */
	bool is_promise() const  { return type_ == T_PROMISE; }
	///	Make this cell a promise.
	void set_promise(pointer code, pointer env) { 
		type_ = T_PROMISE; 
		flag_ = 0;
		set_cell(code, env); 
	}

	///	Return true if this cell is an environment cell.
	bool is_environment() const { return type_ == T_ENVIRONMENT; }
	///	Make this cell an environment cell.
	void set_environment(pointer n, pointer env) { 
		type_ = T_ENVIRONMENT;
		flag_ = F_IMMUTABLE; 
		set_cell(n, env); 
	}


	///	Returns true if the ALINK flag is on.
	/**	The mark algorithm needs to store a single bit in a cell
	 *	while it is marking the structure.  It can use the atom bit, which
	 *	it will restore later.  Since there are spare bits in the flag, we
	 *	have it use a separate flag bit for clarity.
	 */
	bool flag_alink() const { return !!(flag_ & F_ALINK); }
	///	Set the ALINK flag.
	void set_alink() { flag_ |= F_ALINK; }
	///	Clear the ALLINK flag.
	void clr_alink() { flag_ &= M_CLRALINK; }

	///	Return true if the ATOM flag is on.
	/**	If the ATOM flag is on, then the cell does not contain
	 *	cell pointers that have to be marked.
	 */
	bool flag_atom() const { return !!(flag_ & F_ATOM); }

	///	Return true if the cell is marked.
	/**	The mark algorithm follows all pointer links from the known cells,
	 *	marking each one.  If the marked cell is not an atom, it marks
	 *	everything that it points to also.  After marking everything that is
	 *	reachable, it garbage collects everything else.
	 */
	bool flag_mark() const { return !!(flag_ & F_MARK); }
	///	Set the mark.
	void setmark() { flag_ |= F_MARK; }
	///	Clear the mark.
	void clrmark() { flag_ &= M_UNMARK; }

	///	Return true if the cell is immutable.
	bool is_immutable() const { return !!(flag_ & F_IMMUTABLE); }
	///	Make the cell immutable.
	/**	Once the cell is marked immutable, there is no mechanism to 
	 *	make it mutable (except to overwrite it completely).
	 */
	void set_immutable() { flag_ |= F_IMMUTABLE; }

	pointer caar() const	 { return car()->car(); }
	pointer cadr() const	 { return cdr()->car(); }
	pointer cdar() const     { return car()->cdr(); }
	pointer cddr() const     { return cdr()->cdr(); }
	pointer cadar() const    { return car()->cdr()->car(); }
	pointer caddr() const    { return cdr()->cdr()->car(); }
	pointer cadaar() const   { return car()->car()->cdr()->car(); }
	pointer cadddr() const   { return cdr()->cdr()->cdr()->car(); }
	pointer cddddr() const   { return cdr()->cdr()->cdr()->cdr(); }

protected:
	//	disable copy constructor
	cell(const cell& val);

	//	By packing to a 2 byte boundary, the cell is 10 bytes.
#pragma pack(push)
#pragma pack(2)
  unsigned char type_;
  unsigned char flag_;
  union {
    struct {
      cell* car_;
      cell* cdr_;
    } cons;
	struct {
		scheme_string* svalue_;
	} string;
	struct {
		scheme_vector* vec_;
	} vector;
	long ivalue_;
	double rvalue_;
	long procvalue_;
	int cvalue_;
    port *port_;
	stack_frame* frame_;
    foreign_func ff_;
#if MEMBLK
	struct {
	  int*   mvalue_;
	  int     mlen_;
	} memblock;
#endif
  } object;
#pragma pack(pop)
};

class scheme;

//-///////////////
//-/  storage  ///
//-///////////////

///	Manages the storage (alloocation, deallocation, and garbage collection)
///	  of cells.
/**	This knows how to allocate and deallocate cells and the data structures such
 *	as strings, ports, and frames that are associated with them.
 *	<br>
 *	Cells are allocaed in segments, which are blocks of storage that hold only cells.
 *	Segments are allocated as needed.  
 *	The address of each segment is kept in an array, sorted in increasing address.
 *	Cells that are not in use are kept in a free list, which is kept sorted in address order.
 */
class storage {
public:
	///	Construct storage.
	storage(int max_cell_seg);
	///	Destruct storage.
	~storage();
	///	Initialize.
	/**	Returns true if it was impossible to allocate the
	 *	requested memory, false otherwise.
	 */
	bool init(scheme& sc, int initial_segs);

	void mark(pointer a);

	typedef bool (*test_predicate)(pointer p);
	typedef void (*marker)(storage&, pointer);
	class mark_entry {
	public:
		mark_entry() : test_(0), mark_(0) {}
		mark_entry(char test, marker mark) : test_(test), mark_(mark) {}
		char test_;
		marker mark_;
	};
	void extend_marker(mark_entry* mark_table);

	typedef void (*finalizer)(storage&, pointer);
	class final_entry {
	public:
		final_entry() : test_(0), finalize_(0) {}
		final_entry(char test, finalizer finalize) : 
			test_(test), 
			finalize_(finalize) {}
		char test_;
		finalizer finalize_;
	};
	void extend_finalizer(final_entry* final_table);

	/// Special cell for failed memory allocation.
	static const pointer sink_;			 

	///	Return true if the system has run out of memory.
	bool no_memory() const { return no_memory_; }
	///	Return true if gc prints out each time it runs.
	bool gc_verbose() const { return gc_verbose_; }
	///	Sets the gc verbose flag.
	void set_gc_verbose(bool verbose) { gc_verbose_ = verbose; }
	///	This is called to signal that there is no memory left.
	void set_no_memory() { no_memory_ = true; }
	///	Allocate additional memory to store cells.
	/**	The number of additional segments is given.
	 *	Returns the number of cell segments actually allocated. 
	 *	If the number iallocated was less than the number requested,
	 *	then a no_memory indication will be set.
	 */
	int alloc_cells(
		///	The number of additional segments needed.
		int n_seg = 1	
		);
	///	Get more cell storage.
	/**	If there are free cells, then this does nothing.
	 *	If not, it does a garbage collection.
	 *	If that does not make enough cells, it allocates more.
	 */
	void get_more_cells(pointer a, pointer b);
	///	Returns a free cell.
	/**	Gets a free cell from the free list.
	 *	If the free list is empty, returns 0.
	 */
	pointer get_free_cell();
	///	Returns a free cell.
	/**	Gets a cell from the free list.
	 *	If there is none in the free list, collect garbage.
	 *	If that does not give any free cells, allocate more,
	 *	<br>
	 *	The pointers passed as parameters are marked, along with all
	 *	the other cells accessible from the global environment.
	 */
	pointer find_cell(pointer a = 0, pointer b = 0);
	///	Perform garbage collection.
	/**	This is done in the following steps:
	 *	<br>
	 *	1.  Mark all accessible cells, including the ones passed in.
	 *	<br>
	 *	2.  Deallocate all other cells and put them on the free list.
	 */
	int gc(pointer a = 0, pointer b = 0);

	///	Create a pair.
	/**	If immutable is true, then the IMMUTABLE flag is set. */
	pointer cons(pointer a, pointer b, bool immutable = false);
	///	Create an immutable pair.
	pointer immutable_cons(pointer a, pointer b) { return cons(a, b, true); }
	///	Create a port cell, with the given port.
	pointer mk_port(port* p);					// cch
	pointer mk_port(FILE* f, int prop, bool closeit = false);
	pointer mk_port(const std::string& filename, int prop);
	pointer mk_port(char* start, char* past_the_end, int prop);
	///	Create a foreign function cell.
	pointer mk_foreign_func(foreign_func f);
	///	Create an environment cell.
	pointer mk_environment(pointer new_frame, pointer old_env);
	///	Create a character cell.
	pointer mk_character(int c);
	///	Create an integer cell.
	pointer mk_integer(long n);
	///	Create a real cell.
	pointer mk_real(double d);
	///	Create an integer or real cell.
	pointer mk_number(num n);
	///	Create a symbol cell.
	/**	Note that this does not enter it in the symbol table. */
	pointer mk_symbol(const std::string& name);
	///	Create a string cell.
	pointer mk_string(const char* str, int len);
	///	Create a string cell.
	pointer mk_string(const char* str);
	///	Create a string cell.
	pointer mk_string(const std::string& str);
	///	Create a string cell.
	pointer mk_string(int len, char fill);
	///	Create a frame cell.
	pointer mk_frame(Op op, pointer args, pointer envir, pointer code);	
	///	Create a vector cell
	pointer mk_vector(int len);
	///	Create a closure cell.
	pointer mk_closure(pointer c, pointer e);
	///	Create a promise cell.
	pointer mk_promise(pointer c, pointer e);
	///	Create a continuation cell.
	pointer mk_continuation(pointer d);
	///	Create a proc cell.
	pointer mk_proc(Op op);

	///	Returns the cell to the free list.
	/**	If it is a string, port, frame, or vector, return the
	 *	associated storage as well.
	 */
	void finalize_cell(pointer a);


private:
	int collect_garbage();		// collects all non-marked calls onto the free list

	void add_free_cells(pointer first, pointer last);
	void insert_new_seg(pointer newp);

	static cell cellSink_;		// used when out of storage, and a function must return a cell pointer

	//	cell allocation constants and variables
	enum { CELL_SEGSIZE  =  5000 };  // # of cells in one segment 
	pointer* cell_seg_;			//	Array of pointers to cell blocks, in address order
	int last_cell_seg_;			//	Index of last segment.
	int max_cell_seg_;			//	Total number of segments allowed.

	pointer free_cell_;       // pointer to top of free cells 
	long    fcells_;          // # of free cells 
	bool    no_memory_;       // if true, mem. alloc. has failed 
	bool    gc_verbose_;      // if true, print gc status 

	//	Some cell types need special mark functions.
	//	This maps the cell type to the function to call.
	typedef std::map<char, marker> mark_map_type;
	typedef std::pair<char, marker> mark_map_key;
	mark_map_type mark_map_;

	//	Some cell types need special finalize functions.
	//	This maps the cell type to the function to call.
	typedef std::map<char, finalizer> final_map_type;
	typedef std::pair<char, finalizer> final_map_key;
	final_map_type final_map_;

	scheme* sc_;			// this is used only to print gc messages
};

//-///////////////////
//-/  environment  ///
//-///////////////////

///	Manages the environment
class environment {
public:
	///	Construct an environment object.
	environment(storage* store) : store_(store), envir_(0) {
#if CHECK
		assert(store != 0);
#endif
	}

	///	Return the current environment.
	/**	This is suitable for storing and setting later. */
	pointer curr() const { return envir_; }
	///	Set the environment to a previous state.
	void set_curr(pointer envir) { envir_ = envir; }

	///	Mark all cells reachable from this environment.
	void mark() { store_->mark(envir_); }

	///	Push a new, empty, environment.
	void push(pointer old_env = 0);
	///	Look up and return the variable in the given environment.
	pointer get_binding(pointer variable, pointer envir) const;
	///	Look up and return the variable in the current environment.
	pointer get_binding(pointer variable) const;

	///	Store the variable in the given environment.
	pointer define_new_binding(pointer variable, pointer value, pointer env);
	///	Store the variable in the current environment.
	/**	The variable must not already be bound in the environmet */
	pointer define_new_binding(pointer variable, pointer value);

	//	Store the variable in the given environment.
	bool define_binding(pointer variable, pointer value, pointer envir, 
						bool top, bool add);
	///	Store the variable in the current environment.
	/**	If the variable is already bound, then update it.
	 *	Otherwise, (if add is true) add it to the current environment.
	 *	If top is true, use the top environment instead of the current.
	*/
	bool define_binding(pointer variable, pointer value,
						bool top, bool add);
private:
	// The interaction-environment has about 300 variables in it.
	enum { base_frame_size = 461 }; 

	static pointer base(pointer env, pointer variable);
	static pointer search_top_binding(pointer env, pointer variable);
	static pointer search_binding(pointer env, pointer variable);

	//	Return the initial environment.
	pointer new_frame(pointer old_env);

	//	variables
	pointer envir_;		// the current environment

	storage* store_;	// used to get storage cells
};

//-//////////////
//-/  oblist  ///
//-//////////////

///	Stores the symbol table.
/**	Every symbol that is encountered is stored here.
 *	No information is associated with the symbol, except its
 *	address in the store.  Any two symbols with the same name
 *	will result in the same pointer.
 *
 *	This implementation uses a vector as a hash table, with a list at each
 *	  position,  Only the list at the hash location is searched for the symbol
 */
class oblist {
public:
	oblist(storage* store) : store_(store) {}
	///	Initialize the symbol table.
	void init();
	///	Mark all cells accessible through the symbol table.
	void mark() { store_->mark(oblist_); }
	///	Reset the symbol table.
	void clear();

	///	Adds the given symbol to the symbol table.
	pointer add_symbol(
		///	The symbol name to add..
		const std::string& name);
	///	Finds the given symbol in the symbol table.
	/**	Returns a pointer to the symbol cell.
	 *	If not found, returns 0.
	 */
	pointer find_symbol(
		///	The symbol to look for.
		const std::string& name);
	///	Returns a list of all symbols in the symbol table.
	/**	In one implementationm, this returns the internal list used to
	 *	store the symbol table.  In the other, it must build this list
	 *	specifically to satisfy the request. 
	 */
	pointer all_symbols();

private:

	//	The size of the symbol table hash table
	enum {oblist_size = 461 };	// probably should be bigger

	//	Search for the name in the given list.
	static pointer search_symbol(const std::string& name, pointer p);

	//	Varibles
	pointer oblist_;	//	The symbol table.
	//	Used to get the storage for the symbol table.
	storage* store_;
};


//-//////////////
//-/          ///
//-/  scheme  ///
//-/          ///
//-//////////////
//	These execute s_actions and then return some value
//	They are used within opexe_ routines to re-dispatch another operation
#define s_error0(str)		{ s_error_action((str), 0); return; }
#define s_error1(str, a)	{ s_error_action((str), (a)); return; }
#define s_goto(a)			{ s_goto_action(a); return; }
#define s_return(a)			{ s_return_action(a); return; } 
#define s_retbool(tf)		{ s_return_action((tf) ? T_ : F_); return; }

///	This provides an instance of a self-contained scheme interpreter.
/**	Instantiate one of these to create a scheme interpreter.
 *	Different instances give entirely separate interpreters, which cannot
 *	interact with one another.
 */
class scheme {
public:
	///	Construct an interpreter.
	scheme(int first_cell_seg,		///<	Initial number of cell segments to allocate.	
		int max_cell_seg			///<	Maximum number of cell segments to allocate.
		);
	///	Destroy interpreter.
	virtual ~scheme();
	///	Returns true if the interpreter is valid.
	/**	If the interpreter runs out of memory or cannot function for
	 *	some other reason, this returns false.
	 */
	bool is_valid() { return !store_.no_memory(); }

	///	Set the input port to the given FILE*.
	void set_input_port(
		FILE *fin						///< The input file.
		);
	///	Set the input port to the given string.
	void set_input_port(
		char *start,					///<	The input string.
		char *past_the_end				///<	Points one past the string end.
		);
	///	Set the otput port to the given FILE*.
	void set_output_port(
		FILE *fin						//<	The output file.
		);		
	///	Set the output port to the given string.
	void set_output_port(
		char *start,					///<	The output string
		char *past_the_end				///<	Points one past the string end.
		);
	///	Load from the given FILE*.
	//*	Returns the retcode */
	int load(
		FILE *fin						///<	The file to load from.
		);
	///	Load from the given string.
	/**	Returns the retcode. */
	int load(
		char *cmd						///<	The string to load.
		);
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
	void define(
		pointer symbol,			///<	The symbol to deifine.
		pointer value,			///<	The symbol's value.
		pointer env = 0			///<	Te environment in which to make the definition.
		);

	///	Create a pair cell.
	/**	Obtain a free cell, make it a pair cell, and tores the given
	 *	pointers in its cons and cdr.
	 *	If immutable is true, then set the immutable flag on the cell.
	 *	Returns the cons cell.
	 */
	pointer cons(
		pointer a,				///<	The pair car.
		pointer b,				///<	The pair cdr.
		bool immutable = false	///<	If true, make immutable.
		) {
		return store_.cons(a, b, immutable); 
	}
;
	///	Create an immutable pair.
	pointer immutable_cons(
		pointer a,			///<	The pair car.
		pointer b			///<	The pair cdr.
		) { return cons(a, b, true); }

	///	Create an integer cell with the given value.
	pointer mk_integer(
		long num			///<	The value of the integer.
		) {
			return store_.mk_integer(num);
		}
	///	Create a real cell with the given value.
	pointer mk_real(
		double num			///<	The value of the real.
		) {
			return store_.mk_real(num);
		}
	///	Create a new symbol cell.  
	/**	This does not put the symbol in the symbol table. */	
	pointer mk_symbol(
		const std::string& name		///<	The symbol name.
		) {
			return store_.mk_symbol(name);
		}
	/// Make symbol or number, depending on input.
	pointer mk_atom(
		const char *str,				///<	The name or number.
		bool must_be_number = false,	///<	If true, it must be a number.
		int radix = 10					///<	If a number, the number base.
		);
	/// Make symbol or number, depending on input.
	pointer mk_atom(
		const std::string& str,					///<	The name or number
		bool must_be_number = false,	///<	If true, it must be a number.
		int radix = 10					///<	If a number, the number base.
		);
	///	Create a new internal symbol.
	/**	Returns the symbol cell. */
	pointer gensym();
	///	Create a string cell that holds (a copy of) the given string.
	/**	The string must be null terminated. */
	pointer mk_string(
		const std::string& str			///<	Copy initial value from this string.
		) {
			return store_.mk_string(str);
		}
	///	Create a string cell that holds (a copy of) the given string.
	/**	The string must be null terminated. */
	pointer mk_string(
		const char* str		///<	Copy this character string.
		) {
			return store_.mk_string(str);
		}
	///	Create a string cell that holds (a copy of) the given string.
	/**	The string is defined by its start and length, and may contain nulls. */
	pointer mk_string(
		const char* str,		///<	Copy this character string.
		int len					///<	The string length.  Use if not null terminated.
		) {
			return store_.mk_string(str, len);
		}
	///	Create a string cell that holds a string of the given length.
	/**	The string is initially filled with the fill character. */
	pointer mk_string(
		int len,				///<	The string length.
		char fill = ' '			///<	Fill it with this character.
		) {
			return store_.mk_string(len, fill);
		}
	///	Create a character cell.
	pointer mk_character(
		int c					///<	The character value.
		) {
			return store_.mk_character(c);
		}

	///	Create a vector cell, with the given length.
	pointer mk_vector(
		int len					///<	The vector length.
		) {
			return store_.mk_vector(len);
		}
	///	Create a foreign function cell.
	pointer mk_foreign_func(
		foreign_func f			///<	The foreign function pointer.
		) {
			return store_.mk_foreign_func(f);
		}
	///	Create a new environment frame, linking with the old environent.
	pointer mk_environment(
		pointer new_frame,		///<	The new environment.
		pointer old_env			///<	The encolsing environment.
		) {
			return store_.mk_environment(new_frame, old_env);
		}

	///	Create or look up a symbol.
	/**	If the symbol is in the symbol table, return it.
	 *	If not, create a symbol and add it to the symbol table.
	 *	In either case, return a symbol cell.
	 */
	pointer define_symbol(
		const std::string& name			///<	The symbol name (as string).
		);

	///	Returns true if the pointer is not #f.
	static bool is_true(pointer p) { return p != F_; }
	///	Returns true if the pointer is #f.
	/** () is #t in R5RS. */
	static bool is_false(pointer p) { return p == F_; }

	///	Returns the global environment.
	pointer global_env() const { return global_env_; }
	///	Returns the current environment.
	pointer envir() const { return env_.curr(); }
	///	Returns the retcode.
	/**	retcode values: <br>
	 *	   0:	normal <br>
	 *	  -1:	called (error) <br>
	 *	   1:	load error <br>
	 *	other:	called (quit n) for some n
	 */
	int retcode() const { return retcode_; }

	///	Return the reverse if a list.
	pointer reverse(
		pointer a			///<	The list to reverse.
		);
	///	Returns the list a appended with b.
	pointer append(
		pointer a,			///<	The first list.
		pointer b			///<	This is appended to the end of a.
		);

	///	Create a symbol for this name, add it to the environmet, and set its syntax flag.
	/**	This is used to associate a name with a special form.
	 *	Return the symbol cell.
	 */
	pointer assign_syntax(
		const char* name		///<	The symbol name.
		);
	///	Create a symbol for the name, add it to the environment, and associate
	///	it with the Op.
	/**	Return the resulting proc cell. */
	pointer assign_proc(
		Op op,				///<	The op code to execute for this proc.
		const std::string& name	///<	The name of the proc.
		);

	static const pointer NIL_;             ///< special cell representing empty cell 
	static const pointer T_;               ///<special cell representing #t 
	static const pointer F_;               ///< special cell representing #f 
	static const pointer VOID_;            ///< special cell representing void 
	static const pointer EOF_OBJ_;         ///< special cell representing end-of-file object 

protected:

	enum symbol_type { type_symbol, type_syntax };

	//	helper functions for ops
	inline pointer arg0() const	{ return args_->car(); }
	inline pointer arg_tail() const { return args_->cdr(); }
	inline pointer arg1() const	{ return args_->cadr(); }
	inline pointer arg2() const	{ return args_->caddr(); }
	inline bool arg1_nil() const	{ return args_->cdr() == NIL_; }
	inline pointer arg0_code() const { return code_->car(); }
	inline pointer code_tail() const { return code_->cdr(); }
	inline pointer arg1_code() const { return code_->cadr(); }


	//	Some tokens for the parser
	enum {
		TOK_EOF =    (-1),
		TOK_LPAREN = 0,
		TOK_RPAREN,
		TOK_DOT,
		TOK_ATOM,
		TOK_QUOTE,
		TOK_COMMENT,
		TOK_DQUOTE,
		TOK_BQUOTE,
		TOK_COMMA,
		TOK_ATMARK,
		TOK_SHARP,
		TOK_SHARP_CONST,
		TOK_VEC
	};


	//	internal functions
	bool ok_abbrev(pointer p) { return  p->is_pair() && p->cdr() == NIL_; }

	pointer list_star(pointer d);

	bool file_push(const std::string& fname);
	void file_pop();
	bool file_interactive();

	pointer mk_number(num n) {
		return store_.mk_number(n);
	}


	pointer mk_sharp_const(const char *name);

	//	port
	pointer mk_port(port *p);
	void free_port(port* port);
	pointer mk_port(FILE* file, int prop);
	pointer mk_port(const std::string& filename, int prop);
	pointer mk_port(char *start, char *past_the_end, int prop);

	//	gc
	void mark_for_gc(pointer a, pointer b);

	//	input/output
	void putstr(const char *s) {	//	put out a string
	  outport_->portvalue()->putstr(s);
	}
	void putchars(const char *s, int len) {	// put out a counted string
	  outport_->portvalue()->putstr(s, len);
	}
	void putcharacter(int c) {		//	put out a character
	  outport_->portvalue()->putcharacter(c);
	}
	int inchar();				//	read a character
	void backchar(int c) {		// put back character to input buffer 
	  if(c == EOF) return;
	  inport_->portvalue()->backchar(c);
	}
	char* readstr_upto(char *delim);
	pointer readstrexp(char* buff, int buf_len);
	void skipspace();
	int token();
	void record_output() {if (file_interactive()) interactive_out_ = true; }

	const char* atom2str(pointer l, bool f, int radix = 10);
	void printatom(pointer l, bool f);

	//	control structures
	pointer mk_proc(Op op) {
		return store_.mk_proc(op);
	}
	pointer mk_closure(pointer c, pointer e) {
		return store_.mk_closure(c, e);
	}
	pointer mk_promise(pointer c, pointer e) {
		return store_.mk_promise(c, e);
	}
	pointer mk_continuation(pointer d) {
		return store_.mk_continuation(d);
	}

	//	op functions
	void s_error_action(const char *s, pointer a = 0);
	void s_goto_action(Op op) { op_ = op; }
	void s_return_action(pointer a);
	void s_save(Op op, pointer args = NIL_, pointer code = NIL_);

	//	tests the number of args to builtins
	//	sets op_ if there is a problem
	void test_builtin();

	void write_common();
	bool check_nesting();

	void stop_execution() { break_ = false; }
	void reset_execution() { break_ = true; }
	bool continue_execution() const { return break_; }

	//	Individual operations
	//	Each of these corresponds to an Op code
	//    and is executed to perform that operation.
	//	The dispatch_table describes the common argument
	//	  checking that takes place before executing the 
	//	  operation.
	friend op_code_table;
public:
	void op_load(scheme& sc);
	void op_t0lvl(scheme& sc);
	void op_t1lvl(scheme& sc);
	void op_read_internal(scheme& sc);
	void op_gensym(scheme& sc);
	void op_valueprint(scheme& sc);
	void op_eval(scheme& sc);
	void op_real_eval(scheme& sc);
	void op_e0args(scheme& sc);
	void op_e1args(scheme& sc);
	void op_apply(scheme& sc);
	void op_real_apply(scheme& sc);
#if USE_TRACING
	void op_tracing(scheme& sc);
#endif
	void op_domacro(scheme& sc);
	void op_lambda(scheme& sc);
	void op_mkclosure(scheme& sc);
	void op_quote(scheme& sc);
	void op_def0(scheme& sc);
	void op_def1(scheme& sc);
	void op_defp(scheme& sc);
	void op_begin(scheme& sc);
	void op_if0(scheme& sc);
	void op_if1(scheme& sc);
	void op_set0(scheme& sc);
	void op_set1(scheme& sc);
	void op_let0(scheme& sc);
	void op_let1(scheme& sc);
	void op_let2(scheme& sc);
	void op_let0star(scheme& sc);
	void op_let1star(scheme& sc);
	void op_let2star(scheme& sc);
	void op_let0rec(scheme& sc);
	void op_let1rec(scheme& sc);
	void op_let2rec(scheme& sc);
	void op_cond0(scheme& sc);
	void op_cond1(scheme& sc);
	void op_delay(scheme& sc);
	void op_and0(scheme& sc);
	void op_and1(scheme& sc);
	void op_or0(scheme& sc);
	void op_or1(scheme& sc);
	void op_c0stream(scheme& sc);
	void op_c1stream(scheme& sc);
	void op_macro0(scheme& sc);
	void op_macro1(scheme& sc);
	void op_case0(scheme& sc);
	void op_case1(scheme& sc);
	void op_case2(scheme& sc);
	void op_peval(scheme& sc);
	void op_papply(scheme& sc);
	void op_continuation(scheme& sc);
	void op_inex2ex(scheme& sc);
	void op_ex2inex(scheme& sc);
	void op_exp(scheme& sc);
	void op_log(scheme& sc);
	void op_sin(scheme& sc);
	void op_cos(scheme& sc);
	void op_tan(scheme& sc);
	void op_asin(scheme& sc);
	void op_acos(scheme& sc);
	void op_atan(scheme& sc);
	void op_sqrt(scheme& sc);
	void op_expt(scheme& sc);
	void op_floor(scheme& sc);
	void op_ceiling(scheme& sc);
	void op_truncate(scheme& sc);
	void op_round(scheme& sc);
	void op_exact(scheme& sc);
	void op_inexact(scheme& sc);
	void op_odd(scheme& sc);
	void op_even(scheme& sc);
	void op_zero(scheme& sc);
	void op_positive(scheme& sc);
	void op_negative(scheme& sc);
	void op_add(scheme& sc);
	void op_sub(scheme& sc);
	void op_mul(scheme& sc);
	void op_div(scheme& sc);
	void op_intdiv(scheme& sc);
	void op_rem(scheme& sc);
	void op_mod(scheme& sc);
	void op_max(scheme& sc);
	void op_min(scheme& sc);
	void op_abs(scheme& sc);
	void op_gcd(scheme& sc);
	void op_lcm(scheme& sc);
	void op_car(scheme& sc);
	void op_cdr(scheme& sc);
	void op_cons(scheme& sc);
	void op_setcar(scheme& sc);
	void op_setcdr(scheme& sc);
	void op_char2int(scheme& sc);
	void op_int2char(scheme& sc);
	void op_charupcase(scheme& sc);
	void op_chardncase(scheme& sc);
	void op_sym2str(scheme& sc);
	void op_atom2str(scheme& sc);
	void op_str2sym(scheme& sc);
	void op_str2atom(scheme& sc);
	void op_mkstring(scheme& sc);
	void op_strlen(scheme& sc);
	void op_strref(scheme& sc);
	void op_strset(scheme& sc);
	void op_substr(scheme& sc);
	void op_strappend(scheme& sc);
	void op_str2list(scheme& sc);
	void op_string(scheme& sc);
	void op_list2str(scheme& sc);
	void op_str2num(scheme& sc);
	void op_num2str(scheme& sc);
	void op_strfill(scheme& sc);
	void op_list(scheme& sc);
	void op_listtail(scheme& sc);
	void op_listref(scheme& sc);
	void op_memq(scheme& sc);
	void op_memv(scheme& sc);
	void op_member(scheme& sc);
	void op_assq(scheme& sc);
	void op_assv(scheme& sc);
	void op_assoc(scheme& sc);
	void op_vector(scheme& sc);
	void op_mkvector(scheme& sc);
	void op_veclen(scheme& sc);
	void op_vecref(scheme& sc);
	void op_vecset(scheme& sc);
	void op_list2vec(scheme& sc);
	void op_vecfill(scheme& sc);
	void op_vec2list(scheme& sc);
	void op_not(scheme& sc);
	void op_boolp(scheme& sc);
	void op_eofobjp(scheme& sc);
	void op_nullp(scheme& sc);
	void op_comp(scheme& sc);
	void op_symbolp(scheme& sc);
	void op_numberp(scheme& sc);
	void op_stringp(scheme& sc);
	void op_integerp(scheme& sc);
	void op_realp(scheme& sc);
	void op_charp(scheme& sc);
	void op_charap(scheme& sc);
	void op_charnp(scheme& sc);
	void op_charwp(scheme& sc);
	void op_charup(scheme& sc);
	void op_charlp(scheme& sc);
	void op_portp(scheme& sc);
	void op_inportp(scheme& sc);
	void op_outportp(scheme& sc);
	void op_procp(scheme& sc);
	void op_pairp(scheme& sc);
	void op_listp(scheme& sc);
	void op_envp(scheme& sc);
	void op_vectorp(scheme& sc);
	void op_eq(scheme& sc);
	void op_eqv(scheme& sc);
	void op_equal(scheme& sc);
	void op_vecequal(scheme& sc);
	void op_force(scheme& sc);
	void op_save_forced(scheme& sc);
	void op_write(scheme& sc);
	void op_write_char(scheme& sc);
	void op_display(scheme& sc);
	void op_newline(scheme& sc);
	void op_err0(scheme& sc);
	void op_err1(scheme& sc);
	void op_reverse(scheme& sc);
	void op_list_star(scheme& sc);
	void op_append(scheme& sc);
	void op_put(scheme& sc);
	void op_get(scheme& sc);
	void op_quit(scheme& sc);
	void op_gc(scheme& sc);
	void op_gcverb(scheme& sc);
	void op_newsegment(scheme& sc);
	void op_oblist(scheme& sc);
	void op_curr_inport(scheme& sc);
	void op_curr_outport(scheme& sc);
	void op_open_infile(scheme& sc);
	void op_open_outfile(scheme& sc);
	void op_open_inoutfile(scheme& sc);
	void op_open_instring(scheme& sc);
	void op_open_outstring(scheme& sc);
	void op_open_inoutstring(scheme& sc);
	void op_close_inport(scheme& sc);
	void op_close_outport(scheme& sc);
	void op_int_env(scheme& sc);
	void op_curr_env(scheme& sc);
	void op_read(scheme& sc);
	void op_read_char(scheme& sc);	// also peek-char
	void op_char_ready(scheme& sc);
	void op_set_inport(scheme& sc);
	void op_set_outport(scheme& sc);
	void op_rdsexpr(scheme& sc);
	void op_rdlist(scheme& sc);
	void op_rddot(scheme& sc);
	void op_rdquote(scheme& sc);
	void op_rdqquote(scheme& sc);
	void op_rdqquotevec(scheme& sc);
	void op_rdunquote(scheme& sc);
	void op_rduqtsp(scheme& sc);
	void op_rdvec(scheme& sc);
	void op_p0list(scheme& sc);
	void op_p1list(scheme& sc);
	void op_pvecfrom(scheme& sc);
	void op_list_length(scheme& sc);
	void op_get_closure(scheme& sc);
	void op_closurep(scheme& sc);
	void op_macrop(scheme& sc);
	void op_time(scheme& sc);
	void op_illegal(scheme& sc);

	const char* procname(pointer x);

	//	test table
	typedef bool (*test_predicate)(pointer);
	class test_entry {
	public:
		test_entry() : test_(0), fct_(0), kind_(0) {}
		test_entry(char test, test_predicate fct, const char* kind) :
		  test_(test),
		  fct_(fct),
		  kind_(kind) {}
		static bool check();

		char test_;
		test_predicate fct_;
		const char* kind_;
	};

	bool extend_tests(test_entry* tests);
	const test_entry& test(int i) const {
		return tests_[i];
	}

	//	print table
	typedef const char* (*print_fn)(scheme&, pointer, bool, int, char*);
	class print_entry {
	public:
		print_entry() : test_(0), print_(0) {}
		print_entry(char test, print_fn print) :
		    test_(test),
			print_(print) {}
		char test_;
		print_fn print_;
	};
	void extend_print(print_entry* print);

	//	syntax table
	void extend_syntax(char** stntax_table);

protected:
	//	returns the retcode
	int load();

	int eval_cycle_top_level();
	int eval_cycle_eval();
	int eval_cycle_apply();

	//	returns number of opcodes executed
	int eval_cycle(Op op);

	//////////////////
	//// variables ///
	//////////////////

	int retcode_;		// return code
	bool break_;		// if true, breaks out of eval loop

	//	flags that control tracing
	enum {
		trace_none = 0,
		trace_eval = 0x0001,
		trace_op =	 0x0002
	};
	int tracing_;		// if not zero, tracing is on

	//	test table
	std::vector<test_entry> tests_;

	//	print table
	//	Each cell type has a print function.
	//	This maps the cell type to the function to call.
	typedef std::map<char, print_fn> print_map_type;
	typedef std::pair<char, print_fn> print_map_key;
	print_map_type print_map_;

	//	this could be passed in if two interpreters wanted
	//	  to share a store
	storage store_;

	// Main interpreter state
	pointer args_;            // arguments of function 
	pointer code_;            // current code 
	dump_stack dump_;		  // stack of continuatiuons
	environment env_;		  // execution enviroment

	op_code_table op_table_;

	bool interactive_repl_;    // are we in an interactive REPL? 

	static cell cellNIL_;
	static cell cellT_;
	static cell cellF_;
	static cell cellVOID_;
	static cell cellEOF_OBJ_;

	oblist oblist_;			  //  symbol table
	pointer global_env_;      // pointer to global environment 

	// global pointers to special symbols 
	pointer LAMBDA_;          // pointer to syntax lambda 
	pointer QUOTE_;           // pointer to syntax quote 
	pointer QQUOTE_;          // pointer to symbol quasiquote 
	pointer UNQUOTE_;         // pointer to symbol unquote 
	pointer UNQUOTESP_;       // pointer to symbol unquote-splicing 
	pointer FEED_TO_;         // => 
	pointer COLON_HOOK_;      // *colon-hook* 
	pointer ERROR_HOOK_;      // *error-hook* 
	pointer SHARP_HOOK_;		 // *sharp-hook* 

	pointer inport_;		  // the input port
	pointer outport_;		  // the output port
	pointer save_inport_;	  // saved input port
	pointer loadport_;		  // load port

	//	stack of load ports, for nested loading
	class port_stack {
	public:
		port_stack() : 
		  i_(0) {
			nesting_[0] = 0;
		  }
		port* reset() { i_ = 0; return &load_[i_]; }
		bool base() { return i_ == 0; }
		port* curr_port() { return &load_[i_]; }
		int curr_nesting() { return nesting_[i_]; }
		FILE* curr_file() { return load_[i_].file(); }
		//	these keep track of load nesting
		port* push() {
			i_++;
			nesting_[i_] = 0;
			return &load_[i_];
		}
		port* pop() {
			i_--;
			return &load_[i_];
		}
		//	these keep track of parenthes nesting
		void nesting_inc() { nesting_[i_]++; }
		void nesting_dec() { nesting_[i_]--; }
	private:
		enum { MAXFIL = 64 };
		int i_;
		port load_[MAXFIL];
		int nesting_[MAXFIL];
	};
	port_stack load_;		//	stack of load state
	int nesting_;			// global parenthesis nesting level

	char    strbuff_[256];		// temporary used for forming messages

	int tok_;				//  last token read -- used in reading S expressions
	bool print_flag_;
	pointer value_;
	Op op_;					//	the machine operation currently executing

	void *ext_data_;			// for the benefit of foreign functions 
	long gensym_cnt_;			// stores the state of the gensym

	bool interactive_out_;		// true if someting has been produced on stdout;

	friend storage;
	friend dump_stack_safe;
};

//	declaration of helper functions
const char* procname(pointer x);

inline bool is_any(pointer p)			{ return true;}
inline bool is_num_integer(pointer p)	{ return p->is_number(); }
inline bool is_nonneg(pointer p)		{ return p->is_nonneg_integer(); }
inline char type(pointer p)				{ return p->type(); }
inline bool is_type(pointer p, char t)	{ return p->is_type(t); }
inline bool is_atom(pointer p)			{ return p->is_atom(); }
inline bool mark_fun(pointer p)			{ return p->mark_fun(); }
inline bool final_fun(pointer p)		{ return p->final_fun(); }
inline bool is_string(pointer p)		{ return p->is_string(); }
inline std::string& strvalue(pointer p)	{ return p->strvalue(); }
inline int strlength(pointer p)			{ return p->strlength(); }
inline bool is_vector(pointer p)		{ return p->is_vector(); }
inline long veclen(pointer p)			{ return p->veclen(); }
inline pointer vector_elem(pointer p, int ielem) { return p->vector_elem(ielem); }
inline pointer set_vector_elem(pointer p, int ielem, pointer a) {
	return p->set_vector_elem(ielem, a); 
}
inline bool is_number(pointer p)		{ return p->is_number(); }
inline bool is_integer(pointer p)		{ return p->is_integer(); }
inline long ivalue(pointer p)			{ return p->ivalue(); }
inline void set_integer(pointer p, long ivalue) { p->set_integer(ivalue); }
inline bool is_real(pointer p)			{ return p->is_real(); }
inline double rvalue(pointer p)			{ return p->rvalue(); }
inline void set_real(pointer p, double rvalue)	 { p->set_real(rvalue); }
inline bool is_character(pointer p)		{ return p->is_character(); }
inline char charvalue(pointer p)		{ return p->charvalue(); }

inline bool is_port(pointer p)			{ return p->is_port(); }
inline bool is_inport(pointer p)		{ return p->is_inport(); }
inline bool is_outport(pointer p)		{ return p->is_outport(); }

inline bool is_pair(pointer p)			{ return p->is_pair(); }
inline bool is_pair_or_nil(pointer p)	{ return p->is_pair() || p == scheme::NIL_; }
inline pointer car(pointer p)			{ return p->car(); }
inline pointer cdr(pointer p)			{ return p->cdr(); }
inline pointer set_car(pointer p, pointer q) { return p->set_car(q); }
inline pointer set_cdr(pointer p, pointer q) { return p->set_cdr(q); }
inline void set_cell(pointer p, pointer car, pointer cdr) {
	p->set_cell(car, cdr);
}

inline bool is_symbol(pointer p)		{ return p->is_symbol(); }
inline const std::string& symname(pointer p)	{ return p->symname(); }
inline bool hasprop(pointer p)			{ return p->hasprop(); }
inline pointer symprop(pointer p)		{ return p->symprop(); }
inline pointer set_symprop(pointer p, pointer q) { return p->set_symprop(q); }

inline bool is_syntax(pointer p)		{ return p->is_syntax(); }
inline bool is_proc(pointer p)			{ return p->is_proc(); }
inline bool is_foreign(pointer p)		{ return p->is_foreign(); }

inline bool is_closure(pointer p)		{ return p->is_closure(); }
inline bool is_macro(pointer p)			{ return p->is_macro(); }
inline bool is_promise(pointer p)		{ return p->is_promise(); }
inline bool is_environment(pointer p)	{ return p->is_environment(); }
inline bool is_continuation(pointer p)	{ return p->is_continuation(); }

inline bool is_immutable(pointer p)		{ return p->is_immutable(); }
inline void set_immutable(pointer p)	{ p->set_immutable(); }

inline pointer caar(pointer p)			{ return  p->caar(); }
inline pointer cadr(pointer p)			{ return  p->cadr(); }
inline pointer cdar(pointer p)			{ return  p->cdar(); }
inline pointer cddr(pointer p)			{ return  p->cddr(); }
inline pointer cadar(pointer p)			{ return  p->cadar(); }
inline pointer caddr(pointer p)			{ return  p->caddr(); }
inline pointer cadaar(pointer p)		{ return  p->cadaar(); }
inline pointer cadddr(pointer p)		{ return  p->cadddr(); }
inline pointer cddddr(pointer p)		{ return  p->cddddr(); }

pointer reverse_in_place(pointer list, pointer term = 0);

#endif

