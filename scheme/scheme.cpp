//
//	C++ version by Charles Hayden (chayden@comcast.net)
//
/* T I N Y S C H E M E    1 . 3 1
 *   Dimitrios Souflis (dsouflis@acm.org)
 *   Based on MiniScheme (original credits follow)
 * (MINISCM)               coded by Atsushi Moriwaki (11/5/1989)
 * (MINISCM)           E-MAIL :  moriwaki@kurims.kurims.kyoto-u.ac.jp
 * (MINISCM) This version has been modified by R.C. Secrist.
 * (MINISCM)
 * (MINISCM) Mini-Scheme is now maintained by Akira KIDA.
 * (MINISCM)
 * (MINISCM) This is a revised and modified version by Akira KIDA.
 * (MINISCM)	current version is 0.85k4 (15 May 1994)
 *
 */
#include "scheme.h"

#include <math.h>

#include <limits.h>
#include <float.h>
#include <ctype.h>
#include <time.h>


# define BACKQUOTE '`'

#include <string.h>
#include <stdlib.h>

#ifndef PROMPT
# define PROMPT "> "
#endif

bool streq(const char* a, const char*b) { return strcmp(a, b) == 0; }
bool streq_nc(const char* a, const char*b) { return _stricmp(a, b) == 0; }

static num num_zero_;
static num num_one_;


cell scheme::cellNIL_;
const pointer scheme::NIL_(&cellNIL_);
cell scheme::cellT_;
const pointer scheme::T_(&cellT_);
cell scheme::cellF_;
const pointer scheme::F_(&cellF_);
cell scheme::cellVOID_;
const pointer scheme::VOID_(&cellVOID_);
cell scheme::cellEOF_OBJ_;
const pointer scheme::EOF_OBJ_(&cellEOF_OBJ_);
cell storage::cellSink_;
const pointer storage::sink_(&cellSink_);

//	These are convenience definitions
//	  for when the "object" notation is inconvenient

Op procvalue(pointer p)				{ return p->procvalue(); }
pointer closure_code(pointer p)		{ return p->closure_code(); }
pointer closure_env(pointer p)		{ return p->closure_env(); }
pointer continuation(pointer p)		{ return p->continuation(); }

bool flag_mark(pointer p)			{ return p->flag_mark(); }

bool is_frame(pointer p)			{ return p->is_frame(); }

pointer free_cdr(pointer p)			{ return p->free_cdr(); }


#if USE_CHAR_CLASSIFIERS
static int Cisalpha(int c) { return isascii(c) && isalpha(c); }
static int Cisdigit(int c) { return isascii(c) && isdigit(c); }
static int Cisspace(int c) { return isascii(c) && isspace(c); }
static int Cisupper(int c) { return isascii(c) && isupper(c); }
static int Cislower(int c) { return isascii(c) && islower(c); }
#endif

#if USE_ASCII_NAMES
static const char* charnames[32]={
	 "nul",
	 "soh",
	 "stx",
	 "etx",
	 "eot",
	 "enq",
	 "ack",
	 "bel",
	 "bs",
	 "ht",
	 "lf",
	 "vt",
	 "ff",
	 "cr",
	 "so",
	 "si",
	 "dle",
	 "dc1",
	 "dc2",
	 "dc3",
	 "dc4",
	 "nak",
	 "syn",
	 "etb",
	 "can",
	 "em",
	 "sub",
	 "esc",
	 "fs",
	 "gs",
	 "rs",
	 "us"
};
static int size_charnames = sizeof(charnames)/sizeof(char*);

static bool is_ascii_name(const char* name, int *pc) {
  for(int i = 0; i < size_charnames; i++) {
     if(streq_nc(name, charnames[i])) {
          *pc = i;
          return true;
     }
  }
  if(streq_nc(name, "del")) {
     *pc = 127;
     return true;
  }
  return false;
}

#endif

//	forward
static Op syntaxnum(pointer p);

////////////////////////
///  num arithmetic  ///
////////////////////////
inline static long num_ivalue(num n) { 
	return n.as_integer(); 
}

inline static double num_rvalue(num n) { 
	return n.as_real();
}

num::num(pointer p) {
	if (p->is_integer()) {
		set(p->ivalue());
	} else {
		set(p->rvalue());
	}
}

//	true if both are fixed
static bool both_int(const num& a, const num& b) {
	return a.is_integer() && b.is_integer();
}

//	define arithmetic operators for integer and real operations
typedef long (*i_binop)(long, long);
typedef double (*d_binop)(double, double);

//	Combine two num values using the supplied operators.
//	If both operands are integers, use the integer op,
//	  otherwise convert to real and use the real op.
num op(num a, num b, i_binop iop, d_binop dop) {
	if (both_int(a, b)) {
		return num(iop(a.ivalue(), b.ivalue()));
	} else {
		return num(dop(a.as_real(), b.as_real()));
	}
}

//	add
inline static long iadd(long a, long b) { return a+b; }
inline static double dadd(double a, double b) { return a+b; }
num num::add(num b) {
	return op(*this, b, iadd, dadd);
}

//	sub
inline static long isub(long a, long b) { return a-b; }
inline static double dsub(double a, double b) { return a-b; }
num num::sub(num b) {
	return op(*this, b, isub, dsub);
}

//	mul
inline static long imul(long a, long b) { return a*b; }
inline static double dmul(double a, double b) { return a*b; }
num num::mul(num b) {
	return op(*this, b, imul, dmul);
}

num num::div(num b) {
	if(both_int(*this, b) && (ivalue() % b.ivalue()) == 0) {
		return num(ivalue() / b.ivalue());
	} else {
		return num(as_real() / b.as_real());
	}
}

//	div
inline static long idiv(long a, long b) { return a/b; }
inline static double ddiv(double a, double b) { return a/b; }
num num::intdiv(num b) {
	return op(*this, b, idiv, ddiv);
}

num num::rem(num b) {
	if (both_int(*this, b)) {
		long e1 = ivalue();
		long e2 = b.ivalue();
		long res = e1 % e2;
		if(res*e1 < 0) {    // remainder should have same sign as first operand 
			e2 = labs(e2);
			if(res > 0) {
			  res -= e2;
			} else {
			  res += e2;
			}
		}
		return num(res);
	} else {
		return num(0.0);
	}
}

num num::mod(num b) {
	if (both_int(*this, b)) {
		long e1 = ivalue();
		long e2 = b.ivalue();
		long res = e1 % e2;
		if(res*e2 < 0) {    // modulo should have same sign as second operand 
			e2 = labs(e2);
			if(res > 0) {
			  res -= e2;
			} else {
			  res += e2;
			}
		}
		return num(res);
	} else {
		return num(0.0);
	}
}


//	max
inline static long imax(long a, long b) { return (a > b) ? a : b; }
inline static double dmax(double a, double b) { return (a > b) ? a : b; }
num num::max(num b) {
	return op(*this, b, imax, dmax);
}

//	min
inline static long imin(long a, long b) { return (a > b) ? b : a; }
inline static double dmin(double a, double b) { return (a > b) ? b : a; }
num num::min(num b) {
	return op(*this, b, imin, dmin);
}

inline static long iabs(long a, long b) { return (a >= 0) ? a : -a; }
inline static double dabs(double a, double b) { return (a >0.0) ? a : -a; }
num num::abs() {
	return op(*this, num_zero_, iabs, dabs);
}

//	define predicates for integer and real operations
typedef bool (*i_comp)(long,long);
typedef bool (*d_comp)(double,double);

//	Compare two num values using the supplied operators.
//	If both operands are integers, use the integer op,
//	  otherwise convert to real and use the real op.
bool comp(num a, num b, i_comp icomp, d_comp dcomp) {
	if (both_int(a, b)) {
		return icomp(a.ivalue(), b.ivalue());
	} else {
		return dcomp(a.as_real(), b.as_real());
	}
}

//	eq
inline static bool ieq(long a, long b) { return a == b; }
inline static bool deq(double a, double b) { return a == b; }
bool num::eq(num b) {
	return comp(*this, b, ieq, deq);
}

//	gt
inline static bool igt(long a, long b) { return a > b; }
inline static bool dgt(double a, double b) { return a > b; }
bool num::gt(num b) {
	return comp(*this, b, igt, dgt);
}

//	gt
inline static bool ilt(long a, long b) { return a < b; }
inline static bool dlt(double a, double b) { return a < b; }
bool num::lt(num b) {
	return comp(*this, b, ilt, dlt);
}

inline static bool ipositive(long a, long b) { return a > 0; }
inline static bool dpositive(double a, double b) { return a > 0.0; }
bool num::positive() {
	return comp(*this, num_zero_, ipositive, dpositive);
}

inline static bool inegative(long a, long b) { return a < 0; }
inline static bool dnegative(double a, double b) { return a < 0.0; }
bool num::negative() {
	return comp(*this, num_zero_, inegative, dnegative);
}

inline static bool izero(long a, long b) { return a == 0; }
inline static bool dzero(double a, double b) { return a == 0.0; }
bool num::zero() {
	return comp(*this, num_zero_, izero, dzero);
}

//	some two-arg versions for convenience
inline static bool num_eq(num a, num b) { return a.eq(b); }
inline static bool num_gt(num a, num b) { return a.gt(b); }
inline static bool num_ge(num a, num b) { return !a.lt(b); }
inline static bool num_lt(num a, num b) { return a.lt(b); }
inline static bool num_le(num a, num b) { return !num_gt(a,b); }

//	return the appropriate compare function for the given operator
typedef bool (*comp_func)(num, num);
static comp_func get_comp(Op op) {
    comp_func cf;
	switch(op.as_long()) {
	   case Op::NUMEQ: cf = num_eq; break;
	   case Op::LESS:  cf = num_lt; break;
	   case Op::GRE:   cf = num_gt; break;
	   case Op::LEQ:   cf = num_le; break;
	   case Op::GEQ:   cf = num_ge; break;
	   default:		   cf = 0;
	}
	return cf;
}

// Round to nearest. Round to even if midway 
static double round_per_R5RS(double x) {
	double fl = floor(x);	// return one of fl, ce
	double ce = ceil(x);
	double dfl = x-fl;		// dist to floor
	double dce = ce-x;		// dist to ceil
	if(dfl > dce) {
		return ce;
	} else if(dfl < dce) {
		return fl;
	} else {				// dfl == dce
		if(fmod(fl, 2.0) == 0.0) {       
			return fl;		// if even, return it
		} else {
			return ce;
		}
	}
}

inline static int is_zero_double(double x) {
 return x < DBL_MIN && x > -DBL_MIN;
}

/////////////////
///  storage  ///
/////////////////

storage::storage(int max_cell_seg) : 
	last_cell_seg_(-1),
	max_cell_seg_(max_cell_seg),
	free_cell_(scheme::NIL_),
	fcells_(0),
	no_memory_(true),
	gc_verbose_(false) {
		cell_seg_ = new pointer[max_cell_seg_];
	}

storage::~storage() {
	for (int i = 0; i <= last_cell_seg_; i++) {
		pointer p = cell_seg_[i];
		while (p < cell_seg_[i] + CELL_SEGSIZE) {
			if (! p->is_clear()) {
				finalize_cell(p);
			}
			p++;
		}
		delete cell_seg_[i];
	}
	delete [] cell_seg_;
}

bool storage::init(scheme& sc, int initial_segs) {
	sc_ = &sc;
	int n = alloc_cells(initial_segs);
#if CHECK
	assert(n == initial_segs);
#endif
	no_memory_ = (n < initial_segs);
	return no_memory_;
}

void storage::add_free_cells(pointer first, pointer last) {
	// insert new cells in address order on free list 
	if (free_cell_ == scheme::NIL_ || last < free_cell_) {
		//	new cells get put at front of free_cell list
		last->clear(0, free_cell_);
		free_cell_ = first;
	} else {
		// splice new list into free list
		pointer p = free_cell_;
		while (free_cdr(p) != scheme::NIL_ && first > free_cdr(p))
		   p = free_cdr(p);
		//	link new cells into free_cell list in the middle
		last->set_free_cdr(free_cdr(p));
		p->set_free_cdr(first);
	}
}

//	clear all cells in the block and link them together
static void clear_free_cells(pointer first, pointer last) {
    pointer p;	
	for (p = first; p < last; p++) {
		p->clear(scheme::NIL_, p+1);
	}
	//	also clear cell pointed to by last
	p->clear(scheme::NIL_, scheme::NIL_);
}

void storage::insert_new_seg(pointer newp) {
	long i = ++last_cell_seg_;

	//  Insert new segment in address order 
	//	Each cell_seg entry holds the address of a contiguous
	//	  block of cell memory.
	//	Existing blocks are sorted in ascending memory address value.
	cell_seg_[i] = newp;
	while (i > 0 && cell_seg_[i - 1] > cell_seg_[i]) {
		//	new block gets moved into proper position
        pointer p = cell_seg_[i];
		cell_seg_[i] = cell_seg_[i - 1];
        cell_seg_[i - 1] = p;
		--i;
	}
	fcells_ += CELL_SEGSIZE;
}

//	Allocate n new cell segments.
//	Each segment holds CELL_SEGSIZE cells.
//	Returns the number of segments actually allocated.
int storage::alloc_cells(int n_seg) {
     for (int k = 0; k < n_seg; k++) {
		 if (last_cell_seg_ > max_cell_seg_) {
			return k;		// out of segments
		 }
		//	Allocate space for CELL_SEGSIZE cells.
        pointer newp = new cell[CELL_SEGSIZE];
        if (newp == 0) return k;		// out of memory
		//	Insert the new cell segment into cell_seg
		insert_new_seg(newp);
		pointer last = newp + CELL_SEGSIZE - 1;		// points to last cell
		//	Make the memory block into cells and link them together.
		clear_free_cells(newp, last);
		//	Add the new cells to the free list.
		add_free_cells(newp, last);
     }
     return n_seg;
}

//	Get more cell storage.
//	a and b are protected in case a gc is needed.
//	Either free_cell_ != NIL or no_memory_ is true;
void storage::get_more_cells(pointer a, pointer b) {
	if(no_memory_) return;

	if (free_cell_ == scheme::NIL_) {
		gc(a, b);
		if (fcells_ < last_cell_seg_*8 || free_cell_ == scheme::NIL_) {
			//  if only a few recovered, get more to avoid fruitless gc's 
			//	this might fail to allocate more cells, in which case get_free_cell
			//    will return 0 the next time it is called.
			alloc_cells();
		}
	}
}

//	Get a free cell from the free list.
//	If there is none, return 0.
//	The cell is cleared, even though it should already be clear in the free list.
pointer storage::get_free_cell() {
	if (free_cell_ != scheme::NIL_) {
		pointer x = free_cell_;
		free_cell_ = free_cdr(x);
		--fcells_;
		x->clear();
		return x;
	}
	return 0;
}

//	Finds a free cell.
//	Uses one of the in the free_cell list, or calls get_more_cells
//	  to allocate or free more cells.
//	a and b are protected if gc is needed.
//	If it cannot find a cel, it returns a dummy pointer.
pointer storage::find_cell(pointer a, pointer b) {
	pointer x = get_free_cell();
	if (x) return x;
	get_more_cells(a, b);
	x = get_free_cell();
	if (x) return x;
	//	If there are no more cells, set flags that will stop things.
	no_memory_ = true;
	sc_->stop_execution();
#if CHECK
	assert(false);
#endif
	//	Return a valid pointer so callers don't have to check.
	return sink_;
}

//////////////////////
///  symbol table  ///
//////////////////////

///	Return the hash function of the string.
int hash_fn(const std::string& str, int table_size) { 
	unsigned int hashed = 0; 
	int bits_per_int = sizeof(unsigned int)*8; 

	for (unsigned int i = 0; i < str.size(); i++) {
		/* letters have about 5 bits in them */ 
		hashed = (hashed<<5) | (hashed>>(bits_per_int-5)); 
		hashed ^= str.at(i); 
	} 
	return hashed % table_size; 
} 

pointer oblist::search_symbol(const std::string& name, pointer p) {
	for (pointer x = p; x != scheme::NIL_; x = cdr(x)) {
		if(streq_nc(name.c_str(), symname(car(x)).c_str())) {
			return car(x);
		}
	}
	return 0;
}


void oblist::init() { 
  oblist_ = store_->mk_vector(oblist_size);  
} 

// returns the new symbol 
pointer oblist::add_symbol(const std::string& name) { 
	pointer x = store_->mk_symbol(name);
	int location = hash_fn(name, veclen(oblist_)); 
	set_vector_elem(oblist_, location, 
				  store_->immutable_cons(x, vector_elem(oblist_, location))); 
	return x; 
} 

pointer oblist::find_symbol(const std::string& name) { 
	int location = hash_fn(name, veclen(oblist_)); 
	return search_symbol(name, vector_elem(oblist_, location));
} 

pointer oblist::all_symbols() { 
	pointer ob_list = scheme::NIL_; 

  for (int i = 0; i < veclen(oblist_); i++) { 
	  for (pointer x  = vector_elem(oblist_, i); x != scheme::NIL_; x = cdr(x)) { 
      ob_list = store_->cons(x, ob_list); 
    } 
  } 
  return ob_list; 
} 

void oblist::clear() { oblist_ = scheme::NIL_; }

////////////////////////////////////
///  storage cell creation  ///
////////////////////////////////////

static const char* file_rw(int prop) {
  if(prop == (port::port_input|port::port_output)) {
    return "a+";
  } 
  if(prop == port::port_output) {
    return "w";
  } 
  return "r";
}

// return new cons cell 
pointer storage::cons(pointer a, pointer b, bool immutable) {
	pointer x = find_cell(a, b);
	x->set_pair(a, b, immutable);
	return x;
}

//	make a new port cell
pointer storage::mk_port(port* p) {
	if (p == 0) return scheme::NIL_;
	pointer x = find_cell();
	x->set_port(p);
	return x;
}

//	make a new port cell
pointer storage::mk_port(FILE* f, int prop, bool closeit) {
	if (f == 0) return 0;
	pointer x = find_cell();
	port* pt = port::make_port(f, prop, closeit);
	x->set_port(pt);
	return x;
}

//	make a new port cell
pointer storage::mk_port(const std::string& filename, int prop) {
  const char* rw = file_rw(prop);
  FILE* f = fopen(filename.c_str(), rw);
  return mk_port(f, prop, true);
}

//	make a new port cell
pointer storage::mk_port(char* start, char* past_the_end, int prop) {
	pointer x = find_cell();
	port* pt = port::make_port(start, past_the_end, prop);
	x->set_port(pt);
	return x;
}

//	make a new foreign function cell
pointer storage::mk_foreign_func(foreign_func f) {
	pointer x = find_cell();
	x->set_foreign(f);
	return x;
}

//	make a new environment cell
pointer storage::mk_environment(pointer new_frame, pointer old_env) {
	pointer e = find_cell(new_frame, old_env);
	e->set_environment(new_frame, old_env);
	return e;
}

//	make a new character cell
pointer storage::mk_character(int c) {
	pointer x = find_cell();
	x->set_character(c);
	return x;
}

// make a new integer cell
pointer storage::mk_integer(long n) {
	pointer x = find_cell();
	x->set_integer(n);
	return x;
}

//	make a new real cell
pointer storage::mk_real(double d) {
	pointer x = find_cell();
	x->set_real(d);
	return x;
}

//	make a new number cell
pointer storage::mk_number(num n) {
	if(n.is_integer()) {
		return mk_integer(n.ivalue());
	} else {
		return mk_real(n.rvalue());
	}
}

//	make a new symbol cell
pointer storage::mk_symbol(const std::string& name) {
	pointer x = mk_string(name);
	set_immutable(x);
	pointer y = find_cell(x, scheme::NIL_);
	y->set_symbol(x, scheme::NIL_);
	return y;
}

//	make new string cell using given length
pointer storage::mk_string(const char* str, int len) {
     pointer x = find_cell();
	 x->set_string(str, len);
     return x;
}

// make new string cell using char*
pointer storage::mk_string(const char* str) {
     pointer x = find_cell();
	 x->set_string(str);
     return x;
}

// make new string cell using string
pointer storage::mk_string(const std::string& str) {
     pointer x = find_cell();
	 x->set_string(str);
     return x;
}

//	make a new string cell using fill and length
pointer storage::mk_string(int len, char fill) {
     pointer x = find_cell();
	 x->set_string(len, fill);
     return x;
}

//	make a new frame cell
pointer storage::mk_frame(Op op, pointer args, pointer envir, pointer code) {
	stack_frame* f = stack_frame::make_frame(op, args, envir, code);
	if (f == 0) return scheme::NIL_;
	pointer x = find_cell();
	x->set_frame(f);
	return x;
}

//	make a new vector cell
pointer storage::mk_vector(int len) {
	pointer x = find_cell();
	x->set_vector(len);
	x->fill_vector(scheme::NIL_);
	return x;
}

// make closure cell: c is code. e is environment 
pointer storage::mk_closure(pointer c, pointer e) {
     pointer x = find_cell(c, e);
     x->set_closure(c, e);
     return x;
}

// make promise cell: c is code. e is environment 
pointer storage::mk_promise(pointer c, pointer e) {
     pointer x = find_cell(c, e);
     x->set_promise(c, e);
     return x;
}

// make continuation cell 
pointer storage::mk_continuation(pointer d) {
     pointer x = find_cell(d);
     x->set_continuation(d);
     return x;
}

//	Make a proc cell with the given op
pointer storage::mk_proc(Op op) {
     pointer y = find_cell();
     y->set_proc(op);
     return y;
}

//////////////////////////
///  frame operations  ///
//////////////////////////

//	The use of the FRAME cell type is an optimization.
//	The frame stores the whole state in one structure, rather than consing
//	  it together into a list.  
//	The mark algorithm needs to know that it needs to follow all the pointers
//	  in the frame when marking this kind of cell.
//	Storage is allocated for frames from the same allocator, as with ports and strings.
//	The frame storage is deallocated along with the cell that refers to it.

//	Mark everything reachable from the current frame.
void stack_frame::mark(storage& st) {
	st.mark(args_);
	st.mark(envir_);
	st.mark(code_);
}

stack_frame* stack_frame::free_list_ = 0;

//	Define these so that the constructor can be used to
//	  initialize the stack_frame.
void* operator new(size_t, stack_frame* location) { return location; }
//	This never gets called, but is required to shut VC++ up.
void operator delete(void *pMem, stack_frame *rBuffer) {}


//	Factory for stack_frame
stack_frame* stack_frame::make_frame(Op op, pointer args, pointer envir, pointer code,
									 stack_frame* curr) { 
	if (free_list_ == 0) {
		//	allocate a new one
		return new stack_frame(op, args, envir, code, curr);
	} else {
		//	get one from the free list
		stack_frame* next = free_list_->link_;
		stack_frame* frame = new (free_list_) stack_frame(op, args, envir, code, curr);
		free_list_ = next;
		return frame;
	}
}

static void cleanup_list(stack_frame* frame) {
	stack_frame* next;
	for (stack_frame* curr = frame; curr != 0; curr = next) {
		next = curr->link();
		delete curr;
	}
}

static void free_list(stack_frame* frame) {
	stack_frame* next;
	for (stack_frame* curr = frame; curr != 0; curr = next) {
		next = curr->link();
		curr->free();
	}
}

void stack_frame::free() { 
	link_ = free_list_;
	free_list_ = this; 
}

void stack_frame::cleanup() {
	cleanup_list(free_list_);
}

////////////////
///  vector  ///
////////////////

//	Fills the vector with the given object.
//	The existing vector length is used.
void cell::fill_vector(pointer obj) {
     int num = veclen();
	 scheme_vector& vec = vecvalue();
     for(int i = 0; i < num; i++) {
          vec[i] = obj;
     }
}

//	gets a vector element by finding the cell and 
//	  extracting the appropriate part of it
pointer cell::vector_elem(int ielem) {
	return vecvalue()[ielem];
}

pointer cell::set_vector_elem(int ielem, pointer a) {
	vecvalue()[ielem] = a;
	return a;
}

////////////////
///  symbol  ///
////////////////

//	Look up symbol in symbol table.
//	If not found, enter it.
//	Return the symbol.
pointer scheme::define_symbol(const std::string& name) {
	pointer x = oblist_.find_symbol(name);
	if (x != 0) return x;
	return oblist_.add_symbol(name);
}

pointer scheme::gensym() {
     for(; gensym_cnt_ < LONG_MAX; gensym_cnt_++) {
	      char name_buf[40];
          sprintf(name_buf, "gensym-%ld", gensym_cnt_);
		  std::string name(name_buf);

          // first check oblist 
	 	  pointer x = oblist_.find_symbol(name);

          if (x != 0) {
               continue;		// duplicate - try again
          } else {
			   gensym_cnt_++;	// start at next one next time
			   return oblist_.add_symbol(name);
          }
     }

     return NIL_;
}

// make symbol or number atom from string 
pointer scheme::mk_atom(const char* str, bool must_be_number, int radix) {
     char *p;
     bool has_dec_point = false;
     bool has_fp_exp = false;

	 //	we need a lower case version of q
	 char* q = strbuff_;
	 strcpy(q, str);
	 _strlwr(q); 

#if USE_COLON_HOOK
     if((p = strstr(q, "::")) != 0) {
          *p = 0;
          return cons(COLON_HOOK_,
              cons(
                  cons(
                       QUOTE_,
                       cons(mk_atom(p+2, must_be_number, radix), NIL_)),
                  cons(define_symbol(q), NIL_)));
     }
#endif

     p = q;
	 //	this tests the input enough to know what kind
	 //	  of atom to make, and then makes it
     char c = *p++; 
     if ((c == '+') || (c == '-')) { 
       c = *p++; 
       if (c == '.') { 
         has_dec_point = true; 
		 c = *p++; 
       }
	   //	it might start with + or - and perhaps period
	   //	  but is not a number
       if (!isdigit(c)) { 
		if (must_be_number) return F_;
		return define_symbol(q); 
       } 
     } else if (c == '.') { 
       has_dec_point = true; 
       c = *p++; 
	   //	it starts with period, but is not a number
       if (!isdigit(c)) { 
		if (must_be_number) return F_;
		return define_symbol(q); 
       } 
     } else if (!isdigit(c)) { 
		 //	it does not look anything like a number
		if (must_be_number) return F_;
		return define_symbol(q); 
     }

	 //	ok -- it might be a number
	 //	see if consists solely of numbers, or looks like a real constant
     for (; (c = *p) != 0; ++p) {
          if (!isdigit(c)) {
               if(c=='.') {
                    if(!has_dec_point) {
                         has_dec_point = true;
                         continue;
                    }
               } else if (c == 'e') {
                       if(!has_fp_exp) {
                          has_dec_point = true; // decimal point illegal from now on 
                          p++;
                          if ((*p == '-') || (*p == '+') || isdigit(*p)) {
                             continue;
                          }
                       }  
               }
			   // it started like a number, but it is not
			   if (must_be_number) return F_;
               return define_symbol(q);
          }
     }
	 //	the whole thing looks like a nunber
     if(has_dec_point) {
		 //	has a decimal point -- make a float
		 if (radix != 10) return F_;
         return mk_real(strtod(q, 0));
     }
	 //	no decimal point -- make an integer
	 long ll = strtol(q, 0, radix);
	 if (ll == LONG_MAX || ll == LONG_MIN) {
		 //	too big to fit into an integer -- make a float
		 if (radix != 10) return F_;
         return mk_real(strtod(q, 0));
	 } else {
		return mk_integer(ll);
	 }
}

pointer scheme::mk_atom(const std::string& str, bool must_be_number, int radix) {
	return mk_atom(str.c_str(), must_be_number, radix);
}

inline static long decode_octal(const char* str) {
	return strtol(str, 0, 8);
}

inline static long decode_decimal(const char* str) {
	return strtol(str, 0, 10);
}

inline static long decode_hex(const char* str) {
	return strtol(str, 0, 16);
}

//	as long as s consists of 0s and 1s, built binary number
static long decode_binary(const char* str) {
	return strtol(str, 0, 2);
}


// make constant 
pointer scheme::mk_sharp_const(const char* name) {
     if (streq(name, "t"))
          return T_;
     else if (streq(name, "f"))
          return F_;
     else if (*name == 'o') {// #o (octal) 
          return mk_integer(decode_octal(name+1));
     } else if (*name == 'd') {    // #d (decimal) 
          return mk_integer(decode_decimal(name+1));
     } else if (*name == 'x') {    // #x (hex) 
          return mk_integer(decode_hex(name+1));
     } else if (*name == 'b') {    // #b (binary) 
          return mk_integer(decode_binary(name+1));
     } else if (*name == '\\') { // #\w (character) 
          int c = 0;
          if(streq_nc(name+1, "space")) {
               c=' ';
          } else if(streq_nc(name+1, "newline")) {
               c='\n';
          } else if(streq_nc(name+1, "return")) {
               c='\r';
          } else if(streq_nc(name+1, "tab")) {
               c='\t';
		  } else if(name[1]=='x' && name[2] != 0) {
			int c1 = 0;
			if(sscanf(name+2, "%x", &c1) == 1 && c1 < 256) {
				c = c1;
			} else {
				return NIL_;
			}
#if USE_ASCII_NAMES
          } else if(is_ascii_name(name+1, &c)) {
               // nothing 
#endif               
          } else if(name[2] == 0) {
               c = name[1];
          } else {
               return NIL_;
          }
          return mk_character(c);
     } else
          return NIL_;
}

// ========== garbage collector ========== 

static void mark_vector(storage& st, pointer p) {
	const scheme_vector& vec = p->vecvalue();
	for (unsigned int i = 0; i < vec.size(); i++) {
		st.mark(vec[i]);
	}
}

static void mark_frame(storage& st, pointer p) {
	p->framevalue()->mark(st);
}

static storage::mark_entry static_markers[] = {
	storage::mark_entry(cell::T_VECTOR, mark_vector),
	storage::mark_entry(cell::T_FRAME, mark_frame),
	storage::mark_entry(0, 0),
};

void storage::extend_marker(mark_entry* mark_table) {
	for (int i = 0; mark_table[i].test_ != 0; i++) {
		mark_entry& me = mark_table[i];
		mark_map_.insert(mark_map_key(me.test_, me.mark_));
	}
}

void storage::extend_finalizer(final_entry* final_table) {
	for (int i = 0; final_table[i].test_ != 0; i++) {
		final_entry& fe = final_table[i];
		final_map_.insert(final_map_key(fe.test_, fe.finalize_));
	}
}

/*--
 *  We use algorithm E (Knuth, The Art of Computer Programming Vol.1,
 *  sec 2.3.5, algorithm E) for marking.
 *	For clarity, a separate bit is used to store whether we went
 *	  down an alink.
 */
void storage::mark(pointer a) {
     pointer q;
     pointer t = 0;
     pointer p = a;
E2:  p->setmark();
	 if (p->mark_fun()) {
		//	Recursively mark ponters in vectors and frames.
		mark_map_type::iterator i = mark_map_.find(p->type());
		if (i != mark_map_.end()) {
			(*((*i).second))(*this, p);
		}
	 }
//E3:
     if (p->flag_atom())
          goto E6;
//E4:
     q = p->car();
     if (q && !q->flag_mark()) {
          p->set_alink();
          p->set_car(t);
          t = p;
          p = q;
          goto E2;
     }
E5:  q = p->cdr();
     if (q && !q->flag_mark()) {
          p->set_cdr(t);
          t = p;
          p = q;
          goto E2;
     }
E6:  if (!t)
          return;
     q = t;
     if (q->flag_alink()) {
          q->clr_alink();
          t = q->car();
          q->set_car(p);
          p = q;
          goto E5;
     } else {
          t = q->cdr();
          q->set_cdr(p);
          p = q;
          goto E6;
     }
}

//	Protect all the cells in use, prior to garbage collection
void scheme::mark_for_gc(pointer a, pointer b) {
	// mark system globals 
	oblist_.mark();
	store_.mark(global_env_);

	// mark current registers 
	store_.mark(args_);
	env_.mark();
	store_.mark(code_);
	dump_.mark(store_);
	store_.mark(value_);
	store_.mark(inport_);
	store_.mark(save_inport_);
	store_.mark(outport_);
	store_.mark(loadport_);

	// mark variables a, b 
	if (a) store_.mark(a);
	if (b) store_.mark(b);
}

//	garbage collection core algorithm
int storage::collect_garbage() {
	fcells_ = 0;
	free_cell_ = scheme::NIL_;
	/* Free-list is kept sorted by address . Here we scan the cells
	 * (which are also kept sorted by address) downwards to build the
	 * free-list in sorted order.
	 */
	for (int i = last_cell_seg_; i >= 0; i--) {
		pointer p = cell_seg_[i] + CELL_SEGSIZE;
		while (--p >= cell_seg_[i]) {
			if (p->flag_mark()) {
				p->clrmark();
			} else {
				if (! p->is_clear()) {
				// reclaim cell 
				finalize_cell(p);
				++fcells_;
			}
			set_cell(p, scheme::NIL_, free_cell_);
			free_cell_ = p;
			}
		}
	}
	return fcells_;
}

// garbage collection. 
//	parameter a, b are protected if set
int storage::gc(pointer a, pointer b) {
	if(gc_verbose_) {
		sc_->putstr("gc...");
	}

	sc_->mark_for_gc(a, b);

	//	garbage collect
	scheme::NIL_->clrmark();
	int fcells = collect_garbage();

	if (gc_verbose_) {
		char buf[100];
		sprintf(buf, "done: %ld cells were recovered.\n", fcells);
		sc_->putstr(buf);
	}
	return fcells;
}

static void finalize_string(storage& st, pointer a) { 
	a->delete_string();
}

static void finalize_port(storage& st, pointer a) { 
	port* pt = a->portvalue();
	pt->close(port::port_input|port::port_output);
	pt->free();
}

static void finalize_frame(storage& st, pointer a) { 
	a->framevalue()->free();
}

static void finalize_vector(storage& st, pointer a) { 
	a->delete_vector();
}

static storage::final_entry static_finalizers[] = {
	storage::final_entry(cell::T_STRING, finalize_string),
	storage::final_entry(cell::T_PORT, finalize_port),
	storage::final_entry(cell::T_FRAME, finalize_frame),
	storage::final_entry(cell::T_VECTOR, finalize_vector),
	storage::final_entry(0, 0),
};

void storage::finalize_cell(pointer a) {
	if (a->final_fun()) {
		final_map_type::iterator i = final_map_.find(a->type());
		if (i != final_map_.end()) {
			(*((*i).second))(*this, a);
		}
	}
	a->clear();
}

// ========== Routines for Reading ========== 

void port_close(pointer p, int flag) {
	port* pt = p->portvalue();
	pt->closeit(flag);
}

// check if c is in chars 
static bool is_one_of(const char* s, int c) {
     if(c == EOF) return true;
     while (*s) {
          if (*s++ == c)
               return true;
	 }
     return false;
}

bool scheme::file_push(const std::string& fname) {
	FILE* fin = fopen(fname.c_str(), "r");
	if(fin != 0) {
		port* load = load_.push();
		load->set(port::port_input, fin, true);
		loadport_->set_port(load);
	}
	return fin != 0;
}

void scheme::file_pop() {
	nesting_ = load_.curr_nesting();
	if(! load_.base()) {
		port_close(loadport_, port::port_input);
		port* load = load_.pop();
		loadport_->set_port(load);
		if(file_interactive()) {
			putstr(PROMPT);
		}
	}
}

bool scheme::file_interactive() {
	return load_.base() && load_.curr_file() == stdin &&
		inport_->portvalue()->is_file();
}

///////////////////////
///  port creation  ///
///////////////////////

//	make a new port cell from existing port
pointer scheme::mk_port(port* p) {
	return store_.mk_port(p);
}

//	create a port cell from a FILE*
pointer scheme::mk_port(FILE* f, int prop) {
  return store_.mk_port(f, prop, false);
}

//	create a port cell from a filename
pointer scheme::mk_port(const std::string& filename, int prop) {
  return store_.mk_port(filename, prop);
}

//	create a port cel from a strig
pointer scheme::mk_port(char* start, char* past_the_end, int prop) {
  return store_.mk_port(start, past_the_end, prop);
}

void scheme::set_input_port(FILE* fin) {
  inport_ = mk_port(fin, port::port_input);
}

void scheme::set_input_port(char* start, char* past_the_end) {
  inport_ = mk_port(start, past_the_end, port::port_input);
}

void scheme::set_output_port(FILE* fout) {
  outport_ = mk_port(fout, port::port_output);
}

void scheme::set_output_port(char* start, char* past_the_end) {
  outport_ = mk_port(start, past_the_end, port::port_output);
}
//////////////////////////
///  port operations   ///
//////////////////////////

//	Define these so that the constructor can be used to
//	  initialize the port.
void* operator new(size_t, port* location) { return location; }
//	This never gets called, but is required to shut VC++ up.
void operator delete(void *pMem, port *rBuffer) {}

port* port::free_list_ = 0;

//	Allocate a new port or get one from the free list.
//	In either case, initialize it with the constructor.
port* port::make_port(FILE* f, int prop, bool closeit) { 
	if (free_list_ == 0) {
		//	create a new port
		return new port(f, prop, closeit); 
	} else {
		//	get a port off the free list
		port* curr = free_list_;
		free_list_ = curr->link();
		return new (curr) port(f, prop, closeit);
	}
}

port* port::make_port(char* start, char* past_the_end, int prop) { 
	if (free_list_ == 0) {
		//	create a new port
		return new port(start, past_the_end, prop); 
	} else {
		//	get a port off the free list
		port* curr = free_list_;
		free_list_ = curr->link();
		return new (curr) port(start, past_the_end, prop); 
	}
}

static void cleanup_list(port* p) {
	port* next;
	for (port* curr = p; curr != 0; curr = next) {
		next = curr->link();
		delete curr;
	}
}

void port::free() { 
	port* curr = free_list_;
	free_list_ = this;
	kind_ = port_free;
	rep.list.link_ = curr;
}

void port::cleanup() {
	cleanup_list(free_list_);
}

//	get a character from an input port
int port::inchar() {
	if(is_file()) {
		return fgetc(file());
	} else {
		if(*curr() == 0 || curr() == past_the_end()) {
			return EOF;
		} else {
			return *next();
		}
	}
}

//	push a character back into a port
void port::backchar(int c) {
	if(is_file()) {
		ungetc(c, file());
	} else {
		if(curr() != start()) {
			prev();
		}
	}
}

//	put a string to a port
void port::putstr(const char* s) {
	if(is_file()) {
		fputs(s, file());
	} else {
		for( ; *s != 0; s++) {
			if(curr() != past_the_end()) {
				append(*s++);
			}
		}
	}
}

//	put a counted string to a port
void port::putstr(const char* s, int len) {
	if(is_file()) {
		fwrite(s, 1, len, file());
	} else {
		for( ; len != 0; len--) {
			if(curr() != past_the_end()) {
				append(*s++);
			}
		}
	}
}

//	put a character to a port
void port::putcharacter(int c) {
	if(is_file()) {
		fputc(c, file());
	} else {
		if(curr() != past_the_end()) {
			append(c);
		}
	}
}

// cch todo	understand the closeit stuff

//	close a port
void port::closeit(int flag) {
	mask(flag);
	if(!is_io()) {
		if(is_file()) {
			fclose(file());
		}
	}
}

//	close a port
void port::close(int flag) {
	if(is_file() && closeit()) {
		closeit(port::port_input|port::port_output);
	}
}

// get new character from input file 
int scheme::inchar() {
 again:
	port* pt = inport_->portvalue();
	int c = pt->inchar();
	if(c == EOF && inport_ == loadport_ && !load_.base()) {
		file_pop();
		if(nesting_ != 0) {
			return EOF;
		}
		goto again;
	}
	return c;
}

// read chacters up to delimiter, but cater to character constants 
char* scheme::readstr_upto(char* delim) {
	char* p = strbuff_;
	while (!is_one_of(delim, (*p++ = inchar())));
	if(p == strbuff_+2 && p[-2]=='\\') {
		*p=0;
	} else {
		backchar(p[-1]);
		*--p = '\0';
	}
	return strbuff_;
}

// read string expression "xxx...xxx" 
pointer scheme::readstrexp(char* buff, int buf_len) {
	char* p = buff;
	int c1 = 0;
	//	after \ :
	//	  one of the special characters n t r "
	//	  or x followed by two hex chars 
	//	  otherwise next char is itself
	enum { 
		st_ok, 
		st_bsl,		// seen \ 
		st_x1,		// seen \[xX]
		st_x2		// seen \[Xx][0-9a-f] 
	} state = st_ok;
  
	for (;;) {
		int c = inchar();
		if(c == EOF || p-buff > buf_len-1) {
		  return F_;
		}
		switch(state) {
		case st_ok:
			switch(c) {
			case '\\':
				state = st_bsl;
				break;
			case '"':
				*p=0;
				return mk_string(buff, p-buff);
			default:
				*p++=c;
				break;
			}
			break;
		case st_bsl:
			switch(c) {
			case 'x':
			case 'X':
				state = st_x1;
				c1 = 0;
				break;
			case 'n':
				*p++ = '\n';
				state = st_ok;
				break;
			case 't':
				*p++ = '\t';
				state = st_ok;
				break;
			case 'r':
				*p++ = '\r';
				state = st_ok;
				break;
			case '"':
				*p++ = '"';
				state = st_ok;
				break;
			default:
				*p++ = c;
				state = st_ok;
				break;
			}
			break;
		case st_x1:
		case st_x2: {
			c = toupper(c);
			if(c >= '0' && c <= 'F') {
				if(c <= '9') {
				  c1 = (c1<<4)+c-'0';
				} else {
				  c1 = (c1<<4)+c-'A'+10;
				}
				if(state == st_x1) {
				  state = st_x2;
				} else {
				  *p++ = c1;
				  state=st_ok;
				}
			} else {
				return F_;
			}
			break;
			}
		}
	}
}

// skip white characters 
void scheme::skipspace() {
     int c;
     while (isspace(c = inchar()))
          ;
     if(c != EOF) {
          backchar(c);
     }
}

// get token 
int scheme::token() {
     skipspace();
	 int c = inchar();
     switch (c) {
     case EOF:
          return TOK_EOF;
     case '(':
          return TOK_LPAREN;
     case ')':
          return TOK_RPAREN;
     case '.':
          c = inchar();
          if(is_one_of(" \n\t", c)) {
				return TOK_DOT;
          } else {
				backchar(c);
				backchar('.');
				return TOK_ATOM;
          }
     case '\'':
          return TOK_QUOTE;
     case ';':
          return TOK_COMMENT;
     case '"':
          return TOK_DQUOTE;
     case BACKQUOTE:
          return TOK_BQUOTE;
     case ',':
		  c = inchar();
          if (c == '@')
               return TOK_ATMARK;
          else {
               backchar(c);
               return TOK_COMMA;
          }
     case '#':
          c = inchar();
          if (c == '(') {
               return TOK_VEC;
          } else if(c == '!') {
               return TOK_COMMENT;
          } else {
               backchar(c);
               if(is_one_of(" tfodxb\\", c)) {
                    return TOK_SHARP_CONST;
               } else {
                    return TOK_SHARP;
               }
          }
     default:
          backchar(c);
          return TOK_ATOM;
     }
}

//////////////////
///  printing  ///
//////////////////

static char hex_char(int d) {
	if(d < 10) {
		return d+'0';
	} else {
		return d-10+'A';
	}
}

// print all kinds of atoms 
void scheme::printatom(pointer l, bool f) {
	const char* p = atom2str(l, f);
	int len = strlen(p);
	putchars(p, len);
	if (len > 0) record_output();
}

//	Print p enclosed in double quotes
//    using backslash escapes as necessary.
char* make_slashstring(const std::string& p, int len, char* buf) {
	int pos = 0;
	buf[pos++] = '"';
	for (int i = 0; i < len; i++) {
		unsigned char c = p.at(i);
		if(c == 0xff || c == '"' || c < ' ' || c == '\\') {
			buf[pos++] = '\\';
			switch(c) {
			case '"':
			case '\\':
				buf[pos++] = c;
				break;
			case '\n':
				buf[pos++] = 'n';
				break;
			case '\t':
				buf[pos++] = 't';
				break;
			case '\r':
				buf[pos++] = 'r';
				break;
			default:
				buf[pos++] = 'x';
				buf[pos++] = hex_char(c / 16);
				buf[pos++] = hex_char(c % 16);
		  }
		} else {
			buf[pos++] = c;
		}
	}
	buf[pos++] = '"';
	buf[pos++] = 0;
	return buf;
}

//	These print depending on the cell type flag.
static const char* print_atom(scheme&, pointer p, bool, int, char*) {
	if (p == scheme::NIL_) return "()";
	if (p == scheme::T_) return "#t";
	if (p == scheme::F_) return "#f";
	if (p == scheme::VOID_) return "";
	if (p == scheme::EOF_OBJ_) return "#<EOF>";
	return "Unknown atom"; 
}

static const char* print_port(scheme&, pointer, bool, int, char*) { 
	return "#<PORT>"; 
}

static const char* print_integer(scheme&, pointer p, bool, int radix, char* buf) { 
	_itoa(ivalue(p), buf, radix);
	return buf; 
}

static const char* print_real(scheme&, pointer p, bool, int, char* buf) { 
	_gcvt(rvalue(p), 17, buf);
	return buf; 
}

static const char* print_string(scheme&, pointer p, bool f, int, char* buf) { 
	if (!f) {
		return strvalue(p).c_str();
	} else { 
		return make_slashstring(strvalue(p), strlength(p), buf);
	}
}

static const char* print_character(scheme&, pointer p, bool f, int, char* buf) { 
	char c = charvalue(p);
	if (!f) {
	   buf[0] = c;
	   buf[1] = 0;
	   return buf;
	} 
	switch(c) {
	case ' ':
		return "#\\space"; 
	case '\n':
		return "#\\newline"; 
	case '\r':
		return "#\\return"; 
	case '\t':
		return "#\\tab"; 
	default:
#if USE_ASCII_NAMES
		if(c == 127) {
			 return "#\\del"; 
		} 
		if(c < 32) {
			 strcpy(buf, "#\\"); 
			 strcat(buf, charnames[c]); 
			 return buf;
		}
#else
		if(c < 32) {
		  sprintf(buf, "#\\x%x", c); 
		  return buf;
		}
#endif
		sprintf(buf, "#\\%c", c); 
		return buf;
	}
}

static const char* print_symbol(scheme&, pointer p, bool, int, char*) { 
	return symname(p).c_str(); 
}

static const char* print_proc(scheme& sc, pointer p, bool, int, char* buf) { 
	sprintf(buf, "%s", sc.procname(p));
	return buf; 
}

static const char* print_macro(scheme&, pointer, bool, int, char*) { 
	return "#<MACRO>"; 
}

static const char* print_closure(scheme&, pointer, bool, int, char*) { 
	return "#<CLOSURE>"; 
}

static const char* print_promise(scheme&, pointer, bool, int, char*) { 
	return "#<PROMISE>"; 
}

static const char* print_foreign(scheme&, pointer p, bool, int, char* buf) { 
	sprintf(buf, "#<FOREIGN PROCEDURE %ld>", procvalue(p).as_long());
	return buf; 
}

static const char* print_continuation(scheme&, pointer, bool, int, char*) { 
	return "#<CONTINUATION>"; 
}

static const char* print_error(scheme&, pointer, bool, int, char*) { 
	return "#<ERROR>"; 
}

//	This table associates the cell type with the print function.
static scheme::print_entry static_print[] = {
	scheme::print_entry(cell::T_ATOM, print_atom),
	scheme::print_entry(cell::T_PORT, print_port),
	scheme::print_entry(cell::T_INTEGER, print_integer),
	scheme::print_entry(cell::T_REAL, print_real),
	scheme::print_entry(cell::T_STRING, print_string),
	scheme::print_entry(cell::T_CHARACTER, print_character),
	scheme::print_entry(cell::T_SYMBOL, print_symbol),
	scheme::print_entry(cell::T_PROC, print_proc),
	scheme::print_entry(cell::T_MACRO, print_macro),
	scheme::print_entry(cell::T_CLOSURE, print_closure),
	scheme::print_entry(cell::T_PROMISE, print_promise),
	scheme::print_entry(cell::T_FOREIGN, print_foreign),
	scheme::print_entry(cell::T_CONTINUATION, print_continuation),
	scheme::print_entry(0, 0),
};

//	Adds the information in a print table to the interpreter.
void scheme::extend_print(print_entry* print_table) {
	for (int i = 0; print_table[i].test_ != 0; i++) {
		print_entry& pe = print_table[i];
		print_map_.insert(print_map_key(pe.test_, pe.print_));
	}
}

const char* scheme::atom2str(pointer l, bool f, int radix) {
	print_map_type::iterator i = print_map_.find(l->type());
	if (i != print_map_.end()) {
		return (*((*i).second))(*this, l, f, radix, strbuff_);
	}
    return print_error(*this, l, f, radix, strbuff_);
}

pointer scheme::list_star(pointer d) {
	if(cdr(d) == NIL_) {
		return car(d);
	}
	pointer p = cons(car(d),cdr(d));
	pointer q = p;
	while(cdr(cdr(p)) != NIL_) {
		d = cons(car(p), cdr(p));
		if(cdr(cdr(p)) != NIL_) {
			p = cdr(d);
		}
	}
	set_cdr(p, car(cdr(p)));
	return q;
}

// reverse list -- produce new list 
pointer scheme::reverse(pointer a) {
	// a must be checked by gc 
     pointer p;	
     for (p = NIL_; is_pair(a); a = cdr(a)) {
          p = cons(car(a), p);
     }
     return p;
}

// reverse list --- in-place 
pointer reverse_in_place(pointer list, pointer term) {
     pointer p = list;
	 pointer result = (term == 0) ? scheme::NIL_ : term;

     while (p != scheme::NIL_) {
          pointer q = cdr(p);
          set_cdr(p, result);
          result = p;
          p = q;
     }
     return result;
}

// append list -- produce new list 
pointer scheme::append(pointer a, pointer b) {
     pointer p = b;
     if (a != NIL_) {
          a = reverse(a);
          while (a != NIL_) {
               pointer q = cdr(a);
               set_cdr(a, p);
               p = a;
               a = q;
          }
     }
     return p;
}

//	Returns the number of cells in a list.
//	If the first element is not a pair, or if it does not 
//	  end in a pair, then returns -1;
static int list_length(pointer a) {
     register int v = 0;
     register pointer x;	 
     for (x = a; is_pair(x); x = cdr(x)) {
          v++;
     }
     if(x == scheme::NIL_) {
          return v;
     }
     return -1;
}

//	equality of pointers
static bool eq(pointer a, pointer b) {
	return a == b;
}

// equivalence of atoms 
static bool eqv(pointer a, pointer b) {
     if (is_string(a)) {
          if (is_string(b))
               return a == b;
          else
               return false;
     } else if (is_number(a)) {
          if (is_number(b))
               return num_eq(num(a), num(b));
          else
               return false;
     } else if (is_character(a)) {
          if (is_character(b))
               return charvalue(a) == charvalue(b);
          else
               return false;
     } else if (is_port(a)) {
          if (is_port(b))
               return a == b;
          else
               return false;
     } else if (is_proc(a)) {
          if (is_proc(b))
               return a->procvalue() == b->procvalue();
          else
               return false;
     } else {
          return (a == b);
     }
}

static bool equal(pointer a, pointer b);	// forward

static bool vec_equal(pointer a, pointer b) {
	if (!is_vector(a) || !is_vector(b)) return false;
	int len = veclen(a);
	if (len != veclen(b)) return false;
	for (int i = 0; i < len; i++) {
		if (! equal(a->vector_elem(i), b->vector_elem(i))) return false;
	}
	return true;
}

static bool equal(pointer a, pointer b) {
	if (is_pair(a)) {
		if (is_pair(b)) {
			if (! equal(car(a), car(b))) return false;
			return equal(cdr(a), cdr(b));
		} else
			return false;
	} else if (is_vector(a)) {
		if (is_vector(b))
			return vec_equal(a, b);
		else
			return false;
	} else if(is_string(a)) {
		if (is_string(b))
			return strvalue(a) == strvalue(b);
		else
			return false;
	} else {
		return eqv(a, b);
	}
}

/////////////////////
///  environment  ///
/////////////////////

//	A new frame is either (1) a hash table, or (2) an empty list.
pointer environment::new_frame(pointer old_env) {
	if (old_env == scheme::NIL_) { 
		//	Create a vector for the base environment.
		return store_->mk_vector(base_frame_size); 
	} else { 
		//	Create a simple list for the others.
		return scheme::NIL_; 
	} 
}


//	A variable/value pair is consed together and (1) inserted into the
//	  hash table, or (2) linked up to the head of the list.
pointer environment::define_new_binding(pointer variable, pointer value, pointer env) { 
	pointer slot = store_->immutable_cons(variable, value); 
	if (is_vector(car(env))) { 
		int location = hash_fn(symname(variable), veclen(car(env))); 
		set_vector_elem(car(env), location, 
						store_->immutable_cons(slot, vector_elem(car(env), location))); 
	} else { 
		set_car(env, store_->immutable_cons(slot, car(env))); 
	}
	return slot;
} 

//	The base of the environment is (1) the list at the hash table location, 
//	  or (2) the car of the environment.
pointer environment::base(pointer env, pointer variable) {
    if (is_vector(car(env))) { 
		int location = hash_fn(symname(variable), veclen(car(env))); 
		return vector_elem(car(env), location); 
    } else { 
		return car(env); 
    } 
}

//	Make a new frame in the environment.
//	This maintains a list of environments starting at the current
//	  one, which is in envir_.
//	The only way environments are "poppped" is for the
//	  environment to be set to some previous value.
//	Environments pointers are stored in the dump stack.
void environment::push(pointer old_env) { 
	if (old_env == 0) old_env = envir_;
	envir_ = store_->mk_environment(new_frame(old_env), old_env);
} 

//	Search the top environment for a slot matching variable.
//	A slot consists of a (variable . value) pair.
//	This returns the pair where the variable matches the one given.
//	If no such variable is found, returns 0.
pointer environment::search_top_binding(pointer env, pointer variable) { 
	for (pointer y = base(env, variable); y != scheme::NIL_; y = cdr(y)) {
		pointer z = car(y);
		if (car(z) == variable) {
		   return z;		// found it
		}
	}
	return 0;
} 

//	Search the given environment for a slot matching variable.
//	A slot consists of a (variable . value) pair.
//	This returns the slot where the variable matches the one given.
//	If no such variable is found, returns 0.
pointer environment::search_binding(pointer env, pointer variable) { 
	//	search down environment list frame by frame
    for (pointer x = env; x != scheme::NIL_; x = cdr(x)) {
		 // search down slots in this frame
		for (pointer y = base(x, variable); y != scheme::NIL_; y = cdr(y)) {
 			pointer z = car(y);
			if (car(z) == variable) {
				return z;		// found it
			}
         }
    }
    return 0;				// did not find it
}


//	Create a (variable pointer) pair and store in the environment
//	  in the current frame.
//	Return the slot used to store.
pointer environment::define_new_binding(pointer variable, pointer value) { 
	return define_new_binding(variable, value, envir_); 
} 

//	Look up the variable in the given environment and return
//	  its value.
//	If the variable is not found, return 0.
pointer environment::get_binding(pointer variable, pointer envir) const {
	pointer slot = search_binding(envir, variable);
	if (slot) return cdr(slot);
	else return 0;
}

//	Look up the variable in the current environment and return
//	  its value.
//	If the variable is not found, return 0.
pointer environment::get_binding(pointer variable) const {
	return get_binding(variable, envir_);
}

//	Store the variable with the given value in the given environment.
//	If top, then consider only the top environment.
//	If add, then add the binding if it is not already present.
//	Returns true if the varible was stored.
bool environment::define_binding(pointer variable, pointer value, pointer envir, 
						bool top, bool add) {
	pointer slot = (top) ? search_top_binding(envir, variable) : 
		search_binding(envir, variable);;
	if (slot != 0) {
        slot->set_cdr(value);
	} else {
		if (!add) return false;
		define_new_binding(variable, value); 
	}
	return true;
}

//	Store the variable with the given value in the current environment.
bool environment::define_binding(pointer variable, pointer value,
						bool top, bool add) {
	return define_binding(variable, value, envir_, top, add);
}


////////////////////
///  Dump Stack  ///
////////////////////

////////////////////
///  fast stack  ///
////////////////////

dump_stack_fast::~dump_stack_fast() {
	cleanup_list(dump_);
}

void dump_stack_fast::reset() { 
	free_list(dump_);
	dump_ = 0; 
}

//	Mark everything reachable from the stack.
void dump_stack_fast::mark(storage& st) { 
	//	mark everything in the stack
	for (stack_frame* curr = dump_; curr != 0; curr = curr->link()) {
		curr->mark(st);
	}
} 

//	Fast stack is a list of frames.
//	Allocate new frames as needed.
//	Store released frames on a free list.
void dump_stack_fast::push_frame(Op op, 
							pointer args, 
							pointer envir, 
							pointer code) {
	//	get a frame and store registers in it
	dump_ = stack_frame::make_frame(op, args, envir, code, dump_);
}

//	Pop by discarding top element.
bool dump_stack_fast::pop_frame(Op& op, pointer& args, pointer& envir, pointer& code) {
	if (dump_ == 0) return true;
	//	restore registers
	dump_->get(op, args, envir, code);
	//	pop top frame off stack
	stack_frame* curr = dump_;
	dump_ = dump_->link();
	curr->free();
	return false;
}

////////////////////
///  safe stack  ///
////////////////////

//	When switching from fast to safe, we need to copy the fast
//	  stack to a safe stack.
//	This function knows that the fast stack can be traversed
//	  starting at top() by calling next() until it returns 0.
void dump_stack_safe::copy(dump_stack_fast& orig) {
	reset();
	for (stack_frame* curr = orig.top(); 
		 curr != 0; 
		 curr = orig.next(curr)) {
		Op op;
		pointer args;
		pointer envir;
		pointer code;
		curr->get(op, args, envir, code);
		push_frame(op, args, envir, code);
	}
	dump_ = reverse_in_place(dump_);
}

bool dump_stack_safe::is_reset() const { 
	return dump_ == scheme::NIL_; 
}

void dump_stack_safe::reset() { 
	dump_ = scheme::NIL_; 
}

void dump_stack_safe::mark(storage& st) { 
	st.mark(dump_); 
}

void dump_stack_safe::push_frame(Op op, pointer args, pointer envir, pointer code) {
	dump_ = store_->cons(
		store_->mk_frame(op, args, envir, code), 
		dump_);
}

bool dump_stack_safe::pop_frame(Op& op, pointer& args, pointer& envir, pointer& code) {
    if(dump_ == scheme::NIL_) return true; 
	stack_frame* frame = car(dump_)->framevalue();
	frame->get(op, args, envir, code);
	dump_ = cdr(dump_);
    return false; 
}

void scheme::s_return_action(pointer a) { 
    value_ = a;
	pointer envir;
	bool err =  dump_.pop_frame(op_, args_, envir, code_);
	if (err) stop_execution();
	env_.set_curr(envir);
} 

void scheme::s_save(Op op, pointer args, pointer code) { 
	dump_.push_frame(op, args, env_.curr(), code);
} 


//////////////////////////
///  evaluation cycle  ///
//////////////////////////



void scheme::s_error_action(const char* s, pointer a) {
#if USE_ERROR_HOOK
     pointer hook = env_.get_binding(ERROR_HOOK_);
     if (hook) {
         if(a != 0) {
               code_ = cons(cons(QUOTE_, cons(a, NIL_)), NIL_);
         } else {
               code_ = NIL_;
         }
         code_ = cons(mk_string(s), code_);
         set_immutable(car(code_));
         code_ = cons(hook, code_);
         s_goto(Op::EVAL);
    }
#endif

    if(a != 0) {
          args_ = cons(a, NIL_);
    } else {
          args_ = NIL_;
    }
    args_ = cons(mk_string(s), args_);
    set_immutable(car(args_));
    s_goto(Op::ERR0);
}

//////////////////////////
///                    ///
///  opcode execution  ///
///                    ///
//////////////////////////

//	load filename
void scheme::op_load(scheme& sc) {
	scheme_string& fn = scheme_string();
	if (is_string(arg0())) {
		fn = strvalue(arg0());
	} else if (is_symbol(arg0())) {
		fn = symname(arg0());
	} else {
		s_error0("load: arg must be string or symbol");
	}
	if(file_interactive()) {
		sprintf(strbuff_, "Loading %s\n", fn.c_str());
		putstr(strbuff_);
	}
	if (!file_push(fn)) {
	   s_error1("unable to open", arg0());
	}
	s_goto(Op::T0LVL);
}

//	top-level
void scheme::op_t0lvl(scheme& sc) {
	if(file_interactive() && interactive_out_) {
	   putstr("\n");
	}
	interactive_out_ = false;
	nesting_  = 0;
	dump_.reset(); 
	env_.set_curr(global_env_);
	save_inport_ = inport_;
	inport_ = loadport_;
	s_save(Op::T0LVL);
	s_save(Op::VALUEPRINT);
	s_save(Op::T1LVL);
	if (file_interactive()) {
	  putstr(PROMPT);
	}
	s_goto(Op::READ_INTERNAL);
}

//	top level
void scheme::op_t1lvl(scheme& sc) {
	code_ = value_;
	inport_ = save_inport_;
	s_goto(Op::EVAL);
}

void scheme::op_read_internal(scheme& sc) {
	tok_ = token();
	if(tok_ ==TOK_EOF) {
	   if(inport_ == loadport_) {
			args_ = NIL_;
			s_goto(Op::QUIT);
	   } else {
			s_return(EOF_OBJ_);
	   }
	}
	s_goto(Op::RDSEXPR);
}

//	(gensym)
void scheme::op_gensym(scheme& sc) {
	s_return(gensym());
}

//	print evaluation result
void scheme::op_valueprint(scheme& sc) {
	  /* Op::VALUEPRINT is always pushed, because when changing from
		 non-interactive to interactive mode, it needs to be
		 already on the stack */
	if((tracing_ & trace_eval) != 0) {
		putstr("\nGives: ");
	}
	if(file_interactive()) {
		print_flag_ = true;
		args_ = value_;
		s_goto(Op::P0LIST);
	} else {
		s_return(value_);
	}

}

//	tracing part of eval
void scheme::op_eval(scheme& sc) {
#if USE_TRACING
	if((tracing_ & trace_eval) != 0) {
		//s_save(Op::VALUEPRINT);
		s_save(Op::REAL_EVAL, args_, code_);
		args_ = code_;
		putstr("\nEval: ");
		s_goto(Op::P0LIST);
	} else {
		s_goto(Op::REAL_EVAL);
	}
#else
	s_goto(Op::REAL_EVAL);
#endif
}

//	the real eval opcode
void scheme::op_real_eval(scheme& sc) {
	if (is_symbol(code_)) {    // symbol 
		pointer hook = env_.get_binding(code_);
		if (hook) {
			s_return(hook);
		} else {
			s_error1("eval: unbound variable:", code_);
		}
	} else if (is_pair(code_)) {
		pointer x = arg0_code();
		if (is_syntax(x)) {     
			// SYNTAX 
			code_ = code_tail();
			s_goto(syntaxnum(x));
		} else {
			// first, eval top element and eval arguments 
			s_save(Op::E0ARGS, NIL_, code_);
			// If no macros => s_save(Op::E1ARGS, NIL_, code_tail());
			code_ = arg0_code();
			s_goto(Op::EVAL);
		}
	} else {
	   s_return(code_);
	}
}

//`evaluate arguments
void scheme::op_e0args(scheme& sc) {
	if (is_macro(value_)) {    // macro expansion 
		s_save(Op::DOMACRO);
		args_ = cons(code_, NIL_);
		code_ = value_;
		s_goto(Op::APPLY);
	} else {
		code_ = code_tail();
		s_goto(Op::E1ARGS);
	}
}

void scheme::op_e1args(scheme& sc) {
	args_ = cons(value_, args_);
	if (is_pair(code_)) { // continue 
	   s_save(Op::E1ARGS, args_, code_tail());
	   code_ = arg0_code();
	   args_ = NIL_;
	   s_goto(Op::EVAL);
	} else {  // end 
	   args_ = reverse_in_place(args_);
	   code_ = arg0();
	   args_ = arg_tail();
	   s_goto(Op::APPLY);
	}

}

void scheme::op_apply(scheme& sc) {
#if USE_TRACING
	if((tracing_ & trace_eval) != 0) {
		s_save(Op::REAL_APPLY, args_, code_);
		print_flag_ = true;
		//	 args_ = cons(code_, args_);
		putstr("\nApply ");  
		if (is_proc(code_)) putstr(procname(code_)); 
		putstr(" to: ");
		s_goto(Op::P0LIST);
	} else {
		s_goto(Op::REAL_APPLY);
	}
#else
	s_goto(Op::REAL_APPLY);
#endif
}

void scheme::op_real_apply(scheme& sc) {
	if (is_proc(code_)) {
	   s_goto(procvalue(code_));   // PROCEDURE 
	} else if (is_foreign(code_)) {
	   s_return(code_->call_foreign(*this, args_));
	} else if (is_closure(code_) || is_macro(code_) || is_promise(code_)) { // CLOSURE 
	   // Should not accept promise 
	   // make environment 
	   env_.push(closure_env(code_));
	   pointer x, y;
	   for (x = car(closure_code(code_)), y = args_;
			is_pair(x); x = cdr(x), y = cdr(y)) {
			if (y == NIL_) {
				 s_error0("not enough arguments");
			} else {
                 env_.define_new_binding(car(x), car(y)); 
			}
	   }
	   if (x == NIL_) {
			/*--
			 * if (y != NIL_) {
			 *   s_error0("too many arguments");
			 * }
			 */
	   } else if (is_symbol(x))
            env_.define_new_binding(x, y); 
	   else {
			s_error1("syntax error in closure: not a symbol", x);
	   }
	   code_ = cdr(closure_code(code_));
	   args_ = NIL_;
	   s_goto(Op::BEGIN);
	} else if (is_continuation(code_)) { // CONTINUATION 
		dump_.set_curr(continuation(code_));
		s_return(args_ != NIL_ ? car(args_) : NIL_);
	} else {
	   s_error0("illegal function");
	}
}

#if USE_TRACING
//	(tracing flag)
void scheme::op_tracing(scheme& sc) {
	int flag = tracing_;
	tracing_ = ivalue(arg0());
	s_return(mk_integer(flag));
}
#endif

void scheme::op_domacro(scheme& sc) {
	code_ = value_;
	s_goto(Op::EVAL);
}

void scheme::op_lambda(scheme& sc) {
	s_return(mk_closure(code_, env_.curr()));
}

//	(make-closure)
//	(make-closure envir)
void scheme::op_mkclosure(scheme& sc) {
	pointer x = arg0();
	if(car(x) == LAMBDA_) {
		x = cdr(x);
	}
	pointer y = (arg1_nil()) ? env_.curr() : arg1();
	s_return(mk_closure(x, y));
}

//	(quote <datum>)
void scheme::op_quote(scheme& sc) {
	s_return(arg0_code());
}

void scheme::op_def0(scheme& sc) {
	pointer x;
	if (is_pair(car(code_))) {
		x = caar(code_);
		code_ = cons(LAMBDA_, cons(cdar(code_), cdr(code_)));
	} else {
		x = arg0_code();
		code_ = arg1_code();
	}
	if (!is_symbol(x)) {
	   s_error0("variable is not a symbol");
	}
	s_save(Op::DEF1, NIL_, x);
	s_goto(Op::EVAL);
}

void scheme::op_def1(scheme& sc) {
	env_.define_binding(code_, value_, true, true);
	s_return(VOID_);
}

void scheme::op_defp(scheme& sc) {
	pointer x = (arg1_nil()) ? env_.curr() : arg1();
	s_retbool(env_.get_binding(arg0(), x) != 0);
}

void scheme::op_begin(scheme& sc) {
	if (!is_pair(code_)) {
	   s_return(code_);
	}
	if (code_tail() != NIL_) {
	   s_save(Op::BEGIN, NIL_, code_tail());
	}
	code_ = arg0_code();
	s_goto(Op::EVAL);
}

void scheme::op_if0(scheme& sc) {
	s_save(Op::IF1, NIL_, code_tail());
	code_ = arg0_code();
	s_goto(Op::EVAL);
}

void scheme::op_if1(scheme& sc) {
	if (is_true(value_))
	   code_ = arg0_code();
	else if (code_tail() == NIL_)
		code_ = VOID_;
	else
	   code_ = arg1_code();  // (if #f 1) ==> () because car(NIL_) = NIL_ 
	s_goto(Op::EVAL);
}

void scheme::op_set0(scheme& sc) {
	s_save(Op::SET1, NIL_, arg0_code());
	code_ = arg1_code();
	s_goto(Op::EVAL);
}

void scheme::op_set1(scheme& sc) {
	bool updated = env_.define_binding(code_, value_, false, false);
	if (! updated) {
	   s_error1("set!: unbound variable:", code_);
	}
	s_return(VOID_);
}

void scheme::op_let0(scheme& sc) {
	args_ = NIL_;
	value_ = code_;
	code_ = is_symbol(arg0_code()) ? arg1_code() : arg0_code();
	s_goto(Op::LET1);
}

void scheme::op_let1(scheme& sc) {
	args_ = cons(value_, args_);
	if (is_pair(code_)) {		// continue 
	   s_save(Op::LET1, args_, code_tail());
	   code_ = cadar(code_);
	   args_ = NIL_;
	   s_goto(Op::EVAL);
	} else {					// end 
	   args_ = reverse_in_place(args_);
	   code_ = arg0();
	   args_ = arg_tail();
	   s_goto(Op::LET2);
	}
}

void scheme::op_let2(scheme& sc) {
    env_.push(); 
	pointer x, y;
	for (x = is_symbol(arg0_code()) ? arg1_code() : arg0_code(), y = args_;
		 y != NIL_; 
		 x = cdr(x), y = cdr(y)) {
			env_.define_new_binding(caar(x), car(y)); 
	}
	if (is_symbol(arg0_code())) {    // named let 
	   for (x = arg1_code(), args_ = NIL_; x != NIL_; x = cdr(x)) {
			args_ = cons(caar(x), args_);
	   }
	   x = mk_closure(cons(reverse_in_place(args_), cddr(code_)), env_.curr());
       env_.define_new_binding(car(code_), x); 
	   code_ = cddr(code_);
	   args_ = NIL_;
	} else {
	   code_ = code_tail();
	   args_ = NIL_;
	}
	s_goto(Op::BEGIN);
}

//	let*
void scheme::op_let0star(scheme& sc) {
	if (car(code_) == NIL_) {
        env_.push(); 
		code_ = code_tail();
		s_goto(Op::BEGIN);
	}
	s_save(Op::LET1STAR, cdr(code_), car(code_));
	code_ = cadaar(code_);
	s_goto(Op::EVAL);
}

//	let* (make new frame)
void scheme::op_let1star(scheme& sc) {
    env_.push(); 
	s_goto(Op::LET2STAR);
}

// let* (caluculate parameters) 
void scheme::op_let2star(scheme& sc) {
    env_.define_new_binding(caar(code_), value_); 
	code_ = code_tail();
	if (is_pair(code_)) { // continue 
	   s_save(Op::LET2STAR, args_, code_);
	   code_ = cadar(code_);
	   args_ = NIL_;
	   s_goto(Op::EVAL);
	} else {  // end 
	   code_ = args_;
	   args_ = NIL_;
	   s_goto(Op::BEGIN);
	}
}

//	letrec
void scheme::op_let0rec(scheme& sc) {
    env_.push(); 
	args_ = NIL_;
	value_ = code_;
	code_ = car(code_);
	s_goto(Op::LET1REC);
}

// letrec (caluculate parameters) 
void scheme::op_let1rec(scheme& sc) {
	args_ = cons(value_, args_);
	if (is_pair(code_)) { // continue 
	   s_save(Op::LET1REC, args_, code_tail());
	   code_ = cadar(code_);
	   args_ = NIL_;
	   s_goto(Op::EVAL);
	} else {  // end 
	   args_ = reverse_in_place(args_);
	   code_ = car(args_);
	   args_ = cdr(args_);
	   s_goto(Op::LET2REC);
	}
}

void scheme::op_let2rec(scheme& sc) {
	pointer x, y;
	for (x = car(code_), y = args_; y != NIL_; x = cdr(x), y = cdr(y)) {
       env_.define_new_binding(caar(x), car(y)); 
	}
	code_ = code_tail();
	args_ = NIL_;
	s_goto(Op::BEGIN);
}

void scheme::op_cond0(scheme& sc) {
	if (!is_pair(code_)) {
	   s_error0("syntax error in cond");
	}
	s_save(Op::COND1, NIL_, code_);
	code_ = caar(code_);
	s_goto(Op::EVAL);
}

void scheme::op_cond1(scheme& sc) {
	if (is_true(value_)) {
	   if ((code_ = cdar(code_)) == NIL_) {
			s_return(value_);
	   }
	   if(car(code_) == FEED_TO_) {
			if(!is_pair(cdr(code_))) {
				 s_error0("syntax error in cond");
			}
			pointer x = cons(QUOTE_, cons(value_, NIL_));
			code_ =cons(cadr(code_),cons(x,NIL_));
			s_goto(Op::EVAL);
	   }
	   s_goto(Op::BEGIN);
	} else {
	   if ((code_ = code_tail()) == NIL_) {
			s_return(VOID_);
	   } else {
			s_save(Op::COND1, NIL_, code_);
			code_ = caar(code_);
			s_goto(Op::EVAL);
	   }
	}
}

void scheme::op_delay(scheme& sc) {
	pointer x = mk_promise(cons(NIL_, code_), env_.curr());
	s_return(x);
}

//	(and <test1> ...)
//	startup
void scheme::op_and0(scheme& sc) {
	if (code_ == NIL_) {
	   s_return(T_);
	}
	s_save(Op::AND1, NIL_, code_tail());
	code_ = arg0_code();
	s_goto(Op::EVAL);
}

//	(and <test1> ...)
//	after each <test> evaluation
void scheme::op_and1(scheme& sc) {
	if (is_false(value_)) {
	   s_return(value_);
	} else if (code_ == NIL_) {
	   s_return(value_);
	} else {
	   s_save(Op::AND1, NIL_, code_tail());
	   code_ = arg0_code();
	   s_goto(Op::EVAL);
	}
}

//	(or <test1> ...)
//	startup
void scheme::op_or0(scheme& sc) {
	if (code_ == NIL_) {
	   s_return(F_);
	}
	s_save(Op::OR1, NIL_, code_tail());
	code_ = arg0_code();
	s_goto(Op::EVAL);
}

//	(or <test1> ...)
//	after each <test> evaluation
void scheme::op_or1(scheme& sc) {
	if (is_true(value_)) {
	   s_return(value_);
	} else if (code_ == NIL_) {
	   s_return(value_);
	} else {
	   s_save(Op::OR1, NIL_, code_tail());
	   code_ = arg0_code();
	   s_goto(Op::EVAL);
	}
}

 // cons-stream 
void scheme::op_c0stream(scheme& sc) {
	s_save(Op::C1STREAM, NIL_, code_tail());
	code_ = car(code_);
	s_goto(Op::EVAL);
}

// cons-stream 
void scheme::op_c1stream(scheme& sc) {
	args_ = value_;  // save value_ to register args_ for gc 
	pointer x = mk_promise(cons(NIL_, code_), env_.curr());
	s_return(cons(args_, x));
}

void scheme::op_macro0(scheme& sc) {
	pointer x;
	if (is_pair(car(code_))) {
	   x = caar(code_);
	   code_ = cons(LAMBDA_, cons(cdar(code_), cdr(code_)));
	} else {
	   x = car(code_);
	   code_ = cadr(code_);
	}
	if (!is_symbol(x)) {
	   s_error0("variable is not a symbol");
	}
	s_save(Op::MACRO1, NIL_, x);
	s_goto(Op::EVAL);
}

void scheme::op_macro1(scheme& sc) {
	value_->set_macro();
	env_.define_binding(code_, value_, true, true);
	s_return(VOID_);
}

void scheme::op_case0(scheme& sc) {
	s_save(Op::CASE1, NIL_, code_tail());
	code_ = car(code_);
	s_goto(Op::EVAL);
}

void scheme::op_case1(scheme& sc) {
	pointer x, y;
	for (x = code_; x != NIL_; x = cdr(x)) {
	   if (!is_pair(y = caar(x))) {
			break;
	   }
	   for ( ; y != NIL_; y = cdr(y)) {
			if (eqv(car(y), value_)) {
				 break;
			}
	   }
	   if (y != NIL_) {
			break;
	   }
	}
	if (x != NIL_) {
	   if (is_pair(caar(x))) {
			code_ = cdar(x);
			s_goto(Op::BEGIN);
	   } else {// else 
			s_save(Op::CASE2, NIL_, cdar(x));
			code_ = caar(x);
			s_goto(Op::EVAL);
	   }
	} else {
	   s_return(VOID_);
	}
}

void scheme::op_case2(scheme& sc) {
	if (is_true(value_)) {
	   s_goto(Op::BEGIN);
	} else {
	   s_return(NIL_);
	}
}

void scheme::op_peval(scheme& sc) {
	if(cdr(args_) != NIL_) {
	   env_.set_curr(arg1());
	}
	code_ = arg0();
	s_goto(Op::EVAL);
}

void scheme::op_papply(scheme& sc) {
	code_ = arg0();
	args_ = list_star(cdr(args_));
	//args_ = cadr(args_);
	s_goto(Op::APPLY);
}

void scheme::op_continuation(scheme& sc) {
	code_ = arg0();
	args_ = cons(mk_continuation(dump_.curr()), NIL_);
	s_goto(Op::APPLY);

}


//	(ineaxct->exact z)
void scheme::op_inex2ex(scheme& sc) {
	pointer z = arg0();
	if(is_integer(z)) {
	   s_return(z);
	} else {
		long v = ivalue(z);
		double dd;
		modf(rvalue(z), &dd);
		if (fabs((double)v - dd) > 1.0) {
			s_error1("inexact->exact: not convertable:", z);
		} else {
			s_return(mk_integer(v));
		}
	}
}

//	(exact->inexact z)
void scheme::op_ex2inex(scheme& sc) {
	pointer z = arg0();
	if(is_real(z)) {
	   s_return(z);
	} else {
	   s_return(mk_real(ivalue(z)));
	}
}

//	(exp z)
void scheme::op_exp(scheme& sc) {
	s_return(mk_real(exp(rvalue(arg0()))));
}

//	(log z)
void scheme::op_log(scheme& sc) {
	s_return(mk_real(log(rvalue(arg0()))));
}

//	(sin z)
void scheme::op_sin(scheme& sc) {
	s_return(mk_real(sin(rvalue(arg0()))));
}

//	(cos z)
void scheme::op_cos(scheme& sc) {
	s_return(mk_real(cos(rvalue(arg0()))));
}

//	(tan z)
void scheme::op_tan(scheme& sc) {
	s_return(mk_real(tan(rvalue(arg0()))));
}

//	(asin z)
void scheme::op_asin(scheme& sc) {
	s_return(mk_real(asin(rvalue(arg0()))));
}

//	(acos z)
void scheme::op_acos(scheme& sc) {
	s_return(mk_real(acos(rvalue(arg0()))));
}

//	(atan z)
//	(atan y x)
void scheme::op_atan(scheme& sc) {
	if(arg1_nil()) {
	   s_return(mk_real(atan(rvalue(arg0()))));
	} else {
	   s_return(mk_real(atan2(rvalue(arg0()), rvalue(arg1()))));
	}
}

//	(sqrt z)
void scheme::op_sqrt(scheme& sc) {
	s_return(mk_real(sqrt(rvalue(arg0()))));
}

//	(expt z1 z2)
void scheme::op_expt(scheme& sc) {
	s_return(mk_real(pow(rvalue(arg0()),rvalue(arg1()))));
}

//	(floor x)
void scheme::op_floor(scheme& sc) {
	s_return(mk_real(floor(rvalue(arg0()))));
}

//	(ceiling x)
void scheme::op_ceiling(scheme& sc) {
	s_return(mk_real(ceil(rvalue(arg0()))));
}

//	(truncate x)
void scheme::op_truncate(scheme& sc) {
	double rvalue_of_x = rvalue(arg0()) ;
	if (rvalue_of_x > 0) {
		s_return(mk_real(floor(rvalue_of_x)));
	} else {
		s_return(mk_real(ceil(rvalue_of_x)));
	}
}

//	(round x)
void scheme::op_round(scheme& sc) {
	pointer x = arg0();
	s_return(mk_real(round_per_R5RS(rvalue(x))));
}

//	(exact? z)
void scheme::op_exact(scheme& sc) {
	s_retbool(is_integer(arg0()));
}

//	(inexact? z)
void scheme::op_inexact(scheme& sc) {
	s_retbool(is_real(arg0()));
}

//	(odd? n)
void scheme::op_odd(scheme& sc) {
	s_retbool(ivalue(arg0()) % 2 != 0);
}

//	(even? n)
void scheme::op_even(scheme& sc) {
	s_retbool(ivalue(arg0()) % 2 == 0);
}

//	(zero? z)
void scheme::op_zero(scheme& sc) {
	s_retbool(num(arg0()).zero());
}

//	(positive? z)
void scheme::op_positive(scheme& sc) {
	s_retbool(num(arg0()).positive());
}

//	(negative? z)
void scheme::op_negative(scheme& sc) {
	s_retbool(num(arg0()).negative());
}

//	(+ z1 ...)
void scheme::op_add(scheme& sc) {
	num v = num_zero_;
	for (pointer x = args_; x != NIL_; x = cdr(x)) {
		v = v.add(num(car(x)));
	}
	s_return(mk_number(v));
}

//	(- z)
//	(- z1 z2)
//	(- z1 z2 ...)
void scheme::op_sub(scheme& sc) {
	pointer x;
	num v;
	if(arg1_nil()) {
		x = args_;
		v = num_zero_;
	} else {
		x = cdr(args_);
		v = num(arg0());
	}
	for (; x != NIL_; x = cdr(x)) {
		v = v.sub(num(car(x)));
	}
	s_return(mk_number(v));
}

//	(* z1 ...)
void scheme::op_mul(scheme& sc) {
	num v = num_one_;
	for (pointer x = args_; x != NIL_; x = cdr(x)) {
		v = v.mul(num(car(x)));
	}
	s_return(mk_number(v));
}

//	(/ z)
//	(/ z1 z2)
//	(/ z1 z2 ...)
void scheme::op_div(scheme& sc) {
	pointer x;
	num v;
	if(arg1_nil()) {
		x = args_;
		v = num_one_;
	} else {
		x = cdr(args_);
		v = num(arg0());
	}
	for (; x != NIL_; x = cdr(x)) {
		if (is_zero_double(rvalue(car(x)))) {
			s_error0("/: division by zero");
		} else {
			v = v.div(num(car(x)));
		}
	}
	s_return(mk_number(v));
}

//	(quotient n)
//	(quotient n1 n2)
//	(quotient n1 n2 ...)
void scheme::op_intdiv(scheme& sc) {
	pointer x;
	num v;
	if(arg1_nil()) {
	   x = args_;
	   v = num_one_;
	} else {
	   x = cdr(args_);
	   v = num(arg0());
	}
	for (; x != NIL_; x = cdr(x)) {
		if (ivalue(car(x)) == 0) {
			s_error0("quotient: division by zero");
		} else {
			v = v.intdiv(num(car(x)));
		}
	}
	s_return(mk_number(v));
}

//	(remainder n1 n2)
void scheme::op_rem(scheme& sc) {
	num v(arg0());
	if (ivalue(arg1()) == 0) {
	   s_error0("remainder: division by zero");
	} else {
	   v = v.rem(num(arg1()));
	}
	s_return(mk_number(v));
}

//	(modulo n1 n2)
void scheme::op_mod(scheme& sc) {
	num v(arg0());
	if (ivalue(arg1()) == 0) {
	   s_error0("modulo: division by zero");
	} else {
	   v = v.mod(num(arg1()));
	}
	s_return(mk_number(v));
}

//	(max x1 x2 ...)
void scheme::op_max(scheme& sc) {
	num v(arg0());
	for (pointer x = cdr(args_); x != NIL_; x = cdr(x)) {
		v = v.max(num(car(x)));
	}
	s_return(mk_number(v));
}

//	(min x1 x2 ...)
void scheme::op_min(scheme& sc) {
	num v(arg0());
	for (pointer x = cdr(args_); x != NIL_; x = cdr(x)) {
		v = v.min(num(car(x)));
	}
	s_return(mk_number(v));
}

//	(abs x)
void scheme::op_abs(scheme& sc) {
	num v(arg0());
	s_return(mk_number(v.abs()));
}

//	gcd of two positive args
static long gcd(long a, long b) {
	while (b != 0) {
		long temp = b;
		b = a % b;
		a = temp;
	}
	return a;
}

//	lcm of two positive args
static long lcm(long a, long b) {
	return (a/gcd(a, b)) * b;
}

//	(gcd n1 ...)
void scheme::op_gcd(scheme& sc) {
	long a = 0;
	for (pointer x = args_; x != NIL_; x = cdr(x)) {
		long b = abs(ivalue(car(x)));
		a = gcd(a, b);
	}
	s_return(mk_integer(a));
}

//	(lcm n1 ...)
void scheme::op_lcm(scheme& sc) {
	long a = 1;
	for (pointer x = args_; x != NIL_; x = cdr(x)) {
		long b = abs(ivalue(car(x)));
		a = lcm(a, b);
	}
	s_return(mk_integer(a));
}

//	(car pair)
void scheme::op_car(scheme& sc) {
	s_return(car(arg0()));
}

//	(cdr pair)
void scheme::op_cdr(scheme& sc) {
	s_return(cdr(arg0()));
}

//	(cons obj1 obj2)
void scheme::op_cons(scheme& sc) {
	set_cdr(args_, arg1());
	s_return(args_);
}

//	(set-car! pair obj)
void scheme::op_setcar(scheme& sc) {
	if(is_immutable(arg0())) {
		s_error0("set-car!: unable to alter immutable pair");
	} else {
		set_car(arg0(), arg1());
		s_return(VOID_);
	}
}

//	(set-cdr! pair obj)
void scheme::op_setcdr(scheme& sc) {
	if(is_immutable(arg0())) {
		s_error0("set-cdr!: unable to alter immutable pair");
	} else {
		set_cdr(arg0(), arg1());
		s_return(VOID_);
	}
}

//	(char->integer char)
void scheme::op_char2int(scheme& sc) {
	s_return(mk_integer((char)charvalue(arg0())));
}

//	(integer->char n)
void scheme::op_int2char(scheme& sc) {
	s_return(mk_character((char)ivalue(arg0())));
}

//	(char upcase char)
void scheme::op_charupcase(scheme& sc) {
	s_return(mk_character(toupper((char)charvalue(arg0()))));
}

//	(char-downcase char)
void scheme::op_chardncase(scheme& sc) {
	s_return(mk_character(tolower((char)charvalue(arg0()))));
}

//	(symbol->string symbol)
void scheme::op_sym2str(scheme& sc) {
	pointer x = mk_string(symname(arg0()));
	set_immutable(x);
	s_return(x);
}

//	(atom->string atom)
void scheme::op_atom2str(scheme& sc) {
	pointer a = arg0();
	if(is_number(a) || is_character(a) || is_string(a) || is_symbol(a)) {
		const char* p = atom2str(a, false);
		s_return(mk_string(p));
	} else {
		s_error1("atom->string: not an atom:", a);
	}
}

//	(string->symbol string)
void scheme::op_str2sym(scheme& sc) {
	s_return(define_symbol(strvalue(arg0())));
}

//	(string->atom string) 
void scheme::op_str2atom(scheme& sc) {
	const std::string& str = strvalue(arg0());
	if(str.at(0) == '#') {
		s_return(mk_sharp_const(str.c_str()+1));
	} else {
		s_return(mk_atom(str));
	}
}

//	(make-string k)
//	(make-string k char)
void scheme::op_mkstring(scheme& sc) {
	int len = ivalue(arg0());
	char fill = (arg1_nil()) ? ' ': charvalue(arg1());
	s_return(mk_string(len, fill));
}

//	(string-length string)
void scheme::op_strlen(scheme& sc) {
	s_return(mk_integer(strlength(arg0())));
}

//	(string-ref string k)
void scheme::op_strref(scheme& sc) {
	const std::string& str = strvalue(arg0());
	int k = ivalue(arg1());

	if(k >= strlength(arg0())) {
	   s_error1("string-ref: out of bounds:", arg1());
	}
	s_return(mk_character(str.at(k)));
}

//	(string-set! string k char)
void scheme::op_strset(scheme& sc) {
	if(is_immutable(arg0())) {
		s_error1("string-set!: unable to alter immutable string:", arg0());
	}

	pointer x = arg0();
	int k = ivalue(arg1());
	if(k >= strlength(x)) {
		s_error1("string-set!: out of bounds:", arg1());
	}

	x->set_string_elem(k, charvalue(arg2()));
	s_return(VOID_);
}

//	(substring string start)
//	(substring string start end)
void scheme::op_substr(scheme& sc) {
	int start = ivalue(arg1());

	if(start > strlength(arg0())) {
	   s_error1("substring: start out of bounds:", arg1());
	}

	int end;
	if(cddr(args_) != NIL_) {
	   end = ivalue(arg2());
	   if(end > strlength(arg0()) || end < start) {
			s_error1("substring: end out of bounds:", arg2());
	   }
	} else {
	   end = strlength(arg0());
	   //	need not check for end < start -- is must be OK
	}

	const char* str = strvalue(arg0()).c_str();
	pointer x = mk_string(str+start, end - start);

	s_return(x);
}

//	(string-append string ...)
void scheme::op_strappend(scheme& sc) {
	int len = list_length(args_);
	if(len < 0) {
	   s_error1("strappend: not a proper list:", args_);
	}
	pointer x = mk_string(0);
	std::string& str = strvalue(x);
	for (pointer y = args_; is_pair(y); y = cdr(y)) {
		str.append(strvalue(car(y)));
	}
	s_return(x);
}

//	(string->list string)
void scheme::op_str2list(scheme& sc) {
	int slen = strlength(arg0());
	const std::string& str = strvalue(arg0());
	pointer x = NIL_;
	for (int i = 0; i < slen; i++) {
		x = cons(mk_character(str.at(i)), x);
	}
	s_return(reverse_in_place(x));
}

//	(string char ...)
void scheme::op_string(scheme& sc) {
	int len = list_length(args_);
	pointer s = mk_string(len);
	std::string& str = strvalue(s);
	int i = 0;
	for (pointer y = args_; is_pair(y); y = cdr(y)) {
		if (!is_character(car(y))) {
			s_error0("string: not character");
		}
//		str.append(1, charvalue(car(y)));
		str.at(i++) = charvalue(car(y));
	}
	s_return(s);
}

//	(list->string chars)
void scheme::op_list2str(scheme& sc) {
	int len = list_length(arg0());
	if(len < 0) {
	   s_error1("list2str: not a proper list:", arg0());
	}
	pointer x = mk_string(len);
	std::string& str = strvalue(x);
	int i = 0;
	for (pointer y = arg0(); is_pair(y); y = cdr(y)) {
		if (!is_character(car(y))) {
			s_error0("list2str: not character");
		}
		str.append(1, charvalue(car(y)));
	}
	s_return(x);
}

//	(string->number string)
//	(string->number string radix)
void scheme::op_str2num(scheme& sc) {
	int radix = (arg_tail() != NIL_) ? ivalue(arg1()): 10;
	if (radix != 2 && radix != 8 && radix != 10 && radix != 16) {
		s_error1("str2num: illegal radix: ", arg1());
	}
	s_return(mk_atom(strvalue(arg0()), true, radix));
}

//	(number->string number)
//	(number->string number radix)
void scheme::op_num2str(scheme& sc) {
	pointer number = arg0();
	int radix = 10;
	if (arg_tail() != NIL_) {
		radix = ivalue(arg1());
	}
	if(is_number(number)) {
		const char* p = atom2str(number, false, radix);
		s_return(mk_string(p));
	} else {
		s_error1("num->string: not a number:", number);
	}
}

//	(string-fill! string char)
void scheme::op_strfill(scheme& sc) {
	pointer s = arg0();
	int slen = strlength(s);
	std::string& str = strvalue(s);
	int fill = charvalue(arg1());
	str.append(slen, fill);
	s_return(s);
}

//	(list obj ...)
void scheme::op_list(scheme& sc) {
	s_return(args_);
}

//	(list-tail list k)
void scheme::op_listtail(scheme& sc) {
	int k = ivalue(arg1());
	if (k <0) s_error1("list-ref: too short:", arg1());
	for (pointer x = arg0(); x != NIL_; x = cdr(x)) {
		if (k == 0) s_return(x);
		k--;
	}
	//	It is unclear how this should behave in these cases.
	//	Should it return NIL or give an error ?
	//	Which cases are considered an error ?
	if (k == 0) s_return(NIL_);
	s_error1("list-tail: too short:", arg1());
}

//	(list-ref list k)
void scheme::op_listref(scheme& sc) {
	int k = ivalue(arg1());
	if (k <0) s_error1("list-ref: too short:", arg1());
	for (pointer x = arg0(); x != NIL_; x = cdr(x)) {
		if (k == 0) s_return(car(x));
		k--;
	}
	s_error1("list-ref: too short:", arg1());
}

typedef bool (*test_predicate2)(pointer a, pointer b);
//	member test using given compare function
pointer memtest(pointer a, pointer b, test_predicate2 tst) {
	for (pointer x = b; x != scheme::NIL_; x = cdr(x)) {
		if ((*tst)(car(x), a)) return x;;
	}
	return scheme::F_;
}

//	(memq obj list)
void scheme::op_memq(scheme& sc) {
	s_return(memtest(arg0(), arg1(), eq));
}

//	(memv obj list)
void scheme::op_memv(scheme& sc) {
	s_return(memtest(arg0(), arg1(), eqv));
}

//	(member obj list)
void scheme::op_member(scheme& sc) {
	s_return(memtest(arg0(), arg1(), equal));
}

//	assoc test using given compare function
pointer asstest(pointer a, pointer b, test_predicate2 tst) {
	pointer x = a;
    pointer y;
	for (y = b; is_pair(y); y = cdr(y)) {
	   if (!is_pair(car(y))) {
			return 0;
	   }
	   if (tst(x, caar(y)))
			break;
	}
	if (is_pair(y)) {
	   return car(y);
	} else {
		return scheme::F_;
	}
}


//	(assq obj alist)
void scheme::op_assq(scheme& sc) {
	pointer x = asstest(arg0(), arg1(), eq);
	if (x == 0) 
		s_error0("unable to handle non pair element");
	s_return(x);
}

//	(assv obj alist)
void scheme::op_assv(scheme& sc) {
	pointer x = asstest(arg0(), arg1(), eqv);
	if (x == 0) 
		s_error0("unable to handle non pair element");
	s_return(x);
}

//	(assoc obj alist)
void scheme::op_assoc(scheme& sc) {
	pointer x = asstest(arg0(), arg1(), equal);
	if (x == 0) 
		s_error0("unable to handle non pair element");
	s_return(x);
}


//	(vector obj ...)
void scheme::op_vector(scheme& sc) {
	int len = list_length(args_);
	if(len < 0) {
	   s_error1("vector: not a proper list:", args_);
	}
	pointer vec = mk_vector(len);
	int i = 0;
	for (pointer x = args_; is_pair(x); x = cdr(x)) {
	   set_vector_elem(vec, i++, car(x));
	}
	s_return(vec);
}

//	(make-vector k)
//	(make-vector k fill)
void scheme::op_mkvector(scheme& sc) {
	pointer vec = mk_vector(ivalue(arg0()));

	if(!arg1_nil()) {
	   vec->fill_vector(arg1());
	}
	s_return(vec);
}

//	(vector-length vector)
void scheme::op_veclen(scheme& sc) {
	s_return(mk_integer(veclen(arg0())));
}

//	(vector-ref vector k)
void scheme::op_vecref(scheme& sc) {
	int k = ivalue(arg1());
	if(k >= veclen(arg0())) {
	   s_error1("vector-ref: out of bounds:", arg1());
	}

	s_return(vector_elem(arg0(), k));
}

//	(vector-set! vector k obj)
void scheme::op_vecset(scheme& sc) {
	if(is_immutable(arg0())) {
	   s_error1("vector-set!: unable to alter immutable vector:", arg0());
	}

	int k = ivalue(arg1());
	if(k >= veclen(arg0())) {
	   s_error1("vector-set!: out of bounds:", arg1());
	}

	set_vector_elem(arg0(), k, arg2());
	s_return(VOID_);
}

//	(list->vector list)
void scheme::op_list2vec(scheme& sc) {
	int len = list_length(arg0());
	if(len < 0) {
	   s_error1("vector->list: not a proper list:", arg0());
	}
	pointer vec = mk_vector(len);
	int i = 0;
	for (pointer x = arg0(); is_pair(x); x = cdr(x)) {
	   set_vector_elem(vec, i++, car(x));
	}
	s_return(vec);
}

//	(vector-fill! vector fill)
void scheme::op_vecfill(scheme& sc) {
	int len = veclen(arg0());
	for (int i = 0; i < len; i++) {
		set_vector_elem(arg0(), i, arg1());
	}
	s_return(VOID_);
}

//	(vector->list vector)
void scheme::op_vec2list(scheme& sc) {
	int len = veclen(arg0());
	pointer x = NIL_;
	for (int i = 0; i < len; i++) {
		x = cons(vector_elem(arg0(), i), x);
	}
	s_return(reverse_in_place(x));
}

//	(not obj)
void scheme::op_not(scheme& sc) {
	s_retbool(is_false(arg0()));
}

//	(boolean? obj)
void scheme::op_boolp(scheme& sc) {
	s_retbool(arg0() == F_ || arg0() == T_);
}

//	(oef-object? obj)
void scheme::op_eofobjp(scheme& sc) {
	s_retbool(arg0() == EOF_OBJ_);
}

//	(null? obj)
void scheme::op_nullp(scheme& sc) {
	s_retbool(arg0() == NIL_);
}

//	(= z1 z2 z3 ...)
//	(< z1 z2 z3 ...)
//	(> z1 z2 z3 ...)
//	(<= z1 z2 z3 ...)
//	(>= z1 z2 z3 ...)
void scheme::op_comp(scheme& sc) {
	comp_func cf = get_comp(op_);
	pointer x = args_;
	num z(car(x));
	x = cdr(x);

	for (; x != NIL_; x = cdr(x)) {
	   if(!cf(z, num(car(x)))) {
			s_retbool(false);
	   }
	   z = num(car(x));
	}
	s_retbool(true);
}

//	(symbol? obj)
void scheme::op_symbolp(scheme& sc) {
	s_retbool(is_symbol(arg0()));
}

//	(number? obj)
void scheme::op_numberp(scheme& sc) {
	s_retbool(is_number(arg0()));
}

//	(string? obj)
void scheme::op_stringp(scheme& sc) {
	s_retbool(is_string(arg0()));
}

//	(integer? obj)
void scheme::op_integerp(scheme& sc) {
	s_retbool(is_integer(arg0()));
}

//	(real? obj)
void scheme::op_realp(scheme& sc) {
	s_retbool(is_number(arg0()));
}

//	(char? obj)
void scheme::op_charp(scheme& sc) {
	s_retbool(is_character(arg0()));
}

#if USE_CHAR_CLASSIFIERS
//	(char-alphabetic? char)
void scheme::op_charap(scheme& sc) {
	s_retbool(Cisalpha(charvalue(arg0())));
}

//	(char-numeric? char)
void scheme::op_charnp(scheme& sc) {
	s_retbool(Cisdigit(charvalue(arg0())));
}

//	(char-whitespace? char)
void scheme::op_charwp(scheme& sc) {
	s_retbool(Cisspace(charvalue(arg0())));
}

//	(char-upper-case? letter)
void scheme::op_charup(scheme& sc) {
	s_retbool(Cisupper(charvalue(arg0())));
}

//	(char-lower-case? letter)
void scheme::op_charlp(scheme& sc) {
	s_retbool(Cislower(charvalue(arg0())));
}

#endif

// port? 
void scheme::op_portp(scheme& sc) {
	s_retbool(is_port(arg0()));
}

//	(input-port? obj)
void scheme::op_inportp(scheme& sc) {
	s_retbool(is_inport(arg0()));
}

//	(output-port? obj)
void scheme::op_outportp(scheme& sc) {
	s_retbool(is_outport(arg0()));
}

//	(procedure? obj)
void scheme::op_procp(scheme& sc) {
	 /*--
	  * continuation should be procedure by the example
	  * (call-with-current-continuation procedure?) ==> #t
		 * in R^3 report sec. 6.9
	  */
	s_retbool(is_proc(arg0()) || is_closure(arg0())
		 || is_continuation(arg0()) || is_foreign(arg0()));
}

//	(pair? obj)
void scheme::op_pairp(scheme& sc) {
	s_retbool(is_pair(arg0()));
}

//	(list? obj)
void scheme::op_listp(scheme& sc) {
	pointer slow, fast;
	slow = fast = arg0();
	while (1) {
		if (!is_pair(fast)) s_retbool(fast == NIL_);
		fast = cdr(fast);
		if (!is_pair(fast)) s_retbool(fast == NIL_);
		fast = cdr(fast);
		slow = cdr(slow);
		if (fast == slow) {
		  // The fast pointer has looped back around and caught up
		  //   with the slow pointer, hence the structure is circular,
		  //   not of finite length, and therefore not a list
		  s_retbool(false);
		}
	}
}

//	(environment? obj)
void scheme::op_envp(scheme& sc) {
	s_retbool(is_environment(arg0()));
}

//	(vector? obj)
void scheme::op_vectorp(scheme& sc) {
	s_retbool(is_vector(arg0()));
}

//	(eq? obj1 obj2)
void scheme::op_eq(scheme& sc) {
	s_retbool(arg0() == arg1());
}

//	(eqv? obj1 obj2)
void scheme::op_eqv(scheme& sc) {
	s_retbool(eqv(arg0(), arg1()));
}

//	(equal? obj1 obj2)
void scheme::op_equal(scheme& sc) {
	s_retbool(equal(arg0(), arg1()));
}

//	(vector-equal? obj1 obj2)
void scheme::op_vecequal(scheme& sc) {
	s_retbool(vec_equal(arg0(), arg1()));
}

//	(force promise)
void scheme::op_force(scheme& sc) {
	code_ = arg0();
	if (is_promise(code_)) {
	   // Should change type to closure here 
	   s_save(Op::SAVE_FORCED, NIL_, code_);
	   args_ = NIL_;
	   s_goto(Op::APPLY);
	} else {
	   s_return(code_);
	}
}

// Save forced value replacing promise 
void scheme::op_save_forced(scheme& sc) {
	if (is_promise(code_)) {
		*code_ = *value_;
		s_return(value_);
	} else {
	   s_return(code_);
	}
}

void scheme::write_common() {
	if(is_pair(cdr(args_))) {
	   if(arg1() != outport_) {
			pointer x = cons(outport_, NIL_);
			s_save(Op::SET_OUTPORT, x);
			outport_ = arg1();
	   }
	}
	args_ = arg0();
}

//	(write obj)
//	(write obj port)
void scheme::op_write(scheme& sc) {
	write_common();
	print_flag_ = true;
	s_goto(Op::P0LIST);
}

//	(wite-char char)
//	(write-char char port)
void scheme::op_write_char(scheme& sc) {
	write_common();
	print_flag_ = false;
	s_goto(Op::P0LIST);
}

//	(display obj)
//	(display obj port)
void scheme::op_display(scheme& sc) {
	write_common();
	print_flag_ = false;
	s_goto(Op::P0LIST);
}

//	(newline)
//  (newline port)
void scheme::op_newline(scheme& sc) {
	if(is_pair(args_)) {
	   if(arg0() != outport_) {
			pointer x = cons(outport_, NIL_);
			s_save(Op::SET_OUTPORT, x);
			outport_ = arg0();
	   }
	}
	putstr("\n");
	record_output();
	s_return(VOID_);
}

// error 
void scheme::op_err0(scheme& sc) {
	retcode_ = -1;
	if (!is_string(arg0())) {
	   args_ = cons(mk_string(" -- "),args_);
	   set_immutable(arg0());
	}
	putstr("s_error: ");
	putstr(strvalue(arg0()).c_str());
	args_ = cdr(args_);
	s_goto(Op::ERR1);
}

// error 
void scheme::op_err1(scheme& sc) {
	putstr(" ");
	if (args_ != NIL_) {
	   s_save(Op::ERR1, cdr(args_));
	   args_ = arg0();
	   print_flag_ = true;
	   s_goto(Op::P0LIST);
	} else {
	   putstr("\n");
	   if(interactive_repl_) {
			s_goto(Op::T0LVL);
	   } else {
			stop_execution();
			return;
	   }
	}
}

//	(reverse list)
void scheme::op_reverse(scheme& sc) {
	s_return(reverse(arg0()));
}

// list* 
void scheme::op_list_star(scheme& sc) {
	s_return(list_star(args_));
}

//	(append list)
void scheme::op_append(scheme& sc) {
	if(args_ == NIL_) {
	   s_return(NIL_);
	}
	pointer x = arg0();
	if(arg1_nil()) {
		s_return(args_);
	}
	for (pointer y = cdr(args_); y != NIL_; y = cdr(y)) {
	   x = append(x, car(y));
	}
	s_return(x);
}

// put 
void scheme::op_put(scheme& sc) {
	if (!hasprop(arg0()) || !hasprop(arg1())) {
	   s_error0("illegal use of put");
	}
	pointer x, y;
	for (x = symprop(arg0()), y = arg1(); x != NIL_; x = cdr(x)) {
	   if (caar(x) == y) {
			break;
	   }
	}
	if (x != NIL_) {
	   set_cdr(car(x), arg2());
	} else {
	   set_symprop(arg0(), cons(cons(y, arg2()),
						symprop(arg0())));
	}
	s_return(VOID_);

}

// get 
void scheme::op_get(scheme& sc) {
	if (!hasprop(arg0()) || !hasprop(arg1())) {
	   s_error0("illegal use of get");
	}
	pointer x, y;
	for (x = symprop(arg0()), y = arg1(); x != NIL_; x = cdr(x)) {
	   if (caar(x) == y) {
			break;
	   }
	}
	if (x != NIL_) {
	   s_return(cdar(x));
	} else {
	   s_return(NIL_);
	}
}

//	(quit)
//	(quit retcode)
void scheme::op_quit(scheme& sc) {
	if(is_pair(args_)) {
	   retcode_ = ivalue(arg0());
	}
	stop_execution();
}

//	(gc) 
void scheme::op_gc(scheme& sc) {
	store_.gc();
	s_return(T_);
}

//	(gc-verbose flag)
void scheme::op_gcverb(scheme& sc) {
	bool was = store_.gc_verbose();
	store_.set_gc_verbose(arg0() != F_);
	s_retbool(was);
}

//	(new-segment)
//	(new-segment n)
//	Allocates n (default to 1) new cell segments.
void scheme::op_newsegment(scheme& sc) {
	int n;
	if (!is_pair(args_)) {
		n = store_.alloc_cells(1);
	} else if (!is_number(arg0())) {
	   s_error0("new-segment: argument must be a number");
	} else {
		n = store_.alloc_cells(ivalue(arg0()));
	}
	s_return(mk_integer(n));
}

//	(oblist)
//	Returns list of all symbols.
void scheme::op_oblist(scheme& sc) {
	s_return(oblist_.all_symbols());
}

//	(current-input-port) 
void scheme::op_curr_inport(scheme& sc) {
	s_return(inport_);
}

//	(current-output-port) 
void scheme::op_curr_outport(scheme& sc) {
	s_return(outport_);
}

//	(open-input-file filename)
void scheme::op_open_infile(scheme& sc) {
	pointer p = mk_port(strvalue(arg0()), port::port_input);
	if (p == NIL_) { s_return(F_); }
	else { s_return(p); }
}

//	(open-output-file filename)
void scheme::op_open_outfile(scheme& sc) {
	pointer p = mk_port(strvalue(arg0()), port::port_output);
	if (p == NIL_) { s_return(F_); }
	else { s_return(p); }
}

//	(open-input-output-file  filename)
void scheme::op_open_inoutfile(scheme& sc) {
	pointer p = mk_port(strvalue(arg0()), port::port_input|port::port_output);
	if (p == NIL_) { s_return(F_); }
	else { s_return(p); }
}

//	(open-input-string string string)
void scheme::op_open_instring(scheme& sc) {
	const std::string& str = strvalue(arg0());
	int len = strlength(arg0());
	pointer p = mk_port(str, port::port_input);
	s_retbool(p != NIL_);
}

//	(open-output-string string)
void scheme::op_open_outstring(scheme& sc) {
	const std::string& str = strvalue(arg0());
	int len = strlength(arg0());
	pointer p = mk_port(str, port::port_output);
	s_retbool(p != NIL_);
}

//	(open-input-output-string string)
void scheme::op_open_inoutstring(scheme& sc) {
	const std::string& str = strvalue(arg0());
	int len = strlength(arg0());
	pointer p = mk_port(str, port::port_input|port::port_output);
	s_retbool(p != NIL_);
}

//	(close-input-port port)
void scheme::op_close_inport(scheme& sc) {
	port_close(arg0(), port::port_input);
	s_return(VOID_);
}

//	(close-output-port port)
void scheme::op_close_outport(scheme& sc) {
	port_close(arg0(), port::port_output);
	s_return(VOID_);
}

//	(interaction-environment)
void scheme::op_int_env(scheme& sc) {
	s_return(global_env_);
}

//	(current-environment) 
void scheme::op_curr_env(scheme& sc) {
	s_return(env_.curr());
}

bool scheme::check_nesting() {
	if(nesting_ != 0) {
		int n = nesting_;
		nesting_ = 0;
		retcode_ = -1;
		s_error_action("unmatched parentheses:", mk_integer(n));
		return true;
	}
	return false;
}

//	(read)
//	(read port)
void scheme::op_read(scheme& sc) {
	if (check_nesting()) return;
	if(!is_pair(args_)) {
	   s_goto(Op::READ_INTERNAL);
	}
	if(!is_inport(arg0())) {
	   s_error1("read: not an input port:", arg0());
	}
	if(arg0() == inport_) {
	   s_goto(Op::READ_INTERNAL);
	}
	pointer x = inport_;
	inport_ = arg0();
	x = cons(x, NIL_);
	s_save(Op::SET_INPORT, x);
	s_goto(Op::READ_INTERNAL);
}

//	(read-char)
//	(read-char port)
//	(peek-char)
//	(peek-char port)
void scheme::op_read_char(scheme& sc) {
	if (check_nesting()) return;
	if(is_pair(args_)) {
	   if(arg0() != inport_) {
			pointer x = inport_;
			x = cons(x, NIL_);
			s_save(Op::SET_INPORT, x);
			inport_ = arg0();
	   }
	}
	int c = inchar();
	if(c == EOF) {
	   s_return(EOF_OBJ_);
	}
	if(op_ == Op::PEEK_CHAR) {
	   backchar(c);
	}
	s_return(mk_character(c));
}

//	(char-ready?)
//	(char-ready? port)
void scheme::op_char_ready(scheme& sc) {
	if (check_nesting()) return;
	pointer p = (is_pair(args_)) ? arg0(): inport_;
	int res = p->portvalue()->is_string();
	s_retbool(res != 0);
}

//	(set-input-port! port)
void scheme::op_set_inport(scheme& sc) {
	if (check_nesting()) return;
	inport_ = arg0();
	s_return(value_);
}

//	(set-output-port! port)
void scheme::op_set_outport(scheme& sc) {
	if (check_nesting()) return;
	outport_ = arg0();
	s_return(value_);
}

//	primitive
//	read s-expression
void scheme::op_rdsexpr(scheme& sc) {
	if (check_nesting()) return;
	pointer x;
	switch (tok_) {
	case TOK_EOF:
	   if(inport_ == loadport_) {
			args_ = NIL_;
			s_goto(Op::QUIT);
	   } else {
			s_return(EOF_OBJ_);
	   }
	case TOK_COMMENT: {
	   int c;
	   while ((c=inchar()) != '\n' && c!=EOF)
			;
	   tok_ = token();
	   s_goto(Op::RDSEXPR);
	}
	case TOK_VEC:
	   s_save(Op::RDVEC);
	   // fallthrough 
	case TOK_LPAREN:
	   tok_ = token();
	   if (tok_ == TOK_RPAREN) {
			s_return(NIL_);
	   } else if (tok_ == TOK_DOT) {
			s_error0("syntax error: illegal dot expression");
	   } else {
			load_.nesting_inc();
			s_save(Op::RDLIST);
			s_goto(Op::RDSEXPR);
	   }
	case TOK_QUOTE:
	   s_save(Op::RDQUOTE);
	   tok_ = token();
	   s_goto(Op::RDSEXPR);
	case TOK_BQUOTE:
	   tok_ = token();
	   if(tok_ ==TOK_VEC) {
			s_save(Op::RDQQUOTEVEC);
			tok_ =TOK_LPAREN;
			s_goto(Op::RDSEXPR);
		} else {
			s_save(Op::RDQQUOTE);
		}
	   s_goto(Op::RDSEXPR);
	case TOK_COMMA:
	   s_save(Op::RDUNQUOTE);
	   tok_ = token();
	   s_goto(Op::RDSEXPR);
	case TOK_ATMARK:
	   s_save(Op::RDUQTSP);
	   tok_ = token();
	   s_goto(Op::RDSEXPR);
	case TOK_ATOM:
	   s_return(mk_atom(readstr_upto("();\t\n\r ")));
	case TOK_DQUOTE:
	   x = readstrexp(strbuff_, sizeof(strbuff_));
	   if(x == F_) {
			s_error0("error reading string");
	   }
	   set_immutable(x);
	   s_return(x);
	case TOK_SHARP: {
		   pointer hook = env_.get_binding(SHARP_HOOK_);
		   if(hook) {
				code_ = cons(hook, NIL_); 
				s_goto(Op::EVAL);
		   } else {
				s_error0("undefined sharp expression");
		   }
		}
	case TOK_SHARP_CONST:
	   if ((x = mk_sharp_const(readstr_upto("();\t\n\r "))) == NIL_) {
			s_error0("undefined sharp expression");
	   } else {
			s_return(x);
	   }
	default:
	   s_error0("syntax error: illegal token");
	}
}

//	primitive
//	read list
void scheme::op_rdlist(scheme& sc) {
	if (check_nesting()) return;
	args_ = cons(value_, args_);
	tok_ = token();
	if (tok_ == TOK_COMMENT) {
	   int c;
	   while ((c=inchar()) != '\n' && c!=EOF)
			;
	   tok_ = token();
	}
	if (tok_ == TOK_RPAREN) {
	   load_.nesting_dec();		// do this before inchar in case it is EOF
	   int c = inchar();
	   if (c != '\n') backchar(c);
	   s_return(reverse_in_place(args_));
	} else if (tok_ == TOK_DOT) {
	   s_save(Op::RDDOT, args_);
	   tok_ = token();
	   s_goto(Op::RDSEXPR);
	} else {
	   s_save(Op::RDLIST, args_);;
	   s_goto(Op::RDSEXPR);
	}
}

//	primitive
void scheme::op_rddot(scheme& sc) {
	if (check_nesting()) return;
	if (token() != TOK_RPAREN) {
	   s_error0("syntax error: illegal dot expression");
	} else {
	   load_.nesting_dec();
	   s_return(reverse_in_place(args_, value_));
	}
}

//	primitive
void scheme::op_rdquote(scheme& sc) {
	if (check_nesting()) return;
	s_return(cons(QUOTE_, cons(value_, NIL_)));
}

//	primitive
void scheme::op_rdqquote(scheme& sc) {
	if (check_nesting()) return;
	s_return(cons(QQUOTE_, cons(value_, NIL_)));
}

//	primitive
void scheme::op_rdqquotevec(scheme& sc) {
	if (check_nesting()) return;
	s_return(cons(define_symbol("apply"),
		cons(define_symbol("vector"), 
			 cons(cons(QQUOTE_, 
			  cons(value_,NIL_)),
			  NIL_))));
}

//	primitive
void scheme::op_rdunquote(scheme& sc) {
	if (check_nesting()) return;
	s_return(cons(UNQUOTE_, cons(value_, NIL_)));
}

//	primitive
void scheme::op_rduqtsp(scheme& sc) {
	if (check_nesting()) return;
	s_return(cons(UNQUOTESP_, cons(value_, NIL_)));
}

//	primitive
void scheme::op_rdvec(scheme& sc) {
	if (check_nesting()) return;
	/*code_ = cons(mk_proc(Op::VECTOR),value_);
	s_goto(Op::EVAL); Cannot be quoted*/
	/*x = cons(mk_proc(Op::VECTOR),value_);
	s_return(x); Cannot be part of pairs*/
	/*code_ = mk_proc(Op::VECTOR);
	args_ = value_;
	s_goto(Op::APPLY);*/
	args_ =value_;
	s_goto(Op::VECTOR);
}

//	primitive
//	Print a list stored in args_
void scheme::op_p0list(scheme& sc) {
	if (check_nesting()) return;
	if(is_vector(args_)) {
	   putstr("#(");
	   record_output();
	   args_ = cons(args_,mk_integer(0));
	   s_goto(Op::PVECFROM);
	} else if(is_environment(args_)) {
	   putstr("#<ENVIRONMENT>");
	   record_output();
	   s_return(VOID_);
	} else if (!is_pair(args_)) {
	   printatom(args_, print_flag_);
	   s_return(VOID_);
	} else if (arg0() == QUOTE_ && ok_abbrev(arg_tail())) {
	   putstr("'");
	   args_ = arg1();
	   s_goto(Op::P0LIST);
	} else if (arg0() == QQUOTE_ && ok_abbrev(arg_tail())) {
	   putstr("`");
	   args_ = arg1();
	   s_goto(Op::P0LIST);
	} else if (arg0() == UNQUOTE_ && ok_abbrev(arg_tail())) {
	   putstr(", ");
	   args_ = arg1();
	   s_goto(Op::P0LIST);
	} else if (arg0() == UNQUOTESP_ && ok_abbrev(arg_tail())) {
	   putstr(",@");
	   args_ = arg1();
	   s_goto(Op::P0LIST);
	} else {
	   putstr("(");
	   record_output();
	   s_save(Op::P1LIST, arg_tail());
	   args_ = arg0();
	   s_goto(Op::P0LIST);
	}
}

//	primitive
//	Print a list where args_ points to the tail of a list
//	  that has already started printing in p0list.
void scheme::op_p1list(scheme& sc) {
	if (check_nesting()) return;
	if (is_pair(args_)) {
		s_save(Op::P1LIST, arg_tail());
		putstr(" ");
		args_ = arg0();
		s_goto(Op::P0LIST);
	} else if(is_vector(args_)) {
		s_save(Op::P1LIST);
		putstr(" . ");
		s_goto(Op::P0LIST);
	} else {
		if (args_ != NIL_) {
		  putstr(" . ");
		  printatom(args_, print_flag_);
		}
		putstr(")");
		s_return(VOID_);
	}
}

//	primitive
void scheme::op_pvecfrom(scheme& sc) {
	if (check_nesting()) return;
	int i = ivalue(cdr(args_));
	pointer vec = arg0();
	int len = veclen(vec);
	if(i == len) {
	   putstr(")");
	   s_return(VOID_);
	} else {
	   pointer elem = vector_elem(vec, i);
	   set_integer(cdr(args_), i+1);
	   s_save(Op::PVECFROM, args_);
	   args_ = elem;
	   putstr(" ");
	   s_goto(Op::P0LIST);
	}
}

// (length list)
void scheme::op_list_length(scheme& sc) {
	long len = list_length(arg0());
	if(len < 0) {
	   s_error1("length: not a list:", arg0());
	}
	s_return(mk_integer(len));
}

//	(get-closure-code  xx)
void scheme::op_get_closure(scheme& sc) {
	args_ = arg0();
	if (args_ == NIL_) {
	   s_return(F_);
	} else if (is_closure(args_)) {
	   s_return(cons(LAMBDA_, closure_code(value_)));
	} else if (is_macro(args_)) {
	   s_return(cons(LAMBDA_, closure_code(value_)));
	} else {
	   s_return(F_);
	}
}

//	(closure? obj)
void scheme::op_closurep(scheme& sc) {
	/*
	 * Note, macro object is also a closure.
	 * Therefore, (closure? <#MACRO>) ==> #t
	 */
	s_retbool(is_closure(arg0()));
}

//	(macro? obj)
void scheme::op_macrop(scheme& sc) {
	s_retbool(is_macro(arg0()));
}

void scheme::op_time(scheme& sc) {
	time_t t;
	time(&t);
	s_return(mk_integer(t));
}

void scheme::op_illegal(scheme& sc) {
	s_error0("illegal op");
}

//////////////////////////////
///  end of op_ functions  ///
//////////////////////////////

//	This table is used to build the test table.
//	The first element is used as an index in the test table for the
//	  given test. 
//	The next element is the function to call to perform the test.
//	The last element is used to form the error message.
static scheme::test_entry static_tests[] = {
  scheme::test_entry(TST_NONE,		is_any,			0),
  scheme::test_entry(TST_ANY,		is_any,			0),
  scheme::test_entry(TST_STRING,	is_string,		"string"),
  scheme::test_entry(TST_SYMBOL,	is_symbol,		"symbol"),
  scheme::test_entry(TST_PORT,		is_port,		"port"),
  scheme::test_entry(TST_INPORT,	is_inport,		"input port"),
  scheme::test_entry(TST_OUTPORT,	is_outport,		"output_port"),
  scheme::test_entry(TST_ENVIRONMENT, is_environment, "environment"),
  scheme::test_entry(TST_PAIR,		is_pair,		"pair"),
  scheme::test_entry(TST_LIST,		is_pair_or_nil, "pair or '()"),
  scheme::test_entry(TST_CHAR,		is_character,	"character"),
  scheme::test_entry(TST_VECTOR,	is_vector,		"vector"),
  scheme::test_entry(TST_NUMBER,	is_number,		"number"),
  scheme::test_entry(TST_INTEGER,	is_num_integer, "integer"),
  scheme::test_entry(TST_NATURAL,	is_nonneg,		"non-negative integer"),
  scheme::test_entry(-1, 0, 0)
};

//	Build the test table.
bool scheme::extend_tests(test_entry* tests) {
	int i;
	//	Determine the highest test index.
	int max_test = 0;
	for (i = 0; tests[i].test_ != -1; i++) {
		if (tests[i].test_ > max_test) max_test = tests[i].test_; 
	}
	max_test++;
	//	Make sure the test table is big enough.
	tests_.resize(max_test);
	//	Assign the test table entries.
	for (i = 0; tests[i].test_ != -1; i++) {
		int dest = tests[i].test_;
		if (tests_[dest].fct_ != 0) return false;	// err
		tests_[dest] = tests[i];
	}
	return true;
}

//	Tests to make sure the op union overlay is working as expected.
bool Op::check() {
	Op test;
	test.rep.ss.tbl_ = 1;
	test.rep.ss.op_ = 2;
	if (test.rep.l.oper_ == 0x00010002) return true;
	return false;
}

//	The entries in this table are indexed by the Op::OP enum values,
//	  so they must correspond exactly.  This means that 
//	  dispatch_table[Op::LOAD] is the entry for Op::LOAD.
//	This means that not only that the entries must be in the same order as
//	  in the enum OP, but also the #ifdef directives must match as well.
static op_code_ent<scheme> base_op_code_table[] = {
#define PROC_ENTRY0(op, func, name, minargs, maxargs) \
	op_code_ent<scheme>(Op::op, op_code_entry::proc, scheme::func, name, minargs, maxargs, TST_NONE),
#define PROC_ENTRY1(op, func, name, minargs, maxargs, type1) \
	op_code_ent<scheme>(Op::op, op_code_entry::proc, scheme::func, name, minargs, maxargs, type1),
#define PROC_ENTRY2(op, func, name, minargs, maxargs, type1, type2) \
	op_code_ent<scheme>(Op::op, op_code_entry::proc, scheme::func, name, minargs, maxargs, type1, type2),
#define PROC_ENTRY3(op, func, name, minargs, maxargs, type1, type2, type3) \
	op_code_ent<scheme>(Op::op, op_code_entry::proc, scheme::func, name, minargs, maxargs, type1, type2, type3),
#define PRIM_ENTRY(op, func, name) \
	op_code_ent<scheme>(Op::op, op_code_entry::prim, scheme::func, name),
#include "op.h"
};

op_code_table::op_code_table() :
	dispatch_table_(),
	size_(),
	sc_(0) {}

//	All the callback functor were allocated dynamically
//	  by newCallback, so they must be destroyed.
op_code_table::~op_code_table() {
	for (unsigned int i = 0; i < dispatch_table_.size(); i++) {
		for (int j = 0; j < size_[i]; j++) {
			dispatch_table_[i][j].func_->destroy();
		}
		delete [] dispatch_table_[i];
	}
}

//	If the op_code_entry table is out of order, this takes care
//	  of putting it into the right place.
int op_code_table::extend(op_code_entry* base, scheme* sc) {
	int size = get_size(base);
	op_code_entry* ent = new op_code_entry[size];
	int n = dispatch_table_.size();
	dispatch_table_.push_back(ent);
	size_.push_back(size);
	for (int i = 0; i < size; i++) {
		int dest = base[i].op();
		ent[dest] = base[i];
	}
	init_proc(n, size);
	return n;
}

//	Look up the name of a built-in proc
const char* scheme::procname(pointer x) {
	Op op = x->procvalue();
	const std::string& name = op_table_.lookup_name(op);
	if(name.size() == 0) {
		return "ILLEGAL!";
	}
	return name.c_str();
}

//	Given a string cell, look up the opcode that implements the proc.
//	This is similar to a perfect hash.  
//	Used only for syntax.	  
//	Remember to rewrite if more are added! 
static Op syntaxnum(pointer p) {
	 const std::string& s = strvalue(car(p));
     switch(s.size()) {
     case 2:										// xx
          if(s.at(0) == 'i') return Op::IF0;        // if 
          else return Op::OR0;						// or 
		  
     case 3:										// xxx
          if(s.at(0) == 'a') return Op::AND0;		// and 
          else return Op::LET0;						// let 
     case 4:
          switch(s.at(3)) {							// xxxx
          case 'e': return Op::CASE0;				// case 
          case 'd': return Op::COND0;				// cond 
          case '*': return Op::LET0STAR;			// let* 
          default: return Op::SET0;					// set!           
          }
     case 5:
          switch(s.at(2)) {							// xxxxx
          case 'g': return Op::BEGIN;				// begin 
          case 'l': return Op::DELAY;				// delay 
          case 'c': return Op::MACRO0;				// macro 
          default: return Op::QUOTE;				// quote 
          }
     case 6:
          switch(s.at(2)) {							// xxxxxx
          case 'm': return Op::LAMBDA;				// lambda 
          case 'f': return Op::DEF0;				// define 
          default: return Op::LET0REC;				// letrec 
          }
     default:
          return Op::C0STREAM;						// cons-stream 
     }
}

///////////////////////
///  op code table  ///
///////////////////////

//	Test each of the args to a built-in proc to verify that it is
//	  of the correct type.
//	If return is false, args are ok.
//	If true, there was an error and there is a message in errbuf.
//	Arg_tests_encoding is a string of chars describing the
//	  required type of each arg, in order.
bool op_code_entry::test_arg_types(const scheme& sc, pointer arglist, char* errbuf) const {
	const char* t = arg_tests_encoding_;
	if (*t == 0) return false;
	for(int i = 0; arglist != scheme::NIL_; i++) {
		pointer arg = car(arglist);
		int j = (int)*t;
		const scheme::test_entry& test = sc.test(j);
		if(!test.fct_(arg)) {
			sprintf(errbuf, "%s: argument %d must be: %s",
				name_.c_str(),
				i+1,
				test.kind_);
			return true;
		}

		if(*(t+1)) t++;				// next char, unless we have run out
		arglist = cdr(arglist);		//	next arg
	} 
	return false;
}

//	Tests the arg list for all built-in procs.
//	Checks that the right number of args are passed, and if
//	  so, calls test_arg_types to test their types.
bool op_code_entry::test_args(const scheme& sc, pointer args, char* errbuf) const {
	if (type_ != proc) return false;

	// Check number of arguments 
	int n = list_length(args);
	if(n < min_arity_) {
		sprintf(errbuf, "%s: needs%s %d argument(s)",
			name_.c_str(),
			min_arity_ == max_arity_ ? "" : " at least",
			min_arity_);
		return true;
	}
	if(n > max_arity_) {
		sprintf(errbuf, "%s: needs%s %d argument(s)",
			name_.c_str(),
			min_arity_ == max_arity_ ? "" : " at most",
			max_arity_);
		return true;
	}
#if FAST
	if(arg_tests_encoding_[0] != 0) {		// try to avoid a function call
		return test_arg_types(sc, args, errbuf);
	}
#else
	return test_arg_types(args, errbuf);
#endif
	return false;
}

//	For each builtin (entries with proc flag), 
//	  create a proc with that name and put it in the
//	  global environment.
//	Primitives also have names, but these are for trracing
//	  only and are not placed in the environment.
void op_code_table::init_proc(int n_table, int size) {
  for(int i = 0; i < size; i++) {
	  Op op(n_table, i);
	  if(op_code_table::is_proc(op)) {
		sc_->assign_proc(op, lookup_name(op));
    }
  }
}

//	Test args_ to make sure it is compatible with op_.
//	Return true if there is something wrong.
void scheme::test_builtin() {
#if FAST
	const op_code_entry* ent = &op_table_.dispatch_table_[op_.rep.ss.tbl_][op_.rep.ss.op_];
	if (ent->type_ != op_code_entry::proc)		// avoid call to test-args
		return;
#else
	const op_code_entry* ent = op_table_.lookup_op(op_);
	if (ent == 0) {
		s_error_action("unknown op");
		return;
	}
#endif

	if(ent->test_args(*this, args_, strbuff_)) {
		s_error_action(strbuff_);		// sets op_
	}
}

//	Kernel of the intepreter.
int scheme::eval_cycle(Op oper) {
	op_ = oper;
	reset_execution();

#if FAST
    int count;
	for (count = 0; break_; count++) {
#else
	for (int count = 0; continue_execution(); count++) {
#endif
		test_builtin();		// if test fails, sets op_
#if USE_TRACING
		if (tracing_ & trace_op) {
			if (op_table_.is_proc(op_)) putstr("*"); else putstr(" ");
			putstr(op_table_.lookup_name(op_).c_str()); 
			putstr("\n");
		}
#endif
#if FAST
		//	This code is faster, at least in debug mode, because it skips
		//	all the function calls.
		(*op_table_.dispatch_table_[op_.rep.ss.tbl_][op_.rep.ss.op_].func_)(*this);
#else
		//	This code is cleaner
		op_table_.call(op_);
#endif
	}
	//	no memory will cause a break
	if(store_.no_memory()) {
	  fprintf(stderr, "No memory!\n");
	}
	return count;
}

//	These are meant to be called by the application to
//	  initiate some operation in scheme.

//	Run interactive top level.
int scheme::eval_cycle_top_level() {
	dump_.reset();
	env_.set_curr(global_env_);
	retcode_ = 0;
	return eval_cycle(Op::T0LVL);
}

//	Evaluate the symbol at code_.
int scheme::eval_cycle_eval() {
	dump_.reset();
	env_.set_curr(global_env_);
	retcode_ = 0;
	return eval_cycle(Op::EVAL);
}

//	Apply the function at code_ to the args in args_.
int scheme::eval_cycle_apply() {
	dump_.reset();
	env_.set_curr(global_env_);
	retcode_ = 0;
	return eval_cycle(Op::APPLY);
}

// ========== Initialization of internal keywords ========== 

pointer scheme::assign_syntax(const char* name) {
	pointer x = oblist_.add_symbol(name);
	x->set_syntax();
	return x;
}

//	Make a proc cell associating the name with the op
//	  and put it in the global environment.
pointer scheme::assign_proc(Op op, const std::string& name) {
     return env_.define_new_binding(define_symbol(name), mk_proc(op)); 
}

static char* static_syntax_table[] = {
	"lambda",
	"quote",
	"define",
	"if",
	"begin",
	"set!",
	"let",
	"let*",
	"letrec",
	"cond",
	"delay",
	"and",
	"or",
	"cons-stream",
	"macro",
	"case",
	0,
};

void scheme::extend_syntax(char** syntax_table) {
	for (int i = 0; syntax_table[i] != 0; i++) {
		assign_syntax(syntax_table[i]);
	}
}

// initialization of TinyScheme 
scheme::scheme(int first_cell_seg, 
			   int max_cell_seg) :
	store_(max_cell_seg),
	gensym_cnt_(0),
	inport_(NIL_),
	outport_(NIL_),
	save_inport_(NIL_),
	loadport_(NIL_),
	nesting_(0),
	op_table_(),
	interactive_repl_(false),
	args_(NIL_),
	code_(NIL_),
	break_(true),
	tracing_(trace_none),
	env_(&store_),
	oblist_(&store_),
	dump_(&store_) {
	op_table_.init(this);
	//	this should be first, because no cells can be allocated before
	//	it is done
	if (store_.init(*this, first_cell_seg)) return;
	//	  and symbol table
	oblist_.init();

	//	set 0 and 1
	num_zero_.set(0L);
	num_one_.set(1L);

	// init NIL_ 
	NIL_->set_marked_atom(NIL_, NIL_);
	// init T_ 
	T_->set_marked_atom(T_, T_);
	// init F_ 
	F_->set_marked_atom(F_, F_);
	// init VOID_ 
	VOID_->set_marked_atom(VOID_, VOID_);
	// init EOF_OBJ 
	EOF_OBJ_->set_marked_atom(EOF_OBJ_, EOF_OBJ_);
	// init global_env 
	env_.push(NIL_); 
	global_env_ = env_.curr();
	// init else 
	env_.define_new_binding(define_symbol("else"), VOID_); 

	bool ok = extend_tests(static_tests);
	if (! ok) {
		fprintf(stderr, "Static_tests incorrect\n");
		store_.set_no_memory();
		return;
	}

	ok = Op::check();
	if (! ok) {
		fprintf(stderr, "Op union incorrect\n");
		store_.set_no_memory();
		return;
	}

	//	Add the base procs into the oblist.
	op_table_.extend(base_op_code_table, this);

	//	Build marker table
	store_.extend_marker(static_markers);

	//	Build finalizer table.
	store_.extend_finalizer(static_finalizers);

	//	Build print table.
	extend_print(static_print);

	//	Build syntax table.
	extend_syntax(static_syntax_table);

	// intialization of global pointers to special symbols 
	LAMBDA_		= define_symbol("lambda");
	QUOTE_		= define_symbol("quote");
	QQUOTE_		= define_symbol("quasiquote");
	UNQUOTE_	= define_symbol("unquote");
	UNQUOTESP_	= define_symbol("unquote-splicing");
	FEED_TO_	= define_symbol("=>");
	COLON_HOOK_ = define_symbol("*colon-hook*");
	ERROR_HOOK_ = define_symbol("*error-hook*");
	SHARP_HOOK_ = define_symbol("*sharp-hook*");
}

scheme::~scheme() {
	if (store_.no_memory()) return;

	oblist_.clear();
	env_.set_curr(NIL_);
	code_ = NIL_;
	args_ = NIL_;
	value_ = NIL_;
	if(is_port(inport_)) {
		inport_->set_atom();
	}
	inport_ = NIL_;
	outport_ = NIL_;
	if(is_port(save_inport_)) {
		save_inport_->set_atom();
	}
	save_inport_ = NIL_;
	if(is_port(loadport_)) {
		loadport_->set_atom();
	}
	loadport_ = NIL_;
	store_.set_gc_verbose(false);
	//	This will finalize all cells, freeing ports, frames, strings, etc.
	store_.gc();
	stack_frame::cleanup();
	port::cleanup();
}


void scheme::set_external_data(void *p) {
	ext_data_ = p;
}

//	read from loadport and evaluate expressions
int scheme::load() {
#if CHECK
	assert(is_valid());
#endif
	port* port = load_.reset();
	loadport_ = mk_port(port);
	inport_ = loadport_;
	eval_cycle_top_level();
	loadport_->set_atom();
	if(retcode_ == 0) {
		retcode_ = (nesting_ != 0);
	}
	return retcode_;
}

//	load from a FILE
//	entry point if standalone
int scheme::load(FILE* fin) {
#if CHECK
	assert(is_valid());
#endif
	port* port = load_.reset();
	port->set(port::port_input, fin, true);
	if(fin == stdin) {
		interactive_repl_ = true;
	}
	return load();
}

//	load from a string
//	entry point if standalone
int scheme::load(char* cmd) {
#if CHECK
	assert(is_valid());
#endif
	port* port = load_.reset();
	port->set(port::port_input, cmd, cmd+strlen(cmd), cmd);
	interactive_repl_ = false;
	return load();
}

//	Define the symbol in the given environment.
//	If envir == 0 then define in global environment.
//	If already defined, replace the definition.
void scheme::define(pointer symbol, pointer value, pointer envir) {
#if CHECK
	assert(is_valid());
#endif
	if (envir == 0) envir = global_env_;
	env_.define_binding(symbol, value, envir, true, true);
}

//	Apply the proc "procname" with no args.
//	Entry point if not standalone
int scheme::apply0(const char* procname) {
#if CHECK
	assert(is_valid());
#endif
	code_ = cons(define_symbol(procname), NIL_);
	interactive_repl_ = false;
	eval_cycle_eval();
	return retcode_;
}

//	Apply the given function with the given args.
int scheme::call(pointer func, pointer args) { 
 #if CHECK
	assert(is_valid());
#endif
	code_ = func; 
	args_ = args; 
	interactive_repl_ = false; 
	return eval_cycle_apply(); 
} 

