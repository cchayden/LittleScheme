
#include "glob.h"

#include "scheme_memblk.h"

bool is_memblock(pointer p) { return ((mem_pointer)p)->is_memblock(); }
int* memvalue(pointer p) { return ((mem_pointer)p)->memvalue(); }
int memlen(pointer p) { return ((mem_pointer)p)->memlen(); }

//	test table
static scheme::test_entry static_tests[] = {
  scheme::test_entry(TST_BLOCK, is_memblock, "block"),
  scheme::test_entry(-1, 0, 0)
};

static void mark_memblock(storage& st, pointer p) {
	int* mvalue = ((mem_pointer)p)->memvalue();
	int mlen = ((mem_pointer)p)->memlen();
	for (int i = 0; i < mlen; i++) {
		st.mark((pointer)mvalue[i]);
	}
}

//	mark table
static storage::mark_entry static_markers[] = {
	storage::mark_entry(cell_memblock::T_MEMBLOCK, mark_memblock),
	storage::mark_entry(0, 0),
};

static void finalize_memblock(storage& st, pointer p) { 
	delete [] ((mem_pointer)p)->memvalue();
}

//	finalize table
static storage::final_entry static_finalizers[] = {
	storage::final_entry(cell_memblock::T_MEMBLOCK, finalize_memblock),
	storage::final_entry(0, 0),
};

//	op code table
static op_code_ent<scheme_memblock> static_entry[] = {
	op_code_ent<scheme_memblock>(-1, op_code_entry::proc, scheme_memblock::op_mkblock, 
	   "make-block", 1, 2, TST_NATURAL, TST_INTEGER),
	op_code_ent<scheme_memblock>(-1, op_code_entry::proc, scheme_memblock::op_blocklen, 
	   "block-len", 1, 1, TST_BLOCK),
	op_code_ent<scheme_memblock>(-1, op_code_entry::proc, scheme_memblock::op_blockref, 
	   "block-ref", 2, 2, TST_BLOCK, TST_NATURAL),
	op_code_ent<scheme_memblock>(-1, op_code_entry::proc, scheme_memblock::op_blockset, 
	   "block-set!", 3, 3, TST_BLOCK, TST_NATURAL, TST_INTEGER),
	op_code_ent<scheme_memblock>(-1, op_code_entry::proc, scheme_memblock::op_blockp, 
	   "block?", 1, 1, TST_ANY),
	op_code_ent<scheme_memblock>(Op::ILLEGAL, op_code_entry::prim, 0, "")
};

static const char* print_memblock(scheme&, pointer, bool, int, char*) { 
	return "#<MEMBLK>"; 
}

//	print table
static scheme::print_entry static_print[] = {
	scheme::print_entry(cell_memblock::T_MEMBLOCK, print_memblock),
	scheme::print_entry(0, 0),
};

//	constructor
scheme_memblock::scheme_memblock(int first_cell_seg, 
			 int max_cell_seg): 
	scheme(first_cell_seg, max_cell_seg) {

	bool ok = extend_tests(static_tests);

	int n = op_table_.extend(static_entry, this);

	store_.extend_marker(static_markers);
	store_.extend_finalizer(static_finalizers);

	extend_print(static_print);
}

//	make memblock
pointer scheme_memblock::mk_memblock(int len, int fill) {
	mem_pointer x = (mem_pointer)store_.find_cell();
	int* p = new int[len];
	for (int i = 0; i < len; i++) p[i] = fill;
	x->set_memblock(p, len);
	return x;
}

//
//	op codes
//
void scheme_memblock::op_mkblock(scheme& sc) {
	int len = ivalue(arg0());	// already verified nonneg integer
	int fill = 0;
	if (arg_tail() != NIL_) {
		fill = ivalue(arg1());
	}
	s_return(mk_memblock(len, fill));
}

void scheme_memblock::op_blocklen(scheme& sc) {
	s_return(mk_integer(memlen(arg0())));
}

void scheme_memblock::op_blockref(scheme& sc) {
	int* mem = memvalue(arg0());  // guaranteed it is a block
	int index = ivalue(arg1());    // guaranteed it is nonneg integer    
	if(index >= memlen(arg0())) {
		s_error1("block-ref: out of bounds:", arg1());
	}
	s_return(mk_integer(mem[index]));
}

void scheme_memblock::op_blockset(scheme& sc) {
	if(is_immutable(arg0())) {
		s_error1("block-set!: unable to alter immutable memory block:", arg0());
	}
	int* mem = memvalue(arg0());	// guaranteed to be block
	int index=ivalue(arg1());    // guaranteed to be nonneg
	if(index >= memlen(arg0())) {
		s_error1("block-set: out of bounds:", arg1());
	}
	mem[index] = ivalue(arg2());   // guaranteed to be integer
	s_return(VOID_);
}

void scheme_memblock::op_blockp(scheme& sc) {
	s_retbool(is_memblock(arg0()));
}

