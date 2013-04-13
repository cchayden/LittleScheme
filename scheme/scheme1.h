#ifndef _SCHEME1_H
#define _SCHEME1_H

#include "scheme.h"

class scheme1 : public scheme {
public:
	scheme1(int first_cell_seg, 
				 int max_cell_seg);

	//	testing
	void test1(scheme& sc);
	void test2a(scheme& sc);
	void test2b(scheme& sc);
	void test3(scheme& sc);
};



//	testing
void scheme1::test1(scheme& sc) {
	s_return(mk_string("this is test1"));
}

void scheme1::test2a(scheme& sc) {
	s_return(mk_string("this is test2a"));
}

void scheme1::test2b(scheme& sc) {
	s_return(mk_string("this is test2b"));
}

void scheme1::test3(scheme& sc) {
	s_return(mk_string("this is test3"));
}

static op_code_ent<scheme1> static_entry[] = {
	op_code_ent<scheme1>(0, op_code_entry::proc, scheme1::test3, "test3"),
	op_code_ent<scheme1>(Op::ILLEGAL, op_code_entry::prim, 0, "")
};

scheme1::scheme1(int first_cell_seg, 
				 int max_cell_seg): 
	scheme(first_cell_seg, max_cell_seg) {

	// testing

		//	Create op_code_entry table, and extend with it.
	op_code_entry entry1[2];
	entry1[0] = op_code_entry(0, op_code_entry::proc, 
		newCallback((op_code_entry::op_callback)0, this, &scheme1::test1), 
		"test1");
	entry1[1] = op_code_entry(Op::ILLEGAL, op_code_entry::prim, 0, "illegal");
	op_table_.extend(entry1, this);

	//	Create op_code_ent<scheme1> and extend with it
	op_code_ent<scheme1> entry2[3];
	entry2[0] = op_code_ent<scheme1>(-1, op_code_entry::proc, scheme1::test2a, "test2a");
	entry2[1] = op_code_ent<scheme1>(0, op_code_entry::proc, scheme1::test2b, "test2b");
	entry2[2] = op_code_ent<scheme1>(Op::ILLEGAL, op_code_entry::prim, 0, "illegal");
	op_table_.extend(entry2, this);

	//	use static table
	op_table_.extend(static_entry, this);
}

#endif
