#ifndef SCHEME_MEMBLK
#define SCHEME_MEMBLK

#include "scheme.h"


class cell_memblock : public cell {
public:
	enum {
		T_MEMBLOCK =	    'b',	///<	stores a block
	};

    bool is_memblock() const  { return type_ == T_MEMBLOCK; }
	int* memvalue() const { 
#if CHECK
		assert(is_memblock());
#endif
		return object.memblock.mvalue_; 
	}

	int memlen() const { 
#if CHECK
		assert(is_memblock());
#endif
		return object.memblock.mlen_; 
	}
	void set_memblock(int* mvalue, int mlen) { 
		type_ = T_MEMBLOCK;
		flag_ = F_ATOM | F_MARKFUN | F_FINALFUN;
		object.memblock.mvalue_ = mvalue;	
		object.memblock.mlen_ = mlen;
	}
};

typedef cell_memblock* mem_pointer;

#define TST_BLOCK   15

class scheme_memblock : public scheme {
public:
	 scheme_memblock(int first_cell_seg, int max_cell_seg);
	 void op_mkblock(scheme& sc);
	 void op_blocklen(scheme& sc);
	 void op_blockref(scheme& sc);
	 void op_blockset(scheme& sc);
	 void op_blockp(scheme& sc);
private:
	pointer mk_memblock(int len, int fill);
};



#endif
