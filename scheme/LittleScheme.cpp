
#include "scheme.h"
#include "LittleScheme.h"

scheme* scheme_factory::create(int first_cell_seg, int max_cell_seg) {
	return new scheme(first_cell_seg, max_cell_seg);
}

LittleScheme::LittleScheme(scheme_factory& sf, int first_cell_seg, int max_cell_seg) {
	sc_ = sf.create(first_cell_seg, max_cell_seg);
}

LittleScheme::~LittleScheme() {
	delete sc_;
}

bool LittleScheme::is_valid() { return sc_ != 0 && sc_->is_valid(); }

void LittleScheme::set_input_port(FILE *fin) { sc_->set_input_port(fin); }
void LittleScheme::set_input_port(char *start,	///<	The input string.
	char *past_the_end				///<	Points one past the string end.
	) { sc_->set_input_port(start, past_the_end);}
void LittleScheme::set_output_port(FILE *fin) { sc_->set_output_port(fin); }	
void LittleScheme::set_output_port(char *start,	///<	The output string
	char *past_the_end				///<	Points one past the string end.
	) { sc_->set_output_port(start, past_the_end); }
int LittleScheme::load(FILE *fin) { return sc_->load(fin); }
int LittleScheme::load(char *cmd) { return sc_->load(cmd); }
int LittleScheme::apply0(const char *procname) { return sc_->apply0(procname); }
int LittleScheme::call(pointer func, pointer args) { return sc_->call(func, args); }
void LittleScheme::set_external_data(void* p) { sc_->set_external_data(p); }
void LittleScheme::define(pointer symbol, pointer value, pointer env) {
	sc_->define(symbol, value, env);
}

pointer LittleScheme::cons(pointer a, pointer b, bool immutable) {
	return sc_->cons(a, b, immutable);
}
pointer LittleScheme::immutable_cons(pointer a, pointer b) { return sc_->cons(a, b, true); }

pointer LittleScheme::mk_integer(long num) { return sc_->mk_integer(num); }
pointer LittleScheme::mk_real(double num) { return sc_->mk_real(num); }
pointer LittleScheme::mk_atom(const char *str, bool must_be_number, int radix) { 
	return sc_->mk_atom(str, must_be_number, radix); 
}
pointer LittleScheme::gensym() { return sc_->gensym(); }
pointer LittleScheme::mk_string(const char* str) { return sc_->mk_string(str); }
pointer LittleScheme::mk_string(const char* str, int len) { return sc_->mk_string(str, len); }
pointer LittleScheme::mk_string(int len, char fill) { return sc_->mk_string(len, fill); }
pointer LittleScheme::mk_character(int c) { return sc_->mk_character(c); }
pointer LittleScheme::mk_vector(int len) { return sc_->mk_vector(len); }
pointer LittleScheme::mk_foreign_func(foreign_func f) { return sc_->mk_foreign_func(f); }
pointer LittleScheme::mk_environment(pointer new_frame, pointer old_env) {
	return sc_->mk_environment(new_frame, old_env); 
}
pointer LittleScheme::define_symbol(const char* name) { return sc_->define_symbol(name); }

pointer LittleScheme::NIL() { return scheme::NIL_; }
pointer LittleScheme::T() { return scheme::T_; }
pointer LittleScheme::F() { return scheme::F_; }
bool LittleScheme::is_true(pointer p) { return scheme::is_true(p); }
bool LittleScheme::is_false(pointer p) { return scheme::is_false(p); }

pointer LittleScheme::global_env() const { return sc_->global_env(); }
int LittleScheme::retcode() const { return sc_->retcode(); }

pointer LittleScheme::reverse(pointer a) { return sc_->reverse(a); }
pointer LittleScheme::append(pointer a, pointer b) { return sc_->append(a, b); }

pointer LittleScheme::assign_syntax(const char* name) { return sc_->assign_syntax(name); }

