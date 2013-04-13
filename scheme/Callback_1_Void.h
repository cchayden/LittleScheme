#ifndef CALLBACK_1_VOID_H
#define CALLBACK_1_VOID_H

//////////////////////////////////////////////////////
//                                                  //
// NOTICE:  THIS IS A GENERATED FILE -- DO NOT EDIT //
//                                                  //
//////////////////////////////////////////////////////


#include "CallbackFlags.h"


/**
 * @class Callback_1_Void<Arg1>
 *
 * @brief Base class for Callbacks that take 1 argument(s) and return void.
 */

template<class Arg1>
class Callback_1_Void {
public:
	virtual void	operator () (Arg1 a1) = 0;
	void	destroy() { delete this; }
	// Dtor should be protected to prevent folks from making one on
	// the stack, but then we need a way to delete memory.
	virtual ~Callback_1_Void() { }
};

template<class Function, class Arg1>
class Callback_1_Void_Function : public Callback_1_Void<Arg1> {
public:
	Callback_1_Void_Function(Function fn):
		function_(fn)
	{
	}

	virtual void operator () (Arg1 a1) {
		(*function_)(a1);
	}
private:
	Function function_;
};

template<class Function, class Arg1>
Callback_1_Void_Function<Function, Arg1>*
newCallback(const Callback_1_Void<Arg1>*, Function fn) {
	return new Callback_1_Void_Function<Function, Arg1>(fn);
}

template<class Object, class Member, class Arg1>
class Callback_1_Void_Member : public Callback_1_Void<Arg1> {
public:
	Callback_1_Void_Member(Object ot, Member mr):
		object_(ot),
		member_(mr)
	{
	}

	virtual void operator () (Arg1 a1) {
		(object_->*member_)(a1);
	}
private:
	Object object_;
	Member member_;
};

template<class Object, class Member, class Arg1>
Callback_1_Void_Member<Object, Member, Arg1>*
newCallback(const Callback_1_Void<Arg1>*, Object ot, Member mr) {
	return new Callback_1_Void_Member<Object, Member, Arg1>(ot, mr);
}

template<class Object, class Member, class Arg1>
class Callback_1_Void_Tmp_Member : public Callback_1_Void<Arg1> {
public:
	Callback_1_Void_Tmp_Member(Object ot, Member mr):
		object_(ot),
		member_(mr)
	{
	}
	~Callback_1_Void_Tmp_Member() {
		delete object_;
		object_ = 0;
	}

	virtual void operator () (Arg1 a1) {
		(object_->*member_)(a1);
	}
private:
	Object object_;
	Member member_;
};

template<class Object, class Member, class Arg1>
Callback_1_Void_Tmp_Member<Object, Member, Arg1>*
newCallback(const Callback_1_Void<Arg1>*, CallbackFlags::Value, Object ot, Member mr) {
	return new Callback_1_Void_Tmp_Member<Object, Member, Arg1>(ot, mr);
}

#endif /* CALLBACK_1_VOID_H */
