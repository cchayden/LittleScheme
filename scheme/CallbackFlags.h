#ifndef CALLBACK_FLAGS_H
#define CALLBACK_FLAGS_H

/**
 * @class CallbackFlags
 *
 * @brief This class exists to give a scope to the Flags enumerations used by
 * the derived callback classes to pass flags, such as the
 * DELETE_OBJECT_WHEN_DONE flag.
 */

class CallbackFlags {
public:
	enum Value {
		DELETE_OBJECT_WHEN_DONE
	};
};

#endif /* CALLBACK_FLAGS_H */
