/*
 * GetSet.h
 *
 *  Created on: Jun 20, 2019
 *      Author: cotrik
 */

#ifndef LIBCOTRIK_SRC_GETSET_H_
#define LIBCOTRIK_SRC_GETSET_H_

#define GETSET(Type, MemberName) \
    Type Get##MemberName() const { \
        return MemberName; \
    }; \
    void Set##MemberName(Type value) { \
        MemberName = value; \
    }

#define GETSETR(Type, MemberName) \
    const Type &Get##MemberName() const { \
        return MemberName; \
    }; \
    void Set##MemberName(const Type &value) { \
        MemberName = value; \
    }

#define GET(Type, MemberName) \
    Type Get##MemberName() const { \
        return MemberName; \
    }

#define GETR(Type, MemberName) \
    const Type &Get##MemberName() const { \
        return MemberName; \
    }

#define GETRNC(Type, MemberName) \
    Type &Get##MemberName() { \
        return MemberName; \
    }

#define SET(Type, MemberName) \
    void Set##MemberName(const Type &value) { \
        MemberName = value; \
    }

#endif /* LIBCOTRIK_SRC_GETSET_H_ */
