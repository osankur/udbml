#ifndef UDBM_STUBS_H_
#define UDBM_STUBS_H_

#include <vector>
#include <dbm/dbm.h>
using namespace dbm;
using namespace base;
typedef std::vector<int> carray_t;
#define get_cvector(x) ((carray_t*)Data_custom_val(x))


/* This is the main type for Dbms.
		The encapsulation is needed since
		we free dbm_t* manually, while
		the dbm_wrap_t will be freed by the gc.
*/
typedef struct dbm_wrap_t {
  dbm_t * d;
} dbm_wrap_t;

typedef struct fed_wrap_t {
  fed_t * f;
} fed_wrap_t;

typedef struct fed_it_wrap_t {
	fed_t::iterator * d;
} fed_it_wrap_t;

typedef struct bitvector_wrap_t {
	BitString * b;
} bitvector_wrap_t;


//#define get_dbm_ptr(x) static_cast<dbm_t*>(Data_custom_val(x))
#define get_dbm_ptr(x) (((dbm_wrap_t*)Data_custom_val(x))->d)

// forward declaration
namespace dbm {
    class pdbm_t;
}

// smart inclusion test for unpriced zones
bool
dbm_closure_leq(const raw_t * const dr1, const raw_t * const dr2, cindex_t dim,
                const std::vector<int> &lbounds, const std::vector<int> &ubounds);

// smart inclusion test for priced zones
bool
pdbm_square_inclusion_exp(const dbm::pdbm_t &, const dbm::pdbm_t &, const std::vector<int> &);

#endif  // UDBM_STUBS_H_
