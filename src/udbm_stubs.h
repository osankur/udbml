#ifndef UDBM_STUBS_H_
#define UDBM_STUBS_H_

#include <vector>

typedef std::vector<int> carray_t;
#define get_cvector(x) ((carray_t*)Data_custom_val(x))

#define get_dbm_ptr(x) static_cast<dbm_t*>(Data_custom_val(x))

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
