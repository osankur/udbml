#ifndef UDBM_STUBS_H_
#define UDBM_STUBS_H_

#include <vector>

typedef std::vector<int> carray_t;
#define get_cvector(x) ((carray_t*)Data_custom_val(x))

#define get_dbm_ptr(x) static_cast<dbm_t*>(Data_custom_val(x))

#endif  // UDBM_STUBS_H_
