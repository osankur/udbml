extern "C" {
#include <stdio.h>
#include <stdlib.h>
#include <caml/mlvalues.h>
#include <caml/memory.h>
#include <caml/alloc.h>
#include <caml/callback.h>
#include <caml/custom.h>
}
#include <dbm/fed.h>
#include <dbm/pfed.h>
#include <base/bitstring.h>

#include <unordered_map>

#include "udbm_stubs.h"

#define get_pdbm_ptr(x) (static_cast<pdbm_wrap_t*>(Data_custom_val(x)))
#define get_pfed_ptr(x) (static_cast<pfed_wrap_t*>(Data_custom_val(x)))
#define get_pfed_it_ptr(x) (static_cast<pfed_it_wrap_t*>(Data_custom_val(x)))

using namespace dbm;

typedef pdbm_t pdbm_wrap_t;
typedef pfed_t pfed_wrap_t;
typedef pfed_t::iterator pfed_it_wrap_t;

extern "C" void finalize_pdbm(value v){
    get_pdbm_ptr(v)->~pdbm_t();
}

extern "C" void finalize_pfed(value v){
    get_pfed_ptr(v)->~pfed_t();
}

extern "C" void finalize_pfed_it(value v){
    typedef pfed_t::iterator pfed_it;
    get_pfed_it_ptr(v)->~pfed_it();
}

extern "C" int compare_pdbm(value v1, value v2){
    const pdbm_t * i1 = get_pdbm_ptr(v1);
    const pdbm_t * i2 = get_pdbm_ptr(v2);
    relation_t rel = pdbm_relation(*i1, *i2, i1->getDimension());
    if (rel & base_EQUAL) return 0;
    if (i1 < i2) return -1;
    else if (i1 > i2) return 1;
    return 0;
}

extern "C" int compare_pfed(value v1, value v2){
    const pfed_t * f1 = get_pfed_ptr(v1);
    const pfed_t * f2 = get_pfed_ptr(v2);
    relation_t rel;
    if (f1->size() == 1 && f2->size() == 1)
        rel = f1->exactRelation(*f2);
    else
        rel = f1->relation(*f2);
    if (rel & base_EQUAL) return 0;
    if (f1 < f2) return -1;
    else if (f1 > f2) return 1;
    return 0;
}

extern "C" int compare_pfed_it(value v1, value v2){
    const pfed_it_wrap_t * i1 = get_pfed_it_ptr(v1);
    const pfed_it_wrap_t * i2 = get_pfed_it_ptr(v2);
    if (v1 < v2) return -1;
    if (v1 > v2) return 1;
    return 0;
}

extern "C" long hash_pdbm(value v){
    return (long)get_pdbm_ptr(v)->hash();
}

extern "C" long hash_pfed(value v){
    return (long)get_pfed_ptr(v)->hash(0);
}

extern "C" long hash_pfed_it(value v){
    return (long)get_pfed_it_ptr(v);
}

static struct custom_operations custom_ops_pdbm = {
    .identifier         = (char*)"pdbm_wrap_t handling",
    .finalize           = finalize_pdbm,
    .compare            = compare_pdbm,
    .hash               = hash_pdbm,
    .serialize          = custom_serialize_default,
    .deserialize        = custom_deserialize_default
};

static struct custom_operations custom_ops_pfed = {
    .identifier         = (char*)"pfed_wrap_t handling",
    .finalize           = finalize_pfed,
    .compare            = compare_pfed,
    .hash               = hash_pfed,
    .serialize          = custom_serialize_default,
    .deserialize        = custom_deserialize_default
};

static struct custom_operations custom_ops_pfed_it = {
    .identifier         = (char*)"pfed_it_wrap_t handling",
    .finalize           = finalize_pfed_it,
    .compare            = compare_pfed_it,
    .hash               = hash_pfed_it,
    .serialize          = custom_serialize_default,
    .deserialize        = custom_deserialize_default
};

namespace std
{
    template<>
    struct hash<pdbm_t>
    {
        size_t operator()(const pdbm_t &p) const
        {
            return p.hash();
        }
    };
}

// The pdbm interface
extern "C" CAMLprim value
stub_pdbm_create(value size)
{
    CAMLparam1(size);
    CAMLlocal1(res);
    cindex_t dim = Int_val(size);
    res = caml_alloc_custom(&custom_ops_pdbm, sizeof(pdbm_wrap_t), 0, 1);
    new (Data_custom_val(res)) pdbm_wrap_t(dim);
    CAMLreturn(res);
}

extern "C" CAMLprim value
stub_pdbm_dimension(value pdbm)
{
    return Val_int(get_pdbm_ptr(pdbm)->getDimension());
}

extern "C" CAMLprim value
stub_pdbm_copy(value v)
{
    CAMLparam1(v);
    CAMLlocal1(res);
    res = caml_alloc_custom(&custom_ops_pdbm, sizeof(pdbm_wrap_t), 0, 1);
    new (Data_custom_val(res)) pdbm_wrap_t(*get_pdbm_ptr(v));
    CAMLreturn(res);
}

extern "C" CAMLprim value
stub_pdbm_set_init(value v)
{
    pdbm_t * d = get_pdbm_ptr(v);
    pdbm_init(*d, d->getDimension());
    return Val_unit;
}

extern "C" CAMLprim value
stub_pdbm_set_zero(value v)
{
    pdbm_t * d = get_pdbm_ptr(v);
    pdbm_zero(*d, d->getDimension());
    return Val_unit;
}

extern "C" CAMLprim value
stub_pdbm_is_empty(value v)
{
    pdbm_wrap_t * pdbm = get_pdbm_ptr(v);
    return Val_bool(pdbm_isEmpty(*pdbm, pdbm->getDimension()));
}

extern "C" CAMLprim value
stub_pdbm_is_unbounded(value v)
{
    pdbm_wrap_t * pdbm = get_pdbm_ptr(v);
    return Val_bool(pdbm_isUnbounded(*pdbm, pdbm->getDimension()));
}

extern "C" CAMLprim value
stub_pdbm_hash(value v)
{
    pdbm_wrap_t * pdbm = get_pdbm_ptr(v);
    return Val_int(pdbm->hash());
}

extern "C" CAMLprim value
stub_pdbm_at(value v, value i, value j)
{
    CAMLparam3(v, i, j);
    CAMLlocal1(res);
    const pdbm_wrap_t & pdbm = *get_pdbm_ptr(v);
    raw_t r = pdbm(Int_val(i), Int_val(j));
    res = caml_alloc(2, 0);
    Store_field(res, 0, Val_int(dbm_raw2bound(r)));
    Store_field(res, 1, Val_int(dbm_raw2strict(r)));
    CAMLreturn(res);
}

extern "C" CAMLprim value
stub_pdbm_at_bound(value v, value i, value j)
{
    const pdbm_wrap_t & pdbm = *get_pdbm_ptr(v);
    raw_t r = pdbm(Int_val(i), Int_val(j));
    return Val_int(dbm_raw2bound(r));
}

extern "C" CAMLprim value
stub_pdbm_equal(value v1, value v2)
{
    const pdbm_wrap_t &p1 = *get_pdbm_ptr(v1);
    const pdbm_wrap_t &p2 = *get_pdbm_ptr(v2);
    relation_t rel = pdbm_relation(p1, p2, p1.getDimension());
    if (rel & base_EQUAL) return Val_bool(true);
    return Val_bool(false);
}

extern "C" CAMLprim value
stub_pdbm_not_equal(value v1, value v2)
{
    const pdbm_wrap_t &p1 = *get_pdbm_ptr(v1);
    const pdbm_wrap_t &p2 = *get_pdbm_ptr(v2);
    relation_t rel = pdbm_relation(p1, p2, p1.getDimension());
    if (rel & base_EQUAL) return Val_bool(false);
    return Val_bool(true);
}

extern "C" CAMLprim value
stub_pdbm_lt(value v1, value v2)
{
    const pdbm_wrap_t &p1 = *get_pdbm_ptr(v1);
    const pdbm_wrap_t &p2 = *get_pdbm_ptr(v2);
    relation_t rel = pdbm_relation(p1, p2, p1.getDimension());
    if (rel == base_LESS) return Val_bool(true);
    return Val_bool(false);
}

extern "C" CAMLprim value
stub_pdbm_gt(value v1, value v2)
{
    const pdbm_wrap_t &p1 = *get_pdbm_ptr(v1);
    const pdbm_wrap_t &p2 = *get_pdbm_ptr(v2);
    relation_t rel = pdbm_relation(p1, p2, p1.getDimension());
    if (rel == base_GREATER) return Val_bool(true);
    return Val_bool(false);
}

extern "C" CAMLprim value
stub_pdbm_leq(value v1, value v2)
{
    const pdbm_wrap_t &p1 = *get_pdbm_ptr(v1);
    const pdbm_wrap_t &p2 = *get_pdbm_ptr(v2);
    relation_t rel = pdbm_relation(p1, p2, p1.getDimension());
    if (rel & base_LESS) return Val_bool(true);
    return Val_bool(false);
}

extern "C" CAMLprim value
stub_pdbm_geq(value v1, value v2)
{
    const pdbm_wrap_t &p1 = *get_pdbm_ptr(v1);
    const pdbm_wrap_t &p2 = *get_pdbm_ptr(v2);
    relation_t rel = pdbm_relation(p1, p2, p1.getDimension());
    if (rel & base_GREATER) return Val_bool(true);
    return Val_bool(false);
}

// a class that stores the preorder on clocks for a given zone
class clock_po_t
{
    friend class clock_po_iterator;
public:
    explicit clock_po_t() {}

    size_t
    size() const
    {
        return _order.size();
    }

    const std::vector<std::vector<int> > &
    operator[](int i) const
    {
        return _order[i];
    }

    void
    push_back(const std::vector<std::vector<int>> &v)
    {
        _order.push_back(v);
    }

    void
    push_back(std::vector<std::vector<int>> &&v)
    {
        _order.push_back(v);
    }

    void
    add_below(int i)
    {
        _below.push_back(i);
    }

    void
    add_above(int i)
    {}

    std::vector<std::vector<std::vector<int> > >::iterator
    begin()
    {
        return _order.begin();
    }

    std::vector<std::vector<std::vector<int> > >::iterator
    end()
    {
        return _order.begin();
    }

    std::vector<std::vector<std::vector<int> > >::const_iterator
    begin() const
    {
        return _order.begin();
    }

    std::vector<std::vector<std::vector<int> > >::const_iterator
    end() const
    {
        return _order.begin();
    }

private:
    std::vector<int> _below;
    // there is fact no need to store clocks that are above their bound
    std::vector<std::vector<std::vector<int> > > _order;
};

class clock_po_iterator
{
public:
    explicit clock_po_iterator(const clock_po_t &o, int dim)
    : _order(o)
    , _state(o.size())
    , _is_in(dim)
    , _size(0)
    {
        for (auto &i : _order._below)
        {
            _is_in[i] = true;
            ++_size;
        }
    }

    bool
    is_in(int i) const
    {
        return _is_in[i];
    }

    cindex_t
    count() const
    {
        return _size;
    }

    bool done() const { return _size == _order._below.size(); }

    void
    next()
    {
        for (size_t i = 0; i < _order.size(); ++i)
        {
            if (_state[i] == _order[i].size())
            {
                _state[i] = 0;
                // remove indices from _is_in
                for (size_t j = 0; j < _order[i].size(); ++j)
                {
                    for (auto &v : _order[i][j])
                    {
                        _is_in[v] = false;
                        --_size;
                    }
                }
            }
            else
            {
                // add indices to is_in
                for (auto &v : _order[i][_state[i]])
                {
                    assert(v < _is_in.size());
                    _is_in[v] = true;
                    ++_size;
                }
                ++(_state[i]);
                return;
            }
        }
    }

    void
    print_preorder() const
    {
        printf("preorder is\n");
        for (size_t i = 0; i < _order.size(); ++i)
        {
            printf("chain:\n");
            for (size_t j = 0; j < _order[i].size(); ++j)
            {
                printf("cycle:");
                for (auto &v : _order[i][j])
                {
                    printf("%d,", v);
                }
                printf("\n");
            }
        }
        printf("\n\n");
    }


private:
    const clock_po_t & _order;
    std::vector<int> _state; // same size as _order
    std::vector<bool> _is_in; // same size as dim
    cindex_t _size; // the size of the current subset
};

clock_po_t
build_clock_preorder(const pdbm_t &z, const std::vector<int> &mbounds)
{
    int dim = z.getDimension();
    // build the preorder on clocks wrt z
    clock_po_t result;
    for (int i = 1; i < dim; ++i)
    {
        if (z(i,0) <= dbm_boundbool2raw(mbounds[i], false))
        {
            result.add_below(i);
            continue;
        }
        if (z(0,i) <= dbm_boundbool2raw(-mbounds[i], true))
        {
            result.add_above(i);
            continue;
        }

        bool to_insert = true;
        for (auto it = result.begin(); to_insert && it != result.end(); ++it)
        {
            int k = 0;
            while (to_insert)
            {
                if (k == it->size())
                {
                    it->insert(it->begin()+k, std::vector<int>(1,i));
                    to_insert = false;
                }
                else
                {
                    int j = (*it)[k][0];
                    bool i_preceq_j;
                    bool j_preceq_i;
                    i_preceq_j = (z(i,j) <= dbm_boundbool2raw(mbounds[i] - mbounds[j], false));
                    j_preceq_i = (z(j,i) <= dbm_boundbool2raw(mbounds[j] - mbounds[i], false));

                    if (i_preceq_j && j_preceq_i)
                    {
                        (*it)[k].push_back(i);
                        to_insert = false;
                    }
                    else if (i_preceq_j)
                    {
                        it->insert(it->begin()+k, std::vector<int>(1,i));
                        to_insert = false;
                    }
                    else if (j_preceq_i)
                    {
                        ++k;
                    }
                    else
                    {
                        break;
                    }
                }
            }

        }

        if (to_insert)
        {
            result.push_back(std::vector<std::vector<int> >(1, std::vector<int>(1,i)));
        }
    }
    return result;
}

// a function that builds the product priced zone for a given subset of clocks Y
// the corresponding function has rates
//   - r(x)     if x is in Y
//   r'(x)      if x' is in Y'
//   - r(x)     if x is in X-Y
//   - r'(x)    if x' is in X'-Y'
pdbm_t
pdbm_build_product(const pdbm_t &z1,
                   const pdbm_t &z2,
                   const std::vector<int> &mbounds,
                   const clock_po_iterator &y)
{
    assert(z1.getDimension() == z2.getDimension());
    int dim = z1.getDimension();
    int ndim = 2*dim-1-y.count();
    // create a new dbm of size ndim (the reference clock and those in Y need not be duplicated)
    pdbm_t result = pdbm_t(ndim);
    // initialize
    pdbm_init(result, ndim);
    // constrain it
    // TODO     consider using pdbm_constrainN instead
    //          or build the dbm as a matrix of raw_t and canonize once at the end
    // WARNING  when building the corresponding cost function
    //          negate the function, so that its infimum corresponds to the supremum we want
    //          (the API computes the infimum rather than the supremum)
    for (int i = 0, ki = dim; i < dim; ++i)
    {
        int ii = (i && !y.is_in(i)) ? ki++ : i;
        if (i > 0)
        {
            // constrain it with x <= M for x in Y, x > M otherwise
            if (y.is_in(i))
            {
                pdbm_constrain1(result, ndim, i, 0, dbm_boundbool2raw(mbounds[i], false));

                // set rates
                pdbm_setRate(result, ndim, i, pdbm_getRate(z1, dim, i) - pdbm_getRate(z2, dim, i));
            }
            else
            {
                pdbm_constrain1(result, ndim, 0, i, dbm_boundbool2raw(-mbounds[i], true));
                pdbm_constrain1(result, ndim, 0, ii, dbm_boundbool2raw(-mbounds[i], true));

                // set rates
                pdbm_setRate(result, ndim, i, pdbm_getRate(z1, dim, i));
                pdbm_setRate(result, ndim, ii, pdbm_getRate(z2, dim, i));
            }
        }

        for (int j = 0, kj = dim; j < dim; ++j)
        {
            int jj = (j && !y.is_in(j)) ? kj++ : j;
            // constrain from z1
            pdbm_constrain1(result, ndim, i, j, z1(i,j));
            // constrain from z2
            pdbm_constrain1(result, ndim, ii, jj, z2(i,j));
        }
    }

    // return the result
    return result;
}

// check whether a priced zone dominates another
// this is based on the exploration of all subsets of clocks
bool
pdbm_square_inclusion_exp(const pdbm_t &z1, const pdbm_t &z2, const std::vector<int> &mbounds)
{
    // check whether all valuations in Z have an equivalent in Z'
    // it is the same as testing the inclusion of Z in the closure of Z'
    if (! dbm_closure_leq(z1.const_dbm(), z2.const_dbm(), z1.getDimension(), mbounds, mbounds))
        return false;

    // a cache for preorders on clocks
    // TODO use dbm_t (and not pdbm_t) as key for the map
    static std::unordered_map<pdbm_t, clock_po_t> order_cache;
    clock_po_t preorder;
    {
        auto it = order_cache.find(z1);
        if (it == order_cache.end())
        {
            // compute the preorder on clocks
            preorder = build_clock_preorder(z1, mbounds);
            // and store it
            order_cache[z1] = preorder;
        }
        else
        {
            preorder = it->second;
        }
    }

    // get the dimension
    int dim = z1.getDimension();
    // initialize Y to emptyset
    clock_po_iterator currentY(preorder, dim);

    // the main loop
    do
    {
        // first check whether Z_Y is empty
        pdbm_t zy1 = z1;
        bool zy1_not_empty = true;
        for (int i = 1; zy1_not_empty && i < dim; ++i)
        {
            if (currentY.is_in(i))
            {
                zy1_not_empty = pdbm_constrain1(zy1, dim, i, 0, dbm_boundbool2raw(mbounds[i], false));
            } else {
                zy1_not_empty = pdbm_constrain1(zy1, dim, 0, i, dbm_boundbool2raw(-mbounds[i], true));
            }
        }

        // TODO do not build zy1, build the product directly
        //      zy1 empty iff the product is empty

        if (zy1_not_empty)
        {
            // build the product for the current Y
            pdbm_t prody = pdbm_build_product(z1, z2, mbounds, currentY);

            // get the valuation where the infimum of the product zone is reached
            int32_t * inf_val = new int32_t[prody.getDimension()];
            // TODO do not build this array every time, make it static somehow
            bool * free_clocks = new bool[prody.getDimension()];
            for (int i = 0; i < prody.getDimension(); ++i)
            {
                free_clocks[i] = true;
            }
            pdbm_getInfimumValuation(prody, prody.getDimension(), inf_val, free_clocks);

            // use this valuation to evaluate the searched sup for the current Y
            // inf_val = (v0,v0') and local_sup = z2(v0') - z1(v0)
            int32_t * v0p = new int32_t[dim];
            v0p[0] = inf_val[0];
            for (int i = 1, ki = dim; i < dim; ++i)
            {
                int ii = currentY.is_in(i) ? i : ki++;
                v0p[i] = inf_val[ii];
            }
            int32_t local_sup =   pdbm_getCostOfValuation(z2, dim, v0p)
                                - pdbm_getCostOfValuation(z1, dim, inf_val);

            // free
            delete [] v0p;
            delete [] free_clocks;
            delete [] inf_val;

            // if positive, return early
            if (local_sup > 0)
                return false;
        }

        // go to next Y
        currentY.next();

    } while (! currentY.done());

    return true;
}

extern "C" CAMLprim value
stub_pdbm_square_inclusion_exp(value t1, value t2, value mvec)
{
    bool result = pdbm_square_inclusion_exp(*get_pdbm_ptr(t1), *get_pdbm_ptr(t2), *get_cvector(mvec));
    return Val_bool(result);
}

extern "C" CAMLprim value
stub_pdbm_constrain(value t, value ct)
{
    // a constraint is (i,j,(b,ineq))
    int i = Int_val(Field(ct,0));
    int j = Int_val(Field(ct,1));
    int b = Int_val(Field(Field(ct,2),0));
    int ineq = Int_val(Field(Field(ct,2),1));
    pdbm_t * d = get_pdbm_ptr(t);
    raw_t r = dbm_boundbool2raw(b, ineq == 0);
    pdbm_constrain1(*d, d->getDimension(), i, j, r);
    return Val_unit;
}

extern "C" CAMLprim value
stub_pdbm_intersect(value v1, value v2)
{
    pdbm_t * d1 = get_pdbm_ptr(v1);
    const dbm_t & d2 = *get_dbm_ptr(v2);
    raw_t * rdbm = pdbm_getMutableMatrix(*d1, d1->getDimension());
    // TODO d1 and d2 should be non-empty
    // dbm_intersection returns true iff d1 is non-empty afterwards
    dbm_intersection(rdbm, d2.const_dbm(), d1->getDimension());
    pdbm_close(*d1, d1->getDimension());
    return Val_unit;
}

extern "C" CAMLprim value
stub_pdbm_infimum(value v)
{
    pdbm_t * d = get_pdbm_ptr(v);
    return Val_int(pdbm_getInfimum(*d, d->getDimension()));
}

extern "C" CAMLprim value
stub_pdbm_to_string(value v)
{
    CAMLparam1(v);
    pdbm_wrap_t * d = get_pdbm_ptr(v);
    std::stringstream os;
    pdbm_print(os, *d, d->getDimension());
    CAMLreturn(caml_copy_string(os.str().c_str()));
}

// The pfed interface
extern "C" CAMLprim value
stub_pfed_create(value size)
{
    CAMLparam1(size);
    CAMLlocal1(res);
    cindex_t dim = Int_val(size);
    res = caml_alloc_custom(&custom_ops_pfed, sizeof(pfed_wrap_t), 0, 1);
    new (Data_custom_val(res)) pfed_wrap_t(dim);
    CAMLreturn(res);
}

extern "C" CAMLprim value
stub_pfed_add_dbm(value v, value d)
{
    pfed_t & fed = *get_pfed_ptr(v);
    const pdbm_t & dbm = *get_pdbm_ptr(d);
    fed |= dbm;
    return Val_unit;
}

// returns true iff d is non-empty afterwards
bool
pdbm_intersect_dbm(pdbm_t & d, const dbm_t & dbm)
{
    cindex_t dim = d.getDimension();
    pdbm_t old_d = d;
    raw_t * rdbm = pdbm_getMutableMatrix(d, dim);
    // TODO d1 and d2 should be non-empty
    // dbm_intersection returns true iff d1 is non-empty afterwards
    bool res = dbm_intersection(rdbm, dbm.const_dbm(), dim);
    if (!res) return false;

    // otherwise
    // TODO should the pdbm be closed?
    // pdbm_close(d, dim);
    int32_t * offset = new int32_t[dim];
    offset[0] = 0;
    pdbm_getOffset(d, dim, offset);
    int32_t cost = pdbm_getCostOfValuation(old_d, dim, offset);
    pdbm_setCostAtOffset(d, dim, cost);
    delete[] offset;
    return true;
}

extern "C" CAMLprim value
stub_pfed_intersect_dbm(value v, value d)
{
    pfed_t & fed = *get_pfed_ptr(v);
    const dbm_t & dbm = *get_dbm_ptr(d);
    for (pfed_t::iterator it = fed.beginMutable(); it != fed.endMutable(); ++it)
    {
        bool not_empty = pdbm_intersect_dbm(*it, dbm);
        if (!not_empty)
        {
            it = fed.erase(it);
        }
    }
    return Val_unit;
}

extern "C" CAMLprim value
stub_pfed_hash(value v)
{
    CAMLparam1(v);
    CAMLreturn(Val_long(get_pfed_ptr(v)->hash(0)));
}

extern "C" CAMLprim value
stub_pfed_up(value v)
{
    CAMLparam1(v);
    get_pfed_ptr(v)->up(1);
    CAMLreturn(Val_unit);
}

extern "C" CAMLprim value
stub_pfed_update_value(value t, value c, value b)
{
    CAMLparam3(t,c,b);
    get_pfed_ptr(t)->updateValue(Int_val(c), Int_val(b));
    CAMLreturn(Val_unit);
}

extern "C" CAMLprim value
stub_pfed_iterator_get(value v)
{
    CAMLparam1(v);
    CAMLlocal1(res);
    res = caml_alloc_custom(&custom_ops_pdbm, sizeof(pdbm_wrap_t), 0, 1);
    new (Data_custom_val(res)) pdbm_wrap_t(**get_pfed_it_ptr(v));
    CAMLreturn(res);
}

extern "C" CAMLprim value
stub_pfed_iterator_notequal(value v1, value v2)
{
    CAMLparam2(v1, v2);
    const pfed_it_wrap_t & i1 = *get_pfed_it_ptr(v1);
    const pfed_it_wrap_t & i2 = *get_pfed_it_ptr(v2);
    CAMLreturn(Val_bool(i1 != i2));
}

extern "C" CAMLprim value
stub_pfed_iterator_incr(value v)
{
    CAMLparam1(v);
    pfed_it_wrap_t & i = *get_pfed_it_ptr(v);
    ++i;
    CAMLreturn(Val_unit);
}

extern "C" CAMLprim value
stub_pfed_iterator_begin(value v)
{
    CAMLparam1(v);
    CAMLlocal1(res);
    res = caml_alloc_custom(&custom_ops_pfed_it, sizeof(pfed_it_wrap_t), 0, 1);
    new (Data_custom_val(res)) pfed_it_wrap_t(get_pfed_ptr(v)->beginMutable());
    CAMLreturn(res);
}

extern "C" CAMLprim value
stub_pfed_iterator_end(value v)
{
    CAMLparam1(v);
    CAMLlocal1(res);
    res = caml_alloc_custom(&custom_ops_pfed_it, sizeof(pfed_it_wrap_t), 0, 1);
    new (Data_custom_val(res)) pfed_it_wrap_t(get_pfed_ptr(v)->endMutable());
    CAMLreturn(res);
}
