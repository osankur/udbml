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
stub_pdbm_set_empty(value v)
{
    CAMLparam1(v);
    get_pdbm_ptr(v)->~pdbm_t();
    new (Data_custom_val(v)) pdbm_wrap_t();
    CAMLreturn(Val_unit);
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
    // build the preorder on clocks wrt z
    explicit clock_po_t(const pdbm_t &z, const std::vector<int> &mbounds)
    {
        int dim = z.getDimension();

        for (int i = 1; i < dim; ++i)
        {
            if (z(i,0) <= dbm_boundbool2raw(mbounds[i], false))
            {
                _below.push_back(i);
                continue;
            }
            if (z(0,i) <= dbm_boundbool2raw(-mbounds[i], true))
            {
                // nothing to do, as clocks above their bound are not stored
                continue;
            }

            bool to_insert = true;
            for (auto it = _order.begin(); to_insert && it != _order.end(); ++it)
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
                _order.push_back(std::vector<std::vector<int> >(1, std::vector<int>(1,i)));
            }
        }
    }

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

    std::string
    toString() const
    {
        std::stringstream res;
        res << "Y = { ";
        for (int i = 0; i < _is_in.size(); ++i)
        {
            if (_is_in[i])
                res << i << ",";
        }
        res << "}";
        return res.str();
    }


private:
    const clock_po_t & _order;
    std::vector<int> _state; // same size as _order
    std::vector<bool> _is_in; // same size as dim
    cindex_t _size; // the size of the current subset
};

// forward declaration
bool
pdbm_intersect_dbm(pdbm_t & d, const raw_t * const dbm);

// check whether a priced zone dominates another
// this is based on the exploration of all subsets of clocks
bool
pdbm_square_inclusion_exp(const pdbm_t &z1, const pdbm_t &z2, const std::vector<int> &mbounds)
{
    // check whether all valuations in Z have an equivalent in Z'
    // it is the same as testing the inclusion of Z in the closure of Z'
    if (! dbm_closure_leq(z1.const_dbm(), z2.const_dbm(), z1.getDimension(), mbounds, mbounds))
        return false;

    // get the dimension
    cindex_t dim = z1.getDimension();

    // the clock order for the lhs zone
    clock_po_t preorder(z1, mbounds);
    // initialize Y to emptyset
    clock_po_iterator currentY(preorder, dim);

    // the main loop
    // TODO make it parallel, e.g. one thread for each Y
    do
    {
        // first build Z_Y and Z_Y'
        pdbm_t zy1 = z1;
        pdbm_t zy2 = z2;
        bool zy1_not_empty = true;
        for (int i = 1; zy1_not_empty && i != dim; ++i)
        {
            if (currentY.is_in(i))
            {
                zy1_not_empty = pdbm_constrain1(zy1, dim, i, 0, dbm_boundbool2raw(mbounds[i], false));
                pdbm_constrain1(zy2, dim, i, 0, dbm_boundbool2raw(mbounds[i], false));
            }
            else
            {
                zy1_not_empty = pdbm_constrain1(zy1, dim, 0, i, dbm_boundbool2raw(-mbounds[i], true));
                pdbm_constrain1(zy2, dim, 0, i, dbm_boundbool2raw(-mbounds[i], true));
            }
        }

        // if Z_Y is not empty
        if (zy1_not_empty)
        {
            // relax, i.e. take the topological closure
            pdbm_relax(zy1, dim);
            pdbm_relax(zy2, dim);

            // use feds for which all the needed operations are already implemented
            pfed_t facets1(zy1, dim), facets2(zy2, dim);

            // for each clock not in Y, do a projection
            for (cindex_t cl = 1; cl != dim; ++cl)
            {
                if (!currentY.is_in(cl))
                {
                    facets1.updateValue(cl, 0);
                    facets2.updateValue(cl, 0);
                }
            }

            // if a facet of Z is not subsumed by those of Z', we loose
            for (pdbm_t f1 : facets1)
            {
                fed_t cover(dim);
                for (pdbm_t f2 : facets2)
                {
                    pdbm_t f12 = f1;
                    if (pdbm_intersect_dbm(f12, f2.const_dbm()))
                    {
                        relation_t rel = pdbm_relation(f12, f2, dim);
                        if (rel == base_SUBSET || rel == base_EQUAL)
                        {
                            cover |= dbm_t(f12.const_dbm(), dim);
                        }
                    }
                }
                if (!cover.ge(f1.const_dbm(), dim))
                {
                    return false;
                }
            }
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

extern "C" CAMLprim value
stub_pfed_has(value t, value z)
{
    CAMLparam2(t, z);
    bool res = false;
    for (pfed_t::const_iterator it = get_pfed_ptr(t)->begin();
         it != get_pfed_ptr(t)->end(); ++it)
    {
        if (*it == *get_pdbm_ptr(z))
        {
            res = true;
            break;
        }
    }
    CAMLreturn(Val_bool(res));
}

// returns true iff d is non-empty afterwards
bool
pdbm_intersect_dbm(pdbm_t & d, const raw_t * const dbm)
{
    cindex_t dim = d.getDimension();
    pdbm_t old_d = d;
    raw_t * rdbm = pdbm_getMutableMatrix(d, dim);
    // TODO d1 and d2 should be non-empty
    // dbm_intersection returns true iff d1 is non-empty afterwards
    bool res = dbm_intersection(rdbm, dbm, dim);
    if (!res) return false;

    // otherwise
    // TODO should the pdbm be closed?
    // pdbm_close(d, dim);
    int32_t * offset = new int32_t[dim];
    offset[0] = 0;
    pdbm_getOffset(d, dim, offset);
    int32_t cost = pdbm_getCostOfVertex(old_d, dim, offset);
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
        bool not_empty = pdbm_intersect_dbm(*it, dbm.const_dbm());
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
    return Val_long(get_pfed_ptr(v)->hash(0));
}

extern "C" CAMLprim value
stub_pfed_is_empty(value v)
{
    return Val_bool(get_pfed_ptr(v)->isEmpty());
}

extern "C" CAMLprim value
stub_pfed_set_empty(value t)
{
    CAMLparam1(t);
    pfed_t * f = get_pfed_ptr(t);
    f->setEmpty();
    CAMLreturn(Val_unit);
}

extern "C" CAMLprim value
stub_pfed_up(value v, value rate)
{
    get_pfed_ptr(v)->up(Int_val(rate));
    return Val_unit;
}
extern "C" CAMLprim value
stub_pfed_dimension(value pf)
{
    return Val_int(get_pfed_ptr(pf)->getDimension());
}

extern "C" CAMLprim value
stub_pfed_update_value(value t, value c, value b)
{
    get_pfed_ptr(t)->updateValue(Int_val(c), Int_val(b));
    return Val_unit;
}

extern "C" CAMLprim value
stub_pfed_free_clock(value t, value cl)
{
	CAMLparam2(t,cl);
    get_pfed_ptr(t)->freeClock(Int_val(cl));
	CAMLreturn(Val_unit);
}

extern "C" CAMLprim value
stub_pdbm_free_clock(value t, value cl)
{
	CAMLparam2(t,cl);
    get_pdbm_ptr(t)->freeClock(Int_val(cl));
	CAMLreturn(Val_unit);
}


extern "C" CAMLprim value
stub_pfed_iterate(value t, value f)
{
    CAMLparam2(t,f);
    CAMLlocal1(z);

    for (pfed_t::iterator it = get_pfed_ptr(t)->beginMutable(); it != get_pfed_ptr(t)->endMutable(); ++it)
    {
        // initialize z
        z = caml_alloc_custom(&custom_ops_pdbm, sizeof(pdbm_wrap_t), 0, 1);
        new (Data_custom_val(z)) pdbm_wrap_t(*it);

        // callback
        caml_callback(f, z);
        // in case f has changed the value of z, update it with the new value of z
        *it = *get_pdbm_ptr(z);
    }

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
    const pfed_it_wrap_t & i1 = *get_pfed_it_ptr(v1);
    const pfed_it_wrap_t & i2 = *get_pfed_it_ptr(v2);
    return Val_bool(i1 != i2);
}

extern "C" CAMLprim value
stub_pfed_iterator_incr(value v)
{
    pfed_it_wrap_t & i = *get_pfed_it_ptr(v);
    ++i;
    return Val_unit;
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
