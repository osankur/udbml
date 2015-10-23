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
stub_pfed_intersect_dbm(value v, value d)
{
    pfed_t & fed = *get_pfed_ptr(v);
    const dbm_t & dbm = *get_dbm_ptr(d);
    for (pfed_t::iterator it = fed.beginMutable(); it != fed.endMutable(); ++it)
    {
        pdbm_t & d = *it;
        raw_t * rdbm = pdbm_getMutableMatrix(d, d.getDimension());
        // TODO d1 and d2 should be non-empty
        // dbm_intersection returns true iff d1 is non-empty afterwards
        dbm_intersection(rdbm, dbm.const_dbm(), d.getDimension());
        pdbm_close(d, d.getDimension());
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
