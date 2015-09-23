extern "C" {
#include <stdio.h>
#include <stdlib.h>
#include <caml/mlvalues.h>
#include <caml/memory.h>
#include <caml/alloc.h>
#include <caml/callback.h>
#include <caml/custom.h>
}
#include <dbm/constraints.h>
#include <dbm/fed.h>
#include <base/bitstring.h>
#include <iostream>
#include <sstream>
#include <vector>

#define get_dbm_ptr(x) ((dbm_wrap_t*)Data_custom_val(x))
#define get_fed_tp(x) ((fed_t*)((fed_wrap_t*)Data_custom_val(x))->d)
#define get_fed_it_tp(x) ((fed_t::iterator*)((fed_it_wrap_t*)Data_custom_val(x))->d)
#define get_bitvector_tp(x) ((BitString*)((bitvector_wrap_t*)Data_custom_val(x))->b)
using namespace dbm;
using namespace base;
/* This is the main type for Dbms.
		The encapsulation is needed since
		we free dbm_t* manually, while
		the dbm_wrap_t will be freed by the gc.
*/
typedef dbm_t dbm_wrap_t;
typedef struct fed_wrap_t {
	fed_t * d;
} fed_wrap_t;
typedef struct fed_it_wrap_t {
	fed_t::iterator * d;
} fed_it_wrap_t;

typedef struct bitvector_wrap_t {
	BitString * b;
} bitvector_wrap_t;


/*
	 Warning/Reminder:
	 The following functions are for the dbm_wrap_t which 
	 has Custom_tag in the Caml's internal representation.
	 No call to GC should occur here, so do not use
	 CAMLparam, CAMLreturn or any allocation function inside!
 */
extern "C" void finalize_dbm(value v){
    ((dbm_wrap_t *)Data_custom_val(v))->~dbm_t();
}
extern "C" void finalize_fed(value v){
	fed_wrap_t * d = (fed_wrap_t *)Data_custom_val(v);
	delete(d->d);
}
extern "C" void finalize_fed_it(value v){
	fed_it_wrap_t * d = (fed_it_wrap_t *)Data_custom_val(v);
	delete(d->d);
}
extern "C" void finalize_bitvector(value v){
	BitString * b = get_bitvector_tp(v);
	delete(b);
}


extern "C" int compare_dbm(value v1, value v2){
	const dbm_t * i1 = get_dbm_ptr(v1);
	const dbm_t * i2 = get_dbm_ptr(v2);
	if ( i1->sameAs(*i2) ) return 0;
	if ( i1 < i2 ) return -1;
	else if (i1 > i2) return 1;
	return 0;
}

extern "C" int compare_fed(value v1, value v2){
	fed_t * i1 = get_fed_tp(v1);
	fed_t * i2 = get_fed_tp(v2);
	if ( i1->sameAs(*i2) ) return 0;
	if ( i1 < i2 ) return -1;
	else if (i1 > i2) return 1;
	return 0;
}
extern "C" int compare_fed_it(value v1, value v2){
	fed_t::iterator * i1 = get_fed_it_tp(v1);
	fed_t::iterator * i2 = get_fed_it_tp(v2);
	if ( i1 < i2 ) return -1;
	else if (i1 > i2) return 1;
	return 0;
}
extern "C" int compare_bitvector(value v1, value v2){
	BitString * i1 = get_bitvector_tp(v1);
	BitString * i2 = get_bitvector_tp(v2);
	if ( i1 < i2 ) return -1;
	else if (i1 > i2) return 1;
	return 0;
}

extern "C" long hash_dbm(value v){
	return (long)get_dbm_ptr(v)->hash();
}
extern "C" long hash_fed(value v){
	fed_t * f = get_fed_tp(v);
	return (long)f->hash();
}
extern "C" long hash_bitvector(value v){
	return (long)get_bitvector_tp(v);
}
extern "C" long hash_fed_it(value v){
	return (long)get_fed_it_tp(v);
}

static struct custom_operations custom_ops_dbm = {
    .identifier     = (char*)"dbm_wrap_t handling",
    .finalize       = finalize_dbm,
    .compare        = compare_dbm,
    .hash           = hash_dbm,
    .serialize      = custom_serialize_default,
    .deserialize    = custom_deserialize_default
};
static struct custom_operations custom_ops_fed = {
    .identifier     = (char*)"fed_wrap_t handling",
    .finalize       = finalize_fed,
    .compare        = compare_fed,
    .hash           = hash_fed,
    .serialize      = custom_serialize_default,
    .deserialize    = custom_deserialize_default
};
static struct custom_operations custom_ops_fed_it = {
    .identifier     = (char*)"fed_it_wrap_t handling",
    .finalize       = finalize_fed_it,
    .compare        = compare_fed_it,
    .hash           = hash_fed_it,
    .serialize      = custom_serialize_default,
    .deserialize    = custom_deserialize_default
};

static struct custom_operations custom_ops_bitvector = {
    .identifier     = (char*)"bitvector_wrap_t handling",
    .finalize       = finalize_bitvector,
    .compare        = compare_bitvector,
    .hash           = hash_bitvector,
    .serialize      = custom_serialize_default,
    .deserialize    = custom_deserialize_default
};


// The bitvector interface
extern "C" CAMLprim value
stub_bitvector_create(value vsize)
{
	CAMLparam1(vsize);
	CAMLlocal1(bw);
	BitString * b = new BitString(Int_val(vsize));
	bw = caml_alloc_custom(&custom_ops_bitvector, sizeof(bitvector_wrap_t), 0, 1);
  ((bitvector_wrap_t*)Data_custom_val(bw))->b = b;
	CAMLreturn(bw);
}

extern "C" CAMLprim value
stub_bitvector_copy(value bw)
{
	CAMLparam1(bw);
	CAMLlocal1(bwnew);
	BitString * b = get_bitvector_tp(bw);
	BitString * bnew = new BitString(*b);
	bwnew = caml_alloc_custom(&custom_ops_bitvector, sizeof(bitvector_wrap_t), 0, 1);
  ((bitvector_wrap_t*)Data_custom_val(bwnew))->b = bnew;
	CAMLreturn(bw);
}

extern "C" CAMLprim value
stub_bitvector_assign_bit(value bw, value vind, value vv)
{
	CAMLparam3(bw, vind, vv);
	BitString * b = get_bitvector_tp(bw);
	bit_t bit = Bool_val(vv)? ONE : ZERO;
	b->assignBit(Int_val(vind), bit);
	printf("bitcount: %d\n", (int)b->count());
	CAMLreturn(Val_unit);
}
extern "C" CAMLprim value
stub_bitvector_get_bit(value bw, value vind)
{
	CAMLparam2(bw, vind);
	BitString * b = get_bitvector_tp(bw);
	CAMLreturn(Val_int(b->getBit(Bool_val(vind))));
}
extern "C" CAMLprim value
stub_bitvector_to_string(value bw)
{
	CAMLparam1(bw);
	BitString * b = get_bitvector_tp(bw);
	std::stringstream buf;
	b->prettyPrint(buf);
	CAMLreturn(caml_copy_string(buf.str().c_str()));
}

extern "C" CAMLprim value
stub_bitvector_count(value bw)
{
	CAMLparam1(bw);
	BitString * b = get_bitvector_tp(bw);
	CAMLreturn(Val_int(b->count()));
}

/// The C array interface
#define get_cvector(x) ((carray_t*)Data_custom_val(x))

typedef std::vector<int> carray_t;

extern "C" void finalize_carray(value v) {
    ((carray_t*)Data_custom_val(v))->~vector<int>();
}
extern "C" int compare_carray(value v1, value v2) {
    const std::vector<int> &vec1 = *get_cvector(v1);
    const std::vector<int> &vec2 = *get_cvector(v2);
    if (v1 < v2)
        return -1;
    else if (v1 > v2)
        return 1;
    return 0;
}
extern "C" long hash_carray(value v){
    const std::vector<int> &vv = *get_cvector(v);
    long res = 0;
    for (std::vector<int>::const_iterator it = vv.begin();
         it != vv.end(); ++it) {
        res ^= (*it) + 0x9e3779b9 + (res << 6) + (res >> 2);
    }
    return res;
}

static struct custom_operations custom_ops_carray = {
    .identifier     = (char*)"carray handling",
    .finalize       = finalize_carray,
    .compare        = compare_carray,
    .hash           = hash_carray,
    .serialize      = custom_serialize_default,
    .deserialize    = custom_deserialize_default,
};

extern "C" CAMLprim value
stub_carray_to_c(value v, value size)
{
    CAMLparam2(v, size);
    CAMLlocal1(res);
    int dim = Int_val(size);
    res = caml_alloc_custom(&custom_ops_carray, sizeof(carray_t), 0, 1);
    new (Data_custom_val(res)) carray_t(dim);
    std::vector<int> & d = *get_cvector(res);
    for (int i = 0; i < dim; ++i) {
        d[i] = Int_val(Field(v,i));
    }
    CAMLreturn(res);
}




// The dbm interface
extern "C" CAMLprim value
stub_dbm_create(value size)
{
	CAMLparam1(size);
	CAMLlocal1(dw);
	cindex_t dim = Int_val(size);
	dw = caml_alloc_custom(&custom_ops_dbm, sizeof(dbm_wrap_t), 0, 1);
    new (Data_custom_val(dw)) dbm_wrap_t(dim);
	CAMLreturn(dw);
}

extern "C" CAMLprim value
stub_dbm_copy(value dw)
{
	CAMLparam1(dw);
	CAMLlocal1(new_dw);
	const dbm_t & d = *get_dbm_ptr(dw);
	new_dw = caml_alloc_custom(&custom_ops_dbm, sizeof(dbm_wrap_t), 0, 1);
    new (Data_custom_val(new_dw)) dbm_wrap_t(d);
	CAMLreturn(new_dw);
}

extern "C" CAMLprim value
stub_dbm_dimension(value t)
{
    return Val_int(get_dbm_ptr(t)->getDimension());
}

extern "C" CAMLprim value
stub_dbm_intern(value t)
{
	CAMLparam1(t);
    get_dbm_ptr(t)->intern();
	CAMLreturn(Val_unit);
}
extern "C" CAMLprim value
stub_dbm_is_empty(value t)
{
	CAMLparam1(t);
    int ret = get_dbm_ptr(t)->isEmpty();
	CAMLreturn(Val_bool(ret));
}

extern "C" CAMLprim value
stub_dbm_has_zero(value t)
{
	CAMLparam1(t);
    int ret = get_dbm_ptr(t)->hasZero();
	CAMLreturn(Val_bool(ret));
}
extern "C" CAMLprim value
stub_dbm_hash(value t)
{
	CAMLparam1(t);
	const dbm_t * d = get_dbm_ptr(t);
	unsigned ret = d->hash();
	CAMLreturn(Val_int(ret));
}

extern "C" CAMLprim value
stub_dbm_set_init(value t)
{
    get_dbm_ptr(t)->setInit();
	return Val_unit;
}

extern "C"
CAMLprim value
stub_get_infty(value unit){
	return Val_int(dbm_INFINITY);
}


extern "C" CAMLprim value
stub_dbm_set_zero(value t)
{
	CAMLparam1(t);
    get_dbm_ptr(t)->setZero();
	CAMLreturn(Val_unit);
}
extern "C" CAMLprim value
stub_dbm_is_zero(value t)
{
	CAMLparam1(t);
	const dbm_t * d = get_dbm_ptr(t);
	CAMLreturn(Val_bool(d->isZero()));
}
extern "C" CAMLprim value
stub_dbm_is_init(value t)
{
	CAMLparam1(t);
	const dbm_t * d = get_dbm_ptr(t);
	CAMLreturn(Val_bool(d->isInit()));
}

/*
CAMLprim value
make_constraint(int bound, int ineq)
{
	CAMLlocal1(ret);
	ret = caml_alloc(2,0);
	Store_field(ret, 0, Val_int(bound));
	Store_field(ret, 1, Val_int(ineq));
	CAMLreturn(ret);
}
*/

extern "C" CAMLprim value
stub_dbm_at(value t, value i, value j){
	CAMLparam3(t,i,j);
	CAMLlocal1(ret);
	const dbm_t & d = *get_dbm_ptr(t);
	raw_t r = d(Int_val(i),Int_val(j));
	printf("Getting value: %d\n", dbm_raw2bound(r));
	// Make pair
	ret = caml_alloc(2,0);
	Store_field(ret, 0, Val_int(dbm_raw2bound(r)));
	Store_field(ret, 1, Val_int(dbm_raw2strict(r)));
	CAMLreturn(ret);
}
extern "C" CAMLprim value
stub_dbm_at_bound(value t, value i, value j){
	const dbm_t & d = *get_dbm_ptr(t);
	return Val_int(dbm_raw2bound(d(Int_val(i),Int_val(j))));
}

extern "C" CAMLprim value
stub_dbm_equal(value t, value u)
{
	CAMLparam2(t,u);
	const dbm_t & dt = *get_dbm_ptr(t);
	const dbm_t & du = *get_dbm_ptr(u);
	int ret = (dt == du);
	CAMLreturn(Val_bool(ret));
}
extern "C" CAMLprim value
stub_dbm_notequal(value t, value u)
{
	CAMLparam2(t,u);
    const dbm_t & dt = *get_dbm_ptr(t);
    const dbm_t & du = *get_dbm_ptr(u);
	int ret = (dt != du);
	CAMLreturn(Val_bool(ret));
}
extern "C" CAMLprim value
stub_dbm_lt(value t, value u)
{
	CAMLparam2(t,u);
    const dbm_t & dt = *get_dbm_ptr(t);
    const dbm_t & du = *get_dbm_ptr(u);
	int ret = (dt < du);
	CAMLreturn(Val_bool(ret));
}
extern "C" CAMLprim value
stub_dbm_gt(value t, value u)
{
	CAMLparam2(t,u);
    const dbm_t & dt = *get_dbm_ptr(t);
    const dbm_t & du = *get_dbm_ptr(u);
	int ret = (dt > du);
	CAMLreturn(Val_bool(ret));
}
extern "C" CAMLprim value
stub_dbm_leq(value t, value u)
{
	CAMLparam2(t,u);
    const dbm_t & dt = *get_dbm_ptr(t);
    const dbm_t & du = *get_dbm_ptr(u);
	int ret = (dt <= du);
	CAMLreturn(Val_bool(ret));
}
extern "C" CAMLprim value
stub_dbm_geq(value t, value u)
{
	CAMLparam2(t,u);
    const dbm_t & dt = *get_dbm_ptr(t);
    const dbm_t & du = *get_dbm_ptr(u);
	int ret = (dt >= du);
	CAMLreturn(Val_bool(ret));
}

extern "C" CAMLprim value
stub_dbm_constrain(value t, value ct)
{
	CAMLparam2(t,ct);
	CAMLlocal4(i,j,b,ineq);
	// a constraint is (i,j,(b,ineq))
	i = Field(ct,0);
	j = Field(ct,1);
	b = Field(Field(ct,2),0);
	ineq = Field(Field(ct,2),1);
	dbm_t * d = get_dbm_ptr(t);
	raw_t r = dbm_boundbool2raw(Int_val(b), Int_val(ineq) == 0);
	d->constrain(Int_val(i), Int_val(j), r);
	CAMLreturn(Val_unit);
}

extern "C" CAMLprim value
stub_dbm_copy_to(value t)
{
	CAMLparam1(t);
	CAMLlocal2(ar,tmp);
	const dbm_t & d = *get_dbm_ptr(t);
	int dim = d.getDimension();
	// Create and initialize the array
	// before the next call to caml_alloc
	ar = caml_alloc(dim*dim,0);
	for (int i = 0; i < dim; i++){
		for (int j = 0; j < dim; j++){
			Store_field(ar, i*dim + j, Val_unit);
		}
	}
	// Fill in the array
	for (int i = 0; i < dim; i++){
		for (int j = 0; j < dim; j++){
			tmp = caml_alloc(2,0);
			Store_field(tmp, 0, Val_int(dbm_raw2bound(d(i,j))));
			Store_field(tmp, 1, Val_int(dbm_raw2strict(d(i,j))));
			Store_field(ar, i*dim + j, tmp);
		}
	}
	CAMLreturn(ar);
}

// CHECK all semantics operations, intersect, up, down, free..
extern "C" CAMLprim value
stub_dbm_intersect(value t, value u)
{
    dbm_t & dt = *get_dbm_ptr(t);
    const dbm_t & du = *get_dbm_ptr(u);
	dt &= du;
	return Val_unit;
}
extern "C" CAMLprim value
stub_dbm_intersects(value t, value u)
{
	CAMLparam2(t,u);
    const dbm_t & dt = *get_dbm_ptr(t);
    const dbm_t & du = *get_dbm_ptr(u);
	CAMLreturn(Val_bool(dt.intersects(du)));
}

extern "C" CAMLprim value
stub_dbm_up(value t)
{
	CAMLparam1(t);
	dbm_t * d = get_dbm_ptr(t);
	d->up();
	CAMLreturn(Val_unit);
}
extern "C" CAMLprim value
stub_dbm_down(value t)
{
	CAMLparam1(t);
	dbm_t * d = get_dbm_ptr(t);
	d->down();
	CAMLreturn(Val_unit);
}

extern "C" CAMLprim value
stub_dbm_free_clock(value t, value cl)
{
	CAMLparam2(t,cl);
	dbm_t * d = get_dbm_ptr(t);
	cindex_t icl = Int_val(cl);
	d->freeClock(icl);
	CAMLreturn(Val_unit);
}

extern "C" CAMLprim value
stub_dbm_free_up(value t, value cl)
{
	CAMLparam2(t,cl);
	dbm_t * d = get_dbm_ptr(t);
	cindex_t icl = Int_val(cl);
	d->freeUp(icl);
	CAMLreturn(Val_unit);
}
extern "C" CAMLprim value
stub_dbm_free_down(value t, value cl)
{
	CAMLparam2(t,cl);
	dbm_t * d = get_dbm_ptr(t);
	cindex_t icl = Int_val(cl);
	d->freeDown(icl);
	CAMLreturn(Val_unit);
}

extern "C" CAMLprim value
stub_dbm_free_all_up(value t)
{
	CAMLparam1(t);
	dbm_t * d = get_dbm_ptr(t);
	d->freeAllUp();
	CAMLreturn(Val_unit);
}
extern "C" CAMLprim value
stub_dbm_free_all_down(value t)
{
	CAMLparam1(t);
	dbm_t * d = get_dbm_ptr(t);
	d->freeAllDown();
	CAMLreturn(Val_unit);
}

extern "C" CAMLprim value
stub_dbm_is_unbounded(value t)
{
	CAMLparam1(t);
	const dbm_t * d = get_dbm_ptr(t);
	int ret = d->isUnbounded();
	CAMLreturn(Val_bool(ret));
}

// CHECK
extern "C" CAMLprim value
stub_dbm_copy_from(value ar, value vdim)
{
	CAMLparam1(ar);
    CAMLlocal1(dw);
	int dim = Int_val(vdim);
	dbm_t d = dbm_t(dim);
	raw_t * r = new raw_t[dim*dim];
	int b, ineq;
	for(int i = 0; i < dim*dim; i++){
		b = Int_val(Field(Field(ar, i), 0));
		ineq = Int_val(Field(Field(ar,i),1));
		r[i] = dbm_boundbool2raw(b,ineq==0);
	}
    dw = caml_alloc_custom(&custom_ops_dbm, sizeof(dbm_wrap_t), 0, 1);
    new (Data_custom_val(dw)) dbm_wrap_t(dim);
    get_dbm_ptr(dw)->copyFrom(r, dim);
	CAMLreturn(dw);
}

extern "C" CAMLprim value
stub_dbm_convex_add(value t, value u)
{
	CAMLparam2(t,u);
	dbm_t & dt = *get_dbm_ptr(t);
	const dbm_t & du = *get_dbm_ptr(u);
	dt += du;
	CAMLreturn(Val_unit);
}

extern "C" CAMLprim value
stub_dbm_update_value(value t, value vx, value vb)
{
	CAMLparam3(t,vx,vb);
	dbm_t * d = get_dbm_ptr(t);
	d->updateValue(Int_val(vx), Int_val(vb));
	CAMLreturn(Val_unit);
}

extern "C" CAMLprim value
stub_dbm_update_clock(value t, value vx, value vy)
{
	CAMLparam3(t,vx,vy);
	dbm_t * d = get_dbm_ptr(t);
	d->updateClock(Int_val(vx), Int_val(vy));
	CAMLreturn(Val_unit);
}
extern "C" CAMLprim value
stub_dbm_update_increment(value t, value vx, value vb)
{
	CAMLparam3(t,vx,vb);
	dbm_t * d = get_dbm_ptr(t);
	d->updateIncrement(Int_val(vx), Int_val(vb));
	CAMLreturn(Val_unit);
}
extern "C" CAMLprim value
stub_dbm_update(value t, value vx, value vy, value vb)
{
	CAMLparam4(t,vx,vy,vb);
	dbm_t * d = get_dbm_ptr(t);
	d->update(Int_val(vx), Int_val(vy), Int_val(vb));
	CAMLreturn(Val_unit);
}

extern "C" CAMLprim value
stub_dbm_satisfies(value t, value ct)
{
	CAMLparam2(t,ct);
	// a constraint is (i,j,(b,ineq))
	int i = Field(ct, 0);
	int j = Field(ct, 1);
	raw_t r = dbm_boundbool2raw(Field(Field(ct,2),0),
															Field(Field(ct,2),1) ==0 );
	const dbm_t * d = get_dbm_ptr(t);
	CAMLreturn(Val_bool(d->satisfies(i,j,r)));
}

extern "C" bool
dbm_closure_leq(const dbm_t &d1, const dbm_t &d2,
                const std::vector<int> &lbounds, const std::vector<int> &ubounds)
{
    // forall x,y
    //      Z_{x,0} < (<=, -U_x)
    // or   Z'_{x,y} + (<, -L_y) >= Z_{x,0}
    // or   Z'_{x,y} >= Z_{x,y}
    int dim = d1.getDimension();
    const raw_t * dr1 = d1.const_dbm();
    const raw_t * dr2 = d2.const_dbm();
    int n = 0;
    while (n < dim*dim-1) {
        const raw_t & zx0 = dr1[n];
        if (zx0 >= dbm_boundbool2raw(-ubounds[n/dim], false)) {
            const raw_t & zpxy = dr2[n];
            if (zpxy < dr1[n] && dbm_addRawRaw(zpxy, dbm_boundbool2raw(-lbounds[n % dim], true)) < zx0) {
                return false;
            } else {
                ++n;
            }
        } else {
            n = ((n/dim)+1)*dim;
        }
    }
    return true;
}

extern "C" CAMLprim value
stub_dbm_closure_leq(value vlbounds, value vubounds, value t1, value t2)
{
    const dbm_t & d1 = *get_dbm_ptr(t1);
    const dbm_t & d2 = *get_dbm_ptr(t2);
    int dim = d1.getDimension();
    assert(dim == d2.getDimension());
    return Val_bool(dbm_closure_leq(d1, d2, *get_cvector(vlbounds), *get_cvector(vubounds)));
}

extern "C" CAMLprim value
stub_dbm_extrapolate_max_bounds(value t, value vbounds)
{
	CAMLparam2(t,vbounds);
	dbm_t * d = get_dbm_ptr(t);
	int * bounds = get_cvector(vbounds)->data();
	d->extrapolateMaxBounds(bounds);
	CAMLreturn(Val_unit);
}
extern "C" CAMLprim value
stub_dbm_diagonal_extrapolate_max_bounds(value t, value vbounds)
{
	CAMLparam2(t,vbounds);
	dbm_t * d = get_dbm_ptr(t);
    int * bounds = get_cvector(vbounds)->data();
	d->diagonalExtrapolateMaxBounds(bounds);
	CAMLreturn(Val_unit);
}

extern "C" CAMLprim value
stub_dbm_extrapolate_lu_bounds(value t, value vlbounds, value vubounds)
{
	dbm_t * d = get_dbm_ptr(t);
	int * lbounds = get_cvector(vlbounds)->data();
	int * ubounds = get_cvector(vubounds)->data();
	d->extrapolateLUBounds(lbounds, ubounds);
	return Val_unit;
}

extern "C" CAMLprim value
stub_dbm_diagonal_extrapolate_lu_bounds(value t, value vlbounds, value vubounds)
{
	CAMLparam3(t,vlbounds, vubounds);
	dbm_t * d = get_dbm_ptr(t);
    int * lbounds = get_cvector(vlbounds)->data();
    int * ubounds = get_cvector(vubounds)->data();
	d->diagonalExtrapolateLUBounds(lbounds, ubounds);
	CAMLreturn(Val_unit);
}



extern "C" CAMLprim value
stub_dbm_is_subtraction_empty(value t, value u)
{
	CAMLparam2(t,u);
	const dbm_t & dt = *get_dbm_ptr(t);
	const dbm_t & du = *get_dbm_ptr(u);
	CAMLreturn(Val_bool(dt.isSubtractionEmpty(du)));
}

extern "C" CAMLprim value
stub_dbm_resize(value t, value vbvsrc, value vbvdst)
{
	CAMLparam3(t,vbvsrc, vbvdst);
	CAMLlocal1(retar);
	dbm_t * dt = get_dbm_ptr(t);
	BitString & src = *get_bitvector_tp(vbvsrc);
	BitString & dst = *get_bitvector_tp(vbvdst);
	size_t bitSize = src.intSize();
	cindex_t * table = new cindex_t[bitSize*32];
	memset(table, 0, bitSize*32);
	dt->resize(src(), dst(), bitSize, table);
	retar = caml_alloc(bitSize*32, 0);
	for(int i = 0; i < bitSize*32; i++){
		Field(retar, i) = Val_int(table[i]);
	}
	delete(table);
	CAMLreturn(retar);
}
 
extern "C" CAMLprim value
stub_dbm_max_dim(value v){
	return Val_int(dbm_t::MAX_DIM);
}
extern "C" CAMLprim value
stub_dbm_max_dim_power(value v){
	return Val_int(dbm_t::MAX_DIM_POWER);
}

extern "C" CAMLprim value
stub_dbm__internal_addr(value t)
{
	CAMLparam1(t);
	dbm_t * d = get_dbm_ptr(t);
	CAMLreturn(Val_long((long)d));
}


// The Fed interface
extern "C" CAMLprim value
stub_fed_create(value size)
{
	CAMLparam1(size);
	CAMLlocal1(fw);
	cindex_t dim = Int_val(size);
	fed_t * f = new fed_t(dim);
	fw = caml_alloc_custom(&custom_ops_fed, sizeof(fed_wrap_t), 0, 1);
  ((fed_wrap_t*)Data_custom_val(fw))->d = f;
	CAMLreturn(fw);
}

extern "C" CAMLprim value
stub_fed_copy(value fw)
{
	CAMLparam1(fw);
	CAMLlocal1(new_fw);
	fed_t * f = get_fed_tp(fw);
	fed_t * new_f = new fed_t(*f);
	new_fw = caml_alloc_custom(&custom_ops_fed, sizeof(fed_wrap_t), 0, 1);
  ((fed_wrap_t*)Data_custom_val(new_fw))->d = new_f;
	CAMLreturn(new_fw);
}

extern "C" CAMLprim value
stub_fed_dimension(value t)
{
	CAMLparam1(t);
	CAMLreturn(Val_int(get_fed_tp(t)->getDimension()));
}
extern "C" CAMLprim value
stub_fed_intern(value t)
{
	CAMLparam1(t);
	fed_t * f = get_fed_tp(t);
	f->intern();
	CAMLreturn(Val_unit);
}
extern "C" CAMLprim value
stub_fed_is_empty(value t)
{
	CAMLparam1(t);
	fed_t * f = get_fed_tp(t);
	int ret = f->isEmpty();
	CAMLreturn(Val_bool(ret));
}

extern "C" CAMLprim value
stub_fed_has_zero(value t)
{
	CAMLparam1(t);
	fed_t * f = get_fed_tp(t);
	int ret = f->hasZero();
	CAMLreturn(Val_bool(ret));
}
extern "C" CAMLprim value
stub_fed_hash(value t)
{
	CAMLparam1(t);
	fed_t * f = get_fed_tp(t);
	unsigned ret = f->hash();
	CAMLreturn(Val_int(ret));
}

extern "C" CAMLprim value
stub_fed_approx_equal(value t, value u)
{
	CAMLparam2(t,u);
	fed_t * dt = get_fed_tp(t);
	fed_t * du = get_fed_tp(u);
	CAMLreturn(Val_bool((*dt) == (*du)));
}
extern "C" CAMLprim value
stub_fed_approx_notequal(value t, value u)
{
	CAMLparam2(t,u);
	fed_t * dt = get_fed_tp(t);
	fed_t * du = get_fed_tp(u);
	CAMLreturn(Val_bool((*dt) != (*du)));
}

extern "C" CAMLprim value
stub_fed_approx_lt(value t, value u)
{
	CAMLparam2(t,u);
	fed_t * dt = get_fed_tp(t);
	fed_t * du = get_fed_tp(u);
	CAMLreturn(Val_bool((*dt) < (*du)));
}
extern "C" CAMLprim value
stub_fed_approx_leq(value t, value u)
{
	CAMLparam2(t,u);
	fed_t * dt = get_fed_tp(t);
	fed_t * du = get_fed_tp(u);
	CAMLreturn(Val_bool((*dt) <= (*du)));
}
extern "C" CAMLprim value
stub_fed_approx_gt(value t, value u)
{
	CAMLparam2(t,u);
	fed_t * dt = get_fed_tp(t);
	fed_t * du = get_fed_tp(u);
	CAMLreturn(Val_bool((*dt) > (*du)));
}
extern "C" CAMLprim value
stub_fed_approx_geq(value t, value u)
{
	CAMLparam2(t,u);
	fed_t * dt = get_fed_tp(t);
	fed_t * du = get_fed_tp(u);
	CAMLreturn(Val_bool((*dt) >= (*du)));
}

extern "C" CAMLprim value
stub_fed_approx_equal_dbm(value t, value u)
{
	CAMLparam2(t,u);
	fed_t * dt = get_fed_tp(t);
	const dbm_t & du = *get_dbm_ptr(u);
	CAMLreturn(Val_bool((*dt) == du));
}
extern "C" CAMLprim value
stub_fed_approx_notequal_dbm(value t, value u)
{
	CAMLparam2(t,u);
	fed_t * dt = get_fed_tp(t);
	const dbm_t & du = *get_dbm_ptr(u);
	CAMLreturn(Val_bool((*dt) != du));
}
extern "C" CAMLprim value
stub_fed_approx_lt_dbm(value t, value u)
{
	CAMLparam2(t,u);
	fed_t * dt = get_fed_tp(t);
	const dbm_t & du = *get_dbm_ptr(u);
	CAMLreturn(Val_bool((*dt) < du));
}

extern "C" CAMLprim value
stub_fed_approx_gt_dbm(value t, value u)
{
	CAMLparam2(t,u);
	fed_t * dt = get_fed_tp(t);
	const dbm_t & du = *get_dbm_ptr(u);
	CAMLreturn(Val_bool((*dt) > du));
}

extern "C" CAMLprim value
stub_fed_approx_geq_dbm(value t, value u)
{
	CAMLparam2(t,u);
	fed_t * dt = get_fed_tp(t);
	const dbm_t & du = *get_dbm_ptr(u);
	CAMLreturn(Val_bool((*dt) >= du));
}

extern "C" CAMLprim value
stub_fed_exact_equal(value t, value u)
{
	CAMLparam2(t,u);
	fed_t * dt = get_fed_tp(t);
	fed_t * du = get_fed_tp(u);
	CAMLreturn(Val_bool(dt->eq(*du)));
}
extern "C" CAMLprim value
stub_fed_exact_lt(value t, value u)
{
	CAMLparam2(t,u);
	fed_t * dt = get_fed_tp(t);
	fed_t * du = get_fed_tp(u);
	CAMLreturn(Val_bool(dt->lt(*du)));
}
extern "C" CAMLprim value
stub_fed_exact_gt(value t, value u)
{
	CAMLparam2(t,u);
	fed_t * dt = get_fed_tp(t);
	fed_t * du = get_fed_tp(u);
	CAMLreturn(Val_bool(dt->gt(*du)));
}
extern "C" CAMLprim value
stub_fed_exact_leq(value t, value u)
{
	CAMLparam2(t,u);
	fed_t * dt = get_fed_tp(t);
	fed_t * du = get_fed_tp(u);
	CAMLreturn(Val_bool(dt->le(*du)));
}
extern "C" CAMLprim value
stub_fed_exact_geq(value t, value u)
{
	CAMLparam2(t,u);
	fed_t * dt = get_fed_tp(t);
	fed_t * du = get_fed_tp(u);
	CAMLreturn(Val_bool(dt->ge(*du)));
}

extern "C" CAMLprim value
stub_fed_exact_equal_dbm(value t, value u)
{
	CAMLparam2(t,u);
	fed_t * dt = get_fed_tp(t);
	const dbm_t & du = *get_dbm_ptr(u);
	CAMLreturn(Val_bool(dt->ge(du)));
}
extern "C" CAMLprim value
stub_fed_exact_leq_dbm(value t, value u)
{
	CAMLparam2(t,u);
	fed_t * dt = get_fed_tp(t);
	const dbm_t & du = *get_dbm_ptr(u);
	CAMLreturn(Val_bool(dt->le(du)));
}
extern "C" CAMLprim value
stub_fed_exact_lt_dbm(value t, value u)
{
	CAMLparam2(t,u);
	fed_t * dt = get_fed_tp(t);
	const dbm_t & du = *get_dbm_ptr(u);
	CAMLreturn(Val_bool(dt->lt(du)));
}
extern "C" CAMLprim value
stub_fed_exact_geq_dbm(value t, value u)
{
	CAMLparam2(t,u);
	fed_t * dt = get_fed_tp(t);
	const dbm_t & du = *get_dbm_ptr(u);
	CAMLreturn(Val_bool(dt->ge(du)));
}
extern "C" CAMLprim value
stub_fed_exact_gt_dbm(value t, value u)
{
	CAMLparam2(t,u);
	fed_t * dt = get_fed_tp(t);
	const dbm_t & du = *get_dbm_ptr(u);
	CAMLreturn(Val_bool(dt->gt(du)));
}
extern "C" CAMLprim value
stub_fed_set_init(value t)
{
	CAMLparam1(t);
	fed_t * d = get_fed_tp(t);
	d->setInit();
	CAMLreturn(Val_unit);
}
extern "C" CAMLprim value
stub_fed_convex_hull(value t)
{
	CAMLparam1(t);
	fed_t * d = get_fed_tp(t);
	d->convexHull();
	CAMLreturn(Val_unit);
}
extern "C" CAMLprim value
stub_fed_union(value t, value u)
{
	CAMLparam2(t,u);
	fed_t * dt = get_fed_tp(t);
	fed_t * du = get_fed_tp(u);
	(*dt) |= *du;
	CAMLreturn(Val_unit);
}
extern "C" CAMLprim value
stub_fed_add(value t, value u)
{
	CAMLparam2(t,u);
	fed_t * dt = get_fed_tp(t);
	fed_t * du = get_fed_tp(u);
	dt->add(*du);
	CAMLreturn(Val_unit);
}
extern "C" CAMLprim value
stub_fed_add_dbm(value t, value u)
{
	CAMLparam2(t,u);
	fed_t * dt = get_fed_tp(t);
	const dbm_t & du = *get_dbm_ptr(u);
	dt->add(du);
	CAMLreturn(Val_unit);
}


extern "C" CAMLprim value
stub_fed_append(value t, value u)
{
	CAMLparam2(t,u);
	fed_t * dt = get_fed_tp(t);
	fed_t * du = get_fed_tp(u);
	dt->append(*du);
	CAMLreturn(Val_unit);
}
extern "C" CAMLprim value
stub_fed_append_end(value t, value u)
{
	CAMLparam2(t,u);
	fed_t * dt = get_fed_tp(t);
	fed_t * du = get_fed_tp(u);
	dt->appendEnd(*du);
	CAMLreturn(Val_unit);
}
extern "C" CAMLprim value
stub_fed_append_begin(value t, value u)
{
	CAMLparam2(t,u);
	fed_t * dt = get_fed_tp(t);
	fed_t * du = get_fed_tp(u);
	dt->appendBegin(*du);
	CAMLreturn(Val_unit);
}
extern "C" CAMLprim value
stub_fed_steal(value t, value u)
{
	CAMLparam2(t,u);
	fed_t * dt = get_fed_tp(t);
	fed_t * du = get_fed_tp(u);
	dt->steal(*du);
	CAMLreturn(Val_unit);
}
extern "C" CAMLprim value
stub_fed_convex_union(value t, value u)
{
	CAMLparam2(t,u);
	fed_t * dt = get_fed_tp(t);
	fed_t * du = get_fed_tp(u);
	(*dt) += (*du);
	CAMLreturn(Val_unit);
}
extern "C" CAMLprim value
stub_fed_convex_union_dbm(value t, value u)
{
	CAMLparam2(t,u);
	fed_t * dt = get_fed_tp(t);
	const dbm_t & du = *get_dbm_ptr(u);
	(*dt) += du;
	CAMLreturn(Val_unit);
}
extern "C" CAMLprim value
stub_fed_intersect(value t, value u)
{
	CAMLparam2(t,u);
	fed_t * dt = get_fed_tp(t);
	fed_t * du = get_fed_tp(u);
	(*dt) &= (*du);
	CAMLreturn(Val_unit);
}
extern "C" CAMLprim value
stub_fed_intersect_dbm(value t, value u)
{
	CAMLparam2(t,u);
	fed_t * dt = get_fed_tp(t);
	const dbm_t & du = *get_dbm_ptr(u);
	(*dt) &= du;
	CAMLreturn(Val_unit);
}
extern "C" CAMLprim value
stub_fed_subtract(value t, value u)
{
	CAMLparam2(t,u);
	fed_t * dt = get_fed_tp(t);
	fed_t * du = get_fed_tp(u);
	(*dt) -= (*du);
	CAMLreturn(Val_unit);
}
extern "C" CAMLprim value
stub_fed_subtract_dbm(value t, value u)
{
	CAMLparam2(t,u);
	fed_t * dt = get_fed_tp(t);
	const dbm_t & du = *get_dbm_ptr(u);
	(*dt) -= du;
	CAMLreturn(Val_unit);
}
extern "C" CAMLprim value
stub_fed_up(value t)
{
	CAMLparam1(t);
	fed_t * d = get_fed_tp(t);
	d->up();
	CAMLreturn(Val_unit);
}
extern "C" CAMLprim value
stub_fed_down(value t)
{
	CAMLparam1(t);
	fed_t * d = get_fed_tp(t);
	d->down();
	CAMLreturn(Val_unit);
}
extern "C" CAMLprim value
stub_fed_free_clock(value t, value cl)
{
	CAMLparam2(t,cl);
	fed_t * d = get_fed_tp(t);
	cindex_t icl = Int_val(cl);
	d->freeClock(icl);
	CAMLreturn(Val_unit);
}
extern "C" CAMLprim value
stub_fed_constrain(value t, value ct)
{
	CAMLparam2(t,ct);
	CAMLlocal4(i,j,b,ineq);
	// a constraint is (i,j,(b,ineq))
	i = Field(ct,0);
	j = Field(ct,1);
	b = Field(Field(ct,2),0);
	ineq = Field(Field(ct,2),1);
	fed_t * d = get_fed_tp(t);
	raw_t r = dbm_boundbool2raw(Int_val(b), Int_val(ineq) == 0);
	d->constrain(Int_val(i), Int_val(j), r);
	CAMLreturn(Val_unit);
}

extern "C" CAMLprim value
stub_fed_reduce(value t)
{
	CAMLparam1(t);
	fed_t * d = get_fed_tp(t);
	d->reduce();
	CAMLreturn(Val_unit);
}
extern "C" CAMLprim value
stub_fed_expensive_reduce(value t)
{
	CAMLparam1(t);
	fed_t * d = get_fed_tp(t);
	d->expensiveReduce();
	CAMLreturn(Val_unit);
}
extern "C" CAMLprim value
stub_fed_merge_reduce(value t, value vskip, value vlevel)
{
	CAMLparam1(t);
	fed_t * d = get_fed_tp(t);
	d->mergeReduce(Int_val(vskip), Int_val(vlevel));
	CAMLreturn(Val_unit);
}
extern "C" CAMLprim value
stub_fed_convex_reduce(value t)
{
	CAMLparam1(t);
	fed_t * d = get_fed_tp(t);
	d->convexReduce();
	CAMLreturn(Val_unit);
}
extern "C" CAMLprim value
stub_fed_expensive_convex_reduce(value t)
{
	CAMLparam1(t);
	fed_t * d = get_fed_tp(t);
	d->expensiveConvexReduce();
	CAMLreturn(Val_unit);
}
extern "C" CAMLprim value
stub_fed_partition_reduce(value t)
{
	CAMLparam1(t);
	fed_t * d = get_fed_tp(t);
	d->partitionReduce();
	CAMLreturn(Val_unit);
}

extern "C" CAMLprim value
stub_fed_predt(value t, value vbad, value vrestrict)
{
	CAMLparam1(t);
	fed_t * d = get_fed_tp(t);
	fed_t * bad = get_fed_tp(vbad);
	const raw_t * restrict = (Val_int(vrestrict) == 0)? NULL : (get_dbm_ptr(Field(vrestrict,0)))->const_dbm();
	d->predt(*bad, restrict);
	CAMLreturn(Val_unit);
}
extern "C" CAMLprim value
stub_fed__internal_addr(value t)
{
	CAMLparam1(t);
	fed_t * d = get_fed_tp(t);
	CAMLreturn(Val_long((long)d));
}

extern "C" CAMLprim value
stub_fed_is_included_in_predt(value t, value vgood, value vbad)
{
	CAMLparam1(t);
	fed_t * d = get_fed_tp(t);
	fed_t * bad = get_fed_tp(vbad);
	fed_t * good = get_fed_tp(vgood);
	CAMLreturn(Val_bool(d->isIncludedInPredt(*good, *bad)));
}
extern "C" CAMLprim value
stub_fed_begin_it(value t)
{
	CAMLparam1(t);
	CAMLlocal1(fitw);
	fed_t * d = get_fed_tp(t);
	fitw = caml_alloc_custom(&custom_ops_fed_it, sizeof(fed_it_wrap_t), 0, 1);
	fed_t::iterator * it = new fed_t::iterator();
	*it = d->beginMutable();
  ((fed_it_wrap_t*)Data_custom_val(fitw))->d = it;
	CAMLreturn(fitw);
}

// Fed.Iterator interface
extern "C" CAMLprim value
stub_fed_iterator_get(value t)
{
	CAMLparam1(t);
	CAMLlocal1(ret);
	fed_t::iterator * it = get_fed_it_tp(t);
	ret = caml_alloc_custom(&custom_ops_dbm, sizeof(dbm_wrap_t), 0, 1);
    new (Data_custom_val(ret)) dbm_wrap_t(*(it->operator->()));
	CAMLreturn(ret);
}

extern "C" CAMLprim value
stub_fed_iterator_incr(value t)
{
	CAMLparam1(t);
	fed_t::iterator * it = get_fed_it_tp(t);
	(*it).operator++();
	CAMLreturn(Val_unit);
}

extern "C" CAMLprim value
stub_fed_iterator_is_null(value t)
{
	CAMLparam1(t);
	fed_t::iterator * it = get_fed_it_tp(t);
	CAMLreturn(	Val_bool((*it).null()));
}

extern "C" CAMLprim value
stub_fed_iterator_has_next(value t)
{
	CAMLparam1(t);
	fed_t::iterator * it = get_fed_it_tp(t);
	CAMLreturn(	Val_bool((*it).hasNext()));
}

extern "C" CAMLprim value
stub_fed_iterator_remove(value t)
{
	CAMLparam1(t);
	fed_t::iterator * it = get_fed_it_tp(t);
	(*it).remove();
	CAMLreturn(Val_unit);
}
extern "C" CAMLprim value
stub_fed_iterator_insert(value t, value dbm)
{
	CAMLparam2(t, dbm);
	fed_t::iterator * it = get_fed_it_tp(t);
	it->insert(fdbm_t::create(*get_dbm_ptr(dbm), NULL));
	CAMLreturn(Val_unit);
}
