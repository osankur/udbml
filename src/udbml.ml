
module type BASIC_TYPES =
sig
 type bound_t
  type inequality_t = DBM_STRICT | DBM_WEAK
  type raw_t = bound_t * inequality_t
  type cindex_t
  type clock_constraint_t = cindex_t * cindex_t * raw_t
end

module Basic_types =
struct
  (** Abstract type for UDBM's bitvector used for resizing DBMs *)
  type bitvector
  (** Type for the constants that appear in DBMs' components *)
  type bound_t = int
  type inequality_t = DBM_STRICT | DBM_WEAK
  type raw_t = bound_t * inequality_t
  type cindex_t = int
  type clock_constraint_t = cindex_t * cindex_t * raw_t
  external get_infty : unit -> int = "stub_get_infty" [@@noalloc]
  let infty = get_infty()

  let string_of_raw = function
    | (b,_) when b = infty -> "<INF"
    | (b,DBM_WEAK) -> Printf.sprintf "<=%d" b
    | (b,DBM_STRICT) -> Printf.sprintf "<%d" b

  let add_raw (b,ineq) (b',ineq') =
    let ineq'' =
      if ineq = DBM_STRICT || ineq' = DBM_STRICT then DBM_STRICT
      else DBM_WEAK in
    if b = infty || b' = infty then(
      (infty, DBM_STRICT)
    ) else (b+b', ineq'')

  let invert_inequality = function
    | DBM_WEAK -> DBM_STRICT
    | DBM_STRICT -> DBM_WEAK

  let invert_raw (k,ineq) = (-k, invert_inequality ineq)

  let negate_constraint (x,y,(b,ineq)) =
    (y,x,invert_raw (b,ineq))

end

module Carray =
struct
  type t

  external _register_carray : unit -> unit = "caml_udbml_register_carray";;
  let _ = _register_carray ()

  external to_c : int array -> int -> t = "stub_carray_to_c";;
  external to_string : t -> string = "stub_carray_to_string";;
end


