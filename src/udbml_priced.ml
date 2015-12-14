open Udbml_unpriced

module PDbm =
struct
  type t

  include Udbml.Basic_types

  let max_dim_power = Dbm.max_dim_power
  let max_dim = Dbm.max_dim

  external create : int -> t = "stub_pdbm_create";;
  external dimension : t -> int = "stub_pdbm_dimension" "noalloc";;
  external copy : t -> t = "stub_pdbm_copy";;
  external is_empty : t -> bool = "stub_pdbm_is_empty" "noalloc";;
  external hash : t -> int = "stub_pdbm_hash" "noalloc";;
  external at : t -> int -> int -> raw_t = "stub_pdbm_at";;
  external at_bound : t -> int -> int -> int = "stub_pdbm_at_bound" "noalloc";;

  external equal : t -> t -> bool = "stub_pdbm_equal" "noalloc";;
  external notequal : t -> t -> bool = "stub_pdbm_not_equal" "noalloc";;
  external lt : t -> t -> bool = "stub_pdbm_lt" "noalloc";;
  external gt : t -> t -> bool = "stub_pdbm_gt" "noalloc";;
  external leq : t -> t -> bool = "stub_pdbm_leq" "noalloc";;
  external geq : t -> t -> bool = "stub_pdbm_geq" "noalloc";;

  external square_inclusion_exp : t -> t -> Udbml.Carray.t -> bool = "stub_pdbm_square_inclusion_exp" "noalloc";;

  external constrain : t -> clock_constraint_t -> unit =
    "stub_pdbm_constrain" "noalloc";;

  external set_init : t -> unit = "stub_pdbm_set_init" "noalloc";;
  external set_zero : t -> unit = "stub_pdbm_set_zero" "noalloc";;

  external intersect : t -> Udbml_unpriced.Dbm.t -> unit = "stub_pdbm_intersect" "noalloc";;

  external is_unbounded : t -> bool = "stub_pdbm_is_unbounded" "noalloc";;

  external infimum : t -> int = "stub_pdbm_infimum" "noalloc";;

  external to_string : t -> string = "stub_pdbm_to_string";;
end

module PFed =
struct
  type t

  include Udbml.Basic_types

  external create : int -> t = "stub_pfed_create";;
  external hash : t -> int = "stub_pfed_hash" "noalloc";;

  external is_empty : t -> bool = "stub_pfed_is_empty" "noalloc";;
  external set_empty : t -> unit = "stub_pfed_set_empty";;
  external add_dbm : t -> PDbm.t -> unit = "stub_pfed_add_dbm" "noalloc";;

  external up : t -> unit = "stub_pfed_up" "noalloc";;
  external update_value : t -> cindex_t -> bound_t -> unit = "stub_pfed_update_value" "noalloc";;
  external intersect_dbm : t -> Udbml_unpriced.Dbm.t -> unit = "stub_pfed_intersect_dbm" "noalloc";;

  (* TODO should be hidden *)
  module Iterator =
  struct
    type t
    
    external get : t -> PDbm.t = "stub_pfed_iterator_get";;
    external incr : t -> unit = "stub_pfed_iterator_incr" "noalloc";;
    external not_equal : t -> t -> bool = "stub_pfed_iterator_notequal" "noalloc";;
(*
    external remove : t -> unit = "stub_pfed_iterator_remove";;
    external insert : t -> PDbm.t -> unit = "stub_pfed_iterator_insert";;
*)
  end

  type iterator_t = t * Iterator.t

  (* TODO should be hidden *)
  external _begin_it : t -> Iterator.t = "stub_pfed_iterator_begin";;
  external _end_it : t -> Iterator.t = "stub_pfed_iterator_end";;

  let begin_it p = (p, _begin_it p)
  let end_it p = (p, _end_it p)

  let get_it (_,t) = Iterator.get t

  let incr_it (_,t) = Iterator.incr t

  let neq_it (a,t1) (a,t2) = Iterator.not_equal t1 t2

  let iter t f =
    let it = begin_it t in
    let ite = end_it t in
    while (neq_it it ite) do
      f (get_it it);
      incr_it it
    done

  let from_dbm a =
    let f = create (PDbm.dimension a) in
    add_dbm f a;
    f
end

