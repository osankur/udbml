(**
   Notes
   
   1) All 32 bits integers are represented by int in this module.
   We should use Int32 in principle although we don't want to slow things down.
   Similarly all unsigned integers are int in this module. 
   This is not a problem for the cindex_t type and neither should it be for hash values.

   2) When loading from or copying to a raw_t array, dbm[i,j] = array[i*dim + j]
   
   3) TODO LIST
   - use "noalloc" when possible
   - ClockAccessor
   - Pretty printing
   - Print to a channel given in argument, not only to stdout
*)

module type BASIC_TYPES =
sig
  type bound_t
  type inequality_t = DBM_STRICT | DBM_WEAK
  type raw_t = bound_t * inequality_t
  type cindex_t
  type clock_constraint_t = cindex_t * cindex_t * raw_t
end

module Carray =
struct
  type t

  external to_c : int array -> int -> t = "stub_carray_to_c";;
end

module type IDBM =
sig
  type t
  
  include BASIC_TYPES

  val infty : int
  val max_dim_power : int
  val max_dim : int

  val create : int -> t
  val copy : t -> t
  val copy_to : t -> raw_t array
  val copy_from : raw_t array -> int -> t
  val is_empty : t -> bool
  val has_zero : t -> bool
  val hash : t -> int
  val dimension : t -> int
  val at : t -> int -> int -> raw_t
  val at_bound : t -> int -> int -> int

  val equal : t -> t -> bool
  val notequal : t -> t -> bool
  val lt : t -> t -> bool
  val gt : t -> t -> bool
  val leq : t -> t -> bool
  val geq : t -> t -> bool
  val closure_leq : Carray.t -> Carray.t -> t -> t -> bool

  val constrain : t -> clock_constraint_t -> unit

  val set_init : t -> unit
  val set_zero : t -> unit
  val is_init : t -> bool
  val is_zero : t -> bool

  val intersect : t -> t -> unit
  val intersects : t -> t -> bool
  val up : t -> unit
  val down : t -> unit

  val free_clock : t -> cindex_t -> unit
  val free_up : t -> cindex_t -> unit
  val free_down : t -> cindex_t -> unit
  val free_all_up : t -> unit
  val free_all_down : t -> unit

  val update_value : t -> cindex_t -> bound_t -> unit

  val is_unbounded : t -> bool

  val extrapolate_max_bounds : t -> Carray.t -> unit
  val diagonal_extrapolate_max_bounds : t -> Carray.t -> unit
  val extrapolate_lu_bounds : t -> Carray.t -> Carray.t -> unit
  val diagonal_extrapolate_lu_bounds : t -> Carray.t -> Carray.t -> unit

  val to_string : t -> string
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
  external get_infty : unit -> int = "stub_get_infty" "noalloc"
  let infty = get_infty()

  let string_of_raw = function
    | (b,_) when b = infty -> "(INF,<)"
    | (b,DBM_WEAK) -> Printf.sprintf "(%d,<=)" b
    | (b,DBM_STRICT) -> Printf.sprintf "(%d,<)" b
end

module Bit_vector = 
struct
  type t
  external create : int -> t = "stub_bitvector_create";;
  external copy : t -> t = "stub_bitvector_copy";;
  external assign_bit : t -> int -> bool -> unit = "stub_bitvector_assign_bit";;
  external get_bit : t  -> int -> bool = "stub_bitvector_get_bit";;
  external to_string : t -> string = "stub_bitvector_to_string";;
  external count : t -> int = "stub_bitvector_count";;
end

module ClockAccessor = 
struct
  (* TODO Fill in *)
end


module Dbm = 
struct
  type t
  include Basic_types
  external _max_dim_power : unit -> int = "stub_dbm_max_dim_power" "noalloc";;
  external _max_dim : unit -> int = "stub_dbm_max_dim" "noalloc";;
  
  let max_dim_power = _max_dim_power()
  let max_dim = _max_dim()

  external create : int -> t = "stub_dbm_create";;
  external dimension : t -> int = "stub_dbm_dimension" "noalloc";;
  external copy : t -> t = "stub_dbm_copy";;
  external is_empty : t -> bool = "stub_dbm_is_empty";;
  external has_zero : t -> bool = "stub_dbm_has_zero";;
  external hash : t -> int = "stub_dbm_hash";;
  external intern : t -> unit = "stub_dbm_intern";;

  external at : t -> int -> int -> raw_t = "stub_dbm_at";;
  external at_bound : t -> int -> int -> int = "stub_dbm_at_bound" "noalloc";;
  external equal : t -> t -> bool = "stub_dbm_equal";;
  external notequal : t -> t -> bool = "stub_dbm_notequal";;
  external lt : t -> t -> bool = "stub_dbm_lt";;
  external gt : t -> t -> bool = "stub_dbm_gt";;
  external leq : t -> t -> bool = "stub_dbm_leq";;
  external geq : t -> t -> bool = "stub_dbm_geq";;

  external closure_leq : Carray.t -> Carray.t -> t -> t -> bool =
    "stub_dbm_closure_leq" "noalloc";;

  external constrain : t -> clock_constraint_t -> unit = "stub_dbm_constrain";;
  external copy_to : t -> raw_t array = "stub_dbm_copy_to";;
  external set_init : t -> unit = "stub_dbm_set_init" "noalloc";;
  external set_zero : t -> unit = "stub_dbm_set_zero";;
  external is_zero : t -> bool = "stub_dbm_is_zero";;
  external is_init : t -> bool = "stub_dbm_is_init";;
  external intersect : t -> t -> unit = "stub_dbm_intersect" "noalloc";;
  external intersects : t -> t -> bool = "stub_dbm_intersects";;
  external up : t -> unit = "stub_dbm_up";;
  external down : t -> unit = "stub_dbm_down";;
  external free_clock : t -> cindex_t -> unit = "stub_dbm_free_clock";;

  external free_up : t -> cindex_t -> unit = "stub_dbm_free_up";;
  external free_down : t -> cindex_t -> unit = "stub_dbm_free_down";;
  external free_all_up : t -> unit = "stub_dbm_free_all_up";;
  external free_all_down : t -> unit = "stub_dbm_free_all_down";;

  external is_unbounded : t -> bool = "stub_dbm_is_unbounded";;
  
  let to_string t =
    if (is_empty t) then "false"
    else
      (
        let buf = Buffer.create 256 in
        let ar = copy_to t in
        let dim = dimension t in
        for i = 0 to (dim - 1) do
          for j = 0 to (dim - 1) do
            Buffer.add_string buf (Printf.sprintf "%s " (string_of_raw ar.(i*dim +j)))
          done;
          Buffer.add_string buf "\n"
        done;
        Buffer.contents buf
      )

  let print t = 
    let ar = copy_to t in
    let dim = dimension t in
    for i = 0 to (dim - 1) do
      for j = 0 to (dim - 1) do
        Printf.printf "%s " (string_of_raw ar.(i*dim +j))
      done;
      Printf.printf "\n"
    done
    
  external copy_from : raw_t array -> int -> t = "stub_dbm_copy_from";;
  external convex_add : t -> t -> unit = "stub_dbm_convex_add";;

  external update_value : t -> cindex_t -> bound_t -> unit = "stub_dbm_update_value";;
  external update_clock : t -> cindex_t -> cindex_t -> unit = "stub_dbm_update_clock";;
  external update_increment : t -> cindex_t -> bound_t -> unit = "stub_dbm_update_increment";;
  external update : t -> cindex_t -> cindex_t -> bound_t -> unit = "stub_dbm_update";;

  external satisfies : t -> clock_constraint_t -> bool = "stub_dbm_satisfies";;
  external is_subtraction_empty : t -> t -> bool = "stub_dbm_is_subtraction_empty";;

  
  external extrapolate_max_bounds : t -> Carray.t -> unit = "stub_dbm_extrapolate_max_bounds";;
  external diagonal_extrapolate_max_bounds : t -> Carray.t -> unit = "stub_dbm_diagonal_extrapolate_max_bounds";;
  external extrapolate_lu_bounds : t -> Carray.t -> Carray.t -> unit = "stub_dbm_extrapolate_lu_bounds" "noalloc";;
  external diagonal_extrapolate_lu_bounds : t -> Carray.t -> Carray.t -> unit = "stub_dbm_diagonal_extrapolate_lu_bounds";;
  
  external resize : t -> Bit_vector.t ->  Bit_vector.t -> int array = "stub_dbm_resize";;

  external _internal_addr : t -> int = "stub_dbm__internal_addr";;
end

module Fed =
struct
  type t
  include Basic_types
      
  external create : int -> t = "stub_fed_create";;
  external copy : t -> t = "stub_fed_copy";;
  external dimension : t -> int = "stub_fed_dimension";;
  external intern : t -> unit = "stub_fed_intern";;
  external is_empty : t -> bool = "stub_fed_is_empty";;
  external has_zero : t -> bool = "stub_fed_has_zero";;
  external hash : t -> int = "stub_fed_hash";;
  external set_init : t -> unit = "stub_fed_set_init";;

  external approx_equal : t -> t -> bool = "stub_fed_approx_equal";;
  external approx_notequal : t -> t -> bool = "stub_fed_approx_notequal";;
  external approx_lt : t -> t -> bool = "stub_fed_approx_lt";;
  external approx_leq : t -> t -> bool = "stub_fed_approx_leq";;
  external approx_gt : t -> t -> bool = "stub_fed_approx_gt";;
  external approx_geq : t -> t -> bool = "stub_fed_approx_geq";;

  external approx_equal_dbm : t -> Dbm.t -> bool = "stub_fed_approx_equal_dbm";;
  external approx_notequal_dbm : t -> Dbm.t -> bool = "stub_fed_approx_notequal_dbm";;
  external approx_lt_dbm : t -> Dbm.t -> bool = "stub_fed_approx_lt_dbm";;
  external approx_gt_dbm : t -> Dbm.t -> bool = "stub_fed_approx_gt_dbm";;
  external approx_geq_dbm : t -> Dbm.t -> bool = "stub_fed_approx_geq_dbm";;

  external exact_equal : t -> t -> bool = "stub_fed_exact_equal";;
  external exact_lt : t -> t -> bool = "stub_fed_exact_lt";;
  external exact_leq : t -> t -> bool = "stub_fed_exact_leq";;
  external exact_gt : t -> t -> bool = "stub_fed_exact_gt";;
  external exact_geq : t -> t -> bool = "stub_fed_exact_geq";;

  external exact_equal_dbm : t -> Dbm.t -> bool = "stub_fed_exact_equal_dbm";;
  external exact_leq_dbm : t -> Dbm.t -> bool = "stub_fed_exact_leq_dbm";;
  external exact_lt_dbm : t -> Dbm.t -> bool = "stub_fed_exact_leq_dbm";;
  external exact_gt_dbm : t -> Dbm.t -> bool = "stub_fed_exact_gt_dbm";;
  external exact_geq_dbm : t -> Dbm.t -> bool = "stub_fed_exact_geq_dbm";;

  external convex_hull : t -> t = "stub_fed_convex_hull";;
  external union : t -> t -> unit = "stub_fed_union";;
  external add : t -> t -> unit = "stub_fed_add";;
  external add_dbm : t -> Dbm.t -> unit = "stub_fed_add_dbm";;
  external append : t -> t -> unit = "stub_fed_append";;
  external append_end : t -> t -> unit = "stub_fed_append_end";;
  external append_begin : t -> t -> unit = "stub_fed_append_begin";;
  external steal : t -> t -> unit = "stub_fed_steal";;

  external convex_union : t -> t -> unit = "stub_fed_convex_union";;
  external convex_union_dbm : t -> Dbm.t -> unit = "stub_fed_convex_union_dbm";;
  external intersect : t -> t -> unit = "stub_fed_intersect";;
  external intersect_dbm : t -> Dbm.t -> unit = "stub_fed_intersect_dbm";;
  external subtract : t -> t -> unit = "stub_fed_subtract";;
  external subtract_dbm : t -> Dbm.t -> unit = "stub_fed_subtract_dbm";;
  external up : t -> unit = "stub_fed_up";;
  external down : t -> unit = "stub_fed_down";;
  external free_clock : t -> cindex_t -> unit = "stub_fed_free_clock";;
  external constrain : t -> clock_constraint_t -> unit = "stub_fed_constrain";;
  external reduce : t -> unit = "stub_fed_reduce";;
  external expensive_reduce : t -> unit = "stub_fed_expensive_reduce";;
  external merge_reduce : t -> int -> int -> unit = "stub_fed_merge_reduce";;
  external convex_reduce : t -> unit = "stub_fed_convex_reduce";;
  external expensive_convex_reduce : t -> unit = "stub_fed_convex_reduce";;
  external partition_reduce : t -> unit = "stub_fed_convex_reduce";;

  external predt : t -> t -> Dbm.t option -> unit = "stub_fed_predt";;
  external is_included_in_predt : t -> t -> t -> bool = "stub_fed_is_included_in_predt";;

  external _internal_addr : t -> int = "stub_fed__internal_addr";;
  
  let from_dbm a = 
    let b = create (Dbm.dimension a) in
    add_dbm b a;
    b
    

  module Iterator =
  struct 
    type t
    external get : t -> Dbm.t = "stub_fed_iterator_get";;
    external incr : t -> unit = "stub_fed_iterator_incr";;
    external is_null : t -> bool = "stub_fed_iterator_is_null";;
    external has_next : t -> bool = "stub_fed_iterator_has_next";;
    external remove : t -> unit = "stub_fed_iterator_remove";;
    external insert : t -> Dbm.t -> unit = "stub_fed_iterator_insert";;
  end 
  
	external begin_it : t -> Iterator.t  = "stub_fed_begin_it";;
  let iter t f = 
    let it = begin_it t in
    while not(Iterator.is_null it) do
      let dbm = Iterator.get it in
      f dbm;
      Iterator.incr it;
    done

end


