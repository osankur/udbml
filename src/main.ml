open Udbml;;
open Basic_types;;

let test_cmp () = 
  let a = Dbm.create 3 in
  Dbm.set_init a;
  Dbm.constrain a (1,0,(6,DBM_WEAK));
  Dbm.constrain a (0,1,(-2,DBM_WEAK));
  let b = Dbm.create 3 in
  Dbm.set_init b;
  Dbm.constrain b (1,0,(6,DBM_WEAK));
  Dbm.constrain b (0,1,(-2,DBM_WEAK));
  Printf.printf "eq: %b\n" (Dbm.equal a b);
  Printf.printf "%b %b\n" (Dbm.leq a b) (not (Dbm.lt a b));
  Printf.printf "%b %b\n" (Dbm.geq b a) (not (Dbm.gt b a));
  Dbm.constrain b (2,0,(3,DBM_WEAK));
  Printf.printf "%b \n" (Dbm.lt b a);
  Printf.printf "%b \n" (Dbm.gt a b);
  Printf.printf "noteq: %b\n" (Dbm.notequal a b);
  ()

let test_ops() = 
  let a = Dbm.create 3 in
  Dbm.set_init a;
  Dbm.constrain a (1,0,(6,DBM_WEAK));
  Dbm.constrain a (0,1,(-2,DBM_WEAK));
  Dbm.constrain a (2,0,(2,DBM_WEAK));
  Dbm.constrain a (0,2,(-1,DBM_WEAK));
  let b = Dbm.create 3 in
  Dbm.set_init b;
  Dbm.constrain b (1,0,(8,DBM_WEAK));
  Dbm.constrain b (0,1,(-8,DBM_WEAK));
  Dbm.constrain b (2,0,(1,DBM_WEAK));
  Dbm.constrain b (0,2,(-1,DBM_WEAK));
  Printf.printf "a:\n%s\n" (Dbm.to_string a);
  Printf.printf "b:\n%s\n" (Dbm.to_string b);
  let c = Dbm.copy b in
  Dbm.convex_add c a;
  Printf.printf "Convex add:\n%s\n" (Dbm.to_string c);
  Printf.printf "Does (a+b) contain %b?\n" (Dbm.leq a c );
  (*
  Printf.printf "Intersects: %b\n" (Dbm.intersects a b);
  Printf.printf "a:\n%s\n\n" (Dbm.to_string a);
  Dbm.intersect a b;
  Printf.printf "Intersection:\n%s\n\n" (Dbm.to_string a);
  Dbm.free_all_up a;
  Printf.printf "Free_clock :\n%s\n\n" (Dbm.to_string a);
     *)
  ()

let test_fed() = 
  let a = Dbm.create 3 in
  Dbm.set_init a;
	(*
  Dbm.constrain a (1,0,(6,DBM_WEAK));
  Dbm.constrain a (0,1,(-2,DBM_WEAK));
  Dbm.constrain a (2,0,(2,DBM_WEAK));
  Dbm.constrain a (0,2,(-1,DBM_WEAK));
	 *)
	(*
  let b = Dbm.create 3 in
  Dbm.set_init b;
  Dbm.constrain b (1,0,(8,DBM_WEAK));
  Dbm.constrain b (0,1,(-7,DBM_WEAK));
  Dbm.constrain b (2,0,(2,DBM_WEAK));
  Dbm.constrain b (0,2,(-1,DBM_WEAK));
	 *)
  let c = Fed.create 3 in
  Fed.add_dbm c a;
(*  Fed.add_dbm c b;
  Printf.printf "%b\n" (Fed.exact_leq (Fed.from_dbm a) c);
	Printf.printf "Iterating and printing all dbms of the fed\n";
 *)
	let it = Fed.begin_it c in
	Printf.printf "Just got the iterator\n";
	flush stdout;
	Gc.full_major();
	let dbm = Fed.Iterator.get it in
	Printf.printf "Just applied get\n";
	flush stdout;
	Gc.full_major();
	(*
	let count = ref 0 in
	while (not (Fed.Iterator.is_null it)) do
		let dbm = Fed.Iterator.get it in
		Printf.printf "Counting: %d\n" !count;
		incr(count);
		flush stdout;

									(*
		Printf.printf "%s\n\n" (Dbm.to_string dbm)
		Fed.Iterator.incr it
									 *)
	done;
	 *)
  ()
  
let test_resize() = 
  let dim = 4 in
  let a = Dbm.create dim in
  Dbm.set_init a;
  Dbm.constrain a (1,0,(6,DBM_WEAK));
  Dbm.constrain a (0,1,(-2,DBM_WEAK));
  Dbm.constrain a (2,0,(2,DBM_WEAK));
  Dbm.constrain a (0,2,(-1,DBM_WEAK));
  let subset = Bit_vector.create dim in
  Bit_vector.assign_bit subset 0 true;
  Bit_vector.assign_bit subset 1 true;
  Bit_vector.assign_bit subset 2 false;
  Bit_vector.assign_bit subset 3 true;
  let allclocks = Bit_vector.create dim in
  Bit_vector.assign_bit allclocks 0 true;
  Bit_vector.assign_bit allclocks 1 true;
  Bit_vector.assign_bit allclocks 2 true;
  Bit_vector.assign_bit allclocks 3 true;
  Printf.printf "subset:\t\t %s\nallclocks:\t%s\n" 
    (Bit_vector.to_string subset) (Bit_vector.to_string allclocks);
  let table = Dbm.resize a allclocks subset in
  for i = 0 to dim-1 do
    Printf.printf "%d: %d\n" i table.(i)
  done;
  ()
  
let test_mem() =
  let n = 200000 in
  let dim = 50 in
  for i = 0 to n do
    let a = Dbm.create dim in
    Dbm.set_init a
  done;
  ()
(*
  let a = Array.create n 0 in
  for i = 0 to n-1 do
    a.(i) <- a.(i) + 1
  done
  *)
let test_mem2() =
  let n = 1000000 in
  let dim = 50 in
  for i = 0 to n do
    let a = Array.create (dim*dim) 0 in
    a.(0) <- 1
  done;
  ()

let _ = 
  test_fed();
  ()
  
