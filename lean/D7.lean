/-
  Exact integer checks for "The Sphere Packing Problem in Dimension 7".
  The lattice L = {x ∈ E8 : x₁ + x₂ = 0} is a copy of E7; Λ = √2·L.
  Doubled coordinates make every vector integral.
-/

abbrev Vec := List Int

def dot (x y : Vec) : Int := (List.zipWith (· * ·) x y).foldl (· + ·) 0
def sub (x y : Vec) : Vec := List.zipWith (· - ·) x y
def smul (c : Int) (x : Vec) : Vec := x.map (c * ·)

/-- all sign vectors of length n with entries ±1 -/
def signs : Nat → List Vec
  | 0 => [[]]
  | n+1 => (signs n).flatMap fun s => [1 :: s, (-1) :: s]

def unit8 (i : Nat) (c : Int) : Vec := (List.range 8).map fun k => if k = i then c else 0

/-- the 240 roots of E8, doubled -/
def rootsE8x2 : List Vec :=
  let integral := (List.range 8).flatMap fun i => (List.range 8).flatMap fun j =>
    if i < j then [(1:Int),-1].flatMap fun a => [(1:Int),-1].map fun b =>
      List.zipWith (· + ·) (unit8 i (2*a)) (unit8 j (2*b)) else []
  let half := (signs 8).filter fun s => (s.filter (· < 0)).length % 2 = 0
  integral ++ half

/-- doubled roots of L: those with x₁ + x₂ = 0 -/
def rootsLx2 : List Vec := rootsE8x2.filter fun r => r[0]! + r[1]! = 0

/-- the 56 deep holes of L, doubled: ½(0,0,ε) with even number of +, and (±½,∓½,±e_i) -/
def deepLx2 : List Vec :=
  let a := ((signs 6).filter fun s => (s.filter (· > 0)).length % 2 = 0).map fun s => 0 :: 0 :: s
  let b := [(1:Int),-1].flatMap fun e => (List.range 6).flatMap fun i => [(1:Int),-1].map fun t =>
      e :: (-e) :: ((List.range 6).map fun k => if k = i then 2*t else 0)
  a ++ b

theorem card_roots : rootsLx2.length = 126 := by native_decide
theorem roots_norm : rootsLx2.all (fun r => dot r r = 8) = true := by native_decide
theorem card_deep : deepLx2.length = 56 := by native_decide
theorem deep_norm : deepLx2.all (fun v => dot v v = 6) = true := by native_decide

/-- each deep hole lies on the facets ⟨w,α⟩ ≤ 1 (doubled: ⟨2w,2α⟩ ≤ 4), with equality for exactly 27 roots -/
theorem deep_in_cell : deepLx2.all (fun v => rootsLx2.all fun r => dot v r ≤ 4) = true := by native_decide
theorem deep_27 : deepLx2.all (fun v => (rootsLx2.filter fun r => dot v r = 4).length = 27) = true := by
  native_decide

/-- extra centre y = (4/3) w.  With Y = 6·(2w) = 12 y... use 3·(2y) = 4·(2w): compare 9·|2y-2α|² ≥ 9·8.
    In the L-picture the packing condition |y - α|² ≥ 2 reads |4(2w) - 3(2α)|² ≥ 72. -/
theorem extra_far : deepLx2.all (fun v => rootsLx2.all fun r =>
    dot (sub (smul 4 v) (smul 3 r)) (sub (smul 4 v) (smul 3 r)) ≥ 72) = true := by native_decide
theorem extra_touch27 : deepLx2.all (fun v => (rootsLx2.filter fun r =>
    dot (sub (smul 4 v) (smul 3 r)) (sub (smul 4 v) (smul 3 r)) = 72).length = 27) = true := by native_decide
/-- |y| ≥ √2 in L, i.e. |4(2w)|² ≥ 9·8 -/
theorem extra_far_origin : deepLx2.all (fun v => dot (smul 4 v) (smul 4 v) ≥ 72) = true := by native_decide
/-- pairwise distance of extra centres ≥ (8/3)/√2 in L: |4(2w)-4(2w')|² ≥ 9·4·(32/9)=128 -/
theorem extra_pairwise : deepLx2.all (fun v => deepLx2.all fun w =>
    v == w || dot (sub (smul 4 v) (smul 4 w)) (sub (smul 4 v) (smul 4 w)) ≥ 128) = true := by native_decide

/-- the cap formula at the touching radius: 16 - 56/105 = 232/15 and 1/105 = 1/1680 of 16 -/
theorem arith : 16 * 105 * 15 - 56 * 15 = 232 * 105 ∧ 16 * 105 = 1680 := by decide

/-- Lemma 4.11: distinct deep holes have inner products 1/2, -1/2 or -3/2 (doubled: 2, -2, -6) -/
theorem deep_gram : deepLx2.all (fun v => deepLx2.all fun w =>
    v == w || dot v w = 2 || dot v w = -2 || dot v w = -6) = true := by native_decide

/-- the deep holes are pairwise distinct -/
theorem deep_nodup : deepLx2.Nodup := by native_decide

/-- the roots are pairwise distinct -/
theorem roots_nodup : rootsLx2.Nodup := by native_decide
