(* Content-type: application/vnd.wolfram.cdf.text *)

(*** Wolfram CDF File ***)
(* http://www.wolfram.com/cdf *)

(* CreatedBy='Wolfram 14.3' *)

(***************************************************************************)
(*                                                                         *)
(*                                                                         *)
(*  Under the Wolfram FreeCDF terms of use, this file and its content are  *)
(*  bound by the Creative Commons BY-SA Attribution-ShareAlike license.    *)
(*                                                                         *)
(*        For additional information concerning CDF licensing, see:        *)
(*                                                                         *)
(*         www.wolfram.com/cdf/adopting-cdf/licensing-options.html         *)
(*                                                                         *)
(*                                                                         *)
(***************************************************************************)

(*CacheID: 234*)
(* Internal cache information:
NotebookFileLineBreakTest
NotebookFileLineBreakTest
NotebookDataPosition[      1084,         20]
NotebookDataLength[    175927,       2783]
NotebookOptionsPosition[    175274,       2757]
NotebookOutlinePosition[    175715,       2774]
CellTagsIndexPosition[    175672,       2771]
WindowTitle->Platonic Solids -- Interactive Viewer
WindowFrame->Normal*)

(* Beginning of Notebook Content *)
Notebook[{

Cell[CellGroupData[{
Cell["Platonic Self-Dual Solitons \[LongDash] Interactive Viewer", "Title",ExpressionUUID->"555a3f2a-a95e-411c-a6e2-0b9b2c90c275"],

Cell["\<\
All 5 Platonic solids, U and V, live range control. Giuseppe Braz da Silva \
Marcelino \[LongDash] 2026\
\>", "Subtitle",ExpressionUUID->"434989e0-36bc-463a-86a0-6b89300c1aa1"],

Cell[CellGroupData[{

Cell["1. Setup (evaluate this once)", "Section",ExpressionUUID->"2e426786-0546-4dd8-9700-c5976fa3b105"],

Cell[CellGroupData[{

Cell["\<\
(* ============================================================
   Platonic self-dual pre-potentials -- setup
   Vertex sets, rotation, U and V builders for all 5 solids
   ============================================================ *)

rotAxis = Normalize[{1, 1, 0}];
rotAngle = Pi/6;

RotMatrix[axis_, theta_] := Module[{n = axis, K},
  K = {{0, -n[[3]], n[[2]]}, {n[[3]], 0, -n[[1]]}, {-n[[2]], n[[1]], 0}};
  IdentityMatrix[3] Cos[theta] + (1 - Cos[theta]) Outer[Times, n, n] + \
Sin[theta] K
];
Rmat = RotMatrix[rotAxis, rotAngle];
Rot[v_] := Rmat . v;

tau = N[(1 + Sqrt[5])/2];

(* ---------- Vertex sets (unit sphere) ---------- *)
VT = Normalize /@ N[{{1, 1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, 1}}];

VOraw = N[{{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, \
-1}}];
VO = Rot /@ VOraw;

VCraw = Normalize /@ N[Tuples[{1, -1}, 3]];
VC = Rot /@ VCraw;

VIraw = Normalize /@ DeleteDuplicates[N[Join[
    Tuples[{{0}, {1, -1}, {tau, -tau}}],
    Tuples[{{1, -1}, {tau, -tau}, {0}}],
    Tuples[{{tau, -tau}, {0}, {1, -1}}]
    ]]];
VI = Rot /@ VIraw;

VDraw = N[Join[
   Tuples[{1, -1}, 3]/Sqrt[3],
   Tuples[{{0}, {tau, -tau}, {1/tau, -1/tau}}]/Sqrt[3],
   Tuples[{{1/tau, -1/tau}, {0}, {tau, -tau}}]/Sqrt[3],
   Tuples[{{tau, -tau}, {1/tau, -1/tau}, {0}}]/Sqrt[3]
   ]];
VD = Rot /@ VDraw;

(* ---------- Pre-potential U, symbolic in phi1,phi2 ---------- *)
UExpr[v_, phi1_, phi2_] := Module[{r2, p, Cc},
  r2 = phi1^2 + phi2^2;
  p = {2 phi1, 2 phi2, r2 - 1}/(1 + r2);
  Cc = Times @@ (1/(1 - #[[3]]) & /@ v);
  Cc*Times @@ ((1 - #.p) & /@ v)
];

(* Build compiled U and V. *)
BuildUV[v_] := Module[{phi1, phi2, Uexpr, dU, etaScalar, etaInvMat, Vexpr, \
Ucomp, Vcomp},
  Uexpr = UExpr[v, phi1, phi2];
  dU = {D[Uexpr, phi1], D[Uexpr, phi2]};
  etaScalar = (1 + phi1^2 + phi2^2)^2/4;
  etaInvMat = etaScalar*IdentityMatrix[2]; (* eta^{-1}_{ab}, diagonal for \
this conformal metric *)
  Vexpr = (1/2) dU . etaInvMat . dU; (* full a,b sum: (1/2) \
Sum[etaInvMat[[a,b]] dU[[a]] dU[[b]], {a,2}, {b,2}] *)
  Ucomp = Compile[{{phi1, _Real}, {phi2, _Real}}, Evaluate[Uexpr]];
  Vcomp = Compile[{{phi1, _Real}, {phi2, _Real}}, Evaluate[Vexpr]];
  {Ucomp, Vcomp}
];

{UTc, VTc} = BuildUV[VT];
{UOc, VOc} = BuildUV[VO];
{UCc, VCc} = BuildUV[VC];
{UIc, VIc} = BuildUV[VI];
{UDc, VDc} = BuildUV[VD];

(*  *)
solidData = <|
   \"Tetrahedron\" -> <|\"U\" -> UTc, \"V\" -> VTc, \"Umax\" -> 4/3|>,
   \"Octahedron\" -> <|\"U\" -> UOc, \"V\" -> VOc, \"Umax\" -> 1.548|>,
   \"Cube\" -> <|\"U\" -> UCc, \"V\" -> VCc, \"Umax\" -> 2.023|>,
   \"Icosahedron\" -> <|\"U\" -> UIc, \"V\" -> VIc, \"Umax\" -> 1.376|>,
   \"Dodecahedron\" -> <|\"U\" -> UDc, \"V\" -> VDc, \"Umax\" -> 1.938|>
|>;

(* V has no simple closed form, so its normalization is found numerically *)
Print[\"Computing V normalization for each solid (a few seconds)...\"];
vMaxFixed = Association[
   # -> Max[Flatten[Table[solidData[#][\"V\"][x, y], {x, -5., 5., 0.05}, {y, \
-5., 5., 0.05}]]] & /@ Keys[solidData]
];
Print[\"Done. V normalization constants: \", vMaxFixed];\
\>", "Input",
 CellChangeTimes->{{3.994254938598755*^9, 
  3.99425497453087*^9}},ExpressionUUID->"f89fdac5-27df-4551-8df8-\
980dd68db759"],

Cell[CellGroupData[{

Cell[BoxData["\<\"Computing V normalization for each solid (a few seconds)...\
\"\>"], "Print",
 CellChangeTimes->{3.994254881528*^9},
 CellLabel->
  "During evaluation of \
In[1]:=",ExpressionUUID->"39f1c8ed-aed1-42a7-9610-98bd6745a5cc"],

Cell[BoxData[
 InterpretationBox[
  RowBox[{"\<\"Done. V normalization constants: \"\>", "\[InvisibleSpace]", 
   RowBox[{"\[LeftAssociation]", 
    RowBox[{
     RowBox[{"\<\"Tetrahedron\"\>", "\[Rule]", "1.5690432921983544`"}], ",", 
     RowBox[{"\<\"Octahedron\"\>", "\[Rule]", "3.586920069617875`"}], ",", 
     RowBox[{"\<\"Cube\"\>", "\[Rule]", "5.816698900306769`"}], ",", 
     RowBox[{"\<\"Icosahedron\"\>", "\[Rule]", "6.156279288428893`"}], ",", 
     RowBox[{"\<\"Dodecahedron\"\>", "\[Rule]", "10.625843871752878`"}]}], 
    "\[RightAssociation]"}]}],
  SequenceForm[
  "Done. V normalization constants: ", <|
   "Tetrahedron" -> 1.5690432921983544`, "Octahedron" -> 3.586920069617875, 
    "Cube" -> 5.816698900306769, "Icosahedron" -> 6.156279288428893, 
    "Dodecahedron" -> 10.625843871752878`|>],
  Editable->False]], "Print",
 CellChangeTimes->{3.994254881935417*^9},
 CellLabel->
  "During evaluation of \
In[1]:=",ExpressionUUID->"deaa8895-7a87-4ab6-9b84-849b8228df8b"]
}, Open  ]]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell["2. Interactive viewer", "Section",ExpressionUUID->"562eeede-2301-442c-8f88-b77a8055cde6"],

Cell[CellGroupData[{

Cell["\<\
Manipulate[
 Module[{data = solidData[solidName], fn, norm, label, pts},
  If[quantity == \"U\",
   fn = data[\"U\"]; norm = data[\"Umax\"];
   label = \"U(\\!\\(\\*SubscriptBox[\\(\[Phi]\\), \
\\(1\\)]\\),\\!\\(\\*SubscriptBox[\\(\[Phi]\\), \\(2\\)]\\))  \[LongDash]  \" \
<> solidName,
   fn = data[\"V\"]; norm = vMaxFixed[solidName];
   label = \"V(\\!\\(\\*SubscriptBox[\\(\[Phi]\\), \
\\(1\\)]\\),\\!\\(\\*SubscriptBox[\\(\[Phi]\\), \\(2\\)]\\))  \[LongDash]  \" \
<> solidName
   ];
  pts = Round[35 + 35 (15 - Min[r, 15])/14];
  Plot3D[fn[x, y]/norm, {x, -r, r}, {y, -r, r},
   MeshFunctions -> {#3 &}, Mesh -> 50,
   ColorFunction -> (ColorData[\"TemperatureMap\"][#3] &),
   ColorFunctionScaling -> False,
   PlotRange -> {0, 1.05}, PlotPoints -> pts,
   AxesLabel -> {\"\\!\\(\\*SubscriptBox[\\(\[Phi]\\), \\(1\\)]\\)\", \
\"\\!\\(\\*SubscriptBox[\\(\[Phi]\\), \\(2\\)]\\)\", \"\"},
   PlotLabel -> label, ImageSize -> 520,
   PerformanceGoal -> \"Quality\"]
  ],
 {{solidName, \"Dodecahedron\", \"Solid\"}, {\"Tetrahedron\", \"Octahedron\", \
\"Cube\", \"Icosahedron\", \"Dodecahedron\"}},
 {{quantity, \"U\", \"Quantity\"}, {\"U\", \"V\"}},
 {{r, 3, \"range \[PlusMinus]\"}, 1, 15, 0.5, Appearance -> \"Labeled\"},
 ControlPlacement -> Top,
 SaveDefinitions -> True,
 TrackedSymbols :> {solidName, quantity, r}
]\
\>", "Input",
 CellChangeTimes->{{3.9942551008684397`*^9, 3.994255120453586*^9}, {
  3.994255256495343*^9, 3.994255281837858*^9}},
 CellLabel->"In[49]:=",ExpressionUUID->"b27cfa94-ad6b-4664-a3c4-0f7566674086"],

Cell[BoxData[
 TagBox[
  StyleBox[
   DynamicModuleBox[{$CellContext`quantity$$ = "U", $CellContext`r$$ = 
    3, $CellContext`solidName$$ = "Cube", Typeset`show$$ = True, 
    Typeset`bookmarkList$$ = {}, Typeset`bookmarkMode$$ = "Menu", 
    Typeset`animator$$, Typeset`animvar$$ = 1, Typeset`name$$ = 
    "\"untitled\"", Typeset`specs$$ = {{{
       Hold[$CellContext`solidName$$], "Dodecahedron", "Solid"}, {
      "Tetrahedron", "Octahedron", "Cube", "Icosahedron", "Dodecahedron"}, 
      ControlType -> SetterBar}, {{
       Hold[$CellContext`quantity$$], "U", "Quantity"}, {"U", "V"}, 
      ControlType -> SetterBar}, {{
       Hold[$CellContext`r$$], 3, "range \[PlusMinus]"}, 1, 15, 0.5, 
      ControlType -> Manipulator}}, Typeset`size$$ = {
    520., {207.20478000624374`, 212.82978000624374`}}, Typeset`update$$ = 0, 
    Typeset`initDone$$, Typeset`skipInitDone$$ = False, 
    Typeset`keyframeActionsQ$$ = False, Typeset`keyframeList$$ = {}}, 
    DynamicBox[Manipulate`ManipulateBoxes[
     1, StandardForm, 
      "Variables" :> {$CellContext`quantity$$ = "U", $CellContext`r$$ = 
        3, $CellContext`solidName$$ = "Dodecahedron"}, 
      "ControllerVariables" :> {}, 
      "OtherVariables" :> {
       Typeset`show$$, Typeset`bookmarkList$$, Typeset`bookmarkMode$$, 
        Typeset`animator$$, Typeset`animvar$$, Typeset`name$$, 
        Typeset`specs$$, Typeset`size$$, Typeset`update$$, Typeset`initDone$$,
         Typeset`skipInitDone$$, Typeset`keyframeActionsQ$$, 
        Typeset`keyframeList$$}, "Body" :> 
      Module[{$CellContext`data$ = \
$CellContext`solidData[$CellContext`solidName$$], $CellContext`fn$, \
$CellContext`norm$, $CellContext`label$, $CellContext`pts$}, 
        If[$CellContext`quantity$$ == 
          "U", $CellContext`fn$ = $CellContext`data$[
            "U"]; $CellContext`norm$ = $CellContext`data$[
            "Umax"]; $CellContext`label$ = 
           "U(\!\(\*SubscriptBox[\(\[Phi]\), \(1\)]\),\!\(\*SubscriptBox[\(\
\[Phi]\), \(2\)]\))  \[LongDash]  " <> $CellContext`solidName$$, \
$CellContext`fn$ = $CellContext`data$[
            "V"]; $CellContext`norm$ = \
$CellContext`vMaxFixed[$CellContext`solidName$$]; $CellContext`label$ = 
           "V(\!\(\*SubscriptBox[\(\[Phi]\), \(1\)]\),\!\(\*SubscriptBox[\(\
\[Phi]\), \(2\)]\))  \[LongDash]  " <> $CellContext`solidName$$]; \
$CellContext`pts$ = Round[35 + 35 ((15 - Min[$CellContext`r$$, 15])/14)]; 
        Plot3D[$CellContext`fn$[$CellContext`x, \
$CellContext`y]/$CellContext`norm$, {$CellContext`x, -$CellContext`r$$, \
$CellContext`r$$}, {$CellContext`y, -$CellContext`r$$, $CellContext`r$$}, 
          MeshFunctions -> {#3& }, Mesh -> 50, 
          ColorFunction -> (ColorData["TemperatureMap"][#3]& ), 
          ColorFunctionScaling -> False, PlotRange -> {0, 1.05}, 
          PlotPoints -> $CellContext`pts$, 
          AxesLabel -> {
           "\!\(\*SubscriptBox[\(\[Phi]\), \(1\)]\)", 
            "\!\(\*SubscriptBox[\(\[Phi]\), \(2\)]\)", ""}, 
          PlotLabel -> $CellContext`label$, ImageSize -> 520, PerformanceGoal -> 
          "Quality"]], 
      "Specifications" :> {{{$CellContext`solidName$$, "Dodecahedron", 
          "Solid"}, {
         "Tetrahedron", "Octahedron", "Cube", "Icosahedron", 
          "Dodecahedron"}}, {{$CellContext`quantity$$, "U", "Quantity"}, {
         "U", "V"}}, {{$CellContext`r$$, 3, "range \[PlusMinus]"}, 1, 15, 0.5,
          Appearance -> "Labeled"}}, 
      "Options" :> {
       ControlPlacement -> Top, 
        TrackedSymbols :> {$CellContext`solidName$$, $CellContext`quantity$$, \
$CellContext`r$$}}, "DefaultOptions" :> {}],
     ImageSizeCache->{
      565.0999999999999, {276.15868625624375`, 281.78368625624375`}},
     SingleEvaluation->True],
    Deinitialization:>None,
    DynamicModuleValues:>{},
    Initialization:>({$CellContext`solidData = <|
        "Tetrahedron" -> <|"U" -> CompiledFunction[{11, 14.3, 5470}, {
              Blank[Real], 
              Blank[Real]}, {{3, 0, 0}, {3, 0, 1}, {3, 0, 
             22}}, {{0.5773502691896258, {3, 0, 15}}, {-1, {2, 0, 1}}, {
              2.2500000000000013`, {3, 0, 17}}, {-0.5773502691896258, {3, 0, 
               7}}, {1.1547005383792517`, {3, 0, 9}}, {
              1, {2, 0, 0}}, {-1.1547005383792517`, {3, 0, 11}}}, {0, 2, 23, 
             0, 0}, {{40, 56, 3, 0, 0, 3, 0, 2}, {40, 56, 3, 0, 1, 3, 0, 3}, {
              10, 0, 4}, {13, 4, 2, 3, 4}, {40, 60, 3, 0, 4, 3, 0, 5}, {10, 1,
               6}, {13, 6, 2, 3, 6}, {16, 7, 6, 5, 8}, {16, 9, 0, 5, 10}, {16,
               11, 1, 5, 12}, {16, 11, 0, 5, 13}, {16, 9, 1, 5, 14}, {16, 15, 
              6, 5, 16}, {10, 0, 18}, {13, 18, 13, 12, 8, 18}, {10, 0, 19}, {
              13, 19, 10, 14, 8, 19}, {10, 0, 20}, {13, 20, 10, 12, 16, 20}, {
              10, 0, 21}, {13, 21, 13, 14, 16, 21}, {16, 17, 18, 19, 20, 21, 
              22}, {1}}, 
             Function[{$CellContext`phi1$22384, $CellContext`phi2$22384}, 
              
              Block[{Compile`$2, Compile`$3, Compile`$5, Compile`$6, 
                Compile`$9, Compile`$10, Compile`$12, Compile`$8, Compile`$7, 
                Compile`$13, Compile`$15}, 
               Compile`$2 = $CellContext`phi1$22384^2; 
               Compile`$3 = $CellContext`phi2$22384^2; 
               Compile`$5 = 1 + Compile`$2 + Compile`$3; 
               Compile`$6 = Compile`$5^(-1); 
               Compile`$9 = -1 + Compile`$2 + Compile`$3; 
               Compile`$10 = (-0.5773502691896258) Compile`$9 Compile`$6; 
               Compile`$12 = 
                1.1547005383792517` $CellContext`phi1$22384 Compile`$6; 
               Compile`$8 = (-1.1547005383792517`) $CellContext`phi2$22384 
                 Compile`$6; 
               Compile`$7 = (-1.1547005383792517`) $CellContext`phi1$22384 
                 Compile`$6; 
               Compile`$13 = 
                1.1547005383792517` $CellContext`phi2$22384 Compile`$6; 
               Compile`$15 = 0.5773502691896258 Compile`$9 Compile`$6; 
               2.2500000000000013` (1 + Compile`$7 + Compile`$8 + 
                 Compile`$10) (1 + Compile`$12 + Compile`$13 + Compile`$10) (
                 1 + Compile`$12 + Compile`$8 + Compile`$15) (1 + Compile`$7 + 
                 Compile`$13 + Compile`$15)]], Evaluate], "V" -> 
           CompiledFunction[{11, 14.3, 5470}, {
              Blank[Real], 
              Blank[Real]}, {{3, 0, 0}, {3, 0, 1}, {3, 0, 
             43}}, {{0.5773502691896258, {3, 0, 19}}, {
              0.5, {3, 0, 37}}, {-1, {2, 0, 1}}, {
              2.3094010767585034`, {3, 0, 27}}, {
              2.2500000000000013`, {3, 0, 39}}, {-0.5773502691896258, {3, 0, 
               10}}, {1.1547005383792517`, {3, 0, 8}}, {
              1, {2, 0, 0}}, {-2.3094010767585034`, {3, 0, 
               21}}, {-1.1547005383792517`, {3, 0, 12}}, {
              0.25, {3, 0, 38}}}, {0, 2, 46, 0, 
             0}, {{40, 56, 3, 0, 0, 3, 0, 2}, {40, 56, 3, 0, 1, 3, 0, 3}, {10,
               0, 4}, {13, 4, 2, 3, 4}, {40, 56, 3, 0, 4, 3, 0, 5}, {40, 60, 
              3, 0, 5, 3, 0, 6}, {40, 60, 3, 0, 4, 3, 0, 5}, {10, 1, 7}, {13, 
              7, 2, 3, 7}, {16, 8, 0, 5, 9}, {16, 10, 7, 5, 11}, {16, 12, 1, 
              5, 13}, {16, 12, 0, 7, 6, 14}, {16, 12, 0, 5, 15}, {10, 0, 
              16}, {13, 16, 15, 13, 11, 16}, {16, 8, 1, 5, 17}, {10, 0, 18}, {
              13, 18, 9, 17, 11, 18}, {16, 19, 7, 5, 20}, {16, 21, 2, 6, 
              22}, {16, 21, 0, 1, 6, 23}, {16, 8, 5, 24}, {10, 0, 25}, {13, 
              25, 9, 13, 20, 25}, {10, 0, 26}, {13, 26, 15, 17, 20, 26}, {16, 
              27, 2, 6, 28}, {16, 27, 0, 1, 6, 29}, {16, 8, 0, 7, 6, 30}, {16,
               12, 5, 31}, {40, 56, 3, 0, 4, 3, 0, 32}, {16, 12, 1, 7, 6, 
              33}, {16, 21, 3, 6, 34}, {16, 27, 3, 6, 35}, {16, 8, 1, 7, 6, 
              36}, {13, 28, 23, 14, 31, 9, 40}, {16, 39, 40, 16, 18, 25, 
              41}, {13, 22, 29, 14, 24, 9, 40}, {16, 39, 40, 16, 18, 26, 
              42}, {13, 22, 23, 30, 24, 15, 40}, {16, 39, 40, 16, 25, 26, 
              43}, {13, 28, 29, 30, 31, 15, 40}, {16, 39, 40, 18, 25, 26, 
              44}, {13, 41, 42, 43, 44, 41}, {40, 56, 3, 0, 41, 3, 0, 42}, {
              16, 38, 32, 42, 41}, {13, 29, 34, 33, 24, 17, 42}, {16, 39, 42, 
              16, 18, 25, 43}, {13, 23, 35, 33, 31, 17, 42}, {16, 39, 42, 16, 
              18, 26, 44}, {13, 23, 34, 36, 24, 13, 42}, {16, 39, 42, 16, 25, 
              26, 40}, {13, 29, 35, 36, 31, 13, 42}, {16, 39, 42, 18, 25, 26, 
              45}, {13, 43, 44, 40, 45, 43}, {40, 56, 3, 0, 43, 3, 0, 44}, {
              16, 38, 32, 44, 43}, {13, 41, 43, 41}, {16, 37, 41, 43}, {1}}, 
             Function[{$CellContext`phi1$22384, $CellContext`phi2$22384}, 
              
              Block[{Compile`$11, Compile`$14, Compile`$16, Compile`$18, 
                Compile`$23, Compile`$21, Compile`$25, Compile`$29, 
                Compile`$28, Compile`$22, Compile`$27, Compile`$30, 
                Compile`$31, Compile`$32, Compile`$33, Compile`$36, 
                Compile`$20, Compile`$38, Compile`$34, Compile`$40, 
                Compile`$19, Compile`$37, Compile`$42, Compile`$24, 
                Compile`$17, Compile`$51, Compile`$50, Compile`$54, 
                Compile`$57}, Compile`$11 = $CellContext`phi1$22384^2; 
               Compile`$14 = $CellContext`phi2$22384^2; 
               Compile`$16 = 1 + Compile`$11 + Compile`$14; 
               Compile`$18 = Compile`$16^(-2); Compile`$23 = Compile`$16^(-1); 
               Compile`$21 = -1 + Compile`$11 + Compile`$14; 
               Compile`$25 = 
                1.1547005383792517` $CellContext`phi1$22384 Compile`$23; 
               Compile`$29 = (-0.5773502691896258) Compile`$21 Compile`$23; 
               Compile`$28 = (-1.1547005383792517`) $CellContext`phi2$22384 
                 Compile`$23; 
               Compile`$22 = (-1.1547005383792517`) $CellContext`phi1$22384 
                 Compile`$21 Compile`$18; 
               Compile`$27 = (-1.1547005383792517`) $CellContext`phi1$22384 
                 Compile`$23; 
               Compile`$30 = 1 + Compile`$27 + Compile`$28 + Compile`$29; 
               Compile`$31 = 
                1.1547005383792517` $CellContext`phi2$22384 Compile`$23; 
               Compile`$32 = 1 + Compile`$25 + Compile`$31 + Compile`$29; 
               Compile`$33 = 0.5773502691896258 Compile`$21 Compile`$23; 
               Compile`$36 = (-2.3094010767585034`) Compile`$11 Compile`$18; 
               Compile`$20 = (-2.3094010767585034`) $CellContext`phi1$22384 \
$CellContext`phi2$22384 Compile`$18; 
               Compile`$38 = 1.1547005383792517` Compile`$23; 
               Compile`$34 = 1 + Compile`$25 + Compile`$28 + Compile`$33; 
               Compile`$40 = 1 + Compile`$27 + Compile`$31 + Compile`$33; 
               Compile`$19 = 2.3094010767585034` Compile`$11 Compile`$18; 
               Compile`$37 = 
                2.3094010767585034` $CellContext`phi1$22384 \
$CellContext`phi2$22384 Compile`$18; 
               Compile`$42 = 
                1.1547005383792517` $CellContext`phi1$22384 Compile`$21 
                 Compile`$18; 
               Compile`$24 = (-1.1547005383792517`) Compile`$23; 
               Compile`$17 = Compile`$16^2; 
               Compile`$51 = (-1.1547005383792517`) $CellContext`phi2$22384 
                 Compile`$21 Compile`$18; 
               Compile`$50 = (-2.3094010767585034`) Compile`$14 Compile`$18; 
               Compile`$54 = 2.3094010767585034` Compile`$14 Compile`$18; 
               Compile`$57 = 
                1.1547005383792517` $CellContext`phi2$22384 Compile`$21 
                 Compile`$18; ((
                   Compile`$17 (
                    2.2500000000000013` (Compile`$19 + Compile`$20 + 
                    Compile`$22 + Compile`$24 + Compile`$25) Compile`$30 
                    Compile`$32 Compile`$34 + 
                    2.2500000000000013` (Compile`$36 + Compile`$37 + 
                    Compile`$22 + Compile`$38 + Compile`$25) Compile`$30 
                    Compile`$32 Compile`$40 + 
                    2.2500000000000013` (Compile`$36 + Compile`$20 + 
                    Compile`$42 + Compile`$38 + Compile`$27) Compile`$30 
                    Compile`$34 Compile`$40 + 
                    2.2500000000000013` (Compile`$19 + Compile`$37 + 
                    Compile`$42 + Compile`$24 + Compile`$27) Compile`$32 
                    Compile`$34 Compile`$40)^2)/
                  4 + (Compile`$17 (
                    2.2500000000000013` (Compile`$37 + Compile`$50 + 
                    Compile`$51 + Compile`$38 + Compile`$31) Compile`$30 
                    Compile`$32 Compile`$34 + 
                    2.2500000000000013` (Compile`$20 + Compile`$54 + 
                    Compile`$51 + Compile`$24 + Compile`$31) Compile`$30 
                    Compile`$32 Compile`$40 + 
                    2.2500000000000013` (Compile`$20 + Compile`$50 + 
                    Compile`$57 + Compile`$38 + Compile`$28) Compile`$30 
                    Compile`$34 Compile`$40 + 
                    2.2500000000000013` (Compile`$37 + Compile`$54 + 
                    Compile`$57 + Compile`$24 + Compile`$28) Compile`$32 
                    Compile`$34 Compile`$40)^2)/4)/2]], Evaluate], "Umax" -> 
           4/3|>, "Octahedron" -> <|"U" -> CompiledFunction[{11, 14.3, 5470}, {
              Blank[Real], 
              Blank[Real]}, {{3, 0, 0}, {3, 0, 1}, {3, 0, 
             13}}, {{1.8660254037844386`, {3, 0, 
               21}}, {-1.8660254037844386`, {3, 0, 20}}, {-1, {2, 0, 
               1}}, {-0.7071067811865475, {3, 0, 12}}, {
              5.224489795918365, {3, 0, 11}}, {-0.8660254037844386, {3, 0, 
               16}}, {-0.35355339059327373`, {3, 0, 7}}, {1, {2, 0, 0}}, {
              0.8660254037844386, {3, 0, 26}}, {-0.1339745962155614, {3, 0, 
               19}}, {0.35355339059327373`, {3, 0, 9}}, {
              0.7071067811865475, {3, 0, 14}}, {
              0.1339745962155614, {3, 0, 22}}}, {0, 2, 29, 0, 
             0}, {{40, 56, 3, 0, 0, 3, 0, 2}, {40, 56, 3, 0, 1, 3, 0, 3}, {10,
               0, 4}, {13, 4, 2, 3, 4}, {40, 60, 3, 0, 4, 3, 0, 5}, {10, 1, 
              6}, {13, 6, 2, 3, 6}, {16, 7, 6, 5, 8}, {16, 9, 6, 5, 10}, {16, 
              12, 0, 5, 13}, {16, 14, 1, 5, 15}, {16, 16, 6, 5, 17}, {10, 0, 
              18}, {13, 18, 13, 15, 17, 18}, {16, 19, 0, 5, 15}, {16, 20, 1, 
              5, 13}, {10, 0, 17}, {13, 17, 15, 13, 8, 17}, {16, 21, 0, 5, 
              13}, {16, 22, 1, 5, 15}, {10, 0, 23}, {13, 23, 13, 15, 8, 23}, {
              16, 20, 0, 5, 13}, {16, 19, 1, 5, 15}, {10, 0, 24}, {13, 24, 13,
               15, 10, 24}, {16, 22, 0, 5, 13}, {16, 21, 1, 5, 15}, {10, 0, 
              25}, {13, 25, 13, 15, 10, 25}, {16, 14, 0, 5, 13}, {16, 12, 1, 
              5, 15}, {16, 26, 6, 5, 27}, {10, 0, 28}, {13, 28, 13, 15, 27, 
              28}, {16, 11, 18, 17, 23, 24, 25, 28, 13}, {1}}, 
             Function[{$CellContext`phi1$22387, $CellContext`phi2$22387}, 
              
              Block[{Compile`$26, Compile`$35, Compile`$39, Compile`$41, 
                Compile`$45, Compile`$52, Compile`$61}, 
               Compile`$26 = $CellContext`phi1$22387^2; 
               Compile`$35 = $CellContext`phi2$22387^2; 
               Compile`$39 = 1 + Compile`$26 + Compile`$35; 
               Compile`$41 = Compile`$39^(-1); 
               Compile`$45 = -1 + Compile`$26 + Compile`$35; 
               Compile`$52 = (-0.35355339059327373`) Compile`$45 Compile`$41; 
               Compile`$61 = 0.35355339059327373` Compile`$45 Compile`$41; 
               5.224489795918365 (1 - 
                 0.7071067811865475 $CellContext`phi1$22387 Compile`$41 + 
                 0.7071067811865475 $CellContext`phi2$22387 Compile`$41 - 
                 0.8660254037844386 Compile`$45 Compile`$41) (1 - 
                 0.1339745962155614 $CellContext`phi1$22387 Compile`$41 - 
                 1.8660254037844386` $CellContext`phi2$22387 Compile`$41 + 
                 Compile`$52) (1 + 
                 1.8660254037844386` $CellContext`phi1$22387 Compile`$41 + 
                 0.1339745962155614 $CellContext`phi2$22387 Compile`$41 + 
                 Compile`$52) (1 - 
                 1.8660254037844386` $CellContext`phi1$22387 Compile`$41 - 
                 0.1339745962155614 $CellContext`phi2$22387 Compile`$41 + 
                 Compile`$61) (1 + 
                 0.1339745962155614 $CellContext`phi1$22387 Compile`$41 + 
                 1.8660254037844386` $CellContext`phi2$22387 Compile`$41 + 
                 Compile`$61) (1 + 
                 0.7071067811865475 $CellContext`phi1$22387 Compile`$41 - 
                 0.7071067811865475 $CellContext`phi2$22387 Compile`$41 + 
                 0.8660254037844386 Compile`$45 Compile`$41)]], Evaluate], 
           "V" -> CompiledFunction[{11, 14.3, 5470}, {
              Blank[Real], 
              Blank[Real]}, {{3, 0, 0}, {3, 0, 1}, {3, 0, 
             70}}, {{0.5, {3, 0, 63}}, {1.8660254037844386`, {3, 0, 24}}, {
              0.2679491924311228, {3, 0, 56}}, {-1.8660254037844386`, {3, 0, 
               21}}, {-1, {2, 0, 1}}, {-1.414213562373095, {3, 0, 
               43}}, {-0.2679491924311228, {3, 0, 
               46}}, {-1.7320508075688772`, {3, 0, 
               67}}, {-0.7071067811865475, {3, 0, 12}}, {
              5.224489795918365, {3, 0, 65}}, {
              1.7320508075688772`, {3, 0, 69}}, {-0.8660254037844386, {3, 0, 
               16}}, {-0.35355339059327373`, {3, 0, 8}}, {1, {2, 0, 0}}, {
              0.8660254037844386, {3, 0, 38}}, {-0.1339745962155614, {3, 0, 
               19}}, {0.35355339059327373`, {3, 0, 
               10}}, {-3.732050807568877, {3, 0, 53}}, {
              0.7071067811865475, {3, 0, 14}}, {
              1.414213562373095, {3, 0, 60}}, {
              3.732050807568877, {3, 0, 49}}, {
              0.1339745962155614, {3, 0, 26}}, {0.25, {3, 0, 64}}}, {0, 2, 78,
              0, 0}, {{40, 56, 3, 0, 0, 3, 0, 2}, {40, 56, 3, 0, 1, 3, 0, 
              3}, {10, 0, 4}, {13, 4, 2, 3, 4}, {40, 56, 3, 0, 4, 3, 0, 5}, {
              40, 60, 3, 0, 5, 3, 0, 6}, {40, 60, 3, 0, 4, 3, 0, 5}, {10, 1, 
              7}, {13, 7, 2, 3, 7}, {16, 8, 7, 5, 9}, {16, 10, 7, 5, 11}, {16,
               12, 0, 5, 13}, {16, 14, 1, 5, 15}, {16, 16, 7, 5, 17}, {10, 0, 
              18}, {13, 18, 13, 15, 17, 18}, {16, 19, 0, 5, 20}, {16, 21, 1, 
              5, 22}, {10, 0, 23}, {13, 23, 20, 22, 9, 23}, {16, 24, 0, 5, 
              25}, {16, 26, 1, 5, 27}, {10, 0, 28}, {13, 28, 25, 27, 9, 28}, {
              16, 21, 0, 5, 29}, {16, 19, 1, 5, 30}, {10, 0, 31}, {13, 31, 29,
               30, 11, 31}, {16, 14, 0, 5, 32}, {16, 12, 0, 7, 6, 33}, {16, 
              26, 0, 5, 34}, {16, 24, 1, 5, 35}, {10, 0, 36}, {13, 36, 34, 35,
               11, 36}, {16, 12, 1, 5, 37}, {16, 38, 7, 5, 39}, {10, 0, 40}, {
              13, 40, 32, 37, 39, 40}, {16, 14, 0, 7, 6, 41}, {40, 56, 3, 0, 
              4, 3, 0, 42}, {16, 43, 0, 1, 6, 44}, {16, 12, 5, 45}, {16, 46, 
              0, 1, 6, 47}, {16, 24, 5, 48}, {16, 49, 0, 1, 6, 50}, {16, 12, 
              1, 7, 6, 51}, {16, 19, 5, 52}, {16, 53, 0, 1, 6, 54}, {16, 26, 
              5, 55}, {16, 56, 0, 1, 6, 57}, {16, 14, 1, 7, 6, 58}, {16, 21, 
              5, 59}, {16, 60, 0, 1, 6, 61}, {16, 14, 5, 62}, {16, 43, 2, 6, 
              66}, {16, 67, 0, 7, 6, 68}, {16, 69, 0, 5, 70}, {13, 66, 61, 68,
               62, 70, 66}, {16, 65, 66, 18, 23, 28, 31, 36, 68}, {16, 46, 2, 
              6, 66}, {13, 66, 54, 33, 55, 32, 66}, {16, 65, 66, 18, 23, 28, 
              31, 40, 70}, {16, 49, 2, 6, 66}, {13, 66, 57, 33, 59, 32, 66}, {
              16, 65, 66, 18, 23, 28, 36, 40, 71}, {16, 53, 2, 6, 66}, {13, 
              66, 47, 41, 48, 13, 66}, {16, 65, 66, 18, 23, 31, 36, 40, 72}, {
              16, 56, 2, 6, 66}, {13, 66, 50, 41, 52, 13, 66}, {16, 65, 66, 
              18, 28, 31, 36, 40, 73}, {16, 60, 2, 6, 66}, {16, 69, 0, 7, 6, 
              74}, {16, 67, 0, 5, 75}, {13, 66, 44, 74, 45, 75, 66}, {16, 65, 
              66, 23, 28, 31, 36, 40, 74}, {13, 68, 70, 71, 72, 73, 74, 68}, {
              40, 56, 3, 0, 68, 3, 0, 70}, {16, 64, 42, 70, 68}, {16, 60, 3, 
              6, 70}, {16, 67, 1, 7, 6, 71}, {16, 69, 1, 5, 72}, {13, 44, 70, 
              71, 45, 72, 73}, {16, 65, 73, 18, 23, 28, 31, 36, 70}, {16, 53, 
              3, 6, 73}, {13, 47, 73, 51, 48, 15, 71}, {16, 65, 71, 18, 23, 
              28, 31, 40, 73}, {16, 56, 3, 6, 71}, {13, 50, 71, 51, 52, 15, 
              72}, {16, 65, 72, 18, 23, 28, 36, 40, 71}, {16, 46, 3, 6, 72}, {
              13, 54, 72, 58, 55, 37, 74}, {16, 65, 74, 18, 23, 31, 36, 40, 
              72}, {16, 49, 3, 6, 74}, {13, 57, 74, 58, 59, 37, 66}, {16, 65, 
              66, 18, 28, 31, 36, 40, 74}, {16, 43, 3, 6, 66}, {16, 69, 1, 7, 
              6, 75}, {16, 67, 1, 5, 76}, {13, 61, 66, 75, 62, 76, 77}, {16, 
              65, 77, 23, 28, 31, 36, 40, 66}, {13, 70, 73, 71, 72, 74, 66, 
              70}, {40, 56, 3, 0, 70, 3, 0, 73}, {16, 64, 42, 73, 70}, {13, 
              68, 70, 68}, {16, 63, 68, 70}, {1}}, 
             Function[{$CellContext`phi1$22387, $CellContext`phi2$22387}, 
              
              Block[{Compile`$43, Compile`$44, Compile`$46, Compile`$48, 
                Compile`$58, Compile`$55, Compile`$69, Compile`$76, 
                Compile`$63, Compile`$64, Compile`$65, Compile`$66, 
                Compile`$67, Compile`$68, Compile`$70, Compile`$71, 
                Compile`$72, Compile`$73, Compile`$74, Compile`$75, 
                Compile`$77, Compile`$86, Compile`$84, Compile`$78, 
                Compile`$79, Compile`$80, Compile`$88, Compile`$89, 
                Compile`$90, Compile`$99, Compile`$47, Compile`$109, 
                Compile`$111, Compile`$98, Compile`$100, Compile`$104, 
                Compile`$124, Compile`$105, Compile`$83, Compile`$85, 
                Compile`$93, Compile`$131, Compile`$94, Compile`$53, 
                Compile`$59}, Compile`$43 = $CellContext`phi1$22387^2; 
               Compile`$44 = $CellContext`phi2$22387^2; 
               Compile`$46 = 1 + Compile`$43 + Compile`$44; 
               Compile`$48 = Compile`$46^(-2); Compile`$58 = Compile`$46^(-1); 
               Compile`$55 = -1 + Compile`$43 + Compile`$44; 
               Compile`$69 = (-0.35355339059327373`) Compile`$55 Compile`$58; 
               Compile`$76 = 0.35355339059327373` Compile`$55 Compile`$58; 
               Compile`$63 = (-0.7071067811865475) $CellContext`phi1$22387 
                 Compile`$58; 
               Compile`$64 = 
                0.7071067811865475 $CellContext`phi2$22387 Compile`$58; 
               Compile`$65 = (-0.8660254037844386) Compile`$55 Compile`$58; 
               Compile`$66 = 1 + Compile`$63 + Compile`$64 + Compile`$65; 
               Compile`$67 = (-0.1339745962155614) $CellContext`phi1$22387 
                 Compile`$58; 
               Compile`$68 = (-1.8660254037844386`) $CellContext`phi2$22387 
                 Compile`$58; 
               Compile`$70 = 1 + Compile`$67 + Compile`$68 + Compile`$69; 
               Compile`$71 = 
                1.8660254037844386` $CellContext`phi1$22387 Compile`$58; 
               Compile`$72 = 
                0.1339745962155614 $CellContext`phi2$22387 Compile`$58; 
               Compile`$73 = 1 + Compile`$71 + Compile`$72 + Compile`$69; 
               Compile`$74 = (-1.8660254037844386`) $CellContext`phi1$22387 
                 Compile`$58; 
               Compile`$75 = (-0.1339745962155614) $CellContext`phi2$22387 
                 Compile`$58; 
               Compile`$77 = 1 + Compile`$74 + Compile`$75 + Compile`$76; 
               Compile`$86 = 
                0.7071067811865475 $CellContext`phi1$22387 Compile`$58; 
               Compile`$84 = (-0.7071067811865475) $CellContext`phi1$22387 
                 Compile`$55 Compile`$48; 
               Compile`$78 = 
                0.1339745962155614 $CellContext`phi1$22387 Compile`$58; 
               Compile`$79 = 
                1.8660254037844386` $CellContext`phi2$22387 Compile`$58; 
               Compile`$80 = 1 + Compile`$78 + Compile`$79 + Compile`$76; 
               Compile`$88 = (-0.7071067811865475) $CellContext`phi2$22387 
                 Compile`$58; 
               Compile`$89 = 0.8660254037844386 Compile`$55 Compile`$58; 
               Compile`$90 = 1 + Compile`$86 + Compile`$88 + Compile`$89; 
               Compile`$99 = 
                0.7071067811865475 $CellContext`phi1$22387 Compile`$55 
                 Compile`$48; Compile`$47 = Compile`$46^2; 
               Compile`$109 = (-1.414213562373095) $CellContext`phi1$22387 \
$CellContext`phi2$22387 Compile`$48; 
               Compile`$111 = (-0.7071067811865475) Compile`$58; 
               Compile`$98 = (-0.2679491924311228) $CellContext`phi1$22387 \
$CellContext`phi2$22387 Compile`$48; 
               Compile`$100 = 1.8660254037844386` Compile`$58; 
               Compile`$104 = 
                3.732050807568877 $CellContext`phi1$22387 \
$CellContext`phi2$22387 Compile`$48; 
               Compile`$124 = (-0.7071067811865475) $CellContext`phi2$22387 
                 Compile`$55 Compile`$48; 
               Compile`$105 = (-0.1339745962155614) Compile`$58; 
               Compile`$83 = (-3.732050807568877) $CellContext`phi1$22387 \
$CellContext`phi2$22387 Compile`$48; 
               Compile`$85 = 0.1339745962155614 Compile`$58; 
               Compile`$93 = 
                0.2679491924311228 $CellContext`phi1$22387 \
$CellContext`phi2$22387 Compile`$48; 
               Compile`$131 = 
                0.7071067811865475 $CellContext`phi2$22387 Compile`$55 
                 Compile`$48; 
               Compile`$94 = (-1.8660254037844386`) Compile`$58; 
               Compile`$53 = 
                1.414213562373095 $CellContext`phi1$22387 \
$CellContext`phi2$22387 Compile`$48; 
               Compile`$59 = 
                0.7071067811865475 
                 Compile`$58; ((
                   Compile`$47 (
                    5.224489795918365 ((-1.414213562373095) Compile`$43 
                    Compile`$48 + Compile`$53 - 
                    1.7320508075688772` $CellContext`phi1$22387 Compile`$55 
                    Compile`$48 + Compile`$59 + 
                    1.7320508075688772` $CellContext`phi1$22387 Compile`$58) 
                    Compile`$66 Compile`$70 Compile`$73 Compile`$77 
                    Compile`$80 + 
                    5.224489795918365 ((-0.2679491924311228) Compile`$43 
                    Compile`$48 + Compile`$83 + Compile`$84 + Compile`$85 + 
                    Compile`$86) Compile`$66 Compile`$70 Compile`$73 
                    Compile`$77 Compile`$90 + 
                    5.224489795918365 (
                    3.732050807568877 Compile`$43 Compile`$48 + Compile`$93 + 
                    Compile`$84 + Compile`$94 + Compile`$86) Compile`$66 
                    Compile`$70 Compile`$73 Compile`$80 Compile`$90 + 
                    5.224489795918365 ((-3.732050807568877) Compile`$43 
                    Compile`$48 + Compile`$98 + Compile`$99 + Compile`$100 + 
                    Compile`$63) Compile`$66 Compile`$70 Compile`$77 
                    Compile`$80 Compile`$90 + 
                    5.224489795918365 (
                    0.2679491924311228 Compile`$43 Compile`$48 + Compile`$104 + 
                    Compile`$99 + Compile`$105 + Compile`$63) Compile`$66 
                    Compile`$73 Compile`$77 Compile`$80 Compile`$90 + 
                    5.224489795918365 (
                    1.414213562373095 Compile`$43 Compile`$48 + Compile`$109 + 
                    1.7320508075688772` $CellContext`phi1$22387 Compile`$55 
                    Compile`$48 + Compile`$111 - 
                    1.7320508075688772` $CellContext`phi1$22387 Compile`$58) 
                    Compile`$70 Compile`$73 Compile`$77 Compile`$80 
                    Compile`$90)^2)/
                  4 + (Compile`$47 (
                    5.224489795918365 (Compile`$109 + 
                    1.414213562373095 Compile`$44 Compile`$48 - 
                    1.7320508075688772` $CellContext`phi2$22387 Compile`$55 
                    Compile`$48 + Compile`$111 + 
                    1.7320508075688772` $CellContext`phi2$22387 Compile`$58) 
                    Compile`$66 Compile`$70 Compile`$73 Compile`$77 
                    Compile`$80 + 
                    5.224489795918365 (Compile`$98 - 3.732050807568877 
                    Compile`$44 Compile`$48 + Compile`$124 + Compile`$100 + 
                    Compile`$64) Compile`$66 Compile`$70 Compile`$73 
                    Compile`$77 Compile`$90 + 
                    5.224489795918365 (Compile`$104 + 
                    0.2679491924311228 Compile`$44 Compile`$48 + Compile`$124 + 
                    Compile`$105 + Compile`$64) Compile`$66 Compile`$70 
                    Compile`$73 Compile`$80 Compile`$90 + 
                    5.224489795918365 (Compile`$83 - 0.2679491924311228 
                    Compile`$44 Compile`$48 + Compile`$131 + Compile`$85 + 
                    Compile`$88) Compile`$66 Compile`$70 Compile`$77 
                    Compile`$80 Compile`$90 + 
                    5.224489795918365 (Compile`$93 + 
                    3.732050807568877 Compile`$44 Compile`$48 + Compile`$131 + 
                    Compile`$94 + Compile`$88) Compile`$66 Compile`$73 
                    Compile`$77 Compile`$80 Compile`$90 + 
                    5.224489795918365 (Compile`$53 - 1.414213562373095 
                    Compile`$44 Compile`$48 + 
                    1.7320508075688772` $CellContext`phi2$22387 Compile`$55 
                    Compile`$48 + Compile`$59 - 
                    1.7320508075688772` $CellContext`phi2$22387 Compile`$58) 
                    Compile`$70 Compile`$73 Compile`$77 Compile`$80 
                    Compile`$90)^2)/4)/2]], Evaluate], "Umax" -> 1.548|>, 
         "Cube" -> <|"U" -> CompiledFunction[{11, 14.3, 5470}, {
              Blank[Real], 
              Blank[Real]}, {{3, 0, 0}, {3, 0, 1}, {3, 0, 
             15}}, {{10.239999999999998`, {3, 0, 11}}, {
              0.5, {3, 0, 9}}, {-0.5917517095361372, {3, 0, 14}}, {-1, {2, 0, 
               1}}, {-0.5, {3, 0, 7}}, {-1.5629488288431146`, {3, 0, 19}}, {
              0.5917517095361372, {3, 0, 12}}, {-0.09175170953613698, {3, 0, 
               26}}, {-0.908248290463863, {3, 0, 16}}, {
              1.4082482904638631`, {3, 0, 25}}, {
              0.7464522479153887, {3, 0, 21}}, {
              1, {2, 0, 0}}, {-1.4082482904638631`, {3, 0, 24}}, {
              0.908248290463863, {3, 0, 32}}, {-0.7464522479153887, {3, 0, 
               20}}, {1.5629488288431146`, {3, 0, 22}}, {
              0.09175170953613698, {3, 0, 29}}}, {0, 2, 35, 0, 
             0}, {{40, 56, 3, 0, 0, 3, 0, 2}, {40, 56, 3, 0, 1, 3, 0, 3}, {10,
               0, 4}, {13, 4, 2, 3, 4}, {40, 60, 3, 0, 4, 3, 0, 5}, {10, 1, 
              6}, {13, 6, 2, 3, 6}, {16, 7, 6, 5, 8}, {16, 9, 6, 5, 10}, {16, 
              12, 0, 5, 13}, {16, 14, 1, 5, 15}, {16, 16, 6, 5, 17}, {10, 0, 
              18}, {13, 18, 13, 15, 17, 18}, {16, 19, 0, 5, 15}, {16, 20, 1, 
              5, 13}, {10, 0, 17}, {13, 17, 15, 13, 8, 17}, {16, 21, 0, 5, 
              13}, {16, 22, 1, 5, 15}, {10, 0, 23}, {13, 23, 13, 15, 8, 23}, {
              16, 24, 0, 5, 15}, {16, 25, 1, 5, 13}, {16, 26, 6, 5, 27}, {10, 
              0, 28}, {13, 28, 15, 13, 27, 28}, {16, 25, 0, 5, 15}, {16, 24, 
              1, 5, 13}, {16, 29, 6, 5, 27}, {10, 0, 30}, {13, 30, 15, 13, 27,
               30}, {16, 20, 0, 5, 15}, {16, 19, 1, 5, 13}, {10, 0, 27}, {13, 
              27, 15, 13, 10, 27}, {16, 22, 0, 5, 15}, {16, 21, 1, 5, 13}, {
              10, 0, 31}, {13, 31, 15, 13, 10, 31}, {16, 14, 0, 5, 15}, {16, 
              12, 1, 5, 13}, {16, 32, 6, 5, 33}, {10, 0, 34}, {13, 34, 15, 13,
               33, 34}, {16, 11, 18, 17, 23, 28, 30, 27, 31, 34, 15}, {1}}, 
             Function[{$CellContext`phi1$22390, $CellContext`phi2$22390}, 
              
              Block[{Compile`$49, Compile`$56, Compile`$60, Compile`$62, 
                Compile`$87, Compile`$97, Compile`$119}, 
               Compile`$49 = $CellContext`phi1$22390^2; 
               Compile`$56 = $CellContext`phi2$22390^2; 
               Compile`$60 = 1 + Compile`$49 + Compile`$56; 
               Compile`$62 = Compile`$60^(-1); 
               Compile`$87 = -1 + Compile`$49 + Compile`$56; 
               Compile`$97 = (-0.5) Compile`$87 Compile`$62; 
               Compile`$119 = 0.5 Compile`$87 Compile`$62; 
               10.239999999999998` (1 + 
                 0.5917517095361372 $CellContext`phi1$22390 Compile`$62 - 
                 0.5917517095361372 $CellContext`phi2$22390 Compile`$62 - 
                 0.908248290463863 Compile`$87 Compile`$62) (1 - 
                 1.5629488288431146` $CellContext`phi1$22390 Compile`$62 - 
                 0.7464522479153887 $CellContext`phi2$22390 Compile`$62 + 
                 Compile`$97) (1 + 
                 0.7464522479153887 $CellContext`phi1$22390 Compile`$62 + 
                 1.5629488288431146` $CellContext`phi2$22390 Compile`$62 + 
                 Compile`$97) (1 - 
                 1.4082482904638631` $CellContext`phi1$22390 Compile`$62 + 
                 1.4082482904638631` $CellContext`phi2$22390 Compile`$62 - 
                 0.09175170953613698 Compile`$87 Compile`$62) (1 + 
                 1.4082482904638631` $CellContext`phi1$22390 Compile`$62 - 
                 1.4082482904638631` $CellContext`phi2$22390 Compile`$62 + 
                 0.09175170953613698 Compile`$87 Compile`$62) (1 - 
                 0.7464522479153887 $CellContext`phi1$22390 Compile`$62 - 
                 1.5629488288431146` $CellContext`phi2$22390 Compile`$62 + 
                 Compile`$119) (1 + 
                 1.5629488288431146` $CellContext`phi1$22390 Compile`$62 + 
                 0.7464522479153887 $CellContext`phi2$22390 Compile`$62 + 
                 Compile`$119) (1 - 
                 0.5917517095361372 $CellContext`phi1$22390 Compile`$62 + 
                 0.5917517095361372 $CellContext`phi2$22390 Compile`$62 + 
                 0.908248290463863 Compile`$87 Compile`$62)]], Evaluate], "V" -> 
           CompiledFunction[{11, 14.3, 5470}, {
              Blank[Real], 
              Blank[Real]}, {{3, 0, 0}, {3, 0, 1}, {3, 0, 
             93}}, {{10.239999999999998`, {3, 0, 88}}, {
              0.5, {3, 0, 10}}, {-0.5917517095361372, {3, 0, 14}}, {
              3.125897657686229, {3, 0, 79}}, {-0.18350341907227397`, {3, 0, 
               95}}, {1.1835034190722744`, {3, 0, 59}}, {
              2.8164965809277263`, {3, 0, 73}}, {-3.125897657686229, {3, 0, 
               62}}, {-1, {2, 0, 1}}, {-1., {3, 0, 44}}, {
              1.4929044958307773`, {3, 0, 65}}, {
              0.18350341907227397`, {3, 0, 97}}, {-0.5, {3, 0, 
               8}}, {-1.5629488288431146`, {3, 0, 19}}, {
              0.5917517095361372, {3, 0, 12}}, {-0.09175170953613698, {3, 0, 
               33}}, {1.816496580927726, {3, 0, 92}}, {-0.908248290463863, {3,
                0, 16}}, {
              1.4082482904638631`, {3, 0, 31}}, {-2.8164965809277263`, {3, 0, 
               70}}, {0.7464522479153887, {3, 0, 24}}, {1, {2, 0, 0}}, {
              1., {3, 0, 46}}, {-1.4082482904638631`, {3, 0, 
               29}}, {-1.4929044958307773`, {3, 0, 76}}, {
              0.908248290463863, {3, 0, 53}}, {-0.7464522479153887, {3, 0, 
               21}}, {1.5629488288431146`, {3, 0, 26}}, {
              0.09175170953613698, {3, 0, 38}}, {
              0.25, {3, 0, 87}}, {-1.816496580927726, {3, 0, 
               90}}, {-1.1835034190722744`, {3, 0, 84}}}, {0, 2, 105, 0, 
             0}, {{40, 56, 3, 0, 0, 3, 0, 2}, {40, 56, 3, 0, 1, 3, 0, 3}, {10,
               0, 4}, {13, 4, 2, 3, 4}, {40, 56, 3, 0, 4, 3, 0, 5}, {40, 60, 
              3, 0, 5, 3, 0, 6}, {40, 60, 3, 0, 4, 3, 0, 5}, {10, 1, 7}, {13, 
              7, 2, 3, 7}, {16, 8, 7, 5, 9}, {16, 10, 7, 5, 11}, {16, 12, 0, 
              5, 13}, {16, 14, 1, 5, 15}, {16, 16, 7, 5, 17}, {10, 0, 18}, {
              13, 18, 13, 15, 17, 18}, {16, 19, 0, 5, 20}, {16, 21, 1, 5, 
              22}, {10, 0, 23}, {13, 23, 20, 22, 9, 23}, {16, 24, 0, 5, 25}, {
              16, 26, 1, 5, 27}, {10, 0, 28}, {13, 28, 25, 27, 9, 28}, {16, 
              29, 0, 5, 30}, {16, 31, 1, 5, 32}, {16, 33, 7, 5, 34}, {10, 0, 
              35}, {13, 35, 30, 32, 34, 35}, {16, 31, 0, 5, 36}, {16, 29, 1, 
              5, 37}, {16, 38, 7, 5, 39}, {10, 0, 40}, {13, 40, 36, 37, 39, 
              40}, {16, 21, 0, 5, 41}, {16, 19, 1, 5, 42}, {10, 0, 43}, {13, 
              43, 41, 42, 11, 43}, {16, 44, 0, 7, 6, 45}, {16, 46, 0, 5, 
              47}, {16, 26, 0, 5, 48}, {16, 24, 1, 5, 49}, {10, 0, 50}, {13, 
              50, 48, 49, 11, 50}, {16, 14, 0, 5, 51}, {16, 12, 1, 5, 52}, {
              16, 53, 7, 5, 54}, {10, 0, 55}, {13, 55, 51, 52, 54, 55}, {16, 
              46, 0, 7, 6, 56}, {16, 44, 0, 5, 57}, {40, 56, 3, 0, 4, 3, 0, 
              58}, {16, 59, 0, 1, 6, 60}, {16, 12, 5, 61}, {16, 62, 0, 1, 6, 
              63}, {16, 24, 5, 64}, {16, 65, 0, 1, 6, 66}, {16, 44, 1, 7, 6, 
              67}, {16, 19, 5, 68}, {16, 46, 1, 5, 69}, {16, 70, 0, 1, 6, 
              71}, {16, 29, 5, 72}, {16, 73, 0, 1, 6, 74}, {16, 31, 5, 75}, {
              16, 76, 0, 1, 6, 77}, {16, 26, 5, 78}, {16, 79, 0, 1, 6, 80}, {
              16, 46, 1, 7, 6, 81}, {16, 21, 5, 82}, {16, 44, 1, 5, 83}, {16, 
              84, 0, 1, 6, 85}, {16, 14, 5, 86}, {16, 59, 2, 6, 89}, {16, 90, 
              0, 7, 6, 91}, {16, 92, 0, 5, 93}, {13, 89, 85, 91, 86, 93, 
              89}, {16, 88, 89, 18, 23, 28, 35, 40, 43, 50, 91}, {16, 62, 2, 
              6, 89}, {13, 89, 77, 45, 78, 47, 89}, {16, 88, 89, 18, 23, 28, 
              35, 40, 43, 55, 93}, {16, 65, 2, 6, 89}, {13, 89, 80, 45, 82, 
              47, 89}, {16, 88, 89, 18, 23, 28, 35, 40, 50, 55, 94}, {16, 70, 
              2, 6, 89}, {16, 95, 0, 7, 6, 96}, {16, 97, 0, 5, 98}, {13, 89, 
              74, 96, 75, 98, 89}, {16, 88, 89, 18, 23, 28, 35, 43, 50, 55, 
              96}, {16, 73, 2, 6, 89}, {16, 97, 0, 7, 6, 98}, {16, 95, 0, 5, 
              99}, {13, 89, 71, 98, 72, 99, 89}, {16, 88, 89, 18, 23, 28, 40, 
              43, 50, 55, 98}, {16, 76, 2, 6, 89}, {13, 89, 63, 56, 64, 57, 
              89}, {16, 88, 89, 18, 23, 35, 40, 43, 50, 55, 99}, {16, 79, 2, 
              6, 89}, {13, 89, 66, 56, 68, 57, 89}, {16, 88, 89, 18, 28, 35, 
              40, 43, 50, 55, 100}, {16, 84, 2, 6, 89}, {16, 92, 0, 7, 6, 
              101}, {16, 90, 0, 5, 102}, {13, 89, 60, 101, 61, 102, 89}, {16, 
              88, 89, 23, 28, 35, 40, 43, 50, 55, 101}, {13, 91, 93, 94, 96, 
              98, 99, 100, 101, 91}, {40, 56, 3, 0, 91, 3, 0, 93}, {16, 87, 
              58, 93, 91}, {16, 84, 3, 6, 93}, {16, 90, 1, 7, 6, 94}, {16, 92,
               1, 5, 96}, {13, 60, 93, 94, 61, 96, 98}, {16, 88, 98, 18, 23, 
              28, 35, 40, 43, 50, 93}, {16, 76, 3, 6, 98}, {13, 63, 98, 67, 
              64, 69, 94}, {16, 88, 94, 18, 23, 28, 35, 40, 43, 55, 98}, {16, 
              79, 3, 6, 94}, {13, 66, 94, 67, 68, 69, 96}, {16, 88, 96, 18, 
              23, 28, 35, 40, 50, 55, 94}, {16, 73, 3, 6, 96}, {16, 95, 1, 7, 
              6, 99}, {16, 97, 1, 5, 100}, {13, 71, 96, 99, 72, 100, 101}, {
              16, 88, 101, 18, 23, 28, 35, 43, 50, 55, 96}, {16, 70, 3, 6, 
              101}, {16, 97, 1, 7, 6, 99}, {16, 95, 1, 5, 100}, {13, 74, 101, 
              99, 75, 100, 89}, {16, 88, 89, 18, 23, 28, 40, 43, 50, 55, 
              101}, {16, 62, 3, 6, 89}, {13, 77, 89, 81, 78, 83, 99}, {16, 88,
               99, 18, 23, 35, 40, 43, 50, 55, 89}, {16, 65, 3, 6, 99}, {13, 
              80, 99, 81, 82, 83, 100}, {16, 88, 100, 18, 28, 35, 40, 43, 50, 
              55, 99}, {16, 59, 3, 6, 100}, {16, 92, 1, 7, 6, 102}, {16, 90, 
              1, 5, 103}, {13, 85, 100, 102, 86, 103, 104}, {16, 88, 104, 23, 
              28, 35, 40, 43, 50, 55, 100}, {13, 93, 98, 94, 96, 101, 89, 99, 
              100, 93}, {40, 56, 3, 0, 93, 3, 0, 98}, {16, 87, 58, 98, 93}, {
              13, 91, 93, 91}, {16, 10, 91, 93}, {1}}, 
             Function[{$CellContext`phi1$22390, $CellContext`phi2$22390}, 
              
              Block[{Compile`$81, Compile`$82, Compile`$91, Compile`$95, 
                Compile`$106, Compile`$102, Compile`$118, Compile`$136, 
                Compile`$112, Compile`$113, Compile`$114, Compile`$115, 
                Compile`$116, Compile`$117, Compile`$120, Compile`$121, 
                Compile`$122, Compile`$123, Compile`$125, Compile`$126, 
                Compile`$127, Compile`$128, Compile`$129, Compile`$130, 
                Compile`$132, Compile`$133, Compile`$134, Compile`$135, 
                Compile`$137, Compile`$144, Compile`$146, Compile`$138, 
                Compile`$139, Compile`$140, Compile`$148, Compile`$149, 
                Compile`$150, Compile`$151, Compile`$174, Compile`$176, 
                Compile`$92, Compile`$185, Compile`$187, Compile`$173, 
                Compile`$175, Compile`$180, Compile`$200, Compile`$181, 
                Compile`$201, Compile`$166, Compile`$168, Compile`$159, 
                Compile`$161, Compile`$143, Compile`$145, Compile`$154, 
                Compile`$218, Compile`$155, Compile`$219, Compile`$101, 
                Compile`$107}, Compile`$81 = $CellContext`phi1$22390^2; 
               Compile`$82 = $CellContext`phi2$22390^2; 
               Compile`$91 = 1 + Compile`$81 + Compile`$82; 
               Compile`$95 = Compile`$91^(-2); 
               Compile`$106 = Compile`$91^(-1); 
               Compile`$102 = -1 + Compile`$81 + Compile`$82; 
               Compile`$118 = (-0.5) Compile`$102 Compile`$106; 
               Compile`$136 = 0.5 Compile`$102 Compile`$106; 
               Compile`$112 = 
                0.5917517095361372 $CellContext`phi1$22390 Compile`$106; 
               Compile`$113 = (-0.5917517095361372) $CellContext`phi2$22390 
                 Compile`$106; 
               Compile`$114 = (-0.908248290463863) Compile`$102 Compile`$106; 
               Compile`$115 = 1 + Compile`$112 + Compile`$113 + Compile`$114; 
               Compile`$116 = (-1.5629488288431146`) $CellContext`phi1$22390 
                 Compile`$106; 
               Compile`$117 = (-0.7464522479153887) $CellContext`phi2$22390 
                 Compile`$106; 
               Compile`$120 = 1 + Compile`$116 + Compile`$117 + Compile`$118; 
               Compile`$121 = 
                0.7464522479153887 $CellContext`phi1$22390 Compile`$106; 
               Compile`$122 = 
                1.5629488288431146` $CellContext`phi2$22390 Compile`$106; 
               Compile`$123 = 1 + Compile`$121 + Compile`$122 + Compile`$118; 
               Compile`$125 = (-1.4082482904638631`) $CellContext`phi1$22390 
                 Compile`$106; 
               Compile`$126 = 
                1.4082482904638631` $CellContext`phi2$22390 Compile`$106; 
               Compile`$127 = (-0.09175170953613698) Compile`$102 
                 Compile`$106; 
               Compile`$128 = 1 + Compile`$125 + Compile`$126 + Compile`$127; 
               Compile`$129 = 
                1.4082482904638631` $CellContext`phi1$22390 Compile`$106; 
               Compile`$130 = (-1.4082482904638631`) $CellContext`phi2$22390 
                 Compile`$106; 
               Compile`$132 = 0.09175170953613698 Compile`$102 Compile`$106; 
               Compile`$133 = 1 + Compile`$129 + Compile`$130 + Compile`$132; 
               Compile`$134 = (-0.7464522479153887) $CellContext`phi1$22390 
                 Compile`$106; 
               Compile`$135 = (-1.5629488288431146`) $CellContext`phi2$22390 
                 Compile`$106; 
               Compile`$137 = 1 + Compile`$134 + Compile`$135 + Compile`$136; 
               Compile`$144 = -$CellContext`phi1$22390 Compile`$102 
                 Compile`$95; 
               Compile`$146 = 1. $CellContext`phi1$22390 Compile`$106; 
               Compile`$138 = 
                1.5629488288431146` $CellContext`phi1$22390 Compile`$106; 
               Compile`$139 = 
                0.7464522479153887 $CellContext`phi2$22390 Compile`$106; 
               Compile`$140 = 1 + Compile`$138 + Compile`$139 + Compile`$136; 
               Compile`$148 = (-0.5917517095361372) $CellContext`phi1$22390 
                 Compile`$106; 
               Compile`$149 = 
                0.5917517095361372 $CellContext`phi2$22390 Compile`$106; 
               Compile`$150 = 0.908248290463863 Compile`$102 Compile`$106; 
               Compile`$151 = 1 + Compile`$148 + Compile`$149 + Compile`$150; 
               Compile`$174 = 
                1. $CellContext`phi1$22390 Compile`$102 Compile`$95; 
               Compile`$176 = -$CellContext`phi1$22390 Compile`$106; 
               Compile`$92 = Compile`$91^2; 
               Compile`$185 = 
                1.1835034190722744` $CellContext`phi1$22390 \
$CellContext`phi2$22390 Compile`$95; 
               Compile`$187 = 0.5917517095361372 Compile`$106; 
               Compile`$173 = (-3.125897657686229) $CellContext`phi1$22390 \
$CellContext`phi2$22390 Compile`$95; 
               Compile`$175 = 0.7464522479153887 Compile`$106; 
               Compile`$180 = 
                1.4929044958307773` $CellContext`phi1$22390 \
$CellContext`phi2$22390 Compile`$95; 
               Compile`$200 = -$CellContext`phi2$22390 Compile`$102 
                 Compile`$95; 
               Compile`$181 = (-1.5629488288431146`) Compile`$106; 
               Compile`$201 = 1. $CellContext`phi2$22390 Compile`$106; 
               Compile`$166 = (-2.8164965809277263`) $CellContext`phi1$22390 \
$CellContext`phi2$22390 Compile`$95; 
               Compile`$168 = (-1.4082482904638631`) Compile`$106; 
               Compile`$159 = 
                2.8164965809277263` $CellContext`phi1$22390 \
$CellContext`phi2$22390 Compile`$95; 
               Compile`$161 = 1.4082482904638631` Compile`$106; 
               Compile`$143 = (-1.4929044958307773`) $CellContext`phi1$22390 \
$CellContext`phi2$22390 Compile`$95; 
               Compile`$145 = 1.5629488288431146` Compile`$106; 
               Compile`$154 = 
                3.125897657686229 $CellContext`phi1$22390 \
$CellContext`phi2$22390 Compile`$95; 
               Compile`$218 = 
                1. $CellContext`phi2$22390 Compile`$102 Compile`$95; 
               Compile`$155 = (-0.7464522479153887) Compile`$106; 
               Compile`$219 = -$CellContext`phi2$22390 Compile`$106; 
               Compile`$101 = (-1.1835034190722744`) $CellContext`phi1$22390 \
$CellContext`phi2$22390 Compile`$95; 
               Compile`$107 = (-0.5917517095361372) 
                 Compile`$106; ((
                   Compile`$92 (
                    10.239999999999998` (
                    1.1835034190722744` Compile`$81 Compile`$95 + 
                    Compile`$101 - 1.816496580927726 $CellContext`phi1$22390 
                    Compile`$102 Compile`$95 + Compile`$107 + 
                    1.816496580927726 $CellContext`phi1$22390 Compile`$106) 
                    Compile`$115 Compile`$120 Compile`$123 Compile`$128 
                    Compile`$133 Compile`$137 Compile`$140 + 
                    10.239999999999998` ((-3.125897657686229) Compile`$81 
                    Compile`$95 + Compile`$143 + Compile`$144 + Compile`$145 + 
                    Compile`$146) Compile`$115 Compile`$120 Compile`$123 
                    Compile`$128 Compile`$133 Compile`$137 Compile`$151 + 
                    10.239999999999998` (
                    1.4929044958307773` Compile`$81 Compile`$95 + 
                    Compile`$154 + Compile`$144 + Compile`$155 + Compile`$146)
                     Compile`$115 Compile`$120 Compile`$123 Compile`$128 
                    Compile`$133 Compile`$140 Compile`$151 + 
                    10.239999999999998` ((-2.8164965809277263`) Compile`$81 
                    Compile`$95 + Compile`$159 - 
                    0.18350341907227397` $CellContext`phi1$22390 Compile`$102 
                    Compile`$95 + Compile`$161 + 
                    0.18350341907227397` $CellContext`phi1$22390 Compile`$106)
                     Compile`$115 Compile`$120 Compile`$123 Compile`$128 
                    Compile`$137 Compile`$140 Compile`$151 + 
                    10.239999999999998` (
                    2.8164965809277263` Compile`$81 Compile`$95 + 
                    Compile`$166 + 
                    0.18350341907227397` $CellContext`phi1$22390 Compile`$102 
                    Compile`$95 + Compile`$168 - 
                    0.18350341907227397` $CellContext`phi1$22390 Compile`$106)
                     Compile`$115 Compile`$120 Compile`$123 Compile`$133 
                    Compile`$137 Compile`$140 Compile`$151 + 
                    10.239999999999998` ((-1.4929044958307773`) Compile`$81 
                    Compile`$95 + Compile`$173 + Compile`$174 + Compile`$175 + 
                    Compile`$176) Compile`$115 Compile`$120 Compile`$128 
                    Compile`$133 Compile`$137 Compile`$140 Compile`$151 + 
                    10.239999999999998` (
                    3.125897657686229 Compile`$81 Compile`$95 + Compile`$180 + 
                    Compile`$174 + Compile`$181 + Compile`$176) Compile`$115 
                    Compile`$123 Compile`$128 Compile`$133 Compile`$137 
                    Compile`$140 Compile`$151 + 
                    10.239999999999998` ((-1.1835034190722744`) Compile`$81 
                    Compile`$95 + Compile`$185 + 
                    1.816496580927726 $CellContext`phi1$22390 Compile`$102 
                    Compile`$95 + Compile`$187 - 
                    1.816496580927726 $CellContext`phi1$22390 Compile`$106) 
                    Compile`$120 Compile`$123 Compile`$128 Compile`$133 
                    Compile`$137 Compile`$140 Compile`$151)^2)/
                  4 + (Compile`$92 (
                    10.239999999999998` (Compile`$185 - 1.1835034190722744` 
                    Compile`$82 Compile`$95 - 
                    1.816496580927726 $CellContext`phi2$22390 Compile`$102 
                    Compile`$95 + Compile`$187 + 
                    1.816496580927726 $CellContext`phi2$22390 Compile`$106) 
                    Compile`$115 Compile`$120 Compile`$123 Compile`$128 
                    Compile`$133 Compile`$137 Compile`$140 + 
                    10.239999999999998` (Compile`$173 - 1.4929044958307773` 
                    Compile`$82 Compile`$95 + Compile`$200 + Compile`$175 + 
                    Compile`$201) Compile`$115 Compile`$120 Compile`$123 
                    Compile`$128 Compile`$133 Compile`$137 Compile`$151 + 
                    10.239999999999998` (Compile`$180 + 
                    3.125897657686229 Compile`$82 Compile`$95 + Compile`$200 + 
                    Compile`$181 + Compile`$201) Compile`$115 Compile`$120 
                    Compile`$123 Compile`$128 Compile`$133 Compile`$140 
                    Compile`$151 + 
                    10.239999999999998` (Compile`$166 + 
                    2.8164965809277263` Compile`$82 Compile`$95 - 
                    0.18350341907227397` $CellContext`phi2$22390 Compile`$102 
                    Compile`$95 + Compile`$168 + 
                    0.18350341907227397` $CellContext`phi2$22390 Compile`$106)
                     Compile`$115 Compile`$120 Compile`$123 Compile`$128 
                    Compile`$137 Compile`$140 Compile`$151 + 
                    10.239999999999998` (Compile`$159 - 2.8164965809277263` 
                    Compile`$82 Compile`$95 + 
                    0.18350341907227397` $CellContext`phi2$22390 Compile`$102 
                    Compile`$95 + Compile`$161 - 
                    0.18350341907227397` $CellContext`phi2$22390 Compile`$106)
                     Compile`$115 Compile`$120 Compile`$123 Compile`$133 
                    Compile`$137 Compile`$140 Compile`$151 + 
                    10.239999999999998` (Compile`$143 - 3.125897657686229 
                    Compile`$82 Compile`$95 + Compile`$218 + Compile`$145 + 
                    Compile`$219) Compile`$115 Compile`$120 Compile`$128 
                    Compile`$133 Compile`$137 Compile`$140 Compile`$151 + 
                    10.239999999999998` (Compile`$154 + 
                    1.4929044958307773` Compile`$82 Compile`$95 + 
                    Compile`$218 + Compile`$155 + Compile`$219) Compile`$115 
                    Compile`$123 Compile`$128 Compile`$133 Compile`$137 
                    Compile`$140 Compile`$151 + 
                    10.239999999999998` (Compile`$101 + 
                    1.1835034190722744` Compile`$82 Compile`$95 + 
                    1.816496580927726 $CellContext`phi2$22390 Compile`$102 
                    Compile`$95 + Compile`$107 - 
                    1.816496580927726 $CellContext`phi2$22390 Compile`$106) 
                    Compile`$120 Compile`$123 Compile`$128 Compile`$133 
                    Compile`$137 Compile`$140 Compile`$151)^2)/4)/2]], 
             Evaluate], "Umax" -> 2.023|>, 
         "Icosahedron" -> <|"U" -> CompiledFunction[{11, 14.3, 5470}, {
              Blank[Real], 
              Blank[Real]}, {{3, 0, 0}, {3, 0, 1}, {3, 0, 
             11}}, {{-0.257782435890779, {3, 0, 40}}, {
              0.15454602115124183`, {3, 0, 41}}, {-0.11487646027368059`, {3, 
               0, 33}}, {-0.5508111925525423, {3, 0, 
               21}}, {-0.4857136330295899, {3, 0, 52}}, {-1, {2, 0, 1}}, {
              0.257782435890779, {3, 0, 28}}, {
              0.5508111925525423, {3, 0, 49}}, {
              0.4857136330295899, {3, 0, 16}}, {
              1.959084052594859, {3, 0, 39}}, {-1.21558798367449, {3, 0, 
               51}}, {0.7560469761587875, {3, 0, 53}}, {-0.7560469761587875, {
               3, 0, 17}}, {-0.15454602115124183`, {3, 0, 29}}, {
              0.5310663415434267, {3, 0, 47}}, {
              1.21558798367449, {3, 0, 15}}, {-1.6577706315987935`, {3, 0, 
               32}}, {1.6577706315987935`, {3, 0, 36}}, {
              0.6719355684716646, {3, 0, 55}}, {-0.6719355684716646, {3, 0, 
               8}}, {1, {2, 0, 0}}, {0.8670620122047428, {3, 0, 23}}, {
              0.48662449473386504`, {3, 0, 45}}, {-0.3795266557666026, {3, 0, 
               10}}, {-1.5825285657816939`, {3, 0, 
               48}}, {-0.48662449473386504`, {3, 0, 25}}, {
              1.5169014046705556`, {3, 0, 44}}, {-0.9225592270127267, {3, 0, 
               12}}, {-0.8670620122047428, {3, 0, 43}}, {
              0.11487646027368059`, {3, 0, 37}}, {
              30.616258720840584`, {3, 0, 7}}, {
              1.5825285657816939`, {3, 0, 20}}, {-1.0949932093435537`, {3, 0, 
               31}}, {-1.959084052594859, {3, 0, 27}}, {
              0.9225592270127267, {3, 0, 57}}, {-1.5169014046705556`, {3, 0, 
               24}}, {0.3795266557666026, {3, 0, 56}}, {
              1.0949932093435537`, {3, 0, 35}}, {-0.5310663415434267, {3, 0, 
               19}}}, {0, 2, 59, 0, 
             0}, {{40, 56, 3, 0, 0, 3, 0, 2}, {40, 56, 3, 0, 1, 3, 0, 3}, {10,
               0, 4}, {13, 4, 2, 3, 4}, {40, 60, 3, 0, 4, 3, 0, 5}, {10, 1, 
              6}, {13, 6, 2, 3, 6}, {16, 8, 0, 5, 9}, {16, 10, 1, 5, 11}, {16,
               12, 6, 5, 13}, {10, 0, 14}, {13, 14, 9, 11, 13, 14}, {16, 15, 
              0, 5, 11}, {16, 16, 1, 5, 9}, {16, 17, 6, 5, 13}, {10, 0, 18}, {
              13, 18, 11, 9, 13, 18}, {16, 19, 0, 5, 9}, {16, 20, 1, 5, 11}, {
              16, 21, 6, 5, 13}, {10, 0, 22}, {13, 22, 9, 11, 13, 22}, {16, 
              23, 0, 5, 11}, {16, 24, 1, 5, 9}, {16, 25, 6, 5, 13}, {10, 0, 
              26}, {13, 26, 11, 9, 13, 26}, {16, 27, 0, 5, 9}, {16, 28, 1, 5, 
              11}, {16, 29, 6, 5, 13}, {10, 0, 30}, {13, 30, 9, 11, 13, 30}, {
              16, 31, 0, 5, 11}, {16, 32, 1, 5, 9}, {16, 33, 6, 5, 13}, {10, 
              0, 34}, {13, 34, 11, 9, 13, 34}, {16, 35, 0, 5, 9}, {16, 36, 1, 
              5, 11}, {16, 37, 6, 5, 13}, {10, 0, 38}, {13, 38, 9, 11, 13, 
              38}, {16, 39, 0, 5, 11}, {16, 40, 1, 5, 9}, {16, 41, 6, 5, 
              13}, {10, 0, 42}, {13, 42, 11, 9, 13, 42}, {16, 43, 0, 5, 9}, {
              16, 44, 1, 5, 11}, {16, 45, 6, 5, 13}, {10, 0, 46}, {13, 46, 9, 
              11, 13, 46}, {16, 47, 0, 5, 11}, {16, 48, 1, 5, 9}, {16, 49, 6, 
              5, 13}, {10, 0, 50}, {13, 50, 11, 9, 13, 50}, {16, 51, 0, 5, 
              9}, {16, 52, 1, 5, 11}, {16, 53, 6, 5, 13}, {10, 0, 54}, {13, 
              54, 9, 11, 13, 54}, {16, 55, 0, 5, 11}, {16, 56, 1, 5, 9}, {16, 
              57, 6, 5, 13}, {10, 0, 58}, {13, 58, 11, 9, 13, 58}, {16, 7, 14,
               18, 22, 26, 30, 34, 38, 42, 46, 50, 54, 58, 11}, {1}}, 
             Function[{$CellContext`phi1$22393, $CellContext`phi2$22393}, 
              
              Block[{Compile`$96, Compile`$103, Compile`$108, Compile`$110, 
                Compile`$147}, Compile`$96 = $CellContext`phi1$22393^2; 
               Compile`$103 = $CellContext`phi2$22393^2; 
               Compile`$108 = 1 + Compile`$96 + Compile`$103; 
               Compile`$110 = Compile`$108^(-1); 
               Compile`$147 = -1 + Compile`$96 + Compile`$103; 
               30.616258720840584` (1 - 
                 0.6719355684716646 $CellContext`phi1$22393 Compile`$110 - 
                 0.3795266557666026 $CellContext`phi2$22393 Compile`$110 - 
                 0.9225592270127267 Compile`$147 Compile`$110) (1 + 
                 1.21558798367449 $CellContext`phi1$22393 Compile`$110 + 
                 0.4857136330295899 $CellContext`phi2$22393 Compile`$110 - 
                 0.7560469761587875 Compile`$147 Compile`$110) (1 - 
                 0.5310663415434267 $CellContext`phi1$22393 Compile`$110 + 
                 1.5825285657816939` $CellContext`phi2$22393 Compile`$110 - 
                 0.5508111925525423 Compile`$147 Compile`$110) (1 + 
                 0.8670620122047428 $CellContext`phi1$22393 Compile`$110 - 
                 1.5169014046705556` $CellContext`phi2$22393 Compile`$110 - 
                 0.48662449473386504` Compile`$147 Compile`$110) (1 - 
                 1.959084052594859 $CellContext`phi1$22393 Compile`$110 + 
                 0.257782435890779 $CellContext`phi2$22393 Compile`$110 - 
                 0.15454602115124183` Compile`$147 Compile`$110) (1 - 
                 1.0949932093435537` $CellContext`phi1$22393 Compile`$110 - 
                 1.6577706315987935` $CellContext`phi2$22393 Compile`$110 - 
                 0.11487646027368059` Compile`$147 Compile`$110) (1 + 
                 1.0949932093435537` $CellContext`phi1$22393 Compile`$110 + 
                 1.6577706315987935` $CellContext`phi2$22393 Compile`$110 + 
                 0.11487646027368059` Compile`$147 Compile`$110) (1 + 
                 1.959084052594859 $CellContext`phi1$22393 Compile`$110 - 
                 0.257782435890779 $CellContext`phi2$22393 Compile`$110 + 
                 0.15454602115124183` Compile`$147 Compile`$110) (1 - 
                 0.8670620122047428 $CellContext`phi1$22393 Compile`$110 + 
                 1.5169014046705556` $CellContext`phi2$22393 Compile`$110 + 
                 0.48662449473386504` Compile`$147 Compile`$110) (1 + 
                 0.5310663415434267 $CellContext`phi1$22393 Compile`$110 - 
                 1.5825285657816939` $CellContext`phi2$22393 Compile`$110 + 
                 0.5508111925525423 Compile`$147 Compile`$110) (1 - 
                 1.21558798367449 $CellContext`phi1$22393 Compile`$110 - 
                 0.4857136330295899 $CellContext`phi2$22393 Compile`$110 + 
                 0.7560469761587875 Compile`$147 Compile`$110) (1 + 
                 0.6719355684716646 $CellContext`phi1$22393 Compile`$110 + 
                 0.3795266557666026 $CellContext`phi2$22393 Compile`$110 + 
                 0.9225592270127267 Compile`$147 Compile`$110)]], Evaluate], 
           "V" -> CompiledFunction[{11, 14.3, 5470}, {
              Blank[Real], 
              Blank[Real]}, {{3, 0, 0}, {3, 0, 1}, {3, 0, 
             109}}, {{-1.6577706315987935`, {3, 0, 45}}, {
              2.43117596734898, {3, 0, 105}}, {
              0.6719355684716646, {3, 0, 85}}, {-0.9714272660591798, {3, 0, 
               143}}, {1.7341240244094855`, {3, 0, 115}}, {
              0.9225592270127267, {3, 0, 89}}, {-1.21558798367449, {3, 0, 
               78}}, {-0.9225592270127267, {3, 0, 12}}, {
              3.033802809341111, {3, 0, 137}}, {
              0.4857136330295899, {3, 0, 17}}, {-0.515564871781558, {3, 0, 
               134}}, {-1.3438711369433292`, {3, 0, 96}}, {
              3.315541263197587, {3, 0, 131}}, {
              0.257782435890779, {3, 0, 38}}, {
              0.22975292054736118`, {3, 0, 128}}, {-3.033802809341111, {3, 0, 
               116}}, {-0.257782435890779, {3, 0, 
               59}}, {-0.22975292054736118`, {3, 0, 127}}, {
              1.3438711369433292`, {3, 0, 145}}, {
              0.7560469761587875, {3, 0, 82}}, {-0.11487646027368059`, {3, 0, 
               47}}, {-1.5825285657816939`, {3, 0, 73}}, {
              0.15454602115124183`, {3, 0, 61}}, {
              1.5169014046705556`, {3, 0, 66}}, {-3.315541263197587, {3, 0, 
               126}}, {-0.6719355684716646, {3, 0, 8}}, {
              0.515564871781558, {3, 0, 121}}, {-0.3795266557666026, {3, 0, 
               10}}, {3.918168105189718, {3, 0, 133}}, {
              0.9714272660591798, {3, 0, 106}}, {
              1.21558798367449, {3, 0, 15}}, {-0.48662449473386504`, {3, 0, 
               33}}, {0.11487646027368059`, {3, 0, 
               54}}, {-0.7590533115332052, {3, 0, 98}}, {
              0.9732489894677301, {3, 0, 118}}, {
              1.0621326830868534`, {3, 0, 139}}, {-1.8451184540254535`, {3, 0,
                100}}, {-0.5508111925525423, {3, 0, 26}}, {
              0.48662449473386504`, {3, 0, 68}}, {-3.918168105189718, {3, 0, 
               120}}, {-1.5169014046705556`, {3, 0, 31}}, {
              1.1016223851050846`, {3, 0, 113}}, {
              2.1899864186871074`, {3, 0, 130}}, {
              1.512093952317575, {3, 0, 108}}, {
              0.30909204230248366`, {3, 0, 123}}, {
              0.7590533115332052, {3, 0, 146}}, {
              1.5825285657816939`, {3, 0, 24}}, {-1, {2, 0, 1}}, {
              1.0949932093435537`, {3, 0, 50}}, {
              1.6577706315987935`, {3, 0, 52}}, {
              1.8451184540254535`, {3, 0, 103}}, {
              30.616258720840584`, {3, 0, 95}}, {-1.0621326830868534`, {3, 0, 
               110}}, {
              3.1650571315633877`, {3, 0, 111}}, {-2.1899864186871074`, {3, 0,
                125}}, {-1.1016223851050846`, {3, 0, 
               112}}, {-1.512093952317575, {3, 0, 
               107}}, {-0.9732489894677301, {3, 0, 117}}, {
              1, {2, 0, 0}}, {-0.15454602115124183`, {3, 0, 
               40}}, {-0.4857136330295899, {3, 0, 
               80}}, {-0.30909204230248366`, {3, 0, 
               122}}, {-0.5310663415434267, {3, 0, 22}}, {
              1.959084052594859, {3, 0, 57}}, {-1.0949932093435537`, {3, 0, 
               43}}, {0.25, {3, 0, 94}}, {
              0.5310663415434267, {3, 0, 71}}, {-1.959084052594859, {3, 0, 
               36}}, {-1.7341240244094855`, {3, 0, 
               136}}, {-2.43117596734898, {3, 0, 
               142}}, {-3.1650571315633877`, {3, 0, 140}}, {
              0.5508111925525423, {3, 0, 75}}, {-0.7560469761587875, {3, 0, 
               19}}, {0.8670620122047428, {3, 0, 29}}, {0.5, {3, 0, 93}}, {
              0.3795266557666026, {3, 0, 87}}, {-0.8670620122047428, {3, 0, 
               64}}}, {0, 2, 149, 0, 
             0}, {{40, 56, 3, 0, 0, 3, 0, 2}, {40, 56, 3, 0, 1, 3, 0, 3}, {10,
               0, 4}, {13, 4, 2, 3, 4}, {40, 56, 3, 0, 4, 3, 0, 5}, {40, 60, 
              3, 0, 5, 3, 0, 6}, {40, 60, 3, 0, 4, 3, 0, 5}, {10, 1, 7}, {13, 
              7, 2, 3, 7}, {16, 8, 0, 5, 9}, {16, 10, 1, 5, 11}, {16, 12, 7, 
              5, 13}, {10, 0, 14}, {13, 14, 9, 11, 13, 14}, {16, 15, 0, 5, 
              16}, {16, 17, 1, 5, 18}, {16, 19, 7, 5, 20}, {10, 0, 21}, {13, 
              21, 16, 18, 20, 21}, {16, 22, 0, 5, 23}, {16, 24, 1, 5, 25}, {
              16, 26, 7, 5, 27}, {10, 0, 28}, {13, 28, 23, 25, 27, 28}, {16, 
              29, 0, 5, 30}, {16, 31, 1, 5, 32}, {16, 33, 7, 5, 34}, {10, 0, 
              35}, {13, 35, 30, 32, 34, 35}, {16, 36, 0, 5, 37}, {16, 38, 1, 
              5, 39}, {16, 40, 7, 5, 41}, {10, 0, 42}, {13, 42, 37, 39, 41, 
              42}, {16, 43, 0, 5, 44}, {16, 45, 1, 5, 46}, {16, 47, 7, 5, 
              48}, {10, 0, 49}, {13, 49, 44, 46, 48, 49}, {16, 50, 0, 5, 
              51}, {16, 52, 1, 5, 53}, {16, 54, 7, 5, 55}, {10, 0, 56}, {13, 
              56, 51, 53, 55, 56}, {16, 57, 0, 5, 58}, {16, 59, 1, 5, 60}, {
              16, 61, 7, 5, 62}, {10, 0, 63}, {13, 63, 58, 60, 62, 63}, {16, 
              64, 0, 5, 65}, {16, 66, 1, 5, 67}, {16, 68, 7, 5, 69}, {10, 0, 
              70}, {13, 70, 65, 67, 69, 70}, {16, 71, 0, 5, 72}, {16, 73, 1, 
              5, 74}, {16, 75, 7, 5, 76}, {10, 0, 77}, {13, 77, 72, 74, 76, 
              77}, {16, 78, 0, 5, 79}, {16, 80, 1, 5, 81}, {16, 82, 7, 5, 
              83}, {10, 0, 84}, {13, 84, 79, 81, 83, 84}, {16, 85, 0, 5, 
              86}, {16, 87, 1, 5, 88}, {16, 89, 7, 5, 90}, {10, 0, 91}, {13, 
              91, 86, 88, 90, 91}, {40, 56, 3, 0, 4, 3, 0, 92}, {16, 96, 2, 6,
               97}, {16, 98, 0, 1, 6, 99}, {16, 100, 0, 7, 6, 101}, {16, 85, 
              5, 102}, {16, 103, 0, 5, 104}, {13, 97, 99, 101, 102, 104, 
              97}, {16, 95, 97, 14, 21, 28, 35, 42, 49, 56, 63, 70, 77, 84, 
              99}, {16, 105, 2, 6, 101}, {16, 106, 0, 1, 6, 104}, {16, 107, 0,
               7, 6, 102}, {16, 78, 5, 97}, {16, 108, 0, 5, 109}, {13, 101, 
              104, 102, 97, 109, 101}, {16, 95, 101, 14, 21, 28, 35, 42, 49, 
              56, 63, 70, 77, 91, 104}, {16, 110, 2, 6, 102}, {16, 111, 0, 1, 
              6, 109}, {16, 112, 0, 7, 6, 97}, {16, 71, 5, 101}, {16, 113, 0, 
              5, 114}, {13, 102, 109, 97, 101, 114, 102}, {16, 95, 102, 14, 
              21, 28, 35, 42, 49, 56, 63, 70, 84, 91, 109}, {16, 115, 2, 6, 
              97}, {16, 116, 0, 1, 6, 114}, {16, 117, 0, 7, 6, 101}, {16, 64, 
              5, 102}, {16, 118, 0, 5, 119}, {13, 97, 114, 101, 102, 119, 
              97}, {16, 95, 97, 14, 21, 28, 35, 42, 49, 56, 63, 77, 84, 91, 
              114}, {16, 120, 2, 6, 101}, {16, 121, 0, 1, 6, 119}, {16, 122, 
              0, 7, 6, 102}, {16, 57, 5, 97}, {16, 123, 0, 5, 124}, {13, 101, 
              119, 102, 97, 124, 101}, {16, 95, 101, 14, 21, 28, 35, 42, 49, 
              56, 70, 77, 84, 91, 119}, {16, 125, 2, 6, 102}, {16, 126, 0, 1, 
              6, 124}, {16, 127, 0, 7, 6, 97}, {16, 50, 5, 101}, {16, 128, 0, 
              5, 129}, {13, 102, 124, 97, 101, 129, 102}, {16, 95, 102, 14, 
              21, 28, 35, 42, 49, 63, 70, 77, 84, 91, 124}, {16, 130, 2, 6, 
              97}, {16, 131, 0, 1, 6, 129}, {16, 128, 0, 7, 6, 102}, {16, 43, 
              5, 101}, {16, 127, 0, 5, 132}, {13, 97, 129, 102, 101, 132, 
              97}, {16, 95, 97, 14, 21, 28, 35, 42, 56, 63, 70, 77, 84, 91, 
              129}, {16, 133, 2, 6, 102}, {16, 134, 0, 1, 6, 132}, {16, 123, 
              0, 7, 6, 97}, {16, 36, 5, 101}, {16, 122, 0, 5, 135}, {13, 102, 
              132, 97, 101, 135, 102}, {16, 95, 102, 14, 21, 28, 35, 49, 56, 
              63, 70, 77, 84, 91, 132}, {16, 136, 2, 6, 97}, {16, 137, 0, 1, 
              6, 135}, {16, 118, 0, 7, 6, 102}, {16, 29, 5, 101}, {16, 117, 0,
               5, 138}, {13, 97, 135, 102, 101, 138, 97}, {16, 95, 97, 14, 21,
               28, 42, 49, 56, 63, 70, 77, 84, 91, 135}, {16, 139, 2, 6, 
              102}, {16, 140, 0, 1, 6, 138}, {16, 113, 0, 7, 6, 97}, {16, 22, 
              5, 101}, {16, 112, 0, 5, 141}, {13, 102, 138, 97, 101, 141, 
              102}, {16, 95, 102, 14, 21, 35, 42, 49, 56, 63, 70, 77, 84, 91, 
              138}, {16, 142, 2, 6, 97}, {16, 143, 0, 1, 6, 141}, {16, 108, 0,
               7, 6, 102}, {16, 15, 5, 101}, {16, 107, 0, 5, 144}, {13, 97, 
              141, 102, 101, 144, 97}, {16, 95, 97, 14, 28, 35, 42, 49, 56, 
              63, 70, 77, 84, 91, 141}, {16, 145, 2, 6, 102}, {16, 146, 0, 1, 
              6, 144}, {16, 103, 0, 7, 6, 97}, {16, 8, 5, 101}, {16, 100, 0, 
              5, 147}, {13, 102, 144, 97, 101, 147, 102}, {16, 95, 102, 21, 
              28, 35, 42, 49, 56, 63, 70, 77, 84, 91, 144}, {13, 99, 104, 109,
               114, 119, 124, 129, 132, 135, 138, 141, 144, 99}, {40, 56, 3, 
              0, 99, 3, 0, 104}, {16, 94, 92, 104, 99}, {16, 96, 0, 1, 6, 
              104}, {16, 98, 3, 6, 109}, {16, 100, 1, 7, 6, 114}, {16, 87, 5, 
              119}, {16, 103, 1, 5, 124}, {13, 104, 109, 114, 119, 124, 
              104}, {16, 95, 104, 14, 21, 28, 35, 42, 49, 56, 63, 70, 77, 84, 
              109}, {16, 105, 0, 1, 6, 104}, {16, 106, 3, 6, 114}, {16, 107, 
              1, 7, 6, 119}, {16, 80, 5, 124}, {16, 108, 1, 5, 129}, {13, 104,
               114, 119, 124, 129, 104}, {16, 95, 104, 14, 21, 28, 35, 42, 49,
               56, 63, 70, 77, 91, 114}, {16, 110, 0, 1, 6, 104}, {16, 111, 3,
               6, 119}, {16, 112, 1, 7, 6, 124}, {16, 73, 5, 129}, {16, 113, 
              1, 5, 132}, {13, 104, 119, 124, 129, 132, 104}, {16, 95, 104, 
              14, 21, 28, 35, 42, 49, 56, 63, 70, 84, 91, 119}, {16, 115, 0, 
              1, 6, 104}, {16, 116, 3, 6, 124}, {16, 117, 1, 7, 6, 129}, {16, 
              66, 5, 132}, {16, 118, 1, 5, 135}, {13, 104, 124, 129, 132, 135,
               104}, {16, 95, 104, 14, 21, 28, 35, 42, 49, 56, 63, 77, 84, 91,
               124}, {16, 120, 0, 1, 6, 104}, {16, 121, 3, 6, 129}, {16, 122, 
              1, 7, 6, 132}, {16, 59, 5, 135}, {16, 123, 1, 5, 138}, {13, 104,
               129, 132, 135, 138, 104}, {16, 95, 104, 14, 21, 28, 35, 42, 49,
               56, 70, 77, 84, 91, 129}, {16, 125, 0, 1, 6, 104}, {16, 126, 3,
               6, 132}, {16, 127, 1, 7, 6, 135}, {16, 52, 5, 138}, {16, 128, 
              1, 5, 141}, {13, 104, 132, 135, 138, 141, 104}, {16, 95, 104, 
              14, 21, 28, 35, 42, 49, 63, 70, 77, 84, 91, 132}, {16, 130, 0, 
              1, 6, 104}, {16, 131, 3, 6, 135}, {16, 128, 1, 7, 6, 138}, {16, 
              45, 5, 141}, {16, 127, 1, 5, 144}, {13, 104, 135, 138, 141, 144,
               104}, {16, 95, 104, 14, 21, 28, 35, 42, 56, 63, 70, 77, 84, 91,
               135}, {16, 133, 0, 1, 6, 104}, {16, 134, 3, 6, 138}, {16, 123, 
              1, 7, 6, 141}, {16, 38, 5, 144}, {16, 122, 1, 5, 102}, {13, 104,
               138, 141, 144, 102, 104}, {16, 95, 104, 14, 21, 28, 35, 49, 56,
               63, 70, 77, 84, 91, 138}, {16, 136, 0, 1, 6, 104}, {16, 137, 3,
               6, 141}, {16, 118, 1, 7, 6, 144}, {16, 31, 5, 102}, {16, 117, 
              1, 5, 97}, {13, 104, 141, 144, 102, 97, 104}, {16, 95, 104, 14, 
              21, 28, 42, 49, 56, 63, 70, 77, 84, 91, 141}, {16, 139, 0, 1, 6,
               104}, {16, 140, 3, 6, 144}, {16, 113, 1, 7, 6, 102}, {16, 24, 
              5, 97}, {16, 112, 1, 5, 101}, {13, 104, 144, 102, 97, 101, 
              104}, {16, 95, 104, 14, 21, 35, 42, 49, 56, 63, 70, 77, 84, 91, 
              144}, {16, 142, 0, 1, 6, 104}, {16, 143, 3, 6, 102}, {16, 108, 
              1, 7, 6, 97}, {16, 17, 5, 101}, {16, 107, 1, 5, 147}, {13, 104, 
              102, 97, 101, 147, 104}, {16, 95, 104, 14, 28, 35, 42, 49, 56, 
              63, 70, 77, 84, 91, 102}, {16, 145, 0, 1, 6, 104}, {16, 146, 3, 
              6, 97}, {16, 103, 1, 7, 6, 101}, {16, 10, 5, 147}, {16, 100, 1, 
              5, 148}, {13, 104, 97, 101, 147, 148, 104}, {16, 95, 104, 21, 
              28, 35, 42, 49, 56, 63, 70, 77, 84, 91, 97}, {13, 109, 114, 119,
               124, 129, 132, 135, 138, 141, 144, 102, 97, 109}, {40, 56, 3, 
              0, 109, 3, 0, 114}, {16, 94, 92, 114, 109}, {13, 99, 109, 99}, {
              16, 93, 99, 109}, {1}}, 
             Function[{$CellContext`phi1$22393, $CellContext`phi2$22393}, 
              
              Block[{Compile`$141, Compile`$142, Compile`$152, Compile`$156, 
                Compile`$163, Compile`$160, Compile`$169, Compile`$170, 
                Compile`$171, Compile`$172, Compile`$177, Compile`$178, 
                Compile`$179, Compile`$182, Compile`$183, Compile`$184, 
                Compile`$186, Compile`$188, Compile`$189, Compile`$190, 
                Compile`$191, Compile`$192, Compile`$193, Compile`$194, 
                Compile`$195, Compile`$196, Compile`$197, Compile`$198, 
                Compile`$199, Compile`$202, Compile`$203, Compile`$204, 
                Compile`$205, Compile`$206, Compile`$207, Compile`$208, 
                Compile`$209, Compile`$210, Compile`$211, Compile`$212, 
                Compile`$213, Compile`$214, Compile`$215, Compile`$216, 
                Compile`$217, Compile`$220, Compile`$221, Compile`$222, 
                Compile`$223, Compile`$224, Compile`$232, Compile`$233, 
                Compile`$234, Compile`$235, Compile`$153}, 
               Compile`$141 = $CellContext`phi1$22393^2; 
               Compile`$142 = $CellContext`phi2$22393^2; 
               Compile`$152 = 1 + Compile`$141 + Compile`$142; 
               Compile`$156 = Compile`$152^(-2); 
               Compile`$163 = Compile`$152^(-1); 
               Compile`$160 = -1 + Compile`$141 + Compile`$142; 
               Compile`$169 = (-0.6719355684716646) $CellContext`phi1$22393 
                 Compile`$163; 
               Compile`$170 = (-0.3795266557666026) $CellContext`phi2$22393 
                 Compile`$163; 
               Compile`$171 = (-0.9225592270127267) Compile`$160 Compile`$163; 
               Compile`$172 = 1 + Compile`$169 + Compile`$170 + Compile`$171; 
               Compile`$177 = 
                1.21558798367449 $CellContext`phi1$22393 Compile`$163; 
               Compile`$178 = 
                0.4857136330295899 $CellContext`phi2$22393 Compile`$163; 
               Compile`$179 = (-0.7560469761587875) Compile`$160 Compile`$163; 
               Compile`$182 = 1 + Compile`$177 + Compile`$178 + Compile`$179; 
               Compile`$183 = (-0.5310663415434267) $CellContext`phi1$22393 
                 Compile`$163; 
               Compile`$184 = 
                1.5825285657816939` $CellContext`phi2$22393 Compile`$163; 
               Compile`$186 = (-0.5508111925525423) Compile`$160 Compile`$163; 
               Compile`$188 = 1 + Compile`$183 + Compile`$184 + Compile`$186; 
               Compile`$189 = 
                0.8670620122047428 $CellContext`phi1$22393 Compile`$163; 
               Compile`$190 = (-1.5169014046705556`) $CellContext`phi2$22393 
                 Compile`$163; 
               Compile`$191 = (-0.48662449473386504`) Compile`$160 
                 Compile`$163; 
               Compile`$192 = 1 + Compile`$189 + Compile`$190 + Compile`$191; 
               Compile`$193 = (-1.959084052594859) $CellContext`phi1$22393 
                 Compile`$163; 
               Compile`$194 = 
                0.257782435890779 $CellContext`phi2$22393 Compile`$163; 
               Compile`$195 = (-0.15454602115124183`) Compile`$160 
                 Compile`$163; 
               Compile`$196 = 1 + Compile`$193 + Compile`$194 + Compile`$195; 
               Compile`$197 = (-1.0949932093435537`) $CellContext`phi1$22393 
                 Compile`$163; 
               Compile`$198 = (-1.6577706315987935`) $CellContext`phi2$22393 
                 Compile`$163; 
               Compile`$199 = (-0.11487646027368059`) Compile`$160 
                 Compile`$163; 
               Compile`$202 = 1 + Compile`$197 + Compile`$198 + Compile`$199; 
               Compile`$203 = 
                1.0949932093435537` $CellContext`phi1$22393 Compile`$163; 
               Compile`$204 = 
                1.6577706315987935` $CellContext`phi2$22393 Compile`$163; 
               Compile`$205 = 0.11487646027368059` Compile`$160 Compile`$163; 
               Compile`$206 = 1 + Compile`$203 + Compile`$204 + Compile`$205; 
               Compile`$207 = 
                1.959084052594859 $CellContext`phi1$22393 Compile`$163; 
               Compile`$208 = (-0.257782435890779) $CellContext`phi2$22393 
                 Compile`$163; 
               Compile`$209 = 0.15454602115124183` Compile`$160 Compile`$163; 
               Compile`$210 = 1 + Compile`$207 + Compile`$208 + Compile`$209; 
               Compile`$211 = (-0.8670620122047428) $CellContext`phi1$22393 
                 Compile`$163; 
               Compile`$212 = 
                1.5169014046705556` $CellContext`phi2$22393 Compile`$163; 
               Compile`$213 = 0.48662449473386504` Compile`$160 Compile`$163; 
               Compile`$214 = 1 + Compile`$211 + Compile`$212 + Compile`$213; 
               Compile`$215 = 
                0.5310663415434267 $CellContext`phi1$22393 Compile`$163; 
               Compile`$216 = (-1.5825285657816939`) $CellContext`phi2$22393 
                 Compile`$163; 
               Compile`$217 = 0.5508111925525423 Compile`$160 Compile`$163; 
               Compile`$220 = 1 + Compile`$215 + Compile`$216 + Compile`$217; 
               Compile`$221 = (-1.21558798367449) $CellContext`phi1$22393 
                 Compile`$163; 
               Compile`$222 = (-0.4857136330295899) $CellContext`phi2$22393 
                 Compile`$163; 
               Compile`$223 = 0.7560469761587875 Compile`$160 Compile`$163; 
               Compile`$224 = 1 + Compile`$221 + Compile`$222 + Compile`$223; 
               Compile`$232 = 
                0.6719355684716646 $CellContext`phi1$22393 Compile`$163; 
               Compile`$233 = 
                0.3795266557666026 $CellContext`phi2$22393 Compile`$163; 
               Compile`$234 = 0.9225592270127267 Compile`$160 Compile`$163; 
               Compile`$235 = 1 + Compile`$232 + Compile`$233 + Compile`$234; 
               Compile`$153 = 
                Compile`$152^2; ((
                   Compile`$153 (
                    30.616258720840584` ((-1.3438711369433292`) Compile`$141 
                    Compile`$156 - 
                    0.7590533115332052 $CellContext`phi1$22393 \
$CellContext`phi2$22393 Compile`$156 - 
                    1.8451184540254535` $CellContext`phi1$22393 Compile`$160 
                    Compile`$156 + 0.6719355684716646 Compile`$163 + 
                    1.8451184540254535` $CellContext`phi1$22393 Compile`$163) 
                    Compile`$172 Compile`$182 Compile`$188 Compile`$192 
                    Compile`$196 Compile`$202 Compile`$206 Compile`$210 
                    Compile`$214 Compile`$220 Compile`$224 + 
                    30.616258720840584` (
                    2.43117596734898 Compile`$141 Compile`$156 + 
                    0.9714272660591798 $CellContext`phi1$22393 \
$CellContext`phi2$22393 Compile`$156 - 
                    1.512093952317575 $CellContext`phi1$22393 Compile`$160 
                    Compile`$156 - 1.21558798367449 Compile`$163 + 
                    1.512093952317575 $CellContext`phi1$22393 Compile`$163) 
                    Compile`$172 Compile`$182 Compile`$188 Compile`$192 
                    Compile`$196 Compile`$202 Compile`$206 Compile`$210 
                    Compile`$214 Compile`$220 Compile`$235 + 
                    30.616258720840584` ((-1.0621326830868534`) Compile`$141 
                    Compile`$156 + 
                    3.1650571315633877` $CellContext`phi1$22393 \
$CellContext`phi2$22393 Compile`$156 - 
                    1.1016223851050846` $CellContext`phi1$22393 Compile`$160 
                    Compile`$156 + 0.5310663415434267 Compile`$163 + 
                    1.1016223851050846` $CellContext`phi1$22393 Compile`$163) 
                    Compile`$172 Compile`$182 Compile`$188 Compile`$192 
                    Compile`$196 Compile`$202 Compile`$206 Compile`$210 
                    Compile`$214 Compile`$224 Compile`$235 + 
                    30.616258720840584` (
                    1.7341240244094855` Compile`$141 Compile`$156 - 
                    3.033802809341111 $CellContext`phi1$22393 \
$CellContext`phi2$22393 Compile`$156 - 
                    0.9732489894677301 $CellContext`phi1$22393 Compile`$160 
                    Compile`$156 - 0.8670620122047428 Compile`$163 + 
                    0.9732489894677301 $CellContext`phi1$22393 Compile`$163) 
                    Compile`$172 Compile`$182 Compile`$188 Compile`$192 
                    Compile`$196 Compile`$202 Compile`$206 Compile`$210 
                    Compile`$220 Compile`$224 Compile`$235 + 
                    30.616258720840584` ((-3.918168105189718) Compile`$141 
                    Compile`$156 + 
                    0.515564871781558 $CellContext`phi1$22393 \
$CellContext`phi2$22393 Compile`$156 - 
                    0.30909204230248366` $CellContext`phi1$22393 Compile`$160 
                    Compile`$156 + 1.959084052594859 Compile`$163 + 
                    0.30909204230248366` $CellContext`phi1$22393 Compile`$163)
                     Compile`$172 Compile`$182 Compile`$188 Compile`$192 
                    Compile`$196 Compile`$202 Compile`$206 Compile`$214 
                    Compile`$220 Compile`$224 Compile`$235 + 
                    30.616258720840584` ((-2.1899864186871074`) Compile`$141 
                    Compile`$156 - 
                    3.315541263197587 $CellContext`phi1$22393 \
$CellContext`phi2$22393 Compile`$156 - 
                    0.22975292054736118` $CellContext`phi1$22393 Compile`$160 
                    Compile`$156 + 1.0949932093435537` Compile`$163 + 
                    0.22975292054736118` $CellContext`phi1$22393 Compile`$163)
                     Compile`$172 Compile`$182 Compile`$188 Compile`$192 
                    Compile`$196 Compile`$202 Compile`$210 Compile`$214 
                    Compile`$220 Compile`$224 Compile`$235 + 
                    30.616258720840584` (
                    2.1899864186871074` Compile`$141 Compile`$156 + 
                    3.315541263197587 $CellContext`phi1$22393 \
$CellContext`phi2$22393 Compile`$156 + 
                    0.22975292054736118` $CellContext`phi1$22393 Compile`$160 
                    Compile`$156 - 1.0949932093435537` Compile`$163 - 
                    0.22975292054736118` $CellContext`phi1$22393 Compile`$163)
                     Compile`$172 Compile`$182 Compile`$188 Compile`$192 
                    Compile`$196 Compile`$206 Compile`$210 Compile`$214 
                    Compile`$220 Compile`$224 Compile`$235 + 
                    30.616258720840584` (
                    3.918168105189718 Compile`$141 Compile`$156 - 
                    0.515564871781558 $CellContext`phi1$22393 \
$CellContext`phi2$22393 Compile`$156 + 
                    0.30909204230248366` $CellContext`phi1$22393 Compile`$160 
                    Compile`$156 - 1.959084052594859 Compile`$163 - 
                    0.30909204230248366` $CellContext`phi1$22393 Compile`$163)
                     Compile`$172 Compile`$182 Compile`$188 Compile`$192 
                    Compile`$202 Compile`$206 Compile`$210 Compile`$214 
                    Compile`$220 Compile`$224 Compile`$235 + 
                    30.616258720840584` ((-1.7341240244094855`) Compile`$141 
                    Compile`$156 + 
                    3.033802809341111 $CellContext`phi1$22393 \
$CellContext`phi2$22393 Compile`$156 + 
                    0.9732489894677301 $CellContext`phi1$22393 Compile`$160 
                    Compile`$156 + 0.8670620122047428 Compile`$163 - 
                    0.9732489894677301 $CellContext`phi1$22393 Compile`$163) 
                    Compile`$172 Compile`$182 Compile`$188 Compile`$196 
                    Compile`$202 Compile`$206 Compile`$210 Compile`$214 
                    Compile`$220 Compile`$224 Compile`$235 + 
                    30.616258720840584` (
                    1.0621326830868534` Compile`$141 Compile`$156 - 
                    3.1650571315633877` $CellContext`phi1$22393 \
$CellContext`phi2$22393 Compile`$156 + 
                    1.1016223851050846` $CellContext`phi1$22393 Compile`$160 
                    Compile`$156 - 0.5310663415434267 Compile`$163 - 
                    1.1016223851050846` $CellContext`phi1$22393 Compile`$163) 
                    Compile`$172 Compile`$182 Compile`$192 Compile`$196 
                    Compile`$202 Compile`$206 Compile`$210 Compile`$214 
                    Compile`$220 Compile`$224 Compile`$235 + 
                    30.616258720840584` ((-2.43117596734898) Compile`$141 
                    Compile`$156 - 
                    0.9714272660591798 $CellContext`phi1$22393 \
$CellContext`phi2$22393 Compile`$156 + 
                    1.512093952317575 $CellContext`phi1$22393 Compile`$160 
                    Compile`$156 + 1.21558798367449 Compile`$163 - 
                    1.512093952317575 $CellContext`phi1$22393 Compile`$163) 
                    Compile`$172 Compile`$188 Compile`$192 Compile`$196 
                    Compile`$202 Compile`$206 Compile`$210 Compile`$214 
                    Compile`$220 Compile`$224 Compile`$235 + 
                    30.616258720840584` (
                    1.3438711369433292` Compile`$141 Compile`$156 + 
                    0.7590533115332052 $CellContext`phi1$22393 \
$CellContext`phi2$22393 Compile`$156 + 
                    1.8451184540254535` $CellContext`phi1$22393 Compile`$160 
                    Compile`$156 - 0.6719355684716646 Compile`$163 - 
                    1.8451184540254535` $CellContext`phi1$22393 Compile`$163) 
                    Compile`$182 Compile`$188 Compile`$192 Compile`$196 
                    Compile`$202 Compile`$206 Compile`$210 Compile`$214 
                    Compile`$220 Compile`$224 Compile`$235)^2)/
                  4 + (Compile`$153 (
                    30.616258720840584` ((-1.3438711369433292`) \
$CellContext`phi1$22393 $CellContext`phi2$22393 Compile`$156 - 
                    0.7590533115332052 Compile`$142 Compile`$156 - 
                    1.8451184540254535` $CellContext`phi2$22393 Compile`$160 
                    Compile`$156 + 0.3795266557666026 Compile`$163 + 
                    1.8451184540254535` $CellContext`phi2$22393 Compile`$163) 
                    Compile`$172 Compile`$182 Compile`$188 Compile`$192 
                    Compile`$196 Compile`$202 Compile`$206 Compile`$210 
                    Compile`$214 Compile`$220 Compile`$224 + 
                    30.616258720840584` (
                    2.43117596734898 $CellContext`phi1$22393 \
$CellContext`phi2$22393 Compile`$156 + 
                    0.9714272660591798 Compile`$142 Compile`$156 - 
                    1.512093952317575 $CellContext`phi2$22393 Compile`$160 
                    Compile`$156 - 0.4857136330295899 Compile`$163 + 
                    1.512093952317575 $CellContext`phi2$22393 Compile`$163) 
                    Compile`$172 Compile`$182 Compile`$188 Compile`$192 
                    Compile`$196 Compile`$202 Compile`$206 Compile`$210 
                    Compile`$214 Compile`$220 Compile`$235 + 
                    30.616258720840584` ((-1.0621326830868534`) \
$CellContext`phi1$22393 $CellContext`phi2$22393 Compile`$156 + 
                    3.1650571315633877` Compile`$142 Compile`$156 - 
                    1.1016223851050846` $CellContext`phi2$22393 Compile`$160 
                    Compile`$156 - 1.5825285657816939` Compile`$163 + 
                    1.1016223851050846` $CellContext`phi2$22393 Compile`$163) 
                    Compile`$172 Compile`$182 Compile`$188 Compile`$192 
                    Compile`$196 Compile`$202 Compile`$206 Compile`$210 
                    Compile`$214 Compile`$224 Compile`$235 + 
                    30.616258720840584` (
                    1.7341240244094855` $CellContext`phi1$22393 \
$CellContext`phi2$22393 Compile`$156 - 3.033802809341111 Compile`$142 
                    Compile`$156 - 0.9732489894677301 $CellContext`phi2$22393 
                    Compile`$160 Compile`$156 + 
                    1.5169014046705556` Compile`$163 + 
                    0.9732489894677301 $CellContext`phi2$22393 Compile`$163) 
                    Compile`$172 Compile`$182 Compile`$188 Compile`$192 
                    Compile`$196 Compile`$202 Compile`$206 Compile`$210 
                    Compile`$220 Compile`$224 Compile`$235 + 
                    30.616258720840584` ((-3.918168105189718) \
$CellContext`phi1$22393 $CellContext`phi2$22393 Compile`$156 + 
                    0.515564871781558 Compile`$142 Compile`$156 - 
                    0.30909204230248366` $CellContext`phi2$22393 Compile`$160 
                    Compile`$156 - 0.257782435890779 Compile`$163 + 
                    0.30909204230248366` $CellContext`phi2$22393 Compile`$163)
                     Compile`$172 Compile`$182 Compile`$188 Compile`$192 
                    Compile`$196 Compile`$202 Compile`$206 Compile`$214 
                    Compile`$220 Compile`$224 Compile`$235 + 
                    30.616258720840584` ((-2.1899864186871074`) \
$CellContext`phi1$22393 $CellContext`phi2$22393 Compile`$156 - 
                    3.315541263197587 Compile`$142 Compile`$156 - 
                    0.22975292054736118` $CellContext`phi2$22393 Compile`$160 
                    Compile`$156 + 1.6577706315987935` Compile`$163 + 
                    0.22975292054736118` $CellContext`phi2$22393 Compile`$163)
                     Compile`$172 Compile`$182 Compile`$188 Compile`$192 
                    Compile`$196 Compile`$202 Compile`$210 Compile`$214 
                    Compile`$220 Compile`$224 Compile`$235 + 
                    30.616258720840584` (
                    2.1899864186871074` $CellContext`phi1$22393 \
$CellContext`phi2$22393 Compile`$156 + 
                    3.315541263197587 Compile`$142 Compile`$156 + 
                    0.22975292054736118` $CellContext`phi2$22393 Compile`$160 
                    Compile`$156 - 1.6577706315987935` Compile`$163 - 
                    0.22975292054736118` $CellContext`phi2$22393 Compile`$163)
                     Compile`$172 Compile`$182 Compile`$188 Compile`$192 
                    Compile`$196 Compile`$206 Compile`$210 Compile`$214 
                    Compile`$220 Compile`$224 Compile`$235 + 
                    30.616258720840584` (
                    3.918168105189718 $CellContext`phi1$22393 \
$CellContext`phi2$22393 Compile`$156 - 0.515564871781558 Compile`$142 
                    Compile`$156 + 
                    0.30909204230248366` $CellContext`phi2$22393 Compile`$160 
                    Compile`$156 + 0.257782435890779 Compile`$163 - 
                    0.30909204230248366` $CellContext`phi2$22393 Compile`$163)
                     Compile`$172 Compile`$182 Compile`$188 Compile`$192 
                    Compile`$202 Compile`$206 Compile`$210 Compile`$214 
                    Compile`$220 Compile`$224 Compile`$235 + 
                    30.616258720840584` ((-1.7341240244094855`) \
$CellContext`phi1$22393 $CellContext`phi2$22393 Compile`$156 + 
                    3.033802809341111 Compile`$142 Compile`$156 + 
                    0.9732489894677301 $CellContext`phi2$22393 Compile`$160 
                    Compile`$156 - 1.5169014046705556` Compile`$163 - 
                    0.9732489894677301 $CellContext`phi2$22393 Compile`$163) 
                    Compile`$172 Compile`$182 Compile`$188 Compile`$196 
                    Compile`$202 Compile`$206 Compile`$210 Compile`$214 
                    Compile`$220 Compile`$224 Compile`$235 + 
                    30.616258720840584` (
                    1.0621326830868534` $CellContext`phi1$22393 \
$CellContext`phi2$22393 Compile`$156 - 3.1650571315633877` Compile`$142 
                    Compile`$156 + 
                    1.1016223851050846` $CellContext`phi2$22393 Compile`$160 
                    Compile`$156 + 1.5825285657816939` Compile`$163 - 
                    1.1016223851050846` $CellContext`phi2$22393 Compile`$163) 
                    Compile`$172 Compile`$182 Compile`$192 Compile`$196 
                    Compile`$202 Compile`$206 Compile`$210 Compile`$214 
                    Compile`$220 Compile`$224 Compile`$235 + 
                    30.616258720840584` ((-2.43117596734898) \
$CellContext`phi1$22393 $CellContext`phi2$22393 Compile`$156 - 
                    0.9714272660591798 Compile`$142 Compile`$156 + 
                    1.512093952317575 $CellContext`phi2$22393 Compile`$160 
                    Compile`$156 + 0.4857136330295899 Compile`$163 - 
                    1.512093952317575 $CellContext`phi2$22393 Compile`$163) 
                    Compile`$172 Compile`$188 Compile`$192 Compile`$196 
                    Compile`$202 Compile`$206 Compile`$210 Compile`$214 
                    Compile`$220 Compile`$224 Compile`$235 + 
                    30.616258720840584` (
                    1.3438711369433292` $CellContext`phi1$22393 \
$CellContext`phi2$22393 Compile`$156 + 
                    0.7590533115332052 Compile`$142 Compile`$156 + 
                    1.8451184540254535` $CellContext`phi2$22393 Compile`$160 
                    Compile`$156 - 0.3795266557666026 Compile`$163 - 
                    1.8451184540254535` $CellContext`phi2$22393 Compile`$163) 
                    Compile`$182 Compile`$188 Compile`$192 Compile`$196 
                    Compile`$202 Compile`$206 Compile`$210 Compile`$214 
                    Compile`$220 Compile`$224 Compile`$235)^2)/4)/2]], 
             Evaluate], "Umax" -> 1.376|>, 
         "Dodecahedron" -> <|"U" -> CompiledFunction[{11, 14.3, 5470}, {
              Blank[Real], 
              Blank[Real]}, {{3, 0, 0}, {3, 0, 1}, {3, 0, 
             15}}, {{1.9955006726933704`, {3, 0, 54}}, {
              0.09175170953613698, {3, 0, 57}}, {-0.6392967992847324, {3, 0, 
               29}}, {-0.005279474328467426, {3, 0, 79}}, {
              341.07848212762275`, {3, 0, 11}}, {
              1.6953842579395209`, {3, 0, 36}}, {-1.4908780339819565`, {3, 0, 
               28}}, {-0.935172654052801, {3, 0, 
               16}}, {-0.021262810534837584`, {3, 0, 51}}, {
              0.20412414523193154`, {3, 0, 61}}, {-0.7083647052177124, {3, 0, 
               80}}, {-1.3263986939676073`, {3, 0, 23}}, {
              0.021262810534837584`, {3, 0, 55}}, {-0.5917517095361372, {3, 0,
                20}}, {-0.6828613346970941, {3, 0, 
               25}}, {-1.9955006726933704`, {3, 0, 
               50}}, {-0.20412414523193154`, {3, 0, 43}}, {
              0.005279474328467426, {3, 0, 12}}, {-0.5, {3, 0, 7}}, {
              0.6392967992847324, {3, 0, 70}}, {
              0.908248290463863, {3, 0, 77}}, {
              0.4564354645876384, {3, 0, 65}}, {
              1.5629488288431146`, {3, 0, 34}}, {
              0.7909944487358057, {3, 0, 42}}, {-0.09175170953613698, {3, 0, 
               47}}, {0.7464522479153887, {3, 0, 33}}, {
              0.5406837195602691, {3, 0, 64}}, {
              0.6828613346970941, {3, 0, 75}}, {
              1.790994448735806, {3, 0, 41}}, {-1.790994448735806, {3, 0, 
               59}}, {-1.6953842579395209`, {3, 0, 63}}, {
              0.5917517095361372, {3, 0, 19}}, {
              1.4908780339819565`, {3, 0, 69}}, {
              1.3263986939676073`, {3, 0, 73}}, {-1, {2, 0, 1}}, {
              1.4082482904638631`, {3, 0, 46}}, {-0.6127545144214275, {3, 0, 
               74}}, {-0.7464522479153887, {3, 0, 32}}, {
              0.7083647052177124, {3, 0, 14}}, {
              1, {2, 0, 0}}, {-0.12715595476793853`, {3, 0, 
               53}}, {-1.4082482904638631`, {3, 0, 45}}, {
              0.6127545144214275, {3, 0, 24}}, {-0.4564354645876384, {3, 0, 
               38}}, {-0.5406837195602691, {3, 0, 37}}, {
              0.37746668394347527`, {3, 0, 68}}, {
              0.12715595476793853`, {3, 0, 49}}, {-1.5629488288431146`, {3, 0,
                31}}, {
              0.5, {3, 0, 9}}, {-0.37746668394347527`, {3, 0, 
               27}}, {-0.7909944487358057, {3, 0, 60}}, {-0.908248290463863, {
               3, 0, 21}}, {0.935172654052801, {3, 0, 81}}}, {0, 2, 83, 0, 
             0}, {{40, 56, 3, 0, 0, 3, 0, 2}, {40, 56, 3, 0, 1, 3, 0, 3}, {10,
               0, 4}, {13, 4, 2, 3, 4}, {40, 60, 3, 0, 4, 3, 0, 5}, {10, 1, 
              6}, {13, 6, 2, 3, 6}, {16, 7, 6, 5, 8}, {16, 9, 6, 5, 10}, {16, 
              12, 0, 5, 13}, {16, 14, 1, 5, 15}, {16, 16, 6, 5, 17}, {10, 0, 
              18}, {13, 18, 13, 15, 17, 18}, {16, 19, 0, 5, 15}, {16, 20, 1, 
              5, 13}, {16, 21, 6, 5, 17}, {10, 0, 22}, {13, 22, 15, 13, 17, 
              22}, {16, 23, 0, 5, 13}, {16, 24, 1, 5, 15}, {16, 25, 6, 5, 
              17}, {10, 0, 26}, {13, 26, 13, 15, 17, 26}, {16, 27, 0, 5, 
              15}, {16, 28, 1, 5, 13}, {16, 29, 6, 5, 17}, {10, 0, 30}, {13, 
              30, 15, 13, 17, 30}, {16, 31, 0, 5, 13}, {16, 32, 1, 5, 15}, {
              10, 0, 17}, {13, 17, 13, 15, 8, 17}, {16, 33, 0, 5, 15}, {16, 
              34, 1, 5, 13}, {10, 0, 35}, {13, 35, 15, 13, 8, 35}, {16, 36, 0,
               5, 13}, {16, 37, 1, 5, 15}, {16, 38, 6, 5, 39}, {10, 0, 40}, {
              13, 40, 13, 15, 39, 40}, {16, 41, 0, 5, 15}, {16, 42, 1, 5, 
              13}, {16, 43, 6, 5, 39}, {10, 0, 44}, {13, 44, 15, 13, 39, 
              44}, {16, 45, 0, 5, 13}, {16, 46, 1, 5, 15}, {16, 47, 6, 5, 
              39}, {10, 0, 48}, {13, 48, 13, 15, 39, 48}, {16, 49, 0, 5, 
              15}, {16, 50, 1, 5, 13}, {16, 51, 6, 5, 39}, {10, 0, 52}, {13, 
              52, 15, 13, 39, 52}, {16, 53, 0, 5, 13}, {16, 54, 1, 5, 15}, {
              16, 55, 6, 5, 39}, {10, 0, 56}, {13, 56, 13, 15, 39, 56}, {16, 
              46, 0, 5, 13}, {16, 45, 1, 5, 15}, {16, 57, 6, 5, 39}, {10, 0, 
              58}, {13, 58, 13, 15, 39, 58}, {16, 59, 0, 5, 15}, {16, 60, 1, 
              5, 13}, {16, 61, 6, 5, 39}, {10, 0, 62}, {13, 62, 15, 13, 39, 
              62}, {16, 63, 0, 5, 13}, {16, 64, 1, 5, 15}, {16, 65, 6, 5, 
              39}, {10, 0, 66}, {13, 66, 13, 15, 39, 66}, {16, 32, 0, 5, 
              13}, {16, 31, 1, 5, 15}, {10, 0, 39}, {13, 39, 13, 15, 10, 
              39}, {16, 34, 0, 5, 13}, {16, 33, 1, 5, 15}, {10, 0, 67}, {13, 
              67, 13, 15, 10, 67}, {16, 68, 0, 5, 15}, {16, 69, 1, 5, 13}, {
              16, 70, 6, 5, 71}, {10, 0, 72}, {13, 72, 15, 13, 71, 72}, {16, 
              73, 0, 5, 13}, {16, 74, 1, 5, 15}, {16, 75, 6, 5, 71}, {10, 0, 
              76}, {13, 76, 13, 15, 71, 76}, {16, 20, 0, 5, 13}, {16, 19, 1, 
              5, 15}, {16, 77, 6, 5, 71}, {10, 0, 78}, {13, 78, 13, 15, 71, 
              78}, {16, 79, 0, 5, 15}, {16, 80, 1, 5, 13}, {16, 81, 6, 5, 
              71}, {10, 0, 82}, {13, 82, 15, 13, 71, 82}, {16, 11, 18, 22, 26,
               30, 17, 35, 40, 44, 48, 52, 56, 58, 62, 66, 39, 67, 72, 76, 78,
               82, 15}, {1}}, 
             Function[{$CellContext`phi1$22396, $CellContext`phi2$22396}, 
              
              Block[{Compile`$157, Compile`$158, Compile`$162, Compile`$164, 
                Compile`$225, Compile`$246, Compile`$285}, 
               Compile`$157 = $CellContext`phi1$22396^2; 
               Compile`$158 = $CellContext`phi2$22396^2; 
               Compile`$162 = 1 + Compile`$157 + Compile`$158; 
               Compile`$164 = Compile`$162^(-1); 
               Compile`$225 = -1 + Compile`$157 + Compile`$158; 
               Compile`$246 = (-0.5) Compile`$225 Compile`$164; 
               Compile`$285 = 0.5 Compile`$225 Compile`$164; 
               341.07848212762275` (1 + 
                 0.005279474328467426 $CellContext`phi1$22396 Compile`$164 + 
                 0.7083647052177124 $CellContext`phi2$22396 Compile`$164 - 
                 0.935172654052801 Compile`$225 Compile`$164) (1 + 
                 0.5917517095361372 $CellContext`phi1$22396 Compile`$164 - 
                 0.5917517095361372 $CellContext`phi2$22396 Compile`$164 - 
                 0.908248290463863 Compile`$225 Compile`$164) (1 - 
                 1.3263986939676073` $CellContext`phi1$22396 Compile`$164 + 
                 0.6127545144214275 $CellContext`phi2$22396 Compile`$164 - 
                 0.6828613346970941 Compile`$225 Compile`$164) (1 - 
                 0.37746668394347527` $CellContext`phi1$22396 Compile`$164 - 
                 1.4908780339819565` $CellContext`phi2$22396 Compile`$164 - 
                 0.6392967992847324 Compile`$225 Compile`$164) (1 - 
                 1.5629488288431146` $CellContext`phi1$22396 Compile`$164 - 
                 0.7464522479153887 $CellContext`phi2$22396 Compile`$164 + 
                 Compile`$246) (1 + 
                 0.7464522479153887 $CellContext`phi1$22396 Compile`$164 + 
                 1.5629488288431146` $CellContext`phi2$22396 Compile`$164 + 
                 Compile`$246) (1 + 
                 1.6953842579395209` $CellContext`phi1$22396 Compile`$164 - 
                 0.5406837195602691 $CellContext`phi2$22396 Compile`$164 - 
                 0.4564354645876384 Compile`$225 Compile`$164) (1 + 
                 1.790994448735806 $CellContext`phi1$22396 Compile`$164 + 
                 0.7909944487358057 $CellContext`phi2$22396 Compile`$164 - 
                 0.20412414523193154` Compile`$225 Compile`$164) (1 - 
                 1.4082482904638631` $CellContext`phi1$22396 Compile`$164 + 
                 1.4082482904638631` $CellContext`phi2$22396 Compile`$164 - 
                 0.09175170953613698 Compile`$225 Compile`$164) (1 + 
                 0.12715595476793853` $CellContext`phi1$22396 Compile`$164 - 
                 1.9955006726933704` $CellContext`phi2$22396 Compile`$164 - 
                 0.021262810534837584` Compile`$225 Compile`$164) (1 - 
                 0.12715595476793853` $CellContext`phi1$22396 Compile`$164 + 
                 1.9955006726933704` $CellContext`phi2$22396 Compile`$164 + 
                 0.021262810534837584` Compile`$225 Compile`$164) (1 + 
                 1.4082482904638631` $CellContext`phi1$22396 Compile`$164 - 
                 1.4082482904638631` $CellContext`phi2$22396 Compile`$164 + 
                 0.09175170953613698 Compile`$225 Compile`$164) (1 - 
                 1.790994448735806 $CellContext`phi1$22396 Compile`$164 - 
                 0.7909944487358057 $CellContext`phi2$22396 Compile`$164 + 
                 0.20412414523193154` Compile`$225 Compile`$164) (1 - 
                 1.6953842579395209` $CellContext`phi1$22396 Compile`$164 + 
                 0.5406837195602691 $CellContext`phi2$22396 Compile`$164 + 
                 0.4564354645876384 Compile`$225 Compile`$164) (1 - 
                 0.7464522479153887 $CellContext`phi1$22396 Compile`$164 - 
                 1.5629488288431146` $CellContext`phi2$22396 Compile`$164 + 
                 Compile`$285) (1 + 
                 1.5629488288431146` $CellContext`phi1$22396 Compile`$164 + 
                 0.7464522479153887 $CellContext`phi2$22396 Compile`$164 + 
                 Compile`$285) (1 + 
                 0.37746668394347527` $CellContext`phi1$22396 Compile`$164 + 
                 1.4908780339819565` $CellContext`phi2$22396 Compile`$164 + 
                 0.6392967992847324 Compile`$225 Compile`$164) (1 + 
                 1.3263986939676073` $CellContext`phi1$22396 Compile`$164 - 
                 0.6127545144214275 $CellContext`phi2$22396 Compile`$164 + 
                 0.6828613346970941 Compile`$225 Compile`$164) (1 - 
                 0.5917517095361372 $CellContext`phi1$22396 Compile`$164 + 
                 0.5917517095361372 $CellContext`phi2$22396 Compile`$164 + 
                 0.908248290463863 Compile`$225 Compile`$164) (1 - 
                 0.005279474328467426 $CellContext`phi1$22396 Compile`$164 - 
                 0.7083647052177124 $CellContext`phi2$22396 Compile`$164 + 
                 0.935172654052801 Compile`$225 Compile`$164)]], Evaluate], 
           "V" -> CompiledFunction[{11, 14.3, 5470}, {
              Blank[Real], 
              Blank[Real]}, {{3, 0, 0}, {3, 0, 1}, {3, 0, 
             174}}, {{-2.981756067963913, {3, 0, 191}}, {
              1.4929044958307773`, {3, 0, 149}}, {-0.9128709291752768, {3, 0, 
               197}}, {-0.5406837195602691, {3, 0, 52}}, {
              0.935172654052801, {3, 0, 133}}, {-0.12715595476793853`, {3, 0, 
               78}}, {0.010558948656934852`, {3, 0, 173}}, {
              1.0813674391205381`, {3, 0, 223}}, {
              1.2785935985694648`, {3, 0, 193}}, {-1.9955006726933704`, {3, 0,
                73}}, {1.5629488288431146`, {3, 0, 47}}, {1., {3, 0, 138}}, {
              2.981756067963913, {3, 0, 226}}, {
              0.5917517095361372, {3, 0, 19}}, {-0.7083647052177124, {3, 0, 
               131}}, {1.4082482904638631`, {3, 0, 66}}, {
              1.4167294104354249`, {3, 0, 175}}, {
              3.125897657686229, {3, 0, 163}}, {
              1.816496580927726, {3, 0, 183}}, {
              0.09175170953613698, {3, 0, 87}}, {-0.010558948656934852`, {3, 
               0, 233}}, {1.225509028842855, {3, 0, 185}}, {
              1.5819888974716114`, {3, 0, 203}}, {-0.6392967992847324, {3, 0, 
               37}}, {0.04252562106967517, {3, 0, 213}}, {
              2.6527973879352147`, {3, 0, 230}}, {-1.5819888974716114`, {3, 0,
                219}}, {-1.816496580927726, {3, 0, 
               182}}, {-1.1835034190722744`, {3, 0, 
               168}}, {-0.5917517095361372, {3, 0, 21}}, {
              341.07848212762275`, {3, 0, 172}}, {1, {2, 0, 0}}, {
              1.9955006726933704`, {3, 0, 80}}, {-0.7909944487358057, {3, 0, 
               92}}, {0.021262810534837584`, {3, 0, 82}}, {-0.5, {3, 0, 8}}, {
              0.7083647052177124, {3, 0, 14}}, {-0.6127545144214275, {3, 0, 
               119}}, {-0.20412414523193154`, {3, 0, 
               61}}, {-1.3263986939676073`, {3, 0, 26}}, {
              1.3263986939676073`, {3, 0, 117}}, {-0.005279474328467426, {3, 
               0, 129}}, {-1.225509028842855, {3, 0, 
               231}}, {-0.7549333678869505, {3, 0, 
               190}}, {-1.790994448735806, {3, 0, 90}}, {
              0.6828613346970941, {3, 0, 121}}, {-0.37746668394347527`, {3, 0,
                33}}, {-0.4082482904638631, {3, 0, 204}}, {
              0.4564354645876384, {3, 0, 101}}, {-0.4564354645876384, {3, 0, 
               54}}, {0.7464522479153887, {3, 0, 45}}, {
              0.37746668394347527`, {3, 0, 110}}, {
              0.18350341907227397`, {3, 0, 208}}, {
              0.25431190953587707`, {3, 0, 209}}, {-1.3657226693941882`, {3, 
               0, 186}}, {0.5, {3, 0, 10}}, {
              0.5406837195602691, {3, 0, 99}}, {
              3.581988897471612, {3, 0, 202}}, {
              1.790994448735806, {3, 0, 57}}, {
              0.005279474328467426, {3, 0, 12}}, {
              0.12715595476793853`, {3, 0, 71}}, {-0.935172654052801, {3, 0, 
               16}}, {0.9128709291752768, {3, 0, 
               200}}, {-0.25431190953587707`, {3, 0, 
               215}}, {-0.18350341907227397`, {3, 0, 
               207}}, {-1.4908780339819565`, {3, 0, 35}}, {
              3.3907685158790417`, {3, 0, 195}}, {
              0.20412414523193154`, {3, 0, 94}}, {
              0.6127545144214275, {3, 0, 28}}, {
              1.1835034190722744`, {3, 0, 143}}, {-1.870345308105602, {3, 0, 
               177}}, {-1.2785935985694648`, {3, 0, 192}}, {
              1.3657226693941882`, {3, 0, 188}}, {
              0.4082482904638631, {3, 0, 205}}, {-3.125897657686229, {3, 0, 
               146}}, {
              1.870345308105602, {3, 0, 180}}, {-1, {2, 0, 
               1}}, {-0.09175170953613698, {3, 0, 68}}, {
              0.6392967992847324, {3, 0, 114}}, {
              3.991001345386741, {3, 0, 216}}, {-2.8164965809277263`, {3, 0, 
               154}}, {-1.4929044958307773`, {3, 0, 160}}, {
              0.7549333678869505, {3, 0, 225}}, {-2.6527973879352147`, {3, 0, 
               184}}, {
              0.7909944487358057, {3, 0, 59}}, {-0.7464522479153887, {3, 0, 
               42}}, {-1.0813674391205381`, {3, 0, 
               196}}, {-1.6953842579395209`, {3, 0, 
               97}}, {-0.908248290463863, {3, 0, 23}}, {
              1.4908780339819565`, {3, 0, 112}}, {-1., {3, 0, 
               136}}, {-0.021262810534837584`, {3, 0, 75}}, {
              0.908248290463863, {3, 0, 126}}, {-3.581988897471612, {3, 0, 
               218}}, {2.8164965809277263`, {3, 0, 157}}, {
              0.25, {3, 0, 171}}, {-1.4082482904638631`, {3, 0, 
               64}}, {-0.6828613346970941, {3, 0, 30}}, {-3.991001345386741, {
               3, 0, 210}}, {
              1.6953842579395209`, {3, 0, 50}}, {-0.04252562106967517, {3, 0, 
               211}}, {-3.3907685158790417`, {3, 0, 
               222}}, {-1.5629488288431146`, {3, 0, 
               40}}, {-1.4167294104354249`, {3, 0, 234}}}, {0, 2, 238, 0, 
             0}, {{40, 56, 3, 0, 0, 3, 0, 2}, {40, 56, 3, 0, 1, 3, 0, 3}, {10,
               0, 4}, {13, 4, 2, 3, 4}, {40, 56, 3, 0, 4, 3, 0, 5}, {40, 60, 
              3, 0, 5, 3, 0, 6}, {40, 60, 3, 0, 4, 3, 0, 5}, {10, 1, 7}, {13, 
              7, 2, 3, 7}, {16, 8, 7, 5, 9}, {16, 10, 7, 5, 11}, {16, 12, 0, 
              5, 13}, {16, 14, 1, 5, 15}, {16, 16, 7, 5, 17}, {10, 0, 18}, {
              13, 18, 13, 15, 17, 18}, {16, 19, 0, 5, 20}, {16, 21, 1, 5, 
              22}, {16, 23, 7, 5, 24}, {10, 0, 25}, {13, 25, 20, 22, 24, 
              25}, {16, 26, 0, 5, 27}, {16, 28, 1, 5, 29}, {16, 30, 7, 5, 
              31}, {10, 0, 32}, {13, 32, 27, 29, 31, 32}, {16, 33, 0, 5, 
              34}, {16, 35, 1, 5, 36}, {16, 37, 7, 5, 38}, {10, 0, 39}, {13, 
              39, 34, 36, 38, 39}, {16, 40, 0, 5, 41}, {16, 42, 1, 5, 43}, {
              10, 0, 44}, {13, 44, 41, 43, 9, 44}, {16, 45, 0, 5, 46}, {16, 
              47, 1, 5, 48}, {10, 0, 49}, {13, 49, 46, 48, 9, 49}, {16, 50, 0,
               5, 51}, {16, 52, 1, 5, 53}, {16, 54, 7, 5, 55}, {10, 0, 56}, {
              13, 56, 51, 53, 55, 56}, {16, 57, 0, 5, 58}, {16, 59, 1, 5, 
              60}, {16, 61, 7, 5, 62}, {10, 0, 63}, {13, 63, 58, 60, 62, 
              63}, {16, 64, 0, 5, 65}, {16, 66, 1, 5, 67}, {16, 68, 7, 5, 
              69}, {10, 0, 70}, {13, 70, 65, 67, 69, 70}, {16, 71, 0, 5, 
              72}, {16, 73, 1, 5, 74}, {16, 75, 7, 5, 76}, {10, 0, 77}, {13, 
              77, 72, 74, 76, 77}, {16, 78, 0, 5, 79}, {16, 80, 1, 5, 81}, {
              16, 82, 7, 5, 83}, {10, 0, 84}, {13, 84, 79, 81, 83, 84}, {16, 
              66, 0, 5, 85}, {16, 64, 1, 5, 86}, {16, 87, 7, 5, 88}, {10, 0, 
              89}, {13, 89, 85, 86, 88, 89}, {16, 90, 0, 5, 91}, {16, 92, 1, 
              5, 93}, {16, 94, 7, 5, 95}, {10, 0, 96}, {13, 96, 91, 93, 95, 
              96}, {16, 97, 0, 5, 98}, {16, 99, 1, 5, 100}, {16, 101, 7, 5, 
              102}, {10, 0, 103}, {13, 103, 98, 100, 102, 103}, {16, 42, 0, 5,
               104}, {16, 40, 1, 5, 105}, {10, 0, 106}, {13, 106, 104, 105, 
              11, 106}, {16, 47, 0, 5, 107}, {16, 45, 1, 5, 108}, {10, 0, 
              109}, {13, 109, 107, 108, 11, 109}, {16, 110, 0, 5, 111}, {16, 
              112, 1, 5, 113}, {16, 114, 7, 5, 115}, {10, 0, 116}, {13, 116, 
              111, 113, 115, 116}, {16, 117, 0, 5, 118}, {16, 119, 1, 5, 
              120}, {16, 121, 7, 5, 122}, {10, 0, 123}, {13, 123, 118, 120, 
              122, 123}, {16, 21, 0, 5, 124}, {16, 19, 1, 5, 125}, {16, 126, 
              7, 5, 127}, {10, 0, 128}, {13, 128, 124, 125, 127, 128}, {16, 
              129, 0, 5, 130}, {16, 131, 1, 5, 132}, {16, 133, 7, 5, 134}, {
              10, 0, 135}, {13, 135, 130, 132, 134, 135}, {16, 136, 0, 7, 6, 
              137}, {16, 138, 0, 5, 139}, {16, 138, 0, 7, 6, 140}, {16, 136, 
              0, 5, 141}, {40, 56, 3, 0, 4, 3, 0, 142}, {16, 143, 0, 1, 6, 
              144}, {16, 19, 5, 145}, {16, 146, 0, 1, 6, 147}, {16, 45, 5, 
              148}, {16, 149, 0, 1, 6, 150}, {16, 136, 1, 7, 6, 151}, {16, 40,
               5, 152}, {16, 138, 1, 5, 153}, {16, 154, 0, 1, 6, 155}, {16, 
              64, 5, 156}, {16, 157, 0, 1, 6, 158}, {16, 66, 5, 159}, {16, 
              160, 0, 1, 6, 161}, {16, 47, 5, 162}, {16, 163, 0, 1, 6, 164}, {
              16, 138, 1, 7, 6, 165}, {16, 42, 5, 166}, {16, 136, 1, 5, 
              167}, {16, 168, 0, 1, 6, 169}, {16, 21, 5, 170}, {16, 173, 2, 6,
               174}, {16, 175, 0, 1, 6, 176}, {16, 177, 0, 7, 6, 178}, {16, 
              129, 5, 179}, {16, 180, 0, 5, 181}, {13, 174, 176, 178, 179, 
              181, 174}, {16, 172, 174, 18, 25, 32, 39, 44, 49, 56, 63, 70, 
              77, 84, 89, 96, 103, 106, 109, 116, 123, 128, 176}, {16, 143, 2,
               6, 174}, {16, 182, 0, 7, 6, 179}, {16, 183, 0, 5, 178}, {13, 
              174, 169, 179, 170, 178, 174}, {16, 172, 174, 18, 25, 32, 39, 
              44, 49, 56, 63, 70, 77, 84, 89, 96, 103, 106, 109, 116, 123, 
              135, 179}, {16, 184, 2, 6, 178}, {16, 185, 0, 1, 6, 174}, {16, 
              186, 0, 7, 6, 181}, {16, 117, 5, 187}, {16, 188, 0, 5, 189}, {
              13, 178, 174, 181, 187, 189, 178}, {16, 172, 178, 18, 25, 32, 
              39, 44, 49, 56, 63, 70, 77, 84, 89, 96, 103, 106, 109, 116, 128,
               135, 174}, {16, 190, 2, 6, 181}, {16, 191, 0, 1, 6, 189}, {16, 
              192, 0, 7, 6, 187}, {16, 110, 5, 178}, {16, 193, 0, 5, 194}, {
              13, 181, 189, 187, 178, 194, 181}, {16, 172, 181, 18, 25, 32, 
              39, 44, 49, 56, 63, 70, 77, 84, 89, 96, 103, 106, 109, 123, 128,
               135, 189}, {16, 146, 2, 6, 181}, {13, 181, 161, 137, 162, 139, 
              181}, {16, 172, 181, 18, 25, 32, 39, 44, 49, 56, 63, 70, 77, 84,
               89, 96, 103, 106, 116, 123, 128, 135, 187}, {16, 149, 2, 6, 
              181}, {13, 181, 164, 137, 166, 139, 181}, {16, 172, 181, 18, 25,
               32, 39, 44, 49, 56, 63, 70, 77, 84, 89, 96, 103, 109, 116, 123,
               128, 135, 178}, {16, 195, 2, 6, 194}, {16, 196, 0, 1, 6, 
              181}, {16, 197, 0, 7, 6, 198}, {16, 97, 5, 199}, {16, 200, 0, 5,
               201}, {13, 194, 181, 198, 199, 201, 194}, {16, 172, 194, 18, 
              25, 32, 39, 44, 49, 56, 63, 70, 77, 84, 89, 96, 106, 109, 116, 
              123, 128, 135, 181}, {16, 202, 2, 6, 198}, {16, 203, 0, 1, 6, 
              201}, {16, 204, 0, 7, 6, 199}, {16, 90, 5, 194}, {16, 205, 0, 5,
               206}, {13, 198, 201, 199, 194, 206, 198}, {16, 172, 198, 18, 
              25, 32, 39, 44, 49, 56, 63, 70, 77, 84, 89, 103, 106, 109, 116, 
              123, 128, 135, 201}, {16, 154, 2, 6, 198}, {16, 207, 0, 7, 6, 
              194}, {16, 208, 0, 5, 199}, {13, 198, 158, 194, 159, 199, 
              198}, {16, 172, 198, 18, 25, 32, 39, 44, 49, 56, 63, 70, 77, 84,
               96, 103, 106, 109, 116, 123, 128, 135, 194}, {16, 209, 2, 6, 
              199}, {16, 210, 0, 1, 6, 198}, {16, 211, 0, 7, 6, 206}, {16, 78,
               5, 212}, {16, 213, 0, 5, 214}, {13, 199, 198, 206, 212, 214, 
              199}, {16, 172, 199, 18, 25, 32, 39, 44, 49, 56, 63, 70, 77, 89,
               96, 103, 106, 109, 116, 123, 128, 135, 198}, {16, 215, 2, 6, 
              206}, {16, 216, 0, 1, 6, 214}, {16, 213, 0, 7, 6, 199}, {16, 71,
               5, 212}, {16, 211, 0, 5, 217}, {13, 206, 214, 199, 212, 217, 
              206}, {16, 172, 206, 18, 25, 32, 39, 44, 49, 56, 63, 70, 84, 89,
               96, 103, 106, 109, 116, 123, 128, 135, 214}, {16, 157, 2, 6, 
              206}, {16, 208, 0, 7, 6, 199}, {16, 207, 0, 5, 212}, {13, 206, 
              155, 199, 156, 212, 206}, {16, 172, 206, 18, 25, 32, 39, 44, 49,
               56, 63, 77, 84, 89, 96, 103, 106, 109, 116, 123, 128, 135, 
              199}, {16, 218, 2, 6, 212}, {16, 219, 0, 1, 6, 206}, {16, 205, 
              0, 7, 6, 217}, {16, 57, 5, 220}, {16, 204, 0, 5, 221}, {13, 212,
               206, 217, 220, 221, 212}, {16, 172, 212, 18, 25, 32, 39, 44, 
              49, 56, 70, 77, 84, 89, 96, 103, 106, 109, 116, 123, 128, 135, 
              206}, {16, 222, 2, 6, 217}, {16, 223, 0, 1, 6, 221}, {16, 200, 
              0, 7, 6, 212}, {16, 50, 5, 220}, {16, 197, 0, 5, 224}, {13, 217,
               221, 212, 220, 224, 217}, {16, 172, 217, 18, 25, 32, 39, 44, 
              49, 63, 70, 77, 84, 89, 96, 103, 106, 109, 116, 123, 128, 135, 
              221}, {16, 160, 2, 6, 217}, {13, 217, 147, 140, 148, 141, 
              217}, {16, 172, 217, 18, 25, 32, 39, 44, 56, 63, 70, 77, 84, 89,
               96, 103, 106, 109, 116, 123, 128, 135, 212}, {16, 163, 2, 6, 
              217}, {13, 217, 150, 140, 152, 141, 217}, {16, 172, 217, 18, 25,
               32, 39, 49, 56, 63, 70, 77, 84, 89, 96, 103, 106, 109, 116, 
              123, 128, 135, 220}, {16, 225, 2, 6, 224}, {16, 226, 0, 1, 6, 
              217}, {16, 193, 0, 7, 6, 227}, {16, 33, 5, 228}, {16, 192, 0, 5,
               229}, {13, 224, 217, 227, 228, 229, 224}, {16, 172, 224, 18, 
              25, 32, 44, 49, 56, 63, 70, 77, 84, 89, 96, 103, 106, 109, 116, 
              123, 128, 135, 217}, {16, 230, 2, 6, 227}, {16, 231, 0, 1, 6, 
              229}, {16, 188, 0, 7, 6, 224}, {16, 26, 5, 228}, {16, 186, 0, 5,
               232}, {13, 227, 229, 224, 228, 232, 227}, {16, 172, 227, 18, 
              25, 39, 44, 49, 56, 63, 70, 77, 84, 89, 96, 103, 106, 109, 116, 
              123, 128, 135, 229}, {16, 168, 2, 6, 227}, {16, 183, 0, 7, 6, 
              224}, {16, 182, 0, 5, 228}, {13, 227, 144, 224, 145, 228, 
              227}, {16, 172, 227, 18, 32, 39, 44, 49, 56, 63, 70, 77, 84, 89,
               96, 103, 106, 109, 116, 123, 128, 135, 224}, {16, 233, 2, 6, 
              228}, {16, 234, 0, 1, 6, 227}, {16, 180, 0, 7, 6, 232}, {16, 12,
               5, 235}, {16, 177, 0, 5, 236}, {13, 228, 227, 232, 235, 236, 
              228}, {16, 172, 228, 25, 32, 39, 44, 49, 56, 63, 70, 77, 84, 89,
               96, 103, 106, 109, 116, 123, 128, 135, 227}, {13, 176, 179, 
              174, 189, 187, 178, 181, 201, 194, 198, 214, 199, 206, 221, 212,
               220, 217, 229, 224, 227, 176}, {40, 56, 3, 0, 176, 3, 0, 
              179}, {16, 171, 142, 179, 176}, {16, 173, 0, 1, 6, 179}, {16, 
              175, 3, 6, 174}, {16, 177, 1, 7, 6, 189}, {16, 131, 5, 187}, {
              16, 180, 1, 5, 178}, {13, 179, 174, 189, 187, 178, 179}, {16, 
              172, 179, 18, 25, 32, 39, 44, 49, 56, 63, 70, 77, 84, 89, 96, 
              103, 106, 109, 116, 123, 128, 174}, {16, 168, 3, 6, 179}, {16, 
              182, 1, 7, 6, 189}, {16, 183, 1, 5, 187}, {13, 144, 179, 189, 
              145, 187, 178}, {16, 172, 178, 18, 25, 32, 39, 44, 49, 56, 63, 
              70, 77, 84, 89, 96, 103, 106, 109, 116, 123, 135, 179}, {16, 
              184, 0, 1, 6, 178}, {16, 185, 3, 6, 189}, {16, 186, 1, 7, 6, 
              187}, {16, 119, 5, 181}, {16, 188, 1, 5, 201}, {13, 178, 189, 
              187, 181, 201, 178}, {16, 172, 178, 18, 25, 32, 39, 44, 49, 56, 
              63, 70, 77, 84, 89, 96, 103, 106, 109, 116, 128, 135, 189}, {16,
               190, 0, 1, 6, 178}, {16, 191, 3, 6, 187}, {16, 192, 1, 7, 6, 
              181}, {16, 112, 5, 201}, {16, 193, 1, 5, 194}, {13, 178, 187, 
              181, 201, 194, 178}, {16, 172, 178, 18, 25, 32, 39, 44, 49, 56, 
              63, 70, 77, 84, 89, 96, 103, 106, 109, 123, 128, 135, 187}, {16,
               160, 3, 6, 178}, {13, 147, 178, 151, 148, 153, 181}, {16, 172, 
              181, 18, 25, 32, 39, 44, 49, 56, 63, 70, 77, 84, 89, 96, 103, 
              106, 116, 123, 128, 135, 178}, {16, 163, 3, 6, 181}, {13, 150, 
              181, 151, 152, 153, 201}, {16, 172, 201, 18, 25, 32, 39, 44, 49,
               56, 63, 70, 77, 84, 89, 96, 103, 109, 116, 123, 128, 135, 
              181}, {16, 195, 0, 1, 6, 201}, {16, 196, 3, 6, 194}, {16, 197, 
              1, 7, 6, 198}, {16, 99, 5, 214}, {16, 200, 1, 5, 199}, {13, 201,
               194, 198, 214, 199, 201}, {16, 172, 201, 18, 25, 32, 39, 44, 
              49, 56, 63, 70, 77, 84, 89, 96, 106, 109, 116, 123, 128, 135, 
              194}, {16, 202, 0, 1, 6, 201}, {16, 203, 3, 6, 198}, {16, 204, 
              1, 7, 6, 214}, {16, 92, 5, 199}, {16, 205, 1, 5, 206}, {13, 201,
               198, 214, 199, 206, 201}, {16, 172, 201, 18, 25, 32, 39, 44, 
              49, 56, 63, 70, 77, 84, 89, 103, 106, 109, 116, 123, 128, 135, 
              198}, {16, 157, 3, 6, 201}, {16, 207, 1, 7, 6, 214}, {16, 208, 
              1, 5, 199}, {13, 155, 201, 214, 156, 199, 206}, {16, 172, 206, 
              18, 25, 32, 39, 44, 49, 56, 63, 70, 77, 84, 96, 103, 106, 109, 
              116, 123, 128, 135, 201}, {16, 209, 0, 1, 6, 206}, {16, 210, 3, 
              6, 214}, {16, 211, 1, 7, 6, 199}, {16, 80, 5, 221}, {16, 213, 1,
               5, 212}, {13, 206, 214, 199, 221, 212, 206}, {16, 172, 206, 18,
               25, 32, 39, 44, 49, 56, 63, 70, 77, 89, 96, 103, 106, 109, 116,
               123, 128, 135, 214}, {16, 215, 0, 1, 6, 206}, {16, 216, 3, 6, 
              199}, {16, 213, 1, 7, 6, 221}, {16, 73, 5, 212}, {16, 211, 1, 5,
               220}, {13, 206, 199, 221, 212, 220, 206}, {16, 172, 206, 18, 
              25, 32, 39, 44, 49, 56, 63, 70, 84, 89, 96, 103, 106, 109, 116, 
              123, 128, 135, 199}, {16, 154, 3, 6, 206}, {16, 208, 1, 7, 6, 
              221}, {16, 207, 1, 5, 212}, {13, 158, 206, 221, 159, 212, 
              220}, {16, 172, 220, 18, 25, 32, 39, 44, 49, 56, 63, 77, 84, 89,
               96, 103, 106, 109, 116, 123, 128, 135, 206}, {16, 218, 0, 1, 6,
               220}, {16, 219, 3, 6, 221}, {16, 205, 1, 7, 6, 212}, {16, 59, 
              5, 217}, {16, 204, 1, 5, 229}, {13, 220, 221, 212, 217, 229, 
              220}, {16, 172, 220, 18, 25, 32, 39, 44, 49, 56, 70, 77, 84, 89,
               96, 103, 106, 109, 116, 123, 128, 135, 221}, {16, 222, 0, 1, 6,
               220}, {16, 223, 3, 6, 212}, {16, 200, 1, 7, 6, 217}, {16, 52, 
              5, 229}, {16, 197, 1, 5, 224}, {13, 220, 212, 217, 229, 224, 
              220}, {16, 172, 220, 18, 25, 32, 39, 44, 49, 63, 70, 77, 84, 89,
               96, 103, 106, 109, 116, 123, 128, 135, 212}, {16, 146, 3, 6, 
              220}, {13, 161, 220, 165, 162, 167, 217}, {16, 172, 217, 18, 25,
               32, 39, 44, 56, 63, 70, 77, 84, 89, 96, 103, 106, 109, 116, 
              123, 128, 135, 220}, {16, 149, 3, 6, 217}, {13, 164, 217, 165, 
              166, 167, 229}, {16, 172, 229, 18, 25, 32, 39, 49, 56, 63, 70, 
              77, 84, 89, 96, 103, 106, 109, 116, 123, 128, 135, 217}, {16, 
              225, 0, 1, 6, 229}, {16, 226, 3, 6, 224}, {16, 193, 1, 7, 6, 
              227}, {16, 35, 5, 228}, {16, 192, 1, 5, 232}, {13, 229, 224, 
              227, 228, 232, 229}, {16, 172, 229, 18, 25, 32, 44, 49, 56, 63, 
              70, 77, 84, 89, 96, 103, 106, 109, 116, 123, 128, 135, 224}, {
              16, 230, 0, 1, 6, 229}, {16, 231, 3, 6, 227}, {16, 188, 1, 7, 6,
               228}, {16, 28, 5, 232}, {16, 186, 1, 5, 235}, {13, 229, 227, 
              228, 232, 235, 229}, {16, 172, 229, 18, 25, 39, 44, 49, 56, 63, 
              70, 77, 84, 89, 96, 103, 106, 109, 116, 123, 128, 135, 227}, {
              16, 143, 3, 6, 229}, {16, 183, 1, 7, 6, 228}, {16, 182, 1, 5, 
              232}, {13, 169, 229, 228, 170, 232, 235}, {16, 172, 235, 18, 32,
               39, 44, 49, 56, 63, 70, 77, 84, 89, 96, 103, 106, 109, 116, 
              123, 128, 135, 229}, {16, 233, 0, 1, 6, 235}, {16, 234, 3, 6, 
              228}, {16, 180, 1, 7, 6, 232}, {16, 14, 5, 236}, {16, 177, 1, 5,
               237}, {13, 235, 228, 232, 236, 237, 235}, {16, 172, 235, 25, 
              32, 39, 44, 49, 56, 63, 70, 77, 84, 89, 96, 103, 106, 109, 116, 
              123, 128, 135, 228}, {13, 174, 179, 189, 187, 178, 181, 194, 
              198, 201, 214, 199, 206, 221, 212, 220, 217, 224, 227, 229, 228,
               174}, {40, 56, 3, 0, 174, 3, 0, 179}, {16, 171, 142, 179, 
              174}, {13, 176, 174, 176}, {16, 10, 176, 174}, {1}}, 
             Function[{$CellContext`phi1$22396, $CellContext`phi2$22396}, 
              
              Block[{Compile`$165, Compile`$167, Compile`$226, Compile`$228, 
                Compile`$237, Compile`$231, Compile`$260, Compile`$300, 
                Compile`$241, Compile`$242, Compile`$243, Compile`$244, 
                Compile`$245, Compile`$247, Compile`$248, Compile`$249, 
                Compile`$250, Compile`$251, Compile`$252, Compile`$253, 
                Compile`$254, Compile`$255, Compile`$256, Compile`$257, 
                Compile`$258, Compile`$259, Compile`$261, Compile`$262, 
                Compile`$263, Compile`$264, Compile`$265, Compile`$266, 
                Compile`$267, Compile`$268, Compile`$269, Compile`$270, 
                Compile`$271, Compile`$272, Compile`$273, Compile`$274, 
                Compile`$275, Compile`$276, Compile`$277, Compile`$278, 
                Compile`$279, Compile`$280, Compile`$281, Compile`$282, 
                Compile`$283, Compile`$284, Compile`$286, Compile`$287, 
                Compile`$288, Compile`$289, Compile`$290, Compile`$291, 
                Compile`$292, Compile`$293, Compile`$294, Compile`$295, 
                Compile`$296, Compile`$297, Compile`$298, Compile`$299, 
                Compile`$301, Compile`$302, Compile`$303, Compile`$304, 
                Compile`$305, Compile`$306, Compile`$307, Compile`$308, 
                Compile`$309, Compile`$310, Compile`$311, Compile`$312, 
                Compile`$313, Compile`$314, Compile`$315, Compile`$316, 
                Compile`$324, Compile`$325, Compile`$326, Compile`$327, 
                Compile`$345, Compile`$347, Compile`$413, Compile`$415, 
                Compile`$227, Compile`$438, Compile`$440, Compile`$412, 
                Compile`$414, Compile`$419, Compile`$481, Compile`$420, 
                Compile`$482, Compile`$391, Compile`$393, Compile`$370, 
                Compile`$372, Compile`$344, Compile`$346, Compile`$351, 
                Compile`$541, Compile`$352, Compile`$542, Compile`$319, 
                Compile`$321}, Compile`$165 = $CellContext`phi1$22396^2; 
               Compile`$167 = $CellContext`phi2$22396^2; 
               Compile`$226 = 1 + Compile`$165 + Compile`$167; 
               Compile`$228 = Compile`$226^(-2); 
               Compile`$237 = Compile`$226^(-1); 
               Compile`$231 = -1 + Compile`$165 + Compile`$167; 
               Compile`$260 = (-0.5) Compile`$231 Compile`$237; 
               Compile`$300 = 0.5 Compile`$231 Compile`$237; 
               Compile`$241 = 
                0.005279474328467426 $CellContext`phi1$22396 Compile`$237; 
               Compile`$242 = 
                0.7083647052177124 $CellContext`phi2$22396 Compile`$237; 
               Compile`$243 = (-0.935172654052801) Compile`$231 Compile`$237; 
               Compile`$244 = 1 + Compile`$241 + Compile`$242 + Compile`$243; 
               Compile`$245 = 
                0.5917517095361372 $CellContext`phi1$22396 Compile`$237; 
               Compile`$247 = (-0.5917517095361372) $CellContext`phi2$22396 
                 Compile`$237; 
               Compile`$248 = (-0.908248290463863) Compile`$231 Compile`$237; 
               Compile`$249 = 1 + Compile`$245 + Compile`$247 + Compile`$248; 
               Compile`$250 = (-1.3263986939676073`) $CellContext`phi1$22396 
                 Compile`$237; 
               Compile`$251 = 
                0.6127545144214275 $CellContext`phi2$22396 Compile`$237; 
               Compile`$252 = (-0.6828613346970941) Compile`$231 Compile`$237; 
               Compile`$253 = 1 + Compile`$250 + Compile`$251 + Compile`$252; 
               Compile`$254 = (-0.37746668394347527`) $CellContext`phi1$22396 
                 Compile`$237; 
               Compile`$255 = (-1.4908780339819565`) $CellContext`phi2$22396 
                 Compile`$237; 
               Compile`$256 = (-0.6392967992847324) Compile`$231 Compile`$237; 
               Compile`$257 = 1 + Compile`$254 + Compile`$255 + Compile`$256; 
               Compile`$258 = (-1.5629488288431146`) $CellContext`phi1$22396 
                 Compile`$237; 
               Compile`$259 = (-0.7464522479153887) $CellContext`phi2$22396 
                 Compile`$237; 
               Compile`$261 = 1 + Compile`$258 + Compile`$259 + Compile`$260; 
               Compile`$262 = 
                0.7464522479153887 $CellContext`phi1$22396 Compile`$237; 
               Compile`$263 = 
                1.5629488288431146` $CellContext`phi2$22396 Compile`$237; 
               Compile`$264 = 1 + Compile`$262 + Compile`$263 + Compile`$260; 
               Compile`$265 = 
                1.6953842579395209` $CellContext`phi1$22396 Compile`$237; 
               Compile`$266 = (-0.5406837195602691) $CellContext`phi2$22396 
                 Compile`$237; 
               Compile`$267 = (-0.4564354645876384) Compile`$231 Compile`$237; 
               Compile`$268 = 1 + Compile`$265 + Compile`$266 + Compile`$267; 
               Compile`$269 = 
                1.790994448735806 $CellContext`phi1$22396 Compile`$237; 
               Compile`$270 = 
                0.7909944487358057 $CellContext`phi2$22396 Compile`$237; 
               Compile`$271 = (-0.20412414523193154`) Compile`$231 
                 Compile`$237; 
               Compile`$272 = 1 + Compile`$269 + Compile`$270 + Compile`$271; 
               Compile`$273 = (-1.4082482904638631`) $CellContext`phi1$22396 
                 Compile`$237; 
               Compile`$274 = 
                1.4082482904638631` $CellContext`phi2$22396 Compile`$237; 
               Compile`$275 = (-0.09175170953613698) Compile`$231 
                 Compile`$237; 
               Compile`$276 = 1 + Compile`$273 + Compile`$274 + Compile`$275; 
               Compile`$277 = 
                0.12715595476793853` $CellContext`phi1$22396 Compile`$237; 
               Compile`$278 = (-1.9955006726933704`) $CellContext`phi2$22396 
                 Compile`$237; 
               Compile`$279 = (-0.021262810534837584`) Compile`$231 
                 Compile`$237; 
               Compile`$280 = 1 + Compile`$277 + Compile`$278 + Compile`$279; 
               Compile`$281 = (-0.12715595476793853`) $CellContext`phi1$22396 
                 Compile`$237; 
               Compile`$282 = 
                1.9955006726933704` $CellContext`phi2$22396 Compile`$237; 
               Compile`$283 = 0.021262810534837584` Compile`$231 Compile`$237; 
               Compile`$284 = 1 + Compile`$281 + Compile`$282 + Compile`$283; 
               Compile`$286 = 
                1.4082482904638631` $CellContext`phi1$22396 Compile`$237; 
               Compile`$287 = (-1.4082482904638631`) $CellContext`phi2$22396 
                 Compile`$237; 
               Compile`$288 = 0.09175170953613698 Compile`$231 Compile`$237; 
               Compile`$289 = 1 + Compile`$286 + Compile`$287 + Compile`$288; 
               Compile`$290 = (-1.790994448735806) $CellContext`phi1$22396 
                 Compile`$237; 
               Compile`$291 = (-0.7909944487358057) $CellContext`phi2$22396 
                 Compile`$237; 
               Compile`$292 = 0.20412414523193154` Compile`$231 Compile`$237; 
               Compile`$293 = 1 + Compile`$290 + Compile`$291 + Compile`$292; 
               Compile`$294 = (-1.6953842579395209`) $CellContext`phi1$22396 
                 Compile`$237; 
               Compile`$295 = 
                0.5406837195602691 $CellContext`phi2$22396 Compile`$237; 
               Compile`$296 = 0.4564354645876384 Compile`$231 Compile`$237; 
               Compile`$297 = 1 + Compile`$294 + Compile`$295 + Compile`$296; 
               Compile`$298 = (-0.7464522479153887) $CellContext`phi1$22396 
                 Compile`$237; 
               Compile`$299 = (-1.5629488288431146`) $CellContext`phi2$22396 
                 Compile`$237; 
               Compile`$301 = 1 + Compile`$298 + Compile`$299 + Compile`$300; 
               Compile`$302 = 
                1.5629488288431146` $CellContext`phi1$22396 Compile`$237; 
               Compile`$303 = 
                0.7464522479153887 $CellContext`phi2$22396 Compile`$237; 
               Compile`$304 = 1 + Compile`$302 + Compile`$303 + Compile`$300; 
               Compile`$305 = 
                0.37746668394347527` $CellContext`phi1$22396 Compile`$237; 
               Compile`$306 = 
                1.4908780339819565` $CellContext`phi2$22396 Compile`$237; 
               Compile`$307 = 0.6392967992847324 Compile`$231 Compile`$237; 
               Compile`$308 = 1 + Compile`$305 + Compile`$306 + Compile`$307; 
               Compile`$309 = 
                1.3263986939676073` $CellContext`phi1$22396 Compile`$237; 
               Compile`$310 = (-0.6127545144214275) $CellContext`phi2$22396 
                 Compile`$237; 
               Compile`$311 = 0.6828613346970941 Compile`$231 Compile`$237; 
               Compile`$312 = 1 + Compile`$309 + Compile`$310 + Compile`$311; 
               Compile`$313 = (-0.5917517095361372) $CellContext`phi1$22396 
                 Compile`$237; 
               Compile`$314 = 
                0.5917517095361372 $CellContext`phi2$22396 Compile`$237; 
               Compile`$315 = 0.908248290463863 Compile`$231 Compile`$237; 
               Compile`$316 = 1 + Compile`$313 + Compile`$314 + Compile`$315; 
               Compile`$324 = (-0.005279474328467426) $CellContext`phi1$22396 
                 Compile`$237; 
               Compile`$325 = (-0.7083647052177124) $CellContext`phi2$22396 
                 Compile`$237; 
               Compile`$326 = 0.935172654052801 Compile`$231 Compile`$237; 
               Compile`$327 = 1 + Compile`$324 + Compile`$325 + Compile`$326; 
               Compile`$345 = -$CellContext`phi1$22396 Compile`$231 
                 Compile`$228; 
               Compile`$347 = 1. $CellContext`phi1$22396 Compile`$237; 
               Compile`$413 = 
                1. $CellContext`phi1$22396 Compile`$231 Compile`$228; 
               Compile`$415 = -$CellContext`phi1$22396 Compile`$237; 
               Compile`$227 = Compile`$226^2; 
               Compile`$438 = 
                1.1835034190722744` $CellContext`phi1$22396 \
$CellContext`phi2$22396 Compile`$228; 
               Compile`$440 = 0.5917517095361372 Compile`$237; 
               Compile`$412 = (-3.125897657686229) $CellContext`phi1$22396 \
$CellContext`phi2$22396 Compile`$228; 
               Compile`$414 = 0.7464522479153887 Compile`$237; 
               Compile`$419 = 
                1.4929044958307773` $CellContext`phi1$22396 \
$CellContext`phi2$22396 Compile`$228; 
               Compile`$481 = -$CellContext`phi2$22396 Compile`$231 
                 Compile`$228; 
               Compile`$420 = (-1.5629488288431146`) Compile`$237; 
               Compile`$482 = 1. $CellContext`phi2$22396 Compile`$237; 
               Compile`$391 = (-2.8164965809277263`) $CellContext`phi1$22396 \
$CellContext`phi2$22396 Compile`$228; 
               Compile`$393 = (-1.4082482904638631`) Compile`$237; 
               Compile`$370 = 
                2.8164965809277263` $CellContext`phi1$22396 \
$CellContext`phi2$22396 Compile`$228; 
               Compile`$372 = 1.4082482904638631` Compile`$237; 
               Compile`$344 = (-1.4929044958307773`) $CellContext`phi1$22396 \
$CellContext`phi2$22396 Compile`$228; 
               Compile`$346 = 1.5629488288431146` Compile`$237; 
               Compile`$351 = 
                3.125897657686229 $CellContext`phi1$22396 \
$CellContext`phi2$22396 Compile`$228; 
               Compile`$541 = 
                1. $CellContext`phi2$22396 Compile`$231 Compile`$228; 
               Compile`$352 = (-0.7464522479153887) Compile`$237; 
               Compile`$542 = -$CellContext`phi2$22396 Compile`$237; 
               Compile`$319 = (-1.1835034190722744`) $CellContext`phi1$22396 \
$CellContext`phi2$22396 Compile`$228; 
               Compile`$321 = (-0.5917517095361372) 
                 Compile`$237; ((
                   Compile`$227 (
                    341.07848212762275` (
                    0.010558948656934852` Compile`$165 Compile`$228 + 
                    1.4167294104354249` $CellContext`phi1$22396 \
$CellContext`phi2$22396 Compile`$228 - 
                    1.870345308105602 $CellContext`phi1$22396 Compile`$231 
                    Compile`$228 - 0.005279474328467426 Compile`$237 + 
                    1.870345308105602 $CellContext`phi1$22396 Compile`$237) 
                    Compile`$244 Compile`$249 Compile`$253 Compile`$257 
                    Compile`$261 Compile`$264 Compile`$268 Compile`$272 
                    Compile`$276 Compile`$280 Compile`$284 Compile`$289 
                    Compile`$293 Compile`$297 Compile`$301 Compile`$304 
                    Compile`$308 Compile`$312 Compile`$316 + 
                    341.07848212762275` (
                    1.1835034190722744` Compile`$165 Compile`$228 + 
                    Compile`$319 - 1.816496580927726 $CellContext`phi1$22396 
                    Compile`$231 Compile`$228 + Compile`$321 + 
                    1.816496580927726 $CellContext`phi1$22396 Compile`$237) 
                    Compile`$244 Compile`$249 Compile`$253 Compile`$257 
                    Compile`$261 Compile`$264 Compile`$268 Compile`$272 
                    Compile`$276 Compile`$280 Compile`$284 Compile`$289 
                    Compile`$293 Compile`$297 Compile`$301 Compile`$304 
                    Compile`$308 Compile`$312 Compile`$327 + 
                    341.07848212762275` ((-2.6527973879352147`) Compile`$165 
                    Compile`$228 + 
                    1.225509028842855 $CellContext`phi1$22396 \
$CellContext`phi2$22396 Compile`$228 - 
                    1.3657226693941882` $CellContext`phi1$22396 Compile`$231 
                    Compile`$228 + 1.3263986939676073` Compile`$237 + 
                    1.3657226693941882` $CellContext`phi1$22396 Compile`$237) 
                    Compile`$244 Compile`$249 Compile`$253 Compile`$257 
                    Compile`$261 Compile`$264 Compile`$268 Compile`$272 
                    Compile`$276 Compile`$280 Compile`$284 Compile`$289 
                    Compile`$293 Compile`$297 Compile`$301 Compile`$304 
                    Compile`$308 Compile`$316 Compile`$327 + 
                    341.07848212762275` ((-0.7549333678869505) Compile`$165 
                    Compile`$228 - 
                    2.981756067963913 $CellContext`phi1$22396 \
$CellContext`phi2$22396 Compile`$228 - 
                    1.2785935985694648` $CellContext`phi1$22396 Compile`$231 
                    Compile`$228 + 0.37746668394347527` Compile`$237 + 
                    1.2785935985694648` $CellContext`phi1$22396 Compile`$237) 
                    Compile`$244 Compile`$249 Compile`$253 Compile`$257 
                    Compile`$261 Compile`$264 Compile`$268 Compile`$272 
                    Compile`$276 Compile`$280 Compile`$284 Compile`$289 
                    Compile`$293 Compile`$297 Compile`$301 Compile`$304 
                    Compile`$312 Compile`$316 Compile`$327 + 
                    341.07848212762275` ((-3.125897657686229) Compile`$165 
                    Compile`$228 + Compile`$344 + Compile`$345 + Compile`$346 + 
                    Compile`$347) Compile`$244 Compile`$249 Compile`$253 
                    Compile`$257 Compile`$261 Compile`$264 Compile`$268 
                    Compile`$272 Compile`$276 Compile`$280 Compile`$284 
                    Compile`$289 Compile`$293 Compile`$297 Compile`$301 
                    Compile`$308 Compile`$312 Compile`$316 Compile`$327 + 
                    341.07848212762275` (
                    1.4929044958307773` Compile`$165 Compile`$228 + 
                    Compile`$351 + Compile`$345 + Compile`$352 + Compile`$347)
                     Compile`$244 Compile`$249 Compile`$253 Compile`$257 
                    Compile`$261 Compile`$264 Compile`$268 Compile`$272 
                    Compile`$276 Compile`$280 Compile`$284 Compile`$289 
                    Compile`$293 Compile`$297 Compile`$304 Compile`$308 
                    Compile`$312 Compile`$316 Compile`$327 + 
                    341.07848212762275` (
                    3.3907685158790417` Compile`$165 Compile`$228 - 
                    1.0813674391205381` $CellContext`phi1$22396 \
$CellContext`phi2$22396 Compile`$228 - 
                    0.9128709291752768 $CellContext`phi1$22396 Compile`$231 
                    Compile`$228 - 1.6953842579395209` Compile`$237 + 
                    0.9128709291752768 $CellContext`phi1$22396 Compile`$237) 
                    Compile`$244 Compile`$249 Compile`$253 Compile`$257 
                    Compile`$261 Compile`$264 Compile`$268 Compile`$272 
                    Compile`$276 Compile`$280 Compile`$284 Compile`$289 
                    Compile`$293 Compile`$301 Compile`$304 Compile`$308 
                    Compile`$312 Compile`$316 Compile`$327 + 
                    341.07848212762275` (
                    3.581988897471612 Compile`$165 Compile`$228 + 
                    1.5819888974716114` $CellContext`phi1$22396 \
$CellContext`phi2$22396 Compile`$228 - 
                    0.4082482904638631 $CellContext`phi1$22396 Compile`$231 
                    Compile`$228 - 1.790994448735806 Compile`$237 + 
                    0.4082482904638631 $CellContext`phi1$22396 Compile`$237) 
                    Compile`$244 Compile`$249 Compile`$253 Compile`$257 
                    Compile`$261 Compile`$264 Compile`$268 Compile`$272 
                    Compile`$276 Compile`$280 Compile`$284 Compile`$289 
                    Compile`$297 Compile`$301 Compile`$304 Compile`$308 
                    Compile`$312 Compile`$316 Compile`$327 + 
                    341.07848212762275` ((-2.8164965809277263`) Compile`$165 
                    Compile`$228 + Compile`$370 - 
                    0.18350341907227397` $CellContext`phi1$22396 Compile`$231 
                    Compile`$228 + Compile`$372 + 
                    0.18350341907227397` $CellContext`phi1$22396 Compile`$237)
                     Compile`$244 Compile`$249 Compile`$253 Compile`$257 
                    Compile`$261 Compile`$264 Compile`$268 Compile`$272 
                    Compile`$276 Compile`$280 Compile`$284 Compile`$293 
                    Compile`$297 Compile`$301 Compile`$304 Compile`$308 
                    Compile`$312 Compile`$316 Compile`$327 + 
                    341.07848212762275` (
                    0.25431190953587707` Compile`$165 Compile`$228 - 
                    3.991001345386741 $CellContext`phi1$22396 \
$CellContext`phi2$22396 Compile`$228 - 
                    0.04252562106967517 $CellContext`phi1$22396 Compile`$231 
                    Compile`$228 - 0.12715595476793853` Compile`$237 + 
                    0.04252562106967517 $CellContext`phi1$22396 Compile`$237) 
                    Compile`$244 Compile`$249 Compile`$253 Compile`$257 
                    Compile`$261 Compile`$264 Compile`$268 Compile`$272 
                    Compile`$276 Compile`$280 Compile`$289 Compile`$293 
                    Compile`$297 Compile`$301 Compile`$304 Compile`$308 
                    Compile`$312 Compile`$316 Compile`$327 + 
                    341.07848212762275` ((-0.25431190953587707`) Compile`$165 
                    Compile`$228 + 
                    3.991001345386741 $CellContext`phi1$22396 \
$CellContext`phi2$22396 Compile`$228 + 
                    0.04252562106967517 $CellContext`phi1$22396 Compile`$231 
                    Compile`$228 + 0.12715595476793853` Compile`$237 - 
                    0.04252562106967517 $CellContext`phi1$22396 Compile`$237) 
                    Compile`$244 Compile`$249 Compile`$253 Compile`$257 
                    Compile`$261 Compile`$264 Compile`$268 Compile`$272 
                    Compile`$276 Compile`$284 Compile`$289 Compile`$293 
                    Compile`$297 Compile`$301 Compile`$304 Compile`$308 
                    Compile`$312 Compile`$316 Compile`$327 + 
                    341.07848212762275` (
                    2.8164965809277263` Compile`$165 Compile`$228 + 
                    Compile`$391 + 
                    0.18350341907227397` $CellContext`phi1$22396 Compile`$231 
                    Compile`$228 + Compile`$393 - 
                    0.18350341907227397` $CellContext`phi1$22396 Compile`$237)
                     Compile`$244 Compile`$249 Compile`$253 Compile`$257 
                    Compile`$261 Compile`$264 Compile`$268 Compile`$272 
                    Compile`$280 Compile`$284 Compile`$289 Compile`$293 
                    Compile`$297 Compile`$301 Compile`$304 Compile`$308 
                    Compile`$312 Compile`$316 Compile`$327 + 
                    341.07848212762275` ((-3.581988897471612) Compile`$165 
                    Compile`$228 - 
                    1.5819888974716114` $CellContext`phi1$22396 \
$CellContext`phi2$22396 Compile`$228 + 
                    0.4082482904638631 $CellContext`phi1$22396 Compile`$231 
                    Compile`$228 + 1.790994448735806 Compile`$237 - 
                    0.4082482904638631 $CellContext`phi1$22396 Compile`$237) 
                    Compile`$244 Compile`$249 Compile`$253 Compile`$257 
                    Compile`$261 Compile`$264 Compile`$268 Compile`$276 
                    Compile`$280 Compile`$284 Compile`$289 Compile`$293 
                    Compile`$297 Compile`$301 Compile`$304 Compile`$308 
                    Compile`$312 Compile`$316 Compile`$327 + 
                    341.07848212762275` ((-3.3907685158790417`) Compile`$165 
                    Compile`$228 + 
                    1.0813674391205381` $CellContext`phi1$22396 \
$CellContext`phi2$22396 Compile`$228 + 
                    0.9128709291752768 $CellContext`phi1$22396 Compile`$231 
                    Compile`$228 + 1.6953842579395209` Compile`$237 - 
                    0.9128709291752768 $CellContext`phi1$22396 Compile`$237) 
                    Compile`$244 Compile`$249 Compile`$253 Compile`$257 
                    Compile`$261 Compile`$264 Compile`$272 Compile`$276 
                    Compile`$280 Compile`$284 Compile`$289 Compile`$293 
                    Compile`$297 Compile`$301 Compile`$304 Compile`$308 
                    Compile`$312 Compile`$316 Compile`$327 + 
                    341.07848212762275` ((-1.4929044958307773`) Compile`$165 
                    Compile`$228 + Compile`$412 + Compile`$413 + Compile`$414 + 
                    Compile`$415) Compile`$244 Compile`$249 Compile`$253 
                    Compile`$257 Compile`$261 Compile`$268 Compile`$272 
                    Compile`$276 Compile`$280 Compile`$284 Compile`$289 
                    Compile`$293 Compile`$297 Compile`$301 Compile`$304 
                    Compile`$308 Compile`$312 Compile`$316 Compile`$327 + 
                    341.07848212762275` (
                    3.125897657686229 Compile`$165 Compile`$228 + 
                    Compile`$419 + Compile`$413 + Compile`$420 + Compile`$415)
                     Compile`$244 Compile`$249 Compile`$253 Compile`$257 
                    Compile`$264 Compile`$268 Compile`$272 Compile`$276 
                    Compile`$280 Compile`$284 Compile`$289 Compile`$293 
                    Compile`$297 Compile`$301 Compile`$304 Compile`$308 
                    Compile`$312 Compile`$316 Compile`$327 + 
                    341.07848212762275` (
                    0.7549333678869505 Compile`$165 Compile`$228 + 
                    2.981756067963913 $CellContext`phi1$22396 \
$CellContext`phi2$22396 Compile`$228 + 
                    1.2785935985694648` $CellContext`phi1$22396 Compile`$231 
                    Compile`$228 - 0.37746668394347527` Compile`$237 - 
                    1.2785935985694648` $CellContext`phi1$22396 Compile`$237) 
                    Compile`$244 Compile`$249 Compile`$253 Compile`$261 
                    Compile`$264 Compile`$268 Compile`$272 Compile`$276 
                    Compile`$280 Compile`$284 Compile`$289 Compile`$293 
                    Compile`$297 Compile`$301 Compile`$304 Compile`$308 
                    Compile`$312 Compile`$316 Compile`$327 + 
                    341.07848212762275` (
                    2.6527973879352147` Compile`$165 Compile`$228 - 
                    1.225509028842855 $CellContext`phi1$22396 \
$CellContext`phi2$22396 Compile`$228 + 
                    1.3657226693941882` $CellContext`phi1$22396 Compile`$231 
                    Compile`$228 - 1.3263986939676073` Compile`$237 - 
                    1.3657226693941882` $CellContext`phi1$22396 Compile`$237) 
                    Compile`$244 Compile`$249 Compile`$257 Compile`$261 
                    Compile`$264 Compile`$268 Compile`$272 Compile`$276 
                    Compile`$280 Compile`$284 Compile`$289 Compile`$293 
                    Compile`$297 Compile`$301 Compile`$304 Compile`$308 
                    Compile`$312 Compile`$316 Compile`$327 + 
                    341.07848212762275` ((-1.1835034190722744`) Compile`$165 
                    Compile`$228 + Compile`$438 + 
                    1.816496580927726 $CellContext`phi1$22396 Compile`$231 
                    Compile`$228 + Compile`$440 - 
                    1.816496580927726 $CellContext`phi1$22396 Compile`$237) 
                    Compile`$244 Compile`$253 Compile`$257 Compile`$261 
                    Compile`$264 Compile`$268 Compile`$272 Compile`$276 
                    Compile`$280 Compile`$284 Compile`$289 Compile`$293 
                    Compile`$297 Compile`$301 Compile`$304 Compile`$308 
                    Compile`$312 Compile`$316 Compile`$327 + 
                    341.07848212762275` ((-0.010558948656934852`) 
                    Compile`$165 Compile`$228 - 
                    1.4167294104354249` $CellContext`phi1$22396 \
$CellContext`phi2$22396 Compile`$228 + 
                    1.870345308105602 $CellContext`phi1$22396 Compile`$231 
                    Compile`$228 + 0.005279474328467426 Compile`$237 - 
                    1.870345308105602 $CellContext`phi1$22396 Compile`$237) 
                    Compile`$249 Compile`$253 Compile`$257 Compile`$261 
                    Compile`$264 Compile`$268 Compile`$272 Compile`$276 
                    Compile`$280 Compile`$284 Compile`$289 Compile`$293 
                    Compile`$297 Compile`$301 Compile`$304 Compile`$308 
                    Compile`$312 Compile`$316 Compile`$327)^2)/
                  4 + (Compile`$227 (
                    341.07848212762275` (
                    0.010558948656934852` $CellContext`phi1$22396 \
$CellContext`phi2$22396 Compile`$228 + 
                    1.4167294104354249` Compile`$167 Compile`$228 - 
                    1.870345308105602 $CellContext`phi2$22396 Compile`$231 
                    Compile`$228 - 0.7083647052177124 Compile`$237 + 
                    1.870345308105602 $CellContext`phi2$22396 Compile`$237) 
                    Compile`$244 Compile`$249 Compile`$253 Compile`$257 
                    Compile`$261 Compile`$264 Compile`$268 Compile`$272 
                    Compile`$276 Compile`$280 Compile`$284 Compile`$289 
                    Compile`$293 Compile`$297 Compile`$301 Compile`$304 
                    Compile`$308 Compile`$312 Compile`$316 + 
                    341.07848212762275` (Compile`$438 - 1.1835034190722744` 
                    Compile`$167 Compile`$228 - 
                    1.816496580927726 $CellContext`phi2$22396 Compile`$231 
                    Compile`$228 + Compile`$440 + 
                    1.816496580927726 $CellContext`phi2$22396 Compile`$237) 
                    Compile`$244 Compile`$249 Compile`$253 Compile`$257 
                    Compile`$261 Compile`$264 Compile`$268 Compile`$272 
                    Compile`$276 Compile`$280 Compile`$284 Compile`$289 
                    Compile`$293 Compile`$297 Compile`$301 Compile`$304 
                    Compile`$308 Compile`$312 Compile`$327 + 
                    341.07848212762275` ((-2.6527973879352147`) \
$CellContext`phi1$22396 $CellContext`phi2$22396 Compile`$228 + 
                    1.225509028842855 Compile`$167 Compile`$228 - 
                    1.3657226693941882` $CellContext`phi2$22396 Compile`$231 
                    Compile`$228 - 0.6127545144214275 Compile`$237 + 
                    1.3657226693941882` $CellContext`phi2$22396 Compile`$237) 
                    Compile`$244 Compile`$249 Compile`$253 Compile`$257 
                    Compile`$261 Compile`$264 Compile`$268 Compile`$272 
                    Compile`$276 Compile`$280 Compile`$284 Compile`$289 
                    Compile`$293 Compile`$297 Compile`$301 Compile`$304 
                    Compile`$308 Compile`$316 Compile`$327 + 
                    341.07848212762275` ((-0.7549333678869505) \
$CellContext`phi1$22396 $CellContext`phi2$22396 Compile`$228 - 
                    2.981756067963913 Compile`$167 Compile`$228 - 
                    1.2785935985694648` $CellContext`phi2$22396 Compile`$231 
                    Compile`$228 + 1.4908780339819565` Compile`$237 + 
                    1.2785935985694648` $CellContext`phi2$22396 Compile`$237) 
                    Compile`$244 Compile`$249 Compile`$253 Compile`$257 
                    Compile`$261 Compile`$264 Compile`$268 Compile`$272 
                    Compile`$276 Compile`$280 Compile`$284 Compile`$289 
                    Compile`$293 Compile`$297 Compile`$301 Compile`$304 
                    Compile`$312 Compile`$316 Compile`$327 + 
                    341.07848212762275` (Compile`$412 - 1.4929044958307773` 
                    Compile`$167 Compile`$228 + Compile`$481 + Compile`$414 + 
                    Compile`$482) Compile`$244 Compile`$249 Compile`$253 
                    Compile`$257 Compile`$261 Compile`$264 Compile`$268 
                    Compile`$272 Compile`$276 Compile`$280 Compile`$284 
                    Compile`$289 Compile`$293 Compile`$297 Compile`$301 
                    Compile`$308 Compile`$312 Compile`$316 Compile`$327 + 
                    341.07848212762275` (Compile`$419 + 
                    3.125897657686229 Compile`$167 Compile`$228 + 
                    Compile`$481 + Compile`$420 + Compile`$482) Compile`$244 
                    Compile`$249 Compile`$253 Compile`$257 Compile`$261 
                    Compile`$264 Compile`$268 Compile`$272 Compile`$276 
                    Compile`$280 Compile`$284 Compile`$289 Compile`$293 
                    Compile`$297 Compile`$304 Compile`$308 Compile`$312 
                    Compile`$316 Compile`$327 + 
                    341.07848212762275` (
                    3.3907685158790417` $CellContext`phi1$22396 \
$CellContext`phi2$22396 Compile`$228 - 1.0813674391205381` Compile`$167 
                    Compile`$228 - 0.9128709291752768 $CellContext`phi2$22396 
                    Compile`$231 Compile`$228 + 
                    0.5406837195602691 Compile`$237 + 
                    0.9128709291752768 $CellContext`phi2$22396 Compile`$237) 
                    Compile`$244 Compile`$249 Compile`$253 Compile`$257 
                    Compile`$261 Compile`$264 Compile`$268 Compile`$272 
                    Compile`$276 Compile`$280 Compile`$284 Compile`$289 
                    Compile`$293 Compile`$301 Compile`$304 Compile`$308 
                    Compile`$312 Compile`$316 Compile`$327 + 
                    341.07848212762275` (
                    3.581988897471612 $CellContext`phi1$22396 \
$CellContext`phi2$22396 Compile`$228 + 
                    1.5819888974716114` Compile`$167 Compile`$228 - 
                    0.4082482904638631 $CellContext`phi2$22396 Compile`$231 
                    Compile`$228 - 0.7909944487358057 Compile`$237 + 
                    0.4082482904638631 $CellContext`phi2$22396 Compile`$237) 
                    Compile`$244 Compile`$249 Compile`$253 Compile`$257 
                    Compile`$261 Compile`$264 Compile`$268 Compile`$272 
                    Compile`$276 Compile`$280 Compile`$284 Compile`$289 
                    Compile`$297 Compile`$301 Compile`$304 Compile`$308 
                    Compile`$312 Compile`$316 Compile`$327 + 
                    341.07848212762275` (Compile`$391 + 
                    2.8164965809277263` Compile`$167 Compile`$228 - 
                    0.18350341907227397` $CellContext`phi2$22396 Compile`$231 
                    Compile`$228 + Compile`$393 + 
                    0.18350341907227397` $CellContext`phi2$22396 Compile`$237)
                     Compile`$244 Compile`$249 Compile`$253 Compile`$257 
                    Compile`$261 Compile`$264 Compile`$268 Compile`$272 
                    Compile`$276 Compile`$280 Compile`$284 Compile`$293 
                    Compile`$297 Compile`$301 Compile`$304 Compile`$308 
                    Compile`$312 Compile`$316 Compile`$327 + 
                    341.07848212762275` (
                    0.25431190953587707` $CellContext`phi1$22396 \
$CellContext`phi2$22396 Compile`$228 - 3.991001345386741 Compile`$167 
                    Compile`$228 - 
                    0.04252562106967517 $CellContext`phi2$22396 Compile`$231 
                    Compile`$228 + 1.9955006726933704` Compile`$237 + 
                    0.04252562106967517 $CellContext`phi2$22396 Compile`$237) 
                    Compile`$244 Compile`$249 Compile`$253 Compile`$257 
                    Compile`$261 Compile`$264 Compile`$268 Compile`$272 
                    Compile`$276 Compile`$280 Compile`$289 Compile`$293 
                    Compile`$297 Compile`$301 Compile`$304 Compile`$308 
                    Compile`$312 Compile`$316 Compile`$327 + 
                    341.07848212762275` ((-0.25431190953587707`) \
$CellContext`phi1$22396 $CellContext`phi2$22396 Compile`$228 + 
                    3.991001345386741 Compile`$167 Compile`$228 + 
                    0.04252562106967517 $CellContext`phi2$22396 Compile`$231 
                    Compile`$228 - 1.9955006726933704` Compile`$237 - 
                    0.04252562106967517 $CellContext`phi2$22396 Compile`$237) 
                    Compile`$244 Compile`$249 Compile`$253 Compile`$257 
                    Compile`$261 Compile`$264 Compile`$268 Compile`$272 
                    Compile`$276 Compile`$284 Compile`$289 Compile`$293 
                    Compile`$297 Compile`$301 Compile`$304 Compile`$308 
                    Compile`$312 Compile`$316 Compile`$327 + 
                    341.07848212762275` (Compile`$370 - 2.8164965809277263` 
                    Compile`$167 Compile`$228 + 
                    0.18350341907227397` $CellContext`phi2$22396 Compile`$231 
                    Compile`$228 + Compile`$372 - 
                    0.18350341907227397` $CellContext`phi2$22396 Compile`$237)
                     Compile`$244 Compile`$249 Compile`$253 Compile`$257 
                    Compile`$261 Compile`$264 Compile`$268 Compile`$272 
                    Compile`$280 Compile`$284 Compile`$289 Compile`$293 
                    Compile`$297 Compile`$301 Compile`$304 Compile`$308 
                    Compile`$312 Compile`$316 Compile`$327 + 
                    341.07848212762275` ((-3.581988897471612) \
$CellContext`phi1$22396 $CellContext`phi2$22396 Compile`$228 - 
                    1.5819888974716114` Compile`$167 Compile`$228 + 
                    0.4082482904638631 $CellContext`phi2$22396 Compile`$231 
                    Compile`$228 + 0.7909944487358057 Compile`$237 - 
                    0.4082482904638631 $CellContext`phi2$22396 Compile`$237) 
                    Compile`$244 Compile`$249 Compile`$253 Compile`$257 
                    Compile`$261 Compile`$264 Compile`$268 Compile`$276 
                    Compile`$280 Compile`$284 Compile`$289 Compile`$293 
                    Compile`$297 Compile`$301 Compile`$304 Compile`$308 
                    Compile`$312 Compile`$316 Compile`$327 + 
                    341.07848212762275` ((-3.3907685158790417`) \
$CellContext`phi1$22396 $CellContext`phi2$22396 Compile`$228 + 
                    1.0813674391205381` Compile`$167 Compile`$228 + 
                    0.9128709291752768 $CellContext`phi2$22396 Compile`$231 
                    Compile`$228 - 0.5406837195602691 Compile`$237 - 
                    0.9128709291752768 $CellContext`phi2$22396 Compile`$237) 
                    Compile`$244 Compile`$249 Compile`$253 Compile`$257 
                    Compile`$261 Compile`$264 Compile`$272 Compile`$276 
                    Compile`$280 Compile`$284 Compile`$289 Compile`$293 
                    Compile`$297 Compile`$301 Compile`$304 Compile`$308 
                    Compile`$312 Compile`$316 Compile`$327 + 
                    341.07848212762275` (Compile`$344 - 3.125897657686229 
                    Compile`$167 Compile`$228 + Compile`$541 + Compile`$346 + 
                    Compile`$542) Compile`$244 Compile`$249 Compile`$253 
                    Compile`$257 Compile`$261 Compile`$268 Compile`$272 
                    Compile`$276 Compile`$280 Compile`$284 Compile`$289 
                    Compile`$293 Compile`$297 Compile`$301 Compile`$304 
                    Compile`$308 Compile`$312 Compile`$316 Compile`$327 + 
                    341.07848212762275` (Compile`$351 + 
                    1.4929044958307773` Compile`$167 Compile`$228 + 
                    Compile`$541 + Compile`$352 + Compile`$542) Compile`$244 
                    Compile`$249 Compile`$253 Compile`$257 Compile`$264 
                    Compile`$268 Compile`$272 Compile`$276 Compile`$280 
                    Compile`$284 Compile`$289 Compile`$293 Compile`$297 
                    Compile`$301 Compile`$304 Compile`$308 Compile`$312 
                    Compile`$316 Compile`$327 + 
                    341.07848212762275` (
                    0.7549333678869505 $CellContext`phi1$22396 \
$CellContext`phi2$22396 Compile`$228 + 
                    2.981756067963913 Compile`$167 Compile`$228 + 
                    1.2785935985694648` $CellContext`phi2$22396 Compile`$231 
                    Compile`$228 - 1.4908780339819565` Compile`$237 - 
                    1.2785935985694648` $CellContext`phi2$22396 Compile`$237) 
                    Compile`$244 Compile`$249 Compile`$253 Compile`$261 
                    Compile`$264 Compile`$268 Compile`$272 Compile`$276 
                    Compile`$280 Compile`$284 Compile`$289 Compile`$293 
                    Compile`$297 Compile`$301 Compile`$304 Compile`$308 
                    Compile`$312 Compile`$316 Compile`$327 + 
                    341.07848212762275` (
                    2.6527973879352147` $CellContext`phi1$22396 \
$CellContext`phi2$22396 Compile`$228 - 1.225509028842855 Compile`$167 
                    Compile`$228 + 
                    1.3657226693941882` $CellContext`phi2$22396 Compile`$231 
                    Compile`$228 + 0.6127545144214275 Compile`$237 - 
                    1.3657226693941882` $CellContext`phi2$22396 Compile`$237) 
                    Compile`$244 Compile`$249 Compile`$257 Compile`$261 
                    Compile`$264 Compile`$268 Compile`$272 Compile`$276 
                    Compile`$280 Compile`$284 Compile`$289 Compile`$293 
                    Compile`$297 Compile`$301 Compile`$304 Compile`$308 
                    Compile`$312 Compile`$316 Compile`$327 + 
                    341.07848212762275` (Compile`$319 + 
                    1.1835034190722744` Compile`$167 Compile`$228 + 
                    1.816496580927726 $CellContext`phi2$22396 Compile`$231 
                    Compile`$228 + Compile`$321 - 
                    1.816496580927726 $CellContext`phi2$22396 Compile`$237) 
                    Compile`$244 Compile`$253 Compile`$257 Compile`$261 
                    Compile`$264 Compile`$268 Compile`$272 Compile`$276 
                    Compile`$280 Compile`$284 Compile`$289 Compile`$293 
                    Compile`$297 Compile`$301 Compile`$304 Compile`$308 
                    Compile`$312 Compile`$316 Compile`$327 + 
                    341.07848212762275` ((-0.010558948656934852`) \
$CellContext`phi1$22396 $CellContext`phi2$22396 Compile`$228 - 
                    1.4167294104354249` Compile`$167 Compile`$228 + 
                    1.870345308105602 $CellContext`phi2$22396 Compile`$231 
                    Compile`$228 + 0.7083647052177124 Compile`$237 - 
                    1.870345308105602 $CellContext`phi2$22396 Compile`$237) 
                    Compile`$249 Compile`$253 Compile`$257 Compile`$261 
                    Compile`$264 Compile`$268 Compile`$272 Compile`$276 
                    Compile`$280 Compile`$284 Compile`$289 Compile`$293 
                    Compile`$297 Compile`$301 Compile`$304 Compile`$308 
                    Compile`$312 Compile`$316 Compile`$327)^2)/4)/2]], 
             Evaluate], "Umax" -> 1.938|>|>, 
       Attributes[$CellContext`phi1$22384] = {Temporary}, 
       Attributes[$CellContext`phi2$22384] = {Temporary}, 
       Attributes[$CellContext`phi1$22387] = {Temporary}, 
       Attributes[$CellContext`phi2$22387] = {Temporary}, 
       Attributes[$CellContext`phi1$22390] = {Temporary}, 
       Attributes[$CellContext`phi2$22390] = {Temporary}, 
       Attributes[$CellContext`phi1$22393] = {Temporary}, 
       Attributes[$CellContext`phi2$22393] = {Temporary}, 
       Attributes[$CellContext`phi1$22396] = {Temporary}, 
       Attributes[$CellContext`phi2$22396] = {
        Temporary}, $CellContext`vMaxFixed = <|
        "Tetrahedron" -> 1.5690432921983544`, "Octahedron" -> 
         3.586920069617875, "Cube" -> 5.816698900306769, "Icosahedron" -> 
         6.156279288428893, "Dodecahedron" -> 10.625843871752878`|>}; 
     Typeset`initDone$$ = True),
    SynchronousInitialization->True,
    UndoTrackedVariables:>{Typeset`show$$, Typeset`bookmarkMode$$},
    UnsavedVariables:>{Typeset`initDone$$},
    UntrackedVariables:>{Typeset`size$$}], "Manipulate",
   Deployed->True,
   StripOnInput->False],
  Manipulate`InterpretManipulate[1]]], "Output",
 CellChangeTimes->{
  3.9942549823800097`*^9, {3.994255105067416*^9, 3.994255121635496*^9}, 
   3.994255292938599*^9},
 CellLabel->"Out[49]=",ExpressionUUID->"b7648be2-3a44-4dfa-820b-b8662ed8a4bf"]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell["", "Section",
 CellChangeTimes->{{3.994259022934743*^9, 
  3.9942590238058567`*^9}},ExpressionUUID->"e2f56fb6-901c-44e0-84f7-\
f3ece20bfe96"],

Cell["", "Text",
 CellChangeTimes->{{3.9942590201442537`*^9, 
  3.994259020494626*^9}},ExpressionUUID->"3e656aeb-f5cb-4b47-ae86-\
9cfb643cdbe6"]
}, Open  ]]
}, Open  ]]
},
WindowSize->{960, 546},
WindowMargins->{{0, Automatic}, {Automatic, 0}},
WindowTitle->"Platonic Solids -- Interactive Viewer",
FrontEndVersion->"14.3 for Linux x86 (64-bit) (July 8, 2025)",
StyleDefinitions->"Default.nb",
ExpressionUUID->"20a85d55-444e-436a-ab5a-3f6cd00947cd"
]
(* End of Notebook Content *)

(* Internal cache information *)
(*CellTagsOutline
CellTagsIndex->{}
*)
(*CellTagsIndex
CellTagsIndex->{}
*)
(*NotebookFileOutline
Notebook[{
Cell[CellGroupData[{
Cell[1557, 36, 130, 0, 96, "Title",ExpressionUUID->"555a3f2a-a95e-411c-a6e2-0b9b2c90c275"],
Cell[1690, 38, 186, 3, 52, "Subtitle",ExpressionUUID->"434989e0-36bc-463a-86a0-6b89300c1aa1"],
Cell[CellGroupData[{
Cell[1901, 45, 103, 0, 65, "Section",ExpressionUUID->"2e426786-0546-4dd8-9700-c5976fa3b105"],
Cell[CellGroupData[{
Cell[2029, 49, 3203, 92, 1408, "Input",ExpressionUUID->"f89fdac5-27df-4551-8df8-980dd68db759"],
Cell[CellGroupData[{
Cell[5257, 145, 238, 5, 23, "Print",ExpressionUUID->"39f1c8ed-aed1-42a7-9610-98bd6745a5cc"],
Cell[5498, 152, 992, 20, 41, "Print",ExpressionUUID->"deaa8895-7a87-4ab6-9b84-849b8228df8b"]
}, Open  ]]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[6551, 179, 95, 0, 65, "Section",ExpressionUUID->"562eeede-2301-442c-8f88-b77a8055cde6"],
Cell[CellGroupData[{
Cell[6671, 183, 1546, 35, 437, "Input",ExpressionUUID->"b27cfa94-ad6b-4664-a3c4-0f7566674086"],
Cell[8220, 220, 166683, 2519, 574, "Output",ExpressionUUID->"b7648be2-3a44-4dfa-820b-b8662ed8a4bf"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[174952, 2745, 147, 3, 65, "Section",ExpressionUUID->"e2f56fb6-901c-44e0-84f7-f3ece20bfe96"],
Cell[175102, 2750, 144, 3, 33, "Text",ExpressionUUID->"3e656aeb-f5cb-4b47-ae86-9cfb643cdbe6"]
}, Open  ]]
}, Open  ]]
}
]
*)

(* End of internal cache information *)

(* NotebookSignature vuDZ8MTGD4lvSBK#UxXSp6f# *)
