/*
 * ThermalNewtonianIncompressibleMultigrid3D.java
 */

import com.comsol.model.*;
import com.comsol.model.util.*;

/** Model exported on Dec 3 2024, 15:27 by COMSOL 6.1.0.357. */
public class ThermalNewtonianIncompressibleMultigrid3D {

  public static Model run() {
    Model model = ModelUtil.create("Model");

    model.modelPath("/work/sm1100/SLS");

    model.label("ThermalNewtonianIncompressible3D.mph");

    model.param().set("MaterialQuantities", "0", "___________________________________________________");
    model.param().set("rhomat", "1000 [kg/m^3]", "Material density");
    model.param().set("rhoair", "1 [kg/m^3]", "Air density");
    model.param().set("mumatsolid", "1e6 [Pa*s]", "Material viscosity in solid phase");
    model.param().set("mumatliq", "1 [Pa*s]", "Material viscosity in liquid phase");
    model.param().set("muair", "1 [Pa*s]", "Air viscosity");
    model.param().set("TemperatureQuantities", "0", "___________________________________________________");
    model.param().set("kair", "0.03 [W/(m*K)]");
    model.param().set("kmat", "0.6 [W/(m*K)]");
    model.param().set("cpair", "1 [J/(kg*K)]");
    model.param().set("cpmatsolid", "2098 [J/(kg*K)]");
    model.param().set("cpmatliq", "1484 [J/(kg*K)]");
    model.param().set("Hmat", "4000");
    model.param().set("Q", "1e8 [W/m^3]");
    model.param().set("deltat", "3 [K]");
    model.param().set("T0", "310[K]");
    model.param().set("Tm", "320 [K]");
    model.param().set("epsilon", "5e-6", "interface thickness");
    model.param().set("g", "9.81", "Gravitational acceleration");
    model.param().set("sigma", "0.05", "Surface Tension");

    model.component().create("comp1", true);

    model.component("comp1").geom().create("geom1", 3);

    model.result().table().create("evl3", "Table");

    model.component("comp1").func().create("step1", "Step");
    model.component("comp1").func().create("gp1", "GaussianPulse");
    model.component("comp1").func().create("step2", "Step");
    model.component("comp1").func().create("an1", "Analytic");
    model.component("comp1").func("step1").label("Viscosity");
    model.component("comp1").func("step1").set("funcname", "mumat");
    model.component("comp1").func("step1").set("location", "Tm");
    model.component("comp1").func("step1").set("from", "mumatsolid");
    model.component("comp1").func("step1").set("to", "mumatliq");
    model.component("comp1").func("step1").set("smooth", "deltat");
    model.component("comp1").func("gp1").label("Hjump");
    model.component("comp1").func("gp1").set("funcname", "Hjump");
    model.component("comp1").func("gp1").set("location", "Tm");
    model.component("comp1").func("gp1").set("sigma", "deltat/2");
    model.component("comp1").func("gp1").set("normalization", "peak");
    model.component("comp1").func("gp1").set("peakvalue", "Hmat");
    model.component("comp1").func("step2").label("Hstep");
    model.component("comp1").func("step2").set("funcname", "Hstep");
    model.component("comp1").func("step2").set("location", "Tm");
    model.component("comp1").func("step2").set("from", "cpmatsolid");
    model.component("comp1").func("step2").set("to", "cpmatliq");
    model.component("comp1").func("an1").label("Specific Heat");
    model.component("comp1").func("an1").set("funcname", "Cpmat");
    model.component("comp1").func("an1").set("expr", "Hjump(T) + Hstep(T)");
    model.component("comp1").func("an1").set("args", new String[]{"T"});
    model.component("comp1").func("an1").set("argunit", new String[]{""});
    model.component("comp1").func("an1").set("plotargs", new String[][]{{"T", "300", "400"}});

    model.component("comp1").mesh().create("mesh1");

    model.component("comp1").geom("geom1").lengthUnit("\u00b5m");
    model.component("comp1").geom("geom1").create("blk1", "Block");
    model.component("comp1").geom("geom1").feature("blk1").set("size", new int[]{100, 100, 100});
    model.component("comp1").geom("geom1").create("sph1", "Sphere");
    model.component("comp1").geom("geom1").feature("sph1").set("pos", new int[]{50, 50, 10});
    model.component("comp1").geom("geom1").feature("sph1").set("r", 10);
    model.component("comp1").geom("geom1").run();

    model.component("comp1").physics().create("spf", "LaminarFlow", "geom1");
    model.component("comp1").physics("spf").create("out1", "OutletBoundary", 2);
    model.component("comp1").physics("spf").feature("out1").selection().set(4);
    model.component("comp1").physics("spf").create("constr1", "PointwiseConstraint", 2);
    model.component("comp1").physics("spf").feature("constr1").selection().set(4);
    model.component("comp1").physics().create("pf", "PhaseField", "geom1");
    model.component("comp1").physics("pf").feature("initfluid2").selection().set(2);
    model.component("comp1").physics("pf").create("out1", "Outlet", 2);
    model.component("comp1").physics("pf").feature("out1").selection().set(4);
    model.component("comp1").physics().create("ht", "HeatTransferInFluids", "geom1");
    model.component("comp1").physics("ht").create("ofl1", "ConvectiveOutflow", 2);
    model.component("comp1").physics("ht").feature("ofl1").selection().set(4);
    model.component("comp1").physics("ht").create("hs1", "HeatSource", 3);
    model.component("comp1").physics("ht").feature("hs1").selection().set(2);
    model.component("comp1").physics().create("viscode", "DomainODE", "geom1");
    model.component("comp1").physics("viscode").identifier("viscode");
    model.component("comp1").physics("viscode").field("dimensionless").field("mufield");
    model.component("comp1").physics("viscode").field("dimensionless").component(new String[]{"mu2"});
    model.component("comp1").physics("viscode").prop("Units").set("DependentVariableQuantity", "none");
    model.component("comp1").physics("viscode").prop("Units").set("CustomDependentVariableUnit", "Pa*s");
    model.component("comp1").physics("viscode").create("aleq1", "AlgebraicEquation", 3);
    model.component("comp1").physics("viscode").feature("aleq1").selection().all();
    model.component("comp1").physics().create("cpode", "DomainODE", "geom1");
    model.component("comp1").physics("cpode").identifier("cpode");
    model.component("comp1").physics("cpode").field("dimensionless").field("cpfield");
    model.component("comp1").physics("cpode").field("dimensionless").component(new String[]{"Cp2"});
    model.component("comp1").physics("cpode").prop("Units").set("DependentVariableQuantity", "none");
    model.component("comp1").physics("cpode").prop("Units").set("CustomDependentVariableUnit", "J/(kg*K)");
    model.component("comp1").physics("cpode").create("aleq1", "AlgebraicEquation", 3);
    model.component("comp1").physics("cpode").feature("aleq1").selection().all();

    model.component("comp1").multiphysics().create("tpf1", "TwoPhaseFlowPhaseField", 3);
    model.component("comp1").multiphysics("tpf1").selection().all();

    model.component("comp1").mesh("mesh1").create("size1", "Size");
    model.component("comp1").mesh("mesh1").create("cr1", "CornerRefinement");
    model.component("comp1").mesh("mesh1").create("ftet1", "FreeTet");
    model.component("comp1").mesh("mesh1").create("bl1", "BndLayer");
    model.component("comp1").mesh("mesh1").feature("size1").selection().geom("geom1", 2);
    model.component("comp1").mesh("mesh1").feature("size1").selection().set(1, 2, 3, 5, 14);
    model.component("comp1").mesh("mesh1").feature("cr1").selection().geom("geom1", 3);
    model.component("comp1").mesh("mesh1").feature("cr1").selection().set(1, 2);
    model.component("comp1").mesh("mesh1").feature("ftet1").create("dis1", "Distribution");
    model.component("comp1").mesh("mesh1").feature("ftet1").feature("dis1").selection()
         .set(1, 2, 3, 4, 5, 6, 7, 8, 21, 22, 23, 24);
    model.component("comp1").mesh("mesh1").feature("bl1").selection().geom("geom1", 3);
    model.component("comp1").mesh("mesh1").feature("bl1").selection().set(1, 2);
    model.component("comp1").mesh("mesh1").feature("bl1").create("blp1", "BndLayerProp");
    model.component("comp1").mesh("mesh1").feature("bl1").feature("blp1").selection().set(1, 2, 3, 5, 14);

    model.result().table("evl3").label("Evaluation 3D");
    model.result().table("evl3").comments("Interactive 3D values");

    model.component("comp1").physics("spf").prop("PhysicalModelProperty").set("IncludeGravity", true);
    model.component("comp1").physics("spf").feature("init1").set("p_init", "1 [atm]");
    model.component("comp1").physics("spf").feature("out1").set("p0", "1[atm]");
    model.component("comp1").physics("spf").feature("out1").set("SuppressBackflow", false);
    model.component("comp1").physics("spf").feature("constr1").set("constraintExpression", "p-(1[atm])");
    model.component("comp1").physics("pf").prop("ShapeProperty").set("order_phasefield", 2);
    model.component("comp1").physics("pf").feature("pfm1").set("epsilon_pf", "epsilon");
    model.component("comp1").physics("pf").feature("pfm1").set("chiOption", "velocity");
    model.component("comp1").physics("pf").feature("pfm1").set("chi", "sqrt(u^2 + v^2 + w^2)");
    model.component("comp1").physics("pf").feature("pfm1").set("U", "sqrt(u^2 + v^2 + w^2)");
    model.component("comp1").physics("ht").feature("fluid1").set("Cp_mat", "userdef");
    model.component("comp1").physics("ht").feature("fluid1").set("Cp", "cpair*(1-phipf)/2 + Cp2*(1+phipf)/2");
    model.component("comp1").physics("ht").feature("fluid1").set("rho_mat", "userdef");
    model.component("comp1").physics("ht").feature("fluid1").set("rho", "rhoair*(1-phipf)/2 + rhomat*(1+phipf)/2");
    model.component("comp1").physics("ht").feature("fluid1").set("k_mat", "userdef");
    model.component("comp1").physics("ht").feature("fluid1")
         .set("k", new String[][]{{"kair*(1-phipf)/2 + kmat*(1+phipf)/2"}, {"0"}, {"0"}, {"0"}, {"kair*(1-phipf)/2 + kmat*(1+phipf)/2"}, {"0"}, {"0"}, {"0"}, {"kair*(1-phipf)/2 + kmat*(1+phipf)/2"}});
    model.component("comp1").physics("ht").feature("fluid1").set("u", new String[][]{{"u"}, {"v"}, {"w"}});
    model.component("comp1").physics("ht").feature("init1").set("Tinit", "T0");
    model.component("comp1").physics("ht").feature("hs1").set("Q0", "Q");
    model.component("comp1").physics("ht").feature("hs1").set("materialType", "nonSolid");
    model.component("comp1").physics("viscode").label("Viscosity");
    model.component("comp1").physics("viscode").prop("EquationForm").set("form", "Automatic");
    model.component("comp1").physics("viscode").prop("Units").set("CustomSourceTermUnit", "Pa*s");
    model.component("comp1").physics("viscode").feature("init1").set("mu2", "mumat(T0)");
    model.component("comp1").physics("viscode").feature("init1").set("mu2t", "d(mumat(T0), T0) * d(T, t)");
    model.component("comp1").physics("viscode").feature("aleq1").set("f", "mu2 - mumat(T)");
    model.component("comp1").physics("cpode").label("Specific Heat");
    model.component("comp1").physics("cpode").prop("EquationForm").set("form", "Automatic");
    model.component("comp1").physics("cpode").prop("Units").set("CustomSourceTermUnit", "J/(kg*K)");
    model.component("comp1").physics("cpode").feature("init1").set("Cp2", "Cpmat(T0)");
    model.component("comp1").physics("cpode").feature("init1").set("Cp2t", "d(Cpmat(T0), T0) * d(T, t)");
    model.component("comp1").physics("cpode").feature("aleq1").set("f", "Cp2 - Cpmat(T)");

    model.component("comp1").multiphysics("tpf1").set("rho1_mat", "userdef");
    model.component("comp1").multiphysics("tpf1").set("rho1", "rhoair");
    model.component("comp1").multiphysics("tpf1").set("mu1_mat", "userdef");
    model.component("comp1").multiphysics("tpf1").set("mu1", "muair");
    model.component("comp1").multiphysics("tpf1").set("rho2_mat", "userdef");
    model.component("comp1").multiphysics("tpf1").set("rho2", "rhomat");
    model.component("comp1").multiphysics("tpf1").set("mu2_mat", "userdef");
    model.component("comp1").multiphysics("tpf1").set("mu2", "mu2");
    model.component("comp1").multiphysics("tpf1").set("SurfaceTensionCoefficient", "userdef");
    model.component("comp1").multiphysics("tpf1").set("sigma", "sigma");

    model.component("comp1").mesh("mesh1").feature("size").set("hauto", 3);
    model.component("comp1").mesh("mesh1").feature("size").set("table", "cfd");
    model.component("comp1").mesh("mesh1").feature("size1").set("hauto", 3);
    model.component("comp1").mesh("mesh1").feature("size1").set("table", "cfd");
    model.component("comp1").mesh("mesh1").feature("cr1").selection("boundary").set(1, 2, 3, 5, 14);
    model.component("comp1").mesh("mesh1").feature("ftet1").feature("dis1").set("numelem", 200);
    model.component("comp1").mesh("mesh1").feature("bl1").set("sharpcorners", "trim");
    model.component("comp1").mesh("mesh1").feature("bl1").feature("blp1").set("blnlayers", 2);
    model.component("comp1").mesh("mesh1").feature("bl1").feature("blp1").set("blhminfact", 5);
    model.component("comp1").mesh("mesh1").run();

    model.study().create("std1");
    model.study("std1").create("phasei", "PhaseInitialization");
    model.study("std1").create("time", "Transient");
    model.study("std1").feature("phasei")
         .set("activate", new String[]{"spf", "off", "pf", "on", "ht", "on", "viscode", "on", "cpode", "on", 
         "frame:spatial1", "on", "frame:material1", "on"});

    model.sol().create("sol4");
    model.sol("sol4").study("std1");
    model.sol("sol4").create("st1", "StudyStep");
    model.sol("sol4").create("v1", "Variables");
    model.sol("sol4").create("s1", "Stationary");
    model.sol("sol4").create("su1", "StoreSolution");
    model.sol("sol4").create("st2", "StudyStep");
    model.sol("sol4").create("v2", "Variables");
    model.sol("sol4").create("t1", "Time");
    model.sol("sol4").feature("s1").create("fc1", "FullyCoupled");
    model.sol("sol4").feature("s1").create("i1", "Iterative");
    model.sol("sol4").feature("s1").create("d1", "Direct");
    model.sol("sol4").feature("s1").feature("i1").create("mg1", "Multigrid");
    model.sol("sol4").feature("s1").feature("i1").feature("mg1").feature("pr").create("sl1", "SORLine");
    model.sol("sol4").feature("s1").feature("i1").feature("mg1").feature("po").create("sl1", "SORLine");
    model.sol("sol4").feature("s1").feature("i1").feature("mg1").feature("cs").create("d1", "Direct");
    model.sol("sol4").feature("s1").feature().remove("fcDef");
    model.sol("sol4").feature("t1").create("se1", "Segregated");
    model.sol("sol4").feature("t1").create("d1", "Direct");
    model.sol("sol4").feature("t1").create("i1", "Iterative");
    model.sol("sol4").feature("t1").create("i2", "Iterative");
    model.sol("sol4").feature("t1").create("i3", "Iterative");
    model.sol("sol4").feature("t1").create("d2", "Direct");
    model.sol("sol4").feature("t1").create("d3", "Direct");
    model.sol("sol4").feature("t1").feature("se1").create("ss1", "SegregatedStep");
    model.sol("sol4").feature("t1").feature("se1").create("ss2", "SegregatedStep");
    model.sol("sol4").feature("t1").feature("se1").create("ss3", "SegregatedStep");
    model.sol("sol4").feature("t1").feature("se1").create("ss4", "SegregatedStep");
    model.sol("sol4").feature("t1").feature("se1").create("ss5", "SegregatedStep");
    model.sol("sol4").feature("t1").feature("se1").create("ll1", "LowerLimit");
    model.sol("sol4").feature("t1").feature("se1").feature().remove("ssDef");
    model.sol("sol4").feature("t1").feature("i1").create("mg1", "Multigrid");
    model.sol("sol4").feature("t1").feature("i1").feature("mg1").feature("pr").create("sc1", "SCGS");
    model.sol("sol4").feature("t1").feature("i1").feature("mg1").feature("po").create("sc1", "SCGS");
    model.sol("sol4").feature("t1").feature("i1").feature("mg1").feature("cs").create("d1", "Direct");
    model.sol("sol4").feature("t1").feature("i2").create("mg1", "Multigrid");
    model.sol("sol4").feature("t1").feature("i2").feature("mg1").feature("pr").create("sc1", "SCGS");
    model.sol("sol4").feature("t1").feature("i2").feature("mg1").feature("po").create("sc1", "SCGS");
    model.sol("sol4").feature("t1").feature("i2").feature("mg1").feature("cs").create("d1", "Direct");
    model.sol("sol4").feature("t1").feature("i3").create("mg1", "Multigrid");
    model.sol("sol4").feature("t1").feature("i3").feature("mg1").feature("pr").create("so1", "SOR");
    model.sol("sol4").feature("t1").feature("i3").feature("mg1").feature("po").create("so1", "SOR");
    model.sol("sol4").feature("t1").feature("i3").feature("mg1").feature("cs").create("d1", "Direct");
    model.sol("sol4").feature("t1").feature().remove("fcDef");

    model.result().dataset().create("dset3", "Solution");
    model.result().dataset("dset2").set("solution", "sol4");
    model.result().dataset("dset3").set("solution", "sol5");
    model.result().dataset().remove("dset1");

    model.nodeGroup().create("grp1", "Definitions", "comp1");
    model.nodeGroup("grp1").set("type", "func");
    model.nodeGroup("grp1").placeAfter(null);

    model.study("std1").feature("time").set("tlist", "range(0,0.001,0.5)");
    model.study("std1").feature("time").set("useinitsol", true);
    model.study("std1").feature("time").set("initstudy", "std1");
    model.study("std1").feature("time").set("solnum", "auto");

    model.sol("sol4").feature("st1").label("Compile Equations: Phase Initialization");
    model.sol("sol4").feature("v1").label("Dependent Variables 1.1");
    model.sol("sol4").feature("v1").set("clist", new String[]{"5.0E-4[s]"});
    model.sol("sol4").feature("s1").label("Stationary Solver 1.1");
    model.sol("sol4").feature("s1").feature("dDef").label("Direct 2");
    model.sol("sol4").feature("s1").feature("aDef").label("Advanced 1");
    model.sol("sol4").feature("s1").feature("aDef").set("cachepattern", true);
    model.sol("sol4").feature("s1").feature("fc1").label("Fully Coupled 1.1");
    model.sol("sol4").feature("s1").feature("fc1").set("linsolver", "i1");
    model.sol("sol4").feature("s1").feature("fc1").set("initstep", 0.01);
    model.sol("sol4").feature("s1").feature("fc1").set("minstep", 1.0E-6);
    model.sol("sol4").feature("s1").feature("fc1").set("maxiter", 50);
    model.sol("sol4").feature("s1").feature("i1").label("AMG, interface distance (pf)");
    model.sol("sol4").feature("s1").feature("i1").set("nlinnormuse", true);
    model.sol("sol4").feature("s1").feature("i1").set("maxlinit", 1000);
    model.sol("sol4").feature("s1").feature("i1").feature("ilDef").label("Incomplete LU 1");
    model.sol("sol4").feature("s1").feature("i1").feature("mg1").label("Multigrid 1.1");
    model.sol("sol4").feature("s1").feature("i1").feature("mg1").set("prefun", "saamg");
    model.sol("sol4").feature("s1").feature("i1").feature("mg1").set("maxcoarsedof", 50000);
    model.sol("sol4").feature("s1").feature("i1").feature("mg1").set("saamgcompwise", true);
    model.sol("sol4").feature("s1").feature("i1").feature("mg1").set("usesmooth", false);
    model.sol("sol4").feature("s1").feature("i1").feature("mg1").feature("pr").label("Presmoother 1");
    model.sol("sol4").feature("s1").feature("i1").feature("mg1").feature("pr").feature("soDef").label("SOR 1");
    model.sol("sol4").feature("s1").feature("i1").feature("mg1").feature("pr").feature("sl1").label("SOR Line 1.1");
    model.sol("sol4").feature("s1").feature("i1").feature("mg1").feature("pr").feature("sl1")
         .set("linesweeptype", "ssor");
    model.sol("sol4").feature("s1").feature("i1").feature("mg1").feature("pr").feature("sl1").set("iter", 1);
    model.sol("sol4").feature("s1").feature("i1").feature("mg1").feature("pr").feature("sl1").set("linerelax", 0.7);
    model.sol("sol4").feature("s1").feature("i1").feature("mg1").feature("pr").feature("sl1")
         .set("linemethod", "uncoupled");
    model.sol("sol4").feature("s1").feature("i1").feature("mg1").feature("pr").feature("sl1").set("relax", 0.5);
    model.sol("sol4").feature("s1").feature("i1").feature("mg1").feature("po").label("Postsmoother 1");
    model.sol("sol4").feature("s1").feature("i1").feature("mg1").feature("po").feature("soDef").label("SOR 1");
    model.sol("sol4").feature("s1").feature("i1").feature("mg1").feature("po").feature("sl1").label("SOR Line 1.1");
    model.sol("sol4").feature("s1").feature("i1").feature("mg1").feature("po").feature("sl1")
         .set("linesweeptype", "ssor");
    model.sol("sol4").feature("s1").feature("i1").feature("mg1").feature("po").feature("sl1").set("iter", 1);
    model.sol("sol4").feature("s1").feature("i1").feature("mg1").feature("po").feature("sl1").set("linerelax", 0.7);
    model.sol("sol4").feature("s1").feature("i1").feature("mg1").feature("po").feature("sl1")
         .set("linemethod", "uncoupled");
    model.sol("sol4").feature("s1").feature("i1").feature("mg1").feature("po").feature("sl1").set("relax", 0.5);
    model.sol("sol4").feature("s1").feature("i1").feature("mg1").feature("cs").label("Coarse Solver 1");
    model.sol("sol4").feature("s1").feature("i1").feature("mg1").feature("cs").feature("dDef").label("Direct 2");
    model.sol("sol4").feature("s1").feature("i1").feature("mg1").feature("cs").feature("d1").label("Direct 1.1");
    model.sol("sol4").feature("s1").feature("i1").feature("mg1").feature("cs").feature("d1")
         .set("linsolver", "pardiso");
    model.sol("sol4").feature("s1").feature("i1").feature("mg1").feature("cs").feature("d1")
         .set("pivotperturb", 1.0E-13);
    model.sol("sol4").feature("s1").feature("d1").label("Direct, interface distance (pf)");
    model.sol("sol4").feature("s1").feature("d1").set("linsolver", "pardiso");
    model.sol("sol4").feature("s1").feature("d1").set("pivotperturb", 1.0E-13);
    model.sol("sol4").feature("su1").label("Solution Store 2");
    model.sol("sol4").feature("st2").label("Compile Equations: Time Dependent");
    model.sol("sol4").feature("st2").set("studystep", "time");
    model.sol("sol4").feature("v2").label("Dependent Variables 2.1");
    model.sol("sol4").feature("v2").set("initsol", "sol4");
    model.sol("sol4").feature("v2").set("solnum", "auto");
    model.sol("sol4").feature("v2").set("resscalemethod", "manual");
    model.sol("sol4").feature("v2").set("notsolmethod", "sol");
    model.sol("sol4").feature("v2").set("notsol", "sol4");
    model.sol("sol4").feature("v2").set("notsolnum", "auto");
    model.sol("sol4").feature("v2").set("clist", new String[]{"range(0,0.001,0.5)", "5.0E-4[s]"});
    model.sol("sol4").feature("v2").feature("comp1_phipf").set("scalemethod", "manual");
    model.sol("sol4").feature("v2").feature("comp1_phipf").set("scaleval", 1);
    model.sol("sol4").feature("t1").label("Time-Dependent Solver 1.1");
    model.sol("sol4").feature("t1").set("control", "time");
    model.sol("sol4").feature("t1").set("tlist", "range(0,0.001,0.5)");
    model.sol("sol4").feature("t1").set("rtol", 0.005);
    model.sol("sol4").feature("t1").set("atolglobalfactor", 0.05);
    model.sol("sol4").feature("t1")
         .set("atolmethod", new String[]{"comp1_Cp2", "global", "comp1_GI", "global", "comp1_mu2", "global", "comp1_p", "scaled", "comp1_phipf", "global", 
         "comp1_psi", "global", "comp1_T", "global", "comp1_u", "global"});
    model.sol("sol4").feature("t1")
         .set("atolfactor", new String[]{"comp1_Cp2", "0.1", "comp1_GI", "0.1", "comp1_mu2", "0.1", "comp1_p", "1", "comp1_phipf", "0.1", 
         "comp1_psi", "0.1", "comp1_T", "0.1", "comp1_u", "0.1"});
    model.sol("sol4").feature("t1").set("maxorder", 2);
    model.sol("sol4").feature("t1").set("stabcntrl", true);
    model.sol("sol4").feature("t1").set("bwinitstepfrac", 0.01);
    model.sol("sol4").feature("t1").set("estrat", "exclude");
    model.sol("sol4").feature("t1").set("rescaleafterinitbw", true);
    model.sol("sol4").feature("t1").feature("dDef").label("Direct 4");
    model.sol("sol4").feature("t1").feature("aDef").label("Advanced 1");
    model.sol("sol4").feature("t1").feature("aDef").set("cachepattern", true);
    model.sol("sol4").feature("t1").feature("se1").label("Segregated 1.1");
    model.sol("sol4").feature("t1").feature("se1").set("ntolfact", 0.5);
    model.sol("sol4").feature("t1").feature("se1").set("segstabacc", "segaacc");
    model.sol("sol4").feature("t1").feature("se1").set("segaaccdim", 5);
    model.sol("sol4").feature("t1").feature("se1").set("segaaccmix", 0.9);
    model.sol("sol4").feature("t1").feature("se1").feature("ss1").label("Specific Heat");
    model.sol("sol4").feature("t1").feature("se1").feature("ss1").set("segvar", new String[]{"comp1_Cp2"});
    model.sol("sol4").feature("t1").feature("se1").feature("ss2").label("Temperature");
    model.sol("sol4").feature("t1").feature("se1").feature("ss2").set("segvar", new String[]{"comp1_T"});
    model.sol("sol4").feature("t1").feature("se1").feature("ss2").set("linsolver", "d1");
    model.sol("sol4").feature("t1").feature("se1").feature("ss2").set("subdamp", "0.8");
    model.sol("sol4").feature("t1").feature("se1").feature("ss2").set("subjtech", "once");
    model.sol("sol4").feature("t1").feature("se1").feature("ss3").label("Velocity u, Pressure p");
    model.sol("sol4").feature("t1").feature("se1").feature("ss3").set("segvar", new String[]{"comp1_u", "comp1_p"});
    model.sol("sol4").feature("t1").feature("se1").feature("ss3").set("linsolver", "i1");
    model.sol("sol4").feature("t1").feature("se1").feature("ss3").set("subdamp", "0.8");
    model.sol("sol4").feature("t1").feature("se1").feature("ss3").set("subjtech", "once");
    model.sol("sol4").feature("t1").feature("se1").feature("ss4").label("Phase field variables");
    model.sol("sol4").feature("t1").feature("se1").feature("ss4").set("linsolver", "i2");
    model.sol("sol4").feature("t1").feature("se1").feature("ss4").set("subdamp", "0.8");
    model.sol("sol4").feature("t1").feature("se1").feature("ss4").set("subjtech", "once");
    model.sol("sol4").feature("t1").feature("se1").feature("ss5").label("Viscosity");
    model.sol("sol4").feature("t1").feature("se1").feature("ss5").set("segvar", new String[]{"comp1_mu2"});
    model.sol("sol4").feature("t1").feature("se1").feature("ll1").label("Lower Limit 1.1");
    model.sol("sol4").feature("t1").feature("se1").feature("ll1").set("lowerlimit", "comp1.T 0 ");
    model.sol("sol4").feature("t1").feature("se1").feature().remove("ht1");
    model.sol("sol4").feature("t1").feature("d1").label("Direct, heat transfer variables (ht)");
    model.sol("sol4").feature("t1").feature("d1").set("linsolver", "pardiso");
    model.sol("sol4").feature("t1").feature("d1").set("pivotperturb", 1.0E-13);
    model.sol("sol4").feature("t1").feature("i1").label("AMG, fluid flow variables (spf)");
    model.sol("sol4").feature("t1").feature("i1").set("maxlinit", 100);
    model.sol("sol4").feature("t1").feature("i1").set("rhob", 20);
    model.sol("sol4").feature("t1").feature("i1").feature("ilDef").label("Incomplete LU 1");
    model.sol("sol4").feature("t1").feature("i1").feature("mg1").label("Multigrid 1.1");
    model.sol("sol4").feature("t1").feature("i1").feature("mg1").set("prefun", "saamg");
    model.sol("sol4").feature("t1").feature("i1").feature("mg1").set("maxcoarsedof", 80000);
    model.sol("sol4").feature("t1").feature("i1").feature("mg1").set("strconn", 0.02);
    model.sol("sol4").feature("t1").feature("i1").feature("mg1").set("saamgcompwise", true);
    model.sol("sol4").feature("t1").feature("i1").feature("mg1").set("usesmooth", false);
    model.sol("sol4").feature("t1").feature("i1").feature("mg1").feature("pr").label("Presmoother 1");
    model.sol("sol4").feature("t1").feature("i1").feature("mg1").feature("pr").feature("soDef").label("SOR 1");
    model.sol("sol4").feature("t1").feature("i1").feature("mg1").feature("pr").feature("sc1").label("SCGS 1.1");
    model.sol("sol4").feature("t1").feature("i1").feature("mg1").feature("pr").feature("sc1")
         .set("linesweeptype", "ssor");
    model.sol("sol4").feature("t1").feature("i1").feature("mg1").feature("pr").feature("sc1").set("iter", 0);
    model.sol("sol4").feature("t1").feature("i1").feature("mg1").feature("pr").feature("sc1").set("approxscgs", true);
    model.sol("sol4").feature("t1").feature("i1").feature("mg1").feature("pr").feature("sc1")
         .set("scgsdirectmaxsize", 1000);
    model.sol("sol4").feature("t1").feature("i1").feature("mg1").feature("po").label("Postsmoother 1");
    model.sol("sol4").feature("t1").feature("i1").feature("mg1").feature("po").feature("soDef").label("SOR 1");
    model.sol("sol4").feature("t1").feature("i1").feature("mg1").feature("po").feature("sc1").label("SCGS 1.1");
    model.sol("sol4").feature("t1").feature("i1").feature("mg1").feature("po").feature("sc1")
         .set("linesweeptype", "ssor");
    model.sol("sol4").feature("t1").feature("i1").feature("mg1").feature("po").feature("sc1").set("iter", 1);
    model.sol("sol4").feature("t1").feature("i1").feature("mg1").feature("po").feature("sc1").set("approxscgs", true);
    model.sol("sol4").feature("t1").feature("i1").feature("mg1").feature("po").feature("sc1")
         .set("scgsdirectmaxsize", 1000);
    model.sol("sol4").feature("t1").feature("i1").feature("mg1").feature("cs").label("Coarse Solver 1");
    model.sol("sol4").feature("t1").feature("i1").feature("mg1").feature("cs").feature("dDef").label("Direct 2");
    model.sol("sol4").feature("t1").feature("i1").feature("mg1").feature("cs").feature("d1").label("Direct 1.1");
    model.sol("sol4").feature("t1").feature("i1").feature("mg1").feature("cs").feature("d1")
         .set("linsolver", "pardiso");
    model.sol("sol4").feature("t1").feature("i1").feature("mg1").feature("cs").feature("d1")
         .set("pivotperturb", 1.0E-13);
    model.sol("sol4").feature("t1").feature("i2").label("AMG, phase field variables (pf)");
    model.sol("sol4").feature("t1").feature("i2").set("maxlinit", 50);
    model.sol("sol4").feature("t1").feature("i2").set("rhob", 20);
    model.sol("sol4").feature("t1").feature("i2").feature("ilDef").label("Incomplete LU 1");
    model.sol("sol4").feature("t1").feature("i2").feature("mg1").label("Multigrid 1.1");
    model.sol("sol4").feature("t1").feature("i2").feature("mg1").set("prefun", "saamg");
    model.sol("sol4").feature("t1").feature("i2").feature("mg1").set("maxcoarsedof", 50000);
    model.sol("sol4").feature("t1").feature("i2").feature("mg1").set("usesmooth", false);
    model.sol("sol4").feature("t1").feature("i2").feature("mg1").feature("pr").label("Presmoother 1");
    model.sol("sol4").feature("t1").feature("i2").feature("mg1").feature("pr").feature("soDef").label("SOR 1");
    model.sol("sol4").feature("t1").feature("i2").feature("mg1").feature("pr").feature("sc1").label("SCGS 1.1");
    model.sol("sol4").feature("t1").feature("i2").feature("mg1").feature("pr").feature("sc1")
         .set("linesweeptype", "ssor");
    model.sol("sol4").feature("t1").feature("i2").feature("mg1").feature("pr").feature("sc1").set("iter", 0);
    model.sol("sol4").feature("t1").feature("i2").feature("mg1").feature("pr").feature("sc1").set("approxscgs", true);
    model.sol("sol4").feature("t1").feature("i2").feature("mg1").feature("pr").feature("sc1")
         .set("scgsdirectmaxsize", 1000);
    model.sol("sol4").feature("t1").feature("i2").feature("mg1").feature("po").label("Postsmoother 1");
    model.sol("sol4").feature("t1").feature("i2").feature("mg1").feature("po").feature("soDef").label("SOR 1");
    model.sol("sol4").feature("t1").feature("i2").feature("mg1").feature("po").feature("sc1").label("SCGS 1.1");
    model.sol("sol4").feature("t1").feature("i2").feature("mg1").feature("po").feature("sc1")
         .set("linesweeptype", "ssor");
    model.sol("sol4").feature("t1").feature("i2").feature("mg1").feature("po").feature("sc1").set("iter", 1);
    model.sol("sol4").feature("t1").feature("i2").feature("mg1").feature("po").feature("sc1").set("approxscgs", true);
    model.sol("sol4").feature("t1").feature("i2").feature("mg1").feature("po").feature("sc1")
         .set("scgsdirectmaxsize", 1000);
    model.sol("sol4").feature("t1").feature("i2").feature("mg1").feature("cs").label("Coarse Solver 1");
    model.sol("sol4").feature("t1").feature("i2").feature("mg1").feature("cs").feature("dDef").label("Direct 2");
    model.sol("sol4").feature("t1").feature("i2").feature("mg1").feature("cs").feature("d1").label("Direct 1.1");
    model.sol("sol4").feature("t1").feature("i2").feature("mg1").feature("cs").feature("d1")
         .set("linsolver", "pardiso");

    return model;
  }

  public static Model run2(Model model) {
    model.sol("sol4").feature("t1").feature("i2").feature("mg1").feature("cs").feature("d1")
         .set("pivotperturb", 1.0E-13);
    model.sol("sol4").feature("t1").feature("i3").label("AMG, heat transfer variables (ht)");
    model.sol("sol4").feature("t1").feature("i3").set("rhob", 20);
    model.sol("sol4").feature("t1").feature("i3").feature("ilDef").label("Incomplete LU 1");
    model.sol("sol4").feature("t1").feature("i3").feature("mg1").label("Multigrid 1.1");
    model.sol("sol4").feature("t1").feature("i3").feature("mg1").set("prefun", "saamg");
    model.sol("sol4").feature("t1").feature("i3").feature("mg1").set("maxcoarsedof", 50000);
    model.sol("sol4").feature("t1").feature("i3").feature("mg1").set("saamgcompwise", true);
    model.sol("sol4").feature("t1").feature("i3").feature("mg1").set("usesmooth", false);
    model.sol("sol4").feature("t1").feature("i3").feature("mg1").feature("pr").label("Presmoother 1");
    model.sol("sol4").feature("t1").feature("i3").feature("mg1").feature("pr").feature("soDef").label("SOR 2");
    model.sol("sol4").feature("t1").feature("i3").feature("mg1").feature("pr").feature("so1").label("SOR 1.1");
    model.sol("sol4").feature("t1").feature("i3").feature("mg1").feature("pr").feature("so1").set("relax", 0.9);
    model.sol("sol4").feature("t1").feature("i3").feature("mg1").feature("po").label("Postsmoother 1");
    model.sol("sol4").feature("t1").feature("i3").feature("mg1").feature("po").feature("soDef").label("SOR 2");
    model.sol("sol4").feature("t1").feature("i3").feature("mg1").feature("po").feature("so1").label("SOR 1.1");
    model.sol("sol4").feature("t1").feature("i3").feature("mg1").feature("po").feature("so1").set("relax", 0.9);
    model.sol("sol4").feature("t1").feature("i3").feature("mg1").feature("cs").label("Coarse Solver 1");
    model.sol("sol4").feature("t1").feature("i3").feature("mg1").feature("cs").feature("dDef").label("Direct 2");
    model.sol("sol4").feature("t1").feature("i3").feature("mg1").feature("cs").feature("d1").label("Direct 1.1");
    model.sol("sol4").feature("t1").feature("i3").feature("mg1").feature("cs").feature("d1")
         .set("linsolver", "pardiso");
    model.sol("sol4").feature("t1").feature("i3").feature("mg1").feature("cs").feature("d1")
         .set("pivotperturb", 1.0E-13);
    model.sol("sol4").feature("t1").feature("d2").label("Direct, fluid flow variables (spf)");
    model.sol("sol4").feature("t1").feature("d2").set("linsolver", "pardiso");
    model.sol("sol4").feature("t1").feature("d2").set("pivotperturb", 1.0E-13);
    model.sol("sol4").feature("t1").feature("d3").label("Direct, phase field variables (pf)");
    model.sol("sol4").feature("t1").feature("d3").set("linsolver", "pardiso");
    model.sol("sol4").feature("t1").feature("d3").set("pivotperturb", 1.0E-13);

    return model;
  }

  public static void main(String[] args) {
    Model model = run();
    run2(model);
  }

}
