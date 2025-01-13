/*
 * ThermalNewtonianIncompressible2D.java
 */

import com.comsol.model.*;
import com.comsol.model.util.*;

/** Model exported on Jan 13 2025, 14:26 by COMSOL 6.1.0.357. */
public class ThermalNewtonianIncompressible2D {

  public static Model run() {
    Model model = ModelUtil.create("Model");

    model.modelPath("/nas/longleaf/home/mshourya/workspace/SLS/Particle");

    model.component()
         .insert("/nas/longleaf/home/mshourya/workspace/SelectiveLaserSintering/Particle/ClusterNewtonianIncompressibleODE.mph", new String[]{"comp1"}, new String[]{});

    model.param().remove(new String[]{"mumatsolid", "Tm", "deltat", "mumatliq"});

    model.label("ThermalNewtonianIncompressibleODE.mph");

    model.param().set("MaterialQuantities", "0");
    model.param().descr("MaterialQuantities", "___________________________________________________");
    model.param().set("rhomat", "1000 [kg/m^3]");
    model.param().descr("rhomat", "Material density");
    model.param().set("rhoair", "1 [kg/m^3]");
    model.param().descr("rhoair", "Air density");
    model.param().set("mumatsolid", "1e6 [Pa*s]");
    model.param().descr("mumatsolid", "Material viscosity in solid phase");
    model.param().set("mumatliq", "1 [Pa*s]");
    model.param().descr("mumatliq", "Material viscosity in liquid phase");
    model.param().set("muair", "1 [Pa*s]");
    model.param().descr("muair", "Air viscosity");
    model.param().set("TemperatureQuantities", "0");
    model.param().descr("TemperatureQuantities", "___________________________________________________");
    model.param().set("deltat", "50 [K]");
    model.param().descr("deltat", "");
    model.param().set("T0", "400[K]");
    model.param().descr("T0", "");
    model.param().set("alpha", "0[K/s]");
    model.param().descr("alpha", "Heating rate");
    model.param().set("Tm", "320 [K]");
    model.param().descr("Tm", "");
    model.param().set("OtherQuantities", "0");
    model.param().descr("OtherQuantities", "___________________________________________________");
    model.param().set("epsilon", "1e-6");
    model.param().descr("epsilon", "interface thickness");
    model.param().set("g", "9.81");
    model.param().descr("g", "Gravitational acceleration");
    model.param().set("sigma", "0.05");
    model.param().descr("sigma", "Surface Tension");

    model.component("comp1").physics().create("ht", "HeatTransferInFluids", "geom2");
    model.component("comp1").physics("ht").prop("InconsistentStabilization").set("HeatIsotropicDiffusion", true);
    model.component("comp1").physics().remove("ge");
    model.component("comp1").physics("ht").field("temperature").field("T");
    model.component("comp1").physics("ht").feature("fluid1").set("u", new String[]{"u", "v", "0"});
    model.component("comp1").physics("ht").feature("fluid1").set("k_mat", "userdef");
    model.component("comp1").physics("ht").feature("fluid1")
         .set("k", new String[]{"kair*(1-phi)/2 + kmat*(1+phi)/2", "0", "0", "0", "kair*(1-phi)/2 + kmat*(1+phi)/2", "0", "0", "0", "kair*(1-phi)/2 + kmat*(1+phi)/2"});
    model.component("comp1").physics("ht").feature("fluid1").set("fluidType", "gasLiquid");
    model.component("comp1").physics("ht").feature("fluid1").set("rho_mat", "userdef");
    model.component("comp1").physics("ht").feature("fluid1").set("rho", "rhoair*(1-phi)/2 + rhomat*(1+phi)/2");
    model.component("comp1").physics("ht").feature("fluid1").set("Cp_mat", "userdef");
    model.component("comp1").physics("ht").feature("fluid1").set("Cp", "cpair*(1-phi)/2 + rhomat*(1+phi)/2");
    model.component("comp1").physics("ht").feature("fluid1").set("gamma_not_IG_mat", "auto");
    model.component("comp1").physics("ht").feature("fluid1")
         .set("k", new String[]{"kair*(1-phipf)/2 + kmat*(1+phipf)/2", "0", "0", "0", "kair*(1-phipf)/2 + kmat*(1+phipf)/2", "0", "0", "0", "kair*(1-phipf)/2 + kmat*(1+phipf)/2"});
    model.component("comp1").physics("ht").feature("fluid1").set("rho", "rhoair*(1-phipf)/2 + rhomat*(1+phipf)/2");
    model.component("comp1").physics("ht").feature("fluid1").set("Cp", "cpair*(1-phipf)/2 + rhomat*(1+phipf)/2");
    model.component("comp1").physics("ht").feature("init1").set("Tinit", "T0");

    model.study().create("std1");
    model.study("std1").create("phasei", "PhaseInitialization");
    model.study("std1").feature("phasei").set("solnum", "auto");
    model.study("std1").feature("phasei").set("notsolnum", "auto");
    model.study("std1").feature("phasei").set("ngenAUX", "1");
    model.study("std1").feature("phasei").set("goalngenAUX", "1");
    model.study("std1").feature("phasei").set("ngenAUX", "1");
    model.study("std1").feature("phasei").set("goalngenAUX", "1");
    model.study("std1").feature("phasei").setSolveFor("/physics/spf", true);
    model.study("std1").feature("phasei").setSolveFor("/physics/pf", true);
    model.study("std1").feature("phasei").setSolveFor("/physics/viscode", true);
    model.study("std1").feature("phasei").setSolveFor("/physics/ht", true);
    model.study("std1").feature("phasei").setSolveFor("/multiphysics/tpf1", true);
    model.study("std1").feature("phasei").setSolveFor("/physics/spf", false);
    model.study("std1").create("time", "Transient");
    model.study("std1").feature("time").set("initstudy", "std1");
    model.study("std1").feature("time").set("notstudy", "std1");
    model.study("std1").feature("time").set("initialtime", "0");
    model.study("std1").feature("time").set("useinitsol", "on");
    model.study("std1").feature("time").set("notsolmethod", "sol");
    model.study("std1").feature("time").setSolveFor("/physics/spf", true);
    model.study("std1").feature("time").setSolveFor("/physics/pf", true);
    model.study("std1").feature("time").setSolveFor("/physics/viscode", true);
    model.study("std1").feature("time").setSolveFor("/physics/ht", true);
    model.study("std1").feature("time").setSolveFor("/multiphysics/tpf1", true);
    model.study("std1").feature("time").set("tlist", "range(0,0.01,1)");

    model.component("comp1").physics("viscode").feature("init1").set("mu2t", "d(mumat(T0), T0)*d(T, t)");

    model.sol().create("sol1");

    model.component("comp1").mesh("mesh2").stat().selection().geom(2);
    model.component("comp1").mesh("mesh2").stat().selection().set(1, 2);
    model.component("comp1").mesh("mesh2").stat().selection().geom(2);
    model.component("comp1").mesh("mesh2").stat().selection().set(1, 2);
    model.component("comp1").mesh("mesh2").stat().selection().geom(2);
    model.component("comp1").mesh("mesh2").stat().selection().set(1, 2);
    model.component("comp1").mesh("mesh2").stat().selection().geom(2);
    model.component("comp1").mesh("mesh2").stat().selection().set(1, 2);

    model.sol("sol1").study("std1");

    model.study("std1").feature("phasei").set("notlistsolnum", 1);
    model.study("std1").feature("phasei").set("notsolnum", "auto");
    model.study("std1").feature("phasei").set("listsolnum", 1);
    model.study("std1").feature("phasei").set("solnum", "auto");
    model.study("std1").feature("time").set("notlistsolnum", 1);
    model.study("std1").feature("time").set("notsolnum", "auto");
    model.study("std1").feature("time").set("listsolnum", 1);
    model.study("std1").feature("time").set("solnum", "auto");

    model.sol("sol1").create("st1", "StudyStep");
    model.sol("sol1").feature("st1").set("study", "std1");
    model.sol("sol1").feature("st1").set("studystep", "phasei");
    model.sol("sol1").create("v1", "Variables");
    model.sol("sol1").feature("v1").set("control", "phasei");
    model.sol("sol1").create("s1", "Stationary");
    model.sol("sol1").feature("s1").feature("aDef").set("cachepattern", true);
    model.sol("sol1").feature("s1").create("fc1", "FullyCoupled");
    model.sol("sol1").feature("s1").feature("fc1").set("dtech", "auto");
    model.sol("sol1").feature("s1").feature("fc1").set("initstep", 0.01);
    model.sol("sol1").feature("s1").feature("fc1").set("minstep", 1.0E-6);
    model.sol("sol1").feature("s1").feature("fc1").set("maxiter", 50);
    model.sol("sol1").feature("s1").create("d1", "Direct");
    model.sol("sol1").feature("s1").feature("d1").set("linsolver", "pardiso");
    model.sol("sol1").feature("s1").feature("d1").set("pivotperturb", 1.0E-13);
    model.sol("sol1").feature("s1").feature("d1").label("Direct, interface distance (pf)");
    model.sol("sol1").feature("s1").create("i1", "Iterative");
    model.sol("sol1").feature("s1").feature("i1").set("linsolver", "gmres");
    model.sol("sol1").feature("s1").feature("i1").set("prefuntype", "left");
    model.sol("sol1").feature("s1").feature("i1").set("itrestart", 50);
    model.sol("sol1").feature("s1").feature("i1").set("rhob", 400);
    model.sol("sol1").feature("s1").feature("i1").set("maxlinit", 1000);
    model.sol("sol1").feature("s1").feature("i1").set("nlinnormuse", "on");
    model.sol("sol1").feature("s1").feature("i1").label("AMG, interface distance (pf)");
    model.sol("sol1").feature("s1").feature("i1").create("mg1", "Multigrid");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").set("prefun", "saamg");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").set("mgcycle", "v");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").set("maxcoarsedof", 50000);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").set("strconn", 0.01);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").set("nullspace", "constant");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").set("usesmooth", false);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").set("saamgcompwise", true);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").set("loweramg", true);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").set("compactaggregation", false);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr").create("sl1", "SORLine");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr").feature("sl1")
         .set("linesweeptype", "ssor");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr").feature("sl1").set("iter", 1);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr").feature("sl1").set("linerelax", 0.7);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr").feature("sl1")
         .set("linealgorithm", "mesh");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr").feature("sl1")
         .set("linemethod", "uncoupled");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr").feature("sl1").set("seconditer", 1);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr").feature("sl1").set("relax", 0.5);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po").create("sl1", "SORLine");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po").feature("sl1")
         .set("linesweeptype", "ssor");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po").feature("sl1").set("iter", 1);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po").feature("sl1").set("linerelax", 0.7);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po").feature("sl1")
         .set("linealgorithm", "mesh");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po").feature("sl1")
         .set("linemethod", "uncoupled");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po").feature("sl1").set("seconditer", 1);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po").feature("sl1").set("relax", 0.5);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("cs").create("d1", "Direct");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("cs").feature("d1")
         .set("linsolver", "pardiso");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("cs").feature("d1")
         .set("pivotperturb", 1.0E-13);
    model.sol("sol1").feature("s1").feature("fc1").set("linsolver", "d1");
    model.sol("sol1").feature("s1").feature("fc1").set("dtech", "auto");
    model.sol("sol1").feature("s1").feature("fc1").set("initstep", 0.01);
    model.sol("sol1").feature("s1").feature("fc1").set("minstep", 1.0E-6);
    model.sol("sol1").feature("s1").feature("fc1").set("maxiter", 50);
    model.sol("sol1").feature("s1").feature().remove("fcDef");
    model.sol("sol1").create("su1", "StoreSolution");
    model.sol("sol1").create("st2", "StudyStep");
    model.sol("sol1").feature("st2").set("study", "std1");
    model.sol("sol1").feature("st2").set("studystep", "time");
    model.sol("sol1").create("v2", "Variables");
    model.sol("sol1").feature("v2").feature("comp1_phipf").set("scalemethod", "manual");
    model.sol("sol1").feature("v2").feature("comp1_phipf").set("scaleval", "1");
    model.sol("sol1").feature("v2").set("initmethod", "sol");
    model.sol("sol1").feature("v2").set("initsol", "sol1");
    model.sol("sol1").feature("v2").set("initsoluse", "su1");
    model.sol("sol1").feature("v2").set("notsolmethod", "sol");
    model.sol("sol1").feature("v2").set("notsol", "sol1");
    model.sol("sol1").feature("v2").set("notsoluse", "su1");
    model.sol("sol1").feature("v2").set("control", "time");
    model.sol("sol1").create("t1", "Time");
    model.sol("sol1").feature("t1").set("tlist", "range(0,0.01,1)");
    model.sol("sol1").feature("t1").set("plot", "off");
    model.sol("sol1").feature("t1").set("plotgroup", "Default");
    model.sol("sol1").feature("t1").set("plotfreq", "tout");
    model.sol("sol1").feature("t1").set("probesel", "all");
    model.sol("sol1").feature("t1").set("probes", new String[]{});
    model.sol("sol1").feature("t1").set("probefreq", "tsteps");
    model.sol("sol1").feature("t1").set("rtol", 0.005);
    model.sol("sol1").feature("t1").set("atolglobalmethod", "scaled");
    model.sol("sol1").feature("t1").set("atolglobalfactor", 0.05);
    model.sol("sol1").feature("t1").set("atolglobalvaluemethod", "factor");
    model.sol("sol1").feature("t1")
         .set("atolmethod", new String[]{"comp1_GI", "global", "comp1_mu2", "global", "comp1_p", "scaled", "comp1_phipf", "global", "comp1_psi", "global", 
         "comp1_T", "global", "comp1_u", "global"});
    model.sol("sol1").feature("t1")
         .set("atol", new String[]{"comp1_GI", "1e-3", "comp1_mu2", "1e-3", "comp1_p", "1e-3", "comp1_phipf", "1e-3", "comp1_psi", "1e-3", 
         "comp1_T", "1e-3", "comp1_u", "1e-3"});
    model.sol("sol1").feature("t1")
         .set("atolvaluemethod", new String[]{"comp1_GI", "factor", "comp1_mu2", "factor", "comp1_p", "factor", "comp1_phipf", "factor", "comp1_psi", "factor", 
         "comp1_T", "factor", "comp1_u", "factor"});
    model.sol("sol1").feature("t1")
         .set("atolfactor", new String[]{"comp1_GI", "0.1", "comp1_mu2", "0.1", "comp1_p", "1", "comp1_phipf", "0.1", "comp1_psi", "0.1", 
         "comp1_T", "0.1", "comp1_u", "0.1"});
    model.sol("sol1").feature("t1").set("reacf", true);
    model.sol("sol1").feature("t1").set("storeudot", true);
    model.sol("sol1").feature("t1").set("endtimeinterpolation", true);
    model.sol("sol1").feature("t1").set("estrat", "exclude");
    model.sol("sol1").feature("t1").set("rhoinf", 0.5);
    model.sol("sol1").feature("t1").set("predictor", "constant");
    model.sol("sol1").feature("t1").set("maxorder", 2);
    model.sol("sol1").feature("t1").set("stabcntrl", true);
    model.sol("sol1").feature("t1").set("rescaleafterinitbw", true);
    model.sol("sol1").feature("t1").set("bwinitstepfrac", "0.01");
    model.sol("sol1").feature("t1").set("control", "time");
    model.sol("sol1").feature("t1").feature("aDef").set("cachepattern", true);
    model.sol("sol1").feature("t1").create("seDef", "Segregated");
    model.sol("sol1").feature("t1").create("se1", "Segregated");
    model.sol("sol1").feature("t1").feature("se1").feature().remove("ssDef");
    model.sol("sol1").feature("t1").feature("se1").create("ss1", "SegregatedStep");
    model.sol("sol1").feature("t1").feature("se1").feature("ss1").set("segvar", new String[]{"comp1_T", "comp1_mu2"});
    model.sol("sol1").feature("t1").feature("se1").feature("ss1").set("subdamp", 0.8);
    model.sol("sol1").feature("t1").feature("se1").feature("ss1").set("subjtech", "once");
    model.sol("sol1").feature("t1").create("d1", "Direct");
    model.sol("sol1").feature("t1").feature("d1").set("linsolver", "mumps");
    model.sol("sol1").feature("t1").feature("d1").set("pivotperturb", 1.0E-13);
    model.sol("sol1").feature("t1").feature("d1").label("Direct, heat transfer variables (ht) (merged)");
    model.sol("sol1").feature("t1").feature("se1").feature("ss1").set("linsolver", "d1");
    model.sol("sol1").feature("t1").feature("se1").feature("ss1").label("Merged variables");
    model.sol("sol1").feature("t1").feature("se1").create("ss2", "SegregatedStep");
    model.sol("sol1").feature("t1").feature("se1").feature("ss2").set("segvar", new String[]{"comp1_u", "comp1_p"});
    model.sol("sol1").feature("t1").feature("se1").feature("ss2").set("subdamp", 0.8);
    model.sol("sol1").feature("t1").feature("se1").feature("ss2").set("subjtech", "once");
    model.sol("sol1").feature("t1").create("d2", "Direct");
    model.sol("sol1").feature("t1").feature("d2").set("linsolver", "pardiso");
    model.sol("sol1").feature("t1").feature("d2").set("pivotperturb", 1.0E-13);
    model.sol("sol1").feature("t1").feature("d2").label("Direct, fluid flow variables (spf)");
    model.sol("sol1").feature("t1").feature("se1").feature("ss2").set("linsolver", "d2");
    model.sol("sol1").feature("t1").feature("se1").feature("ss2").label("Velocity u, Pressure p");
    model.sol("sol1").feature("t1").feature("se1").create("ss3", "SegregatedStep");
    model.sol("sol1").feature("t1").feature("se1").feature("ss3")
         .set("segvar", new String[]{"comp1_phipf", "comp1_psi"});
    model.sol("sol1").feature("t1").feature("se1").feature("ss3").set("subdamp", 0.8);
    model.sol("sol1").feature("t1").feature("se1").feature("ss3").set("subjtech", "once");
    model.sol("sol1").feature("t1").create("d3", "Direct");
    model.sol("sol1").feature("t1").feature("d3").set("linsolver", "pardiso");
    model.sol("sol1").feature("t1").feature("d3").set("pivotperturb", 1.0E-13);
    model.sol("sol1").feature("t1").feature("d3").label("Direct, phase field variables (pf)");
    model.sol("sol1").feature("t1").feature("se1").feature("ss3").set("linsolver", "d3");
    model.sol("sol1").feature("t1").feature("se1").feature("ss3").label("Phase field variables");
    model.sol("sol1").feature("t1").feature("se1").set("segstabacc", "segaacc");
    model.sol("sol1").feature("t1").feature("se1").set("segaaccdim", 5);
    model.sol("sol1").feature("t1").feature("se1").set("segaaccmix", 0.9);
    model.sol("sol1").feature("t1").feature("se1").set("segaaccdelay", 0);
    model.sol("sol1").feature("t1").feature("se1").set("ntolfact", 0.5);
    model.sol("sol1").feature("t1").feature("se1").set("maxsegiter", 10);
    model.sol("sol1").feature("t1").feature("se1").create("ll1", "LowerLimit");
    model.sol("sol1").feature("t1").feature("se1").feature("ll1").set("lowerlimit", "comp1.T 0 ");
    model.sol("sol1").feature("t1").create("i1", "Iterative");
    model.sol("sol1").feature("t1").feature("i1").set("linsolver", "gmres");
    model.sol("sol1").feature("t1").feature("i1").set("prefuntype", "left");
    model.sol("sol1").feature("t1").feature("i1").set("itrestart", 50);
    model.sol("sol1").feature("t1").feature("i1").set("rhob", 20);
    model.sol("sol1").feature("t1").feature("i1").set("maxlinit", 10000);
    model.sol("sol1").feature("t1").feature("i1").set("nlinnormuse", "on");
    model.sol("sol1").feature("t1").feature("i1").label("AMG, heat transfer variables (ht)");
    model.sol("sol1").feature("t1").feature("i1").create("mg1", "Multigrid");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").set("prefun", "saamg");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").set("mgcycle", "v");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").set("maxcoarsedof", 50000);
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").set("strconn", 0.01);
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").set("nullspace", "constant");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").set("usesmooth", false);
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").set("saamgcompwise", true);
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").set("loweramg", true);
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("pr").create("so1", "SOR");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("pr").feature("so1").set("iter", 2);
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("pr").feature("so1").set("relax", 0.9);
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("po").create("so1", "SOR");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("po").feature("so1").set("iter", 2);
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("po").feature("so1").set("relax", 0.9);
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("cs").create("d1", "Direct");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("cs").feature("d1")
         .set("linsolver", "pardiso");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("cs").feature("d1")
         .set("pivotperturb", 1.0E-13);
    model.sol("sol1").feature("t1").create("i2", "Iterative");
    model.sol("sol1").feature("t1").feature("i2").set("linsolver", "gmres");
    model.sol("sol1").feature("t1").feature("i2").set("prefuntype", "left");
    model.sol("sol1").feature("t1").feature("i2").set("itrestart", 50);
    model.sol("sol1").feature("t1").feature("i2").set("rhob", 20);
    model.sol("sol1").feature("t1").feature("i2").set("maxlinit", 100);
    model.sol("sol1").feature("t1").feature("i2").set("nlinnormuse", "on");
    model.sol("sol1").feature("t1").feature("i2").label("AMG, fluid flow variables (spf)");
    model.sol("sol1").feature("t1").feature("i2").create("mg1", "Multigrid");
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").set("prefun", "saamg");
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").set("mgcycle", "v");
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").set("maxcoarsedof", 80000);
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").set("strconn", 0.02);
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").set("nullspace", "constant");
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").set("usesmooth", false);
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").set("saamgcompwise", true);
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").set("loweramg", true);
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").set("compactaggregation", false);
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("pr").create("sc1", "SCGS");
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("pr").feature("sc1")
         .set("linesweeptype", "ssor");
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("pr").feature("sc1").set("iter", 0);
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("pr").feature("sc1").set("scgsrelax", 0.7);
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("pr").feature("sc1")
         .set("scgsmethod", "lines");
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("pr").feature("sc1")
         .set("scgsvertexrelax", 0.7);
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("pr").feature("sc1").set("relax", 0.5);
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("pr").feature("sc1")
         .set("scgssolv", "stored");
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("pr").feature("sc1").set("approxscgs", true);
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("pr").feature("sc1")
         .set("scgsdirectmaxsize", 1000);
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("po").create("sc1", "SCGS");
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("po").feature("sc1")
         .set("linesweeptype", "ssor");
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("po").feature("sc1").set("iter", 1);
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("po").feature("sc1").set("scgsrelax", 0.7);
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("po").feature("sc1")
         .set("scgsmethod", "lines");
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("po").feature("sc1")
         .set("scgsvertexrelax", 0.7);
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("po").feature("sc1").set("relax", 0.5);
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("po").feature("sc1")
         .set("scgssolv", "stored");
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("po").feature("sc1").set("approxscgs", true);
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("po").feature("sc1")
         .set("scgsdirectmaxsize", 1000);
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("cs").create("d1", "Direct");
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("cs").feature("d1")
         .set("linsolver", "pardiso");
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("cs").feature("d1")
         .set("pivotperturb", 1.0E-13);
    model.sol("sol1").feature("t1").create("i3", "Iterative");
    model.sol("sol1").feature("t1").feature("i3").set("linsolver", "gmres");
    model.sol("sol1").feature("t1").feature("i3").set("prefuntype", "left");
    model.sol("sol1").feature("t1").feature("i3").set("itrestart", 50);
    model.sol("sol1").feature("t1").feature("i3").set("rhob", 20);
    model.sol("sol1").feature("t1").feature("i3").set("maxlinit", 50);
    model.sol("sol1").feature("t1").feature("i3").set("nlinnormuse", "on");
    model.sol("sol1").feature("t1").feature("i3").label("AMG, phase field variables (pf)");
    model.sol("sol1").feature("t1").feature("i3").create("mg1", "Multigrid");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").set("prefun", "saamg");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").set("mgcycle", "v");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").set("maxcoarsedof", 50000);
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").set("strconn", 0.01);
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").set("nullspace", "constant");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").set("usesmooth", false);
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").set("saamgcompwise", false);
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").set("loweramg", true);
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").set("compactaggregation", false);
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("pr").create("sc1", "SCGS");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("pr").feature("sc1")
         .set("linesweeptype", "ssor");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("pr").feature("sc1").set("iter", 0);
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("pr").feature("sc1").set("scgsrelax", 0.7);
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("pr").feature("sc1")
         .set("scgsmethod", "lines");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("pr").feature("sc1")
         .set("scgsvertexrelax", 0.7);
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("pr").feature("sc1").set("relax", 0.5);
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("pr").feature("sc1")
         .set("scgssolv", "stored");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("pr").feature("sc1").set("approxscgs", true);
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("pr").feature("sc1")
         .set("scgsdirectmaxsize", 1000);
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("po").create("sc1", "SCGS");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("po").feature("sc1")
         .set("linesweeptype", "ssor");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("po").feature("sc1").set("iter", 1);
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("po").feature("sc1").set("scgsrelax", 0.7);
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("po").feature("sc1")
         .set("scgsmethod", "lines");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("po").feature("sc1")
         .set("scgsvertexrelax", 0.7);
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("po").feature("sc1").set("relax", 0.5);
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("po").feature("sc1")
         .set("scgssolv", "stored");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("po").feature("sc1").set("approxscgs", true);
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("po").feature("sc1")
         .set("scgsdirectmaxsize", 1000);
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("cs").create("d1", "Direct");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("cs").feature("d1")
         .set("linsolver", "pardiso");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("cs").feature("d1")
         .set("pivotperturb", 1.0E-13);
    model.sol("sol1").feature("t1").feature().remove("fcDef");
    model.sol("sol1").feature("t1").feature().remove("seDef");
    model.sol("sol1").feature("v2").set("solnum", "auto");
    model.sol("sol1").feature("v2").set("solvertype", "solnum");
    model.sol("sol1").feature("v2").set("listsolnum", new String[]{"1"});
    model.sol("sol1").feature("v2").set("solnum", "auto");
    model.sol("sol1").feature("v2").set("listsolnum", new String[]{"1"});
    model.sol("sol1").feature("v2").set("solnum", "auto");
    model.sol("sol1").feature("v2").set("control", "time");
    model.sol("sol1").attach("std1");

    model.result().dataset("dset1").set("geom", "geom2");
    model.result().create("pg1", "PlotGroup2D");
    model.result("pg1").label("Velocity (spf)");
    model.result("pg1").set("frametype", "spatial");
    model.result("pg1").set("defaultPlotID", "ResultDefaults_SinglePhaseFlow/icom1/pdef1/pcond1/pg1");
    model.result("pg1").feature().create("surf1", "Surface");
    model.result("pg1").feature("surf1").label("Surface");
    model.result("pg1").feature("surf1").set("showsolutionparams", "on");
    model.result("pg1").feature("surf1").set("smooth", "internal");
    model.result("pg1").feature("surf1").set("showsolutionparams", "on");
    model.result("pg1").feature("surf1").set("data", "parent");
    model.result().create("pg2", "PlotGroup2D");
    model.result("pg2").label("Pressure (spf)");
    model.result("pg2").set("frametype", "spatial");
    model.result("pg2").set("defaultPlotID", "ResultDefaults_SinglePhaseFlow/icom1/pdef1/pcond1/pg2");
    model.result("pg2").feature().create("con1", "Contour");
    model.result("pg2").feature("con1").label("Contour");
    model.result("pg2").feature("con1").set("showsolutionparams", "on");
    model.result("pg2").feature("con1").set("expr", "p");
    model.result("pg2").feature("con1").set("number", 40);
    model.result("pg2").feature("con1").set("levelrounding", false);
    model.result("pg2").feature("con1").set("smooth", "internal");
    model.result("pg2").feature("con1").set("showsolutionparams", "on");
    model.result("pg2").feature("con1").set("data", "parent");
    model.result().create("pg3", "PlotGroup2D");
    model.result("pg3").label("Volume Fraction of Fluid 1 (pf)");

    return model;
  }

  public static Model run2(Model model) {
    model.result("pg3").set("frametype", "spatial");
    model.result("pg3").set("defaultPlotID", "StandaloneTwoPhaseFlowPhysicsInterfaces/icom1/pdef1/pcond1/pg1");
    model.result("pg3").feature().create("surf1", "Surface");
    model.result("pg3").feature("surf1").set("expr", "pf.Vf1");
    model.result("pg3").feature("surf1").set("smooth", "internal");
    model.result("pg3").feature("surf1").set("data", "parent");
    model.result("pg3").feature().create("con1", "Contour");
    model.result("pg3").feature("con1").set("expr", "pf.Vf1");
    model.result("pg3").feature("con1").set("levelmethod", "levels");
    model.result("pg3").feature("con1").set("levels", "0.5");
    model.result("pg3").feature("con1").set("coloring", "uniform");
    model.result("pg3").feature("con1").set("colorlegend", false);
    model.result("pg3").feature("con1").set("color", "gray");
    model.result("pg3").feature("con1").set("smooth", "none");
    model.result("pg3").feature("con1").set("data", "parent");
    model.result().create("pg4", "PlotGroup2D");
    model.result("pg4").set("data", "dset1");
    model.result("pg4").create("surf1", "Surface");
    model.result("pg4").label("Viscosity");
    model.result("pg4").feature("surf1").set("expr", "mu2");
    model.result().create("pg5", "PlotGroup2D");
    model.result("pg5").label("Temperature (ht)");
    model.result("pg5").set("defaultPlotID", "ht/HT_PhysicsInterfaces/icom8/pdef1/pcond2/pcond1/pcond2/pg1");
    model.result("pg5").feature().create("surf1", "Surface");
    model.result("pg5").feature("surf1").label("Surface");
    model.result("pg5").feature("surf1").set("expr", "T");
    model.result("pg5").feature("surf1").set("colortable", "HeatCameraLight");
    model.result("pg5").feature("surf1").set("data", "parent");
    model.result().create("pg6", "PlotGroup2D");
    model.result("pg6").label("Isothermal Contours (ht)");
    model.result("pg6").set("defaultPlotID", "ht/HT_PhysicsInterfaces/icom8/pdef1/pcond2/pg1");
    model.result("pg6").feature().create("con1", "Contour");
    model.result("pg6").feature("con1").label("Contour");
    model.result("pg6").feature("con1").set("expr", "T");
    model.result("pg6").feature("con1").set("levelrounding", false);
    model.result("pg6").feature("con1").set("colortable", "HeatCameraLight");
    model.result("pg6").feature("con1").set("smooth", "internal");
    model.result("pg6").feature("con1").set("data", "parent");
    model.result().remove("pg2");
    model.result().remove("pg1");
    model.result().remove("pg4");
    model.result().remove("pg3");
    model.result().remove("pg6");
    model.result().remove("pg5");

    model.sol("sol1").feature("t1").set("initialstepbdfactive", true);
    model.sol("sol1").feature("t1").set("initialstepbdf", "1e-10");

    model.result().dataset("dset1").set("geom", "geom2");
    model.result().create("pg1", "PlotGroup2D");
    model.result("pg1").label("Velocity (spf)");
    model.result("pg1").set("frametype", "spatial");
    model.result("pg1").set("defaultPlotID", "ResultDefaults_SinglePhaseFlow/icom1/pdef1/pcond1/pg1");
    model.result("pg1").feature().create("surf1", "Surface");
    model.result("pg1").feature("surf1").label("Surface");
    model.result("pg1").feature("surf1").set("showsolutionparams", "on");
    model.result("pg1").feature("surf1").set("smooth", "internal");
    model.result("pg1").feature("surf1").set("showsolutionparams", "on");
    model.result("pg1").feature("surf1").set("data", "parent");
    model.result().create("pg2", "PlotGroup2D");
    model.result("pg2").label("Pressure (spf)");
    model.result("pg2").set("frametype", "spatial");
    model.result("pg2").set("defaultPlotID", "ResultDefaults_SinglePhaseFlow/icom1/pdef1/pcond1/pg2");
    model.result("pg2").feature().create("con1", "Contour");
    model.result("pg2").feature("con1").label("Contour");
    model.result("pg2").feature("con1").set("showsolutionparams", "on");
    model.result("pg2").feature("con1").set("expr", "p");
    model.result("pg2").feature("con1").set("number", 40);
    model.result("pg2").feature("con1").set("levelrounding", false);
    model.result("pg2").feature("con1").set("smooth", "internal");
    model.result("pg2").feature("con1").set("showsolutionparams", "on");
    model.result("pg2").feature("con1").set("data", "parent");
    model.result().create("pg3", "PlotGroup2D");
    model.result("pg3").label("Volume Fraction of Fluid 1 (pf)");
    model.result("pg3").set("frametype", "spatial");
    model.result("pg3").set("defaultPlotID", "StandaloneTwoPhaseFlowPhysicsInterfaces/icom1/pdef1/pcond1/pg1");
    model.result("pg3").feature().create("surf1", "Surface");
    model.result("pg3").feature("surf1").set("expr", "pf.Vf1");
    model.result("pg3").feature("surf1").set("smooth", "internal");
    model.result("pg3").feature("surf1").set("data", "parent");
    model.result("pg3").feature().create("con1", "Contour");
    model.result("pg3").feature("con1").set("expr", "pf.Vf1");
    model.result("pg3").feature("con1").set("levelmethod", "levels");
    model.result("pg3").feature("con1").set("levels", "0.5");
    model.result("pg3").feature("con1").set("coloring", "uniform");
    model.result("pg3").feature("con1").set("colorlegend", false);
    model.result("pg3").feature("con1").set("color", "gray");
    model.result("pg3").feature("con1").set("smooth", "none");
    model.result("pg3").feature("con1").set("data", "parent");
    model.result().create("pg4", "PlotGroup2D");
    model.result("pg4").set("data", "dset1");
    model.result("pg4").create("surf1", "Surface");
    model.result("pg4").label("Viscosity");
    model.result("pg4").feature("surf1").set("expr", "mu2");
    model.result().create("pg5", "PlotGroup2D");
    model.result("pg5").label("Temperature (ht)");
    model.result("pg5").set("defaultPlotID", "ht/HT_PhysicsInterfaces/icom8/pdef1/pcond2/pcond1/pcond2/pg1");
    model.result("pg5").feature().create("surf1", "Surface");
    model.result("pg5").feature("surf1").label("Surface");
    model.result("pg5").feature("surf1").set("expr", "T");
    model.result("pg5").feature("surf1").set("colortable", "HeatCameraLight");
    model.result("pg5").feature("surf1").set("data", "parent");
    model.result().create("pg6", "PlotGroup2D");
    model.result("pg6").label("Isothermal Contours (ht)");
    model.result("pg6").set("defaultPlotID", "ht/HT_PhysicsInterfaces/icom8/pdef1/pcond2/pg1");
    model.result("pg6").feature().create("con1", "Contour");
    model.result("pg6").feature("con1").label("Contour");
    model.result("pg6").feature("con1").set("expr", "T");
    model.result("pg6").feature("con1").set("levelrounding", false);
    model.result("pg6").feature("con1").set("colortable", "HeatCameraLight");
    model.result("pg6").feature("con1").set("smooth", "internal");
    model.result("pg6").feature("con1").set("data", "parent");
    model.result().remove("pg2");
    model.result().remove("pg1");
    model.result().remove("pg4");
    model.result().remove("pg3");
    model.result().remove("pg6");
    model.result().remove("pg5");

    model.param().set("kair", "1");
    model.param().set("kmat", "1");
    model.param().set("cpair", "1");
    model.param().set("cpmat", "1");

    model.result().dataset("dset1").set("geom", "geom2");
    model.result().create("pg1", "PlotGroup2D");
    model.result("pg1").label("Velocity (spf)");
    model.result("pg1").set("frametype", "spatial");
    model.result("pg1").set("defaultPlotID", "ResultDefaults_SinglePhaseFlow/icom1/pdef1/pcond1/pg1");
    model.result("pg1").feature().create("surf1", "Surface");
    model.result("pg1").feature("surf1").label("Surface");
    model.result("pg1").feature("surf1").set("showsolutionparams", "on");
    model.result("pg1").feature("surf1").set("smooth", "internal");
    model.result("pg1").feature("surf1").set("showsolutionparams", "on");
    model.result("pg1").feature("surf1").set("data", "parent");
    model.result().create("pg2", "PlotGroup2D");
    model.result("pg2").label("Pressure (spf)");
    model.result("pg2").set("frametype", "spatial");
    model.result("pg2").set("defaultPlotID", "ResultDefaults_SinglePhaseFlow/icom1/pdef1/pcond1/pg2");
    model.result("pg2").feature().create("con1", "Contour");
    model.result("pg2").feature("con1").label("Contour");
    model.result("pg2").feature("con1").set("showsolutionparams", "on");
    model.result("pg2").feature("con1").set("expr", "p");
    model.result("pg2").feature("con1").set("number", 40);
    model.result("pg2").feature("con1").set("levelrounding", false);
    model.result("pg2").feature("con1").set("smooth", "internal");
    model.result("pg2").feature("con1").set("showsolutionparams", "on");
    model.result("pg2").feature("con1").set("data", "parent");
    model.result().create("pg3", "PlotGroup2D");
    model.result("pg3").label("Volume Fraction of Fluid 1 (pf)");
    model.result("pg3").set("frametype", "spatial");
    model.result("pg3").set("defaultPlotID", "StandaloneTwoPhaseFlowPhysicsInterfaces/icom1/pdef1/pcond1/pg1");
    model.result("pg3").feature().create("surf1", "Surface");
    model.result("pg3").feature("surf1").set("expr", "pf.Vf1");
    model.result("pg3").feature("surf1").set("smooth", "internal");
    model.result("pg3").feature("surf1").set("data", "parent");
    model.result("pg3").feature().create("con1", "Contour");
    model.result("pg3").feature("con1").set("expr", "pf.Vf1");
    model.result("pg3").feature("con1").set("levelmethod", "levels");
    model.result("pg3").feature("con1").set("levels", "0.5");
    model.result("pg3").feature("con1").set("coloring", "uniform");
    model.result("pg3").feature("con1").set("colorlegend", false);
    model.result("pg3").feature("con1").set("color", "gray");
    model.result("pg3").feature("con1").set("smooth", "none");
    model.result("pg3").feature("con1").set("data", "parent");
    model.result().create("pg4", "PlotGroup2D");
    model.result("pg4").set("data", "dset1");
    model.result("pg4").create("surf1", "Surface");
    model.result("pg4").label("Viscosity");
    model.result("pg4").feature("surf1").set("expr", "mu2");
    model.result().create("pg5", "PlotGroup2D");
    model.result("pg5").label("Temperature (ht)");
    model.result("pg5").set("defaultPlotID", "ht/HT_PhysicsInterfaces/icom8/pdef1/pcond2/pcond1/pcond2/pg1");
    model.result("pg5").feature().create("surf1", "Surface");
    model.result("pg5").feature("surf1").label("Surface");
    model.result("pg5").feature("surf1").set("expr", "T");
    model.result("pg5").feature("surf1").set("colortable", "HeatCameraLight");
    model.result("pg5").feature("surf1").set("data", "parent");
    model.result().create("pg6", "PlotGroup2D");
    model.result("pg6").label("Isothermal Contours (ht)");
    model.result("pg6").set("defaultPlotID", "ht/HT_PhysicsInterfaces/icom8/pdef1/pcond2/pg1");
    model.result("pg6").feature().create("con1", "Contour");
    model.result("pg6").feature("con1").label("Contour");
    model.result("pg6").feature("con1").set("expr", "T");
    model.result("pg6").feature("con1").set("levelrounding", false);
    model.result("pg6").feature("con1").set("colortable", "HeatCameraLight");
    model.result("pg6").feature("con1").set("smooth", "internal");
    model.result("pg6").feature("con1").set("data", "parent");
    model.result().remove("pg2");
    model.result().remove("pg1");
    model.result().remove("pg4");
    model.result().remove("pg3");
    model.result().remove("pg6");
    model.result().remove("pg5");

    model.sol("sol1").feature("t1").feature("d1").set("linsolver", "pardiso");

    model.result().dataset("dset1").set("geom", "geom2");
    model.result().create("pg1", "PlotGroup2D");
    model.result("pg1").label("Velocity (spf)");
    model.result("pg1").set("frametype", "spatial");
    model.result("pg1").set("defaultPlotID", "ResultDefaults_SinglePhaseFlow/icom1/pdef1/pcond1/pg1");
    model.result("pg1").feature().create("surf1", "Surface");
    model.result("pg1").feature("surf1").label("Surface");
    model.result("pg1").feature("surf1").set("showsolutionparams", "on");
    model.result("pg1").feature("surf1").set("smooth", "internal");
    model.result("pg1").feature("surf1").set("showsolutionparams", "on");
    model.result("pg1").feature("surf1").set("data", "parent");
    model.result().create("pg2", "PlotGroup2D");
    model.result("pg2").label("Pressure (spf)");
    model.result("pg2").set("frametype", "spatial");
    model.result("pg2").set("defaultPlotID", "ResultDefaults_SinglePhaseFlow/icom1/pdef1/pcond1/pg2");
    model.result("pg2").feature().create("con1", "Contour");
    model.result("pg2").feature("con1").label("Contour");
    model.result("pg2").feature("con1").set("showsolutionparams", "on");
    model.result("pg2").feature("con1").set("expr", "p");
    model.result("pg2").feature("con1").set("number", 40);
    model.result("pg2").feature("con1").set("levelrounding", false);
    model.result("pg2").feature("con1").set("smooth", "internal");
    model.result("pg2").feature("con1").set("showsolutionparams", "on");
    model.result("pg2").feature("con1").set("data", "parent");
    model.result().create("pg3", "PlotGroup2D");
    model.result("pg3").label("Volume Fraction of Fluid 1 (pf)");
    model.result("pg3").set("frametype", "spatial");
    model.result("pg3").set("defaultPlotID", "StandaloneTwoPhaseFlowPhysicsInterfaces/icom1/pdef1/pcond1/pg1");
    model.result("pg3").feature().create("surf1", "Surface");
    model.result("pg3").feature("surf1").set("expr", "pf.Vf1");
    model.result("pg3").feature("surf1").set("smooth", "internal");
    model.result("pg3").feature("surf1").set("data", "parent");
    model.result("pg3").feature().create("con1", "Contour");
    model.result("pg3").feature("con1").set("expr", "pf.Vf1");
    model.result("pg3").feature("con1").set("levelmethod", "levels");
    model.result("pg3").feature("con1").set("levels", "0.5");
    model.result("pg3").feature("con1").set("coloring", "uniform");
    model.result("pg3").feature("con1").set("colorlegend", false);
    model.result("pg3").feature("con1").set("color", "gray");
    model.result("pg3").feature("con1").set("smooth", "none");
    model.result("pg3").feature("con1").set("data", "parent");
    model.result().create("pg4", "PlotGroup2D");
    model.result("pg4").set("data", "dset1");
    model.result("pg4").create("surf1", "Surface");
    model.result("pg4").label("Viscosity");
    model.result("pg4").feature("surf1").set("expr", "mu2");
    model.result().create("pg5", "PlotGroup2D");
    model.result("pg5").label("Temperature (ht)");
    model.result("pg5").set("defaultPlotID", "ht/HT_PhysicsInterfaces/icom8/pdef1/pcond2/pcond1/pcond2/pg1");
    model.result("pg5").feature().create("surf1", "Surface");
    model.result("pg5").feature("surf1").label("Surface");
    model.result("pg5").feature("surf1").set("expr", "T");
    model.result("pg5").feature("surf1").set("colortable", "HeatCameraLight");
    model.result("pg5").feature("surf1").set("data", "parent");
    model.result().create("pg6", "PlotGroup2D");
    model.result("pg6").label("Isothermal Contours (ht)");
    model.result("pg6").set("defaultPlotID", "ht/HT_PhysicsInterfaces/icom8/pdef1/pcond2/pg1");
    model.result("pg6").feature().create("con1", "Contour");
    model.result("pg6").feature("con1").label("Contour");
    model.result("pg6").feature("con1").set("expr", "T");
    model.result("pg6").feature("con1").set("levelrounding", false);
    model.result("pg6").feature("con1").set("colortable", "HeatCameraLight");
    model.result("pg6").feature("con1").set("smooth", "internal");
    model.result("pg6").feature("con1").set("data", "parent");
    model.result().remove("pg2");
    model.result().remove("pg1");
    model.result().remove("pg4");
    model.result().remove("pg3");
    model.result().remove("pg6");
    model.result().remove("pg5");
    model.result().create("pg1", "PlotGroup2D");
    model.result("pg1").run();
    model.result("pg1").label("Velocity and Volume Fraction");
    model.result("pg1").create("surf1", "Surface");
    model.result("pg1").run();
    model.result("pg1").label("Temperature and Volume Fraction");
    model.result("pg1").run();
    model.result("pg1").feature("surf1").set("expr", "ht.T");
    model.result("pg1").run();
    model.result("pg1").feature("surf1").set("colortable", "Magma");
    model.result("pg1").run();
    model.result("pg1").feature("surf1").set("rangecoloractive", true);
    model.result("pg1").feature("surf1").set("rangecolormin", 300);
    model.result("pg1").run();
    model.result("pg1").run();
    model.result("pg1").create("con1", "Contour");
    model.result("pg1").feature("con1").set("expr", "pf.Vf2");
    model.result("pg1").run();
    model.result("pg1").feature("con1").set("coloring", "uniform");
    model.result("pg1").feature("con1").set("color", "black");
    model.result("pg1").feature("con1").set("number", 10);
    model.result("pg1").run();
    model.result("pg1").run();
    model.result("pg1").stepFirst(0);
    model.result("pg1").run();
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result().create("pg2", "PlotGroup2D");
    model.result("pg2").run();
    model.result("pg2").label("Pressure and Volume Fraction");
    model.result("pg2").create("surf1", "Surface");
    model.result("pg2").feature("surf1").set("expr", "pf.Vf2");
    model.result("pg2").run();
    model.result("pg2").feature("surf1").set("rangecoloractive", true);
    model.result("pg2").feature("surf1").set("coloring", "gradient");
    model.result("pg2").feature("surf1").set("topcolor", "blue");
    model.result("pg2").feature("surf1").set("bottomcolor", "white");
    model.result("pg2").run();
    model.result("pg2").run();
    model.result("pg2").create("con1", "Contour");
    model.result("pg2").feature("con1").set("expr", "spf.p");
    model.result("pg2").feature("con1").set("unit", "atm");
    model.result("pg2").feature("con1").set("expr", "p");
    model.result("pg2").run();
    model.result("pg2").feature("con1").set("coloring", "gradient");
    model.result("pg2").feature("con1").set("topcolor", "green");
    model.result("pg2").run();
    model.result("pg2").feature("con1").stepFirst(0);
    model.result("pg2").run();
    model.result("pg2").feature("con1").stepNext(0);
    model.result("pg2").run();
    model.result("pg2").feature("con1").stepNext(0);
    model.result("pg2").run();

    model.component("comp1").physics().create("ewbe", "ElectromagneticWavesBeamEnvelopes", "geom2");

    model.study("std1").feature("phasei").setSolveFor("/physics/ewbe", true);
    model.study("std1").feature("time").setSolveFor("/physics/ewbe", true);

    model.component("comp1").physics().create("ht2", "HeatTransfer", "geom2");

    model.study("std1").feature("phasei").setSolveFor("/physics/ht2", true);
    model.study("std1").feature("time").setSolveFor("/physics/ht2", true);

    model.component("comp1").physics("ht2").prop("PhysicalModelProperty").set("dz", "1[m]");

    model.component("comp1").multiphysics().create("emh1", "ElectromagneticHeating", 2);

    model.study("std1").feature("phasei").setSolveFor("/multiphysics/emh1", true);
    model.study("std1").feature("time").setSolveFor("/multiphysics/emh1", true);

    model.component("comp1").multiphysics("emh1").set("EMHeat_physics", "ewbe");
    model.component("comp1").multiphysics("emh1").set("Heat_physics", "ht2");
    model.component("comp1").multiphysics("emh1").selection().all();

    model.component("comp1").physics().remove("ewbe");
    model.component("comp1").physics().remove("ht2");

    model.label("ThermalNewtonianIncompressibleODE.mph");

    model.component("comp1").physics("ht").create("temp1", "TemperatureBoundary", 1);
    model.component("comp1").physics("ht").feature("temp1").selection().set(3);
    model.component("comp1").physics("ht").feature("temp1").set("T0", "T0+alpha*t");
    model.component("comp1").physics("ht").create("ofl1", "ConvectiveOutflow", 1);
    model.component("comp1").physics("ht").feature("ofl1").selection().set(3);

    model.result("pg2").run();
    model.result("pg2").stepNext(0);
    model.result("pg2").run();

    model.param().set("T0", "300[K]");
    model.param().set("alpha", "1[K/s]");

    model.sol("sol1").runAll();

    model.result("pg1").run();
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result("pg1").run();
    model.result("pg1").feature("con1").set("color", "green");
    model.result("pg1").feature("con1").stepNext(0);
    model.result("pg1").run();
    model.result("pg1").feature("con1").stepNext(0);
    model.result("pg1").run();
    model.result("pg1").feature("con1").stepNext(0);
    model.result("pg1").run();
    model.result("pg1").feature("con1").stepLast(0);
    model.result("pg1").run();
    model.result().table().create("evl2", "Table");
    model.result().table("evl2").comments("Interactive 2D values");
    model.result().table("evl2").label("Evaluation 2D");
    model.result().table("evl2")
         .addRow(new double[]{-0.5731449127197266, 23.2696475982666, 300}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-2.5218353271484375, 23.6135311126709, 300}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-3.7827510833740234, 24.301305770874023, 300}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-4.470521926879883, 24.75982093811035, 300}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-4.814409255981445, 24.75982093811035, 300}, new double[]{0, 0, 0});

    model.component("comp1").physics("ht").feature("temp1").selection().set(3);
    model.component("comp1").physics("ht").feature().remove("temp1");
    model.component("comp1").physics("ht").create("hs1", "HeatSource", 2);
    model.component("comp1").physics("ht").feature("hs1").selection().set(2);
    model.component("comp1").physics("ht").feature("hs1").set("Q0", 1);

    model.sol("sol1").runAll();

    model.result("pg1").run();
    model.result().table("evl2")
         .addRow(new double[]{-2.407205581665039, -21.550214767456055, 300.000001317784}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-4.929037094116211, -20.17466926574707, 300.000001317784}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{0.5731430053710938, -15.474888801574707, 300.000001317784}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{0.22925758361816406, -19.601524353027344, 300.000001317784}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-6.419212341308594, -17.652835845947266, 300.000001317784}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-9.284933090209961, -13.984713554382324, 300.000001317784}, new double[]{0, 0, 0});
    model.result("pg1").stepFirst(0);
    model.result("pg1").run();
    model.result().table("evl2")
         .addRow(new double[]{1.3755455017089844, -23.957420349121094, 300}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{0.8024005889892578, -21.32095718383789, 300}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-0.9170303344726562, -20.6331844329834, 300}, new double[]{0, 0, 0});
    model.result("pg1").stepNext(0);
    model.result("pg1").run();

    model.component("comp1").physics("ht").feature("hs1").setIndex("materialType", "nonSolid", 0);
    model.component("comp1").physics("ht").feature("hs1").set("Q0", 1000);

    model.study("std1").feature("time").set("tlist", "range(0,0.01,5)");

    model.result("pg1").run();
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result("pg1").stepLast(0);
    model.result("pg1").run();
    model.result().table("evl2")
         .addRow(new double[]{1.1462879180908203, -17.4235782623291, 300.0047849753266}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-1.7194328308105469, -16.391918182373047, 300.0047849757099}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-4.126638412475586, -17.996723175048828, 300.0047849758177}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-4.585151672363281, -19.143009185791016, 300.0047849758211}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-2.407205581665039, -23.498905181884766, 300.00478497591877}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{4.126636505126953, -21.77947235107422, 300.0047849758566}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{4.585151672363281, -20.862442016601562, 300.00478497584106}, new double[]{0, 0, 0});

    model.component("comp1").physics("ht").feature("fluid1").set("Cp", "cpair*(1-phipf)/2 + cpmat*(1+phipf)/2");
    model.component("comp1").physics("ht").feature("hs1").set("Q0", 1);

    model.result("pg1").run();
    model.result("pg1").stepFirst(0);
    model.result("pg1").run();
    model.result().table("evl2")
         .addRow(new double[]{-11.118993759155273, -17.76746368408203, 300}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{17.423580169677734, -14.787116050720215, 300}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{2.177947998046875, -20.28929901123047, 300}, new double[]{0, 0, 0});
    model.result("pg1").stepLast(0);
    model.result("pg1").run();
    model.result().table("evl2")
         .addRow(new double[]{-2.9803504943847656, -20.17466926574707, 300.0020347213692}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-2.9803504943847656, -20.060039520263672, 300.0020347213692}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-5.158294677734375, -19.83078384399414, 300.0020347213692}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-5.158294677734375, -19.83078384399414, 300.0020347213692}, new double[]{0, 0, 0});

    model.component("comp1").physics("ht").feature("hs1").set("Q0", 100);

    model.result("pg1").run();
    model.result().table("evl2")
         .addRow(new double[]{7.221614837646484, -21.206329345703125, 300.00477634114634}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{1.948690414428711, -21.43558692932129, 300.0047763412899}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-1.1462879180908203, -17.88209342956543, 300.00477634124735}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-3.3242359161376953, -17.76746368408203, 300.0047763412189}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-6.648469924926758, -15.360260009765625, 300.0047763411153}, new double[]{0, 0, 0});
    model.result("pg2").run();
    model.result("pg2").stepPrevious(0);
    model.result("pg2").run();
    model.result("pg2").stepPrevious(0);
    model.result("pg2").run();
    model.result("pg2").stepPrevious(0);
    model.result("pg2").run();
    model.result("pg2").stepPrevious(0);
    model.result("pg2").run();
    model.result("pg2").stepPrevious(0);
    model.result("pg2").run();
    model.result("pg2").stepNext(0);
    model.result("pg2").run();
    model.result("pg2").stepNext(0);
    model.result("pg2").run();
    model.result("pg2").stepNext(0);
    model.result("pg2").run();
    model.result("pg2").stepNext(0);
    model.result("pg2").run();
    model.result("pg2").stepNext(0);
    model.result("pg2").run();
    model.result("pg1").run();

    model.component("comp1").physics("ht").feature("hs1").selection().set(1, 2);
    model.component("comp1").physics("ht").feature("hs1").set("heatSourceType", "LinearSource");
    model.component("comp1").physics("ht").feature().remove("hs1");
    model.component("comp1").physics("ht").create("temp1", "TemperatureBoundary", 1);
    model.component("comp1").physics("ht").feature().remove("temp1");
    model.component("comp1").physics("ht").create("bhs1", "BoundaryHeatSource", 1);
    model.component("comp1").physics("ht").feature("bhs1").selection().set();
    model.component("comp1").physics("ht").feature().remove("bhs1");
    model.component("comp1").physics("ht").create("temp1", "TemperatureBoundary", 1);
    model.component("comp1").physics("ht").feature("temp1").selection().set(6, 7, 8, 9);
    model.component("comp1").physics("ht").feature("temp1").set("T0", "T0+alpha*t");
    model.component("comp1").physics("ht").feature().remove("temp1");

    model.result("pg1").run();
    model.result().create("pg3", "PlotGroup2D");
    model.result("pg3").run();
    model.result("pg3").create("surf1", "Surface");
    model.result("pg3").run();
    model.result("pg3").run();
    model.result("pg3").run();
    model.result("pg3").feature("surf1").set("colortable", "Magma");
    model.result("pg3").run();
    model.result("pg3").run();
    model.result("pg3").feature("surf1").set("rangecoloractive", true);
    model.result("pg3").feature("surf1").set("rangecolormax", "1e-4");
    model.result("pg3").run();
    model.result("pg3").run();
    model.result("pg3").create("con1", "Contour");
    model.result("pg3").feature("con1").set("number", 10);
    model.result("pg3").feature("con1").set("coloring", "uniform");
    model.result("pg3").feature("con1").set("color", "green");
    model.result("pg3").feature("con1").set("expr", "pf.Vf2");
    model.result("pg3").run();
    model.result("pg2").run();
    model.result("pg3").run();
    model.result("pg3").label("Velocity and Volume Fraction");

    model.component("comp1").physics("ht").create("hs1", "HeatSource", 2);
    model.component("comp1").physics("ht").feature("hs1").selection().set(1, 2);
    model.component("comp1").physics("ht").feature("hs1").set("Q0", "1e10");

    model.result("pg2").run();

    model.component("comp1").physics("ht").feature("hs1").set("Q0", "1e4");

    model.result("pg3").run();
    model.result("pg1").run();
    model.result().table("evl2")
         .addRow(new double[]{7.33624267578125, -6.877727508544922, 342.6304751894849}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-5.387552261352539, -9.39956283569336, 342.63047388707963}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-8.597160339355469, -13.067683219909668, 342.63047296304313}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-19.83078384399414, -16.506547927856445, 342.63047458855334}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-15.245631217956543, -20.747814178466797, 342.6304731399391}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-9.858076095581055, -20.518556594848633, 342.6304711341336}, new double[]{0, 0, 0});
    model.result("pg2").run();
    model.result("pg2").stepPrevious(0);
    model.result("pg2").run();

    model.param().set("deltat", "5 [K]");

    model.component("comp1").physics("ht").feature("hs1").selection().set(2);
    model.component("comp1").physics("ht").feature("hs1").set("Q0", "1e6");

    model.result("pg2").run();
    model.result("pg1").run();
    model.result("pg1").stepPrevious(0);
    model.result("pg1").run();
    model.result("pg1").stepPrevious(0);
    model.result("pg1").run();
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result().table("evl2")
         .addRow(new double[]{3.897378921508789, -16.506547927856445, 319.1053641910389}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{1.948690414428711, -15.131002426147461, 319.1053641232544}, new double[]{0, 0, 0});

    return model;
  }

  public static Model run3(Model model) {
    model.result().table("evl2")
         .addRow(new double[]{-2.7510929107666016, -11.577508926391602, 319.10536357887105}, new double[]{0, 0, 0});
    model.result("pg2").run();
    model.result("pg3").run();
    model.result("pg3").run();
    model.result().table("evl2")
         .addRow(new double[]{-1.6048030853271484, -0.22925758361816406, 6.434224627873536E-10}, new double[]{0, 0, 0});
    model.result("pg3").stepPrevious(0);
    model.result("pg3").run();
    model.result("pg3").stepFirst(0);
    model.result("pg3").run();
    model.result("pg3").run();
    model.result("pg3").feature("surf1").set("rangecoloractive", false);
    model.result("pg1").run();

    model.param().set("kair", "1e-3");

    model.component("comp1").func().create("gp1", "GaussianPulse");
    model.component("comp1").func("gp1").label("diracdelta");
    model.component("comp1").func("gp1").set("funcname", "diracdelta");
    model.component("comp1").func("gp1").set("location", "Tm");
    model.component("comp1").func("gp1").set("sigma", "deltat/2");
    model.component("comp1").func("gp1").set("normalization", "peak");
    model.component("comp1").func("gp1").label("jumpfunc");
    model.component("comp1").func("gp1").set("funcname", "jumpfunc");

    model.param().rename("cpmat", "cpmats");
    model.param().rename("cpmats", "cpmatsolid");
    model.param().set("cpmatliq", "2");
    model.param().set("Hmat", "100");

    model.component("comp1").func("gp1").set("peakvalue", "Hmat");
    model.component("comp1").func().create("step2", "Step");
    model.component("comp1").func("step2").label("stepfunc");
    model.component("comp1").func("step2").set("location", "Tm");
    model.component("comp1").func("step2").set("from", "cpmatsolid");
    model.component("comp1").func("step2").set("to", "cpmatliq");
    model.component("comp1").func("step2").set("smooth", "deltat");
    model.component("comp1").func().create("an1", "Analytic");
    model.component("comp1").func("an1").label("Specific Heat");
    model.component("comp1").func("an1").set("funcname", "cpmat");
    model.component("comp1").func("an1").set("args", "T");
    model.component("comp1").func("an1").set("expr", "jumpfunc(T) + stepfunc(T)");
    model.component("comp1").func("step2").set("funcname", "stepfunc");
    model.component("comp1").func("an1").setIndex("plotargs", 300, 0, 1);
    model.component("comp1").func("an1").setIndex("plotargs", 400, 0, 2);

    model.component("comp1").physics().create("dode", "DomainODE", new String[][]{{"u2"}});

    model.study("std1").feature("phasei").setSolveFor("/physics/dode", true);
    model.study("std1").feature("time").setSolveFor("/physics/dode", true);

    model.component("comp1").physics("dode").prop("EquationForm").set("form", "Automatic");
    model.component("comp1").physics("dode").label("Specific Heat");
    model.component("comp1").physics("dode").tag("cpode");
    model.component("comp1").physics("cpode").feature("init1").set("u2", "cpmat(T0)");
    model.component("comp1").physics("cpode").feature("init1").set("u2t", "d(cpmat(T0), T0)*d(T, t)");
    model.component("comp1").physics("cpode").field("dimensionless").field("cpfield");
    model.component("comp1").physics("cpode").field("dimensionless").component(1, "Cp2");
    model.component("comp1").physics("cpode").create("aleq1", "AlgebraicEquation", 2);
    model.component("comp1").physics("cpode").feature("aleq1").setIndex("f", "Cp2 - cpmat(T)", 0);
    model.component("comp1").physics("cpode").feature("aleq1").selection().set(1, 2);
    model.component("comp1").physics("ht").feature("fluid1").set("Cp", "cpair*(1-phipf)/2 + Cp2*(1+phipf)/2");

    model.result("pg2").run();
    model.result("pg2").stepLast(0);
    model.result("pg2").run();
    model.result("pg3").run();
    model.result("pg1").run();
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result().table("evl2")
         .addRow(new double[]{-7.794757843017578, -19.716154098510742, 315.6366938259374}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-6.648469924926758, -20.403926849365234, 315.6366940286563}, new double[]{0, 0, 0});
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result().table("evl2")
         .addRow(new double[]{-5.043666839599609, -17.88209342956543, 316.12358383063975}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-5.5021820068359375, -18.22597885131836, 316.12358377021775}, new double[]{0, 0, 0});
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result().table("evl2")
         .addRow(new double[]{-4.126638412475586, -17.194320678710938, 316.8317332443574}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-5.387552261352539, -17.4235782623291, 316.83173300975056}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-5.846067428588867, -18.799123764038086, 316.83173303619446}, new double[]{0, 0, 0});
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result().table("evl2")
         .addRow(new double[]{-7.794757843017578, -18.799123764038086, 317.3867792815993}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-6.304582595825195, -19.257638931274414, 317.38677940906337}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-5.616809844970703, -20.6331844329834, 317.38677953764915}, new double[]{0, 0, 0});
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result().table("evl2")
         .addRow(new double[]{-4.355894088745117, -22.008729934692383, 318.2215301707762}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-5.846067428588867, -20.518556594848633, 318.2215298717971}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-5.272924423217773, -18.111351013183594, 318.22152986626145}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-2.5218353271484375, -15.818775177001953, 318.22153001299966}, new double[]{0, 0, 0});
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result().table("evl2")
         .addRow(new double[]{-4.699779510498047, -21.550214767456055, 318.87324039608626}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-5.846067428588867, -20.518556594848633, 318.87324014081406}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-4.126638412475586, -17.5382080078125, 318.8732403159556}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-4.126638412475586, -17.5382080078125, 318.8732403159556}, new double[]{0, 0, 0});
    model.result("pg1").stepNext(0);
    model.result("pg1").run();

    model.param().set("Hmat", "1000");

    model.component("comp1").func("gp1").set("sigma", "deltat");

    model.sol("sol1").runAll();

    model.result("pg1").run();
    model.result("pg1").stepLast(0);
    model.result("pg1").run();
    model.result().table("evl2")
         .addRow(new double[]{-6.648469924926758, -18.799123764038086, 319.8329447245952}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-7.794757843017578, -18.569866180419922, 319.8329447000703}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-6.533840179443359, -15.131002426147461, 319.83294468411674}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-6.419212341308594, -15.016373634338379, 319.8329446842116}, new double[]{0, 0, 0});
    model.result("pg2").run();
    model.result("pg3").run();
    model.result("pg3").stepLast(0);
    model.result("pg3").run();
    model.result("pg2").run();
    model.result("pg1").run();
    model.result().create("pg4", "PlotGroup2D");
    model.result("pg4").run();
    model.result("pg4").label("Viscosity");
    model.result("pg4").create("surf1", "Surface");
    model.result("pg4").feature("surf1").set("expr", "comp1.mu2");
    model.result("pg4").run();
    model.result().table("evl2")
         .addRow(new double[]{2.3863658905029297, 3.409090042114258, 562307.1884163866}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-6.477272033691406, 7.386363983154297, 562307.2360032706}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-9.318180084228516, 14.09090805053711, 562307.2879593191}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-9.545454025268555, 15.454544067382812, 562307.2951086677}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-9.545454025268555, 15.454544067382812, 562307.2951086677}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{0.22727394104003906, -19.772727966308594, 562307.1665807336}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-3.863636016845703, -20, 562306.8213687898}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-3.1818180084228516, -22.727272033691406, 562306.7465962179}, new double[]{0, 0, 0});
    model.result("pg4").run();
    model.result("pg4").run();
    model.result("pg4").feature("surf1").set("expr", "muair*(1-phipf)/2 + comp1.mu2*(1+phipf)/2");
    model.result("pg4").run();
    model.result().table("evl2")
         .addRow(new double[]{-3.5227279663085938, -23.75, 244698.01028505026}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-4.204545974731445, -22.954545974731445, 249207.9966030543}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-5.113636016845703, -21.136363983154297, 232420.3932244964}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-4.886363983154297, -18.522727966308594, 257211.19499716788}, new double[]{0, 0, 0});
    model.result("pg4").feature("surf1").set("coloring", "gradient");
    model.result("pg4").feature("surf1").set("topcolor", "blue");
    model.result("pg4").feature("surf1").set("bottomcolor", "white");
    model.result("pg4").feature("surf1").set("topcolor", "magenta");
    model.result("pg4").run();
    model.result("pg4").create("con1", "Contour");
    model.result("pg4").feature("con1").set("expr", "pf.Vf2");
    model.result("pg4").run();
    model.result("pg4").feature("con1").stepPrevious(0);
    model.result("pg4").run();
    model.result("pg4").feature("con1").stepFirst(0);
    model.result("pg4").run();
    model.result("pg1").run();
    model.result("pg4").run();

    model.component("comp1").physics("ht").feature("hs1").set("Q0", "1e7");

    model.result("pg3").run();
    model.result("pg2").run();
    model.result("pg1").run();
    model.result("pg2").run();
    model.result("pg2").stepLast(0);
    model.result("pg2").run();
    model.result("pg2").stepPrevious(0);
    model.result("pg2").run();
    model.result("pg2").stepPrevious(0);
    model.result("pg2").run();
    model.result("pg1").run();
    model.result("pg1").stepFirst(0);
    model.result("pg1").run();
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result("pg4").run();
    model.result("pg4").stepLast(0);
    model.result("pg4").run();
    model.result("pg1").run();
    model.result().table("evl2")
         .addRow(new double[]{-5.846067428588867, -17.76746368408203, 310.1111595181136}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-4.241266250610352, -16.048032760620117, 310.1111598160689}, new double[]{0, 0, 0});
    model.result("pg1").stepLast(0);
    model.result("pg1").run();
    model.result().table("evl2")
         .addRow(new double[]{-5.846067428588867, -20.060039520263672, 322.5057767524697}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-6.18995475769043, -19.601524353027344, 322.50578427934477}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-3.5534934997558594, -18.455238342285156, 322.5058590974268}, new double[]{0, 0, 0});
    model.result("pg4").run();
    model.result("pg4").stepFirst(0);
    model.result("pg4").run();
    model.result("pg4").stepNext(0);
    model.result("pg4").run();
    model.result("pg4").stepNext(0);
    model.result("pg4").run();
    model.result("pg4").stepNext(0);
    model.result("pg4").run();
    model.result("pg4").stepNext(0);
    model.result("pg4").run();
    model.result("pg4").stepNext(0);
    model.result("pg4").run();

    model.study("std1").feature("time").set("tlist", "range(0,0.01,1)");

    model.label("ThermalNewtonianIncompressibleODE.mph");

    model.study("std1").run(false);

    model.label("ThermalNewtonianIncompressibleODE.mph");

    model.result("pg2").run();
    model.result("pg2").stepLast(0);
    model.result("pg2").run();
    model.result("pg1").run();
    model.result("pg1").stepLast(0);
    model.result("pg1").run();
    model.result("pg1").stepFirst(0);
    model.result("pg1").run();
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result().export().create("tbl1", "Table");
    model.result().export().remove("tbl1");
    model.result("pg2").run();
    model.result().export().create("anim1", "Animation");
    model.result().export("anim1").set("fontsize", "9");
    model.result().export("anim1").set("colortheme", "globaltheme");
    model.result().export("anim1").set("customcolor", new double[]{1, 1, 1});
    model.result().export("anim1").set("background", "color");
    model.result().export("anim1").set("gltfincludelines", "on");
    model.result().export("anim1").set("title1d", "on");
    model.result().export("anim1").set("legend1d", "on");
    model.result().export("anim1").set("logo1d", "on");
    model.result().export("anim1").set("options1d", "on");
    model.result().export("anim1").set("title2d", "on");
    model.result().export("anim1").set("legend2d", "on");
    model.result().export("anim1").set("logo2d", "on");
    model.result().export("anim1").set("options2d", "on");
    model.result().export("anim1").set("title3d", "on");
    model.result().export("anim1").set("legend3d", "on");
    model.result().export("anim1").set("logo3d", "on");
    model.result().export("anim1").set("options3d", "off");
    model.result().export("anim1").set("axisorientation", "on");
    model.result().export("anim1").set("grid", "on");
    model.result().export("anim1").set("axes1d", "on");
    model.result().export("anim1").set("axes2d", "on");
    model.result().export("anim1").set("showgrid", "on");
    model.result().export("anim1").label("T and Vf2");
    model.result("pg1").run();
    model.result("pg1").set("titletype", "manual");
    model.result("pg1").set("title", "Temperature (K) with Volume Fraction Contours");
    model.result("pg1").run();
    model.result("pg1").set("showlegendsmaxmin", true);
    model.result("pg1").set("showlegendsunit", true);
    model.result("pg1").run();
    model.result().export("anim1").set("framesel", "all");
    model.result().export().duplicate("anim2", "anim1");
    model.result().export("anim2").label("p and Vf2");
    model.result().export("anim2").set("plotgroup", "pg2");
    model.result().export().duplicate("anim3", "anim2");
    model.result().export("anim3").label("U and Vf2");
    model.result().export("anim3").set("plotgroup", "pg3");
    model.result().export().duplicate("anim4", "anim3");
    model.result().export("anim4").label("mu and Vf2");
    model.result().export("anim4")
         .set("giffilename", "/nas/longleaf/home/mshourya/workspace/SelectiveLaserSintering/Particle/Results/TemperatureDependent/mu+Vf2.gif");
    model.result().export("anim3")
         .set("giffilename", "/nas/longleaf/home/mshourya/workspace/SelectiveLaserSintering/Particle/Results/TemperatureDependent/U+Vf2.gif");
    model.result().export("anim2")
         .set("giffilename", "/nas/longleaf/home/mshourya/workspace/SelectiveLaserSintering/Particle/Results/TemperatureDependent/Vf2+p.gif");
    model.result().export("anim1")
         .set("giffilename", "/nas/longleaf/home/mshourya/workspace/SelectiveLaserSintering/Particle/Results/TemperatureDependent/T+Vf2.gif");
    model.result("pg1").run();
    model.result("pg2").run();
    model.result("pg2").set("titletype", "manual");
    model.result("pg2").set("title", "Volume Fraction with Pressure (Pa) contours");
    model.result("pg2").run();
    model.result("pg3").run();
    model.result("pg3").set("titletype", "manual");
    model.result("pg3").set("title", "Velocity magnitude (m/s) with Volume Fraction contours");
    model.result("pg3").run();
    model.result("pg4").run();
    model.result("pg4").set("titletype", "manual");
    model.result("pg4").set("title", "Fluid viscosity (Pa*s) with Volume Fraction contours");
    model.result("pg4").run();
    model.result("pg1").run();
    model.result("pg2").run();
    model.result("pg1").run();
    model.result("pg1").stepLast(0);
    model.result("pg1").run();
    model.result("pg2").run();
    model.result("pg4").run();
    model.result("pg4").stepLast(0);
    model.result("pg4").run();
    model.result("pg4").stepFirst(0);
    model.result("pg4").run();
    model.result("pg4").run();
    model.result("pg4").feature("con1").set("number", 10);
    model.result("pg4").run();
    model.result("pg2").run();
    model.result("pg2").run();
    model.result("pg1").run();
    model.result("pg2").run();
    model.result("pg3").run();
    model.result("pg4").run();
    model.result("pg1").run();
    model.result().export("anim1").run();
    model.result().export("anim2").run();
    model.result().export("anim3").run();
    model.result().export("anim4").set("plotgroup", "pg4");
    model.result().export("anim4").run();
    model.result("pg4").run();

    model.component("comp1").func("step1").createPlot("pg5");

    model.result("pg5").run();
    model.result("pg5").label("Material Viscosity");
    model.result("pg5").set("title", "Viscosity of Material (Pa*s)");
    model.result("pg5").set("ylabelactive", true);
    model.result("pg5").set("ylabel", "Pa*s");
    model.result("pg5").set("xlabelactive", true);
    model.result("pg5").set("xlabel", "T (K)");
    model.result("pg5").run();

    model.component("comp1").func("an1").createPlot("pg6");

    model.result("pg6").run();
    model.result("pg6").label("Specific Heat of Material");
    model.result("pg6").set("title", "Specific Heat of Material ");
    model.result("pg6").run();
    model.result("pg6").set("title", "Specific Heat of Material (J/(kg*K))");
    model.result("pg6").set("xlabelactive", true);
    model.result("pg6").set("xlabel", "T (K)");
    model.result("pg6").set("ylabel", "J/(kg*K)");
    model.result("pg6").run();

    model.param().set("cpair", "1 [J/(kg*K)]");
    model.param().set("cpmatsolid", "1 [J/(kg*K)]");
    model.param().set("cpmatliq", "2 [J/(kg*K)]");

    model.component("comp1").physics("cpode").prop("Units").set("CustomDependentVariableUnit", "1");
    model.component("comp1").physics("cpode").prop("Units").set("DependentVariableQuantity", "none");
    model.component("comp1").physics("cpode").prop("Units").setIndex("CustomDependentVariableUnit", "J/(kg*K)", 0, 0);

    model.param().set("kair", "1e-3 [W/(m*K)]");
    model.param().set("kmat", "0.6 [W/(m*K)]");
    model.param().set("kair", "0.03 [W/(m*K)]");
    model.param().set("cpair", "1 [J/(kg*K)]");
    model.param().set("cpmatsolid", "2098 [J/(kg*K)]");
    model.param().set("cpmatliq", "1.484 [J/(kg*K)]");
    model.param().set("deltat", "1 [K]");

    model.component("comp1").physics("ht").feature("hs1").set("Q0", "Q");

    model.param().set("Q", "1e8");
    model.param().set("Hmat", "4000");

    model.result("pg2").run();
    model.result("pg2").stepLast(0);
    model.result("pg2").run();
    model.result("pg2").stepPrevious(0);
    model.result("pg2").run();
    model.result("pg2").stepLast(0);
    model.result("pg2").run();
    model.result("pg2").set("paramindicator", "");
    model.result("pg2").run();
    model.result("pg3").run();
    model.result("pg1").run();
    model.result("pg1").set("paramindicator", "t=eval(t,s)");
    model.result("pg1").run();
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result("pg2").run();
    model.result("pg2").set("paramindicator", "t=eval(t,s)");
    model.result("pg2").run();
    model.result("pg3").run();
    model.result("pg3").set("paramindicator", "t=eval(t)");
    model.result("pg3").run();
    model.result("pg3").stepLast(0);
    model.result("pg3").run();
    model.result("pg4").run();
    model.result("pg4").set("paramindicator", "t=eval(t,s)");
    model.result("pg4").run();
    model.result("pg4").set("paramindicator", "t=eval(t,s) (s)");
    model.result("pg3").run();
    model.result("pg3").set("paramindicator", "t=eval(t) (s)");
    model.result("pg2").run();
    model.result("pg2").set("paramindicator", "t=eval(t,s) (s)");
    model.result("pg1").run();
    model.result("pg1").set("paramindicator", "t=eval(t,s) (s)");

    model.study("std1").feature("time").set("tlist", "range(0,0.001,0.5)");

    model.label("ThermalNewtonianIncompressibleODE.mph");

    model.study("std1").run(false);

    model.label("ThermalNewtonianIncompressibleODE.mph");

    model.result("pg4").run();
    model.result("pg4").stepNext(0);
    model.result("pg4").run();
    model.result("pg4").stepNext(0);
    model.result("pg4").run();
    model.result("pg4").stepNext(0);
    model.result("pg4").run();
    model.result("pg4").stepLast(0);
    model.result("pg4").run();
    model.result("pg4").stepFirst(0);
    model.result("pg4").run();
    model.result("pg3").run();
    model.result("pg3").stepNext(0);
    model.result("pg3").run();
    model.result("pg3").setIndex("looplevel", 42, 0);
    model.result("pg3").stepNext(0);
    model.result("pg3").run();
    model.result("pg3").stepLast(0);
    model.result("pg3").run();
    model.result().export("anim4")
         .set("giffilename", "/nas/longleaf/home/mshourya/workspace/SelectiveLaserSintering/Particle/Results/PhysicalTemperatureDependent/mu+Vf2.gif");
    model.result().export("anim3")
         .set("giffilename", "/nas/longleaf/home/mshourya/workspace/SelectiveLaserSintering/Particle/Results/PhysicalTemperatureDependent/U+Vf2.gif");
    model.result().export("anim2")
         .set("giffilename", "/nas/longleaf/home/mshourya/workspace/SelectiveLaserSintering/Particle/Results/PhysicalTemperatureDependent/Vf2+p.gif");
    model.result().export("anim1")
         .set("giffilename", "/nas/longleaf/home/mshourya/workspace/SelectiveLaserSintering/Particle/Results/PhysicalTemperatureDependent/T+Vf2.gif");
    model.result().export("anim1").run();
    model.result().export("anim4").run();

    model.label("ThermalNewtonianIncompressibleODE.mph");

    model.result("pg1").run();
    model.result("pg2").run();
    model.result("pg1").run();
    model.result("pg1").stepNext(0);
    model.result("pg1").run();
    model.result("pg1").stepNext(0);
    model.result("pg1").run();

    model.component("comp1").multiphysics().remove("emh1");

    model.result().param().clear();

    model.component("comp1").probe().create("dom1", "Domain");
    model.component("comp1").probe("dom1").set("intsurface", true);
    model.component("comp1").probe("dom1").set("intvolume", true);
    model.component("comp1").probe("dom1").selection().set(2);
    model.component("comp1").probe().remove("dom1");

    model.component("comp1").geom("geom2").run("c2");
    model.component("comp1").geom("geom2").create("pt1", "Point");
    model.component("comp1").geom("geom2").feature("pt1").set("p", new double[]{0, -20});
    model.component("comp1").geom("geom2").run("pt1");

    model.component("comp1").probe().create("point1", "Point");

    model.component("comp1").geom("geom2").run();

    model.component("comp1").probe("point1").selection().set(5);
    model.component("comp1").probe("point1").set("expr", "ht.T");
    model.component("comp1").probe("point1").genResult("sol1");

    model.result().numerical("pev1").set("table", "tbl1");
    model.result().numerical("pev1").set("innerinput", "all");
    model.result().numerical("pev1").set("outerinput", "all");
    model.result().numerical("pev1").setResult();
    model.result("pg7").feature("tblp1").set("plotcolumns", new int[]{2});
    model.result("pg7").feature("tblp1").set("xaxisdata", "auto");
    model.result("pg7").run();

    model.label("ThermalNewtonianIncompressibleODE.mph");

    model.param().set("cpmatliq", "1484 [J/(kg*K)]");

    model.label("ThermalNewtonianIncompressibleODE.mph");

    model.study("std1").run(false);

    model.label("ThermalNewtonianIncompressibleODE.mph");

    model.result("pg1").run();
    model.result("pg3").run();
    model.result().export("anim2").set("framesel", "number");
    model.result().export("anim2").set("maxframes", 201);
    model.result().export("anim2").run();
    model.result().export("anim1").set("framesel", "number");
    model.result().export("anim1").set("maxframes", 201);
    model.result().export("anim1").run();
    model.result().export("anim4").set("framesel", "number");
    model.result().export("anim4").set("maxframes", 201);
    model.result().export("anim4").run();

    model.param().set("alpha", "0");
    model.param().remove("alpha");
    model.param().set("Q", "1e8 [W/m^3]");
    model.param().remove("OtherQuantities");

    model.result("pg6").run();
    model.result("pg7").set("window", "window2");
    model.result("pg7").run();
    model.result("pg7").set("window", "window2");
    model.result("pg7").run();
    model.result("pg7").set("window", "window2");
    model.result("pg7").run();
    model.result("pg7").set("xlabelactive", true);
    model.result("pg7").set("xlabel", "t (s)");
    model.result("pg7").set("ylabelactive", true);
    model.result("pg7").set("ylabel", "T (K)");
    model.result("pg7").set("window", "window2");
    model.result("pg7").run();
    model.result("pg7").set("showlegends", false);
    model.result("pg7").set("window", "window2");
    model.result("pg7").run();
    model.result("pg7").set("titletype", "manual");
    model.result("pg7").set("title", "Temperature Probe");
    model.result("pg7").set("window", "window2");
    model.result("pg7").run();
    model.result("pg7").label("Temperature Probe");

    model.label("ThermalNewtonianIncompressibleODE.mph");

    model.result("pg7").set("window", "window2");
    model.result("pg7").run();

    model.component("comp1").mesh("mesh2").feature("ftri2").feature("size1").set("hmax", 1);
    model.component("comp1").mesh("mesh2").run();

    model.component("comp1").probe("point1").genResult("none");

    model.sol("sol1").runAll();

    model.result("pg1").run();
    model.result("pg1").stepLast(0);
    model.result("pg1").run();
    model.result("pg3").run();
    model.result("pg3").stepPrevious(0);
    model.result("pg3").run();
    model.result("pg3").stepPrevious(0);
    model.result("pg3").run();
    model.result("pg3").stepPrevious(0);
    model.result("pg3").run();
    model.result("pg3").stepPrevious(0);
    model.result("pg3").run();
    model.result("pg3").stepPrevious(0);
    model.result("pg3").run();
    model.result("pg3").stepPrevious(0);
    model.result("pg3").run();
    model.result("pg3").stepPrevious(0);
    model.result("pg3").run();
    model.result("pg3").stepPrevious(0);
    model.result("pg3").run();
    model.result("pg4").run();
    model.result("pg4").stepLast(0);
    model.result("pg4").run();
    model.result().export("anim4").set("target", "player");
    model.result().export("anim4").set("stopped", true);
    model.result().export("anim4").run();

    model.label("ThermalNewtonianIncompressibleODE.mph");

    model.result().export("anim4").showFrame();

    model.label("ThermalNewtonianIncompressibleODE.mph");

    model.nodeGroup().create("grp1", "Definitions", "comp1");
    model.nodeGroup("grp1").set("type", "func");
    model.nodeGroup("grp1").add("func", "step1");
    model.nodeGroup("grp1").add("func", "gp1");
    model.nodeGroup("grp1").remove("func", "gp1", false);
    model.nodeGroup("grp1").add("func", "gp1");
    model.nodeGroup("grp1").remove("func", "gp1", false);
    model.nodeGroup("grp1").add("func", "an1");
    model.nodeGroup("grp1").remove("func", "an1", false);
    model.nodeGroup("grp1").label("Materials Properties");

    model.label("ThermalNewtonianIncompressibleODE.mph");

    model.component("comp1").probe("point1").genResult("none");

    model.component("comp1").physics("spf").create("constr2", "PointwiseConstraint", 2);
    model.component("comp1").physics("spf").feature("constr2").selection().set(1, 2);
    model.component("comp1").physics("spf").feature("constr2").set("constraintExpression", "u*(1-phipf)/2");

    model.component("comp1").func().create("step3", "Step");
    model.component("comp1").func("step3").label("step");
    model.component("comp1").func("step3").set("funcname", "step");
    model.component("comp1").func("step3").set("from", 1);
    model.component("comp1").func("step3").set("to", 0);
    model.component("comp1").func("step3").set("location", "Tm");
    model.component("comp1").func("step3").set("smooth", "deltat");

    model.component("comp1").physics("spf").feature("constr2").set("constraintExpression", "u*step(T)");
    model.component("comp1").physics("spf").feature("constr2").label("Solid Constraint");

    model.component("comp1").probe("point1").genResult("none");

    model.component("comp1").physics("spf").feature("constr2").set("constraintExpression", "u*step(T)*(1+phi)/2");

    model.component("comp1").probe("point1").genResult("none");

    model.component("comp1").physics("spf").feature("constr2").set("constraintExpression", "u*step(T)*(1+phipf)/2");

    model.component("comp1").probe("point1").genResult("none");

    model.component("comp1").physics("spf").feature("constr2").set("constraintType", "symmetricConstraint");
    model.component("comp1").physics("spf").feature().remove("constr2");
    model.component("comp1").physics("spf").feature().remove("weak2");
    model.component("comp1").physics("spf").feature().remove("weak1");
    model.component("comp1").physics("ht").feature("hs1").selection().set(1, 2);

    return model;
  }

  public static Model run4(Model model) {
    model.component("comp1").physics("ht").feature("hs1").set("Q0", "Q*(y+25)");

    model.component("comp1").probe("point1").genResult("none");

    model.param().set("Q", "1 [W/m^3]");

    model.component("comp1").probe("point1").genResult("none");

    model.result("pg1").run();
    model.result("pg1").stepFirst(0);
    model.result("pg1").run();
    model.result().table("evl2")
         .addRow(new double[]{-7.1707305908203125, 11.65243911743164, 300.00000000002126}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{1.664632797241211, 6.146343231201172, 300.00000000002126}, new double[]{0, 0, 0});
    model.result("pg1").stepLast(0);
    model.result("pg1").run();
    model.result().table("evl2")
         .addRow(new double[]{-3.2012195587158203, -0.38414573669433594, 300.00016684979227}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-6.530487060546875, 1.7926826477050781, 300.0001669124788}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{-9.987804412841797, 3.8414649963378906, 300.00016696500666}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{5.76219367980957, 6.0182952880859375, 300.00016699429585}, new double[]{0, 0, 0});
    model.result().table("evl2")
         .addRow(new double[]{9.347558975219727, 1.280487060546875, 300.0001669120045}, new double[]{0, 0, 0});

    model.component("comp1").physics("ht").feature("hs1").set("heatSourceType", "GeneralSource");

    model.label("ThermalNewtonianIncompressibleODE.mph");

    model.result("pg7").set("window", "window2");
    model.result("pg7").run();
    model.result("pg6").run();
    model.result("pg7").set("window", "window2");
    model.result("pg7").run();
    model.result("pg5").run();
    model.result("pg4").run();
    model.result("pg4").stepFirst(0);
    model.result("pg4").run();
    model.result("pg3").run();
    model.result("pg3").set("data", "dset1");
    model.result("pg3").stepFirst(0);
    model.result("pg3").run();
    model.result("pg4").run();
    model.result("pg4").stepLast(0);
    model.result("pg4").run();

    model.component("comp1").probe("point1").genResult("none");

    model.sol("sol1").runAll();

    model.result("pg1").run();
    model.result("pg3").run();
    model.result("pg3").stepLast(0);
    model.result("pg3").run();
    model.result("pg3").stepFirst(0);
    model.result("pg3").run();
    model.result("pg3").stepLast(0);
    model.result("pg3").run();

    model.component("comp1").physics("ht").feature("hs1").selection().set(2);
    model.component("comp1").physics("ht").feature("hs1").set("Q0", "Q");

    model.component("comp1").probe("point1").genResult("none");

    model.param().set("Q", "10 [W/m^3]");

    model.component("comp1").probe("point1").genResult("none");

    model.param().set("Q", "1000 [W/m^3]");

    model.component("comp1").probe("point1").genResult("none");

    model.param().set("Q", "1e9 [W/m^3]");

    model.component("comp1").probe("point1").genResult("none");

    model.label("ThermalNewtonianIncompressibleODE.mph");
    model.label("ThermalNewtonianIncompressibleODE.mph");

    return model;
  }

  public static void main(String[] args) {
    Model model = run();
    model = run2(model);
    model = run3(model);
    run4(model);
  }

}
