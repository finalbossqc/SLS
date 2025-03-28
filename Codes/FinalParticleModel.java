/*
 * FinalParticleModel.java
 */

import com.comsol.model.*;
import com.comsol.model.util.*;

/** Model exported on Feb 11 2025, 22:50 by COMSOL 6.1.0.357. */
public class FinalParticleModel {

  public static Model run() {
    Model model = ModelUtil.create("Model");

    model.modelPath("/nas/longleaf/home/mshourya/workspace/SLS/Particle");

    model.label("FinalParticleModel.mph");

    model.param().set("MaterialQuantities", "0", "___________________________________________________");
    model.param().set("rhomat", "1000 [kg/m^3]", "Material density");
    model.param().set("mumatsolid", "1e5 [Pa*s]", "Material viscosity in solid phase");
    model.param().set("mumatliq", "1 [Pa*s]", "Material viscosity in liquid phase");
    model.param().set("TemperatureQuantities", "0", "___________________________________________________");
    model.param().set("kmat", "0.6 [W/(m*K)]", "Material thermal conductivity");
    model.param().set("cpair", "1 [J/(kg*K)]", "Air specific heat");
    model.param().set("cpmatsolid", "2098 [J/(kg*K)]", "Material specific heat solid");
    model.param().set("cpmatliq", "1484 [J/(kg*K)]", "Material specific heat liquid");
    model.param().set("Hmat", "4000", "Peak Cp");
    model.param().set("deltat", "1 [K]", "Phase transition regime size");
    model.param().set("alpha", "1e-6 [1/K]", "Thermal expansion");
    model.param().set("T0", "300[K]", "Initial temperature");
    model.param().set("Tm", "320 [K]", "Melting temperature");
    model.param().set("Tamb", "300 [K]", "Ambient temperature");
    model.param().set("B", "0.1");
    model.param().set("epsilon", "1.5e-5 [m]", "interface thickness");
    model.param().set("g", "9.81", "Gravitational acceleration");
    model.param().set("sigma", "0.1 [N/m]");
    model.param().set("LaserParameters", "0", "_______________________");
    model.param().set("vlaser", "3 [mm/s]", "Laser speed");
    model.param().set("Ep", "0.1 [J/s]", "Laser power");
    model.param().set("Pw", "100 [um]");
    model.param().set("D", "100 [um]");
    model.param().set("A", "0.12");
    model.param().set("xr", "-50 [um]");
    model.param().set("yr", "200 [um]");
    model.param().set("xd", "100 [um]");
    model.param().set("ontime", "0.02 [s]");
    model.param().set("offtime", "0.16 [s]");

    model.component().create("comp1", true);

    model.component("comp1").geom().create("geom2", 2);

    model.result().table().create("evl2", "Table");
    model.result().table().create("tbl1", "Table");

    model.component("comp1").func().create("an1", "Analytic");
    model.component("comp1").func().create("gp1", "GaussianPulse");
    model.component("comp1").func().create("step2", "Step");
    model.component("comp1").func().create("step3", "Step");
    model.component("comp1").func().create("rect1", "Rectangle");
    model.component("comp1").func().create("an2", "Analytic");
    model.component("comp1").func("an1").label("Specific Heat");
    model.component("comp1").func("an1").set("funcname", "cpmat");
    model.component("comp1").func("an1").set("expr", "jumpfunc(T) + stepfunc(T)");
    model.component("comp1").func("an1").set("args", new String[]{"T"});
    model.component("comp1").func("an1").set("fununit", "J/(kg*K)");
    model.component("comp1").func("an1").set("argunit", new String[]{"K"});
    model.component("comp1").func("an1").set("plotargs", new String[][]{{"T", "300", "400"}});
    model.component("comp1").func("gp1").label("jumpfunc");
    model.component("comp1").func("gp1").set("funcname", "jumpfunc");
    model.component("comp1").func("gp1").set("location", "Tm");
    model.component("comp1").func("gp1").set("sigma", "deltat");
    model.component("comp1").func("gp1").set("normalization", "peak");
    model.component("comp1").func("gp1").set("peakvalue", "Hmat");
    model.component("comp1").func("step2").label("stepfunc");
    model.component("comp1").func("step2").set("funcname", "stepfunc");
    model.component("comp1").func("step2").set("location", "Tm");
    model.component("comp1").func("step2").set("from", "cpmatsolid");
    model.component("comp1").func("step2").set("to", "cpmatliq");
    model.component("comp1").func("step2").set("smooth", "deltat");
    model.component("comp1").func("step3").label("step");
    model.component("comp1").func("step3").set("funcname", "step");
    model.component("comp1").func("step3").set("location", "Tm");
    model.component("comp1").func("step3").set("from", 1);
    model.component("comp1").func("step3").set("to", 0);
    model.component("comp1").func("step3").set("smooth", "deltat");
    model.component("comp1").func("rect1").label("onoff");
    model.component("comp1").func("rect1").set("funcname", "onoff");
    model.component("comp1").func("rect1").set("lower", "ontime");
    model.component("comp1").func("rect1").set("upper", "offtime");
    model.component("comp1").func("rect1").set("smooth", 0.001);
    model.component("comp1").func("an2").label("mumat");
    model.component("comp1").func("an2").set("funcname", "mumat");
    model.component("comp1").func("an2").set("expr", "mumatliq+exp(-1*(T-Tm))");
    model.component("comp1").func("an2").set("args", new String[]{"T"});
    model.component("comp1").func("an2").set("argunit", new String[]{"K"});
    model.component("comp1").func("an2").set("plotargs", new String[][]{{"T", "323", "324"}});

    model.component("comp1").mesh().create("mesh3");

    model.component("comp1").geom("geom2").label("Geometry 1");
    model.component("comp1").geom("geom2").lengthUnit("\u00b5m");
    model.component("comp1").geom("geom2").create("c1", "Circle");
    model.component("comp1").geom("geom2").feature("c1").set("pos", new int[]{-45, -15});
    model.component("comp1").geom("geom2").feature("c1").set("r", 50);
    model.component("comp1").geom("geom2").create("c2", "Circle");
    model.component("comp1").geom("geom2").feature("c2").active(false);
    model.component("comp1").geom("geom2").feature("c2").set("pos", new int[]{-10, -15});
    model.component("comp1").geom("geom2").feature("c2").set("r", 5);
    model.component("comp1").geom("geom2").create("pt1", "Point");
    model.component("comp1").geom("geom2").feature("pt1").set("p", new int[]{0, -20});
    model.component("comp1").geom("geom2").create("arr1", "Array");
    model.component("comp1").geom("geom2").feature("arr1").set("fullsize", new int[]{10, 3});
    model.component("comp1").geom("geom2").feature("arr1").set("displ", new int[]{100, 100});
    model.component("comp1").geom("geom2").feature("arr1").selection("input").set("c1");
    model.component("comp1").geom("geom2").create("r1", "Rectangle");
    model.component("comp1").geom("geom2").feature("r1").set("pos", new int[]{-100, -70});
    model.component("comp1").geom("geom2").feature("r1").set("size", new int[]{1010, 380});
    model.component("comp1").geom("geom2").run();

    model.component("comp1").variable().create("var1");
    model.component("comp1").variable("var1").set("Ed", "Ep / (Pw*(pi*0.25*D^2))");
    model.component("comp1").variable("var1")
         .set("G_space", "exp(-1* ( (x  - (xr+vlaser*t) )^2 + (y - yr)^2 ) / (2 * xd^2 ) )");
    model.component("comp1").variable("var1").set("Qlaser", "A*Ed*G_space*onoff(t) * (1 + phipf ) / 2");
    model.component("comp1").variable("var1").set("Qsb", "B*(Tamb^4 - T^4)");
    model.component("comp1").variable("var1").set("mu2", "mumat(T)", "Dynamic viscosity of fluid 2");
    model.component("comp1").variable("var1").set("U", "sqrt(u^2 + v^2)", "Velocity magnitude");
    model.component("comp1").variable("var1").set("k", "mat1.def.k_iso*(1-phipf)/2 + kmat*(1+phipf)/2");
    model.component("comp1").variable("var1")
         .set("rho", "mat1.def.rho(p, T)*(1-phipf)/2 + rhomat*(1+phipf)/2", "Density");
    model.component("comp1").variable("var1")
         .set("Cp", "mat1.def.Cp(T)*(1-phipf)/2 + cpmat(T)*(1+phipf)/2", "Heat capacity at constant pressure");
    model.component("comp1").variable("var1")
         .set("mu", "mat1.def.mu*(1-phipf)/2 + mu2*(1+phipf)/2", "Dynamic viscosity");
    model.component("comp1").variable("var1").set("lambda", "3*epsilon*sigma/sqrt(8)");
    model.component("comp1").variable("var1").set("Fstx", "lambda/(epsilon^2) * psi * phipfx");
    model.component("comp1").variable("var1").set("Fsty", "lambda/(epsilon^2) * psi * phipfy");

    model.component("comp1").material().create("mat1", "Common");
    model.component("comp1").material("mat1").propertyGroup("def").func().create("eta", "Piecewise");
    model.component("comp1").material("mat1").propertyGroup("def").func().create("Cp", "Piecewise");
    model.component("comp1").material("mat1").propertyGroup("def").func().create("rho", "Analytic");
    model.component("comp1").material("mat1").propertyGroup("def").func().create("k", "Piecewise");
    model.component("comp1").material("mat1").propertyGroup("def").func().create("cs", "Analytic");
    model.component("comp1").material("mat1").propertyGroup("def").func().create("an1", "Analytic");
    model.component("comp1").material("mat1").propertyGroup("def").func().create("an2", "Analytic");
    model.component("comp1").material("mat1").propertyGroup().create("NonlinearModel", "Nonlinear model");
    model.component("comp1").material("mat1").propertyGroup().create("idealGas", "Ideal gas");
    model.component("comp1").material("mat1").propertyGroup("idealGas").func().create("Cp", "Piecewise");

    model.component("comp1").physics().create("spf", "LaminarFlow", "geom2");
    model.component("comp1").physics("spf").create("out1", "OutletBoundary", 1);
    model.component("comp1").physics("spf").feature("out1").selection().set(3);
    model.component("comp1").physics("spf").create("vf1", "VolumeForce", 2);
    model.component("comp1").physics("spf").feature("vf1").selection().all();
    model.component("comp1").physics("spf").create("vf2", "VolumeForce", 2);
    model.component("comp1").physics("spf").feature("vf2").selection().all();
    model.component("comp1").physics("spf").create("inl1", "InletBoundary", 1);
    model.component("comp1").physics("spf").feature("inl1").selection().set(1, 4);
    model.component("comp1").physics().create("pf", "PhaseField", "geom2");
    model.component("comp1").physics("pf").feature("initfluid2").selection()
         .set(2, 3, 4, 7, 8, 9, 12, 13, 14, 17, 18, 19, 22, 23, 24, 27, 28, 29, 32, 33, 34, 37, 38, 39, 42, 43, 44, 47, 48, 49);
    model.component("comp1").physics().create("ht", "HeatTransferInFluids", "geom2");
    model.component("comp1").physics("ht").create("hs1", "HeatSource", 2);
    model.component("comp1").physics("ht").feature("hs1").selection()
         .set(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49);
    model.component("comp1").physics("ht").create("hs2", "HeatSource", 2);
    model.component("comp1").physics("ht").feature("hs2").selection()
         .set(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49);

    model.component("comp1").mesh("mesh3").create("ftri1", "FreeTri");
    model.component("comp1").mesh("mesh3").feature("ftri1").create("size1", "Size");

    model.component("comp1").probe().create("point1", "Point");
    model.component("comp1").probe().create("point2", "Point");
    model.component("comp1").probe().create("point3", "Point");
    model.component("comp1").probe("point1").selection().set(27);
    model.component("comp1").probe("point2").selection().set(26);
    model.component("comp1").probe("point3").selection().set(34);

    model.result().table("evl2").label("Evaluation 2D");
    model.result().table("evl2").comments("Interactive 2D values");
    model.result().table("tbl1").label("Probe Table 1");

    model.component("comp1").view("view1").axis().set("xmin", -150.5);
    model.component("comp1").view("view1").axis().set("xmax", 960.5);
    model.component("comp1").view("view1").axis().set("ymin", -185.43170166015625);
    model.component("comp1").view("view1").axis().set("ymax", 425.43170166015625);

    model.component("comp1").material("mat1").label("Air");
    model.component("comp1").material("mat1").set("family", "air");
    model.component("comp1").material("mat1").propertyGroup("def").func("eta").set("arg", "T");
    model.component("comp1").material("mat1").propertyGroup("def").func("eta")
         .set("pieces", new String[][]{{"200.0", "1600.0", "-8.38278E-7+8.35717342E-8*T^1-7.69429583E-11*T^2+4.6437266E-14*T^3-1.06585607E-17*T^4"}});
    model.component("comp1").material("mat1").propertyGroup("def").func("eta").set("argunit", "K");
    model.component("comp1").material("mat1").propertyGroup("def").func("eta").set("fununit", "Pa*s");
    model.component("comp1").material("mat1").propertyGroup("def").func("Cp").set("arg", "T");
    model.component("comp1").material("mat1").propertyGroup("def").func("Cp")
         .set("pieces", new String[][]{{"200.0", "1600.0", "1047.63657-0.372589265*T^1+9.45304214E-4*T^2-6.02409443E-7*T^3+1.2858961E-10*T^4"}});
    model.component("comp1").material("mat1").propertyGroup("def").func("Cp").set("argunit", "K");
    model.component("comp1").material("mat1").propertyGroup("def").func("Cp").set("fununit", "J/(kg*K)");
    model.component("comp1").material("mat1").propertyGroup("def").func("rho")
         .set("expr", "pA*0.02897/R_const[K*mol/J]/T");
    model.component("comp1").material("mat1").propertyGroup("def").func("rho").set("args", new String[]{"pA", "T"});
    model.component("comp1").material("mat1").propertyGroup("def").func("rho").set("dermethod", "manual");
    model.component("comp1").material("mat1").propertyGroup("def").func("rho")
         .set("argders", new String[][]{{"pA", "d(pA*0.02897/R_const/T,pA)"}, {"T", "d(pA*0.02897/R_const/T,T)"}});
    model.component("comp1").material("mat1").propertyGroup("def").func("rho").set("fununit", "kg/m^3");
    model.component("comp1").material("mat1").propertyGroup("def").func("rho")
         .set("argunit", new String[]{"Pa", "K"});
    model.component("comp1").material("mat1").propertyGroup("def").func("rho")
         .set("plotargs", new String[][]{{"pA", "101325", "101325"}, {"T", "273.15", "293.15"}});
    model.component("comp1").material("mat1").propertyGroup("def").func("k").set("arg", "T");
    model.component("comp1").material("mat1").propertyGroup("def").func("k")
         .set("pieces", new String[][]{{"200.0", "1600.0", "-0.00227583562+1.15480022E-4*T^1-7.90252856E-8*T^2+4.11702505E-11*T^3-7.43864331E-15*T^4"}});
    model.component("comp1").material("mat1").propertyGroup("def").func("k").set("argunit", "K");
    model.component("comp1").material("mat1").propertyGroup("def").func("k").set("fununit", "W/(m*K)");
    model.component("comp1").material("mat1").propertyGroup("def").func("cs")
         .set("expr", "sqrt(1.4*R_const[K*mol/J]/0.02897*T)");
    model.component("comp1").material("mat1").propertyGroup("def").func("cs").set("args", new String[]{"T"});
    model.component("comp1").material("mat1").propertyGroup("def").func("cs").set("dermethod", "manual");
    model.component("comp1").material("mat1").propertyGroup("def").func("cs").set("fununit", "m/s");
    model.component("comp1").material("mat1").propertyGroup("def").func("cs").set("argunit", new String[]{"K"});
    model.component("comp1").material("mat1").propertyGroup("def").func("cs")
         .set("plotargs", new String[][]{{"T", "273.15", "373.15"}});
    model.component("comp1").material("mat1").propertyGroup("def").func("an1").label("Analytic ");
    model.component("comp1").material("mat1").propertyGroup("def").func("an1").set("funcname", "alpha_p");
    model.component("comp1").material("mat1").propertyGroup("def").func("an1")
         .set("expr", "-1/rho(pA,T)*d(rho(pA,T),T)");
    model.component("comp1").material("mat1").propertyGroup("def").func("an1").set("args", new String[]{"pA", "T"});
    model.component("comp1").material("mat1").propertyGroup("def").func("an1").set("fununit", "1/K");
    model.component("comp1").material("mat1").propertyGroup("def").func("an1")
         .set("argunit", new String[]{"Pa", "K"});
    model.component("comp1").material("mat1").propertyGroup("def").func("an1")
         .set("plotargs", new String[][]{{"pA", "101325", "101325"}, {"T", "273.15", "373.15"}});
    model.component("comp1").material("mat1").propertyGroup("def").func("an2").set("funcname", "muB");
    model.component("comp1").material("mat1").propertyGroup("def").func("an2").set("expr", "0.6*eta(T)");
    model.component("comp1").material("mat1").propertyGroup("def").func("an2").set("args", new String[]{"T"});
    model.component("comp1").material("mat1").propertyGroup("def").func("an2").set("fununit", "Pa*s");
    model.component("comp1").material("mat1").propertyGroup("def").func("an2").set("argunit", new String[]{"K"});
    model.component("comp1").material("mat1").propertyGroup("def").func("an2")
         .set("plotargs", new String[][]{{"T", "200", "1600"}});
    model.component("comp1").material("mat1").propertyGroup("def").set("thermalexpansioncoefficient", "");
    model.component("comp1").material("mat1").propertyGroup("def").set("molarmass", "");
    model.component("comp1").material("mat1").propertyGroup("def").set("bulkviscosity", "");
    model.component("comp1").material("mat1").propertyGroup("def")
         .set("thermalexpansioncoefficient", new String[]{"alpha_p(pA,T)", "0", "0", "0", "alpha_p(pA,T)", "0", "0", "0", "alpha_p(pA,T)"});
    model.component("comp1").material("mat1").propertyGroup("def").set("molarmass", "0.02897[kg/mol]");
    model.component("comp1").material("mat1").propertyGroup("def").set("bulkviscosity", "muB(T)");
    model.component("comp1").material("mat1").propertyGroup("def").set("dynamicviscosity", "eta(T)");
    model.component("comp1").material("mat1").propertyGroup("def").set("ratioofspecificheat", "1.4");
    model.component("comp1").material("mat1").propertyGroup("def")
         .set("electricconductivity", new String[]{"0[S/m]", "0", "0", "0", "0[S/m]", "0", "0", "0", "0[S/m]"});
    model.component("comp1").material("mat1").propertyGroup("def").set("heatcapacity", "Cp(T)");
    model.component("comp1").material("mat1").propertyGroup("def").set("density", "rho(pA,T)");
    model.component("comp1").material("mat1").propertyGroup("def")
         .set("thermalconductivity", new String[]{"k(T)", "0", "0", "0", "k(T)", "0", "0", "0", "k(T)"});
    model.component("comp1").material("mat1").propertyGroup("def").set("soundspeed", "cs(T)");
    model.component("comp1").material("mat1").propertyGroup("def").addInput("temperature");
    model.component("comp1").material("mat1").propertyGroup("def").addInput("pressure");
    model.component("comp1").material("mat1").propertyGroup("NonlinearModel").set("BA", "(def.gamma+1)/2");
    model.component("comp1").material("mat1").propertyGroup("idealGas").func("Cp").label("Piecewise 2");
    model.component("comp1").material("mat1").propertyGroup("idealGas").func("Cp").set("arg", "T");
    model.component("comp1").material("mat1").propertyGroup("idealGas").func("Cp")
         .set("pieces", new String[][]{{"200.0", "1600.0", "1047.63657-0.372589265*T^1+9.45304214E-4*T^2-6.02409443E-7*T^3+1.2858961E-10*T^4"}});
    model.component("comp1").material("mat1").propertyGroup("idealGas").func("Cp").set("argunit", "K");
    model.component("comp1").material("mat1").propertyGroup("idealGas").func("Cp").set("fununit", "J/(kg*K)");
    model.component("comp1").material("mat1").propertyGroup("idealGas").set("Rs", "R_const/Mn");
    model.component("comp1").material("mat1").propertyGroup("idealGas").set("heatcapacity", "Cp(T)");
    model.component("comp1").material("mat1").propertyGroup("idealGas").set("ratioofspecificheat", "1.4");
    model.component("comp1").material("mat1").propertyGroup("idealGas").set("molarmass", "0.02897");
    model.component("comp1").material("mat1").propertyGroup("idealGas").addInput("temperature");
    model.component("comp1").material("mat1").propertyGroup("idealGas").addInput("pressure");
    model.component("comp1").material("mat1").materialType("nonSolid");

    model.component("comp1").coordSystem("sys1").set("name", "sys2");

    model.component("comp1").physics("spf").prop("ShapeProperty").set("order_fluid", 3);
    model.component("comp1").physics("spf").prop("PhysicalModelProperty").set("IncludeGravity", true);
    model.component("comp1").physics("spf").prop("PhysicalModelProperty").set("UseReducedPressure", true);
    model.component("comp1").physics("spf").prop("PhysicalModelProperty").set("Tref", "T0");
    model.component("comp1").physics("spf").prop("TurbulenceModelProperty").set("TurbulenceModel", "AlgebraicYplus");
    model.component("comp1").physics("spf").prop("InconsistentStabilization").set("IsotropicDiffusion", true);
    model.component("comp1").physics("spf").prop("AdvancedSettingProperty").set("useBNS", true);
    model.component("comp1").physics("spf").feature("fp1").set("rho_mat", "userdef");
    model.component("comp1").physics("spf").feature("fp1").set("rho", "rho");
    model.component("comp1").physics("spf").feature("fp1").set("minput_temperature_src", "userdef");
    model.component("comp1").physics("spf").feature("fp1").set("minput_temperature", "T");
    model.component("comp1").physics("spf").feature("fp1").set("mu_mat", "userdef");
    model.component("comp1").physics("spf").feature("fp1").set("mu", "mu");
    model.component("comp1").physics("spf").feature("init1")
         .set("u_init", new String[][]{{"1e-12"}, {"1e-12"}, {"0"}});
    model.component("comp1").physics("spf").feature("init1").set("p_init", "1[atm]");
    model.component("comp1").physics("spf").feature("init1").set("CompensateForHydrostaticPressure", false);
    model.component("comp1").physics("spf").feature("out1").set("p0", "1 [atm]");
    model.component("comp1").physics("spf").feature("vf1")
         .set("F", new String[][]{{"0"}, {"rho * g_const * alpha * (T - 293.15 [K])"}, {"0"}});
    model.component("comp1").physics("spf").feature("vf2").set("F", new String[][]{{"Fstx"}, {"Fsty"}, {"0"}});
    model.component("comp1").physics("pf").prop("ShapeProperty").set("order_phasefield", 3);
    model.component("comp1").physics("pf").feature("pfm1").set("epsilon_pf", "epsilon");
    model.component("comp1").physics("pf").feature("pfm1").set("chiOption", "velocity");
    model.component("comp1").physics("pf").feature("pfm1").set("chi", "0.001 * epsilon");
    model.component("comp1").physics("pf").feature("pfm1").set("U", "U");
    model.component("comp1").physics("pf").feature("pfm1").set("u", new String[][]{{"u"}, {"v"}, {"0"}});
    model.component("comp1").physics("pf").feature("pfm1").set("sigma", "sigma");
    model.component("comp1").physics("pf").feature("ww1").set("surfaceE1", "0.001[J/m^2]");
    model.component("comp1").physics("pf").feature("ww1").set("surfaceE2", "0.001[J/m^2]");
    model.component("comp1").physics("ht").prop("ShapeProperty").set("order_temperature", 3);
    model.component("comp1").physics("ht").prop("PhysicalModelProperty").set("Tref", 300);
    model.component("comp1").physics("ht").prop("InconsistentStabilization").set("HeatIsotropicDiffusion", true);
    model.component("comp1").physics("ht").feature("fluid1").set("Cp_mat", "userdef");
    model.component("comp1").physics("ht").feature("fluid1").set("Cp", "Cp");
    model.component("comp1").physics("ht").feature("fluid1").set("rho_mat", "userdef");
    model.component("comp1").physics("ht").feature("fluid1").set("rho", "rho");
    model.component("comp1").physics("ht").feature("fluid1").set("k_mat", "userdef");
    model.component("comp1").physics("ht").feature("fluid1")
         .set("k", new String[][]{{"k"}, {"0"}, {"0"}, {"0"}, {"k"}, {"0"}, {"0"}, {"0"}, {"k"}});
    model.component("comp1").physics("ht").feature("fluid1").set("u", new String[][]{{"u"}, {"v"}, {"0"}});
    model.component("comp1").physics("ht").feature("fluid1").set("minput_pressure_src", "root.comp1.spf.pA");
    model.component("comp1").physics("ht").feature("fluid1").set("minput_pressure", "p");
    model.component("comp1").physics("ht").feature("init1").set("Tinit", "T0");
    model.component("comp1").physics("ht").feature("hs1").set("Q0", "Qlaser");
    model.component("comp1").physics("ht").feature("hs2").set("Q0", "Qsb");

    model.component("comp1").mesh("mesh3").label("Mesh 2");
    model.component("comp1").mesh("mesh3").feature("size").set("hauto", 6);
    model.component("comp1").mesh("mesh3").feature("size").set("custom", "on");
    model.component("comp1").mesh("mesh3").feature("size").set("table", "cfd");
    model.component("comp1").mesh("mesh3").feature("size").set("hmax", 10);
    model.component("comp1").mesh("mesh3").feature("size").set("hmin", 8);
    model.component("comp1").mesh("mesh3").feature("size").set("hgrad", 2);
    model.component("comp1").mesh("mesh3").feature("ftri1").set("smoothmaxiter", 10);
    model.component("comp1").mesh("mesh3").feature("ftri1").set("smoothmaxdepth", 8);
    model.component("comp1").mesh("mesh3").feature("ftri1").set("method", "af");
    model.component("comp1").mesh("mesh3").feature("ftri1").feature("size1").set("hauto", 3);
    model.component("comp1").mesh("mesh3").feature("ftri1").feature("size1").set("custom", "on");
    model.component("comp1").mesh("mesh3").feature("ftri1").feature("size1").set("hmax", 16);
    model.component("comp1").mesh("mesh3").feature("ftri1").feature("size1").set("hmaxactive", true);
    model.component("comp1").mesh("mesh3").feature("ftri1").feature("size1").set("hmin", 15);
    model.component("comp1").mesh("mesh3").feature("ftri1").feature("size1").set("hminactive", true);
    model.component("comp1").mesh("mesh3").feature("ftri1").feature("size1").set("hcurve", 1);
    model.component("comp1").mesh("mesh3").feature("ftri1").feature("size1").set("hcurveactive", true);
    model.component("comp1").mesh("mesh3").feature("ftri1").feature("size1").set("hnarrow", 0.1);
    model.component("comp1").mesh("mesh3").feature("ftri1").feature("size1").set("hnarrowactive", true);
    model.component("comp1").mesh("mesh3").feature("ftri1").feature("size1").set("hgrad", 5);
    model.component("comp1").mesh("mesh3").feature("ftri1").feature("size1").set("hgradactive", true);
    model.component("comp1").mesh("mesh3").run();

    model.component("comp1").probe("point1").set("expr", "ht.T");
    model.component("comp1").probe("point1").set("unit", "K");
    model.component("comp1").probe("point1").set("descr", "Temperature");
    model.component("comp1").probe("point1").set("table", "tbl1");
    model.component("comp1").probe("point1").set("window", "window2");
    model.component("comp1").probe("point2").set("expr", "ht.T");
    model.component("comp1").probe("point2").set("unit", "K");
    model.component("comp1").probe("point2").set("descr", "Temperature");
    model.component("comp1").probe("point2").set("table", "tbl1");
    model.component("comp1").probe("point2").set("window", "window2");
    model.component("comp1").probe("point3").set("expr", "comp1.mu2");
    model.component("comp1").probe("point3").set("unit", "Pa*s");
    model.component("comp1").probe("point3").set("descr", "Dependent variable mu2");
    model.component("comp1").probe("point3").set("table", "tbl1");
    model.component("comp1").probe("point3").set("window", "window2");

    model.component("comp1").physics("ht").feature("fluid1").set("minput_pressure_src", "root.comp1.spf.pA");

    model.study().create("std1");
    model.study("std1").create("param", "Parametric");
    model.study("std1").create("phasei", "PhaseInitialization");
    model.study("std1").create("time", "Transient");
    model.study("std1").feature("phasei")
         .set("activate", new String[]{"spf", "off", "pf", "on", "ht", "on", "frame:spatial1", "on", "frame:material1", "on"});

    model.sol().create("sol1");
    model.sol("sol1").study("std1");
    model.sol("sol1").create("st1", "StudyStep");
    model.sol("sol1").create("v1", "Variables");
    model.sol("sol1").create("s1", "Stationary");
    model.sol("sol1").create("su1", "StoreSolution");
    model.sol("sol1").create("st2", "StudyStep");
    model.sol("sol1").create("v2", "Variables");
    model.sol("sol1").create("t1", "Time");
    model.sol("sol1").feature("s1").create("fc1", "FullyCoupled");
    model.sol("sol1").feature("s1").create("d1", "Direct");
    model.sol("sol1").feature("s1").create("i1", "Iterative");
    model.sol("sol1").feature("s1").feature("i1").create("mg1", "Multigrid");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr").create("sl1", "SORLine");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po").create("sl1", "SORLine");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("cs").create("d1", "Direct");
    model.sol("sol1").feature("s1").feature().remove("fcDef");
    model.sol("sol1").feature("t1").create("fc1", "FullyCoupled");
    model.sol("sol1").feature("t1").create("d1", "Direct");
    model.sol("sol1").feature("t1").create("i1", "Iterative");
    model.sol("sol1").feature("t1").create("i2", "Iterative");
    model.sol("sol1").feature("t1").create("i3", "Iterative");
    model.sol("sol1").feature("t1").feature("i1").create("mg1", "Multigrid");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("pr").create("sc1", "SCGS");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("po").create("sc1", "SCGS");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("cs").create("d1", "Direct");
    model.sol("sol1").feature("t1").feature("i2").create("bns1", "BlockNavierStokes");
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("vs").create("mg1", "Multigrid");
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1").feature("pr")
         .create("sl1", "SORLine");
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1").feature("po")
         .create("sl1", "SORLine");
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1").feature("cs")
         .create("d1", "Direct");
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("ps").create("mg1", "Multigrid");
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1").feature("pr")
         .create("sl1", "SORLine");
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1").feature("po")
         .create("sl1", "SORLine");
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1").feature("cs")
         .create("d1", "Direct");
    model.sol("sol1").feature("t1").feature("i3").create("mg1", "Multigrid");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("pr").create("so1", "SOR");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("po").create("so1", "SOR");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("cs").create("d1", "Direct");
    model.sol("sol1").feature("t1").feature().remove("fcDef");
    model.sol().create("sol3");
    model.sol("sol3").study("std1");
    model.sol().create("sol4");

    model.result().dataset().remove("dset4");
    model.result().dataset("dset3").set("solution", "sol4");

    model.sol("sol4").study("std1");
    model.sol("sol4").attach("std1");
    model.sol("sol4").create("st1", "StudyStep");
    model.sol("sol4").create("v1", "Variables");
    model.sol("sol4").create("s1", "Stationary");
    model.sol("sol4").create("su1", "StoreSolution");
    model.sol("sol4").create("st2", "StudyStep");
    model.sol("sol4").create("v2", "Variables");
    model.sol("sol4").create("t1", "Time");
    model.sol("sol4").feature("s1").create("fc1", "FullyCoupled");
    model.sol("sol4").feature("s1").create("d1", "Direct");
    model.sol("sol4").feature("s1").create("i1", "Iterative");
    model.sol("sol4").feature("s1").feature("i1").create("mg1", "Multigrid");
    model.sol("sol4").feature("s1").feature("i1").feature("mg1").feature("pr").create("sl1", "SORLine");
    model.sol("sol4").feature("s1").feature("i1").feature("mg1").feature("po").create("sl1", "SORLine");
    model.sol("sol4").feature("s1").feature("i1").feature("mg1").feature("cs").create("d1", "Direct");

    return model;
  }

  public static Model run2(Model model) {
    model.sol("sol4").feature("s1").feature().remove("fcDef");
    model.sol("sol4").feature("t1").create("fc1", "FullyCoupled");
    model.sol("sol4").feature("t1").create("d1", "Direct");
    model.sol("sol4").feature("t1").create("i1", "Iterative");
    model.sol("sol4").feature("t1").create("i2", "Iterative");
    model.sol("sol4").feature("t1").create("i3", "Iterative");
    model.sol("sol4").feature("t1").feature("i1").create("mg1", "Multigrid");
    model.sol("sol4").feature("t1").feature("i1").feature("mg1").feature("pr").create("sc1", "SCGS");
    model.sol("sol4").feature("t1").feature("i1").feature("mg1").feature("po").create("sc1", "SCGS");
    model.sol("sol4").feature("t1").feature("i1").feature("mg1").feature("cs").create("d1", "Direct");
    model.sol("sol4").feature("t1").feature("i2").create("bns1", "BlockNavierStokes");
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("vs").create("mg1", "Multigrid");
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1").feature("pr")
         .create("sl1", "SORLine");
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1").feature("po")
         .create("sl1", "SORLine");
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1").feature("cs")
         .create("d1", "Direct");
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("ps").create("mg1", "Multigrid");
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1").feature("pr")
         .create("sl1", "SORLine");
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1").feature("po")
         .create("sl1", "SORLine");
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1").feature("cs")
         .create("d1", "Direct");
    model.sol("sol4").feature("t1").feature("i3").create("mg1", "Multigrid");
    model.sol("sol4").feature("t1").feature("i3").feature("mg1").feature("pr").create("so1", "SOR");
    model.sol("sol4").feature("t1").feature("i3").feature("mg1").feature("po").create("so1", "SOR");
    model.sol("sol4").feature("t1").feature("i3").feature("mg1").feature("cs").create("d1", "Direct");
    model.sol("sol4").feature("t1").feature().remove("fcDef");

    model.result().dataset().create("an1_ds1", "Grid1D");
    model.result().dataset().create("avh1", "Average");
    model.result().dataset().create("avh2", "Average");
    model.result().dataset().create("avh3", "Average");
    model.result().dataset().create("dset5", "Solution");
    model.result().dataset().create("dset6", "Solution");
    model.result().dataset("an1_ds1").set("data", "none");
    model.result().dataset("dset3").set("probetag", "point3");
    model.result().dataset("avh1").set("probetag", "point1");
    model.result().dataset("avh1").set("data", "dset3");
    model.result().dataset("avh1").selection().geom("geom2", 0);
    model.result().dataset("avh1").selection().set(27);
    model.result().dataset("avh2").set("probetag", "point2");
    model.result().dataset("avh2").set("data", "dset3");
    model.result().dataset("avh2").selection().geom("geom2", 0);
    model.result().dataset("avh2").selection().set(26);
    model.result().dataset("avh3").set("probetag", "point3");
    model.result().dataset("avh3").set("data", "dset3");
    model.result().dataset("avh3").selection().geom("geom2", 0);
    model.result().dataset("avh3").selection().set(34);
    model.result().dataset("dset4").set("solution", "sol3");
    model.result().dataset("dset5").set("solution", "sol4");
    model.result().dataset("dset6").set("solution", "sol5");
    model.result().numerical().create("pev1", "EvalPoint");
    model.result().numerical().create("pev2", "EvalPoint");
    model.result().numerical().create("pev3", "EvalPoint");
    model.result().numerical("pev1").set("probetag", "point1");
    model.result().numerical("pev2").set("probetag", "point2");
    model.result().numerical("pev3").set("probetag", "point3");
    model.result().create("pg1", "PlotGroup2D");
    model.result().create("pg10", "PlotGroup2D");
    model.result().create("pg2", "PlotGroup2D");
    model.result().create("pg3", "PlotGroup2D");
    model.result().create("pg4", "PlotGroup2D");
    model.result().create("pg6", "PlotGroup1D");
    model.result().create("pg7", "PlotGroup1D");
    model.result().create("pg8", "PlotGroup2D");
    model.result().create("pg9", "PlotGroup2D");
    model.result("pg1").set("data", "dset5");
    model.result("pg1").create("surf1", "Surface");
    model.result("pg1").create("con1", "Contour");
    model.result("pg1").feature("surf1").set("expr", "Qlaser");
    model.result("pg1").feature("con1").set("expr", "pf.Vf2");
    model.result("pg10").set("data", "dset5");
    model.result("pg10").create("surf1", "Surface");
    model.result("pg10").create("con1", "Contour");
    model.result("pg10").feature("surf1").set("expr", "T");
    model.result("pg10").feature("con1").set("expr", "pf.Vf2");
    model.result("pg2").set("data", "dset5");
    model.result("pg2").create("surf1", "Surface");
    model.result("pg2").create("con1", "Contour");
    model.result("pg2").feature("surf1").set("expr", "pf.Vf2");
    model.result("pg2").feature("con1").set("expr", "p");
    model.result("pg3").set("data", "dset5");
    model.result("pg3").create("surf1", "Surface");
    model.result("pg3").create("con1", "Contour");
    model.result("pg3").create("arwl1", "ArrowLine");
    model.result("pg3").feature("con1").set("expr", "pf.Vf2");
    model.result("pg4").set("data", "dset5");
    model.result("pg4").create("surf1", "Surface");
    model.result("pg4").create("con1", "Contour");
    model.result("pg4").feature("surf1").set("expr", "mu");
    model.result("pg4").feature("con1").set("expr", "pf.Vf2");
    model.result("pg6").create("plot1", "Function");
    model.result("pg6").feature("plot1").set("expr", "comp1.cpmat(T)");
    model.result("pg7").set("probetag", "window2_default");
    model.result("pg7").create("tblp1", "Table");
    model.result("pg7").feature("tblp1").set("probetag", "point1,point2,point3");
    model.result("pg8").create("surf1", "Surface");
    model.result("pg8").feature("surf1").set("expr", "rho");
    model.result("pg9").create("surf1", "Surface");
    model.result("pg9").feature("surf1").set("expr", "Qlaser");
    model.result().export().create("anim1", "Animation");
    model.result().export().create("anim2", "Animation");
    model.result().export().create("anim3", "Animation");
    model.result().export().create("anim4", "Animation");

    model.component("comp1").probe("point1").genResult(null);
    model.component("comp1").probe("point2").genResult(null);
    model.component("comp1").probe("point3").genResult(null);

    model.study("std1").feature("param").set("sweeptype", "filled");
    model.study("std1").feature("param").set("pname", new String[]{"vlaser", "Ep"});
    model.study("std1").feature("param")
         .set("plistarr", new String[]{"2 [mm/s], 3 [mm/s], 5 [mm/s], 10 [mm/s]", "0.1 [W], 1 [W], 10 [W]"});
    model.study("std1").feature("param").set("punit", new String[]{"m/s", "W"});
    model.study("std1").feature("time").set("tlist", "range(0,0.001,1)");
    model.study("std1").feature("time").set("useinitsol", true);
    model.study("std1").feature("time").set("initstudy", "std1");
    model.study("std1").feature("time").set("solnum", "auto");

    model.sol("sol1").feature("st1").label("Compile Equations: Phase Initialization");
    model.sol("sol1").feature("v1").label("Dependent Variables 1.1");
    model.sol("sol1").feature("s1").label("Stationary Solver 1.1");
    model.sol("sol1").feature("s1").feature("dDef").label("Direct 2");
    model.sol("sol1").feature("s1").feature("aDef").label("Advanced 1");
    model.sol("sol1").feature("s1").feature("aDef").set("cachepattern", true);
    model.sol("sol1").feature("s1").feature("fc1").label("Fully Coupled 1.1");
    model.sol("sol1").feature("s1").feature("fc1").set("linsolver", "d1");
    model.sol("sol1").feature("s1").feature("fc1").set("initstep", 0.01);
    model.sol("sol1").feature("s1").feature("fc1").set("minstep", 1.0E-6);
    model.sol("sol1").feature("s1").feature("fc1").set("maxiter", 50);
    model.sol("sol1").feature("s1").feature("d1").label("Direct, interface distance (pf)");
    model.sol("sol1").feature("s1").feature("d1").set("linsolver", "pardiso");
    model.sol("sol1").feature("s1").feature("d1").set("pivotperturb", 1.0E-13);
    model.sol("sol1").feature("s1").feature("i1").label("AMG, interface distance (pf)");
    model.sol("sol1").feature("s1").feature("i1").set("nlinnormuse", true);
    model.sol("sol1").feature("s1").feature("i1").set("maxlinit", 1000);
    model.sol("sol1").feature("s1").feature("i1").feature("ilDef").label("Incomplete LU 1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").label("Multigrid 1.1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").set("prefun", "saamg");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").set("maxcoarsedof", 50000);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").set("saamgcompwise", true);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").set("usesmooth", false);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr").label("Presmoother 1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr").feature("soDef").label("SOR 1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr").feature("sl1").label("SOR Line 1.1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr").feature("sl1")
         .set("linesweeptype", "ssor");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr").feature("sl1").set("iter", 1);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr").feature("sl1").set("linerelax", 0.7);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr").feature("sl1")
         .set("linemethod", "uncoupled");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr").feature("sl1").set("relax", 0.5);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po").label("Postsmoother 1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po").feature("soDef").label("SOR 1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po").feature("sl1").label("SOR Line 1.1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po").feature("sl1")
         .set("linesweeptype", "ssor");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po").feature("sl1").set("iter", 1);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po").feature("sl1").set("linerelax", 0.7);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po").feature("sl1")
         .set("linemethod", "uncoupled");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po").feature("sl1").set("relax", 0.5);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("cs").label("Coarse Solver 1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("cs").feature("dDef").label("Direct 2");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("cs").feature("d1").label("Direct 1.1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("cs").feature("d1")
         .set("linsolver", "pardiso");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("cs").feature("d1")
         .set("pivotperturb", 1.0E-13);
    model.sol("sol1").feature("su1").label("Solution Store 1.1");
    model.sol("sol1").feature("st2").label("Compile Equations: Time Dependent");
    model.sol("sol1").feature("st2").set("studystep", "time");
    model.sol("sol1").feature("v2").label("Dependent Variables 2.1");
    model.sol("sol1").feature("v2").set("initsol", "sol1");
    model.sol("sol1").feature("v2").set("solnum", "auto");
    model.sol("sol1").feature("v2").set("resscalemethod", "manual");
    model.sol("sol1").feature("v2").set("notsolmethod", "sol");
    model.sol("sol1").feature("v2").set("notsol", "sol1");
    model.sol("sol1").feature("v2").set("notsolnum", "auto");
    model.sol("sol1").feature("v2").set("clist", new String[]{"range(0,0.001,1)", "0.001[s]"});
    model.sol("sol1").feature("v2").feature("comp1_phipf").set("scalemethod", "manual");
    model.sol("sol1").feature("v2").feature("comp1_phipf").set("scaleval", 1);
    model.sol("sol1").feature("t1").label("Time-Dependent Solver 1.1");
    model.sol("sol1").feature("t1").set("control", "time");
    model.sol("sol1").feature("t1").set("tlist", "range(0,0.001,1)");
    model.sol("sol1").feature("t1").set("rtol", 0.005);
    model.sol("sol1").feature("t1").set("atolglobalfactor", 0.05);
    model.sol("sol1").feature("t1")
         .set("atolmethod", new String[]{"comp1_GI", "global", "comp1_p", "scaled", "comp1_phipf", "global", "comp1_psi", "global", "comp1_T", "global", 
         "comp1_u", "global"});
    model.sol("sol1").feature("t1")
         .set("atolfactor", new String[]{"comp1_GI", "0.1", "comp1_p", "1", "comp1_phipf", "0.1", "comp1_psi", "0.1", "comp1_T", "0.1", 
         "comp1_u", "0.1"});
    model.sol("sol1").feature("t1").set("maxorder", 2);
    model.sol("sol1").feature("t1").set("stabcntrl", true);
    model.sol("sol1").feature("t1").set("bwinitstepfrac", 0.01);
    model.sol("sol1").feature("t1").set("estrat", "exclude");
    model.sol("sol1").feature("t1").feature("dDef").label("Direct 2");
    model.sol("sol1").feature("t1").feature("aDef").label("Advanced 1");
    model.sol("sol1").feature("t1").feature("aDef").set("cachepattern", true);
    model.sol("sol1").feature("t1").feature("fc1").label("Fully Coupled 1.1");
    model.sol("sol1").feature("t1").feature("fc1").set("linsolver", "d1");
    model.sol("sol1").feature("t1").feature("fc1").set("maxiter", 8);
    model.sol("sol1").feature("t1").feature("fc1").set("ntolfact", 0.5);
    model.sol("sol1").feature("t1").feature("fc1").set("damp", "0.9");
    model.sol("sol1").feature("t1").feature("fc1").set("jtech", "once");
    model.sol("sol1").feature("t1").feature("fc1").set("stabacc", "aacc");
    model.sol("sol1").feature("t1").feature("fc1").set("aaccdim", 5);
    model.sol("sol1").feature("t1").feature("fc1").set("aaccmix", 0.9);
    model.sol("sol1").feature("t1").feature("fc1").set("aaccdelay", 1);
    model.sol("sol1").feature("t1").feature("d1").label("Direct, heat transfer variables (ht) (merged)");
    model.sol("sol1").feature("t1").feature("d1").set("linsolver", "pardiso");
    model.sol("sol1").feature("t1").feature("d1").set("pivotperturb", 1.0E-13);
    model.sol("sol1").feature("t1").feature("i1").label("AMG, phase field variables (pf)");
    model.sol("sol1").feature("t1").feature("i1").set("maxlinit", 50);
    model.sol("sol1").feature("t1").feature("i1").set("rhob", 20);
    model.sol("sol1").feature("t1").feature("i1").feature("ilDef").label("Incomplete LU 1");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").label("Multigrid 1.1");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").set("prefun", "saamg");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").set("maxcoarsedof", 50000);
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").set("usesmooth", false);
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("pr").label("Presmoother 1");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("pr").feature("soDef").label("SOR 1");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("pr").feature("sc1").label("SCGS 1.1");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("pr").feature("sc1")
         .set("linesweeptype", "ssor");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("pr").feature("sc1").set("iter", 0);
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("pr").feature("sc1").set("approxscgs", true);
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("pr").feature("sc1")
         .set("scgsdirectmaxsize", 1000);
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("po").label("Postsmoother 1");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("po").feature("soDef").label("SOR 1");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("po").feature("sc1").label("SCGS 1.1");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("po").feature("sc1")
         .set("linesweeptype", "ssor");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("po").feature("sc1").set("iter", 1);
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("po").feature("sc1").set("approxscgs", true);
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("po").feature("sc1")
         .set("scgsdirectmaxsize", 1000);
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("cs").label("Coarse Solver 1");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("cs").feature("dDef").label("Direct 2");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("cs").feature("d1").label("Direct 1.1");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("cs").feature("d1")
         .set("linsolver", "pardiso");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("cs").feature("d1")
         .set("pivotperturb", 1.0E-13);
    model.sol("sol1").feature("t1").feature("i2").label("Block Navier-Stokes, fluid flow variables (spf)");
    model.sol("sol1").feature("t1").feature("i2").set("maxlinit", 100);
    model.sol("sol1").feature("t1").feature("i2").set("rhob", 20);
    model.sol("sol1").feature("t1").feature("i2").feature("ilDef").label("Incomplete LU 1");
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").label("Block Navier-Stokes 1.1");
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").set("schurcomplementapproximation", "abssumf");
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").set("velocityvars", new String[]{"comp1_u"});
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").set("pressurevars", new String[]{"comp1_p"});
    model.sol("sol1").feature("t1").feature("i2").feature("bns1")
         .set("hybridvar", new String[]{"comp1_u", "comp1_p"});
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("vs").label("Velocity Solver 1");
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("vs").feature("dDef").label("Direct 1");
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1").label("Multigrid 1.1");
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1").set("prefun", "saamg");
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1").set("iter", 1);
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1")
         .set("maxcoarsedof", 50000);
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1").set("strconn", 0.02);
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1")
         .set("saamgcompwise", true);
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1")
         .set("usesmooth", false);
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1").feature("pr")
         .label("Presmoother 1");
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1").feature("pr")
         .feature("soDef").label("SOR 1");
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1").feature("pr")
         .feature("sl1").label("SOR Line 1.1");
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1").feature("pr")
         .feature("sl1").set("iter", 1);
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1").feature("po")
         .label("Postsmoother 1");
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1").feature("po")
         .feature("soDef").label("SOR 1");
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1").feature("po")
         .feature("sl1").label("SOR Line 1.1");
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1").feature("po")
         .feature("sl1").set("iter", 1);
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1").feature("cs")
         .label("Coarse Solver 1");
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1").feature("cs")
         .feature("dDef").label("Direct 2");
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1").feature("cs")
         .feature("d1").label("Direct 1.1");
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1").feature("cs")
         .feature("d1").set("linsolver", "pardiso");
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1").feature("cs")
         .feature("d1").set("pivotperturb", 1.0E-13);
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("ps").label("Pressure Solver 1");
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("ps").feature("dDef").label("Direct 1");
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1").label("Multigrid 1.1");
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1").set("prefun", "saamg");
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1")
         .set("maxcoarsedof", 50000);
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1").set("strconn", 0.02);
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1")
         .set("usesmooth", false);
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1").feature("pr")
         .label("Presmoother 1");
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1").feature("pr")
         .feature("soDef").label("SOR 1");
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1").feature("pr")
         .feature("sl1").label("SOR Line 1.1");
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1").feature("pr")
         .feature("sl1").set("iter", 1);
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1").feature("pr")
         .feature("sl1").set("linemethod", "uncoupled");
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1").feature("po")
         .label("Postsmoother 1");
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1").feature("po")
         .feature("soDef").label("SOR 1");
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1").feature("po")
         .feature("sl1").label("SOR Line 1.1");
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1").feature("po")
         .feature("sl1").set("iter", 1);
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1").feature("po")
         .feature("sl1").set("linemethod", "uncoupled");
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1").feature("cs")
         .label("Coarse Solver 1");
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1").feature("cs")
         .feature("dDef").label("Direct 2");
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1").feature("cs")
         .feature("d1").label("Direct 1.1");
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1").feature("cs")
         .feature("d1").set("linsolver", "pardiso");
    model.sol("sol1").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1").feature("cs")
         .feature("d1").set("pivotperturb", 1.0E-13);
    model.sol("sol1").feature("t1").feature("i3").label("AMG, heat transfer variables (ht)");
    model.sol("sol1").feature("t1").feature("i3").set("rhob", 20);
    model.sol("sol1").feature("t1").feature("i3").feature("ilDef").label("Incomplete LU 1");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").label("Multigrid 1.1");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").set("prefun", "saamg");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").set("maxcoarsedof", 50000);
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").set("saamgcompwise", true);
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").set("usesmooth", false);
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("pr").label("Presmoother 1");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("pr").feature("soDef").label("SOR 2");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("pr").feature("so1").label("SOR 1.1");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("pr").feature("so1").set("relax", 0.9);
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("po").label("Postsmoother 1");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("po").feature("soDef").label("SOR 2");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("po").feature("so1").label("SOR 1.1");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("po").feature("so1").set("relax", 0.9);
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("cs").label("Coarse Solver 1");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("cs").feature("dDef").label("Direct 2");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("cs").feature("d1").label("Direct 1.1");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("cs").feature("d1")
         .set("linsolver", "pardiso");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("cs").feature("d1")
         .set("pivotperturb", 1.0E-13);
    model.sol("sol1").runAll();
    model.sol("sol3").label("Refined Mesh Solution 1");
    model.sol("sol4").attach("std1");
    model.sol("sol4").feature("st1").label("Compile Equations: Phase Initialization");
    model.sol("sol4").feature("v1").label("Dependent Variables 1.1");
    model.sol("sol4").feature("s1").label("Stationary Solver 1.1");
    model.sol("sol4").feature("s1").feature("dDef").label("Direct 2");
    model.sol("sol4").feature("s1").feature("aDef").label("Advanced 1");
    model.sol("sol4").feature("s1").feature("aDef").set("cachepattern", true);
    model.sol("sol4").feature("s1").feature("fc1").label("Fully Coupled 1.1");
    model.sol("sol4").feature("s1").feature("fc1").set("linsolver", "d1");
    model.sol("sol4").feature("s1").feature("fc1").set("initstep", 0.01);
    model.sol("sol4").feature("s1").feature("fc1").set("minstep", 1.0E-6);
    model.sol("sol4").feature("s1").feature("fc1").set("maxiter", 50);
    model.sol("sol4").feature("s1").feature("d1").label("Direct, interface distance (pf)");
    model.sol("sol4").feature("s1").feature("d1").set("linsolver", "pardiso");
    model.sol("sol4").feature("s1").feature("d1").set("pivotperturb", 1.0E-13);
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
    model.sol("sol4").feature("su1").label("Solution Store 2.1");
    model.sol("sol4").feature("st2").label("Compile Equations: Time Dependent");
    model.sol("sol4").feature("st2").set("studystep", "time");
    model.sol("sol4").feature("v2").label("Dependent Variables 2.1");
    model.sol("sol4").feature("v2").set("initsol", "sol4");
    model.sol("sol4").feature("v2").set("solnum", "auto");
    model.sol("sol4").feature("v2").set("resscalemethod", "manual");
    model.sol("sol4").feature("v2").set("notsolmethod", "sol");
    model.sol("sol4").feature("v2").set("notsol", "sol4");
    model.sol("sol4").feature("v2").set("notsoluse", "sol5");
    model.sol("sol4").feature("v2").set("notsolnum", "auto");
    model.sol("sol4").feature("v2").set("clist", new String[]{"range(0,0.001,1)", "0.001[s]"});
    model.sol("sol4").feature("v2").feature("comp1_phipf").set("scalemethod", "manual");
    model.sol("sol4").feature("v2").feature("comp1_phipf").set("scaleval", 1);

    return model;
  }

  public static Model run3(Model model) {
    model.sol("sol4").feature("t1").label("Time-Dependent Solver 1.1");
    model.sol("sol4").feature("t1").set("control", "time");
    model.sol("sol4").feature("t1").set("tlist", "range(0,0.001,1)");
    model.sol("sol4").feature("t1").set("rtol", 0.005);
    model.sol("sol4").feature("t1").set("atolglobalfactor", 0.05);
    model.sol("sol4").feature("t1")
         .set("atolmethod", new String[]{"comp1_GI", "global", "comp1_p", "scaled", "comp1_phipf", "global", "comp1_psi", "global", "comp1_T", "global", 
         "comp1_u", "global"});
    model.sol("sol4").feature("t1")
         .set("atolfactor", new String[]{"comp1_GI", "0.1", "comp1_p", "1", "comp1_phipf", "0.1", "comp1_psi", "0.1", "comp1_T", "0.1", 
         "comp1_u", "0.1"});
    model.sol("sol4").feature("t1").set("maxorder", 2);
    model.sol("sol4").feature("t1").set("stabcntrl", true);
    model.sol("sol4").feature("t1").set("bwinitstepfrac", 0.01);
    model.sol("sol4").feature("t1").set("estrat", "exclude");
    model.sol("sol4").feature("t1").feature("dDef").label("Direct 2");
    model.sol("sol4").feature("t1").feature("aDef").label("Advanced 1");
    model.sol("sol4").feature("t1").feature("aDef").set("cachepattern", true);
    model.sol("sol4").feature("t1").feature("fc1").label("Fully Coupled 1.1");
    model.sol("sol4").feature("t1").feature("fc1").set("linsolver", "d1");
    model.sol("sol4").feature("t1").feature("fc1").set("maxiter", 8);
    model.sol("sol4").feature("t1").feature("fc1").set("ntolfact", 0.5);
    model.sol("sol4").feature("t1").feature("fc1").set("damp", "0.9");
    model.sol("sol4").feature("t1").feature("fc1").set("jtech", "once");
    model.sol("sol4").feature("t1").feature("fc1").set("stabacc", "aacc");
    model.sol("sol4").feature("t1").feature("fc1").set("aaccdim", 5);
    model.sol("sol4").feature("t1").feature("fc1").set("aaccmix", 0.9);
    model.sol("sol4").feature("t1").feature("fc1").set("aaccdelay", 1);
    model.sol("sol4").feature("t1").feature("d1").label("Direct, heat transfer variables (ht) (merged)");
    model.sol("sol4").feature("t1").feature("d1").set("linsolver", "pardiso");
    model.sol("sol4").feature("t1").feature("d1").set("pivotperturb", 1.0E-13);
    model.sol("sol4").feature("t1").feature("i1").label("AMG, phase field variables (pf)");
    model.sol("sol4").feature("t1").feature("i1").set("maxlinit", 50);
    model.sol("sol4").feature("t1").feature("i1").set("rhob", 20);
    model.sol("sol4").feature("t1").feature("i1").feature("ilDef").label("Incomplete LU 1");
    model.sol("sol4").feature("t1").feature("i1").feature("mg1").label("Multigrid 1.1");
    model.sol("sol4").feature("t1").feature("i1").feature("mg1").set("prefun", "saamg");
    model.sol("sol4").feature("t1").feature("i1").feature("mg1").set("maxcoarsedof", 50000);
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
    model.sol("sol4").feature("t1").feature("i2").label("Block Navier-Stokes, fluid flow variables (spf)");
    model.sol("sol4").feature("t1").feature("i2").set("maxlinit", 100);
    model.sol("sol4").feature("t1").feature("i2").set("rhob", 20);
    model.sol("sol4").feature("t1").feature("i2").feature("ilDef").label("Incomplete LU 1");
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").label("Block Navier-Stokes 1.1");
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").set("schurcomplementapproximation", "abssumf");
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").set("velocityvars", new String[]{"comp1_u"});
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").set("pressurevars", new String[]{"comp1_p"});
    model.sol("sol4").feature("t1").feature("i2").feature("bns1")
         .set("hybridvar", new String[]{"comp1_u", "comp1_p"});
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("vs").label("Velocity Solver 1");
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("vs").feature("dDef").label("Direct 1");
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1").label("Multigrid 1.1");
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1").set("prefun", "saamg");
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1").set("iter", 1);
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1")
         .set("maxcoarsedof", 50000);
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1").set("strconn", 0.02);
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1")
         .set("saamgcompwise", true);
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1")
         .set("usesmooth", false);
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1").feature("pr")
         .label("Presmoother 1");
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1").feature("pr")
         .feature("soDef").label("SOR 1");
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1").feature("pr")
         .feature("sl1").label("SOR Line 1.1");
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1").feature("pr")
         .feature("sl1").set("iter", 1);
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1").feature("po")
         .label("Postsmoother 1");
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1").feature("po")
         .feature("soDef").label("SOR 1");
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1").feature("po")
         .feature("sl1").label("SOR Line 1.1");
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1").feature("po")
         .feature("sl1").set("iter", 1);
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1").feature("cs")
         .label("Coarse Solver 1");
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1").feature("cs")
         .feature("dDef").label("Direct 2");
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1").feature("cs")
         .feature("d1").label("Direct 1.1");
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1").feature("cs")
         .feature("d1").set("linsolver", "pardiso");
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("vs").feature("mg1").feature("cs")
         .feature("d1").set("pivotperturb", 1.0E-13);
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("ps").label("Pressure Solver 1");
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("ps").feature("dDef").label("Direct 1");
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1").label("Multigrid 1.1");
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1").set("prefun", "saamg");
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1")
         .set("maxcoarsedof", 50000);
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1").set("strconn", 0.02);
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1")
         .set("usesmooth", false);
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1").feature("pr")
         .label("Presmoother 1");
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1").feature("pr")
         .feature("soDef").label("SOR 1");
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1").feature("pr")
         .feature("sl1").label("SOR Line 1.1");
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1").feature("pr")
         .feature("sl1").set("iter", 1);
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1").feature("pr")
         .feature("sl1").set("linemethod", "uncoupled");
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1").feature("po")
         .label("Postsmoother 1");
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1").feature("po")
         .feature("soDef").label("SOR 1");
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1").feature("po")
         .feature("sl1").label("SOR Line 1.1");
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1").feature("po")
         .feature("sl1").set("iter", 1);
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1").feature("po")
         .feature("sl1").set("linemethod", "uncoupled");
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1").feature("cs")
         .label("Coarse Solver 1");
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1").feature("cs")
         .feature("dDef").label("Direct 2");
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1").feature("cs")
         .feature("d1").label("Direct 1.1");
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1").feature("cs")
         .feature("d1").set("linsolver", "pardiso");
    model.sol("sol4").feature("t1").feature("i2").feature("bns1").feature("ps").feature("mg1").feature("cs")
         .feature("d1").set("pivotperturb", 1.0E-13);
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
    model.sol("sol4").runAll();

    model.result().dataset("an1_ds1").label("Grid 1D 1a");
    model.result().dataset("an1_ds1").set("function", "all");
    model.result().dataset("an1_ds1").set("par1", "T");
    model.result().dataset("an1_ds1").set("parmin1", 300);
    model.result().dataset("an1_ds1").set("parmax1", 400);
    model.result().dataset("an1_ds1").set("res1", 10000);
    model.result().dataset("an1_ds1").set("distribution", "mixed");
    model.result().dataset("dset3").label("Probe Solution 3");
    model.result().numerical("pev1").set("descr", new String[]{"", "Temperature", ""});
    model.result().numerical("pev3").set("unit", new String[]{""});
    model.result().numerical("pev1").setResult();
    model.result("pg1").label("Heating and Volume Fraction");
    model.result("pg1").set("looplevel", new int[]{201});
    model.result("pg1").set("titletype", "manual");
    model.result("pg1").set("title", "Heating (W/m^3) with Volume Fraction Contours");
    model.result("pg1").set("paramindicator", "t=eval(t,s) (s)");
    model.result("pg1").set("edges", false);
    model.result("pg1").set("showlegendsmaxmin", true);
    model.result("pg1").set("showlegendsunit", true);
    model.result("pg1").feature("surf1").set("coloring", "gradient");
    model.result("pg1").feature("surf1").set("topcolor", "red");
    model.result("pg1").feature("surf1").set("bottomcolor", "white");
    model.result("pg1").feature("surf1").set("resolution", "normal");
    model.result("pg1").feature("con1").set("coloring", "gradient");
    model.result("pg1").feature("con1").set("topcolor", "black");
    model.result("pg1").feature("con1").set("bottomcolor", "white");
    model.result("pg1").feature("con1").set("resolution", "normal");
    model.result("pg10").label("Temperature and Volume Fraction");
    model.result("pg10").set("looplevel", new int[]{48});
    model.result("pg10").set("titletype", "manual");
    model.result("pg10").set("title", "Temperature (K) with Volume Fraction Contours");
    model.result("pg10").set("paramindicator", "t=eval(t,s) (s)");
    model.result("pg10").set("edges", false);
    model.result("pg10").set("showlegendsmaxmin", true);
    model.result("pg10").set("showlegendsunit", true);
    model.result("pg10").feature("surf1").set("coloring", "gradient");
    model.result("pg10").feature("surf1").set("topcolor", "red");
    model.result("pg10").feature("surf1").set("bottomcolor", "white");
    model.result("pg10").feature("surf1").set("resolution", "normal");
    model.result("pg10").feature("con1").set("coloring", "gradient");
    model.result("pg10").feature("con1").set("topcolor", "black");
    model.result("pg10").feature("con1").set("bottomcolor", "white");
    model.result("pg10").feature("con1").set("resolution", "normal");
    model.result("pg2").label("Pressure and Volume Fraction");
    model.result("pg2").set("titletype", "manual");
    model.result("pg2").set("title", "Volume Fraction with Pressure (Pa) contours");
    model.result("pg2").set("paramindicator", "t=eval(t,s) (s)");
    model.result("pg2").feature("surf1").set("rangecoloractive", true);
    model.result("pg2").feature("surf1").set("rangecolormax", 1);
    model.result("pg2").feature("surf1").set("coloring", "gradient");
    model.result("pg2").feature("surf1").set("topcolor", "blue");
    model.result("pg2").feature("surf1").set("bottomcolor", "white");
    model.result("pg2").feature("surf1").set("resolution", "normal");
    model.result("pg2").feature("con1").set("coloring", "gradient");
    model.result("pg2").feature("con1").set("topcolor", "green");
    model.result("pg2").feature("con1").set("resolution", "normal");
    model.result("pg3").label("Velocity and Volume Fraction");
    model.result("pg3").set("looplevel", new int[]{199});
    model.result("pg3").set("titletype", "manual");
    model.result("pg3").set("title", "Velocity magnitude (m/s) with Volume Fraction contours");
    model.result("pg3").set("paramindicator", "t=eval(t) (s)");
    model.result("pg3").feature("surf1").set("colortable", "Magma");
    model.result("pg3").feature("surf1").set("resolution", "normal");
    model.result("pg3").feature("con1").set("coloring", "uniform");
    model.result("pg3").feature("con1").set("color", "green");
    model.result("pg3").feature("con1").set("resolution", "normal");
    model.result("pg3").feature("arwl1").active(false);
    model.result("pg3").feature("arwl1").set("scale", 4.224745667893173E7);
    model.result("pg3").feature("arwl1").set("scaleactive", false);
    model.result("pg4").label("Viscosity");
    model.result("pg4").set("looplevel", new int[]{195});
    model.result("pg4").set("titletype", "manual");
    model.result("pg4").set("title", "Fluid viscosity (Pa*s) with Volume Fraction contours");
    model.result("pg4").set("paramindicator", "t=eval(t,s) (s)");
    model.result("pg4").feature("surf1").set("descr", "mu");
    model.result("pg4").feature("surf1").set("coloring", "gradient");
    model.result("pg4").feature("surf1").set("topcolor", "magenta");
    model.result("pg4").feature("surf1").set("bottomcolor", "white");
    model.result("pg4").feature("surf1").set("resolution", "normal");
    model.result("pg4").feature("con1").set("resolution", "normal");
    model.result("pg6").label("Specific Heat of Material");
    model.result("pg6").set("data", "an1_ds1");
    model.result("pg6").set("solrepresentation", "solnum");
    model.result("pg6").set("titletype", "manual");
    model.result("pg6").set("title", "Specific Heat of Material (J/(kg*K))");
    model.result("pg6").set("xlabel", "T (K)");
    model.result("pg6").set("xlabelactive", true);
    model.result("pg6").set("ylabel", "J/(kg*K)");
    model.result("pg6").set("ylabelactive", true);
    model.result("pg6").feature("plot1").set("solrepresentation", "solnum");
    model.result("pg6").feature("plot1").set("unit", "");
    model.result("pg6").feature("plot1").set("descr", "cpmat(T)");
    model.result("pg6").feature("plot1").set("xdataexpr", "T");
    model.result("pg6").feature("plot1").set("xdataunit", "m");
    model.result("pg6").feature("plot1").set("xdatadescractive", true);
    model.result("pg6").feature("plot1").set("xdatadescr", "");
    model.result("pg6").feature("plot1").set("lowerbound", 300);
    model.result("pg6").feature("plot1").set("upperbound", 400);
    model.result("pg6").feature("plot1").set("linewidth", "preference");
    model.result("pg7").label("Temperature Probe");
    model.result("pg7").set("titletype", "manual");
    model.result("pg7").set("title", "Temperature Probe");
    model.result("pg7").set("xlabel", "t (s)");
    model.result("pg7").set("xlabelactive", true);
    model.result("pg7").set("ylabel", "T (K)");
    model.result("pg7").set("ylabelactive", true);
    model.result("pg7").set("showlegends", false);
    model.result("pg7").set("windowtitle", "Probe Plot 2");
    model.result("pg7").create("tblp1", "Table");
    model.result("pg7").feature("tblp1").label("Probe Table Graph 1");
    model.result("pg7").feature("tblp1").set("table", "tbl1");
    model.result("pg7").feature("tblp1").set("plotcolumninput", "manual");
    model.result("pg7").feature("tblp1").set("legend", true);
    model.result("pg7").feature().remove("tblp2");
    model.result("pg8").label("Density");
    model.result("pg8").feature("surf1").set("resolution", "normal");
    model.result("pg9").label("Heating");
    model.result("pg9").feature("surf1").set("resolution", "normal");
    model.result().export("anim1").label("T and Vf2");
    model.result().export("anim1")
         .set("giffilename", "/nas/longleaf/home/mshourya/workspace/SLS/Particle/Results/ParticleModel/heating+Vf2.gif");
    model.result().export("anim1").set("fps", 6);
    model.result().export("anim1").set("framesel", "all");
    model.result().export("anim1").set("size", "current");
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
    model.result().export("anim2").label("p and Vf2");
    model.result().export("anim2").set("plotgroup", "pg2");
    model.result().export("anim2")
         .set("giffilename", "/nas/longleaf/home/mshourya/workspace/SLS/Particle/Results/ParticleModel/Vf2+p.gif");
    model.result().export("anim2").set("framesel", "all");
    model.result().export("anim2").set("size", "current");
    model.result().export("anim2").set("fontsize", "9");
    model.result().export("anim2").set("colortheme", "globaltheme");
    model.result().export("anim2").set("customcolor", new double[]{1, 1, 1});
    model.result().export("anim2").set("background", "color");
    model.result().export("anim2").set("gltfincludelines", "on");
    model.result().export("anim2").set("title1d", "on");
    model.result().export("anim2").set("legend1d", "on");
    model.result().export("anim2").set("logo1d", "on");
    model.result().export("anim2").set("options1d", "on");
    model.result().export("anim2").set("title2d", "on");
    model.result().export("anim2").set("legend2d", "on");
    model.result().export("anim2").set("logo2d", "on");
    model.result().export("anim2").set("options2d", "on");
    model.result().export("anim2").set("title3d", "on");
    model.result().export("anim2").set("legend3d", "on");
    model.result().export("anim2").set("logo3d", "on");
    model.result().export("anim2").set("options3d", "off");
    model.result().export("anim2").set("axisorientation", "on");
    model.result().export("anim2").set("grid", "on");
    model.result().export("anim2").set("axes1d", "on");
    model.result().export("anim2").set("axes2d", "on");
    model.result().export("anim2").set("showgrid", "on");
    model.result().export("anim3").label("U and Vf2");
    model.result().export("anim3").set("plotgroup", "pg3");
    model.result().export("anim3")
         .set("giffilename", "/nas/longleaf/home/mshourya/workspace/SLS/Particle/Results/ParticleModel/U+Vf2.gif");
    model.result().export("anim3").set("framesel", "all");
    model.result().export("anim3").set("fontsize", "9");
    model.result().export("anim3").set("colortheme", "globaltheme");
    model.result().export("anim3").set("customcolor", new double[]{1, 1, 1});
    model.result().export("anim3").set("background", "color");
    model.result().export("anim3").set("gltfincludelines", "on");
    model.result().export("anim3").set("title1d", "on");
    model.result().export("anim3").set("legend1d", "on");
    model.result().export("anim3").set("logo1d", "on");
    model.result().export("anim3").set("options1d", "on");
    model.result().export("anim3").set("title2d", "on");
    model.result().export("anim3").set("legend2d", "on");
    model.result().export("anim3").set("logo2d", "on");
    model.result().export("anim3").set("options2d", "on");
    model.result().export("anim3").set("title3d", "on");
    model.result().export("anim3").set("legend3d", "on");
    model.result().export("anim3").set("logo3d", "on");
    model.result().export("anim3").set("options3d", "off");
    model.result().export("anim3").set("axisorientation", "on");
    model.result().export("anim3").set("grid", "on");
    model.result().export("anim3").set("axes1d", "on");
    model.result().export("anim3").set("axes2d", "on");
    model.result().export("anim3").set("showgrid", "on");
    model.result().export("anim4").label("mu and Vf2");
    model.result().export("anim4").set("plotgroup", "pg4");
    model.result().export("anim4")
         .set("giffilename", "/nas/longleaf/home/mshourya/workspace/SLS/Particle/Results/ParticleModel/mu+Vf2.gif");
    model.result().export("anim4").set("framesel", "all");
    model.result().export("anim4").set("size", "current");
    model.result().export("anim4").set("fontsize", "9");
    model.result().export("anim4").set("colortheme", "globaltheme");
    model.result().export("anim4").set("customcolor", new double[]{1, 1, 1});
    model.result().export("anim4").set("background", "color");
    model.result().export("anim4").set("gltfincludelines", "on");
    model.result().export("anim4").set("title1d", "on");
    model.result().export("anim4").set("legend1d", "on");
    model.result().export("anim4").set("logo1d", "on");
    model.result().export("anim4").set("options1d", "on");
    model.result().export("anim4").set("title2d", "on");
    model.result().export("anim4").set("legend2d", "on");
    model.result().export("anim4").set("logo2d", "on");
    model.result().export("anim4").set("options2d", "on");
    model.result().export("anim4").set("title3d", "on");
    model.result().export("anim4").set("legend3d", "on");
    model.result().export("anim4").set("logo3d", "on");
    model.result().export("anim4").set("options3d", "off");
    model.result().export("anim4").set("axisorientation", "on");
    model.result().export("anim4").set("grid", "on");
    model.result().export("anim4").set("axes1d", "on");
    model.result().export("anim4").set("axes2d", "on");
    model.result().export("anim4").set("showgrid", "on");

    return model;
  }

  public static void main(String[] args) {
    Model model = run();
    model = run2(model);
    run3(model);
  }

}
