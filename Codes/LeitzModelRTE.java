/*
 * LeitzModelRTE.java
 */

import com.comsol.model.*;
import com.comsol.model.util.*;

/** Model exported on Mar 1 2025, 12:33 by COMSOL 6.1.0.357. */
public class LeitzModelRTE {

  public static Model run() {
    Model model = ModelUtil.create("Model");

    model.modelPath("/work/users/m/s/mshourya/SLS/Particle");

    model.label("LeitzModelRTE.mph");

    model.param().set("R", "8.31446261815324", "Gas constant");
    model.param().set("V", "300 [um]*1000 [um] * 1 [m]");
    model.param().set("TemperatureQuantities", "0", "___________________________________________________");
    model.param().set("Hf", "(33.6 / 95.5) [kJ / g]", "Heat of fusion");
    model.param().set("deltat", "2 [K]", "Phase transition regime size");
    model.param().set("Dp1", "1/(3*kappa)");
    model.param().set("n", "1");
    model.param().set("kappa", "1/pd");
    model.param().set("sigma_b", "5.67e-8 [W/(m^2*K^4)]");
    model.param().set("T0", "300[K]", "Initial temperature");
    model.param().set("Tm", "2896 [K]", "Melting temperature");
    model.param().set("Tamb", "300 [K]", "Ambient temperature");
    model.param().set("B", "1e-6");
    model.param().set("epsilon", "1.5e-5 [m]", "interface thickness");
    model.param().set("LaserParameters", "0", "_______________________");
    model.param().set("vlaser", "500 [mm/s]", "Laser speed");
    model.param().set("Ep", "100 [W]", "Laser power");
    model.param().set("Pw", "100 [um]");
    model.param().set("D", "50 [um]");
    model.param().set("A", "0.12");
    model.param().set("xr", "-10 [mm]");
    model.param().set("yr", "210 [um]");
    model.param().set("xd", "50 [um]");
    model.param().set("ontime", "0.02 [s]");
    model.param().set("offtime", "0.021 [s]");
    model.param().set("pd", "80 [um]");

    model.component().create("comp1", true);

    model.component("comp1").geom().create("geom2", 2);

    model.result().table().create("evl2", "Table");

    model.component("comp1").func().create("gp1", "GaussianPulse");
    model.component("comp1").func().create("rect1", "Rectangle");
    model.component("comp1").func().create("pw1", "Piecewise");
    model.component("comp1").func("gp1").label("delta");
    model.component("comp1").func("gp1").set("funcname", "delta");
    model.component("comp1").func("gp1").set("sigma", "deltat");
    model.component("comp1").func("rect1").label("onoff");
    model.component("comp1").func("rect1").set("funcname", "onoff");
    model.component("comp1").func("rect1").set("lower", "ontime");
    model.component("comp1").func("rect1").set("upper", "offtime");
    model.component("comp1").func("rect1").set("smooth", "0.0001");
    model.component("comp1").func("pw1").label("step");
    model.component("comp1").func("pw1").set("funcname", "step");
    model.component("comp1").func("pw1").set("smooth", "contd2");
    model.component("comp1").func("pw1").set("smoothzone", "0.000001");
    model.component("comp1").func("pw1")
         .set("pieces", new String[][]{{"0", "ontime-0.001", "0.001"}, {"ontime-0.001", "offtime+0.001", "1e-5"}, {"offtime+0.001", "10000", "10"}});
    model.component("comp1").func("pw1").set("argunit", "s");
    model.component("comp1").func("pw1").set("fununit", "s");

    model.component("comp1").mesh().create("mesh3");

    model.component("comp1").geom("geom2").label("Geometry 1");
    model.component("comp1").geom("geom2").lengthUnit("\u00b5m");
    model.component("comp1").geom("geom2").create("c1", "Circle");
    model.component("comp1").geom("geom2").feature("c1").set("pos", new int[]{-45, 85});
    model.component("comp1").geom("geom2").feature("c1").set("r", 50);
    model.component("comp1").geom("geom2").create("c2", "Circle");
    model.component("comp1").geom("geom2").feature("c2").active(false);
    model.component("comp1").geom("geom2").feature("c2").set("pos", new int[]{-10, -15});
    model.component("comp1").geom("geom2").feature("c2").set("r", 5);
    model.component("comp1").geom("geom2").create("arr1", "Array");
    model.component("comp1").geom("geom2").feature("arr1").set("fullsize", new int[]{10, 2});
    model.component("comp1").geom("geom2").feature("arr1").set("displ", new int[]{100, 100});
    model.component("comp1").geom("geom2").feature("arr1").selection("input").set("c1");
    model.component("comp1").geom("geom2").create("r1", "Rectangle");
    model.component("comp1").geom("geom2").feature("r1").set("pos", new int[]{-100, -70});
    model.component("comp1").geom("geom2").feature("r1").set("size", new int[]{1010, 400});
    model.component("comp1").geom("geom2").create("r2", "Rectangle");
    model.component("comp1").geom("geom2").feature("r2").set("pos", new int[]{-100, -70});
    model.component("comp1").geom("geom2").feature("r2").set("size", new int[]{1010, 100});
    model.component("comp1").geom("geom2").run();

    model.component("comp1").variable().create("var1");
    model.component("comp1").variable("var1").set("Tref", "T0", "Reference temperature");
    model.component("comp1").variable("var1").set("Ed", "Ep / (Pw*(pi*0.25*D^2))");
    model.component("comp1").variable("var1").set("G_space", "exp(-1* ( (x  - (xr+vlaser*t) )^2) / (2 * xd^2 ) )");
    model.component("comp1").variable("var1").set("BL", "exp((y - yr)/pd)");
    model.component("comp1").variable("var1").set("Qlaser", "A*Ed*G_space*BL*onoff(t) * (1 + phipf ) / 2");
    model.component("comp1").variable("var1").set("Qsb", "B*(Tamb^4 - T^4)");
    model.component("comp1").variable("var1").set("U", "sqrt(u^2 + v^2)", "Velocity magnitude");
    model.component("comp1").variable("var1")
         .set("alpha", "mat1.def.alpha_p(1[atm],T0)*(1 - phipf)/2 + mat2.def.alpha_iso*(1+phipf)/2");
    model.component("comp1").variable("var1").set("kheat", "mat1.def.k(T0)*(1-phipf)/2 + mat2.def.k_iso*(1+phipf)/2");
    model.component("comp1").variable("var1")
         .set("rho", "mat1.def.rho(1[atm], T0)*(1-phipf)/2 + mat2.def.rho*(1+phipf)/2", "Density");
    model.component("comp1").variable("var1")
         .set("Cp", "mat1.def.Cp(T0)*(1-phipf)/2 + mat2.def.Cp(T)*(1+phipf)/2", "Heat capacity at constant pressure");
    model.component("comp1").variable("var1")
         .set("mu", "mat1.def.eta(T0)*(1-phipf)/2 + mat2.def.mu(T)*(1+phipf)/2", "Dynamic viscosity");
    model.component("comp1").variable("var1").set("sigma", "mat2.def.sigma(T)", "Surface tension coefficient");
    model.component("comp1").variable("var1").set("lambda", "3*epsilon*sigma/sqrt(8)");
    model.component("comp1").variable("var1").set("Fstx", "lambda/(epsilon^2) * psi * phipfx");
    model.component("comp1").variable("var1").set("Fsty", "lambda/(epsilon^2) * psi * phipfy");
    model.component("comp1").variable("var1").set("Vmat", "integratedomain( (1 + phipf) / 2 ) * 1 [m]");
    model.component("comp1").variable("var1").set("porosity", "Vmat/V");
    model.component("comp1").variable("var1").set("Qr", "kappa*(Grad - 4*pi*Ib)");
    model.component("comp1").variable("var1").set("Ib", "n^2*sigma_b*T^4/pi");

    model.component("comp1").material().create("mat1", "Common");
    model.component("comp1").material().create("mat2", "Common");
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
    model.component("comp1").material("mat2").propertyGroup("def").func().create("an1", "Analytic");
    model.component("comp1").material("mat2").propertyGroup("def").func().create("an2", "Analytic");
    model.component("comp1").material("mat2").propertyGroup("def").func().create("an3", "Analytic");
    model.component("comp1").material("mat2").propertyGroup().create("Enu", "Young's modulus and Poisson's ratio");
    model.component("comp1").material("mat2").propertyGroup().create("Murnaghan", "Murnaghan");

    model.component("comp1").cpl().create("intop1", "Integration");
    model.component("comp1").cpl("intop1").selection().all();

    model.component("comp1").physics().create("spf", "LaminarFlow", "geom2");
    model.component("comp1").physics("spf").create("vf1", "VolumeForce", 2);
    model.component("comp1").physics("spf").feature("vf1").selection().all();
    model.component("comp1").physics("spf").create("vf2", "VolumeForce", 2);
    model.component("comp1").physics("spf").feature("vf2").selection().all();
    model.component("comp1").physics("spf").create("out1", "OutletBoundary", 1);
    model.component("comp1").physics("spf").feature("out1").selection().set(5);
    model.component("comp1").physics("spf").create("constr1", "PointwiseConstraint", 1);
    model.component("comp1").physics("spf").feature("constr1").selection().set(5);
    model.component("comp1").physics().create("pf", "PhaseField", "geom2");
    model.component("comp1").physics("pf").feature("initfluid2").selection()
         .set(1, 3, 4, 6, 7, 9, 10, 12, 13, 15, 16, 18, 19, 21, 22, 24, 25, 27, 28, 30, 31);
    model.component("comp1").physics("pf").create("out1", "Outlet", 1);
    model.component("comp1").physics("pf").feature("out1").selection().set(5);
    model.component("comp1").physics().create("ht", "HeatTransferInFluids", "geom2");
    model.component("comp1").physics("ht").create("hs1", "HeatSource", 2);
    model.component("comp1").physics("ht").feature("hs1").selection()
         .set(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31);
    model.component("comp1").physics("ht").create("hs2", "HeatSource", 2);
    model.component("comp1").physics("ht").feature("hs2").selection()
         .set(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31);
    model.component("comp1").physics().create("rteeq", "PoissonEquation", "geom2");
    model.component("comp1").physics("rteeq").identifier("rteeq");
    model.component("comp1").physics("rteeq").field("dimensionless").field("Grad");
    model.component("comp1").physics("rteeq").prop("Units").set("DependentVariableQuantity", "none");
    model.component("comp1").physics("rteeq").prop("Units").set("CustomDependentVariableUnit", "W*m^-2");
    model.component("comp1").physics("rteeq").create("flux1", "FluxBoundary", 1);
    model.component("comp1").physics("rteeq").feature("flux1").selection().set(5);

    model.component("comp1").mesh("mesh3").create("ftri1", "FreeTri");
    model.component("comp1").mesh("mesh3").feature("ftri1").create("size1", "Size");
    model.component("comp1").mesh("mesh3").feature("ftri1").create("dis2", "Distribution");
    model.component("comp1").mesh("mesh3").feature("ftri1").create("dis3", "Distribution");
    model.component("comp1").mesh("mesh3").feature("ftri1").feature("dis2").selection().set(2);
    model.component("comp1").mesh("mesh3").feature("ftri1").feature("dis3").selection().set(5);

    model.component("comp1").probe().create("point1", "Point");
    model.component("comp1").probe().create("point2", "Point");
    model.component("comp1").probe().create("point3", "Point");
    model.component("comp1").probe().create("point4", "Point");
    model.component("comp1").probe().create("point5", "Point");
    model.component("comp1").probe().create("point6", "Point");
    model.component("comp1").probe().create("point7", "Point");
    model.component("comp1").probe().create("point8", "Point");
    model.component("comp1").probe().create("point9", "Point");
    model.component("comp1").probe().create("point10", "Point");
    model.component("comp1").probe("point1").selection().set(8);
    model.component("comp1").probe("point2").selection().set(13);
    model.component("comp1").probe("point3").selection().set(18);
    model.component("comp1").probe("point4").selection().set(23);
    model.component("comp1").probe("point5").selection().set(28);
    model.component("comp1").probe("point6").selection().set(33);
    model.component("comp1").probe("point7").selection().set(38);
    model.component("comp1").probe("point8").selection().set(43);
    model.component("comp1").probe("point9").selection().set(48);
    model.component("comp1").probe("point10").selection().set(53);

    model.result().table("evl2").label("Evaluation 2D");
    model.result().table("evl2").comments("Interactive 2D values");

    model.component("comp1").view("view1").axis().set("xmin", -299.34442138671875);
    model.component("comp1").view("view1").axis().set("xmax", 1208.1097412109375);
    model.component("comp1").view("view1").axis().set("ymin", -367.45831298828125);
    model.component("comp1").view("view1").axis().set("ymax", 630.0618896484375);

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
    model.component("comp1").material("mat2").label("Molybdenum");
    model.component("comp1").material("mat2").set("family", "custom");
    model.component("comp1").material("mat2").set("customspecular", new double[]{0.7843137254901961, 1, 1});
    model.component("comp1").material("mat2")
         .set("customdiffuse", new double[]{0.7843137254901961, 0.7843137254901961, 0.7843137254901961});
    model.component("comp1").material("mat2").set("fresnel", 0.3);
    model.component("comp1").material("mat2").set("roughness", 0.1);
    model.component("comp1").material("mat2").set("metallic", 0);
    model.component("comp1").material("mat2").set("pearl", 0);
    model.component("comp1").material("mat2").set("diffusewrap", 0);
    model.component("comp1").material("mat2").set("clearcoat", 0);
    model.component("comp1").material("mat2").set("reflectance", 0);
    model.component("comp1").material("mat2").propertyGroup("def").func("an1").label("Cp");
    model.component("comp1").material("mat2").propertyGroup("def").func("an1").set("funcname", "Cp");
    model.component("comp1").material("mat2").propertyGroup("def").func("an1")
         .set("expr", "(34.2 + 1.13 *10^(-3) * (T - Tm ) ) * ( 1 / 95.5 ) + Hf * delta( T - Tm )");
    model.component("comp1").material("mat2").propertyGroup("def").func("an1").set("args", new String[]{"T"});
    model.component("comp1").material("mat2").propertyGroup("def").func("an1").set("fununit", "J/(g*K)");
    model.component("comp1").material("mat2").propertyGroup("def").func("an1").set("argunit", new String[]{"K"});
    model.component("comp1").material("mat2").propertyGroup("def").func("an1")
         .set("plotargs", new String[][]{{"T", "300", "350"}});
    model.component("comp1").material("mat2").propertyGroup("def").func("an2").label("sigma");
    model.component("comp1").material("mat2").propertyGroup("def").func("an2").set("funcname", "sigma");
    model.component("comp1").material("mat2").propertyGroup("def").func("an2")
         .set("expr", "(2.29 * 10^3 - 0.26 * (T - Tm) ) * 0.001");
    model.component("comp1").material("mat2").propertyGroup("def").func("an2").set("args", new String[]{"T"});
    model.component("comp1").material("mat2").propertyGroup("def").func("an2").set("fununit", "N/m");
    model.component("comp1").material("mat2").propertyGroup("def").func("an2").set("argunit", new String[]{"K"});
    model.component("comp1").material("mat2").propertyGroup("def").func("an2")
         .set("plotargs", new String[][]{{"T", "0", "3000"}});
    model.component("comp1").material("mat2").propertyGroup("def").func("an3").label("mu");
    model.component("comp1").material("mat2").propertyGroup("def").func("an3").set("funcname", "mu");
    model.component("comp1").material("mat2").propertyGroup("def").func("an3")
         .set("expr", "0.27 * exp(73e3 / (R*T) ) * 0.001");
    model.component("comp1").material("mat2").propertyGroup("def").func("an3").set("args", new String[]{"T"});
    model.component("comp1").material("mat2").propertyGroup("def").func("an3").set("fununit", "Pa*s");
    model.component("comp1").material("mat2").propertyGroup("def").func("an3").set("argunit", new String[]{"K"});
    model.component("comp1").material("mat2").propertyGroup("def").func("an3")
         .set("plotargs", new String[][]{{"T", "1000", "3000"}});
    model.component("comp1").material("mat2").propertyGroup("def")
         .set("thermalexpansioncoefficient", new String[]{"5.3e-5[1/K]", "0", "0", "0", "5.3e-5[1/K]", "0", "0", "0", "5.3e-5[1/K]"});
    model.component("comp1").material("mat2").propertyGroup("def").set("density", "10200[kg/m^3]");
    model.component("comp1").material("mat2").propertyGroup("def").set("heatcapacity", "Cp(T)");
    model.component("comp1").material("mat2").propertyGroup("def")
         .set("thermalconductivity", new String[]{"138[W/(m*K)]", "0", "0", "0", "138[W/(m*K)]", "0", "0", "0", "138[W/(m*K)]"});
    model.component("comp1").material("mat2").propertyGroup("Enu").set("E", "312[GPa]");
    model.component("comp1").material("mat2").propertyGroup("Enu").set("nu", "0.31");
    model.component("comp1").material("mat2").propertyGroup("Murnaghan").set("l", "-300[GPa]");
    model.component("comp1").material("mat2").propertyGroup("Murnaghan").set("m", "-850[GPa]");
    model.component("comp1").material("mat2").propertyGroup("Murnaghan").set("n", "-910[GPa]");

    model.component("comp1").cpl("intop1").label("integrate");
    model.component("comp1").cpl("intop1").set("opname", "integratedomain");

    model.component("comp1").coordSystem("sys1").set("name", "sys2");

    model.component("comp1").physics("spf").prop("ShapeProperty").set("order_fluid", 3);
    model.component("comp1").physics("spf").prop("PhysicalModelProperty")
         .set("Compressibility", "CompressibleMALT03");
    model.component("comp1").physics("spf").prop("PhysicalModelProperty").set("IncludeGravity", true);
    model.component("comp1").physics("spf").prop("PhysicalModelProperty").set("UseReducedPressure", true);
    model.component("comp1").physics("spf").prop("PhysicalModelProperty").set("Tref", "T0");
    model.component("comp1").physics("spf").prop("TurbulenceModelProperty").set("TurbulenceModel", "SST");
    model.component("comp1").physics("spf").prop("InconsistentStabilization").set("IsotropicDiffusion", true);
    model.component("comp1").physics("spf").prop("InconsistentStabilization").set("delid", 0.1);
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
    model.component("comp1").physics("spf").feature("vf1")
         .set("F", new String[][]{{"0"}, {"rho * g_const * alpha * (T - Tref)"}, {"0"}});
    model.component("comp1").physics("spf").feature("vf1").active(false);
    model.component("comp1").physics("spf").feature("vf2").set("F", new String[][]{{"Fstx"}, {"Fsty"}, {"0"}});
    model.component("comp1").physics("spf").feature("out1").set("p0", "1 [atm]");
    model.component("comp1").physics("spf").feature("out1").set("AverageTotalPressure", false);
    model.component("comp1").physics("spf").feature("out1").set("PressureType", "TotalPressure");
    model.component("comp1").physics("spf").feature("constr1").set("constraintExpression", "u");
    model.component("comp1").physics("spf").feature("constr1").set("constraintMethod", "nodal");
    model.component("comp1").physics("pf").prop("ShapeProperty").set("order_phasefield", 3);
    model.component("comp1").physics("pf").feature("pfm1").set("epsilon_pf", "epsilon");
    model.component("comp1").physics("pf").feature("pfm1").set("chiOption", "velocity");
    model.component("comp1").physics("pf").feature("pfm1").set("chi", "0.001 * epsilon");
    model.component("comp1").physics("pf").feature("pfm1").set("U", "1e-12");
    model.component("comp1").physics("pf").feature("pfm1").set("u", new String[][]{{"u"}, {"v"}, {"0"}});
    model.component("comp1").physics("pf").feature("pfm1").set("sigma", "sigma");
    model.component("comp1").physics("pf").feature("ww1").set("surfaceE1", "0.001[J/m^2]");
    model.component("comp1").physics("pf").feature("ww1").set("surfaceE2", "0.001[J/m^2]");
    model.component("comp1").physics("ht").prop("ShapeProperty").set("order_temperature", 3);
    model.component("comp1").physics("ht").prop("PhysicalModelProperty").set("dz", "100[um]");
    model.component("comp1").physics("ht").prop("PhysicalModelProperty").set("Tref", 300);
    model.component("comp1").physics("ht").prop("InconsistentStabilization").set("HeatIsotropicDiffusion", true);
    model.component("comp1").physics("ht").feature("fluid1").set("Cp_mat", "userdef");
    model.component("comp1").physics("ht").feature("fluid1").set("Cp", "Cp");
    model.component("comp1").physics("ht").feature("fluid1").set("rho_mat", "userdef");
    model.component("comp1").physics("ht").feature("fluid1").set("rho", "rho");
    model.component("comp1").physics("ht").feature("fluid1").set("k_mat", "userdef");
    model.component("comp1").physics("ht").feature("fluid1")
         .set("k", new String[][]{{"kheat"}, {"0"}, {"0"}, {"0"}, {"kheat"}, {"0"}, {"0"}, {"0"}, {"kheat"}});
    model.component("comp1").physics("ht").feature("fluid1").set("u", new String[][]{{"u"}, {"v"}, {"0"}});
    model.component("comp1").physics("ht").feature("fluid1").set("minput_pressure_src", "root.comp1.spf.pA");
    model.component("comp1").physics("ht").feature("fluid1").set("minput_pressure", "p");
    model.component("comp1").physics("ht").feature("init1").set("Tinit", "T0");
    model.component("comp1").physics("ht").feature("hs1").set("Q0", "Qlaser");
    model.component("comp1").physics("ht").feature("hs1").set("materialType", "from_mat");
    model.component("comp1").physics("ht").feature("hs2").set("Q0", "Qr");
    model.component("comp1").physics("ht").feature("hs2").set("materialType", "nonSolid");
    model.component("comp1").physics("rteeq").label("RTE");
    model.component("comp1").physics("rteeq").prop("EquationForm").set("form", "Automatic");
    model.component("comp1").physics("rteeq").prop("Units").set("CustomSourceTermUnit", "W*m^-3");
    model.component("comp1").physics("rteeq").feature("peq1").set("f", "-Qr");
    model.component("comp1").physics("rteeq").feature("peq1").set("c", new String[][]{{"Dp1", "0", "0", "Dp1"}});
    model.component("comp1").physics("rteeq").feature("flux1").set("g", "0.5*(4*pi*Ib - Grad)");

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

    return model;
  }

  public static Model run2(Model model) {
    model.component("comp1").mesh("mesh3").feature("ftri1").feature("size1").set("hmin", 15);
    model.component("comp1").mesh("mesh3").feature("ftri1").feature("size1").set("hminactive", true);
    model.component("comp1").mesh("mesh3").feature("ftri1").feature("size1").set("hcurve", 1);
    model.component("comp1").mesh("mesh3").feature("ftri1").feature("size1").set("hcurveactive", true);
    model.component("comp1").mesh("mesh3").feature("ftri1").feature("size1").set("hnarrow", 0.1);
    model.component("comp1").mesh("mesh3").feature("ftri1").feature("size1").set("hnarrowactive", true);
    model.component("comp1").mesh("mesh3").feature("ftri1").feature("size1").set("hgrad", 5);
    model.component("comp1").mesh("mesh3").feature("ftri1").feature("size1").set("hgradactive", true);
    model.component("comp1").mesh("mesh3").feature("ftri1").feature("dis2").active(false);
    model.component("comp1").mesh("mesh3").feature("ftri1").feature("dis2").set("numelem", 20);
    model.component("comp1").mesh("mesh3").feature("ftri1").feature("dis3").active(false);
    model.component("comp1").mesh("mesh3").feature("ftri1").feature("dis3").set("numelem", 20);
    model.component("comp1").mesh("mesh3").run();

    model.component("comp1").physics("ht").feature("fluid1").set("minput_pressure_src", "root.comp1.spf.pA");

    model.study().create("std1");
    model.study("std1").create("phasei", "PhaseInitialization");
    model.study("std1").create("time", "Transient");
    model.study("std1").feature("phasei")
         .set("activate", new String[]{"spf", "off", "pf", "on", "ht", "on", "rteeq", "on", "frame:spatial1", "on", 
         "frame:material1", "on"});

    model.sol().create("sol4");
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
    model.sol("sol4").feature("s1").feature().remove("fcDef");
    model.sol("sol4").feature("t1").create("fc1", "FullyCoupled");
    model.sol("sol4").feature("t1").create("d1", "Direct");
    model.sol("sol4").feature("t1").create("i1", "Iterative");
    model.sol("sol4").feature("t1").create("i2", "Iterative");
    model.sol("sol4").feature("t1").create("i3", "Iterative");
    model.sol("sol4").feature("t1").create("se1", "Segregated");
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
    model.sol("sol4").feature("t1").feature("se1").create("ss1", "SegregatedStep");
    model.sol("sol4").feature("t1").feature("se1").create("ss2", "SegregatedStep");
    model.sol("sol4").feature("t1").feature("se1").create("rteeq1", "SegregatedStep");
    model.sol("sol4").feature("t1").feature().remove("fcDef");

    model.result().dataset().create("dset4", "Solution");
    model.result().dataset().create("an1_ds1", "Grid1D");
    model.result().dataset().create("dset5", "Solution");
    model.result().dataset().create("dset6", "Solution");
    model.result().dataset("dset1").set("solution", "none");
    model.result().dataset("dset2").set("solution", "none");
    model.result().dataset("dset4").set("solution", "none");
    model.result().dataset("an1_ds1").set("data", "none");
    model.result().dataset("dset6").set("solution", "sol5");
    model.result().numerical().create("gev1", "EvalGlobal");
    model.result().numerical("gev1").set("probetag", "none");
    model.result().create("pg1", "PlotGroup2D");
    model.result().create("pg10", "PlotGroup2D");
    model.result().create("pg2", "PlotGroup2D");
    model.result().create("pg3", "PlotGroup2D");
    model.result().create("pg4", "PlotGroup2D");
    model.result().create("pg6", "PlotGroup1D");
    model.result().create("pg8", "PlotGroup2D");
    model.result().create("pg11", "PlotGroup1D");
    model.result().create("pg12", "PlotGroup2D");
    model.result("pg1").set("data", "dset5");
    model.result("pg1").create("surf1", "Surface");
    model.result("pg1").create("con1", "Contour");
    model.result("pg1").feature("surf1").set("expr", "Qr");
    model.result("pg1").feature("con1").set("expr", "pf.Vf2");
    model.result("pg10").set("data", "dset5");
    model.result("pg10").create("surf1", "Surface");
    model.result("pg10").create("con1", "Contour");
    model.result("pg10").feature("surf1").set("expr", "T");
    model.result("pg10").feature("con1").set("expr", "pf.Vf2");
    model.result("pg2").set("data", "dset5");
    model.result("pg2").create("surf1", "Surface");
    model.result("pg2").create("con1", "Contour");
    model.result("pg2").create("arws2", "ArrowSurface");
    model.result("pg2").create("arws3", "ArrowSurface");
    model.result("pg2").feature("surf1").set("expr", "comp1.phipf");
    model.result("pg2").feature("con1").set("expr", "pf.Vf2");
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
    model.result("pg8").set("data", "dset5");
    model.result("pg8").create("surf1", "Surface");
    model.result("pg8").feature("surf1").set("expr", "rho");
    model.result("pg11").set("data", "dset5");
    model.result("pg11").create("glob1", "Global");
    model.result("pg12").set("data", "dset5");
    model.result("pg12").create("surf1", "Surface");
    model.result("pg12").create("con1", "Contour");
    model.result("pg12").create("arws2", "ArrowSurface");
    model.result("pg12").create("arws3", "ArrowSurface");
    model.result("pg12").feature("surf1").set("expr", "comp1.phipf");
    model.result("pg12").feature("con1").set("expr", "pf.Vf2");
    model.result().export().create("anim1", "Animation");
    model.result().export().create("anim2", "Animation");
    model.result().export().create("anim3", "Animation");
    model.result().export().create("anim4", "Animation");
    model.result().export().create("anim5", "Animation");

    model.component("comp1").probe("point1").genResult(null);
    model.component("comp1").probe("point2").genResult(null);
    model.component("comp1").probe("point3").genResult(null);
    model.component("comp1").probe("point4").genResult(null);
    model.component("comp1").probe("point5").genResult(null);
    model.component("comp1").probe("point6").genResult(null);
    model.component("comp1").probe("point7").genResult(null);
    model.component("comp1").probe("point8").genResult(null);
    model.component("comp1").probe("point9").genResult(null);
    model.component("comp1").probe("point10").genResult(null);

    model.study("std1").feature("time").set("tlist", "range(0,0.0001,0.2)");
    model.study("std1").feature("time").set("useinitsol", true);
    model.study("std1").feature("time").set("initstudy", "std1");
    model.study("std1").feature("time").set("solnum", "auto");

    model.sol("sol4").attach("std1");
    model.sol("sol4").feature("st1").label("Compile Equations: Phase Initialization");
    model.sol("sol4").feature("v1").label("Dependent Variables 1.1");
    model.sol("sol4").feature("v1").set("clist", new String[]{"4.0E-5[s]"});
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
    model.sol("sol4").feature("su1").label("Solution Store 1.1");
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
    model.sol("sol4").feature("v2").set("clist", new String[]{"range(0,0.001,0.04)", "4.0E-5[s]"});
    model.sol("sol4").feature("v2").feature("comp1_phipf").set("scalemethod", "manual");
    model.sol("sol4").feature("v2").feature("comp1_phipf").set("scaleval", 1);
    model.sol("sol4").feature("v2").feature("comp1_psi").set("resscalemethod", "manual");
    model.sol("sol4").feature("v2").feature("comp1_psi").set("resscaleval", 10);
    model.sol("sol4").feature("v2").feature("comp1_T").set("scalemethod", "manual");
    model.sol("sol4").feature("v2").feature("comp1_T").set("scaleval", 5000);
    model.sol("sol4").feature("v2").feature("comp1_u").set("scalemethod", "auto");
    model.sol("sol4").feature("t1").label("Time-Dependent Solver 2.1");
    model.sol("sol4").feature("t1").set("tlist", "range(0,0.001,0.04)");
    model.sol("sol4").feature("t1").set("atolglobalfactor", 0.01);
    model.sol("sol4").feature("t1")
         .set("atolmethod", new String[]{"comp1_GI", "global", "comp1_p", "scaled", "comp1_phipf", "global", "comp1_psi", "global", "comp1_T", "global", 
         "comp1_u", "global", "comp1_Grad", "global"});
    model.sol("sol4").feature("t1")
         .set("atolfactor", new String[]{"comp1_GI", "0.1", "comp1_p", "1", "comp1_phipf", "0.1", "comp1_psi", "0.1", "comp1_T", "0.1", 
         "comp1_u", "0.1", "comp1_Grad", "0.1"});
    model.sol("sol4").feature("t1").set("initialstepbdf", "1e-15");
    model.sol("sol4").feature("t1").set("maxstepconstraintbdf", "expr");
    model.sol("sol4").feature("t1").set("maxstepexpressionbdf", "comp1.step(t)");
    model.sol("sol4").feature("t1").set("maxorder", 1);
    model.sol("sol4").feature("t1").set("stabcntrl", true);
    model.sol("sol4").feature("t1").set("bwinitstepfrac", 0.01);
    model.sol("sol4").feature("t1").set("bwinitfactor", 200);
    model.sol("sol4").feature("t1").set("estrat", "exclude");
    model.sol("sol4").feature("t1").set("rescaleafterinitbw", true);
    model.sol("sol4").feature("t1").set("reacf", false);
    model.sol("sol4").feature("t1").set("storeudot", false);
    model.sol("sol4").feature("t1").feature("dDef").label("Direct 2");
    model.sol("sol4").feature("t1").feature("aDef").label("Advanced 1");
    model.sol("sol4").feature("t1").feature("aDef").set("blocksize", 500);
    model.sol("sol4").feature("t1").feature("aDef").set("blocksizeactive", true);
    model.sol("sol4").feature("t1").feature("aDef").set("assemloc", false);
    model.sol("sol4").feature("t1").feature("fc1").active(true);
    model.sol("sol4").feature("t1").feature("fc1").label("Fully Coupled 1.1");
    model.sol("sol4").feature("t1").feature("fc1").set("linsolver", "d1");
    model.sol("sol4").feature("t1").feature("fc1").set("maxiter", 5);
    model.sol("sol4").feature("t1").feature("fc1").set("ntolfact", 0.1);
    model.sol("sol4").feature("t1").feature("fc1").set("termonres", "both");
    model.sol("sol4").feature("t1").feature("fc1").set("damp", "0.9");
    model.sol("sol4").feature("t1").feature("fc1").set("jtech", "once");
    model.sol("sol4").feature("t1").feature("fc1").set("stabacc", "aacc");
    model.sol("sol4").feature("t1").feature("fc1").set("aaccdim", 5);
    model.sol("sol4").feature("t1").feature("fc1").set("aaccmix", 0.9);
    model.sol("sol4").feature("t1").feature("fc1").set("aaccdelay", 1);
    model.sol("sol4").feature("t1").feature("d1").label("Direct, heat transfer variables (ht) (merged)");
    model.sol("sol4").feature("t1").feature("d1").set("linsolver", "pardiso");
    model.sol("sol4").feature("t1").feature("d1").set("pivotperturb", 1.0E-13);
    model.sol("sol4").feature("t1").feature("d1").set("nliniterrefine", true);
    model.sol("sol4").feature("t1").feature("i1").label("AMG");
    model.sol("sol4").feature("t1").feature("i1").set("irestol", 0.1);
    model.sol("sol4").feature("t1").feature("i1").set("maxlinit", 100);
    model.sol("sol4").feature("t1").feature("i1").set("rhob", 20);
    model.sol("sol4").feature("t1").feature("i1").set("maxilinit", 1000);
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
    model.sol("sol4").feature("t1").feature("i1").feature("mg1").feature("cs").feature("d1")
         .set("pardreorder", "ndmt");
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
    model.sol("sol4").feature("t1").feature("se1").label("Segregated 1.1");
    model.sol("sol4").feature("t1").feature("se1").set("maxsegiter", 15);
    model.sol("sol4").feature("t1").feature("se1").set("segtermonres", "both");
    model.sol("sol4").feature("t1").feature("se1").feature("ssDef").label("Fluid Variables");
    model.sol("sol4").feature("t1").feature("se1").feature("ssDef").set("segvar", new String[]{"comp1_p", "comp1_u"});
    model.sol("sol4").feature("t1").feature("se1").feature("ssDef").set("linsolver", "i1");
    model.sol("sol4").feature("t1").feature("se1").feature("ssDef").set("subdamp", "0.9");
    model.sol("sol4").feature("t1").feature("se1").feature("ssDef").set("subjtech", "once");
    model.sol("sol4").feature("t1").feature("se1").feature("ss1").label("Phase Field");
    model.sol("sol4").feature("t1").feature("se1").feature("ss1")
         .set("segvar", new String[]{"comp1_GI", "comp1_phipf", "comp1_psi"});
    model.sol("sol4").feature("t1").feature("se1").feature("ss1").set("linsolver", "i1");
    model.sol("sol4").feature("t1").feature("se1").feature("ss1").set("subdamp", "0.9");
    model.sol("sol4").feature("t1").feature("se1").feature("ss1").set("subjtech", "once");
    model.sol("sol4").feature("t1").feature("se1").feature("ss2").label("Temperature");

    return model;
  }

  public static Model run3(Model model) {
    model.sol("sol4").feature("t1").feature("se1").feature("ss2").set("segvar", new String[]{"comp1_T"});
    model.sol("sol4").feature("t1").feature("se1").feature("ss2").set("linsolver", "i1");
    model.sol("sol4").feature("t1").feature("se1").feature("ss2").set("subdamp", "0.9");
    model.sol("sol4").feature("t1").feature("se1").feature("ss2").set("subjtech", "once");
    model.sol("sol4").feature("t1").feature("se1").feature("rteeq1").label("Segregated Step 1");
    model.sol("sol4").feature("t1").feature("se1").feature("rteeq1").set("segvar", new String[]{"comp1_Grad"});
    model.sol("sol4").feature("t1").feature("se1").feature("rteeq1").set("linsolver", "d1");
    model.sol("sol4").runAll();

    model.result().dataset("an1_ds1").label("Grid 1D 1a");
    model.result().dataset("an1_ds1").set("function", "all");
    model.result().dataset("an1_ds1").set("par1", "T");
    model.result().dataset("an1_ds1").set("parmin1", 300);
    model.result().dataset("an1_ds1").set("parmax1", 400);
    model.result().dataset("an1_ds1").set("res1", 10000);
    model.result().dataset("an1_ds1").set("distribution", "mixed");
    model.result().dataset().remove("dset7");
    model.result().numerical("gev1").set("expr", new String[]{"porosity"});
    model.result().numerical("gev1").set("unit", new String[]{"1"});
    model.result().numerical("gev1").set("descr", new String[]{""});
    model.result("pg1").label("Heating and Volume Fraction");
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
    model.result("pg2").label("Pressure and Volume Fraction with Arrows for Surface Tension");
    model.result("pg2").set("titletype", "manual");
    model.result("pg2").set("title", "Volume Fraction with Pressure (Pa) contours");
    model.result("pg2").set("paramindicator", "t=eval(t,s) (s)");
    model.result("pg2").feature("surf1").set("coloring", "gradient");
    model.result("pg2").feature("surf1").set("topcolor", "blue");
    model.result("pg2").feature("surf1").set("bottomcolor", "white");
    model.result("pg2").feature("surf1").set("resolution", "normal");
    model.result("pg2").feature("con1").set("coloring", "uniform");
    model.result("pg2").feature("con1").set("color", "magenta");
    model.result("pg2").feature("con1").set("resolution", "normal");
    model.result("pg2").feature("arws2").set("expr", new String[]{"Fstx", "Fsty"});
    model.result("pg2").feature("arws2").set("descr", "");
    model.result("pg2").feature("arws2").set("scale", 3.1394639721073104E-9);
    model.result("pg2").feature("arws2").set("scaleactive", true);
    model.result("pg2").feature("arws2").set("color", "green");
    model.result("pg2").feature("arws3").active(false);
    model.result("pg2").feature("arws3").set("expr", new String[]{"0", "rho * g_const * alpha * (T - Tref)"});
    model.result("pg2").feature("arws3").set("descr", "");
    model.result("pg2").feature("arws3").set("scale", 2.1582861549393718E12);
    model.result("pg2").feature("arws3").set("color", "cyan");
    model.result("pg2").feature("arws3").set("scaleactive", false);
    model.result("pg3").label("Velocity and Volume Fraction");
    model.result("pg3").set("titletype", "manual");
    model.result("pg3").set("title", "Velocity magnitude (m/s) with Volume Fraction contours");
    model.result("pg3").set("paramindicator", "t=eval(t) (s)");
    model.result("pg3").feature("surf1").set("colortable", "Magma");
    model.result("pg3").feature("surf1").set("resolution", "normal");
    model.result("pg3").feature("con1").set("coloring", "uniform");
    model.result("pg3").feature("con1").set("color", "green");
    model.result("pg3").feature("con1").set("resolution", "normal");
    model.result("pg3").feature("arwl1").set("scale", 15335.485005874825);
    model.result("pg3").feature("arwl1").set("scaleactive", false);
    model.result("pg4").label("Viscosity");
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
    model.result("pg6").set("titletype", "manual");
    model.result("pg6").set("title", "Specific Heat of Material (J/(kg*K))");
    model.result("pg6").set("xlabel", "T (K)");
    model.result("pg6").set("xlabelactive", true);
    model.result("pg6").set("ylabel", "J/(kg*K)");
    model.result("pg6").set("ylabelactive", true);
    model.result("pg6").feature("plot1").set("solrepresentation", "solnum");
    model.result("pg6").feature("plot1").set("descr", "cpmat(T)");
    model.result("pg6").feature("plot1").set("xdataexpr", "T");
    model.result("pg6").feature("plot1").set("xdataunit", "m");
    model.result("pg6").feature("plot1").set("xdatadescractive", true);
    model.result("pg6").feature("plot1").set("xdatadescr", "");
    model.result("pg6").feature("plot1").set("lowerbound", 300);
    model.result("pg6").feature("plot1").set("upperbound", 400);
    model.result("pg6").feature("plot1").set("linewidth", "preference");
    model.result("pg8").label("Density");
    model.result("pg8").feature("surf1").set("colortable", "Twilight");
    model.result("pg8").feature("surf1").set("resolution", "normal");
    model.result("pg11").label("Porosity");
    model.result("pg11").set("xlabel", "Time (s)");
    model.result("pg11").set("xlabelactive", false);
    model.result("pg11").feature("glob1").set("expr", new String[]{"Vmat/V"});
    model.result("pg11").feature("glob1").set("unit", new String[]{"1"});
    model.result("pg11").feature("glob1").set("descr", new String[]{""});
    model.result("pg11").feature("glob1").set("linewidth", "preference");
    model.result("pg12").label("Pressure and Volume Fraction with Arrows for Thermal Expansion Force");
    model.result("pg12").set("titletype", "manual");
    model.result("pg12").set("title", "Volume Fraction with Pressure (Pa) contours");
    model.result("pg12").set("paramindicator", "t=eval(t,s) (s)");
    model.result("pg12").feature("surf1").set("coloring", "gradient");
    model.result("pg12").feature("surf1").set("topcolor", "blue");
    model.result("pg12").feature("surf1").set("bottomcolor", "white");
    model.result("pg12").feature("surf1").set("resolution", "normal");
    model.result("pg12").feature("con1").set("coloring", "uniform");
    model.result("pg12").feature("con1").set("color", "magenta");
    model.result("pg12").feature("con1").set("resolution", "normal");
    model.result("pg12").feature("arws2").set("expr", new String[]{"0", "rho * g_const * alpha * (T - Tref)"});
    model.result("pg12").feature("arws2").set("descr", "");
    model.result("pg12").feature("arws2").set("scale", 7.0E-4);
    model.result("pg12").feature("arws2").set("scaleactive", true);
    model.result("pg12").feature("arws2").set("color", "green");
    model.result("pg12").feature("arws3").active(false);
    model.result("pg12").feature("arws3").set("expr", new String[]{"0", "rho * g_const * alpha * (T - Tref)"});
    model.result("pg12").feature("arws3").set("descr", "");
    model.result("pg12").feature("arws3").set("scale", 2.1582861549393718E12);
    model.result("pg12").feature("arws3").set("color", "cyan");
    model.result("pg12").feature("arws3").set("scaleactive", false);
    model.result().remove("pg13");
    model.result().export("anim1").label("Heating and Vf2");
    model.result().export("anim1")
         .set("giffilename", "/work/users/m/s/mshourya/SLS/Particle/Results/Heating+Vf2.gif");
    model.result().export("anim1").set("maxframes", 41);
    model.result().export("anim1").set("size", "current");
    model.result().export("anim1").set("repeat", "forever");
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
    model.result().export("anim2").label("Temperature and Vf2");
    model.result().export("anim2").set("plotgroup", "pg10");
    model.result().export("anim2").set("giffilename", "/work/users/m/s/mshourya/SLS/Particle/Results/T+Vf2.gif");
    model.result().export("anim2").set("maxframes", 41);
    model.result().export("anim2").set("size", "current");
    model.result().export("anim2").set("repeat", "forever");
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
    model.result().export("anim3").label("Velocity and Vf2");
    model.result().export("anim3").set("plotgroup", "pg3");
    model.result().export("anim3").set("giffilename", "/work/users/m/s/mshourya/SLS/Particle/Results/U+Vf2.gif");
    model.result().export("anim3").set("maxframes", 41);
    model.result().export("anim3").set("size", "current");
    model.result().export("anim3").set("repeat", "forever");
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
    model.result().export("anim4").label("Volume Fraction and Arrows for Surface Tension");
    model.result().export("anim4").set("plotgroup", "pg2");
    model.result().export("anim4")
         .set("giffilename", "/work/users/m/s/mshourya/SLS/Particle/Results/p+Vf2+sigma.gif");
    model.result().export("anim4").set("maxframes", 41);
    model.result().export("anim4").set("size", "current");
    model.result().export("anim4").set("repeat", "forever");
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
    model.result().export("anim5").label("Volume Fraction and Arrows for Thermal Expansion Force");
    model.result().export("anim5").set("plotgroup", "pg12");
    model.result().export("anim5")
         .set("giffilename", "/work/users/m/s/mshourya/SLS/Particle/Results/p+Vf2+Falpha.gif");
    model.result().export("anim5").set("maxframes", 41);
    model.result().export("anim5").set("size", "current");
    model.result().export("anim5").set("repeat", "forever");
    model.result().export("anim5").set("fontsize", "9");
    model.result().export("anim5").set("colortheme", "globaltheme");
    model.result().export("anim5").set("customcolor", new double[]{1, 1, 1});
    model.result().export("anim5").set("background", "color");
    model.result().export("anim5").set("gltfincludelines", "on");
    model.result().export("anim5").set("title1d", "on");
    model.result().export("anim5").set("legend1d", "on");
    model.result().export("anim5").set("logo1d", "on");
    model.result().export("anim5").set("options1d", "on");
    model.result().export("anim5").set("title2d", "on");
    model.result().export("anim5").set("legend2d", "on");
    model.result().export("anim5").set("logo2d", "on");
    model.result().export("anim5").set("options2d", "on");
    model.result().export("anim5").set("title3d", "on");
    model.result().export("anim5").set("legend3d", "on");
    model.result().export("anim5").set("logo3d", "on");
    model.result().export("anim5").set("options3d", "off");
    model.result().export("anim5").set("axisorientation", "on");
    model.result().export("anim5").set("grid", "on");
    model.result().export("anim5").set("axes1d", "on");
    model.result().export("anim5").set("axes2d", "on");
    model.result().export("anim5").set("showgrid", "on");

    return model;
  }

  public static void main(String[] args) {
    Model model = run();
    model = run2(model);
    run3(model);
  }

}
