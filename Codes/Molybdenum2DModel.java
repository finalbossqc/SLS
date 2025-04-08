/*
 * Molybdenum2DModel.java
 */

import com.comsol.model.*;
import com.comsol.model.util.*;

/** Model exported on Apr 8 2025, 13:12 by COMSOL 6.1.0.357. */
public class Molybdenum2DModel {

  public static Model run() {
    Model model = ModelUtil.create("Model");

    model.modelPath("/users/m/s/mshourya/SLS/Complete/Tests");

    model.label("Molybdenum2DModel.mph");

    model.param().set("R", "8.31446261815324 [J/(mol*K)]", "Gas constant");
    model.param().set("T0", "50 [degC]", "Initial temperature");
    model.param().set("Tamb", "T0", "Ambient temperature");
    model.param().set("dwidth", "100 [um]", "Thickness value for boundary conduction");
    model.param().set("epsilon", "0.15*D", "interface thickness");
    model.param().set("D", "100 [um]", "Particle size");
    model.param().set("Nx", "10", "Number of particles in x direction");
    model.param().set("Ny", "2", "Number of particles in the y direction");
    model.param().set("Lheight", "100 [um]", "Height of previous layer");
    model.param().set("vlaser", "500 [mm/s]", "Laser speed");
    model.param().set("Ep", "300 [W]", "Laser power");
    model.param().set("laserpenetration", "D", "Laser penetration depth");
    model.param().set("A", "0.5", "Laser power absorption factor");
    model.param().set("rlaser", "50 [um]", "Radial distance of laser width");
    model.param().set("xr", "onlocation - vlaser * ontime");
    model.param().set("yr", "Lheight + D*Ny");
    model.param().set("onlocation", "2*D + D/2");
    model.param().set("ontime", "0.02 [s]");
    model.param().set("offtime", "0.16 [s]");
    model.param().set("Dp1", "1/(3*kappa)", "Radiation factor assuming no scattering");
    model.param().set("n", "1");
    model.param().set("kappa", "1/pd", "Absorption coefficient over all wavelengths of light");
    model.param().set("sigma_b", "5.67e-8 [W/(m^2*K^4)]", "Stefan Boltzmann constant");
    model.param().set("pd", "80 [um]", "Penetration depth of light");

    model.component().create("comp1", true);

    model.component("comp1").geom().create("geom2", 2);

    model.result().table().create("evl2", "Table");
    model.result().table().create("tbl1", "Table");
    model.result().table().create("tbl2", "Table");
    model.result().table().create("tbl4", "Table");

    model.component("comp1").func().create("rect1", "Rectangle");
    model.component("comp1").func().create("gp1", "GaussianPulse");
    model.component("comp1").func("rect1").label("onoff");
    model.component("comp1").func("rect1").set("funcname", "onoff");
    model.component("comp1").func("rect1").set("lower", "ontime");
    model.component("comp1").func("rect1").set("upper", "offtime");
    model.component("comp1").func("rect1").set("smooth", 0.001);
    model.component("comp1").func("gp1").label("delta");
    model.component("comp1").func("gp1").set("funcname", "delta");
    model.component("comp1").func("gp1").set("sigma", 5);

    model.component("comp1").mesh().create("mesh3");

    model.component("comp1").geom("geom2").label("Geometry 1");
    model.component("comp1").geom("geom2").lengthUnit("\u00b5m");
    model.component("comp1").geom("geom2").create("c1", "Circle");
    model.component("comp1").geom("geom2").feature("c1").set("pos", new String[]{"D/2 + D", "D/2 + Lheight"});
    model.component("comp1").geom("geom2").feature("c1").set("r", "D/2");
    model.component("comp1").geom("geom2").create("arr1", "Array");
    model.component("comp1").geom("geom2").feature("arr1").set("fullsize", new String[]{"Nx", "Ny"});
    model.component("comp1").geom("geom2").feature("arr1").set("displ", new String[]{"D", "D"});
    model.component("comp1").geom("geom2").feature("arr1").selection("input").set("c1");
    model.component("comp1").geom("geom2").create("r1", "Rectangle");
    model.component("comp1").geom("geom2").feature("r1").set("pos", new int[]{0, 0});
    model.component("comp1").geom("geom2").feature("r1").set("size", new String[]{"D*(Nx+2)", "Lheight + D*Ny*2"});
    model.component("comp1").geom("geom2").create("r2", "Rectangle");
    model.component("comp1").geom("geom2").feature("r2").set("pos", new int[]{0, 0});
    model.component("comp1").geom("geom2").feature("r2").set("size", new String[]{"D*(Nx+2)", "Lheight"});
    model.component("comp1").geom("geom2").run();

    model.component("comp1").variable().create("var1");
    model.component("comp1").variable("var1").set("Ed", "Ep / (D*pi*rlaser^2)");
    model.component("comp1").variable("var1")
         .set("G_space", "exp(-1* ( (x  - (xr+vlaser*t) )^2 ) / (2 * rlaser^2 ) )");
    model.component("comp1").variable("var1").set("BeerLambert", "exp( (y - yr) / laserpenetration)");
    model.component("comp1").variable("var1").set("Qlaser", "A*Ed*G_space*onoff(t) * BeerLambert  * pf.Vf2");
    model.component("comp1").variable("var1").set("k", "mat6.def.k(T0)*pf.Vf1 + mat7.def.k(T)*pf.Vf2");
    model.component("comp1").variable("var1")
         .set("rho", "mat6.def.rho(1 [atm], T0)*pf.Vf1 + mat7.def.rho(T0)*pf.Vf2", "Density");
    model.component("comp1").variable("var1")
         .set("Cp", "mat6.def.Cp(T0)*pf.Vf1 + mat7.def.Cp(T)*pf.Vf2", "Heat capacity at constant pressure");
    model.component("comp1").variable("var1")
         .set("mu", "mat6.def.eta(T0)*pf.Vf1 + mat7.def.mu(T)*pf.Vf2 + 0.005", "Dynamic viscosity");
    model.component("comp1").variable("var1").set("sigma", "mat7.def.sigma(T)");
    model.component("comp1").variable("var1").set("lambda", "3*epsilon*sigma/sqrt(8)");
    model.component("comp1").variable("var1").set("Fstx", "lambda/(epsilon^2) * psi * phipfx");
    model.component("comp1").variable("var1").set("Fsty", "lambda/(epsilon^2) * psi * phipfy");
    model.component("comp1").variable("var1").set("Qr", "kappa*(Grad - 4*pi*Ib)");
    model.component("comp1").variable("var1").set("Ib", "n^2*sigma_b*T^4/pi * pf.Vf2");

    model.component("comp1").material().create("mat6", "Common");
    model.component("comp1").material().create("mat7", "Common");
    model.component("comp1").material("mat6").selection().set();
    model.component("comp1").material("mat6").propertyGroup("def").func().create("eta", "Piecewise");
    model.component("comp1").material("mat6").propertyGroup("def").func().create("Cp", "Piecewise");
    model.component("comp1").material("mat6").propertyGroup("def").func().create("rho", "Analytic");
    model.component("comp1").material("mat6").propertyGroup("def").func().create("k", "Piecewise");
    model.component("comp1").material("mat6").propertyGroup().create("idealGas", "Ideal gas");
    model.component("comp1").material("mat6").propertyGroup("idealGas").func().create("Cp", "Piecewise");
    model.component("comp1").material("mat7").propertyGroup("def").func().create("an1", "Analytic");
    model.component("comp1").material("mat7").propertyGroup("def").func().create("an2", "Analytic");
    model.component("comp1").material("mat7").propertyGroup("def").func().create("an4", "Analytic");
    model.component("comp1").material("mat7").propertyGroup("def").func().create("an5", "Analytic");
    model.component("comp1").material("mat7").propertyGroup("def").func().create("pw1", "Piecewise");
    model.component("comp1").material("mat7").propertyGroup().create("Enu", "Young's modulus and Poisson's ratio");
    model.component("comp1").material("mat7").propertyGroup().create("Murnaghan", "Murnaghan");

    model.component("comp1").physics().create("spf", "LaminarFlow", "geom2");
    model.component("comp1").physics("spf").create("vf2", "VolumeForce", 2);
    model.component("comp1").physics("spf").feature("vf2").selection().all();
    model.component("comp1").physics("spf").create("out1", "OutletBoundary", 1);
    model.component("comp1").physics("spf").feature("out1").selection().set(5);
    model.component("comp1").physics().create("pf", "PhaseField", "geom2");
    model.component("comp1").physics("pf").feature("initfluid2").selection()
         .set(1, 12, 13, 15, 16, 18, 19, 21, 22, 24, 25, 27, 28, 30, 31, 33, 34, 36, 37, 39, 40);
    model.component("comp1").physics("pf").create("out1", "Outlet", 1);
    model.component("comp1").physics("pf").feature("out1").selection().set(5);
    model.component("comp1").physics().create("ht", "HeatTransferInFluids", "geom2");
    model.component("comp1").physics("ht").create("hs1", "HeatSource", 2);
    model.component("comp1").physics("ht").feature("hs1").selection()
         .set(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40);
    model.component("comp1").physics("ht").create("hs2", "HeatSource", 2);
    model.component("comp1").physics("ht").feature("hs2").selection()
         .set(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40);
    model.component("comp1").physics("ht").create("hf1", "HeatFluxBoundary", 1);
    model.component("comp1").physics("ht").feature("hf1").selection().all();
    model.component("comp1").physics().create("rteeq", "PoissonEquation", "geom2");
    model.component("comp1").physics("rteeq").identifier("rteeq");
    model.component("comp1").physics("rteeq").field("dimensionless").field("Grad");
    model.component("comp1").physics("rteeq").prop("Units").set("DependentVariableQuantity", "none");
    model.component("comp1").physics("rteeq").prop("Units").set("CustomDependentVariableUnit", "W*m^-2");
    model.component("comp1").physics("rteeq").create("flux1", "FluxBoundary", 1);
    model.component("comp1").physics("rteeq").feature("flux1").selection().set(1, 2, 3, 5, 16, 17);

    model.component("comp1").mesh("mesh3").create("ftri1", "FreeTri");
    model.component("comp1").mesh("mesh3").feature("ftri1").create("size1", "Size");

    model.component("comp1").probe().create("point3", "Point");
    model.component("comp1").probe().create("point4", "Point");
    model.component("comp1").probe("point3").selection()
         .set(4, 6, 7, 9, 11, 12, 14, 16, 17, 19, 21, 22, 24, 26, 27, 29, 31, 32, 34, 36, 37, 39, 41, 42, 44, 46, 47, 49, 51, 52, 54);
    model.component("comp1").probe("point4").selection().set(6, 11, 16, 21, 26, 31, 36, 41, 46, 51);

    model.result().table("evl2").label("Evaluation 2D");
    model.result().table("evl2").comments("Interactive 2D values");
    model.result().table("tbl1").label("Probe Table 1");
    model.result().table("tbl2").label("Probe Table 2");
    model.result().table("tbl4").label("Objective Table 4");

    model.component("comp1").view("view1").axis().set("xmin", -30.000059127807617);
    model.component("comp1").view("view1").axis().set("xmax", 1230);
    model.component("comp1").view("view1").axis().set("ymin", -233.37283325195312);
    model.component("comp1").view("view1").axis().set("ymax", 733.372802734375);

    model.component("comp1").material("mat6").label("Oxygen");
    model.component("comp1").material("mat6").set("family", "air");
    model.component("comp1").material("mat6").propertyGroup("def").func("eta").set("arg", "T");
    model.component("comp1").material("mat6").propertyGroup("def").func("eta")
         .set("pieces", new String[][]{{"150.0", "600.0", "-5.55818182E-7+9.24202797E-8*T^1-8.71841492E-11*T^2+4.82983683E-14*T^3"}});
    model.component("comp1").material("mat6").propertyGroup("def").func("eta").set("argunit", "K");
    model.component("comp1").material("mat6").propertyGroup("def").func("eta").set("fununit", "Pa*s");
    model.component("comp1").material("mat6").propertyGroup("def").func("Cp").set("arg", "T");
    model.component("comp1").material("mat6").propertyGroup("def").func("Cp")
         .set("pieces", new String[][]{{"150.0", "600.0", "959.514545-0.416383077*T^1+7.63158508E-4*T^2+1.46018648E-6*T^3-3.24009324E-9*T^4+1.6E-12*T^5"}});
    model.component("comp1").material("mat6").propertyGroup("def").func("Cp").set("argunit", "K");
    model.component("comp1").material("mat6").propertyGroup("def").func("Cp").set("fununit", "J/(kg*K)");
    model.component("comp1").material("mat6").propertyGroup("def").func("rho")
         .set("expr", "pA*0.032/R_const[K*mol/J]/T");
    model.component("comp1").material("mat6").propertyGroup("def").func("rho").set("args", new String[]{"pA", "T"});
    model.component("comp1").material("mat6").propertyGroup("def").func("rho").set("dermethod", "manual");
    model.component("comp1").material("mat6").propertyGroup("def").func("rho")
         .set("argders", new String[][]{{"pA", "d(pA*0.032/R_const/T,pA)"}, {"T", "d(pA*0.032/R_const/T,T)"}});
    model.component("comp1").material("mat6").propertyGroup("def").func("rho").set("fununit", "kg/m^3");
    model.component("comp1").material("mat6").propertyGroup("def").func("rho")
         .set("argunit", new String[]{"Pa", "K"});
    model.component("comp1").material("mat6").propertyGroup("def").func("rho")
         .set("plotargs", new String[][]{{"pA", "101325", "101325"}, {"T", "273.15", "293.15"}});
    model.component("comp1").material("mat6").propertyGroup("def").func("k").set("arg", "T");
    model.component("comp1").material("mat6").propertyGroup("def").func("k")
         .set("pieces", new String[][]{{"150.0", "600.0", "-0.0070110303+1.688723E-4*T^1-2.28911422E-7*T^2+1.6991453E-10*T^3"}});
    model.component("comp1").material("mat6").propertyGroup("def").func("k").set("argunit", "K");
    model.component("comp1").material("mat6").propertyGroup("def").func("k").set("fununit", "W/(m*K)");
    model.component("comp1").material("mat6").propertyGroup("def").set("dynamicviscosity", "eta(T)");
    model.component("comp1").material("mat6").propertyGroup("def").set("ratioofspecificheat", "1.4");
    model.component("comp1").material("mat6").propertyGroup("def").set("heatcapacity", "Cp(T)");
    model.component("comp1").material("mat6").propertyGroup("def").set("density", "rho(pA,T)");
    model.component("comp1").material("mat6").propertyGroup("def")
         .set("thermalconductivity", new String[]{"k(T)", "0", "0", "0", "k(T)", "0", "0", "0", "k(T)"});
    model.component("comp1").material("mat6").propertyGroup("def").addInput("temperature");
    model.component("comp1").material("mat6").propertyGroup("def").addInput("pressure");
    model.component("comp1").material("mat6").propertyGroup("idealGas").func("Cp").label("Piecewise 2");
    model.component("comp1").material("mat6").propertyGroup("idealGas").func("Cp").set("arg", "T");
    model.component("comp1").material("mat6").propertyGroup("idealGas").func("Cp")
         .set("pieces", new String[][]{{"150.0", "600.0", "959.514545-0.416383077*T^1+7.63158508E-4*T^2+1.46018648E-6*T^3-3.24009324E-9*T^4+1.6E-12*T^5"}});
    model.component("comp1").material("mat6").propertyGroup("idealGas").func("Cp").set("argunit", "K");
    model.component("comp1").material("mat6").propertyGroup("idealGas").func("Cp").set("fununit", "J/(kg*K)");
    model.component("comp1").material("mat6").propertyGroup("idealGas").set("Rs", "R_const/Mn");
    model.component("comp1").material("mat6").propertyGroup("idealGas").set("heatcapacity", "Cp(T)");
    model.component("comp1").material("mat6").propertyGroup("idealGas").set("ratioofspecificheat", "1.4");
    model.component("comp1").material("mat6").propertyGroup("idealGas").set("molarmass", "0.03199");
    model.component("comp1").material("mat6").propertyGroup("idealGas").addInput("temperature");
    model.component("comp1").material("mat6").propertyGroup("idealGas").addInput("pressure");
    model.component("comp1").material("mat7").label("Molybdenum");
    model.component("comp1").material("mat7").set("family", "custom");
    model.component("comp1").material("mat7").set("customspecular", new double[]{0.7843137254901961, 1, 1});
    model.component("comp1").material("mat7")
         .set("customdiffuse", new double[]{0.7843137254901961, 0.7843137254901961, 0.7843137254901961});
    model.component("comp1").material("mat7").set("fresnel", 0.3);
    model.component("comp1").material("mat7").set("roughness", 0.1);
    model.component("comp1").material("mat7").set("metallic", 0);
    model.component("comp1").material("mat7").set("pearl", 0);
    model.component("comp1").material("mat7").set("diffusewrap", 0);
    model.component("comp1").material("mat7").set("clearcoat", 0);
    model.component("comp1").material("mat7").set("reflectance", 0);
    model.component("comp1").material("mat7").propertyGroup("def").func("an1").label("Cp");
    model.component("comp1").material("mat7").propertyGroup("def").func("an1").set("funcname", "Cp");
    model.component("comp1").material("mat7").propertyGroup("def").func("an1")
         .set("expr", "(34.2 + 1.13 *10^(-3) * (T - 2896 [K] ) ) * ( 1 / 95.5 ) + 3.5183e5 * delta( T - 2896 [K] )");
    model.component("comp1").material("mat7").propertyGroup("def").func("an1").set("args", new String[]{"T"});
    model.component("comp1").material("mat7").propertyGroup("def").func("an1").set("fununit", "J/(g*K)");
    model.component("comp1").material("mat7").propertyGroup("def").func("an1").set("argunit", new String[]{"K"});
    model.component("comp1").material("mat7").propertyGroup("def").func("an1")
         .set("plotargs", new String[][]{{"T", "2000", "2550"}});
    model.component("comp1").material("mat7").propertyGroup("def").func("an2").label("sigma");
    model.component("comp1").material("mat7").propertyGroup("def").func("an2").set("funcname", "sigma");
    model.component("comp1").material("mat7").propertyGroup("def").func("an2")
         .set("expr", "(2.29 * 10^3 - 0.26 * (T - 2896 [K]) ) * 0.001");
    model.component("comp1").material("mat7").propertyGroup("def").func("an2").set("args", new String[]{"T"});
    model.component("comp1").material("mat7").propertyGroup("def").func("an2").set("fununit", "N/m");
    model.component("comp1").material("mat7").propertyGroup("def").func("an2").set("argunit", new String[]{"K"});
    model.component("comp1").material("mat7").propertyGroup("def").func("an2")
         .set("plotargs", new String[][]{{"T", "0", "3000"}});
    model.component("comp1").material("mat7").propertyGroup("def").func("an4").label("rho");
    model.component("comp1").material("mat7").propertyGroup("def").func("an4").set("funcname", "rho");
    model.component("comp1").material("mat7").propertyGroup("def").func("an4")
         .set("expr", "10.2*exp(12.270764915751679e-6*(273- T) + 0.008718337229355673e-6*(273^2 - T^2)/2 - 7.799703040318278e-12*(273^3 - T^3)/3 + 5.77704867949318e-15*(273^4 - T^4)/4)");
    model.component("comp1").material("mat7").propertyGroup("def").func("an4").set("args", new String[]{"T"});
    model.component("comp1").material("mat7").propertyGroup("def").func("an4").set("fununit", "g/cm^3");
    model.component("comp1").material("mat7").propertyGroup("def").func("an4").set("argunit", new String[]{"K"});
    model.component("comp1").material("mat7").propertyGroup("def").func("an4")
         .set("plotargs", new String[][]{{"T", "300", "3000"}});
    model.component("comp1").material("mat7").propertyGroup("def").func("an5").label("k");
    model.component("comp1").material("mat7").propertyGroup("def").func("an5").set("funcname", "k");
    model.component("comp1").material("mat7").propertyGroup("def").func("an5").set("expr", "138[W/(m*K)]");
    model.component("comp1").material("mat7").propertyGroup("def").func("an5").set("args", new String[]{"T"});
    model.component("comp1").material("mat7").propertyGroup("def").func("an5").set("fununit", "W/(m*K)");
    model.component("comp1").material("mat7").propertyGroup("def").func("an5").set("argunit", new String[]{"K"});
    model.component("comp1").material("mat7").propertyGroup("def").func("an5")
         .set("plotargs", new String[][]{{"T", "0", "1"}});
    model.component("comp1").material("mat7").propertyGroup("def").func("pw1").label("mu");
    model.component("comp1").material("mat7").propertyGroup("def").func("pw1").set("funcname", "mu");
    model.component("comp1").material("mat7").propertyGroup("def").func("pw1").set("arg", "T");
    model.component("comp1").material("mat7").propertyGroup("def").func("pw1").set("smooth", "contd2");
    model.component("comp1").material("mat7").propertyGroup("def").func("pw1")
         .set("pieces", new String[][]{{"200", "500", "120"}, {"500", "3000", "0.27 * exp(73e3 / (R*T) ) * 0.001"}});
    model.component("comp1").material("mat7").propertyGroup("def").func("pw1").set("argunit", "K");
    model.component("comp1").material("mat7").propertyGroup("def").func("pw1").set("fununit", "Pa*s");
    model.component("comp1").material("mat7").propertyGroup("def")
         .set("thermalexpansioncoefficient", new String[]{"5.3e-5[1/K]", "0", "0", "0", "5.3e-5[1/K]", "0", "0", "0", "5.3e-5[1/K]"});
    model.component("comp1").material("mat7").propertyGroup("def").set("density", "10200[kg/m^3]");
    model.component("comp1").material("mat7").propertyGroup("def").set("heatcapacity", "Cp(T)");
    model.component("comp1").material("mat7").propertyGroup("def")
         .set("thermalconductivity", new String[]{"138[W/(m*K)]", "0", "0", "0", "138[W/(m*K)]", "0", "0", "0", "138[W/(m*K)]"});
    model.component("comp1").material("mat7").propertyGroup("Enu").set("E", "312[GPa]");
    model.component("comp1").material("mat7").propertyGroup("Enu").set("nu", "0.31");
    model.component("comp1").material("mat7").propertyGroup("Murnaghan").set("l", "-300[GPa]");
    model.component("comp1").material("mat7").propertyGroup("Murnaghan").set("m", "-850[GPa]");
    model.component("comp1").material("mat7").propertyGroup("Murnaghan").set("n", "-910[GPa]");

    model.component("comp1").coordSystem("sys1").set("name", "sys2");

    model.component("comp1").physics("spf").prop("ShapeProperty").set("order_fluid", 3);
    model.component("comp1").physics("spf").prop("EquationForm").set("form", "Transient");
    model.component("comp1").physics("spf").prop("PhysicalModelProperty").set("IncludeGravity", true);
    model.component("comp1").physics("spf").prop("PhysicalModelProperty").set("Tref", "T0");
    model.component("comp1").physics("spf").prop("PhysicalModelProperty").set("pref", "1 [atm]");
    model.component("comp1").physics("spf").prop("TurbulenceModelProperty").set("TurbulenceModel", "AlgebraicYplus");
    model.component("comp1").physics("spf").prop("InconsistentStabilization").set("IsotropicDiffusion", true);
    model.component("comp1").physics("spf").prop("AdvancedSettingProperty").set("useBNS", true);
    model.component("comp1").physics("spf").feature("fp1").set("rho_mat", "userdef");
    model.component("comp1").physics("spf").feature("fp1").set("rho", "rho");
    model.component("comp1").physics("spf").feature("fp1").set("minput_temperature_src", "userdef");
    model.component("comp1").physics("spf").feature("fp1").set("minput_temperature", "T");
    model.component("comp1").physics("spf").feature("fp1").set("mu_mat", "userdef");
    model.component("comp1").physics("spf").feature("fp1").set("mu", "mu");
    model.component("comp1").physics("spf").feature("init1").set("p_init", "1[atm]");
    model.component("comp1").physics("spf").feature("init1")
         .set("CompensateForHydrostaticPressureApproximation", false);
    model.component("comp1").physics("spf").feature("init1").set("CompensateForHydrostaticPressure", false);
    model.component("comp1").physics("spf").feature("wallbc1").set("beta_factor", 1);
    model.component("comp1").physics("spf").feature("vf2").set("F", new String[][]{{"Fstx"}, {"Fsty"}, {"0"}});
    model.component("comp1").physics("spf").feature("vf2").label("Surface Tension");
    model.component("comp1").physics("spf").feature("out1").set("p0", "1 [atm]");
    model.component("comp1").physics("spf").feature("out1").set("NormalFlow", true);
    model.component("comp1").physics("spf").feature("out1").set("CompensateForHydrostaticPressure", false);
    model.component("comp1").physics("pf").prop("ShapeProperty").set("order_phasefield", 2);
    model.component("comp1").physics("pf").feature("pfm1").set("epsilon_pf", "epsilon");
    model.component("comp1").physics("pf").feature("pfm1").set("chi", 0);
    model.component("comp1").physics("pf").feature("pfm1").set("U", "1e-12");
    model.component("comp1").physics("pf").feature("pfm1").set("u", new String[][]{{"u"}, {"v"}, {"0"}});
    model.component("comp1").physics("pf").feature("pfm1").set("sigma", "sigma");
    model.component("comp1").physics("pf").feature("ww1").set("surfaceE1", "0.001[J/m^2]");
    model.component("comp1").physics("pf").feature("ww1").set("surfaceE2", "0.001[J/m^2]");
    model.component("comp1").physics("ht").prop("PhysicalModelProperty").set("dz", "dwidth");
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
    model.component("comp1").physics("ht").feature("hs1").set("materialType", "nonSolid");
    model.component("comp1").physics("ht").feature("hs2").set("Q0", "Qr");
    model.component("comp1").physics("ht").feature("hs2").set("materialType", "nonSolid");
    model.component("comp1").physics("ht").feature("hf1").set("q0_input", "k*(Tamb - T) / dwidth");
    model.component("comp1").physics("rteeq").label("RTE");
    model.component("comp1").physics("rteeq").prop("ShapeProperty").set("order", 1);
    model.component("comp1").physics("rteeq").prop("ShapeProperty").set("valueType", "real");
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
    model.component("comp1").mesh("mesh3").feature("ftri1").set("method", "del");
    model.component("comp1").mesh("mesh3").feature("ftri1").feature("size1").set("hauto", 3);
    model.component("comp1").mesh("mesh3").feature("ftri1").feature("size1").set("custom", "on");
    model.component("comp1").mesh("mesh3").feature("ftri1").feature("size1").set("hmax", "0.3*D");
    model.component("comp1").mesh("mesh3").feature("ftri1").feature("size1").set("hmaxactive", true);
    model.component("comp1").mesh("mesh3").feature("ftri1").feature("size1").set("hmin", "0.15*D");
    model.component("comp1").mesh("mesh3").feature("ftri1").feature("size1").set("hminactive", true);
    model.component("comp1").mesh("mesh3").feature("ftri1").feature("size1").set("hcurve", 0.5);
    model.component("comp1").mesh("mesh3").feature("ftri1").feature("size1").set("hcurveactive", true);
    model.component("comp1").mesh("mesh3").feature("ftri1").feature("size1").set("hnarrow", 0.1);
    model.component("comp1").mesh("mesh3").feature("ftri1").feature("size1").set("hnarrowactive", true);
    model.component("comp1").mesh("mesh3").feature("ftri1").feature("size1").set("hgrad", 5);
    model.component("comp1").mesh("mesh3").feature("ftri1").feature("size1").set("hgradactive", false);
    model.component("comp1").mesh("mesh3").run();

    model.component("comp1").probe("point3").label("Density Probe");
    model.component("comp1").probe("point3").set("probename", "rhoprobe");
    model.component("comp1").probe("point3").set("expr", "rho");
    model.component("comp1").probe("point3").set("unit", "kg/m^3");
    model.component("comp1").probe("point3").set("descr", "Density");
    model.component("comp1").probe("point3").set("table", "tbl2");
    model.component("comp1").probe("point3").set("window", "window2");
    model.component("comp1").probe("point4").label("Porosity Probe");
    model.component("comp1").probe("point4").set("probename", "porosity");
    model.component("comp1").probe("point4").set("expr", "pf.Vf2");
    model.component("comp1").probe("point4").set("unit", "1");
    model.component("comp1").probe("point4").set("descr", "Volume fraction of fluid 2");
    model.component("comp1").probe("point4").set("table", "tbl2");
    model.component("comp1").probe("point4").set("window", "window2");

    model.component("comp1").physics("ht").feature("fluid1").set("minput_pressure_src", "root.comp1.spf.pA");

    model.study().create("std1");
    model.study("std1").create("phasei", "PhaseInitialization");
    model.study("std1").create("time", "Transient");
    model.study("std1").feature("phasei")
         .set("activate", new String[]{"spf", "off", "pf", "on", "ht", "on", "rteeq", "off", "frame:spatial1", "on", 
         "frame:material1", "on"});

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

    return model;
  }

  public static Model run2(Model model) {
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
    model.sol("sol4").feature("s1").feature().remove("fcDef");
    model.sol("sol4").feature("t1").create("fc1", "FullyCoupled");
    model.sol("sol4").feature("t1").create("d1", "Direct");
    model.sol("sol4").feature("t1").create("i1", "Iterative");
    model.sol("sol4").feature("t1").create("i2", "Iterative");
    model.sol("sol4").feature("t1").create("i3", "Iterative");
    model.sol("sol4").feature("t1").create("tp1", "TimeParametric");
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
    model.sol().create("sol6");
    model.sol("sol6").study("std1");
    model.sol("sol6").label("Parametric Solutions 1");
    model.sol().create("sol19");
    model.sol("sol19").study("std1");
    model.sol("sol19").label("Parametric Solutions 2");
    model.sol().create("sol33");
    model.sol("sol33").study("std1");
    model.sol("sol33").label("Parametric Solutions 3");
    model.sol().create("sol35");
    model.sol("sol35").study("std1");
    model.sol("sol35").label("Parametric Solutions 4");

    model.result().dataset().create("an1_ds1", "Grid1D");
    model.result().dataset().create("avh4", "Average");
    model.result().dataset().create("avh5", "Average");
    model.result().dataset().create("dset9", "Solution");
    model.result().dataset().create("dset10", "Solution");
    model.result().dataset().create("pw1_ds1", "Grid1D");
    model.result().dataset("dset3").set("probetag", "point4");
    model.result().dataset("dset4").set("solution", "sol3");
    model.result().dataset("an1_ds1").set("data", "none");
    model.result().dataset("dset5").set("solution", "sol4");
    model.result().dataset("dset6").set("solution", "sol5");
    model.result().dataset("dset7").set("solution", "sol6");
    model.result().dataset("avh4").set("probetag", "point3");
    model.result().dataset("avh4").set("data", "dset3");
    model.result().dataset("avh4").selection().geom("geom2", 0);
    model.result().dataset("avh4").selection()
         .set(4, 6, 7, 9, 11, 12, 14, 16, 17, 19, 21, 22, 24, 26, 27, 29, 31, 32, 34, 36, 37, 39, 41, 42, 44, 46, 47, 49, 51, 52, 54);
    model.result().dataset("avh5").set("probetag", "point4");
    model.result().dataset("avh5").set("data", "dset3");
    model.result().dataset("avh5").selection().geom("geom2", 0);
    model.result().dataset("avh5").selection().set(6, 11, 16, 21, 26, 31, 36, 41, 46, 51);
    model.result().dataset("dset8").set("solution", "sol19");
    model.result().dataset("dset9").set("solution", "sol33");
    model.result().dataset("dset10").set("solution", "sol35");
    model.result().dataset("pw1_ds1").set("data", "none");
    model.result().numerical().create("pev4", "EvalPoint");
    model.result().numerical().create("pev5", "EvalPoint");
    model.result().numerical("pev4").set("probetag", "point3");
    model.result().numerical("pev5").set("probetag", "point4");
    model.result().create("pg1", "PlotGroup2D");
    model.result().create("pg10", "PlotGroup2D");
    model.result().create("pg2", "PlotGroup2D");
    model.result().create("pg3", "PlotGroup2D");
    model.result().create("pg4", "PlotGroup2D");
    model.result().create("pg7", "PlotGroup1D");
    model.result().create("pg8", "PlotGroup2D");
    model.result().create("pg11", "PlotGroup1D");
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
    model.result("pg3").set("data", "dset10");
    model.result("pg3").create("surf1", "Surface");
    model.result("pg3").create("con1", "Contour");
    model.result("pg3").create("arwl1", "ArrowLine");
    model.result("pg3").feature("con1").set("expr", "pf.Vf2");
    model.result("pg4").set("data", "dset5");
    model.result("pg4").create("surf1", "Surface");
    model.result("pg4").create("con1", "Contour");
    model.result("pg4").feature("surf1").set("expr", "mu");
    model.result("pg4").feature("con1").set("expr", "pf.Vf2");
    model.result("pg7").set("probetag", "window2_default");
    model.result("pg7").create("tblp3", "Table");
    model.result("pg7").create("tblp4", "Table");
    model.result("pg7").feature("tblp3").set("probetag", "dom1,point3,point4");
    model.result("pg7").feature("tblp4").set("probetag", "point1,point2");
    model.result("pg8").set("data", "dset5");
    model.result("pg8").create("surf1", "Surface");
    model.result("pg8").feature("surf1").set("expr", "rho");
    model.result("pg11").create("plot1", "Function");
    model.result("pg11").feature("plot1").set("expr", "comp1.mat7.def.mu(T)");
    model.result().export().create("anim1", "Animation");
    model.result().export().create("anim2", "Animation");
    model.result().export().create("anim3", "Animation");
    model.result().export().create("data1", "Data");

    model.component("comp1").probe("point3").genResult(null);
    model.component("comp1").probe("point4").genResult(null);

    model.study("std1").feature("time").set("tlist", "range(0,0.0001,0.04)");
    model.study("std1").feature("time").set("usertol", true);
    model.study("std1").feature("time").set("rtol", 0.1);
    model.study("std1").feature("time").set("plot", true);
    model.study("std1").feature("time").set("plotgroup", "pg10");
    model.study("std1").feature("time").set("plotfreq", "tsteps");
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
    model.sol("sol1").feature("v2").set("control", "user");
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
    model.sol("sol1").feature("t1").set("tlist", "range(0,0.001,1)");
    model.sol("sol1").feature("t1").set("rtol", 0.005);
    model.sol("sol1").feature("t1").set("atolglobalfactor", 0.05);
    model.sol("sol1").feature("t1")
         .set("atolmethod", new String[]{"comp1_GI", "global", "comp1_p", "scaled", "comp1_phipf", "global", "comp1_psi", "global", "comp1_T", "global", 
         "comp1_u", "global", "comp1_Grad", "global"});
    model.sol("sol1").feature("t1")
         .set("atolfactor", new String[]{"comp1_GI", "0.1", "comp1_p", "1", "comp1_phipf", "0.1", "comp1_psi", "0.1", "comp1_T", "0.1", 
         "comp1_u", "0.1", "comp1_Grad", "0.1"});
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
    model.sol("sol3").label("Refined Mesh Solution 1");
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

    return model;
  }

  public static Model run3(Model model) {
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
    model.sol("sol4").feature("su1").label("Solution Store 1");
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
    model.sol("sol4").feature("v2").set("clist", new String[]{"range(0,0.0001,0.04)", "4.0E-5[s]"});
    model.sol("sol4").feature("v2").feature("comp1_phipf").set("scalemethod", "manual");
    model.sol("sol4").feature("v2").feature("comp1_phipf").set("scaleval", 1);
    model.sol("sol4").feature("v2").feature("comp1_T").set("scalemethod", "manual");
    model.sol("sol4").feature("v2").feature("comp1_T").set("scaleval", 10);
    model.sol("sol4").feature("t1").label("Time-Dependent Solver 1.1");
    model.sol("sol4").feature("t1").set("control", "time");
    model.sol("sol4").feature("t1").set("tlist", "range(0,0.0001,0.04)");
    model.sol("sol4").feature("t1").set("rtol", 0.1);
    model.sol("sol4").feature("t1").set("atolglobalfactor", 1);
    model.sol("sol4").feature("t1")
         .set("atolmethod", new String[]{"comp1_GI", "global", "comp1_p", "scaled", "comp1_phipf", "global", "comp1_psi", "global", "comp1_T", "global", 
         "comp1_u", "global", "comp1_Grad", "global"});
    model.sol("sol4").feature("t1")
         .set("atolfactor", new String[]{"comp1_GI", "0.1", "comp1_p", "1", "comp1_phipf", "0.1", "comp1_psi", "0.1", "comp1_T", "0.1", 
         "comp1_u", "0.1", "comp1_Grad", "0.1"});
    model.sol("sol4").feature("t1").set("maxstepconstraintbdf", "const");
    model.sol("sol4").feature("t1").set("maxstepbdf", 0.001);
    model.sol("sol4").feature("t1").set("maxorder", 1);
    model.sol("sol4").feature("t1").set("stabcntrl", true);
    model.sol("sol4").feature("t1").set("bwinitstepfrac", 0.01);
    model.sol("sol4").feature("t1").set("estrat", "exclude");
    model.sol("sol4").feature("t1").set("plot", true);
    model.sol("sol4").feature("t1").set("plotgroup", "pg10");
    model.sol("sol4").feature("t1").set("plotfreq", "tsteps");
    model.sol("sol4").feature("t1").feature("dDef").label("Direct 2");
    model.sol("sol4").feature("t1").feature("aDef").label("Advanced 1");
    model.sol("sol4").feature("t1").feature("aDef").set("cachepattern", true);
    model.sol("sol4").feature("t1").feature("fc1").label("Fully Coupled 1.1");
    model.sol("sol4").feature("t1").feature("fc1").set("linsolver", "d1");
    model.sol("sol4").feature("t1").feature("fc1").set("dtech", "auto");
    model.sol("sol4").feature("t1").feature("fc1").set("maxiter", 15);
    model.sol("sol4").feature("t1").feature("fc1").set("ntolfact", 0.05);
    model.sol("sol4").feature("t1").feature("fc1").set("termonres", "both");
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
    model.sol("sol4").feature("t1").feature("tp1").active(false);
    model.sol("sol4").feature("t1").feature("tp1").label("Time Parametric 1.1");
    model.sol("sol4").feature("t1").feature("tp1").set("control", "time");
    model.sol("sol4").runAll();

    model.result().dataset("dset3").label("Probe Solution 3");
    model.result().dataset("an1_ds1").label("Grid 1D 1a 1");
    model.result().dataset("an1_ds1").set("function", "all");
    model.result().dataset("an1_ds1").set("par1", "T");
    model.result().dataset("an1_ds1").set("parmin1", 300);
    model.result().dataset("an1_ds1").set("parmax1", 400);
    model.result().dataset("an1_ds1").set("res1", 10000);
    model.result().dataset("an1_ds1").set("distribution", "mixed");
    model.result().dataset("pw1_ds1").set("functionlist", "material/mat7/def");
    model.result().dataset("pw1_ds1").label("Grid 1D 1");
    model.result().dataset("pw1_ds1").set("function", "pw1");
    model.result().dataset("pw1_ds1").set("par1", "T");
    model.result().dataset("pw1_ds1").set("parmin1", -80);
    model.result().dataset("pw1_ds1").set("parmax1", 3280);
    model.result().dataset("pw1_ds1").set("res1", 10000);
    model.result().dataset("pw1_ds1").set("distribution", "mixed");
    model.result().numerical("pev4").setResult();
    model.result("pg1").label("Heating and Volume Fraction");
    model.result("pg1").set("looplevel", new int[]{218});
    model.result("pg1").set("titletype", "label");
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
    model.result("pg10").set("looplevel", new int[]{394});
    model.result("pg10").set("titletype", "label");
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
    model.result("pg2").set("looplevel", new int[]{1});
    model.result("pg2").set("titletype", "label");
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
    model.result("pg3").set("titletype", "label");
    model.result("pg3").feature("surf1").set("colortable", "Magma");
    model.result("pg3").feature("surf1").set("resolution", "normal");
    model.result("pg3").feature("con1").set("coloring", "uniform");
    model.result("pg3").feature("con1").set("color", "green");
    model.result("pg3").feature("con1").set("resolution", "normal");
    model.result("pg3").feature("arwl1").active(false);
    model.result("pg3").feature("arwl1").set("scale", 4.224745667893173E7);
    model.result("pg3").feature("arwl1").set("scaleactive", false);
    model.result("pg4").label("Viscosity");
    model.result("pg4").set("looplevel", new int[]{218});
    model.result("pg4").set("titletype", "manual");
    model.result("pg4").set("title", "Fluid viscosity (Pa*s) with Volume Fraction contours");
    model.result("pg4").set("paramindicator", "t=eval(t,s) (s)");
    model.result("pg4").set("showlegendsmaxmin", true);
    model.result("pg4").feature("surf1").set("descr", "mu");
    model.result("pg4").feature("surf1").set("coloring", "gradient");
    model.result("pg4").feature("surf1").set("topcolor", "magenta");
    model.result("pg4").feature("surf1").set("bottomcolor", "white");
    model.result("pg4").feature("surf1").set("resolution", "normal");
    model.result("pg4").feature("con1").set("resolution", "normal");
    model.result("pg7").label("Density Probe Measurement");
    model.result("pg7").set("xlabel", "t (s)");
    model.result("pg7").set("xlabelactive", true);
    model.result("pg7").set("ylabel", "Density");
    model.result("pg7").set("ylabelactive", true);
    model.result("pg7").set("showlegends", false);
    model.result("pg7").set("windowtitle", "Probe Plot 2");
    model.result("pg7").feature("tblp3").label("Probe Table Graph 3");
    model.result("pg7").feature("tblp3").set("table", "tbl2");
    model.result("pg7").feature("tblp3").set("plotcolumninput", "manual");
    model.result("pg7").feature("tblp3").set("legend", true);
    model.result("pg7").feature("tblp4").label("Probe Table Graph 4");
    model.result("pg7").feature("tblp4").set("table", "tbl1");
    model.result("pg7").feature("tblp4").set("plotcolumninput", "manual");
    model.result("pg7").feature("tblp4").set("legend", true);
    model.result("pg7").feature().remove("tblp5");
    model.result("pg8").label("Density");
    model.result("pg8").set("looplevel", new int[]{218});
    model.result("pg8").set("showlegendsmaxmin", true);
    model.result("pg8").feature("surf1").set("colortable", "Twilight");
    model.result("pg8").feature("surf1").set("resolution", "normal");
    model.result("pg11").label("Viscosity Plot");
    model.result("pg11").set("data", "pw1_ds1");
    model.result("pg11").set("solrepresentation", "solnum");
    model.result("pg11").set("titletype", "manual");
    model.result("pg11").set("title", "Temperature Dependent Viscosity of Molybdenum");
    model.result("pg11").set("xlabel", "Temperature (K)");
    model.result("pg11").set("xlabelactive", true);
    model.result("pg11").set("ylabel", "Viscosity (Pa s)");
    model.result("pg11").set("ylabelactive", true);
    model.result("pg11").feature("plot1").set("solrepresentation", "solnum");
    model.result("pg11").feature("plot1").set("unit", "");
    model.result("pg11").feature("plot1").set("descr", "mu(T)");
    model.result("pg11").feature("plot1").set("xdataexpr", "T");
    model.result("pg11").feature("plot1").set("xdataunit", "m");
    model.result("pg11").feature("plot1").set("xdatadescractive", true);
    model.result("pg11").feature("plot1").set("xdatadescr", "");
    model.result("pg11").feature("plot1").set("lowerbound", 200);
    model.result("pg11").feature("plot1").set("upperbound", 3000);
    model.result("pg11").feature("plot1").set("extrapolation", "leftright");
    model.result("pg11").feature("plot1").set("linewidth", "preference");
    model.result().export("anim1").label("Temperature and Vf2");
    model.result().export("anim1").set("plotgroup", "pg10");
    model.result().export("anim1").set("giffilename", "/users/m/s/mshourya/SLS/Complete/Results/T+Vf2.gif");
    model.result().export("anim1").set("fps", 50);
    model.result().export("anim1").set("framesel", "all");
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
    model.result().export("anim2").label("Pressure and Vf2");
    model.result().export("anim2").set("plotgroup", "pg2");
    model.result().export("anim2").set("giffilename", "/users/m/s/mshourya/SLS/Complete/Results/p+Vf2.gif");
    model.result().export("anim2").set("fps", 50);
    model.result().export("anim2").set("framesel", "all");
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
    model.result().export("anim3").label("Heating and Vf2");
    model.result().export("anim3").set("giffilename", "/users/m/s/mshourya/SLS/Complete/Results/Heating+Vf2.gif");
    model.result().export("anim3").set("fps", 50);
    model.result().export("anim3").set("framesel", "all");
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
    model.result().export("data1").label("Average Density");
    model.result().export("data1").set("data", "avh4");
    model.result().export("data1").set("looplevelinput", new String[]{"all"});
    model.result().export("data1").set("filename", "/users/m/s/mshourya/SLS/Complete/Results/densityprobe.csv");
    model.result().export("data1").set("smooth", "internal");
    model.result().export("data1").set("separator", ",");

    return model;
  }

  public static void main(String[] args) {
    Model model = run();
    model = run2(model);
    run3(model);
  }

}
