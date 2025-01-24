/*
 * ThermalNewtonianIncompressible2D.java
 */

import java.io.*;
import java.util.*;
import java.util.regex.*;
import com.comsol.model.*;
import com.comsol.model.util.*;

/** Model exported on Jan 23 2025, 19:50 by COMSOL 6.1.0.357. */
public class ThermalNewtonianIncompressible2D {
  private static int count = 1;
  public static void drawCircle(GeomSequence geom, List<Double> pos, Double r) {
    String name = "c" + count;
    count = count + 1;

    double[] x = {pos.get(0), pos.get(1)};
    geom.create(name, "Circle");
    geom.feature(name).set("pos", x);
    geom.feature(name).set("r", r);
  }

  public static void drawCircles(GeomSequence geom) {
	List<List<Double>> xList = Arrays.asList(
            Arrays.asList(342.60509461 / 100.0, 368.57425392 / 100.0),
            Arrays.asList(527.68322253 / 100.0, 447.69537508 / 100.0),
            Arrays.asList(401.98493574 / 100.0, 486.1268805 / 100.0),
            Arrays.asList(326.87080771 / 100.0, 425.80642869 / 100.0),
            Arrays.asList(621.27828372 / 100.0, 472.97660615 / 100.0),
            Arrays.asList(168.87637113 / 100.0, 413.05250695 / 100.0),
            Arrays.asList(215.16501926 / 100.0, 455.68443012 / 100.0),
            Arrays.asList(358.55500712 / 100.0, 466.7750568 / 100.0),
            Arrays.asList(151.32939835 / 100.0, 327.50969774 / 100.0),
            Arrays.asList(516.40931341 / 100.0, 377.50483994 / 100.0),
            Arrays.asList( 63.55115054 / 100.0, 478.12699404 / 100.0),
            Arrays.asList(499.43753436 / 100.0, 479.86696999 / 100.0),
            Arrays.asList(427.27502307 / 100.0, 351.77238457 / 100.0),
            Arrays.asList( 93.14295786 / 100.0, 438.87094157 / 100.0),
            Arrays.asList(109.26704982 / 100.0, 372.58293678 / 100.0),
            Arrays.asList(653.01602005 / 100.0, 414.22848953 / 100.0),
            Arrays.asList(972.81756493 / 100.0, 376.91207467 / 100.0),
            Arrays.asList(848.28136044 / 100.0, 386.93775474 / 100.0),
            Arrays.asList(123.83211792 / 100.0, 477.73049012 / 100.0),
            Arrays.asList(448.47145171 / 100.0, 467.69591183 / 100.0),
            Arrays.asList(455.05757232 / 100.0, 409.88400151 / 100.0),
            Arrays.asList(591.84927683 / 100.0, 412.68155169 / 100.0),
            Arrays.asList(813.8106275  / 100.0, 421.37202391 / 100.0),
            Arrays.asList(906.36321242 / 100.0, 383.1462884  / 100.0),
            Arrays.asList(937.32827635 / 100.0, 314.72749465 / 100.0),
            Arrays.asList(807.34117127 / 100.0, 332.79567952 / 100.0),
            Arrays.asList(857.84506074 / 100.0, 456.89392998 / 100.0),
            Arrays.asList(733.57431107 / 100.0, 344.85030257 / 100.0),
            Arrays.asList(159.22051177 / 100.0, 485.93852568 / 100.0),
            Arrays.asList(782.69231581 / 100.0, 455.50958479 / 100.0),
            Arrays.asList(271.1338555  / 100.0, 321.8032231  / 100.0),
            Arrays.asList(397.81446251 / 100.0, 415.29845396 / 100.0),
            Arrays.asList( 43.87068476 / 100.0, 387.42025727 / 100.0),
            Arrays.asList(701.44438514 / 100.0, 451.01528442 / 100.0),
            Arrays.asList(760.76098973 / 100.0, 391.70769196 / 100.0),
            Arrays.asList(236.07764836 / 100.0, 479.02014185 / 100.0),
            Arrays.asList( 84.21812426 / 100.0, 327.05129047 / 100.0),
            Arrays.asList(814.38053898 / 100.0, 489.08861811 / 100.0),
            Arrays.asList(673.77379004 / 100.0, 363.119294   / 100.0),
            Arrays.asList(254.52803474 / 100.0, 401.96797293 / 100.0),
            Arrays.asList(950.36045561 / 100.0, 450.36045561 / 100.0),
            Arrays.asList(247.28347438 / 100.0, 449.45866109 / 100.0),
            Arrays.asList(193.81776867 / 100.0, 478.72762112 / 100.0),
            Arrays.asList(403.99073978 / 100.0, 459.37293395 / 100.0),
            Arrays.asList(640.18056501 / 100.0, 436.6869517  / 100.0),
            Arrays.asList(591.15576067 / 100.0, 324.99652179 / 100.0),
            Arrays.asList( 20.89622292 / 100.0, 479.03548992 / 100.0),
            Arrays.asList(290.31964477 / 100.0, 464.93684542 / 100.0),
            Arrays.asList(489.77342951 / 100.0, 443.61399967 / 100.0),
            Arrays.asList(568.35088881 / 100.0, 474.10426229 / 100.0)
        );

        // Hardcoded rList (radii), each divided by 100
	List<Double> rList = Arrays.asList(
    0.4080955327, 0.2266256335, 0.1383163727,
    0.1854029626, 0.2702339385, 0.4893113965,
    0.1019478032, 0.332249432, 0.3833570294,
    0.4349262408, 0.2177574285, 0.2013303001,
    0.3846742863, 0.2737788862, 0.2325520194,
    0.1189864944, 0.2717765618, 0.2751194836,
    0.2226950988, 0.3230408817, 0.2594091681,
    0.3986006238, 0.1350373183, 0.3071011124,
    0.4439083312, 0.4032255269, 0.4310607002,
    0.1936826243, 0.1406147432, 0.3228791386,
    0.4454769059, 0.3154658211, 0.4387068476,
    0.4898471558, 0.3484577145, 0.2097985815,
    0.2866764212, 0.1088092941, 0.4313462958,
    0.3736966282, 0.4963954439, 0.1066815949,
    0.2127237888, 0.1291849708, 0.1392516616,
    0.4776305276, 0.2086811149, 0.3505173665,
    0.155021577, 0.2589573771
);
    	for (int i = 0; i < xList.size(); i++) {
	    drawCircle(geom, xList.get(i), rList.get(i));
    	}
  }
  public static Model run() {
    Model model = ModelUtil.create("Model");

    model.modelPath("/nas/longleaf/home/mshourya/workspace/SLS/Particle");

    model.label("ThermalNewtonianIncompressibleODE.mph");

    model.param().set("MaterialQuantities", "0", "___________________________________________________");
    model.param().set("rhomat", "1000 [kg/m^3]", "Material density");
    model.param().set("rhoair", "1.204 [kg/m^3]", "Air density");
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
    model.param().set("Q", "0 [W/m^3]");
    model.param().set("deltat", "1 [K]");
    model.param().set("T0", "300[K]");
    model.param().set("Tm", "320 [K]");
    model.param().set("epsilon", "1.5e-6", "interface thickness");
    model.param().set("g", "9.81", "Gravitational acceleration");
    model.param().set("sigma", "0.05", "Surface Tension");

    model.component().create("comp1", true);

    model.component("comp1").geom().create("geom2", 2);

    model.result().table().create("evl2", "Table");
    model.result().table().create("tbl1", "Table");

    model.component("comp1").func().create("step1", "Step");
    model.component("comp1").func().create("an1", "Analytic");
    model.component("comp1").func().create("gp1", "GaussianPulse");
    model.component("comp1").func().create("step2", "Step");
    model.component("comp1").func().create("step3", "Step");
    model.component("comp1").func("step1").label("Viscosity");
    model.component("comp1").func("step1").set("funcname", "mumat");
    model.component("comp1").func("step1").set("location", "Tm");
    model.component("comp1").func("step1").set("from", "mumatsolid");
    model.component("comp1").func("step1").set("to", "mumatliq");
    model.component("comp1").func("step1").set("smooth", "deltat");
    model.component("comp1").func("an1").label("Specific Heat");
    model.component("comp1").func("an1").set("funcname", "cpmat");
    model.component("comp1").func("an1").set("expr", "jumpfunc(T) + stepfunc(T)");
    model.component("comp1").func("an1").set("args", new String[]{"T"});
    model.component("comp1").func("an1").set("argunit", new String[]{""});
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

    model.component("comp1").mesh().create("mesh2");

    model.component("comp1").geom("geom2").label("Geometry 1");
    model.component("comp1").geom("geom2").lengthUnit("\u00b5m");
    model.component("comp1").geom("geom2").create("sq1", "Square");
    model.component("comp1").geom("geom2").feature("sq1").set("pos", new int[]{0, 0});
    model.component("comp1").geom("geom2").feature("sq1").set("base", "corner");
    model.component("comp1").geom("geom2").feature("sq1").set("size", 100);
    drawCircles(model.component("comp1").geom("geom2"));
    model.component("comp1").geom("geom2").run();

    model.component("comp1").variable().create("var1");
    model.component("comp1").variable("var1").set("chi", "pf.hmax*Uscale/(3*sqrt(2)*sigma*epsilon)");
    model.component("comp1").variable("var1").set("lambda", "3*epsilon*sigma/(sqrt(8))");
    model.component("comp1").variable("var1").set("gamma", "chi/(epsilon^2)");

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
    model.component("comp1").physics("spf").create("vf1", "VolumeForce", 2);
    model.component("comp1").physics("spf").feature("vf1").selection().set(1, 2);
    model.component("comp1").physics("spf").create("constr1", "PointwiseConstraint", 0);
    model.component("comp1").physics("spf").feature("constr1").selection().set(2, 8);
    model.component("comp1").physics("spf").create("disc1", "Discretization", -1);
    model.component("comp1").physics("spf").create("out1", "OutletBoundary", 1);
    model.component("comp1").physics("spf").feature("out1").selection().set(3);
    model.component("comp1").physics().create("pf", "PhaseField", "geom2");
    model.component("comp1").physics("pf").feature("initfluid2").selection().set(2);
    model.component("comp1").physics("pf").create("inl1", "InletBoundary", 1);
    model.component("comp1").physics("pf").feature("inl1").selection().set(1, 5);
    model.component("comp1").physics("pf").create("constr1", "PointwiseConstraint", 2);
    model.component("comp1").physics("pf").feature("constr1").selection().set(1, 2);
    model.component("comp1").physics("pf").create("constr2", "PointwiseConstraint", 1);
    model.component("comp1").physics("pf").feature("constr2").selection().set(2, 4);
    model.component("comp1").physics("pf").create("disc1", "Discretization", -1);
    model.component("comp1").physics().create("viscode", "DomainODE", "geom2");
    model.component("comp1").physics("viscode").identifier("viscode");
    model.component("comp1").physics("viscode").field("dimensionless").field("muode");
    model.component("comp1").physics("viscode").field("dimensionless").component(new String[]{"mu2"});
    model.component("comp1").physics("viscode").prop("Units").set("DependentVariableQuantity", "none");
    model.component("comp1").physics("viscode").prop("Units").set("CustomDependentVariableUnit", "Pa*s");
    model.component("comp1").physics("viscode").create("aleq1", "AlgebraicEquation", 2);
    model.component("comp1").physics("viscode").feature("aleq1").selection().set(1, 2);
    model.component("comp1").physics().create("ht", "HeatTransferInFluids", "geom2");
    model.component("comp1").physics("ht").create("ofl1", "ConvectiveOutflow", 1);
    model.component("comp1").physics("ht").feature("ofl1").selection().set(3);
    model.component("comp1").physics("ht").create("hs1", "HeatSource", 2);
    model.component("comp1").physics("ht").feature("hs1").selection().set(2);
    model.component("comp1").physics().create("cpode", "DomainODE", "geom2");
    model.component("comp1").physics("cpode").identifier("cpode");
    model.component("comp1").physics("cpode").field("dimensionless").field("cpfield");
    model.component("comp1").physics("cpode").field("dimensionless").component(new String[]{"Cp2"});
    model.component("comp1").physics("cpode").prop("Units").set("DependentVariableQuantity", "none");
    model.component("comp1").physics("cpode").prop("Units").set("CustomDependentVariableUnit", "J/(kg*K)");
    model.component("comp1").physics("cpode").create("aleq1", "AlgebraicEquation", 2);
    model.component("comp1").physics("cpode").feature("aleq1").selection().set(1, 2);

    model.component("comp1").multiphysics().create("tpf1", "TwoPhaseFlowPhaseField", 2);
    model.component("comp1").multiphysics("tpf1").selection().all();

    model.component("comp1").mesh("mesh2").create("ftri2", "FreeTri");
    model.component("comp1").mesh("mesh2").feature("ftri2").create("dis1", "Distribution");
    model.component("comp1").mesh("mesh2").feature("ftri2").create("size1", "Size");
    model.component("comp1").mesh("mesh2").feature("ftri2").feature("dis1").selection().set(6, 7, 8, 9);

    model.component("comp1").probe().create("point1", "Point");
    model.component("comp1").probe("point1").selection().set();

    model.result().table("evl2").label("Evaluation 2D");
    model.result().table("evl2").comments("Interactive 2D values");
    model.result().table("tbl1").label("Probe Table 1");

    model.component("comp1").view("view1").axis().set("xmin", -65.11363983154297);
    model.component("comp1").view("view1").axis().set("xmax", 65.11363983154297);
    model.component("comp1").view("view1").axis().set("ymin", -22.5);
    model.component("comp1").view("view1").axis().set("ymax", 82.5);

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

    model.component("comp1").physics("spf").prop("PhysicalModelProperty").set("UseReducedPressure", true);
    model.component("comp1").physics("spf").prop("TurbulenceModelProperty").set("TurbulenceModel", "AlgebraicYplus");
    model.component("comp1").physics("spf").prop("InconsistentStabilization").set("IsotropicDiffusion", true);
    model.component("comp1").physics("spf").feature("init1").set("p_init", "1[atm]");
    model.component("comp1").physics("spf").feature("init1").set("CompensateForHydrostaticPressure", false);
    model.component("comp1").physics("spf").feature("vf1")
         .set("F", new String[][]{{"0"}, {"-g_const*(1+phipf)/2 * rhomat"}, {"0"}});
    model.component("comp1").physics("spf").feature("constr1").set("constraintExpression", "p - (1 [atm])");
    model.component("comp1").physics("spf").feature("constr1").set("constraintMethod", "nodal");
    model.component("comp1").physics("spf").feature("disc1").set("order_fluid", 3);
    model.component("comp1").physics("spf").feature("out1").set("p0", "1 [atm]");
    model.component("comp1").physics("pf").feature("pfm1").set("epsilon_pf", "epsilon");
    model.component("comp1").physics("pf").feature("pfm1").set("chiOption", "velocity");
    model.component("comp1").physics("pf").feature("pfm1").set("chi", "chi");
    model.component("comp1").physics("pf").feature("pfm1").set("U", "sqrt(u^2+v^2)");
    model.component("comp1").physics("pf").feature("ww1").set("surfaceE1", "0.001[J/m^2]");
    model.component("comp1").physics("pf").feature("ww1").set("surfaceE2", "0.001[J/m^2]");
    model.component("comp1").physics("pf").feature("inl1").active(false);
    model.component("comp1").physics("pf").feature("constr1").set("constraintExpression", "psixx+psiyy");
    model.component("comp1").physics("pf").feature("constr1").active(false);
    model.component("comp1").physics("pf").feature("constr2").set("constraintExpression", "unx*phix+uny*phiy");
    model.component("comp1").physics("pf").feature("constr2").active(false);
    model.component("comp1").physics("pf").feature("disc1").set("order_phasefield", 3);
    model.component("comp1").physics("viscode").label("Viscosity");
    model.component("comp1").physics("viscode").prop("EquationForm").set("form", "Automatic");
    model.component("comp1").physics("viscode").prop("Units").set("CustomSourceTermUnit", "Pa");
    model.component("comp1").physics("viscode").feature("init1").set("mu2", "mumat(T0)");
    model.component("comp1").physics("viscode").feature("init1").set("mu2t", "d(mumat(T0), T0)*d(T, t)");
    model.component("comp1").physics("viscode").feature("aleq1").set("f", "mu2 - mumat(T)");
    model.component("comp1").physics("ht").prop("InconsistentStabilization").set("HeatIsotropicDiffusion", true);
    model.component("comp1").physics("ht").feature("fluid1").set("Cp_mat", "userdef");
    model.component("comp1").physics("ht").feature("fluid1").set("Cp", "cpair*(1-phipf)/2 + Cp2*(1+phipf)/2");
    model.component("comp1").physics("ht").feature("fluid1").set("rho_mat", "userdef");
    model.component("comp1").physics("ht").feature("fluid1").set("rho", "rhoair*(1-phipf)/2 + rhomat*(1+phipf)/2");
    model.component("comp1").physics("ht").feature("fluid1").set("k_mat", "userdef");
    model.component("comp1").physics("ht").feature("fluid1")
         .set("k", new String[][]{{"kair*(1-phipf)/2 + kmat*(1+phipf)/2"}, {"0"}, {"0"}, {"0"}, {"kair*(1-phipf)/2 + kmat*(1+phipf)/2"}, {"0"}, {"0"}, {"0"}, {"kair*(1-phipf)/2 + kmat*(1+phipf)/2"}});
    model.component("comp1").physics("ht").feature("fluid1").set("u", new String[][]{{"u"}, {"v"}, {"0"}});
    model.component("comp1").physics("ht").feature("init1").set("Tinit", "T0");
    model.component("comp1").physics("ht").feature("hs1").set("Q0", "Q");
    model.component("comp1").physics("cpode").label("Specific Heat");
    model.component("comp1").physics("cpode").prop("EquationForm").set("form", "Automatic");
    model.component("comp1").physics("cpode").feature("init1").set("Cp2", "cpmat(T0)");
    model.component("comp1").physics("cpode").feature("init1").set("Cp2t", "d(cpmat(T0), T0)*d(T, t)");
    model.component("comp1").physics("cpode").feature("aleq1").set("f", "Cp2 - cpmat(T)");

    model.component("comp1").multiphysics("tpf1").set("Fluid1", "mat1");
    model.component("comp1").multiphysics("tpf1").set("rho1_mat", "userdef");
    model.component("comp1").multiphysics("tpf1").set("rho1", "rhoair");
    model.component("comp1").multiphysics("tpf1").set("mu1_mat", "userdef");
    model.component("comp1").multiphysics("tpf1").set("mu1", "muair");
    model.component("comp1").multiphysics("tpf1").set("rho2_mat", "userdef");
    model.component("comp1").multiphysics("tpf1").set("rho2", "rhomat");
    model.component("comp1").multiphysics("tpf1").set("mu2_mat", "userdef");
    model.component("comp1").multiphysics("tpf1").set("mu2", "mu2");
    model.component("comp1").multiphysics("tpf1").set("IncludeSurfaceTensionGradientEffect", true);
    model.component("comp1").multiphysics("tpf1").set("SurfaceTensionCoefficient", "userdef");
    model.component("comp1").multiphysics("tpf1").set("sigma", "sigma");

    model.component("comp1").mesh("mesh2").label("Mesh 1");
    model.component("comp1").mesh("mesh2").feature("size").set("hauto", 3);
    model.component("comp1").mesh("mesh2").feature("ftri2").feature("dis1").active(false);
    model.component("comp1").mesh("mesh2").feature("ftri2").feature("dis1").set("numelem", 20);
    model.component("comp1").mesh("mesh2").feature("ftri2").feature("size1").set("hauto", 2);
    model.component("comp1").mesh("mesh2").feature("ftri2").feature("size1").set("custom", "on");
    model.component("comp1").mesh("mesh2").feature("ftri2").feature("size1").set("hmax", 5);
    model.component("comp1").mesh("mesh2").feature("ftri2").feature("size1").set("hmaxactive", true);
    model.component("comp1").mesh("mesh2").feature("ftri2").feature("size1").set("hmin", 0.00375);
    model.component("comp1").mesh("mesh2").feature("ftri2").feature("size1").set("hminactive", false);
    model.component("comp1").mesh("mesh2").run();

    model.component("comp1").probe("point1").set("expr", "ht.T");
    model.component("comp1").probe("point1").set("unit", "K");
    model.component("comp1").probe("point1").set("descr", "Temperature");
    model.component("comp1").probe("point1").set("table", "tbl1");
    model.component("comp1").probe("point1").set("window", "window2");

    model.study().create("std1");
    model.study("std1").create("phasei", "PhaseInitialization");
    model.study("std1").create("time", "Transient");
    model.study("std1").feature("phasei")
         .set("activate", new String[]{"spf", "off", "pf", "on", "viscode", "on", "ht", "on", "cpode", "on", 
         "frame:spatial1", "on", "frame:material1", "on"});

    model.sol().create("sol1");
    model.sol("sol1").study("std1");
    model.sol("sol1").attach("std1");
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
    model.sol("sol1").feature("t1").create("se1", "Segregated");
    model.sol("sol1").feature("t1").create("d1", "Direct");
    model.sol("sol1").feature("t1").create("d2", "Direct");
    model.sol("sol1").feature("t1").create("d3", "Direct");
    model.sol("sol1").feature("t1").create("i1", "Iterative");
    model.sol("sol1").feature("t1").create("i2", "Iterative");
    model.sol("sol1").feature("t1").create("i3", "Iterative");
    model.sol("sol1").feature("t1").feature("se1").create("ss1", "SegregatedStep");
    model.sol("sol1").feature("t1").feature("se1").create("ss2", "SegregatedStep");
    model.sol("sol1").feature("t1").feature("se1").create("ss3", "SegregatedStep");
    model.sol("sol1").feature("t1").feature("se1").create("ll1", "LowerLimit");
    model.sol("sol1").feature("t1").feature("se1").create("cpode1", "SegregatedStep");
    model.sol("sol1").feature("t1").feature("se1").feature().remove("ssDef");
    model.sol("sol1").feature("t1").feature("i1").create("mg1", "Multigrid");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("pr").create("so1", "SOR");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("po").create("so1", "SOR");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("cs").create("d1", "Direct");
    model.sol("sol1").feature("t1").feature("i2").create("mg1", "Multigrid");
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("pr").create("sc1", "SCGS");
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("po").create("sc1", "SCGS");
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("cs").create("d1", "Direct");
    model.sol("sol1").feature("t1").feature("i3").create("mg1", "Multigrid");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("pr").create("sc1", "SCGS");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("po").create("sc1", "SCGS");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("cs").create("d1", "Direct");
    model.sol("sol1").feature("t1").feature().remove("fcDef");

    model.result().dataset().create("step1_ds1", "Grid1D");
    model.result().dataset().create("an1_ds1", "Grid1D");
    model.result().dataset().create("dset3", "Solution");
    model.result().dataset().create("avh1", "Average");
    model.result().dataset("step1_ds1").set("data", "none");
    model.result().dataset("an1_ds1").set("data", "none");
    model.result().dataset("dset3").set("probetag", "point1");
    model.result().dataset("avh1").set("probetag", "point1");
    model.result().dataset("avh1").set("data", "dset3");
    model.result().dataset("avh1").selection().geom("geom2", 0);
    model.result().numerical().create("pev1", "EvalPoint");
    model.result().numerical("pev1").set("probetag", "point1");
    model.result().create("pg1", "PlotGroup2D");
    model.result().create("pg2", "PlotGroup2D");
    model.result().create("pg3", "PlotGroup2D");
    model.result().create("pg4", "PlotGroup2D");
    model.result().create("pg5", "PlotGroup1D");
    model.result().create("pg6", "PlotGroup1D");
    model.result().create("pg7", "PlotGroup1D");
    model.result().create("pg8", "PlotGroup2D");
    model.result("pg1").create("surf1", "Surface");
    model.result("pg1").create("con1", "Contour");
    model.result("pg1").feature("surf1").set("expr", "ht.T");
    model.result("pg1").feature("con1").set("expr", "pf.Vf2");
    model.result("pg2").create("surf1", "Surface");
    model.result("pg2").create("con1", "Contour");
    model.result("pg2").feature("surf1").set("expr", "pf.Vf2");
    model.result("pg2").feature("con1").set("expr", "p");
    model.result("pg3").create("surf1", "Surface");
    model.result("pg3").create("con1", "Contour");
    model.result("pg3").feature("con1").set("expr", "pf.Vf2");
    model.result("pg4").create("surf1", "Surface");
    model.result("pg4").create("con1", "Contour");
    model.result("pg4").feature("surf1").set("expr", "muair*(1-phipf)/2 + comp1.mu2*(1+phipf)/2");
    model.result("pg4").feature("con1").set("expr", "pf.Vf2");
    model.result("pg5").create("plot1", "Function");
    model.result("pg5").feature("plot1").set("expr", "comp1.mumat(x[1/m][K])");
    model.result("pg6").create("plot1", "Function");
    model.result("pg6").feature("plot1").set("expr", "comp1.cpmat(T)");
    model.result("pg7").set("probetag", "window2_default");
    model.result("pg7").create("tblp1", "Table");
    model.result("pg7").feature("tblp1").set("probetag", "point1");
    model.result("pg8").create("surf1", "Surface");
    model.result("pg8").feature("surf1").set("expr", "-g_const*(1+phipf)/2 * rhomat");
    model.result().export().create("anim1", "Animation");
    model.result().export().create("anim2", "Animation");
    model.result().export().create("anim3", "Animation");
    model.result().export().create("anim4", "Animation");

    model.component("comp1").probe("point1").genResult(null);

    model.result("pg9").tag("pg7");

    model.nodeGroup().create("grp1", "Definitions", "comp1");
    model.nodeGroup("grp1").set("type", "func");
    model.nodeGroup("grp1").placeAfter(null);

    return model;
  }

  public static Model run2(Model model) {

    model.study("std1").feature("time").set("tlist", "range(0,0.01,5)");
    model.study("std1").feature("time").set("useinitsol", true);
    model.study("std1").feature("time").set("initstudy", "std1");
    model.study("std1").feature("time").set("solnum", "auto");

    model.sol("sol1").attach("std1");
    model.sol("sol1").feature("st1").label("Compile Equations: Phase Initialization");
    model.sol("sol1").feature("v1").label("Dependent Variables 1.1");
    model.sol("sol1").feature("v1").set("clist", new String[]{"1e-10[s]"});
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
    model.sol("sol1").feature("v2").set("notsoluse", "sol2");
    model.sol("sol1").feature("v2").set("notsolnum", "auto");
    model.sol("sol1").feature("v2").set("clist", new String[]{"range(0,0.01,5)", "1e-10[s]"});
    model.sol("sol1").feature("v2").feature("comp1_phipf").set("scalemethod", "manual");
    model.sol("sol1").feature("v2").feature("comp1_phipf").set("scaleval", 1);
    model.sol("sol1").feature("t1").label("Time-Dependent Solver 1.1");
    model.sol("sol1").feature("t1").set("control", "time");
    model.sol("sol1").feature("t1").set("tlist", "range(0,0.01,5)");
    model.sol("sol1").feature("t1").set("rtol", 0.005);
    model.sol("sol1").feature("t1").set("atolglobalfactor", 0.05);
    model.sol("sol1").feature("t1")
         .set("atolmethod", new String[]{"comp1_GI", "global", "comp1_mu2", "global", "comp1_p", "scaled", "comp1_phipf", "global", "comp1_psi", "global", 
         "comp1_T", "global", "comp1_u", "global", "comp1_Cp2", "global"});
    model.sol("sol1").feature("t1")
         .set("atolfactor", new String[]{"comp1_GI", "0.1", "comp1_mu2", "0.1", "comp1_p", "1", "comp1_phipf", "0.1", "comp1_psi", "0.1", 
         "comp1_T", "0.1", "comp1_u", "0.1", "comp1_Cp2", "0.1"});
    model.sol("sol1").feature("t1").set("initialstepbdf", "1e-10");
    model.sol("sol1").feature("t1").set("initialstepbdfactive", true);
    model.sol("sol1").feature("t1").set("maxorder", 2);
    model.sol("sol1").feature("t1").set("stabcntrl", true);
    model.sol("sol1").feature("t1").set("bwinitstepfrac", 0.01);
    model.sol("sol1").feature("t1").set("estrat", "exclude");
    model.sol("sol1").feature("t1").set("rescaleafterinitbw", true);
    model.sol("sol1").feature("t1").feature("dDef").label("Direct 4");
    model.sol("sol1").feature("t1").feature("aDef").label("Advanced 1");
    model.sol("sol1").feature("t1").feature("aDef").set("cachepattern", true);
    model.sol("sol1").feature("t1").feature("se1").label("Segregated 1.1");
    model.sol("sol1").feature("t1").feature("se1").set("ntolfact", 0.5);
    model.sol("sol1").feature("t1").feature("se1").set("segstabacc", "segaacc");
    model.sol("sol1").feature("t1").feature("se1").set("segaaccdim", 5);
    model.sol("sol1").feature("t1").feature("se1").set("segaaccmix", 0.9);
    model.sol("sol1").feature("t1").feature("se1").feature("ss1").label("Merged variables");
    model.sol("sol1").feature("t1").feature("se1").feature("ss1").set("segvar", new String[]{"comp1_T", "comp1_mu2"});
    model.sol("sol1").feature("t1").feature("se1").feature("ss1").set("linsolver", "d1");
    model.sol("sol1").feature("t1").feature("se1").feature("ss1").set("subdamp", "0.8");
    model.sol("sol1").feature("t1").feature("se1").feature("ss1").set("subjtech", "once");
    model.sol("sol1").feature("t1").feature("se1").feature("ss2").label("Velocity u, Pressure p");
    model.sol("sol1").feature("t1").feature("se1").feature("ss2").set("segvar", new String[]{"comp1_u", "comp1_p"});
    model.sol("sol1").feature("t1").feature("se1").feature("ss2").set("linsolver", "d2");
    model.sol("sol1").feature("t1").feature("se1").feature("ss2").set("subdamp", "0.8");
    model.sol("sol1").feature("t1").feature("se1").feature("ss2").set("subjtech", "once");
    model.sol("sol1").feature("t1").feature("se1").feature("ss3").label("Phase field variables");
    model.sol("sol1").feature("t1").feature("se1").feature("ss3")
         .set("segvar", new String[]{"comp1_phipf", "comp1_psi"});
    model.sol("sol1").feature("t1").feature("se1").feature("ss3").set("linsolver", "d3");
    model.sol("sol1").feature("t1").feature("se1").feature("ss3").set("subdamp", "0.8");
    model.sol("sol1").feature("t1").feature("se1").feature("ss3").set("subjtech", "once");
    model.sol("sol1").feature("t1").feature("se1").feature("ll1").label("Lower Limit 1.1");
    model.sol("sol1").feature("t1").feature("se1").feature("ll1").set("lowerlimit", "comp1.T 0 ");
    model.sol("sol1").feature("t1").feature("se1").feature("cpode1").label("Segregated Step 1");
    model.sol("sol1").feature("t1").feature("se1").feature("cpode1").set("segvar", new String[]{"comp1_Cp2"});
    model.sol("sol1").feature("t1").feature("se1").feature("cpode1").set("linsolver", "d3");
    model.sol("sol1").feature("t1").feature("se1").feature().remove("spf1");
    model.sol("sol1").feature("t1").feature("se1").feature().remove("ht1");
    model.sol("sol1").feature("t1").feature("d1").label("Direct, heat transfer variables (ht) (merged)");
    model.sol("sol1").feature("t1").feature("d1").set("linsolver", "pardiso");
    model.sol("sol1").feature("t1").feature("d1").set("pivotperturb", 1.0E-13);
    model.sol("sol1").feature("t1").feature("d2").label("Direct, fluid flow variables (spf)");
    model.sol("sol1").feature("t1").feature("d2").set("linsolver", "pardiso");
    model.sol("sol1").feature("t1").feature("d2").set("pivotperturb", 1.0E-13);
    model.sol("sol1").feature("t1").feature("d3").label("Direct, phase field variables (pf)");
    model.sol("sol1").feature("t1").feature("d3").set("linsolver", "pardiso");
    model.sol("sol1").feature("t1").feature("d3").set("pivotperturb", 1.0E-13);
    model.sol("sol1").feature("t1").feature("i1").label("AMG, heat transfer variables (ht)");
    model.sol("sol1").feature("t1").feature("i1").set("rhob", 20);
    model.sol("sol1").feature("t1").feature("i1").feature("ilDef").label("Incomplete LU 1");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").label("Multigrid 1.1");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").set("prefun", "saamg");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").set("maxcoarsedof", 50000);
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").set("saamgcompwise", true);
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").set("usesmooth", false);
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("pr").label("Presmoother 1");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("pr").feature("soDef").label("SOR 2");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("pr").feature("so1").label("SOR 1.1");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("pr").feature("so1").set("relax", 0.9);
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("po").label("Postsmoother 1");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("po").feature("soDef").label("SOR 2");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("po").feature("so1").label("SOR 1.1");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("po").feature("so1").set("relax", 0.9);
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("cs").label("Coarse Solver 1");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("cs").feature("dDef").label("Direct 2");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("cs").feature("d1").label("Direct 1.1");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("cs").feature("d1")
         .set("linsolver", "pardiso");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("cs").feature("d1")
         .set("pivotperturb", 1.0E-13);
    model.sol("sol1").feature("t1").feature("i2").label("AMG, fluid flow variables (spf)");
    model.sol("sol1").feature("t1").feature("i2").set("maxlinit", 100);
    model.sol("sol1").feature("t1").feature("i2").set("rhob", 20);
    model.sol("sol1").feature("t1").feature("i2").feature("ilDef").label("Incomplete LU 1");
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").label("Multigrid 1.1");
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").set("prefun", "saamg");
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").set("maxcoarsedof", 80000);
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").set("strconn", 0.02);
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").set("saamgcompwise", true);
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").set("usesmooth", false);
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("pr").label("Presmoother 1");
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("pr").feature("soDef").label("SOR 1");
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("pr").feature("sc1").label("SCGS 1.1");
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("pr").feature("sc1")
         .set("linesweeptype", "ssor");
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("pr").feature("sc1").set("iter", 0);
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("pr").feature("sc1").set("approxscgs", true);
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("pr").feature("sc1")
         .set("scgsdirectmaxsize", 1000);
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("po").label("Postsmoother 1");
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("po").feature("soDef").label("SOR 1");
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("po").feature("sc1").label("SCGS 1.1");
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("po").feature("sc1")
         .set("linesweeptype", "ssor");
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("po").feature("sc1").set("iter", 1);
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("po").feature("sc1").set("approxscgs", true);
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("po").feature("sc1")
         .set("scgsdirectmaxsize", 1000);
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("cs").label("Coarse Solver 1");
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("cs").feature("dDef").label("Direct 2");
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("cs").feature("d1").label("Direct 1.1");
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("cs").feature("d1")
         .set("linsolver", "pardiso");
    model.sol("sol1").feature("t1").feature("i2").feature("mg1").feature("cs").feature("d1")
         .set("pivotperturb", 1.0E-13);
    model.sol("sol1").feature("t1").feature("i3").label("AMG, phase field variables (pf)");
    model.sol("sol1").feature("t1").feature("i3").set("maxlinit", 50);
    model.sol("sol1").feature("t1").feature("i3").set("rhob", 20);
    model.sol("sol1").feature("t1").feature("i3").feature("ilDef").label("Incomplete LU 1");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").label("Multigrid 1.1");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").set("prefun", "saamg");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").set("maxcoarsedof", 50000);
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").set("usesmooth", false);
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("pr").label("Presmoother 1");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("pr").feature("soDef").label("SOR 1");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("pr").feature("sc1").label("SCGS 1.1");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("pr").feature("sc1")
         .set("linesweeptype", "ssor");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("pr").feature("sc1").set("iter", 0);
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("pr").feature("sc1").set("approxscgs", true);
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("pr").feature("sc1")
         .set("scgsdirectmaxsize", 1000);
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("po").label("Postsmoother 1");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("po").feature("soDef").label("SOR 1");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("po").feature("sc1").label("SCGS 1.1");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("po").feature("sc1")
         .set("linesweeptype", "ssor");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("po").feature("sc1").set("iter", 1);
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("po").feature("sc1").set("approxscgs", true);
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("po").feature("sc1")
         .set("scgsdirectmaxsize", 1000);
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("cs").label("Coarse Solver 1");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("cs").feature("dDef").label("Direct 2");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("cs").feature("d1").label("Direct 1.1");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("cs").feature("d1")
         .set("linsolver", "pardiso");
    model.sol("sol1").feature("t1").feature("i3").feature("mg1").feature("cs").feature("d1")
         .set("pivotperturb", 1.0E-13);
    model.sol("sol1").runAll();

    model.result().dataset("step1_ds1").set("function", "step1");
    model.result().dataset("step1_ds1").set("parmin1", 312.5);
    model.result().dataset("step1_ds1").set("parmax1", 327.5);
    model.result().dataset("an1_ds1").set("function", "all");
    model.result().dataset("an1_ds1").set("par1", "T");
    model.result().dataset("an1_ds1").set("parmin1", 300);
    model.result().dataset("an1_ds1").set("parmax1", 400);
    model.result().dataset("an1_ds1").set("res1", 10000);
    model.result().dataset("an1_ds1").set("distribution", "mixed");
    model.result().dataset("dset3").label("Probe Solution 3");
    model.result().numerical("pev1").setResult();
    model.result("pg1").label("Temperature and Volume Fraction");
    model.result("pg1").set("looplevel", new int[]{3});
    model.result("pg1").set("titletype", "manual");
    model.result("pg1").set("title", "Temperature (K) with Volume Fraction Contours");
    model.result("pg1").set("paramindicator", "t=eval(t,s) (s)");
    model.result("pg1").set("showlegendsmaxmin", true);
    model.result("pg1").set("showlegendsunit", true);
    model.result("pg1").feature("surf1").set("rangecoloractive", true);
    model.result("pg1").feature("surf1").set("rangecolormin", 300);
    model.result("pg1").feature("surf1").set("rangecolormax", 399.9999999999998);
    model.result("pg1").feature("surf1").set("colortable", "Magma");
    model.result("pg1").feature("surf1").set("resolution", "normal");
    model.result("pg1").feature("con1").set("number", 10);
    model.result("pg1").feature("con1").set("coloring", "uniform");
    model.result("pg1").feature("con1").set("color", "green");
    model.result("pg1").feature("con1").set("resolution", "normal");
    model.result("pg2").label("Pressure and Volume Fraction");
    model.result("pg2").set("looplevel", new int[]{2});
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
    model.result("pg3").set("looplevel", new int[]{2});
    model.result("pg3").set("titletype", "manual");
    model.result("pg3").set("title", "Velocity magnitude (m/s) with Volume Fraction contours");
    model.result("pg3").set("paramindicator", "t=eval(t) (s)");
    model.result("pg3").feature("surf1").set("colortable", "Magma");
    model.result("pg3").feature("surf1").set("resolution", "normal");
    model.result("pg3").feature("con1").set("number", 10);
    model.result("pg3").feature("con1").set("coloring", "uniform");
    model.result("pg3").feature("con1").set("color", "green");
    model.result("pg3").feature("con1").set("resolution", "normal");
    model.result("pg4").label("Viscosity");
    model.result("pg4").set("looplevel", new int[]{2});
    model.result("pg4").set("titletype", "manual");
    model.result("pg4").set("title", "Fluid viscosity (Pa*s) with Volume Fraction contours");
    model.result("pg4").set("paramindicator", "t=eval(t,s) (s)");
    model.result("pg4").feature("surf1").set("coloring", "gradient");
    model.result("pg4").feature("surf1").set("topcolor", "magenta");
    model.result("pg4").feature("surf1").set("bottomcolor", "white");
    model.result("pg4").feature("surf1").set("resolution", "normal");
    model.result("pg4").feature("con1").set("number", 10);
    model.result("pg4").feature("con1").set("resolution", "normal");
    model.result("pg5").label("Material Viscosity");
    model.result("pg5").set("data", "step1_ds1");
    model.result("pg5").set("solrepresentation", "solnum");
    model.result("pg5").set("titletype", "manual");
    model.result("pg5").set("title", "Viscosity of Material (Pa*s)");
    model.result("pg5").set("xlabel", "T (K)");
    model.result("pg5").set("xlabelactive", true);
    model.result("pg5").set("ylabel", "Pa*s");
    model.result("pg5").set("ylabelactive", true);
    model.result("pg5").feature("plot1").set("solrepresentation", "solnum");
    model.result("pg5").feature("plot1").set("descr", "mumat(x)");
    model.result("pg5").feature("plot1").set("xdataexpr", "x");
    model.result("pg5").feature("plot1").set("xdataunit", "m");
    model.result("pg5").feature("plot1").set("xdatadescractive", true);
    model.result("pg5").feature("plot1").set("xdatadescr", "x (K)");
    model.result("pg5").feature("plot1").set("lowerbound", 312.5);
    model.result("pg5").feature("plot1").set("upperbound", 327.5);
    model.result("pg5").feature("plot1").set("linewidth", "preference");
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
    model.result("pg8").label("Weight");
    model.result("pg8").set("looplevel", new int[]{2});
    model.result("pg8").feature("surf1").set("resolution", "normal");
    model.result().export("anim1").label("T and Vf2");
    model.result().export("anim1")
         .set("giffilename", "/nas/longleaf/home/mshourya/workspace/SelectiveLaserSintering/Particle/Results/PhysicalTemperatureDependent/T+Vf2.gif");
    model.result().export("anim1").set("maxframes", 201);
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
         .set("giffilename", "/nas/longleaf/home/mshourya/workspace/SelectiveLaserSintering/Particle/Results/PhysicalTemperatureDependent/Vf2+p.gif");
    model.result().export("anim2").set("maxframes", 201);
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
    model.result().export("anim3").set("target", "player");
    model.result().export("anim3")
         .set("giffilename", "/nas/longleaf/home/mshourya/workspace/SelectiveLaserSintering/Particle/Results/PhysicalTemperatureDependent/U+Vf2.gif");
    model.result().export("anim3").set("framesel", "all");
    model.result().export("anim3").set("showframe", 222);
    model.result().export("anim3").set("shownparameter", "2.21");
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
    model.result().export("anim4").set("target", "player");
    model.result().export("anim4")
         .set("giffilename", "/nas/longleaf/home/mshourya/workspace/SelectiveLaserSintering/Particle/Results/PhysicalTemperatureDependent/mu+Vf2.gif");
    model.result().export("anim4").set("maxframes", 201);
    model.result().export("anim4").set("showframe", 194);
    model.result().export("anim4").set("shownparameter", "0.482");
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

    return model;
  }

  public static Model run3(Model model) {
    model.result().export("anim4").set("logo3d", "on");
    model.result().export("anim4").set("options3d", "off");
    model.result().export("anim4").set("axisorientation", "on");
    model.result().export("anim4").set("grid", "on");
    model.result().export("anim4").set("axes1d", "on");
    model.result().export("anim4").set("axes2d", "on");
    model.result().export("anim4").set("showgrid", "on");

    model.nodeGroup("grp1").label("Materials Properties");
    model.nodeGroup("grp1").add("func", "step1");

    return model;
  }

  public static void main(String[] args) {
    Model model = run();
    model = run2(model);
    run3(model);
  }

}
