package Reactor;

import UserInterface.MainClass;
import NumericalMethod.Function;
import Reactions.PrimaryReaction;
import Reactions.Reaction;
import Reactions.SecondaryReaction;
import UserInterface.FileReader;

public class PFR implements Function {
	public static final double R = 8.314;// Gas law constant (J/mol·K)
	public static final double[] Cp298 = { 0, 0, 0, 0, 0 };// EDC,NH3,EDA,HCl,VC,NH4Cl
	public static final double a = 4.602;// (m2/m3) surface
	public static final double mc = 55;// mass flow of the coolant
	// public static final double Ta = 298;// heat of the
	public static final double Cpc = 31.5;// heat capacity of the coolant
	public static final double U = 500;// overall heat transfer coefficient (J/m2.s.K)
	private final PrimaryReaction R1;
	private final SecondaryReaction R2;
	private final double[] g_opConditions;

	// Constructor
	public PFR(PrimaryReaction R1, SecondaryReaction R2) {
		this.R1 = R1;
		this.R2 = R2;
		this.g_opConditions = FileReader.g_opConditions;
		// deep copying the FileReader g_opConditions array to the private instance
		// variable array g_opConditions in this class
	}

	// Copy constructor
	public PFR(PFR source) {

		if (source == null)
			throw new IllegalArgumentException("Source PFR cannot be null");
		if (source.R1 == null)
			throw new IllegalArgumentException("The reaction object is null");
		this.R1 = source.R1;
		if (source.R2 == null)
			throw new IllegalArgumentException("The reaction object is null");
		this.R2 = source.R2;
		if (source.g_opConditions == null)
			throw new IllegalArgumentException("The g_opConditions array object is null");
		this.g_opConditions = new double[source.g_opConditions.length];
		System.arraycopy(source.g_opConditions, 0, this.g_opConditions, 0, source.g_opConditions.length);
	}

	// This method is overridden from the implemented from the interface Function
	// The equations of each ODE are coded here allowing for these to be calculated
	// at every iteration
	public double[] calculateValue(double x, double[] y, int index) {
	    // Extract current state variables
	    double T = y[0]; // Temperature (K)
	    double P = Math.max(y[1], 0.1); // Pressure (bar), prevent going below 0.1 bar
	    double F_EDC = Math.max(y[2], 0); // EDC flow (kmol/h)
	    double F_NH3 = Math.max(y[3], 0); // NH3 flow (kmol/h)
	    double F_EDA = Math.max(y[4], 0); // EDA flow (kmol/h)
	    double F_HCl = Math.max(y[5], 0); // HCl flow (kmol/h)
	    double F_VC = Math.max(y[6], 0); // VC flow (kmol/h)
	    double F_NH4Cl = Math.max(y[7], 0); // NH4Cl flow (kmol/h)
	    double Ta = y[8]; // Coolant temperature (K)
	    
	    // Calculate total molar flow rate
	    double F_T = F_EDC + F_NH3 + F_EDA + F_HCl + F_VC + F_NH4Cl;
	    if (F_T < 1e-6) F_T = 1e-6; // Prevent division by zero
	    
	    // Convert pressure to Pa for ideal gas calculations
	    double P_Pa = P * 100000;
	    
	    // Calculate concentrations as mol fractions * (P/RT)
	    double C_EDC = (F_EDC / F_T) * (P_Pa / (R * T));
	    double C_NH3 = (F_NH3 / F_T) * (P_Pa / (R * T));
	    
	    // Debug output
	    System.out.printf("Debug - V: %.4f, T: %.2f, P: %.2f, C_EDC: %.6f, C_NH3: %.6f%n", 
	                     x, T, P_Pa, C_EDC, C_NH3);
	    
	    // Calculate reaction rates with error handling
	    double r1 = 0, r2 = 0;
	    try {
	        r1 = -R1.calculateReactionRate(T, C_EDC, C_NH3);
	        r2 = -R2.calculateReactionRate(T, C_EDC, C_NH3);
	    } catch (Exception e) {
	        System.err.println("Warning at x=" + x + ": " + e.getMessage());
	        // Use very small rates if calculation fails
	        r1 = -1e-10;
	        r2 = -1e-10;
	    }
	    
	    // Calculate component reaction rates
	    double rEDC = r1 + r2;
	    double rNH3 = 2 * r1 + r2; // 2 mol NH3 per mol EDC in primary reaction
	    double rEDA = -r1;
	    double rHCl = -2 * r1; // 2 mol HCl per mol EDC
	    double rVC = -r2;
	    double rNH4Cl = -r2;
	    
	    // Calculate heat effects
	    double[] HRX = {R1.calculateReactionEnthalpy(T), R2.calculateReactionEnthalpy(T)};
	    
	    // Heat generation
	    double heat_generation = r1 * HRX[0] + r2 * HRX[1];
	    
	    // Heat removal
	    double heat_removal = U * a * (T - Ta);
	    
	    // Calculate heat capacity sum
	    double[] Fi = {F_EDC, F_NH3, F_EDA, F_HCl, F_VC, F_NH4Cl};
	    double sumFjCpj = 0;
	    for (int i = 0; i < Fi.length; i++) {
	        sumFjCpj += Fi[i] * Reaction.calculateHeatCapacity(T)[i];
	    }
	    if (sumFjCpj < 1e-6) sumFjCpj = 1e-6; // Prevent division by zero
	    
	    // Temperature differential
	    double dTdV = (heat_generation - heat_removal) / sumFjCpj;
	    
	    // Limit temperature changes to prevent instability
	    double maxTempChange = 10.0; // K/m³
	    if (Math.abs(dTdV) > maxTempChange) {
	        dTdV = Math.signum(dTdV) * maxTempChange;
	    }
	    
	    // Coolant temperature differential
	    double dTadV = U * a * (T - Ta) / (mc * Cpc);
	    
	    // Pressure differential - gentle pressure drop
	    double X = (MainClass.F_EDC0 - F_EDC) / MainClass.F_EDC0;
	    if (X > 1.0) X = 1.0; // Cap conversion at 100%
	    
	    double dPdV = 0;
	    if (P > 0.1) {
	        // Gentler pressure drop model
	        dPdV = -0.01 * P * (1 + 0.05 * X);
	    }
	    
	    // Assemble differential equation results
	    double[] rhs = new double[9];
	    rhs[0] = dTdV;      // dT/dV
	    rhs[1] = dPdV;      // dP/dV
	    rhs[2] = rEDC;      // dF_EDC/dV
	    rhs[3] = rNH3;      // dF_NH3/dV
	    rhs[4] = rEDA;      // dF_EDA/dV
	    rhs[5] = rHCl;      // dF_HCl/dV
	    rhs[6] = rVC;       // dF_VC/dV
	    rhs[7] = rNH4Cl;    // dF_NH4Cl/dV
	    rhs[8] = dTadV;     // dTa/dV
	    
	    return rhs;
	}

	// method that returns the results calculated in the ODE solver class
	// (RungeKutta 4)
	public double[] integrateODES(double x_0, double x_f, double[] y_0, double delx, int maxIt, Function f) {
		return NumericalMethod.RungeKutta4.rungeKutta4Integrate(x_0, x_f, y_0, delx, maxIt, this);
	}
	// equals method is inherited from the parent class PBR and is simply called on

}
