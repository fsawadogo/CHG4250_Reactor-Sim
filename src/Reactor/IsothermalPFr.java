
package Reactor;

import NumericalMethod.Function;
import Reactions.PrimaryReaction;
import Reactions.Reaction;
import Reactions.SecondaryReaction;
import UserInterface.FileReader;
import UserInterface.MainClass;

public class IsothermalPFr implements Function {
    public static final double R = 8.314;//Gas law constant (J/mol·K)
    public static final double [] Cp298 ={0,0,0,0,0};//EDC,NH3,EDA,HCl,VC,NH4Cl
    public static final double a= 4.602;// (m2/m3) surface
    public static final double mc = 55;// mass flow of the coolant
    //public static final double Ta = 298;// heat of the
    public static final double Cpc = 31.5;// heat capacity of the coolant
    public static final double U = 500;// overall heat transfer coefficient (J/m2.s.K)
    private final PrimaryReaction R1;
    private final SecondaryReaction R2;
    private final double[] g_opConditions;
    private final double v0 = 61.067; //m3/h

    //Constructor
    public IsothermalPFr(PrimaryReaction R1, SecondaryReaction R2) {
        this.R1 = R1;
        this.R2 = R2;
        this.g_opConditions = FileReader.g_opConditions;
        // deep copying the FileReader g_opConditions array to the private instance variable array g_opConditions in this class
    }
    //Copy constructor
    public IsothermalPFr(IsothermalPFr source){

        if (source == null) throw new IllegalArgumentException ("Source PFR cannot be null");
        if(source.R1 == null) throw new IllegalArgumentException("The reaction object is null");
        this.R1=source.R1;
        if(source.R2==null) throw new IllegalArgumentException("The reaction object is null");
        this.R2=source.R2;
        if (source.g_opConditions==null) throw new IllegalArgumentException("The g_opConditions array object is null");
        this.g_opConditions = new double [source.g_opConditions.length];
        System.arraycopy(source.g_opConditions, 0, this.g_opConditions, 0, source.g_opConditions.length);
    }

    // This method is overridden from the implemented from the interface Function
    // The equations of each ODE are coded here allowing for these to be calculated at every iteration
    public double[] calculateValue(double x, double[] y, int index){
        // These are the values calculated for each parameter in the y array to clarify which parameter is associated with which element
        double  [] HRX = new double[] {R1.calculateReactionEnthalpy(y[0]), R2.calculateReactionEnthalpy(y[0])};
        double T = y[0]; // Temperature (K)
        double P = y[1]; // Pressure (atm)
        double F_EDC = y[2]; // Flow rate of EDC (kmol/h)
        double F_NH3 = y[3]; // Flow rate of Ammonia (kmol/h)
        double F_EDA = y[4]; // Flow rate of EDA (kmol/h)
        double F_HCl = y[5]; // Flow rate of HCl(kmol/h)
        double F_VC = y[6]; // Flow rate of HCl(kmol/h)
        double F_NH4Cl = y[7]; // Flow rate of HCl(kmol/h)
        double Ta = y[8];

        double F_T = F_EDC + F_NH3 + F_EDA + F_HCl+ F_VC+ F_NH4Cl; // Total molar flow rate
        double C_EDC= (F_EDC/v0);
        double C_NH3 = (F_NH3/v0);
        double r1 = -R1.calculateReactionRate(T,C_EDC, C_NH3);//change the values to match the variables that give the proper value
        double r2 = -R2.calculateReactionRate(T,C_EDC,C_NH3); // Reaction rate (mol/kg·s)
        double rEDC = r1+r2;
        double rNH3 = 1./2*r1+r2;
        double rEDA = -r1;
        double rHCl = -r1;
        double rVC = -r2;
        double rNH4Cl = -r2;

        double[] rhs = new double[9];

        double [] Fi = {y[2],y[3],y[4],y[5],y[6],y[7]};
        double sumFjCpj = 0;
        for(int i=0;i<Fi.length;i++){
            sumFjCpj += Fi[i]* Reaction.calculateHeatCapacity(T)[i];// denominator of the differential equation dT/dV
            }
        double X =(MainClass.F_EDC0-F_EDC)/F_EDC;//conversion



        //double a = A/V;
        double term1 = rEDC*(HRX[0]+HRX[1])+rNH3*(HRX[0]+HRX[1])+rNH3*(HRX[0])+rNH3*(HRX[0])+rNH3*(HRX[0])+rNH3*(HRX[0]);
        double term2 = U*a*(T-Ta);
        //double sumFjCj += Fi[i]*Reaction.Cp[i]// denominator of the differential equation dT/dV

        rhs[0] = 0; // dT/dV, temperature change with multiple reaction
        rhs[1] =  U*a*(Ta-y[0])/(mc*Cpc);//coolant temperature Ua(Ta-T)/mcCpc counter-current heat exchange
        rhs[2] = 0; //dP/dV
        rhs[3] = rEDC; // dF_EDC/dV
        rhs[4] = rNH3; // dF_NH3/dV
        rhs[5] = rEDA; // dF_EDA/dV
        rhs[6] = rHCl; // dF_HCl/dV
        rhs[7] =rVC;// dF_VC/dV
        rhs[8] = rNH4Cl;// dF_NH4Cl/dV

        System.out.println("The molar flow rate of EDC is " + y[2] + "kmol/h"+" ,the conversion is " +X+ "the volume of the reactor is"+ x );
        return rhs;

    }

    // method that returns the results calculated in the ODE solver class (RungeKutta 4)
    public double[] integrateODES(double x_0, double x_f, double[] y_0, double delx, int maxIt, Function f) {
        return NumericalMethod.RungeKutta4.rungeKutta4Integrate(x_0, x_f, y_0, delx, maxIt, this);
    }
    // equals method is inherited from the parent class PBR and is simply called on


    }





