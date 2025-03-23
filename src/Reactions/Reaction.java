package Reactions;
public abstract class Reaction {
    public final static double R = 8.314;
    public final static double TR = 298; // reference temperature in K
    //Parameters to calculate the heat capacities of each component as a function of temperature EDC, NH3, EDA, HCl, VC, NH4Cl
    public final static double[] A = {2,19.99563,2, 32.12392,2,2101.470};//(J/mol.K) NH3 and HCl
    public final static double[] B = {2,49.77119,2,-13.45805,2,-11194.40};//(J/mol.K) NH3 and HCl
    public final static double[] C = {2,-15.37599,2, 19.86852,2,23342.30};//(J/mol.K) NH3 and HCl
    public final static double[] D = {2,1.921168,2, -6.853936,2,-16996.70};//(J/mol.K) NH3 and HCl
    public final static double[] E = {2,0.189174,2,-0.0496772,2,-26.83530};// (J/mol.K) NH3 and HCl
    public final static double[] CP298 = {122.2, 178.78,87.69,86.45}; //Heat capacities at 298 K(J/mol.K) for EDC and EDA, NH4Cl
    public final static double [] Ho_f ={-167.2,-45.90,-63.01,-92.31,22.,-314.55};//kJ/mol heat of reaction at the reference temperature Tr of 298K

    public static double[] calculateHeatCapacity (double T){
        double[] Cp = new double[6];
        Cp[0] = Reaction.CP298[0]; //EDC heat capacity in function of temperature
        Cp[1] = Reaction.A[0] + Reaction.B[0] * T + Reaction.C[0] * Math.pow(T, 2) + Reaction.D[0] * Math.pow(T, 3) + Reaction.E[0] * (1 / Math.pow(T, 2.));//NH3 heat capacity in function of temperature
        Cp[2] = Reaction.CP298[1];//VC heat capacity in function of temperature
        Cp[3] = Reaction.A[5] + Reaction.B[5] * T + Reaction.C[5] * Math.pow(T, 2) + Reaction.D[5] * Math.pow(T, 3) + Reaction.E[5] * (1 / Math.pow(T, 2.));//NH4Cl heat capacity in function of temperature
        Cp[4] = Reaction.A[5] + Reaction.B[5] * T + Reaction.C[5] * Math.pow(T, 2) + Reaction.D[5] * Math.pow(T, 3) + Reaction.E[5] * (1 / Math.pow(T, 2.));//
        Cp[5] = Reaction.A[5] + Reaction.B[5] * T + Reaction.C[5] * Math.pow(T, 2) + Reaction.D[5] * Math.pow(T, 3) + Reaction.E[5] * (1 / Math.pow(T, 2.));//
        return Cp;
    };

    public abstract double calculateReactionRate (double T, double C_EDC, double C_NH3);
    public abstract double calculateReactionEnthalpy(double T);


}
    /*public static final double K433 = 10.; // Reaction rate constant at 433 K
    private double k;// Reaction rate constant dependent of temperature

    // Standard Constructor
    public Reaction (double k){
        if (k<0.0) throw new IllegalArgumentException("The value of k is negative");
        this.k=k;
    }
    //Copy Constructor
    public Reaction (Reaction source){
        if (source== null) throw new IllegalArgumentException("The source passed is null");
        this.k=source.k;
    }

    //Clone Method
    @Override
    public abstract Reaction clone();


    //Accessor Method
    public double getK(){
        return this.k;
    }

    //Mutator Method
    public boolean setK(){
        if (k<0.0) return false; //boundary set; k value should be non-negative
        this.k=k;
        return true;
    }

    //Method for calculating reaction rate
    public  abstract double calculateReactionRate(double T, double C_EDC, double C_NH3);
    public abstract double[] calculateHeatCapacities (double T, double[] A, double[] B, double []C, double []D);
    public abstract double calculateReactionEnthalpy(double T );

    //Equals Method
    public boolean equals (Object comparator){
        if (comparator==null) return false;
        if (this.getClass()!=comparator.getClass())return false;
        if(this.k!=((PrimaryReaction)comparator).k) return false;
        return true;
    }
}

     */
