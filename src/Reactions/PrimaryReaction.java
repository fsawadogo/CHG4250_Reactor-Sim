package Reactions;
public class PrimaryReaction extends Reaction {
    public static final double k1_0 = 1.12;//s-1
    private double k1;

        public PrimaryReaction(double k1) {
            if (k1 < 0.) throw new IllegalArgumentException("The rate constant cannot be negative");
            this.k1 = k1;
        }

        public PrimaryReaction(Reactions.PrimaryReaction source) {
            if (source == null)
                throw new IllegalArgumentException("The PrimaryReaction object sent to the copy constructor is null");
            this.k1 = source.k1;
        }

        public double getk1() {
            return this.k1;
        }

        public boolean setk1(double k1) {
            if (k1 < 0.)
                return false;
            this.k1 = k1;
            return true;
        }



    public double calculateReactionRate(double T, double C_EDC, double C_NH3) {
        if (T <= 0 || C_EDC < 0 || C_NH3 < 0) {
            throw new IllegalArgumentException("Invalid Temperature or concentration");
        }
        //Derivative of Arrhenius equation
        k1 = Math.exp(14.288 - 53775 / (R * T));
        return k1 * C_EDC * Math.pow(C_NH3, 2);

    }

    public double calculateReactionEnthalpy(double T) {
        double H0_RX = Reaction.Ho_f[2] + Reaction.Ho_f[3] - Reaction.Ho_f[0] - 2 * Reaction.Ho_f[1];
        double Cp_prod = Reaction.CP298[0] * (T - Reaction.TR) + 1. / 2 * (Reaction.A[0] * (T - Reaction.TR))
                + 1. / 2 * (Reaction.B[0] * (Math.pow(T - Reaction.TR, 2))) + 1. / 3 * (Reaction.C[0] * Math.pow(T - Reaction.TR, 3))
                + 1. / 4 * (Reaction.D[0] * Math.pow(T - Reaction.TR, 4)) + (Reaction.E[0] * 1 / (T - Reaction.TR));
        double Cp_rea = Reaction.CP298[1] * (T - Reaction.TR) +2*( 1. / 2 * (Reaction.A[1] * (T - Reaction.TR))
                + 1. / 2 * (Reaction.B[1] * (Math.pow(T - Reaction.TR, 2))) + 1. / 3 * (Reaction.C[1] * Math.pow(T - Reaction.TR, 3))
                + 1. / 4 * (Reaction.D[1] * Math.pow(T - Reaction.TR, 4)) + (Reaction.E[1] * 1 / (T - Reaction.TR)));
        double delta_CP = Cp_prod - Cp_rea;
        double HRX;
        HRX = H0_RX + delta_CP;
        return HRX;
    }
}
    /*// default constructor
    public PrimaryReaction(double k) {
        super(k);
    }

    //copy constructor
    public PrimaryReaction(PrimaryReaction source) {
        super(source);
    }

    // clone
    public PrimaryReaction clone() {
        return new PrimaryReaction(this);
    }

    @Override
    public double calculateReactionRate(double T, double C_EDC, double C_NH3) {
        if (T <= 0 || C_EDC < 0 || C_NH3 < 0) //boundary set to ensure that the values of concentration and temperature are physically attainable before running the simulation
        {
            throw new IllegalArgumentException("Invalid temperature or concentration.");
        }
        return -this.getK() * C_EDC * Math.pow(C_NH3, 2);
    }

    @Override
    public double[] calculateHeatCapacities(double T, double[] A, double[] B, double[] C, double[] D) {
        double[] Cp = new double[6];
        Cp[0] = MainClass.CP298[0]; //heat capacity in function of temperature
        Cp[1] = MainClass.A[0] + MainClass.B[0] * T + MainClass.C[0] * Math.pow(T, 2) + MainClass.D[0] * Math.pow(T, 3) + MainClass.E[0] * (1 / Math.pow(T, 2.));//heat capacity in function of temperature
        Cp[2] = MainClass.CP298[1];//heat capacity in function of temperature
        Cp[3] = MainClass.A[1] + MainClass.B[1] * T + MainClass.C[1] * Math.pow(T, 2) + MainClass.D[1] * Math.pow(T, 3) + MainClass.E[1] * (1 / Math.pow(T, 2.));//heat capacity in function of temperature
        return Cp;



        @Override
        public boolean equals (Object comparator){
            return super.equals(comparator);
        }
    }

    @Override
    public double calculateReactionEnthalpy(double T) {
        double H0_RX = MainClass.Ho_f[2] + MainClass.Ho_f[3] - MainClass.Ho_f[0] - 2 * MainClass.Ho_f[1];
        double Cp_prod = MainClass.CP298[0] * (T - MainClass.TR) + 1. / 2 * (MainClass.A[0] * (T - MainClass.TR)) + 1. / 2 * (MainClass.B[0] * (Math.pow(T - MainClass.TR, 2))) + 1. / 3 * (MainClass.C[0] * Math.pow(T - MainClass.TR, 3)) +
                1. / 4 * (MainClass.D[0] * Math.pow(T - MainClass.TR, 4)) + (MainClass.E[0] * 1 / (T - MainClass.TR));
        double Cp_rea = MainClass.CP298[1] * (T - MainClass.TR) + 1. / 2 * (MainClass.A[1] * (T - MainClass.TR)) + 1. / 2 * (MainClass.B[1] * (Math.pow(T - MainClass.TR, 2))) + 1. / 3 * (MainClass.C[1] * Math.pow(T - MainClass.TR, 3)) +
                1. / 4 * (MainClass.D[1] * Math.pow(T - MainClass.TR, 4)) + (MainClass.E[1] * 1 / (T - MainClass.TR));
        double delta_CP = Cp_prod - Cp_rea;
        double HRX= H0_RX + delta_CP;
        return HRX;
    }
}


*/
