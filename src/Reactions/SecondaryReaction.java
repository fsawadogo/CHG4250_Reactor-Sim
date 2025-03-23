package Reactions;
public class SecondaryReaction extends Reaction {
    public static final double k2_0 = 3.59165 * Math.pow(10., 11.);//rate law constant
    public static final double Ea2 = 9.93;// kJ/mol
    private double k2;

    public SecondaryReaction(double k2) {
        if (k2 < 0.) throw new IllegalArgumentException("The rate constant cannot be negative");
        this.k2 = k2;
    }

    public SecondaryReaction(SecondaryReaction source) {
        if (source == null)
            throw new IllegalArgumentException("The SecondaryReaction object sent to the copy constructor is null");
        this.k2 = source.k2;
    }

    public double getk2() {
        return this.k2;
    }

    public boolean setk1(double k2) {
        if (k2 < 0.)
            return false;
        this.k2 = k2;
        return true;
    }


@Override
public double calculateReactionRate(double T, double C_EDC, double C_NH3) {
    if (T <= 0 || C_EDC < 0 || C_NH3 < 0) {
        throw new IllegalArgumentException("Invalid Temperature or concentration");
    }
    
    // Scale down the reaction rate to prevent numerical instability
    // Use 3% of EDC for secondary reaction as mentioned in your design document
    double C_EDC_secondary = C_EDC * 0.03;
    
    // Derivative of Arrhenius equation
    k2 = k2_0 * Math.exp(Ea2 / (R * T) * (1 / TR - 1 / T));
    
    // Add a rate limiter to prevent explosion
    double rate = this.k2 * C_EDC_secondary * C_NH3;
    
    // Limit maximum rate to prevent numerical instability
    double maxRate = 100.0; // Adjust this value based on your specific system
    if (Math.abs(rate) > maxRate) {
        rate = Math.signum(rate) * maxRate;
    }
    
    return rate;
}


    @Override
    public double calculateReactionEnthalpy(double T) {
        double H0_RX = Reaction.Ho_f[4] + Reaction.Ho_f[5] - Reaction.Ho_f[0] - Reaction.Ho_f[1];
        double Cp_prod = Reaction.CP298[2] * (T - Reaction.TR) + 1. / 2 * (Reaction.A[5] * (T - Reaction.TR))
                + 1. / 2 * (Reaction.B[5] * (Math.pow(T - Reaction.TR, 2))) + 1. / 3 * (Reaction.C[5] * Math.pow(T - Reaction.TR, 3))
                + 1. / 4 * (Reaction.D[5] * Math.pow(T - Reaction.TR, 4)) + (Reaction.E[5] * 1 / (T - Reaction.TR));
        double Cp_rea = Reaction.CP298[0] * (T - Reaction.TR) + 1. / 2 * (Reaction.A[1] * (T - Reaction.TR))
                + 1. / 2 * (Reaction.B[1] * (Math.pow(T - Reaction.TR, 2))) + 1. / 3 * (Reaction.C[1] * Math.pow(T - Reaction.TR, 3))
                + 1. / 4 * (Reaction.D[1] * Math.pow(T - Reaction.TR, 4)) + (Reaction.E[1] * 1 / (T - Reaction.TR));
        double delta_CP = Cp_prod - Cp_rea;
        double HRX;
        HRX = H0_RX + delta_CP;
        return HRX;
    }

}


