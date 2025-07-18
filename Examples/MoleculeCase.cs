using CintSharp.DataStructures;

namespace Examples
{
    public class MoleculeCase(string name, Atom[] atoms, 
        string basisName,
        double rhfEnergy, double mp2Energy, 
        double ccdEnergy, double ccsdEnergy,
        string auxBasisName,
        double dfRHFEnergy)
    {
        public string Name { get; } = name;

        public Atom[] Atoms { get; } = atoms;

        public string BasisName { get; } = basisName;

        public double RhfEnergy { get; } = rhfEnergy;

        public double Mp2Energy { get; } = mp2Energy;

        public double CCDEnergy { get; } = ccdEnergy;

        public double CCSDEnergy { get; } = ccsdEnergy;

        public string AuxBasisName { get; } = auxBasisName;

        public double DfRHFEnergy { get; } = dfRHFEnergy;

        public static MoleculeCase CH3OH => new(
            "CH3OH",
            [
                new("C", 0.664617, 0.033231, 0.000000),
                new("H", 1.338596, -1.873145, 0.000000),
                new("H", 1.338631, 0.986406, -1.650963),
                new("H", -1.357391, 0.033256, 0.000000),
                new("O", 1.565402, 1.307100, 2.206427),
                new("H", 0.961143, 3.017646, 2.207215),
            ],
            "sto-3g",
            -113.54406552303,
            -113.629390527013,
            -113.54406552303 - 0.114891844396,
            -113.54406552303 - 0.115430350863,
            "def2-universal-jkfit",
            -113.54426993866
        );

        public static MoleculeCase CH3OH_st_ => new(
            "CH3OH",
            [
                new("C", 0.664617, 0.033231, 0.000000),
                new("H", 1.338596, -1.873145, 0.000000),
                new("H", 1.338631, 0.986406, -1.650963),
                new("H", -1.357391, 0.033256, 0.000000),
                new("O", 1.565402, 1.307100, 2.206427),
                new("H", 0.961143, 3.017646, 2.207215),
            ],
            "6-31g_st_",
            -115.032381671115,
            -115.34376394189,
            -115.032381671115 - 0.3314715986909997,
            -115.032381671115 - 0.3332518402504928,
            "def2-universal-jkfit",
            0.0
        );
    }
}
