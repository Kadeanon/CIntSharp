using CintSharp.DataStructures;

namespace Examples
{
    internal class MoleculeCase(string name, Atom[] atoms, string basisName, 
        double rhfEnergy, double mp2Energy, 
        double ccdEnergy, double ccsdEnergy)
    {
        public string Name { get; } = name;

        public Atom[] Atoms { get; } = atoms;

        public string BasisName { get; } = basisName;

        public double RhfEnergy { get; } = rhfEnergy;

        public double Mp2Energy { get; } = mp2Energy;

        public double CCDEnergy { get; } = ccdEnergy;

        public double CCSDEnergy { get; } = ccsdEnergy;

        public static MoleculeCase CH3OH => new MoleculeCase(
            "CH3OH",
            new Atom[]
            {
                new("C", 0.664617, 0.033231, 0.000000),
                new("H", 1.338596, -1.873145, 0.000000),
                new("H", 1.338631, 0.986406, -1.650963),
                new("H", -1.357391, 0.033256, 0.000000),
                new("O", 1.565402, 1.307100, 2.206427),
                new("H", 0.961143, 3.017646, 2.207215),
            },
            "sto-3g",
            -113.54406552303,
            -113.629390517639,
            -113.54406552303 - 0.114891666546,
            -113.54406552303 - 0.115430350863
        );
    }
}
