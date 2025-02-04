using CintSharp.DataStructures;
using CintSharp.Intor;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics.Tensors;
using System.Text;
using System.Threading.Tasks;

namespace CintSharp.Test
{
    [TestClass]
    public class TestCintor
    {

        [TestMethod]
        public void BuildCintor1e() 
        {
            List<Atom> list = [
                new Atom("O", 0.0, 0.0, 0.0), 
                new Atom("H", 0.0, 0.0, 1.0), 
                new Atom("H", 0.0, 1.0, 0.0)
                ];
            CIntEnvs cintor = CIntEnvs.Create(list, atm => "sto-3g");
            IntorBase ovlp = cintor.CreateIntor("int1e_ovlp");
            Tensor<double> result = ovlp.Invoke();
            double[][] standard = 
                [[1.00000002,0.23670394,.0         , .0         ,.0         ,0.14130744, 0.14130744],
                [0.23670394,1.00000003,.0         , .0         ,.0         ,0.77505773,0.77505773],
                [.0         ,.0         ,1.00000002,.0         ,.0         ,.0         ,.0         ],
                [.0         ,.0         ,.0         ,1.00000002,.0         ,.0         ,0.45204608],
                [.0         ,.0         ,.0         ,.0         ,1.00000002,0.45204608,.0         ],
                [0.14130744,0.77505773,.0         , .0         ,0.45204608,0.99999999,0.65439938],
                [0.14130744,0.77505773,.0         ,0.45204608,.0         ,0.65439938,0.99999999]];
            Console.WriteLine(result);

            for (int i = 0; i < result.Lengths[0]; i++)
            {
                for (int j = i + 1; j < result.Lengths[1]; j++)
                {
                    Assert.AreEqual(result[i, j], standard[i][j], 1e-7);
                }
            }
        }
    }
}
