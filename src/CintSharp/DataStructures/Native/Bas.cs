using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CIntSharp.DataStructures.Native
{
    public struct Bas
    {
        /// <summary>
        /// 0-based index of corresponding atom
        /// </summary>
        public int atomOf;

        /// <summary>
        /// Angular momentum
        /// </summary>
        public int angleOf;

        /// <summary>
        /// Number of primitive GTO
        /// </summary>
        public int numOfPrim;

        /// <summary>
        /// Number of contracted GTO
        ///</summary>
        public int numOfCont;

        /// <summary>
        /// Kappa for spinor GTO
        /// < 0 the basis ~ j = l + 1/2.
        /// > 0 the basis ~j = l - 1/2.
        /// = 0 the basis includes both j = l + 1/2 and j = l - 1 / 2
        /// </summary>
        public int kappaForSpinor;

        /// <summary>
        /// Ptr of exponents, 
        /// storing the offset of the primitive GTO exponents in env
        /// </summary>
        public int pointerOfExps;

        /// <summary>
        /// Ptr of coefficients,
        /// storing the offset of the column-major contracted GTO coefficients in env
        /// </summary>
        public int pointerOfCoeff;

        /// <summary>
        /// Unused
        /// </summary>
        public int reserved;
    }
}
