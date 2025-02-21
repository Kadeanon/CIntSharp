using CintSharp.DataStructures;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics.Tensors;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace CintSharp.Intor
{
    public static partial class IntorUtils
    {
        public static readonly Dictionary<string, (int comp4Scalar, int comp4Spinor)> intorDict = new()
        {
            {"int1e_ovlp", (1, 1)},
            {"int1e_nuc", (1, 1)},
            {"int1e_kin", (1, 1)},
            {"int1e_ia01p", (3, 3)},
            {"int1e_giao_irjxp", (3, 3)},
            {"int1e_cg_irxp", (3, 3)},
            {"int1e_giao_a11part", (9, 9)},
            {"int1e_cg_a11part", (9, 9)},
            {"int1e_a01gp", (9, 9)},
            {"int1e_igkin", (3, 3)},
            {"int1e_igovlp", (3, 3)},
            {"int1e_ignuc", (3, 3)},
            {"int1e_pnucp", (1, 1)},
            {"int1e_z", (1, 1)},
            {"int1e_zz", (1, 1)},
            {"int1e_r", (3, 3)},
            {"int1e_r2", (1, 1)},
            {"int1e_r4", (1, 1)},
            {"int1e_rr", (9, 9)},
            {"int1e_rrr", (27, 27)},
            {"int1e_rrrr", (81, 81)},
            {"int1e_z_origj", (1, 1)},
            {"int1e_zz_origj", (1, 1)},
            {"int1e_r_origj", (3, 3)},
            {"int1e_rr_origj", (9, 9)},
            {"int1e_r2_origj", (1, 1)},
            {"int1e_r4_origj", (1, 1)},
            {"int1e_p4", (1, 1)},
            {"int1e_prinvp", (1, 1)},
            {"int1e_prinvxp", (3, 3)},
            {"int1e_pnucxp", (3, 3)},
            {"int1e_irp", (9, 9)},
            {"int1e_irrp", (27, 27)},
            {"int1e_irpr", (27, 27)},
            {"int1e_ggovlp", (9, 9)},
            {"int1e_ggkin", (9, 9)},
            {"int1e_ggnuc", (9, 9)},
            {"int1e_grjxp", (9, 9)},
            {"int2e", (1, 1)},
            {"int2e_ig1", (3, 3)},
            {"int2e_gg1", (9, 9)},
            {"int2e_g1g2", (9, 9)},
            {"int2e_ip1v_rc1", (9, 9)},
            {"int2e_ip1v_r1", (9, 9)},
            {"int2e_ipvg1_xp1", (9, 9)},
            {"int2e_ipvg2_xp1", (9, 9)},
            {"int2e_p1vxp1", (3, 3)},
            {"int1e_inuc_rcxp", (3, 3)},
            {"int1e_inuc_rxp", (3, 3)},
            {"int1e_sigma", (12, 3)},
            {"int1e_spsigmasp", (12, 3)},
            {"int1e_srsr", (4, 1)},
            {"int1e_sr", (4, 1)},
            {"int1e_srsp", (4, 1)},
            {"int1e_spsp", (4, 1)},
            {"int1e_sp", (4, 1)},
            {"int1e_spnucsp", (4, 1)},
            {"int1e_sprinvsp", (4, 1)},
            {"int1e_srnucsr", (4, 1)},
            {"int1e_sprsp", (12, 3)},
            {"int1e_govlp", (3, 3)},
            {"int1e_gnuc", (3, 3)},
            {"int1e_cg_sa10sa01", (36, 9)},
            {"int1e_cg_sa10sp", (12, 3)},
            {"int1e_cg_sa10nucsp", (12, 3)},
            {"int1e_giao_sa10sa01", (36, 9)},
            {"int1e_giao_sa10sp", (12, 3)},
            {"int1e_giao_sa10nucsp", (12, 3)},
            {"int1e_sa01sp", (12, 3)},
            {"int1e_spgsp", (12, 3)},
            {"int1e_spgnucsp", (12, 3)},
            {"int1e_spgsa01", (36, 9)},
            {"int2e_spsp1", (4, 1)},
            {"int2e_spsp1spsp2", (16, 1)},
            {"int2e_srsr1", (4, 1)},
            {"int2e_srsr1srsr2", (16, 1)},
            {"int2e_cg_sa10sp1", (12, 3)},
            {"int2e_cg_sa10sp1spsp2", (48, 3)},
            {"int2e_giao_sa10sp1", (12, 3)},
            {"int2e_giao_sa10sp1spsp2", (48, 3)},
            {"int2e_g1", (12, 3)},
            {"int2e_spgsp1", (12, 3)},
            {"int2e_g1spsp2", (12, 3)},
            {"int2e_spgsp1spsp2", (48, 3)},
            {"int2e_pp1", (1, 1)},
            {"int2e_pp2", (1, 1)},
            {"int2e_pp1pp2", (1, 1)},
            {"int1e_spspsp", (4, 1)},
            {"int1e_spnuc", (4, 1)},
            {"int2e_spv1", (4, 1)},
            {"int2e_vsp1", (4, 1)},
            {"int2e_spsp2", (4, 1)},
            {"int2e_spv1spv2", (16, 1)},
            {"int2e_vsp1spv2", (16, 1)},
            {"int2e_spv1vsp2", (16, 1)},
            {"int2e_vsp1vsp2", (16, 1)},
            {"int2e_spv1spsp2", (16, 1)},
            {"int2e_vsp1spsp2", (16, 1)},
            {"int1e_ipovlp", (3, 3)},
            {"int1e_ipkin", (3, 3)},
            {"int1e_ipnuc", (3, 3)},
            {"int1e_iprinv", (3, 3)},
            {"int1e_rinv", (1, 1)},
            {"int1e_ipspnucsp", (12, 3)},
            {"int1e_ipsprinvsp", (12, 3)},
            {"int1e_ippnucp", (3, 3)},
            {"int1e_ipprinvp", (3, 3)},
            {"int2e_ip1", (3, 3)},
            {"int2e_ip2", (3, 3)},
            {"int2e_ipspsp1", (12, 3)},
            {"int2e_ip1spsp2", (12, 3)},
            {"int2e_ipspsp1spsp2", (48, 3)},
            {"int2e_ipsrsr1", (12, 3)},
            {"int2e_ip1srsr2", (12, 3)},
            {"int2e_ipsrsr1srsr2", (48, 3)},
            {"int2e_ssp1ssp2", (16, 1)},
            {"int2e_ssp1sps2", (16, 1)},
            {"int2e_sps1ssp2", (16, 1)},
            {"int2e_sps1sps2", (16, 1)},
            {"int2e_cg_ssa10ssp2", (48, 3)},
            {"int2e_giao_ssa10ssp2", (18, 3)},
            {"int2e_gssp1ssp2", (18, 3)},
            {"int2e_gauge_r1_ssp1ssp2", (0, 1)},
            {"int2e_gauge_r1_ssp1sps2", (0, 1)},
            {"int2e_gauge_r1_sps1ssp2", (0, 1)},
            {"int2e_gauge_r1_sps1sps2", (0, 1)},
            {"int2e_gauge_r2_ssp1ssp2", (0, 1)},
            {"int2e_gauge_r2_ssp1sps2", (0, 1)},
            {"int2e_gauge_r2_sps1ssp2", (0, 1)},
            {"int2e_gauge_r2_sps1sps2", (0, 1)},
            {"int1e_ipipovlp", (9, 9)},
            {"int1e_ipovlpip", (9, 9)},
            {"int1e_ipipkin", (9, 9)},
            {"int1e_ipkinip", (9, 9)},
            {"int1e_ipipnuc", (9, 9)},
            {"int1e_ipnucip", (9, 9)},
            {"int1e_ipiprinv", (9, 9)},
            {"int1e_iprinvip", (9, 9)},
            {"int2e_ipip1", (9, 9)},
            {"int2e_ipvip1", (9, 9)},
            {"int2e_ip1ip2", (9, 9)},
            {"int1e_ipippnucp", (9, 9)},
            {"int1e_ippnucpip", (9, 9)},
            {"int1e_ipipprinvp", (9, 9)},
            {"int1e_ipprinvpip", (9, 9)},
            {"int1e_ipipspnucsp", (36, 9)},
            {"int1e_ipspnucspip", (36, 9)},
            {"int1e_ipipsprinvsp", (36, 9)},
            {"int1e_ipsprinvspip", (36, 9)},
            {"int3c2e", (1, 1)},
            {"int3c2e_ip1", (3, 3)},
            {"int3c2e_ip2", (3, 3)},
            {"int3c2e_pvp1", (1, 1)},
            {"int3c2e_pvxp1", (3, 3)},
            {"int2c2e_ip1", (3, 3)},
            {"int2c2e_ip2", (3, 3)},
            {"int3c2e_ig1", (3, 3)},
            {"int3c2e_spsp1", (4, 1)},
            {"int3c2e_ipspsp1", (12, 3)},
            {"int3c2e_spsp1ip2", (12, 3)},
            {"int3c2e_ipip1", (9, 9)},
            {"int3c2e_ipip2", (9, 9)},
            {"int3c2e_ipvip1", (9, 9)},
            {"int3c2e_ip1ip2", (9, 9)},
            {"int2c2e_ip1ip2", (9, 9)},
            {"int2c2e_ipip1", (9, 9)},
            {"int3c1e", (1, 1)},
            {"int3c1e_ip1", (3, 3)},
            {"int3c1e_p2", (1, 1)},
            {"int3c1e_iprinv", (3, 3)},
            {"int2c2e", (1, 1)},
            {"int2e_yp", (1, 1)},
            {"int2e_stg", (1, 1)},
            {"int2e_coulerf", (1, 1)},
            {"int1e_grids", (1, 1)},
            {"int1e_grids_ip", (3, 3)},
            {"int1e_grids_spvsp", (4, 1)},
            {"ECPscalar", (1, 0)},
            {"ECPscalar_ipnuc", (3, 0)},
            {"ECPscalar_iprinv", (3, 0)},
            {"ECPscalar_ignuc", (3, 0)},
            {"ECPscalar_iprinvip", (9, 0)},
            {"ECPso", (3, 1)}
        };

        public enum IntorType
        {
            Cartesian,
            Spheric,
            Spinor
        }

        public static int GetIntorComp(ref string intorName, out IntorType type)
        {
            Regex regex = SuffixRegex();
            Match match = regex.Match(intorName);
            if (!match.Success)
            {
                throw new ArgumentException($"Invalid intor name: {intorName}");
            }
            intorName = match.Groups[1].Value;
            type = match.Groups[2].Value switch
            {
                "_cart" => IntorType.Cartesian,
                "_spin" => IntorType.Spinor,
                _ => IntorType.Spheric
            };
            if (!intorDict.TryGetValue(intorName, out var value))
            {
                throw new ArgumentException($"Invalid intor name: {intorName}");
            }
            var (comp4Scalar, comp4Spinor) = value;
            return type == IntorType.Spinor ? comp4Spinor : comp4Scalar;
        }

        [GeneratedRegex(@"(.*?)(_cart|_sph|_spin)?$")]
        private static partial Regex SuffixRegex();

        #region Intors


        public static Tensor<double> InvokeIntor(this CIntEnvs envs, string intorName)
        {
            using var intor = envs.CreateIntor(intorName);
            return intor.Invoke();
        }

        public static Tensor<double> GetOvlp(this CIntEnvs envs) => envs.InvokeIntor("int1e_ovlp");

        public static Tensor<double> GetKin(this CIntEnvs envs) => envs.InvokeIntor("int1e_kin");

        public static Tensor<double> GetNuc(this CIntEnvs envs) => envs.InvokeIntor("int1e_nuc");

        public static Tensor<double> GetERI(this CIntEnvs envs) => envs.InvokeIntor("int2e");

        public static Tensor<double> GetHCore(this CIntEnvs envs) =>
            Tensor.Add(envs.GetKin().AsReadOnlyTensorSpan(), envs.GetNuc());

        #endregion
    }
}
