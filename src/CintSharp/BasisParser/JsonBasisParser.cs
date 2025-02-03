using CintSharp.DataStructures;
using CintSharp.Native.Libcint;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Text.Json;
using System.Text.Json.Nodes;
using System.Threading.Tasks;

namespace CintSharp.BasisParser
{
    public class JsonBasisParser(JsonBasisParser.BasisStreamSlover? CustomBasisStreamSlover = null) : IBasisParser
    {

        public delegate Stream BasisStreamSlover(string name);

        public JsonObject? CurrentJObj { get; private set; }

        public string BasisName { get; private set; } = "";

        public BasisStreamSlover? CustomBasisStreamSlover { get; set; }
            = CustomBasisStreamSlover;

        void LoadBasisFile(string basisName)
        {
            if(basisName == "")
            {
                throw new ArgumentException("Basis name cannot be empty");
            }

            if(BasisName == basisName)
            {
                return;
            }
            JsonNode? jsonNode;
            Stream stream = FindBasisFileStream(basisName);
            if (stream == null)
            {
                throw new ArgumentException($"Failed to open basis {basisName}: basis file not found.");
            }
            try
            {
                jsonNode = JsonNode.Parse(stream);
            }
            catch (Exception e)
            {
                throw new ArgumentException($"Failed to open basis {basisName}: {e.Message}");
            }

            jsonNode = jsonNode?["elements"]??throw new ArgumentException($"Failed to open basis {basisName}: elements not found");
            
            if (jsonNode is not JsonObject jsonObj)
            {
                throw new ArgumentException($"Failed to open basis {basisName}: elements is not an object");
            }
            CurrentJObj = jsonObj;
            BasisName = basisName;
        }

        Stream FindBasisFileStream(string name)
        {
            name = name.ToLower().Replace("*", "_st_");
            return
                Assembly.GetExecutingAssembly().GetManifestResourceStream($"CintSharp.Basis.{name}.json") ??
                CustomBasisStreamSlover?.Invoke(name) ??
                throw new ArgumentException($"Basis file {name} not found");
        }

        public Shells ParseBasis(BasisDict.BasisKey key, List<double> envs)
        {
            var basisName = key.BasisName;
            var atomNumber = key.AtomNumber;
            LoadBasisFile(basisName);
            var atomShells = CurrentJObj?[atomNumber.ToString()]["electron_shells"] as JsonArray
                ?? throw new ArgumentException
                ($"Failed to open basis {basisName}: electron_shells not found for atom {atomNumber}");
            var shells = new Shells();

            for (int layerIndex = 0; layerIndex < atomShells.Count; layerIndex++)//电子层,一个电子层可以包含libcint中的多个shell，例如sp层包括一个s轨道和一组p轨道（这组p轨道轨道共同构成一个libcint中的一个shell，简称shl）
            {
                var layer = atomShells[layerIndex];
                var mulShellCoefsJArray = layer["coefficients"] as JsonArray;//每个角动量对应一套系数，同时对应一个shl

                var expsJArray = layer["exponents"] as JsonArray;//同一电子层的不同shell共用一套指数

                double[] exps = new double[expsJArray.Count];
                for (int i = 0; i < exps.Length; i++)
                {
                    exps[i] = Convert.ToDouble(expsJArray[i].GetValue<string>());
                }

                for (int shellIndex = 0; shellIndex < mulShellCoefsJArray.Count; shellIndex++)//shell
                {
                    JsonArray coefsJArray = mulShellCoefsJArray[shellIndex] as JsonArray;
                    double[] coefs = new double[coefsJArray.Count];
                    JsonArray angJArray = layer["angular_momentum"] as JsonArray;
                    int angMomentum = (int)angJArray[shellIndex];//壳层劈裂时每个系数具有自己对应的角动量

                    for (int i = 0; i < coefs.Length; i++)
                    {
                        double coef = Convert.ToDouble(coefsJArray[i].GetValue<string>());
                        double normalized = LibcintHandler.CINTgto_norm(angMomentum, exps[i]);//当场归一化，方便后续使用
                        coefs[i] = coef * normalized;
                    }

                    Shell shell = new(exps, coefs, angMomentum, envs);
                    shells.Add(shell);
                }
            }

            return shells;
        }
    }
}
