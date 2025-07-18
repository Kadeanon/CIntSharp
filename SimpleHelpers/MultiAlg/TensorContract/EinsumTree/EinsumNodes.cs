using SimpleHelpers.Indices;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Numerics;
using System.Numerics.Tensors;
using System.Text;
using System.Threading.Tasks;

namespace SimpleHelpers.MultiAlg.TensorContract.EinsumTree
{
    public abstract class EinsumNode(string name)
    {
        public string Name { get; set; } = name;

        protected internal abstract NDArray Kernel(ReadOnlySpan<NDArray> inputs);

        protected internal abstract void Check(ReadOnlySpan<NDArray> inputs);
    }

    public class InputNode(string symbol, int index,
        SingleIndice[] indices) : EinsumNode(symbol)
    {
        public SingleIndice[] Indices { get; set; } = indices;

        public int Index { get; set; } = index;

        public override string ToString()
        {
            StringBuilder sb = new();
            var rank = Name.Length;
            sb.Append(Name);
            sb.Append(' ');
            for (int i = 0; i < rank; i++)
            {
                var c = Name[i];
                var ind = Indices[i];
                sb.Append($"{c}[{ind.Length}:{ind.Stride}]");
            }
            return sb.ToString();
        }

        protected internal override NDArray Kernel(ReadOnlySpan<NDArray> inputs)
        {
            return inputs[Index];
        }

        protected internal override void Check(ReadOnlySpan<NDArray> inputs)
        {
            if (inputs.Length <= Index)
            {
                throw new ArgumentOutOfRangeException(nameof(inputs),
                    $"Input index {Index} is out of range for inputs of length {inputs.Length}.");
            }
        }
    }

    [DebuggerDisplay($"{{{nameof(Expression)}}}")]
    public class ContractedNode(string symbol,
        EinsumNode left, EinsumNode right) : EinsumNode(symbol)
    {
        public EinsumNode NodeLeft { get; set; } = left;

        public EinsumNode NodeRight { get; set; } = right;

        public string Expression
            => $"{NodeLeft.Name},{NodeRight.Name}->{Name}";

        protected internal override NDArray Kernel(ReadOnlySpan<NDArray> inputs)
        {
            var left = NodeLeft.Kernel(inputs);
            var right = NodeRight.Kernel(inputs);
            return NDArray.Contract(Expression,
                left, right);
        }

        protected internal override void Check(ReadOnlySpan<NDArray> inputs)
        {
            NodeLeft.Check(inputs);
            NodeRight.Check(inputs);
        }
    }

    public class OutputNode : EinsumNode
    {
        public EinsumNode Child { get; set; }

        public OutputNode(string symbol, EinsumNode child) : base(symbol)
        {
            if (child.Name != symbol)
            {
                if (child is ContractedNode)
                {
                    child.Name = symbol;
                }
            }
            Child = child;
        }

        public override string ToString()
        {
            return $"{Child.Name}->{Name}";
        }

        protected internal override NDArray Kernel(ReadOnlySpan<NDArray> inputs)
        {
            return Child.Kernel(inputs);
        }

        public NDArray Invoke(params ReadOnlySpan<NDArray> inputs)
        {
            if (inputs.Length == 0)
            {
                throw new ArgumentException("No inputs provided for the output node.");
            }
            Check(inputs);
            var result = Kernel(inputs);
            return result;
        }

        protected internal override void Check(ReadOnlySpan<NDArray> inputs)
        {
            Child.Check(inputs);
        }
    }
}
