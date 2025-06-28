using SimpleHelpers.Indices;
using SimpleHelpers.Utilities.Pools;
using System.Buffers;
using System.Text;

namespace SimpleHelpers.MultiAlg.Helpers
{
    public class ContinuousInfo
    {
        public int NumLayer;

        public List<ContinuousLayer> Layers;

        public ContinuousInfo(ReadOnlySpan<nint> head, bool shouldSort = true)
        {
            if (head.Length == 0)
            {
                NumLayer = 0;
                Layers = [];
                return;
            }
            else
            {
                int rank = head.Length / 2;
                Span<SingleIndice> dims = stackalloc SingleIndice[rank];
                for (int i = 0; i < rank; i++)
                {
                    dims[i] = new SingleIndice(head[i], head[rank + i]);
                }
                NumLayer = rank;
                Layers = new List<ContinuousLayer>(NumLayer * 2);
                var buffer = ArrayPool<SingleIndice>.Shared.Rent(NumLayer);
                buffer.AsSpan(0, rank).Clear();
                for (int i = 0; i < rank; i++)
                {
                    buffer[i] = new(head[i], head[rank + i]);
                }
                if (shouldSort)
                {
                    // sort by stride descending
                    Array.Sort(buffer, (a, b) => b.Stride.CompareTo(a.Stride));
                }
                int dimIndex = 0;
                int layerIndex = 0;
                int headIndex = 0;
                var dim = buffer[dimIndex];
                var layer = ContinuousLayer.HeadFromDiminfo(dim, layerIndex);
                var dimOrigIndex = 0;
                Layers.Add(layer);//Placehold
                Layers.Add(ContinuousLayer.FromDiminfo(dim, dimOrigIndex));
                dimIndex++;
                while (dimIndex < NumLayer)
                {
                    dim = buffer[dimIndex];

                    if (dim.Stride == layer.Size)
                    {
                        layer.Length *= dim.Length;
                        layer.Size *= dim.Length;
                    }
                    else // should create new layer
                    {
                        Layers[headIndex] = layer;//write back
                        headIndex = Layers.Count;// new head index
                        layerIndex++;
                        layer = ContinuousLayer.HeadFromDiminfo(dim, layerIndex);
                        Layers.Add(layer);
                    }
                    dimOrigIndex++;
                    Layers.Add(ContinuousLayer.FromDiminfo(dim, dimOrigIndex));
                    dimIndex++;
                }
                Layers[headIndex] = layer;//write back
                NumLayer = layerIndex + 1;
                buffer.AsSpan(0, NumLayer).Clear();
                ArrayPool<SingleIndice>.Shared.Return(buffer);
            }

        }

        public void ToString(StringBuilder sb)
        {
            var head = default(ContinuousLayer);
            bool first = true;
            foreach (var layer in Layers)
            {
                if (layer.IsHead)
                {
                    if (!first)
                    {
                        sb.Append($"totalStep={head.Stride}  ")
                        .Append($"totalLength={head.Length}  ")
                        .Append($"totalTotal={head.Size}")
                        .AppendLine();
                    }
                    head = layer;
                    sb.AppendLine($"Layer {layer.Index}: ");
                }
                else
                {
                    sb.Append($">>DimIndex: {layer.Index} -> ")
                        .Append($"step={layer.Stride}  ")
                        .Append($"length={layer.Length}  ")
                        .Append($"total={layer.Size}")
                        .AppendLine();
                }
            }
            sb.Append($"totalStep={head.Stride}  ")
            .Append($"totalLength={head.Length}  ")
            .Append($"totalSize={head.Size}");
        }

        public override string ToString()
        {
            using (StringBuilderPool.Borrow(out var sb))
            {
                ToString(sb);
                return sb.ToString();
            }
        }
    }
}
