using SimpleHelpers.LinearAlg;
using SimpleHelpers.MultiAlg;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics.Tensors;
using System.Text;
using System.Threading.Tasks;

namespace SimpleHelpers
{
    public static class Other2SNTensor
    {
        public static Tensor<double> AsSNTensor(this Matrix mat)
        {
            if(mat.ColStride > mat.RowStride)
            {
                return System.Numerics.Tensors.Tensor.Create(mat.Data, (int)mat.Offset,
                    [mat.Cols, mat.Rows], [mat.ColStride, mat.RowStride]).
                    PermuteDimensions([1, 0]);
            }
            else
            {
                return System.Numerics.Tensors.Tensor.Create(mat.Data, (int)mat.Offset,
                    [mat.Rows, mat.Cols], [mat.RowStride, mat.ColStride]);
            }
        }

        public static Vector AsVector(this System.Numerics.Tensors.Tensor<double> tensor)
        {
            if (tensor.Rank != 1)
                throw new ArgumentException("Tensor must be 1-dimensional to convert to Vector.");
            Vector result = Vector.Create(tensor.Lengths[0], uninited: true);
            for (int i = 0; i < tensor.Lengths[0]; i++)
            {
                result[i] = tensor[i];
            }
            return result;
        }

        public static System.Numerics.Tensors.Tensor<double> AsSNTensor(
            this Vector vec)
        {
            return System.Numerics.Tensors.Tensor.Create
                (vec.Data, (int)vec.Offset, [vec.Length], [vec.Stride]);
        }

        public static NDArray AsNDArray(this Matrix mat)
        {
            return new NDArray(mat.Data, mat.Offset,
                [mat.Rows, mat.Cols], [mat.RowStride, mat.ColStride]);
        }

        public static NDArray AsNDArray(this Vector vec)
        {
            return new NDArray(vec.Data, vec.Offset,
                [vec.Length], [vec.Stride]);
        }

        public static Matrix AsMatrix(this NDArray array)
        {
            if (array.Rank != 2)
                throw new ArgumentException("NDArray must be 2-dimensional to convert to Matrix.");
            return new Matrix(array.Data, (int)array.Offset,
                array.Lengths[0], array.Lengths[1],
                array.Strides[0], array.Strides[1]);
        }

        public static Vector AsVector(this NDArray array)
        {
            if (array.Rank != 1)
                throw new ArgumentException("NDArray must be 1-dimensional to convert to Vector.");
            return new Vector(array.Data, (int)array.Offset,
                array.Lengths[0], array.Strides[0]);
        }
    }
}
