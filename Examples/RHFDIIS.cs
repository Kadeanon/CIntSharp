using SimpleHelpers;
using SimpleHelpers.LinearAlg;
using SimpleHelpers.MultiAlg;

namespace Examples
{
    internal class RHFDIIS
    {
        const int max_diis = 6;
        readonly Matrix ovlp_;
        readonly Matrix error_ovlp_;
        readonly NDArray fock_list_;
        readonly Matrix error_list_;
        int current_size_;
        int current_index_;

        public RHFDIIS(RHF rhf)
        {
            nint nao = rhf.Ovlp.Rows;
            ovlp_ = rhf.Ovlp;
            error_ovlp_ = Matrix.Create(nao + 1, nao + 1);
            fock_list_ = NDArray.Create(max_diis, nao, nao);
            error_list_ = Matrix.Create(max_diis, nao * nao);
            current_size_ = 0;
            current_index_ = -1;
        }

        public Matrix Update(Matrix fock, Matrix density)
        {
            Matrix B = fock * density * ovlp_ - ovlp_ * density * fock;
            Vector error = B.Flatten();
            current_index_ = (current_index_ + 1) % max_diis;
            current_size_ = Math.Clamp(current_size_ + 1, 0, max_diis);
            fock_list_[current_index_, .., ..] = fock.AsNDArray();   
            error_list_[current_index_, ..] = error; 

            if (current_size_ <= 1)
                return fock;

            var workSize = current_size_ + 1;
            var errors = error_ovlp_[..workSize, ..workSize];
            for (int i = 0; i < current_size_; i++)
            {
                var error_i = error_list_[i, ..];
                var ovlp = error * error_i;
                errors[i, current_size_] = -1;
                errors[i, current_index_] = ovlp;
                errors[current_index_, i] = ovlp;
                errors[current_size_, i] = -1;
            }
            errors[current_size_, current_size_] = 0;
            var b = Vector.Create(workSize);
            b[current_size_] = -1;
            var lu = errors.LU();

            var coeffs = lu.Solve(b);
            coeffs[current_index_] -= 1;
            NDArray.Axpy(coeffs[..current_size_], 
                fock_list_[..current_size_, .., ..], fock.AsNDArray());

            return fock;
        }
    }
}
