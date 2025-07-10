//@HEADER
// ************************************************************************
//
//                     		       Pressio
//                             Copyright 2019
//    National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Chris Wentland (crwentl@sandia.gov)
//
// ************************************************************************
//@HEADER


/*
    These are largely copied from pressio-tutorials, simplified for my own reading comprehension
*/

#ifndef PRESSIODEMOAPPS_SCHWARZ_ROMUTILS_
#define PRESSIODEMOAPPS_SCHWARZ_ROMUTILS_

#include <string>
#include <fstream>
#include <iostream>
#include "Eigen/Dense"


namespace pschwarz {

void checkfile(const std::string & fileIn){
    std::ifstream infile(fileIn);
    if (infile.good() == 0) {
        throw std::runtime_error("Cannot find file: " + fileIn);
    }
}

template<class ScalarType>
auto read_matrix_from_binary(const std::string & fileName, int numColsToRead) {

    using matrix_type = Eigen::Matrix<ScalarType, -1, -1, Eigen::ColMajor>;
    using sc_t  = typename matrix_type::Scalar;

    checkfile(fileName);

    matrix_type M;

    std::ifstream fin(fileName, std::ios::in | std::ios::binary);
    fin.exceptions(std::ifstream::failbit | std::ifstream::badbit);

    // read 2 8-byte integer header, size matrix accordingly
    std::size_t rows = {};
    std::size_t cols = {};
    fin.read((char*) (&rows), sizeof(std::size_t));
    fin.read((char*) (&cols), sizeof(std::size_t));
    M.resize(rows, numColsToRead);

    // read matrix
    // assumed that binary file has same record length as intended matrix (float, double, etc.)
    const auto nBytes = rows * numColsToRead * sizeof(sc_t);
    fin.read( (char *) M.data(), nBytes );

    if (!fin){
        std::cerr << strerror(errno) << std::endl;
        throw std::runtime_error("ERROR READING binary file");
    }
    else{
        std::cout << fin.gcount() << " bytes read\n";
    }
    fin.close();
    return M;

}

template<class ScalarType>
auto read_vector_from_binary(const std::string & fileName) {
    auto Vmat = read_matrix_from_binary<ScalarType>(fileName, 1);

    using vector_type = Eigen::Vector<ScalarType, -1>;
    vector_type V(Vmat.rows());
    for (std::size_t i = 0; i < Vmat.rows(); ++i) {
        V(i) = Vmat(i, 0);
    }

    return V;
}

// TODO: adjust I/O naming conventions to reflect return type
template<class ScalarType>
auto read_vector_from_ascii(const std::string & fileName)
{

    checkfile(fileName);

    std::ifstream source;
    source.open(fileName, std::ios_base::in);
    std::string line, value;
    std::vector<ScalarType> v;
    while (std::getline(source, line) ) {
        std::istringstream in(line);
        in >> value;
        if (std::is_floating_point<ScalarType>::value) {
            v.push_back(std::atof(value.c_str()));
        }
        else {
            v.push_back(std::atoi(value.c_str()));
        }
    }
    source.close();

    return v;
}

template<class MatType>
void write_matrix_to_binary(const std::string & outfile, MatType outmat)
{
    std::ofstream out(outfile.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
    typename MatType::Index rows=outmat.rows(), cols=outmat.cols();
    out.write((char*) (&rows), sizeof(typename MatType::Index));
    out.write((char*) (&cols), sizeof(typename MatType::Index));
    out.write((char*) outmat.data(), rows*cols*sizeof(typename MatType::Scalar) );
    out.close();
}

// pretty much everything below this is ripped directly from pressio-tutorials
template <typename T = int32_t>
auto create_cell_gids_vector_and_fill_from_ascii(const std::string & fileName)
{
    using vector_type = Eigen::Matrix<T, -1, 1>;

    const auto v = read_vector_from_ascii<T>(fileName);
    vector_type result(v.size());
    for (std::size_t i = 0; i < v.size(); ++i) {
        result[i] = v[i];
    }
    return result;
}

template<class IntType>
IntType find_index(const Eigen::Matrix<IntType, -1, 1> & vector, IntType target)
{
    for (std::size_t j = 0; j < vector.size(); ++j) {
        if (vector[j] == target) {
            return j;
        }
    }

    return std::numeric_limits<IntType>::max();
}

// class required to pass to LSPG hyper-reduction problem
// provides a sort of BLAS interface for vector and matrix addition
// also supplies mapping from sample mesh indices to stencil mesh indices
template<class ScalarType>
struct HypRedUpdater
{
    using vec_operand_type = Eigen::Matrix<ScalarType, -1, 1>;
    using mat_ll_operand_type = Eigen::Matrix<ScalarType, -1, -1>;
    std::vector<int> indices_;
    int numDofsPerCell_ = {};

    explicit HypRedUpdater(const int numDofsPerCell, const std::string & stfile, const std::string & safile)
        : numDofsPerCell_(numDofsPerCell)
    {
        const auto stencilMeshGids = create_cell_gids_vector_and_fill_from_ascii(stfile);
        const auto sampleMeshGids  = create_cell_gids_vector_and_fill_from_ascii(safile);

        indices_.resize(sampleMeshGids.size());
        for (std::size_t i = 0; i < indices_.size(); ++i) {
            const auto index = find_index<int>(stencilMeshGids, sampleMeshGids[i]);
            assert(index != std::numeric_limits<int>::max());
            indices_[i] = index;
        }
    }

    void updateSampleMeshOperandWithStencilMeshOne(
        vec_operand_type & a,
        ScalarType alpha,
        const vec_operand_type & b,
        ScalarType beta) const
    {
        for (std::size_t i = 0; i < indices_.size(); ++i) {
            const std::size_t r = i * numDofsPerCell_;
            const std::size_t g = indices_[i] * numDofsPerCell_;
            for (std::size_t k = 0; k < numDofsPerCell_; ++k) {
                a(r+k) = alpha * a(r+k) + beta * b(g+k);
            }
        }
    }

    void updateSampleMeshOperandWithStencilMeshOne(
        mat_ll_operand_type & a,
        ScalarType alpha,
        const mat_ll_operand_type & b,
        ScalarType beta) const
    {
        for (std::size_t j = 0; j < b.cols(); ++j) {
            for (std::size_t i = 0; i < indices_.size(); ++i) {
                const std::size_t r = i * numDofsPerCell_;
                const std::size_t g = indices_[i] * numDofsPerCell_;
                for (std::size_t k = 0; k < numDofsPerCell_; ++k) {
                    a(r+k, j) = alpha * a(r+k, j) + beta * b(g+k, j);
                }
            }
        }
    }
};

template<class mesh_t>
auto create_hyper_updater(
    const int numDofsPerCell,
    const std::string & stfile,
    const std::string & safile)
{
    using scalar_type = typename mesh_t::scalar_t;

    checkfile(stfile);
    checkfile(safile);

    using return_type = HypRedUpdater<scalar_type>;
    return return_type(numDofsPerCell, stfile, safile);

}

// extract stencil mesh values from full-order matrix
// rows are unrolled state vector, columns are basis/snapshot/etc. vectors
// OperandType should really be an Eigen::Matrix
template<class OperandType, class CellGidsVectorType>
auto reduce_matrix_on_stencil_mesh(
    const OperandType & operand,
    const CellGidsVectorType & stencilMeshGids,
    const int numDofsPerCell)
{

    const auto totStencilDofs = stencilMeshGids.size() * numDofsPerCell;
    OperandType result(totStencilDofs, operand.cols());
    for (int i = 0; i < stencilMeshGids.size(); ++i) {
        for (int k = 0; k < numDofsPerCell; ++k){
            const int row = i * numDofsPerCell + k;
            const int ind = stencilMeshGids[i] * numDofsPerCell + k;
            for (int j = 0; j < operand.cols(); ++j){
                result(row, j) = operand(ind, j);
            }
        }
    }
    return result;
}

// extract stencil mesh values from full-order vector
// OperandType should really be an Eigen::Vector
template<class OperandType, class CellGidsVectorType>
auto reduce_vector_on_stencil_mesh(
    const OperandType & operand,
    const CellGidsVectorType & stencilMeshGids,
    const int numDofsPerCell)
{

    const auto totStencilDofs = stencilMeshGids.size() * numDofsPerCell;
    OperandType result(totStencilDofs);
    for (int i = 0; i < stencilMeshGids.size(); ++i) {
        for (int k = 0; k < numDofsPerCell; ++k) {
            const int row = i * numDofsPerCell + k;
            const int ind = stencilMeshGids[i] * numDofsPerCell + k;
            result(row) = operand(ind);
        }
    }
    return result;
}

// Gappy POD weighting operator passed to create_gauss_newton_solver()
// Computes operator [ S * Psi ]^+,
//      where S is the sample matrix and Psi is the gappy POD regressor matrix of choice
// `nmodes` can be different than the number of modes in the trial basis
template<class scalar_t>
class Weigher {

    using matrix_type = Eigen::Matrix<scalar_t, -1, -1, Eigen::ColMajor>;

public:

    Weigher(
        const std::string & weigher_type,
        const std::string & basisfile,
        const std::string & samplefile,
        const int nmodes,
        const int numDofsPerCell)
    {

        m_weigher_type = weigher_type;

        const auto sampleGids = pschwarz::create_cell_gids_vector_and_fill_from_ascii(samplefile);

        if (weigher_type == "identity") {
            leading_dim = sampleGids.rows();
        }
        else if (weigher_type == "gappy_pod") {
            leading_dim = nmodes;

            // compute Z * Phi
            auto basis_gpod = pschwarz::read_matrix_from_binary<scalar_t>(basisfile, nmodes);
            auto basis_sample = pschwarz::reduce_matrix_on_stencil_mesh(basis_gpod, sampleGids, numDofsPerCell);

            // size matrices
            std::size_t numsamps = sampleGids.rows();
            m_gpod_operator.resize(nmodes, numsamps);

            // compute A = pinv(Z * Phi)
            m_gpod_operator = basis_sample.completeOrthogonalDecomposition().pseudoInverse();

        }
        else {
            throw std::runtime_error("Invalid weigher_type: " + weigher_type);
        }
    }

    int leadingDim() const { return leading_dim; }

    // operator on residual
    void operator()(const Eigen::Matrix<scalar_t, -1, 1> & operand,
                    Eigen::Matrix<scalar_t, -1, 1> & result) const
    {
        if (m_weigher_type == "identity") {
            // copy
            result = operand;
        }
        else if (m_weigher_type == "gappy_pod") {
            // multiply weighting operator
            pressio::ops::product(
                ::pressio::nontranspose(),
                1., m_gpod_operator, operand,
                0., result
            );
        }
    }

    // operator on Jacobian
    void operator()(const Eigen::Matrix<scalar_t, -1, -1> & operand,
                    Eigen::Matrix<scalar_t, -1, -1> & result) const
    {
        if (m_weigher_type == "identity") {
            // copy
            result = operand;
        }
        else if (m_weigher_type == "gappy_pod") {
            // multiply weighting operator
            pressio::ops::product(
                ::pressio::nontranspose(), ::pressio::nontranspose(),
                1., m_gpod_operator, operand,
                0., result
            );
        }
    }

private:

    int leading_dim;
    std::string m_weigher_type;
    matrix_type m_gpod_operator;

};

}

#endif
