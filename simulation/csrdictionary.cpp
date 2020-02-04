#include <algorithm>
#include <stdexcept>
#include <iostream>
#include "csrdictionary.h"

icy::CSRDictionary::CSRDictionary()
{
    // these are typical values for medium-sized FEM geometry
    int initial_size = 50000;
    int reserve = 100000;

    rows_staticNeighbors.reserve(reserve);
    rows_staticNeighbors.resize(initial_size);

    rows_dynamicNeighbors.reserve(reserve);
    rows_dynamicNeighbors.resize(initial_size);

    rows_pcsr.reserve(reserve);
    rows_sortedAllNeighbors.reserve(reserve);

    csr_cols = csr_rows = nullptr;
}

void icy::CSRDictionary::updateMaxRowIndex(int rowIndex)
{
    if (rowIndex <= maxRowIndex) return;
    maxRowIndex = rowIndex; // update the row count
    rows_staticNeighbors.resize(maxRowIndex+1);
    rows_dynamicNeighbors.resize(maxRowIndex+1);
}

void icy::CSRDictionary::ClearStatic()
{
    maxRowIndex = -1;
    for(auto &rowSet : rows_staticNeighbors) rowSet.clear();
}

void icy::CSRDictionary::ClearDynamic()
{
    for(auto &rowSet : rows_dynamicNeighbors) rowSet.clear();
}

void icy::CSRDictionary::AddDynamic(int row, int column)
{
    rows_dynamicNeighbors[row].insert(column);
}

void icy::CSRDictionary::AddStatic(int row, int column)
{
    if (row > maxRowIndex) updateMaxRowIndex(row);
    rows_staticNeighbors[row].insert(column);
}

void icy::CSRDictionary::CreateStructure()
{
    // matrix size
    N = maxRowIndex + 1;

    rows_sortedAllNeighbors.resize(N);
    rows_pcsr.resize(N);

    // unite and sort neighbors
#pragma omp parallel for
    for(int i=0;i<N;i++)
    {
        rows_pcsr[i].clear();
        std::vector<int> & sortedAllNeighbors = rows_sortedAllNeighbors[i];
        sortedAllNeighbors.clear();
        sortedAllNeighbors.reserve(30); // typical value (maybe overkill) for tetrahedral meshes

        std::set<int> &staticNeighbors = rows_staticNeighbors[i];
        std::set<int> &dynamicNeighbors = rows_dynamicNeighbors[i];

//        if(dynamicNeighbors.size() != 0) std::cout << "row " << i << "; st " << staticNeighbors.size() << "; dy " << dynamicNeighbors.size();

        std::set_union(staticNeighbors.begin(),staticNeighbors.end(),
                       dynamicNeighbors.begin(),dynamicNeighbors.end(),
                       std::back_inserter(sortedAllNeighbors));

//        if(dynamicNeighbors.size() != 0) std::cout << "; union " << sortedAllNeighbors.size() << std::endl;
    }

    // count non-zero entries
    nnz = 0;
    for(auto const &vec : rows_sortedAllNeighbors) nnz+=vec.size();

    // allocate structure arrays
    if(csr_rows_size < N+1) {
        csr_rows_size = N+1;
        delete csr_rows;
        csr_rows = new int[csr_rows_size];
    }

    if(csr_cols_size < nnz) {
        csr_cols_size = nnz*2;
        delete csr_cols;
        csr_cols = new int[csr_cols_size];
    }

    csr_rows[N] = nnz;

    // enumerate entries
    for(int i=0,count=0;i<N;i++)
    {
        csr_rows[i] = count;
        for(int const &local_column : rows_sortedAllNeighbors[i])
        {
            rows_pcsr[i][local_column] = count;
            csr_cols[count] = local_column;
            count++;
        }
    }
}

int icy::CSRDictionary::offset(int row, int column)
{
    return rows_pcsr[row][column];
}

void icy::CSRDictionary::Assert()
{
    if(csr_rows[0] != 0) throw std::runtime_error("rows[0] != 0");
    if(csr_rows[N] != nnz) throw std::runtime_error("rows[N] != nnz");
    for (int i = 1; i < N + 1; i++)
        if (csr_rows[i] <= csr_rows[i - 1]) throw std::runtime_error("rows[i] is not increasing");

    // verify columns array, upper triangular
    for (int i = 0; i < N; i++)
    {
        if (csr_cols[csr_rows[i]] != i) throw std::runtime_error("structure not UT");
        for (int j = csr_rows[i]; j < csr_rows[i + 1] - 1; j++)
            if (csr_cols[j + 1] <= csr_cols[j]) throw std::runtime_error("cols in same row not increasing");
    }
}
