#include "csrdictionary.h"
#include <algorithm>

icy::CSRDictionary::CSRDictionary()
{
    int initial_size = 50000;
    int reserve = 100000;

    rows_staticNeighbors.reserve(reserve);
    rows_staticNeighbors.resize(initial_size);

    rows_dynamicNeighbors.reserve(reserve);
    rows_dynamicNeighbors.resize(initial_size);

    rows_pcsr.reserve(reserve);
    rows_pcsr.resize(initial_size);

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
        sortedAllNeighbors.reserve(30);

        std::set<int> &staticNeighbors = rows_staticNeighbors[i];
        std::set<int> &dynamicNeighbors = rows_dynamicNeighbors[i];

        std::set_union(staticNeighbors.begin(),staticNeighbors.end(),
                       dynamicNeighbors.begin(),dynamicNeighbors.end(),
                       std::back_inserter(sortedAllNeighbors));
    }

    // count non-zero entries
    nnz = 0;
    for(auto const &vec : rows_sortedAllNeighbors) nnz+=vec.size();


}

int icy::CSRDictionary::offset(int row, int column)
{
    return rows_pcsr[row][column];
}

void icy::CSRDictionary::Assert(bool nonSymmetric)
{

}
