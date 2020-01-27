#ifndef CSRDICTIONARY_H
#define CSRDICTIONARY_H

//#include <utility>
#include <map>
#include <vector>
#include <set>

namespace icy {
class CSRDictionary;
}

class icy::CSRDictionary
{
public:

    std::vector<std::set<int>> rows_staticNeighbors;
    std::vector<std::set<int>> rows_allNeighbors;
    int *csr_rows;
    int *csr_cols;            // structure arrays of the sparse matrix
    int N = 0;   // number of variables
    int nnz = 0; // number of non-zero entries

    CSRDictionary();
    void updateMaxRowIndex(int rowIndex);   // expand "rows_" vectors as needed
    void ClearStatic();
    void ClearDynamic();
    void AddDynamic(int row, int column);   // insert dynamic non-zero element
    void AddStatic(int row, int column);
    void CreateStructure();
    int offset(int row, int column);    // returns the id of an i-j pair (pairs are consecutively numbered)
    void Assert(bool nonSymmetric); // check matrix consistency

private:
    int maxRowIndex = -1;   // maximum row index, zero-based, equal to resulting matrix height minus 1
    std::vector<std::map<int,int>> rows_pcsr;   // per row mappings between columns and offset in "values"

};

#endif // CSRDICTIONARY_H
