#include <metis.h>
#include <iostream>

using namespace std;
int main()
{

    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_DBGLVL] = METIS_DBG_INFO; // Enable debug info

    idx_t num_elements = 8;
    idx_t num_nodes = 10;
    idx_t eptr[] = {0, 3, 6, 9, 12, 15, 18, 21, 25};
    int eind[] = {0, 1, 2,
                  1, 3, 2,
                  1, 4, 3,
                  4, 5, 3,
                  5, 6, 3,
                  3, 6, 7,
                  2, 3, 7,
                  0, 2, 8, 9};

    idx_t ncommon = 2;   // Minimum number of nodes shared by two elements (2 since it's a 2D mesh)
    idx_t num_parts = 2; // Number of desired partitions
    idx_t objval{-1};

    idx_t *epart = new idx_t[num_elements];
    idx_t *npart = new idx_t[num_nodes];

    int ret_val = METIS_PartMeshDual(&num_elements, &num_nodes, eptr, eind,
                                     NULL, NULL, &ncommon, &num_parts, NULL, options,
                                     &objval, epart, npart);

    /*Print result*/

    cout << "Status: " << ret_val << endl;

    for (int i = 0; i < num_elements; ++i)
        cout << "Element " << i << " is in Partition " << epart[i] << endl;

    cout << endl;

    for (int i = 0; i < num_nodes; ++i)
        cout << "Node " << i << " is in Partition " << npart[i] << endl;

    cout << "objval: " << objval << endl;
}
