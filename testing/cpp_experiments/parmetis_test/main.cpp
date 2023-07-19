
#include <iostream>
#include <metis.h>
#include <map>

using namespace std;

const map<rstatus_et, string> statuses = {{METIS_OK, "METIS_OK"},
                                          {METIS_ERROR_INPUT, "METIS_ERROR_INPUT"},
                                          {METIS_ERROR_MEMORY, "METIS_ERROR_MEMORY"},
                                          {METIS_ERROR, "METIS_ERROR"}};

int main(int argc, char *argv[])
{
    cout << "HALLA\n";

    idx_t eptr[] = {0, 3, 6, 9, 12};
    idx_t eind[] = {0, 1, 2,
                    1, 3, 4,
                    1, 4, 2,
                    2, 4, 5};

    idx_t ne = 4;
    idx_t nn = 6;
    idx_t ncommon = 2;
    idx_t nparts = 3;
    idx_t objval;
    idx_t *epart = new idx_t[ne];
    idx_t *npart = new idx_t[nn];

    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_DBGLVL] = METIS_DBG_INFO; // Enable debug info

    cerr << "starting metis\n";
    int status = METIS_PartMeshDual(&ne, &nn, eptr, eind, NULL, NULL,
                                    &ncommon, &nparts, NULL, options, &objval, epart, npart);

    cout << "Status: " << status << ", " << statuses.at(static_cast<rstatus_et>(status)) << endl;

    for (int i{0}; i < ne; i++)
    {
        cout << "epart[" << i << "] = " << epart[i] << endl;
    }
    cout << "\n\n";
    for (int i{0}; i < ne; i++)
    {
        cout << "npart[" << i << "] = " << npart[i] << endl;
    }

    cout << "objval " << objval << endl;
}