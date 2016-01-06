//! [snippet1]

#include <clustering_1.hpp>

using namespace std;

/********************************************************************************/

// We define some constant strings for names of command line parameters
static const char* STR_FOO = "-foo";

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
clustering_1::clustering_1 ()  : Tool ("clustering_1")
{
    // We add some custom arguments for command line interface
    getParser()->push_front (new OptionOneParam (STR_FOO, "my option",  false, "1"));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void clustering_1::execute ()
{
    // We can do here anything we want.
    // For further information about the Tool class, please have a look
    // on the ToyTool snippet  (http://gatb-core.gforge.inria.fr/snippets_tools.html)

    // We gather some statistics.
    getInfo()->add (1, "input");
    getInfo()->add (2, STR_FOO,  "%d",  getInput()->getInt(STR_FOO));
    getInfo()->add (1, &LibraryInfo::getInfo());
}