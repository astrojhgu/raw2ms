#include <cassert>
#include <mscreate.hpp>
#include <ms/MeasurementSets.h>
#include <tables/DataMan/IncrementalStMan.h>
#include <tables/DataMan/StandardStMan.h>
#include <tables/DataMan/TiledColumnStMan.h>
#include <tables/DataMan/TiledStManAccessor.h>
#include <tables/Tables/SetupNewTab.h>
#include <tables/Tables/TableDesc.h>
#include <tables/Tables/ArrColDesc.h>
#include <tables/Tables/ScaColDesc.h>
#include <tables/Tables/TableRecord.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Cube.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Containers/Block.h>
#include <casa/Containers/Record.h>
#include <measures/Measures/MPosition.h>
#include <measures/Measures/MBaseline.h>
#include <measures/Measures/Muvw.h>
#include <measures/Measures/MeasTable.h>
#include <measures/Measures/Stokes.h>
#include <measures/Measures/MCBaseline.h>
#include <measures/Measures/MeasConvert.h>
#include <casa/Quanta/MVEpoch.h>
#include <casa/Quanta/MVDirection.h>
#include <casa/Quanta/MVPosition.h>
#include <casa/Quanta/MVBaseline.h>
#include <casa/OS/Time.h>
#include <casa/OS/SymLink.h>
#include <casa/BasicSL/Constants.h>
#include <casa/Utilities/Assert.h>
#include <casa/Exceptions/Error.h>
#include <casa/Arrays/Slicer.h>
#include <casa/Arrays/Slice.h>

using namespace casacore;


int main (int argc, char *argv[])
{
    if (argc != 2)
        {
            std::cerr << "Usage:" << argv[0] << " <input file>" << std::endl;
            std::cerr << "Input file: a 4-space-separated-columns file, with the columns of "
                         "antenna name, WCS84 x, y, and z for each antenna station"
                      << std::endl;
            return -1;
        }

    TableDesc td ("ANTENNA", TableDesc::New);
    td.addColumn (ScalarColumnDesc<double> ("DISH_DIAMETER", "Physical diameter of dish",
                                            "StandardStMan", "StandardStMan"));
    td.addColumn (ScalarColumnDesc<bool> ("FLAG_ROW", "Flag for this row", "StandardStMan", "StandardStMan"));
    td.addColumn (ArrayColumnDesc<double> ("PHASE_REFERENCE", "Beamformer phase reference position",
                                           "StandardStMan", "StandardStMan", 1));
    td.addColumn (ScalarColumnDesc<int> ("STATION_ID", "ID in STATION table", "StandardStMan", "StandardStMan"));
    td.addColumn (ScalarColumnDesc<String> ("MOUNT", "Mount type e.g. alt-az, equatorial, etc.",
                                            "StandardStMan", "StandardStMan"));
    td.addColumn (ScalarColumnDesc<String> ("NAME", "Antenna name", "StandardStMan", "StandardStMan"));
    td.addColumn (ArrayColumnDesc<double> ("OFFSET", "Axes offset of mount to FEED REFERENCE point",
                                           "StandardStMan", "StandardStMan", 1));
    td.addColumn (ArrayColumnDesc<double> ("POSITION", "Antenna X,Y,Z phase reference position",
                                           "StandardStMan", "StandardStMan", 1));
    td.addColumn (ScalarColumnDesc<String> ("STATION", "Station (antenna pad) name", "StandardStMan", "StandardStMan"));
    td.addColumn (ScalarColumnDesc<String> ("TYPE", "Antenna type (e.g., SPACE-BASED",
                                            "StandardStMan", "StandardStMan"));

    SetupNewTable snt ("ANTENNA", td, Table::New);
    Table anttab (snt);
    std::ifstream antenna_pos (argv[1]);
    for (int i = 0;; ++i)
        {
            std::string station_name;
            double x, y, z;
            antenna_pos >> station_name >> x >> y >> z;
            if (!antenna_pos.good ())
                {
                    std::cerr << i << " antennas written" << std::endl;
                    break;
                }
            std::cerr << station_name << " " << x << " " << y << " " << z << endl;
            anttab.addRow ();
            ScalarColumn<double> (anttab, "DISH_DIAMETER").put (i, 30);
            ScalarColumn<bool> (anttab, "FLAG_ROW").put (i, false);
            Array<double> phref (IPosition (1, 3));
            phref (IPosition (1, 0)) = 0;
            phref (IPosition (1, 1)) = 0;
            phref (IPosition (1, 2)) = 0;
            ArrayColumn<double> (anttab, "PHASE_REFERENCE").put (i, phref);
            ScalarColumn<int> (anttab, "STATION_ID").put (i, i);
            ScalarColumn<String> (anttab, "MOUNT").put (i, "alt-az");
            ScalarColumn<String> (anttab, "NAME").put (i, station_name);
            Array<double> offset (IPosition (1, 3));
            offset (IPosition (1, 0)) = 0;
            offset (IPosition (1, 1)) = 0;
            offset (IPosition (1, 2)) = 0;
            ArrayColumn<double> (anttab, "OFFSET").put (i, offset);
            Array<double> position (IPosition (1, 3));
            position (IPosition (1, 0)) = x;
            position (IPosition (1, 1)) = y;
            position (IPosition (1, 2)) = z;
            ArrayColumn<double> (anttab, "POSITION").put (i, position);
            ScalarColumn<String> (anttab, "STATION").put (i, "21CMA");
            ScalarColumn<String> (anttab, "TYPE").put (i, "GROUND-BASED");
        }
    anttab.flush ();
}
