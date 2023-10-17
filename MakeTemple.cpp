//It needs to make an temple of config for calculatings

#include <fstream>


int main()
{
	std::ofstream out;
    out.open("TEMPLE_CONFIG.txt");
    char *text = 
                "POWER\t\tVOLUME\n"
            "150\t\t\t1350295.792\n\n"

            "HEIGHT\t\tRADIUS\n"
            "130.0\t\t57.5\n\n"

            "TEMPERATURE_OF_MODERATOR\n"
            "571\n\n"

            "NUMBER_OF_ITERATIONS\n"
            "2\n\n"

            "EFFECTIVE_DAYS_STEP\n"
            "200\n\n"

            "RODS_ON_OR_OF\n"
            "1\n\n"

            "B-10\t\tB-11\n"
            "20\t\t\t80\n\n"

            "ELEMENTS\t\t\tCONCENTRATION\t\t\tTEMPERATURE\n"
            "B10\t\t\t\t\t0\t\t\t\t\t\t0\n"
            "B11\t\t\t\t\t0\t\t\t\t\t\t0\n"
            "C\t\t\t\t\t0\t\t\t\t\t\t0\n"
            "Th232\t\t\t\t0\t\t\t\t\t\t0\n"
            "U233\t\t\t\t0\t\t\t\t\t\t0\n"
            "U234\t\t\t\t0\t\t\t\t\t\t0\n"
            "U235\t\t\t\t0.00049732042433\t\t650\n"
            "U236\t\t\t\t0\t\t\t\t\t\t0\n"
            "U238\t\t\t\t0.00217644529786\t\t650\n"
            "Pu239\t\t\t\t0\t\t\t\t\t\t0\n"
            "Pu240\t\t\t\t0\t\t\t\t\t\t0\n"
            "Pu241\t\t\t\t0\t\t\t\t\t\t0\n"
            "Pu242\t\t\t\t0\t\t\t\t\t\t0\n"
            "FF3\t\t\t\t\t0\t\t\t\t\t\t0\n"
            "FF5\t\t\t\t\t0\t\t\t\t\t\t0\n"
            "FF9\t\t\t\t\t0\t\t\t\t\t\t0\n"
            "H\t\t\t\t\t0.02341756645166\t\t571\n"
            "O\t\t\t\t\t0.01705631467019\t\t571\n"
            "Zr\t\t\t\t\t0.00826734217661\t\t571\n"
            "Si\t\t\t\t\t0.00079634440313\t\t571\n"
            "Nb\t\t\t\t\t0.00008199630668\t\t571\n"
            "Al \t\t\t\t\t0.00746048811801\t\t571";
            
            out << text;
}
