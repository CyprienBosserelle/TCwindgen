#include "TCwindgen.h"


std::string findparameter(std::string parameterstr, std::string line)
{
	std::size_t found, Numberstart, Numberend;
	std::string parameternumber;
	found = line.find(parameterstr);
	if (found != std::string::npos) // found a line that has Lonmin
	{
		//std::cout <<"found LonMin at : "<< found << std::endl;
		Numberstart = found + parameterstr.length();
		found = line.find(";");
		if (found != std::string::npos) // found a line that has Lonmin
		{
			Numberend = found;
		}
		else
		{
			Numberend = line.length();
		}
		parameternumber = line.substr(Numberstart, Numberend - Numberstart);
		//std::cout << parameternumber << std::endl;

	}
	return parameternumber;
}

param readparamstr(std::string line, param grid)
{


	std::string parameterstr, parametervalue;

	//
	parameterstr = "LonMin =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		grid.LonMin = std::stod(parametervalue);
	}

	//
	parameterstr = "LonMax =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		grid.LonMax = std::stod(parametervalue);
	}

	//
	parameterstr = "LatMin =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		grid.LatMin = std::stod(parametervalue);
	}

	//
	parameterstr = "LatMax =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		grid.LatMax = std::stod(parametervalue);
	}

	//
	parameterstr = "dlon =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		grid.dlon = std::stod(parametervalue);
	}

	//
	parameterstr = "dlat =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		grid.dlat = std::stod(parametervalue);
	}

	//
	parameterstr = "trackfile =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		grid.Trackfile = parametervalue;
	}

	//
	parameterstr = "ncoutputfile =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		grid.Outputncfile = parametervalue;
		//std::cout << grid.Outputfile << std::endl;
	}
	return grid;
}
