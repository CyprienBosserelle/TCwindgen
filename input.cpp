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
		// Need to remove leading or trailing spaces

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

std::vector<TCparam> readBSHfile(std::string BSHfilename)
{
	//
	std::vector<TCparam> testTCparamlist;

	TCparam TCparinline;

	TCparinline.TClat = -18.0;//Latitude of TC centre
	TCparinline.TClon = 178.0;//Longitude of TC centre

	TCparinline.cP = 900.0; //central pressure hpa
	TCparinline.eP = 1013.0; //Env pressure hpa
	TCparinline.rMax = 40.0; // Radius of maximum wind (km)
	TCparinline.vFm = 15.0; //Foward speed of the storm(m / s)
	TCparinline.thetaFm = 180.0; //Forward direction of the storm(geographic bearing, positive clockwise);
	TCparinline.beta = 1.30;
	TCparinline.rho = 1.15;

	std::ifstream fs(BSHfilename);

	if (fs.fail()){
		std::cerr << BSHfilename <<" file could not be opened" << std::endl;
		exit(1);
	}

	std::string line;
	while (std::getline(fs, line))
	{
		//std::cout << line << std::endl;

		// skip empty lines
		if (!line.empty())
		{
			//Data should be in teh format :
			//BASIN,CY,YYYYMMDDHH,TECHNUM/MIN,TECH,TAU,LatN/S,LonE/W,VMAX,MSLP,TY,RAD,WINDCODE,RAD1,RAD2,RAD3,RAD4,RADP,RRP,MRD,GUSTS,EYE,SUBREGION,MAXSEAS,INITIALS,DIR,SPEED,STORMNAME,DEPTH,SEAS,SEASCODE,SEAS1,SEAS2,SEAS3,SEAS4,USERDEFINED,userdata

			TCparinline = readBSHline(line);
			testTCparamlist.push_back(TCparinline);
			//std::cout << line << std::endl;
		}

	}
	fs.close();

	return testTCparamlist;


}


TCparam readBSHline(std::string line)
{
	TCparam TCparinline;
	std::size_t found, Numberstart, Numberend;
	//Data should be in teh format :
	//BASIN,CY,YYYYMMDDHH,TECHNUM/MIN,TECH,TAU,LatN/S,LonE/W,VMAX,MSLP,TY,RAD,WINDCODE,RAD1,RAD2,RAD3,RAD4,RADP,RRP,MRD,GUSTS,EYE,SUBREGION,MAXSEAS,INITIALS,DIR,SPEED,STORMNAME,DEPTH,SEAS,SEASCODE,SEAS1,SEAS2,SEAS3,SEAS4,USERDEFINED,userdata
	
	// skip the first 3 items
	std::vector<std::string> x = split(line, ',');

	//std::cout << x[1] << std::endl;

	TCparinline.TClat = -18.0;//Latitude of TC centre
	TCparinline.TClon = 178.0;//Longitude of TC centre

	TCparinline.cP = 900.0; //central pressure hpa
	TCparinline.eP = 1013.0; //Env pressure hpa
	TCparinline.rMax = 40.0; // Radius of maximum wind (km)
	TCparinline.vFm = 15.0; //Foward speed of the storm(m / s)
	TCparinline.thetaFm = 180.0; //Forward direction of the storm(geographic bearing, positive clockwise);
	TCparinline.beta = 1.30;
	TCparinline.rho = 1.15;

	
	//Read Latitude of TC centre
	std::string Latstring=x[6];
	std::string Latstrnum = Latstring.substr(0, Latstring.length()-1);
	found = Latstring.find("S");
	if (found != std::string::npos) // found a line that has Lonmin
	{
		TCparinline.TClat = std::stod(Latstrnum) / 10.0 * -1.0;//
	}
	else
	{
		TCparinline.TClat = std::stod(Latstrnum) / 10.0;
	}
		
	//read Longitude of TC centre
	std::string Lonstring = x[7];
	std::string Lonstrnum = Lonstring.substr(0, Lonstring.length() - 1);
	found = Lonstring.find("W");
	if (found != std::string::npos) // found a line that has Lonmin
	{
		TCparinline.TClon = std::stod(Lonstrnum) / 10.0 * -1.0;//
	}
	else
	{
		TCparinline.TClon = std::stod(Lonstrnum) / 10.0 *-1.0;
	}
	

	TCparinline.cP = std::stod(x[9]) ; //central pressure hpa
	TCparinline.eP = 1013.0; //Env pressure hpa
	TCparinline.rMax = std::stod(x[19])*1.852; // Radius of maximum wind (km)
	TCparinline.vFm = std::stod(x[26])*0.5144444; //Foward speed of the storm(m / s)
	TCparinline.thetaFm = std::stod(x[25])*1.0; //Forward direction of the storm(geographic bearing, positive clockwise);
	//TCparinline.beta = 1.30;
	//TCparinline.rho = 1.15;

	
	return TCparinline;

}

void split(const std::string &s, char delim, std::vector<std::string> &elems) {
	std::stringstream ss;
	ss.str(s);
	std::string item;
	while (std::getline(ss, item, delim)) {
		elems.push_back(item);
	}
}


std::vector<std::string> split(const std::string &s, char delim) {
	std::vector<std::string> elems;
	split(s, delim, elems);
	return elems;
}