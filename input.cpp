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
	return trim(parameternumber," ");
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

	//grid.datestart
	parameterstr = "datestart =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		//need check that string is 15 character long 
		grid.datestart.tm_year=std::stoi(parametervalue.substr(0,4))-1900; // starting from 1900
		grid.datestart.tm_mon = std::stoi(parametervalue.substr(4, 2)) - 1; //0 to 11
		grid.datestart.tm_mday = std::stoi(parametervalue.substr(6, 2));
		grid.datestart.tm_hour = std::stoi(parametervalue.substr(9, 2));
		grid.datestart.tm_min = std::stoi(parametervalue.substr(11, 2));

		grid.datestart.tm_sec = std::stoi(parametervalue.substr(13, 2));




		//std::cout << grid.Outputfile << std::endl;
	}
	parameterstr = "dateend =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		//need to check that string is 15 character long 
		grid.dateend.tm_year = std::stoi(parametervalue.substr(0, 4)) - 1900; // starting from 1900
		grid.dateend.tm_mon = std::stoi(parametervalue.substr(4, 2)) - 1; //0 to 11
		grid.dateend.tm_mday = std::stoi(parametervalue.substr(6, 2));
		grid.dateend.tm_hour = std::stoi(parametervalue.substr(9, 2));
		grid.dateend.tm_min = std::stoi(parametervalue.substr(11, 2));

		grid.dateend.tm_sec = std::stoi(parametervalue.substr(13, 2));




		//std::cout << grid.Outputfile << std::endl;
	}

	parameterstr = "dt =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		//need to check that string is 15 character long 
		grid.dt = std::stod(parametervalue);
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
	TCparinline.rMax = 47.0; // Radius of maximum wind (km)
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
	//TCparinline.datetime.sec
	std::string parametervalue = trim(x[2]," ");

	//std::cout << x[2] << std::endl;
	TCparinline.datetime.tm_year = std::stoi(parametervalue.substr(0, 4)) - 1900; // starting from 1900
	TCparinline.datetime.tm_mon = std::stoi(parametervalue.substr(4, 2)) - 1; //0 to 11
	TCparinline.datetime.tm_mday = std::stoi(parametervalue.substr(6, 2));
	TCparinline.datetime.tm_hour = std::stoi(parametervalue.substr(8, 2));
	TCparinline.datetime.tm_min = 0;

	TCparinline.datetime.tm_sec = 0;
	




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
		TCparinline.TClon = 360.0 - std::stod(Lonstrnum) / 10.0;//
	}
	else
	{
		TCparinline.TClon = std::stod(Lonstrnum) / 10.0;
	}
	
	



	TCparinline.cP = std::stod(x[9]) ; //central pressure hpa
	TCparinline.eP = 1013.0; //Env pressure hpa
	TCparinline.rMax = std::stod(x[19])*1.852; // Radius of maximum wind (km)
	if (TCparinline.rMax < 0.1)
	{
		TCparinline.rMax = 47.0; //default is set to the mean value estimated by :
		//S. A. Hsu and Zhongde Yana (Spring 1998). "A Note on the Radius of Maximum Winds for Hurricanes". Journal of Coastal Research. Coastal Education & Research Foundation, Inc. 12 (2): 667–668. JSTOR 4298820.
	}
	TCparinline.vFm = std::stod(x[26])*0.5144444; //Foward speed of the storm(m / s)
	// if the forward speed is 0.0 it should be calculated again
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

std::string trim(const std::string& str, const std::string& whitespace)
{
	const auto strBegin = str.find_first_not_of(whitespace);
	if (strBegin == std::string::npos)
		return ""; // no content

	const auto strEnd = str.find_last_not_of(whitespace);
	const auto strRange = strEnd - strBegin + 1;

	return str.substr(strBegin, strRange);
}

template <class T> const T& max(const T& a, const T& b) {
	return (a<b) ? b : a;     // or: return comp(a,b)?b:a; for version (2)
}

template <class T> const T& min(const T& a, const T& b) {
	return !(b<a) ? a : b;     // or: return !comp(b,a)?a:b; for version (2)
}

double interptime(double next, double prev, double timenext, double time)
{
	return prev + (time)/(timenext)*(next - prev) ;
}
//windv = windvold + (totaltime - rtwind)*(windvnew - windvold) / (windtime - rtwind);