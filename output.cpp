#include "TCwindgen.h"


void write_text_to_log_file(const std::string &text)
{
	std::ofstream log_file(
		"TCwindgen_log.txt", std::ios_base::out | std::ios_base::app);
	log_file << text << std::endl;
	log_file.close(); //destructoir implicitly does it
}

extern "C" void creatncfile(std::string ncfileout, int nx, int ny, float totaltime, float * xval, float * yval, float *R, float *V, float *Z)
{
	int status;
	int ncid, xx_dim, yy_dim, time_dim, p_dim;
	size_t nxx, nyy, nnpart;
	int  var_dimids[3], var_dimzb[2];

	int R_id, time_id, xx_id, yy_id, V_id, Z_id;

	nxx = nx;
	nyy = ny;
	//nnpart=npart;

	const char *outfile = ncfileout.c_str();

	static size_t start[] = { 0, 0, 0 }; // start at first value 
	static size_t count[] = { 1, ny, nx };
	static size_t zbstart[] = { 0, 0 }; // start at first value 
	static size_t zbcount[] = { ny, nx };
	//static size_t pstart[] = {0, 0}; // start at first value 
	// 	static size_t pcount[] = {1, npart};
	static size_t tst[] = { 0 };
	static size_t xstart[] = { 0 }; // start at first value 
	static size_t xcount[] = { nx };

	static size_t ystart[] = { 0 }; // start at first value 
	static size_t ycount[] = { ny };



	//create the netcdf dataset
	status = nc_create(outfile, NC_NOCLOBBER | NC_NETCDF4, &ncid);

	//Define dimensions: Name and length

	status = nc_def_dim(ncid, "x", nxx, &xx_dim);
	status = nc_def_dim(ncid, "y", nyy, &yy_dim);
	//status = nc_def_dim(ncid, "npart",nnpart,&p_dim);
	status = nc_def_dim(ncid, "time", NC_UNLIMITED, &time_dim);
	

	int tdim[] = { time_dim };
	int xdim[] = { xx_dim };
	int ydim[] = { yy_dim };
	//int pdim[2];
	//pdim[0]=time_dim;
	//pdim[1]=p_dim;
	//define variables: Name, Type,...
	var_dimids[0] = time_dim;
	var_dimids[1] = yy_dim;
	var_dimids[2] = xx_dim;
	var_dimzb[0] = yy_dim;
	var_dimzb[1] = xx_dim;

	status = nc_def_var(ncid, "time", NC_FLOAT, 1, tdim, &time_id);
	status = nc_def_var(ncid, "x", NC_FLOAT, 1, xdim, &xx_id);
	status = nc_def_var(ncid, "y", NC_FLOAT, 1, ydim, &yy_id);

	static char txtname[] = "time";
	status = nc_put_att_text(ncid, time_id, "standard_name", strlen(txtname), txtname);
	std::string timestr = "seconds";
	const char* timeunit = timestr.c_str();

	status = nc_put_att_text(ncid, time_id, "units", strlen(timeunit), timeunit);

	std::string xaxis = "X";
	std::string yaxis = "Y";
	status = nc_put_att_text(ncid, xx_id, "axis", xaxis.size(), xaxis.c_str());
	status = nc_put_att_text(ncid, yy_id, "axis", yaxis.size(), yaxis.c_str());

	int shuffle = 1;
	int deflate = 1;        // This switches compression on (1) or off (0).
	int deflate_level = 9;  // This is the compression level in range 1 (less) - 9 (more).

	status = nc_def_var(ncid, "P", NC_FLOAT, 3, var_dimids, &R_id);
	std::string Pstdname = "pressure";
	std::string Plongname = "Atmospheric pressure";
	std::string Punit = "pa";

	status = nc_put_att_text(ncid, R_id, "standard_name", Pstdname.size(), Pstdname.c_str());
	status = nc_put_att_text(ncid, R_id, "long_name", Plongname.size(), Plongname.c_str());
	status = nc_put_att_text(ncid, R_id, "units", Punit.size(), Punit.c_str());
	nc_def_var_deflate(ncid, R_id, shuffle, deflate, deflate_level);


	status = nc_def_var(ncid, "U", NC_FLOAT, 3, var_dimids, &V_id);
	std::string Ustdname = "uwind";
	std::string Ulongname = "U Wind";
	std::string Uunit = "m/s";

	status = nc_put_att_text(ncid, V_id, "standard_name", Ustdname.size(), Ustdname.c_str());
	status = nc_put_att_text(ncid, V_id, "long_name", Ulongname.size(), Ulongname.c_str());
	status = nc_put_att_text(ncid, V_id, "units", Uunit.size(), Uunit.c_str());
	nc_def_var_deflate(ncid, V_id, shuffle, deflate, deflate_level);

	status = nc_def_var(ncid, "V", NC_FLOAT, 3, var_dimids, &Z_id);
	std::string Vstdname = "vwind";
	std::string Vlongname = "V Wind";
	std::string Vunit = "m/s";

	status = nc_put_att_text(ncid, Z_id, "standard_name", Vstdname.size(), Vstdname.c_str());
	status = nc_put_att_text(ncid, Z_id, "long_name", Vlongname.size(), Vlongname.c_str());
	status = nc_put_att_text(ncid, Z_id, "units", Vunit.size(), Vunit.c_str());
	nc_def_var_deflate(ncid, Z_id, shuffle, deflate, deflate_level);
	//put attriute: assign attibute values
	//nc_put_att

	//End definitions: leave define mode
	status = nc_enddef(ncid);

	//Provide values for variables
	status = nc_put_var1_float(ncid, time_id, tst, &totaltime);
	status = nc_put_vara_float(ncid, xx_id, xstart, xcount, xval);
	status = nc_put_vara_float(ncid, yy_id, ystart, ycount, yval);


	status = nc_put_vara_float(ncid, R_id, start, count, R);
	status = nc_put_vara_float(ncid, V_id, start, count, V);
	status = nc_put_vara_float(ncid, Z_id, start, count, Z);
	// U, V, P are also needed

	//close and save new file
	status = nc_close(ncid);
}


extern "C" void writestep2nc(std::string ncfileout, int nx, int ny, float totaltime, float *R, float *V, float *Z)
{
	int status;
	int ncid, time_dim, recid;
	size_t nxx, nyy;
	int time_id, R_id, V_id, Z_id;
	static size_t start[] = { 0, 0, 0 }; // start at first value 
	static size_t count[] = { 1, ny, nx };
	//static size_t pstart[] = {0, 0}; // start at first value 
	//	static size_t pcount[] = {1, npart};
	static size_t tst[] = { 0 };

	nxx = nx;
	nyy = ny;

	const char *outfile = ncfileout.c_str();

	static size_t nrec;
	status = nc_open(outfile, NC_WRITE, &ncid);

	//read id from time dimension
	status = nc_inq_unlimdim(ncid, &recid);
	status = nc_inq_dimlen(ncid, recid, &nrec);
	//printf("nrec=%d\n", nrec);

	//read file for variable ids
	status = nc_inq_varid(ncid, "time", &time_id);

	status = nc_inq_varid(ncid, "P", &R_id);
	status = nc_inq_varid(ncid, "U", &V_id);
	status = nc_inq_varid(ncid, "V", &Z_id);

	//status = nc_inq_varid(ncid, "xxp", &xxp_id);
	//status = nc_inq_varid(ncid, "yyp", &yyp_id);


	start[0] = nrec;
	//pstart[0] = nrec;    
	tst[0] = nrec;

	//Provide values for variables
	status = nc_put_var1_float(ncid, time_id, tst, &totaltime);

	status = nc_put_vara_float(ncid, R_id, start, count, R);
	status = nc_put_vara_float(ncid, V_id, start, count, V);
	status = nc_put_vara_float(ncid, Z_id, start, count, Z);




	//close and save
	status = nc_close(ncid);


}

void createSWANwindfileFPF(std::string SWANfileOut, int nx, int ny, float *U, float *V)
{
	FILE * ofile;
	const char * outfile = SWANfileOut.c_str();

	ofile = fopen(outfile, "w");

	for (int j = 0; j < ny; j++)
	{
		for (int i = 0; i < nx; i++)
		{
			//fs.precision(2);
			fprintf(ofile, "%.2f\t", U[i + j*nx]);
		}
		fprintf(ofile, "\n");
	}
	for (int j = 0; j < ny; j++)
	{
		for (int i = 0; i < nx; i++)
		{
			//fs.precision(2);
			fprintf(ofile, "%.2f\t", V[i + j*nx]);;
		}
		fprintf(ofile, "\n");
	}
}
void createSWANwindfile(std::string SWANfileOut, int nx, int ny,float *U, float *V)
{
	//
	std::ofstream fs(SWANfileOut);

	if (fs.fail()){
		std::cerr << SWANfileOut << "file could not be opened" << std::endl;
		exit(1);
	}
	
	
	for (int j = 0; j < ny; j++)
	{
		for (int i = 0; i < nx; i++)
		{
			//fs.precision(2);
			fs << std::fixed << std::setprecision(2) << U[i + j*nx] << "\t";
		}
		fs << "\n";
	}
	for (int j = 0; j < ny; j++)
	{
		for (int i = 0; i < nx; i++)
		{
			//fs.precision(2);
			fs << std::fixed << std::setprecision(2) << V[i + j*nx] << "\t";
		}
		fs << "\n";
	}
	fs.close();
}

void writeSWANWindstepFPF(std::string SWANfileOut, int nx, int ny, float *U, float *V)
{
	//
	FILE * ofile;
	const char * outfile = SWANfileOut.c_str();

	ofile = fopen(outfile, "a");

	for (int j = 0; j < ny; j++)
	{
		for (int i = 0; i < nx; i++)
		{
			//fs.precision(2);
			fprintf(ofile, "%.2f\t", U[i + j*nx]);
		}
		fprintf(ofile, "\n");
	}
	for (int j = 0; j < ny; j++)
	{
		for (int i = 0; i < nx; i++)
		{
			//fs.precision(2);
			fprintf(ofile, "%.2f\t", V[i + j*nx]);;
		}
		fprintf(ofile, "\n");
	}
}
void writeSWANWindstep(std::string SWANfileOut, int nx, int ny, float *U, float *V)
{
	//
	std::ofstream fs;
	fs.open(SWANfileOut, std::ios::out | std::ios::app);

	if (fs.fail()){
		std::cerr << SWANfileOut << "file could not be opened" << std::endl;
		exit(1);
	}
	
	for (int j = 0; j < ny; j++)
	{
		for (int i = 0; i < nx; i++)
		{
			
			fs << std::fixed << std::setprecision(2) << U[i + j*nx] << "\t";
		}
		fs << "\n";
	}
	for (int j = 0; j < ny; j++)
	{
		for (int i = 0; i < nx; i++)
		{
			//fs.precision(2);
			fs << std::fixed << std::setprecision(2) << V[i + j*nx] << "\t";
		}
		fs << "\n";
	}
	fs.close();
}
void createD3DAtmfileFPF(std::string D3dfileOut, int nx, int ny, double lonmin, double latmin, double dlon, double dlat, tm datestart, float *P, float *U, float *V)
{
	FILE * fsUU, *fsVV, *fsPP;
	std::string Ufile, Vfile, Pfile;

	Ufile = D3dfileOut + ".amu";
	Vfile = D3dfileOut + ".amv";
	Pfile = D3dfileOut + ".amp";

	//const char * outfile = U.c_str();

	fsUU = fopen(Ufile.c_str(), "w");
	fsVV = fopen(Vfile.c_str(), "w");
	fsPP = fopen(Pfile.c_str(), "w");

	time_t rawtime, dstart;
	struct tm * timeinfo;
	char buffer[80];

	time(&rawtime);
	timeinfo = localtime(&rawtime);

	strftime(buffer, 80, "%d-%m-%Y ", timeinfo);
	std::string strtimenow(buffer);

	strftime(buffer, 80, "%Y-%m-%d %H:%M:%S", &datestart);
	std::string strdatestart(buffer);

	//std::cout << str;

	// write header for U file
	fprintf(fsUU ,"### START OF HEADER\n");
	fprintf(fsUU, "### Created by TCWindgen on %s\n" , strtimenow.c_str() );
	fprintf(fsUU, "### Put other header info Here\n");
	fprintf(fsUU, "Fileversion\t=\t1.0\n");
	fprintf(fsUU, "filetype\t=\tmeteo_on_equidistant_grid\n");
	fprintf(fsUU, "NODATA_value\t=\t-9999.0\n");
	fprintf(fsUU , "n_cols\t=\t%i\n" , nx );
	fprintf(fsUU , "n_rows\t=\t%i\n", ny);
	fprintf(fsUU , "grid_unit\t=\tdegree\n");
	fprintf(fsUU , "x_llcenter\t=\t%f\n" ,lonmin );
	fprintf(fsUU , "y_llcenter\t=\t%f\n" , latmin );
	fprintf(fsUU , "dx\t=\t%f\n" , dlon );
	fprintf(fsUU , "dy\t=\t%f\n" , dlat );
	fprintf(fsUU , "n_quantity\t=\t1\n");
	fprintf(fsUU , "quantity1\t=\tx_wind\n");
	fprintf(fsUU , "unit1\t=\tm s-1\n");
	fprintf(fsUU , "### END OF HEADER\n");
	fprintf(fsUU , "TIME = 0 seconds since %s+00:00\n" , strdatestart.c_str() );

	for (int j = 0; j < ny; j++)
	{
		for (int i = 0; i < nx; i++)
		{
			//fs.precision(2);
			fprintf(fsUU, "%.2f\t" , U[i + j*nx] );
		}
		fprintf(fsUU, "\n");
	}
	// write header for V file
	fprintf(fsVV, "### START OF HEADER\n");
	fprintf(fsVV, "### Created by TCWindgen on %s\n", strtimenow.c_str());
	fprintf(fsVV, "### Put other header info Here\n");
	fprintf(fsVV, "Fileversion\t=\t1.0\n");
	fprintf(fsVV, "filetype\t=\tmeteo_on_equidistant_grid\n");
	fprintf(fsVV, "NODATA_value\t=\t-9999.0\n");
	fprintf(fsVV, "n_cols\t=\t%i\n", nx);
	fprintf(fsVV, "n_rows\t=\t%i\n", ny);
	fprintf(fsVV, "grid_unit\t=\tdegree\n");
	fprintf(fsVV, "x_llcenter\t=\t%f\n", lonmin);
	fprintf(fsVV, "y_llcenter\t=\t%f\n", latmin);
	fprintf(fsVV, "dx\t=\t%f\n", dlon);
	fprintf(fsVV, "dy\t=\t%f\n", dlat);
	fprintf(fsVV, "n_quantity\t=\t1\n");
	fprintf(fsVV, "quantity1\t=\ty_wind\n");
	fprintf(fsVV, "unit1\t=\tm s-1\n");
	fprintf(fsVV, "### END OF HEADER\n");
	fprintf(fsVV, "TIME = 0 seconds since %s+00:00\n", strdatestart.c_str());

	for (int j = 0; j < ny; j++)
	{
		for (int i = 0; i < nx; i++)
		{
			//fs.precision(2);
			fprintf(fsVV, "%.2f\t", V[i + j*nx]);
		}
		fprintf(fsVV, "\n");
	}
	// write header for P file
	fprintf(fsPP, "### START OF HEADER\n");
	fprintf(fsPP, "### Created by TCWindgen on %s\n", strtimenow.c_str());
	fprintf(fsPP, "### Put other header info Here\n");
	fprintf(fsPP, "Fileversion\t=\t1.0\n");
	fprintf(fsPP, "filetype\t=\tmeteo_on_equidistant_grid\n");
	fprintf(fsPP, "NODATA_value\t=\t-9999.0\n");
	fprintf(fsPP, "n_cols\t=\t%i\n", nx);
	fprintf(fsPP, "n_rows\t=\t%i\n", ny);
	fprintf(fsPP, "grid_unit\t=\tdegree\n");
	fprintf(fsPP, "x_llcenter\t=\t%f\n", lonmin);
	fprintf(fsPP, "y_llcenter\t=\t%f\n", latmin);
	fprintf(fsPP, "dx\t=\t%f\n", dlon);
	fprintf(fsPP, "dy\t=\t%f\n", dlat);
	fprintf(fsPP, "n_quantity\t=\t1\n");
	fprintf(fsPP, "quantity1\t=\tair_pressure\n");
	fprintf(fsPP, "unit1\t=\tPa\n");
	fprintf(fsPP, "### END OF HEADER\n");
	fprintf(fsPP, "TIME = 0 seconds since %s+00:00\n", strdatestart.c_str());

	for (int j = 0; j < ny; j++)
	{
		for (int i = 0; i < nx; i++)
		{
			//fs.precision(2);
			fprintf(fsPP, "%.2f\t", P[i + j*nx]);
		}
		fprintf(fsPP, "\n");
	}
	fclose(fsUU);
	fclose(fsVV);
	fclose(fsPP);
}
void createD3DAtmfile(std::string D3dfileOut, int nx, int ny, double lonmin, double latmin, double dlon, double dlat, tm datestart, float *P, float *U, float *V)
{

	std::ofstream fsUU, fsVV, fsPP;

	time_t rawtime,dstart;
	struct tm * timeinfo;
	char buffer[80];

	time(&rawtime);
	timeinfo = localtime(&rawtime);

	strftime(buffer, 80, "%d-%m-%Y ", timeinfo);
	std::string strtimenow(buffer);

	strftime(buffer, 80, "%Y-%m-%d %H:%M:%S", &datestart); 
	std::string strdatestart(buffer);

	//std::cout << str;


	fsUU.open(D3dfileOut + "_U.amu", std::ios::out | std::ios::trunc);
	fsVV.open(D3dfileOut + "_V.amv", std::ios::out | std::ios::trunc);
	fsPP.open(D3dfileOut + "_P.amp", std::ios::out | std::ios::trunc);
	//
	if (fsUU.fail() || fsVV.fail() || fsPP.fail()){
		std::cerr << D3dfileOut << "files could not be opened" << std::endl;
		exit(1);
	}

	// write header for U file
	fsUU << "### START OF HEADER\n";
	fsUU << "### Created by TCWindgen on " << strtimenow << "\n";
	fsUU << "### Put other header info Here\n";
	fsUU << "Fileversion\t=\t1.0\n";
	fsUU << "filetype\t=\tmeteo_on_equidistant_grid\n";
	fsUU << "NODATA_value\t=\t-9999.0\n";
	fsUU << "n_cols\t=\t"<< nx <<"\n";
	fsUU << "n_rows\t=\t"<< ny <<"\n";
	fsUU << "grid_unit\t=\tdegree\n";
	fsUU << "x_llcenter\t=\t" << std::fixed << std::setprecision(7) << lonmin << "\n";
	fsUU << "y_llcenter\t=\t" << std::fixed << std::setprecision(7) << latmin << "\n";
	fsUU << "dx\t=\t" << std::fixed << std::setprecision(7) << dlon << "\n";
	fsUU << "dy\t=\t" << std::fixed << std::setprecision(7) << dlat << "\n";
	fsUU << "n_quantity\t=\t1\n";
	fsUU << "quantity1\t=\tx_wind\n";
	fsUU << "unit1\t=\tm s-1\n";
	fsUU << "### END OF HEADER\n";
	fsUU << "TIME = 0 seconds since " << strdatestart << "+00:00\n";
	for (int j = 0; j < ny; j++)
	{
		for (int i = 0; i < nx; i++)
		{
			//fs.precision(2);
			fsUU << std::fixed << std::setprecision(2) << U[i + j*nx] << "\t";
		}
		fsUU << "\n";
	}

	// write header for V file
	fsVV << "### START OF HEADER\n";
	fsVV << "### Created by TCWindgen on " << strtimenow << "\n";
	fsVV << "### Put other header info Here\n";
	fsVV << "Fileversion\t=\t1.0\n";
	fsVV << "filetype\t=\tmeteo_on_equidistant_grid\n";
	fsVV << "NODATA_value\t=\t-9999.0\n";
	fsVV << "n_cols\t=\t" << nx << "\n";
	fsVV << "n_rows\t=\t" << ny << "\n";
	fsVV << "grid_unit\t=\tdegree\n";
	fsVV << "x_llcenter\t=\t" << std::fixed << std::setprecision(7) << lonmin << "\n";
	fsVV << "y_llcenter\t=\t" << std::fixed << std::setprecision(7) << latmin << "\n";
	fsVV << "dx\t=\t" << std::fixed << std::setprecision(7) << dlon << "\n";
	fsVV << "dy\t=\t" << std::fixed << std::setprecision(7) << dlat << "\n";
	fsVV << "n_quantity\t=\t1\n";
	fsVV << "quantity1\t=\ty_wind\n";
	fsVV << "unit1\t=\tm s-1\n";
	fsVV << "### END OF HEADER\n";
	fsVV << "TIME = 0 seconds since " << strdatestart << "+00:00\n";
	for (int j = 0; j < ny; j++)
	{
		for (int i = 0; i < nx; i++)
		{
			//fs.precision(2);
			fsVV << std::fixed << std::setprecision(2) << V[i + j*nx] << "\t";
		}
		fsVV << "\n";
	}

	// write header for P file
	fsPP << "### START OF HEADER\n";
	fsPP << "### Created by TCWindgen on " << strtimenow << "\n";
	fsPP << "### Put other header info Here\n";
	fsPP << "Fileversion\t=\t1.0\n";
	fsPP << "filetype\t=\tmeteo_on_equidistant_grid\n";
	fsPP << "NODATA_value\t=\t-9999.0\n";
	fsPP << "n_cols\t=\t" << nx << "\n";
	fsPP << "n_rows\t=\t" << ny << "\n";
	fsPP << "grid_unit\t=\tdegree\n";
	fsPP << "x_llcenter\t=\t" << std::fixed << std::setprecision(7) << lonmin << "\n";
	fsPP << "y_llcenter\t=\t" << std::fixed << std::setprecision(7) << latmin << "\n";
	fsPP << "dx\t=\t" << std::fixed << std::setprecision(7) << dlon << "\n";
	fsPP << "dy\t=\t" << std::fixed << std::setprecision(7) << dlat << "\n";
	fsPP << "n_quantity\t=\t1\n";
	fsPP << "quantity1\t=\tair_pressure\n";
	fsPP << "unit1\t=\tPa\n";
	fsPP << "### END OF HEADER\n";
	fsPP << "TIME = 0 seconds since " << strdatestart << "+00:00\n";
	for (int j = 0; j < ny; j++)
	{
		for (int i = 0; i < nx; i++)
		{
			//fs.precision(2);
			fsPP << std::fixed << std::setprecision(2) << P[i + j*nx] << "\t";
		}
		fsPP << "\n";
	}

	fsUU.close();
	fsVV.close();
	fsPP.close();



}
void writeD3DAtmstepFPF(std::string D3dfileOut, int nx, int ny, tm datestart, double totaltime, float *P, float *U, float *V)
{
	FILE * fsUU, *fsVV, *fsPP;
	char buffer[80];

	strftime(buffer, 80, "%Y-%m-%d %H:%M:%S", &datestart);
	std::string strdatestart(buffer);

	std::string Ufile, Vfile, Pfile;

	Ufile = D3dfileOut + ".amu";
	Vfile = D3dfileOut + ".amv";
	Pfile = D3dfileOut + ".amp";

	//const char * outfile = U.c_str();

	fsUU = fopen(Ufile.c_str(), "a");
	fsVV = fopen(Vfile.c_str(), "a");
	fsPP = fopen(Pfile.c_str(), "a");

	fprintf(fsUU ,"TIME = %i seconds since %s+00:00\n" ,(int) round(totaltime) , strdatestart.c_str() );
	fprintf(fsVV, "TIME = %i seconds since %s+00:00\n", (int)round(totaltime), strdatestart.c_str());
	fprintf(fsPP, "TIME = %i seconds since %s+00:00\n", (int)round(totaltime), strdatestart.c_str());

	for (int j = 0; j < ny; j++)
	{
		for (int i = 0; i < nx; i++)
		{
			//fs.precision(2);
			fprintf(fsUU, "%.2f\t", U[i + j*nx]);
		}
		fprintf(fsUU, "\n");
	}

	for (int j = 0; j < ny; j++)
	{
		for (int i = 0; i < nx; i++)
		{
			//fs.precision(2);
			fprintf(fsPP, "%.2f\t", P[i + j*nx]);
		}
		fprintf(fsPP, "\n");
	}

	for (int j = 0; j < ny; j++)
	{
		for (int i = 0; i < nx; i++)
		{
			//fs.precision(2);
			fprintf(fsVV, "%.2f\t", V[i + j*nx]);
		}
		fprintf(fsVV, "\n");
	}
	fclose(fsUU);
	fclose(fsVV);
	fclose(fsPP);

}

void writeD3DAtmstep(std::string D3dfileOut, int nx, int ny, tm datestart,double totaltime, float *P, float *U, float *V)
{
	std::ofstream fsUU, fsVV, fsPP;
	char buffer[80];

	strftime(buffer, 80, "%Y-%m-%d %H:%M:%S", &datestart);
	std::string strdatestart(buffer);

	fsUU.open(D3dfileOut + "_U.amu", std::ios::out | std::ios::app);
	fsVV.open(D3dfileOut + "_V.amv", std::ios::out | std::ios::app);
	fsPP.open(D3dfileOut + "_P.amp", std::ios::out | std::ios::app);
	//
	fsUU << "TIME = "<< (int) round(totaltime) <<" seconds since " << strdatestart << "+00:00\n";
	for (int j = 0; j < ny; j++)
	{
		for (int i = 0; i < nx; i++)
		{
			//fs.precision(2);
			fsUU << std::fixed << std::setprecision(2) << U[i + j*nx] << "\t";
		}
		fsUU << "\n";
	}

	fsVV << "TIME = " << (int)round(totaltime) << " seconds since " << strdatestart << "+00:00\n";
	for (int j = 0; j < ny; j++)
	{
		for (int i = 0; i < nx; i++)
		{
			//fs.precision(2);
			fsVV << std::fixed << std::setprecision(2) << V[i + j*nx] << "\t";
		}
		fsVV << "\n";
	}

	fsPP << "TIME = " << (int)round(totaltime) << " seconds since " << strdatestart << "+00:00\n";
	for (int j = 0; j < ny; j++)
	{
		for (int i = 0; i < nx; i++)
		{
			//fs.precision(2);
			fsPP << std::fixed << std::setprecision(2) << P[i + j*nx] << "\t";
		}
		fsPP << "\n";
	}

}