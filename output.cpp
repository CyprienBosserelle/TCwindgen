#include "TCwindgen.h"

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
	status = nc_create(outfile, NC_NOCLOBBER, &ncid);

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



	status = nc_def_var(ncid, "P", NC_FLOAT, 3, var_dimids, &R_id);
	status = nc_def_var(ncid, "U", NC_FLOAT, 3, var_dimids, &V_id);
	status = nc_def_var(ncid, "V", NC_FLOAT, 3, var_dimids, &Z_id);

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
	printf("nrec=%d\n", nrec);

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