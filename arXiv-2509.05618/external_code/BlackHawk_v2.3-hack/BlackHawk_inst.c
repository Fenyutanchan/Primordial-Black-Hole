// Copyright (c) 2025 Quan-feng WU <wuquanfeng@ihep.ac.cn>
// 
// This software is released under the MIT License.
// https://opensource.org/licenses/MIT
// 
// Modified from `BlackHawk_inst.c` in BlackHawk v2.3
// Date: 2025-05-20
// ============================================================================
// 
// BlackHawk_inst is a program aimed at computing the Hawking
// instantaneous emission of a distribution of black holes.
// Last modification: 07 April 2024
// Authors: Alexandre Arbey alexandre.arbey@ens-lyon.fr & Jérémy Auffinger jeremy.auffinger.blackhawk@protonmail.com

#include "src/include.h"

// &&&&&&&&&&&&&&&&&&&&& MAIN &&&&&&&&&&&&&&&&&&&&&&&&
// This is the BlackHawk_inst main routine. It executes all intermediate steps
// by calling routines defined in './src'. It returns 0 if the computation was
// completed, 1 otherwise.
int main(int argc,char **argv)
{
	/*
	printf("\n                                                       `-    o/                                     \n");
	printf("                                               --     sy   -N+                                      \n");
	printf("                                               /m`   :My  :Nm` -h/                                  \n");
	printf("                                               +M/  -NM:.hMd.-yNs                                   \n");
	printf("                                               hM+ /NMdoNMh-sNN/  `:`                               \n");
	printf("                                              :MM/sMMMMMMhsNMh. -sd/                                \n");
	printf("                                             `mMMmMMMMMMNNMNo-+dNh-                                 \n");
	printf("          `-   .ho   ``                      yMMMMMMMMMMMMmshNMd/``-+.                              \n");
	printf("         `ds  .mM/  /m:                     /MMMMMMMMMMMMMNMMd+:+hmd+`                              \n");
	printf("    `.   yM/ `mMm `yMh  -o`                `mMMMMMMMMMMMMMMMddmNNy/.`                               \n");
	printf("    yh  /MM` yMM/`yMM-`sNh  :.             oMMMMMMMMMMMMMMMMMMMNddmmd/                              \n");
	printf("   `MM``mMd :MMm`yMMy-dMm`-hm``/`          mMMMMMMMMMMMMMMMMMMMMMMMMs:                              \n");
	printf("   :MM/+MMh`mMMyhMMNoNMN-oNN-+mN`.`      s/MMMMMMMMMMMMMMMMMMMMMMMMMd/                              \n");
	printf("   :MMhmMMNoMMMNMMMNMMMohMN+hMMhsm/`     MNMMMMMMMMMMMMMMMMMMMMMMMMMd/                              \n");
	printf("   .MMMMMMMMMMMMMMMMMMMNMMmNMMMMMMhd-    dMMMMMMMMMMMMMMMMMMMMMMMMMMd.                              \n");
	printf("    mMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM:``  :MMMMMMMMMMMMMMMMMMMMMMMMMMd`                              \n");
	printf("    oMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMdms   hMMMMMMMMMMMMMMMMMMMMMMMMMN/                              \n");
	printf("    .NMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMy.. .NMMMMMMMMMMMMMMMMMMMMMMMMMh                              \n");
	printf("     oMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMNNNy/dMMMMMMMMMMMMMMMMMMMMMMMMMN`                             \n");
	printf("     `dMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM-                             \n");
	printf("      -NMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMy                             \n");
	printf("       /MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM+                             \n");
	printf("        +MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMh                              \n");	
	printf("      +dohMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMN/                              \n");
	printf("       :mMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMo                      -//-     \n");
	printf("         :yNMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMm:                  `/sdNMMMMN.   \n");
	printf("           `/yNMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMNo               `-+ymMMMMMMMMMN    \n");
	printf("              `:sdNMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMNho+:`            ./sdNMMMMMMMMMMMMy.    \n");
	printf("                  ``-/+shmMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMms/-:---........-/ohNMMMMMMMMMMMMMMMN+      \n");
	printf("                         `.+dMMMMMMMMMMMMMMMMMMMMMMMMMMMMNNNMMMNNNNNNNNNMMMMMMMMMMMMMMMMMMm+`       \n");
	printf("                            `/dMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMh/`         \n");
	printf("                              `/NMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMNs-            \n");
	printf("                                oMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMd+.              \n");
	printf("                               /NMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMNy:`                \n");
	printf("                             -hNMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMmdhdNMMMMMMMMNh+.                   \n");
	printf("                           :hNMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMdy+-.` `.-/oyhhy+.`                     \n");
	printf("                           NMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMNs`                                       \n");
	printf("                          -MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMN.                                         \n");
	printf("                          +MMMMMMMmy+//+NMMMMMMNo+NMMMMMMd                                          \n");
	printf("                          -MMMMMm+`     yMMMMMMd  oMMMMMMs                                          \n");
	printf("                           sMMMh.       sMMMMMMM: `MMMMMMy                                          \n");
	printf("                            sNN.        dMMMMMMMd  mMMMMN.                                          \n");
	printf("                             ..         NMMMMMMmo .MMMNo:                                           \n");
	printf("                                        NMMMMMdy `hMMh:                                             \n");
	printf("                                        NMMMNo:msdMd/                                               \n");
	printf("                                 `-`   `MMMm+ydMMMy`                                                \n");
	printf("                                 `my`  sMMh+yys++smy/`                                              \n");
	printf("                             :o+.`oNd+hMMs        `ohd.                                             \n");
	printf("                             :+dmmhhMMMMN/-::-`      -                                              \n");
	printf("                               -/sdMNy+/+oyhs/s                                                     \n");
	printf("                              syss+-");
	printf("\n\n");*/
	printf("\t ############################\n");
	printf("\t #      BLACKHAWK v2.3      #\n");
	printf("\t #     HAWKING SPECTRUM     #\n");
	printf("\t #    COMPUTATION DEVICE    #\n");
	printf("\t ############################\n\n");
	fflush(stdout);
	
	if(argc < 2){
		printf("\t [main] : ERROR THIS PROGRAM NEEDS ONE PARAMETER\n");
		printf("\t [main] : PLEASE GIVE THE NAME OF THE PARAMETERS FILE\n\n");
		printf("[main] : END OF EXECUTION\n");
		fflush(stdout);
		exit(1);
	};
	
	printf("[main] : STARTING EXECUTION...\n");
	
	// ******READING THE RUN PARAMETERS*******
	struct param parameters;
	char param_name[128]; // name of the parameters file
	sprintf(param_name,"%s",argv[1]);
	
	printf("[main] : READING THE RUN PARAMETERS IN '%s'...",param_name);
	fflush(stdout);
	
	if(!read_params(&parameters,param_name,1)){ // if there was an error in the input parameters, we end the execution
		printf("[main] : END OF EXECUTION\n");
		fflush(stdout);
		return 1;
	}
	
	if(parameters.full_output){
		printf("\t DONE\n");
	}
	else{
		printf("\t\t DONE\n");
	}
	
	// *******ESTIMATION OF THE MEMORY USE*******
	// printf("[main] : ESTIMATION OF THE MEMORY USE...");
	// fflush(stdout);
	
	// if(!memory_estimation(&parameters,1)){ // if the user refuses the run
	// 	if(parameters.full_output){
	// 		printf("\n\t DONE\n");
	// 	}
	// 	else{
	// 		printf("\t\t DONE\n");
	// 	}
	// 	printf("[main] : END OF EXECUTION\n");
	// 	fflush(stdout);
	// 	return 1;
	// }
	
	if(parameters.full_output){
		printf("\t DONE\n");
	}
	else{
		printf("\t\t\t\t DONE\n");
	}
	
	// *******SAVING RUN PARAMETERS*******
	printf("[main] : SAVING RUN PARAMETERS...");
	fflush(stdout);
	
	DIR *directory;
	char dir_name[500];
	sprintf(dir_name,"results/%s",parameters.destination_folder);
	directory = opendir(dir_name);
	char create_dir[500];
	int overwrite = -1;
	if(directory){ // if the output directory selected already exists
		if(parameters.full_output){
			printf("\n\n\t Destination folder '%s/' already exists.\n",parameters.destination_folder);
			printf("\n\t Do you want to overwrite data? (type y or n) ");
			fflush(stdout);
			char answer;
			scanf(" %c",&answer);
			if(answer == 'y'){ // in this case data will be overwritten
				overwrite = 1;
				printf("\n");
				fflush(stdout);
			}
			else if(answer == 'n'){ // in this case we don't overwrite, execution will be stopped
				overwrite = 0;
				printf("\n\t DONE\n");
				fflush(stdout);
			}
			else{
				printf("\n[main] : ERROR WRONG ANSWER !\n");
				printf("[main] : END OF EXECUTION\n");
				fflush(stdout);
				exit(1);
			}
		}
		else{
			overwrite = 1;
		}
		closedir(directory);
	}
	
	switch(overwrite){
		case -1:{ // the output folder didn't exist, we create it
			chdir("results");
			sprintf(create_dir,"mkdir %s",parameters.destination_folder);
			system(create_dir);
			chdir("..");
			break;
		}
		case 0:{ // the output folder already existed, we don't overwrite data
			printf("[main] : END OF EXECUTION\n");
			fflush(stdout);
			return 1;
			break;
		}
		case 1:{ // the output folder already existed, we overwrite data
			if(parameters.full_output){
				printf("\t Data in '%s/' will be overwritten.\n\n",parameters.destination_folder);
			}
			else{
				printf("\n\n\t\t Data in '%s/' will be overwritten.\n",parameters.destination_folder);
			}
			fflush(stdout);
			break;
		}
		default :{
			printf("[main] : ERROR WRONG WRITING CHOICE !\n");
			fflush(stdout);
			exit(1);
			break;
		}
	}
	
	char new_parameters[500]; // we copy the parameters into the output directory for later interpretation
	sprintf(new_parameters,"cp %s ./results/%s/%s.txt",param_name,parameters.destination_folder,parameters.destination_folder);
	system(new_parameters);
	
	if(parameters.full_output){
		printf("\t DONE\n");
	}
	else{
		if(overwrite == -1){
			printf("\t\t\t\t\t DONE\n");
		}
		else{
			printf("\t\t\t\t\t\t\t\t\t DONE\n");
		}
	}
	fflush(stdout);
	
	// *******COMPUTING THE BH INSTANTANEOUS DISTRIBUTION*******
	double *init_masses = (double *)malloc(parameters.BH_number*sizeof(double)); // contains the masses of BHs
	double *init_params = (double *)malloc(parameters.param_number*sizeof(double)); // contains the spins of BHs
	double **spec_table = (double **)malloc(parameters.BH_number*sizeof(double *)); // contains the corresponding number codensity
	for(int i = 0;i<parameters.BH_number;i++){
		spec_table[i] = (double *)malloc(parameters.param_number*sizeof(double));
	}
	
	if(parameters.spectrum_choice != -1){ // in this case we use pre-built distributions
		printf("[main] : COMPUTING THE INITIAL DISTRIBUTION OF BLACK HOLES...");
		fflush(stdout);
		
		spectrum(init_masses,init_params,spec_table,&parameters);
		
		printf("\t\t DONE\n");
	}
	else{ // in this case we use the User's distribution
		printf("[main] : READING USER'S BLACK HOLES DISTRIBUTION TABLE...");
		fflush(stdout);
		
		read_users_table(init_masses,init_params,spec_table,&parameters);
		
		printf("\t\t DONE\n");
	}
	
	// *******WRITING INTO FILE "BH_spectrum.txt"*******
	printf("[main] : WRITING INTO FILE 'BH_spectrum.txt'...");
	fflush(stdout);
	
	write_spectrum(init_masses,init_params,spec_table,&parameters);
	
	printf("\t\t\t\t DONE\n");
	
	// *******READING GAMMA TABLES*******
	printf("[main] : READING GAMMA TABLES...");
	fflush(stdout);
	
	double ***gammas = (double ***)malloc(parameters.nb_gamma_spins*sizeof(double **)); // contains the tabulated greybody factors Gamma
	double *gamma_param = (double *)malloc(parameters.nb_gamma_param*sizeof(double)); // contains the tabulated energies
	double *gamma_x = (double *)malloc(parameters.nb_gamma_x*sizeof(double)); // contains the tabulated masses
	for(int i = 0;i<parameters.nb_gamma_spins;i++){
		gammas[i] = (double **)malloc(parameters.nb_gamma_param*sizeof(double *));
		for(int j = 0;j<parameters.nb_gamma_param;j++){
			gammas[i][j] = (double *)malloc(parameters.nb_gamma_x*sizeof(double));
		}
	}
	
	read_gamma_tables(gammas,gamma_param,gamma_x,&parameters);
	
	if(parameters.full_output){
		printf("\t DONE\n");
	}
	else{
		printf("\t\t\t\t\t DONE\n");
	}
	fflush(stdout);
	
	printf("[main] : READING FIT TABLES...");
	fflush(stdout);
	
	double ***fits = (double ***)malloc(parameters.nb_gamma_spins*sizeof(double **)); // contains the tabulated asymptotical Kerr BH emissivities
	for(int i = 0;i<parameters.nb_gamma_spins;i++){
		fits[i] = (double **)malloc(parameters.nb_gamma_param*sizeof(double *));
		for(int j = 0;j<parameters.nb_gamma_param;j++){
			fits[i][j] = (double *)malloc(parameters.nb_gamma_fits*sizeof(double));
		}
	}
	
	read_asymp_fits(fits,&parameters);
	
	if(parameters.full_output){
		printf("\t DONE\n");
	}
	else{
		printf("\t\t\t\t\t\t DONE\n");
	}
	fflush(stdout);
	
	// *******COMPUTING THE BH PRIMARY HAWKING SPECTRA*******
	printf("[main] : COMPUTING PRIMARY SPECTRA...");
	fflush(stdout);
	
	// Defining particle information
	double *dof = (double *)malloc((parameters.particle_number+parameters.grav+parameters.add_DM)*sizeof(double)); // contains the number of degrees of freedom of each particle
	dof[0] = 2.; // photons (2 helicity states)
	dof[5] = 6.; // neutrinos (3 particle species)*(1 helicity state)*(2 for antiparticles)
	dof[6] = 4.; // electrons (2 helicity states)*(2 for antiparticles)
	dof[7] = 4.; // muons (2 helicity states)*(2 for antiparticles)
	if(parameters.hadronization_choice == 3){ // in the Hazma case we consider that only photons, muons, electrons, neutrinos and pions are directly Hawking radiated at low energy
		dof[1] = 0.; // no gluons at low energy
		dof[2] = 0.; // no higgs boson at low energy
		dof[3] = 0.; // no W+- bosons at low energy
		dof[4] = 0.; // no Z0 boson at low energy
		dof[8] = 0.; // no taus at low energy
		dof[9] = 0.; // no up quark at low energy
		dof[10] = 0.; // no down quark at low energy
		dof[11] = 0.; // no charm quark at low energy
		dof[12] = 0.; // no strange quark at low energy
		dof[13] = 0.; // no top quark at low energy
		dof[14] = 0.; // no bottom quark at low energy
	}
	else{
		dof[1] = 16.; // gluons (8 particle species)*(2 helicity states)
		dof[2] = 1.; // higgs boson (1 helicity state)
		dof[3] = 6.; // W+- bosons (2 particle species)*(3 helicity states)
		dof[4] = 3.; // Z0 boson (3 helicity states)
		dof[8] = 4.; // taus (2 helicity states)*(2 for antiparticles)
		dof[9] = 12.; // up quark (2 helicity states)*(3 color states)*(2 for antiparticles)
		dof[10] = 12.; // down quark (2 helicity states)*(3 color states)*(2 for antiparticles)
		dof[11] = 12.; // charm quark (2 helicity states)*(3 color states)*(2 for antiparticles)
		dof[12] = 12.; // strange quark (2 helicity states)*(3 color states)*(2 for antiparticles)
		dof[13] = 12.; // top quark (2 helicity states)*(3 color states)*(2 for antiparticles)
		dof[14] = 12.; // bottom quark (2 helicity states)*(3 color states)*(2 for antiparticles)
	}
	if(parameters.grav){
		dof[15] = 2.; // graviton (2 polarization states)
		if(parameters.add_DM){
			dof[16] = parameters.dof_DM;
		}
	}
	else{
		if(parameters.add_DM){
			dof[15] = parameters.dof_DM;
		}
	}
		
	double *spins = (double *)malloc((parameters.particle_number+parameters.grav+parameters.add_DM)*sizeof(double)); // contains the spin of each particle
	spins[0] = 1.;
	spins[1] = 1.;
	spins[2] = 0.;
	spins[3] = 1.;
	spins[4] = 1.;
	spins[5] = 0.5;
	spins[6] = 0.5;
	spins[7] = 0.5;
	spins[8] = 0.5;
	spins[9] = 0.5;
	spins[10] = 0.5;
	spins[11] = 0.5;
	spins[12] = 0.5;
	spins[13] = 0.5;
	spins[14] = 0.5;
	if(parameters.grav){
		spins[15] = 2.;
		if(parameters.add_DM){
			spins[16] = parameters.spin_DM;
		}
	}
	else{
		if(parameters.add_DM){
			spins[15] = parameters.spin_DM;
		}
	}
	
	double *masses_primary = (double *)malloc((parameters.particle_number+parameters.grav+parameters.add_DM)*sizeof(double)); // contains the mass of each particle
	masses_primary[0] = m_photon;
	masses_primary[1] = m_gluon;
	masses_primary[2] = m_higgs;
	masses_primary[3] = m_wpm;
	masses_primary[4] = m_z0;
	masses_primary[5] = m_neutrino;
	masses_primary[6] = m_electron;
	masses_primary[7] = m_muon;
	masses_primary[8] = m_tau;
	masses_primary[9] = m_up;
	masses_primary[10] = m_down;
	masses_primary[11] = m_charm;
	masses_primary[12] = m_strange;
	masses_primary[13] = m_top;
	masses_primary[14] = m_bottom;
	if(parameters.grav){
		masses_primary[15] = m_graviton;
		if(parameters.add_DM){
			masses_primary[16] = parameters.m_DM;
		}
	}
	else{
		if(parameters.add_DM){
			masses_primary[15] = parameters.m_DM;
		}
	}
	
	double **instantaneous_primary_spectra = (double **)malloc((parameters.particle_number+parameters.grav+parameters.add_DM)*sizeof(double *)); // contains the Hawking spectra of each primary particle
	for(int i = 0;i<parameters.particle_number+parameters.grav+parameters.add_DM;i++){
		instantaneous_primary_spectra[i] = (double *)malloc(parameters.E_number*sizeof(double));
		for(int j = 0;j<parameters.E_number;j++){
			instantaneous_primary_spectra[i][j] = 0.;
		}
	}
	
	double *energies = (double *)malloc(parameters.E_number*sizeof(double)); // contains the initial energies
	for(int i = 0;i<parameters.E_number;i++){ // logarithmic distribution
		energies[i] = pow(10.,log10(parameters.Emin) + i*(log10(parameters.Emax) - log10(parameters.Emin))/(parameters.E_number-1));
	}
	
	double **BH_masses = (double **)malloc(parameters.BH_number*sizeof(double *)); // contains the values of instantaneous BH masses
	double **BH_params = (double **)malloc(parameters.BH_number*sizeof(double *)); // contains the values of instantaneous BH spins
	for(int i = 0;i<parameters.BH_number;i++){
		BH_masses[i] = (double *)malloc(parameters.param_number*sizeof(double));
		BH_params[i] = (double *)malloc(parameters.param_number*sizeof(double));
		for(int j = 0;j<parameters.param_number;j++){
			BH_masses[i][j] = init_masses[i];
			BH_params[i][j] = init_params[j];
		}
	}
	
	int *compute_primary = (int *)malloc((parameters.particle_number+parameters.grav+parameters.add_DM)*sizeof(int)); // this table decides whether the corresponding spectra are computed or not
	compute_primary[0] = 1; // photons
	compute_primary[1] = 1; // gluons
	compute_primary[2] = 1; // higgs boson
	compute_primary[3] = 1; // W+- boson
	compute_primary[4] = 1; // Z0 boson
	compute_primary[5] = 1; // neutrinos
	compute_primary[6] = 1; // electron
	compute_primary[7] = 1; // muon
	compute_primary[8] = 1; // tau
	compute_primary[9] = 1; // up quark
	compute_primary[10] = 1; // down quark
	compute_primary[11] = 1; // charm quark
	compute_primary[12] = 1; // strange quark
	compute_primary[13] = 1; // top quark
	compute_primary[14] = 1; // bottom quark
	if(parameters.grav){
		compute_primary[15] = 1; // graviton
		if(parameters.add_DM){
			compute_primary[16] = 1; // DM
		}
	}
	else if(parameters.add_DM){
		compute_primary[15] = 1; // DM
	}
	
	int **counters_param = (int **)malloc(parameters.BH_number*sizeof(int *)); // contains the running position in the tabulated parameters
	int ***counters_x = (int ***)malloc(parameters.BH_number*sizeof(int **)); // contains the running position in the tabulated energies
	for(int i = 0;i<parameters.BH_number;i++){
		counters_param[i] = (int *)malloc(parameters.param_number*sizeof(int));
		counters_x[i] = (int **)malloc(parameters.param_number*sizeof(int *));
		for(int j = 0;j<parameters.param_number;j++){
			counters_x[i][j] = (int *)malloc(parameters.E_number*sizeof(int));
		}
	}
	
	instantaneous_primary_spectrum(instantaneous_primary_spectra,BH_masses,BH_params,spec_table,energies,gammas,gamma_param,gamma_x,fits,dof,spins,masses_primary,counters_param,counters_x,compute_primary,&parameters);
	
	if(parameters.full_output){
		printf("\t DONE\n");
	}
	else{
		printf("\t\t\t\t\t DONE\n");
	}
	
	free2D_double(BH_masses,parameters.BH_number);
	free2D_double(BH_params,parameters.BH_number);
	free2D_double(spec_table,parameters.BH_number);
	free3D_double(gammas,parameters.nb_gamma_spins,parameters.nb_gamma_param);
	free1D_double(gamma_param);
	free1D_double(gamma_x);
	free3D_double(fits,parameters.nb_gamma_spins,parameters.nb_gamma_param);
	free1D_double(dof); // freeing particle informations
	free1D_double(spins);
	free1D_int(compute_primary);
	free2D_int(counters_param,parameters.BH_number);
	free3D_int(counters_x,parameters.BH_number,parameters.param_number);
	
	// *******WRITING INTO FILE "instantaneous_primary_spectra.txt"*******
	printf("[main] : WRITING INTO FILE 'instantaneous_primary_spectra.txt'...");
	fflush(stdout);
	
	write_instantaneous_primary_spectra(instantaneous_primary_spectra,energies,&parameters);
	
	printf("\t DONE\n");
	fflush(stdout);
	
	if(parameters.primary_only == 0){
		
		// *******READING THE HADRONIZATION TABLES IN FOLDERS "*_tables/"*******
		printf("[main] : READING HADRONIZATION TABLES...");
		fflush(stdout);
		
		double ****tables = (double ****)malloc(parameters.nb_fin_part*sizeof(double ***)); // contains the branching ratios dN(E',E) of the hadronization tables (PYTHIA, HERWIG) of dN/dE(E',E) (HAZMA)
		for(int i = 0;i<parameters.nb_fin_part;i++){
			tables[i] = (double ***)malloc(parameters.nb_init_en*sizeof(double **));
			for(int j = 0;j<parameters.nb_init_en;j++){
				tables[i][j] = (double **)malloc(parameters.nb_fin_en*sizeof(double *));
				for(int k = 0;k<parameters.nb_fin_en;k++){
					tables[i][j][k] = (double *)malloc(parameters.nb_init_part*sizeof(double));
				}
			}
		}
		double *initial_energies = (double *)malloc(parameters.nb_init_en*sizeof(double)); // contains the tabulated initial energies
		double *final_energies = (double *)malloc(parameters.nb_fin_en*sizeof(double)); // contains the tabulated final energies
		
		read_hadronization_tables(tables,initial_energies,final_energies,&parameters);
		
		if(parameters.full_output){
			printf("\t DONE\n");
		}
		else{
			printf("\t\t\t\t DONE\n");
		}
		
		// *******HADRONIZING THE INSTANTANEOUS PRIMARY SPECTRA*******			
		printf("[main] : HADRONIZING PARTICLES...");
		fflush(stdout);
		
		double *masses_secondary = (double *)malloc(parameters.nb_fin_part*sizeof(double)); // secondary particles informations
		switch(parameters.hadronization_choice){
			case 0: case 1:{ // BBN epoch particles
				masses_secondary[0] = m_photon; // photon
				masses_secondary[1] = m_electron; // electron
				masses_secondary[2] = m_muon; // muon 
				masses_secondary[3] = m_neutrino; // nu_e
				masses_secondary[4] = m_neutrino; // nu_mu 
				masses_secondary[5] = m_neutrino; // nu_tau
				masses_secondary[6] = m_pipm; // pi+-
				masses_secondary[7] = m_K0L; // K0L
				masses_secondary[8] = m_Kpm; // K+-
				masses_secondary[9] = m_proton; // proton
				masses_secondary[10] = m_neutron; // neutron
				break;
			}
			case 2: case 4:{ // present epoch particles
				masses_secondary[0] = m_photon; // photon
				masses_secondary[1] = m_electron; // electron
				masses_secondary[2] = m_neutrino; // nu_e
				masses_secondary[3] = m_neutrino; // nu_mu
				masses_secondary[4] = m_neutrino; // nu_tau
				masses_secondary[5] = m_proton; // proton
				break;
			}
			case 3:{ // Hazma particles
				masses_secondary[0] = m_photon; // photon
				masses_secondary[1] = m_electron; // electron
				break;
			}
			default:{
				printf("[main] : ERROR WRONG HADRONIZATION CHOICE !\n");
				fflush(stdout);
				exit(1);
				break;
			}
		}
		
		double ***instantaneous_hadronized_spectra = (double ***)malloc(parameters.nb_fin_part*sizeof(double **)); // contains the primary spectra multiplied by the branching ratios for each primary particle
		for(int i = 0;i<parameters.nb_fin_part;i++){
			instantaneous_hadronized_spectra[i] = (double **)malloc(parameters.E_number*sizeof(double *));
			for(int j = 0;j<parameters.E_number;j++){
				instantaneous_hadronized_spectra[i][j] = (double *)malloc(parameters.nb_fin_en*sizeof(double));
				for(int k = 0;k<parameters.nb_fin_en;k++){
					instantaneous_hadronized_spectra[i][j][k] = 0.;
				}
			}
		}
		
		int *compute = (int *)malloc(parameters.nb_fin_part*sizeof(int)); // decides whether the particle contribution is computed or not
		if(parameters.hadronization_choice == 0 || parameters.hadronization_choice == 1){ // BBN epoch secondary particles
			compute[0] = 1; // photon
			compute[1] = 1; // electron
			compute[2] = 1; // muon
			compute[3] = 1; // nu_e
			compute[4] = 1; // nu_mu
			compute[5] = 1; // nu_tau
			compute[6] = 1; // pi+-
			compute[7] = 1; // K0 long
			compute[8] = 1; // K+-
			compute[9] = 1; // proton
			compute[10] = 1; // neutron
		}
		else if(parameters.hadronization_choice == 2 || parameters.hadronization_choice == 4){ // today epoch secondary particles
			compute[0] = 1; // photon
			compute[1] = 1; // electron
			compute[2] = 1; // nu_e
			compute[3] = 1; // nu_mu
			compute[4] = 1; // nu_tau
			compute[5] = 1; // proton
		}
		else{ // Hazma secondary particles
			compute[0] = 1; // photon
			compute[1] = 1; // electron
		}
		
		hadronize_instantaneous(instantaneous_hadronized_spectra,tables,initial_energies,final_energies,instantaneous_primary_spectra,energies,masses_secondary,compute,&parameters);
		
		if(parameters.full_output){
			printf("\t DONE\n");
		}
		else{
			printf("\t\t\t\t\t DONE\n");
		}
		
		free4D_double(tables,parameters.nb_fin_part,parameters.nb_init_en,parameters.nb_fin_en); // freeing arrays not used anymore
		free1D_double(masses_secondary);
		
		// *******INTEGRATING THE INITIAL ENERGIES*******
		printf("[main] : INTEGRATING OVER INITIAL ENERGIES...");
		fflush(stdout);
					
		double **instantaneous_integrated_hadronized_spectra = (double **)malloc(parameters.nb_fin_part*sizeof(double *)); // contains the secondary spectra for each secondary particle
		for(int i = 0;i<parameters.nb_fin_part;i++){
			instantaneous_integrated_hadronized_spectra[i] = (double *)malloc(parameters.nb_fin_en*sizeof(double));
		}
		
		integrate_initial_energies_instantaneous(instantaneous_hadronized_spectra,instantaneous_integrated_hadronized_spectra,energies,final_energies,compute,&parameters);
		
		if(compute[0]){
			add_photons_instantaneous(instantaneous_primary_spectra,instantaneous_integrated_hadronized_spectra,energies,final_energies,&parameters);
			if(parameters.hadronization_choice == 3){
				add_FSR_instantaneous(instantaneous_primary_spectra,instantaneous_integrated_hadronized_spectra,energies,final_energies,masses_primary,&parameters);
			}
		}
		if(((parameters.hadronization_choice == 0 || parameters.hadronization_choice == 1) && compute[3]+compute[4]+compute[5]) || ((parameters.hadronization_choice == 2 || parameters.hadronization_choice == 4) && compute[2]+compute[3]+compute[4])){
			add_neutrinos_instantaneous(instantaneous_primary_spectra,instantaneous_integrated_hadronized_spectra,energies,final_energies,&parameters);
		}
		if(compute[1]){
			add_electrons_instantaneous(instantaneous_primary_spectra,instantaneous_integrated_hadronized_spectra,energies,final_energies,&parameters);
		}
		
		if(parameters.full_output){
			printf("\t DONE\n");
		}
		else{
			printf("\t\t\t\t DONE\n");
		}
		fflush(stdout);
		
		free2D_double(instantaneous_primary_spectra,parameters.particle_number+parameters.grav+parameters.add_DM);
		free1D_double(initial_energies); // freeing arrays not used anymore
		free3D_double(instantaneous_hadronized_spectra,parameters.nb_fin_part,parameters.E_number);
		free1D_double(energies);
		free1D_double(masses_primary);
		free1D_int(compute);
		
		// *******WRITING INTO FILE "instantaneous_secondary_spectra.txt"*******
		printf("[main] : WRITING INTO FILE 'instantaneous_secondary_spectra.txt'...");
		fflush(stdout);
		
		write_instantaneous_secondary_spectra(instantaneous_integrated_hadronized_spectra,final_energies,&parameters);
		
		printf("\t DONE\n");
		fflush(stdout);
		
		free1D_double(final_energies); // freeing last arrays
		free2D_double(instantaneous_integrated_hadronized_spectra,parameters.nb_fin_part);
	}
	
	printf("[main] : END OF EXECUTION\n");
	fflush(stdout);
	
	return 0; // this is the only possibility of returning 0: it means that all steps have succeeded
}
// &&&&&&&&&&&&&&&&&&&&&&&&& END OF MAIN &&&&&&&&&&&&&&&&&&&&&&&&&&&&&
