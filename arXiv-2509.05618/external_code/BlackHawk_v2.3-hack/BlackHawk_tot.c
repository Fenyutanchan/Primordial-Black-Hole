// Copyright (c) 2025 Quan-feng WU <wuquanfeng@ihep.ac.cn>
// 
// This software is released under the MIT License.
// https://opensource.org/licenses/MIT
// 
// Modified from `BlackHawk_tot.c` in BlackHawk v2.3
// Date: 2025-05-20
// ============================================================================
// 
// BlackHawk_tot is a program aimed at computing the Hawking
// emission of a distribution of black holes through time.
// Last modification: 15 September 2023
// Authors: Jérémy Auffinger jeremy.auffinger.blackhawk@protonmail.com & Alexandre Arbey alexandre.arbey@ens-lyon.fr

#include "src/include.h"

// &&&&&&&&&&&&&&&&&&&&& MAIN &&&&&&&&&&&&&&&&&&&&&&&&
// This is the BlackHawk_tot main routine. It executes all intermediate steps
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
		exit(0);
	};
	
	printf("[main] : STARTING EXECUTION...\n");
	
	// ******READING THE RUN PARAMETERS*******
	struct param parameters;
	char param_name[128]; // name of the parameters file
	sprintf(param_name,"%s",argv[1]);
	
	printf("[main] : READING THE RUN PARAMETERS IN '%s'...",param_name);
	fflush(stdout);
	
	if(!read_params(&parameters,param_name,0)){ // if there was an error in the input parameters, we end the execution
		printf("[main] : END OF EXECUTION\n");
		fflush(stdout);
		return 1;
	}
	
	printf("\t DONE\n");
	
	// *******ESTIMATION OF THE MEMORY USE*******
	// printf("[main] : ESTIMATION OF THE MEMORY USE...");
	// fflush(stdout);
	
	// if(!memory_estimation(&parameters,0)){ // if the user refuses the run
	// 	printf("\n\t DONE\n");
	// 	printf("[main] : END OF EXECUTION\n");
	// 	fflush(stdout);
	// 	return 1;
	// }
	
	if(parameters.full_output){
		printf("\t DONE\n");
	}
	else{
		printf("\t\t\t DONE\n");
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
			char answer = 'y';
			// scanf(" %c",&answer);
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
				exit(0);
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
				printf("\n\n\t Data in '%s/' will be overwritten.\n\n",parameters.destination_folder);
			}
			fflush(stdout);
			break;
		}
		default :{
			printf("[main] : ERROR WRONG WRITING CHOICE !\n");
			fflush(stdout);
			exit(0);
			break;
		}
	}
	
	char new_parameters[256]; // we copy the parameters into the output directory for later interpretation
	sprintf(new_parameters,"cp %s ./results/%s/%s.txt",param_name,parameters.destination_folder,parameters.destination_folder);
	system(new_parameters);
	
	if(parameters.full_output){
		printf("\t DONE\n");
	}
	else{
		if(overwrite == -1){
			printf("\t\t\t\t DONE\n");
		}
		else{
			printf("\t\t\t\t\t\t\t\t DONE\n");
		}
	}
	fflush(stdout);
	
	// *******COMPUTATION OR READING OF THE MASS SPECTRUM OF BH*******
	double *init_masses = (double *)malloc(parameters.BH_number*sizeof(double)); // contains the list of BH initial masses
	double *init_params = (double *)malloc(parameters.param_number*sizeof(double)); // contains the list of BH initial spins
	double **spec_table = (double **)malloc(parameters.BH_number*sizeof(double *)); // contains the comoving number density for each mass
	for(int i = 0;i<parameters.BH_number;i++){
		spec_table[i] = (double *)malloc(parameters.param_number*sizeof(double));
	}
	
	if(parameters.spectrum_choice != -1){ // in this case we use pre-built distributions
		printf("[main] : COMPUTING THE INITIAL DISTRIBUTION OF BLACK HOLES...");
		fflush(stdout);
		
		spectrum(init_masses,init_params,spec_table,&parameters);
	}
	else{ // in this case we use the User's distribution
		printf("[main] : READING USER'S BLACK HOLES DISTRIBUTION TABLE...");
		fflush(stdout);
		
		read_users_table(init_masses,init_params,spec_table,&parameters);
	}
	
	printf("\t DONE\n");
	
	// *******WRITING INTO FILE "BH_spectrum.txt"*******
	printf("[main] : WRITING INTO FILE 'BH_spectrum.txt'...");
	fflush(stdout);
	
	write_spectrum(init_masses,init_params,spec_table,&parameters);
	
	printf("\t\t\t DONE\n");
	
	// *******READING f/g TABLES in files './fM_tables/*.txt'*******
	printf("[main] : READING EVOLUTION TABLES...");
	fflush(stdout);
	
	double **fM_table = (double **)malloc(parameters.nb_fM_masses*sizeof(double *)); // contains the f(M,param) values
	double **gM_table;
	double *fM_masses = (double *)malloc(parameters.nb_fM_masses*sizeof(double)); // contains the tabulated masses
	double *fM_param = (double *)malloc(parameters.nb_fM_param*sizeof(double)); // contains the tabulated spins
	if(parameters.metric == 0){
		gM_table = (double **)malloc(parameters.nb_fM_masses*sizeof(double *)); // contains the g(M,a) values
	}
	for(int i = 0;i<parameters.nb_fM_masses;i++){
		fM_table[i] = (double *)malloc(parameters.nb_fM_param*sizeof(double));
		if(parameters.metric == 0){
			gM_table[i] = (double *)malloc(parameters.nb_fM_param*sizeof(double));
		}
	}
	
	read_fM_table(fM_table,fM_masses,fM_param,&parameters);
	
	if(parameters.metric == 0){
		read_gM_table(gM_table,fM_masses,fM_param,&parameters);
	}
	
	if(parameters.full_output){
		printf("\t DONE\n");
	}
	else{
		printf("\t\t\t\t DONE\n");
	}
	fflush(stdout);
	
	// *******COMPUTING THE EVOLUTION OF PBHs*******
	printf("[main] : COMPUTING THE EVOLUTION OF BLACK HOLES...");
	fflush(stdout);
	
	double **evol_times;
	double *sorted_times;
	int **rank = (int **)malloc(parameters.BH_number*sizeof(int *)); // contains the rank of BHs in the sorted lifetimes table
	for(int i = 0;i<parameters.BH_number;i++){
		rank[i] = (int *)malloc(parameters.param_number*sizeof(int));
	}
	evol_times = (double **)malloc(parameters.BH_number*sizeof(double *)); // contains the BH lifetimes
	for(int i = 0;i<parameters.BH_number;i++){
		evol_times[i] = (double *)malloc(parameters.param_number*sizeof(double));
	}

	evolution_times(init_masses,init_params,evol_times,fM_table,gM_table,fM_masses,fM_param,&parameters);

	sorted_times = (double *)malloc(parameters.BH_number*parameters.param_number*sizeof(double)); // contains the sorted BH lifetimes
	
	sort_lifetimes(evol_times,sorted_times,rank,&parameters);
	
	double ***life_masses = (double ***)malloc(parameters.BH_number*sizeof(double **)); // contains the list of BH masses as a function of time
	double ***life_params = (double ***)malloc(parameters.BH_number*sizeof(double **)); // contains the list of BH spins as a function of time
	double *life_times = (double *)malloc(parameters.limit*sizeof(double)); // contains the corresponding time steps
	double *dts = (double *)malloc(parameters.limit*sizeof(double)); // contains the explicit time intervals to avoid null differences between close time steps
	for(int i = 0;i<parameters.BH_number;i++){
		life_masses[i] = (double **)malloc(parameters.param_number*sizeof(double *));
		life_params[i] = (double **)malloc(parameters.param_number*sizeof(double *));
		for(int j = 0;j<parameters.param_number;j++){
			life_masses[i][j] = (double *)malloc(parameters.limit*sizeof(double));
			life_params[i][j] = (double *)malloc(parameters.limit*sizeof(double));
		}
	}
	
	life_evolution(life_masses,life_params,life_times,dts,init_masses,init_params,rank,fM_table,gM_table,fM_masses,fM_param,&parameters);
	
	if(parameters.full_output){
		printf("\t DONE\n");
	}
	else{
		printf("\t\t DONE\n");
	}
	
	free1D_double(init_masses); // freeing arrays not used anymore
	free1D_double(init_params);
	free2D_double(fM_table,parameters.nb_fM_masses);
	free1D_double(fM_param);
	if(parameters.metric == 0){
		free1D_double(sorted_times);
		free2D_double(gM_table,parameters.nb_fM_masses);
		free2D_double(evol_times,parameters.BH_number);
	}
	free1D_double(fM_masses);
	free2D_int(rank,parameters.BH_number);
	
	// *******WRITING INTO FILE "life_evolutions.txt"*******
	printf("[main] : WRITING INTO FILE 'life_evolutions.txt'...");
	fflush(stdout);
	
	write_life_evolutions(life_masses,life_params,life_times,dts,&parameters);
	
	free1D_double(dts);
	
	printf("\t\t DONE\n");
	fflush(stdout);
	
	// *******READING GAMMA TABLES*******
	printf("[main] : READING GAMMA TABLES...");
	fflush(stdout);
	
	double ***gammas = (double ***)malloc(parameters.nb_gamma_spins*sizeof(double **)); // contains the tabulated greybody factors Gamma
	double *gamma_param = (double *)malloc(parameters.nb_gamma_param*sizeof(double)); // contains the tabulated parameters
	double *gamma_x = (double *)malloc(parameters.nb_gamma_x*sizeof(double)); // contains the tabulated energies
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
		printf("\t\t\t\t DONE\n");
	}
	fflush(stdout);
	
	printf("[main] : READING FIT TABLES...");
	fflush(stdout);
	
	double ***fits = (double ***)malloc(parameters.nb_gamma_spins*sizeof(double **)); // contains the tabulated asymptotical fits for Kerr BH emissivities
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
		printf("\t\t\t\t\t DONE\n");
	}
	fflush(stdout);
	
	// Defining particle information
	double *dof = (double *)malloc((parameters.particle_number+parameters.grav+parameters.add_DM)*sizeof(double)); // contains the number of degrees of freedom for each particle
	dof[0] = 2.; // photons (2 helicity states)
	dof[1] = 16.; // gluons (8 particle species)*(2 helicity states)
	dof[2] = 1.; // higgs boson (1 helicity state)
	dof[3] = 6.; // W+- bosons (2 particle species)*(3 helicity states)
	dof[4] = 3.; // Z0 bosons (3 helicity states)
	dof[5] = 6.; // neutrinos (3 particle species)*(1 helicity state)*(2 for antiparticles)
	dof[6] = 4.; // electrons (2 helicity states)*(2 for antiparticles)
	dof[7] = 4.; // muons (2 helicity states)*(2 for antiparticles)
	dof[8] = 4.; // taus (2 helicity states)*(2 for antiparticles)
	dof[9] = 12.; // up quark (2 helicity states)*(3 color states)*(2 for antiparticles)
	dof[10] = 12.; // down quark (2 helicity states)*(3 color states)*(2 for antiparticles)
	dof[11] = 12.; // charm quark (2 helicity states)*(3 color states)*(2 for antiparticles)
	dof[12] = 12.; // strange quark (2 helicity states)*(3 color states)*(2 for antiparticles)
	dof[13] = 12.; // top quark (2 helicity states)*(3 color states)*(2 for antiparticles)
	dof[14] = 12.; // bottom quark (2 helicity states)*(3 color states)*(2 for antiparticles)
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
	
	double *times = (double *)malloc(parameters.nb_fin_times*sizeof(double)); // contains the final evolution times
	for(int i = 0;i<parameters.nb_fin_times;i++){
		times[i] = life_times[i];
	}
	free1D_double(life_times);
				
	double *energies = (double *)malloc(parameters.E_number*sizeof(double)); // contains the initial energies
	for(int i = 0;i<parameters.E_number;i++){
		energies[i] = pow(10.,log10(parameters.Emin) + i*(log10(parameters.Emax) - log10(parameters.Emin))/(parameters.E_number-1)); // initial energies are distributed logarithmically
	}
	
	double ****tables; // contains the tabulated branching ratios dN(E',E) of the hadronization tables (PYTHIA, HERWIG) of dN/dE(E',E) (HAZMA)
	double *initial_energies; // contains the tabulated initial energies
	double *final_energies; // contains the tabulated final energies
	
	if(parameters.primary_only == 0){
		// *******READING THE HADRONIZATION TABLES*******
		printf("[main] : READING HADRONIZATION TABLES...");
		fflush(stdout);
		
		tables = (double ****)malloc(parameters.nb_fin_part*sizeof(double ***));
		for(int i = 0;i<parameters.nb_fin_part;i++){
			tables[i] = (double ***)malloc(parameters.nb_init_en*sizeof(double **));
			for(int j = 0;j<parameters.nb_init_en;j++){
				tables[i][j] = (double **)malloc(parameters.nb_fin_en*sizeof(double *));
				for(int k = 0;k<parameters.nb_fin_en;k++){
					tables[i][j][k] = (double *)malloc(parameters.nb_init_part*sizeof(double));
				}
			}
		}
		initial_energies = (double *)malloc(parameters.nb_init_en*sizeof(double));
		final_energies = (double *)malloc(parameters.nb_fin_en*sizeof(double));
		
		read_hadronization_tables(tables,initial_energies,final_energies,&parameters);
		
		if(parameters.full_output){
			printf("\t DONE\n");
		}
		else{
			printf("\t\t\t DONE\n");
		}
	}
		
	// *******COMPUTATION OF PRIMARY AND SECONDARY SPECTRA*******
	printf("[main] : COMPUTING SPECTRA...");
	fflush(stdout);
	
	double **partial_primary_spectra = (double **)malloc((parameters.particle_number + parameters.grav+parameters.add_DM)*sizeof(double *)); // contains the snapshot of the primary spectra at each time
	for(int i = 0;i<parameters.particle_number+parameters.grav+parameters.add_DM;i++){
		partial_primary_spectra[i] = (double *)malloc(parameters.E_number*sizeof(double));
	}
	
	double ***partial_hadronized_spectra; // contains the snapshot of the primary spectra multiplied by the branching ratios
	double **partial_integrated_hadronized_spectra; // contains the snapshot of the integrated hadronized spectra at each time
	double *masses_secondary = (double *)malloc(parameters.nb_fin_part*sizeof(double)); // secondary particles informations
	if(parameters.primary_only == 0){
		partial_hadronized_spectra = (double ***)malloc(parameters.nb_fin_part*sizeof(double **));
		for(int i = 0;i<parameters.nb_fin_part;i++){
			partial_hadronized_spectra[i] = (double **)malloc(parameters.E_number*sizeof(double *));
			for(int j = 0;j<parameters.E_number;j++){
				partial_hadronized_spectra[i][j] = (double *)malloc(parameters.nb_fin_en*sizeof(double));
				for(int k = 0;k<parameters.nb_fin_en;k++){
					partial_hadronized_spectra[i][j][k] = 0.;
				}
			}
		}
		partial_integrated_hadronized_spectra = (double **)malloc(parameters.nb_fin_part*sizeof(double *));
		for(int i = 0;i<parameters.nb_fin_part;i++){
			partial_integrated_hadronized_spectra[i] = (double *)malloc(parameters.nb_fin_en*sizeof(double));
		}
		
		if(parameters.hadronization_choice == 0 || parameters.hadronization_choice == 1){ // BBN epoch particles
			masses_secondary[0] = m_photon;
			masses_secondary[1] = m_electron;
			masses_secondary[2] = m_muon;
			masses_secondary[3] = m_neutrino;
			masses_secondary[4] = m_neutrino;
			masses_secondary[5] = m_neutrino;
			masses_secondary[6] = m_pipm;
			masses_secondary[7] = m_K0L;
			masses_secondary[8] = m_Kpm;
			masses_secondary[9] = m_proton;
			masses_secondary[10] = m_neutron;
		}
		else if(parameters.hadronization_choice == 2 || parameters.hadronization_choice == 4){ // present epoch particles
			masses_secondary[0] = m_photon;
			masses_secondary[1] = m_electron;
			masses_secondary[2] = m_neutrino;
			masses_secondary[3] = m_neutrino;
			masses_secondary[4] = m_neutrino;
			masses_secondary[5] = m_proton;
		}
		else{ // Hazma particles
			masses_secondary[0] = m_photon;
			masses_secondary[1] = m_electron;
		}
	}
	
	double **BH_masses = (double **)malloc(parameters.BH_number*sizeof(double *)); // contains the instantaneous values of the BH masses
	double **BH_params = (double **)malloc(parameters.BH_number*sizeof(double *)); // contains the instantaneous values of the BH spins
	for(int i = 0;i<parameters.BH_number;i++){
		BH_masses[i] = (double *)malloc(parameters.param_number*sizeof(double));
		BH_params[i] = (double *)malloc(parameters.param_number*sizeof(double));
	}
	
	total_spectra(partial_hadronized_spectra,partial_primary_spectra,partial_integrated_hadronized_spectra,tables,initial_energies,final_energies,spec_table,times,life_masses,BH_masses,life_params,BH_params,energies,masses_primary,spins,dof,masses_secondary,gammas,gamma_param,gamma_x,fits,&parameters);
	
	if(parameters.full_output){
		printf("\t DONE\n");
	}
	else{
		printf("\t\t\t\t\t DONE\n");
	}
	
	free3D_double(gammas,5,parameters.nb_gamma_param); // freeing last arrays
	free1D_double(gamma_param);
	free1D_double(gamma_x);
	free3D_double(fits,5,parameters.nb_gamma_param);
	free1D_double(times);
	free1D_double(energies);
	free2D_double(BH_masses,parameters.BH_number);
	free2D_double(BH_params,parameters.BH_number);
	free1D_double(masses_primary);
	free1D_double(spins);
	free1D_double(dof);
	free2D_double(partial_primary_spectra,parameters.particle_number+parameters.grav+parameters.add_DM);
	if(parameters.primary_only == 0){
		free4D_double(tables,parameters.nb_fin_part,parameters.nb_init_en,parameters.nb_fin_en);
		free1D_double(final_energies);
		free1D_double(initial_energies);
		free2D_double(partial_integrated_hadronized_spectra,parameters.nb_fin_part);
		free3D_double(partial_hadronized_spectra,parameters.nb_fin_part,parameters.E_number);
		free1D_double(masses_secondary);
	}
	
	printf("[main] : END OF EXECUTION\n");
	fflush(stdout);
	
	return 0; // this is the only possibility of returning 0: it means that all steps have succeeded
}
// &&&&&&&&&&&&&&&&&&&&&&&&& END OF MAIN &&&&&&&&&&&&&&&&&&&&&&&&&&&&&
