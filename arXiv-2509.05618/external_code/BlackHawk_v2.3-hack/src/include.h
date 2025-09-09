// Header file for all functions in BlackHawk.c
// Last modification: 07 April 2024
// Author: Alexandre Arbey alexandre.arbey@ens-lyon.fr & Jérémy Auffinger jeremy.auffinger.blackhawk@protonmail.com

#include <stdbool.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <dirent.h>
#include <sys/stat.h>
#include <complex.h>
#include <strings.h>
#include <unistd.h>

//#define HARDTABLES

// unit conversions, CGS to GeV
#define mass_conversion  5.60958884e+23 // g to GeV
#define time_conversion  1.519267407e+24 // s to GeV^(-1)
#define leng_conversion  5.06773058e+13 // cm to GeV^(-1)
#define rate_conversion  (1.e-100*time_conversion*pow(leng_conversion,3.)) // cm^(-3)s^(-1)GeV^(-1) to GeV^3, plus the correction of 1.e-100
#define dens_conversion  (1.e-100*pow(leng_conversion,3.)) // cm^(-3) to GeV^3, plus the correction of 1.e-100

// some constants
#define pi  				3.141592653589795
#define machine_precision  	1.e-10 // this is the precision used to decide whether #define numbers are equal or not
#define G  					(6.67408e-11*pow(leng_conversion*100.,3.)/(mass_conversion*1000.)/pow(time_conversion,2.)) // Newton constant in GeV
#define Mp  				1.221e+19 // Planck mass in GeV				
#define alpha_em			1./137. // dimensionless electroweak coupling constant

// masses of the Standard Model particles (PDG 2017) in GeV
#define m_electron  0.5109989461e-3
#define m_neutrino  0. // simplifying hypothesis, considering that even the heaviest primordial black hole has a temperature >> m_neutrinos
#define m_muon  	105.6583745e-3
#define m_tau  		1.77686
#define m_up 		2.2e-3
#define m_down  	4.7e-3
#define m_strange  	96e-3
#define m_charm  	1.27
#define m_top  		173.21
#define m_bottom  	4.18
#define m_higgs  	125.03
#define m_photon   	0.
#define m_gluon   	0.2 // effective mass
#define m_wpm  		80.403
#define m_z0  		91.1876
#define m_graviton  0.

#define m_pi0		134.9766e-3
#define m_pipm  	139.57018e-3
#define m_K0L   	497.7e-3
#define m_Kpm   	493.7e-3
#define m_proton    938.272e-3
#define m_neutron   939.5654e-3

// structure
struct param{ // this structure will contain all the parameters of the run
	char destination_folder[200]; // destination folder of the output
	int full_output; // 0: reduced output, 1: extensive output
	int interpolation_method; // 0: we perform linear interpolations, 1: we perform logarithmic interpolations
	
	int metric; // 0: Kerr BHs, 1: polymerized BHs, 2: charged BHs, 3: higher-dimensional BHs
	
	int BH_number; // number of initial black hole masses
	double Mmin; // minimum initial mass of black holes
	double Mmax; // maximum initial mass of black holes
	int param_number; // number of black holes parameters
	double amin; // minimum initial black hole spin
	double amax; // maximum initial black hole spin
	double Qmin; // minimum initial black hole charge
	double Qmax; // maximum initial black hole charge
	double epsilon_LQG; // the epsilon parameter for the polymerized metric (epsilon < 0.794)
	double a0_LQG; // the a0 parameter for the polymerized metric (a0 = 0. or 0.11)
	double n; // the number of additional spatial dimensions
	double M_star; // modified Planck mass in extra dimension models
	int spectrum_choice; // chooses the initial black holes mass distribution. 0: Dirac, 1: lognormal distribution (mass), 11: lognormal distribution (number), 2: power-law distribution, 
							// 3: critical collapse distribution, 4: peak-theory distribution, 5:uniform, -1: User's distribution
	int spectrum_choice_param; // chooses the initial black holes spin distribution. 0: Dirac, 1: uniform, 2: gaussian.
	double amplitude; // amplitude of the distribution
	double stand_dev; // standard deviation of the mass distribution
	double crit_mass; // characteristic mass
	double eqstate; // equation of state parameter P = w*rho
	double stand_dev_param; // standard deviation of the parameter distribution
	double mean_param; // mean of the parameter distribution
	char table[32]; // name of the User's mass/spin distribution file
	
	int tmin_manual; // 0: automatic tmin for Dirac distribution, 1: manual tmin
	double tmin; // time of black hole formation
	int limit; // maximum iteration limit when computing the life evolution of black holes
	int nb_fin_times; // contains the number of integration times
	int BH_remnant; // 0: total evaporation after Planck mass, 1: stop the evaporation at M_relic
	double M_remnant; // mass of the BH remnant
	
	int E_number; // number of initial energies
	double Emin; // minimum initial energy
	double Emax; // maximum initial energy
	int particle_number; // number of particle species to be emitted
	int grav; // 0: no gravitons emitted, 1: gravitons emitted
	int add_DM; // 0: no DM, 1: one primary DM particle
	double m_DM; // DM mass in GeV
	double spin_DM; // spin of DM among 0., 1., 2., 0.5, 1.5
	double dof_DM; // number of DM degrees of freedom
	
	int primary_only; // if set to 1, the code will skip the hadronization part, if set to 0, the secondary spectra are computed
	int hadronization_choice; // choice of hadronization tables (0: PYTHIA BBN, 1: HERWIG BBN, 2: PYTHIA today, 3: HAZMA today, 4: HDMSpectra today)
		
	double Mmin_fM; // minimum BH mass in the f / g tables
	double Mmax_fM; // maximum BH mass in the f / g tables
	int nb_fM_masses; // number of black hole masses in the f / g tables
	int nb_fM_param; // number of black hole spins in the f(M,a) / g(M,a) tables, or BH espilon in the f(M,epsilon) tables
	double param_min; // minimum BH secondary parameter in the f,g and gamma tables
	double param_max; // maximum BH secondary parameter in the f,g and gamma tables
	int nb_gamma_param; // number of particle spins in the gamma(E,M,a) tables, or BH epsilon in the gamma(E,M,epsilon,a0) tables
	int nb_gamma_x; // number of the BH masses in the gamma tables
	int nb_gamma_spins; // number of spins (0, 1, 2, 1/2, [3/2]) tabulated
	int nb_gamma_fits; // number of fitting coefficients in the gamma fits tables
	double Emin_hadro; // minimum initial energy in the hadronization tables
	double Emax_hadro; // maximum initial energy in the hadronization tables
	int nb_init_en; // number of initial energies in the hadronization tables
	int nb_fin_en; // number of final energies in the hadronization tables
	int nb_init_part; // number of initial particles in the hadronization tables
	int nb_fin_part; // number of final particles in the hadronization tables
};

// int main();

// general functions in general.c
int read_params(struct param *parameters, char name[],int session);
int read_fM_infos(struct param *parameters);
int read_gamma_infos(struct param *parameters);
int read_hadronization_infos(struct param *parameters);
int memory_estimation(struct param *parameters,int session);

// functions computing the initial mass spectrum of PBH in spectrum.c
void read_users_table(double *init_masses,double *init_params,double **spec_table,struct param *parameters);
void spectrum(double *init_masses,double *init_params,double **spec_table,struct param *parameters);
double nu(double M);
double M_dist(double M,struct param *parameters);
double param_dist(double a,struct param *parameters);
void write_spectrum(double *init_masses,double *init_params,double **spec_table,struct param *parameters);

// functions computing the time evolution of PBH in evolution.c
double rplus_Kerr(double M,double a);
double temp_Kerr(double M,double a);
double P_LQG(double epsilon);
double m_LQG(double M,double epsilon);
double temp_LQG(double M,double epsilon,double a0);
double rplus_charged(double M,double Q);
double rminus_charged(double M,double Q);
double temp_charged(double M, double Q);
double rH_higher(double M,double n,double M_star);
double temp_higher(double M,double n,double M_star);
void read_fM_table(double **fM_table,double *fM_masses,double *fM_param,struct param *parameters);
void read_gM_table(double **gM_table,double *fM_masses,double *fM_param,struct param *parameters);
double loss_rate_M(double M,double param,double **fM_table,double *fM_masses,double *fM_param,int counter_M,int counter_param,struct param *parameters);
double loss_rate_param(double M,double param,double **fM_table,double **gM_table,double *fM_masses,double *fM_param,int counter_M,int counter_param,
	struct param *parameters);
void evolution_times(double *init_masses,double *init_params,double **evol_times,double **fM_table,double **gM_table,double *fM_masses,
	double *fM_param,struct param *parameters);
void sort_lifetimes(double **evol_times,double *sorted_times,int **rank,struct param *parameters);
void life_evolution(double ***life_masses,double ***life_params,double *life_times,double *dts,double *init_masses,
	double *init_params,int **rank,double **fM_table,double **gM_table,double *fM_masses,double *fM_param,struct param *parameters);
void write_life_evolutions(double ***life_masses,double ***life_params,double *life_times,double *dts,struct param *parameters);

// functions computing the primary emission spectra of PBH in primary.c
void read_gamma_tables(double ***gammas,double *gamma_param,double *gamma_x,struct param *parameters);
void read_asymp_fits(double ***fits,struct param *parameters);
double dNdtdE(double E,double M,double param,int particle_index,double ***gammas,double *gamma_param,double *gamma_x,double ***fits,double *dof,
	double *spins,double *masses_primary,int counter_param,int counter_x,struct param *parameters);
void instantaneous_primary_spectrum(double **instantaneous_primary_spectra,double **BH_masses,double **BH_params,double **spec_table,double *energies,
	double ***gammas,double *gamma_param,double *gamma_x,double ***fits,double *dof,double *spins,double *masses,int **counters_param,int ***counters_x,
	int *compute,struct param *parameters);
void write_instantaneous_primary_spectra(double **instantaneous_primary_spectra,double *energies,struct param *parameters);

// functions realizing the hadronization to get the secondary spectra of PBH in secondary.c
void read_hadronization_tables(double ****tables,double *initial_energies,double *final_energies,struct param *parameters);
void convert_hadronization_tables(double ****tables,double *initial_energies,double *final_energies,struct param *parameters);
void total_spectra(double ***partial_hadronized_spectra,double **partial_primary_spectra,double **partial_integrated_hadronized_spectra,
	double ****tables,double *initial_energies,double *final_energies,double **spec_table,double *times,double ***life_masses,double **BH_masses,
	double ***life_params,double **BH_params,double *energies,double *masses_primary,double *spins,double *dof,double *masses_secondary,
	double ***gammas,double *gamma_param,double *gamma_x,double ***fits,struct param *parameters);

void read_writing_instructions(int *write_primary,int *write_secondary,struct param *parameters);
void write_lines(char **file_names_primary,char **file_names_secondary,double **partial_primary_spectra,double **partial_integrated_hadronized_spectra,
	int *write_primary,int *write_secondary,double time,struct param *parameters);
void hadronize_instantaneous(double ***instantaneous_hadronized_spectra,double ****tables,double *initial_energies,double *final_energies,
	double **instantaneous_primary_spectra,double *energies,double *masses_secondary,int *compute,struct param *parameters);
double contribution_instantaneous(int j,int counter,int k,double **instantaneous_primary_spectra,double ****tables,double *initial_energies,
	double *final_energies,int particle_type,struct param *parameters);
void integrate_initial_energies_instantaneous(double ***instantaneous_hadronized_spectra,double **instantaneous_integrated_hadronized_spectra,
	double *energies,double *final_energies,int *compute,struct param *parameters);
void write_instantaneous_secondary_spectra(double **instantaneous_integrated_hadronized_spectra,double *hadronized_energies,struct param *parameters);
void write_lowE_photon_spectra(double **photons,double *final_energies,struct param *parameters);
void add_photons_instantaneous(double **instantaneous_primary_spectra,double **instantaneous_integrated_hadronized_spectra,double *energies,
	double *final_energies,struct param *parameters);
void add_neutrinos_instantaneous(double **instantaneous_primary_spectra,double **instantaneous_integrated_hadronized_spectra,double *energies,
	double *final_energies,struct param *parameters);
void add_electrons_instantaneous(double **instantaneous_primary_spectra,double **instantaneous_integrated_hadronized_spectra,double *energies,
	double *final_energies,struct param *parameters);
void add_FSR_instantaneous(double **instantaneous_primary_spectra,double **instantaneous_integrated_hadronized_spectra,
	double *energies,double *final_energies,double *masses_primary,struct param *parameters);

// function reading the PYTHIA hadronization table in hadro_pythia.c
void read_hadronization_pythia(double ****tables,double *initial_energies,double *final_energies,struct param *parameters);

// function reading the HERWIG hadronization table in hadro_herwig.c
void read_hadronization_herwig(double ****tables,double *initial_energies,double *final_energies,struct param *parameters);

// function reading the new PYTHIA hadronization table in hadro_pythianew.c
void read_hadronization_pythianew(double ****tables,double *initial_energies,double *final_energies,struct param *parameters);

// function reading the HAZMA hadronization table in hadro_hazma.c
void read_hadronization_hazma(double ****tables,double *initial_energies,double *final_energies,struct param *parameters);

// function reading the HDMSpectra hadronization table in hadro_hdmspectra.c
void read_hadronization_hdmspectra(double ****tables,double *initial_energies,double *final_energies,struct param *parameters);

// useful methods in technical.c
double exp_adapt(double x);
double trapeze(double x1,double x2,double y1,double y2);
void free1D_double(double *array);
void free1D_int(int *array);
void free2D_int(int **array,int l_1stD);
void free2D_char(char **array,int l_1stD);
void free2D_double(double **array,int l_1stD);
void free3D_int(int ***array,int l_1stD,int l_2ndD);
void free3D_double(double ***array,int l_1stD,int l_2ndD);
void free4D_double(double ****array,int l_1stD,int l_2ndD,int l_3rdD);
int ind_max(double *table,int llength);
void fusion(double *table,int start1,int end1,int end2);
void sort_fusion_bis(double *table,int start,int end);
void sort_fusion(double *table,int llength);
