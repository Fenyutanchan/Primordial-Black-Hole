using Artifacts

mutable struct BlackHawkParameters
    destination_folder::String              # name of the output folder in results/
    full_output::Int                        # quantity of information displayed (0=less, 1=more)
    interpolation_method::Int               # interpolation in the numerical tables (0=linear, 1=logarithmic)

    metric::Int                             # BH metric: 0=Kerr, 1=polymerized, 2=charged, 3=higher-dimensional

    BH_number::Int                          # number of BH masses (should be the number of tabulated masses if spectrum_choice=5)
    Mmin::Float64                           # lowest BH mass in g (larger than the Planck mass)
    Mmax::Float64                           # highest BH mass in g (larger than the Planck mass)
    param_number::Int                       # number of Kerr spins
    amin::Float64                           # lowest Kerr spin
    amax::Float64                           # highest Kerr spin
    Qmin::Float64                           # lowest Reissner-Norström charge
    Qmax::Float64                           # highest Reissner-Norström charge

    epsilon_LQG::Float64                    # dimensionless epsilon parameter for the polymerized metric
    a0_LQG::Float64                         # minimal area for the polymerized metric in GeV^(-2)
    n::Int                                  # number of extra spatial dimensions in higher-dimensional metric

    spectrum_choice::Int                    # form of the BH distribution: 0=Dirac, 1=log-normal for the mass, 11: log-normal for the number, 2=power-law, 3=critical collapse, 4=peak theory, 5=uniform -1=user-defined
    spectrum_choice_param::Int              # form of the spin distribution for each mass: 0=Dirac, 1=uniform, 2=Gaussian

    amplitude_lognormal::Float64            # amplitude of the log-normal (mass density) distribution in g.cm^-3
    amplitude_lognormal2::Float64           # amplitude of the log-normal (number density) distribution in cm^-3
    stand_dev_lognormal::Float64            # dimensionless variance of the log-normal distribution
    crit_mass_lognormal::Float64            # characteristic mass of the log-normal distribution in g

    amplitude_powerlaw::Float64             # amplitude of the power-law distribution in g^(gamma-1).cm^-3
    eqstate_powerlaw::Float64               # equation of state of the Universe at the BH formation time P = w.rho

    amplitude_critical_collapse::Float64    # amplitude of the critical collapse distribution in g^(-2.85).cm^-3
    crit_mass_critical_collapse::Float64    # characteristic mass of the critical collapse distribution in g

    amplitude_uniform::Float64              # amplitude of the uniform mass distribution in cm^(-3)

    stand_dev_param_Gaussian::Float64       # standard deviation of the Gaussian spin distribution
    mean_param_Gaussian::Float64            # mean of the Gaussian spin distribution

    table::String                           # table containing the User's BH distribution

    tmin_manual::Int                        # 1: user-defined tmin, 0:automatically set tmin
    tmin::Float64                           # initial integration time of the evolution of BH in s
    limit::Int                              # iteration limit when computing the time evolution of a single BH
    BH_remnant::Int                         # 0: total evaporation, 1: BH relic at mass M_relic
    M_remnant::Float64                      # BH relic mass in g

    E_number::Int                           # number of primary particles energies to be simulated
    Emin::Float64                           # minimal energy in GeV of the primary particles
    Emax::Float64                           # maximal energy in GeV of the primary particles

    grav::Int                               # 0=no graviton, 1=emission of gravitons
    add_DM::Int                             # 0=no DM added, 1=one DM particle
    m_DM::Float64                           # DM mass in GeV
    spin_DM::Float64                        # DM spin
    dof_DM::Float64                         # number of DM degrees of freedom

    primary_only::Int                       # 1=no secondary spectrum, 0=secondary spectrum computed

    hadronization_choice::Int               # 0=PYTHIA at the BBN epoch, 1=HERWIG at the BBN epoch, 2=PYTHIA (new) at the present epoch, 3=HAZMA at the present epoch, 4=HDMSpectra at the present epoch

    function BlackHawkParameters(
        ; verbose::Bool=false
    )                                       # default parameter setting
        verbose && @info "Initializing BlackHawkParameters with default values!"

        return new(
            "test",                         # name of the output folder in results/
            1,                              # quantity of information displayed (0=less, 1=more)
            0,                              # interpolation in the numerical tables (0=linear, 1=logarithmic)
            0,                              # BH metric: 0=Kerr, 1=polymerized, 2=charged, 3=higher-dimensional
            1,                              # number of BH masses (should be the number of tabulated masses if spectrum_choice=5)
            1e9,                            # lowest BH mass in g (larger than the Planck mass)
            1e16,                           # highest BH mass in g (larger than the Planck mass)
            1,                              # number of Kerr spins
            0.,                             # lowest Kerr spin
            .5,                             # highest Kerr spin
            0.,                             # lowest Reissner-Norström charge
            .7,                             # highest Reissner-Norström charge
            1.5,                            # dimensionless epsilon parameter for the polymerized metric
            0.,                             # minimal area for the polymerized metric in GeV^(-2)
            0,                              # number of extra spatial dimensions in higher-dimensional metric
            0,                              # form of the BH distribution: 0=Dirac, 1=log-normal for the mass, 11: log-normal for the number, 2=power-law, 3=critical collapse, 4=peak theory, 5=uniform -1=user-defined
            0,                              # form of the spin distribution for each mass: 0=Dirac, 1=uniform, 2=Gaussian
            1.,                             # amplitude of the log-normal (mass density) distribution in g.cm^-3
            1.,                             # amplitude of the log-normal (number density) distribution in cm^-3
            1.,                             # dimensionless variance of the log-normal distribution
            1.,                             # characteristic mass of the log-normal distribution in g
            1.,                             # amplitude of the power-law distribution in g^(gamma-1).cm^-3
            1/3,                            # equation of state of the Universe at the BH formation time P = w.rho
            1.,                             # amplitude of the critical collapse distribution in g^(-2.85).cm^-3
            1.,                             # characteristic mass of the critical collapse distribution in g
            1.,                             # amplitude of the uniform mass distribution in cm^(-3)
            1.,                             # standard deviation of the Gaussian spin distribution
            1.,                             # mean of the Gaussian spin distribution
            "spin_distribution_BH.txt",     # table containing the User's BH distribution
            0,                              # 1: user-defined tmin, 0:automatically set tmin
            1e-30,                          # initial integration time of the evolution of BH in s
            5000,                           # iteration limit when computing the time evolution of a single BH
            0,                              # 0: total evaporation, 1: BH relic at mass M_relic
            1e-4,                           # BH relic mass in g
            10,                             # number of primary particles energies to be simulated
            5,                              # minimal energy in GeV of the primary particles
            1e5,                            # maximal energy in GeV of the primary particles
            1,                              # 0=no graviton, 1=emission of gravitons
            1,                              # 0=no DM added, 1=one DM particle
            1.,                             # DM mass in GeV
            0.,                             # DM spin
            1.,                             # number of DM degrees of freedom
            0,                              # 1=no secondary spectrum, 0=secondary spectrum computed
            2                               # 0=PYTHIA at the BBN epoch, 1=HERWIG at the BBN epoch, 2=PYTHIA (new) at the present epoch, 3=HAZMA at the present epoch, 4=HDMSpectra at the present epoch
        )
    end
end

write_parameters(params::BlackHawkParameters) =
    open(joinpath(BlackHawk_directory, "$(params.destination_folder).in"), "w+") do io
        write(io, """
        destination_folder = $(params.destination_folder)
        full_output = $(params.full_output)
        interpolation_method = $(params.interpolation_method)
        metric = $(params.metric)
        BH_number = $(params.BH_number)
        Mmin = $(params.Mmin)
        Mmax = $(params.Mmax)
        param_number = $(params.param_number)
        amin = $(params.amin)
        amax = $(params.amax)
        Qmin = $(params.Qmin)
        Qmax = $(params.Qmax)
        epsilon_LQG = $(params.epsilon_LQG)
        a0_LQG = $(params.a0_LQG)
        n = $(params.n)
        spectrum_choice = $(params.spectrum_choice)
        spectrum_choice_param = $(params.spectrum_choice_param)
        amplitude_lognormal = $(params.amplitude_lognormal)
        amplitude_lognormal2 = $(params.amplitude_lognormal2)
        stand_dev_lognormal = $(params.stand_dev_lognormal)
        crit_mass_lognormal = $(params.crit_mass_lognormal)
        amplitude_powerlaw = $(params.amplitude_powerlaw)
        eqstate_powerlaw = $(params.eqstate_powerlaw)
        amplitude_critical_collapse = $(params.amplitude_critical_collapse)
        crit_mass_critical_collapse = $(params.crit_mass_critical_collapse)
        amplitude_uniform = $(params.amplitude_uniform)
        stand_dev_param_Gaussian = $(params.stand_dev_param_Gaussian)
        mean_param_Gaussian = $(params.mean_param_Gaussian)
        table = $(params.table)
        tmin_manual = $(params.tmin_manual)
        tmin = $(params.tmin)
        limit = $(params.limit)
        BH_remnant = $(params.BH_remnant)
        M_remnant = $(params.M_remnant)
        E_number = $(params.E_number)
        Emin = $(params.Emin)
        Emax = $(params.Emax)
        grav = $(params.grav)
        add_DM = $(params.add_DM)
        m_DM = $(params.m_DM)
        spin_DM = $(params.spin_DM)
        dof_DM = $(params.dof_DM)
        primary_only = $(params.primary_only)
        hadronization_choice = $(params.hadronization_choice)
        """)
    end

function set_BlackHawk_name!(params::BlackHawkParameters, name::String)
    params.destination_folder = name
    return params
end

function set_BlackHawk_mass!(params::BlackHawkParameters, mass::EnergyUnit; NU::Union{Missing, NaturalUnit}=missing)
    if ismissing(NU)
        NU = (NaturalUnit ∘ typeof)(mass)
    end

    @check_EU_dimension mass 1
    @check_positive_value mass

    # @assert mass > NU.m_Pl "Mass must be larger than the Planck mass."

    params.BH_number = 1
    params.Mmin = mass / NU.g
    params.Mmax = 10 * params.Mmin

    return params
end

function set_BlackHawk_time!(params::BlackHawkParameters,
    t::EnergyUnit, iteration_limit::Int;
    NU::Union{Missing, NaturalUnit}=missing
)
    if ismissing(NU)
        NU = (NaturalUnit ∘ typeof)(t)
    end

    @check_EU_dimension t -1
    @check_positive_value t

    params.tmin_manual = 1
    params.tmin = t / NU.s
    params.limit = iteration_limit

    return params
end

function set_BlackHawk_energy_spectrum!(params::BlackHawkParameters,
    Eₘᵢₙ::EnergyUnit, Eₘₐₓ::EnergyUnit, number_of_energies::Int
)
    @check_EU_dimension Eₘᵢₙ 1
    @check_EU_dimension Eₘₐₓ 1
    @check_positive_value Eₘᵢₙ
    @check_positive_value Eₘₐₓ

    @assert Eₘᵢₙ < Eₘₐₓ "Eₘᵢₙ must be smaller than Eₘₐₓ."

    params.E_number = number_of_energies
    params.Emin = EUval(GeV, Eₘᵢₙ)
    params.Emax = EUval(GeV, Eₘₐₓ)

    return params
end

function close_BlackHawk_dark_matter!(params::BlackHawkParameters)
    params.add_DM = 0
    return params
end

function set_BlackHawk_BBN!(params)
    params.hadronization_choice = 0
    return params
end

function run_BlackHawk_tot(params::BlackHawkParameters; verbose::Bool=false)
    write_parameters(params)
    working_directory = pwd()
    cd(BlackHawk_directory)
    verbose ? run(`./BlackHawk_tot.x $(params.destination_folder).in`) :
        (run ∘ pipeline)(`./BlackHawk_tot.x $(params.destination_folder).in`, devnull)
    cd(working_directory)
    rm(joinpath(BlackHawk_directory, "$(params.destination_folder).in"); force=true, recursive=true)
    verbose && @info "`BlackHawk_tot.x`: Completed successfully."

    return nothing
end

function run_BlackHawk_inst(params::BlackHawkParameters; verbose::Bool=false)
    write_parameters(params)
    working_directory = pwd()
    cd(BlackHawk_directory)
    verbose ? run(`./BlackHawk_inst.x $(params.destination_folder).in`) :
        (run ∘ pipeline)(`./BlackHawk_inst.x $(params.destination_folder).in`, devnull)
    cd(working_directory)
    rm(joinpath(BlackHawk_directory, "$(params.destination_folder).in"); force=true, recursive=true)
    verbose && @info "`BlackHawk_inst.x`: Completed successfully."

    return nothing
end

function read_BlackHawk_instantaneous_spectra(BlackHawk_result_directory::String)
    @assert isdir(BlackHawk_result_directory) "The directory $(BlackHawk_result_directory) does not exist."

    primary_spectra_file = joinpath(BlackHawk_result_directory, "instantaneous_primary_spectra.txt")
    secondary_spectra_file = joinpath(BlackHawk_result_directory, "instantaneous_secondary_spectra.txt")

    @assert isfile(primary_spectra_file) "The file $(primary_spectra_file) does not exist."
    @assert isfile(secondary_spectra_file) "The file $(secondary_spectra_file) does not exist."

    primary_spectra_data, primary_spectra_header = readdlm(primary_spectra_file; header=true, skipstart=1)
    @assert first(primary_spectra_header) == "energy/particle"
    primary_spectra = Dict{String, Vector{Float64}}()
    primary_spectra["energy"] = primary_spectra_data[:, 1]
    for ii ∈ 2:size(primary_spectra_data, 2)
        primary_spectra[primary_spectra_header[ii]] = primary_spectra_data[:, ii]
    end

    secondary_spectra_data, secondary_spectra_header = readdlm(secondary_spectra_file; header=true, skipstart=1)
    @assert first(secondary_spectra_header) == "energy/particle"
    secondary_spectra = Dict{String, Vector{Float64}}()
    secondary_spectra["energy"] = secondary_spectra_data[:, 1]
    for ii ∈ 2:size(secondary_spectra_data, 2)
        secondary_spectra[secondary_spectra_header[ii]] = secondary_spectra_data[:, ii]
    end

    return primary_spectra, secondary_spectra
end

function prepare_BlackHawk(; overwrite_flag::Bool=true)
    overwrite_flag && rm(BlackHawk_directory; force=true, recursive=true)
    if isdir(BlackHawk_directory)
        @warn "`BlackHawk` directory already exists. Skipping initialization..."
        return nothing
    end
    rm(BlackHawk_directory; force=true, recursive=true)
    cp(joinpath(artifact"BlackHawk", "blackhawk_v2.3"), BlackHawk_directory)

    for (root, dirs, files) ∈ walkdir(BlackHawk_hack_directory)
        target_dir = BlackHawk_directory * replace(root, BlackHawk_hack_directory => "")
        for file ∈ files
            source_file = joinpath(root, file)
            target_file = joinpath(target_dir, file)
            cp(source_file, target_file; force=true)
        end
    end

    for (root, dirs, files) ∈ (walkdir ∘ joinpath)(BlackHawk_directory, "src", "tables")
        chmod.(joinpath.(root, files), 0o644)
    end

    working_directory = pwd()
    cd(BlackHawk_directory)
    run(`make BlackHawk_inst.c BlackHawk_tot.c`)
    cd(working_directory)
    @info "BlackHawk directory initialized successfully."

    return nothing
end
