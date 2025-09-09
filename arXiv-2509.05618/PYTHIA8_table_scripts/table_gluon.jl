using Distributed, JLD2

@info """
Added $(nworkers()) workers for $(basename(@__FILE__)).
Please run the following command to start the workers if the number of workers is not your desired one:
`julia --procs=<number of workers> $(basename(@__FILE__))`
"""

@everywhere import PYTHIA8

@everywhere include("tool_script-count_energy.jl")
@everywhere include("tool_script-geomspace.jl")
include("tool_script-directories.jl")

number_of_initial_energies = 100
Emin_initial = 1e3
Emax_initial = 1e5 # Above this energy, Pythia8 works too slowly.
@everywhere number_of_final_energies = 10000
@everywhere Emin_final = 1e-1
@everywhere Emax_final = 1e5
@everywhere number_of_events = 10^5

E_initial_list = geomspace(Emin_initial, Emax_initial, number_of_initial_energies)
@everywhere E_final_list = geomspace(Emin_final, Emax_final, number_of_final_energies)

@everywhere pythia = begin
    pythia = PYTHIA8.Pythia("", false)
    PYTHIA8.readString(pythia, "Beams:idA = 11")
    PYTHIA8.readString(pythia, "Beams:idB = -11")
    PYTHIA8.readString(pythia, "PDF:lepton = off")
    PYTHIA8.readString(pythia, "PhaseSpace:pTHatMin = 0.")
    PYTHIA8.readString(pythia, "HiggsSM:ffbar2H = on")
    PYTHIA8.readString(pythia, "PartonLevel:ISR = off")
    PYTHIA8.readString(pythia, "25:onMode = off")
    PYTHIA8.readString(pythia, "25:onIfAny = 21")
    PYTHIA8.readString(pythia, "Print:quiet = on")

    pythia
end

@info "Pythia 8 initialized."

@everywhere log_dir = joinpath((first ∘ splitext)(@__FILE__) * "_logs")
rm(log_dir; recursive=true, force=true)
mkdir(log_dir)
@info """
Please run `watch 'ls $log_dir/*.tmp | wc -l'` to monitor the number of files.
"""

@everywhere function magic(E_initial::Real)
    photon = zeros(Int, number_of_final_energies)
    electron = zeros(Int, number_of_final_energies)
    muon = zeros(Int, number_of_final_energies)
    nu_e = zeros(Int, number_of_final_energies)
    nu_mu = zeros(Int, number_of_final_energies)
    nu_tau = zeros(Int, number_of_final_energies)
    pipm = zeros(Int, number_of_final_energies)
    K0L = zeros(Int, number_of_final_energies)
    Kpm = zeros(Int, number_of_final_energies)
    proton = zeros(Int, number_of_final_energies)
    neutron = zeros(Int, number_of_final_energies)

    # actual_number_of_events = number_of_events

    PYTHIA8.readString(pythia, "Beams:eCM = $(2 * E_initial)")
    PYTHIA8.init(pythia)
    event = PYTHIA8.event(pythia)

    ii_event = 1

    while ii_event ≤ number_of_events
        while !PYTHIA8.next(pythia)
            PYTHIA8.next(pythia)
        end

        event_continue_flag = false
        event_size = PYTHIA8.size(event)
        for ii_particle = 1:event_size
            particle = event[ii_particle]

            if PYTHIA8.status(particle) == -23
                PDG_ID = PYTHIA8.id(particle)
                if PDG_ID ≠ 21
                    event_continue_flag = true
                    break
                end
                continue
            end

            PYTHIA8.isFinal(particle) || continue

            PDG_ID = PYTHIA8.id(particle)
            energy = PYTHIA8.e(particle)

            if PDG_ID == 22
                count_energy!(energy, E_final_list, photon)
            elseif PDG_ID == 11 || PDG_ID == -11
                count_energy!(energy, E_final_list, electron)
            elseif PDG_ID == 13 || PDG_ID == -13
                count_energy!(energy, E_final_list, muon)
            elseif PDG_ID == 12 || PDG_ID == -12
                count_energy!(energy, E_final_list, nu_e)
            elseif PDG_ID == 14 || PDG_ID == -14
                count_energy!(energy, E_final_list, nu_mu)
            elseif PDG_ID == 16 || PDG_ID == -16
                count_energy!(energy, E_final_list, nu_tau)
            elseif PDG_ID == 211 || PDG_ID == -211
                count_energy!(energy, E_final_list, pipm)
            elseif PDG_ID == 130
                count_energy!(energy, E_final_list, K0L)
            elseif PDG_ID == 321 || PDG_ID == -321
                count_energy!(energy, E_final_list, Kpm)
            elseif PDG_ID == 2212 || PDG_ID == -2212
                count_energy!(energy, E_final_list, proton)
            elseif PDG_ID == 2112 || PDG_ID == -2112
                count_energy!(energy, E_final_list, neutron)
            end
        end

        event_continue_flag && continue
        ii_event += 1
    end

    write(joinpath(log_dir, "$(E_initial).tmp"), "")

    return Dict(
        "photon" => photon,
        "electron" => electron,
        "muon" => muon,
        "nu_e" => nu_e,
        "nu_mu" => nu_mu,
        "nu_tau" => nu_tau,
        "pipm" => pipm,
        "K0L" => K0L,
        "Kpm" => Kpm,
        "proton" => proton,
        "neutron" => neutron
    )
end

results = pmap(magic, E_initial_list)
# rm(log_dir; recursive=true, force=true)

data_path = joinpath(output_data_directory, (first ∘ splitext ∘ basename)(@__FILE__) * ".jld2")
jldopen(data_path, "w"; compress=true) do jld
    @info "Saving results to `$data_path`..."
    jld["number of events"] = number_of_events
    jld["E [GeV] (initial)"] = E_initial_list
    jld["E [GeV] (final)"] = E_final_list
    
    photon_data = zeros(Int, number_of_initial_energies, number_of_final_energies)
    electron_data = zeros(Int, number_of_initial_energies, number_of_final_energies)
    muon_data = zeros(Int, number_of_initial_energies, number_of_final_energies)
    nu_e_data = zeros(Int, number_of_initial_energies, number_of_final_energies)
    nu_mu_data = zeros(Int, number_of_initial_energies, number_of_final_energies)
    nu_tau_data = zeros(Int, number_of_initial_energies, number_of_final_energies)
    pipm_data = zeros(Int, number_of_initial_energies, number_of_final_energies)
    K0L_data = zeros(Int, number_of_initial_energies, number_of_final_energies)
    Kpm_data = zeros(Int, number_of_initial_energies, number_of_final_energies)
    proton_data = zeros(Int, number_of_initial_energies, number_of_final_energies)
    neutron_data = zeros(Int, number_of_initial_energies, number_of_final_energies)
    for (ii, E_initial) ∈ enumerate(E_initial_list)
        result = results[ii]
        photon_data[ii, :] = result["photon"]
        electron_data[ii, :] = result["electron"]
        muon_data[ii, :] = result["muon"]
        nu_e_data[ii, :] = result["nu_e"]
        nu_mu_data[ii, :] = result["nu_mu"]
        nu_tau_data[ii, :] = result["nu_tau"]
        pipm_data[ii, :] = result["pipm"]
        K0L_data[ii, :] = result["K0L"]
        Kpm_data[ii, :] = result["Kpm"]
        proton_data[ii, :] = result["proton"]
        neutron_data[ii, :] = result["neutron"]
    end
    jld["photon"] = photon_data
    jld["electron"] = electron_data
    jld["muon"] = muon_data
    jld["nu_e"] = nu_e_data
    jld["nu_mu"] = nu_mu_data
    jld["nu_tau"] = nu_tau_data
    jld["pipm"] = pipm_data
    jld["K0L"] = K0L_data
    jld["Kpm"] = Kpm_data
    jld["proton"] = proton_data
    jld["neutron"] = neutron_data

    @info "Results saved successfully."
end
