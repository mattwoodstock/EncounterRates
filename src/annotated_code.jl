using LinearAlgebra, Statistics, DataFrames, CSV, Tables, Dates
using PlanktonIndividuals, Plots
import JLD2

cd("D:/pz_encounter")

# Function to calculate the optimal movement vector for prey to avoid multiple predators in 3D
function optimal_movement(prey_position, predator_positions, avoidance_distance)
    # Calculate direction away from all predators

    direction_vector = transpose(prey_position) .- predator_positions
    distances = sum(direction_vector.^2, dims=2)
    closest_predator = argmin(distances)
  
    direction_vector = direction_vector[closest_predator, :]
    direction_vector /= sqrt(sum(direction_vector.^2))  # Normalize
    
    # Calculate the move fraction based on predator locations

    # Example: Move fraction decreases as the average distance to predators increases
    average_distance = mean(sqrt.(sum((transpose(prey_position) .- predator_positions).^2, dims=2)))

    move_frac = min(1,1.0 - 0.1 * average_distance)
  
    # Calculate distance to move (avoidance distance)
    move_distance = move_frac * min(avoidance_distance, sqrt(sum(direction_vector.^2)))
  
    # Update prey's position
    optimal_movement = prey_position .+ direction_vector * move_distance
  
    return optimal_movement
end


function plot_model(model::PlanktonModel, z_pause::Vector{Float64}, ΔT::Float64, deactivate::Bool, encounter_diags::Dict{String,Any},densities::Matrix{Int},isubexp::Int64,iteration::Int64; z_pause_ΔT=0.0::Float64, velo1=1e-4::Float64, velo2 = 1e-5::Float64, encounter_tresh=0.001::Float64,prey_vis=1::Float64, pred_vis = 1::Float64)

    # Blank count for prey numbers
    prey_num = []

    # square here as sqrt is not taken in computation of distance
    encounter_tresh_sq = encounter_tresh^2

    # ac stands for active
    ip_ac = model.individuals.phytos.sp1.data.ac .== 1.0
    p_pos = transpose(hcat(Array(model.individuals.phytos.sp1.data.x[ip_ac]), Array(model.individuals.phytos.sp1.data.y[ip_ac]), Array(model.individuals.phytos.sp1.data.z[ip_ac])))
    
    iz_ac = model.individuals.phytos.sp2.data.ac .== 1.0
    z_pos = transpose(hcat(Array(model.individuals.phytos.sp2.data.x[iz_ac]), Array(model.individuals.phytos.sp2.data.y[iz_ac]), Array(model.individuals.phytos.sp2.data.z[iz_ac])))
    
    # convert logical to integer array
    ip_ac_int = findall(ip_ac)
    iz_ac_int = findall(iz_ac)
    
    
    # compute distance
    dist_3d = zeros(Float64,length(ip_ac_int),length(iz_ac_int))

    for i in 1:size(dist_3d,1)
        for j in 1:size(dist_3d,2)
                dist_3d[i,j] = sqrt((p_pos[1,i]-z_pos[1,j])^2+(p_pos[2,i]-z_pos[2,j])^2+(p_pos[3,i]-z_pos[3,j])^2)
        end
    end

    old_dist =dist_3d[1,1]


    ## Prey movement (Does not calculate encounters.)
    if velo2 > 0.0
        vec = Array{Float64}(undef, 3)
        max_movement_len = ΔT * velo2
        # square here as sqrt was not taken in computation of distance
        max_movement_len_sq = max_movement_len

        # find closest and build vector
        for (ip, iz_model) in enumerate(ip_ac_int)
            iz = argmin(dist_3d[ip, :]) #Find closest predator

            #@assert z_pos[1, iz] == model.individuals.phytos.sp2.data.x[iz_model]
            #@assert z_pos[2, iz] == model.individuals.phytos.sp2.data.y[iz_model]
            #@assert z_pos[3, iz] == model.individuals.phytos.sp2.data.z[iz_model]

            prey_visual_range = ((3 / (4 * π)) * prey_vis)^(1/3) #Calculates the radial distance the prey can visualize based on the proportion provided.
            
            #Placeholder
            # check if predators are within range. Need to change so predator finds "optimal location".

            if dist_3d[ip, iz] > prey_visual_range #Animal moves to random vector because closest animal is outside of visual range
                for idim in 1:3
                    vec[idim] = 2* rand() -1 #Animal moves to random vector between -1 and 1

                    #if dist_periodic[idim, ip, iz] #Remove boundary effects?
                    #    println(dist_periodic[idim,ip,iz])
                    #    throw(ErrorException("stop"))
                    #    if vec[idim] > 0.0
                    #        vec[idim] = vec[idim] - 1.0
                    #    else
                    #        vec[idim] = vec[idim] + 1.0
                    #    end
                    #end
                end

                #vec /= sqrt(sum(vec.^2))

                p_pos[:,ip] = p_pos[:,ip] .+ vec .* max_movement_len

            else #Animal makes decision to move to "optimal" location.
                #Prey position
                prey_position = [model.individuals.phytos.sp1.data.x[ip],model.individuals.phytos.sp1.data.y[ip],model.individuals.phytos.sp1.data.z[ip]]


                #Find all preds within visual range
                ind_dists = dist_3d[ip,:]
                sub_dist = ind_dists[ind_dists .< prey_visual_range]

                pred_positions = Matrix{Float64}(undef,size(sub_dist,1),3)
                sort_vector =sort(dist_3d[ip,:])

                if size(sub_dist,1) > 0
                    for i in 1:size(sub_dist,1)

                        pred_dist = sort_vector[i] #Next closest prey
                        index = findall(x -> x == pred_dist,ind_dists)

                        #Gather x,y,z coordinates of next closest predator
                        x = model.individuals.phytos.sp2.data.x[index][1]
                        y = model.individuals.phytos.sp2.data.y[index][1]
                        z = model.individuals.phytos.sp2.data.z[index][1]

                        pred_positions[i,:] .= (x,y,z)
                    end

                    max_distance = max_movement_len_sq #Can complicate this to correspond to a flexible time frame 

                    optimal_location = optimal_movement(prey_position,pred_positions,max_distance)

                    p_pos[:, ip] = optimal_location
                end
            end

            #Update location
            model.individuals.phytos.sp1.data.x[iz_model] = p_pos[1, ip]
            model.individuals.phytos.sp1.data.y[iz_model] = p_pos[2, ip]
            model.individuals.phytos.sp1.data.z[iz_model] = p_pos[3, ip]

        end
    end

    ## Predator movement (Calculates encounters; i.e., predators catch up to prey)
    
    # ac stands for active
    ip_ac = model.individuals.phytos.sp1.data.ac .== 1.0
    p_pos = transpose(hcat(Array(model.individuals.phytos.sp1.data.x[ip_ac]), Array(model.individuals.phytos.sp1.data.y[ip_ac]), Array(model.individuals.phytos.sp1.data.z[ip_ac])))

    iz_ac = model.individuals.phytos.sp2.data.ac .== 1.0
    z_pos = transpose(hcat(Array(model.individuals.phytos.sp2.data.x[iz_ac]), Array(model.individuals.phytos.sp2.data.y[iz_ac]), Array(model.individuals.phytos.sp2.data.z[iz_ac])))


    # convert logical to integer array
    ip_ac_int = findall(ip_ac)
    iz_ac_int = findall(iz_ac)

    # compute distance
    for i in 1:size(dist_3d,1)
        for j in 1:size(dist_3d,2)
            dist_3d[i,j] = sqrt((p_pos[1,i]-z_pos[1,j])^2+(p_pos[2,i]-z_pos[2,j])^2+(p_pos[3,i]-z_pos[3,j])^2)
        end
    end

    if velo1 > 0.0
        vec = Array{Float64}(undef, 3)
        max_movement_len = ΔT * velo1
        # square here as sqrt was not taken in computation of distance
        max_movement_len_sq = max_movement_len

        # find closest and build vector
        for (iz, iz_model) in enumerate(iz_ac_int)
            if z_pause[iz_model] > 0.0
                z_pause[iz_model] = max(0.0, z_pause[iz_model] - ΔT)
                continue
            end
            ip = argmin(dist_3d[:, iz])

            #@assert z_pos[1, iz] == model.individuals.phytos.sp2.data.x[iz_model]
            #@assert z_pos[2, iz] == model.individuals.phytos.sp2.data.y[iz_model]
            #@assert z_pos[3, iz] == model.individuals.phytos.sp2.data.z[iz_model]

            pred_visual_range = ((3 / (4 * π)) * pred_vis)^(1/3) #Calculates the radial distance the prey can visualize based on the proportion provided.

            if iteration == 300  #Want the final result of the prey densities to show the effect
                densities[iz,isubexp] = length(dist_3d[:,iz][dist_3d[:,iz] .< pred_visual_range])
            end
            count = 1
            # check if closest phytoplankton can be reached in this timestep
            if dist_3d[ip, iz] ≤ max_movement_len_sq
                #println("moving zooplankton $(iz_model) directly to phytoplankton $(ip_ac_int[ip])")
                z_pos[:, iz] = p_pos[:, ip]

            elseif dist_3d[ip,iz] ≤ pred_visual_range #Predator cannot catch preys, but it sees them. Moves towards closest prey
                vec .= p_pos[:, ip] .- z_pos[:, iz]
                    #if dist_periodic[idim, ip, iz]
                    #    if vec[idim] > 0.0
                    #        vec[idim] = vec[idim] - 1.0
                    #    else
                    #        vec[idim] = vec[idim] + 1.0
                    #    end
                    #end
                #vec /= sqrt(sum(vec.^2))
                normalized_movement_vector = vec / dist_3d[ip,iz]
                
                z_pos[:, iz] = z_pos[:,iz] .+ normalized_movement_vector .* max_movement_len

            else #Predator moves in random vector
                for idim in 1:3
                    vec[idim] = 2* rand() -1 #Animal moves to random vector between -1 and 1

                    #if dist_periodic[idim, ip, iz] #Remove boundary effects?
                    #    if vec[idim] > 0.0
                    #        vec[idim] = vec[idim] - 1.0
                    #    else
                    #        vec[idim] = vec[idim] + 1.0
                    #    end
                    #end
                end
                #vec /= sqrt(sum(vec.^2))

                z_pos[:, iz] = z_pos[:,iz] .+ vec .* max_movement_len
            end

            model.individuals.phytos.sp2.data.x[iz_model] = z_pos[1, iz]
            model.individuals.phytos.sp2.data.y[iz_model] = z_pos[2, iz]
            model.individuals.phytos.sp2.data.z[iz_model] = z_pos[3, iz]
        end
    else
        ind_z = z_pause .> 0
        z_pause[ind_z] = max.(0.0, z_pause[ind_z] .- ΔT)
    end

    #
    # encounters and create plot
    #
    for i in 1:size(dist_3d,1)
        for j in 1:size(dist_3d,2)
            dist_3d[i,j] = sqrt((p_pos[1,i]-z_pos[1,j])^2+(p_pos[2,i]-z_pos[2,j])^2+(p_pos[3,i]-z_pos[3,j])^2)
        end
    end

    plt = Plots.scatter(p_pos[1, :], p_pos[2, :], p_pos[3,:], ms=4, color="#228833", markerstrokewidth=0, size=(900, 900), legend=:none)
    Plots.scatter!(z_pos[1, :], z_pos[2, :], z_pos[3,:], ms=5, color="#CCBB44", markerstrokewidth=0, legend=:none)

    encounters = findall(dist_3d .< encounter_tresh_sq)

    iz_vec = [ci.I[2] for ci in encounters]
    Plots.scatter!(z_pos[1, iz_vec], z_pos[2, iz_vec], z_pos[3,iz_vec], ms=5, color = "red", legend=:none)

    num_encounters = 0
    for ci in encounters
        # keep paused individuals unable to encounter
        if z_pause[iz_ac_int[ci.I[2]]] > 0.0
            continue
        end
        # ci is of type CartesianIndex
        Plots.scatter!([p_pos[1, ci.I[1]], z_pos[1, ci.I[2]]], [p_pos[2, ci.I[1]], z_pos[2, ci.I[2]]],[p_pos[3, ci.I[1]], z_pos[3, ci.I[2]]], ms=5, color="#EE6677", markerstrokewidth=0, legend=:none)
        Plots.scatter!([z_pos[1, ci.I[2]]], [z_pos[2, ci.I[2]]], [z_pos[3, ci.I[2]]], ms=5, color="#EE6677", markerstrokewidth=0, legend=:none)
        if deactivate
            #println("deactivating phytoplankton $(ip_ac_int[ci.I[1]])")
            model.individuals.phytos.sp1.data.ac[ip_ac_int[ci.I[1]]] = 0.0
        else
            model.individuals.phytos.sp1.data.x[ip_ac_int[ci.I[1]]] = rand()
            model.individuals.phytos.sp1.data.y[ip_ac_int[ci.I[1]]] = rand()
            model.individuals.phytos.sp1.data.z[ip_ac_int[ci.I[1]]] = rand()
        end
        z_pause[iz_ac_int[ci.I[2]]] = z_pause_ΔT
        num_encounters += 1
    end
    push!(encounter_diags["num_encounters"], num_encounters)

    Plots.xlims!(0, 1)
    Plots.ylims!(0, 1)
    Plots.zlims!(0, 1)
	return densities
end

function run_test(N_replicates=10, N_individual=[2^8, 2^3]::Vector{Int64},densities = zeros()::Matrix{Int}; z_pause_ΔT=0.0::Float64, velo1=1e-4::Float64,velo2=1e-5::Float64, encounter_tresh=0.001::Float64, deactivate=false::Bool, arch=CPU(),prey_vis=1.0::Float64,pred_vis = 1.0::Float64,isubexp = 1::Int64,N_iter = 300::Int64)

    ## Can play around with grid. Scale grid size to larger distance to match predator velocity?
    grid = RectilinearGrid(size=(1,1,1),
                           x = (0, 10meters),
                           y = (0, 10meters),
                           z = (0,-10meters),
                           topology = (Periodic, Periodic, Periodic))


    bgc_params = bgc_params_default()
    bgc_params["κhP"] = 1e-6
    bgc_params["κvP"] = 1e-6

    N_species = 2
    mode = CarbonMode()

    deactivate_relocate = "deactivate"
    if ! deactivate
        deactivate_relocate = "relocate"
    end

    # set division and mortality to 0
    phyt_params = phyt_params_default(N_species, mode)
    phyt_params["dvid_P"] = [0.0, 0.0]
    phyt_params["mort_P"] = [0.0, 0.0]

    ΔT = 60.0

    replicate_diags = Dict{String, Any}()
    replicate_diags["sum_encounters"] = Vector{Int64}()
    replicate_diags["mean_encounters"] = Vector{Float64}()

    z_pause = zeros(Float64, N_individual[2])

    preys = DataFrame(TS = [], Max = [], Min = [], Mean = [], SD = [])

    for i = 1:N_replicates
        model = PlanktonModel(arch, grid; N_species = N_species,
                                          N_individual = N_individual,
                                          max_individuals = sum(N_individual),
                                          mode = mode,
                                          bgc_params = bgc_params,
                                          phyt_params = phyt_params)

        #model.individuals.phytos.sp1.data.ac[N_individual[1]+1:end] .= 0.0
        #model.individuals.phytos.sp2.data.ac[N_individual[1]+1:end] .= 0.0

        z_pause[1:end] .= 0.0

        pow = PlanktonOutputWriter(save_diags=true,
                                   save_plankton=true)

        sim = PlanktonSimulation(model, ΔT=ΔT, iterations = 1, output_writer=pow)
        encounter_diags = Dict{String, Any}()
        encounter_diags["num_encounters"] = Vector{Int64}()

        #anim = @animate for i in 1:N_iter
        for i in 1:N_iter
            update!(sim)

            densities = plot_model(model, z_pause, ΔT, deactivate, encounter_diags,densities,isubexp,i; z_pause_ΔT=z_pause_ΔT, encounter_tresh=encounter_tresh, velo1=velo1, velo2 = velo2,prey_vis=prey_vis,pred_vis = pred_vis)

            #plt = plot_model(model, z_pause, ΔT, deactivate, encounter_diags; z_pause_ΔT=z_pause_ΔT, encounter_tresh=encounter_tresh, velo1=velo1, velo2 = velo2,prey_vis=prey_vis,pred_vis = pred_vis)
            #if i < 10
            #    filename = "frames/frame_movement_$(deactivate_relocate)_3d_00$(i).png"
            #    Plots.savefig(plt, filename)
            #end
            #plt
            
        end
        #gif(anim, joinpath("frames", "exp $(isubexp).gif"), fps = 5) # used to be 15 fps

        #println(encounter_diags)
        push!(replicate_diags["sum_encounters"], sum(encounter_diags["num_encounters"]))
        push!(replicate_diags["mean_encounters"], mean(encounter_diags["num_encounters"]))
    end
    println("Prey swimming speed: $(velo2), n_prey = $(N_individual[1]), prey vis: $(prey_vis), ")
    println("number of encounters: $(replicate_diags["sum_encounters"]), mean number of encounters: $(replicate_diags["mean_encounters"])")

    return replicate_diags["sum_encounters"], N_iter, densities
end

function run_experiments()
    res2 = DataFrame(Pred_velo = [], Prey_velo = [], PredN = [], PreyN = [], Speed = [], Prey_vis = [], Pred_vis = [], Encounter = [],Enc_sd = [])
    new_preys = []

    isubexp = 1

    ΔT = 60.0 # NOTE: not passed on, just used in computation for pause

    n_rep = 10 #Number of replicates to run
    encounter_tresh = 0.01  # increased encounter_tresh
    n_z = 100 #Number of predators
    # exponents = 6:12
    n_pz = [25:25:75;100:100:900;]

    n_it = 1
    #velos1 = (0, 3e-6, 1e-5, 3e-5, 1e-4) #Will want this to be data-derived and scaled to model area
    velos1 = 1e-4 #Predator velocity
    pred_visual = 0.5
    velos2 = (1e-6,1e-5,3e-5,5e-5,7e-5,1e-4,3e-4,1e-3) #Order of magnitude less swimming speed. Will want this to be data-derived

    prey_visual = [0:0.0001:0.0001;0.001;0.025:0.025:0.1;0.2:0.2:1;] #Prey visual ranges to test

    z_pause_ΔT = 8.0*ΔT
    res = Dict{String, Dict}()

    expname = "exp_$(isubexp)"
    N_iter = 300
    densities = zeros(Int,n_z,length(n_pz)*length(velos2)*length(prey_visual))

    for velo1 in velos1
        for pred_vis in pred_visual
                for velo2 in velos2
                    for prey_vis in prey_visual
                    #=
                            res[expname] = Dict(
                                "description" => "varying num(Prey) (Pred speed: $(velo1)) (Prey speed: $(velo2))",
                                "num_encounters" => fill(-1, (length(n_pz), n_rep)),
                                "num_prey" => Vector{Int64}(),
                                "num_pred" => n_z,
                                "pred_speed" => velo1,
                                "prey_speed" => velo2,
                                "encounter_tresh" => encounter_tresh,
                                "encounter_type" => "relocate",
                                "z_pause_ΔT" => z_pause_ΔT,
                            )
                            =#

                            for (i, n) in enumerate(n_pz)
                                start = now()

                                n_p = n
                                #push!(res[expname]["num_prey"], n_p)                        
                                tmp, n_it, densities = run_test(n_rep, [n_p,n_z],densities; velo1=velo1, velo2 = velo2, z_pause_ΔT=z_pause_ΔT, encounter_tresh=encounter_tresh,prey_vis = prey_vis,pred_vis = pred_vis,isubexp,N_iter)
                                #res[expname]["num_encounters"][i, :] = tmp
                                #res[expname]["num_iter"] = n_it

                                new_row = Dict("Pred_velo" => velo1, "Prey_velo" => velo2, "PredN" => n_z, "PreyN" => n_pz[i], "Speed" => velo1/velo2, "Prey_vis" => prey_vis, "Pred_vis" => pred_vis, "Encounter" => mean(tmp),"Enc_sd" => std(tmp))
                                push!(res2,new_row)
                                CSV.write("Prey visual range results.csv",res2)
                                CSV.write("Prey Densities.csv",Tables.table(densities), writeheader=false)
                        isubexp += 1
                        stop = now()
                        println(stop-start)
                    end
                end
            end
        end
    end
end
run_experiments()



