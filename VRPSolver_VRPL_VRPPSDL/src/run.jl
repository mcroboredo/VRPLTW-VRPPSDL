using BaPCodVRPSolver, JuMP, ArgParse, CPLEX, GLPKMathProgInterface, Unicode
include("data.jl")
include("model.jl")
include("solution.jl")

function parse_commandline(args_array::Array{String,1}, appfolder::String)
   s = ArgParseSettings(usage="##### VRPSolver #####\n\n"*
	   "  On interactive mode, call main([\"arg1\", ..., \"argn\"])", exit_after_help=false)
   @add_arg_table s begin
      "instance"
         help = "Instance file path"
      "--cfg", "-c"
         help = "Configuration file path"
         default = "$appfolder/../config/VRPL.cfg"
      "--ub","-u"
         help = "Upper bound (primal bound)"
         arg_type = Float64
         default = 10000000.0
      "--vrppsdl","-v"
         help = "Run the specialized VRPPSDL VRPSolver model. In this case, time windows constraints are turned on"
         action = :store_true
      "--tw","-w"
         help = "Time windows constraints"
         action = :store_true
      "--vrp_1", "-b"
         help = "VRPSolver model 1?"
         action = :store_true
      "--out","-o"
         help = "Path to write the solution found"
      "--tikz","-t"
         help = "Path to write the TikZ figure of the solution found."
   end
   return parse_args(args_array, s)
end

function run_cvrp(app::Dict{String,Any})
   println("Application parameters:")
   for (arg,val) in app
      println("  $arg  =>  $(repr(val))")
   end
   flush(stdout)

   instance_name = split(basename(app["instance"]), ".")[1]

   # ----- Get the instance ---------
   data = readData(app)
   #--------------------------

   if app["vrp_1"]
      (model, x,y,t) = build_model_vrp_1(data,app) #VRPSolver model 1
   else
      if app["vrppsdl"]
         (model, x,y,t) = build_model_vrp_2(data,app) #Specialized VRPPSDL VRPSolver model
      else
         (model, x,y,z,t) = build_model_vrp_2(data,app) #VRPSolver model 2
      end
   end
   
   optimizer = VrpOptimizer(model, app["cfg"], instance_name,baptreedot = "bap_tree.dot")
   set_cutoff!(optimizer, app["ub"]+0.1)
   
   (status, solution_found) = optimize!(optimizer)
   # if some solution is found we run a ILP model only to check the assignment
   if solution_found
      #set of lockers visited by some route
      Lockers = [l + data.nc for l=1:data.nl if get_value(optimizer, t[l]) > 0.9]  
      #customers assgined to some locker 
      Customers = []
      for j=1:data.nc
         for e in data.G′.E′
            if e[1] == j && get_value(optimizer, y[e]) > 0.01 && j ∉ Customers
               push!(Customers,j)
            end
         end
      end
      E⁺ = [(i,j) for i in Customers, j in Lockers if (i,j) ∈ data.G′.E′]
      

      m=Model(solver = GLPKSolverMIP())
      @variable(m, y[e in E⁺], Bin)
      @objective(m, Max, 0)
      for j in Customers
         @constraint(m,sum(y[e] for e in E⁺ if e[1] == j) == 1)
      end
      for l ∈ Lockers
         @constraint(m,sum(y[e] for e in E⁺ if e[2] == l) <= data.G′.V′[l+1].dem_cap)
      end
      status=solve(m)

      y_ = getvalue(y)
      E⁻ = [e for e in E⁺ if y_[e] > 0.9]

      #getting the solution (routes and asignments)
      sol = getsolution(data, optimizer, x,E⁻,t,get_objective_value(optimizer), app)
      #solution checker
      stat = checksol(app,data,sol)

      println("###############################################################################################################")
      println("custom_stats: ", instance_name, " ", app["ub"], " ", optimizer.stats[:bcRecRootDb], " ", optimizer.stats[:bcTimeRootEval]/100, " ", optimizer.stats[:bcCountNodeProc], " ",
      optimizer.stats[:bcRecBestDb], " ", optimizer.stats[:bcRecBestInc], " ", optimizer.stats[:bcTimeMain]/100)
      println("#############################################################################################################")

      println("*****************************************************************************")
      print_routes(data,sol)
      println("*****************************************************************************")

      if app["out"] != nothing
         #write the solution in a file
         writesolution(data,app["out"],sol)
      end

      if app["tikz"] != nothing
         #generate a tex file with the solution
         drawsolution(app["tikz"], data, sol)
      end
   else
      println("No solution was found for the instance ", instance_name)
   end   
end

function main(args)
   appfolder = dirname(@__FILE__)
   app = parse_commandline(args, appfolder)
   isnothing(app) && return
   run_cvrp(app)
end

if isempty(ARGS)
   main(["--help"])
else
   main(ARGS)
end
