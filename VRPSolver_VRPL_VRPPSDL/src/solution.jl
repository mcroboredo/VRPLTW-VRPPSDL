mutable struct Solution
   cost::Float64
   routes::Array{Array{Int}}
   assignments::Dict{Int,Array{Int}}
end

function tw_route(app,data,r)
   time = 0.0
   for k=1:length(r)-1
      time += c(app,data,(r[k],r[k+1]))
      if time < data.G′.V′[r[k+1]+1].tw_begin
         time = data.G′.V′[r[k+1]+1].tw_begin
      end
      if time > data.G′.V′[r[k+1]+1].tw_end
         return false
      end
   end
   return true
end

function checksol(app::Dict{String,Any},data::DataVRPL,sol::Solution)

   E = edges(data) # set of edges of the input graph G′
   E′ = data.G′.E′
   nc = data.nc
   nl = data.nl
   L = [i for i in nc+1:nc+nl] # set of vertices of the graphs G′ and G
   V⁺ = [i for i in 1:nc] # set of customers of the input graph G′
   Q = veh_capacity(data)

   travel_cost = 0.0
   assingnment_cost = 0.0
   objetctive_function = 0.0

   for l in L
      for j in sol.assignments[l]
         if app["vrppsdl"] == false
            assingnment_cost += 0.5*c(app,data,(j,l))
         else
            assingnment_cost += 5
         end
      end
   end
   
   visited_c = Dict()
   visited_l = Dict()

   for j in V⁺
      visited_c[j] = 0
   end

   for l in L
      visited_l[l] = 0
   end

   ca = 0
   for l in L
      if length(sol.assignments[l]) > data.G′.V′[l+1].dem_cap
         println("error : locker ", l, " capacity violated ", data.G′.V′[l+1].dem_cap, " ", length(sol.assignments[l]))
         return -1
      end
      for j in sol.assignments[l]
         ca += 5
         if (j,l) ∉ E′
            println("error: the customer ", j, " can not be assigned to locker ", l, " but it is" )
            return -1
         else
            visited_c[j] += 1
         end
      end
   end

   ct = 0
   for r in sol.routes
      time_r = 0.0
      for k=1:length(r) - 1
         travel_cost += c(app,data,(r[k],r[k+1]))
         t = c(app,data,(r[k],r[k+1]))
         ct += t
         if k != length(r) - 1
            if r[k+1] <= nc
               visited_c[r[k+1]] += 1
               if app["vrppsdl"]
                  t += 5
               end
            else
               visited_l[r[k+1]] += 1
               if app["vrppsdl"]
                  t += 10
               end
            end
         end
         if time_r + t < data.G′.V′[r[k+1]+1].tw_begin
            time_r = data.G′.V′[r[k+1]+1].tw_begin
         else
            time_r += t
         end
         if time_r > data.G′.V′[r[k+1]+1].tw_end && (app["vrppsdl"] || app["tw"])
            println("error: vertex ", r[k+1], " tw is violated ", data.G′.V′[r[k+1]+1].tw_end, " ", time_r)
            return -1
         end
      end
   end

   for j in V⁺
      if visited_c[j] != 1
         println("error: customer ", j, " is being visited or assigned ", visited_c[j], " times ")
      end
   end

   for l in L
      if visited_l[l] > 1
         println("error: locker ", l, " is being visited ", visited_l[l], " times ")
      end
   end

   if app["vrppsdl"]
      objetctive_function = travel_cost + assingnment_cost + length(sol.routes)
   else
      objetctive_function = travel_cost + assingnment_cost
   end

   if objetctive_function > sol.cost + 0.001 || objetctive_function < sol.cost - 0.001
      println("error in the objective function ", objetctive_function, " ", sol.cost)
   end

   return 1

end

function getsolution(data::DataVRPL, optimizer::VrpOptimizer, x,E⁻,t,objval, app::Dict{String,Any})
   nc = data.nc
   nl = data.nl

   routes = []
   assignments = Dict()
   E = edges(data) # set of edges of the input graph G′
   E′=data.G′.E′

   travel_cost = 0.0
   assingnment_cost = 0.0

   for e in E
      if get_value(optimizer, x[e]) > 0.9
         travel_cost += get_value(optimizer, x[e])*c(app,data,e)
      end
   end

   for e in E⁻
      if app["vrppsdl"] == false
         assingnment_cost += 0.5*c(app,data,e)
      else
         assingnment_cost += 5
      end
   end

   if app["vrppsdl"] == false
      #@show travel_cost + assingnment_cost
   else
      #@show travel_cost + assingnment_cost + get_number_of_positive_paths(optimizer)
   end

   for path_id in 1:get_number_of_positive_paths(optimizer)
      route_edeges, route, visited, visits =[], [], Dict(), Dict()

      for e in E
         visits[e] = 0
         if get_value(optimizer, x[e], path_id) > 0.9
            push!(route_edeges,e)
            visited[e] = floor(get_value(optimizer, x[e], path_id)+0.5)
         end
      end

      current = 0
      for e in route_edeges
         if e[1] == 0
            visits[e] = visits[e] + 1
            route = [0, e[2]]
            current = e[2]
            break
         end
      end

      while current != 0
         for e in route_edeges
            if visits[e] < visited[e] && (e[1] == current || e[2] == current)
               visits[e] = visits[e] + 1
               if e[1] == current
                  current = e[2]
               else 
                  current = e[1]
               end
               push!(route, current)
            end
         end
      end
      if !tw_route(app,data,route)
         route = reverse(route)
      end
      #@show route
      if length(route_edeges) > 0 
         push!(routes, route)
      end
   end

   for l=data.nc+1:data.nc+data.nl
      assignments[l] = []
   end

   for e in E⁻
      push!(assignments[e[2]],e[1])
   end
   
   return Solution(objval, routes,assignments)
end

function print_routes(data, solution)
   for (i,r) in enumerate(solution.routes)
      print("Route #$i: ") 
      for j in r
         print("$j ")
      end
      println()
   end



   for j=data.nc+1:data.nc+data.nl
      if length(solution.assignments[j]) > 0
         print("Assignments for the locker $j: " )
         for k in solution.assignments[j]
            print("$k ")
         end
         println()
      end
   end
   

   ct = solution.cost
   println("Cost: $(ct)")
end

# write solution in a file
function writesolution(data,solpath, solution)
   open(solpath, "w") do f
      for (i,r) in enumerate(solution.routes)
         write(f, "Route #$i: ")
         for j in r
            write(f, "$j ") 
         end
         write(f, "\n")
      end

      for j=data.nc+1:data.nc+data.nl
         if length(solution.assignments[j]) > 0
            write(f,"Assignments for the locker $j: " )
            for k in solution.assignments[j]
               write(f,"$k ")
            end
            write(f, "\n")
         end
      end

      write(f, "Cost $(solution.cost)\n")
   end
end

# write solution as TikZ figure (.tex) 
function drawsolution(tikzpath, data, solution)
   open(tikzpath, "w") do f
      write(f,"\\documentclass[crop,tikz]{standalone}\n\\begin{document}\n")
      # get limits to draw
      pos_x_vals = [i.pos_x for i in data.G′.V′]
      pos_y_vals = [i.pos_y for i in data.G′.V′]
      scale_fac = 1/(max(maximum(pos_x_vals),maximum(pos_y_vals))/10)
      write(f,"\\begin{tikzpicture}[thick, scale=1, every node/.style={scale=0.3}]\n")
      for i in data.G′.V′
         x_plot = scale_fac*i.pos_x
         y_plot = scale_fac*i.pos_y
         if i.id == 0 # plot depot
            write(f, "\t\\node[draw, line width=0.1mm, rectangle, fill=white, inner sep=0.05cm, scale=1.4] (v$(i.id)) at ($(x_plot),$(y_plot)) {\\footnotesize $(i.id)};\n")
            # Uncomment to plot without vertex id
            #write(f, "\t\\node[draw, rectangle, fill=yellow, scale=1.4] (v$(i.id_vertex)) at ($(x_plot),$(y_plot)) {};\n")
         elseif i.id <= data.nc #plot customer
            write(f, "\t\\node[draw, line width=0.1mm, circle, fill=white, inner sep=0.05cm] (v$(i.id)) at ($(x_plot),$(y_plot)) {\\footnotesize $(i.id)};\n")
            # Uncomment to plot without vertex id
            #write(f, "\t\\node[draw, circle, fill=white] (v$(i.id_vertex)) at ($(x_plot),$(y_plot)) {};\n")
         else #locker
            write(f, "\t\\node[draw, line width=0.1mm, rectangle, fill=yellow, inner sep=0.05cm] (v$(i.id)) at ($(x_plot),$(y_plot)) {\\footnotesize $(i.id)};\n")
            if data.Name[1] == 'r'
               write(f,"\t\\draw[red,dashed,-,line width=0.2pt,opacity=.4] (v$(i.id)) circle ($(scale_fac*data.radius/3));\n")
            else
               write(f,"\t\\draw[red,dashed,-,line width=0.2pt,opacity=.4] (v$(i.id)) circle ($(scale_fac*data.radius));\n")
            end
         end
      end
      
      for r in solution.routes
         for k=1:length(r)-1
            edge_style =  "dashed,-,line width=0.4pt"
            write(f, "\t\\draw[$(edge_style)] (v$(r[k])) -- (v$(r[k+1]));\n")
         end
      end

      for l=data.nc+1:data.nc+data.nl
         for j in solution.assignments[l]
            edge_style =  "dashed,-,line width=0.4pt"
            write(f, "\t\\draw[blue, $(edge_style)] (v$(l)) -- (v$(j));\n")
         end
      end

      
      write(f, "\\end{tikzpicture}\n")
      write(f, "\\end{document}\n")
   end   
end