
function build_model_vrp_1(data::DataVRPL, app::Dict{String,Any})

   E = edges(data) # set of edges of the input graph G′
   E′ = data.G′.E′
   nc = data.nc
   nl = data.nl
   V = [i for i in 0:nc+nl] # set of vertices of the graphs G′ and G
   V⁺ = [i for i in 1:nc] # set of customers of the input graph G′
   Q = veh_capacity(data)

   # Formulation
   cvrp = VrpModel()

   #variables
   @variable(cvrp.formulation, x[e in E], Int) #Number of times that the arc e is traversed
   @variable(cvrp.formulation, t[k=1:nl], Int) # Is the locker k used in the solution? 
   if app["vrppsdl"] == false
      @variable(cvrp.formulation, y[e in E′], Int) #Is the customer e[1] assigned  to the locker e[2]?
   else
      @variable(cvrp.formulation, y[e in E′], Int)
   end

   #Objective function (1a)
   if app["vrppsdl"] == false
      @objective(cvrp.formulation, Min, sum(c(app,data,e) * x[e] for e in E) + 0.5*sum(c(app,data,e) * y[e] for e in E′))
   else
      @objective(cvrp.formulation, Min, sum(0.5*x[(0,i)] for i=1:nc+nl) + sum(c(app,data,e) * x[e] for e in E) + sum(5 * y[e] for e in E′))
   end

   #Constraints
   @constraint(cvrp.formulation, deg[i in V⁺], sum(x[e] for e in δ(data, i)) + 2*sum(y[e] for e in δ′(data, i))== 2.0) #(1b)
   @constraint(cvrp.formulation, deg2[i=nc+1:nc+nl], sum(x[e] for e in δ(data, i))== 2.0*t[i-data.nc]) #(1c)
   @constraint(cvrp.formulation, y3[e in E′], y[e] <= t[e[2]-nc]) #(1d)
   @constraint(cvrp.formulation, deg5[i=nc+1:nc+nl], sum(y[e] for e in E′ if e[2] == i) <= data.G′.V′[i+1].dem_cap) #(1e)
   

   pack_arc = Dict()
   for j in setdiff(V,[0])
      pack_arc[j] = []
   end

   id_to_edge = Dict()
   edge_to_id = Dict()

   vert_locker = Dict()
   for l=nc+1:nc+nl
      vert_locker[l] = []
   end

   V′=[i for i in V]
   pos = V′[end] + 1
   for e in E′
      push!(V′,pos)
      edge_to_id[e] = pos
      id_to_edge[pos] = e
      #push!(pack_vert[e[1]],pos)
      push!(vert_locker[e[2]],pos)
      pos += 1
   end
   
   function build_graph()
      v_source = v_sink = 0

      #Graph
      G = VrpGraph(cvrp, V′, v_source, v_sink, (0, nc))

      #Resources
      cap_res_id = add_resource!(G, main = true) # R = R_M = {cap_res_id}
      if app["tw"] || app["vrppsdl"]
         time_res_id = add_resource!(G, main = true)
      end

      #Bounds
      for i in V′
         set_resource_bounds!(G, i, cap_res_id, 0.0, Float64(Q)) 
         if app["tw"] || app["vrppsdl"]
            if i <= nc + nl
               set_resource_bounds!(G, i, time_res_id, data.G′.V′[i+1].tw_begin, data.G′.V′[i+1].tw_end) 
            else
               set_resource_bounds!(G, i, time_res_id, data.G′.V′[1].tw_begin, data.G′.V′[1].tw_end) 
            end
         end
      end
      
      # Building set of arcs
      for i=0:nc+nl
         for j=0:nc+nl
            if i != j
               e = (i,j)
               if i > j
                  e = (j,i)
               end

               #criacao do arco
               if i > nc
                  arc_id = add_arc!(G, vert_locker[i][end], j)
               else
                  arc_id = add_arc!(G, i, j)
               end

               #packing sets
               if j <= nc + nl 
                  if j != 0
                     push!(pack_arc[j],arc_id)
                  end
               else
                  cust = id_to_edge[j][1] #cliente
                  push!(pack_arc[cust],arc_id)
               end
               #consumos
               if j <= nc
                  set_arc_consumption!(G, arc_id, cap_res_id, d(data,j))
               end
               if app["tw"] || app["vrppsdl"]
                  if app["tw"]
                     set_arc_consumption!(G, arc_id, time_res_id, c(app,data,e))
                  else
                     if j == 0
                        set_arc_consumption!(G, arc_id, time_res_id, c(app,data,e))
                     else
                        if j <= nc
                           set_arc_consumption!(G, arc_id, time_res_id, c(app,data,e)+ 5)
                        else
                           set_arc_consumption!(G, arc_id, time_res_id, c(app,data,e)+10)
                        end
                     end
                  end
               end

               #mapeamentos
               if j <= nc
                  add_arc_var_mapping!(G, arc_id, [x[e]])
               else
                  add_arc_var_mapping!(G, arc_id, [x[e],t[j-nc]])
               end
            end
         end
      end

      for l=nc+1:nc+nl
         j = id_to_edge[vert_locker[l][1]][1] #cliente
         #arco por cima
         arc_id = add_arc!(G,l, vert_locker[l][1])
         set_arc_consumption!(G, arc_id, cap_res_id,1)
         add_arc_var_mapping!(G, arc_id, [y[(j,l)]])
         #arco por baixo
         arc_id = add_arc!(G,l, vert_locker[l][1])

         for k=1:length(vert_locker[l])-1
            id₁ = vert_locker[l][k]
            id₂ = vert_locker[l][k+1]
            j = id_to_edge[vert_locker[l][k+1]][1] #cliente
            #arco por cima
            arc_id = add_arc!(G,id₁, id₂)
            set_arc_consumption!(G, arc_id, cap_res_id,1)
            add_arc_var_mapping!(G, arc_id, [y[(j,l)]])
            #arco por baixo
            arc_id = add_arc!(G,id₁, id₂)
         end
      end


      return G
   end

   Graphs = []
   G = build_graph()
   add_graph!(cvrp, G)
   push!(Graphs,G)
   #println(G)

   set_arc_packing_sets!(cvrp,[[(G,i) for i in pack_arc[k]] for k=1:nc+nl])
 
   set_branching_priority!(cvrp, "x", 1)
   set_branching_priority!(cvrp, "y", 1)
   set_branching_priority!(cvrp, "t", 10)

   return (cvrp, x,y,t)
end

function build_model_vrp_2(data::DataVRPL, app::Dict{String,Any})

   E = edges(data) # set of edges of the input graph G′
   E′ = data.G′.E′
   nc = data.nc
   nl = data.nl
   V = [i for i in 0:nc+nl] # set of vertices of the graphs G′ and G
   V⁺ = [i for i in 1:nc] # set of customers of the input graph G′
   Q = veh_capacity(data)

   demands = [data.G′.V′[k+1].dem_cap for k=1:nc]
   demands = reverse(sort(demands))
   max_dem = Dict()
   for l=nc+1:nc+nl
      demands = []
      max_dem[l] = 0
      for e in E′
         if e[2] == l
            push!(demands, data.G′.V′[e[1]+1].dem_cap) 
         end
      end
      demands = reverse(sort(demands))
      for k=1:5
         if k <= length(demands)
            max_dem[l] += demands[k]
         end
      end
   end

   # Formulation
   cvrp = VrpModel()

   #variables
   @variable(cvrp.formulation, x[e in E], Int) #Number of times that the arc e is traversed
   if app["vrppsdl"] == false
      @variable(cvrp.formulation, y[e in E′] >= 0,Int) #Is the customer e[1] assigned  to the locker e[2]?
   else
      @variable(cvrp.formulation, y[e in E′] >= 0)
   end
   @variable(cvrp.formulation, t[k=1:nl] <=1, Int) # Is the locker k used in the solution? 
   if app["vrppsdl"] == false
      @variable(cvrp.formulation, z[k=1:nl,c=1:max_dem[nc+k]], Int)
   end

   #Objective function
   if app["vrppsdl"] == false
      @objective(cvrp.formulation, Min, sum(c(app,data,e) * x[e] for e in E) + 0.5*sum(c(app,data,e) * y[e] for e in E′))
   else
      @objective(cvrp.formulation, Min, sum(0.5*x[(0,i)] for i=1:nc+nl) + sum(c(app,data,e) * x[e] for e in E) + sum(5 * y[e] for e in E′))
   end

   #Constraints
   @constraint(cvrp.formulation, deg[i in V⁺], sum(x[e] for e in δ(data, i)) + 2*sum(y[e] for e in δ′(data, i))== 2.0) #(1b)
   @constraint(cvrp.formulation, deg2[i=nc+1:nc+nl], sum(x[e] for e in δ(data, i))== 2.0*t[i-data.nc]) #(1c)
   @constraint(cvrp.formulation, y3[e in E′], y[e] <= t[e[2]-nc]) #(1d)
   @constraint(cvrp.formulation, deg5[i=nc+1:nc+nl], sum(y[e] for e in E′ if e[2] == i) <= data.G′.V′[i+1].dem_cap) #(1e)
   if app["vrppsdl"] == false
      @constraint(cvrp.formulation, deg3[i=nc+1:nc+nl], sum(data.G′.V′[e[1]+1].dem_cap*y[e] for e in E′ if e[2] == i)== sum(k*z[i-nc,k] for k=1:max_dem[i])) #(3c)
      @constraint(cvrp.formulation, deg4[i=nc+1:nc+nl], sum(z[i-nc,k] for k=1:max_dem[i]) <= 1) #(3d)
   end
   

   pack_arc = Dict()
   for j in setdiff(V,[0])
      pack_arc[j] = []
   end
   
   #routing path-generator graph
   function build_graph()
      v_source = v_sink = 0
      L = 0
      U = nc
  
      G = VrpGraph(cvrp, V, v_source, v_sink, (L, U))
      cap_res_id = add_resource!(G, main = true) # resource for the capacity constraints
      if app["tw"] || app["vrppsdl"]
         time_res_id = add_resource!(G, main = true) #resource for the time windows constraints
      end

      #Bounds for the resources
      for i in V
         set_resource_bounds!(G, i, cap_res_id, 0.0, Float64(Q)) 
         if app["tw"] || app["vrppsdl"]
            @show i, data.G′.V′[i+1].tw_begin, data.G′.V′[i+1].tw_end
            set_resource_bounds!(G, i, time_res_id, data.G′.V′[i+1].tw_begin, data.G′.V′[i+1].tw_end) 
         end
      end
      
      # Building set of arcs
      for i=0:nc+nl
         for j=0:nc+nl
            if i != j
               e = (i,j)
               if i > j
                  e = (j,i)
               end
               if j <= nc
                  arc_id = add_arc!(G, i, j)
                  if j != 0
                     set_arc_consumption!(G, arc_id, cap_res_id, d(data,j))
                  end
                  if app["tw"] || app["vrppsdl"]
                     if app["vrppsdl"]
                        if j == 0
                           set_arc_consumption!(G, arc_id, time_res_id, c(app,data,e))
                        else
                           if data.Name[1] == 'r'
                              set_arc_consumption!(G, arc_id, time_res_id, c(app,data,e)+5)
                           else
                              set_arc_consumption!(G, arc_id, time_res_id, c(app,data,e)+5)
                           end
                        end
                     else
                        if app["tw"]
                           set_arc_consumption!(G, arc_id, time_res_id, c(app,data,e))
                        end
                     end
                  end
                  add_arc_var_mapping!(G, arc_id, [x[e]])
               else
                  if app["vrppsdl"] == false
                     for k=1:max_dem[j]
                        arc_id = add_arc!(G, i, j)
                        set_arc_consumption!(G, arc_id, cap_res_id, k)
                        if app["tw"] 
                           set_arc_consumption!(G, arc_id, time_res_id, c(app,data,e))
                        end
                        add_arc_var_mapping!(G, arc_id, [x[e],t[j-nc],z[j-nc,k]])
                     end
                  else
                     arc_id = add_arc!(G, i, j)
                     if data.Name[1] == 'r'
                        set_arc_consumption!(G, arc_id, time_res_id, c(app,data,e)+10)
                     else
                        set_arc_consumption!(G, arc_id, time_res_id, c(app,data,e)+10)
                     end
                     add_arc_var_mapping!(G, arc_id, [x[e],t[j-nc]])
                  end
               end
            end
         end
      end

      return G
   end

   #assignment path-generator graph for each locker l
   function build_graph2(l)
      v_source = 0
      v_sink = 2*nc+1
      Vert = [0]
      for j =1:nc
         if (j,l) in E′
            push!(Vert,j)
         end
      end
      for j =1:nc
         push!(Vert,j+nc)
      end
      push!(Vert,v_sink)
      Gra = VrpGraph(cvrp, Vert, v_source, v_sink, (0, 1))

      # cap_res_id = add_resource!(Gra, main = true) # R = R_M = {cap_res_id}
      # for i in Vert
      #    set_resource_bounds!(Gra, i, cap_res_id, 0.0, Float64(data.G′.V′[l+1].dem_cap)) 
      # end
      
      for j=1:nc+1
         if j == 1
            if j != nc+1 && (j,l) in E′
               if j-1 in Vert && j in Vert
                  arc_id = add_arc!(Gra, j-1, j)
                  add_arc_var_mapping!(Gra, arc_id, [y[(j,l)]])
                  #set_arc_consumption!(Gra, arc_id, cap_res_id, 1)
               end
            end
            if j-1 in Vert
               arc_id = add_arc!(Gra, j-1, nc+j)
            end
         else
            if j != nc +1
               # j is being visited
               if j != nc+1 && (j,l) in E′
                  if j-1 in Vert && j in Vert
                     arc_id = add_arc!(Gra, j-1, j)
                     add_arc_var_mapping!(Gra, arc_id, [y[(j,l)]])
                     #set_arc_consumption!(Gra, arc_id, cap_res_id, 1)
                  end
                  if j in Vert
                     arc_id = add_arc!(Gra, nc+j-1, j)
                     add_arc_var_mapping!(Gra, arc_id, [y[(j,l)]])
                     #set_arc_consumption!(Gra, arc_id, cap_res_id, 1)
                  end
               end
            end

            # j is not being visited
            if j-1 in Vert
               arc_id = add_arc!(Gra, j-1, nc+j)
            end
            arc_id = add_arc!(Gra, nc+j-1, nc+j)

         end

      end


      return Gra
      
   end

   existe = Dict()

   Graphs = []
   G = build_graph()
   add_graph!(cvrp, G)
   push!(Graphs,G)
   for j=1:nc
      existe[(G,j)] = true
   end

   if app["vrppsdl"] == false
      for l=nc+1:nc+nl
         Gra = build_graph2(l)
         add_graph!(cvrp, Gra)
         push!(Graphs,Gra)

         for j=1:nc
            if (j,l) in E′
               existe[(Gra,j)] = true
            else
               existe[(Gra,j)] = false
            end
         end
      end
   end
   
   set_vertex_packing_sets!(cvrp, union([[(g,j) for g in Graphs if existe[(g,j)]] for j=1:nc], [[(G,j)] for j=nc+1:nc+nl]) )
   define_elementarity_sets_distance_matrix!(cvrp, G, [[c(app,data, (i, j)) for j in setdiff(V ,[0])] for i in setdiff(V ,[0])])
 
 
   set_branching_priority!(cvrp, "x", 1)
   set_branching_priority!(cvrp, "t", 10)
   if app["vrppsdl"] == false
      set_branching_priority!(cvrp, "y", 1)
      set_branching_priority!(cvrp, "z", 1)
   end
   
   if !app["vrppsdl"]
      return (cvrp, x,y,z,t)
   else
      return (cvrp, x,y,t)
   end
end
