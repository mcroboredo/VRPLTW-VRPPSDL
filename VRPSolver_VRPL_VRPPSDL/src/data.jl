import Unicode

mutable struct Vertex
   id::Int # 0 is the depot
   pos_x::Float64
   pos_y::Float64
   tw_begin::Float64
   tw_end::Float64
   dem_cap::Int # 0 for depot, demand for customers, or cap for lockers
end

# Undirected graph
mutable struct InputGraph
   V′::Array{Vertex} # set of vertices (access with id_vertex + 1)
   E::Array{Tuple{Int64,Int64}} # set of edges
   cost::Dict{Tuple{Int64,Int64},Float64} # cost for each edge
   E′::Array{Tuple{Int64,Int64}} # set of edges involving customers and lockers for assignment
end

mutable struct DataVRPL
   G′::InputGraph
   Q::Float64 # vehicle capacity
   nc::Int #number of customers
   nl::Int #number of lockers
   radius::Float64
   Name::String
end

# Euclidian distance
function distance(app,data::DataVRPL, arc::Tuple{Int64, Int64})
   e = (arc[1] < arc[2]) ? arc : (arc[2],arc[1])
   u, v = arc
   x_sq = (data.G′.V′[v+1].pos_x - data.G′.V′[u+1].pos_x)^2
   y_sq = (data.G′.V′[v+1].pos_y - data.G′.V′[u+1].pos_y)^2
   if app["vrppsdl"]
      return sqrt(x_sq + y_sq)
   end
   return round(sqrt(x_sq + y_sq), digits=2)
end

contains(p, s) = findnext(s, p, 1) != nothing

function readData(app::Dict{String,Any})

   str = Unicode.normalize(read(app["instance"], String); stripcc=true)
   breaks_in = [' '; ':'; '\n']
   aux = split(str, breaks_in; limit=0, keepempty=false)

   G′ = InputGraph([],[],Dict(),[])
   data = DataVRPL(G′,0.0,0,0,0.0,"")

   id = 0
   for s=1:length(app["instance"])
      if app["instance"][s] == '/'
         id = s
      end
   end
   data.Name = app["instance"][id+1:end-5]


   for i=1:length(aux)
      if contains(aux[i], "CUSTOMERS")
         data.nc = parse(Int, aux[i+1])
      elseif contains(aux[i], "LOCKERS")
         data.nl = parse(Int, aux[i+1])
      elseif contains(aux[i], "RADIUS")
         data.radius = parse(Float64, aux[i+1])
      elseif contains(aux[i], "CUSTOMER_SECTION")
         pos = i
         for k=1:data.nc
            push!(data.G′.V′, Vertex(parse(Int, aux[pos+1]), parse(Float64, aux[pos+2]), parse(Float64, aux[pos+3]), parse(Float64, aux[pos+4]), parse(Float64, aux[pos+5]), 1)) #Int(parse(Float64, aux[pos+6]))
            pos += 6
         end
      elseif contains(aux[i], "LOCKER_SECTION")
         pos = i
         for k=1:data.nl+1
            if k == 1
               pushfirst!(data.G′.V′,Vertex(parse(Int, aux[pos+1]), parse(Float64, aux[pos+2]), parse(Float64, aux[pos+3]), parse(Float64, aux[pos+4]), parse(Float64, aux[pos+5]), Int(parse(Float64, aux[pos+6]))))
            else
               push!(data.G′.V′,Vertex(data.nc+parse(Int, aux[pos+1]), parse(Float64, aux[pos+2]), parse(Float64, aux[pos+3]), parse(Float64, aux[pos+4]), parse(Float64, aux[pos+5]), Int(parse(Float64, aux[pos+6]))))
            end
            pos += 6
         end
      end
   end


   data.Q = ceil(data.nc/2)
   if app["vrppsdl"]
      data.Q = data.nc
   end

   for i=0:data.nc+data.nl-1
      for j=i+1:data.nc+data.nl
         push!(data.G′.E,(i,j))
         data.G′.cost[(i,j)] = distance(app,data,(i,j))
      end
   end

   for i=1:data.nc
      for j=data.nc+1:data.nc+data.nl
         if c(app,data,(i,j)) <= data.radius
            push!(data.G′.E′,(i,j))
         end
      end
   end

   return data
end

veh_capacity(data::DataVRPL) = data.Q
twb(data,i) = data.G′.V′[i+1].tw_begin
twe(data,i) = data.G′.V′[i+1].tw_end
lcap(data,i) = data.lcap[i]
d(data::DataVRPL, i) = data.G′.V′[i+1].dem_cap # return demand of i
edges(data::DataVRPL) = data.G′.E # return set of edges

function c(app,data,e) 
   if e[1] == e[2]
      return 0.0
   end
   if e[1] < e[2] 
      if data.Name[1]=='r'
         return 3*data.G′.cost[e] # cost of the edge e
      else
         return data.G′.cost[e] # cost of the edge e
      end
   else
      if data.Name[1]=='r'
         return 3*data.G′.cost[(e[2],e[1])] # cost of the edge e
      else
         return data.G′.cost[(e[2],e[1])]
      end
   end
end

function δ(data,j)
   incidents = []
   for e in data.G′.E
      if e[1] == j || e[2] == j
         push!(incidents,e)
      end
   end
   return incidents
end

function δ′(data,j)
   incidents = []
   for e in data.G′.E′
      if e[1] == j || e[2] == j
         push!(incidents,e)
      end
   end
   return incidents
end


