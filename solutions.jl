### A Pluto.jl notebook ###
# v0.17.2

using Markdown
using InteractiveUtils

# ╔═╡ 821803bc-93b1-436d-9adb-2a759a782034
md"""
# AoC 2021
"""

# ╔═╡ d152b170-8593-4128-a950-ec14a66af2c2
md"""
## Day 1
"""

# ╔═╡ b745b668-52a3-11ec-2566-c1d79a56cfc2
let
	# first part
	var, res_1 = open("./D1-1.txt") do f
		var = parse.(Int,readlines(f))
		var, sum(var[2:end] - var[1:end-1] .> 0)
	end
	
	# second part
	var_2 = map(i-> sum(@view var[i-1:i+1]), 2:length(var)-1)
	res_2 = sum(var_2[2:end] - var_2[1:end-1] .> 0)

	(res_1, res_2)
end

# ╔═╡ c08f1079-ca65-458a-837d-ffbda5046bf3
md"""
## Day 2
"""

# ╔═╡ 1f3a8a64-e78b-407e-8835-f9a0fa419c46
let
	# first part
	input = readlines("./D2-1.txt")
	hor, depth = 0, 0
	for instruction in input
		command, qty = split(instruction," ")
		if isequal(command, "forward")
			hor += parse(Int, qty)
		elseif isequal(command, "down")
			depth += parse(Int, qty)
		elseif isequal(command, "up")
			depth -= parse(Int, qty)
		end
	end
	res_1 = hor * depth
	
	# second part
	hor, depth, aim = 0, 0, 0
	for instruction in input
		command, qty = split(instruction," ")
		if isequal(command, "down")
			aim += parse(Int, qty)
		elseif isequal(command, "up")
			aim -= parse(Int, qty)
		elseif isequal(command, "forward")
			hor += parse(Int, qty)
			depth += aim * parse(Int, qty)
		end
	end
	res_2 = hor * depth 

	(res_1, res_2)
end

# ╔═╡ f6e71aca-953a-4320-b3ab-be377e0aa22c
md"""
## Day 3
"""

# ╔═╡ 942b6adf-426a-432e-9c04-116ece48bd9d
let
	# first part
	input = readlines("./D3-1.txt")
	m = hcat([parse.(Bool,split(line,"")) for line in input]...)
	γ_bin = sum(m, dims=2) .> length(input) / 2
	γ = parse(Int,prod(["$(Int(b))" for b in γ_bin]), base=2)
	ϵ_bin = .! γ_bin 
	ϵ = parse(Int,prod(["$(Int(b))" for b in ϵ_bin]), base=2)
	
	res_1 = γ * ϵ
	# second part
	function parser!(list; majority::Bool=true)
		# start at first bit
		bitpos = 1
		
		# loop
		while bitpos <= size(list, 1)
			# get majority/minority value
			if majority
				ctlval = count(list[bitpos,:]) >= size(list,2)/2 ? 1 : 0
			else
				ctlval = count(list[bitpos,:]) >= size(list,2)/2 ? 0 : 1
			end
			
			# build filter & filter
			list_filter = findall(x->x==ctlval, list[bitpos,:])
			list = list[:, list_filter]
			
			# stop if required
			if size(list,2) == 1
				break
			end

			bitpos += 1
		end
		return bitstring(list)
	end

	O2 = parse(Int,filter(x -> !isspace(x), parser!(deepcopy(m), majority=true)),base=2)
	CO2 = parse(Int,filter(x -> !isspace(x), parser!(deepcopy(m), majority=false)),base=2)
	res_2 = O2 * CO2

	(res_1, res_2)
end

# ╔═╡ 134adeb8-1759-4fec-8a62-3906bb72c451
md"""
## Day 4
"""

# ╔═╡ 5aa396ce-fed7-45fa-a2a5-6fe40a90bf84
let
	function inputparser(path::String="./D4-2.txt")
		input = readlines(path)
		draw = parse.(Int,split(input[1],","))
		N_fields = round(Int, (length(input) - 1)/6)
		fields = Array{Int}(undef, 5,5, N_fields)
		drawn = zeros(Bool, size(fields))
		fieldnbr = 1 # counter for field
		line = 3 # startline
		while fieldnbr <= N_fields
			for l = 1:5
				res = parse.(Int,filter(x->!isequal(x,""),split(input[line + l - 1]," ")))
				fields[l, :, fieldnbr] .= res
			end
			line += 6
			fieldnbr += 1
		end
		
		return draw, fields, drawn
	end

	function checkbingo(drawn, grid_size::Int=5)
		if grid_size ∈ sum(drawn, dims= 1)
			return true
		elseif  grid_size ∈ sum(drawn, dims= 2)
			return true
		else
			return false
		end
	end

	function playbingo(path::String="./D4-2.txt")
		draw, fields, drawn = inputparser(path)
		winningboards = Any[]
		scores = Any[]
		for number in draw
			# add number
			drawn[findall(x->isequal(x, number), fields)] .= true
			# checkbingo
			bingo = checkbingo(drawn)
			if bingo
				# check simultaneous winners
				winners = unique(map(x->x[3], 
					vcat(findall(x->isequal(x, 5),sum(drawn,dims=1)), 
						 findall(x->isequal(x, 5),sum(drawn,dims=2)))))
				# add winners and their score
				push!(winningboards, winners)
				score = map(card -> sum(fields[:,:,card][findall(.! drawn[:,:,card])]) * number, winners)
				push!(scores, score)
				# remove cards from deck after bingo
				fields = fields[:,:,setdiff(1:size(fields,3), winners)]
				drawn = drawn[:,:,setdiff(1:size(drawn,3), winners)]
			end
			
		end
		return (scores[1][1], scores[end][1])
	end

	playbingo("./D4-1.txt")
end

# ╔═╡ b8993212-0328-416b-8c27-6b6ae3aca8fd
md"""
## Day 5
"""

# ╔═╡ ad35bd5a-401a-44b5-9786-1b1744fabe34
let
	input = readlines("./D5-1.txt")
	# part 1
	obs_horver = Dict()
	obs_horverdiag = Dict()
	for line in input
		p1, p2 = split(line," -> ")
		x1,y1 = parse.(Int,split(p1,","))
		x2,y2 = parse.(Int,split(p2,","))
		if x1 == x2 || y1 == y2 # hor/ver lines
			for point in Iterators.product(x1<=x2 ? (x1:x2) : (x1:-1:x2), y1 <= y2 ? (y1:y2) : (y1:-1:y2))
				obs_horver[point] = get!(obs_horver, point,0) + 1
				obs_horverdiag[point] = get!(obs_horverdiag, point,0) + 1
			end 
		else # diag lines
			xrange = x1<=x2 ? (x1:x2) : (x1:-1:x2)
			y = map(x-> round(Int,y1 + (y2-y1)/(x2-x1) * (x-x1)), xrange)
			for point in zip(xrange, y)
				obs_horverdiag[point] = get!(obs_horverdiag, point,0) + 1
			end
		end
	end
	res_1 = count(x-> x[2] >= 2,obs_horver)
	res_2 = count(x-> x[2] >= 2,obs_horverdiag)
	(res_1, res_2)
end

# ╔═╡ 02ffd706-1616-46f8-b620-47322b22aa3f
md"""
## Day 6
"""

# ╔═╡ 452c9fc0-919f-4824-ab7b-cc4dcaf11bd4
let
	input = parse.(Int,split(readlines("./D6-1.txt")[1],","))

	# part 1 - does not scale well...
	function iteratefish!(x::Vector, n::Int=18)
		for _ = 1:n
			z = findall(x->isequal(x, 0), x)
			x = x .- 1
			x[z] .= 6
			for _ = 1:length(z)
				push!(x, 8)
			end
		end
		return length(x)
	end

	t = 80
	res_1 = iteratefish!(deepcopy(input), t)

	# part 2 
	# tried fitting exponential growth model xₜ = x₀(1+r)ᵗ => no succes :-(
	# better: use counter instead of actual 🐟

	# read start config
	out = Dict(x=> 0 for x in 0:8)
	for item in input
		out[item] += 1
	end

	function nextgen(d)
		future = Dict(x=>0 for x in 0:8)
		for x in 1:8
			future[x-1] = d[x]
		end
		future[6] += d[0]
		future[8] += d[0]
		
		return future
	end
	
	for _ = 1:256
		out = nextgen(out)
	end
	
	res_2 = sum(values(out))

	(res_1, res_2)
end

# ╔═╡ 1ad145a2-2554-4255-9a9e-f67723ad1a22
md"""
## Day 7
"""


# ╔═╡ 68396344-01cf-477a-a677-b0d6b434b7a4
let
	input = parse.(Int, split(readlines("./D7-1.txt")[1],","))
	
	distance_1(v::Vector{Int}, x::Int) = sum(abs.(v .- x))
	distance_2(v::Vector{Int}, x::Int) = sum(map(y -> y*(y+1) ÷ 2, abs.(v .- x)))
	
	function solve(v::Vector{Int}; d_fun::Function=distance_1)
		low = typemax(eltype(v))
		for x in minimum(v):maximum(v)
			if d_fun(v,x) < low
				low = d_fun(v,x)
			end
		end
		return low
	end
		
	(solve(input), solve(input, d_fun=distance_2))
end

# ╔═╡ 2f55f484-b0ff-4134-9214-12b81b0c1757
md"""
## Day 8

"""

# ╔═╡ 78a56791-b64a-4fec-b631-4474d5a93f47
let
	input = readlines("./D8-1.txt")
	# part 1
	res_1 = 0
	for line in input
		outnumbers = split(line," | ")[2]
		outnumbers = split(outnumbers," ")
		for item in outnumbers
			if length(item) == 2 || length(item) == 4 || length(item) == 7 || length(item) == 3
				res_1 += 1
			end
		end
	end

	# part 2	
	function findmap(inbr)
		# start inferring parts
		res = Dict{Symbol, Char}()
		one = collect(inbr[findfirst(x-> length(x) == 2,inbr)])
		four = collect(inbr[findfirst(x-> length(x) == 4,inbr)])
		seven = collect(inbr[findfirst(x-> length(x) == 3,inbr)])
		eight = collect(inbr[findfirst(x-> length(x) == 7,inbr)])
		
		# find top bar:
		res[:top] = setdiff(seven, one)[1]
		
		# identify zero/six/nine:
		zero_six_nine = collect.(inbr[findall(x-> length(x) == 6,inbr)])
		res[:bottom] = filter(x->length(x)==1,map(x -> setdiff(x, four, one, res[:top]), zero_six_nine))[1][1]
		res[:bottomleft] = setdiff(filter(x->length(x)==2,map(x -> setdiff(x, four, one, res[:top]), zero_six_nine))[1], res[:bottom])[1]
		
		# identify two/three/five:
		two_three_five = collect.(inbr[findall(x-> length(x) == 5,inbr)])	
		res[:mid] = filter(x->length(x)==1, map(x-> setdiff(x, one, res[:top], res[:bottom]), two_three_five))[1][1]
		res[:topleft] = filter(x->length(x)==1,map(x-> setdiff(x, one, res[:top], res[:bottom], res[:mid], res[:bottomleft]), two_three_five))[1][1]
		
		# identify bottomright/topright
		six_nine = filter(x-> res[:mid]∈x, zero_six_nine)
		res[:bottomright] = filter(x->length(x)==1, map(x->setdiff(x, values(res)...), six_nine))[1][1]
		res[:topright] = setdiff(eight, collect(values(res)))[1]

		# generate string map
		converter = Dict(
			sort([res[:top], res[:topleft], res[:topright], res[:bottomleft], res[:bottomright], res[:bottom]]) => "0",
			sort([res[:topright], res[:bottomright]])=>"1",
			sort([res[:top], res[:topright], res[:mid],res[:bottomleft], res[:bottom]]) => "2",
			sort([res[:top], res[:topright], res[:mid], res[:bottomright], res[:bottom]]) => "3",
			sort([res[:topleft], res[:topright], res[:mid], res[:bottomright]]) => "4",
			sort([res[:top], res[:topleft], res[:mid], res[:bottomright], res[:bottom]]) => "5",
			sort([res[:top], res[:topleft], res[:mid], res[:bottomleft], res[:bottomright], res[:bottom]]) => "6",
			sort([res[:top], res[:topright], res[:bottomright]]) => "7",
			sort([res[:top], res[:topleft], res[:topright], res[:mid], res[:bottomleft], res[:bottomright], res[:bottom]]) => "8",
			sort([res[:top], res[:topleft], res[:topright], res[:mid], res[:bottomright], res[:bottom]]) => "9",
		)
		return converter
	end

	res_2 = 0
	for line in input
		tobreak  = split(split(line, " | ")[1], " ")
		todecode = split(split(line, " | ")[2], " ")
		mapping = findmap(tobreak)
		broken = parse(Int,prod([mapping[sort(collect(val))] for val in todecode]))
		res_2 += broken
	end

	(res_1, res_2)
end

# ╔═╡ 9732dade-dcfc-4c36-bc48-7fa2490c42a8
md"""
## Day 9
"""

# ╔═╡ b7710d66-9a60-41b6-8477-2f7fbb1d1343
let 
	input = permutedims(hcat(map(x->parse.(Int,x),split.(readlines("./D9-1.txt"),""))...))

	# part 1
	# place input in larger matrix
	A = 9 * ones(Int, size(input) .+ 2)
	A[2:end-1,2:end-1] .= input
	# find all neighbors
	inds(i,j) = [CartesianIndex(v) for v in Iterators.product(i-1:i+1,j-1:j+1) if (v[1],v[2])≠(i,j)]
	# check if location is local min
	localmin(A,i,j) = count(x->x<A[i,j], @view A[inds(i,j)]) == 0
	# get the minima
	mins = findall(map(x-> localmin(A,x...), Iterators.product(2:size(A,1)-1, 2:size(A,2)-1)))
	# get the anser
	res_1 = sum(input[mins]) + length(mins)

	# part 2 - using recurrence
	checked = zeros(Bool, size(A))
	checked[A .== 9] .= true # holds the cells
	# find a starting points
	neigcoor(ind) = [CartesianIndex(Tuple(ind) .- diff) for diff in [(-1,0),(1,0),(0,-1),(0,1)]]

	function crawler(checked::Matrix{Bool}, startind::CartesianIndex{2}, comp::Vector{CartesianIndex{2}}=[startind])
		checked[startind] = true
		for neig in neigcoor(startind)
			if !checked[neig]
				push!(comp, neig)
				crawler(checked, neig, comp)
			end
		end

		return comp
	end
	
	basins = Vector{Vector{CartesianIndex{2}}}()
	while !isnothing(findfirst(iszero, checked))
		startpoint = findfirst(iszero, checked)
		minicomp = crawler(checked, startpoint)
		push!(basins, minicomp)
	end
	res_2 = prod(length.(sort(basins,by=length)[end-2:end]))

	(res_1,res_2)

end

# ╔═╡ 1a1d5dd1-39b8-4fae-9cd6-d551e0c18639
md"""
## Day 10
"""

# ╔═╡ 7cc8077e-63b8-4277-89d5-22104d0a3a55
let 
	input = readlines("./D10-2.txt")
	complement = Dict('(' => ')', '<' => '>', '[' => ']', '{' => '}',
						')' => '0', '>' => '0', '}' => '0', ']' => '0')
	scoring_corrupt = Dict(')' => 3, ']' => 57, '}' => 1197, '>' => 25137)
	scoring_complete = Dict('(' => 1, '[' => 2, '{' => 3, '<' => 4)
	
	function checkline(line)
		stack = Char[]
		for i in eachindex(line)
			push!(stack, line[i])
			if length(stack)>1
				if stack[end] == complement[stack[end-1]]
					pop!(stack); pop!(stack)
				end
			end
		end
		return stack
	end

	function scoreme(res, scores::Dict=scoring_complete)
		score = 0
		for char in reverse(res)
			score *= 5
			score += scores[char]
		end
		return score
	end

	# solution
	res_1 = 0
	completionscores = Int[]

	for line in input
		res_c = checkline(line)
		corrupt = any(map(x -> x ∈ unique(res_c), [')';'}';']';'>']))
		if corrupt
			winner = findfirst(x -> x ∈ [')';'}';']';'>'], res_c)
			res_1 += scoring_corrupt[res_c[winner]]
		else
			push!(completionscores, scoreme(res_c))
		end
	end
	res_2 = sort(completionscores)[floor(Int,length(completionscores)/2)+1]
	
	(res_1, res_2)
end

# ╔═╡ 7ed4f9d0-bbd5-4705-8345-626c7f3a9428
md"""
## Day 11
"""

# ╔═╡ c1b34f1b-722a-4016-84f6-f57a1a09c406
let
	input = permutedims(parse.(Int,hcat(collect.(readlines("./D11-2.txt"))...)))
	flashed = zeros(Bool, size(input))

	inds(i,j, n, m) = [CartesianIndex(v) for v in Iterators.product(i-1:i+1,j-1:j+1) if all([all((v[1],v[2])≠(i,j)), v[1] >= 1, v[2] >= 1, v[1] <= n, v[2] <= m])]
	
	function iterate!(A::Matrix{Int}, f::Matrix{Bool}, (n,m)::Tuple{Int,Int}=size(A))
		flashcount = 0
		f .= 0
		A .+= 1
		flashers = findall(x -> x > 9, A)
		while length(flashers) > 0
			selected = pop!(flashers)
			flashcount += 1
			f[selected] = true
			neigh = inds(selected[1], selected[2],n,m)
			for n in neigh
				A[n] += 1
				if A[n] > 9 && !f[n] && n∉flashers
				 	push!(flashers, n)
				end
			end
		end
		A[A .> 9] .= 0
		return A, flashcount
	end

	# part 1
	M = deepcopy(input)
	buffer = similar(M, Bool)
	res_1 = 0
	for i in 1:100
		M, c = iterate!(M, buffer)
		res_1 += c
	end
	# part 2
	res_2 = 100
	while true
		res_2 += 1
		M, _ = iterate!(M, buffer)
		if all(M .== 0)
			break
		end
	end
	
	(res_1, res_2)
end

# ╔═╡ c61e6dea-42d9-45f2-834d-897734987bb3
md"""
## Day 12
"""

# ╔═╡ 4acc98f4-068d-4f4c-b29e-2cd80c604528


# ╔═╡ Cell order:
# ╟─821803bc-93b1-436d-9adb-2a759a782034
# ╟─d152b170-8593-4128-a950-ec14a66af2c2
# ╠═b745b668-52a3-11ec-2566-c1d79a56cfc2
# ╟─c08f1079-ca65-458a-837d-ffbda5046bf3
# ╠═1f3a8a64-e78b-407e-8835-f9a0fa419c46
# ╟─f6e71aca-953a-4320-b3ab-be377e0aa22c
# ╠═942b6adf-426a-432e-9c04-116ece48bd9d
# ╟─134adeb8-1759-4fec-8a62-3906bb72c451
# ╠═5aa396ce-fed7-45fa-a2a5-6fe40a90bf84
# ╟─b8993212-0328-416b-8c27-6b6ae3aca8fd
# ╠═ad35bd5a-401a-44b5-9786-1b1744fabe34
# ╟─02ffd706-1616-46f8-b620-47322b22aa3f
# ╠═452c9fc0-919f-4824-ab7b-cc4dcaf11bd4
# ╟─1ad145a2-2554-4255-9a9e-f67723ad1a22
# ╠═68396344-01cf-477a-a677-b0d6b434b7a4
# ╟─2f55f484-b0ff-4134-9214-12b81b0c1757
# ╠═78a56791-b64a-4fec-b631-4474d5a93f47
# ╟─9732dade-dcfc-4c36-bc48-7fa2490c42a8
# ╠═b7710d66-9a60-41b6-8477-2f7fbb1d1343
# ╟─1a1d5dd1-39b8-4fae-9cd6-d551e0c18639
# ╠═7cc8077e-63b8-4277-89d5-22104d0a3a55
# ╟─7ed4f9d0-bbd5-4705-8345-626c7f3a9428
# ╠═c1b34f1b-722a-4016-84f6-f57a1a09c406
# ╟─c61e6dea-42d9-45f2-834d-897734987bb3
# ╠═4acc98f4-068d-4f4c-b29e-2cd80c604528
