include("./linked_stack.jl")
using .Stack	
using Test

function populate_cliques(A)

	c = LinkedStack{Array{Integer}}();
	for i in 1:size(A,1)
		push!(c,[i]);
	end
	
	results = [];
	while(!isempty(c))
		item = pop!(c);
		largest = item[length(item)];
		for i in largest+1:size(A,1)
			# if is connected to all
			if sum(A[i,item]) == length(item)
				push!(c, vcat(item, i));
			end
		end
		push!(results, item);
	end
	
	return results;
end

function tests()
	testSimple();
	testBigger();
end

function testSimple()
	#arrange
	A = [1 0 0 0 0 0;
	0 1 0 0 0 0;
	0 0 1 0 0 0;
	0 0 0 1 0 0;
	0 0 0 0 1 0;
	0 0 0 0 0 1];
	
	st = LinkedStack{Array{Integer}}();
	push!(st, [1]);
	push!(st, [2]);
	push!(st, [3]);
	push!(st, [4]);
	push!(st, [5]);
	push!(st, [6]);
	
	#act
	res = populate_cliques(A);
	
	#assert
	for i in 1:length(res)
		@assert pop!(st,pop!(res));
	end
	@test isempty(res) && isempty(st);

end

function testBigger()
	#arrange
	A = [1 1 1 0 0 0 0;
		1 1 1 0 0 0 0;
		1 1 1 1 1 1 1;
		0 0 1 1 1 0 0;
		0 0 1 1 1 0 0;
		0 0 1 0 0 1 1;
		0 0 1 0 0 1 1];  
		
	st = LinkedStack{Array{Integer}}();
	push!(st, [1]);
	push!(st, [2]);
	push!(st, [3]);
	push!(st, [4]);
	push!(st, [5]);
	push!(st, [6]);
	push!(st, [7]);
	push!(st, [1,2]);
	push!(st, [2,3]);
	push!(st, [1,3]);
	push!(st, [3,4]);
	push!(st, [3,5]);
	push!(st, [4,5]);
	push!(st, [3,6]);
	push!(st, [6,7]);
	push!(st, [3,7]);
	push!(st, [1,2,3]);
	push!(st, [3,4,5]);
	push!(st, [3,6,7]);
	
	#act
	res = populate_cliques(A);
	#show(res);
	#assert
	for i in 1:length(res)
		@assert pop!(st,pop!(res));
	end
	@test isempty(res) && isempty(st);
end

#@time @testset "FindingCliques" begin tests() end


