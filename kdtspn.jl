#!/usr/bin/env julia

module kDTSPN

using Dubins

struct Individual
	Sequence::Vector{Int32}
	StartingPoints::Vector{Int32}
end

function crossover(ind1:Individual, ind2:Individual, point:Int32)
	if	point > length(ind1.Sequence)
		error("Crossover point outside possible range");
	end
	
	new_ind1 = copy(ind1[1:point]);
	left_ind1 = copy(ind1[point+1:end]);
	for i in point+1:length(ind1)
		
	end
	
	
end

function getNextIndexNotUsed(arr:Vector{Int32})
	for i in 1:length(arr)
		if(arr[i] != 0)
			return i;
		end
	end
end

function doSomething()
	a = Individual([3,4], [3,5]);
	#a.Sequence = [3,4];
	#a.StartingPoints = [3,5];
	println(a.Sequence, a.StartingPoints)
end

doSomething();

end