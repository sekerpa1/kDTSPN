using PyPlot

Circle = matplotlib.patches.Circle

struct Point{T<:AbstractFloat}
	x::T
	y::T
end

Point(p::Tuple{T,T}) where T <: AbstractFloat = Point{T}(p[1], p[2]);

struct Region{T<:AbstractFloat}
	center::Point{T}
	radius::T
end

function plotCircle(xy, radius)
    pyplot();
    circle = Circle(xy, radius, facecolor="yellow", edgecolor="orange",
        linewidth=1, alpha=0.2);
end


minimum = 3;
maximum = 97;


function plot_circle_and_points()

#Finding projected point for : [ 6.71627 , 8.3767] , [ 10.6659, 6.6348] and point
# c : [ 6.40217, 9.49307],
#         and region: [3.41556 , 10.2848 , 0.5]

	r = Region(Point(3.41556,10.2848), 0.5);
	a = Point(6.71627,8.3767);
	b = Point(10.6659,6.6348);
	c = Point(6.40217,9.49307);
	
	#projected = projected_point(a, b, c, r);
	#println(projected);
	
	xy = [r.center.x, r.center.y];
	
	plt.figure(1)
    plt.clf()
	ax = plt.gca();
    circle = Circle(xy, r.radius, facecolor="red", edgecolor="red", 
        linewidth=1, alpha=0.2);
    ax.add_patch(circle);
	plt.scatter([a.x, b.x, c.x],[a.y, b.y, c.y]);
	
	plt.title("Length = ")
    plt.axis("equal")
    plt.tight_layout()
    plt.pause(0.1)

    # save into a file
    #plt.savefig("something.pdf")
end

function generate_polygon_equations(x)
	XS = x[1];
	YS = x[2];
	s = size(x[1],1);
	# only 5 equations at maximum for polygon instead of 6 points + polygon flag
	s2 = size(x[1],2)-2;
	
	equations = [(0.0,0.0,0.0) for i in 1:s, j in 1:s2];
	
	for i in 1:s
		if XS[i,1] == 3
			eq = [get_line_equation(XS[i,j:j+1], YS[i,j:j+1]) for j in 2:4];
			equations[i,1:3] = eq;
		elseif XS[i,1] == 4
			#square
			eq = [get_line_equation(XS[i,j:j+1], YS[i,j:j+1]) for j in 2:5];
			equations[i,1:4] = eq;
		else
			#pentagon
			eq = [get_line_equation(XS[i,j:j+1], YS[i,j:j+1]) for j in 2:6];
			equations[i,1:5] = eq;
		end
	end
	
	return equations;
end

function generate_polygons()

	polygon_count = 100;
	XS = rand(polygon_count, 7);
	YS = rand(polygon_count, 7);
	
	for i in 1:polygon_count
		polygon_type = rand(3:5);
		size = rand(1:4);
		diff_x = rand(minimum:0.001:maximum);
		diff_y = rand(minimum:0.001:maximum);
		xs = zeros(7);
		ys = zeros(7);
		if polygon_type == 3
			#generate triangle
			xs[1:5] = [3.0,0.0,1.0,0.5,0.0];
			ys[1:5] = [3.0,0.0,0.0,1.0,0.0];
			xs[2:5] = xs[2:5] .* size .+ diff_x;
			ys[2:5] = ys[2:5] .* size .+ diff_y;
			xs[1] = 3;
			ys[1] = 3;
		elseif polygon_type == 4
			#generate square
			xs[1:6] = [4.0,0.0,1.0,1.0,0.0,0.0];
			ys[1:6] = [4.0,0.0,0.0,1.0,1.0,0.0];
			xs[2:6] = xs[2:6] .* size .+ diff_x;
			ys[2:6] = ys[2:6] .* size .+ diff_y;
			xs[1] = 4;
			ys[1] = 4;
		else 
			#generate pentagon
			xs[1:7] = [5.0,0.5,1.5,2.0,1.0,0.0,0.5];
			ys[1:7] = [5.0,0.0,0.0,1.0,2.0,1.0,0.0];
			xs[2:7] = xs[2:7] .* size .+ diff_x;
			ys[2:7] = ys[2:7] .* size .+ diff_y;
			xs[1] = 5;
			ys[1] = 5;
		end
		XS[i,1:7] = xs;
		YS[i,1:7] = ys;
	end
	
	return [XS,YS];
end


function triangle(min, max,size)
	xs = [0.0,1.0,0.5,0.0];
	ys = [0.0,0.0,1.0,0.0];
	diff_x = rand(min:0.001:max);
	diff_y = rand(min:0.001:max);
	xs = xs .*size .+ diff_x;
	ys = ys .*size .+ diff_y;
	return collect(zip(xs,ys));
end

function pentagon(min, max,size)
	xs = [0.5,1.5,2.0,1.0,0.0,0.5];
	ys = [0.0,0.0,1.0,2.0,1.0,0.0];
	diff_x = rand(min:0.001:max);
	diff_y = rand(min:0.001:max);
	xs = xs .*size .+ diff_x;
	ys = ys .*size .+ diff_y;
	return collect(zip(xs,ys));
end

function square(min, max,size)
	xs = [0.0,1.0,1.0,0.0,0.0];
	ys = [0.0,0.0,1.0,1.0,0.0];
	diff_x = rand(min:0.001:max);
	diff_y = rand(min:0.001:max);
	xs = xs .*size .+ diff_x;
	ys = ys .*size .+ diff_y;
	return collect(zip(xs,ys));
end

function get_line_equations(xs, ys)
	eq = [get_line_equation(xs[i:i+1],ys[i:i+1]) for i in 1:length(xs)-1];
	return eq;
end

#=
  Returns line equation in form [-a, +1/0, -b]
  -a*x + y - b = 0
=#
function get_line_equation(xs, ys)
	if xs[1] == xs[2]
		return (1.0, 0.0, -xs[1]);
	elseif ys[1] == ys[2]
		return (0.0, 1.0, -ys[1]);
	else
		a = (ys[2] - ys[1]) / (xs[2] - xs[1]);
		b = ys[2] - a*xs[2];
		return (-a, 1.0 ,-b);
	end
end



function random_colors(n)
	c = rand(0:0.01:1, n, 3);
	return c;
end

#=function plot_triangles()
	lines = Any[];
	for i in 1:3
	   push!(lines, triangle(3,97));
	end
	
	c
end=#

function check_inequation(equation, sign, x, y)
	if sign > 0
		return equation[1]*x + equation[2] * y + equation[3] >= 0;
	else
		return equation[1]*x + equation[2] * y + equation[3] <= 0;
	end
end


function plot_polygon(x)
	
	lines = Any[];
	for i in 1:size(x[1],1)
		if x[1][i,1] == 3
			#triangle
			push!(lines, collect(zip(x[1][i,2:5],x[2][i,2:5])));
		elseif x[1][i,1] == 4
			#square
			push!(lines, collect(zip(x[1][i,2:6],x[2][i,2:6])));
		else
			#pentagon
			push!(lines, collect(zip(x[1][i,2:7],x[2][i,2:7])));
		end
	end
	
	c = random_colors(100);
	line_widths = [2 for i in 1:100];
	
	line_segments = matplotlib.collections.LineCollection(lines, linewidths=line_widths, colors=c);
	
	fig = figure("Line Collection Example",figsize=(100,100));
	ax = PyPlot.axes();
	ax.add_collection(line_segments);
	axis("image");
	
end

function plot_stuff()
	
	plot_polygon();
end

function basic_plot()

	xs = [1.0,2.0,5.0,0.0];
	ys = [2.0,4.0,6.0,7.0];
	lines = Any[collect(zip(xs,ys))];
	
	xs = [3.0,4.0];
	ys = [5.0,6.0];
	push!(lines, collect(zip(xs,ys)));
	println(lines);
	
	c = Vector{Int}[[1,0,0],[0,1,0],[0,0,1]];
	
	line_segments = matplotlib.collections.LineCollection(lines,colors=c)
	
	fig = figure("Line Collection Example",figsize=(10,10));
	ax = PyPlot.axes();
	ax.add_collection(line_segments);
	axis("image");

    #=x = [i for i in 1:100];
    y = [3*i for i in 1:100];
    println(x);
    println(y);
    p = plot(x,y)
    xlabel("X")
    ylabel("Y")
    PyPlot.title("Your Title Goes Here")
    grid("on")=#
end

function inside_square()
	xs = [0.0,1.0,1.0,0.0,0.0];
	ys = [0.0,0.0,1.0,1.0,0.0];
	p_x = 0.0;
	p_y = 0.0;
	
	eq1 = get_line_equation(xs[1:2], ys[1:2]);
	eq2 = get_line_equation(xs[2:3], ys[2:3]);
	eq3 = get_line_equation(xs[3:4], ys[3:4]);
	eq4 = get_line_equation(xs[4:5], ys[4:5]);
	println(eq1);
	println(eq2);
	println(eq3);
	println(eq4);
	isin = [check_inequation(eq1, 1, p_x, p_y), check_inequation(eq2, -1, p_x, p_y),
		check_inequation(eq3, -1, p_x, p_y), check_inequation(eq4, 1, p_x, p_y)];
	println(isin);
end

function inside_triangle()
	xs = [0.0,1.0,0.5,0.0];
	ys = [0.0,0.0,1.0,0.0];
	p_x = 0;
	p_y = 2.5;
	
	eq1 = get_line_equation(xs[1:2], ys[1:2]);
	eq2 = get_line_equation(xs[2:3], ys[2:3]);
	eq3 = get_line_equation(xs[3:4], ys[3:4]);
	println(eq1);
	println(eq2);
	println(eq3);
	isin = [check_inequation(eq1, 1, p_x, p_y), check_inequation(eq2, -1, p_x, p_y),
		check_inequation(eq3, -1, p_x, p_y)];
	println(isin);
end

function inside_pentagon()
	xs = [0.5,1.5,2.0,1.0,0.0,0.5];
	ys = [0.0,0.0,1.0,2.0,1.0,0.0];
	p_x = 1.5;
	p_y = 1.5;
	
	eq1 = get_line_equation(xs[1:2], ys[1:2]);
	eq2 = get_line_equation(xs[2:3], ys[2:3]);
	eq3 = get_line_equation(xs[3:4], ys[3:4]);
	eq4 = get_line_equation(xs[4:5], ys[4:5]);
	eq5 = get_line_equation(xs[5:6], ys[5:6]);
	
	println(eq1);
	println(eq2);
	println(eq3);
	println(eq4);
	println(eq5);
	
	isin = [check_inequation(eq1, 1, p_x, p_y), check_inequation(eq2, 1, p_x, p_y),
		check_inequation(eq3, -1, p_x, p_y), check_inequation(eq4, -1, p_x, p_y), 
		check_inequation(eq5, 1, p_x, p_y)];
	println(isin);
	
end

function point_inside_polygon(x, y)
	inside = [];
	for i in 1:size(ex,1)
		if XS[i,1] == 3
			#triangle
			result = [check_inequation(ex[i,1], 1, x, y), check_inequation(ex[i,2], -1, x, y),
					check_inequation(ex[i,3], -1, x, y)];
			result = (sum(result) == 3);
		elseif XS[i,1] == 4
			#square
			result = [check_inequation(ex[i,1], 1, x, y), check_inequation(ex[i,2], -1, x, y),
					check_inequation(ex[i,3], -1, x, y), check_inequation(ex[i,4], 1, x, y)];
			result = (sum(result) == 4);
			
		else
			#pentagon
			result = [check_inequation(ex[i,1], 1, x, y), check_inequation(ex[i,2], 1, x, y),
					check_inequation(ex[i,3], -1, x, y), check_inequation(ex[i,4], -1, x, y),
					check_inequation(ex[i,5], 1, x, y)];
			result = (sum(result) == 5);
		end
		push!(inside, result);
	end
	
	return inside;
end


plot_circle_and_points();

#basic_plot();
#plot_stuff();

#xs = [0.5,1.5,2.0,1.0,0.0,0.5];
#ys = [0.0,0.0,1.0,2.0,1.0,0.0];

#eq = get_line_equations(xs,ys);
#x = 3;
#y = 1;
#check_inequation(eq[1], 1, x, y);

#x = generate_polygons();
#XS = x[1];
#YS = x[2];
#plot_polygon(x);
#ex = generate_polygon_equations(x);
#point_inside_polygon();


#point_inside_polygon(3.0,4.0);
#println(point_inside_polygon(centroid(XS[5,2:4]),centroid(YS[5,2:4])));
#println(point_inside_polygon(centroid(XS[10,2:4]),centroid(YS[10,2:4])));

#println(ex[15,1:end]);
#println(XS[15,1:end]);
#println(YS[15,1:end]);
#println(centroid(XS[15,2:4]));
#println(centroid(YS[15,2:4]));

#println(centroid(YS[15,2:4]));
#println(centroid(XS[15,2:4]));


#a = point_inside_polygon(centroid(XS[20,2:4]),centroid(YS[20,2:4]));
#println(a .> 0);
#println(point_inside_polygon(centroid(XS[25,2:4]),centroid(YS[25,2:4])));
#println(point_inside_polygon(centroid(XS[30,2:4]),centroid(YS[30,2:4])));
#println(point_inside_polygon(centroid(XS[35,2:4]),centroid(YS[35,2:4])));

function centroid(xs)
	return sum(xs)/length(xs);
end

#centroid([0.0,1.0,2.0,2.5]);

#inside_square();
#inside_triangle();
#inside_pentagon();