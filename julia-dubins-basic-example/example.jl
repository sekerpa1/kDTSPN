#!/usr/bin/env julia

module Example

using Printf
using PyPlot
using TimerOutputs
using Dubins
using Random

# Circle from matplotlib for plotting
Circle = matplotlib.patches.Circle

# Create the timer object
to = TimerOutput()

################################################
# Settings
################################################
problem_file = "./problems/gdip-n10.txt"
turning_radius = 0.7
sensing_radius = 1.0

##################################################
# Functions
##################################################
"""
    load_map(filename)

Read config with goal positions
(This is standard way to document function in Julia. 
It can be shown by typing "?Example.load_map" in REPL.)
"""
function load_map(filename)
    goals = Vector{Vector{Float64}}()
    for line in readlines(open(filename))
        label, x, y = split(line)
        point::Vector{Float64} = [parse(Float64,x), parse(Float64,y)]
        push!(goals, point) 
    end
    return goals
end

function plot_points(points; specs = "b")
    x_val = [x[1] for x in points]
    y_val = [x[2] for x in points]
    plt.plot(x_val, y_val, specs)	 
end 

function plot_circle(xy, radius)
    ax = plt.gca()
    circle = Circle(xy, radius, facecolor="yellow", edgecolor="orange", 
        linewidth=1, alpha=0.2)
    ax.add_patch(circle)
end

function dist_euclidean(coord1, coord2)
    (x1, y1) = coord1
    (x2, y2) = coord2
    (dx, dy) = (x2 - x1, y2 - y1)
    return sqrt(dx * dx + dy * dy)
end

function plot_map(goals, sensing_radius)
    plt.figure(1)
    plt.clf()
    for g in goals
        plot_circle(g, sensing_radius)
        plot_points([g]; specs="rx")
    end
end

function plot_trajectory(trajectory::Vector{DubinsPath})
    for t in trajectory
        samples = dubins_path_sample_many(t, turning_radius * 0.1)[2]
        plot_points(samples)
    end
end

function get_trajectory_length(trajectory::Vector{DubinsPath})
    sum = 0.0
    for t in trajectory
        sum += dubins_path_length(t)
    end
    return sum
end

##################################################
# Main code
##################################################
function solve()
    @timeit to "load" goals = load_map(problem_file)
    n = length(goals)
    @show goals
    # clear the figure 1, and plot the problem
    plot_map(goals, sensing_radius)

    # random permutation of goals
    sequence = randperm(n)
    # random heading angles
    headings = 2pi * rand(n)

    # final trajectory is a array of Dubins paths
    trajectory = DubinsPath[]
    for i in 1:n 
        idx1 = sequence[i]
        idx2 = sequence[mod1(i+1, n)]
        q1 = [goals[idx1]..., headings[idx1]]
        q2 = [goals[idx2]..., headings[idx2]]
		println(q1);
		println(q2);
        @timeit to "dubins" begin
            dubins = dubins_shortest_path(q1, q2, turning_radius)[2]
        end
		println(dubins.params);
        push!(trajectory, dubins)
    end

    plot_trajectory(trajectory)

    traj_length = get_trajectory_length(trajectory)

    plt.title("Length = " * string(traj_length))
    plt.axis("equal")
    plt.tight_layout()
    plt.pause(0.1)

    # save into a file
    plt.savefig("trajectory.pdf")
end

@timeit to "first" solve()
@timeit to "second" solve()

# show computational times
@show to

end