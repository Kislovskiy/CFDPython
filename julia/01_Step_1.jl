### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ 31efc376-fb0b-11ea-3855-7dfe8b07c97b
md"
# Step 1: 1-D Linear Convection

The 1-D Linear Convection equation is the simplest, most basic model that can be used to learn something about CFD. It is surprising that this little equation can teach us so much! Here it is:

$$\frac{\partial u}{\partial t} + c \frac{\partial u}{\partial x} = 0$$

With given initial conditions (understood as a *wave*), the equation represents the propagation of that initial *wave* with speed $c$, without change of shape. Let the initial condition be $u(x,0)=u_0(x)$. Then the exact solution of the equation is $u(x,t)=u_0(x-ct)$.

We discretize this equation in both space and time, using the Forward Difference scheme for the time derivative and the Backward Difference scheme for the space derivative. Consider discretizing the spatial coordinate $x$ into points that we index from $i=0$ to $N$, and stepping in discrete time intervals of size $\Delta t$.

From the definition of a derivative (and simply removing the limit), we know that:

$$\frac{\partial u}{\partial x}\approx \frac{u(x+\Delta x)-u(x)}{\Delta x}$$

Our discrete equation, then, is:

$$\frac{u_i^{n+1}-u_i^n}{\Delta t} + c \frac{u_i^n - u_{i-1}^n}{\Delta x} = 0 $$

Where $n$ and $n+1$ are two consecutive steps in time, while $i-1$ and $i$ are two neighboring points of the discretized $x$ coordinate. If there are given initial conditions, then the only unknown in this discretization is $u_i^{n+1}$.  We can solve for our unknown to get an equation that allows us to advance in time, as follows:

$$u_i^{n+1} = u_i^n - c \frac{\Delta t}{\Delta x}(u_i^n-u_{i-1}^n)$$

Now let's try implementing this in Julia.
"

# ╔═╡ 7ff6c55c-fb1b-11ea-0f97-e5173989b913
md"
First, let's define a few variables; we want to define an evenly spaced grid of points within a spatial domain that is 2 units of length wide, i.e., $x_i\in(0,2)$.  We'll define a variable `NUMBER_OF_GRID_POINTS`, which will be the number of grid points we want and $\Delta x$ will be the distance between any pair of adjacent grid points.
"

# ╔═╡ 2576940e-fb0c-11ea-1c5f-9540a8875719
begin
	NUMBER_OF_GRID_POINTS = 81
	X_LEFT, X_RIGHT = .0, 2.0
	NUMBER_OF_TIMESTEPS = 25
	
	domain_length = abs(X_RIGHT - X_LEFT)

	Δx = 2 / (NUMBER_OF_GRID_POINTS - 1)
	Δt = .025
	c = 1
end

# ╔═╡ e77c8d7e-fb1b-11ea-1016-996e3d924a2f
md"
We also need to set up our initial conditions. The initial velocity $u_0$ is given as 
$u = 2$ in the interval $0.5 \leq x \leq 1$  and $u = 1$ everywhere else in $(0,2)$ (i.e., a hat function).

Here, we use the function `ones()` defining an array which is `NUMBER_OF_GRID_POINTS` elements long with every value equal to 1.
"

# ╔═╡ b682ef82-fb0d-11ea-1e1f-c989cf7e1122
begin
	x = range(X_LEFT, X_RIGHT, step=Δx)

	u = ones(NUMBER_OF_GRID_POINTS)
	
	step_left_idx = floor(Int, .5 / Δx)
	step_right_idx = floor(Int, 1 / Δx + 1)
	
	u[step_left_idx:step_right_idx] .= 2
end

# ╔═╡ 8a041a9c-fb14-11ea-197c-ab80ab7e935f
begin
	using Plots
	plot(x, u, 
		title = "Initial conditions", 
		lw = 3,
		xlabel = "x",
		ylabel = "u",
		leg = false)
end

# ╔═╡ 1385d574-fb1c-11ea-17d6-6129086a307f
md"
Now let's take a look at those initial conditions using Plots.  We've imported the `Plots` plotting library and the plotting function is called `plot`, so we'll call `plot`. To learn about the myriad possibilities of Matplotlib, explore the [Gallery](https://docs.juliaplots.org/latest/graphrecipes/examples/) of example plots.

Here, we use the syntax for a simple 2D plot: `plot(x,y)`, where the `x` values are evenly distributed grid points:
"

# ╔═╡ 62d3534a-fb1c-11ea-04fc-dbd4742455cd
md"
Now it's time to implement the discretization of the convection equation using a finite-difference scheme.  

For every element of our array `u`, we need to perform the operation $u_i^{n+1} = u_i^n - c \frac{\Delta t}{\Delta x}(u_i^n-u_{i-1}^n)$

We'll store the result in a new (temporary) array `un`, which will be the solution $u$ for the next time-step.  We will repeat this operation for as many time-steps as we specify and then we can see how far the wave has convected.  

We first initialize our placeholder array `un` to hold the values we calculate for the $n+1$ timestep, using once again the function `ones()`.

Then, we may think we have two iterative operations: one in space and one in time (we'll learn differently later), so we'll start by nesting one loop inside the other. Note the use of the nifty `range()` function. When we write: `for i in range(2, length = NUMBER_OF_GRID_POINTS - 1)` we will iterate through the `u` array, but we'll be skipping the first element.  *Why?*
"

# ╔═╡ b2cce0e0-fb16-11ea-334b-f9204d71b63e
begin
	un = ones(NUMBER_OF_GRID_POINTS)
	
	for n in range(1, length = NUMBER_OF_TIMESTEPS + 1)
		un = copy(u)
	    for i in range(2, length = NUMBER_OF_GRID_POINTS - 1)
	        u[i] = un[i] - c * Δt / Δx * (un[i] - un[i-1])
		end
	end
end

# ╔═╡ b2461e3a-fb1c-11ea-1f65-ef632cae4b32
md"
**Note**—We will learn later that the code as written above is quite inefficient, and there are better ways to write this.

Now let's try plotting our `u` array after advancing in time.
"

# ╔═╡ 1958374a-fb19-11ea-105a-d57b8a4aecaf
plot(x, u, 
	title = "Solution",
	lw = 3,
	xlabel = "x",
	ylabel = "u",
	leg = false)

# ╔═╡ Cell order:
# ╟─31efc376-fb0b-11ea-3855-7dfe8b07c97b
# ╟─7ff6c55c-fb1b-11ea-0f97-e5173989b913
# ╠═2576940e-fb0c-11ea-1c5f-9540a8875719
# ╟─e77c8d7e-fb1b-11ea-1016-996e3d924a2f
# ╠═b682ef82-fb0d-11ea-1e1f-c989cf7e1122
# ╟─1385d574-fb1c-11ea-17d6-6129086a307f
# ╠═8a041a9c-fb14-11ea-197c-ab80ab7e935f
# ╟─62d3534a-fb1c-11ea-04fc-dbd4742455cd
# ╠═b2cce0e0-fb16-11ea-334b-f9204d71b63e
# ╟─b2461e3a-fb1c-11ea-1f65-ef632cae4b32
# ╠═1958374a-fb19-11ea-105a-d57b8a4aecaf
