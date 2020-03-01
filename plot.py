import matplotlib.pyplot as plt

def plot(Riemann, Solver_solution, grid, t):

    #Riemann problem solution
    Pressure = [Riemann.eval_sampling_point(x/t).pressure for x in grid.cell_position]
    Velocity = [Riemann.eval_sampling_point(x/t).velocity for x in grid.cell_position]
    Density = [Riemann.eval_sampling_point(x/t).rho for x in grid.cell_position]

    #Euler solver output
    Solver_density = [Solver_solution.U_final[i].rho for i in range(len(grid.cell_position))]
    Solver_velocity = [Solver_solution.U_final[i].velocity for i in range(len(grid.cell_position))]
    Solver_pressure = [Solver_solution.U_final[i].pressure for i in range(len(grid.cell_position))]

    fig, axs = plt.subplots(3, sharex=True)
    fig.suptitle("Solution of the Riemann problem\nat t = {}s".format(t))
    axs[0].plot(grid.cell_position, Density)
    axs[0].plot(grid.cell_position, Solver_density, "+")
    axs[1].plot(grid.cell_position, Velocity)
    axs[1].plot(grid.cell_position, Solver_velocity, "+")
    axs[2].plot(grid.cell_position, Pressure)
    axs[2].plot(grid.cell_position, Solver_pressure, "+")

    axs[0].grid()
    axs[0].set(ylabel = "Density")
    axs[1].grid()
    axs[1].set(ylabel = "Velocity")
    axs[2].grid()
    axs[2].set(ylabel = "Pressure")

    plt.xlabel("Location x")
