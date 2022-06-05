from localisation import *

if __name__ == '__main__':
    wave=Wave_Packet(
        epsilon=0.0,    # amplitude of disorder
        spacing=0,      # spacing of disorder
        dt=0.25,        # length of each timestep
        x0=5,           # starting position of wave in x_range
        x_range=40,     # size of landscape, symmetrical so 40 means x is from -20 to 20 
        resolution=100, # how many steps the range in discretised into
        sigma0=1.5,     # shape of wave
        k0=3.0          # momentum of wave
        )

    wave.add_barrier(
        "[45:55]",            # location of barrier, in terms of resolution (must be integers)
        curve="abs(sin(x))",  # equation of curve, must use proper Python operators, e.g. "ax" -> "a*x"
        preview=True          # gives a preview and asks for input to verify if correct
        )
    wave.add_barrier("[10:11]", curve="100", preview=False)
    wave.add_barrier("[89:90]", curve="100", preview=False)

    viewer(
        wave, 
        max_steps=300,      # number of time steps to simulate for
        save=True,          # option to save to set path
        PATH=r"./test.gif"  # save path
        )