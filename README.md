# Evolution in disordered quantum systems
 
 #### This project aims to provide the user the ability to visualise the time evolution of a quantum system in a disordered system of arbitrary potential. The user has the ability to add a Gaussian wavepacket with a set momentum and shape to a landscape of a given length and 'resolution'. In its current form, an arbitrary potential can be added to this landscape by the piecewise addition of various curves, allowing for the visualisation of the time evolution of the Gaussian curve in this potential. In addition, one can then perturb the potential by specifying an amount of disorder (and the spacing of this disorder), which allows the user to see the phenomenon of 'Anderson localisation'. The programme is designed to run on a Windows machine with Python 3.10.1 or higher.

The code can be installed by copying the *localisation.py* and *requirements.txt* packages into a single folder. Then the required packages are installed by running the command:
```
pip install -r requirements.txt
```
in the command line, using a suitable code editor.


For a guided example, see the examples.py file. Below are some a few animations which show the potential of this software. In particular we see: a free particle, a particle in infinite square well, a quantum harmonic osciallator, and a perturbed quantum harmonic oscilaltor. **More animations can be seen in the *Animations* folder in the repository.**
 
 **Free particle:**
 
<img src="https://github.com/emgsharp99/evolution-DQS/blob/main/Animations/free_part.gif" data-canonical-src="https://github.com/emgsharp99/evolution-DQS/blob/main/Animations/free_part.gif" width="750" height="750" />
 
 **Disorder free particle:**
 
 <img src="https://github.com/emgsharp99/evolution-DQS/blob/main/Animations/free_part_pert.gif" data-canonical-src="https://github.com/emgsharp99/evolution-DQS/blob/main/Animations/free_part.gif" width="750" height="750" />
  
 **(Quasi-)infinite square well:**
 
 <img src="https://github.com/emgsharp99/evolution-DQS/blob/main/Animations/inf_square_well.gif" data-canonical-src="https://github.com/emgsharp99/evolution-DQS/blob/main/Animations/inf_square_well.gif" width="750" height="750" /> 
 
 **Quantum harmonic oscillator:**
 
<img src="https://github.com/emgsharp99/evolution-DQS/blob/main/Animations/qho.gif" data-canonical-src="https://github.com/emgsharp99/evolution-DQS/blob/main/Animations/qho.gif" width="750" height="750" /> 
 
 **Disordered quantum harmonic oscillator:**
 
 <img src="https://github.com/emgsharp99/evolution-DQS/blob/main/Animations/qho_pert.gif" data-canonical-src="https://github.com/emgsharp99/evolution-DQS/blob/main/Animations/qho_pert.gif" width="750" height="750" /> 


## How to use the code:
Once the required packages are installed, we must import the code by including;
```
from localisation import *
```
in the code editor, or in the command line.

Next, you must create a *Wave_Packet* object:

```
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
```
 
 
This sets out the characteristics of the landscape and Gaussian wave packet that you wish to simulate.

Next, we can add an arbitrary curve to this landscape, which represents a potential of some function. We do this using the *add_barrier* method:
 
```
wave.add_barrier(
       "[70:230]",            # location of barrier, in terms of resolution (must be integers)
       curve="x**2/150",      # equation of curve, must use proper Python operators, e.g. "ax" -> "a*x"
       preview=False          # gives a preview and asks for input to verify if correct
       )
```
 
 
Any 1 dimensional function V(𝑥) can be used in the *curve* argument, and its location is set by the first of argument of the method. If the user wishes to preview their updated potential landscape before running the full simulation, they simply include *preview=True* as an argument.

Once the user is happy with the potential that they have set up, they can finally call the *viewer* function to generate an animation of the time evolution of the wave function:
 
```
viewer(
    wave, 
    max_steps=300,      # number of time steps to simulate for
    save=True,          # option to save to set path
    PATH=r"./PATH.gif"  # save path
    )
``` 
 
The arguments are fairly straightforward. We select the number of steps we wish to simulate, and whether we wish to save the animation. If the user wishes to save the animation, then they must specify a save path. The default save path is ./test.gif. If the animation is not saved, then it will appear in the code editor.
