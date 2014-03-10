Wave1D
======

Experiments with 1 dimensional waves. I am using this project to learn to use
git, and GitHub.

The system is initialized with a square wave. The UI provides several choices
of problems, and discretization schemes. The wave speed is also a user
controlled parameter. FLTK and OpenGL are used for the UI.

# Wave Equations
The program solves the 1D hyperbolic partial differential equation: 
![](http://www.sciweavers.org/upload/Tex2Img_1394419013/render.png)

It also solves the wave advection equation:
![](http://www.sciweavers.org/upload/Tex2Img_1394419106/render.png)

# Discretization Schemes

The following numerical discretization schemes are used to solve the equations.

## Wave Equation:
1. Explicit time, centered space (unstable).
2. Lax scheme

## Advection equation:
1. Explicit time, centered space (unstable).
2. Lax scheme
3. Upwind differencing
4. BFECC
5. Modified BFECC
6. MacCormack
7. Flux limited MacCormack
