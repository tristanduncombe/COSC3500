#import "@preview/polylux:0.3.1": *

#import themes.simple: *

#set text(font: "Inria Sans")

#show: simple-theme.with(
  footer: [COSC3500: Optimising Serial Electron Interaction Simulation],
)

#title-slide[
  = COSC3500: Optimising Electron Interaction Simulation
  #v(2em)

  Tristan Duncombe s470298 \
  July 23
]

#slide[
  == The Project

  Before the project was implemented, the project was specified.

  The project was consistent of simulating electron interactions in a void. 

  Initially this consisted of 100 electrons over the course of 2 seconds (over intervals of 0.1 seconds), but this extended to 1000 electrons over 50 seconds.

  This was visualised using a python script.
]

#slide[
  == Initial Implementation

  To initially implement the simulation:
  - All 100 electrons were randomly placed within a 1m by 1m area
  - For each 0.1 period of time
  - Every electron calculates the force imposed on it.
  - Using trignometry we can calculate use this force to see how it would change the X,Y,Z position.

  
]

#slide[
  == Initial Implementation: Math
  To implement this simulation, we will use Columb's law:

  #align(center)[
    $|F| = k_e (|q_1||q_2|)/(r^2)$
  ]

  Where $k_e$ is Columb's Constant, $q_1$, $q_2$ are the quantities charge, and r is distance between the charge.

  It is important to notice that $|F|$ is inversely proportional to the distance squared.
]

#slide[
  == Initial Implementation: Code
  ```cpp
  int main() {
    for num in range(0, timeSteps) {
      for num in range(0, numElectrons) {
        for num in range(0, numElectrons - 1) {
          calculateForceToElectron()
        }
        ApplyForce()
      }
    }
  }
  ```
]

#slide[
  == Initial Implementation: Visualisation
  
]

#focus-slide[
  _So Let's Optimise!_
]

#slide[
  == Removing Loop Invariant Code

  A simple optimise technique is to remove invariant code:

  #columns(2)[
    #set text(16pt)
      ```cpp
      for (...) {
         const float m = 9.1093837 * pow(10,-31);
          const float k = 8.987 * pow(10, 9);
          const float e = 1.602 * pow(10, -19);
          const float t = 0.01;
      }
    ```
    \
    \
    \
    \
    \
      ```cpp
      const float m = 9.1093837 * pow(10,-31);
      const float k = 8.987 * pow(10, 9);
      const float e = 1.602 * pow(10, -19);
      const float t = 0.01;

      for (...) {
         ...
      }
    ```
  ]
]

#slide[
  == Removing Loop Invariant Code: Performance
  #align(center)[
    #image("remove_loop_chart.png", width: 40%)
  ]
  Impressive, right?
  Well, the difference was only $8 times 10^(-6)$, which is an improvement, but not substancial
]

#slide[
  == Days Spent Fooling Around Wondering Why Data Wasn't Being Updating Correctly but I Slow Down Time which Meant that it Didn't Move Quickly

  #align(center)[2]
  #align(center)[Let's not discuss this any further]
]

#centered-slide[
  _Let's continue optimising!_
]

#slide[
  == Using Spacial Computing to Optimise Electron Interaction Simulation

  To Further Optimise the Simulation, Spacial Computing was used; \
  \
  #align(center)[As $F("applied to electron") prop 1/r^2$] \
  We can ignore some electrons that are further from the currently evaluated electron.
]
#slide[
  == Spacial Computing: Chunking the Simulation

  If we consider the area within our simulation as distinct chunks we can choose to ignore specific areas and generalise across these areas;
]

#slide[
  == Dynamic slide
  Did you know that...

  #pause
  ...you can see the current section at the top of the slide?
]