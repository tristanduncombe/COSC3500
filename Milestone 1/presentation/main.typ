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
  #align(center)[
    #image("animation.gif", width: 40%)
  ]
  #set text(18pt)
  _Note that this is only displaying $1/10$th the frames generated!_
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
  Well, the difference was only $8 mu"s"$ across all loops, which is an improvement, but not substantial
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
  == Exit Loops Early
  By exiting loops early, we can avoid unncessary computations!

  Within the case of the electron simulation, if we create cut-offs when the distance or force is less than a specific amount, then we can avoid computing these cases.

  As the simulation progresses, this will add up!
]

#slide[
  == Exit Loops Early: Distance
  Given that we know that:

  #align(center)[$F("applied to electron") prop 1/r^2$]
  Once the distance is greater than an amount, then we can ignore these.

  The cases, $3, 5, 7$ were considered, if the distance between the electron was greater than one of these, it was ignored.
]

#slide[
  == Exit Loops Early: Distance Performance
  #align(center)[

  ]

  Evidently, choosing to exit at 3, 5, 7 improved performance at minimum, by around 40%, but choosing 3 improved performance by 76%. 
]

#slide[
  == Exit Loops Early: Distance Accuracy
  #align(center)[

  ]

  Hence, the inaccuracy over each iteration of 0.25cm was acceptable given the 20m by 20m volume, indicating perhaps an exit of 1m or 2m may have appropriate.
]


#slide[
  == Exit Loops Early: Force
  Another way in which we could exit the loop early would be to exit if the force was below a certain amount.

  The range of forces from 10N to 0.00001N at multiples of 10.
]

#slide[
  == Exit Loops Early: Force Performance
  #align(center)[

  ]

  All force cutoffs excluding 1N performed near identically as they were within 0.02s of each other, with 1N performing twice as worse.

  This performance makes little sense as it is expected that the time taken decreases as the force cut-off increases (as it would exclude more iterations), however, this performance (and the accuracy) did not function as expected, and hence, this optimisation was ignored.
]

#centered-slide[
  _Let's continue optimising!_
]


#slide[
  == Spacial Computing: Chunking

  To Further Optimise the Simulation, Spacial Computing was used; \
  \
  This is a technique which separtes space into chunks inorder to only consider those.

  Within the context of the electron interaction simulation mutltiple grid dimensions were tested.
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