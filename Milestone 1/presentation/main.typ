#import "@preview/polylux:0.3.1": *

#import themes.simple: *

#set text(font: "Inria Sans")

#show: simple-theme.with(
  footer: [Simple slides],
)

#title-slide[
  = COSC3500: Optimising Electron Interaction Simulation
  #v(2em)

  Tristan Duncombe s470298 \
  July 23
]

#slide[
  == First slide

  #lorem(20)
]

#focus-slide[
  _Focus!_

  This is very important.
]

#centered-slide[
  = Let's start a new section!
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
  == Days Spent Fooling Around Wondering Why Data Wasn't Being Updating Correctly but I Slow Down Time which Meant that it Didn't Move Quickly

  #align(center)[2]
  #align(center)[Let's not discuss this any further]
]

#slide[
  == Dynamic slide
  Did you know that...

  #pause
  ...you can see the current section at the top of the slide?
]