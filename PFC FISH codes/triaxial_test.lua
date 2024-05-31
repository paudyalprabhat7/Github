;-------------------------------------------------------------------
;                Sleeved Triaxial Test of a Bonded 
;                Material
;-------------------------------------------------------------------
model new
model title 'Sleeved Triaxial Test of a Bonded Material'
;make repeatable by setting the random number seed 
model random 10001
;set the model domain for PFC balls and walls
model domain extent -2 2 condition destroy
;largestrain mode must always be on for coupled simulations
model large-strain on
;apply timestep scaling so the PFC timestep will be 1
model mechanical timestep scale

;define relevant parameters for the cylinder
[rad = 1.0]
[height = 3.0]
[segments = 6]
[halfLen = height/2.0]
[freeRegion = height/2.0*0.8]

;create the cylinder with the geometry logic
;first create arcs making a circle of edges
geometry edge create by-arc origin (0,0,[-halfLen]) ...
   start ([rad*(-1)],0,[-halfLen]) end (0,[rad*(-1)],[-halfLen]) ...
   segments [segments]
geometry edge create by-arc origin (0,0,[-halfLen]) ...
   start (0,[rad*(-1)],[-halfLen]) end ([rad],0,[-halfLen]) ...
   segments [segments]
geometry edge create by-arc origin (0,0,[-halfLen]) ...
   start ([rad],0,[-halfLen]) end (0,[rad],[-halfLen]) ...
   segments [segments]
geometry edge create by-arc origin (0,0,[-halfLen]) ...
   start (0,[rad],[-halfLen]) end ([rad*(-1)],0,[-halfLen]) ...
   segments [segments]
;extrude the edges to make a cylinder
geometry generate from-edges extrude (0,0,[height]) segments [segments*2]

;import the structural elements as a shell from the geometry logic
;and assign groups and properties
structure shell import from-geometry 'Default' element-type dkt-cst 
structure node group 'middle' range position-z [-freeRegion] [freeRegion]
structure node group 'top' range position-z [freeRegion] [freeRegion+1]
structure node group 'bot' range position-z [-freeRegion-1] [-freeRegion]
structure shell group 'middle' range position-z [-freeRegion] [freeRegion]
structure shell property isotropic (1e6, 0.0) thick 0.25 density 930.0
;initialize the structural element related data structures 
model cycle 0
;In order to use the STRUCTURE APPLY command 
;set the local system of each structural
;element to be pointing to the center of the triaxial cell
fish define setLocalSystem
    loop foreach local s struct.node.list()
        local p = struct.node.pos(s)
        local nid = struct.node.id.component(s)
        local mvec = vector(0,0,comp.z(p))
        zdir = math.unit(p-mvec)
        ydir = vector(0,0,1)
        command
            structure node system-local z [zdir] y [ydir] ...
                      range component-id [nid]
        endcommand
    endloop
   command
     structure node fix system-local
   endcommand
end
[setLocalSystem]
;make sure that local damping is active for the structural elements
structure damp local
;fix the structural element nodes for specimen creation
structure node fix velocity rotation

;create a wall representation of the structural element
wall-structure create
;parameter to set the platen width relative to the cylinder radius
[rad2 = rad*1.3]
;create the platens
wall generate name 'platenTop' polygon ([-rad2],[-rad2],[halfLen]) ...
                                       ([rad2],[-rad2],[halfLen])  ...
                                       ([rad2],[rad2],[halfLen])   ...
                                       ([-rad2],[rad2],[halfLen])
wall generate name 'platenBottom' polygon ([-rad2],[-rad2],[-halfLen]) ...
                                          ([rad2],[-rad2],[-halfLen])  ...
                                          ([rad2],[rad2],[-halfLen])   ...
                                          ([-rad2],[rad2],[-halfLen])
;set the wall resolution strategy to full and set the cutoff angle for 
;proximity contacts                                          
wall resolution full
wall attribute cutoff-angle 20

;set the ball modulus and generate a cloud of balls with arbitrary overlap
[ballemod = 1.0e8]
ball distribute box [-rad] [rad] [-rad] [rad] [-halfLen] [halfLen] ...
                porosity 0.3 radius 0.05 0.1 ...
                range cylinder end-1 (0,0,[-halfLen]) end-2 (0,0,[height]) ...
                      radius [rad*0.95]

;set the ball attributes
ball attribute density 2600 damp 0.8
;set the default contact behavior - 
;the deformability method sets properties of the 
;linear portion of the contact model
contact cmat default model linearpbond method deformability ...
                   emod [ballemod] kratio 1.0 
;allow the balls to rearrange, 
;nulling the linear and angular velocities every 100 cycles
model cycle 1000 calm 100

;provide some friction to kill additional energy 
;and apply this to all of the contacts
contact cmat default model linearpbond method deformability ...
                   emod [ballemod] kratio 1.0 property fric 0.3
contact cmat apply
;solve to a small average ratio
model solve ratio-average 1e-8

;bond the ball-ball contacts
contact method bond gap 1.0e-2 range contact type 'ball-ball'
;change the existing contact properties 
;so that the linear force is incremental and
;supply strength
contact property lin_mode 1 lin_force 0 0 0 pb_ten 2e5 pb_coh 2e6
;set the stiffness of the parallel-bond portion of the contact model 
contact method pb_deformability emod [ballemod] kratio 1.0 ...
        range contact type 'ball-ball'

;this set of operations removes all ball velocities 
;and all contact forces in the model
;so that the specimen is completely stress free and bonded
model calm
ball fix velocity spin
model cycle 2
ball free velocity spin

;free the nodes in the middle section of the sleeve 
;while keeping the top and bottom edges fixed
structure node free velocity rotation range group 'middle'

;function for calculating stress and strain as the platens are displaced
;these values will be recorded as a history
;first find the top and bottom platens 
[platenTop = wall.find('platenTop')]
[platenBottom = wall.find('platenBottom')]
;define some variables for the calculation
[failureStress = 0]
[currentStress = 0]
[failureStrain = 0]
[area = math.pi()*rad^2.0]
;define the stress FISH function to measure the stress and strain
fish define stress
    local topForce = math.abs(comp.z(wall.force.contact(platenTop))) 
    local botForce = math.abs(comp.z(wall.force.contact(platenBottom)))
    currentStress = 0.5*(topForce+botForce)/area
    stress = currentStress
    strain = (height - (comp.z(wall.pos(platenTop)) - ...
              comp.z(wall.pos(platenBottom))))/height * 100
    if failureStress <= currentStress
        failureStress = currentStress
        failureStrain = strain
    endif
end

;define the halt FISH function to stop cycling 
;as the platens displace and the material fails
fish define halt
    halt = 0
    if currentStress < failureStress * 0.85
        halt = 1
    endif
end

;define a FISH function to increase the isotropic confining pressure 
;in increments so that the bonded material does not break during loading
fish define rampUp(beginIn,ending,increment)
    command
        ball attribute displacement (0,0,0)
        structure node initialize displacement (0,0,0)
    endcommand
    begin = beginIn
    loop while (math.abs(begin) < math.abs(ending))
        begin = begin + increment
        command
            ;apply the confining stress
            structure shell apply [begin] range group 'middle'
            ;apply the same confining stress on the platens
            wall servo force (0,0,[begin*area]) activate true ...
                 range name 'platenTop'
            wall servo force (0,0,[-begin*area]) activate true ...
                 range name 'platenBottom'
model cycle 200
            model calm
        endcommand
    endloop
    command
        ;once the stress state has been installed cycle and turn off the servo 
        ;mechanism on the walls
model cycle 1000
        wall servo activate false
        wall attribute velocity (0,0,0) range name 'platenTop'
        wall attribute velocity (0,0,0) range name 'platenBottom'
    endcommand
end
;set the platen velocity
[platenVel = 0.000003]

;save the model for future use, including these FISH utility function, before
;any confinement has been applied
model save 'beforeApplication'

;-------------------------------------------------------------------
;test 1: perform a UCS test on the specimen
model restore 'beforeApplication'
structure shell delete
wall attribute velocity-z [-platenVel] range name 'platenTop'
wall attribute velocity-z [platenVel] range name 'platenBottom'
ball attribute displacement (0,0,0)
fish history stress
fish history strain
model solve fish-halt halt
model save 'ucs'
[io.out(string(failureStress) + 'Pa ')]
[io.out('at' + string(failureStrain) + '% strain')]

;-------------------------------------------------------------------
;test 2: perform a triaxial test with isotropic confining stress 1e4
model restore 'beforeApplication'
[rampUp(0,-1e4,-1e3)]
model save 'to1e4'
wall attribute velocity-z [-platenVel] range name 'platenTop'
wall attribute velocity-z [platenVel] range name 'platenBottom'
ball attribute displacement (0,0,0)
structure node initialize displacement (0,0,0)
fish history stress
fish history strain
model solve fish-halt halt
model save 'triaxial1e4'
[io.out(string(failureStress) + 'Pa ')]
[io.out('at' + string(failureStrain) + '% strain')]

;-------------------------------------------------------------------
;test 3: perform a triaxial test with isotropic confining stress 5e4
model restore 'to1e4'
[rampUp(-1e4,-5e4,-1e3)]
model save 'to5e4'
wall attribute velocity-z [-platenVel] range name 'platenTop'
wall attribute velocity-z [platenVel] range name 'platenBottom'
ball attribute displacement (0,0,0)
structure node initialize displacement (0,0,0)
fish history stress
fish history strain
model solve fish-halt halt
model save 'triaxial5e4'
[io.out(string(failureStress) + 'Pa ')]
[io.out('at' + string(failureStrain) + '% strain')]

;-------------------------------------------------------------------
;test 4: perform a triaxial test with isotropic confining stress 1e5
model restore 'to5e4'
[rampUp(-5e4,-1e5,-1e3)]
model save 'to1e5'
wall attribute velocity-z [-platenVel] range name 'platenTop'
wall attribute velocity-z [platenVel] range name 'platenBottom'
ball attribute displacement (0,0,0)
structure node initialize displacement (0,0,0)
fish history stress
fish history strain
model solve fish-halt halt
model save 'triaxial1e5'
[io.out(string(failureStress) + 'Pa ')]
[io.out('at' + string(failureStrain) + '% strain')]