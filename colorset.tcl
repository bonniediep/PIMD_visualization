## DESCRIPTION:
## This TCL script is a simple method to visualize path integrals through VMD from both a PDB and NC file, while producing values for centre of mass. The visualization of each molecule can be changed in terms of its colour, representation type, selection method (e.g. residue type, index, element, xyz coordinates, etc.), and material. The calculated centre of masses replaces the xyz coordinates of each first bead of each atom, which is connected visually. Colouring of atoms are based on standard CPK rules.

HI

## PROCEDURES:
## colorset

## EXAMPLE USAGE (WITHOUT VMD MAIN WINDOW & TK CONSOLE ACCESS):
## vmd -f nameofpdbfile.pdb -f nameofncfile.nc -startup colorset.tcl

## EXAMPLE USAGE (WITH VMD MAIN WINDOW & TK CONSOLE ACCESS):
## vmd -f nameofpdbfile.pdb -f nameofncfile.nc
## Then open Tk Console from Extensions dropdown menu and type in
## source colorset.tcl

proc colorset {} { 
    # Colours: 0=blue; 1=red; 2=gray; 3=orange; 4=yellow; 5=tan; 6=silver; 7=green; 8=white; 9=pink; 10=cyan; 11=purple; 12=lime; 13=mauve; 14=ochre; 15=iceblue; 16=black; 17=yellow2; 18=yellow3; 19=green2; 20=green3; 21=cyan2; 22=cyan3; 23=blue2; 24=blue3; 25=violet; 26=violet2; 27=magenta; 28=magenta2; 29=red2; 30=red3; 31=orange2; 32=orange3 

    # Representations: Lines, Bonds, DynamicBonds, HBonds, Points, VDW, CPK, Licorice, Polyhedra, Trace, Tube, Ribbons, NewRibbons, Cartoon, NewCartoon, PaperChain, Twister, QuickSurf, MSMS, Surf, VolumeSlice, IsoSurface, FieldLines, Orbital, Beads, Dotted, Solvent

    # Selection Methods: name, type, backbonetype, residuetype, index, serial, atomicnumber, element, residue, resname, altloc, resid, insertion, chain, segname, segid, fragment, pfrag, nfrag, numbonds, structure, structure, pucker, user, user2, user3, user4, x, y, z, vx, vy, vz, ufx, ufy, ufz, phi, psi, radius, mass, charge, beta, occupancy, volindex0, volindex1, volindex2, volindex3, volindex4, volindex5, volindex6, volindex7, vol0, vol1, vol2, vol3, vol4, vol5, vol6, vol7, interpvol0, interpvol1, interpvol2, interpvol3, interpvol4, interpvol5, interpvol6, interpvol7 

    # Colouring Methods: Name, Type, Element, ResName, ResType, ResID, Chain, SegName, Conformation, Molecule, Secondary Structure, ColorID, Beta, Occupancy, Mass, Charge, Position, Trajectory, Fragment, Index, Backbone, Throb, Volume
    
    # Material: Opaque, Transparent, BrushedMetal, Diffuse, Ghost, Glass1, Glass2, Glass3, Glossy, HardPlastic, MetallicPastel, Steel, Translucent, Edgy, EdgyShiny, EdgyGlass, GoodShell, AOShiny, AOChalky, AOEdgy, BlownGlass, GlassBubble, RTChrome

    # Can change background display colour:
    color Display Background 2

    ## 0, 1, 2 placed after modselect are the representation numbers                       
    ## "top" is the id of the top molecule or molecule number, can often refer to molID instead 

    # First molecule (pdb file) will have molID 0
    # Second molecule (nc file) will have molID 1
    
    # Note: Ensure that pdb file is in alphabetical order

    # The "top" molecule is the one opened last                                                   
    # To assign a different molecule as the "top" one:                                            
    # mol top molID  

    # nbeads = numatoms in nc / numatoms in pdb
    set numpdb [molinfo 0 get numatoms]
    set numnc [molinfo 1 get numatoms]
    
    # Hide pdb file: (not necessary to show pdb file) (Format: mol off molID)
    mol off 0
    
    # Get the number of beads per atom
    set nbeads [expr {$numnc/$numpdb}]
    
    # Set representation for pdb molecule:
    # To adjust size of chosen representation, default settings chosen if no values given after name of representation; otherwise add values for scaling, thickness, etc. (Format: mol representation NameofRepresentation OptionalValues)
    # Change coloring method (Format: mol color ColourMethod)
    # Apply modifications (Format: mol modrep RepresentationNumber MoleculeNumber)
    mol representation Lines
    mol color Name
    mol modrep 0 0
    
    # Select for each element and find quantity of each
    # Format:
    # set VariableName [atomselect molID SelectionText]
    # set VariableName Value
    set carbon [atomselect 0 "element C"]
    set numC [$carbon num]
    
    set chlorine [atomselect 0 "element Cl"]
    set numCl [$chlorine num]
    
    set fluorine [atomselect 0 "element F"]
    set numF [$fluorine num]

    set hydrogen [atomselect 0 "element H"]
    set numH [$hydrogen num]
    
    set potassium [atomselect 0 "element K"]
    set numK [$potassium num]

    set nitrogen [atomselect 0 "element N"]
    set numN [$nitrogen num]

    set oxygen [atomselect 0 "element O"]
    set numO [$oxygen num]
    
    set phosphorus [atomselect 0 "element P"]
    set numP [$phosphorus num]
    
    set sulfur [atomselect 0 "element S"]
    set numS [$sulfur num]

    # Change colouring method and representation of nc molecule (when using ColorID as ColouringMethod, add value of colour afterwards to adjust)
    # Colour the beads according to standard CPK rules
    # Let n be a counter for RepresentationNumber of nc file
    # Let prev be the current number of beads so far
    # Add new representation for each new element

    mol representation Points 10
    set n 0
        
    if {$numC != 0} {
	mol color ColorID 10
	mol modrep 0 1
 	mol modselect $n top serial<=($numC*$nbeads)
 	set prev $numC
    } else {}

    if {$numCl != 0} {
	set n [expr {$n + 1}]
	mol color ColorID 19
	mol addrep 1
	mol modselect $n top (serial>($prev*$nbeads) and serial<=(($prev + $numCl)*$nbeads))
        set prev [expr {$prev + $numCl}]
    } else {}

    if {$numF != 0} {
	set n [expr {$n + 1}]
	mol color ColorID 7
	mol addrep 1
	mol modselect $n top (serial>($prev*$nbeads) and serial<=(($prev + $numF)*$nbeads))
	set prev [expr {$prev + $numF}]
    } else {}

    if {$numH != 0} {
	set n [expr {$n + 1}]
	mol color ColorID 8
	mol addrep 1
	mol modselect $n top (serial>($prev*$nbeads) and serial<=(($prev + $numH)*$nbeads))
	set prev [expr {$prev + $numH}]
    } else {}

    if {$numK != 0} {
	set n [expr {$n + 1}]
	mol color ColorID 25
	mol addrep 1
	mol modselect $n top (serial>($prev*$nbeads) and serial<=(($prev + $numK)*$nbeads))
	set prev [expr {$prev + $numK}]
    } else {}

    if {$numN != 0} {
	set n [expr {$n + 1}]
	mol color ColorID 0
	mol addrep 1
	mol modselect $n top (serial>($prev*$nbeads) and serial<=(($prev + $numN)*$nbeads))
	set prev [expr {$prev + $numN}]
    } else {}

    if {$numO != 0} {
	set n [expr {$n + 1}]
	mol color ColorID 1
	mol addrep 1
	mol modselect $n top (serial>($prev*$nbeads) and serial<=(($prev + $numO)*$nbeads))
	set prev [expr {$prev + $numO}]
    } else {}

    if {$numP != 0} {
	set n [expr {$n + 1}]
	mol color ColorID 3
	mol addrep 1
	mol modselect $n top (serial>($prev*$nbeads) and serial<=(($prev + $numP)*$nbeads))
	set prev [expr {$prev + $numP}]
    } else {}

    if {$numS != 0} {
	set n [expr {$n + 1}]
	mol color ColorID 4
	mol addrep 1
	mol modselect $n top (serial>($prev*$nbeads) and serial<=(($prev + $numS)*$nbeads))
    } else {}

    # List of all atoms
    set alist {}
    for {set i 1} {$i<=($numnc/$nbeads)} {incr i} {
        lappend alist $i}

    # Calculating centre of mass 
    set coordlist ""
    foreach atom $alist {
	set newcoord {}
	set xsum 0
	set ysum 0
	set zsum 0
       	set beads [atomselect top "serial [expr (($atom-1)*$nbeads+1)] to [expr ($atom*$nbeads)]"]
	set xcoord [$beads get {x}]
	set ycoord [$beads get {y}]
	set zcoord [$beads get {z}]
	foreach x $xcoord {set xsum [expr ($xsum + $x)]}
	foreach y $ycoord {set ysum [expr ($ysum + $y)]}
	foreach z $zcoord {set zsum [expr ($zsum + $z)]}
	set comx [expr ($xsum/$nbeads)]
	set comy [expr ($ysum/$nbeads)]
	set comz [expr ($zsum/$nbeads)]
	lappend newcoord $comx
	lappend newcoord $comy
	lappend newcoord $comz
	lappend coordlist $newcoord}
    
    # Selecting for every first bead of each atom, then change its coordinates to the calculated centre of masses
    set counter 0
    foreach coord $coordlist {
	set counter [expr {$counter + 1}]
	set new [atomselect top "index [expr (($counter-1)*$nbeads+1)]"]
	$new moveto $coord
	set s [$new get {x y z}]
	puts $s}
    
    # Making the selection for the representation
    set sel ""
    foreach a $alist {
        lappend sel "serial [expr (($a-1)*$nbeads+1)]" or}
    set sel [linsert $sel 0 "("]
    set sel [regsub -all {\{|\}} $sel ""]
    set sellist [lreplace $sel [expr [llength $sel]-1] [expr [llength $sel]-1]]
    set sellist [linsert $sellist [llength $sellist] ")"]
    
    # Set colouring and representation type in a new representation                                  
    mol color ColorID 12
    mol representation DynamicBonds 1.6 0.1 12.0
    mol material Transparent
    set n [expr {$n + 1}]
    mol addrep 1
    mol modselect $n top $sellist

    # Start animation
    animate forward}

# call function:
colorset   

# Enter animate pause to pause animation
