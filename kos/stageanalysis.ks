// https://www.reddit.com/r/Kos/comments/3fq4is/a_script_to_get_information_on_stages_in_advance/

lock throttle to 0.
wait 0.05.
 
// the ultimate product of the program
// will contain lists of data about each stage
set stages to list().
 
list engines in engineList.
 
local g is 9.80665.
 
// have to activate engines to get their thrusts (for now)
local activeEngines is list().
for eng in engineList
	if eng:ignition = false
		eng:activate.
	else 
		activeEngines:add(eng).
 
// tagging all decouplers as decoupler by looking at part modules
for part in ship:parts
	for module in part:modules
		if part:getModule(module):name = "ModuleDecouple" 
		or part:getModule(module):name = "ModuleAnchoredDecoupler"
			set part:tag to "decoupler".
 
// sections are groups of parts between decouplers
// the roots of sections are the ship root and all decouplers
local sectionRoots is list().
 
sectionRoots:add(ship:rootPart).
 
for decoupler in ship:partsTagged("decoupler")
	sectionRoots:add(decoupler).
 
// lists of (root part, mass, fuelmass, engines, and fuelflow)
local sections is list().
 
// creates a section from the root of each section
for sectionRoot in sectionRoots{
 
	local sectionMass is 0.
	local sectionFuelMass is 0.
	local sectionEngineList is list().
	local fuelFlow is 0.
 
	local sectionParts is list().
	sectionParts:add(sectionRoot).
 
	// add all children down part tree from the section root to
	// list of section parts unless they are a decoupler or launch clamp
	local i is 0.
	until i = sectionParts:length{
		if sectionParts[i]:children:empty = false
			for child in sectionParts[i]:children
				if child:tag = "decoupler" = false and child:name = "LaunchClamp1" = false
					sectionParts:add(child).
		set i to i + 1.
	}
 
	for part in sectionParts{
 
		set sectionMass to sectionMass + part:mass.
 
		// avoiding adding rcs fuel to fuelmass
		local rcsFlag is false.
 
		if part:resources:empty = false{
			for resource in part:resources
				if resource:name = "monopropellant"
					set rcsFlag to true.
 
			if rcsFlag = false
				set sectionFuelMass to sectionFuelMass + part:mass - part:drymass.
		}
 
		if engineList:contains(part)
			sectionEngineList:add(part).
	}
 
	local section is list(sectionRoot,sectionMass,sectionFuelMass,sectionEngineList,fuelFlow).
	sections:add(section).
}
 
local firstStageNum is 0.
for eng in engineList
	if eng:stage > firstStageNum
		set firstStageNum to eng:stage.
 
// counting down from first (highest number) stage
// to stage 0, creating stage data and
// updating mass and fuelmass of the sections as it goes
// essentially "simulating" staging
local i is firstStageNum.
until i = -1 {
 
	// the four things really needed
	local stageMass is 0.
	local stageThrust is 0.
	local stageFuelFlow is 0.
	local stageBurnTime is 987654321. // starts high cause need to find lowest
 
	// other stuff it may as well calculate
	local stageMinA is 0.
	local stageMaxA is 0.
	local stageISP is 0.
	local stageDeltaV is 0.
 
	local curStage is list().
 
	// if the section decoupler activates on this stage remove that section
	// except the first section root (not a decoupler)
	local k is sections:length - 1.
	until k = 0{
 
		if sections[k][0]:stage = i
			sections:remove(k).
				set k to k - 1.
		}
 
	// generating the stage mass, thrust, fuelflow, and burntime 
	// from the sections that make up the stage
	for section in sections{
 
		local sectionMass is section[1].
		local sectionFuelMass is section[2].
		// resetting fuelflow
		set section[4] to 0.
		local sectionBurnTime is 0.
 
		set stageMass to stageMass + sectionMass.
 
		if section[3]:empty = false{
			for engine in section[3]
				if engine:stage  >= i{
					set stageThrust to stageThrust + engine:maxthrustat(0).
					set stageFuelFlow to stageFuelFlow + engine:maxthrustat(0)/engine:visp.
					set section[4] to section[4] + engine:maxthrustat(0)/engine:visp.
				}
		}
		// if it has fuelflow (active engines)
		if section[4] > 0{
			set sectionBurnTime to g * section[2] / section[4].
 
			// if the section will stage next stage
			// or this if this is the last stage
			if section[0]: stage = i - 1 or i = 0
				if sectionBurnTime < stageBurnTime
					set stageBurnTime to sectionBurnTime.
		}
	}
 
	// only possible if there are no active engines this stage (or god help you)
	if stageBurnTime = 987654321
		set stageBurnTime to 0.
 
	// calculating optional goodies
	if stageBurnTime > 0{
		local stageEndMass is stageMass - stageBurnTime * stageFuelFlow / g.
		set stageMinA to stageThrust / stageMass.
		set stageMaxA to stageThrust/ stageEndMass.
		set stageISP to stageThrust / stageFuelFlow.
		set stageDeltaV to stageISP * g * ln(stageMass / stageEndMass).
	}
 
	// take a deep breath
	set curStage to list(stageMass,stageISP,stageThrust,stageMinA,stageMaxA,stageDeltaV,stageBurnTime).
	stages:add(curStage).
 
	// reduce the mass and fuel mass of sections with active engines
	// according to the burn time of the stage
	for section in sections{
		set section[1] to section[1]- stageBurnTime * section[4] / g.
		set section[2] to section[2] - stageBurnTime * section[4] / g.
		}
 
	set i to i - 1.
}
 
// remove stages with no burn time
// comment out if you're curious, should look more like KERs' "show all stages" in VAB
{
	local i is stages:length - 1.
	until i = -1{
		if stages[i][6] = 0
		stages:remove(i).
		set i to i - 1.
		}
}
// shutting engines back down
for eng in engineList
	if activeEngines:contains(eng) = false
		eng:shutdown.