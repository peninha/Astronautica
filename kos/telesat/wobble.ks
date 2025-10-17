parameter sample_time is 30.0.  // sampling time in seconds
parameter sample_interval is 0.1. // time between samples

function calc_wobble {
    parameter sample_time is 30.0.
    parameter sample_interval is 0.1.

    local samples is round(sample_time / sample_interval).
    local dt is sample_interval.

    // Lista para armazenar as amostras de angularvel
    local facing_vecs is list().
    print "Getting samples... wait " + dt*samples + " seconds.".

    local avg is V(0,0,0).
    local i is 0.
    until i = samples {
        WAIT dt.
        facing_vecs:add(ship:facing:vector).
        set avg to avg + facing_vecs[i].
        set i to i + 1.
    }
    set avg to avg / samples.
    print "Average facing vector: " + avg.

    // Calculate wobbleness (average angle between each facing vector and the average facing vector)
    set total_angle to 0.
    for sample IN facing_vecs {
        set angle_diff to vang(avg, sample).
        set total_angle to total_angle + angle_diff.
    }
    set average_angle to total_angle / samples.
    return average_angle.
}

print "Checking stability...               ".
local wobble is calc_wobble(sample_time, sample_interval).
print "wobble: " + wobble + " deg".
unlock all.