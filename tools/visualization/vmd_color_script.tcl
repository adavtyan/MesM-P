set framecount [molinfo top get numframes]
set sel [atomselect top all]

for {set frame_no 0} {$frame_no<$framecount} {incr frame_no} {
       $sel frame $frame_no
       $sel set user  [$sel get vx]
}
for {set frame_no 0} {$frame_no<$framecount} {incr frame_no} {
       $sel frame $frame_no
       $sel set user2  [$sel get vy]
}
for {set frame_no 0} {$frame_no<$framecount} {incr frame_no} {
       $sel frame $frame_no
       $sel set user3  [$sel get vz]
}
$sel delete

