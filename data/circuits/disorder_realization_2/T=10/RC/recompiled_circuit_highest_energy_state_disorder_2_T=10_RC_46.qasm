OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.18345565) q[0];
sx q[0];
rz(2.3975211) q[0];
sx q[0];
rz(10.476713) q[0];
rz(0.19368859) q[1];
sx q[1];
rz(-0.63313484) q[1];
sx q[1];
rz(-2.5685891) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.409789) q[0];
sx q[0];
rz(-0.62792009) q[0];
sx q[0];
rz(-1.464792) q[0];
x q[1];
rz(-0.51197211) q[2];
sx q[2];
rz(-0.7535156) q[2];
sx q[2];
rz(-1.3775502) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1243368) q[1];
sx q[1];
rz(-1.6410429) q[1];
sx q[1];
rz(-0.4976561) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6735733) q[3];
sx q[3];
rz(-2.2647018) q[3];
sx q[3];
rz(-0.68460195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.28213349) q[2];
sx q[2];
rz(-0.30974516) q[2];
sx q[2];
rz(2.0852883) q[2];
rz(-3.1193962) q[3];
sx q[3];
rz(-2.3789417) q[3];
sx q[3];
rz(1.7215884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2317155) q[0];
sx q[0];
rz(-0.48119369) q[0];
sx q[0];
rz(0.43014446) q[0];
rz(3.0138956) q[1];
sx q[1];
rz(-1.1558497) q[1];
sx q[1];
rz(1.7040303) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4081892) q[0];
sx q[0];
rz(-1.3131307) q[0];
sx q[0];
rz(-2.8297421) q[0];
x q[1];
rz(-1.5135852) q[2];
sx q[2];
rz(-1.5199516) q[2];
sx q[2];
rz(1.0982996) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.592652) q[1];
sx q[1];
rz(-0.90528622) q[1];
sx q[1];
rz(1.9926461) q[1];
rz(-pi) q[2];
rz(-2.6019179) q[3];
sx q[3];
rz(-1.7575348) q[3];
sx q[3];
rz(0.58247551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0251856) q[2];
sx q[2];
rz(-1.4613287) q[2];
sx q[2];
rz(-2.6961668) q[2];
rz(-0.40924254) q[3];
sx q[3];
rz(-1.0990812) q[3];
sx q[3];
rz(-2.2621431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5889848) q[0];
sx q[0];
rz(-3.0486076) q[0];
sx q[0];
rz(-2.4267922) q[0];
rz(2.0843166) q[1];
sx q[1];
rz(-0.42066586) q[1];
sx q[1];
rz(0.32726273) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6598073) q[0];
sx q[0];
rz(-0.57087539) q[0];
sx q[0];
rz(1.8828189) q[0];
x q[1];
rz(0.94560854) q[2];
sx q[2];
rz(-0.66788061) q[2];
sx q[2];
rz(-2.5651074) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.79993805) q[1];
sx q[1];
rz(-2.9886878) q[1];
sx q[1];
rz(-1.7064817) q[1];
rz(-pi) q[2];
rz(0.9483665) q[3];
sx q[3];
rz(-2.7223058) q[3];
sx q[3];
rz(-0.47800999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1383692) q[2];
sx q[2];
rz(-0.97347632) q[2];
sx q[2];
rz(-1.6376015) q[2];
rz(-1.9231632) q[3];
sx q[3];
rz(-0.4156433) q[3];
sx q[3];
rz(-0.98201069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0729436) q[0];
sx q[0];
rz(-1.8953841) q[0];
sx q[0];
rz(-0.15596998) q[0];
rz(0.34863696) q[1];
sx q[1];
rz(-2.5376384) q[1];
sx q[1];
rz(-1.2329996) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5208682) q[0];
sx q[0];
rz(-1.3107131) q[0];
sx q[0];
rz(3.0164032) q[0];
x q[1];
rz(-0.92278752) q[2];
sx q[2];
rz(-1.7402667) q[2];
sx q[2];
rz(3.0344935) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.11248842) q[1];
sx q[1];
rz(-1.8524516) q[1];
sx q[1];
rz(-2.2206578) q[1];
rz(-pi) q[2];
rz(-2.8901222) q[3];
sx q[3];
rz(-1.5850388) q[3];
sx q[3];
rz(2.3492658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4431346) q[2];
sx q[2];
rz(-0.036018697) q[2];
sx q[2];
rz(-1.6301463) q[2];
rz(1.1394507) q[3];
sx q[3];
rz(-1.3316493) q[3];
sx q[3];
rz(0.99036923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.776942) q[0];
sx q[0];
rz(-2.7237837) q[0];
sx q[0];
rz(-1.9885709) q[0];
rz(0.54368377) q[1];
sx q[1];
rz(-1.2666356) q[1];
sx q[1];
rz(2.8797454) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7196822) q[0];
sx q[0];
rz(-1.6043264) q[0];
sx q[0];
rz(1.4872585) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0864429) q[2];
sx q[2];
rz(-2.5005955) q[2];
sx q[2];
rz(-2.6015559) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3186372) q[1];
sx q[1];
rz(-1.3483817) q[1];
sx q[1];
rz(-0.59808029) q[1];
x q[2];
rz(-2.35942) q[3];
sx q[3];
rz(-1.7165136) q[3];
sx q[3];
rz(-2.4268933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.61949817) q[2];
sx q[2];
rz(-2.5358989) q[2];
sx q[2];
rz(1.7363133) q[2];
rz(-2.3960579) q[3];
sx q[3];
rz(-1.8757952) q[3];
sx q[3];
rz(-2.1982927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0010506823) q[0];
sx q[0];
rz(-3.0241835) q[0];
sx q[0];
rz(-2.7897799) q[0];
rz(2.70128) q[1];
sx q[1];
rz(-1.3654717) q[1];
sx q[1];
rz(-0.17131677) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7650314) q[0];
sx q[0];
rz(-0.91463551) q[0];
sx q[0];
rz(-1.2622467) q[0];
x q[1];
rz(-1.2662884) q[2];
sx q[2];
rz(-0.84104462) q[2];
sx q[2];
rz(0.27308057) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.379516) q[1];
sx q[1];
rz(-0.87240309) q[1];
sx q[1];
rz(-1.5435436) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.76254179) q[3];
sx q[3];
rz(-2.2665215) q[3];
sx q[3];
rz(0.97555893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0773086) q[2];
sx q[2];
rz(-1.7694387) q[2];
sx q[2];
rz(2.1997931) q[2];
rz(1.9866379) q[3];
sx q[3];
rz(-0.20808163) q[3];
sx q[3];
rz(-1.8010767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44523859) q[0];
sx q[0];
rz(-2.0360763) q[0];
sx q[0];
rz(0.0048986991) q[0];
rz(-2.694963) q[1];
sx q[1];
rz(-2.547867) q[1];
sx q[1];
rz(2.4536224) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94243455) q[0];
sx q[0];
rz(-2.4068038) q[0];
sx q[0];
rz(-0.25811974) q[0];
x q[1];
rz(-0.46779386) q[2];
sx q[2];
rz(-2.4681598) q[2];
sx q[2];
rz(-1.3040257) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6354949) q[1];
sx q[1];
rz(-2.6539972) q[1];
sx q[1];
rz(-2.0388842) q[1];
rz(-0.68617188) q[3];
sx q[3];
rz(-0.99643512) q[3];
sx q[3];
rz(2.2597093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.530431) q[2];
sx q[2];
rz(-2.7335584) q[2];
sx q[2];
rz(-0.92998695) q[2];
rz(-1.2215349) q[3];
sx q[3];
rz(-2.0285716) q[3];
sx q[3];
rz(-0.59337029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68194836) q[0];
sx q[0];
rz(-1.821803) q[0];
sx q[0];
rz(3.0514858) q[0];
rz(-1.7168761) q[1];
sx q[1];
rz(-2.5017891) q[1];
sx q[1];
rz(2.7065014) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6259436) q[0];
sx q[0];
rz(-1.0727779) q[0];
sx q[0];
rz(1.4048715) q[0];
rz(-2.9138759) q[2];
sx q[2];
rz(-0.68074742) q[2];
sx q[2];
rz(1.6344223) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1798517) q[1];
sx q[1];
rz(-1.4423749) q[1];
sx q[1];
rz(2.1895616) q[1];
rz(0.42436142) q[3];
sx q[3];
rz(-1.5715944) q[3];
sx q[3];
rz(-0.80959807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4789751) q[2];
sx q[2];
rz(-2.0765897) q[2];
sx q[2];
rz(0.62620658) q[2];
rz(-0.19691697) q[3];
sx q[3];
rz(-2.1508689) q[3];
sx q[3];
rz(-1.7530493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5597252) q[0];
sx q[0];
rz(-1.4511061) q[0];
sx q[0];
rz(3.0873121) q[0];
rz(0.42298969) q[1];
sx q[1];
rz(-1.3754247) q[1];
sx q[1];
rz(-1.3351701) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88807073) q[0];
sx q[0];
rz(-1.0253064) q[0];
sx q[0];
rz(2.2175199) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1303875) q[2];
sx q[2];
rz(-1.2647243) q[2];
sx q[2];
rz(2.0843506) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9250592) q[1];
sx q[1];
rz(-2.0833942) q[1];
sx q[1];
rz(-0.61116265) q[1];
rz(-pi) q[2];
rz(-1.7087206) q[3];
sx q[3];
rz(-1.5028126) q[3];
sx q[3];
rz(1.2757511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6173031) q[2];
sx q[2];
rz(-1.4681939) q[2];
sx q[2];
rz(-2.7900901) q[2];
rz(1.7678123) q[3];
sx q[3];
rz(-2.6280554) q[3];
sx q[3];
rz(-2.8721749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7518625) q[0];
sx q[0];
rz(-2.8808012) q[0];
sx q[0];
rz(0.75916284) q[0];
rz(-2.0579386) q[1];
sx q[1];
rz(-1.5307129) q[1];
sx q[1];
rz(-1.0677451) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6747492) q[0];
sx q[0];
rz(-2.3567794) q[0];
sx q[0];
rz(2.270257) q[0];
rz(-0.2944417) q[2];
sx q[2];
rz(-1.8368372) q[2];
sx q[2];
rz(-2.5584115) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8651004) q[1];
sx q[1];
rz(-0.70521627) q[1];
sx q[1];
rz(2.221285) q[1];
rz(-pi) q[2];
rz(1.9633406) q[3];
sx q[3];
rz(-2.4053229) q[3];
sx q[3];
rz(-0.53680778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4296253) q[2];
sx q[2];
rz(-1.9316614) q[2];
sx q[2];
rz(-2.9313226) q[2];
rz(2.353904) q[3];
sx q[3];
rz(-1.7588047) q[3];
sx q[3];
rz(-0.53019607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44337153) q[0];
sx q[0];
rz(-2.5826695) q[0];
sx q[0];
rz(-1.2179751) q[0];
rz(-1.5086077) q[1];
sx q[1];
rz(-1.8764381) q[1];
sx q[1];
rz(1.1988342) q[1];
rz(-0.88268139) q[2];
sx q[2];
rz(-1.8919049) q[2];
sx q[2];
rz(1.6406825) q[2];
rz(2.0785594) q[3];
sx q[3];
rz(-2.8294143) q[3];
sx q[3];
rz(-2.1635273) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
