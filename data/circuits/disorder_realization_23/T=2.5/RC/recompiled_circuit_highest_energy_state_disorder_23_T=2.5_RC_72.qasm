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
rz(-1.3462525) q[0];
sx q[0];
rz(-1.4486382) q[0];
sx q[0];
rz(-1.9899272) q[0];
rz(1.0898074) q[1];
sx q[1];
rz(1.6648219) q[1];
sx q[1];
rz(9.2738168) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.359084) q[0];
sx q[0];
rz(-1.623296) q[0];
sx q[0];
rz(1.8655325) q[0];
rz(-pi) q[1];
x q[1];
rz(0.26246983) q[2];
sx q[2];
rz(-0.32400692) q[2];
sx q[2];
rz(0.25036795) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0651567) q[1];
sx q[1];
rz(-1.5744835) q[1];
sx q[1];
rz(1.5559959) q[1];
rz(-pi) q[2];
rz(-1.6116222) q[3];
sx q[3];
rz(-1.4741952) q[3];
sx q[3];
rz(-0.20303328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0207409) q[2];
sx q[2];
rz(-3.1181702) q[2];
sx q[2];
rz(-0.49129301) q[2];
rz(-1.8190207) q[3];
sx q[3];
rz(-1.5343851) q[3];
sx q[3];
rz(-1.3278495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78111929) q[0];
sx q[0];
rz(-0.55447117) q[0];
sx q[0];
rz(-0.69960272) q[0];
rz(-1.5505002) q[1];
sx q[1];
rz(-0.49483776) q[1];
sx q[1];
rz(2.9717305) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79029523) q[0];
sx q[0];
rz(-1.5322432) q[0];
sx q[0];
rz(-0.077705381) q[0];
rz(2.8367859) q[2];
sx q[2];
rz(-1.4423324) q[2];
sx q[2];
rz(0.94513946) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3292698) q[1];
sx q[1];
rz(-2.053875) q[1];
sx q[1];
rz(2.7365351) q[1];
rz(-pi) q[2];
x q[2];
rz(0.079143123) q[3];
sx q[3];
rz(-2.0266266) q[3];
sx q[3];
rz(1.7439708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.26213172) q[2];
sx q[2];
rz(-2.7792271) q[2];
sx q[2];
rz(1.7246838) q[2];
rz(-2.0587685) q[3];
sx q[3];
rz(-0.88734752) q[3];
sx q[3];
rz(2.3143903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5603492) q[0];
sx q[0];
rz(-2.2214948) q[0];
sx q[0];
rz(-1.5592519) q[0];
rz(0.70445591) q[1];
sx q[1];
rz(-1.9943941) q[1];
sx q[1];
rz(-0.53680435) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5671331) q[0];
sx q[0];
rz(-1.6170814) q[0];
sx q[0];
rz(0.68741902) q[0];
x q[1];
rz(1.6073999) q[2];
sx q[2];
rz(-1.6609909) q[2];
sx q[2];
rz(-2.6548403) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4406261) q[1];
sx q[1];
rz(-1.4203209) q[1];
sx q[1];
rz(-1.6789084) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2800673) q[3];
sx q[3];
rz(-0.71746263) q[3];
sx q[3];
rz(2.2253401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.392936) q[2];
sx q[2];
rz(-1.8042678) q[2];
sx q[2];
rz(-2.362747) q[2];
rz(-0.084176453) q[3];
sx q[3];
rz(-0.86936969) q[3];
sx q[3];
rz(1.4028153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8970784) q[0];
sx q[0];
rz(-2.2205413) q[0];
sx q[0];
rz(0.47065863) q[0];
rz(-1.4587519) q[1];
sx q[1];
rz(-0.0048480821) q[1];
sx q[1];
rz(-2.2768314) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7047459) q[0];
sx q[0];
rz(-0.517874) q[0];
sx q[0];
rz(3.0240866) q[0];
rz(-2.607279) q[2];
sx q[2];
rz(-0.25548894) q[2];
sx q[2];
rz(-1.3380255) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4196849) q[1];
sx q[1];
rz(-0.083217155) q[1];
sx q[1];
rz(0.91099847) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8806609) q[3];
sx q[3];
rz(-2.0892077) q[3];
sx q[3];
rz(0.37004091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.027815759) q[2];
sx q[2];
rz(-2.2769589) q[2];
sx q[2];
rz(1.9846385) q[2];
rz(0.35465869) q[3];
sx q[3];
rz(-1.2361453) q[3];
sx q[3];
rz(-0.11053301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.329634) q[0];
sx q[0];
rz(-2.2060153) q[0];
sx q[0];
rz(-2.6079566) q[0];
rz(-1.2136906) q[1];
sx q[1];
rz(-0.023093725) q[1];
sx q[1];
rz(1.1812706) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5438687) q[0];
sx q[0];
rz(-1.4078316) q[0];
sx q[0];
rz(-1.5242432) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6152739) q[2];
sx q[2];
rz(-0.72229702) q[2];
sx q[2];
rz(2.1762951) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5983728) q[1];
sx q[1];
rz(-0.66987619) q[1];
sx q[1];
rz(1.3175439) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3245246) q[3];
sx q[3];
rz(-0.51760537) q[3];
sx q[3];
rz(2.8235265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.033325) q[2];
sx q[2];
rz(-0.51563087) q[2];
sx q[2];
rz(-0.92833129) q[2];
rz(-2.8632274) q[3];
sx q[3];
rz(-1.3185578) q[3];
sx q[3];
rz(3.0729455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.5334897) q[0];
sx q[0];
rz(-0.52284306) q[0];
sx q[0];
rz(-2.9485517) q[0];
rz(-2.4642384) q[1];
sx q[1];
rz(-3.1132128) q[1];
sx q[1];
rz(1.631564) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2264435) q[0];
sx q[0];
rz(-2.9993176) q[0];
sx q[0];
rz(2.9998542) q[0];
rz(3.0966412) q[2];
sx q[2];
rz(-2.0325902) q[2];
sx q[2];
rz(3.1147946) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.88107706) q[1];
sx q[1];
rz(-2.279638) q[1];
sx q[1];
rz(-2.8167732) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0738689) q[3];
sx q[3];
rz(-0.9091223) q[3];
sx q[3];
rz(-0.53406427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3626927) q[2];
sx q[2];
rz(-0.71234667) q[2];
sx q[2];
rz(2.1544797) q[2];
rz(-2.0569862) q[3];
sx q[3];
rz(-0.54607138) q[3];
sx q[3];
rz(0.62091056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5793594) q[0];
sx q[0];
rz(-2.1491304) q[0];
sx q[0];
rz(1.5935422) q[0];
rz(0.69001895) q[1];
sx q[1];
rz(-0.023612173) q[1];
sx q[1];
rz(1.4366368) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48955917) q[0];
sx q[0];
rz(-1.0512182) q[0];
sx q[0];
rz(1.621031) q[0];
x q[1];
rz(2.4568672) q[2];
sx q[2];
rz(-1.5878356) q[2];
sx q[2];
rz(-0.34340252) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1664983) q[1];
sx q[1];
rz(-2.1121032) q[1];
sx q[1];
rz(-0.06947738) q[1];
rz(-pi) q[2];
rz(-1.2728746) q[3];
sx q[3];
rz(-1.1323338) q[3];
sx q[3];
rz(0.39194731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7684266) q[2];
sx q[2];
rz(-1.2370141) q[2];
sx q[2];
rz(0.58488673) q[2];
rz(-1.5960426) q[3];
sx q[3];
rz(-1.4228179) q[3];
sx q[3];
rz(2.2271631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5717413) q[0];
sx q[0];
rz(-0.91747704) q[0];
sx q[0];
rz(1.566994) q[0];
rz(0.51708108) q[1];
sx q[1];
rz(-2.2061429) q[1];
sx q[1];
rz(1.8749974) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1388113) q[0];
sx q[0];
rz(-1.5673593) q[0];
sx q[0];
rz(-3.1365738) q[0];
rz(-pi) q[1];
rz(-2.3364303) q[2];
sx q[2];
rz(-1.9765719) q[2];
sx q[2];
rz(-2.9540887) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.75830981) q[1];
sx q[1];
rz(-0.66557887) q[1];
sx q[1];
rz(2.1583827) q[1];
x q[2];
rz(1.4483767) q[3];
sx q[3];
rz(-1.6922235) q[3];
sx q[3];
rz(1.1894912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7685585) q[2];
sx q[2];
rz(-0.80058241) q[2];
sx q[2];
rz(-1.2726834) q[2];
rz(-0.4396762) q[3];
sx q[3];
rz(-1.1594073) q[3];
sx q[3];
rz(0.064597733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8398297) q[0];
sx q[0];
rz(-0.01627144) q[0];
sx q[0];
rz(2.8453258) q[0];
rz(-0.06427327) q[1];
sx q[1];
rz(-1.6769033) q[1];
sx q[1];
rz(-1.5110663) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8730174) q[0];
sx q[0];
rz(-2.7843256) q[0];
sx q[0];
rz(-1.1910458) q[0];
rz(-pi) q[1];
rz(-1.3675646) q[2];
sx q[2];
rz(-2.2335238) q[2];
sx q[2];
rz(-1.4858147) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.293773) q[1];
sx q[1];
rz(-0.34983954) q[1];
sx q[1];
rz(1.185525) q[1];
rz(-pi) q[2];
rz(-0.011476403) q[3];
sx q[3];
rz(-1.4496798) q[3];
sx q[3];
rz(0.83601421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5955547) q[2];
sx q[2];
rz(-1.7067355) q[2];
sx q[2];
rz(2.0175374) q[2];
rz(-1.3042287) q[3];
sx q[3];
rz(-0.29574695) q[3];
sx q[3];
rz(-0.058163253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8703363) q[0];
sx q[0];
rz(-2.3503292) q[0];
sx q[0];
rz(-1.9747718) q[0];
rz(-1.600949) q[1];
sx q[1];
rz(-2.8438042) q[1];
sx q[1];
rz(1.8150394) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7739511) q[0];
sx q[0];
rz(-2.2164167) q[0];
sx q[0];
rz(1.2626179) q[0];
x q[1];
rz(0.72188702) q[2];
sx q[2];
rz(-2.6195002) q[2];
sx q[2];
rz(2.4185857) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5656783) q[1];
sx q[1];
rz(-0.88972461) q[1];
sx q[1];
rz(1.8450062) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3438026) q[3];
sx q[3];
rz(-0.41068893) q[3];
sx q[3];
rz(2.3410377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.88663179) q[2];
sx q[2];
rz(-0.16356629) q[2];
sx q[2];
rz(1.6939885) q[2];
rz(-0.12668954) q[3];
sx q[3];
rz(-1.4796175) q[3];
sx q[3];
rz(-2.0255585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3348677) q[0];
sx q[0];
rz(-1.8341478) q[0];
sx q[0];
rz(-1.3679416) q[0];
rz(1.548829) q[1];
sx q[1];
rz(-0.82850414) q[1];
sx q[1];
rz(-2.9864476) q[1];
rz(2.7663076) q[2];
sx q[2];
rz(-2.6659758) q[2];
sx q[2];
rz(-2.804826) q[2];
rz(0.98649377) q[3];
sx q[3];
rz(-1.428008) q[3];
sx q[3];
rz(-1.2724439) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
