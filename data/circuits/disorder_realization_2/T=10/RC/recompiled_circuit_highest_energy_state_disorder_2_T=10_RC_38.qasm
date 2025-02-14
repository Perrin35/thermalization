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
rz(-0.7440716) q[0];
sx q[0];
rz(-1.0519354) q[0];
rz(-2.9479041) q[1];
sx q[1];
rz(-2.5084578) q[1];
sx q[1];
rz(2.5685891) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5405226) q[0];
sx q[0];
rz(-2.1946476) q[0];
sx q[0];
rz(-3.0649351) q[0];
x q[1];
rz(-2.4560375) q[2];
sx q[2];
rz(-1.2289882) q[2];
sx q[2];
rz(-0.58219257) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4595653) q[1];
sx q[1];
rz(-2.6394156) q[1];
sx q[1];
rz(2.9952496) q[1];
rz(0.12270452) q[3];
sx q[3];
rz(-0.70022455) q[3];
sx q[3];
rz(2.6168857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.28213349) q[2];
sx q[2];
rz(-0.30974516) q[2];
sx q[2];
rz(-1.0563043) q[2];
rz(-3.1193962) q[3];
sx q[3];
rz(-0.76265097) q[3];
sx q[3];
rz(-1.7215884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9098772) q[0];
sx q[0];
rz(-2.660399) q[0];
sx q[0];
rz(-2.7114482) q[0];
rz(-3.0138956) q[1];
sx q[1];
rz(-1.9857429) q[1];
sx q[1];
rz(1.7040303) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8970222) q[0];
sx q[0];
rz(-1.8720227) q[0];
sx q[0];
rz(-1.3006849) q[0];
rz(-0.050927863) q[2];
sx q[2];
rz(-1.5136592) q[2];
sx q[2];
rz(-0.47540755) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.592652) q[1];
sx q[1];
rz(-0.90528622) q[1];
sx q[1];
rz(1.1489465) q[1];
rz(-pi) q[2];
rz(-0.35234612) q[3];
sx q[3];
rz(-0.56803136) q[3];
sx q[3];
rz(-1.2888792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.11640707) q[2];
sx q[2];
rz(-1.4613287) q[2];
sx q[2];
rz(0.44542584) q[2];
rz(2.7323501) q[3];
sx q[3];
rz(-2.0425115) q[3];
sx q[3];
rz(-0.87944952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55260783) q[0];
sx q[0];
rz(-0.092985066) q[0];
sx q[0];
rz(2.4267922) q[0];
rz(-2.0843166) q[1];
sx q[1];
rz(-2.7209268) q[1];
sx q[1];
rz(0.32726273) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11571685) q[0];
sx q[0];
rz(-2.1109606) q[0];
sx q[0];
rz(2.9469423) q[0];
rz(-pi) q[1];
rz(0.43250044) q[2];
sx q[2];
rz(-2.0969146) q[2];
sx q[2];
rz(1.8217349) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3416546) q[1];
sx q[1];
rz(-0.1529049) q[1];
sx q[1];
rz(-1.4351109) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8873575) q[3];
sx q[3];
rz(-1.9079068) q[3];
sx q[3];
rz(1.1439307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0032234) q[2];
sx q[2];
rz(-2.1681163) q[2];
sx q[2];
rz(-1.6376015) q[2];
rz(-1.2184294) q[3];
sx q[3];
rz(-2.7259493) q[3];
sx q[3];
rz(-0.98201069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0686491) q[0];
sx q[0];
rz(-1.8953841) q[0];
sx q[0];
rz(0.15596998) q[0];
rz(-2.7929557) q[1];
sx q[1];
rz(-2.5376384) q[1];
sx q[1];
rz(1.9085931) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1240163) q[0];
sx q[0];
rz(-1.4498382) q[0];
sx q[0];
rz(1.308754) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2945588) q[2];
sx q[2];
rz(-2.4748908) q[2];
sx q[2];
rz(1.2444956) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0291042) q[1];
sx q[1];
rz(-1.8524516) q[1];
sx q[1];
rz(-2.2206578) q[1];
rz(1.5560914) q[3];
sx q[3];
rz(-1.3193519) q[3];
sx q[3];
rz(-0.78212839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4431346) q[2];
sx q[2];
rz(-3.105574) q[2];
sx q[2];
rz(1.5114463) q[2];
rz(2.002142) q[3];
sx q[3];
rz(-1.8099433) q[3];
sx q[3];
rz(-2.1512234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3646506) q[0];
sx q[0];
rz(-0.41780892) q[0];
sx q[0];
rz(1.9885709) q[0];
rz(-0.54368377) q[1];
sx q[1];
rz(-1.2666356) q[1];
sx q[1];
rz(0.26184729) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9898997) q[0];
sx q[0];
rz(-1.6542871) q[0];
sx q[0];
rz(0.033647353) q[0];
x q[1];
rz(-1.0864429) q[2];
sx q[2];
rz(-2.5005955) q[2];
sx q[2];
rz(-2.6015559) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8229554) q[1];
sx q[1];
rz(-1.793211) q[1];
sx q[1];
rz(2.5435124) q[1];
rz(-pi) q[2];
rz(1.366794) q[3];
sx q[3];
rz(-2.3424934) q[3];
sx q[3];
rz(2.4287756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5220945) q[2];
sx q[2];
rz(-2.5358989) q[2];
sx q[2];
rz(1.4052793) q[2];
rz(0.74553472) q[3];
sx q[3];
rz(-1.8757952) q[3];
sx q[3];
rz(-2.1982927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0010506823) q[0];
sx q[0];
rz(-0.11740919) q[0];
sx q[0];
rz(2.7897799) q[0];
rz(0.44031269) q[1];
sx q[1];
rz(-1.3654717) q[1];
sx q[1];
rz(-2.9702759) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7650314) q[0];
sx q[0];
rz(-0.91463551) q[0];
sx q[0];
rz(1.879346) q[0];
rz(-1.8753042) q[2];
sx q[2];
rz(-0.84104462) q[2];
sx q[2];
rz(-0.27308057) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.76207668) q[1];
sx q[1];
rz(-0.87240309) q[1];
sx q[1];
rz(1.5435436) q[1];
rz(-pi) q[2];
x q[2];
rz(0.76254179) q[3];
sx q[3];
rz(-2.2665215) q[3];
sx q[3];
rz(-0.97555893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.064284023) q[2];
sx q[2];
rz(-1.7694387) q[2];
sx q[2];
rz(-0.94179955) q[2];
rz(1.1549548) q[3];
sx q[3];
rz(-0.20808163) q[3];
sx q[3];
rz(1.8010767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44523859) q[0];
sx q[0];
rz(-2.0360763) q[0];
sx q[0];
rz(-0.0048986991) q[0];
rz(-2.694963) q[1];
sx q[1];
rz(-2.547867) q[1];
sx q[1];
rz(-0.68797025) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82179994) q[0];
sx q[0];
rz(-1.742779) q[0];
sx q[0];
rz(0.71806192) q[0];
rz(-pi) q[1];
x q[1];
rz(0.46779386) q[2];
sx q[2];
rz(-2.4681598) q[2];
sx q[2];
rz(1.3040257) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1157284) q[1];
sx q[1];
rz(-1.1394355) q[1];
sx q[1];
rz(-2.906745) q[1];
rz(-pi) q[2];
rz(0.68617188) q[3];
sx q[3];
rz(-0.99643512) q[3];
sx q[3];
rz(0.8818834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.530431) q[2];
sx q[2];
rz(-0.40803424) q[2];
sx q[2];
rz(-2.2116057) q[2];
rz(-1.9200578) q[3];
sx q[3];
rz(-1.113021) q[3];
sx q[3];
rz(2.5482224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4596443) q[0];
sx q[0];
rz(-1.821803) q[0];
sx q[0];
rz(0.090106877) q[0];
rz(-1.7168761) q[1];
sx q[1];
rz(-2.5017891) q[1];
sx q[1];
rz(2.7065014) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97532192) q[0];
sx q[0];
rz(-1.7164125) q[0];
sx q[0];
rz(-2.6377489) q[0];
rz(0.66799156) q[2];
sx q[2];
rz(-1.4282303) q[2];
sx q[2];
rz(3.0270456) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7104946) q[1];
sx q[1];
rz(-2.5113547) q[1];
sx q[1];
rz(1.3517387) q[1];
rz(-pi) q[2];
rz(3.1396542) q[3];
sx q[3];
rz(-2.7172305) q[3];
sx q[3];
rz(2.3821609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4789751) q[2];
sx q[2];
rz(-1.0650029) q[2];
sx q[2];
rz(0.62620658) q[2];
rz(0.19691697) q[3];
sx q[3];
rz(-2.1508689) q[3];
sx q[3];
rz(-1.3885434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5597252) q[0];
sx q[0];
rz(-1.6904866) q[0];
sx q[0];
rz(0.054280601) q[0];
rz(2.718603) q[1];
sx q[1];
rz(-1.3754247) q[1];
sx q[1];
rz(-1.8064226) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8568118) q[0];
sx q[0];
rz(-0.81997141) q[0];
sx q[0];
rz(2.3597844) q[0];
rz(-pi) q[1];
rz(2.1077431) q[2];
sx q[2];
rz(-2.5116911) q[2];
sx q[2];
rz(0.065187188) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.74354913) q[1];
sx q[1];
rz(-0.7760007) q[1];
sx q[1];
rz(0.77568027) q[1];
rz(0.068633462) q[3];
sx q[3];
rz(-1.4331927) q[3];
sx q[3];
rz(-0.28561628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6173031) q[2];
sx q[2];
rz(-1.6733988) q[2];
sx q[2];
rz(-0.35150251) q[2];
rz(1.7678123) q[3];
sx q[3];
rz(-2.6280554) q[3];
sx q[3];
rz(-2.8721749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
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
rz(-1.6108797) q[1];
sx q[1];
rz(-2.0738475) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7368374) q[0];
sx q[0];
rz(-2.1421297) q[0];
sx q[0];
rz(0.57147567) q[0];
x q[1];
rz(0.2944417) q[2];
sx q[2];
rz(-1.3047555) q[2];
sx q[2];
rz(0.58318116) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.23087654) q[1];
sx q[1];
rz(-1.9741804) q[1];
sx q[1];
rz(0.97539263) q[1];
rz(-pi) q[2];
rz(-2.8078733) q[3];
sx q[3];
rz(-0.90147831) q[3];
sx q[3];
rz(-0.02726083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4296253) q[2];
sx q[2];
rz(-1.2099313) q[2];
sx q[2];
rz(-0.21027002) q[2];
rz(-0.78768864) q[3];
sx q[3];
rz(-1.7588047) q[3];
sx q[3];
rz(-0.53019607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44337153) q[0];
sx q[0];
rz(-2.5826695) q[0];
sx q[0];
rz(-1.2179751) q[0];
rz(-1.632985) q[1];
sx q[1];
rz(-1.2651545) q[1];
sx q[1];
rz(-1.9427585) q[1];
rz(-0.40661033) q[2];
sx q[2];
rz(-2.2175773) q[2];
sx q[2];
rz(-2.8178136) q[2];
rz(-2.9859424) q[3];
sx q[3];
rz(-1.2990824) q[3];
sx q[3];
rz(-1.6344447) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
