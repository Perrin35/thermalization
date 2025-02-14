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
rz(-1.2019914) q[0];
sx q[0];
rz(-2.658598) q[0];
sx q[0];
rz(1.510409) q[0];
rz(-0.13934879) q[1];
sx q[1];
rz(5.723602) q[1];
sx q[1];
rz(8.6948123) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6522927) q[0];
sx q[0];
rz(-2.518939) q[0];
sx q[0];
rz(-1.3176509) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3454622) q[2];
sx q[2];
rz(-0.98774922) q[2];
sx q[2];
rz(-0.096681373) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7429205) q[1];
sx q[1];
rz(-1.6026596) q[1];
sx q[1];
rz(-0.58018654) q[1];
rz(-pi) q[2];
rz(-2.2375305) q[3];
sx q[3];
rz(-2.781467) q[3];
sx q[3];
rz(1.2670497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1175179) q[2];
sx q[2];
rz(-2.8705609) q[2];
sx q[2];
rz(-0.81895858) q[2];
rz(-3.135318) q[3];
sx q[3];
rz(-1.9082021) q[3];
sx q[3];
rz(0.98244572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55945021) q[0];
sx q[0];
rz(-1.5511976) q[0];
sx q[0];
rz(0.5994125) q[0];
rz(-0.86743152) q[1];
sx q[1];
rz(-2.1116833) q[1];
sx q[1];
rz(0.74554602) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.065577995) q[0];
sx q[0];
rz(-1.4434555) q[0];
sx q[0];
rz(-0.29177427) q[0];
rz(2.8023984) q[2];
sx q[2];
rz(-0.789398) q[2];
sx q[2];
rz(1.9292924) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2170495) q[1];
sx q[1];
rz(-2.1991208) q[1];
sx q[1];
rz(-1.7679067) q[1];
rz(0.98394139) q[3];
sx q[3];
rz(-0.7823669) q[3];
sx q[3];
rz(-0.38228211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6856689) q[2];
sx q[2];
rz(-0.29752877) q[2];
sx q[2];
rz(1.4073184) q[2];
rz(0.84434858) q[3];
sx q[3];
rz(-0.83365369) q[3];
sx q[3];
rz(1.0367397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5832962) q[0];
sx q[0];
rz(-1.946227) q[0];
sx q[0];
rz(2.1287647) q[0];
rz(-2.5335675) q[1];
sx q[1];
rz(-1.5589747) q[1];
sx q[1];
rz(-1.8720522) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2484754) q[0];
sx q[0];
rz(-0.60113827) q[0];
sx q[0];
rz(0.84971835) q[0];
x q[1];
rz(-0.89983268) q[2];
sx q[2];
rz(-1.6353893) q[2];
sx q[2];
rz(3.10499) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.041965466) q[1];
sx q[1];
rz(-1.7668952) q[1];
sx q[1];
rz(2.4134497) q[1];
rz(-pi) q[2];
rz(1.9023015) q[3];
sx q[3];
rz(-2.0391782) q[3];
sx q[3];
rz(1.4862332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6262007) q[2];
sx q[2];
rz(-2.0833368) q[2];
sx q[2];
rz(1.2916279) q[2];
rz(0.0083222566) q[3];
sx q[3];
rz(-2.1885927) q[3];
sx q[3];
rz(2.9279809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0030647) q[0];
sx q[0];
rz(-0.55532885) q[0];
sx q[0];
rz(1.3767161) q[0];
rz(1.6142913) q[1];
sx q[1];
rz(-1.2034028) q[1];
sx q[1];
rz(-0.52070224) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10258612) q[0];
sx q[0];
rz(-1.9780567) q[0];
sx q[0];
rz(-0.38577052) q[0];
x q[1];
rz(2.7815656) q[2];
sx q[2];
rz(-1.0485149) q[2];
sx q[2];
rz(-1.479666) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5129052) q[1];
sx q[1];
rz(-0.73343745) q[1];
sx q[1];
rz(2.4705486) q[1];
rz(2.2124452) q[3];
sx q[3];
rz(-0.7693253) q[3];
sx q[3];
rz(1.820946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.5248096) q[2];
sx q[2];
rz(-2.3526134) q[2];
sx q[2];
rz(1.7899803) q[2];
rz(-2.1221519) q[3];
sx q[3];
rz(-0.56988684) q[3];
sx q[3];
rz(1.9701689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.549642) q[0];
sx q[0];
rz(-1.6787981) q[0];
sx q[0];
rz(-0.098966448) q[0];
rz(0.70676604) q[1];
sx q[1];
rz(-2.2711429) q[1];
sx q[1];
rz(-0.4695355) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75824518) q[0];
sx q[0];
rz(-1.4329785) q[0];
sx q[0];
rz(-0.024017781) q[0];
rz(-pi) q[1];
rz(-2.1624998) q[2];
sx q[2];
rz(-2.5212477) q[2];
sx q[2];
rz(-1.9949335) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5597804) q[1];
sx q[1];
rz(-2.1511934) q[1];
sx q[1];
rz(-2.0349166) q[1];
x q[2];
rz(-1.4245701) q[3];
sx q[3];
rz(-2.1220045) q[3];
sx q[3];
rz(0.92155462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.88064757) q[2];
sx q[2];
rz(-1.8069043) q[2];
sx q[2];
rz(-1.8355231) q[2];
rz(2.7169054) q[3];
sx q[3];
rz(-1.7644019) q[3];
sx q[3];
rz(-0.88596058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0223618) q[0];
sx q[0];
rz(-2.5194118) q[0];
sx q[0];
rz(1.3223883) q[0];
rz(2.0139096) q[1];
sx q[1];
rz(-1.6219982) q[1];
sx q[1];
rz(-1.3816396) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2974325) q[0];
sx q[0];
rz(-1.4844262) q[0];
sx q[0];
rz(-0.81120904) q[0];
x q[1];
rz(-0.83397978) q[2];
sx q[2];
rz(-2.4206941) q[2];
sx q[2];
rz(2.3490459) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4037212) q[1];
sx q[1];
rz(-0.779169) q[1];
sx q[1];
rz(2.8761151) q[1];
rz(-pi) q[2];
rz(-0.30140169) q[3];
sx q[3];
rz(-1.5001138) q[3];
sx q[3];
rz(1.4247198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3809001) q[2];
sx q[2];
rz(-1.6221294) q[2];
sx q[2];
rz(-0.48119989) q[2];
rz(0.5091269) q[3];
sx q[3];
rz(-0.23734084) q[3];
sx q[3];
rz(-2.7595162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5612438) q[0];
sx q[0];
rz(-2.6600397) q[0];
sx q[0];
rz(-3.0522108) q[0];
rz(-1.0786169) q[1];
sx q[1];
rz(-1.6856472) q[1];
sx q[1];
rz(-2.1481029) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8315961) q[0];
sx q[0];
rz(-2.352) q[0];
sx q[0];
rz(1.5168651) q[0];
rz(-pi) q[1];
x q[1];
rz(0.87331302) q[2];
sx q[2];
rz(-2.2255579) q[2];
sx q[2];
rz(-2.8755434) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3741039) q[1];
sx q[1];
rz(-1.2084157) q[1];
sx q[1];
rz(0.15717536) q[1];
x q[2];
rz(-2.4782789) q[3];
sx q[3];
rz(-0.74771008) q[3];
sx q[3];
rz(1.472689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3420458) q[2];
sx q[2];
rz(-0.61508238) q[2];
sx q[2];
rz(2.7395524) q[2];
rz(-1.0953995) q[3];
sx q[3];
rz(-1.9270555) q[3];
sx q[3];
rz(-2.5636165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6680229) q[0];
sx q[0];
rz(-0.80900017) q[0];
sx q[0];
rz(1.3336257) q[0];
rz(2.2857621) q[1];
sx q[1];
rz(-1.7355093) q[1];
sx q[1];
rz(2.2241101) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8569023) q[0];
sx q[0];
rz(-1.7698341) q[0];
sx q[0];
rz(3.0375098) q[0];
rz(-pi) q[1];
rz(-2.1945004) q[2];
sx q[2];
rz(-1.4741885) q[2];
sx q[2];
rz(-1.4132702) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.11115796) q[1];
sx q[1];
rz(-1.3492172) q[1];
sx q[1];
rz(-1.2171924) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.6360953) q[3];
sx q[3];
rz(-2.3878752) q[3];
sx q[3];
rz(-1.5223283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.63460073) q[2];
sx q[2];
rz(-1.4054106) q[2];
sx q[2];
rz(1.8514006) q[2];
rz(-0.90977943) q[3];
sx q[3];
rz(-1.4434283) q[3];
sx q[3];
rz(-3.1006052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9057587) q[0];
sx q[0];
rz(-1.1474778) q[0];
sx q[0];
rz(-0.27715096) q[0];
rz(1.2241036) q[1];
sx q[1];
rz(-1.6033019) q[1];
sx q[1];
rz(1.7701497) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7200206) q[0];
sx q[0];
rz(-2.4457481) q[0];
sx q[0];
rz(2.8369342) q[0];
rz(-pi) q[1];
x q[1];
rz(0.54746898) q[2];
sx q[2];
rz(-1.1927422) q[2];
sx q[2];
rz(0.34502703) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.692457) q[1];
sx q[1];
rz(-2.4583011) q[1];
sx q[1];
rz(1.9016674) q[1];
rz(1.3731433) q[3];
sx q[3];
rz(-2.4098793) q[3];
sx q[3];
rz(-0.55303516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.93854967) q[2];
sx q[2];
rz(-0.86157346) q[2];
sx q[2];
rz(2.4533563) q[2];
rz(2.8271683) q[3];
sx q[3];
rz(-0.36948547) q[3];
sx q[3];
rz(-1.2615874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37453434) q[0];
sx q[0];
rz(-0.86807591) q[0];
sx q[0];
rz(-0.57149291) q[0];
rz(0.68069619) q[1];
sx q[1];
rz(-1.9711767) q[1];
sx q[1];
rz(0.64819711) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2637973) q[0];
sx q[0];
rz(-1.5938984) q[0];
sx q[0];
rz(3.1277411) q[0];
rz(-2.4785751) q[2];
sx q[2];
rz(-2.1923991) q[2];
sx q[2];
rz(-0.66689516) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.56420621) q[1];
sx q[1];
rz(-1.5485475) q[1];
sx q[1];
rz(-0.69907) q[1];
rz(0.33721029) q[3];
sx q[3];
rz(-0.54223947) q[3];
sx q[3];
rz(0.48619871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8615243) q[2];
sx q[2];
rz(-2.0678949) q[2];
sx q[2];
rz(-1.6746707) q[2];
rz(2.9649949) q[3];
sx q[3];
rz(-0.64591518) q[3];
sx q[3];
rz(-1.8471898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8661154) q[0];
sx q[0];
rz(-1.3963516) q[0];
sx q[0];
rz(1.8657952) q[0];
rz(-0.38446174) q[1];
sx q[1];
rz(-1.6356331) q[1];
sx q[1];
rz(2.5148139) q[1];
rz(1.2550884) q[2];
sx q[2];
rz(-1.2366625) q[2];
sx q[2];
rz(-1.2319215) q[2];
rz(-0.53388673) q[3];
sx q[3];
rz(-1.0718828) q[3];
sx q[3];
rz(-0.26477118) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
