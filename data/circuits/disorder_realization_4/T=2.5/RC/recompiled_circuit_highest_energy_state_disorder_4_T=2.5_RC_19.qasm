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
rz(-2.0877617) q[0];
sx q[0];
rz(-0.55019903) q[0];
sx q[0];
rz(-1.3681816) q[0];
rz(2.6693681) q[1];
sx q[1];
rz(-0.11657403) q[1];
sx q[1];
rz(-2.9393458) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8919173) q[0];
sx q[0];
rz(-1.7698827) q[0];
sx q[0];
rz(0.21659918) q[0];
rz(-2.3362581) q[2];
sx q[2];
rz(-2.4676358) q[2];
sx q[2];
rz(1.0588624) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3766915) q[1];
sx q[1];
rz(-0.67923949) q[1];
sx q[1];
rz(1.1107854) q[1];
rz(-pi) q[2];
rz(-1.1047045) q[3];
sx q[3];
rz(-1.6045609) q[3];
sx q[3];
rz(-1.0068144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5998485) q[2];
sx q[2];
rz(-2.2969963) q[2];
sx q[2];
rz(0.94103003) q[2];
rz(0.30432501) q[3];
sx q[3];
rz(-2.2797238) q[3];
sx q[3];
rz(-2.8743675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8180654) q[0];
sx q[0];
rz(-0.37779385) q[0];
sx q[0];
rz(-0.7793119) q[0];
rz(-1.3620954) q[1];
sx q[1];
rz(-0.39924386) q[1];
sx q[1];
rz(-2.0077226) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5061653) q[0];
sx q[0];
rz(-0.44596684) q[0];
sx q[0];
rz(-1.3414499) q[0];
rz(-3.0236565) q[2];
sx q[2];
rz(-1.6291766) q[2];
sx q[2];
rz(-2.749325) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0596327) q[1];
sx q[1];
rz(-1.6784918) q[1];
sx q[1];
rz(-1.3038941) q[1];
rz(-pi) q[2];
rz(-1.5690345) q[3];
sx q[3];
rz(-0.50927466) q[3];
sx q[3];
rz(-0.60520303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.69089943) q[2];
sx q[2];
rz(-1.0319812) q[2];
sx q[2];
rz(0.1965941) q[2];
rz(-0.96902668) q[3];
sx q[3];
rz(-1.6220379) q[3];
sx q[3];
rz(1.903418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3900688) q[0];
sx q[0];
rz(-2.5876434) q[0];
sx q[0];
rz(-0.73177904) q[0];
rz(0.36830184) q[1];
sx q[1];
rz(-0.75804561) q[1];
sx q[1];
rz(-0.36645737) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.059341934) q[0];
sx q[0];
rz(-1.8346795) q[0];
sx q[0];
rz(1.3871971) q[0];
rz(-pi) q[1];
rz(-0.24213893) q[2];
sx q[2];
rz(-1.3506827) q[2];
sx q[2];
rz(1.5907703) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7190105) q[1];
sx q[1];
rz(-1.9937464) q[1];
sx q[1];
rz(-0.99653901) q[1];
rz(-2.5440574) q[3];
sx q[3];
rz(-0.96216494) q[3];
sx q[3];
rz(-3.1239005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9590108) q[2];
sx q[2];
rz(-2.3704447) q[2];
sx q[2];
rz(-2.2491573) q[2];
rz(-2.6871032) q[3];
sx q[3];
rz(-2.317704) q[3];
sx q[3];
rz(-2.9096326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7682122) q[0];
sx q[0];
rz(-2.9990271) q[0];
sx q[0];
rz(-0.40661231) q[0];
rz(-0.82798249) q[1];
sx q[1];
rz(-1.7384638) q[1];
sx q[1];
rz(0.83555317) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0428187) q[0];
sx q[0];
rz(-1.6102428) q[0];
sx q[0];
rz(-1.7516319) q[0];
x q[1];
rz(2.9391772) q[2];
sx q[2];
rz(-1.9107585) q[2];
sx q[2];
rz(0.92325231) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9367076) q[1];
sx q[1];
rz(-2.6815273) q[1];
sx q[1];
rz(-0.33554828) q[1];
x q[2];
rz(-1.3341414) q[3];
sx q[3];
rz(-2.7464205) q[3];
sx q[3];
rz(0.69505668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7416209) q[2];
sx q[2];
rz(-0.30210945) q[2];
sx q[2];
rz(2.2201904) q[2];
rz(1.3678) q[3];
sx q[3];
rz(-2.1801345) q[3];
sx q[3];
rz(-3.0000946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88857404) q[0];
sx q[0];
rz(-2.044401) q[0];
sx q[0];
rz(-2.9827523) q[0];
rz(1.5171492) q[1];
sx q[1];
rz(-2.8327063) q[1];
sx q[1];
rz(1.812017) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3388728) q[0];
sx q[0];
rz(-0.021139806) q[0];
sx q[0];
rz(-2.6129524) q[0];
x q[1];
rz(-1.9407523) q[2];
sx q[2];
rz(-2.5479377) q[2];
sx q[2];
rz(-1.4482738) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6318887) q[1];
sx q[1];
rz(-1.77138) q[1];
sx q[1];
rz(-1.9575809) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8908126) q[3];
sx q[3];
rz(-2.5502) q[3];
sx q[3];
rz(2.0689912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.57881957) q[2];
sx q[2];
rz(-2.4476738) q[2];
sx q[2];
rz(-2.5895183) q[2];
rz(2.3479346) q[3];
sx q[3];
rz(-2.5439883) q[3];
sx q[3];
rz(2.1551267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3119222) q[0];
sx q[0];
rz(-2.2770918) q[0];
sx q[0];
rz(-2.435834) q[0];
rz(-3.0041079) q[1];
sx q[1];
rz(-2.4973713) q[1];
sx q[1];
rz(-2.4286043) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5891083) q[0];
sx q[0];
rz(-2.6875067) q[0];
sx q[0];
rz(-0.54291351) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.001686) q[2];
sx q[2];
rz(-2.0451114) q[2];
sx q[2];
rz(-0.16109078) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9116152) q[1];
sx q[1];
rz(-1.4739081) q[1];
sx q[1];
rz(-1.9703034) q[1];
x q[2];
rz(-1.9401523) q[3];
sx q[3];
rz(-0.31285646) q[3];
sx q[3];
rz(-2.2830056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2116427) q[2];
sx q[2];
rz(-0.86923081) q[2];
sx q[2];
rz(2.8575274) q[2];
rz(-0.12187135) q[3];
sx q[3];
rz(-2.8594696) q[3];
sx q[3];
rz(1.4819283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26509869) q[0];
sx q[0];
rz(-0.6186741) q[0];
sx q[0];
rz(0.70736831) q[0];
rz(-2.7012198) q[1];
sx q[1];
rz(-2.423954) q[1];
sx q[1];
rz(1.3622989) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4934568) q[0];
sx q[0];
rz(-1.8168983) q[0];
sx q[0];
rz(-1.8561234) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7168305) q[2];
sx q[2];
rz(-2.2519037) q[2];
sx q[2];
rz(2.8515138) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.41304663) q[1];
sx q[1];
rz(-2.4056681) q[1];
sx q[1];
rz(-0.82637431) q[1];
rz(-pi) q[2];
rz(-1.6860289) q[3];
sx q[3];
rz(-1.1172022) q[3];
sx q[3];
rz(0.99059243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1661487) q[2];
sx q[2];
rz(-0.41836172) q[2];
sx q[2];
rz(-2.5725906) q[2];
rz(1.0373908) q[3];
sx q[3];
rz(-0.85809696) q[3];
sx q[3];
rz(0.4666127) q[3];
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
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55050945) q[0];
sx q[0];
rz(-0.35174462) q[0];
sx q[0];
rz(-1.2048703) q[0];
rz(0.51271802) q[1];
sx q[1];
rz(-2.3725489) q[1];
sx q[1];
rz(-2.9606294) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9009243) q[0];
sx q[0];
rz(-0.5597194) q[0];
sx q[0];
rz(-2.3842035) q[0];
rz(-pi) q[1];
rz(-0.34951194) q[2];
sx q[2];
rz(-1.4415359) q[2];
sx q[2];
rz(-1.2964647) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1026336) q[1];
sx q[1];
rz(-0.60188434) q[1];
sx q[1];
rz(1.9485056) q[1];
rz(-pi) q[2];
x q[2];
rz(0.11060235) q[3];
sx q[3];
rz(-1.6440653) q[3];
sx q[3];
rz(1.6694809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.60244954) q[2];
sx q[2];
rz(-0.64843833) q[2];
sx q[2];
rz(-3.0446206) q[2];
rz(2.5868331) q[3];
sx q[3];
rz(-1.8192889) q[3];
sx q[3];
rz(-2.7831912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48034126) q[0];
sx q[0];
rz(-0.77965176) q[0];
sx q[0];
rz(3.0338147) q[0];
rz(-2.5906471) q[1];
sx q[1];
rz(-0.44083732) q[1];
sx q[1];
rz(1.5239747) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4827037) q[0];
sx q[0];
rz(-1.802717) q[0];
sx q[0];
rz(3.0085682) q[0];
rz(-0.63568398) q[2];
sx q[2];
rz(-0.7506461) q[2];
sx q[2];
rz(-0.25668609) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.35573449) q[1];
sx q[1];
rz(-1.6140892) q[1];
sx q[1];
rz(2.0329352) q[1];
rz(-pi) q[2];
rz(3.0436936) q[3];
sx q[3];
rz(-2.6479841) q[3];
sx q[3];
rz(2.5312587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.23065755) q[2];
sx q[2];
rz(-0.50369889) q[2];
sx q[2];
rz(-1.1304193) q[2];
rz(2.9039827) q[3];
sx q[3];
rz(-0.41046023) q[3];
sx q[3];
rz(-2.3835278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1199353) q[0];
sx q[0];
rz(-0.1788685) q[0];
sx q[0];
rz(-2.2005431) q[0];
rz(-2.7811116) q[1];
sx q[1];
rz(-1.4070114) q[1];
sx q[1];
rz(-1.2233268) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1038641) q[0];
sx q[0];
rz(-2.1489804) q[0];
sx q[0];
rz(-1.3343156) q[0];
rz(0.5816831) q[2];
sx q[2];
rz(-0.8831199) q[2];
sx q[2];
rz(2.4066642) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0666484) q[1];
sx q[1];
rz(-0.42443353) q[1];
sx q[1];
rz(-2.9702999) q[1];
x q[2];
rz(0.19630614) q[3];
sx q[3];
rz(-2.098791) q[3];
sx q[3];
rz(1.4083901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2337522) q[2];
sx q[2];
rz(-1.3774104) q[2];
sx q[2];
rz(-0.09093786) q[2];
rz(-2.9371373) q[3];
sx q[3];
rz(-0.8046059) q[3];
sx q[3];
rz(-2.0596152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37888708) q[0];
sx q[0];
rz(-1.5305516) q[0];
sx q[0];
rz(-1.332921) q[0];
rz(1.7145722) q[1];
sx q[1];
rz(-2.6418229) q[1];
sx q[1];
rz(-1.5508834) q[1];
rz(0.84080055) q[2];
sx q[2];
rz(-1.2111831) q[2];
sx q[2];
rz(2.2704587) q[2];
rz(-1.6283725) q[3];
sx q[3];
rz(-1.3575469) q[3];
sx q[3];
rz(-1.224106) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
