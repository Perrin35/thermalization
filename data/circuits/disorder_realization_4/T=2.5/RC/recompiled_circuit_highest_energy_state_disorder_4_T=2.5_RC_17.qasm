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
rz(1.7734111) q[0];
rz(2.6693681) q[1];
sx q[1];
rz(-0.11657403) q[1];
sx q[1];
rz(-2.9393458) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8919173) q[0];
sx q[0];
rz(-1.37171) q[0];
sx q[0];
rz(-0.21659918) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0482685) q[2];
sx q[2];
rz(-1.1236345) q[2];
sx q[2];
rz(-1.9856404) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3766915) q[1];
sx q[1];
rz(-0.67923949) q[1];
sx q[1];
rz(1.1107854) q[1];
rz(-pi) q[2];
rz(1.1047045) q[3];
sx q[3];
rz(-1.6045609) q[3];
sx q[3];
rz(-2.1347783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5998485) q[2];
sx q[2];
rz(-0.84459633) q[2];
sx q[2];
rz(2.2005626) q[2];
rz(2.8372676) q[3];
sx q[3];
rz(-0.86186886) q[3];
sx q[3];
rz(0.26722515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32352725) q[0];
sx q[0];
rz(-2.7637988) q[0];
sx q[0];
rz(-2.3622808) q[0];
rz(1.7794973) q[1];
sx q[1];
rz(-2.7423488) q[1];
sx q[1];
rz(2.0077226) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72778217) q[0];
sx q[0];
rz(-1.6690133) q[0];
sx q[0];
rz(-2.006524) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0236565) q[2];
sx q[2];
rz(-1.5124161) q[2];
sx q[2];
rz(2.749325) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.459455) q[1];
sx q[1];
rz(-1.3054779) q[1];
sx q[1];
rz(3.0299761) q[1];
x q[2];
rz(-3.1406088) q[3];
sx q[3];
rz(-2.0800701) q[3];
sx q[3];
rz(2.5343717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.69089943) q[2];
sx q[2];
rz(-2.1096114) q[2];
sx q[2];
rz(2.9449985) q[2];
rz(-0.96902668) q[3];
sx q[3];
rz(-1.6220379) q[3];
sx q[3];
rz(-1.2381747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3900688) q[0];
sx q[0];
rz(-0.5539493) q[0];
sx q[0];
rz(-2.4098136) q[0];
rz(0.36830184) q[1];
sx q[1];
rz(-2.383547) q[1];
sx q[1];
rz(0.36645737) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5598504) q[0];
sx q[0];
rz(-1.3936211) q[0];
sx q[0];
rz(-0.26818256) q[0];
rz(-pi) q[1];
rz(0.75080504) q[2];
sx q[2];
rz(-2.8158203) q[2];
sx q[2];
rz(-0.70394403) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4225821) q[1];
sx q[1];
rz(-1.1478462) q[1];
sx q[1];
rz(2.1450536) q[1];
rz(-pi) q[2];
rz(2.5440574) q[3];
sx q[3];
rz(-2.1794277) q[3];
sx q[3];
rz(-3.1239005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1825819) q[2];
sx q[2];
rz(-2.3704447) q[2];
sx q[2];
rz(-2.2491573) q[2];
rz(-2.6871032) q[3];
sx q[3];
rz(-2.317704) q[3];
sx q[3];
rz(0.23196001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7682122) q[0];
sx q[0];
rz(-0.1425655) q[0];
sx q[0];
rz(0.40661231) q[0];
rz(0.82798249) q[1];
sx q[1];
rz(-1.7384638) q[1];
sx q[1];
rz(-0.83555317) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6064049) q[0];
sx q[0];
rz(-1.390103) q[0];
sx q[0];
rz(3.101493) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9172952) q[2];
sx q[2];
rz(-1.76148) q[2];
sx q[2];
rz(2.4257223) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.16627993) q[1];
sx q[1];
rz(-1.1381835) q[1];
sx q[1];
rz(1.7325425) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3341414) q[3];
sx q[3];
rz(-0.39517212) q[3];
sx q[3];
rz(-0.69505668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.39997175) q[2];
sx q[2];
rz(-0.30210945) q[2];
sx q[2];
rz(0.92140222) q[2];
rz(1.7737927) q[3];
sx q[3];
rz(-2.1801345) q[3];
sx q[3];
rz(-0.14149806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2530186) q[0];
sx q[0];
rz(-1.0971917) q[0];
sx q[0];
rz(-2.9827523) q[0];
rz(-1.6244434) q[1];
sx q[1];
rz(-0.30888638) q[1];
sx q[1];
rz(-1.812017) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8027199) q[0];
sx q[0];
rz(-3.1204528) q[0];
sx q[0];
rz(0.52864023) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9407523) q[2];
sx q[2];
rz(-2.5479377) q[2];
sx q[2];
rz(-1.6933189) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5160421) q[1];
sx q[1];
rz(-0.4333638) q[1];
sx q[1];
rz(2.0651555) q[1];
x q[2];
rz(2.1383189) q[3];
sx q[3];
rz(-1.7470932) q[3];
sx q[3];
rz(0.76667537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.57881957) q[2];
sx q[2];
rz(-2.4476738) q[2];
sx q[2];
rz(2.5895183) q[2];
rz(-0.79365802) q[3];
sx q[3];
rz(-0.59760439) q[3];
sx q[3];
rz(-2.1551267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3119222) q[0];
sx q[0];
rz(-2.2770918) q[0];
sx q[0];
rz(0.70575869) q[0];
rz(0.13748473) q[1];
sx q[1];
rz(-2.4973713) q[1];
sx q[1];
rz(-2.4286043) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1437837) q[0];
sx q[0];
rz(-1.1857872) q[0];
sx q[0];
rz(1.3237756) q[0];
rz(-pi) q[1];
rz(2.4586468) q[2];
sx q[2];
rz(-0.6295528) q[2];
sx q[2];
rz(-0.94972875) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.56604993) q[1];
sx q[1];
rz(-2.7311196) q[1];
sx q[1];
rz(1.3259352) q[1];
x q[2];
rz(1.8637795) q[3];
sx q[3];
rz(-1.4594541) q[3];
sx q[3];
rz(0.35929832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2116427) q[2];
sx q[2];
rz(-0.86923081) q[2];
sx q[2];
rz(-2.8575274) q[2];
rz(-3.0197213) q[3];
sx q[3];
rz(-2.8594696) q[3];
sx q[3];
rz(1.6596644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
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
rz(0.44037285) q[1];
sx q[1];
rz(-0.71763867) q[1];
sx q[1];
rz(-1.3622989) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4934568) q[0];
sx q[0];
rz(-1.3246944) q[0];
sx q[0];
rz(-1.8561234) q[0];
x q[1];
rz(-2.2978034) q[2];
sx q[2];
rz(-1.8966881) q[2];
sx q[2];
rz(-1.0032723) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3063115) q[1];
sx q[1];
rz(-2.0871442) q[1];
sx q[1];
rz(2.5912214) q[1];
x q[2];
rz(2.6853722) q[3];
sx q[3];
rz(-1.4672605) q[3];
sx q[3];
rz(0.52952784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1661487) q[2];
sx q[2];
rz(-0.41836172) q[2];
sx q[2];
rz(-0.56900209) q[2];
rz(-1.0373908) q[3];
sx q[3];
rz(-0.85809696) q[3];
sx q[3];
rz(-0.4666127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5910832) q[0];
sx q[0];
rz(-2.789848) q[0];
sx q[0];
rz(-1.2048703) q[0];
rz(-0.51271802) q[1];
sx q[1];
rz(-2.3725489) q[1];
sx q[1];
rz(-0.18096322) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9009243) q[0];
sx q[0];
rz(-0.5597194) q[0];
sx q[0];
rz(2.3842035) q[0];
x q[1];
rz(2.7788073) q[2];
sx q[2];
rz(-0.37172592) q[2];
sx q[2];
rz(-0.065601018) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6540203) q[1];
sx q[1];
rz(-2.1250238) q[1];
sx q[1];
rz(-2.8934862) q[1];
x q[2];
rz(1.6445141) q[3];
sx q[3];
rz(-1.4604919) q[3];
sx q[3];
rz(0.10681399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5391431) q[2];
sx q[2];
rz(-2.4931543) q[2];
sx q[2];
rz(-0.096972018) q[2];
rz(2.5868331) q[3];
sx q[3];
rz(-1.8192889) q[3];
sx q[3];
rz(-2.7831912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48034126) q[0];
sx q[0];
rz(-2.3619409) q[0];
sx q[0];
rz(0.10777792) q[0];
rz(-0.55094552) q[1];
sx q[1];
rz(-0.44083732) q[1];
sx q[1];
rz(1.617618) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4827037) q[0];
sx q[0];
rz(-1.3388757) q[0];
sx q[0];
rz(-3.0085682) q[0];
x q[1];
rz(-2.5059087) q[2];
sx q[2];
rz(-0.7506461) q[2];
sx q[2];
rz(-2.8849066) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.35573449) q[1];
sx q[1];
rz(-1.6140892) q[1];
sx q[1];
rz(-2.0329352) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.518256) q[3];
sx q[3];
rz(-2.0618304) q[3];
sx q[3];
rz(-0.72140104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.23065755) q[2];
sx q[2];
rz(-2.6378938) q[2];
sx q[2];
rz(-2.0111734) q[2];
rz(-0.23760992) q[3];
sx q[3];
rz(-0.41046023) q[3];
sx q[3];
rz(0.75806481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1199353) q[0];
sx q[0];
rz(-0.1788685) q[0];
sx q[0];
rz(-2.2005431) q[0];
rz(0.36048105) q[1];
sx q[1];
rz(-1.4070114) q[1];
sx q[1];
rz(-1.2233268) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0377285) q[0];
sx q[0];
rz(-0.99261221) q[0];
sx q[0];
rz(-1.8072771) q[0];
x q[1];
rz(-2.1603196) q[2];
sx q[2];
rz(-0.86893493) q[2];
sx q[2];
rz(-1.5379932) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4810923) q[1];
sx q[1];
rz(-1.5005439) q[1];
sx q[1];
rz(2.7226647) q[1];
rz(-pi) q[2];
rz(2.1072708) q[3];
sx q[3];
rz(-1.7400898) q[3];
sx q[3];
rz(0.062549755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2337522) q[2];
sx q[2];
rz(-1.7641822) q[2];
sx q[2];
rz(-3.0506548) q[2];
rz(2.9371373) q[3];
sx q[3];
rz(-0.8046059) q[3];
sx q[3];
rz(-1.0819775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7627056) q[0];
sx q[0];
rz(-1.5305516) q[0];
sx q[0];
rz(-1.332921) q[0];
rz(-1.7145722) q[1];
sx q[1];
rz(-0.49976977) q[1];
sx q[1];
rz(1.5907092) q[1];
rz(-2.6743307) q[2];
sx q[2];
rz(-0.89667761) q[2];
sx q[2];
rz(-2.1368334) q[2];
rz(0.25973917) q[3];
sx q[3];
rz(-0.22077122) q[3];
sx q[3];
rz(2.1833899) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
