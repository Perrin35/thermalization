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
rz(1.053831) q[0];
sx q[0];
rz(-2.5913936) q[0];
sx q[0];
rz(1.3681816) q[0];
rz(2.6693681) q[1];
sx q[1];
rz(-0.11657403) q[1];
sx q[1];
rz(-2.9393458) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36461386) q[0];
sx q[0];
rz(-1.783051) q[0];
sx q[0];
rz(-1.3670761) q[0];
x q[1];
rz(2.3362581) q[2];
sx q[2];
rz(-0.67395681) q[2];
sx q[2];
rz(-2.0827302) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7649012) q[1];
sx q[1];
rz(-0.67923949) q[1];
sx q[1];
rz(2.0308073) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6458166) q[3];
sx q[3];
rz(-2.6743691) q[3];
sx q[3];
rz(0.4969767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5998485) q[2];
sx q[2];
rz(-0.84459633) q[2];
sx q[2];
rz(0.94103003) q[2];
rz(-2.8372676) q[3];
sx q[3];
rz(-2.2797238) q[3];
sx q[3];
rz(-2.8743675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32352725) q[0];
sx q[0];
rz(-2.7637988) q[0];
sx q[0];
rz(2.3622808) q[0];
rz(1.7794973) q[1];
sx q[1];
rz(-0.39924386) q[1];
sx q[1];
rz(-2.0077226) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88863605) q[0];
sx q[0];
rz(-1.1373113) q[0];
sx q[0];
rz(-0.10826464) q[0];
x q[1];
rz(1.629584) q[2];
sx q[2];
rz(-1.4530621) q[2];
sx q[2];
rz(1.9699772) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2782005) q[1];
sx q[1];
rz(-0.28732936) q[1];
sx q[1];
rz(-1.1817688) q[1];
x q[2];
rz(-1.5725582) q[3];
sx q[3];
rz(-2.632318) q[3];
sx q[3];
rz(-0.60520303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.69089943) q[2];
sx q[2];
rz(-1.0319812) q[2];
sx q[2];
rz(-0.1965941) q[2];
rz(-0.96902668) q[3];
sx q[3];
rz(-1.6220379) q[3];
sx q[3];
rz(1.903418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75152385) q[0];
sx q[0];
rz(-2.5876434) q[0];
sx q[0];
rz(-2.4098136) q[0];
rz(0.36830184) q[1];
sx q[1];
rz(-0.75804561) q[1];
sx q[1];
rz(-0.36645737) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0822507) q[0];
sx q[0];
rz(-1.8346795) q[0];
sx q[0];
rz(1.7543955) q[0];
rz(-pi) q[1];
x q[1];
rz(0.75080504) q[2];
sx q[2];
rz(-0.32577235) q[2];
sx q[2];
rz(0.70394403) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0302387) q[1];
sx q[1];
rz(-1.0525647) q[1];
sx q[1];
rz(0.49211647) q[1];
rz(-pi) q[2];
rz(-2.2711804) q[3];
sx q[3];
rz(-1.0910209) q[3];
sx q[3];
rz(-1.9241672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1825819) q[2];
sx q[2];
rz(-2.3704447) q[2];
sx q[2];
rz(-0.89243531) q[2];
rz(-2.6871032) q[3];
sx q[3];
rz(-0.82388866) q[3];
sx q[3];
rz(-0.23196001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7682122) q[0];
sx q[0];
rz(-2.9990271) q[0];
sx q[0];
rz(2.7349803) q[0];
rz(-2.3136102) q[1];
sx q[1];
rz(-1.4031289) q[1];
sx q[1];
rz(0.83555317) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53518772) q[0];
sx q[0];
rz(-1.390103) q[0];
sx q[0];
rz(3.101493) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9391772) q[2];
sx q[2];
rz(-1.2308342) q[2];
sx q[2];
rz(-2.2183403) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4728188) q[1];
sx q[1];
rz(-1.4240648) q[1];
sx q[1];
rz(0.43763486) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8074513) q[3];
sx q[3];
rz(-0.39517212) q[3];
sx q[3];
rz(-0.69505668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7416209) q[2];
sx q[2];
rz(-2.8394832) q[2];
sx q[2];
rz(2.2201904) q[2];
rz(-1.3678) q[3];
sx q[3];
rz(-0.96145815) q[3];
sx q[3];
rz(-3.0000946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2530186) q[0];
sx q[0];
rz(-2.044401) q[0];
sx q[0];
rz(0.15884037) q[0];
rz(1.6244434) q[1];
sx q[1];
rz(-0.30888638) q[1];
sx q[1];
rz(1.812017) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3388728) q[0];
sx q[0];
rz(-3.1204528) q[0];
sx q[0];
rz(2.6129524) q[0];
x q[1];
rz(2.9022568) q[2];
sx q[2];
rz(-1.0221326) q[2];
sx q[2];
rz(1.2557097) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.98011866) q[1];
sx q[1];
rz(-1.9494281) q[1];
sx q[1];
rz(0.21610723) q[1];
x q[2];
rz(1.25078) q[3];
sx q[3];
rz(-2.5502) q[3];
sx q[3];
rz(1.0726014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5627731) q[2];
sx q[2];
rz(-2.4476738) q[2];
sx q[2];
rz(2.5895183) q[2];
rz(0.79365802) q[3];
sx q[3];
rz(-2.5439883) q[3];
sx q[3];
rz(-2.1551267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3119222) q[0];
sx q[0];
rz(-2.2770918) q[0];
sx q[0];
rz(2.435834) q[0];
rz(3.0041079) q[1];
sx q[1];
rz(-0.64422137) q[1];
sx q[1];
rz(-2.4286043) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5214382) q[0];
sx q[0];
rz(-1.7993986) q[0];
sx q[0];
rz(-2.7457353) q[0];
rz(-0.6829459) q[2];
sx q[2];
rz(-2.5120398) q[2];
sx q[2];
rz(-2.1918639) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2299774) q[1];
sx q[1];
rz(-1.4739081) q[1];
sx q[1];
rz(1.1712892) q[1];
rz(-pi) q[2];
rz(-1.2778132) q[3];
sx q[3];
rz(-1.4594541) q[3];
sx q[3];
rz(-2.7822943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9299499) q[2];
sx q[2];
rz(-2.2723618) q[2];
sx q[2];
rz(2.8575274) q[2];
rz(-3.0197213) q[3];
sx q[3];
rz(-0.28212306) q[3];
sx q[3];
rz(1.4819283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26509869) q[0];
sx q[0];
rz(-2.5229186) q[0];
sx q[0];
rz(-0.70736831) q[0];
rz(-2.7012198) q[1];
sx q[1];
rz(-2.423954) q[1];
sx q[1];
rz(1.3622989) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61575261) q[0];
sx q[0];
rz(-2.7670015) q[0];
sx q[0];
rz(-0.84217854) q[0];
x q[1];
rz(2.2978034) q[2];
sx q[2];
rz(-1.8966881) q[2];
sx q[2];
rz(-2.1383204) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.728546) q[1];
sx q[1];
rz(-2.4056681) q[1];
sx q[1];
rz(2.3152183) q[1];
rz(1.4555638) q[3];
sx q[3];
rz(-1.1172022) q[3];
sx q[3];
rz(0.99059243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.97544396) q[2];
sx q[2];
rz(-0.41836172) q[2];
sx q[2];
rz(-2.5725906) q[2];
rz(-2.1042018) q[3];
sx q[3];
rz(-0.85809696) q[3];
sx q[3];
rz(0.4666127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55050945) q[0];
sx q[0];
rz(-0.35174462) q[0];
sx q[0];
rz(-1.9367223) q[0];
rz(2.6288746) q[1];
sx q[1];
rz(-0.7690438) q[1];
sx q[1];
rz(0.18096322) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9009243) q[0];
sx q[0];
rz(-2.5818733) q[0];
sx q[0];
rz(-2.3842035) q[0];
rz(-pi) q[1];
rz(2.7920807) q[2];
sx q[2];
rz(-1.7000568) q[2];
sx q[2];
rz(-1.845128) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.21576439) q[1];
sx q[1];
rz(-1.7811532) q[1];
sx q[1];
rz(-1.0025566) q[1];
rz(1.4970786) q[3];
sx q[3];
rz(-1.4604919) q[3];
sx q[3];
rz(-0.10681399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5391431) q[2];
sx q[2];
rz(-0.64843833) q[2];
sx q[2];
rz(3.0446206) q[2];
rz(-2.5868331) q[3];
sx q[3];
rz(-1.8192889) q[3];
sx q[3];
rz(2.7831912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48034126) q[0];
sx q[0];
rz(-0.77965176) q[0];
sx q[0];
rz(0.10777792) q[0];
rz(-2.5906471) q[1];
sx q[1];
rz(-2.7007553) q[1];
sx q[1];
rz(-1.5239747) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11883987) q[0];
sx q[0];
rz(-1.700239) q[0];
sx q[0];
rz(-1.3368827) q[0];
rz(1.0650159) q[2];
sx q[2];
rz(-2.1518101) q[2];
sx q[2];
rz(0.53321028) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7858582) q[1];
sx q[1];
rz(-1.6140892) q[1];
sx q[1];
rz(-2.0329352) q[1];
rz(-pi) q[2];
rz(-1.518256) q[3];
sx q[3];
rz(-1.0797622) q[3];
sx q[3];
rz(-2.4201916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9109351) q[2];
sx q[2];
rz(-0.50369889) q[2];
sx q[2];
rz(-2.0111734) q[2];
rz(-2.9039827) q[3];
sx q[3];
rz(-0.41046023) q[3];
sx q[3];
rz(2.3835278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1199353) q[0];
sx q[0];
rz(-0.1788685) q[0];
sx q[0];
rz(2.2005431) q[0];
rz(-0.36048105) q[1];
sx q[1];
rz(-1.4070114) q[1];
sx q[1];
rz(1.2233268) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0377285) q[0];
sx q[0];
rz(-2.1489804) q[0];
sx q[0];
rz(1.8072771) q[0];
rz(-pi) q[1];
rz(-2.1603196) q[2];
sx q[2];
rz(-2.2726577) q[2];
sx q[2];
rz(-1.6035994) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4810923) q[1];
sx q[1];
rz(-1.5005439) q[1];
sx q[1];
rz(-2.7226647) q[1];
rz(-pi) q[2];
rz(2.9452865) q[3];
sx q[3];
rz(-1.0428016) q[3];
sx q[3];
rz(1.4083901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.90784043) q[2];
sx q[2];
rz(-1.3774104) q[2];
sx q[2];
rz(-0.09093786) q[2];
rz(0.20445538) q[3];
sx q[3];
rz(-2.3369868) q[3];
sx q[3];
rz(-1.0819775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37888708) q[0];
sx q[0];
rz(-1.611041) q[0];
sx q[0];
rz(1.8086717) q[0];
rz(1.4270205) q[1];
sx q[1];
rz(-0.49976977) q[1];
sx q[1];
rz(1.5907092) q[1];
rz(1.0574404) q[2];
sx q[2];
rz(-0.79887894) q[2];
sx q[2];
rz(0.32499921) q[2];
rz(0.21359276) q[3];
sx q[3];
rz(-1.5145258) q[3];
sx q[3];
rz(0.35888844) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
