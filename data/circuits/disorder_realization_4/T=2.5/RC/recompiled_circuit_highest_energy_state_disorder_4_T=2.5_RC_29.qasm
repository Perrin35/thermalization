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
rz(0.20224686) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7305827) q[0];
sx q[0];
rz(-0.29313335) q[0];
sx q[0];
rz(2.3877451) q[0];
x q[1];
rz(1.0482685) q[2];
sx q[2];
rz(-2.0179581) q[2];
sx q[2];
rz(1.9856404) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.80965197) q[1];
sx q[1];
rz(-2.1686848) q[1];
sx q[1];
rz(0.34418587) q[1];
x q[2];
rz(1.1047045) q[3];
sx q[3];
rz(-1.6045609) q[3];
sx q[3];
rz(-2.1347783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.54174417) q[2];
sx q[2];
rz(-2.2969963) q[2];
sx q[2];
rz(-2.2005626) q[2];
rz(2.8372676) q[3];
sx q[3];
rz(-0.86186886) q[3];
sx q[3];
rz(-2.8743675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32352725) q[0];
sx q[0];
rz(-0.37779385) q[0];
sx q[0];
rz(-2.3622808) q[0];
rz(-1.3620954) q[1];
sx q[1];
rz(-0.39924386) q[1];
sx q[1];
rz(1.13387) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2529566) q[0];
sx q[0];
rz(-1.1373113) q[0];
sx q[0];
rz(-3.033328) q[0];
rz(-pi) q[1];
rz(0.46102779) q[2];
sx q[2];
rz(-3.010058) q[2];
sx q[2];
rz(-1.6361089) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.081959978) q[1];
sx q[1];
rz(-1.4631008) q[1];
sx q[1];
rz(-1.8376985) q[1];
rz(-pi) q[2];
rz(1.5725582) q[3];
sx q[3];
rz(-0.50927466) q[3];
sx q[3];
rz(2.5363896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3900688) q[0];
sx q[0];
rz(-0.5539493) q[0];
sx q[0];
rz(-0.73177904) q[0];
rz(2.7732908) q[1];
sx q[1];
rz(-2.383547) q[1];
sx q[1];
rz(-0.36645737) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5817422) q[0];
sx q[0];
rz(-1.7479715) q[0];
sx q[0];
rz(0.26818256) q[0];
rz(-pi) q[1];
rz(2.3907876) q[2];
sx q[2];
rz(-2.8158203) q[2];
sx q[2];
rz(0.70394403) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7190105) q[1];
sx q[1];
rz(-1.1478462) q[1];
sx q[1];
rz(-0.99653901) q[1];
x q[2];
rz(0.59753527) q[3];
sx q[3];
rz(-2.1794277) q[3];
sx q[3];
rz(3.1239005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1825819) q[2];
sx q[2];
rz(-2.3704447) q[2];
sx q[2];
rz(2.2491573) q[2];
rz(2.6871032) q[3];
sx q[3];
rz(-0.82388866) q[3];
sx q[3];
rz(-2.9096326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3733805) q[0];
sx q[0];
rz(-0.1425655) q[0];
sx q[0];
rz(-2.7349803) q[0];
rz(2.3136102) q[1];
sx q[1];
rz(-1.4031289) q[1];
sx q[1];
rz(-0.83555317) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0428187) q[0];
sx q[0];
rz(-1.5313498) q[0];
sx q[0];
rz(1.7516319) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9172952) q[2];
sx q[2];
rz(-1.76148) q[2];
sx q[2];
rz(0.71587038) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4728188) q[1];
sx q[1];
rz(-1.7175279) q[1];
sx q[1];
rz(-0.43763486) q[1];
rz(-1.8074513) q[3];
sx q[3];
rz(-0.39517212) q[3];
sx q[3];
rz(-2.446536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.39997175) q[2];
sx q[2];
rz(-0.30210945) q[2];
sx q[2];
rz(-0.92140222) q[2];
rz(1.3678) q[3];
sx q[3];
rz(-0.96145815) q[3];
sx q[3];
rz(3.0000946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2530186) q[0];
sx q[0];
rz(-2.044401) q[0];
sx q[0];
rz(-2.9827523) q[0];
rz(1.5171492) q[1];
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
rz(0.76046645) q[0];
sx q[0];
rz(-1.5814578) q[0];
sx q[0];
rz(-3.1233379) q[0];
x q[1];
rz(1.0091802) q[2];
sx q[2];
rz(-1.7744641) q[2];
sx q[2];
rz(-0.18850279) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.161474) q[1];
sx q[1];
rz(-1.9494281) q[1];
sx q[1];
rz(2.9254854) q[1];
rz(1.8908126) q[3];
sx q[3];
rz(-2.5502) q[3];
sx q[3];
rz(-1.0726014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5627731) q[2];
sx q[2];
rz(-2.4476738) q[2];
sx q[2];
rz(2.5895183) q[2];
rz(2.3479346) q[3];
sx q[3];
rz(-0.59760439) q[3];
sx q[3];
rz(-2.1551267) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3119222) q[0];
sx q[0];
rz(-0.86450082) q[0];
sx q[0];
rz(-0.70575869) q[0];
rz(-0.13748473) q[1];
sx q[1];
rz(-0.64422137) q[1];
sx q[1];
rz(0.71298832) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55248439) q[0];
sx q[0];
rz(-0.45408598) q[0];
sx q[0];
rz(2.5986791) q[0];
rz(2.4586468) q[2];
sx q[2];
rz(-2.5120398) q[2];
sx q[2];
rz(0.94972875) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9116152) q[1];
sx q[1];
rz(-1.4739081) q[1];
sx q[1];
rz(1.1712892) q[1];
rz(-0.1162545) q[3];
sx q[3];
rz(-1.2796806) q[3];
sx q[3];
rz(-1.2450041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9299499) q[2];
sx q[2];
rz(-0.86923081) q[2];
sx q[2];
rz(0.28406528) q[2];
rz(-3.0197213) q[3];
sx q[3];
rz(-0.28212306) q[3];
sx q[3];
rz(1.4819283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26509869) q[0];
sx q[0];
rz(-2.5229186) q[0];
sx q[0];
rz(-0.70736831) q[0];
rz(2.7012198) q[1];
sx q[1];
rz(-0.71763867) q[1];
sx q[1];
rz(-1.7792938) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.52584) q[0];
sx q[0];
rz(-0.37459117) q[0];
sx q[0];
rz(0.84217854) q[0];
rz(-pi) q[1];
rz(-2.7168305) q[2];
sx q[2];
rz(-0.889689) q[2];
sx q[2];
rz(2.8515138) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.41304663) q[1];
sx q[1];
rz(-0.73592454) q[1];
sx q[1];
rz(2.3152183) q[1];
rz(-pi) q[2];
rz(1.4555638) q[3];
sx q[3];
rz(-1.1172022) q[3];
sx q[3];
rz(-2.1510002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.97544396) q[2];
sx q[2];
rz(-0.41836172) q[2];
sx q[2];
rz(0.56900209) q[2];
rz(-2.1042018) q[3];
sx q[3];
rz(-2.2834957) q[3];
sx q[3];
rz(-0.4666127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5910832) q[0];
sx q[0];
rz(-2.789848) q[0];
sx q[0];
rz(1.9367223) q[0];
rz(0.51271802) q[1];
sx q[1];
rz(-0.7690438) q[1];
sx q[1];
rz(2.9606294) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7962387) q[0];
sx q[0];
rz(-1.1974044) q[0];
sx q[0];
rz(2.7143584) q[0];
rz(-pi) q[1];
rz(2.7920807) q[2];
sx q[2];
rz(-1.4415359) q[2];
sx q[2];
rz(1.845128) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0389591) q[1];
sx q[1];
rz(-0.60188434) q[1];
sx q[1];
rz(-1.9485056) q[1];
rz(-pi) q[2];
rz(-0.58684565) q[3];
sx q[3];
rz(-0.13258697) q[3];
sx q[3];
rz(2.6574893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5391431) q[2];
sx q[2];
rz(-2.4931543) q[2];
sx q[2];
rz(-3.0446206) q[2];
rz(2.5868331) q[3];
sx q[3];
rz(-1.3223038) q[3];
sx q[3];
rz(2.7831912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48034126) q[0];
sx q[0];
rz(-0.77965176) q[0];
sx q[0];
rz(3.0338147) q[0];
rz(2.5906471) q[1];
sx q[1];
rz(-2.7007553) q[1];
sx q[1];
rz(1.5239747) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1861098) q[0];
sx q[0];
rz(-2.874827) q[0];
sx q[0];
rz(-1.0590932) q[0];
rz(-pi) q[1];
rz(-2.4977106) q[2];
sx q[2];
rz(-1.9877627) q[2];
sx q[2];
rz(-1.3326926) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.904976) q[1];
sx q[1];
rz(-1.1091241) q[1];
sx q[1];
rz(-3.0932337) q[1];
rz(-0.49160853) q[3];
sx q[3];
rz(-1.6171241) q[3];
sx q[3];
rz(-2.2674048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.23065755) q[2];
sx q[2];
rz(-2.6378938) q[2];
sx q[2];
rz(-2.0111734) q[2];
rz(-0.23760992) q[3];
sx q[3];
rz(-2.7311324) q[3];
sx q[3];
rz(2.3835278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0216574) q[0];
sx q[0];
rz(-0.1788685) q[0];
sx q[0];
rz(-0.94104952) q[0];
rz(-2.7811116) q[1];
sx q[1];
rz(-1.7345813) q[1];
sx q[1];
rz(1.2233268) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0377285) q[0];
sx q[0];
rz(-0.99261221) q[0];
sx q[0];
rz(1.3343156) q[0];
rz(-pi) q[1];
rz(2.5599096) q[2];
sx q[2];
rz(-2.2584728) q[2];
sx q[2];
rz(2.4066642) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.66050038) q[1];
sx q[1];
rz(-1.5005439) q[1];
sx q[1];
rz(-0.41892799) q[1];
rz(-pi) q[2];
rz(-1.2480601) q[3];
sx q[3];
rz(-2.5815423) q[3];
sx q[3];
rz(-1.7843475) q[3];
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
rz(2.0596152) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37888708) q[0];
sx q[0];
rz(-1.611041) q[0];
sx q[0];
rz(1.8086717) q[0];
rz(1.7145722) q[1];
sx q[1];
rz(-2.6418229) q[1];
sx q[1];
rz(-1.5508834) q[1];
rz(-2.0841523) q[2];
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
