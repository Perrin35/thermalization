OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6381792) q[0];
sx q[0];
rz(-1.8621651) q[0];
sx q[0];
rz(2.3655565) q[0];
rz(-2.4322721) q[1];
sx q[1];
rz(-1.5172989) q[1];
sx q[1];
rz(0.59901839) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6544908) q[0];
sx q[0];
rz(-1.693049) q[0];
sx q[0];
rz(2.046665) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.28223306) q[2];
sx q[2];
rz(-0.90105173) q[2];
sx q[2];
rz(1.9053659) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.30556074) q[1];
sx q[1];
rz(-1.9232043) q[1];
sx q[1];
rz(2.6120017) q[1];
rz(-pi) q[2];
rz(-2.3612411) q[3];
sx q[3];
rz(-1.8858741) q[3];
sx q[3];
rz(-1.9596069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1385931) q[2];
sx q[2];
rz(-2.004576) q[2];
sx q[2];
rz(2.9453759) q[2];
rz(-2.0972706) q[3];
sx q[3];
rz(-2.315867) q[3];
sx q[3];
rz(-2.830982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0394548) q[0];
sx q[0];
rz(-0.65334833) q[0];
sx q[0];
rz(2.2838604) q[0];
rz(-2.6013382) q[1];
sx q[1];
rz(-1.0539571) q[1];
sx q[1];
rz(2.3589755) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40561134) q[0];
sx q[0];
rz(-0.95658871) q[0];
sx q[0];
rz(-3.0354397) q[0];
x q[1];
rz(2.1925794) q[2];
sx q[2];
rz(-1.3781386) q[2];
sx q[2];
rz(-0.72177835) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.63943422) q[1];
sx q[1];
rz(-0.63259387) q[1];
sx q[1];
rz(0.32986395) q[1];
rz(-2.2223118) q[3];
sx q[3];
rz(-0.79943919) q[3];
sx q[3];
rz(-0.26200595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.1648078) q[2];
sx q[2];
rz(-0.79443496) q[2];
sx q[2];
rz(0.59107333) q[2];
rz(1.0278541) q[3];
sx q[3];
rz(-2.5058392) q[3];
sx q[3];
rz(-3.0052321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9117821) q[0];
sx q[0];
rz(-0.52017838) q[0];
sx q[0];
rz(-1.0562563) q[0];
rz(2.1504869) q[1];
sx q[1];
rz(-1.3737563) q[1];
sx q[1];
rz(0.80449218) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64824087) q[0];
sx q[0];
rz(-1.8364779) q[0];
sx q[0];
rz(0.65346395) q[0];
x q[1];
rz(2.1995538) q[2];
sx q[2];
rz(-1.647718) q[2];
sx q[2];
rz(-2.6475737) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.1002297) q[1];
sx q[1];
rz(-1.8112665) q[1];
sx q[1];
rz(2.9036456) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3785754) q[3];
sx q[3];
rz(-2.2092735) q[3];
sx q[3];
rz(-2.9137127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6911917) q[2];
sx q[2];
rz(-0.85019008) q[2];
sx q[2];
rz(2.5737393) q[2];
rz(-1.9495226) q[3];
sx q[3];
rz(-2.6046533) q[3];
sx q[3];
rz(0.15538628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6464226) q[0];
sx q[0];
rz(-1.7429054) q[0];
sx q[0];
rz(-1.9544741) q[0];
rz(2.8861956) q[1];
sx q[1];
rz(-1.7385769) q[1];
sx q[1];
rz(2.752221) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.02006836) q[0];
sx q[0];
rz(-2.8830155) q[0];
sx q[0];
rz(-2.5917087) q[0];
rz(-pi) q[1];
x q[1];
rz(1.537048) q[2];
sx q[2];
rz(-0.72740388) q[2];
sx q[2];
rz(0.45798618) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2411449) q[1];
sx q[1];
rz(-0.68640781) q[1];
sx q[1];
rz(1.2696086) q[1];
rz(1.3741174) q[3];
sx q[3];
rz(-1.4405319) q[3];
sx q[3];
rz(-1.9288587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4046459) q[2];
sx q[2];
rz(-0.32361042) q[2];
sx q[2];
rz(-1.3673937) q[2];
rz(-0.5395475) q[3];
sx q[3];
rz(-1.5893785) q[3];
sx q[3];
rz(1.8168943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0668199) q[0];
sx q[0];
rz(-0.19817752) q[0];
sx q[0];
rz(-2.5812126) q[0];
rz(0.21891521) q[1];
sx q[1];
rz(-1.0812662) q[1];
sx q[1];
rz(-2.970649) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31256277) q[0];
sx q[0];
rz(-2.0382904) q[0];
sx q[0];
rz(2.3947021) q[0];
rz(0.60336171) q[2];
sx q[2];
rz(-0.24626479) q[2];
sx q[2];
rz(2.6880996) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7485908) q[1];
sx q[1];
rz(-0.95470631) q[1];
sx q[1];
rz(-0.93008091) q[1];
rz(1.6690977) q[3];
sx q[3];
rz(-2.519033) q[3];
sx q[3];
rz(-3.0873201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.35974744) q[2];
sx q[2];
rz(-1.1283987) q[2];
sx q[2];
rz(0.95853364) q[2];
rz(0.16252276) q[3];
sx q[3];
rz(-2.2992117) q[3];
sx q[3];
rz(1.5490279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.573134) q[0];
sx q[0];
rz(-0.064366654) q[0];
sx q[0];
rz(-2.3934613) q[0];
rz(2.023078) q[1];
sx q[1];
rz(-0.41901127) q[1];
sx q[1];
rz(-0.31707877) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7359809) q[0];
sx q[0];
rz(-1.6389585) q[0];
sx q[0];
rz(-0.24644417) q[0];
rz(-2.9002764) q[2];
sx q[2];
rz(-2.9394627) q[2];
sx q[2];
rz(-1.2723107) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2851781) q[1];
sx q[1];
rz(-0.76907255) q[1];
sx q[1];
rz(-1.4682707) q[1];
rz(1.6104782) q[3];
sx q[3];
rz(-0.79426605) q[3];
sx q[3];
rz(3.119885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.76724425) q[2];
sx q[2];
rz(-2.2129462) q[2];
sx q[2];
rz(1.7283776) q[2];
rz(-0.080538571) q[3];
sx q[3];
rz(-1.3486515) q[3];
sx q[3];
rz(2.9197689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-0.11248511) q[0];
sx q[0];
rz(-1.8928098) q[0];
sx q[0];
rz(-1.8674194) q[0];
rz(-1.796465) q[1];
sx q[1];
rz(-2.3990264) q[1];
sx q[1];
rz(-1.8003731) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91914908) q[0];
sx q[0];
rz(-2.0608927) q[0];
sx q[0];
rz(-2.5045913) q[0];
rz(-pi) q[1];
rz(2.6244782) q[2];
sx q[2];
rz(-0.54143751) q[2];
sx q[2];
rz(1.1351897) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.53645027) q[1];
sx q[1];
rz(-1.7243694) q[1];
sx q[1];
rz(-1.3773514) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.92067155) q[3];
sx q[3];
rz(-1.4772282) q[3];
sx q[3];
rz(-0.9182932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.93366569) q[2];
sx q[2];
rz(-1.9542481) q[2];
sx q[2];
rz(0.46736091) q[2];
rz(-2.3815239) q[3];
sx q[3];
rz(-1.4773388) q[3];
sx q[3];
rz(1.2490341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46045983) q[0];
sx q[0];
rz(-2.12119) q[0];
sx q[0];
rz(-1.6946174) q[0];
rz(-0.32577062) q[1];
sx q[1];
rz(-2.9278946) q[1];
sx q[1];
rz(1.5209341) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2163154) q[0];
sx q[0];
rz(-1.0641353) q[0];
sx q[0];
rz(1.8761047) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.09163945) q[2];
sx q[2];
rz(-1.5704201) q[2];
sx q[2];
rz(2.3918652) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3586761) q[1];
sx q[1];
rz(-2.7442928) q[1];
sx q[1];
rz(-2.0372169) q[1];
rz(-pi) q[2];
rz(2.2748442) q[3];
sx q[3];
rz(-1.731385) q[3];
sx q[3];
rz(-0.15204568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1763566) q[2];
sx q[2];
rz(-0.74702817) q[2];
sx q[2];
rz(0.61526862) q[2];
rz(-1.8687013) q[3];
sx q[3];
rz(-2.8764909) q[3];
sx q[3];
rz(-2.7191775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0203005) q[0];
sx q[0];
rz(-1.2739807) q[0];
sx q[0];
rz(2.2913388) q[0];
rz(-1.1212768) q[1];
sx q[1];
rz(-1.774615) q[1];
sx q[1];
rz(2.3698295) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92895011) q[0];
sx q[0];
rz(-1.5947026) q[0];
sx q[0];
rz(-0.8769518) q[0];
rz(-pi) q[1];
x q[1];
rz(0.36649165) q[2];
sx q[2];
rz(-0.19333177) q[2];
sx q[2];
rz(-1.1761348) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.38852019) q[1];
sx q[1];
rz(-2.8917312) q[1];
sx q[1];
rz(-1.3461963) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2466627) q[3];
sx q[3];
rz(-1.3132902) q[3];
sx q[3];
rz(0.83453629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0261592) q[2];
sx q[2];
rz(-1.5449646) q[2];
sx q[2];
rz(-0.8405295) q[2];
rz(3.0606411) q[3];
sx q[3];
rz(-0.5778802) q[3];
sx q[3];
rz(-2.8988083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7131272) q[0];
sx q[0];
rz(-2.9436538) q[0];
sx q[0];
rz(1.1567098) q[0];
rz(-1.0276065) q[1];
sx q[1];
rz(-1.7637858) q[1];
sx q[1];
rz(-0.81046945) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3768328) q[0];
sx q[0];
rz(-0.78928052) q[0];
sx q[0];
rz(-0.30662243) q[0];
rz(-1.9058305) q[2];
sx q[2];
rz(-2.0910237) q[2];
sx q[2];
rz(-0.65641415) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.6582074) q[1];
sx q[1];
rz(-1.2772296) q[1];
sx q[1];
rz(0.29694966) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.23328554) q[3];
sx q[3];
rz(-2.5676227) q[3];
sx q[3];
rz(3.1149394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.79601866) q[2];
sx q[2];
rz(-1.8859665) q[2];
sx q[2];
rz(-1.5239117) q[2];
rz(2.0412622) q[3];
sx q[3];
rz(-1.1805781) q[3];
sx q[3];
rz(-2.9617917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0236459) q[0];
sx q[0];
rz(-1.7272341) q[0];
sx q[0];
rz(-0.9247307) q[0];
rz(-3.0991411) q[1];
sx q[1];
rz(-1.1141384) q[1];
sx q[1];
rz(1.3420807) q[1];
rz(-3.0855865) q[2];
sx q[2];
rz(-1.1650892) q[2];
sx q[2];
rz(2.7964609) q[2];
rz(-2.6841738) q[3];
sx q[3];
rz(-1.6030967) q[3];
sx q[3];
rz(0.23489192) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
