OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.71198553) q[0];
sx q[0];
rz(3.5482121) q[0];
sx q[0];
rz(9.1756048) q[0];
rz(-0.063440032) q[1];
sx q[1];
rz(4.1133147) q[1];
sx q[1];
rz(9.9749554) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8367856) q[0];
sx q[0];
rz(-1.7903622) q[0];
sx q[0];
rz(-3.1138793) q[0];
x q[1];
rz(2.2416441) q[2];
sx q[2];
rz(-0.64621011) q[2];
sx q[2];
rz(-2.3044569) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7632335) q[1];
sx q[1];
rz(-0.86508703) q[1];
sx q[1];
rz(0.50904973) q[1];
x q[2];
rz(1.689765) q[3];
sx q[3];
rz(-0.61421466) q[3];
sx q[3];
rz(-3.0005232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.41574079) q[2];
sx q[2];
rz(-0.44885138) q[2];
sx q[2];
rz(2.501781) q[2];
rz(2.2885627) q[3];
sx q[3];
rz(-2.5363688) q[3];
sx q[3];
rz(0.38133347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2663651) q[0];
sx q[0];
rz(-1.1667644) q[0];
sx q[0];
rz(-0.27045989) q[0];
rz(-2.4282783) q[1];
sx q[1];
rz(-1.0353054) q[1];
sx q[1];
rz(1.6289904) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4166116) q[0];
sx q[0];
rz(-0.98326251) q[0];
sx q[0];
rz(-0.38831098) q[0];
rz(-pi) q[1];
rz(3.1299635) q[2];
sx q[2];
rz(-0.28380576) q[2];
sx q[2];
rz(-1.1139882) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1158893) q[1];
sx q[1];
rz(-1.4374226) q[1];
sx q[1];
rz(-0.89110903) q[1];
x q[2];
rz(-1.8566425) q[3];
sx q[3];
rz(-1.3405521) q[3];
sx q[3];
rz(2.2765991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9397395) q[2];
sx q[2];
rz(-0.24213232) q[2];
sx q[2];
rz(-0.80336037) q[2];
rz(2.0837636) q[3];
sx q[3];
rz(-1.4927031) q[3];
sx q[3];
rz(-3.1159744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2597044) q[0];
sx q[0];
rz(-2.0443679) q[0];
sx q[0];
rz(-0.29552466) q[0];
rz(-0.23513901) q[1];
sx q[1];
rz(-1.7087015) q[1];
sx q[1];
rz(2.3957516) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.495372) q[0];
sx q[0];
rz(-1.4797858) q[0];
sx q[0];
rz(-1.3489086) q[0];
rz(1.7350115) q[2];
sx q[2];
rz(-2.0914372) q[2];
sx q[2];
rz(2.8525713) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1675889) q[1];
sx q[1];
rz(-1.0942735) q[1];
sx q[1];
rz(-1.88489) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6882719) q[3];
sx q[3];
rz(-1.3461777) q[3];
sx q[3];
rz(-2.2090705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0535584) q[2];
sx q[2];
rz(-3.1159248) q[2];
sx q[2];
rz(2.4528465) q[2];
rz(0.050343242) q[3];
sx q[3];
rz(-0.91209948) q[3];
sx q[3];
rz(1.5215727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.922309) q[0];
sx q[0];
rz(-1.0204717) q[0];
sx q[0];
rz(-0.10198378) q[0];
rz(-3.02137) q[1];
sx q[1];
rz(-0.52821237) q[1];
sx q[1];
rz(2.8682958) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1647987) q[0];
sx q[0];
rz(-1.8935888) q[0];
sx q[0];
rz(0.093683634) q[0];
rz(-1.0068514) q[2];
sx q[2];
rz(-1.0014921) q[2];
sx q[2];
rz(-3.0727) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8371493) q[1];
sx q[1];
rz(-1.4238365) q[1];
sx q[1];
rz(-1.3346346) q[1];
rz(-pi) q[2];
rz(-0.64951879) q[3];
sx q[3];
rz(-1.7139385) q[3];
sx q[3];
rz(0.79917819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3499202) q[2];
sx q[2];
rz(-0.73799729) q[2];
sx q[2];
rz(2.6376574) q[2];
rz(-3.062011) q[3];
sx q[3];
rz(-1.1527529) q[3];
sx q[3];
rz(2.7664405) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2824317) q[0];
sx q[0];
rz(-1.0520042) q[0];
sx q[0];
rz(-0.082745634) q[0];
rz(-2.4619608) q[1];
sx q[1];
rz(-1.4928763) q[1];
sx q[1];
rz(2.1544429) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3521096) q[0];
sx q[0];
rz(-1.1687359) q[0];
sx q[0];
rz(-0.91870086) q[0];
x q[1];
rz(-0.012533112) q[2];
sx q[2];
rz(-2.4977411) q[2];
sx q[2];
rz(2.8963793) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.6189177) q[1];
sx q[1];
rz(-2.4213311) q[1];
sx q[1];
rz(2.656225) q[1];
rz(0.042111245) q[3];
sx q[3];
rz(-1.1782421) q[3];
sx q[3];
rz(0.3596572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8905028) q[2];
sx q[2];
rz(-2.3513998) q[2];
sx q[2];
rz(-2.9303072) q[2];
rz(0.42090297) q[3];
sx q[3];
rz(-2.588721) q[3];
sx q[3];
rz(0.063407272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6565276) q[0];
sx q[0];
rz(-1.2542897) q[0];
sx q[0];
rz(0.791839) q[0];
rz(-0.99545288) q[1];
sx q[1];
rz(-0.96770006) q[1];
sx q[1];
rz(-1.2794367) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6364481) q[0];
sx q[0];
rz(-0.2821869) q[0];
sx q[0];
rz(-1.6237153) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7369667) q[2];
sx q[2];
rz(-0.64794108) q[2];
sx q[2];
rz(-0.65587703) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.8705604) q[1];
sx q[1];
rz(-2.6189657) q[1];
sx q[1];
rz(-1.51103) q[1];
rz(-pi) q[2];
x q[2];
rz(0.1412973) q[3];
sx q[3];
rz(-1.7501037) q[3];
sx q[3];
rz(0.80207223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8906158) q[2];
sx q[2];
rz(-2.1112517) q[2];
sx q[2];
rz(2.3708169) q[2];
rz(1.4701014) q[3];
sx q[3];
rz(-2.7272868) q[3];
sx q[3];
rz(-2.9582086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32271785) q[0];
sx q[0];
rz(-0.6495496) q[0];
sx q[0];
rz(3.0859257) q[0];
rz(-0.21559134) q[1];
sx q[1];
rz(-0.76342738) q[1];
sx q[1];
rz(0.0035704426) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9650426) q[0];
sx q[0];
rz(-1.1561484) q[0];
sx q[0];
rz(-1.167017) q[0];
rz(-pi) q[1];
rz(0.75242075) q[2];
sx q[2];
rz(-1.8634897) q[2];
sx q[2];
rz(2.9257286) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.057295416) q[1];
sx q[1];
rz(-0.97505403) q[1];
sx q[1];
rz(2.7792395) q[1];
rz(-pi) q[2];
rz(3.0655541) q[3];
sx q[3];
rz(-1.3730611) q[3];
sx q[3];
rz(-0.34464371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5443762) q[2];
sx q[2];
rz(-0.95149136) q[2];
sx q[2];
rz(-2.2214831) q[2];
rz(-0.19872935) q[3];
sx q[3];
rz(-1.2343497) q[3];
sx q[3];
rz(2.7238817) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6253117) q[0];
sx q[0];
rz(-1.5100864) q[0];
sx q[0];
rz(1.7238808) q[0];
rz(-2.7334546) q[1];
sx q[1];
rz(-2.1022271) q[1];
sx q[1];
rz(-0.67869854) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6060651) q[0];
sx q[0];
rz(-1.8341656) q[0];
sx q[0];
rz(2.8019964) q[0];
rz(-2.6762814) q[2];
sx q[2];
rz(-2.2530472) q[2];
sx q[2];
rz(-2.0899783) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.855367) q[1];
sx q[1];
rz(-1.6337506) q[1];
sx q[1];
rz(2.0952058) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7292761) q[3];
sx q[3];
rz(-1.9915808) q[3];
sx q[3];
rz(0.63384038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9644908) q[2];
sx q[2];
rz(-2.047838) q[2];
sx q[2];
rz(0.91782451) q[2];
rz(1.142189) q[3];
sx q[3];
rz(-2.2873788) q[3];
sx q[3];
rz(0.93723047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98638242) q[0];
sx q[0];
rz(-2.4676403) q[0];
sx q[0];
rz(2.881799) q[0];
rz(0.70867509) q[1];
sx q[1];
rz(-2.8627113) q[1];
sx q[1];
rz(-2.6146467) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5848815) q[0];
sx q[0];
rz(-1.7877794) q[0];
sx q[0];
rz(0.35993872) q[0];
rz(-2.2286345) q[2];
sx q[2];
rz(-1.7857988) q[2];
sx q[2];
rz(-1.022033) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2745797) q[1];
sx q[1];
rz(-1.305294) q[1];
sx q[1];
rz(-0.32151476) q[1];
rz(1.8606887) q[3];
sx q[3];
rz(-1.1416417) q[3];
sx q[3];
rz(0.030062519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8156585) q[2];
sx q[2];
rz(-0.84647536) q[2];
sx q[2];
rz(2.3507067) q[2];
rz(-0.38665006) q[3];
sx q[3];
rz(-1.0145885) q[3];
sx q[3];
rz(-2.8044243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-2.7336422) q[0];
sx q[0];
rz(-2.9691417) q[0];
sx q[0];
rz(-2.1561484) q[0];
rz(0.5685637) q[1];
sx q[1];
rz(-1.111258) q[1];
sx q[1];
rz(2.7808166) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0246968) q[0];
sx q[0];
rz(-1.403406) q[0];
sx q[0];
rz(-0.025339729) q[0];
x q[1];
rz(1.0656409) q[2];
sx q[2];
rz(-2.897122) q[2];
sx q[2];
rz(0.62015647) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.39721397) q[1];
sx q[1];
rz(-2.386464) q[1];
sx q[1];
rz(-2.932764) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.71513128) q[3];
sx q[3];
rz(-2.9188041) q[3];
sx q[3];
rz(-2.1905394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0439904) q[2];
sx q[2];
rz(-0.68949914) q[2];
sx q[2];
rz(2.1981751) q[2];
rz(-2.5676981) q[3];
sx q[3];
rz(-0.4807764) q[3];
sx q[3];
rz(-1.3963612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5671134) q[0];
sx q[0];
rz(-1.6945101) q[0];
sx q[0];
rz(2.2816818) q[0];
rz(-1.7815331) q[1];
sx q[1];
rz(-1.0018476) q[1];
sx q[1];
rz(-0.76029653) q[1];
rz(-3.0439539) q[2];
sx q[2];
rz(-1.7751481) q[2];
sx q[2];
rz(-2.912599) q[2];
rz(-0.73137024) q[3];
sx q[3];
rz(-2.2435355) q[3];
sx q[3];
rz(0.13767195) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
