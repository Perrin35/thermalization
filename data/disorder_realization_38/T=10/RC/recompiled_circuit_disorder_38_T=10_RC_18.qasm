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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.881641) q[0];
sx q[0];
rz(-1.5437484) q[0];
sx q[0];
rz(-1.3511488) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7031541) q[2];
sx q[2];
rz(-2.0619832) q[2];
sx q[2];
rz(0.05471281) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.3783592) q[1];
sx q[1];
rz(-0.86508703) q[1];
sx q[1];
rz(-2.6325429) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.95991858) q[3];
sx q[3];
rz(-1.5023408) q[3];
sx q[3];
rz(1.6144891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7258519) q[2];
sx q[2];
rz(-0.44885138) q[2];
sx q[2];
rz(-0.63981167) q[2];
rz(2.2885627) q[3];
sx q[3];
rz(-0.60522389) q[3];
sx q[3];
rz(2.7602592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-1.2663651) q[0];
sx q[0];
rz(-1.1667644) q[0];
sx q[0];
rz(2.8711328) q[0];
rz(-0.71331435) q[1];
sx q[1];
rz(-2.1062873) q[1];
sx q[1];
rz(1.6289904) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.068802527) q[0];
sx q[0];
rz(-1.2501984) q[0];
sx q[0];
rz(2.1945303) q[0];
rz(1.5741882) q[2];
sx q[2];
rz(-1.2870103) q[2];
sx q[2];
rz(-1.1018745) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1158893) q[1];
sx q[1];
rz(-1.70417) q[1];
sx q[1];
rz(0.89110903) q[1];
rz(-pi) q[2];
rz(-2.2643331) q[3];
sx q[3];
rz(-2.7765397) q[3];
sx q[3];
rz(-0.045543268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2018532) q[2];
sx q[2];
rz(-0.24213232) q[2];
sx q[2];
rz(-0.80336037) q[2];
rz(-2.0837636) q[3];
sx q[3];
rz(-1.4927031) q[3];
sx q[3];
rz(-0.025618205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-1.2597044) q[0];
sx q[0];
rz(-2.0443679) q[0];
sx q[0];
rz(2.846068) q[0];
rz(2.9064536) q[1];
sx q[1];
rz(-1.7087015) q[1];
sx q[1];
rz(-0.74584109) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94507664) q[0];
sx q[0];
rz(-1.3498422) q[0];
sx q[0];
rz(3.0483079) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8638641) q[2];
sx q[2];
rz(-2.5979497) q[2];
sx q[2];
rz(-0.03253983) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1675889) q[1];
sx q[1];
rz(-2.0473192) q[1];
sx q[1];
rz(-1.2567026) q[1];
rz(-pi) q[2];
rz(-2.6607473) q[3];
sx q[3];
rz(-2.6391609) q[3];
sx q[3];
rz(-1.06711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0535584) q[2];
sx q[2];
rz(-3.1159248) q[2];
sx q[2];
rz(2.4528465) q[2];
rz(-3.0912494) q[3];
sx q[3];
rz(-2.2294932) q[3];
sx q[3];
rz(1.6200199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
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
rz(2.922309) q[0];
sx q[0];
rz(-2.121121) q[0];
sx q[0];
rz(3.0396089) q[0];
rz(0.12022262) q[1];
sx q[1];
rz(-0.52821237) q[1];
sx q[1];
rz(-0.27329683) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1647987) q[0];
sx q[0];
rz(-1.8935888) q[0];
sx q[0];
rz(0.093683634) q[0];
rz(-pi) q[1];
rz(-0.69584537) q[2];
sx q[2];
rz(-0.77866422) q[2];
sx q[2];
rz(-0.7960745) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8400152) q[1];
sx q[1];
rz(-1.3372278) q[1];
sx q[1];
rz(-0.15109269) q[1];
x q[2];
rz(-1.7498383) q[3];
sx q[3];
rz(-2.2125707) q[3];
sx q[3];
rz(-2.2620576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.79167241) q[2];
sx q[2];
rz(-0.73799729) q[2];
sx q[2];
rz(-0.50393528) q[2];
rz(-0.079581633) q[3];
sx q[3];
rz(-1.9888398) q[3];
sx q[3];
rz(2.7664405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85916096) q[0];
sx q[0];
rz(-1.0520042) q[0];
sx q[0];
rz(-3.058847) q[0];
rz(0.67963183) q[1];
sx q[1];
rz(-1.6487164) q[1];
sx q[1];
rz(-2.1544429) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8866681) q[0];
sx q[0];
rz(-0.75037557) q[0];
sx q[0];
rz(2.1819941) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5613902) q[2];
sx q[2];
rz(-0.9270037) q[2];
sx q[2];
rz(-0.26088342) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5226749) q[1];
sx q[1];
rz(-0.72026157) q[1];
sx q[1];
rz(-0.48536761) q[1];
rz(1.1779285) q[3];
sx q[3];
rz(-1.53189) q[3];
sx q[3];
rz(-1.2272569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2510898) q[2];
sx q[2];
rz(-0.7901929) q[2];
sx q[2];
rz(-0.21128543) q[2];
rz(-2.7206897) q[3];
sx q[3];
rz(-0.55287164) q[3];
sx q[3];
rz(3.0781854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6565276) q[0];
sx q[0];
rz(-1.2542897) q[0];
sx q[0];
rz(0.791839) q[0];
rz(0.99545288) q[1];
sx q[1];
rz(-0.96770006) q[1];
sx q[1];
rz(-1.8621559) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50514454) q[0];
sx q[0];
rz(-0.2821869) q[0];
sx q[0];
rz(1.6237153) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.40462599) q[2];
sx q[2];
rz(-0.64794108) q[2];
sx q[2];
rz(2.4857156) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2020859) q[1];
sx q[1];
rz(-1.0491976) q[1];
sx q[1];
rz(3.107198) q[1];
x q[2];
rz(-3.0002954) q[3];
sx q[3];
rz(-1.7501037) q[3];
sx q[3];
rz(0.80207223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.25097686) q[2];
sx q[2];
rz(-1.0303409) q[2];
sx q[2];
rz(-2.3708169) q[2];
rz(1.4701014) q[3];
sx q[3];
rz(-2.7272868) q[3];
sx q[3];
rz(-2.9582086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32271785) q[0];
sx q[0];
rz(-2.4920431) q[0];
sx q[0];
rz(3.0859257) q[0];
rz(-2.9260013) q[1];
sx q[1];
rz(-2.3781653) q[1];
sx q[1];
rz(-3.1380222) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5033747) q[0];
sx q[0];
rz(-0.57045454) q[0];
sx q[0];
rz(2.4128782) q[0];
rz(-0.41530208) q[2];
sx q[2];
rz(-2.3447782) q[2];
sx q[2];
rz(2.08564) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.057295416) q[1];
sx q[1];
rz(-2.1665386) q[1];
sx q[1];
rz(-0.36235313) q[1];
rz(-pi) q[2];
rz(3.0655541) q[3];
sx q[3];
rz(-1.3730611) q[3];
sx q[3];
rz(2.7969489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.59721649) q[2];
sx q[2];
rz(-0.95149136) q[2];
sx q[2];
rz(-0.92010951) q[2];
rz(-0.19872935) q[3];
sx q[3];
rz(-1.9072429) q[3];
sx q[3];
rz(0.41771093) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51628095) q[0];
sx q[0];
rz(-1.6315062) q[0];
sx q[0];
rz(-1.4177119) q[0];
rz(0.40813804) q[1];
sx q[1];
rz(-1.0393655) q[1];
sx q[1];
rz(0.67869854) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1269826) q[0];
sx q[0];
rz(-1.8982366) q[0];
sx q[0];
rz(1.2922657) q[0];
rz(-pi) q[1];
rz(-2.6762814) q[2];
sx q[2];
rz(-0.88854549) q[2];
sx q[2];
rz(2.0899783) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.24819599) q[1];
sx q[1];
rz(-1.0475323) q[1];
sx q[1];
rz(0.072695331) q[1];
rz(-pi) q[2];
rz(2.7160866) q[3];
sx q[3];
rz(-1.7153499) q[3];
sx q[3];
rz(1.0021462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9644908) q[2];
sx q[2];
rz(-2.047838) q[2];
sx q[2];
rz(-0.91782451) q[2];
rz(-1.142189) q[3];
sx q[3];
rz(-0.85421383) q[3];
sx q[3];
rz(-2.2043622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98638242) q[0];
sx q[0];
rz(-2.4676403) q[0];
sx q[0];
rz(0.25979364) q[0];
rz(2.4329176) q[1];
sx q[1];
rz(-0.27888137) q[1];
sx q[1];
rz(-2.6146467) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5567112) q[0];
sx q[0];
rz(-1.7877794) q[0];
sx q[0];
rz(2.7816539) q[0];
rz(-pi) q[1];
rz(1.9138463) q[2];
sx q[2];
rz(-2.4545049) q[2];
sx q[2];
rz(-0.81817852) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.37104169) q[1];
sx q[1];
rz(-0.41401225) q[1];
sx q[1];
rz(0.71055926) q[1];
rz(-pi) q[2];
rz(-1.8606887) q[3];
sx q[3];
rz(-1.999951) q[3];
sx q[3];
rz(-3.1115301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8156585) q[2];
sx q[2];
rz(-2.2951173) q[2];
sx q[2];
rz(-0.79088598) q[2];
rz(-0.38665006) q[3];
sx q[3];
rz(-1.0145885) q[3];
sx q[3];
rz(-2.8044243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40795046) q[0];
sx q[0];
rz(-0.17245094) q[0];
sx q[0];
rz(-2.1561484) q[0];
rz(0.5685637) q[1];
sx q[1];
rz(-2.0303346) q[1];
sx q[1];
rz(-2.7808166) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8737313) q[0];
sx q[0];
rz(-2.972313) q[0];
sx q[0];
rz(1.7196359) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.785727) q[2];
sx q[2];
rz(-1.4533918) q[2];
sx q[2];
rz(1.4431151) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.39721397) q[1];
sx q[1];
rz(-2.386464) q[1];
sx q[1];
rz(0.20882864) q[1];
rz(1.4233227) q[3];
sx q[3];
rz(-1.4031938) q[3];
sx q[3];
rz(0.22351219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0976022) q[2];
sx q[2];
rz(-0.68949914) q[2];
sx q[2];
rz(-0.94341755) q[2];
rz(-2.5676981) q[3];
sx q[3];
rz(-0.4807764) q[3];
sx q[3];
rz(1.7452314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(1.7815331) q[1];
sx q[1];
rz(-2.139745) q[1];
sx q[1];
rz(2.3812961) q[1];
rz(-3.0439539) q[2];
sx q[2];
rz(-1.7751481) q[2];
sx q[2];
rz(-2.912599) q[2];
rz(-2.2684569) q[3];
sx q[3];
rz(-0.94948873) q[3];
sx q[3];
rz(-0.82617847) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
