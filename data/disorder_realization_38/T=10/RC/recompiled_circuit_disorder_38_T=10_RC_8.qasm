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
rz(3.0781526) q[1];
sx q[1];
rz(-0.97172207) q[1];
sx q[1];
rz(2.5914153) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8367856) q[0];
sx q[0];
rz(-1.3512304) q[0];
sx q[0];
rz(0.027713393) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1044188) q[2];
sx q[2];
rz(-1.1871157) q[2];
sx q[2];
rz(1.2984315) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7632335) q[1];
sx q[1];
rz(-2.2765056) q[1];
sx q[1];
rz(-0.50904973) q[1];
rz(-pi) q[2];
x q[2];
rz(1.689765) q[3];
sx q[3];
rz(-2.527378) q[3];
sx q[3];
rz(3.0005232) q[3];
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
rz(-0.63981167) q[2];
rz(0.85302991) q[3];
sx q[3];
rz(-2.5363688) q[3];
sx q[3];
rz(2.7602592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.8752276) q[0];
sx q[0];
rz(-1.1667644) q[0];
sx q[0];
rz(-2.8711328) q[0];
rz(0.71331435) q[1];
sx q[1];
rz(-1.0353054) q[1];
sx q[1];
rz(-1.5126022) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4166116) q[0];
sx q[0];
rz(-2.1583301) q[0];
sx q[0];
rz(2.7532817) q[0];
x q[1];
rz(-1.5741882) q[2];
sx q[2];
rz(-1.8545824) q[2];
sx q[2];
rz(-1.1018745) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4334129) q[1];
sx q[1];
rz(-2.4509894) q[1];
sx q[1];
rz(1.7811) q[1];
x q[2];
rz(-0.23962044) q[3];
sx q[3];
rz(-1.2926971) q[3];
sx q[3];
rz(-0.77277377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2018532) q[2];
sx q[2];
rz(-0.24213232) q[2];
sx q[2];
rz(-2.3382323) q[2];
rz(1.057829) q[3];
sx q[3];
rz(-1.4927031) q[3];
sx q[3];
rz(3.1159744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8818883) q[0];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94507664) q[0];
sx q[0];
rz(-1.3498422) q[0];
sx q[0];
rz(-0.093284746) q[0];
rz(-pi) q[1];
rz(-0.27772851) q[2];
sx q[2];
rz(-0.54364294) q[2];
sx q[2];
rz(-0.03253983) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.44887603) q[1];
sx q[1];
rz(-1.2926896) q[1];
sx q[1];
rz(-0.49726185) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3219222) q[3];
sx q[3];
rz(-1.1296774) q[3];
sx q[3];
rz(-2.6114024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0535584) q[2];
sx q[2];
rz(-0.025667889) q[2];
sx q[2];
rz(0.68874613) q[2];
rz(-0.050343242) q[3];
sx q[3];
rz(-0.91209948) q[3];
sx q[3];
rz(1.6200199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21928366) q[0];
sx q[0];
rz(-1.0204717) q[0];
sx q[0];
rz(3.0396089) q[0];
rz(-0.12022262) q[1];
sx q[1];
rz(-2.6133803) q[1];
sx q[1];
rz(2.8682958) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7057987) q[0];
sx q[0];
rz(-1.6596284) q[0];
sx q[0];
rz(-1.8949132) q[0];
rz(0.64812135) q[2];
sx q[2];
rz(-2.0378049) q[2];
sx q[2];
rz(-1.8304706) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.71972388) q[1];
sx q[1];
rz(-0.27742741) q[1];
sx q[1];
rz(-1.006702) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9076505) q[3];
sx q[3];
rz(-0.66286874) q[3];
sx q[3];
rz(0.58594221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.79167241) q[2];
sx q[2];
rz(-2.4035954) q[2];
sx q[2];
rz(0.50393528) q[2];
rz(-0.079581633) q[3];
sx q[3];
rz(-1.9888398) q[3];
sx q[3];
rz(-0.3751522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2824317) q[0];
sx q[0];
rz(-2.0895884) q[0];
sx q[0];
rz(-0.082745634) q[0];
rz(-0.67963183) q[1];
sx q[1];
rz(-1.6487164) q[1];
sx q[1];
rz(-0.98714978) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3521096) q[0];
sx q[0];
rz(-1.9728567) q[0];
sx q[0];
rz(0.91870086) q[0];
rz(-pi) q[1];
x q[1];
rz(0.012533112) q[2];
sx q[2];
rz(-0.64385157) q[2];
sx q[2];
rz(2.8963793) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3293216) q[1];
sx q[1];
rz(-1.8835888) q[1];
sx q[1];
rz(-0.66004628) q[1];
rz(3.0994814) q[3];
sx q[3];
rz(-1.9633506) q[3];
sx q[3];
rz(-2.7819355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2510898) q[2];
sx q[2];
rz(-2.3513998) q[2];
sx q[2];
rz(0.21128543) q[2];
rz(0.42090297) q[3];
sx q[3];
rz(-0.55287164) q[3];
sx q[3];
rz(-0.063407272) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6565276) q[0];
sx q[0];
rz(-1.2542897) q[0];
sx q[0];
rz(-2.3497537) q[0];
rz(2.1461398) q[1];
sx q[1];
rz(-2.1738926) q[1];
sx q[1];
rz(-1.8621559) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56023843) q[0];
sx q[0];
rz(-1.2890153) q[0];
sx q[0];
rz(-0.015334107) q[0];
rz(-pi) q[1];
rz(1.860414) q[2];
sx q[2];
rz(-2.1589303) q[2];
sx q[2];
rz(-2.9786125) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4931603) q[1];
sx q[1];
rz(-1.5409768) q[1];
sx q[1];
rz(-1.0489419) q[1];
rz(-pi) q[2];
rz(-1.7518696) q[3];
sx q[3];
rz(-1.7098134) q[3];
sx q[3];
rz(2.3982323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8188748) q[0];
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
rz(-0.63821793) q[0];
sx q[0];
rz(-2.5711381) q[0];
sx q[0];
rz(0.72871448) q[0];
rz(-pi) q[1];
rz(-2.7262906) q[2];
sx q[2];
rz(-0.79681444) q[2];
sx q[2];
rz(2.08564) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.057295416) q[1];
sx q[1];
rz(-0.97505403) q[1];
sx q[1];
rz(-2.7792395) q[1];
rz(-1.2083863) q[3];
sx q[3];
rz(-0.21167314) q[3];
sx q[3];
rz(0.025312245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5443762) q[2];
sx q[2];
rz(-0.95149136) q[2];
sx q[2];
rz(2.2214831) q[2];
rz(-0.19872935) q[3];
sx q[3];
rz(-1.2343497) q[3];
sx q[3];
rz(2.7238817) q[3];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51628095) q[0];
sx q[0];
rz(-1.5100864) q[0];
sx q[0];
rz(-1.4177119) q[0];
rz(0.40813804) q[1];
sx q[1];
rz(-1.0393655) q[1];
sx q[1];
rz(0.67869854) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1269826) q[0];
sx q[0];
rz(-1.243356) q[0];
sx q[0];
rz(1.8493269) q[0];
rz(-pi) q[1];
rz(-2.6762814) q[2];
sx q[2];
rz(-2.2530472) q[2];
sx q[2];
rz(1.0516143) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2862257) q[1];
sx q[1];
rz(-1.6337506) q[1];
sx q[1];
rz(1.0463868) q[1];
x q[2];
rz(0.42550605) q[3];
sx q[3];
rz(-1.4262428) q[3];
sx q[3];
rz(-2.1394465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.17710182) q[2];
sx q[2];
rz(-2.047838) q[2];
sx q[2];
rz(-0.91782451) q[2];
rz(1.9994036) q[3];
sx q[3];
rz(-0.85421383) q[3];
sx q[3];
rz(0.93723047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98638242) q[0];
sx q[0];
rz(-0.6739524) q[0];
sx q[0];
rz(2.881799) q[0];
rz(2.4329176) q[1];
sx q[1];
rz(-0.27888137) q[1];
sx q[1];
rz(0.52694595) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6078867) q[0];
sx q[0];
rz(-0.41782802) q[0];
sx q[0];
rz(2.582344) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9138463) q[2];
sx q[2];
rz(-0.68708778) q[2];
sx q[2];
rz(-2.3234141) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2745797) q[1];
sx q[1];
rz(-1.8362987) q[1];
sx q[1];
rz(-0.32151476) q[1];
x q[2];
rz(0.44550495) q[3];
sx q[3];
rz(-1.8337436) q[3];
sx q[3];
rz(1.477369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.32593411) q[2];
sx q[2];
rz(-0.84647536) q[2];
sx q[2];
rz(0.79088598) q[2];
rz(-0.38665006) q[3];
sx q[3];
rz(-1.0145885) q[3];
sx q[3];
rz(0.33716831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40795046) q[0];
sx q[0];
rz(-0.17245094) q[0];
sx q[0];
rz(0.98544425) q[0];
rz(0.5685637) q[1];
sx q[1];
rz(-2.0303346) q[1];
sx q[1];
rz(-2.7808166) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4496778) q[0];
sx q[0];
rz(-1.5957818) q[0];
sx q[0];
rz(-1.4033532) q[0];
rz(-pi) q[1];
rz(-1.785727) q[2];
sx q[2];
rz(-1.6882009) q[2];
sx q[2];
rz(-1.4431151) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.68041486) q[1];
sx q[1];
rz(-2.3056681) q[1];
sx q[1];
rz(-1.7635036) q[1];
rz(-pi) q[2];
rz(-0.71513128) q[3];
sx q[3];
rz(-2.9188041) q[3];
sx q[3];
rz(0.95105329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0439904) q[2];
sx q[2];
rz(-0.68949914) q[2];
sx q[2];
rz(0.94341755) q[2];
rz(0.57389456) q[3];
sx q[3];
rz(-2.6608163) q[3];
sx q[3];
rz(1.3963612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5744793) q[0];
sx q[0];
rz(-1.4470826) q[0];
sx q[0];
rz(-0.8599109) q[0];
rz(1.3600596) q[1];
sx q[1];
rz(-1.0018476) q[1];
sx q[1];
rz(-0.76029653) q[1];
rz(-1.1311244) q[2];
sx q[2];
rz(-0.22618539) q[2];
sx q[2];
rz(-2.4629081) q[2];
rz(0.75136649) q[3];
sx q[3];
rz(-1.02117) q[3];
sx q[3];
rz(1.1985967) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];