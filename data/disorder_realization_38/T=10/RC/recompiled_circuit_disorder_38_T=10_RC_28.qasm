OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4296071) q[0];
sx q[0];
rz(-0.40661943) q[0];
sx q[0];
rz(-2.8924195) q[0];
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
rz(0.30480706) q[0];
sx q[0];
rz(-1.3512304) q[0];
sx q[0];
rz(3.1138793) q[0];
rz(-2.1044188) q[2];
sx q[2];
rz(-1.1871157) q[2];
sx q[2];
rz(-1.8431611) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5397415) q[1];
sx q[1];
rz(-1.1907693) q[1];
sx q[1];
rz(2.3439581) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1816741) q[3];
sx q[3];
rz(-1.5023408) q[3];
sx q[3];
rz(-1.5271036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7258519) q[2];
sx q[2];
rz(-0.44885138) q[2];
sx q[2];
rz(-2.501781) q[2];
rz(-0.85302991) q[3];
sx q[3];
rz(-2.5363688) q[3];
sx q[3];
rz(0.38133347) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2663651) q[0];
sx q[0];
rz(-1.1667644) q[0];
sx q[0];
rz(-2.8711328) q[0];
rz(2.4282783) q[1];
sx q[1];
rz(-1.0353054) q[1];
sx q[1];
rz(-1.6289904) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0727901) q[0];
sx q[0];
rz(-1.2501984) q[0];
sx q[0];
rz(-0.94706236) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8578051) q[2];
sx q[2];
rz(-1.5675401) q[2];
sx q[2];
rz(2.6736205) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1158893) q[1];
sx q[1];
rz(-1.4374226) q[1];
sx q[1];
rz(-2.2504836) q[1];
x q[2];
rz(2.9019722) q[3];
sx q[3];
rz(-1.2926971) q[3];
sx q[3];
rz(2.3688189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9397395) q[2];
sx q[2];
rz(-2.8994603) q[2];
sx q[2];
rz(0.80336037) q[2];
rz(1.057829) q[3];
sx q[3];
rz(-1.4927031) q[3];
sx q[3];
rz(-0.025618205) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8818883) q[0];
sx q[0];
rz(-2.0443679) q[0];
sx q[0];
rz(0.29552466) q[0];
rz(0.23513901) q[1];
sx q[1];
rz(-1.4328911) q[1];
sx q[1];
rz(-0.74584109) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.495372) q[0];
sx q[0];
rz(-1.4797858) q[0];
sx q[0];
rz(-1.792684) q[0];
rz(0.27772851) q[2];
sx q[2];
rz(-0.54364294) q[2];
sx q[2];
rz(-3.1090528) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5902139) q[1];
sx q[1];
rz(-2.5775902) q[1];
sx q[1];
rz(2.6022634) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3219222) q[3];
sx q[3];
rz(-1.1296774) q[3];
sx q[3];
rz(0.53019023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.088034257) q[2];
sx q[2];
rz(-0.025667889) q[2];
sx q[2];
rz(-0.68874613) q[2];
rz(-0.050343242) q[3];
sx q[3];
rz(-2.2294932) q[3];
sx q[3];
rz(-1.6200199) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.922309) q[0];
sx q[0];
rz(-1.0204717) q[0];
sx q[0];
rz(0.10198378) q[0];
rz(0.12022262) q[1];
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
rz(0.8768339) q[0];
sx q[0];
rz(-0.33565531) q[0];
sx q[0];
rz(-1.8434974) q[0];
x q[1];
rz(-0.69584537) q[2];
sx q[2];
rz(-0.77866422) q[2];
sx q[2];
rz(2.3455182) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8400152) q[1];
sx q[1];
rz(-1.3372278) q[1];
sx q[1];
rz(-2.9905) q[1];
rz(-pi) q[2];
rz(-2.9076505) q[3];
sx q[3];
rz(-2.4787239) q[3];
sx q[3];
rz(-2.5556504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3499202) q[2];
sx q[2];
rz(-2.4035954) q[2];
sx q[2];
rz(2.6376574) q[2];
rz(3.062011) q[3];
sx q[3];
rz(-1.1527529) q[3];
sx q[3];
rz(0.3751522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2824317) q[0];
sx q[0];
rz(-1.0520042) q[0];
sx q[0];
rz(-0.082745634) q[0];
rz(-0.67963183) q[1];
sx q[1];
rz(-1.4928763) q[1];
sx q[1];
rz(0.98714978) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7894831) q[0];
sx q[0];
rz(-1.1687359) q[0];
sx q[0];
rz(-2.2228918) q[0];
rz(0.64381386) q[2];
sx q[2];
rz(-1.5783196) q[2];
sx q[2];
rz(1.8260337) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3293216) q[1];
sx q[1];
rz(-1.2580039) q[1];
sx q[1];
rz(-0.66004628) q[1];
x q[2];
rz(3.0994814) q[3];
sx q[3];
rz(-1.1782421) q[3];
sx q[3];
rz(2.7819355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2510898) q[2];
sx q[2];
rz(-2.3513998) q[2];
sx q[2];
rz(-0.21128543) q[2];
rz(2.7206897) q[3];
sx q[3];
rz(-2.588721) q[3];
sx q[3];
rz(3.0781854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6565276) q[0];
sx q[0];
rz(-1.2542897) q[0];
sx q[0];
rz(-0.791839) q[0];
rz(-2.1461398) q[1];
sx q[1];
rz(-2.1738926) q[1];
sx q[1];
rz(-1.2794367) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5813542) q[0];
sx q[0];
rz(-1.8525774) q[0];
sx q[0];
rz(-3.1262585) q[0];
x q[1];
rz(-2.7369667) q[2];
sx q[2];
rz(-2.4936516) q[2];
sx q[2];
rz(2.4857156) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.8705604) q[1];
sx q[1];
rz(-0.522627) q[1];
sx q[1];
rz(1.51103) q[1];
rz(-pi) q[2];
rz(-1.389723) q[3];
sx q[3];
rz(-1.7098134) q[3];
sx q[3];
rz(-2.3982323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.25097686) q[2];
sx q[2];
rz(-2.1112517) q[2];
sx q[2];
rz(2.3708169) q[2];
rz(-1.4701014) q[3];
sx q[3];
rz(-0.41430587) q[3];
sx q[3];
rz(0.18338403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32271785) q[0];
sx q[0];
rz(-2.4920431) q[0];
sx q[0];
rz(3.0859257) q[0];
rz(0.21559134) q[1];
sx q[1];
rz(-2.3781653) q[1];
sx q[1];
rz(-3.1380222) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5033747) q[0];
sx q[0];
rz(-0.57045454) q[0];
sx q[0];
rz(-0.72871448) q[0];
rz(2.7262906) q[2];
sx q[2];
rz(-0.79681444) q[2];
sx q[2];
rz(1.0559527) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.53686245) q[1];
sx q[1];
rz(-2.4559048) q[1];
sx q[1];
rz(-2.0525949) q[1];
x q[2];
rz(-1.7690897) q[3];
sx q[3];
rz(-1.6453504) q[3];
sx q[3];
rz(-1.9004746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.59721649) q[2];
sx q[2];
rz(-0.95149136) q[2];
sx q[2];
rz(2.2214831) q[2];
rz(-0.19872935) q[3];
sx q[3];
rz(-1.2343497) q[3];
sx q[3];
rz(-0.41771093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-2.6253117) q[0];
sx q[0];
rz(-1.6315062) q[0];
sx q[0];
rz(1.4177119) q[0];
rz(0.40813804) q[1];
sx q[1];
rz(-1.0393655) q[1];
sx q[1];
rz(0.67869854) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6060651) q[0];
sx q[0];
rz(-1.3074271) q[0];
sx q[0];
rz(-2.8019964) q[0];
rz(2.6762814) q[2];
sx q[2];
rz(-2.2530472) q[2];
sx q[2];
rz(-1.0516143) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8933967) q[1];
sx q[1];
rz(-2.0940603) q[1];
sx q[1];
rz(3.0688973) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4123165) q[3];
sx q[3];
rz(-1.9915808) q[3];
sx q[3];
rz(0.63384038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.17710182) q[2];
sx q[2];
rz(-1.0937546) q[2];
sx q[2];
rz(0.91782451) q[2];
rz(-1.9994036) q[3];
sx q[3];
rz(-2.2873788) q[3];
sx q[3];
rz(-2.2043622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98638242) q[0];
sx q[0];
rz(-2.4676403) q[0];
sx q[0];
rz(0.25979364) q[0];
rz(0.70867509) q[1];
sx q[1];
rz(-2.8627113) q[1];
sx q[1];
rz(-2.6146467) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6078867) q[0];
sx q[0];
rz(-0.41782802) q[0];
sx q[0];
rz(-2.582344) q[0];
x q[1];
rz(-0.2692659) q[2];
sx q[2];
rz(-2.2109647) q[2];
sx q[2];
rz(-0.38538853) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.867013) q[1];
sx q[1];
rz(-1.8362987) q[1];
sx q[1];
rz(-2.8200779) q[1];
rz(0.55836375) q[3];
sx q[3];
rz(-2.6287968) q[3];
sx q[3];
rz(-0.59190291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.32593411) q[2];
sx q[2];
rz(-0.84647536) q[2];
sx q[2];
rz(-2.3507067) q[2];
rz(2.7549426) q[3];
sx q[3];
rz(-2.1270042) q[3];
sx q[3];
rz(-0.33716831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40795046) q[0];
sx q[0];
rz(-0.17245094) q[0];
sx q[0];
rz(-0.98544425) q[0];
rz(-2.573029) q[1];
sx q[1];
rz(-2.0303346) q[1];
sx q[1];
rz(-2.7808166) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0246968) q[0];
sx q[0];
rz(-1.7381867) q[0];
sx q[0];
rz(0.025339729) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3558657) q[2];
sx q[2];
rz(-1.4533918) q[2];
sx q[2];
rz(1.6984775) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.39721397) q[1];
sx q[1];
rz(-0.75512868) q[1];
sx q[1];
rz(-2.932764) q[1];
rz(-pi) q[2];
rz(2.9721857) q[3];
sx q[3];
rz(-1.7161887) q[3];
sx q[3];
rz(1.8190847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0439904) q[2];
sx q[2];
rz(-2.4520935) q[2];
sx q[2];
rz(2.1981751) q[2];
rz(-2.5676981) q[3];
sx q[3];
rz(-0.4807764) q[3];
sx q[3];
rz(1.7452314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-1.3654937) q[2];
sx q[2];
rz(-1.4751954) q[2];
sx q[2];
rz(1.8196646) q[2];
rz(2.4102224) q[3];
sx q[3];
rz(-2.2435355) q[3];
sx q[3];
rz(0.13767195) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
