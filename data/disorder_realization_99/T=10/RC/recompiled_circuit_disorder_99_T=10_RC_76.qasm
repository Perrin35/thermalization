OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1113623) q[0];
sx q[0];
rz(-2.4863939) q[0];
sx q[0];
rz(-0.97638786) q[0];
rz(2.2250277) q[1];
sx q[1];
rz(-2.38148) q[1];
sx q[1];
rz(-2.3488933) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3307318) q[0];
sx q[0];
rz(-0.90002492) q[0];
sx q[0];
rz(3.1109111) q[0];
rz(-1.8007457) q[2];
sx q[2];
rz(-2.1030428) q[2];
sx q[2];
rz(1.4872273) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3789062) q[1];
sx q[1];
rz(-1.1211809) q[1];
sx q[1];
rz(3.0607037) q[1];
x q[2];
rz(2.4281265) q[3];
sx q[3];
rz(-0.93868512) q[3];
sx q[3];
rz(1.5330832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4684326) q[2];
sx q[2];
rz(-0.94753733) q[2];
sx q[2];
rz(-2.6426278) q[2];
rz(-2.7178102) q[3];
sx q[3];
rz(-1.9215923) q[3];
sx q[3];
rz(0.012714816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.760261) q[0];
sx q[0];
rz(-1.6853764) q[0];
sx q[0];
rz(2.1495842) q[0];
rz(2.965062) q[1];
sx q[1];
rz(-2.9361528) q[1];
sx q[1];
rz(0.36546779) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.394709) q[0];
sx q[0];
rz(-0.53046662) q[0];
sx q[0];
rz(0.51325004) q[0];
rz(-pi) q[1];
rz(-0.59132691) q[2];
sx q[2];
rz(-2.1543243) q[2];
sx q[2];
rz(0.11596767) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6668015) q[1];
sx q[1];
rz(-0.7102237) q[1];
sx q[1];
rz(0.53479654) q[1];
rz(-pi) q[2];
rz(2.6037381) q[3];
sx q[3];
rz(-0.23306498) q[3];
sx q[3];
rz(-2.2190998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3084597) q[2];
sx q[2];
rz(-0.59811991) q[2];
sx q[2];
rz(-2.0642521) q[2];
rz(-2.9789553) q[3];
sx q[3];
rz(-1.8789623) q[3];
sx q[3];
rz(-1.8796273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66115528) q[0];
sx q[0];
rz(-0.85395956) q[0];
sx q[0];
rz(-2.5307632) q[0];
rz(1.707533) q[1];
sx q[1];
rz(-1.779498) q[1];
sx q[1];
rz(0.45062137) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0656933) q[0];
sx q[0];
rz(-1.4928627) q[0];
sx q[0];
rz(3.1083641) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9446104) q[2];
sx q[2];
rz(-1.8573559) q[2];
sx q[2];
rz(2.8754004) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.33774099) q[1];
sx q[1];
rz(-2.1954814) q[1];
sx q[1];
rz(-1.9546263) q[1];
rz(-pi) q[2];
rz(-1.9335453) q[3];
sx q[3];
rz(-0.72788531) q[3];
sx q[3];
rz(0.7286275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.11015636) q[2];
sx q[2];
rz(-2.0560196) q[2];
sx q[2];
rz(1.9499367) q[2];
rz(0.84447652) q[3];
sx q[3];
rz(-2.8561487) q[3];
sx q[3];
rz(0.58117956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-3.0464756) q[0];
sx q[0];
rz(-1.3866383) q[0];
sx q[0];
rz(2.0003831) q[0];
rz(-1.8449239) q[1];
sx q[1];
rz(-0.43162391) q[1];
sx q[1];
rz(2.8362714) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1353377) q[0];
sx q[0];
rz(-2.4843785) q[0];
sx q[0];
rz(-0.52657907) q[0];
rz(-1.9919473) q[2];
sx q[2];
rz(-1.9990168) q[2];
sx q[2];
rz(-2.2575833) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.893559) q[1];
sx q[1];
rz(-0.78130022) q[1];
sx q[1];
rz(-1.4206738) q[1];
rz(0.78379811) q[3];
sx q[3];
rz(-0.67516203) q[3];
sx q[3];
rz(-0.53962196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9332463) q[2];
sx q[2];
rz(-1.5758005) q[2];
sx q[2];
rz(0.34710458) q[2];
rz(-2.7582205) q[3];
sx q[3];
rz(-0.88630668) q[3];
sx q[3];
rz(-1.0679519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40889302) q[0];
sx q[0];
rz(-2.2786031) q[0];
sx q[0];
rz(2.6026759) q[0];
rz(-1.7319038) q[1];
sx q[1];
rz(-0.46202818) q[1];
sx q[1];
rz(2.7358823) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8680506) q[0];
sx q[0];
rz(-1.2965856) q[0];
sx q[0];
rz(-0.6875086) q[0];
x q[1];
rz(2.876725) q[2];
sx q[2];
rz(-2.1200074) q[2];
sx q[2];
rz(-2.7153646) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.87863641) q[1];
sx q[1];
rz(-0.43081455) q[1];
sx q[1];
rz(-2.2104435) q[1];
rz(-pi) q[2];
x q[2];
rz(2.898786) q[3];
sx q[3];
rz(-0.77951509) q[3];
sx q[3];
rz(-1.8368343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6695909) q[2];
sx q[2];
rz(-2.3386023) q[2];
sx q[2];
rz(-1.6873138) q[2];
rz(0.83868319) q[3];
sx q[3];
rz(-1.8554153) q[3];
sx q[3];
rz(-1.3735166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8307761) q[0];
sx q[0];
rz(-0.61284471) q[0];
sx q[0];
rz(0.40711656) q[0];
rz(-0.85917568) q[1];
sx q[1];
rz(-1.5645063) q[1];
sx q[1];
rz(-1.3173332) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8562397) q[0];
sx q[0];
rz(-1.0278773) q[0];
sx q[0];
rz(-0.026144233) q[0];
rz(-pi) q[1];
rz(2.4738884) q[2];
sx q[2];
rz(-2.0470847) q[2];
sx q[2];
rz(0.27031937) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0504005) q[1];
sx q[1];
rz(-0.7351774) q[1];
sx q[1];
rz(-1.5749732) q[1];
rz(2.1727354) q[3];
sx q[3];
rz(-2.4638306) q[3];
sx q[3];
rz(2.4399151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.09981) q[2];
sx q[2];
rz(-0.87974292) q[2];
sx q[2];
rz(-0.89522925) q[2];
rz(-1.2081395) q[3];
sx q[3];
rz(-2.4619305) q[3];
sx q[3];
rz(-0.7081379) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5642501) q[0];
sx q[0];
rz(-2.0944493) q[0];
sx q[0];
rz(1.1279001) q[0];
rz(0.212542) q[1];
sx q[1];
rz(-1.8022585) q[1];
sx q[1];
rz(1.2316661) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8679778) q[0];
sx q[0];
rz(-2.5036739) q[0];
sx q[0];
rz(-2.5173553) q[0];
x q[1];
rz(2.5008051) q[2];
sx q[2];
rz(-2.1038342) q[2];
sx q[2];
rz(-1.8719045) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5994508) q[1];
sx q[1];
rz(-2.8891924) q[1];
sx q[1];
rz(-0.42599328) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1197458) q[3];
sx q[3];
rz(-2.1751746) q[3];
sx q[3];
rz(1.5165839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3125375) q[2];
sx q[2];
rz(-2.2237491) q[2];
sx q[2];
rz(-2.5404239) q[2];
rz(-2.6230295) q[3];
sx q[3];
rz(-2.8278567) q[3];
sx q[3];
rz(1.7140088) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9940051) q[0];
sx q[0];
rz(-3.1412791) q[0];
sx q[0];
rz(-3.0684379) q[0];
rz(-1.0271429) q[1];
sx q[1];
rz(-1.5348397) q[1];
sx q[1];
rz(-0.59246078) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62646482) q[0];
sx q[0];
rz(-2.160762) q[0];
sx q[0];
rz(-1.5149087) q[0];
rz(0.7166491) q[2];
sx q[2];
rz(-0.71084329) q[2];
sx q[2];
rz(0.58472842) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5432424) q[1];
sx q[1];
rz(-2.228884) q[1];
sx q[1];
rz(-0.42037085) q[1];
rz(-2.153272) q[3];
sx q[3];
rz(-0.47401014) q[3];
sx q[3];
rz(1.43309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.56117326) q[2];
sx q[2];
rz(-1.6757123) q[2];
sx q[2];
rz(0.44211659) q[2];
rz(2.2413065) q[3];
sx q[3];
rz(-1.0628275) q[3];
sx q[3];
rz(3.1282848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028037926) q[0];
sx q[0];
rz(-0.91911626) q[0];
sx q[0];
rz(-1.8369209) q[0];
rz(1.6944983) q[1];
sx q[1];
rz(-1.9686952) q[1];
sx q[1];
rz(0.56217271) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5931372) q[0];
sx q[0];
rz(-0.89567417) q[0];
sx q[0];
rz(-2.0329342) q[0];
rz(0.093900605) q[2];
sx q[2];
rz(-1.0419838) q[2];
sx q[2];
rz(-1.7319861) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.11645424) q[1];
sx q[1];
rz(-1.0944195) q[1];
sx q[1];
rz(-0.78671793) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9713054) q[3];
sx q[3];
rz(-1.9477961) q[3];
sx q[3];
rz(1.4124944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.28329453) q[2];
sx q[2];
rz(-1.7389745) q[2];
sx q[2];
rz(-2.1600058) q[2];
rz(0.0043407241) q[3];
sx q[3];
rz(-2.0948295) q[3];
sx q[3];
rz(2.5481352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8269862) q[0];
sx q[0];
rz(-1.9208603) q[0];
sx q[0];
rz(-3.0944968) q[0];
rz(1.754952) q[1];
sx q[1];
rz(-1.9853233) q[1];
sx q[1];
rz(-3.093739) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76059607) q[0];
sx q[0];
rz(-2.2254235) q[0];
sx q[0];
rz(-3.1050443) q[0];
rz(2.0132952) q[2];
sx q[2];
rz(-1.6167165) q[2];
sx q[2];
rz(1.0525345) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0117482) q[1];
sx q[1];
rz(-2.1455892) q[1];
sx q[1];
rz(0.5188491) q[1];
rz(-pi) q[2];
rz(0.1286653) q[3];
sx q[3];
rz(-2.7358958) q[3];
sx q[3];
rz(-2.2588244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4973267) q[2];
sx q[2];
rz(-1.9571807) q[2];
sx q[2];
rz(0.93635526) q[2];
rz(-2.3406773) q[3];
sx q[3];
rz(-2.622486) q[3];
sx q[3];
rz(-0.67805964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3180278) q[0];
sx q[0];
rz(-2.2021273) q[0];
sx q[0];
rz(0.20903023) q[0];
rz(-0.50280747) q[1];
sx q[1];
rz(-1.3192758) q[1];
sx q[1];
rz(-2.1154138) q[1];
rz(-2.5903493) q[2];
sx q[2];
rz(-1.5894645) q[2];
sx q[2];
rz(-0.20655256) q[2];
rz(2.3021163) q[3];
sx q[3];
rz(-2.6203733) q[3];
sx q[3];
rz(0.93422191) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
