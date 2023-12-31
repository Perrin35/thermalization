OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0302304) q[0];
sx q[0];
rz(5.6279866) q[0];
sx q[0];
rz(10.401166) q[0];
rz(-0.916565) q[1];
sx q[1];
rz(-0.7601127) q[1];
sx q[1];
rz(-0.7926994) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22098955) q[0];
sx q[0];
rz(-1.5467637) q[0];
sx q[0];
rz(-0.8997957) q[0];
rz(-pi) q[1];
rz(1.8007457) q[2];
sx q[2];
rz(-1.0385498) q[2];
sx q[2];
rz(1.4872273) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.9470889) q[1];
sx q[1];
rz(-0.45634404) q[1];
sx q[1];
rz(-1.7366921) q[1];
rz(2.3000852) q[3];
sx q[3];
rz(-2.2268647) q[3];
sx q[3];
rz(-2.5049202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.67316002) q[2];
sx q[2];
rz(-2.1940553) q[2];
sx q[2];
rz(-0.49896487) q[2];
rz(-0.4237825) q[3];
sx q[3];
rz(-1.2200004) q[3];
sx q[3];
rz(-3.1288778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.760261) q[0];
sx q[0];
rz(-1.4562162) q[0];
sx q[0];
rz(0.99200845) q[0];
rz(0.17653067) q[1];
sx q[1];
rz(-2.9361528) q[1];
sx q[1];
rz(2.7761249) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8651543) q[0];
sx q[0];
rz(-1.8218453) q[0];
sx q[0];
rz(2.6692078) q[0];
rz(-pi) q[1];
rz(-0.89895504) q[2];
sx q[2];
rz(-1.0869173) q[2];
sx q[2];
rz(1.8091786) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6234399) q[1];
sx q[1];
rz(-1.9095416) q[1];
sx q[1];
rz(2.5046196) q[1];
rz(2.6037381) q[3];
sx q[3];
rz(-0.23306498) q[3];
sx q[3];
rz(-2.2190998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3084597) q[2];
sx q[2];
rz(-2.5434727) q[2];
sx q[2];
rz(1.0773405) q[2];
rz(0.16263738) q[3];
sx q[3];
rz(-1.2626303) q[3];
sx q[3];
rz(-1.2619654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4804374) q[0];
sx q[0];
rz(-0.85395956) q[0];
sx q[0];
rz(-0.61082947) q[0];
rz(-1.707533) q[1];
sx q[1];
rz(-1.779498) q[1];
sx q[1];
rz(-0.45062137) q[1];
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
rz(-1.9446104) q[2];
sx q[2];
rz(-1.2842368) q[2];
sx q[2];
rz(2.8754004) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0011117) q[1];
sx q[1];
rz(-1.2621659) q[1];
sx q[1];
rz(-2.4806697) q[1];
rz(1.2080473) q[3];
sx q[3];
rz(-0.72788531) q[3];
sx q[3];
rz(-2.4129652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0314363) q[2];
sx q[2];
rz(-2.0560196) q[2];
sx q[2];
rz(1.9499367) q[2];
rz(-2.2971161) q[3];
sx q[3];
rz(-0.28544393) q[3];
sx q[3];
rz(-0.58117956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.095117005) q[0];
sx q[0];
rz(-1.3866383) q[0];
sx q[0];
rz(2.0003831) q[0];
rz(-1.2966688) q[1];
sx q[1];
rz(-2.7099687) q[1];
sx q[1];
rz(2.8362714) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13320623) q[0];
sx q[0];
rz(-1.8828694) q[0];
sx q[0];
rz(2.5532789) q[0];
rz(0.46378739) q[2];
sx q[2];
rz(-1.1897435) q[2];
sx q[2];
rz(-2.2708937) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0381283) q[1];
sx q[1];
rz(-0.80059073) q[1];
sx q[1];
rz(0.14726463) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6257315) q[3];
sx q[3];
rz(-1.1138042) q[3];
sx q[3];
rz(-0.36992321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2083464) q[2];
sx q[2];
rz(-1.5758005) q[2];
sx q[2];
rz(2.7944881) q[2];
rz(2.7582205) q[3];
sx q[3];
rz(-2.255286) q[3];
sx q[3];
rz(-1.0679519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7326996) q[0];
sx q[0];
rz(-0.86298958) q[0];
sx q[0];
rz(-2.6026759) q[0];
rz(-1.4096889) q[1];
sx q[1];
rz(-2.6795645) q[1];
sx q[1];
rz(2.7358823) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51605326) q[0];
sx q[0];
rz(-0.91359432) q[0];
sx q[0];
rz(1.2217193) q[0];
rz(-pi) q[1];
x q[1];
rz(1.97498) q[2];
sx q[2];
rz(-0.6037854) q[2];
sx q[2];
rz(0.052978901) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5647445) q[1];
sx q[1];
rz(-1.2291359) q[1];
sx q[1];
rz(2.8738352) q[1];
rz(-pi) q[2];
rz(1.3375072) q[3];
sx q[3];
rz(-2.3217215) q[3];
sx q[3];
rz(-1.5017205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6695909) q[2];
sx q[2];
rz(-2.3386023) q[2];
sx q[2];
rz(-1.6873138) q[2];
rz(-2.3029095) q[3];
sx q[3];
rz(-1.2861774) q[3];
sx q[3];
rz(-1.7680761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31081653) q[0];
sx q[0];
rz(-2.5287479) q[0];
sx q[0];
rz(0.40711656) q[0];
rz(2.282417) q[1];
sx q[1];
rz(-1.5645063) q[1];
sx q[1];
rz(-1.3173332) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.84264) q[0];
sx q[0];
rz(-1.5931804) q[0];
sx q[0];
rz(-2.1138666) q[0];
rz(0.69465722) q[2];
sx q[2];
rz(-2.3431871) q[2];
sx q[2];
rz(-2.3677804) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6588905) q[1];
sx q[1];
rz(-1.5735978) q[1];
sx q[1];
rz(0.83562327) q[1];
x q[2];
rz(0.42767151) q[3];
sx q[3];
rz(-2.1139522) q[3];
sx q[3];
rz(3.1205408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0417827) q[2];
sx q[2];
rz(-2.2618497) q[2];
sx q[2];
rz(-0.89522925) q[2];
rz(-1.9334531) q[3];
sx q[3];
rz(-0.67966214) q[3];
sx q[3];
rz(-0.7081379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5642501) q[0];
sx q[0];
rz(-1.0471434) q[0];
sx q[0];
rz(1.1279001) q[0];
rz(2.9290507) q[1];
sx q[1];
rz(-1.8022585) q[1];
sx q[1];
rz(1.9099265) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3198277) q[0];
sx q[0];
rz(-1.9263096) q[0];
sx q[0];
rz(-0.54152821) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2052843) q[2];
sx q[2];
rz(-1.0299183) q[2];
sx q[2];
rz(-2.4782431) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.98037887) q[1];
sx q[1];
rz(-1.8002138) q[1];
sx q[1];
rz(-1.4646261) q[1];
rz(-2.1752862) q[3];
sx q[3];
rz(-1.55282) q[3];
sx q[3];
rz(0.041796587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3125375) q[2];
sx q[2];
rz(-2.2237491) q[2];
sx q[2];
rz(2.5404239) q[2];
rz(-0.51856315) q[3];
sx q[3];
rz(-2.8278567) q[3];
sx q[3];
rz(-1.7140088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1475875) q[0];
sx q[0];
rz(-3.1412791) q[0];
sx q[0];
rz(-3.0684379) q[0];
rz(2.1144497) q[1];
sx q[1];
rz(-1.5348397) q[1];
sx q[1];
rz(2.5491319) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1661467) q[0];
sx q[0];
rz(-1.5243634) q[0];
sx q[0];
rz(0.59068824) q[0];
rz(0.7166491) q[2];
sx q[2];
rz(-2.4307494) q[2];
sx q[2];
rz(2.5568642) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.91208273) q[1];
sx q[1];
rz(-0.76369897) q[1];
sx q[1];
rz(2.0565226) q[1];
x q[2];
rz(2.8665364) q[3];
sx q[3];
rz(-1.1797136) q[3];
sx q[3];
rz(1.0712136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5804194) q[2];
sx q[2];
rz(-1.6757123) q[2];
sx q[2];
rz(2.6994761) q[2];
rz(-2.2413065) q[3];
sx q[3];
rz(-1.0628275) q[3];
sx q[3];
rz(-3.1282848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(3.1135547) q[0];
sx q[0];
rz(-2.2224764) q[0];
sx q[0];
rz(1.8369209) q[0];
rz(-1.6944983) q[1];
sx q[1];
rz(-1.1728975) q[1];
sx q[1];
rz(-2.5794199) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5484555) q[0];
sx q[0];
rz(-0.89567417) q[0];
sx q[0];
rz(-1.1086585) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7299037) q[2];
sx q[2];
rz(-2.6052887) q[2];
sx q[2];
rz(1.2250587) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8852946) q[1];
sx q[1];
rz(-0.89026272) q[1];
sx q[1];
rz(-2.2018432) q[1];
rz(-2.7355248) q[3];
sx q[3];
rz(-1.1998402) q[3];
sx q[3];
rz(3.1379116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.28329453) q[2];
sx q[2];
rz(-1.7389745) q[2];
sx q[2];
rz(-0.98158681) q[2];
rz(0.0043407241) q[3];
sx q[3];
rz(-2.0948295) q[3];
sx q[3];
rz(-0.59345746) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3146064) q[0];
sx q[0];
rz(-1.2207323) q[0];
sx q[0];
rz(3.0944968) q[0];
rz(1.3866407) q[1];
sx q[1];
rz(-1.9853233) q[1];
sx q[1];
rz(-0.047853619) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3536516) q[0];
sx q[0];
rz(-1.5997868) q[0];
sx q[0];
rz(-2.2257462) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.677703) q[2];
sx q[2];
rz(-0.44471834) q[2];
sx q[2];
rz(2.5267548) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0117482) q[1];
sx q[1];
rz(-0.99600345) q[1];
sx q[1];
rz(2.6227436) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.1286653) q[3];
sx q[3];
rz(-0.4056969) q[3];
sx q[3];
rz(0.88276827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4973267) q[2];
sx q[2];
rz(-1.9571807) q[2];
sx q[2];
rz(2.2052374) q[2];
rz(0.80091536) q[3];
sx q[3];
rz(-2.622486) q[3];
sx q[3];
rz(2.463533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8235648) q[0];
sx q[0];
rz(-0.93946539) q[0];
sx q[0];
rz(-2.9325624) q[0];
rz(-2.6387852) q[1];
sx q[1];
rz(-1.8223169) q[1];
sx q[1];
rz(1.0261789) q[1];
rz(1.5488831) q[2];
sx q[2];
rz(-2.1219325) q[2];
sx q[2];
rz(1.3757201) q[2];
rz(-2.3021163) q[3];
sx q[3];
rz(-0.52121938) q[3];
sx q[3];
rz(-2.2073707) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
