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
rz(-0.916565) q[1];
sx q[1];
rz(-0.7601127) q[1];
sx q[1];
rz(2.3488933) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9206031) q[0];
sx q[0];
rz(-1.594829) q[0];
sx q[0];
rz(2.241797) q[0];
x q[1];
rz(1.8007457) q[2];
sx q[2];
rz(-1.0385498) q[2];
sx q[2];
rz(-1.6543653) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1945038) q[1];
sx q[1];
rz(-2.6852486) q[1];
sx q[1];
rz(-1.7366921) q[1];
x q[2];
rz(2.3402432) q[3];
sx q[3];
rz(-1.0145463) q[3];
sx q[3];
rz(0.43503161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4684326) q[2];
sx q[2];
rz(-0.94753733) q[2];
sx q[2];
rz(0.49896487) q[2];
rz(-2.7178102) q[3];
sx q[3];
rz(-1.2200004) q[3];
sx q[3];
rz(-0.012714816) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.760261) q[0];
sx q[0];
rz(-1.4562162) q[0];
sx q[0];
rz(-0.99200845) q[0];
rz(2.965062) q[1];
sx q[1];
rz(-0.20543988) q[1];
sx q[1];
rz(-0.36546779) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1680981) q[0];
sx q[0];
rz(-2.0272278) q[0];
sx q[0];
rz(1.2903851) q[0];
rz(-2.2720049) q[2];
sx q[2];
rz(-0.80539942) q[2];
sx q[2];
rz(-0.7676917) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.47479113) q[1];
sx q[1];
rz(-2.431369) q[1];
sx q[1];
rz(-2.6067961) q[1];
x q[2];
rz(-1.4497827) q[3];
sx q[3];
rz(-1.7704718) q[3];
sx q[3];
rz(1.6691085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3084597) q[2];
sx q[2];
rz(-0.59811991) q[2];
sx q[2];
rz(1.0773405) q[2];
rz(-0.16263738) q[3];
sx q[3];
rz(-1.2626303) q[3];
sx q[3];
rz(-1.8796273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4804374) q[0];
sx q[0];
rz(-2.2876331) q[0];
sx q[0];
rz(2.5307632) q[0];
rz(1.707533) q[1];
sx q[1];
rz(-1.3620946) q[1];
sx q[1];
rz(-0.45062137) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6621646) q[0];
sx q[0];
rz(-3.056884) q[0];
sx q[0];
rz(1.1685632) q[0];
rz(-pi) q[1];
rz(0.89183715) q[2];
sx q[2];
rz(-2.6747181) q[2];
sx q[2];
rz(2.4613949) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0011117) q[1];
sx q[1];
rz(-1.8794267) q[1];
sx q[1];
rz(2.4806697) q[1];
rz(-pi) q[2];
rz(0.30626014) q[3];
sx q[3];
rz(-2.2420886) q[3];
sx q[3];
rz(1.9426395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.11015636) q[2];
sx q[2];
rz(-2.0560196) q[2];
sx q[2];
rz(-1.191656) q[2];
rz(2.2971161) q[3];
sx q[3];
rz(-2.8561487) q[3];
sx q[3];
rz(2.5604131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0464756) q[0];
sx q[0];
rz(-1.3866383) q[0];
sx q[0];
rz(2.0003831) q[0];
rz(-1.2966688) q[1];
sx q[1];
rz(-2.7099687) q[1];
sx q[1];
rz(-0.30532125) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5019708) q[0];
sx q[0];
rz(-1.0143711) q[0];
sx q[0];
rz(-1.9407546) q[0];
rz(-pi) q[1];
rz(1.1496454) q[2];
sx q[2];
rz(-1.1425758) q[2];
sx q[2];
rz(2.2575833) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1034643) q[1];
sx q[1];
rz(-2.3410019) q[1];
sx q[1];
rz(2.994328) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.51586113) q[3];
sx q[3];
rz(-1.1138042) q[3];
sx q[3];
rz(2.7716694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2083464) q[2];
sx q[2];
rz(-1.5758005) q[2];
sx q[2];
rz(2.7944881) q[2];
rz(2.7582205) q[3];
sx q[3];
rz(-0.88630668) q[3];
sx q[3];
rz(-2.0736407) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40889302) q[0];
sx q[0];
rz(-2.2786031) q[0];
sx q[0];
rz(-2.6026759) q[0];
rz(1.4096889) q[1];
sx q[1];
rz(-2.6795645) q[1];
sx q[1];
rz(-2.7358823) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1203128) q[0];
sx q[0];
rz(-0.73182523) q[0];
sx q[0];
rz(-0.41723199) q[0];
rz(-1.0056555) q[2];
sx q[2];
rz(-1.7959776) q[2];
sx q[2];
rz(-1.2852247) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5647445) q[1];
sx q[1];
rz(-1.2291359) q[1];
sx q[1];
rz(0.26775743) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3375072) q[3];
sx q[3];
rz(-0.81987112) q[3];
sx q[3];
rz(-1.6398721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6695909) q[2];
sx q[2];
rz(-0.80299032) q[2];
sx q[2];
rz(-1.6873138) q[2];
rz(-0.83868319) q[3];
sx q[3];
rz(-1.8554153) q[3];
sx q[3];
rz(-1.7680761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8307761) q[0];
sx q[0];
rz(-2.5287479) q[0];
sx q[0];
rz(0.40711656) q[0];
rz(2.282417) q[1];
sx q[1];
rz(-1.5645063) q[1];
sx q[1];
rz(-1.3173332) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.84264) q[0];
sx q[0];
rz(-1.5931804) q[0];
sx q[0];
rz(1.0277261) q[0];
x q[1];
rz(2.4738884) q[2];
sx q[2];
rz(-1.0945079) q[2];
sx q[2];
rz(-0.27031937) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6588905) q[1];
sx q[1];
rz(-1.5679949) q[1];
sx q[1];
rz(-2.3059694) q[1];
rz(-pi) q[2];
rz(0.42767151) q[3];
sx q[3];
rz(-2.1139522) q[3];
sx q[3];
rz(3.1205408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0417827) q[2];
sx q[2];
rz(-0.87974292) q[2];
sx q[2];
rz(2.2463634) q[2];
rz(1.2081395) q[3];
sx q[3];
rz(-2.4619305) q[3];
sx q[3];
rz(-2.4334548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5773425) q[0];
sx q[0];
rz(-1.0471434) q[0];
sx q[0];
rz(-1.1279001) q[0];
rz(-2.9290507) q[1];
sx q[1];
rz(-1.8022585) q[1];
sx q[1];
rz(-1.9099265) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5989482) q[0];
sx q[0];
rz(-2.0751187) q[0];
sx q[0];
rz(1.9796611) q[0];
x q[1];
rz(-2.2052843) q[2];
sx q[2];
rz(-1.0299183) q[2];
sx q[2];
rz(-0.66334954) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.54214189) q[1];
sx q[1];
rz(-2.8891924) q[1];
sx q[1];
rz(-0.42599328) q[1];
rz(-pi) q[2];
x q[2];
rz(1.539174) q[3];
sx q[3];
rz(-2.5368689) q[3];
sx q[3];
rz(-1.5865758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3125375) q[2];
sx q[2];
rz(-0.91784358) q[2];
sx q[2];
rz(2.5404239) q[2];
rz(-2.6230295) q[3];
sx q[3];
rz(-0.31373599) q[3];
sx q[3];
rz(-1.7140088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9940051) q[0];
sx q[0];
rz(-3.1412791) q[0];
sx q[0];
rz(-0.073154733) q[0];
rz(1.0271429) q[1];
sx q[1];
rz(-1.5348397) q[1];
sx q[1];
rz(-2.5491319) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1661467) q[0];
sx q[0];
rz(-1.5243634) q[0];
sx q[0];
rz(0.59068824) q[0];
x q[1];
rz(-0.57581298) q[2];
sx q[2];
rz(-2.0137219) q[2];
sx q[2];
rz(-1.5695614) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5432424) q[1];
sx q[1];
rz(-2.228884) q[1];
sx q[1];
rz(-0.42037085) q[1];
rz(2.153272) q[3];
sx q[3];
rz(-2.6675825) q[3];
sx q[3];
rz(-1.7085027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5804194) q[2];
sx q[2];
rz(-1.6757123) q[2];
sx q[2];
rz(0.44211659) q[2];
rz(-0.90028611) q[3];
sx q[3];
rz(-1.0628275) q[3];
sx q[3];
rz(3.1282848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1135547) q[0];
sx q[0];
rz(-0.91911626) q[0];
sx q[0];
rz(-1.3046718) q[0];
rz(-1.4470944) q[1];
sx q[1];
rz(-1.9686952) q[1];
sx q[1];
rz(-2.5794199) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5484555) q[0];
sx q[0];
rz(-0.89567417) q[0];
sx q[0];
rz(2.0329342) q[0];
rz(1.0400585) q[2];
sx q[2];
rz(-1.4897523) q[2];
sx q[2];
rz(-0.20866742) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8852946) q[1];
sx q[1];
rz(-2.2513299) q[1];
sx q[1];
rz(0.93974944) q[1];
rz(-pi) q[2];
rz(0.4060679) q[3];
sx q[3];
rz(-1.1998402) q[3];
sx q[3];
rz(3.1379116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.28329453) q[2];
sx q[2];
rz(-1.4026182) q[2];
sx q[2];
rz(0.98158681) q[2];
rz(-3.1372519) q[3];
sx q[3];
rz(-2.0948295) q[3];
sx q[3];
rz(2.5481352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8269862) q[0];
sx q[0];
rz(-1.2207323) q[0];
sx q[0];
rz(0.04709588) q[0];
rz(-1.3866407) q[1];
sx q[1];
rz(-1.9853233) q[1];
sx q[1];
rz(0.047853619) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3210147) q[0];
sx q[0];
rz(-2.4860959) q[0];
sx q[0];
rz(1.5232248) q[0];
rz(-pi) q[1];
rz(-2.0132952) q[2];
sx q[2];
rz(-1.5248761) q[2];
sx q[2];
rz(-2.0890582) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0117482) q[1];
sx q[1];
rz(-2.1455892) q[1];
sx q[1];
rz(-2.6227436) q[1];
rz(-pi) q[2];
rz(-0.1286653) q[3];
sx q[3];
rz(-0.4056969) q[3];
sx q[3];
rz(-2.2588244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.64426595) q[2];
sx q[2];
rz(-1.1844119) q[2];
sx q[2];
rz(2.2052374) q[2];
rz(-2.3406773) q[3];
sx q[3];
rz(-0.51910669) q[3];
sx q[3];
rz(0.67805964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8235648) q[0];
sx q[0];
rz(-0.93946539) q[0];
sx q[0];
rz(-2.9325624) q[0];
rz(2.6387852) q[1];
sx q[1];
rz(-1.3192758) q[1];
sx q[1];
rz(-2.1154138) q[1];
rz(0.035632523) q[2];
sx q[2];
rz(-2.590066) q[2];
sx q[2];
rz(-1.8077015) q[2];
rz(-1.1669284) q[3];
sx q[3];
rz(-1.9098017) q[3];
sx q[3];
rz(-3.1168934) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
