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
rz(2.1652048) q[0];
rz(2.2250277) q[1];
sx q[1];
rz(-2.38148) q[1];
sx q[1];
rz(-2.3488933) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22098955) q[0];
sx q[0];
rz(-1.5467637) q[0];
sx q[0];
rz(2.241797) q[0];
x q[1];
rz(0.36926271) q[2];
sx q[2];
rz(-2.5662176) q[2];
sx q[2];
rz(-2.0865666) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.9470889) q[1];
sx q[1];
rz(-0.45634404) q[1];
sx q[1];
rz(1.4049005) q[1];
rz(0.80134942) q[3];
sx q[3];
rz(-1.0145463) q[3];
sx q[3];
rz(2.706561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.67316002) q[2];
sx q[2];
rz(-0.94753733) q[2];
sx q[2];
rz(0.49896487) q[2];
rz(-0.4237825) q[3];
sx q[3];
rz(-1.9215923) q[3];
sx q[3];
rz(-0.012714816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38133165) q[0];
sx q[0];
rz(-1.4562162) q[0];
sx q[0];
rz(0.99200845) q[0];
rz(0.17653067) q[1];
sx q[1];
rz(-2.9361528) q[1];
sx q[1];
rz(2.7761249) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7468837) q[0];
sx q[0];
rz(-0.53046662) q[0];
sx q[0];
rz(0.51325004) q[0];
rz(-0.86958779) q[2];
sx q[2];
rz(-2.3361932) q[2];
sx q[2];
rz(-0.7676917) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.47479113) q[1];
sx q[1];
rz(-0.7102237) q[1];
sx q[1];
rz(2.6067961) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6037381) q[3];
sx q[3];
rz(-2.9085277) q[3];
sx q[3];
rz(-0.92249289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3084597) q[2];
sx q[2];
rz(-2.5434727) q[2];
sx q[2];
rz(-1.0773405) q[2];
rz(-0.16263738) q[3];
sx q[3];
rz(-1.8789623) q[3];
sx q[3];
rz(-1.2619654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
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
rz(0.66115528) q[0];
sx q[0];
rz(-0.85395956) q[0];
sx q[0];
rz(-2.5307632) q[0];
rz(-1.707533) q[1];
sx q[1];
rz(-1.3620946) q[1];
sx q[1];
rz(0.45062137) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47942802) q[0];
sx q[0];
rz(-0.084708609) q[0];
sx q[0];
rz(-1.1685632) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1969823) q[2];
sx q[2];
rz(-1.2842368) q[2];
sx q[2];
rz(-2.8754004) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.9420942) q[1];
sx q[1];
rz(-0.71951413) q[1];
sx q[1];
rz(0.47902963) q[1];
x q[2];
rz(-0.30626014) q[3];
sx q[3];
rz(-2.2420886) q[3];
sx q[3];
rz(-1.9426395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.11015636) q[2];
sx q[2];
rz(-2.0560196) q[2];
sx q[2];
rz(-1.191656) q[2];
rz(0.84447652) q[3];
sx q[3];
rz(-0.28544393) q[3];
sx q[3];
rz(-0.58117956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(0.095117005) q[0];
sx q[0];
rz(-1.3866383) q[0];
sx q[0];
rz(-1.1412096) q[0];
rz(1.2966688) q[1];
sx q[1];
rz(-2.7099687) q[1];
sx q[1];
rz(0.30532125) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0083864) q[0];
sx q[0];
rz(-1.2587233) q[0];
sx q[0];
rz(-0.58831373) q[0];
x q[1];
rz(1.1496454) q[2];
sx q[2];
rz(-1.1425758) q[2];
sx q[2];
rz(2.2575833) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2480337) q[1];
sx q[1];
rz(-0.78130022) q[1];
sx q[1];
rz(1.4206738) q[1];
rz(2.6257315) q[3];
sx q[3];
rz(-2.0277884) q[3];
sx q[3];
rz(-2.7716694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9332463) q[2];
sx q[2];
rz(-1.5657921) q[2];
sx q[2];
rz(-2.7944881) q[2];
rz(-2.7582205) q[3];
sx q[3];
rz(-2.255286) q[3];
sx q[3];
rz(-2.0736407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7326996) q[0];
sx q[0];
rz(-2.2786031) q[0];
sx q[0];
rz(-2.6026759) q[0];
rz(1.7319038) q[1];
sx q[1];
rz(-0.46202818) q[1];
sx q[1];
rz(0.40571037) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1203128) q[0];
sx q[0];
rz(-0.73182523) q[0];
sx q[0];
rz(-0.41723199) q[0];
rz(2.876725) q[2];
sx q[2];
rz(-2.1200074) q[2];
sx q[2];
rz(-2.7153646) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5768482) q[1];
sx q[1];
rz(-1.9124568) q[1];
sx q[1];
rz(2.8738352) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3769576) q[3];
sx q[3];
rz(-1.4009762) q[3];
sx q[3];
rz(-0.09165435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6695909) q[2];
sx q[2];
rz(-0.80299032) q[2];
sx q[2];
rz(1.6873138) q[2];
rz(0.83868319) q[3];
sx q[3];
rz(-1.2861774) q[3];
sx q[3];
rz(-1.7680761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8307761) q[0];
sx q[0];
rz(-0.61284471) q[0];
sx q[0];
rz(2.7344761) q[0];
rz(2.282417) q[1];
sx q[1];
rz(-1.5770864) q[1];
sx q[1];
rz(1.3173332) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8562397) q[0];
sx q[0];
rz(-1.0278773) q[0];
sx q[0];
rz(3.1154484) q[0];
x q[1];
rz(-2.1520734) q[2];
sx q[2];
rz(-0.98810722) q[2];
sx q[2];
rz(-1.6473824) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0560318) q[1];
sx q[1];
rz(-0.83562682) q[1];
sx q[1];
rz(-3.1378156) q[1];
rz(2.156593) q[3];
sx q[3];
rz(-1.9337774) q[3];
sx q[3];
rz(1.3604878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.09981) q[2];
sx q[2];
rz(-0.87974292) q[2];
sx q[2];
rz(2.2463634) q[2];
rz(1.9334531) q[3];
sx q[3];
rz(-2.4619305) q[3];
sx q[3];
rz(-0.7081379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5773425) q[0];
sx q[0];
rz(-1.0471434) q[0];
sx q[0];
rz(1.1279001) q[0];
rz(-0.212542) q[1];
sx q[1];
rz(-1.3393341) q[1];
sx q[1];
rz(1.2316661) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8679778) q[0];
sx q[0];
rz(-0.6379188) q[0];
sx q[0];
rz(-0.62423737) q[0];
rz(-pi) q[1];
rz(0.64078757) q[2];
sx q[2];
rz(-1.0377585) q[2];
sx q[2];
rz(-1.8719045) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1612138) q[1];
sx q[1];
rz(-1.3413789) q[1];
sx q[1];
rz(1.6769665) q[1];
rz(1.539174) q[3];
sx q[3];
rz(-2.5368689) q[3];
sx q[3];
rz(-1.5865758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.82905519) q[2];
sx q[2];
rz(-0.91784358) q[2];
sx q[2];
rz(-2.5404239) q[2];
rz(0.51856315) q[3];
sx q[3];
rz(-2.8278567) q[3];
sx q[3];
rz(1.7140088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9940051) q[0];
sx q[0];
rz(-3.1412791) q[0];
sx q[0];
rz(-0.073154733) q[0];
rz(-2.1144497) q[1];
sx q[1];
rz(-1.606753) q[1];
sx q[1];
rz(2.5491319) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6153529) q[0];
sx q[0];
rz(-0.59229367) q[0];
sx q[0];
rz(-0.083239716) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0560889) q[2];
sx q[2];
rz(-1.0564431) q[2];
sx q[2];
rz(-0.27013847) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2295099) q[1];
sx q[1];
rz(-2.3778937) q[1];
sx q[1];
rz(1.0850701) q[1];
rz(-pi) q[2];
rz(0.27505625) q[3];
sx q[3];
rz(-1.1797136) q[3];
sx q[3];
rz(2.070379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5804194) q[2];
sx q[2];
rz(-1.6757123) q[2];
sx q[2];
rz(-2.6994761) q[2];
rz(-0.90028611) q[3];
sx q[3];
rz(-2.0787652) q[3];
sx q[3];
rz(-3.1282848) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1135547) q[0];
sx q[0];
rz(-2.2224764) q[0];
sx q[0];
rz(-1.3046718) q[0];
rz(-1.6944983) q[1];
sx q[1];
rz(-1.1728975) q[1];
sx q[1];
rz(-2.5794199) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5484555) q[0];
sx q[0];
rz(-2.2459185) q[0];
sx q[0];
rz(1.1086585) q[0];
rz(-pi) q[1];
rz(0.093900605) q[2];
sx q[2];
rz(-2.0996089) q[2];
sx q[2];
rz(-1.4096066) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2562981) q[1];
sx q[1];
rz(-0.89026272) q[1];
sx q[1];
rz(-2.2018432) q[1];
rz(-pi) q[2];
x q[2];
rz(0.77769827) q[3];
sx q[3];
rz(-2.5986528) q[3];
sx q[3];
rz(0.87399769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.28329453) q[2];
sx q[2];
rz(-1.4026182) q[2];
sx q[2];
rz(2.1600058) q[2];
rz(-0.0043407241) q[3];
sx q[3];
rz(-1.0467632) q[3];
sx q[3];
rz(-0.59345746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8269862) q[0];
sx q[0];
rz(-1.2207323) q[0];
sx q[0];
rz(-3.0944968) q[0];
rz(-1.3866407) q[1];
sx q[1];
rz(-1.9853233) q[1];
sx q[1];
rz(0.047853619) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78794107) q[0];
sx q[0];
rz(-1.5997868) q[0];
sx q[0];
rz(0.91584648) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1282975) q[2];
sx q[2];
rz(-1.6167165) q[2];
sx q[2];
rz(2.0890582) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1298444) q[1];
sx q[1];
rz(-2.1455892) q[1];
sx q[1];
rz(-2.6227436) q[1];
rz(-pi) q[2];
rz(1.6258532) q[3];
sx q[3];
rz(-1.1686472) q[3];
sx q[3];
rz(1.0226585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.64426595) q[2];
sx q[2];
rz(-1.1844119) q[2];
sx q[2];
rz(-2.2052374) q[2];
rz(0.80091536) q[3];
sx q[3];
rz(-0.51910669) q[3];
sx q[3];
rz(-2.463533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8235648) q[0];
sx q[0];
rz(-2.2021273) q[0];
sx q[0];
rz(0.20903023) q[0];
rz(-2.6387852) q[1];
sx q[1];
rz(-1.8223169) q[1];
sx q[1];
rz(1.0261789) q[1];
rz(-0.035632523) q[2];
sx q[2];
rz(-0.55152668) q[2];
sx q[2];
rz(1.3338911) q[2];
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
