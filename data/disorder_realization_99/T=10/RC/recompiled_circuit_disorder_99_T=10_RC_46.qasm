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
rz(2.38148) q[1];
sx q[1];
rz(10.217477) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3307318) q[0];
sx q[0];
rz(-0.90002492) q[0];
sx q[0];
rz(3.1109111) q[0];
x q[1];
rz(1.8007457) q[2];
sx q[2];
rz(-2.1030428) q[2];
sx q[2];
rz(1.6543653) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.76268643) q[1];
sx q[1];
rz(-1.1211809) q[1];
sx q[1];
rz(-0.080888918) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3000852) q[3];
sx q[3];
rz(-2.2268647) q[3];
sx q[3];
rz(2.5049202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4684326) q[2];
sx q[2];
rz(-2.1940553) q[2];
sx q[2];
rz(0.49896487) q[2];
rz(-2.7178102) q[3];
sx q[3];
rz(-1.9215923) q[3];
sx q[3];
rz(0.012714816) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.760261) q[0];
sx q[0];
rz(-1.4562162) q[0];
sx q[0];
rz(-0.99200845) q[0];
rz(0.17653067) q[1];
sx q[1];
rz(-2.9361528) q[1];
sx q[1];
rz(-0.36546779) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1680981) q[0];
sx q[0];
rz(-2.0272278) q[0];
sx q[0];
rz(-1.8512076) q[0];
rz(2.2720049) q[2];
sx q[2];
rz(-2.3361932) q[2];
sx q[2];
rz(2.373901) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9531627) q[1];
sx q[1];
rz(-0.97524446) q[1];
sx q[1];
rz(-1.1577391) q[1];
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
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.83313292) q[2];
sx q[2];
rz(-2.5434727) q[2];
sx q[2];
rz(2.0642521) q[2];
rz(-0.16263738) q[3];
sx q[3];
rz(-1.2626303) q[3];
sx q[3];
rz(1.2619654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
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
rz(-1.779498) q[1];
sx q[1];
rz(2.6909713) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.492309) q[0];
sx q[0];
rz(-1.603924) q[0];
sx q[0];
rz(1.6487728) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.89183715) q[2];
sx q[2];
rz(-2.6747181) q[2];
sx q[2];
rz(0.6801978) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1994985) q[1];
sx q[1];
rz(-0.71951413) q[1];
sx q[1];
rz(-2.662563) q[1];
x q[2];
rz(-1.2080473) q[3];
sx q[3];
rz(-0.72788531) q[3];
sx q[3];
rz(2.4129652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0314363) q[2];
sx q[2];
rz(-2.0560196) q[2];
sx q[2];
rz(-1.191656) q[2];
rz(-0.84447652) q[3];
sx q[3];
rz(-2.8561487) q[3];
sx q[3];
rz(2.5604131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.095117005) q[0];
sx q[0];
rz(-1.7549544) q[0];
sx q[0];
rz(-1.1412096) q[0];
rz(-1.8449239) q[1];
sx q[1];
rz(-0.43162391) q[1];
sx q[1];
rz(2.8362714) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5019708) q[0];
sx q[0];
rz(-1.0143711) q[0];
sx q[0];
rz(1.200838) q[0];
rz(-pi) q[1];
rz(-1.9919473) q[2];
sx q[2];
rz(-1.1425758) q[2];
sx q[2];
rz(2.2575833) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.42974737) q[1];
sx q[1];
rz(-1.6763121) q[1];
sx q[1];
rz(-0.79515102) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6257315) q[3];
sx q[3];
rz(-2.0277884) q[3];
sx q[3];
rz(-2.7716694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2083464) q[2];
sx q[2];
rz(-1.5657921) q[2];
sx q[2];
rz(-0.34710458) q[2];
rz(0.38337213) q[3];
sx q[3];
rz(-0.88630668) q[3];
sx q[3];
rz(2.0736407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-0.40889302) q[0];
sx q[0];
rz(-2.2786031) q[0];
sx q[0];
rz(0.53891671) q[0];
rz(-1.7319038) q[1];
sx q[1];
rz(-2.6795645) q[1];
sx q[1];
rz(-2.7358823) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51605326) q[0];
sx q[0];
rz(-2.2279983) q[0];
sx q[0];
rz(-1.9198734) q[0];
rz(1.0056555) q[2];
sx q[2];
rz(-1.7959776) q[2];
sx q[2];
rz(1.2852247) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.87863641) q[1];
sx q[1];
rz(-2.7107781) q[1];
sx q[1];
rz(-2.2104435) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3375072) q[3];
sx q[3];
rz(-2.3217215) q[3];
sx q[3];
rz(1.6398721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6695909) q[2];
sx q[2];
rz(-0.80299032) q[2];
sx q[2];
rz(-1.6873138) q[2];
rz(0.83868319) q[3];
sx q[3];
rz(-1.2861774) q[3];
sx q[3];
rz(1.3735166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31081653) q[0];
sx q[0];
rz(-2.5287479) q[0];
sx q[0];
rz(-0.40711656) q[0];
rz(-0.85917568) q[1];
sx q[1];
rz(-1.5645063) q[1];
sx q[1];
rz(-1.3173332) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8562397) q[0];
sx q[0];
rz(-2.1137153) q[0];
sx q[0];
rz(-0.026144233) q[0];
rz(-0.69465722) q[2];
sx q[2];
rz(-2.3431871) q[2];
sx q[2];
rz(2.3677804) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6588905) q[1];
sx q[1];
rz(-1.5735978) q[1];
sx q[1];
rz(2.3059694) q[1];
rz(-2.7139211) q[3];
sx q[3];
rz(-1.0276405) q[3];
sx q[3];
rz(0.021051858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0417827) q[2];
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
x q[2];
rz(pi/2) q[3];
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
rz(1.5642501) q[0];
sx q[0];
rz(-2.0944493) q[0];
sx q[0];
rz(2.0136925) q[0];
rz(-2.9290507) q[1];
sx q[1];
rz(-1.8022585) q[1];
sx q[1];
rz(1.2316661) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5989482) q[0];
sx q[0];
rz(-2.0751187) q[0];
sx q[0];
rz(1.9796611) q[0];
rz(-2.2052843) q[2];
sx q[2];
rz(-1.0299183) q[2];
sx q[2];
rz(-0.66334954) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1612138) q[1];
sx q[1];
rz(-1.8002138) q[1];
sx q[1];
rz(1.4646261) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1197458) q[3];
sx q[3];
rz(-2.1751746) q[3];
sx q[3];
rz(-1.6250087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3125375) q[2];
sx q[2];
rz(-2.2237491) q[2];
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
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1475875) q[0];
sx q[0];
rz(-3.1412791) q[0];
sx q[0];
rz(0.073154733) q[0];
rz(1.0271429) q[1];
sx q[1];
rz(-1.606753) q[1];
sx q[1];
rz(2.5491319) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5151278) q[0];
sx q[0];
rz(-2.160762) q[0];
sx q[0];
rz(1.626684) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.7166491) q[2];
sx q[2];
rz(-0.71084329) q[2];
sx q[2];
rz(-0.58472842) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.29443024) q[1];
sx q[1];
rz(-1.8995598) q[1];
sx q[1];
rz(-0.86818236) q[1];
rz(-pi) q[2];
rz(0.98832066) q[3];
sx q[3];
rz(-0.47401014) q[3];
sx q[3];
rz(-1.7085027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5804194) q[2];
sx q[2];
rz(-1.6757123) q[2];
sx q[2];
rz(-0.44211659) q[2];
rz(-0.90028611) q[3];
sx q[3];
rz(-1.0628275) q[3];
sx q[3];
rz(3.1282848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1135547) q[0];
sx q[0];
rz(-2.2224764) q[0];
sx q[0];
rz(-1.8369209) q[0];
rz(-1.6944983) q[1];
sx q[1];
rz(-1.1728975) q[1];
sx q[1];
rz(-2.5794199) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2660414) q[0];
sx q[0];
rz(-0.79715675) q[0];
sx q[0];
rz(-2.6334727) q[0];
x q[1];
rz(-0.093900605) q[2];
sx q[2];
rz(-2.0996089) q[2];
sx q[2];
rz(-1.7319861) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8852946) q[1];
sx q[1];
rz(-0.89026272) q[1];
sx q[1];
rz(-2.2018432) q[1];
x q[2];
rz(0.77769827) q[3];
sx q[3];
rz(-0.54293984) q[3];
sx q[3];
rz(2.267595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.28329453) q[2];
sx q[2];
rz(-1.4026182) q[2];
sx q[2];
rz(2.1600058) q[2];
rz(-3.1372519) q[3];
sx q[3];
rz(-2.0948295) q[3];
sx q[3];
rz(2.5481352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3146064) q[0];
sx q[0];
rz(-1.9208603) q[0];
sx q[0];
rz(-3.0944968) q[0];
rz(1.754952) q[1];
sx q[1];
rz(-1.9853233) q[1];
sx q[1];
rz(0.047853619) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3536516) q[0];
sx q[0];
rz(-1.5997868) q[0];
sx q[0];
rz(0.91584648) q[0];
rz(1.1282975) q[2];
sx q[2];
rz(-1.6167165) q[2];
sx q[2];
rz(-1.0525345) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.88356) q[1];
sx q[1];
rz(-2.0000534) q[1];
sx q[1];
rz(-0.92991035) q[1];
rz(-1.6258532) q[3];
sx q[3];
rz(-1.1686472) q[3];
sx q[3];
rz(2.1189342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.64426595) q[2];
sx q[2];
rz(-1.9571807) q[2];
sx q[2];
rz(-2.2052374) q[2];
rz(-0.80091536) q[3];
sx q[3];
rz(-2.622486) q[3];
sx q[3];
rz(0.67805964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
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
rz(-1.5927096) q[2];
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
