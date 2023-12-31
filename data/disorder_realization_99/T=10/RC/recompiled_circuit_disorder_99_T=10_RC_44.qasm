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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22098955) q[0];
sx q[0];
rz(-1.5467637) q[0];
sx q[0];
rz(2.241797) q[0];
x q[1];
rz(-1.3408469) q[2];
sx q[2];
rz(-1.0385498) q[2];
sx q[2];
rz(1.4872273) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.76268643) q[1];
sx q[1];
rz(-2.0204117) q[1];
sx q[1];
rz(0.080888918) q[1];
rz(0.84150746) q[3];
sx q[3];
rz(-2.2268647) q[3];
sx q[3];
rz(2.5049202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.67316002) q[2];
sx q[2];
rz(-2.1940553) q[2];
sx q[2];
rz(-2.6426278) q[2];
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
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.760261) q[0];
sx q[0];
rz(-1.4562162) q[0];
sx q[0];
rz(2.1495842) q[0];
rz(0.17653067) q[1];
sx q[1];
rz(-2.9361528) q[1];
sx q[1];
rz(-0.36546779) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.394709) q[0];
sx q[0];
rz(-0.53046662) q[0];
sx q[0];
rz(-0.51325004) q[0];
x q[1];
rz(-2.2720049) q[2];
sx q[2];
rz(-0.80539942) q[2];
sx q[2];
rz(2.373901) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.47479113) q[1];
sx q[1];
rz(-0.7102237) q[1];
sx q[1];
rz(0.53479654) q[1];
x q[2];
rz(-0.5378546) q[3];
sx q[3];
rz(-0.23306498) q[3];
sx q[3];
rz(0.92249289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3084597) q[2];
sx q[2];
rz(-0.59811991) q[2];
sx q[2];
rz(2.0642521) q[2];
rz(0.16263738) q[3];
sx q[3];
rz(-1.2626303) q[3];
sx q[3];
rz(1.8796273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4804374) q[0];
sx q[0];
rz(-0.85395956) q[0];
sx q[0];
rz(0.61082947) q[0];
rz(1.4340596) q[1];
sx q[1];
rz(-1.779498) q[1];
sx q[1];
rz(2.6909713) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47942802) q[0];
sx q[0];
rz(-0.084708609) q[0];
sx q[0];
rz(1.1685632) q[0];
rz(-0.30655105) q[2];
sx q[2];
rz(-1.9286641) q[2];
sx q[2];
rz(1.4150261) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8038517) q[1];
sx q[1];
rz(-0.9461113) q[1];
sx q[1];
rz(1.9546263) q[1];
rz(0.87617989) q[3];
sx q[3];
rz(-1.8091222) q[3];
sx q[3];
rz(0.5660457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0314363) q[2];
sx q[2];
rz(-2.0560196) q[2];
sx q[2];
rz(-1.9499367) q[2];
rz(2.2971161) q[3];
sx q[3];
rz(-0.28544393) q[3];
sx q[3];
rz(-2.5604131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.095117005) q[0];
sx q[0];
rz(-1.7549544) q[0];
sx q[0];
rz(1.1412096) q[0];
rz(1.2966688) q[1];
sx q[1];
rz(-0.43162391) q[1];
sx q[1];
rz(2.8362714) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0083864) q[0];
sx q[0];
rz(-1.2587233) q[0];
sx q[0];
rz(2.5532789) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.46378739) q[2];
sx q[2];
rz(-1.9518491) q[2];
sx q[2];
rz(-2.2708937) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.42974737) q[1];
sx q[1];
rz(-1.6763121) q[1];
sx q[1];
rz(-2.3464416) q[1];
rz(-0.78379811) q[3];
sx q[3];
rz(-0.67516203) q[3];
sx q[3];
rz(0.53962196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2083464) q[2];
sx q[2];
rz(-1.5758005) q[2];
sx q[2];
rz(0.34710458) q[2];
rz(-0.38337213) q[3];
sx q[3];
rz(-0.88630668) q[3];
sx q[3];
rz(1.0679519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40889302) q[0];
sx q[0];
rz(-2.2786031) q[0];
sx q[0];
rz(-0.53891671) q[0];
rz(1.7319038) q[1];
sx q[1];
rz(-2.6795645) q[1];
sx q[1];
rz(-0.40571037) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8680506) q[0];
sx q[0];
rz(-1.8450071) q[0];
sx q[0];
rz(0.6875086) q[0];
rz(2.1359372) q[2];
sx q[2];
rz(-1.3456151) q[2];
sx q[2];
rz(-1.856368) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0438784) q[1];
sx q[1];
rz(-1.8227302) q[1];
sx q[1];
rz(-1.21752) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3769576) q[3];
sx q[3];
rz(-1.7406165) q[3];
sx q[3];
rz(-0.09165435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.47200176) q[2];
sx q[2];
rz(-0.80299032) q[2];
sx q[2];
rz(-1.6873138) q[2];
rz(0.83868319) q[3];
sx q[3];
rz(-1.8554153) q[3];
sx q[3];
rz(1.7680761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8307761) q[0];
sx q[0];
rz(-2.5287479) q[0];
sx q[0];
rz(-2.7344761) q[0];
rz(2.282417) q[1];
sx q[1];
rz(-1.5770864) q[1];
sx q[1];
rz(1.3173332) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8562397) q[0];
sx q[0];
rz(-2.1137153) q[0];
sx q[0];
rz(-0.026144233) q[0];
x q[1];
rz(0.69465722) q[2];
sx q[2];
rz(-0.79840556) q[2];
sx q[2];
rz(2.3677804) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6588905) q[1];
sx q[1];
rz(-1.5735978) q[1];
sx q[1];
rz(-0.83562327) q[1];
rz(-pi) q[2];
rz(2.7139211) q[3];
sx q[3];
rz(-2.1139522) q[3];
sx q[3];
rz(-3.1205408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.09981) q[2];
sx q[2];
rz(-2.2618497) q[2];
sx q[2];
rz(0.89522925) q[2];
rz(-1.2081395) q[3];
sx q[3];
rz(-2.4619305) q[3];
sx q[3];
rz(-0.7081379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(-2.0944493) q[0];
sx q[0];
rz(2.0136925) q[0];
rz(-2.9290507) q[1];
sx q[1];
rz(-1.8022585) q[1];
sx q[1];
rz(-1.9099265) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2736149) q[0];
sx q[0];
rz(-2.5036739) q[0];
sx q[0];
rz(2.5173553) q[0];
rz(-pi) q[1];
x q[1];
rz(0.93630837) q[2];
sx q[2];
rz(-1.0299183) q[2];
sx q[2];
rz(2.4782431) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.54214189) q[1];
sx q[1];
rz(-0.25240024) q[1];
sx q[1];
rz(-2.7155994) q[1];
rz(-pi) q[2];
rz(-0.96630649) q[3];
sx q[3];
rz(-1.55282) q[3];
sx q[3];
rz(-0.041796587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.82905519) q[2];
sx q[2];
rz(-0.91784358) q[2];
sx q[2];
rz(-2.5404239) q[2];
rz(-2.6230295) q[3];
sx q[3];
rz(-0.31373599) q[3];
sx q[3];
rz(1.4275838) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9940051) q[0];
sx q[0];
rz(-0.00031358263) q[0];
sx q[0];
rz(-3.0684379) q[0];
rz(1.0271429) q[1];
sx q[1];
rz(-1.5348397) q[1];
sx q[1];
rz(-2.5491319) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62646482) q[0];
sx q[0];
rz(-2.160762) q[0];
sx q[0];
rz(-1.626684) q[0];
x q[1];
rz(-2.0855037) q[2];
sx q[2];
rz(-2.0851496) q[2];
sx q[2];
rz(-2.8714542) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5432424) q[1];
sx q[1];
rz(-2.228884) q[1];
sx q[1];
rz(-2.7212218) q[1];
rz(-pi) q[2];
x q[2];
rz(2.153272) q[3];
sx q[3];
rz(-0.47401014) q[3];
sx q[3];
rz(-1.43309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5804194) q[2];
sx q[2];
rz(-1.6757123) q[2];
sx q[2];
rz(-0.44211659) q[2];
rz(0.90028611) q[3];
sx q[3];
rz(-2.0787652) q[3];
sx q[3];
rz(3.1282848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1135547) q[0];
sx q[0];
rz(-0.91911626) q[0];
sx q[0];
rz(1.8369209) q[0];
rz(-1.4470944) q[1];
sx q[1];
rz(-1.9686952) q[1];
sx q[1];
rz(0.56217271) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2660414) q[0];
sx q[0];
rz(-2.3444359) q[0];
sx q[0];
rz(-2.6334727) q[0];
x q[1];
rz(-1.0400585) q[2];
sx q[2];
rz(-1.6518403) q[2];
sx q[2];
rz(-0.20866742) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0251384) q[1];
sx q[1];
rz(-2.0471731) q[1];
sx q[1];
rz(-0.78671793) q[1];
x q[2];
rz(1.9713054) q[3];
sx q[3];
rz(-1.9477961) q[3];
sx q[3];
rz(1.4124944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8582981) q[2];
sx q[2];
rz(-1.4026182) q[2];
sx q[2];
rz(2.1600058) q[2];
rz(-3.1372519) q[3];
sx q[3];
rz(-1.0467632) q[3];
sx q[3];
rz(-2.5481352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8269862) q[0];
sx q[0];
rz(-1.2207323) q[0];
sx q[0];
rz(3.0944968) q[0];
rz(-1.754952) q[1];
sx q[1];
rz(-1.1562693) q[1];
sx q[1];
rz(0.047853619) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3536516) q[0];
sx q[0];
rz(-1.5418058) q[0];
sx q[0];
rz(0.91584648) q[0];
rz(1.677703) q[2];
sx q[2];
rz(-2.6968743) q[2];
sx q[2];
rz(2.5267548) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.25803265) q[1];
sx q[1];
rz(-1.1415392) q[1];
sx q[1];
rz(2.2116823) q[1];
rz(-pi) q[2];
x q[2];
rz(2.738897) q[3];
sx q[3];
rz(-1.6214569) q[3];
sx q[3];
rz(0.5697054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4973267) q[2];
sx q[2];
rz(-1.1844119) q[2];
sx q[2];
rz(-0.93635526) q[2];
rz(2.3406773) q[3];
sx q[3];
rz(-0.51910669) q[3];
sx q[3];
rz(2.463533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8235648) q[0];
sx q[0];
rz(-0.93946539) q[0];
sx q[0];
rz(-2.9325624) q[0];
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
