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
rz(-0.65519873) q[0];
sx q[0];
rz(0.97638786) q[0];
rz(2.2250277) q[1];
sx q[1];
rz(-2.38148) q[1];
sx q[1];
rz(-2.3488933) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22098955) q[0];
sx q[0];
rz(-1.5467637) q[0];
sx q[0];
rz(0.8997957) q[0];
rz(-pi) q[1];
rz(1.8007457) q[2];
sx q[2];
rz(-2.1030428) q[2];
sx q[2];
rz(-1.4872273) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3789062) q[1];
sx q[1];
rz(-1.1211809) q[1];
sx q[1];
rz(0.080888918) q[1];
rz(-pi) q[2];
rz(-2.4281265) q[3];
sx q[3];
rz(-0.93868512) q[3];
sx q[3];
rz(1.6085094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.67316002) q[2];
sx q[2];
rz(-2.1940553) q[2];
sx q[2];
rz(0.49896487) q[2];
rz(2.7178102) q[3];
sx q[3];
rz(-1.9215923) q[3];
sx q[3];
rz(-0.012714816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38133165) q[0];
sx q[0];
rz(-1.4562162) q[0];
sx q[0];
rz(0.99200845) q[0];
rz(2.965062) q[1];
sx q[1];
rz(-2.9361528) q[1];
sx q[1];
rz(-2.7761249) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8651543) q[0];
sx q[0];
rz(-1.3197474) q[0];
sx q[0];
rz(-0.47238484) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2426376) q[2];
sx q[2];
rz(-1.0869173) q[2];
sx q[2];
rz(-1.3324141) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9531627) q[1];
sx q[1];
rz(-2.1663482) q[1];
sx q[1];
rz(-1.1577391) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4497827) q[3];
sx q[3];
rz(-1.3711208) q[3];
sx q[3];
rz(1.6691085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.83313292) q[2];
sx q[2];
rz(-2.5434727) q[2];
sx q[2];
rz(1.0773405) q[2];
rz(-0.16263738) q[3];
sx q[3];
rz(-1.2626303) q[3];
sx q[3];
rz(1.2619654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4804374) q[0];
sx q[0];
rz(-2.2876331) q[0];
sx q[0];
rz(-0.61082947) q[0];
rz(1.707533) q[1];
sx q[1];
rz(-1.779498) q[1];
sx q[1];
rz(-2.6909713) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47942802) q[0];
sx q[0];
rz(-0.084708609) q[0];
sx q[0];
rz(1.1685632) q[0];
rz(0.89183715) q[2];
sx q[2];
rz(-2.6747181) q[2];
sx q[2];
rz(-0.6801978) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8038517) q[1];
sx q[1];
rz(-2.1954814) q[1];
sx q[1];
rz(1.1869663) q[1];
x q[2];
rz(0.87617989) q[3];
sx q[3];
rz(-1.8091222) q[3];
sx q[3];
rz(-2.5755469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.11015636) q[2];
sx q[2];
rz(-2.0560196) q[2];
sx q[2];
rz(1.191656) q[2];
rz(2.2971161) q[3];
sx q[3];
rz(-0.28544393) q[3];
sx q[3];
rz(-2.5604131) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-2.7099687) q[1];
sx q[1];
rz(0.30532125) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5019708) q[0];
sx q[0];
rz(-2.1272215) q[0];
sx q[0];
rz(-1.9407546) q[0];
rz(-pi) q[1];
rz(-1.1496454) q[2];
sx q[2];
rz(-1.9990168) q[2];
sx q[2];
rz(2.2575833) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2480337) q[1];
sx q[1];
rz(-2.3602924) q[1];
sx q[1];
rz(-1.7209189) q[1];
x q[2];
rz(2.0852854) q[3];
sx q[3];
rz(-2.0293651) q[3];
sx q[3];
rz(-1.4460627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9332463) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7326996) q[0];
sx q[0];
rz(-0.86298958) q[0];
sx q[0];
rz(2.6026759) q[0];
rz(1.7319038) q[1];
sx q[1];
rz(-2.6795645) q[1];
sx q[1];
rz(-0.40571037) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6255394) q[0];
sx q[0];
rz(-2.2279983) q[0];
sx q[0];
rz(1.2217193) q[0];
rz(1.97498) q[2];
sx q[2];
rz(-2.5378072) q[2];
sx q[2];
rz(-0.052978901) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0438784) q[1];
sx q[1];
rz(-1.3188625) q[1];
sx q[1];
rz(1.9240727) q[1];
x q[2];
rz(-2.3769576) q[3];
sx q[3];
rz(-1.7406165) q[3];
sx q[3];
rz(-0.09165435) q[3];
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
rz(-2.3029095) q[3];
sx q[3];
rz(-1.8554153) q[3];
sx q[3];
rz(1.7680761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-0.61284471) q[0];
sx q[0];
rz(-0.40711656) q[0];
rz(0.85917568) q[1];
sx q[1];
rz(-1.5645063) q[1];
sx q[1];
rz(1.3173332) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9068127) q[0];
sx q[0];
rz(-0.54348511) q[0];
sx q[0];
rz(-1.5275005) q[0];
rz(-pi) q[1];
rz(0.98951927) q[2];
sx q[2];
rz(-2.1534854) q[2];
sx q[2];
rz(-1.4942102) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.091192186) q[1];
sx q[1];
rz(-2.4064153) q[1];
sx q[1];
rz(1.5749732) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.98499961) q[3];
sx q[3];
rz(-1.2078152) q[3];
sx q[3];
rz(-1.3604878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0417827) q[2];
sx q[2];
rz(-2.2618497) q[2];
sx q[2];
rz(-0.89522925) q[2];
rz(-1.9334531) q[3];
sx q[3];
rz(-2.4619305) q[3];
sx q[3];
rz(-2.4334548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-1.5773425) q[0];
sx q[0];
rz(-1.0471434) q[0];
sx q[0];
rz(-2.0136925) q[0];
rz(-2.9290507) q[1];
sx q[1];
rz(-1.3393341) q[1];
sx q[1];
rz(1.9099265) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3198277) q[0];
sx q[0];
rz(-1.2152831) q[0];
sx q[0];
rz(-2.6000644) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2052843) q[2];
sx q[2];
rz(-1.0299183) q[2];
sx q[2];
rz(-2.4782431) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5269446) q[1];
sx q[1];
rz(-1.467418) q[1];
sx q[1];
rz(-2.9109216) q[1];
rz(-pi) q[2];
x q[2];
rz(0.96630649) q[3];
sx q[3];
rz(-1.5887727) q[3];
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
rz(-2.2237491) q[2];
sx q[2];
rz(-0.60116872) q[2];
rz(-2.6230295) q[3];
sx q[3];
rz(-0.31373599) q[3];
sx q[3];
rz(1.4275838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1475875) q[0];
sx q[0];
rz(-0.00031358263) q[0];
sx q[0];
rz(0.073154733) q[0];
rz(-1.0271429) q[1];
sx q[1];
rz(-1.606753) q[1];
sx q[1];
rz(0.59246078) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5151278) q[0];
sx q[0];
rz(-2.160762) q[0];
sx q[0];
rz(1.5149087) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0560889) q[2];
sx q[2];
rz(-1.0564431) q[2];
sx q[2];
rz(-0.27013847) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5983502) q[1];
sx q[1];
rz(-2.228884) q[1];
sx q[1];
rz(2.7212218) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.153272) q[3];
sx q[3];
rz(-0.47401014) q[3];
sx q[3];
rz(1.43309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5804194) q[2];
sx q[2];
rz(-1.4658804) q[2];
sx q[2];
rz(0.44211659) q[2];
rz(2.2413065) q[3];
sx q[3];
rz(-2.0787652) q[3];
sx q[3];
rz(-3.1282848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.028037926) q[0];
sx q[0];
rz(-0.91911626) q[0];
sx q[0];
rz(1.8369209) q[0];
rz(1.4470944) q[1];
sx q[1];
rz(-1.1728975) q[1];
sx q[1];
rz(0.56217271) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87555128) q[0];
sx q[0];
rz(-0.79715675) q[0];
sx q[0];
rz(-0.50811998) q[0];
rz(-pi) q[1];
rz(-1.7299037) q[2];
sx q[2];
rz(-2.6052887) q[2];
sx q[2];
rz(-1.2250587) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.93974944) q[1];
rz(-2.3638944) q[3];
sx q[3];
rz(-0.54293984) q[3];
sx q[3];
rz(-0.87399769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8582981) q[2];
sx q[2];
rz(-1.4026182) q[2];
sx q[2];
rz(-0.98158681) q[2];
rz(0.0043407241) q[3];
sx q[3];
rz(-1.0467632) q[3];
sx q[3];
rz(-2.5481352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3146064) q[0];
sx q[0];
rz(-1.9208603) q[0];
sx q[0];
rz(3.0944968) q[0];
rz(1.754952) q[1];
sx q[1];
rz(-1.9853233) q[1];
sx q[1];
rz(-3.093739) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82057792) q[0];
sx q[0];
rz(-0.65549675) q[0];
sx q[0];
rz(1.5232248) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4638897) q[2];
sx q[2];
rz(-2.6968743) q[2];
sx q[2];
rz(-2.5267548) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.88356) q[1];
sx q[1];
rz(-2.0000534) q[1];
sx q[1];
rz(0.92991035) q[1];
rz(-pi) q[2];
rz(-0.40269561) q[3];
sx q[3];
rz(-1.6214569) q[3];
sx q[3];
rz(0.5697054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4973267) q[2];
sx q[2];
rz(-1.9571807) q[2];
sx q[2];
rz(0.93635526) q[2];
rz(-0.80091536) q[3];
sx q[3];
rz(-0.51910669) q[3];
sx q[3];
rz(-0.67805964) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8235648) q[0];
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
rz(0.83947635) q[3];
sx q[3];
rz(-0.52121938) q[3];
sx q[3];
rz(-2.2073707) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];