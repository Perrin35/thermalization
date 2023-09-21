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
rz(-0.916565) q[1];
sx q[1];
rz(-0.7601127) q[1];
sx q[1];
rz(2.3488933) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9206031) q[0];
sx q[0];
rz(-1.5467637) q[0];
sx q[0];
rz(0.8997957) q[0];
rz(-0.36926271) q[2];
sx q[2];
rz(-2.5662176) q[2];
sx q[2];
rz(-1.0550261) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.9470889) q[1];
sx q[1];
rz(-2.6852486) q[1];
sx q[1];
rz(1.4049005) q[1];
rz(2.4281265) q[3];
sx q[3];
rz(-2.2029075) q[3];
sx q[3];
rz(1.6085094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.67316002) q[2];
sx q[2];
rz(-2.1940553) q[2];
sx q[2];
rz(0.49896487) q[2];
rz(2.7178102) q[3];
sx q[3];
rz(-1.2200004) q[3];
sx q[3];
rz(-3.1288778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.760261) q[0];
sx q[0];
rz(-1.6853764) q[0];
sx q[0];
rz(-0.99200845) q[0];
rz(0.17653067) q[1];
sx q[1];
rz(-2.9361528) q[1];
sx q[1];
rz(-0.36546779) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8651543) q[0];
sx q[0];
rz(-1.3197474) q[0];
sx q[0];
rz(2.6692078) q[0];
x q[1];
rz(-2.5502657) q[2];
sx q[2];
rz(-0.98726833) q[2];
sx q[2];
rz(0.11596767) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6668015) q[1];
sx q[1];
rz(-2.431369) q[1];
sx q[1];
rz(-0.53479654) q[1];
rz(-pi) q[2];
rz(-0.5378546) q[3];
sx q[3];
rz(-2.9085277) q[3];
sx q[3];
rz(-0.92249289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.83313292) q[2];
sx q[2];
rz(-0.59811991) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66115528) q[0];
sx q[0];
rz(-0.85395956) q[0];
sx q[0];
rz(-0.61082947) q[0];
rz(-1.4340596) q[1];
sx q[1];
rz(-1.3620946) q[1];
sx q[1];
rz(-0.45062137) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.492309) q[0];
sx q[0];
rz(-1.5376687) q[0];
sx q[0];
rz(1.6487728) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2497555) q[2];
sx q[2];
rz(-0.46687451) q[2];
sx q[2];
rz(2.4613949) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.9420942) q[1];
sx q[1];
rz(-0.71951413) q[1];
sx q[1];
rz(-0.47902963) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9335453) q[3];
sx q[3];
rz(-2.4137073) q[3];
sx q[3];
rz(2.4129652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0314363) q[2];
sx q[2];
rz(-1.0855731) q[2];
sx q[2];
rz(-1.191656) q[2];
rz(-0.84447652) q[3];
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
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0464756) q[0];
sx q[0];
rz(-1.3866383) q[0];
sx q[0];
rz(2.0003831) q[0];
rz(1.2966688) q[1];
sx q[1];
rz(-0.43162391) q[1];
sx q[1];
rz(2.8362714) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6396219) q[0];
sx q[0];
rz(-2.1272215) q[0];
sx q[0];
rz(-1.200838) q[0];
rz(-1.9919473) q[2];
sx q[2];
rz(-1.9990168) q[2];
sx q[2];
rz(0.8840094) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.893559) q[1];
sx q[1];
rz(-0.78130022) q[1];
sx q[1];
rz(-1.4206738) q[1];
x q[2];
rz(-0.51586113) q[3];
sx q[3];
rz(-2.0277884) q[3];
sx q[3];
rz(-2.7716694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2083464) q[2];
sx q[2];
rz(-1.5657921) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40889302) q[0];
sx q[0];
rz(-0.86298958) q[0];
sx q[0];
rz(-0.53891671) q[0];
rz(-1.7319038) q[1];
sx q[1];
rz(-0.46202818) q[1];
sx q[1];
rz(-0.40571037) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.273542) q[0];
sx q[0];
rz(-1.8450071) q[0];
sx q[0];
rz(-2.4540841) q[0];
x q[1];
rz(-1.1666127) q[2];
sx q[2];
rz(-0.6037854) q[2];
sx q[2];
rz(-3.0886138) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0438784) q[1];
sx q[1];
rz(-1.3188625) q[1];
sx q[1];
rz(-1.21752) q[1];
rz(-pi) q[2];
rz(1.8040854) q[3];
sx q[3];
rz(-0.81987112) q[3];
sx q[3];
rz(-1.5017205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6695909) q[2];
sx q[2];
rz(-2.3386023) q[2];
sx q[2];
rz(1.4542788) q[2];
rz(0.83868319) q[3];
sx q[3];
rz(-1.2861774) q[3];
sx q[3];
rz(1.3735166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31081653) q[0];
sx q[0];
rz(-2.5287479) q[0];
sx q[0];
rz(-0.40711656) q[0];
rz(0.85917568) q[1];
sx q[1];
rz(-1.5645063) q[1];
sx q[1];
rz(-1.8242594) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9068127) q[0];
sx q[0];
rz(-2.5981075) q[0];
sx q[0];
rz(-1.6140922) q[0];
rz(0.69465722) q[2];
sx q[2];
rz(-2.3431871) q[2];
sx q[2];
rz(0.77381221) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.091192186) q[1];
sx q[1];
rz(-0.7351774) q[1];
sx q[1];
rz(-1.5749732) q[1];
x q[2];
rz(0.98499961) q[3];
sx q[3];
rz(-1.2078152) q[3];
sx q[3];
rz(1.3604878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0417827) q[2];
sx q[2];
rz(-0.87974292) q[2];
sx q[2];
rz(-0.89522925) q[2];
rz(-1.2081395) q[3];
sx q[3];
rz(-0.67966214) q[3];
sx q[3];
rz(-2.4334548) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5773425) q[0];
sx q[0];
rz(-2.0944493) q[0];
sx q[0];
rz(1.1279001) q[0];
rz(-2.9290507) q[1];
sx q[1];
rz(-1.8022585) q[1];
sx q[1];
rz(1.2316661) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8679778) q[0];
sx q[0];
rz(-2.5036739) q[0];
sx q[0];
rz(0.62423737) q[0];
rz(-pi) q[1];
rz(0.93630837) q[2];
sx q[2];
rz(-2.1116743) q[2];
sx q[2];
rz(-2.4782431) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.98037887) q[1];
sx q[1];
rz(-1.3413789) q[1];
sx q[1];
rz(1.4646261) q[1];
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
rz(pi/2) q[1];
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
rz(-1.7140088) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1475875) q[0];
sx q[0];
rz(-0.00031358263) q[0];
sx q[0];
rz(3.0684379) q[0];
rz(-1.0271429) q[1];
sx q[1];
rz(-1.5348397) q[1];
sx q[1];
rz(2.5491319) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5151278) q[0];
sx q[0];
rz(-2.160762) q[0];
sx q[0];
rz(1.5149087) q[0];
rz(-0.7166491) q[2];
sx q[2];
rz(-2.4307494) q[2];
sx q[2];
rz(0.58472842) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.91208273) q[1];
sx q[1];
rz(-2.3778937) q[1];
sx q[1];
rz(-1.0850701) q[1];
rz(1.1660277) q[3];
sx q[3];
rz(-1.824607) q[3];
sx q[3];
rz(2.7491731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.56117326) q[2];
sx q[2];
rz(-1.4658804) q[2];
sx q[2];
rz(0.44211659) q[2];
rz(0.90028611) q[3];
sx q[3];
rz(-1.0628275) q[3];
sx q[3];
rz(0.013307868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(3.1135547) q[0];
sx q[0];
rz(-2.2224764) q[0];
sx q[0];
rz(-1.3046718) q[0];
rz(1.6944983) q[1];
sx q[1];
rz(-1.9686952) q[1];
sx q[1];
rz(-2.5794199) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5931372) q[0];
sx q[0];
rz(-0.89567417) q[0];
sx q[0];
rz(-1.1086585) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7299037) q[2];
sx q[2];
rz(-0.53630398) q[2];
sx q[2];
rz(-1.916534) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1162029) q[1];
sx q[1];
rz(-0.89239489) q[1];
sx q[1];
rz(0.62979001) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9713054) q[3];
sx q[3];
rz(-1.1937965) q[3];
sx q[3];
rz(-1.7290982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.28329453) q[2];
sx q[2];
rz(-1.7389745) q[2];
sx q[2];
rz(0.98158681) q[2];
rz(-3.1372519) q[3];
sx q[3];
rz(-1.0467632) q[3];
sx q[3];
rz(0.59345746) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3146064) q[0];
sx q[0];
rz(-1.9208603) q[0];
sx q[0];
rz(-3.0944968) q[0];
rz(1.3866407) q[1];
sx q[1];
rz(-1.1562693) q[1];
sx q[1];
rz(0.047853619) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-1.4638897) q[2];
sx q[2];
rz(-2.6968743) q[2];
sx q[2];
rz(2.5267548) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8216985) q[1];
sx q[1];
rz(-0.7541637) q[1];
sx q[1];
rz(-2.2241705) q[1];
rz(-pi) q[2];
x q[2];
rz(0.1286653) q[3];
sx q[3];
rz(-2.7358958) q[3];
sx q[3];
rz(0.88276827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4973267) q[2];
sx q[2];
rz(-1.1844119) q[2];
sx q[2];
rz(0.93635526) q[2];
rz(-2.3406773) q[3];
sx q[3];
rz(-0.51910669) q[3];
sx q[3];
rz(-2.463533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.8235648) q[0];
sx q[0];
rz(-2.2021273) q[0];
sx q[0];
rz(0.20903023) q[0];
rz(2.6387852) q[1];
sx q[1];
rz(-1.3192758) q[1];
sx q[1];
rz(-2.1154138) q[1];
rz(-1.5488831) q[2];
sx q[2];
rz(-1.0196601) q[2];
sx q[2];
rz(-1.7658726) q[2];
rz(2.7754178) q[3];
sx q[3];
rz(-1.9504642) q[3];
sx q[3];
rz(-1.4049243) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
