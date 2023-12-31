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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3307318) q[0];
sx q[0];
rz(-2.2415677) q[0];
sx q[0];
rz(0.030681507) q[0];
rz(-pi) q[1];
x q[1];
rz(0.36926271) q[2];
sx q[2];
rz(-2.5662176) q[2];
sx q[2];
rz(-2.0865666) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.77289171) q[1];
sx q[1];
rz(-1.643631) q[1];
sx q[1];
rz(1.1198977) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.84150746) q[3];
sx q[3];
rz(-2.2268647) q[3];
sx q[3];
rz(0.63667242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4684326) q[2];
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
rz(-pi/2) q[1];
x q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38133165) q[0];
sx q[0];
rz(-1.4562162) q[0];
sx q[0];
rz(0.99200845) q[0];
rz(-0.17653067) q[1];
sx q[1];
rz(-2.9361528) q[1];
sx q[1];
rz(0.36546779) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1680981) q[0];
sx q[0];
rz(-1.1143648) q[0];
sx q[0];
rz(-1.2903851) q[0];
x q[1];
rz(2.5502657) q[2];
sx q[2];
rz(-2.1543243) q[2];
sx q[2];
rz(-3.025625) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6234399) q[1];
sx q[1];
rz(-1.9095416) q[1];
sx q[1];
rz(0.63697302) q[1];
rz(-pi) q[2];
rz(1.4497827) q[3];
sx q[3];
rz(-1.7704718) q[3];
sx q[3];
rz(1.4724842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.83313292) q[2];
sx q[2];
rz(-0.59811991) q[2];
sx q[2];
rz(1.0773405) q[2];
rz(-2.9789553) q[3];
sx q[3];
rz(-1.8789623) q[3];
sx q[3];
rz(1.2619654) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66115528) q[0];
sx q[0];
rz(-2.2876331) q[0];
sx q[0];
rz(0.61082947) q[0];
rz(1.707533) q[1];
sx q[1];
rz(-1.779498) q[1];
sx q[1];
rz(0.45062137) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.492309) q[0];
sx q[0];
rz(-1.5376687) q[0];
sx q[0];
rz(1.6487728) q[0];
x q[1];
rz(2.2497555) q[2];
sx q[2];
rz(-2.6747181) q[2];
sx q[2];
rz(0.6801978) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.140481) q[1];
sx q[1];
rz(-1.2621659) q[1];
sx q[1];
rz(2.4806697) q[1];
rz(-1.2080473) q[3];
sx q[3];
rz(-2.4137073) q[3];
sx q[3];
rz(0.7286275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0314363) q[2];
sx q[2];
rz(-1.0855731) q[2];
sx q[2];
rz(-1.9499367) q[2];
rz(-2.2971161) q[3];
sx q[3];
rz(-2.8561487) q[3];
sx q[3];
rz(0.58117956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.095117005) q[0];
sx q[0];
rz(-1.3866383) q[0];
sx q[0];
rz(-1.1412096) q[0];
rz(-1.8449239) q[1];
sx q[1];
rz(-2.7099687) q[1];
sx q[1];
rz(0.30532125) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0083864) q[0];
sx q[0];
rz(-1.2587233) q[0];
sx q[0];
rz(-0.58831373) q[0];
rz(-pi) q[1];
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
sx q[0];
rz(pi/2) q[0];
rz(-1.2480337) q[1];
sx q[1];
rz(-2.3602924) q[1];
sx q[1];
rz(1.7209189) q[1];
rz(-0.51586113) q[3];
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
rz(-1.5758005) q[2];
sx q[2];
rz(-0.34710458) q[2];
rz(0.38337213) q[3];
sx q[3];
rz(-2.255286) q[3];
sx q[3];
rz(-2.0736407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7326996) q[0];
sx q[0];
rz(-2.2786031) q[0];
sx q[0];
rz(0.53891671) q[0];
rz(1.4096889) q[1];
sx q[1];
rz(-2.6795645) q[1];
sx q[1];
rz(0.40571037) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51605326) q[0];
sx q[0];
rz(-2.2279983) q[0];
sx q[0];
rz(1.2217193) q[0];
x q[1];
rz(-2.1359372) q[2];
sx q[2];
rz(-1.3456151) q[2];
sx q[2];
rz(-1.2852247) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.87863641) q[1];
sx q[1];
rz(-0.43081455) q[1];
sx q[1];
rz(2.2104435) q[1];
rz(-pi) q[2];
x q[2];
rz(0.24280663) q[3];
sx q[3];
rz(-2.3620776) q[3];
sx q[3];
rz(-1.8368343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.47200176) q[2];
sx q[2];
rz(-2.3386023) q[2];
sx q[2];
rz(1.4542788) q[2];
rz(-2.3029095) q[3];
sx q[3];
rz(-1.2861774) q[3];
sx q[3];
rz(-1.7680761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8307761) q[0];
sx q[0];
rz(-2.5287479) q[0];
sx q[0];
rz(0.40711656) q[0];
rz(-2.282417) q[1];
sx q[1];
rz(-1.5770864) q[1];
sx q[1];
rz(1.8242594) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2989527) q[0];
sx q[0];
rz(-1.5484122) q[0];
sx q[0];
rz(2.1138666) q[0];
rz(-2.1520734) q[2];
sx q[2];
rz(-0.98810722) q[2];
sx q[2];
rz(1.4942102) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.091192186) q[1];
sx q[1];
rz(-2.4064153) q[1];
sx q[1];
rz(-1.5666195) q[1];
rz(-2.7139211) q[3];
sx q[3];
rz(-2.1139522) q[3];
sx q[3];
rz(3.1205408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0417827) q[2];
sx q[2];
rz(-0.87974292) q[2];
sx q[2];
rz(0.89522925) q[2];
rz(1.2081395) q[3];
sx q[3];
rz(-0.67966214) q[3];
sx q[3];
rz(2.4334548) q[3];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5642501) q[0];
sx q[0];
rz(-1.0471434) q[0];
sx q[0];
rz(1.1279001) q[0];
rz(-2.9290507) q[1];
sx q[1];
rz(-1.3393341) q[1];
sx q[1];
rz(1.9099265) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5989482) q[0];
sx q[0];
rz(-1.066474) q[0];
sx q[0];
rz(1.1619316) q[0];
rz(-pi) q[1];
rz(2.3627794) q[2];
sx q[2];
rz(-2.3328569) q[2];
sx q[2];
rz(-2.8445809) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.54214189) q[1];
sx q[1];
rz(-2.8891924) q[1];
sx q[1];
rz(-0.42599328) q[1];
rz(1.6024186) q[3];
sx q[3];
rz(-2.5368689) q[3];
sx q[3];
rz(1.5865758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3125375) q[2];
sx q[2];
rz(-0.91784358) q[2];
sx q[2];
rz(-0.60116872) q[2];
rz(2.6230295) q[3];
sx q[3];
rz(-0.31373599) q[3];
sx q[3];
rz(-1.4275838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9940051) q[0];
sx q[0];
rz(-0.00031358263) q[0];
sx q[0];
rz(0.073154733) q[0];
rz(2.1144497) q[1];
sx q[1];
rz(-1.5348397) q[1];
sx q[1];
rz(2.5491319) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97544599) q[0];
sx q[0];
rz(-1.6172292) q[0];
sx q[0];
rz(-0.59068824) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0855037) q[2];
sx q[2];
rz(-1.0564431) q[2];
sx q[2];
rz(2.8714542) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5432424) q[1];
sx q[1];
rz(-0.91270869) q[1];
sx q[1];
rz(2.7212218) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9755649) q[3];
sx q[3];
rz(-1.824607) q[3];
sx q[3];
rz(-0.39241957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5804194) q[2];
sx q[2];
rz(-1.6757123) q[2];
sx q[2];
rz(0.44211659) q[2];
rz(-2.2413065) q[3];
sx q[3];
rz(-2.0787652) q[3];
sx q[3];
rz(3.1282848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028037926) q[0];
sx q[0];
rz(-0.91911626) q[0];
sx q[0];
rz(-1.8369209) q[0];
rz(-1.4470944) q[1];
sx q[1];
rz(-1.9686952) q[1];
sx q[1];
rz(0.56217271) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2660414) q[0];
sx q[0];
rz(-2.3444359) q[0];
sx q[0];
rz(-0.50811998) q[0];
rz(0.093900605) q[2];
sx q[2];
rz(-2.0996089) q[2];
sx q[2];
rz(-1.4096066) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2562981) q[1];
sx q[1];
rz(-0.89026272) q[1];
sx q[1];
rz(0.93974944) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1702873) q[3];
sx q[3];
rz(-1.1937965) q[3];
sx q[3];
rz(1.4124944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.28329453) q[2];
sx q[2];
rz(-1.7389745) q[2];
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
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3146064) q[0];
sx q[0];
rz(-1.9208603) q[0];
sx q[0];
rz(3.0944968) q[0];
rz(-1.754952) q[1];
sx q[1];
rz(-1.1562693) q[1];
sx q[1];
rz(0.047853619) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3536516) q[0];
sx q[0];
rz(-1.5997868) q[0];
sx q[0];
rz(2.2257462) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0132952) q[2];
sx q[2];
rz(-1.5248761) q[2];
sx q[2];
rz(1.0525345) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1298444) q[1];
sx q[1];
rz(-2.1455892) q[1];
sx q[1];
rz(2.6227436) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5157394) q[3];
sx q[3];
rz(-1.1686472) q[3];
sx q[3];
rz(-2.1189342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.64426595) q[2];
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
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
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
rz(0.50280747) q[1];
sx q[1];
rz(-1.8223169) q[1];
sx q[1];
rz(1.0261789) q[1];
rz(-1.5927096) q[2];
sx q[2];
rz(-2.1219325) q[2];
sx q[2];
rz(1.3757201) q[2];
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
