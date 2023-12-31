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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3800669) q[0];
sx q[0];
rz(-2.4702284) q[0];
sx q[0];
rz(-1.6094366) q[0];
x q[1];
rz(0.36926271) q[2];
sx q[2];
rz(-2.5662176) q[2];
sx q[2];
rz(1.0550261) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1945038) q[1];
sx q[1];
rz(-2.6852486) q[1];
sx q[1];
rz(-1.4049005) q[1];
x q[2];
rz(-2.3000852) q[3];
sx q[3];
rz(-2.2268647) q[3];
sx q[3];
rz(-0.63667242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.67316002) q[2];
sx q[2];
rz(-2.1940553) q[2];
sx q[2];
rz(-0.49896487) q[2];
rz(0.4237825) q[3];
sx q[3];
rz(-1.9215923) q[3];
sx q[3];
rz(0.012714816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38133165) q[0];
sx q[0];
rz(-1.6853764) q[0];
sx q[0];
rz(2.1495842) q[0];
rz(-2.965062) q[1];
sx q[1];
rz(-2.9361528) q[1];
sx q[1];
rz(-0.36546779) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27643833) q[0];
sx q[0];
rz(-1.3197474) q[0];
sx q[0];
rz(2.6692078) q[0];
x q[1];
rz(0.89895504) q[2];
sx q[2];
rz(-2.0546753) q[2];
sx q[2];
rz(1.8091786) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9531627) q[1];
sx q[1];
rz(-0.97524446) q[1];
sx q[1];
rz(1.9838536) q[1];
x q[2];
rz(1.69181) q[3];
sx q[3];
rz(-1.3711208) q[3];
sx q[3];
rz(-1.6691085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3084597) q[2];
sx q[2];
rz(-0.59811991) q[2];
sx q[2];
rz(-2.0642521) q[2];
rz(-2.9789553) q[3];
sx q[3];
rz(-1.2626303) q[3];
sx q[3];
rz(-1.2619654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66115528) q[0];
sx q[0];
rz(-2.2876331) q[0];
sx q[0];
rz(0.61082947) q[0];
rz(-1.707533) q[1];
sx q[1];
rz(-1.3620946) q[1];
sx q[1];
rz(-2.6909713) q[1];
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
x q[1];
rz(-1.9446104) q[2];
sx q[2];
rz(-1.2842368) q[2];
sx q[2];
rz(-0.26619226) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.140481) q[1];
sx q[1];
rz(-1.2621659) q[1];
sx q[1];
rz(2.4806697) q[1];
rz(-0.30626014) q[3];
sx q[3];
rz(-0.89950409) q[3];
sx q[3];
rz(-1.1989532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.11015636) q[2];
sx q[2];
rz(-1.0855731) q[2];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0464756) q[0];
sx q[0];
rz(-1.7549544) q[0];
sx q[0];
rz(-1.1412096) q[0];
rz(-1.2966688) q[1];
sx q[1];
rz(-0.43162391) q[1];
sx q[1];
rz(-2.8362714) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0062549) q[0];
sx q[0];
rz(-2.4843785) q[0];
sx q[0];
rz(-2.6150136) q[0];
rz(2.6778053) q[2];
sx q[2];
rz(-1.1897435) q[2];
sx q[2];
rz(-0.87069893) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2480337) q[1];
sx q[1];
rz(-2.3602924) q[1];
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
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2083464) q[2];
sx q[2];
rz(-1.5657921) q[2];
sx q[2];
rz(0.34710458) q[2];
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
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40889302) q[0];
sx q[0];
rz(-2.2786031) q[0];
sx q[0];
rz(-0.53891671) q[0];
rz(-1.4096889) q[1];
sx q[1];
rz(-2.6795645) q[1];
sx q[1];
rz(2.7358823) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1203128) q[0];
sx q[0];
rz(-0.73182523) q[0];
sx q[0];
rz(2.7243607) q[0];
rz(-1.1666127) q[2];
sx q[2];
rz(-0.6037854) q[2];
sx q[2];
rz(-3.0886138) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.09771422) q[1];
sx q[1];
rz(-1.8227302) q[1];
sx q[1];
rz(1.9240727) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3375072) q[3];
sx q[3];
rz(-0.81987112) q[3];
sx q[3];
rz(1.5017205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6695909) q[2];
sx q[2];
rz(-0.80299032) q[2];
sx q[2];
rz(-1.4542788) q[2];
rz(2.3029095) q[3];
sx q[3];
rz(-1.2861774) q[3];
sx q[3];
rz(-1.3735166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8307761) q[0];
sx q[0];
rz(-0.61284471) q[0];
sx q[0];
rz(0.40711656) q[0];
rz(2.282417) q[1];
sx q[1];
rz(-1.5770864) q[1];
sx q[1];
rz(1.3173332) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9068127) q[0];
sx q[0];
rz(-0.54348511) q[0];
sx q[0];
rz(-1.6140922) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1520734) q[2];
sx q[2];
rz(-2.1534854) q[2];
sx q[2];
rz(-1.4942102) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0560318) q[1];
sx q[1];
rz(-0.83562682) q[1];
sx q[1];
rz(-0.0037770165) q[1];
rz(-0.96885724) q[3];
sx q[3];
rz(-2.4638306) q[3];
sx q[3];
rz(2.4399151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.09981) q[2];
sx q[2];
rz(-2.2618497) q[2];
sx q[2];
rz(2.2463634) q[2];
rz(-1.2081395) q[3];
sx q[3];
rz(-2.4619305) q[3];
sx q[3];
rz(-0.7081379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
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
rz(-1.3393341) q[1];
sx q[1];
rz(1.9099265) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8679778) q[0];
sx q[0];
rz(-0.6379188) q[0];
sx q[0];
rz(-2.5173553) q[0];
x q[1];
rz(-0.77881323) q[2];
sx q[2];
rz(-0.80873571) q[2];
sx q[2];
rz(-0.29701172) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.98037887) q[1];
sx q[1];
rz(-1.3413789) q[1];
sx q[1];
rz(1.4646261) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6024186) q[3];
sx q[3];
rz(-2.5368689) q[3];
sx q[3];
rz(-1.5550169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.82905519) q[2];
sx q[2];
rz(-2.2237491) q[2];
sx q[2];
rz(-2.5404239) q[2];
rz(-0.51856315) q[3];
sx q[3];
rz(-0.31373599) q[3];
sx q[3];
rz(-1.4275838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.9940051) q[0];
sx q[0];
rz(-3.1412791) q[0];
sx q[0];
rz(-3.0684379) q[0];
rz(-2.1144497) q[1];
sx q[1];
rz(-1.5348397) q[1];
sx q[1];
rz(0.59246078) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62646482) q[0];
sx q[0];
rz(-0.98083064) q[0];
sx q[0];
rz(1.626684) q[0];
rz(-pi) q[1];
rz(0.7166491) q[2];
sx q[2];
rz(-0.71084329) q[2];
sx q[2];
rz(-2.5568642) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5432424) q[1];
sx q[1];
rz(-0.91270869) q[1];
sx q[1];
rz(-0.42037085) q[1];
x q[2];
rz(-0.98832066) q[3];
sx q[3];
rz(-0.47401014) q[3];
sx q[3];
rz(1.7085027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.56117326) q[2];
sx q[2];
rz(-1.4658804) q[2];
sx q[2];
rz(2.6994761) q[2];
rz(-2.2413065) q[3];
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
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1135547) q[0];
sx q[0];
rz(-0.91911626) q[0];
sx q[0];
rz(-1.8369209) q[0];
rz(-1.6944983) q[1];
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
rz(2.2660414) q[0];
sx q[0];
rz(-2.3444359) q[0];
sx q[0];
rz(-0.50811998) q[0];
rz(0.093900605) q[2];
sx q[2];
rz(-2.0996089) q[2];
sx q[2];
rz(1.7319861) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0251384) q[1];
sx q[1];
rz(-2.0471731) q[1];
sx q[1];
rz(-0.78671793) q[1];
rz(0.4060679) q[3];
sx q[3];
rz(-1.9417524) q[3];
sx q[3];
rz(0.0036811034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8582981) q[2];
sx q[2];
rz(-1.7389745) q[2];
sx q[2];
rz(-0.98158681) q[2];
rz(3.1372519) q[3];
sx q[3];
rz(-2.0948295) q[3];
sx q[3];
rz(-2.5481352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8269862) q[0];
sx q[0];
rz(-1.2207323) q[0];
sx q[0];
rz(3.0944968) q[0];
rz(1.3866407) q[1];
sx q[1];
rz(-1.1562693) q[1];
sx q[1];
rz(0.047853619) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78794107) q[0];
sx q[0];
rz(-1.5997868) q[0];
sx q[0];
rz(2.2257462) q[0];
x q[1];
rz(0.050806413) q[2];
sx q[2];
rz(-2.0127957) q[2];
sx q[2];
rz(0.49651422) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8216985) q[1];
sx q[1];
rz(-0.7541637) q[1];
sx q[1];
rz(-0.9174222) q[1];
rz(1.6258532) q[3];
sx q[3];
rz(-1.1686472) q[3];
sx q[3];
rz(1.0226585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4973267) q[2];
sx q[2];
rz(-1.9571807) q[2];
sx q[2];
rz(-0.93635526) q[2];
rz(0.80091536) q[3];
sx q[3];
rz(-0.51910669) q[3];
sx q[3];
rz(0.67805964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3180278) q[0];
sx q[0];
rz(-2.2021273) q[0];
sx q[0];
rz(0.20903023) q[0];
rz(2.6387852) q[1];
sx q[1];
rz(-1.3192758) q[1];
sx q[1];
rz(-2.1154138) q[1];
rz(3.1059601) q[2];
sx q[2];
rz(-0.55152668) q[2];
sx q[2];
rz(1.3338911) q[2];
rz(-1.9746642) q[3];
sx q[3];
rz(-1.231791) q[3];
sx q[3];
rz(0.024699208) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
