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
rz(-0.7926994) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3307318) q[0];
sx q[0];
rz(-2.2415677) q[0];
sx q[0];
rz(3.1109111) q[0];
rz(2.5976074) q[2];
sx q[2];
rz(-1.3731125) q[2];
sx q[2];
rz(-0.20180861) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3789062) q[1];
sx q[1];
rz(-1.1211809) q[1];
sx q[1];
rz(-3.0607037) q[1];
x q[2];
rz(0.84150746) q[3];
sx q[3];
rz(-0.91472799) q[3];
sx q[3];
rz(0.63667242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4684326) q[2];
sx q[2];
rz(-2.1940553) q[2];
sx q[2];
rz(0.49896487) q[2];
rz(2.7178102) q[3];
sx q[3];
rz(-1.9215923) q[3];
sx q[3];
rz(3.1288778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-0.17653067) q[1];
sx q[1];
rz(-2.9361528) q[1];
sx q[1];
rz(0.36546779) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1680981) q[0];
sx q[0];
rz(-1.1143648) q[0];
sx q[0];
rz(-1.8512076) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.86958779) q[2];
sx q[2];
rz(-2.3361932) q[2];
sx q[2];
rz(-0.7676917) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5181527) q[1];
sx q[1];
rz(-1.232051) q[1];
sx q[1];
rz(2.5046196) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.69181) q[3];
sx q[3];
rz(-1.7704718) q[3];
sx q[3];
rz(-1.6691085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3084597) q[2];
sx q[2];
rz(-0.59811991) q[2];
sx q[2];
rz(2.0642521) q[2];
rz(0.16263738) q[3];
sx q[3];
rz(-1.8789623) q[3];
sx q[3];
rz(-1.8796273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66115528) q[0];
sx q[0];
rz(-2.2876331) q[0];
sx q[0];
rz(2.5307632) q[0];
rz(1.4340596) q[1];
sx q[1];
rz(-1.779498) q[1];
sx q[1];
rz(-0.45062137) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47942802) q[0];
sx q[0];
rz(-0.084708609) q[0];
sx q[0];
rz(1.9730294) q[0];
rz(-pi) q[1];
rz(1.9446104) q[2];
sx q[2];
rz(-1.2842368) q[2];
sx q[2];
rz(-2.8754004) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.33774099) q[1];
sx q[1];
rz(-2.1954814) q[1];
sx q[1];
rz(-1.1869663) q[1];
x q[2];
rz(2.8353325) q[3];
sx q[3];
rz(-2.2420886) q[3];
sx q[3];
rz(-1.9426395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0314363) q[2];
sx q[2];
rz(-1.0855731) q[2];
sx q[2];
rz(1.191656) q[2];
rz(-2.2971161) q[3];
sx q[3];
rz(-2.8561487) q[3];
sx q[3];
rz(0.58117956) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.095117005) q[0];
sx q[0];
rz(-1.7549544) q[0];
sx q[0];
rz(2.0003831) q[0];
rz(1.8449239) q[1];
sx q[1];
rz(-0.43162391) q[1];
sx q[1];
rz(-2.8362714) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1353377) q[0];
sx q[0];
rz(-2.4843785) q[0];
sx q[0];
rz(-2.6150136) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9919473) q[2];
sx q[2];
rz(-1.9990168) q[2];
sx q[2];
rz(-0.8840094) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(pi/2) q[1];
rz(-1.9332463) q[2];
sx q[2];
rz(-1.5657921) q[2];
sx q[2];
rz(-2.7944881) q[2];
rz(-2.7582205) q[3];
sx q[3];
rz(-0.88630668) q[3];
sx q[3];
rz(-1.0679519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7326996) q[0];
sx q[0];
rz(-2.2786031) q[0];
sx q[0];
rz(-0.53891671) q[0];
rz(-1.7319038) q[1];
sx q[1];
rz(-0.46202818) q[1];
sx q[1];
rz(2.7358823) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8680506) q[0];
sx q[0];
rz(-1.8450071) q[0];
sx q[0];
rz(0.6875086) q[0];
rz(-pi) q[1];
rz(1.0056555) q[2];
sx q[2];
rz(-1.3456151) q[2];
sx q[2];
rz(-1.2852247) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2629562) q[1];
sx q[1];
rz(-0.43081455) q[1];
sx q[1];
rz(2.2104435) q[1];
rz(-0.76463502) q[3];
sx q[3];
rz(-1.4009762) q[3];
sx q[3];
rz(3.0499383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6695909) q[2];
sx q[2];
rz(-2.3386023) q[2];
sx q[2];
rz(1.6873138) q[2];
rz(-2.3029095) q[3];
sx q[3];
rz(-1.2861774) q[3];
sx q[3];
rz(1.3735166) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31081653) q[0];
sx q[0];
rz(-0.61284471) q[0];
sx q[0];
rz(2.7344761) q[0];
rz(-2.282417) q[1];
sx q[1];
rz(-1.5770864) q[1];
sx q[1];
rz(-1.3173332) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9068127) q[0];
sx q[0];
rz(-0.54348511) q[0];
sx q[0];
rz(-1.6140922) q[0];
rz(-pi) q[1];
rz(2.1520734) q[2];
sx q[2];
rz(-2.1534854) q[2];
sx q[2];
rz(-1.6473824) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4827022) q[1];
sx q[1];
rz(-1.5679949) q[1];
sx q[1];
rz(0.83562327) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.42767151) q[3];
sx q[3];
rz(-1.0276405) q[3];
sx q[3];
rz(-0.021051858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.09981) q[2];
sx q[2];
rz(-0.87974292) q[2];
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
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3198277) q[0];
sx q[0];
rz(-1.2152831) q[0];
sx q[0];
rz(0.54152821) q[0];
x q[1];
rz(2.5008051) q[2];
sx q[2];
rz(-2.1038342) q[2];
sx q[2];
rz(1.2696881) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.98037887) q[1];
sx q[1];
rz(-1.8002138) q[1];
sx q[1];
rz(-1.6769665) q[1];
rz(-3.1197458) q[3];
sx q[3];
rz(-0.96641809) q[3];
sx q[3];
rz(1.6250087) q[3];
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
rz(2.6230295) q[3];
sx q[3];
rz(-0.31373599) q[3];
sx q[3];
rz(1.7140088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1475875) q[0];
sx q[0];
rz(-3.1412791) q[0];
sx q[0];
rz(-3.0684379) q[0];
rz(-1.0271429) q[1];
sx q[1];
rz(-1.5348397) q[1];
sx q[1];
rz(2.5491319) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97544599) q[0];
sx q[0];
rz(-1.6172292) q[0];
sx q[0];
rz(2.5509044) q[0];
rz(0.7166491) q[2];
sx q[2];
rz(-2.4307494) q[2];
sx q[2];
rz(2.5568642) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5432424) q[1];
sx q[1];
rz(-2.228884) q[1];
sx q[1];
rz(-2.7212218) q[1];
x q[2];
rz(1.9755649) q[3];
sx q[3];
rz(-1.824607) q[3];
sx q[3];
rz(0.39241957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5804194) q[2];
sx q[2];
rz(-1.4658804) q[2];
sx q[2];
rz(2.6994761) q[2];
rz(2.2413065) q[3];
sx q[3];
rz(-1.0628275) q[3];
sx q[3];
rz(3.1282848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5931372) q[0];
sx q[0];
rz(-0.89567417) q[0];
sx q[0];
rz(2.0329342) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7299037) q[2];
sx q[2];
rz(-0.53630398) q[2];
sx q[2];
rz(-1.916534) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2562981) q[1];
sx q[1];
rz(-0.89026272) q[1];
sx q[1];
rz(-2.2018432) q[1];
x q[2];
rz(-1.1702873) q[3];
sx q[3];
rz(-1.1937965) q[3];
sx q[3];
rz(1.7290982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.28329453) q[2];
sx q[2];
rz(-1.4026182) q[2];
sx q[2];
rz(-2.1600058) q[2];
rz(-0.0043407241) q[3];
sx q[3];
rz(-2.0948295) q[3];
sx q[3];
rz(0.59345746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8269862) q[0];
sx q[0];
rz(-1.2207323) q[0];
sx q[0];
rz(-0.04709588) q[0];
rz(1.754952) q[1];
sx q[1];
rz(-1.1562693) q[1];
sx q[1];
rz(-0.047853619) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3210147) q[0];
sx q[0];
rz(-0.65549675) q[0];
sx q[0];
rz(1.6183678) q[0];
x q[1];
rz(1.677703) q[2];
sx q[2];
rz(-0.44471834) q[2];
sx q[2];
rz(0.61483785) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0117482) q[1];
sx q[1];
rz(-0.99600345) q[1];
sx q[1];
rz(0.5188491) q[1];
rz(-0.1286653) q[3];
sx q[3];
rz(-2.7358958) q[3];
sx q[3];
rz(-0.88276827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4973267) q[2];
sx q[2];
rz(-1.1844119) q[2];
sx q[2];
rz(-0.93635526) q[2];
rz(-0.80091536) q[3];
sx q[3];
rz(-0.51910669) q[3];
sx q[3];
rz(-0.67805964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
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
rz(0.035632523) q[2];
sx q[2];
rz(-2.590066) q[2];
sx q[2];
rz(-1.8077015) q[2];
rz(-0.36617483) q[3];
sx q[3];
rz(-1.9504642) q[3];
sx q[3];
rz(-1.4049243) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
