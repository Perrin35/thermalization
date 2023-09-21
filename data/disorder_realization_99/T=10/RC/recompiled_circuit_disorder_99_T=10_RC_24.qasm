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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9206031) q[0];
sx q[0];
rz(-1.594829) q[0];
sx q[0];
rz(0.8997957) q[0];
rz(-pi) q[1];
rz(2.7723299) q[2];
sx q[2];
rz(-0.57537503) q[2];
sx q[2];
rz(1.0550261) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3789062) q[1];
sx q[1];
rz(-1.1211809) q[1];
sx q[1];
rz(3.0607037) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3402432) q[3];
sx q[3];
rz(-1.0145463) q[3];
sx q[3];
rz(2.706561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4684326) q[2];
sx q[2];
rz(-0.94753733) q[2];
sx q[2];
rz(2.6426278) q[2];
rz(-0.4237825) q[3];
sx q[3];
rz(-1.9215923) q[3];
sx q[3];
rz(-0.012714816) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38133165) q[0];
sx q[0];
rz(-1.6853764) q[0];
sx q[0];
rz(-0.99200845) q[0];
rz(-0.17653067) q[1];
sx q[1];
rz(-2.9361528) q[1];
sx q[1];
rz(-2.7761249) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.394709) q[0];
sx q[0];
rz(-2.611126) q[0];
sx q[0];
rz(-0.51325004) q[0];
rz(-pi) q[1];
rz(-0.59132691) q[2];
sx q[2];
rz(-2.1543243) q[2];
sx q[2];
rz(0.11596767) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6234399) q[1];
sx q[1];
rz(-1.232051) q[1];
sx q[1];
rz(-2.5046196) q[1];
x q[2];
rz(2.9404853) q[3];
sx q[3];
rz(-1.6893941) q[3];
sx q[3];
rz(0.12242854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3084597) q[2];
sx q[2];
rz(-2.5434727) q[2];
sx q[2];
rz(-1.0773405) q[2];
rz(2.9789553) q[3];
sx q[3];
rz(-1.8789623) q[3];
sx q[3];
rz(-1.2619654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.4804374) q[0];
sx q[0];
rz(-0.85395956) q[0];
sx q[0];
rz(-2.5307632) q[0];
rz(1.4340596) q[1];
sx q[1];
rz(-1.3620946) q[1];
sx q[1];
rz(0.45062137) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47942802) q[0];
sx q[0];
rz(-0.084708609) q[0];
sx q[0];
rz(-1.1685632) q[0];
x q[1];
rz(0.89183715) q[2];
sx q[2];
rz(-0.46687451) q[2];
sx q[2];
rz(-2.4613949) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1994985) q[1];
sx q[1];
rz(-2.4220785) q[1];
sx q[1];
rz(0.47902963) q[1];
rz(-pi) q[2];
rz(0.30626014) q[3];
sx q[3];
rz(-2.2420886) q[3];
sx q[3];
rz(-1.1989532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.11015636) q[2];
sx q[2];
rz(-1.0855731) q[2];
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
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.095117005) q[0];
sx q[0];
rz(-1.7549544) q[0];
sx q[0];
rz(1.1412096) q[0];
rz(-1.8449239) q[1];
sx q[1];
rz(-0.43162391) q[1];
sx q[1];
rz(-0.30532125) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5019708) q[0];
sx q[0];
rz(-1.0143711) q[0];
sx q[0];
rz(1.9407546) q[0];
x q[1];
rz(-2.4112169) q[2];
sx q[2];
rz(-0.59130284) q[2];
sx q[2];
rz(3.0808466) q[2];
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
rz(-2.3577945) q[3];
sx q[3];
rz(-0.67516203) q[3];
sx q[3];
rz(-0.53962196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2083464) q[2];
sx q[2];
rz(-1.5657921) q[2];
sx q[2];
rz(-2.7944881) q[2];
rz(0.38337213) q[3];
sx q[3];
rz(-0.88630668) q[3];
sx q[3];
rz(2.0736407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7326996) q[0];
sx q[0];
rz(-2.2786031) q[0];
sx q[0];
rz(-0.53891671) q[0];
rz(-1.4096889) q[1];
sx q[1];
rz(-0.46202818) q[1];
sx q[1];
rz(-2.7358823) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.021279871) q[0];
sx q[0];
rz(-0.73182523) q[0];
sx q[0];
rz(-2.7243607) q[0];
rz(-pi) q[1];
rz(0.26486764) q[2];
sx q[2];
rz(-2.1200074) q[2];
sx q[2];
rz(2.7153646) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.87863641) q[1];
sx q[1];
rz(-0.43081455) q[1];
sx q[1];
rz(-0.93114914) q[1];
rz(-pi) q[2];
x q[2];
rz(2.898786) q[3];
sx q[3];
rz(-2.3620776) q[3];
sx q[3];
rz(1.8368343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.47200176) q[2];
sx q[2];
rz(-0.80299032) q[2];
sx q[2];
rz(-1.4542788) q[2];
rz(-0.83868319) q[3];
sx q[3];
rz(-1.2861774) q[3];
sx q[3];
rz(1.7680761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31081653) q[0];
sx q[0];
rz(-0.61284471) q[0];
sx q[0];
rz(0.40711656) q[0];
rz(2.282417) q[1];
sx q[1];
rz(-1.5645063) q[1];
sx q[1];
rz(-1.3173332) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.84264) q[0];
sx q[0];
rz(-1.5484122) q[0];
sx q[0];
rz(-1.0277261) q[0];
rz(0.66770422) q[2];
sx q[2];
rz(-2.0470847) q[2];
sx q[2];
rz(-0.27031937) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.085560872) q[1];
sx q[1];
rz(-0.83562682) q[1];
sx q[1];
rz(-0.0037770165) q[1];
x q[2];
rz(-0.42767151) q[3];
sx q[3];
rz(-2.1139522) q[3];
sx q[3];
rz(-3.1205408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(2.4334548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5642501) q[0];
sx q[0];
rz(-1.0471434) q[0];
sx q[0];
rz(2.0136925) q[0];
rz(0.212542) q[1];
sx q[1];
rz(-1.3393341) q[1];
sx q[1];
rz(-1.2316661) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3198277) q[0];
sx q[0];
rz(-1.2152831) q[0];
sx q[0];
rz(2.6000644) q[0];
rz(-pi) q[1];
rz(2.3627794) q[2];
sx q[2];
rz(-2.3328569) q[2];
sx q[2];
rz(0.29701172) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5994508) q[1];
sx q[1];
rz(-2.8891924) q[1];
sx q[1];
rz(2.7155994) q[1];
rz(-1.539174) q[3];
sx q[3];
rz(-2.5368689) q[3];
sx q[3];
rz(-1.5550169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.82905519) q[2];
sx q[2];
rz(-2.2237491) q[2];
sx q[2];
rz(0.60116872) q[2];
rz(2.6230295) q[3];
sx q[3];
rz(-2.8278567) q[3];
sx q[3];
rz(1.4275838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9940051) q[0];
sx q[0];
rz(-0.00031358263) q[0];
sx q[0];
rz(3.0684379) q[0];
rz(-2.1144497) q[1];
sx q[1];
rz(-1.5348397) q[1];
sx q[1];
rz(-2.5491319) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52623977) q[0];
sx q[0];
rz(-0.59229367) q[0];
sx q[0];
rz(0.083239716) q[0];
rz(-2.4249435) q[2];
sx q[2];
rz(-0.71084329) q[2];
sx q[2];
rz(0.58472842) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.91208273) q[1];
sx q[1];
rz(-2.3778937) q[1];
sx q[1];
rz(-2.0565226) q[1];
rz(2.8665364) q[3];
sx q[3];
rz(-1.961879) q[3];
sx q[3];
rz(-1.0712136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.56117326) q[2];
sx q[2];
rz(-1.6757123) q[2];
sx q[2];
rz(-2.6994761) q[2];
rz(-0.90028611) q[3];
sx q[3];
rz(-1.0628275) q[3];
sx q[3];
rz(-0.013307868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1135547) q[0];
sx q[0];
rz(-2.2224764) q[0];
sx q[0];
rz(-1.3046718) q[0];
rz(-1.4470944) q[1];
sx q[1];
rz(-1.1728975) q[1];
sx q[1];
rz(2.5794199) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32414831) q[0];
sx q[0];
rz(-1.2153017) q[0];
sx q[0];
rz(0.72974156) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4116889) q[2];
sx q[2];
rz(-2.6052887) q[2];
sx q[2];
rz(-1.2250587) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8852946) q[1];
sx q[1];
rz(-0.89026272) q[1];
sx q[1];
rz(0.93974944) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1702873) q[3];
sx q[3];
rz(-1.9477961) q[3];
sx q[3];
rz(1.7290982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8582981) q[2];
sx q[2];
rz(-1.7389745) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-1.9853233) q[1];
sx q[1];
rz(3.093739) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82057792) q[0];
sx q[0];
rz(-2.4860959) q[0];
sx q[0];
rz(1.6183678) q[0];
rz(1.1282975) q[2];
sx q[2];
rz(-1.6167165) q[2];
sx q[2];
rz(2.0890582) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1298444) q[1];
sx q[1];
rz(-2.1455892) q[1];
sx q[1];
rz(-2.6227436) q[1];
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
rz(-pi) q[1];
rz(2.4973267) q[2];
sx q[2];
rz(-1.1844119) q[2];
sx q[2];
rz(-2.2052374) q[2];
rz(0.80091536) q[3];
sx q[3];
rz(-2.622486) q[3];
sx q[3];
rz(2.463533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(2.5903493) q[2];
sx q[2];
rz(-1.5521282) q[2];
sx q[2];
rz(2.9350401) q[2];
rz(-0.83947635) q[3];
sx q[3];
rz(-2.6203733) q[3];
sx q[3];
rz(0.93422191) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];