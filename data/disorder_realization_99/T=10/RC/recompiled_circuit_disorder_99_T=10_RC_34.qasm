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
rz(-0.7926994) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22098955) q[0];
sx q[0];
rz(-1.594829) q[0];
sx q[0];
rz(-0.8997957) q[0];
rz(1.3408469) q[2];
sx q[2];
rz(-1.0385498) q[2];
sx q[2];
rz(-1.4872273) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.76268643) q[1];
sx q[1];
rz(-2.0204117) q[1];
sx q[1];
rz(-3.0607037) q[1];
rz(-pi) q[2];
rz(-0.80134942) q[3];
sx q[3];
rz(-2.1270463) q[3];
sx q[3];
rz(2.706561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.67316002) q[2];
sx q[2];
rz(-2.1940553) q[2];
sx q[2];
rz(-2.6426278) q[2];
rz(-0.4237825) q[3];
sx q[3];
rz(-1.2200004) q[3];
sx q[3];
rz(0.012714816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38133165) q[0];
sx q[0];
rz(-1.6853764) q[0];
sx q[0];
rz(-0.99200845) q[0];
rz(2.965062) q[1];
sx q[1];
rz(-2.9361528) q[1];
sx q[1];
rz(0.36546779) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27643833) q[0];
sx q[0];
rz(-1.3197474) q[0];
sx q[0];
rz(-2.6692078) q[0];
x q[1];
rz(-0.59132691) q[2];
sx q[2];
rz(-2.1543243) q[2];
sx q[2];
rz(-3.025625) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(2.9404853) q[3];
sx q[3];
rz(-1.4521986) q[3];
sx q[3];
rz(3.0191641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3084597) q[2];
sx q[2];
rz(-0.59811991) q[2];
sx q[2];
rz(-2.0642521) q[2];
rz(-2.9789553) q[3];
sx q[3];
rz(-1.2626303) q[3];
sx q[3];
rz(1.8796273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
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
rz(0.61082947) q[0];
rz(-1.4340596) q[1];
sx q[1];
rz(-1.779498) q[1];
sx q[1];
rz(0.45062137) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0656933) q[0];
sx q[0];
rz(-1.6487299) q[0];
sx q[0];
rz(0.033228544) q[0];
rz(-pi) q[1];
rz(-0.89183715) q[2];
sx q[2];
rz(-0.46687451) q[2];
sx q[2];
rz(2.4613949) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.33774099) q[1];
sx q[1];
rz(-0.9461113) q[1];
sx q[1];
rz(1.9546263) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-3.0314363) q[2];
sx q[2];
rz(-2.0560196) q[2];
sx q[2];
rz(1.191656) q[2];
rz(-0.84447652) q[3];
sx q[3];
rz(-2.8561487) q[3];
sx q[3];
rz(-0.58117956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.095117005) q[0];
sx q[0];
rz(-1.3866383) q[0];
sx q[0];
rz(-1.1412096) q[0];
rz(1.8449239) q[1];
sx q[1];
rz(-0.43162391) q[1];
sx q[1];
rz(0.30532125) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1353377) q[0];
sx q[0];
rz(-2.4843785) q[0];
sx q[0];
rz(0.52657907) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4112169) q[2];
sx q[2];
rz(-0.59130284) q[2];
sx q[2];
rz(-3.0808466) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0381283) q[1];
sx q[1];
rz(-0.80059073) q[1];
sx q[1];
rz(-0.14726463) q[1];
rz(-pi) q[2];
rz(-1.0563072) q[3];
sx q[3];
rz(-2.0293651) q[3];
sx q[3];
rz(1.6955299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9332463) q[2];
sx q[2];
rz(-1.5657921) q[2];
sx q[2];
rz(-0.34710458) q[2];
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
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7326996) q[0];
sx q[0];
rz(-2.2786031) q[0];
sx q[0];
rz(2.6026759) q[0];
rz(-1.4096889) q[1];
sx q[1];
rz(-2.6795645) q[1];
sx q[1];
rz(2.7358823) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.021279871) q[0];
sx q[0];
rz(-0.73182523) q[0];
sx q[0];
rz(-0.41723199) q[0];
rz(-pi) q[1];
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
rz(0.09771422) q[1];
sx q[1];
rz(-1.3188625) q[1];
sx q[1];
rz(1.9240727) q[1];
x q[2];
rz(2.3769576) q[3];
sx q[3];
rz(-1.4009762) q[3];
sx q[3];
rz(3.0499383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(1.3735166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(0.85917568) q[1];
sx q[1];
rz(-1.5770864) q[1];
sx q[1];
rz(1.8242594) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.84264) q[0];
sx q[0];
rz(-1.5931804) q[0];
sx q[0];
rz(2.1138666) q[0];
rz(-pi) q[1];
x q[1];
rz(0.69465722) q[2];
sx q[2];
rz(-2.3431871) q[2];
sx q[2];
rz(-2.3677804) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4827022) q[1];
sx q[1];
rz(-1.5735978) q[1];
sx q[1];
rz(-0.83562327) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7139211) q[3];
sx q[3];
rz(-2.1139522) q[3];
sx q[3];
rz(-3.1205408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0417827) q[2];
sx q[2];
rz(-2.2618497) q[2];
sx q[2];
rz(-2.2463634) q[2];
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
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5642501) q[0];
sx q[0];
rz(-2.0944493) q[0];
sx q[0];
rz(1.1279001) q[0];
rz(2.9290507) q[1];
sx q[1];
rz(-1.3393341) q[1];
sx q[1];
rz(-1.9099265) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82176498) q[0];
sx q[0];
rz(-1.2152831) q[0];
sx q[0];
rz(-0.54152821) q[0];
x q[1];
rz(0.93630837) q[2];
sx q[2];
rz(-1.0299183) q[2];
sx q[2];
rz(2.4782431) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.61464804) q[1];
sx q[1];
rz(-1.6741747) q[1];
sx q[1];
rz(2.9109216) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6024186) q[3];
sx q[3];
rz(-2.5368689) q[3];
sx q[3];
rz(1.5865758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3125375) q[2];
sx q[2];
rz(-2.2237491) q[2];
sx q[2];
rz(2.5404239) q[2];
rz(2.6230295) q[3];
sx q[3];
rz(-0.31373599) q[3];
sx q[3];
rz(-1.4275838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9940051) q[0];
sx q[0];
rz(-3.1412791) q[0];
sx q[0];
rz(0.073154733) q[0];
rz(1.0271429) q[1];
sx q[1];
rz(-1.5348397) q[1];
sx q[1];
rz(0.59246078) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5151278) q[0];
sx q[0];
rz(-0.98083064) q[0];
sx q[0];
rz(-1.5149087) q[0];
rz(-2.5657797) q[2];
sx q[2];
rz(-1.1278707) q[2];
sx q[2];
rz(1.5720313) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.91208273) q[1];
sx q[1];
rz(-2.3778937) q[1];
sx q[1];
rz(1.0850701) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.153272) q[3];
sx q[3];
rz(-2.6675825) q[3];
sx q[3];
rz(1.7085027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.56117326) q[2];
sx q[2];
rz(-1.6757123) q[2];
sx q[2];
rz(-2.6994761) q[2];
rz(0.90028611) q[3];
sx q[3];
rz(-2.0787652) q[3];
sx q[3];
rz(3.1282848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
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
rz(1.6944983) q[1];
sx q[1];
rz(-1.9686952) q[1];
sx q[1];
rz(-2.5794199) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5484555) q[0];
sx q[0];
rz(-2.2459185) q[0];
sx q[0];
rz(1.1086585) q[0];
rz(1.4116889) q[2];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2562981) q[1];
sx q[1];
rz(-2.2513299) q[1];
sx q[1];
rz(0.93974944) q[1];
rz(-0.4060679) q[3];
sx q[3];
rz(-1.9417524) q[3];
sx q[3];
rz(-0.0036811034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.28329453) q[2];
sx q[2];
rz(-1.4026182) q[2];
sx q[2];
rz(0.98158681) q[2];
rz(0.0043407241) q[3];
sx q[3];
rz(-1.0467632) q[3];
sx q[3];
rz(0.59345746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8269862) q[0];
sx q[0];
rz(-1.9208603) q[0];
sx q[0];
rz(3.0944968) q[0];
rz(1.3866407) q[1];
sx q[1];
rz(-1.9853233) q[1];
sx q[1];
rz(3.093739) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3809966) q[0];
sx q[0];
rz(-2.2254235) q[0];
sx q[0];
rz(-3.1050443) q[0];
x q[1];
rz(2.0132952) q[2];
sx q[2];
rz(-1.5248761) q[2];
sx q[2];
rz(-1.0525345) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1298444) q[1];
sx q[1];
rz(-2.1455892) q[1];
sx q[1];
rz(2.6227436) q[1];
rz(-pi) q[2];
rz(-0.1286653) q[3];
sx q[3];
rz(-2.7358958) q[3];
sx q[3];
rz(-0.88276827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.64426595) q[2];
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
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
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
rz(-0.55124333) q[2];
sx q[2];
rz(-1.5521282) q[2];
sx q[2];
rz(2.9350401) q[2];
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
