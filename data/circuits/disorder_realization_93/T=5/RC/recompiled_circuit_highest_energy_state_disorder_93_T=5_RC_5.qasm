OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6858202) q[0];
sx q[0];
rz(-2.0210285) q[0];
sx q[0];
rz(-1.4135345) q[0];
rz(-2.7780374) q[1];
sx q[1];
rz(-1.8283365) q[1];
sx q[1];
rz(0.29120905) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.04714) q[0];
sx q[0];
rz(-1.5908851) q[0];
sx q[0];
rz(-1.8909341) q[0];
x q[1];
rz(1.2313263) q[2];
sx q[2];
rz(-1.4541153) q[2];
sx q[2];
rz(-1.1543903) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8785581) q[1];
sx q[1];
rz(-0.59459762) q[1];
sx q[1];
rz(-2.8697398) q[1];
x q[2];
rz(-2.5194571) q[3];
sx q[3];
rz(-0.91968007) q[3];
sx q[3];
rz(1.5762005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.0061965813) q[2];
sx q[2];
rz(-1.5781032) q[2];
sx q[2];
rz(3.1175933) q[2];
rz(-2.7855347) q[3];
sx q[3];
rz(-0.67432109) q[3];
sx q[3];
rz(-3.1168028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.085670797) q[0];
sx q[0];
rz(-2.4503158) q[0];
sx q[0];
rz(2.5066277) q[0];
rz(2.3786646) q[1];
sx q[1];
rz(-1.678391) q[1];
sx q[1];
rz(-0.91711226) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2208495) q[0];
sx q[0];
rz(-1.5524327) q[0];
sx q[0];
rz(1.5361274) q[0];
rz(-2.6659662) q[2];
sx q[2];
rz(-2.1948819) q[2];
sx q[2];
rz(-1.9204572) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0113127) q[1];
sx q[1];
rz(-0.44751274) q[1];
sx q[1];
rz(1.391973) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2767601) q[3];
sx q[3];
rz(-0.64084478) q[3];
sx q[3];
rz(-1.3230192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.94747535) q[2];
sx q[2];
rz(-0.52299356) q[2];
sx q[2];
rz(0.73327649) q[2];
rz(1.0670079) q[3];
sx q[3];
rz(-1.4209483) q[3];
sx q[3];
rz(0.1799306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5697524) q[0];
sx q[0];
rz(-1.1886007) q[0];
sx q[0];
rz(0.52823129) q[0];
rz(-2.275548) q[1];
sx q[1];
rz(-1.3287611) q[1];
sx q[1];
rz(-0.51308647) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.899174) q[0];
sx q[0];
rz(-1.1326101) q[0];
sx q[0];
rz(3.0946381) q[0];
rz(-pi) q[1];
rz(1.971286) q[2];
sx q[2];
rz(-2.4520564) q[2];
sx q[2];
rz(-0.51338085) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0167488) q[1];
sx q[1];
rz(-1.4433772) q[1];
sx q[1];
rz(-1.9215073) q[1];
rz(-0.052185197) q[3];
sx q[3];
rz(-0.5002678) q[3];
sx q[3];
rz(-2.1310998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.036006007) q[2];
sx q[2];
rz(-2.0522108) q[2];
sx q[2];
rz(-3.1028683) q[2];
rz(0.45768467) q[3];
sx q[3];
rz(-2.3364412) q[3];
sx q[3];
rz(1.800644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.7078581) q[0];
sx q[0];
rz(-0.41446328) q[0];
sx q[0];
rz(-2.0204954) q[0];
rz(2.3838249) q[1];
sx q[1];
rz(-1.517375) q[1];
sx q[1];
rz(2.4574492) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30267495) q[0];
sx q[0];
rz(-0.7644628) q[0];
sx q[0];
rz(1.27728) q[0];
rz(-pi) q[1];
x q[1];
rz(1.563638) q[2];
sx q[2];
rz(-1.560692) q[2];
sx q[2];
rz(1.2105389) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.63057478) q[1];
sx q[1];
rz(-1.2666432) q[1];
sx q[1];
rz(-1.263878) q[1];
rz(-pi) q[2];
rz(0.512796) q[3];
sx q[3];
rz(-1.0362451) q[3];
sx q[3];
rz(-3.0844546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6571558) q[2];
sx q[2];
rz(-1.6146722) q[2];
sx q[2];
rz(-0.26051513) q[2];
rz(-2.303463) q[3];
sx q[3];
rz(-1.1476436) q[3];
sx q[3];
rz(0.014613541) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92358661) q[0];
sx q[0];
rz(-1.3777233) q[0];
sx q[0];
rz(2.6309784) q[0];
rz(1.249373) q[1];
sx q[1];
rz(-1.9793972) q[1];
sx q[1];
rz(-1.6872663) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4076865) q[0];
sx q[0];
rz(-2.9260194) q[0];
sx q[0];
rz(3.0093497) q[0];
rz(-pi) q[1];
rz(0.32938214) q[2];
sx q[2];
rz(-2.9940146) q[2];
sx q[2];
rz(-1.3775228) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0328503) q[1];
sx q[1];
rz(-1.8287484) q[1];
sx q[1];
rz(2.5775419) q[1];
rz(2.2959443) q[3];
sx q[3];
rz(-0.21580869) q[3];
sx q[3];
rz(-0.25419054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5877567) q[2];
sx q[2];
rz(-2.7601056) q[2];
sx q[2];
rz(2.7510263) q[2];
rz(0.25911123) q[3];
sx q[3];
rz(-1.5026389) q[3];
sx q[3];
rz(2.794877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7873586) q[0];
sx q[0];
rz(-2.7095095) q[0];
sx q[0];
rz(-2.4969192) q[0];
rz(1.8395754) q[1];
sx q[1];
rz(-1.527486) q[1];
sx q[1];
rz(2.7929746) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.020418305) q[0];
sx q[0];
rz(-1.114959) q[0];
sx q[0];
rz(-2.0860079) q[0];
x q[1];
rz(-0.63747823) q[2];
sx q[2];
rz(-2.5036252) q[2];
sx q[2];
rz(1.5359985) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.4069566) q[1];
sx q[1];
rz(-2.6183081) q[1];
sx q[1];
rz(2.1914545) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0246234) q[3];
sx q[3];
rz(-1.5364963) q[3];
sx q[3];
rz(-3.0142865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1031441) q[2];
sx q[2];
rz(-0.95547533) q[2];
sx q[2];
rz(3.1032584) q[2];
rz(-0.96758715) q[3];
sx q[3];
rz(-1.7318232) q[3];
sx q[3];
rz(-2.7819395) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2751145) q[0];
sx q[0];
rz(-0.51946467) q[0];
sx q[0];
rz(0.073609784) q[0];
rz(-1.4069125) q[1];
sx q[1];
rz(-2.584447) q[1];
sx q[1];
rz(-0.73307347) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3033111) q[0];
sx q[0];
rz(-1.065863) q[0];
sx q[0];
rz(-0.61793332) q[0];
rz(-pi) q[1];
rz(-0.60150749) q[2];
sx q[2];
rz(-2.2081293) q[2];
sx q[2];
rz(-1.21797) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9736874) q[1];
sx q[1];
rz(-1.2406772) q[1];
sx q[1];
rz(-2.7658505) q[1];
rz(2.9063609) q[3];
sx q[3];
rz(-1.3266272) q[3];
sx q[3];
rz(-1.5926804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0343895) q[2];
sx q[2];
rz(-1.723039) q[2];
sx q[2];
rz(-3.1305195) q[2];
rz(0.30558807) q[3];
sx q[3];
rz(-1.1253076) q[3];
sx q[3];
rz(1.2529469) q[3];
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
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3333862) q[0];
sx q[0];
rz(-2.5337063) q[0];
sx q[0];
rz(-2.407684) q[0];
rz(-0.58087307) q[1];
sx q[1];
rz(-2.1323233) q[1];
sx q[1];
rz(-3.0427921) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89284586) q[0];
sx q[0];
rz(-0.98631421) q[0];
sx q[0];
rz(-0.92846034) q[0];
rz(-pi) q[1];
rz(0.67811857) q[2];
sx q[2];
rz(-1.0762012) q[2];
sx q[2];
rz(1.8910318) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.96896763) q[1];
sx q[1];
rz(-1.8207112) q[1];
sx q[1];
rz(0.45721284) q[1];
rz(0.12998534) q[3];
sx q[3];
rz(-1.49125) q[3];
sx q[3];
rz(1.6787488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.29844478) q[2];
sx q[2];
rz(-1.8567905) q[2];
sx q[2];
rz(1.0996381) q[2];
rz(-3.1067276) q[3];
sx q[3];
rz(-0.74368447) q[3];
sx q[3];
rz(-0.59665027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1140887) q[0];
sx q[0];
rz(-1.9845668) q[0];
sx q[0];
rz(1.1771359) q[0];
rz(1.4741395) q[1];
sx q[1];
rz(-0.93282229) q[1];
sx q[1];
rz(-1.6966381) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6072877) q[0];
sx q[0];
rz(-1.8518847) q[0];
sx q[0];
rz(0.40977733) q[0];
rz(-pi) q[1];
rz(-2.3722052) q[2];
sx q[2];
rz(-2.0035929) q[2];
sx q[2];
rz(-2.3355049) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9517773) q[1];
sx q[1];
rz(-0.51945247) q[1];
sx q[1];
rz(-1.7096667) q[1];
x q[2];
rz(2.8650896) q[3];
sx q[3];
rz(-0.91405896) q[3];
sx q[3];
rz(1.2402676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9745447) q[2];
sx q[2];
rz(-1.5844774) q[2];
sx q[2];
rz(2.9696999) q[2];
rz(2.3520172) q[3];
sx q[3];
rz(-2.004345) q[3];
sx q[3];
rz(1.8062228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0563141) q[0];
sx q[0];
rz(-0.21509898) q[0];
sx q[0];
rz(0.54157448) q[0];
rz(-1.1594634) q[1];
sx q[1];
rz(-1.8340725) q[1];
sx q[1];
rz(0.83311876) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.706382) q[0];
sx q[0];
rz(-2.4738389) q[0];
sx q[0];
rz(-0.45489648) q[0];
rz(-pi) q[1];
rz(-1.2255185) q[2];
sx q[2];
rz(-1.0547027) q[2];
sx q[2];
rz(0.47613019) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.009506735) q[1];
sx q[1];
rz(-0.888538) q[1];
sx q[1];
rz(-2.8481774) q[1];
rz(-pi) q[2];
rz(0.0019440193) q[3];
sx q[3];
rz(-2.2984347) q[3];
sx q[3];
rz(1.9950642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1200166) q[2];
sx q[2];
rz(-0.77777255) q[2];
sx q[2];
rz(-0.94338256) q[2];
rz(1.0776445) q[3];
sx q[3];
rz(-1.5543944) q[3];
sx q[3];
rz(-0.32850346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42644603) q[0];
sx q[0];
rz(-1.0363415) q[0];
sx q[0];
rz(2.8847726) q[0];
rz(2.9019451) q[1];
sx q[1];
rz(-0.9050723) q[1];
sx q[1];
rz(-1.1368652) q[1];
rz(-0.69375615) q[2];
sx q[2];
rz(-1.0508595) q[2];
sx q[2];
rz(-1.6169242) q[2];
rz(-1.8244459) q[3];
sx q[3];
rz(-1.8719049) q[3];
sx q[3];
rz(-0.20897839) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
