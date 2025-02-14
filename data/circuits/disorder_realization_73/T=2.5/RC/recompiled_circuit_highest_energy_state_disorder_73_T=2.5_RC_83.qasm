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
rz(-1.7492548) q[0];
sx q[0];
rz(5.9049913) q[0];
sx q[0];
rz(9.7752934) q[0];
rz(-2.2502083) q[1];
sx q[1];
rz(-0.79217029) q[1];
sx q[1];
rz(1.9686735) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.354686) q[0];
sx q[0];
rz(-1.7607435) q[0];
sx q[0];
rz(0.19330103) q[0];
x q[1];
rz(-2.1831398) q[2];
sx q[2];
rz(-0.26587379) q[2];
sx q[2];
rz(-1.0622417) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8280331) q[1];
sx q[1];
rz(-2.782722) q[1];
sx q[1];
rz(-0.48017217) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8137742) q[3];
sx q[3];
rz(-0.49428764) q[3];
sx q[3];
rz(-2.0959299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9709836) q[2];
sx q[2];
rz(-1.2759408) q[2];
sx q[2];
rz(-0.21273908) q[2];
rz(-3.0563266) q[3];
sx q[3];
rz(-1.250896) q[3];
sx q[3];
rz(1.8703478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79412115) q[0];
sx q[0];
rz(-2.329282) q[0];
sx q[0];
rz(2.0982657) q[0];
rz(0.13283816) q[1];
sx q[1];
rz(-0.21505198) q[1];
sx q[1];
rz(-1.1588233) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9143694) q[0];
sx q[0];
rz(-2.717319) q[0];
sx q[0];
rz(-2.7844564) q[0];
x q[1];
rz(1.6558596) q[2];
sx q[2];
rz(-0.16189215) q[2];
sx q[2];
rz(-2.7281103) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0483425) q[1];
sx q[1];
rz(-0.56386098) q[1];
sx q[1];
rz(-2.6287931) q[1];
x q[2];
rz(-0.32907991) q[3];
sx q[3];
rz(-1.1606154) q[3];
sx q[3];
rz(-0.044959294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9449238) q[2];
sx q[2];
rz(-0.66201869) q[2];
sx q[2];
rz(0.72466737) q[2];
rz(-0.66257462) q[3];
sx q[3];
rz(-2.1734838) q[3];
sx q[3];
rz(-0.29964963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.1835943) q[0];
sx q[0];
rz(-1.2677001) q[0];
sx q[0];
rz(-2.2268353) q[0];
rz(0.58685189) q[1];
sx q[1];
rz(-1.4205168) q[1];
sx q[1];
rz(-2.9920726) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6059593) q[0];
sx q[0];
rz(-2.0894755) q[0];
sx q[0];
rz(-0.16196047) q[0];
rz(-1.9376041) q[2];
sx q[2];
rz(-1.9123532) q[2];
sx q[2];
rz(1.2836518) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.65626493) q[1];
sx q[1];
rz(-1.0228923) q[1];
sx q[1];
rz(0.62959558) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.87823509) q[3];
sx q[3];
rz(-0.98537579) q[3];
sx q[3];
rz(-0.4957605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2933423) q[2];
sx q[2];
rz(-1.9363656) q[2];
sx q[2];
rz(-0.80043522) q[2];
rz(-2.6750001) q[3];
sx q[3];
rz(-1.1060017) q[3];
sx q[3];
rz(2.485937) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7993497) q[0];
sx q[0];
rz(-1.8514587) q[0];
sx q[0];
rz(-0.85860646) q[0];
rz(-0.41060064) q[1];
sx q[1];
rz(-0.20315367) q[1];
sx q[1];
rz(0.98791775) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84774524) q[0];
sx q[0];
rz(-1.721707) q[0];
sx q[0];
rz(-1.9035089) q[0];
rz(-pi) q[1];
rz(-2.06591) q[2];
sx q[2];
rz(-1.4637814) q[2];
sx q[2];
rz(2.1008976) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4639484) q[1];
sx q[1];
rz(-1.299674) q[1];
sx q[1];
rz(1.8810924) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8186422) q[3];
sx q[3];
rz(-1.7569555) q[3];
sx q[3];
rz(-1.3398088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.19114384) q[2];
sx q[2];
rz(-1.7473651) q[2];
sx q[2];
rz(-2.6960755) q[2];
rz(-0.71349239) q[3];
sx q[3];
rz(-2.4092509) q[3];
sx q[3];
rz(-2.2215686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54818654) q[0];
sx q[0];
rz(-2.0596518) q[0];
sx q[0];
rz(-1.5997546) q[0];
rz(-1.2708739) q[1];
sx q[1];
rz(-1.7393232) q[1];
sx q[1];
rz(-2.612203) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3893376) q[0];
sx q[0];
rz(-1.9175954) q[0];
sx q[0];
rz(-2.8777468) q[0];
rz(-pi) q[1];
rz(0.61754333) q[2];
sx q[2];
rz(-1.7360592) q[2];
sx q[2];
rz(-2.3338855) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.84698757) q[1];
sx q[1];
rz(-2.5683476) q[1];
sx q[1];
rz(-2.6662988) q[1];
x q[2];
rz(0.1481109) q[3];
sx q[3];
rz(-1.9198196) q[3];
sx q[3];
rz(2.7452041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9993837) q[2];
sx q[2];
rz(-1.1915519) q[2];
sx q[2];
rz(1.8184398) q[2];
rz(0.38149825) q[3];
sx q[3];
rz(-2.0254717) q[3];
sx q[3];
rz(2.748446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8489654) q[0];
sx q[0];
rz(-0.75858527) q[0];
sx q[0];
rz(2.1204156) q[0];
rz(1.5317597) q[1];
sx q[1];
rz(-0.94667089) q[1];
sx q[1];
rz(-2.3381332) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96413104) q[0];
sx q[0];
rz(-1.877907) q[0];
sx q[0];
rz(1.7384647) q[0];
x q[1];
rz(1.4864462) q[2];
sx q[2];
rz(-1.9816895) q[2];
sx q[2];
rz(-1.3162168) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.66113055) q[1];
sx q[1];
rz(-1.4509038) q[1];
sx q[1];
rz(2.5698623) q[1];
rz(1.6709953) q[3];
sx q[3];
rz(-0.9462983) q[3];
sx q[3];
rz(-1.1794832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.44567406) q[2];
sx q[2];
rz(-0.58738223) q[2];
sx q[2];
rz(-0.43061259) q[2];
rz(0.63498354) q[3];
sx q[3];
rz(-3.1228784) q[3];
sx q[3];
rz(1.2682605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1465313) q[0];
sx q[0];
rz(-0.54556161) q[0];
sx q[0];
rz(-1.5437641) q[0];
rz(0.68430463) q[1];
sx q[1];
rz(-1.6938208) q[1];
sx q[1];
rz(0.79944557) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9198624) q[0];
sx q[0];
rz(-1.5761901) q[0];
sx q[0];
rz(-0.07660596) q[0];
x q[1];
rz(-2.2861023) q[2];
sx q[2];
rz(-1.6281307) q[2];
sx q[2];
rz(-2.3279026) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3441678) q[1];
sx q[1];
rz(-2.0312211) q[1];
sx q[1];
rz(-1.5595705) q[1];
x q[2];
rz(2.0518028) q[3];
sx q[3];
rz(-2.3481927) q[3];
sx q[3];
rz(-0.97360669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.53576175) q[2];
sx q[2];
rz(-2.3039218) q[2];
sx q[2];
rz(-2.3568995) q[2];
rz(0.11624087) q[3];
sx q[3];
rz(-1.616547) q[3];
sx q[3];
rz(1.2522662) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6004953) q[0];
sx q[0];
rz(-1.2299812) q[0];
sx q[0];
rz(2.1121209) q[0];
rz(-3.0211499) q[1];
sx q[1];
rz(-1.8414626) q[1];
sx q[1];
rz(-0.57074237) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.127205) q[0];
sx q[0];
rz(-2.528101) q[0];
sx q[0];
rz(-0.86742371) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.71855259) q[2];
sx q[2];
rz(-1.6552123) q[2];
sx q[2];
rz(-0.29044232) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.51008049) q[1];
sx q[1];
rz(-1.2464379) q[1];
sx q[1];
rz(2.0381169) q[1];
rz(2.9597046) q[3];
sx q[3];
rz(-1.2396253) q[3];
sx q[3];
rz(1.5665552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.268078) q[2];
sx q[2];
rz(-2.2669078) q[2];
sx q[2];
rz(2.6178005) q[2];
rz(-1.1526456) q[3];
sx q[3];
rz(-2.5609784) q[3];
sx q[3];
rz(1.4279648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48083392) q[0];
sx q[0];
rz(-1.9402215) q[0];
sx q[0];
rz(0.42385605) q[0];
rz(-2.2747874) q[1];
sx q[1];
rz(-1.024217) q[1];
sx q[1];
rz(2.9383235) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.568726) q[0];
sx q[0];
rz(-2.4857268) q[0];
sx q[0];
rz(0.076506581) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6248943) q[2];
sx q[2];
rz(-1.2612169) q[2];
sx q[2];
rz(-0.54540173) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2937676) q[1];
sx q[1];
rz(-1.2888288) q[1];
sx q[1];
rz(-1.7196697) q[1];
x q[2];
rz(-0.8747845) q[3];
sx q[3];
rz(-0.38988567) q[3];
sx q[3];
rz(2.3830151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0554492) q[2];
sx q[2];
rz(-1.2148427) q[2];
sx q[2];
rz(0.27935585) q[2];
rz(1.2934359) q[3];
sx q[3];
rz(-1.3118298) q[3];
sx q[3];
rz(-1.2156585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-2.9781037) q[0];
sx q[0];
rz(-2.3390529) q[0];
sx q[0];
rz(-1.7247024) q[0];
rz(-2.2143927) q[1];
sx q[1];
rz(-1.5798774) q[1];
sx q[1];
rz(-1.4097811) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15148189) q[0];
sx q[0];
rz(-1.6674433) q[0];
sx q[0];
rz(-0.90630177) q[0];
rz(-pi) q[1];
rz(1.0877092) q[2];
sx q[2];
rz(-1.2521241) q[2];
sx q[2];
rz(-2.0234194) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.77892196) q[1];
sx q[1];
rz(-0.92423981) q[1];
sx q[1];
rz(-2.126466) q[1];
rz(-pi) q[2];
rz(-0.042155592) q[3];
sx q[3];
rz(-1.32844) q[3];
sx q[3];
rz(-2.4740296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7659144) q[2];
sx q[2];
rz(-1.8199074) q[2];
sx q[2];
rz(-2.1709757) q[2];
rz(-0.64540234) q[3];
sx q[3];
rz(-2.161945) q[3];
sx q[3];
rz(2.1360883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8476625) q[0];
sx q[0];
rz(-1.6456589) q[0];
sx q[0];
rz(1.8346067) q[0];
rz(-0.38181276) q[1];
sx q[1];
rz(-1.3419071) q[1];
sx q[1];
rz(0.76795427) q[1];
rz(-1.0316331) q[2];
sx q[2];
rz(-2.3066386) q[2];
sx q[2];
rz(2.5021449) q[2];
rz(-2.6991424) q[3];
sx q[3];
rz(-1.9922602) q[3];
sx q[3];
rz(-0.69666399) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
