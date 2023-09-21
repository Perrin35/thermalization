OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.22566158) q[0];
sx q[0];
rz(7.1516501) q[0];
sx q[0];
rz(9.2317543) q[0];
rz(-1.9999737) q[1];
sx q[1];
rz(3.5715754) q[1];
sx q[1];
rz(6.9663098) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30259351) q[0];
sx q[0];
rz(-1.6380881) q[0];
sx q[0];
rz(1.6291314) q[0];
rz(-0.24129759) q[2];
sx q[2];
rz(-1.9120875) q[2];
sx q[2];
rz(-1.8482006) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7072304) q[1];
sx q[1];
rz(-1.1303567) q[1];
sx q[1];
rz(0.67727725) q[1];
rz(-pi) q[2];
rz(-1.0115511) q[3];
sx q[3];
rz(-2.4992001) q[3];
sx q[3];
rz(-1.1005644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0156988) q[2];
sx q[2];
rz(-1.761972) q[2];
sx q[2];
rz(1.0985628) q[2];
rz(-1.0788318) q[3];
sx q[3];
rz(-0.94509411) q[3];
sx q[3];
rz(2.0400955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1612448) q[0];
sx q[0];
rz(-2.8256567) q[0];
sx q[0];
rz(-2.9336477) q[0];
rz(0.57693276) q[1];
sx q[1];
rz(-0.88795841) q[1];
sx q[1];
rz(1.6764486) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10509051) q[0];
sx q[0];
rz(-2.4173792) q[0];
sx q[0];
rz(3.1378531) q[0];
rz(-pi) q[1];
rz(1.0826153) q[2];
sx q[2];
rz(-1.3435875) q[2];
sx q[2];
rz(2.8607334) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2966753) q[1];
sx q[1];
rz(-2.1293318) q[1];
sx q[1];
rz(1.0122453) q[1];
rz(-1.2259237) q[3];
sx q[3];
rz(-0.85540918) q[3];
sx q[3];
rz(1.2777002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0186105) q[2];
sx q[2];
rz(-0.98615042) q[2];
sx q[2];
rz(3.0318276) q[2];
rz(-2.5189853) q[3];
sx q[3];
rz(-0.37125769) q[3];
sx q[3];
rz(-1.4573147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(2.1784172) q[0];
sx q[0];
rz(-0.85997471) q[0];
sx q[0];
rz(-0.51613581) q[0];
rz(-0.57488817) q[1];
sx q[1];
rz(-0.92620414) q[1];
sx q[1];
rz(0.80054545) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0101937) q[0];
sx q[0];
rz(-2.7079765) q[0];
sx q[0];
rz(-1.2244768) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2764552) q[2];
sx q[2];
rz(-1.0312928) q[2];
sx q[2];
rz(-1.7973961) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.264818) q[1];
sx q[1];
rz(-0.41300981) q[1];
sx q[1];
rz(-3.1289711) q[1];
x q[2];
rz(2.27182) q[3];
sx q[3];
rz(-1.3372034) q[3];
sx q[3];
rz(-1.1886532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4642554) q[2];
sx q[2];
rz(-2.8223473) q[2];
sx q[2];
rz(-1.7819972) q[2];
rz(-2.1740186) q[3];
sx q[3];
rz(-1.8680957) q[3];
sx q[3];
rz(1.4250071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1490705) q[0];
sx q[0];
rz(-1.2620121) q[0];
sx q[0];
rz(2.676679) q[0];
rz(2.7930296) q[1];
sx q[1];
rz(-2.8788853) q[1];
sx q[1];
rz(-2.0565313) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5269055) q[0];
sx q[0];
rz(-0.58262107) q[0];
sx q[0];
rz(2.7132062) q[0];
rz(-pi) q[1];
rz(-1.2041353) q[2];
sx q[2];
rz(-1.3385834) q[2];
sx q[2];
rz(2.1249078) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9492053) q[1];
sx q[1];
rz(-2.9678223) q[1];
sx q[1];
rz(-0.6178356) q[1];
rz(0.55235858) q[3];
sx q[3];
rz(-0.68867749) q[3];
sx q[3];
rz(-2.3104582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1365635) q[2];
sx q[2];
rz(-2.0932784) q[2];
sx q[2];
rz(0.68112779) q[2];
rz(2.629771) q[3];
sx q[3];
rz(-2.8183283) q[3];
sx q[3];
rz(-0.26369035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8624449) q[0];
sx q[0];
rz(-2.6036766) q[0];
sx q[0];
rz(-1.408668) q[0];
rz(-2.7092343) q[1];
sx q[1];
rz(-2.2996348) q[1];
sx q[1];
rz(-0.98168215) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3571346) q[0];
sx q[0];
rz(-1.1103837) q[0];
sx q[0];
rz(0.61607342) q[0];
x q[1];
rz(-0.15820299) q[2];
sx q[2];
rz(-1.8100097) q[2];
sx q[2];
rz(-1.6895134) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9202068) q[1];
sx q[1];
rz(-2.033794) q[1];
sx q[1];
rz(1.4057926) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0675689) q[3];
sx q[3];
rz(-0.97462666) q[3];
sx q[3];
rz(1.0130458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0354054) q[2];
sx q[2];
rz(-1.4865439) q[2];
sx q[2];
rz(0.48941082) q[2];
rz(-2.126746) q[3];
sx q[3];
rz(-0.5265407) q[3];
sx q[3];
rz(-1.813252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6569825) q[0];
sx q[0];
rz(-2.9841612) q[0];
sx q[0];
rz(2.7691675) q[0];
rz(1.3308446) q[1];
sx q[1];
rz(-2.1061888) q[1];
sx q[1];
rz(-2.9763124) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80711354) q[0];
sx q[0];
rz(-1.4757336) q[0];
sx q[0];
rz(-1.6001742) q[0];
rz(-pi) q[1];
rz(-1.500962) q[2];
sx q[2];
rz(-2.2572821) q[2];
sx q[2];
rz(1.7718466) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.11152553) q[1];
sx q[1];
rz(-1.0377874) q[1];
sx q[1];
rz(-1.0982606) q[1];
x q[2];
rz(1.1214439) q[3];
sx q[3];
rz(-1.7835788) q[3];
sx q[3];
rz(1.3495812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2281987) q[2];
sx q[2];
rz(-2.1134351) q[2];
sx q[2];
rz(1.4432663) q[2];
rz(-2.7741487) q[3];
sx q[3];
rz(-1.2747217) q[3];
sx q[3];
rz(2.929556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.7845602) q[0];
sx q[0];
rz(-2.050188) q[0];
sx q[0];
rz(0.34926397) q[0];
rz(-2.3941984) q[1];
sx q[1];
rz(-0.29577297) q[1];
sx q[1];
rz(0.73648891) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9518785) q[0];
sx q[0];
rz(-1.6197546) q[0];
sx q[0];
rz(0.017107054) q[0];
rz(0.94524224) q[2];
sx q[2];
rz(-0.40996273) q[2];
sx q[2];
rz(-1.2298825) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.10662096) q[1];
sx q[1];
rz(-1.7502516) q[1];
sx q[1];
rz(2.3751395) q[1];
rz(-0.50076671) q[3];
sx q[3];
rz(-1.1874677) q[3];
sx q[3];
rz(-0.70481833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.13654576) q[2];
sx q[2];
rz(-1.2282635) q[2];
sx q[2];
rz(-2.6100256) q[2];
rz(2.452204) q[3];
sx q[3];
rz(-1.4644943) q[3];
sx q[3];
rz(1.846107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(0.078995973) q[0];
sx q[0];
rz(-1.9364708) q[0];
sx q[0];
rz(-1.8435562) q[0];
rz(-2.334306) q[1];
sx q[1];
rz(-1.9629982) q[1];
sx q[1];
rz(-2.2198026) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73496504) q[0];
sx q[0];
rz(-1.4359183) q[0];
sx q[0];
rz(-2.448003) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8838896) q[2];
sx q[2];
rz(-0.76258341) q[2];
sx q[2];
rz(1.2459754) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.31729749) q[1];
sx q[1];
rz(-1.3166787) q[1];
sx q[1];
rz(1.3586678) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2790518) q[3];
sx q[3];
rz(-1.4338014) q[3];
sx q[3];
rz(-2.9555637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2395997) q[2];
sx q[2];
rz(-2.5805876) q[2];
sx q[2];
rz(1.9699338) q[2];
rz(1.7840067) q[3];
sx q[3];
rz(-1.6849018) q[3];
sx q[3];
rz(1.6931036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(2.8885324) q[0];
sx q[0];
rz(-1.9655515) q[0];
sx q[0];
rz(1.4755479) q[0];
rz(-1.5400003) q[1];
sx q[1];
rz(-1.4614636) q[1];
sx q[1];
rz(-0.67970651) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98201671) q[0];
sx q[0];
rz(-1.3681108) q[0];
sx q[0];
rz(2.6152339) q[0];
x q[1];
rz(-0.18025132) q[2];
sx q[2];
rz(-2.2580574) q[2];
sx q[2];
rz(0.3790516) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.45021536) q[1];
sx q[1];
rz(-2.2274349) q[1];
sx q[1];
rz(-2.0037829) q[1];
rz(-pi) q[2];
rz(-0.71984843) q[3];
sx q[3];
rz(-0.82092972) q[3];
sx q[3];
rz(1.6380701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1264964) q[2];
sx q[2];
rz(-1.6588147) q[2];
sx q[2];
rz(-0.91040197) q[2];
rz(0.67534584) q[3];
sx q[3];
rz(-2.2048435) q[3];
sx q[3];
rz(0.85062406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0891721) q[0];
sx q[0];
rz(-1.9071254) q[0];
sx q[0];
rz(-2.4269379) q[0];
rz(2.4275298) q[1];
sx q[1];
rz(-2.1866182) q[1];
sx q[1];
rz(1.9649327) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6274174) q[0];
sx q[0];
rz(-1.6192993) q[0];
sx q[0];
rz(1.4072627) q[0];
rz(-pi) q[1];
rz(-0.93675905) q[2];
sx q[2];
rz(-2.0271795) q[2];
sx q[2];
rz(2.2697743) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9709819) q[1];
sx q[1];
rz(-0.5628399) q[1];
sx q[1];
rz(-1.4303722) q[1];
x q[2];
rz(2.6305466) q[3];
sx q[3];
rz(-1.647445) q[3];
sx q[3];
rz(-0.49728909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8979793) q[2];
sx q[2];
rz(-1.4695797) q[2];
sx q[2];
rz(1.127355) q[2];
rz(-0.35774287) q[3];
sx q[3];
rz(-2.2704411) q[3];
sx q[3];
rz(2.0991142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42416278) q[0];
sx q[0];
rz(-1.9308199) q[0];
sx q[0];
rz(0.45817026) q[0];
rz(-0.39623109) q[1];
sx q[1];
rz(-3.1165262) q[1];
sx q[1];
rz(-2.810626) q[1];
rz(-1.8489807) q[2];
sx q[2];
rz(-2.2879911) q[2];
sx q[2];
rz(0.0030980274) q[2];
rz(2.7908294) q[3];
sx q[3];
rz(-1.5784932) q[3];
sx q[3];
rz(-2.2898883) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
