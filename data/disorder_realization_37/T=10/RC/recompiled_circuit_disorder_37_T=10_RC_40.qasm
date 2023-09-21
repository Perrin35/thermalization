OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.10387575) q[0];
sx q[0];
rz(-1.9394983) q[0];
sx q[0];
rz(1.9934959) q[0];
rz(-1.8885053) q[1];
sx q[1];
rz(-0.94068599) q[1];
sx q[1];
rz(1.747945) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5210261) q[0];
sx q[0];
rz(-1.5728083) q[0];
sx q[0];
rz(0.30962551) q[0];
x q[1];
rz(0.68140985) q[2];
sx q[2];
rz(-1.0076367) q[2];
sx q[2];
rz(-1.5208706) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4268036) q[1];
sx q[1];
rz(-1.1183294) q[1];
sx q[1];
rz(1.8718375) q[1];
rz(-pi) q[2];
rz(2.0066891) q[3];
sx q[3];
rz(-2.5385058) q[3];
sx q[3];
rz(-3.1325504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6538438) q[2];
sx q[2];
rz(-1.8493435) q[2];
sx q[2];
rz(0.1208819) q[2];
rz(-0.17928784) q[3];
sx q[3];
rz(-0.59569734) q[3];
sx q[3];
rz(-2.9860935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0497465) q[0];
sx q[0];
rz(-2.3738528) q[0];
sx q[0];
rz(-3.0088186) q[0];
rz(-1.4615387) q[1];
sx q[1];
rz(-1.5802054) q[1];
sx q[1];
rz(2.9002088) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7624843) q[0];
sx q[0];
rz(-2.603841) q[0];
sx q[0];
rz(2.3178029) q[0];
x q[1];
rz(-2.3953305) q[2];
sx q[2];
rz(-0.66821874) q[2];
sx q[2];
rz(1.2506739) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.0028210359) q[1];
sx q[1];
rz(-1.4797987) q[1];
sx q[1];
rz(-0.9072733) q[1];
rz(-0.026147141) q[3];
sx q[3];
rz(-1.7215683) q[3];
sx q[3];
rz(-1.025841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1277348) q[2];
sx q[2];
rz(-1.3188136) q[2];
sx q[2];
rz(1.1068809) q[2];
rz(1.3876623) q[3];
sx q[3];
rz(-0.51968402) q[3];
sx q[3];
rz(1.2319516) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36104193) q[0];
sx q[0];
rz(-1.1013958) q[0];
sx q[0];
rz(-2.4011491) q[0];
rz(-0.45117798) q[1];
sx q[1];
rz(-1.8809044) q[1];
sx q[1];
rz(-2.0887451) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7171779) q[0];
sx q[0];
rz(-2.4632235) q[0];
sx q[0];
rz(-2.5413187) q[0];
rz(-2.244433) q[2];
sx q[2];
rz(-1.8067915) q[2];
sx q[2];
rz(1.8434075) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4098674) q[1];
sx q[1];
rz(-1.4508529) q[1];
sx q[1];
rz(0.15323318) q[1];
x q[2];
rz(-2.0009082) q[3];
sx q[3];
rz(-1.5617237) q[3];
sx q[3];
rz(-0.7113925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8911002) q[2];
sx q[2];
rz(-1.6001469) q[2];
sx q[2];
rz(-2.5946674) q[2];
rz(0.28918239) q[3];
sx q[3];
rz(-0.5368036) q[3];
sx q[3];
rz(-0.012399013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1716487) q[0];
sx q[0];
rz(-2.8778853) q[0];
sx q[0];
rz(-1.7893715) q[0];
rz(2.9317454) q[1];
sx q[1];
rz(-2.4317957) q[1];
sx q[1];
rz(-0.23637493) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.198092) q[0];
sx q[0];
rz(-2.8236832) q[0];
sx q[0];
rz(-2.148118) q[0];
rz(2.1430074) q[2];
sx q[2];
rz(-2.5479655) q[2];
sx q[2];
rz(1.3015391) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.826556) q[1];
sx q[1];
rz(-0.16780014) q[1];
sx q[1];
rz(1.5312974) q[1];
rz(-2.866719) q[3];
sx q[3];
rz(-2.636424) q[3];
sx q[3];
rz(1.2884017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9005047) q[2];
sx q[2];
rz(-2.1630478) q[2];
sx q[2];
rz(0.55348712) q[2];
rz(2.2262946) q[3];
sx q[3];
rz(-1.2585636) q[3];
sx q[3];
rz(-2.4966911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6699162) q[0];
sx q[0];
rz(-1.8138509) q[0];
sx q[0];
rz(2.4374403) q[0];
rz(2.0856805) q[1];
sx q[1];
rz(-2.6864955) q[1];
sx q[1];
rz(0.59590894) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51673698) q[0];
sx q[0];
rz(-1.7863818) q[0];
sx q[0];
rz(0.09341021) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4679568) q[2];
sx q[2];
rz(-1.4623702) q[2];
sx q[2];
rz(0.078439586) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6223645) q[1];
sx q[1];
rz(-1.16751) q[1];
sx q[1];
rz(1.3695903) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7704352) q[3];
sx q[3];
rz(-1.3988929) q[3];
sx q[3];
rz(2.2558444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.492505) q[2];
sx q[2];
rz(-2.0272144) q[2];
sx q[2];
rz(-2.5904783) q[2];
rz(2.9344432) q[3];
sx q[3];
rz(-1.3144349) q[3];
sx q[3];
rz(-1.5392039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9838487) q[0];
sx q[0];
rz(-2.284323) q[0];
sx q[0];
rz(-0.95170784) q[0];
rz(-0.47479182) q[1];
sx q[1];
rz(-1.1750849) q[1];
sx q[1];
rz(-0.29528433) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8897032) q[0];
sx q[0];
rz(-1.1767052) q[0];
sx q[0];
rz(0.15051145) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4163383) q[2];
sx q[2];
rz(-1.2631577) q[2];
sx q[2];
rz(-2.6442106) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.59776781) q[1];
sx q[1];
rz(-0.73111594) q[1];
sx q[1];
rz(0.75146971) q[1];
x q[2];
rz(3.0764334) q[3];
sx q[3];
rz(-1.061073) q[3];
sx q[3];
rz(-0.62419696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7513912) q[2];
sx q[2];
rz(-1.3662806) q[2];
sx q[2];
rz(-0.93377101) q[2];
rz(-1.4592524) q[3];
sx q[3];
rz(-1.3542342) q[3];
sx q[3];
rz(-1.4565844) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.075832531) q[0];
sx q[0];
rz(-1.1446784) q[0];
sx q[0];
rz(2.899535) q[0];
rz(2.4767955) q[1];
sx q[1];
rz(-1.2882065) q[1];
sx q[1];
rz(0.40329969) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4664073) q[0];
sx q[0];
rz(-1.2012321) q[0];
sx q[0];
rz(-0.36375605) q[0];
x q[1];
rz(2.7935739) q[2];
sx q[2];
rz(-0.70901477) q[2];
sx q[2];
rz(-1.6955171) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.62337263) q[1];
sx q[1];
rz(-1.5689335) q[1];
sx q[1];
rz(1.7809479) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3808448) q[3];
sx q[3];
rz(-1.5459832) q[3];
sx q[3];
rz(0.97770377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.91281259) q[2];
sx q[2];
rz(-2.0704724) q[2];
sx q[2];
rz(-2.7071803) q[2];
rz(-1.0007535) q[3];
sx q[3];
rz(-2.775511) q[3];
sx q[3];
rz(-0.51030695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0077165724) q[0];
sx q[0];
rz(-3.1354597) q[0];
sx q[0];
rz(-0.49466053) q[0];
rz(-1.6330632) q[1];
sx q[1];
rz(-1.6260908) q[1];
sx q[1];
rz(2.5411434) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4279815) q[0];
sx q[0];
rz(-0.62987721) q[0];
sx q[0];
rz(2.4226818) q[0];
rz(-2.0958488) q[2];
sx q[2];
rz(-1.5480435) q[2];
sx q[2];
rz(-2.3843228) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3214896) q[1];
sx q[1];
rz(-0.21807018) q[1];
sx q[1];
rz(0.58184187) q[1];
rz(-pi) q[2];
x q[2];
rz(0.59372254) q[3];
sx q[3];
rz(-0.36168081) q[3];
sx q[3];
rz(2.7152674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1373458) q[2];
sx q[2];
rz(-2.1477951) q[2];
sx q[2];
rz(-3.0252769) q[2];
rz(-2.7159193) q[3];
sx q[3];
rz(-0.95723546) q[3];
sx q[3];
rz(11*pi/12) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9882934) q[0];
sx q[0];
rz(-2.9635552) q[0];
sx q[0];
rz(1.6631888) q[0];
rz(2.2019745) q[1];
sx q[1];
rz(-1.3213108) q[1];
sx q[1];
rz(0.41752648) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4672887) q[0];
sx q[0];
rz(-1.0397362) q[0];
sx q[0];
rz(3.0896316) q[0];
rz(-pi) q[1];
x q[1];
rz(1.19403) q[2];
sx q[2];
rz(-0.28290877) q[2];
sx q[2];
rz(-0.49809581) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6431943) q[1];
sx q[1];
rz(-0.25499757) q[1];
sx q[1];
rz(-0.71868371) q[1];
rz(-pi) q[2];
rz(-2.0017654) q[3];
sx q[3];
rz(-1.7472072) q[3];
sx q[3];
rz(-0.21709066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6614762) q[2];
sx q[2];
rz(-0.63642234) q[2];
sx q[2];
rz(-0.49368668) q[2];
rz(2.4168329) q[3];
sx q[3];
rz(-1.1586435) q[3];
sx q[3];
rz(-0.31989583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17393728) q[0];
sx q[0];
rz(-0.65615654) q[0];
sx q[0];
rz(-2.4560112) q[0];
rz(2.8441692) q[1];
sx q[1];
rz(-2.90459) q[1];
sx q[1];
rz(1.1313653) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3216074) q[0];
sx q[0];
rz(-1.2258343) q[0];
sx q[0];
rz(2.5489775) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4360043) q[2];
sx q[2];
rz(-1.0552647) q[2];
sx q[2];
rz(2.3956092) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.28045052) q[1];
sx q[1];
rz(-0.44499731) q[1];
sx q[1];
rz(-0.54235561) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4990436) q[3];
sx q[3];
rz(-0.71392871) q[3];
sx q[3];
rz(0.16845265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.82548213) q[2];
sx q[2];
rz(-1.1967412) q[2];
sx q[2];
rz(-0.70739174) q[2];
rz(-0.80983821) q[3];
sx q[3];
rz(-0.67088586) q[3];
sx q[3];
rz(-0.18856089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0512882) q[0];
sx q[0];
rz(-1.123883) q[0];
sx q[0];
rz(-0.71258769) q[0];
rz(-2.9293625) q[1];
sx q[1];
rz(-1.4490912) q[1];
sx q[1];
rz(2.6279411) q[1];
rz(-0.99954188) q[2];
sx q[2];
rz(-2.0959601) q[2];
sx q[2];
rz(-3.0053896) q[2];
rz(0.57701941) q[3];
sx q[3];
rz(-0.95303017) q[3];
sx q[3];
rz(-0.87458761) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];