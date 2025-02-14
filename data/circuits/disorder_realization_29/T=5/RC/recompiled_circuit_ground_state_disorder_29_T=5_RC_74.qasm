OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.090341181) q[0];
sx q[0];
rz(-1.0426961) q[0];
sx q[0];
rz(-2.8629177) q[0];
rz(-1.9947808) q[1];
sx q[1];
rz(-2.6169701) q[1];
sx q[1];
rz(-1.8705179) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2575114) q[0];
sx q[0];
rz(-1.3102828) q[0];
sx q[0];
rz(1.5519014) q[0];
x q[1];
rz(1.4393629) q[2];
sx q[2];
rz(-1.4214941) q[2];
sx q[2];
rz(2.1367913) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.179175) q[1];
sx q[1];
rz(-1.6500874) q[1];
sx q[1];
rz(0.85722629) q[1];
rz(0.99156143) q[3];
sx q[3];
rz(-2.09822) q[3];
sx q[3];
rz(2.9468731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.110454) q[2];
sx q[2];
rz(-1.217507) q[2];
sx q[2];
rz(2.8793867) q[2];
rz(2.0860705) q[3];
sx q[3];
rz(-1.2473829) q[3];
sx q[3];
rz(-3.054255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6156886) q[0];
sx q[0];
rz(-2.6380802) q[0];
sx q[0];
rz(-2.9600034) q[0];
rz(-1.4008105) q[1];
sx q[1];
rz(-1.1927651) q[1];
sx q[1];
rz(0.84091944) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5530295) q[0];
sx q[0];
rz(-2.8210777) q[0];
sx q[0];
rz(-3.0856087) q[0];
rz(-pi) q[1];
rz(-2.6662331) q[2];
sx q[2];
rz(-0.98612758) q[2];
sx q[2];
rz(0.66051403) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4136728) q[1];
sx q[1];
rz(-2.9693954) q[1];
sx q[1];
rz(-1.3581469) q[1];
rz(1.0686458) q[3];
sx q[3];
rz(-2.4663894) q[3];
sx q[3];
rz(1.5218376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.38222868) q[2];
sx q[2];
rz(-2.9266734) q[2];
sx q[2];
rz(2.0767029) q[2];
rz(3.0801638) q[3];
sx q[3];
rz(-1.3353142) q[3];
sx q[3];
rz(1.5016865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29838022) q[0];
sx q[0];
rz(-1.1792264) q[0];
sx q[0];
rz(0.19905736) q[0];
rz(2.3152323) q[1];
sx q[1];
rz(-2.2233456) q[1];
sx q[1];
rz(0.77702776) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26829241) q[0];
sx q[0];
rz(-1.7275294) q[0];
sx q[0];
rz(0.066989338) q[0];
x q[1];
rz(1.0269182) q[2];
sx q[2];
rz(-1.9444398) q[2];
sx q[2];
rz(-2.0622562) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4118581) q[1];
sx q[1];
rz(-0.61690205) q[1];
sx q[1];
rz(2.0286948) q[1];
rz(-pi) q[2];
x q[2];
rz(0.2774419) q[3];
sx q[3];
rz(-2.146477) q[3];
sx q[3];
rz(-0.68439181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1021154) q[2];
sx q[2];
rz(-1.9934883) q[2];
sx q[2];
rz(-2.6285505) q[2];
rz(-2.3593864) q[3];
sx q[3];
rz(-1.2057883) q[3];
sx q[3];
rz(1.7647083) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2143329) q[0];
sx q[0];
rz(-1.5434649) q[0];
sx q[0];
rz(1.4904892) q[0];
rz(-0.53228846) q[1];
sx q[1];
rz(-0.63096255) q[1];
sx q[1];
rz(-1.5375563) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.45134) q[0];
sx q[0];
rz(-1.7336539) q[0];
sx q[0];
rz(1.4630408) q[0];
rz(-2.1170148) q[2];
sx q[2];
rz(-1.4909296) q[2];
sx q[2];
rz(-0.090674222) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.73276943) q[1];
sx q[1];
rz(-0.63946001) q[1];
sx q[1];
rz(1.6896405) q[1];
x q[2];
rz(-1.8054078) q[3];
sx q[3];
rz(-1.0554093) q[3];
sx q[3];
rz(-2.0976515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.535546) q[2];
sx q[2];
rz(-2.594785) q[2];
sx q[2];
rz(0.01288506) q[2];
rz(1.2031201) q[3];
sx q[3];
rz(-1.4667526) q[3];
sx q[3];
rz(2.119868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1383698) q[0];
sx q[0];
rz(-0.62507451) q[0];
sx q[0];
rz(-1.7417997) q[0];
rz(-2.7895582) q[1];
sx q[1];
rz(-2.0875918) q[1];
sx q[1];
rz(-2.3037516) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96124803) q[0];
sx q[0];
rz(-2.367428) q[0];
sx q[0];
rz(-1.3526239) q[0];
rz(-pi) q[1];
rz(2.3532392) q[2];
sx q[2];
rz(-0.78561312) q[2];
sx q[2];
rz(-2.2503302) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.20866933) q[1];
sx q[1];
rz(-2.1341548) q[1];
sx q[1];
rz(-1.7209919) q[1];
rz(-pi) q[2];
rz(-2.9772163) q[3];
sx q[3];
rz(-1.419074) q[3];
sx q[3];
rz(0.80899948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0248523) q[2];
sx q[2];
rz(-0.96472538) q[2];
sx q[2];
rz(-1.9579197) q[2];
rz(-0.27070326) q[3];
sx q[3];
rz(-0.94047061) q[3];
sx q[3];
rz(-1.6089599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.026641332) q[0];
sx q[0];
rz(-0.6237492) q[0];
sx q[0];
rz(0.57914105) q[0];
rz(-1.2222414) q[1];
sx q[1];
rz(-1.3074343) q[1];
sx q[1];
rz(-2.0319669) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0298808) q[0];
sx q[0];
rz(-1.7222026) q[0];
sx q[0];
rz(-0.16686186) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8366361) q[2];
sx q[2];
rz(-1.6246319) q[2];
sx q[2];
rz(0.75537813) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5286104) q[1];
sx q[1];
rz(-0.41944606) q[1];
sx q[1];
rz(-0.64639133) q[1];
rz(-2.3253981) q[3];
sx q[3];
rz(-2.4494054) q[3];
sx q[3];
rz(-0.029267197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1960725) q[2];
sx q[2];
rz(-1.2976982) q[2];
sx q[2];
rz(0.35745364) q[2];
rz(1.1899905) q[3];
sx q[3];
rz(-1.9414732) q[3];
sx q[3];
rz(1.3397217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59915197) q[0];
sx q[0];
rz(-2.1245133) q[0];
sx q[0];
rz(2.3495112) q[0];
rz(0.7041086) q[1];
sx q[1];
rz(-2.6857565) q[1];
sx q[1];
rz(-1.5584996) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7575551) q[0];
sx q[0];
rz(-2.677394) q[0];
sx q[0];
rz(-3.1350101) q[0];
rz(-2.3249945) q[2];
sx q[2];
rz(-0.77864051) q[2];
sx q[2];
rz(1.1731847) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6686791) q[1];
sx q[1];
rz(-1.1448507) q[1];
sx q[1];
rz(-1.6086474) q[1];
rz(-2.3901863) q[3];
sx q[3];
rz(-1.4045241) q[3];
sx q[3];
rz(2.5010482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.030297) q[2];
sx q[2];
rz(-0.78812495) q[2];
sx q[2];
rz(-0.0041858717) q[2];
rz(-2.8748416) q[3];
sx q[3];
rz(-1.6175852) q[3];
sx q[3];
rz(-2.4193616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5608212) q[0];
sx q[0];
rz(-0.72951356) q[0];
sx q[0];
rz(-0.15604493) q[0];
rz(1.5849628) q[1];
sx q[1];
rz(-0.86235756) q[1];
sx q[1];
rz(2.5468266) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7177009) q[0];
sx q[0];
rz(-2.0621952) q[0];
sx q[0];
rz(1.8609079) q[0];
rz(-pi) q[1];
rz(-1.8608708) q[2];
sx q[2];
rz(-0.89934394) q[2];
sx q[2];
rz(-2.4460276) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.22560355) q[1];
sx q[1];
rz(-2.2609613) q[1];
sx q[1];
rz(0.76721104) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.90601633) q[3];
sx q[3];
rz(-2.182534) q[3];
sx q[3];
rz(-0.25463984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5742089) q[2];
sx q[2];
rz(-1.5625861) q[2];
sx q[2];
rz(1.0279921) q[2];
rz(1.3422286) q[3];
sx q[3];
rz(-2.4864311) q[3];
sx q[3];
rz(1.3348234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8780355) q[0];
sx q[0];
rz(-1.4565383) q[0];
sx q[0];
rz(-0.79826075) q[0];
rz(-0.80052605) q[1];
sx q[1];
rz(-0.48897484) q[1];
sx q[1];
rz(-2.5953603) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6123264) q[0];
sx q[0];
rz(-0.61175013) q[0];
sx q[0];
rz(1.3806369) q[0];
rz(1.9990218) q[2];
sx q[2];
rz(-2.1167123) q[2];
sx q[2];
rz(0.8949309) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.47910467) q[1];
sx q[1];
rz(-0.71245414) q[1];
sx q[1];
rz(2.3478048) q[1];
rz(-pi) q[2];
rz(-1.4978374) q[3];
sx q[3];
rz(-1.3449418) q[3];
sx q[3];
rz(0.42499229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3937248) q[2];
sx q[2];
rz(-1.0900898) q[2];
sx q[2];
rz(0.063035034) q[2];
rz(-0.6978327) q[3];
sx q[3];
rz(-0.2581667) q[3];
sx q[3];
rz(0.91309083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6460655) q[0];
sx q[0];
rz(-2.263948) q[0];
sx q[0];
rz(-2.7255507) q[0];
rz(1.7795732) q[1];
sx q[1];
rz(-1.1241309) q[1];
sx q[1];
rz(-1.2927885) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1406547) q[0];
sx q[0];
rz(-1.7998371) q[0];
sx q[0];
rz(2.8662837) q[0];
rz(-pi) q[1];
rz(2.6854158) q[2];
sx q[2];
rz(-1.8229228) q[2];
sx q[2];
rz(1.9664362) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.19504612) q[1];
sx q[1];
rz(-2.1113696) q[1];
sx q[1];
rz(-0.54624301) q[1];
x q[2];
rz(-1.9927916) q[3];
sx q[3];
rz(-1.8099603) q[3];
sx q[3];
rz(1.07475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.68299874) q[2];
sx q[2];
rz(-1.8051882) q[2];
sx q[2];
rz(-0.90751737) q[2];
rz(1.2609743) q[3];
sx q[3];
rz(-2.5670299) q[3];
sx q[3];
rz(-1.6921836) q[3];
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
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64611971) q[0];
sx q[0];
rz(-1.4541805) q[0];
sx q[0];
rz(-2.4021586) q[0];
rz(-0.44838913) q[1];
sx q[1];
rz(-1.3403475) q[1];
sx q[1];
rz(1.1833804) q[1];
rz(-1.2585959) q[2];
sx q[2];
rz(-1.9070503) q[2];
sx q[2];
rz(-2.287434) q[2];
rz(-0.84239324) q[3];
sx q[3];
rz(-1.471023) q[3];
sx q[3];
rz(-1.4820549) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
