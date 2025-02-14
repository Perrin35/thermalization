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
rz(0.49864545) q[0];
sx q[0];
rz(-2.5957624) q[0];
sx q[0];
rz(-0.11830615) q[0];
rz(0.91953295) q[1];
sx q[1];
rz(4.4448648) q[1];
sx q[1];
rz(8.2104609) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81033731) q[0];
sx q[0];
rz(-0.48733586) q[0];
sx q[0];
rz(-2.0361855) q[0];
rz(-2.1636398) q[2];
sx q[2];
rz(-1.0314157) q[2];
sx q[2];
rz(-0.025503615) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0561021) q[1];
sx q[1];
rz(-0.63492763) q[1];
sx q[1];
rz(-1.8024428) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6475347) q[3];
sx q[3];
rz(-2.8676448) q[3];
sx q[3];
rz(-0.63989598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.97293568) q[2];
sx q[2];
rz(-1.1309705) q[2];
sx q[2];
rz(-1.0966148) q[2];
rz(-2.9003669) q[3];
sx q[3];
rz(-1.3696407) q[3];
sx q[3];
rz(-0.84478861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.089652561) q[0];
sx q[0];
rz(-1.0900499) q[0];
sx q[0];
rz(0.3845149) q[0];
rz(0.3081201) q[1];
sx q[1];
rz(-0.48079753) q[1];
sx q[1];
rz(-0.62517977) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70474941) q[0];
sx q[0];
rz(-1.6024717) q[0];
sx q[0];
rz(-0.28625536) q[0];
rz(0.7169508) q[2];
sx q[2];
rz(-1.6670818) q[2];
sx q[2];
rz(1.7801628) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.75698419) q[1];
sx q[1];
rz(-0.9238657) q[1];
sx q[1];
rz(-1.4642843) q[1];
rz(-2.6513807) q[3];
sx q[3];
rz(-1.0102838) q[3];
sx q[3];
rz(1.1309159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.008931) q[2];
sx q[2];
rz(-0.93175685) q[2];
sx q[2];
rz(0.30161944) q[2];
rz(-2.610142) q[3];
sx q[3];
rz(-0.88574946) q[3];
sx q[3];
rz(-2.0187995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6279491) q[0];
sx q[0];
rz(-2.1301837) q[0];
sx q[0];
rz(-1.9568141) q[0];
rz(2.9966677) q[1];
sx q[1];
rz(-1.8791608) q[1];
sx q[1];
rz(-1.5941934) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7892706) q[0];
sx q[0];
rz(-2.0579859) q[0];
sx q[0];
rz(2.4643023) q[0];
rz(-pi) q[1];
rz(1.1357665) q[2];
sx q[2];
rz(-0.41450459) q[2];
sx q[2];
rz(-0.5629102) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6559534) q[1];
sx q[1];
rz(-1.9536634) q[1];
sx q[1];
rz(-0.98321557) q[1];
x q[2];
rz(-0.66158847) q[3];
sx q[3];
rz(-0.49822712) q[3];
sx q[3];
rz(-0.029904043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.80921119) q[2];
sx q[2];
rz(-0.96497649) q[2];
sx q[2];
rz(1.3261718) q[2];
rz(2.2281846) q[3];
sx q[3];
rz(-1.923442) q[3];
sx q[3];
rz(-2.6083045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42295414) q[0];
sx q[0];
rz(-2.664743) q[0];
sx q[0];
rz(1.780321) q[0];
rz(2.4286229) q[1];
sx q[1];
rz(-2.0013335) q[1];
sx q[1];
rz(-2.4580809) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8382156) q[0];
sx q[0];
rz(-1.5607173) q[0];
sx q[0];
rz(-2.1440708) q[0];
rz(0.2095378) q[2];
sx q[2];
rz(-1.7895797) q[2];
sx q[2];
rz(-2.0444586) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6660712) q[1];
sx q[1];
rz(-0.88705766) q[1];
sx q[1];
rz(2.4525675) q[1];
rz(1.470742) q[3];
sx q[3];
rz(-1.3276427) q[3];
sx q[3];
rz(0.2167165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6354562) q[2];
sx q[2];
rz(-0.41716245) q[2];
sx q[2];
rz(0.48640856) q[2];
rz(0.5168612) q[3];
sx q[3];
rz(-1.1812527) q[3];
sx q[3];
rz(-1.0335056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0559693) q[0];
sx q[0];
rz(-1.509868) q[0];
sx q[0];
rz(-1.3729209) q[0];
rz(-1.7347451) q[1];
sx q[1];
rz(-0.56766784) q[1];
sx q[1];
rz(2.6645606) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6876887) q[0];
sx q[0];
rz(-2.6756508) q[0];
sx q[0];
rz(2.0201319) q[0];
rz(-2.1978756) q[2];
sx q[2];
rz(-2.5260128) q[2];
sx q[2];
rz(2.1114388) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.43852311) q[1];
sx q[1];
rz(-1.0756936) q[1];
sx q[1];
rz(-2.7994521) q[1];
rz(-pi) q[2];
rz(-2.414235) q[3];
sx q[3];
rz(-2.2148956) q[3];
sx q[3];
rz(-1.3535045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.10739747) q[2];
sx q[2];
rz(-0.98661462) q[2];
sx q[2];
rz(2.6743215) q[2];
rz(0.16921903) q[3];
sx q[3];
rz(-1.7305814) q[3];
sx q[3];
rz(-2.8362078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53448236) q[0];
sx q[0];
rz(-0.87962532) q[0];
sx q[0];
rz(-0.18950263) q[0];
rz(0.27319187) q[1];
sx q[1];
rz(-1.0739645) q[1];
sx q[1];
rz(-0.80064076) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9439745) q[0];
sx q[0];
rz(-1.6139784) q[0];
sx q[0];
rz(1.3373242) q[0];
rz(-pi) q[1];
rz(0.24836274) q[2];
sx q[2];
rz(-1.6585322) q[2];
sx q[2];
rz(1.9875634) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8756518) q[1];
sx q[1];
rz(-1.3752626) q[1];
sx q[1];
rz(0.13130782) q[1];
x q[2];
rz(0.74108776) q[3];
sx q[3];
rz(-2.4588278) q[3];
sx q[3];
rz(0.53214754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2022986) q[2];
sx q[2];
rz(-2.9677128) q[2];
sx q[2];
rz(2.6046216) q[2];
rz(-1.8592853) q[3];
sx q[3];
rz(-1.8351646) q[3];
sx q[3];
rz(1.6622701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.031161664) q[0];
sx q[0];
rz(-2.3261676) q[0];
sx q[0];
rz(-1.331331) q[0];
rz(-1.4415461) q[1];
sx q[1];
rz(-1.6621637) q[1];
sx q[1];
rz(1.2225245) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.545802) q[0];
sx q[0];
rz(-1.8078601) q[0];
sx q[0];
rz(-1.14181) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4961081) q[2];
sx q[2];
rz(-2.1213581) q[2];
sx q[2];
rz(-2.2979743) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4536087) q[1];
sx q[1];
rz(-0.88064945) q[1];
sx q[1];
rz(0.28120561) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8172713) q[3];
sx q[3];
rz(-2.685084) q[3];
sx q[3];
rz(-2.498561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9407201) q[2];
sx q[2];
rz(-2.6437289) q[2];
sx q[2];
rz(0.3328003) q[2];
rz(-2.8424272) q[3];
sx q[3];
rz(-1.9009512) q[3];
sx q[3];
rz(-0.021473264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0126295) q[0];
sx q[0];
rz(-0.7681995) q[0];
sx q[0];
rz(-0.27498883) q[0];
rz(0.081427447) q[1];
sx q[1];
rz(-0.80137253) q[1];
sx q[1];
rz(1.3224695) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7985991) q[0];
sx q[0];
rz(-1.3396153) q[0];
sx q[0];
rz(-0.061132758) q[0];
x q[1];
rz(2.7733388) q[2];
sx q[2];
rz(-0.91742951) q[2];
sx q[2];
rz(-0.39619941) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7863955) q[1];
sx q[1];
rz(-1.2010368) q[1];
sx q[1];
rz(-2.3936097) q[1];
rz(-2.2127719) q[3];
sx q[3];
rz(-0.95212338) q[3];
sx q[3];
rz(-0.88255461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.72121173) q[2];
sx q[2];
rz(-0.12004852) q[2];
sx q[2];
rz(-1.7397286) q[2];
rz(1.8574572) q[3];
sx q[3];
rz(-1.1893585) q[3];
sx q[3];
rz(-3.106626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5597124) q[0];
sx q[0];
rz(-1.548883) q[0];
sx q[0];
rz(-0.48932073) q[0];
rz(-3.1210461) q[1];
sx q[1];
rz(-1.2362044) q[1];
sx q[1];
rz(2.4403341) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0017635) q[0];
sx q[0];
rz(-1.0047303) q[0];
sx q[0];
rz(-0.28621121) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6270774) q[2];
sx q[2];
rz(-2.4372134) q[2];
sx q[2];
rz(-3.0075108) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.32966954) q[1];
sx q[1];
rz(-0.81336248) q[1];
sx q[1];
rz(-1.2438891) q[1];
x q[2];
rz(2.5908628) q[3];
sx q[3];
rz(-2.3756873) q[3];
sx q[3];
rz(2.285241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7265085) q[2];
sx q[2];
rz(-2.5231611) q[2];
sx q[2];
rz(-2.6611967) q[2];
rz(-1.7324309) q[3];
sx q[3];
rz(-1.3536072) q[3];
sx q[3];
rz(2.8044146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3163863) q[0];
sx q[0];
rz(-1.2403064) q[0];
sx q[0];
rz(0.41148841) q[0];
rz(-1.2443789) q[1];
sx q[1];
rz(-1.139541) q[1];
sx q[1];
rz(-2.8172475) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1169057) q[0];
sx q[0];
rz(-1.1984662) q[0];
sx q[0];
rz(-0.056554746) q[0];
rz(-pi) q[1];
rz(-0.80071189) q[2];
sx q[2];
rz(-2.6691648) q[2];
sx q[2];
rz(-2.1683482) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.577475) q[1];
sx q[1];
rz(-1.0077969) q[1];
sx q[1];
rz(-3.0433118) q[1];
rz(-pi) q[2];
rz(1.5839229) q[3];
sx q[3];
rz(-1.6919976) q[3];
sx q[3];
rz(-1.2016202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1856498) q[2];
sx q[2];
rz(-2.0273384) q[2];
sx q[2];
rz(0.11697098) q[2];
rz(-0.95514917) q[3];
sx q[3];
rz(-2.9626466) q[3];
sx q[3];
rz(-0.54168934) q[3];
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
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2627926) q[0];
sx q[0];
rz(-0.87651064) q[0];
sx q[0];
rz(0.87690092) q[0];
rz(-2.0211438) q[1];
sx q[1];
rz(-0.81605492) q[1];
sx q[1];
rz(1.1647404) q[1];
rz(-2.929959) q[2];
sx q[2];
rz(-1.0734049) q[2];
sx q[2];
rz(2.887101) q[2];
rz(1.6616115) q[3];
sx q[3];
rz(-0.27169966) q[3];
sx q[3];
rz(-3.0063831) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
