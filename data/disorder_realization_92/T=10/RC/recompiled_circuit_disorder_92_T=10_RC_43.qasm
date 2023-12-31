OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.6264412) q[0];
sx q[0];
rz(3.175088) q[0];
sx q[0];
rz(7.6498084) q[0];
rz(-1.051149) q[1];
sx q[1];
rz(-1.4895952) q[1];
sx q[1];
rz(1.1319914) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53349797) q[0];
sx q[0];
rz(-2.29288) q[0];
sx q[0];
rz(0.46790926) q[0];
rz(-pi) q[1];
rz(1.7366473) q[2];
sx q[2];
rz(-1.4115141) q[2];
sx q[2];
rz(0.35370358) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1177897) q[1];
sx q[1];
rz(-2.1315247) q[1];
sx q[1];
rz(0.48052131) q[1];
rz(-pi) q[2];
rz(-2.5868086) q[3];
sx q[3];
rz(-1.3265508) q[3];
sx q[3];
rz(2.3158405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.91036096) q[2];
sx q[2];
rz(-1.3277162) q[2];
sx q[2];
rz(-1.1532016) q[2];
rz(0.48405805) q[3];
sx q[3];
rz(-2.3279133) q[3];
sx q[3];
rz(-2.6822065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83950481) q[0];
sx q[0];
rz(-2.8110101) q[0];
sx q[0];
rz(0.50305811) q[0];
rz(1.5867651) q[1];
sx q[1];
rz(-0.70650548) q[1];
sx q[1];
rz(0.15393004) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.126037) q[0];
sx q[0];
rz(-0.50305788) q[0];
sx q[0];
rz(-1.6004827) q[0];
rz(-pi) q[1];
rz(-2.2203127) q[2];
sx q[2];
rz(-1.8378471) q[2];
sx q[2];
rz(2.9608179) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.82984) q[1];
sx q[1];
rz(-1.1097739) q[1];
sx q[1];
rz(-2.6677569) q[1];
rz(-pi) q[2];
x q[2];
rz(0.0014988621) q[3];
sx q[3];
rz(-2.631819) q[3];
sx q[3];
rz(-1.9382167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8872035) q[2];
sx q[2];
rz(-0.35787359) q[2];
sx q[2];
rz(1.3228234) q[2];
rz(-1.6555188) q[3];
sx q[3];
rz(-1.6008987) q[3];
sx q[3];
rz(0.71050182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94674295) q[0];
sx q[0];
rz(-2.0236334) q[0];
sx q[0];
rz(-2.2316566) q[0];
rz(2.3643156) q[1];
sx q[1];
rz(-2.3060019) q[1];
sx q[1];
rz(2.1562703) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6353778) q[0];
sx q[0];
rz(-0.92184508) q[0];
sx q[0];
rz(-2.9406767) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8172242) q[2];
sx q[2];
rz(-2.2019221) q[2];
sx q[2];
rz(2.1567675) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9039771) q[1];
sx q[1];
rz(-1.3217889) q[1];
sx q[1];
rz(-1.5140837) q[1];
rz(2.4661029) q[3];
sx q[3];
rz(-2.1152861) q[3];
sx q[3];
rz(2.4558223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2788006) q[2];
sx q[2];
rz(-1.9868439) q[2];
sx q[2];
rz(-1.8939691) q[2];
rz(0.4425846) q[3];
sx q[3];
rz(-1.3675888) q[3];
sx q[3];
rz(2.8201593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2414395) q[0];
sx q[0];
rz(-2.1223919) q[0];
sx q[0];
rz(2.0181657) q[0];
rz(-2.5627047) q[1];
sx q[1];
rz(-1.6826948) q[1];
sx q[1];
rz(1.3935864) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9959895) q[0];
sx q[0];
rz(-1.1568501) q[0];
sx q[0];
rz(2.4664509) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4030928) q[2];
sx q[2];
rz(-2.0003013) q[2];
sx q[2];
rz(-1.475856) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.77809282) q[1];
sx q[1];
rz(-2.2908195) q[1];
sx q[1];
rz(-2.5219445) q[1];
rz(-1.3136775) q[3];
sx q[3];
rz(-2.3895279) q[3];
sx q[3];
rz(2.5079692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1018155) q[2];
sx q[2];
rz(-1.5856051) q[2];
sx q[2];
rz(3.139479) q[2];
rz(2.5801616) q[3];
sx q[3];
rz(-1.0174454) q[3];
sx q[3];
rz(2.5533365) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.026022) q[0];
sx q[0];
rz(-1.1239115) q[0];
sx q[0];
rz(-1.2258688) q[0];
rz(1.4670124) q[1];
sx q[1];
rz(-1.8672698) q[1];
sx q[1];
rz(-1.7747169) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5850692) q[0];
sx q[0];
rz(-1.1856623) q[0];
sx q[0];
rz(-1.3655846) q[0];
rz(-1.6010124) q[2];
sx q[2];
rz(-2.317791) q[2];
sx q[2];
rz(3.1291762) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4701159) q[1];
sx q[1];
rz(-1.4911545) q[1];
sx q[1];
rz(-0.86298841) q[1];
x q[2];
rz(2.5451238) q[3];
sx q[3];
rz(-0.58613741) q[3];
sx q[3];
rz(-0.49367192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.86429578) q[2];
sx q[2];
rz(-1.0884476) q[2];
sx q[2];
rz(2.7362291) q[2];
rz(-2.6799485) q[3];
sx q[3];
rz(-0.82563892) q[3];
sx q[3];
rz(-1.595114) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6500403) q[0];
sx q[0];
rz(-1.7164282) q[0];
sx q[0];
rz(2.6382085) q[0];
rz(0.21884306) q[1];
sx q[1];
rz(-1.8914521) q[1];
sx q[1];
rz(2.4898081) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9572243) q[0];
sx q[0];
rz(-1.3807218) q[0];
sx q[0];
rz(-0.24380542) q[0];
rz(-pi) q[1];
rz(-2.3197744) q[2];
sx q[2];
rz(-1.6449252) q[2];
sx q[2];
rz(2.8649583) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.3411322) q[1];
sx q[1];
rz(-0.144185) q[1];
sx q[1];
rz(-1.202281) q[1];
rz(-1.6739453) q[3];
sx q[3];
rz(-1.2892937) q[3];
sx q[3];
rz(2.3867949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.51263222) q[2];
sx q[2];
rz(-1.3971389) q[2];
sx q[2];
rz(0.39166489) q[2];
rz(-0.20646778) q[3];
sx q[3];
rz(-0.73409096) q[3];
sx q[3];
rz(-0.68968836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3702635) q[0];
sx q[0];
rz(-1.6917133) q[0];
sx q[0];
rz(2.1719334) q[0];
rz(-2.5580653) q[1];
sx q[1];
rz(-1.1279761) q[1];
sx q[1];
rz(3.0113509) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9758494) q[0];
sx q[0];
rz(-1.9231057) q[0];
sx q[0];
rz(-1.5177112) q[0];
x q[1];
rz(-2.4415605) q[2];
sx q[2];
rz(-2.3398952) q[2];
sx q[2];
rz(-1.1449555) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.817037) q[1];
sx q[1];
rz(-2.5845924) q[1];
sx q[1];
rz(-2.7213098) q[1];
x q[2];
rz(2.7534915) q[3];
sx q[3];
rz(-1.023205) q[3];
sx q[3];
rz(-2.2263118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.51817259) q[2];
sx q[2];
rz(-1.5815846) q[2];
sx q[2];
rz(-2.0557892) q[2];
rz(3.0893677) q[3];
sx q[3];
rz(-1.6970535) q[3];
sx q[3];
rz(1.0857371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6770342) q[0];
sx q[0];
rz(-0.98419404) q[0];
sx q[0];
rz(1.9751256) q[0];
rz(2.4160066) q[1];
sx q[1];
rz(-1.8478994) q[1];
sx q[1];
rz(0.74277791) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22314534) q[0];
sx q[0];
rz(-0.76876516) q[0];
sx q[0];
rz(0.2785152) q[0];
x q[1];
rz(0.62959813) q[2];
sx q[2];
rz(-0.78231914) q[2];
sx q[2];
rz(0.19537374) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.038866) q[1];
sx q[1];
rz(-2.3042149) q[1];
sx q[1];
rz(2.1329761) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.90950305) q[3];
sx q[3];
rz(-1.0854183) q[3];
sx q[3];
rz(0.73736008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.62721884) q[2];
sx q[2];
rz(-1.6122931) q[2];
sx q[2];
rz(2.880704) q[2];
rz(-1.1076814) q[3];
sx q[3];
rz(-1.8543782) q[3];
sx q[3];
rz(2.0598944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.35659197) q[0];
sx q[0];
rz(-2.3605425) q[0];
sx q[0];
rz(0.93604952) q[0];
rz(0.014135663) q[1];
sx q[1];
rz(-1.8363258) q[1];
sx q[1];
rz(-2.4900808) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0148841) q[0];
sx q[0];
rz(-1.5944905) q[0];
sx q[0];
rz(-0.004304927) q[0];
rz(-pi) q[1];
rz(2.3867943) q[2];
sx q[2];
rz(-0.81364606) q[2];
sx q[2];
rz(-0.99572832) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.028101746) q[1];
sx q[1];
rz(-0.3416225) q[1];
sx q[1];
rz(-1.9373059) q[1];
x q[2];
rz(1.2048079) q[3];
sx q[3];
rz(-1.9848739) q[3];
sx q[3];
rz(0.28596349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8462048) q[2];
sx q[2];
rz(-0.69869763) q[2];
sx q[2];
rz(-0.52337581) q[2];
rz(0.30803099) q[3];
sx q[3];
rz(-2.568646) q[3];
sx q[3];
rz(1.3380922) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.734252) q[0];
sx q[0];
rz(-0.14695209) q[0];
sx q[0];
rz(-1.403632) q[0];
rz(-0.47406468) q[1];
sx q[1];
rz(-1.4539376) q[1];
sx q[1];
rz(1.7636991) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5546075) q[0];
sx q[0];
rz(-2.4025612) q[0];
sx q[0];
rz(-0.87297312) q[0];
rz(-pi) q[1];
rz(-2.2950298) q[2];
sx q[2];
rz(-2.2228974) q[2];
sx q[2];
rz(-0.50154274) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.1174406) q[1];
sx q[1];
rz(-1.1785058) q[1];
sx q[1];
rz(1.3082318) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7337448) q[3];
sx q[3];
rz(-2.1221707) q[3];
sx q[3];
rz(-2.5577302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3906117) q[2];
sx q[2];
rz(-1.6335952) q[2];
sx q[2];
rz(-1.5314468) q[2];
rz(-1.4577929) q[3];
sx q[3];
rz(-2.0861574) q[3];
sx q[3];
rz(2.2916601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67125852) q[0];
sx q[0];
rz(-2.8699734) q[0];
sx q[0];
rz(1.097453) q[0];
rz(-2.509027) q[1];
sx q[1];
rz(-1.0610776) q[1];
sx q[1];
rz(-2.9656596) q[1];
rz(0.67062545) q[2];
sx q[2];
rz(-1.309633) q[2];
sx q[2];
rz(1.39967) q[2];
rz(-0.87631638) q[3];
sx q[3];
rz(-0.96517589) q[3];
sx q[3];
rz(2.518566) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
