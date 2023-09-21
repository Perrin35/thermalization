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
rz(-3.1080973) q[0];
sx q[0];
rz(-1.7749696) q[0];
rz(-1.051149) q[1];
sx q[1];
rz(-1.4895952) q[1];
sx q[1];
rz(1.1319914) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7819408) q[0];
sx q[0];
rz(-1.9160761) q[0];
sx q[0];
rz(0.79202534) q[0];
rz(-pi) q[1];
rz(1.7366473) q[2];
sx q[2];
rz(-1.7300786) q[2];
sx q[2];
rz(-0.35370358) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7992026) q[1];
sx q[1];
rz(-2.4203165) q[1];
sx q[1];
rz(-0.93625416) q[1];
rz(-pi) q[2];
rz(-0.55478401) q[3];
sx q[3];
rz(-1.3265508) q[3];
sx q[3];
rz(0.82575219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.91036096) q[2];
sx q[2];
rz(-1.8138764) q[2];
sx q[2];
rz(1.1532016) q[2];
rz(2.6575346) q[3];
sx q[3];
rz(-0.81367937) q[3];
sx q[3];
rz(0.4593862) q[3];
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
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83950481) q[0];
sx q[0];
rz(-0.33058259) q[0];
sx q[0];
rz(2.6385345) q[0];
rz(-1.5867651) q[1];
sx q[1];
rz(-0.70650548) q[1];
sx q[1];
rz(-0.15393004) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4707697) q[0];
sx q[0];
rz(-1.5851067) q[0];
sx q[0];
rz(2.0736681) q[0];
rz(2.8106861) q[2];
sx q[2];
rz(-0.94793301) q[2];
sx q[2];
rz(-1.9493584) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1679722) q[1];
sx q[1];
rz(-0.64860839) q[1];
sx q[1];
rz(-0.82778511) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1400938) q[3];
sx q[3];
rz(-0.50977364) q[3];
sx q[3];
rz(1.203376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8872035) q[2];
sx q[2];
rz(-2.7837191) q[2];
sx q[2];
rz(-1.3228234) q[2];
rz(1.4860738) q[3];
sx q[3];
rz(-1.6008987) q[3];
sx q[3];
rz(-2.4310908) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1948497) q[0];
sx q[0];
rz(-1.1179593) q[0];
sx q[0];
rz(0.90993607) q[0];
rz(-0.77727708) q[1];
sx q[1];
rz(-0.83559075) q[1];
sx q[1];
rz(-2.1562703) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5062149) q[0];
sx q[0];
rz(-0.92184508) q[0];
sx q[0];
rz(-0.20091591) q[0];
rz(-pi) q[1];
rz(0.91395949) q[2];
sx q[2];
rz(-1.8310391) q[2];
sx q[2];
rz(2.7514806) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9039771) q[1];
sx q[1];
rz(-1.8198038) q[1];
sx q[1];
rz(1.6275089) q[1];
rz(-2.2306973) q[3];
sx q[3];
rz(-2.1351372) q[3];
sx q[3];
rz(2.6499555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.862792) q[2];
sx q[2];
rz(-1.1547487) q[2];
sx q[2];
rz(1.8939691) q[2];
rz(2.6990081) q[3];
sx q[3];
rz(-1.3675888) q[3];
sx q[3];
rz(0.32143337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2414395) q[0];
sx q[0];
rz(-2.1223919) q[0];
sx q[0];
rz(-1.1234269) q[0];
rz(0.57888794) q[1];
sx q[1];
rz(-1.6826948) q[1];
sx q[1];
rz(-1.7480063) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9959895) q[0];
sx q[0];
rz(-1.9847426) q[0];
sx q[0];
rz(-2.4664509) q[0];
rz(-pi) q[1];
rz(2.7921177) q[2];
sx q[2];
rz(-2.6824143) q[2];
sx q[2];
rz(-1.0897204) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.77809282) q[1];
sx q[1];
rz(-0.85077319) q[1];
sx q[1];
rz(0.6196482) q[1];
rz(0.83538576) q[3];
sx q[3];
rz(-1.7454034) q[3];
sx q[3];
rz(-1.1268827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0397772) q[2];
sx q[2];
rz(-1.5559876) q[2];
sx q[2];
rz(3.139479) q[2];
rz(-2.5801616) q[3];
sx q[3];
rz(-1.0174454) q[3];
sx q[3];
rz(-2.5533365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11557065) q[0];
sx q[0];
rz(-2.0176812) q[0];
sx q[0];
rz(-1.9157238) q[0];
rz(1.4670124) q[1];
sx q[1];
rz(-1.8672698) q[1];
sx q[1];
rz(1.3668758) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0923094) q[0];
sx q[0];
rz(-1.7607848) q[0];
sx q[0];
rz(2.7490194) q[0];
x q[1];
rz(2.3943704) q[2];
sx q[2];
rz(-1.5929654) q[2];
sx q[2];
rz(-1.5378466) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.83134507) q[1];
sx q[1];
rz(-0.86569769) q[1];
sx q[1];
rz(-0.10465937) q[1];
rz(-pi) q[2];
rz(2.6392691) q[3];
sx q[3];
rz(-1.886743) q[3];
sx q[3];
rz(1.5497006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.86429578) q[2];
sx q[2];
rz(-2.0531451) q[2];
sx q[2];
rz(2.7362291) q[2];
rz(2.6799485) q[3];
sx q[3];
rz(-0.82563892) q[3];
sx q[3];
rz(-1.5464787) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49155238) q[0];
sx q[0];
rz(-1.7164282) q[0];
sx q[0];
rz(-0.50338411) q[0];
rz(0.21884306) q[1];
sx q[1];
rz(-1.2501406) q[1];
sx q[1];
rz(0.65178451) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0359905) q[0];
sx q[0];
rz(-2.833617) q[0];
sx q[0];
rz(2.4686747) q[0];
rz(2.3197744) q[2];
sx q[2];
rz(-1.4966674) q[2];
sx q[2];
rz(-0.2766343) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5946878) q[1];
sx q[1];
rz(-1.5190131) q[1];
sx q[1];
rz(-1.7054218) q[1];
rz(-pi) q[2];
x q[2];
rz(0.28292803) q[3];
sx q[3];
rz(-1.471721) q[3];
sx q[3];
rz(-2.2968452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6289604) q[2];
sx q[2];
rz(-1.3971389) q[2];
sx q[2];
rz(2.7499278) q[2];
rz(-0.20646778) q[3];
sx q[3];
rz(-0.73409096) q[3];
sx q[3];
rz(-0.68968836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7713292) q[0];
sx q[0];
rz(-1.6917133) q[0];
sx q[0];
rz(2.1719334) q[0];
rz(0.58352739) q[1];
sx q[1];
rz(-2.0136166) q[1];
sx q[1];
rz(-3.0113509) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9758494) q[0];
sx q[0];
rz(-1.2184869) q[0];
sx q[0];
rz(-1.5177112) q[0];
rz(-pi) q[1];
rz(-0.70003216) q[2];
sx q[2];
rz(-0.80169741) q[2];
sx q[2];
rz(-1.1449555) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.817037) q[1];
sx q[1];
rz(-0.55700028) q[1];
sx q[1];
rz(-0.42028285) q[1];
rz(-pi) q[2];
rz(-1.0153766) q[3];
sx q[3];
rz(-0.65952276) q[3];
sx q[3];
rz(-1.5606172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.51817259) q[2];
sx q[2];
rz(-1.5815846) q[2];
sx q[2];
rz(2.0557892) q[2];
rz(3.0893677) q[3];
sx q[3];
rz(-1.6970535) q[3];
sx q[3];
rz(-2.0558555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(2.4645585) q[0];
sx q[0];
rz(-2.1573986) q[0];
sx q[0];
rz(1.1664671) q[0];
rz(0.72558609) q[1];
sx q[1];
rz(-1.2936932) q[1];
sx q[1];
rz(-2.3988147) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9184473) q[0];
sx q[0];
rz(-2.3728275) q[0];
sx q[0];
rz(0.2785152) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.62959813) q[2];
sx q[2];
rz(-0.78231914) q[2];
sx q[2];
rz(2.9462189) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.10272664) q[1];
sx q[1];
rz(-2.3042149) q[1];
sx q[1];
rz(1.0086165) q[1];
rz(2.5524213) q[3];
sx q[3];
rz(-0.99654752) q[3];
sx q[3];
rz(-1.1816927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.62721884) q[2];
sx q[2];
rz(-1.5292995) q[2];
sx q[2];
rz(0.26088866) q[2];
rz(2.0339113) q[3];
sx q[3];
rz(-1.8543782) q[3];
sx q[3];
rz(2.0598944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
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
rz(-3.127457) q[1];
sx q[1];
rz(-1.8363258) q[1];
sx q[1];
rz(0.65151185) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.053035887) q[0];
sx q[0];
rz(-3.1175107) q[0];
sx q[0];
rz(1.7504897) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3867943) q[2];
sx q[2];
rz(-0.81364606) q[2];
sx q[2];
rz(-2.1458643) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7266255) q[1];
sx q[1];
rz(-1.8888998) q[1];
sx q[1];
rz(3.0148562) q[1];
rz(2.4581962) q[3];
sx q[3];
rz(-0.545524) q[3];
sx q[3];
rz(-0.47513902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8462048) q[2];
sx q[2];
rz(-2.442895) q[2];
sx q[2];
rz(-2.6182168) q[2];
rz(0.30803099) q[3];
sx q[3];
rz(-0.57294661) q[3];
sx q[3];
rz(1.8035005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(-1.4073407) q[0];
sx q[0];
rz(-2.9946406) q[0];
sx q[0];
rz(1.7379606) q[0];
rz(2.667528) q[1];
sx q[1];
rz(-1.6876551) q[1];
sx q[1];
rz(1.3778936) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4353367) q[0];
sx q[0];
rz(-1.0284817) q[0];
sx q[0];
rz(-0.52973912) q[0];
rz(0.84656287) q[2];
sx q[2];
rz(-2.2228974) q[2];
sx q[2];
rz(2.6400499) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.024152012) q[1];
sx q[1];
rz(-1.1785058) q[1];
sx q[1];
rz(1.8333608) q[1];
rz(-1.4078478) q[3];
sx q[3];
rz(-2.1221707) q[3];
sx q[3];
rz(-0.58386246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3906117) q[2];
sx q[2];
rz(-1.6335952) q[2];
sx q[2];
rz(1.6101458) q[2];
rz(-1.6837998) q[3];
sx q[3];
rz(-2.0861574) q[3];
sx q[3];
rz(-2.2916601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67125852) q[0];
sx q[0];
rz(-0.27161921) q[0];
sx q[0];
rz(-2.0441396) q[0];
rz(-0.63256565) q[1];
sx q[1];
rz(-2.0805151) q[1];
sx q[1];
rz(0.17593304) q[1];
rz(-1.8995646) q[2];
sx q[2];
rz(-0.9267926) q[2];
sx q[2];
rz(-0.37315858) q[2];
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
