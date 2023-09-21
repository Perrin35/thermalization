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
rz(4.7935901) q[1];
sx q[1];
rz(10.556769) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6080947) q[0];
sx q[0];
rz(-2.29288) q[0];
sx q[0];
rz(-2.6736834) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.79905431) q[2];
sx q[2];
rz(-2.91215) q[2];
sx q[2];
rz(0.45861751) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.72343091) q[1];
sx q[1];
rz(-1.9730113) q[1];
sx q[1];
rz(0.9546141) q[1];
rz(-pi) q[2];
x q[2];
rz(0.44191435) q[3];
sx q[3];
rz(-0.60097296) q[3];
sx q[3];
rz(1.1170944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2312317) q[2];
sx q[2];
rz(-1.3277162) q[2];
sx q[2];
rz(-1.1532016) q[2];
rz(-2.6575346) q[3];
sx q[3];
rz(-0.81367937) q[3];
sx q[3];
rz(2.6822065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83950481) q[0];
sx q[0];
rz(-2.8110101) q[0];
sx q[0];
rz(-0.50305811) q[0];
rz(-1.5867651) q[1];
sx q[1];
rz(-2.4350872) q[1];
sx q[1];
rz(0.15393004) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0494369) q[0];
sx q[0];
rz(-1.0679809) q[0];
sx q[0];
rz(-0.016331971) q[0];
x q[1];
rz(-0.33090654) q[2];
sx q[2];
rz(-2.1936596) q[2];
sx q[2];
rz(1.9493584) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.31175266) q[1];
sx q[1];
rz(-1.1097739) q[1];
sx q[1];
rz(0.47383576) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5716343) q[3];
sx q[3];
rz(-2.0805693) q[3];
sx q[3];
rz(1.9364995) q[3];
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
rz(1.6555188) q[3];
sx q[3];
rz(-1.6008987) q[3];
sx q[3];
rz(2.4310908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(-0.94674295) q[0];
sx q[0];
rz(-1.1179593) q[0];
sx q[0];
rz(-0.90993607) q[0];
rz(-0.77727708) q[1];
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
rz(1.9604208) q[0];
sx q[0];
rz(-2.4665678) q[0];
sx q[0];
rz(-1.8280562) q[0];
rz(0.91395949) q[2];
sx q[2];
rz(-1.3105536) q[2];
sx q[2];
rz(-2.7514806) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2376155) q[1];
sx q[1];
rz(-1.3217889) q[1];
sx q[1];
rz(-1.6275089) q[1];
rz(-pi) q[2];
rz(-2.4661029) q[3];
sx q[3];
rz(-1.0263066) q[3];
sx q[3];
rz(-0.68577037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.862792) q[2];
sx q[2];
rz(-1.9868439) q[2];
sx q[2];
rz(1.2476236) q[2];
rz(0.4425846) q[3];
sx q[3];
rz(-1.7740039) q[3];
sx q[3];
rz(0.32143337) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9001532) q[0];
sx q[0];
rz(-2.1223919) q[0];
sx q[0];
rz(-1.1234269) q[0];
rz(-2.5627047) q[1];
sx q[1];
rz(-1.6826948) q[1];
sx q[1];
rz(1.3935864) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2508586) q[0];
sx q[0];
rz(-2.3669741) q[0];
sx q[0];
rz(-2.5289092) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7384999) q[2];
sx q[2];
rz(-1.1412914) q[2];
sx q[2];
rz(1.475856) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.77809282) q[1];
sx q[1];
rz(-2.2908195) q[1];
sx q[1];
rz(-0.6196482) q[1];
x q[2];
rz(0.83538576) q[3];
sx q[3];
rz(-1.3961892) q[3];
sx q[3];
rz(-2.01471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0397772) q[2];
sx q[2];
rz(-1.5856051) q[2];
sx q[2];
rz(-0.0021136443) q[2];
rz(2.5801616) q[3];
sx q[3];
rz(-1.0174454) q[3];
sx q[3];
rz(2.5533365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.026022) q[0];
sx q[0];
rz(-2.0176812) q[0];
sx q[0];
rz(-1.9157238) q[0];
rz(-1.4670124) q[1];
sx q[1];
rz(-1.2743228) q[1];
sx q[1];
rz(1.3668758) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0923094) q[0];
sx q[0];
rz(-1.7607848) q[0];
sx q[0];
rz(-0.39257322) q[0];
rz(-pi) q[1];
x q[1];
rz(0.032614313) q[2];
sx q[2];
rz(-2.3941052) q[2];
sx q[2];
rz(3.0847197) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.67147672) q[1];
sx q[1];
rz(-1.4911545) q[1];
sx q[1];
rz(0.86298841) q[1];
rz(-pi) q[2];
rz(-1.2138052) q[3];
sx q[3];
rz(-2.0461267) q[3];
sx q[3];
rz(2.9514422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2772969) q[2];
sx q[2];
rz(-2.0531451) q[2];
sx q[2];
rz(0.40536353) q[2];
rz(-0.46164414) q[3];
sx q[3];
rz(-0.82563892) q[3];
sx q[3];
rz(-1.5464787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49155238) q[0];
sx q[0];
rz(-1.4251645) q[0];
sx q[0];
rz(-0.50338411) q[0];
rz(0.21884306) q[1];
sx q[1];
rz(-1.2501406) q[1];
sx q[1];
rz(-2.4898081) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0359905) q[0];
sx q[0];
rz(-0.30797568) q[0];
sx q[0];
rz(-2.4686747) q[0];
rz(-pi) q[1];
rz(1.6794372) q[2];
sx q[2];
rz(-2.3896653) q[2];
sx q[2];
rz(1.7679364) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.1106907) q[1];
sx q[1];
rz(-1.4363524) q[1];
sx q[1];
rz(3.0893374) q[1];
rz(0.28292803) q[3];
sx q[3];
rz(-1.471721) q[3];
sx q[3];
rz(-2.2968452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6289604) q[2];
sx q[2];
rz(-1.7444538) q[2];
sx q[2];
rz(0.39166489) q[2];
rz(-0.20646778) q[3];
sx q[3];
rz(-2.4075017) q[3];
sx q[3];
rz(-2.4519043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3702635) q[0];
sx q[0];
rz(-1.6917133) q[0];
sx q[0];
rz(0.96965924) q[0];
rz(0.58352739) q[1];
sx q[1];
rz(-1.1279761) q[1];
sx q[1];
rz(-0.13024174) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1657432) q[0];
sx q[0];
rz(-1.9231057) q[0];
sx q[0];
rz(-1.6238814) q[0];
rz(-0.70003216) q[2];
sx q[2];
rz(-2.3398952) q[2];
sx q[2];
rz(1.1449555) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0252467) q[1];
sx q[1];
rz(-1.7882008) q[1];
sx q[1];
rz(-2.6245963) q[1];
rz(-pi) q[2];
rz(1.0153766) q[3];
sx q[3];
rz(-0.65952276) q[3];
sx q[3];
rz(-1.5809755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.51817259) q[2];
sx q[2];
rz(-1.5815846) q[2];
sx q[2];
rz(-1.0858034) q[2];
rz(-3.0893677) q[3];
sx q[3];
rz(-1.4445392) q[3];
sx q[3];
rz(1.0857371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6770342) q[0];
sx q[0];
rz(-0.98419404) q[0];
sx q[0];
rz(-1.9751256) q[0];
rz(2.4160066) q[1];
sx q[1];
rz(-1.8478994) q[1];
sx q[1];
rz(0.74277791) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1449435) q[0];
sx q[0];
rz(-1.7631233) q[0];
sx q[0];
rz(0.74914519) q[0];
rz(-pi) q[1];
rz(-0.6767512) q[2];
sx q[2];
rz(-1.1427715) q[2];
sx q[2];
rz(-1.2892436) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0689773) q[1];
sx q[1];
rz(-1.9779357) q[1];
sx q[1];
rz(-0.8168656) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5524213) q[3];
sx q[3];
rz(-2.1450451) q[3];
sx q[3];
rz(1.1816927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5143738) q[2];
sx q[2];
rz(-1.5292995) q[2];
sx q[2];
rz(-2.880704) q[2];
rz(1.1076814) q[3];
sx q[3];
rz(-1.2872144) q[3];
sx q[3];
rz(-1.0816983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7850007) q[0];
sx q[0];
rz(-2.3605425) q[0];
sx q[0];
rz(-2.2055431) q[0];
rz(-0.014135663) q[1];
sx q[1];
rz(-1.8363258) q[1];
sx q[1];
rz(-0.65151185) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6974029) q[0];
sx q[0];
rz(-1.5751) q[0];
sx q[0];
rz(1.547102) q[0];
rz(-pi) q[1];
rz(0.65666171) q[2];
sx q[2];
rz(-1.0495532) q[2];
sx q[2];
rz(-1.148664) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9459322) q[1];
sx q[1];
rz(-1.6911427) q[1];
sx q[1];
rz(1.8912998) q[1];
rz(-pi) q[2];
rz(-1.2048079) q[3];
sx q[3];
rz(-1.9848739) q[3];
sx q[3];
rz(-0.28596349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2953879) q[2];
sx q[2];
rz(-2.442895) q[2];
sx q[2];
rz(2.6182168) q[2];
rz(-0.30803099) q[3];
sx q[3];
rz(-0.57294661) q[3];
sx q[3];
rz(-1.8035005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.734252) q[0];
sx q[0];
rz(-2.9946406) q[0];
sx q[0];
rz(-1.7379606) q[0];
rz(2.667528) q[1];
sx q[1];
rz(-1.6876551) q[1];
sx q[1];
rz(-1.7636991) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70625593) q[0];
sx q[0];
rz(-2.113111) q[0];
sx q[0];
rz(-0.52973912) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3465967) q[2];
sx q[2];
rz(-1.0161875) q[2];
sx q[2];
rz(2.5650052) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5046895) q[1];
sx q[1];
rz(-2.6733589) q[1];
sx q[1];
rz(-2.5813028) q[1];
x q[2];
rz(1.4078478) q[3];
sx q[3];
rz(-2.1221707) q[3];
sx q[3];
rz(-2.5577302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3906117) q[2];
sx q[2];
rz(-1.5079974) q[2];
sx q[2];
rz(-1.5314468) q[2];
rz(-1.4577929) q[3];
sx q[3];
rz(-2.0861574) q[3];
sx q[3];
rz(-0.84993258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4703341) q[0];
sx q[0];
rz(-2.8699734) q[0];
sx q[0];
rz(1.097453) q[0];
rz(0.63256565) q[1];
sx q[1];
rz(-1.0610776) q[1];
sx q[1];
rz(-2.9656596) q[1];
rz(1.8995646) q[2];
sx q[2];
rz(-2.2148) q[2];
sx q[2];
rz(2.7684341) q[2];
rz(0.87631638) q[3];
sx q[3];
rz(-2.1764168) q[3];
sx q[3];
rz(-0.62302667) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];