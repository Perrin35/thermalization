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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3596518) q[0];
sx q[0];
rz(-1.2255166) q[0];
sx q[0];
rz(2.3495673) q[0];
x q[1];
rz(-1.7366473) q[2];
sx q[2];
rz(-1.4115141) q[2];
sx q[2];
rz(2.7878891) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7992026) q[1];
sx q[1];
rz(-2.4203165) q[1];
sx q[1];
rz(-2.2053385) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.44191435) q[3];
sx q[3];
rz(-2.5406197) q[3];
sx q[3];
rz(1.1170944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2312317) q[2];
sx q[2];
rz(-1.8138764) q[2];
sx q[2];
rz(1.1532016) q[2];
rz(-2.6575346) q[3];
sx q[3];
rz(-0.81367937) q[3];
sx q[3];
rz(2.6822065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3020878) q[0];
sx q[0];
rz(-2.8110101) q[0];
sx q[0];
rz(-2.6385345) q[0];
rz(-1.5867651) q[1];
sx q[1];
rz(-0.70650548) q[1];
sx q[1];
rz(2.9876626) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.670823) q[0];
sx q[0];
rz(-1.5851067) q[0];
sx q[0];
rz(-2.0736681) q[0];
rz(-pi) q[1];
rz(-0.33090654) q[2];
sx q[2];
rz(-0.94793301) q[2];
sx q[2];
rz(-1.9493584) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.82984) q[1];
sx q[1];
rz(-1.1097739) q[1];
sx q[1];
rz(2.6677569) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1400938) q[3];
sx q[3];
rz(-2.631819) q[3];
sx q[3];
rz(1.9382167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8872035) q[2];
sx q[2];
rz(-0.35787359) q[2];
sx q[2];
rz(-1.3228234) q[2];
rz(1.6555188) q[3];
sx q[3];
rz(-1.5406939) q[3];
sx q[3];
rz(0.71050182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94674295) q[0];
sx q[0];
rz(-2.0236334) q[0];
sx q[0];
rz(2.2316566) q[0];
rz(-0.77727708) q[1];
sx q[1];
rz(-2.3060019) q[1];
sx q[1];
rz(-0.98532239) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5062149) q[0];
sx q[0];
rz(-2.2197476) q[0];
sx q[0];
rz(-2.9406767) q[0];
rz(2.2276332) q[2];
sx q[2];
rz(-1.3105536) q[2];
sx q[2];
rz(2.7514806) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2376155) q[1];
sx q[1];
rz(-1.3217889) q[1];
sx q[1];
rz(-1.5140837) q[1];
x q[2];
rz(-2.4661029) q[3];
sx q[3];
rz(-2.1152861) q[3];
sx q[3];
rz(0.68577037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2788006) q[2];
sx q[2];
rz(-1.1547487) q[2];
sx q[2];
rz(1.2476236) q[2];
rz(2.6990081) q[3];
sx q[3];
rz(-1.3675888) q[3];
sx q[3];
rz(-2.8201593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2414395) q[0];
sx q[0];
rz(-1.0192008) q[0];
sx q[0];
rz(-1.1234269) q[0];
rz(0.57888794) q[1];
sx q[1];
rz(-1.6826948) q[1];
sx q[1];
rz(-1.7480063) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14560315) q[0];
sx q[0];
rz(-1.9847426) q[0];
sx q[0];
rz(0.67514174) q[0];
rz(-pi) q[1];
rz(-1.4030928) q[2];
sx q[2];
rz(-1.1412914) q[2];
sx q[2];
rz(1.475856) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.046603831) q[1];
sx q[1];
rz(-2.2294083) q[1];
sx q[1];
rz(0.98595001) q[1];
x q[2];
rz(1.8279151) q[3];
sx q[3];
rz(-2.3895279) q[3];
sx q[3];
rz(-0.63362345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1018155) q[2];
sx q[2];
rz(-1.5559876) q[2];
sx q[2];
rz(0.0021136443) q[2];
rz(2.5801616) q[3];
sx q[3];
rz(-1.0174454) q[3];
sx q[3];
rz(2.5533365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11557065) q[0];
sx q[0];
rz(-1.1239115) q[0];
sx q[0];
rz(-1.9157238) q[0];
rz(-1.6745802) q[1];
sx q[1];
rz(-1.2743228) q[1];
sx q[1];
rz(1.7747169) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5565235) q[0];
sx q[0];
rz(-1.9559304) q[0];
sx q[0];
rz(1.3655846) q[0];
rz(-pi) q[1];
rz(-1.6010124) q[2];
sx q[2];
rz(-0.82380166) q[2];
sx q[2];
rz(0.012416427) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1495684) q[1];
sx q[1];
rz(-2.4300886) q[1];
sx q[1];
rz(1.4486538) q[1];
rz(1.2138052) q[3];
sx q[3];
rz(-1.095466) q[3];
sx q[3];
rz(2.9514422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-1.2501406) q[1];
sx q[1];
rz(-2.4898081) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1843684) q[0];
sx q[0];
rz(-1.3807218) q[0];
sx q[0];
rz(0.24380542) q[0];
rz(-1.6794372) q[2];
sx q[2];
rz(-2.3896653) q[2];
sx q[2];
rz(1.3736563) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5469049) q[1];
sx q[1];
rz(-1.6225796) q[1];
sx q[1];
rz(1.7054218) q[1];
x q[2];
rz(-0.34206335) q[3];
sx q[3];
rz(-2.8422589) q[3];
sx q[3];
rz(2.743486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.51263222) q[2];
sx q[2];
rz(-1.7444538) q[2];
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
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3702635) q[0];
sx q[0];
rz(-1.4498793) q[0];
sx q[0];
rz(-2.1719334) q[0];
rz(2.5580653) q[1];
sx q[1];
rz(-2.0136166) q[1];
sx q[1];
rz(3.0113509) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42338615) q[0];
sx q[0];
rz(-1.6206181) q[0];
sx q[0];
rz(0.35276619) q[0];
x q[1];
rz(-0.70003216) q[2];
sx q[2];
rz(-0.80169741) q[2];
sx q[2];
rz(1.9966372) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8091734) q[1];
sx q[1];
rz(-1.0671339) q[1];
sx q[1];
rz(-1.3219576) q[1];
rz(-1.0153766) q[3];
sx q[3];
rz(-0.65952276) q[3];
sx q[3];
rz(-1.5606172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6234201) q[2];
sx q[2];
rz(-1.5600081) q[2];
sx q[2];
rz(1.0858034) q[2];
rz(3.0893677) q[3];
sx q[3];
rz(-1.6970535) q[3];
sx q[3];
rz(1.0857371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6770342) q[0];
sx q[0];
rz(-0.98419404) q[0];
sx q[0];
rz(1.9751256) q[0];
rz(-2.4160066) q[1];
sx q[1];
rz(-1.2936932) q[1];
sx q[1];
rz(-2.3988147) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9966492) q[0];
sx q[0];
rz(-1.3784694) q[0];
sx q[0];
rz(0.74914519) q[0];
x q[1];
rz(-2.4648415) q[2];
sx q[2];
rz(-1.1427715) q[2];
sx q[2];
rz(-1.8523491) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0726154) q[1];
sx q[1];
rz(-1.163657) q[1];
sx q[1];
rz(-2.3247271) q[1];
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
rz(-pi/2) q[1];
sx q[1];
rz(0.62721884) q[2];
sx q[2];
rz(-1.6122931) q[2];
sx q[2];
rz(-0.26088866) q[2];
rz(-1.1076814) q[3];
sx q[3];
rz(-1.8543782) q[3];
sx q[3];
rz(2.0598944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35659197) q[0];
sx q[0];
rz(-0.78105015) q[0];
sx q[0];
rz(-0.93604952) q[0];
rz(3.127457) q[1];
sx q[1];
rz(-1.8363258) q[1];
sx q[1];
rz(2.4900808) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0885568) q[0];
sx q[0];
rz(-0.024081973) q[0];
sx q[0];
rz(1.3911029) q[0];
rz(-pi) q[1];
rz(-2.4849309) q[2];
sx q[2];
rz(-2.0920394) q[2];
sx q[2];
rz(-1.9929287) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1956605) q[1];
sx q[1];
rz(-1.6911427) q[1];
sx q[1];
rz(1.2502928) q[1];
rz(1.9367847) q[3];
sx q[3];
rz(-1.9848739) q[3];
sx q[3];
rz(-0.28596349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2953879) q[2];
sx q[2];
rz(-0.69869763) q[2];
sx q[2];
rz(2.6182168) q[2];
rz(2.8335617) q[3];
sx q[3];
rz(-2.568646) q[3];
sx q[3];
rz(-1.3380922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-1.734252) q[0];
sx q[0];
rz(-2.9946406) q[0];
sx q[0];
rz(-1.7379606) q[0];
rz(0.47406468) q[1];
sx q[1];
rz(-1.4539376) q[1];
sx q[1];
rz(-1.7636991) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5869851) q[0];
sx q[0];
rz(-2.4025612) q[0];
sx q[0];
rz(2.2686195) q[0];
rz(-pi) q[1];
x q[1];
rz(0.7147185) q[2];
sx q[2];
rz(-2.208459) q[2];
sx q[2];
rz(-1.4710466) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6973411) q[1];
sx q[1];
rz(-1.8129983) q[1];
sx q[1];
rz(2.7368306) q[1];
rz(1.4078478) q[3];
sx q[3];
rz(-1.019422) q[3];
sx q[3];
rz(-0.58386246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3906117) q[2];
sx q[2];
rz(-1.6335952) q[2];
sx q[2];
rz(-1.5314468) q[2];
rz(1.4577929) q[3];
sx q[3];
rz(-2.0861574) q[3];
sx q[3];
rz(0.84993258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.4703341) q[0];
sx q[0];
rz(-2.8699734) q[0];
sx q[0];
rz(1.097453) q[0];
rz(-0.63256565) q[1];
sx q[1];
rz(-2.0805151) q[1];
sx q[1];
rz(0.17593304) q[1];
rz(1.2420281) q[2];
sx q[2];
rz(-0.9267926) q[2];
sx q[2];
rz(-0.37315858) q[2];
rz(-2.395527) q[3];
sx q[3];
rz(-2.2545771) q[3];
sx q[3];
rz(1.547326) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
