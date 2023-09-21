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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7819408) q[0];
sx q[0];
rz(-1.9160761) q[0];
sx q[0];
rz(0.79202534) q[0];
rz(-2.9801324) q[2];
sx q[2];
rz(-1.4070639) q[2];
sx q[2];
rz(1.2436359) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0238029) q[1];
sx q[1];
rz(-2.1315247) q[1];
sx q[1];
rz(2.6610713) q[1];
rz(-1.2855929) q[3];
sx q[3];
rz(-1.0343026) q[3];
sx q[3];
rz(-0.59629089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2312317) q[2];
sx q[2];
rz(-1.8138764) q[2];
sx q[2];
rz(1.1532016) q[2];
rz(2.6575346) q[3];
sx q[3];
rz(-0.81367937) q[3];
sx q[3];
rz(-2.6822065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83950481) q[0];
sx q[0];
rz(-2.8110101) q[0];
sx q[0];
rz(-2.6385345) q[0];
rz(-1.5548276) q[1];
sx q[1];
rz(-2.4350872) q[1];
sx q[1];
rz(-0.15393004) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0155556) q[0];
sx q[0];
rz(-2.6385348) q[0];
sx q[0];
rz(-1.6004827) q[0];
rz(1.1459848) q[2];
sx q[2];
rz(-2.4467391) q[2];
sx q[2];
rz(1.7244463) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.82984) q[1];
sx q[1];
rz(-2.0318188) q[1];
sx q[1];
rz(0.47383576) q[1];
rz(-1.5699584) q[3];
sx q[3];
rz(-1.0610233) q[3];
sx q[3];
rz(1.2050932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2543891) q[2];
sx q[2];
rz(-0.35787359) q[2];
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
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1948497) q[0];
sx q[0];
rz(-2.0236334) q[0];
sx q[0];
rz(2.2316566) q[0];
rz(-2.3643156) q[1];
sx q[1];
rz(-2.3060019) q[1];
sx q[1];
rz(0.98532239) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9545427) q[0];
sx q[0];
rz(-1.7304725) q[0];
sx q[0];
rz(2.2295582) q[0];
rz(-pi) q[1];
rz(1.1595721) q[2];
sx q[2];
rz(-0.69934884) q[2];
sx q[2];
rz(1.6388091) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4640376) q[1];
sx q[1];
rz(-2.8863393) q[1];
sx q[1];
rz(0.21932253) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3722234) q[3];
sx q[3];
rz(-2.3017075) q[3];
sx q[3];
rz(1.6826671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.862792) q[2];
sx q[2];
rz(-1.1547487) q[2];
sx q[2];
rz(1.8939691) q[2];
rz(-0.4425846) q[3];
sx q[3];
rz(-1.3675888) q[3];
sx q[3];
rz(-2.8201593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2414395) q[0];
sx q[0];
rz(-2.1223919) q[0];
sx q[0];
rz(-2.0181657) q[0];
rz(2.5627047) q[1];
sx q[1];
rz(-1.6826948) q[1];
sx q[1];
rz(1.7480063) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0279562) q[0];
sx q[0];
rz(-0.96158577) q[0];
sx q[0];
rz(2.0834126) q[0];
rz(-2.7921177) q[2];
sx q[2];
rz(-0.45917837) q[2];
sx q[2];
rz(2.0518722) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.046603831) q[1];
sx q[1];
rz(-0.91218439) q[1];
sx q[1];
rz(2.1556426) q[1];
x q[2];
rz(-2.3062069) q[3];
sx q[3];
rz(-1.3961892) q[3];
sx q[3];
rz(1.1268827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1018155) q[2];
sx q[2];
rz(-1.5856051) q[2];
sx q[2];
rz(0.0021136443) q[2];
rz(-0.56143108) q[3];
sx q[3];
rz(-1.0174454) q[3];
sx q[3];
rz(-0.58825618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.026022) q[0];
sx q[0];
rz(-2.0176812) q[0];
sx q[0];
rz(1.2258688) q[0];
rz(1.4670124) q[1];
sx q[1];
rz(-1.8672698) q[1];
sx q[1];
rz(1.3668758) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5850692) q[0];
sx q[0];
rz(-1.9559304) q[0];
sx q[0];
rz(-1.3655846) q[0];
x q[1];
rz(1.6010124) q[2];
sx q[2];
rz(-0.82380166) q[2];
sx q[2];
rz(-0.012416427) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.67147672) q[1];
sx q[1];
rz(-1.6504382) q[1];
sx q[1];
rz(2.2786042) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6392691) q[3];
sx q[3];
rz(-1.886743) q[3];
sx q[3];
rz(-1.5918921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.86429578) q[2];
sx q[2];
rz(-2.0531451) q[2];
sx q[2];
rz(-0.40536353) q[2];
rz(2.6799485) q[3];
sx q[3];
rz(-2.3159537) q[3];
sx q[3];
rz(1.5464787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6500403) q[0];
sx q[0];
rz(-1.4251645) q[0];
sx q[0];
rz(-0.50338411) q[0];
rz(-0.21884306) q[1];
sx q[1];
rz(-1.8914521) q[1];
sx q[1];
rz(-2.4898081) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8021278) q[0];
sx q[0];
rz(-1.3314684) q[0];
sx q[0];
rz(1.375074) q[0];
x q[1];
rz(0.10105614) q[2];
sx q[2];
rz(-0.82436845) q[2];
sx q[2];
rz(1.9161759) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8004605) q[1];
sx q[1];
rz(-0.144185) q[1];
sx q[1];
rz(-1.9393117) q[1];
rz(1.6739453) q[3];
sx q[3];
rz(-1.8522989) q[3];
sx q[3];
rz(2.3867949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.51263222) q[2];
sx q[2];
rz(-1.3971389) q[2];
sx q[2];
rz(2.7499278) q[2];
rz(2.9351249) q[3];
sx q[3];
rz(-2.4075017) q[3];
sx q[3];
rz(0.68968836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-1.3702635) q[0];
sx q[0];
rz(-1.4498793) q[0];
sx q[0];
rz(-2.1719334) q[0];
rz(-2.5580653) q[1];
sx q[1];
rz(-2.0136166) q[1];
sx q[1];
rz(0.13024174) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1657432) q[0];
sx q[0];
rz(-1.2184869) q[0];
sx q[0];
rz(-1.6238814) q[0];
x q[1];
rz(2.1580556) q[2];
sx q[2];
rz(-2.1526255) q[2];
sx q[2];
rz(-2.8772417) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8091734) q[1];
sx q[1];
rz(-1.0671339) q[1];
sx q[1];
rz(1.8196351) q[1];
rz(2.7534915) q[3];
sx q[3];
rz(-2.1183876) q[3];
sx q[3];
rz(2.2263118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.51817259) q[2];
sx q[2];
rz(-1.5600081) q[2];
sx q[2];
rz(-2.0557892) q[2];
rz(0.052224934) q[3];
sx q[3];
rz(-1.6970535) q[3];
sx q[3];
rz(-1.0857371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6770342) q[0];
sx q[0];
rz(-2.1573986) q[0];
sx q[0];
rz(1.1664671) q[0];
rz(2.4160066) q[1];
sx q[1];
rz(-1.8478994) q[1];
sx q[1];
rz(-2.3988147) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9184473) q[0];
sx q[0];
rz(-2.3728275) q[0];
sx q[0];
rz(-2.8630775) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4648415) q[2];
sx q[2];
rz(-1.9988212) q[2];
sx q[2];
rz(-1.8523491) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0726154) q[1];
sx q[1];
rz(-1.9779357) q[1];
sx q[1];
rz(2.3247271) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5524213) q[3];
sx q[3];
rz(-2.1450451) q[3];
sx q[3];
rz(1.9598999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35659197) q[0];
sx q[0];
rz(-0.78105015) q[0];
sx q[0];
rz(-2.2055431) q[0];
rz(-0.014135663) q[1];
sx q[1];
rz(-1.3052669) q[1];
sx q[1];
rz(0.65151185) q[1];
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
x q[1];
rz(-2.3867943) q[2];
sx q[2];
rz(-2.3279466) q[2];
sx q[2];
rz(-0.99572832) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.028101746) q[1];
sx q[1];
rz(-0.3416225) q[1];
sx q[1];
rz(1.9373059) q[1];
rz(-1.2048079) q[3];
sx q[3];
rz(-1.9848739) q[3];
sx q[3];
rz(2.8556292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2953879) q[2];
sx q[2];
rz(-2.442895) q[2];
sx q[2];
rz(0.52337581) q[2];
rz(-0.30803099) q[3];
sx q[3];
rz(-0.57294661) q[3];
sx q[3];
rz(1.3380922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.734252) q[0];
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
rz(2.5705469) q[0];
sx q[0];
rz(-2.018398) q[0];
sx q[0];
rz(0.96121995) q[0];
rz(2.4268742) q[2];
sx q[2];
rz(-0.93313365) q[2];
sx q[2];
rz(-1.4710466) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.1174406) q[1];
sx q[1];
rz(-1.9630868) q[1];
sx q[1];
rz(1.3082318) q[1];
rz(-pi) q[2];
x q[2];
rz(0.25791191) q[3];
sx q[3];
rz(-2.5690418) q[3];
sx q[3];
rz(0.27975988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3906117) q[2];
sx q[2];
rz(-1.6335952) q[2];
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
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-1.8995646) q[2];
sx q[2];
rz(-0.9267926) q[2];
sx q[2];
rz(-0.37315858) q[2];
rz(-0.74606568) q[3];
sx q[3];
rz(-0.88701556) q[3];
sx q[3];
rz(-1.5942667) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
