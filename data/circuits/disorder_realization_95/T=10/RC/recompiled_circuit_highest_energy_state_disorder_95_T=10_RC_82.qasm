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
rz(0.79701841) q[0];
sx q[0];
rz(-2.2198644) q[0];
sx q[0];
rz(-1.8812688) q[0];
rz(-2.4951275) q[1];
sx q[1];
rz(-0.87352455) q[1];
sx q[1];
rz(-2.8096107) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9063569) q[0];
sx q[0];
rz(-0.42354326) q[0];
sx q[0];
rz(1.5111501) q[0];
rz(-pi) q[1];
rz(-1.9464363) q[2];
sx q[2];
rz(-2.2323998) q[2];
sx q[2];
rz(-2.313569) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5378137) q[1];
sx q[1];
rz(-2.3821476) q[1];
sx q[1];
rz(-0.045995637) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.1370955) q[3];
sx q[3];
rz(-2.6957316) q[3];
sx q[3];
rz(2.9112828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.476764) q[2];
sx q[2];
rz(-2.2085184) q[2];
sx q[2];
rz(-2.942371) q[2];
rz(0.18579379) q[3];
sx q[3];
rz(-1.4150861) q[3];
sx q[3];
rz(-1.7004405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81714565) q[0];
sx q[0];
rz(-0.48878336) q[0];
sx q[0];
rz(2.4031438) q[0];
rz(1.6959408) q[1];
sx q[1];
rz(-1.4001458) q[1];
sx q[1];
rz(0.48283985) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8017977) q[0];
sx q[0];
rz(-1.4301824) q[0];
sx q[0];
rz(-1.7231581) q[0];
rz(-pi) q[1];
rz(0.17530967) q[2];
sx q[2];
rz(-1.2273437) q[2];
sx q[2];
rz(0.94409787) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9571217) q[1];
sx q[1];
rz(-1.0252684) q[1];
sx q[1];
rz(-2.1197182) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0121423) q[3];
sx q[3];
rz(-1.4317272) q[3];
sx q[3];
rz(0.46848759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8740251) q[2];
sx q[2];
rz(-1.4635307) q[2];
sx q[2];
rz(-1.2919424) q[2];
rz(1.1235631) q[3];
sx q[3];
rz(-0.70190391) q[3];
sx q[3];
rz(-1.4152214) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32347754) q[0];
sx q[0];
rz(-2.1577142) q[0];
sx q[0];
rz(-1.4418607) q[0];
rz(-1.2280751) q[1];
sx q[1];
rz(-1.9107198) q[1];
sx q[1];
rz(-2.0679811) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.613205) q[0];
sx q[0];
rz(-2.3484485) q[0];
sx q[0];
rz(-2.6499477) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3153438) q[2];
sx q[2];
rz(-2.4702304) q[2];
sx q[2];
rz(-2.3731212) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0328999) q[1];
sx q[1];
rz(-1.5376904) q[1];
sx q[1];
rz(-2.2733403) q[1];
rz(-1.7341679) q[3];
sx q[3];
rz(-2.022812) q[3];
sx q[3];
rz(2.3491835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.29232612) q[2];
sx q[2];
rz(-0.30305114) q[2];
sx q[2];
rz(-2.1185875) q[2];
rz(-0.060221378) q[3];
sx q[3];
rz(-1.9950208) q[3];
sx q[3];
rz(-2.4165418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66134727) q[0];
sx q[0];
rz(-2.9954973) q[0];
sx q[0];
rz(-0.51937854) q[0];
rz(-2.5410779) q[1];
sx q[1];
rz(-1.0888211) q[1];
sx q[1];
rz(0.75469887) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8225427) q[0];
sx q[0];
rz(-0.62915914) q[0];
sx q[0];
rz(2.1436917) q[0];
rz(-pi) q[1];
rz(2.9982996) q[2];
sx q[2];
rz(-0.87277647) q[2];
sx q[2];
rz(2.3868274) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3364611) q[1];
sx q[1];
rz(-1.0829003) q[1];
sx q[1];
rz(-0.35780235) q[1];
rz(2.5537476) q[3];
sx q[3];
rz(-0.64661542) q[3];
sx q[3];
rz(0.69348303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7924013) q[2];
sx q[2];
rz(-2.3442522) q[2];
sx q[2];
rz(0.1758197) q[2];
rz(1.1261806) q[3];
sx q[3];
rz(-1.3577941) q[3];
sx q[3];
rz(-0.064182909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46614161) q[0];
sx q[0];
rz(-1.4956681) q[0];
sx q[0];
rz(-2.5597036) q[0];
rz(-1.1833082) q[1];
sx q[1];
rz(-0.71154037) q[1];
sx q[1];
rz(-1.2540832) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3235856) q[0];
sx q[0];
rz(-0.76420751) q[0];
sx q[0];
rz(-1.0952522) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0462363) q[2];
sx q[2];
rz(-2.0031906) q[2];
sx q[2];
rz(0.43579416) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5820133) q[1];
sx q[1];
rz(-1.8648284) q[1];
sx q[1];
rz(-0.21802417) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3816816) q[3];
sx q[3];
rz(-1.0737891) q[3];
sx q[3];
rz(-0.86732098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9767655) q[2];
sx q[2];
rz(-0.81563121) q[2];
sx q[2];
rz(0.35484472) q[2];
rz(-1.1497078) q[3];
sx q[3];
rz(-1.0747654) q[3];
sx q[3];
rz(2.3000075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.069816291) q[0];
sx q[0];
rz(-0.27158296) q[0];
sx q[0];
rz(3.1014882) q[0];
rz(-1.6768203) q[1];
sx q[1];
rz(-0.96950871) q[1];
sx q[1];
rz(-2.6693595) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.036333648) q[0];
sx q[0];
rz(-1.3091631) q[0];
sx q[0];
rz(-2.9519269) q[0];
x q[1];
rz(1.5750999) q[2];
sx q[2];
rz(-0.65281463) q[2];
sx q[2];
rz(-3.0497361) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.66989723) q[1];
sx q[1];
rz(-1.9006839) q[1];
sx q[1];
rz(2.8389588) q[1];
x q[2];
rz(-0.29904265) q[3];
sx q[3];
rz(-2.3573993) q[3];
sx q[3];
rz(1.2687781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.27452305) q[2];
sx q[2];
rz(-2.3752866) q[2];
sx q[2];
rz(1.6609894) q[2];
rz(-3.1001422) q[3];
sx q[3];
rz(-0.71235123) q[3];
sx q[3];
rz(-2.1672772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.239045) q[0];
sx q[0];
rz(-2.3426265) q[0];
sx q[0];
rz(2.3833185) q[0];
rz(0.43765086) q[1];
sx q[1];
rz(-1.9270908) q[1];
sx q[1];
rz(-0.95380107) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7051288) q[0];
sx q[0];
rz(-2.7565694) q[0];
sx q[0];
rz(1.692311) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1787011) q[2];
sx q[2];
rz(-2.1959119) q[2];
sx q[2];
rz(-1.4639548) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.055102417) q[1];
sx q[1];
rz(-2.2728572) q[1];
sx q[1];
rz(0.081900077) q[1];
rz(-0.99178828) q[3];
sx q[3];
rz(-0.48509337) q[3];
sx q[3];
rz(2.4320784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4147676) q[2];
sx q[2];
rz(-0.32764062) q[2];
sx q[2];
rz(0.3978351) q[2];
rz(-1.2658524) q[3];
sx q[3];
rz(-1.7222907) q[3];
sx q[3];
rz(-0.46441594) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1606814) q[0];
sx q[0];
rz(-2.2137764) q[0];
sx q[0];
rz(0.30447793) q[0];
rz(-1.132384) q[1];
sx q[1];
rz(-2.4696923) q[1];
sx q[1];
rz(-1.0736046) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17504263) q[0];
sx q[0];
rz(-1.5877519) q[0];
sx q[0];
rz(1.5477563) q[0];
rz(-1.2108285) q[2];
sx q[2];
rz(-1.1774633) q[2];
sx q[2];
rz(-0.29811146) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.5636812) q[1];
sx q[1];
rz(-1.328224) q[1];
sx q[1];
rz(-2.1798308) q[1];
rz(-pi) q[2];
rz(-2.497914) q[3];
sx q[3];
rz(-1.5305909) q[3];
sx q[3];
rz(-2.5098367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.73211804) q[2];
sx q[2];
rz(-1.2827337) q[2];
sx q[2];
rz(0.049962433) q[2];
rz(2.2299855) q[3];
sx q[3];
rz(-0.92602366) q[3];
sx q[3];
rz(0.91134206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53495812) q[0];
sx q[0];
rz(-1.7462523) q[0];
sx q[0];
rz(-0.362679) q[0];
rz(0.72049385) q[1];
sx q[1];
rz(-0.81169218) q[1];
sx q[1];
rz(-1.1788751) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76438457) q[0];
sx q[0];
rz(-2.3297455) q[0];
sx q[0];
rz(-0.7332515) q[0];
rz(-pi) q[1];
rz(-2.9083546) q[2];
sx q[2];
rz(-2.812603) q[2];
sx q[2];
rz(1.1057265) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3549344) q[1];
sx q[1];
rz(-1.7410856) q[1];
sx q[1];
rz(3.0436467) q[1];
rz(-pi) q[2];
x q[2];
rz(0.54525996) q[3];
sx q[3];
rz(-1.2125748) q[3];
sx q[3];
rz(-0.48945603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1371586) q[2];
sx q[2];
rz(-0.68778554) q[2];
sx q[2];
rz(1.9384109) q[2];
rz(0.016599003) q[3];
sx q[3];
rz(-1.4998452) q[3];
sx q[3];
rz(-1.871292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5430629) q[0];
sx q[0];
rz(-2.6180584) q[0];
sx q[0];
rz(1.6981) q[0];
rz(0.47830018) q[1];
sx q[1];
rz(-2.4595478) q[1];
sx q[1];
rz(-1.2102478) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87916527) q[0];
sx q[0];
rz(-1.0813923) q[0];
sx q[0];
rz(0.64206815) q[0];
rz(-pi) q[1];
rz(-1.7487594) q[2];
sx q[2];
rz(-2.2134668) q[2];
sx q[2];
rz(-0.27163525) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5961884) q[1];
sx q[1];
rz(-2.8126723) q[1];
sx q[1];
rz(2.9795582) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5555672) q[3];
sx q[3];
rz(-2.7507493) q[3];
sx q[3];
rz(1.217921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9532507) q[2];
sx q[2];
rz(-0.61389273) q[2];
sx q[2];
rz(-2.4284412) q[2];
rz(0.6330511) q[3];
sx q[3];
rz(-2.1745067) q[3];
sx q[3];
rz(0.87575325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1322094) q[0];
sx q[0];
rz(-1.3306946) q[0];
sx q[0];
rz(2.7015986) q[0];
rz(-2.412759) q[1];
sx q[1];
rz(-0.76580096) q[1];
sx q[1];
rz(-2.0892807) q[1];
rz(-1.8447137) q[2];
sx q[2];
rz(-0.40362457) q[2];
sx q[2];
rz(0.85164126) q[2];
rz(-2.3340529) q[3];
sx q[3];
rz(-1.836946) q[3];
sx q[3];
rz(2.436773) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
