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
rz(-2.8831557) q[0];
sx q[0];
rz(-2.8539477) q[0];
sx q[0];
rz(-2.0339461) q[0];
rz(-1.8183174) q[1];
sx q[1];
rz(-0.32185093) q[1];
sx q[1];
rz(-2.5133361) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4391651) q[0];
sx q[0];
rz(-1.747923) q[0];
sx q[0];
rz(0.11985531) q[0];
rz(2.7034862) q[2];
sx q[2];
rz(-0.77229653) q[2];
sx q[2];
rz(0.4060678) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9440556) q[1];
sx q[1];
rz(-1.6070131) q[1];
sx q[1];
rz(2.9128068) q[1];
rz(-pi) q[2];
rz(-2.2443549) q[3];
sx q[3];
rz(-0.64770401) q[3];
sx q[3];
rz(-3.0239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2884752) q[2];
sx q[2];
rz(-1.709781) q[2];
sx q[2];
rz(2.8004069) q[2];
rz(-2.2513385) q[3];
sx q[3];
rz(-1.2672) q[3];
sx q[3];
rz(0.53513479) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5219236) q[0];
sx q[0];
rz(-2.4578019) q[0];
sx q[0];
rz(0.71440119) q[0];
rz(-1.5996541) q[1];
sx q[1];
rz(-0.49595141) q[1];
sx q[1];
rz(0.094706789) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9908087) q[0];
sx q[0];
rz(-1.2512733) q[0];
sx q[0];
rz(0.099415778) q[0];
x q[1];
rz(-2.1723198) q[2];
sx q[2];
rz(-1.8516632) q[2];
sx q[2];
rz(-1.8619249) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4058806) q[1];
sx q[1];
rz(-2.7196214) q[1];
sx q[1];
rz(-2.5338245) q[1];
x q[2];
rz(-2.6323039) q[3];
sx q[3];
rz(-1.4319766) q[3];
sx q[3];
rz(-1.2535439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.16227214) q[2];
sx q[2];
rz(-1.8956192) q[2];
sx q[2];
rz(0.47404131) q[2];
rz(0.76215172) q[3];
sx q[3];
rz(-0.63665974) q[3];
sx q[3];
rz(2.6704085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8629465) q[0];
sx q[0];
rz(-0.11676783) q[0];
sx q[0];
rz(-0.85292029) q[0];
rz(-2.9184753) q[1];
sx q[1];
rz(-0.47839636) q[1];
sx q[1];
rz(-2.6257637) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86403041) q[0];
sx q[0];
rz(-1.6245961) q[0];
sx q[0];
rz(-0.050431099) q[0];
rz(-pi) q[1];
rz(0.25571172) q[2];
sx q[2];
rz(-2.1923794) q[2];
sx q[2];
rz(2.5598124) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8341537) q[1];
sx q[1];
rz(-2.1527228) q[1];
sx q[1];
rz(0.80640275) q[1];
rz(-pi) q[2];
rz(-2.2688341) q[3];
sx q[3];
rz(-1.3036733) q[3];
sx q[3];
rz(1.1760548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3333266) q[2];
sx q[2];
rz(-0.9767248) q[2];
sx q[2];
rz(-0.094956368) q[2];
rz(0.44148463) q[3];
sx q[3];
rz(-1.7850826) q[3];
sx q[3];
rz(-1.9879742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3917291) q[0];
sx q[0];
rz(-0.57498217) q[0];
sx q[0];
rz(-2.0854501) q[0];
rz(-0.9437584) q[1];
sx q[1];
rz(-2.8902003) q[1];
sx q[1];
rz(-1.4518849) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75987923) q[0];
sx q[0];
rz(-0.96764123) q[0];
sx q[0];
rz(-2.2908014) q[0];
x q[1];
rz(0.69047549) q[2];
sx q[2];
rz(-2.568141) q[2];
sx q[2];
rz(-0.8586463) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8274967) q[1];
sx q[1];
rz(-2.8982867) q[1];
sx q[1];
rz(2.9929586) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.11066505) q[3];
sx q[3];
rz(-1.5783211) q[3];
sx q[3];
rz(-1.2570109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4132495) q[2];
sx q[2];
rz(-2.3409833) q[2];
sx q[2];
rz(0.70017868) q[2];
rz(-2.2705196) q[3];
sx q[3];
rz(-2.0756523) q[3];
sx q[3];
rz(1.7843436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0102608) q[0];
sx q[0];
rz(-2.6483674) q[0];
sx q[0];
rz(-0.48466551) q[0];
rz(0.87801814) q[1];
sx q[1];
rz(-1.1596102) q[1];
sx q[1];
rz(2.7427618) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53610401) q[0];
sx q[0];
rz(-1.7057422) q[0];
sx q[0];
rz(0.55079726) q[0];
x q[1];
rz(-2.7473248) q[2];
sx q[2];
rz(-2.3986001) q[2];
sx q[2];
rz(0.34994469) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.64885253) q[1];
sx q[1];
rz(-1.938174) q[1];
sx q[1];
rz(-2.2671764) q[1];
rz(-pi) q[2];
rz(0.54581235) q[3];
sx q[3];
rz(-0.94792507) q[3];
sx q[3];
rz(0.92178492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.47348076) q[2];
sx q[2];
rz(-1.5405737) q[2];
sx q[2];
rz(-0.87781805) q[2];
rz(-0.38464883) q[3];
sx q[3];
rz(-0.84482241) q[3];
sx q[3];
rz(1.9627242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-0.28410742) q[0];
sx q[0];
rz(-0.50898886) q[0];
sx q[0];
rz(0.2704764) q[0];
rz(0.30666223) q[1];
sx q[1];
rz(-0.51391927) q[1];
sx q[1];
rz(2.564863) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3535296) q[0];
sx q[0];
rz(-1.3856141) q[0];
sx q[0];
rz(-0.089140541) q[0];
rz(0.88097303) q[2];
sx q[2];
rz(-1.2861797) q[2];
sx q[2];
rz(1.5280444) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.57177793) q[1];
sx q[1];
rz(-1.0228436) q[1];
sx q[1];
rz(0.46787365) q[1];
rz(-pi) q[2];
rz(0.0402952) q[3];
sx q[3];
rz(-2.8065805) q[3];
sx q[3];
rz(-2.1744436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3660672) q[2];
sx q[2];
rz(-2.0719353) q[2];
sx q[2];
rz(0.97037399) q[2];
rz(0.35106418) q[3];
sx q[3];
rz(-0.23215663) q[3];
sx q[3];
rz(-0.27950132) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2358667) q[0];
sx q[0];
rz(-2.7577363) q[0];
sx q[0];
rz(-0.50091499) q[0];
rz(-0.72055912) q[1];
sx q[1];
rz(-0.84545285) q[1];
sx q[1];
rz(0.93773425) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36762992) q[0];
sx q[0];
rz(-1.7741307) q[0];
sx q[0];
rz(3.0995661) q[0];
rz(-pi) q[1];
x q[1];
rz(0.47368424) q[2];
sx q[2];
rz(-1.9100744) q[2];
sx q[2];
rz(-0.35655278) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5825469) q[1];
sx q[1];
rz(-1.9314879) q[1];
sx q[1];
rz(2.2466895) q[1];
x q[2];
rz(-0.012091919) q[3];
sx q[3];
rz(-0.68688697) q[3];
sx q[3];
rz(1.6364607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.131989) q[2];
sx q[2];
rz(-2.9961573) q[2];
sx q[2];
rz(-0.84849882) q[2];
rz(-2.2961473) q[3];
sx q[3];
rz(-0.65631056) q[3];
sx q[3];
rz(1.9158624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7312412) q[0];
sx q[0];
rz(-0.9587962) q[0];
sx q[0];
rz(-2.248726) q[0];
rz(-0.061575312) q[1];
sx q[1];
rz(-2.4696746) q[1];
sx q[1];
rz(-2.9796013) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6380986) q[0];
sx q[0];
rz(-1.5862521) q[0];
sx q[0];
rz(-1.6729805) q[0];
rz(-pi) q[1];
rz(1.4561557) q[2];
sx q[2];
rz(-0.87213665) q[2];
sx q[2];
rz(-1.6841152) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2714241) q[1];
sx q[1];
rz(-1.4335911) q[1];
sx q[1];
rz(2.8297564) q[1];
rz(-pi) q[2];
rz(1.4668767) q[3];
sx q[3];
rz(-1.6641518) q[3];
sx q[3];
rz(0.68478497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.091247678) q[2];
sx q[2];
rz(-0.42751905) q[2];
sx q[2];
rz(2.4166935) q[2];
rz(2.8643518) q[3];
sx q[3];
rz(-2.1767949) q[3];
sx q[3];
rz(3.0550756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8898833) q[0];
sx q[0];
rz(-0.6530264) q[0];
sx q[0];
rz(-3.1157893) q[0];
rz(-1.7426527) q[1];
sx q[1];
rz(-1.6394697) q[1];
sx q[1];
rz(0.10765156) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.066641971) q[0];
sx q[0];
rz(-2.5558976) q[0];
sx q[0];
rz(1.8138712) q[0];
rz(-pi) q[1];
rz(2.669966) q[2];
sx q[2];
rz(-1.9666858) q[2];
sx q[2];
rz(2.5079923) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3947666) q[1];
sx q[1];
rz(-2.0438384) q[1];
sx q[1];
rz(-1.5302782) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6844382) q[3];
sx q[3];
rz(-1.5530619) q[3];
sx q[3];
rz(2.4271698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4699576) q[2];
sx q[2];
rz(-1.8822972) q[2];
sx q[2];
rz(1.0939481) q[2];
rz(-2.5661902) q[3];
sx q[3];
rz(-2.3016774) q[3];
sx q[3];
rz(1.1168787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6082918) q[0];
sx q[0];
rz(-2.5550483) q[0];
sx q[0];
rz(2.4285512) q[0];
rz(-0.22480045) q[1];
sx q[1];
rz(-0.59200042) q[1];
sx q[1];
rz(2.0483268) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3521234) q[0];
sx q[0];
rz(-0.88406815) q[0];
sx q[0];
rz(-2.5144469) q[0];
x q[1];
rz(1.8453127) q[2];
sx q[2];
rz(-1.7651254) q[2];
sx q[2];
rz(-0.39482612) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0436127) q[1];
sx q[1];
rz(-0.47046767) q[1];
sx q[1];
rz(-2.8941675) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8328463) q[3];
sx q[3];
rz(-0.42632494) q[3];
sx q[3];
rz(0.48708068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.1405868) q[2];
sx q[2];
rz(-2.4457377) q[2];
sx q[2];
rz(-0.33983964) q[2];
rz(-3.0991683) q[3];
sx q[3];
rz(-1.8570615) q[3];
sx q[3];
rz(2.5137918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9545659) q[0];
sx q[0];
rz(-1.5196336) q[0];
sx q[0];
rz(-1.2484311) q[0];
rz(-2.3595702) q[1];
sx q[1];
rz(-2.2047058) q[1];
sx q[1];
rz(2.4155736) q[1];
rz(1.2371042) q[2];
sx q[2];
rz(-1.9241698) q[2];
sx q[2];
rz(-2.4826222) q[2];
rz(-1.9966765) q[3];
sx q[3];
rz(-1.799634) q[3];
sx q[3];
rz(-1.4992731) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
