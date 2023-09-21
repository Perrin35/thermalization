OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2774529) q[0];
sx q[0];
rz(4.6946445) q[0];
sx q[0];
rz(10.932218) q[0];
rz(-1.545067) q[1];
sx q[1];
rz(-2.5453321) q[1];
sx q[1];
rz(2.526386) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9278487) q[0];
sx q[0];
rz(-0.97457492) q[0];
sx q[0];
rz(-2.5736789) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7841714) q[2];
sx q[2];
rz(-1.7454141) q[2];
sx q[2];
rz(-1.071196) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2088036) q[1];
sx q[1];
rz(-2.1425769) q[1];
sx q[1];
rz(-0.51704452) q[1];
x q[2];
rz(2.6184222) q[3];
sx q[3];
rz(-0.35029951) q[3];
sx q[3];
rz(-1.3594128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.4404099) q[2];
sx q[2];
rz(-1.6117233) q[2];
sx q[2];
rz(0.33828503) q[2];
rz(-1.7017378) q[3];
sx q[3];
rz(-0.9153291) q[3];
sx q[3];
rz(2.2556944) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97025362) q[0];
sx q[0];
rz(-2.4304424) q[0];
sx q[0];
rz(0.030348226) q[0];
rz(0.066210315) q[1];
sx q[1];
rz(-0.98774424) q[1];
sx q[1];
rz(1.617584) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0328007) q[0];
sx q[0];
rz(-1.5724702) q[0];
sx q[0];
rz(1.770442) q[0];
rz(-pi) q[1];
rz(1.6127869) q[2];
sx q[2];
rz(-2.6868372) q[2];
sx q[2];
rz(0.08200478) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1387353) q[1];
sx q[1];
rz(-0.5746952) q[1];
sx q[1];
rz(-1.4434622) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3832983) q[3];
sx q[3];
rz(-1.6631931) q[3];
sx q[3];
rz(-0.20860162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7559738) q[2];
sx q[2];
rz(-0.9884584) q[2];
sx q[2];
rz(1.9937817) q[2];
rz(1.8148445) q[3];
sx q[3];
rz(-1.8170522) q[3];
sx q[3];
rz(-0.23708788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0148934) q[0];
sx q[0];
rz(-0.48148695) q[0];
sx q[0];
rz(0.31578627) q[0];
rz(-2.2029927) q[1];
sx q[1];
rz(-1.4626075) q[1];
sx q[1];
rz(2.8895203) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2825851) q[0];
sx q[0];
rz(-0.95882817) q[0];
sx q[0];
rz(-0.99143272) q[0];
rz(-pi) q[1];
rz(-2.3825233) q[2];
sx q[2];
rz(-2.1096863) q[2];
sx q[2];
rz(-0.76314229) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.38110106) q[1];
sx q[1];
rz(-2.4931506) q[1];
sx q[1];
rz(1.8849424) q[1];
rz(-pi) q[2];
rz(-0.73022233) q[3];
sx q[3];
rz(-0.50958868) q[3];
sx q[3];
rz(0.93917055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0198274) q[2];
sx q[2];
rz(-2.6941507) q[2];
sx q[2];
rz(3.1075409) q[2];
rz(0.017459067) q[3];
sx q[3];
rz(-1.7826467) q[3];
sx q[3];
rz(-1.0954558) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62717342) q[0];
sx q[0];
rz(-1.592941) q[0];
sx q[0];
rz(-1.6148286) q[0];
rz(2.0544255) q[1];
sx q[1];
rz(-0.68030578) q[1];
sx q[1];
rz(0.70708752) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93543816) q[0];
sx q[0];
rz(-1.9618789) q[0];
sx q[0];
rz(2.6576256) q[0];
rz(-pi) q[1];
rz(2.884042) q[2];
sx q[2];
rz(-0.45986816) q[2];
sx q[2];
rz(0.79007733) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1781588) q[1];
sx q[1];
rz(-1.4512832) q[1];
sx q[1];
rz(0.91576373) q[1];
x q[2];
rz(-0.96890038) q[3];
sx q[3];
rz(-0.36877353) q[3];
sx q[3];
rz(-2.8440648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.84918555) q[2];
sx q[2];
rz(-1.9533998) q[2];
sx q[2];
rz(2.6814931) q[2];
rz(1.397331) q[3];
sx q[3];
rz(-1.5528691) q[3];
sx q[3];
rz(-2.8989255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3209155) q[0];
sx q[0];
rz(-0.21235947) q[0];
sx q[0];
rz(-1.3943577) q[0];
rz(-2.0460515) q[1];
sx q[1];
rz(-1.54116) q[1];
sx q[1];
rz(-0.25462338) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29339644) q[0];
sx q[0];
rz(-0.080349803) q[0];
sx q[0];
rz(1.3991762) q[0];
x q[1];
rz(0.014572797) q[2];
sx q[2];
rz(-0.99327786) q[2];
sx q[2];
rz(-2.2968959) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6756145) q[1];
sx q[1];
rz(-1.8120159) q[1];
sx q[1];
rz(0.66959186) q[1];
rz(-pi) q[2];
x q[2];
rz(0.43318627) q[3];
sx q[3];
rz(-0.90900366) q[3];
sx q[3];
rz(1.2695241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.022481) q[2];
sx q[2];
rz(-2.9412061) q[2];
sx q[2];
rz(1.7648034) q[2];
rz(-1.4962176) q[3];
sx q[3];
rz(-1.508537) q[3];
sx q[3];
rz(-1.0866603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45143932) q[0];
sx q[0];
rz(-0.7398766) q[0];
sx q[0];
rz(-0.29944637) q[0];
rz(-1.0401789) q[1];
sx q[1];
rz(-1.6957915) q[1];
sx q[1];
rz(-2.9350231) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5305938) q[0];
sx q[0];
rz(-1.8078513) q[0];
sx q[0];
rz(-0.80942746) q[0];
rz(-1.0429522) q[2];
sx q[2];
rz(-1.9163418) q[2];
sx q[2];
rz(0.31574677) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.2628852) q[1];
sx q[1];
rz(-2.4135114) q[1];
sx q[1];
rz(1.3036149) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2264678) q[3];
sx q[3];
rz(-1.3279928) q[3];
sx q[3];
rz(-1.1806928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.841659) q[2];
sx q[2];
rz(-1.724023) q[2];
sx q[2];
rz(2.858813) q[2];
rz(-2.3287866) q[3];
sx q[3];
rz(-2.7225284) q[3];
sx q[3];
rz(-2.9747484) q[3];
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
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2151826) q[0];
sx q[0];
rz(-1.0535425) q[0];
sx q[0];
rz(2.7600631) q[0];
rz(-0.58386699) q[1];
sx q[1];
rz(-0.54324141) q[1];
sx q[1];
rz(-1.3279703) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.049411557) q[0];
sx q[0];
rz(-1.4565399) q[0];
sx q[0];
rz(1.1294424) q[0];
rz(-pi) q[1];
rz(-0.46220025) q[2];
sx q[2];
rz(-1.000058) q[2];
sx q[2];
rz(-2.5525023) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.77818645) q[1];
sx q[1];
rz(-1.4038329) q[1];
sx q[1];
rz(0.88552514) q[1];
rz(-pi) q[2];
rz(2.4626571) q[3];
sx q[3];
rz(-0.62641615) q[3];
sx q[3];
rz(-2.7501447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5232089) q[2];
sx q[2];
rz(-1.0655468) q[2];
sx q[2];
rz(-1.7810812) q[2];
rz(-1.7112188) q[3];
sx q[3];
rz(-1.0083219) q[3];
sx q[3];
rz(0.84806228) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1469864) q[0];
sx q[0];
rz(-1.1598347) q[0];
sx q[0];
rz(-2.9597136) q[0];
rz(2.6673642) q[1];
sx q[1];
rz(-1.0206181) q[1];
sx q[1];
rz(-2.1906733) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5247105) q[0];
sx q[0];
rz(-1.8572154) q[0];
sx q[0];
rz(-0.25047238) q[0];
rz(-pi) q[1];
rz(3.0792564) q[2];
sx q[2];
rz(-1.3051635) q[2];
sx q[2];
rz(-3.0692284) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.89815318) q[1];
sx q[1];
rz(-2.9228518) q[1];
sx q[1];
rz(0.30355138) q[1];
x q[2];
rz(-0.37787921) q[3];
sx q[3];
rz(-1.9482908) q[3];
sx q[3];
rz(-0.32285238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2237079) q[2];
sx q[2];
rz(-2.6612838) q[2];
sx q[2];
rz(0.075909464) q[2];
rz(0.54801303) q[3];
sx q[3];
rz(-1.3242105) q[3];
sx q[3];
rz(-0.63265911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6178745) q[0];
sx q[0];
rz(-2.0848367) q[0];
sx q[0];
rz(1.3762208) q[0];
rz(-0.41704047) q[1];
sx q[1];
rz(-1.7224256) q[1];
sx q[1];
rz(-2.4818647) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83803672) q[0];
sx q[0];
rz(-0.8350026) q[0];
sx q[0];
rz(2.3124218) q[0];
x q[1];
rz(0.21978901) q[2];
sx q[2];
rz(-2.347749) q[2];
sx q[2];
rz(1.1500051) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1253423) q[1];
sx q[1];
rz(-2.2012735) q[1];
sx q[1];
rz(0.36228212) q[1];
x q[2];
rz(-3.0845853) q[3];
sx q[3];
rz(-1.032864) q[3];
sx q[3];
rz(-1.2008592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.62347162) q[2];
sx q[2];
rz(-2.3770964) q[2];
sx q[2];
rz(1.0260322) q[2];
rz(-2.9927411) q[3];
sx q[3];
rz(-2.1089349) q[3];
sx q[3];
rz(-0.025645105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24348564) q[0];
sx q[0];
rz(-2.4004816) q[0];
sx q[0];
rz(-1.3056668) q[0];
rz(1.1765515) q[1];
sx q[1];
rz(-1.2780317) q[1];
sx q[1];
rz(-1.0356888) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0730597) q[0];
sx q[0];
rz(-1.8951299) q[0];
sx q[0];
rz(-2.7102094) q[0];
rz(-0.87558526) q[2];
sx q[2];
rz(-1.1165809) q[2];
sx q[2];
rz(-2.9241965) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.021691572) q[1];
sx q[1];
rz(-2.7530115) q[1];
sx q[1];
rz(-0.67967023) q[1];
rz(1.3003179) q[3];
sx q[3];
rz(-2.3016553) q[3];
sx q[3];
rz(-2.4887816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5130561) q[2];
sx q[2];
rz(-0.81130242) q[2];
sx q[2];
rz(-1.9899842) q[2];
rz(-2.628905) q[3];
sx q[3];
rz(-2.0435464) q[3];
sx q[3];
rz(-2.776896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7619027) q[0];
sx q[0];
rz(-1.3544449) q[0];
sx q[0];
rz(0.83308573) q[0];
rz(1.5079386) q[1];
sx q[1];
rz(-0.58273756) q[1];
sx q[1];
rz(2.6599463) q[1];
rz(-1.528231) q[2];
sx q[2];
rz(-1.1602989) q[2];
sx q[2];
rz(1.5699671) q[2];
rz(0.61836615) q[3];
sx q[3];
rz(-1.5779839) q[3];
sx q[3];
rz(0.18045651) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];