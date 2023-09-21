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
rz(-1.5885408) q[0];
sx q[0];
rz(1.5074402) q[0];
rz(1.5965257) q[1];
sx q[1];
rz(-0.59626055) q[1];
sx q[1];
rz(0.61520666) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0628375) q[0];
sx q[0];
rz(-0.7987928) q[0];
sx q[0];
rz(-0.90057217) q[0];
x q[1];
rz(-2.7841714) q[2];
sx q[2];
rz(-1.3961785) q[2];
sx q[2];
rz(1.071196) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8021009) q[1];
sx q[1];
rz(-1.9994945) q[1];
sx q[1];
rz(-0.93356737) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7513566) q[3];
sx q[3];
rz(-1.2689586) q[3];
sx q[3];
rz(-2.3328822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.4404099) q[2];
sx q[2];
rz(-1.6117233) q[2];
sx q[2];
rz(0.33828503) q[2];
rz(1.7017378) q[3];
sx q[3];
rz(-0.9153291) q[3];
sx q[3];
rz(-2.2556944) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.171339) q[0];
sx q[0];
rz(-0.71115029) q[0];
sx q[0];
rz(3.1112444) q[0];
rz(-3.0753823) q[1];
sx q[1];
rz(-2.1538484) q[1];
sx q[1];
rz(-1.617584) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0328007) q[0];
sx q[0];
rz(-1.5691225) q[0];
sx q[0];
rz(1.3711506) q[0];
x q[1];
rz(-0.020521684) q[2];
sx q[2];
rz(-2.0251209) q[2];
sx q[2];
rz(0.12873912) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6807032) q[1];
sx q[1];
rz(-1.6398805) q[1];
sx q[1];
rz(-2.1417888) q[1];
rz(-pi) q[2];
rz(2.0321235) q[3];
sx q[3];
rz(-0.20878775) q[3];
sx q[3];
rz(1.8148282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.38561884) q[2];
sx q[2];
rz(-2.1531343) q[2];
sx q[2];
rz(-1.1478109) q[2];
rz(-1.3267481) q[3];
sx q[3];
rz(-1.3245405) q[3];
sx q[3];
rz(0.23708788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0148934) q[0];
sx q[0];
rz(-2.6601057) q[0];
sx q[0];
rz(-0.31578627) q[0];
rz(-2.2029927) q[1];
sx q[1];
rz(-1.4626075) q[1];
sx q[1];
rz(2.8895203) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4323498) q[0];
sx q[0];
rz(-0.81626695) q[0];
sx q[0];
rz(-2.4791251) q[0];
x q[1];
rz(-2.4263072) q[2];
sx q[2];
rz(-0.89865696) q[2];
sx q[2];
rz(2.8298024) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7604916) q[1];
sx q[1];
rz(-2.4931506) q[1];
sx q[1];
rz(-1.2566503) q[1];
x q[2];
rz(-0.73022233) q[3];
sx q[3];
rz(-2.632004) q[3];
sx q[3];
rz(2.2024221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0198274) q[2];
sx q[2];
rz(-0.44744197) q[2];
sx q[2];
rz(3.1075409) q[2];
rz(-0.017459067) q[3];
sx q[3];
rz(-1.7826467) q[3];
sx q[3];
rz(-2.0461369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5144192) q[0];
sx q[0];
rz(-1.5486516) q[0];
sx q[0];
rz(-1.5267641) q[0];
rz(-1.0871672) q[1];
sx q[1];
rz(-2.4612869) q[1];
sx q[1];
rz(-0.70708752) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1336466) q[0];
sx q[0];
rz(-0.61230731) q[0];
sx q[0];
rz(-2.4164819) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.25755067) q[2];
sx q[2];
rz(-2.6817245) q[2];
sx q[2];
rz(2.3515153) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.69869631) q[1];
sx q[1];
rz(-2.2203608) q[1];
sx q[1];
rz(-2.991308) q[1];
rz(-1.2624192) q[3];
sx q[3];
rz(-1.7763419) q[3];
sx q[3];
rz(-1.2984848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2924071) q[2];
sx q[2];
rz(-1.9533998) q[2];
sx q[2];
rz(0.46009955) q[2];
rz(-1.397331) q[3];
sx q[3];
rz(-1.5887235) q[3];
sx q[3];
rz(0.24266711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3209155) q[0];
sx q[0];
rz(-0.21235947) q[0];
sx q[0];
rz(1.3943577) q[0];
rz(-1.0955411) q[1];
sx q[1];
rz(-1.54116) q[1];
sx q[1];
rz(-2.8869693) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1063227) q[0];
sx q[0];
rz(-1.5570886) q[0];
sx q[0];
rz(-1.6499707) q[0];
x q[1];
rz(2.1483634) q[2];
sx q[2];
rz(-1.5830056) q[2];
sx q[2];
rz(2.4075367) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9533206) q[1];
sx q[1];
rz(-0.70536648) q[1];
sx q[1];
rz(-2.7642247) q[1];
rz(-2.7084064) q[3];
sx q[3];
rz(-0.90900366) q[3];
sx q[3];
rz(-1.8720686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1191117) q[2];
sx q[2];
rz(-0.20038651) q[2];
sx q[2];
rz(1.7648034) q[2];
rz(1.6453751) q[3];
sx q[3];
rz(-1.6330556) q[3];
sx q[3];
rz(1.0866603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45143932) q[0];
sx q[0];
rz(-2.4017161) q[0];
sx q[0];
rz(-0.29944637) q[0];
rz(-1.0401789) q[1];
sx q[1];
rz(-1.4458011) q[1];
sx q[1];
rz(-0.20656955) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6109989) q[0];
sx q[0];
rz(-1.3337413) q[0];
sx q[0];
rz(2.3321652) q[0];
rz(-2.1913387) q[2];
sx q[2];
rz(-0.62180078) q[2];
sx q[2];
rz(1.3602464) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.088591136) q[1];
sx q[1];
rz(-2.2676761) q[1];
sx q[1];
rz(2.9104396) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2264678) q[3];
sx q[3];
rz(-1.8135999) q[3];
sx q[3];
rz(1.1806928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2999337) q[2];
sx q[2];
rz(-1.724023) q[2];
sx q[2];
rz(-2.858813) q[2];
rz(-0.81280604) q[3];
sx q[3];
rz(-2.7225284) q[3];
sx q[3];
rz(-0.16684428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92641002) q[0];
sx q[0];
rz(-2.0880501) q[0];
sx q[0];
rz(-2.7600631) q[0];
rz(-2.5577257) q[1];
sx q[1];
rz(-0.54324141) q[1];
sx q[1];
rz(1.3279703) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5751942) q[0];
sx q[0];
rz(-1.1325206) q[0];
sx q[0];
rz(-3.0153494) q[0];
rz(0.94857256) q[2];
sx q[2];
rz(-1.186139) q[2];
sx q[2];
rz(-2.4228061) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.59203397) q[1];
sx q[1];
rz(-2.4394819) q[1];
sx q[1];
rz(-1.3105427) q[1];
x q[2];
rz(-2.4626571) q[3];
sx q[3];
rz(-0.62641615) q[3];
sx q[3];
rz(2.7501447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.61838377) q[2];
sx q[2];
rz(-2.0760459) q[2];
sx q[2];
rz(1.3605114) q[2];
rz(-1.4303738) q[3];
sx q[3];
rz(-1.0083219) q[3];
sx q[3];
rz(-0.84806228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1469864) q[0];
sx q[0];
rz(-1.9817579) q[0];
sx q[0];
rz(0.18187901) q[0];
rz(2.6673642) q[1];
sx q[1];
rz(-1.0206181) q[1];
sx q[1];
rz(0.95091933) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5247105) q[0];
sx q[0];
rz(-1.8572154) q[0];
sx q[0];
rz(2.8911203) q[0];
rz(-1.7958926) q[2];
sx q[2];
rz(-2.8689119) q[2];
sx q[2];
rz(-2.8358012) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2086522) q[1];
sx q[1];
rz(-1.3622074) q[1];
sx q[1];
rz(-1.5044466) q[1];
rz(-pi) q[2];
rz(1.9740281) q[3];
sx q[3];
rz(-1.9208761) q[3];
sx q[3];
rz(-2.0389327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9178847) q[2];
sx q[2];
rz(-0.48030883) q[2];
sx q[2];
rz(-3.0656832) q[2];
rz(-2.5935796) q[3];
sx q[3];
rz(-1.3242105) q[3];
sx q[3];
rz(-0.63265911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52371812) q[0];
sx q[0];
rz(-2.0848367) q[0];
sx q[0];
rz(1.7653718) q[0];
rz(-2.7245522) q[1];
sx q[1];
rz(-1.7224256) q[1];
sx q[1];
rz(-0.65972796) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3035559) q[0];
sx q[0];
rz(-2.3065901) q[0];
sx q[0];
rz(-2.3124218) q[0];
x q[1];
rz(-1.3525891) q[2];
sx q[2];
rz(-2.3404684) q[2];
sx q[2];
rz(-1.4584691) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0162504) q[1];
sx q[1];
rz(-0.94031912) q[1];
sx q[1];
rz(-0.36228212) q[1];
x q[2];
rz(1.6660059) q[3];
sx q[3];
rz(-2.6009437) q[3];
sx q[3];
rz(1.8298061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.518121) q[2];
sx q[2];
rz(-2.3770964) q[2];
sx q[2];
rz(-1.0260322) q[2];
rz(0.14885151) q[3];
sx q[3];
rz(-2.1089349) q[3];
sx q[3];
rz(-0.025645105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24348564) q[0];
sx q[0];
rz(-0.74111104) q[0];
sx q[0];
rz(-1.8359258) q[0];
rz(1.9650412) q[1];
sx q[1];
rz(-1.8635609) q[1];
sx q[1];
rz(-1.0356888) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0730597) q[0];
sx q[0];
rz(-1.8951299) q[0];
sx q[0];
rz(2.7102094) q[0];
rz(-pi) q[1];
rz(2.5752441) q[2];
sx q[2];
rz(-0.95745917) q[2];
sx q[2];
rz(-1.7042421) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2346238) q[1];
sx q[1];
rz(-1.3303489) q[1];
sx q[1];
rz(-2.8333227) q[1];
rz(-pi) q[2];
rz(0.28966784) q[3];
sx q[3];
rz(-2.3710459) q[3];
sx q[3];
rz(0.25911261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5130561) q[2];
sx q[2];
rz(-2.3302902) q[2];
sx q[2];
rz(1.9899842) q[2];
rz(0.51268762) q[3];
sx q[3];
rz(-2.0435464) q[3];
sx q[3];
rz(-2.776896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7619027) q[0];
sx q[0];
rz(-1.3544449) q[0];
sx q[0];
rz(0.83308573) q[0];
rz(-1.6336541) q[1];
sx q[1];
rz(-0.58273756) q[1];
sx q[1];
rz(2.6599463) q[1];
rz(-2.7307636) q[2];
sx q[2];
rz(-1.6098235) q[2];
sx q[2];
rz(3.1237684) q[2];
rz(-3.1291943) q[3];
sx q[3];
rz(-2.5231902) q[3];
sx q[3];
rz(-1.4004422) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];