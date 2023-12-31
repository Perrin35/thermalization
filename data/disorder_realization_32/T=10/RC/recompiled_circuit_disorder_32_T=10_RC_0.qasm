OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.10652868) q[0];
sx q[0];
rz(-1.0892692) q[0];
sx q[0];
rz(0.16103345) q[0];
rz(1.610202) q[1];
sx q[1];
rz(-0.47709563) q[1];
sx q[1];
rz(0.49638003) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3408605) q[0];
sx q[0];
rz(-0.98156089) q[0];
sx q[0];
rz(-0.13829921) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7726937) q[2];
sx q[2];
rz(-1.2374094) q[2];
sx q[2];
rz(0.29793973) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.1094184) q[1];
sx q[1];
rz(-2.2757581) q[1];
sx q[1];
rz(-1.6900307) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6707889) q[3];
sx q[3];
rz(-1.9766207) q[3];
sx q[3];
rz(2.0506746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6266142) q[2];
sx q[2];
rz(-0.86003059) q[2];
sx q[2];
rz(-0.28764763) q[2];
rz(1.3927762) q[3];
sx q[3];
rz(-0.98373047) q[3];
sx q[3];
rz(-1.1028517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2709133) q[0];
sx q[0];
rz(-2.3864855) q[0];
sx q[0];
rz(-0.57759181) q[0];
rz(1.5354935) q[1];
sx q[1];
rz(-2.0237193) q[1];
sx q[1];
rz(1.577852) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6710885) q[0];
sx q[0];
rz(-1.7894063) q[0];
sx q[0];
rz(1.6862306) q[0];
rz(-2.4096476) q[2];
sx q[2];
rz(-0.39790301) q[2];
sx q[2];
rz(-1.0674455) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6985804) q[1];
sx q[1];
rz(-2.6786782) q[1];
sx q[1];
rz(1.7673311) q[1];
rz(-pi) q[2];
rz(2.9868449) q[3];
sx q[3];
rz(-1.7259211) q[3];
sx q[3];
rz(2.250092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.10721283) q[2];
sx q[2];
rz(-2.1652174) q[2];
sx q[2];
rz(-1.0822901) q[2];
rz(-1.1632495) q[3];
sx q[3];
rz(-1.9415559) q[3];
sx q[3];
rz(2.5015586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0796233) q[0];
sx q[0];
rz(-2.847023) q[0];
sx q[0];
rz(2.9586155) q[0];
rz(-3.0984763) q[1];
sx q[1];
rz(-2.1886107) q[1];
sx q[1];
rz(-1.9972237) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2079175) q[0];
sx q[0];
rz(-1.1753923) q[0];
sx q[0];
rz(1.0008706) q[0];
rz(-pi) q[1];
x q[1];
rz(1.794533) q[2];
sx q[2];
rz(-1.1615331) q[2];
sx q[2];
rz(-2.4350016) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.548) q[1];
sx q[1];
rz(-1.9152194) q[1];
sx q[1];
rz(-2.6542632) q[1];
rz(-pi) q[2];
rz(0.010300962) q[3];
sx q[3];
rz(-2.8539538) q[3];
sx q[3];
rz(-1.079263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.019505067) q[2];
sx q[2];
rz(-1.1818037) q[2];
sx q[2];
rz(0.90467492) q[2];
rz(-2.4217862) q[3];
sx q[3];
rz(-2.7147839) q[3];
sx q[3];
rz(-2.3864746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4633789) q[0];
sx q[0];
rz(-2.1713874) q[0];
sx q[0];
rz(-2.4257207) q[0];
rz(-0.32575682) q[1];
sx q[1];
rz(-1.0595067) q[1];
sx q[1];
rz(-0.99266565) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0621322) q[0];
sx q[0];
rz(-1.8347164) q[0];
sx q[0];
rz(0.53702766) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5771192) q[2];
sx q[2];
rz(-1.0130922) q[2];
sx q[2];
rz(-2.0828473) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0996453) q[1];
sx q[1];
rz(-1.4775659) q[1];
sx q[1];
rz(-0.16361841) q[1];
x q[2];
rz(2.1487001) q[3];
sx q[3];
rz(-1.9979457) q[3];
sx q[3];
rz(0.072673365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0830393) q[2];
sx q[2];
rz(-0.71648592) q[2];
sx q[2];
rz(-1.7129718) q[2];
rz(0.84609091) q[3];
sx q[3];
rz(-1.6276136) q[3];
sx q[3];
rz(-1.2231474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59584004) q[0];
sx q[0];
rz(-1.4056453) q[0];
sx q[0];
rz(-0.76675057) q[0];
rz(1.2402361) q[1];
sx q[1];
rz(-2.2940472) q[1];
sx q[1];
rz(-0.3516745) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7809217) q[0];
sx q[0];
rz(-0.70510266) q[0];
sx q[0];
rz(-0.4522001) q[0];
rz(-pi) q[1];
rz(2.2064477) q[2];
sx q[2];
rz(-1.1543373) q[2];
sx q[2];
rz(-2.2098429) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1406447) q[1];
sx q[1];
rz(-1.5854156) q[1];
sx q[1];
rz(-2.7468202) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0293803) q[3];
sx q[3];
rz(-0.57313985) q[3];
sx q[3];
rz(-2.074632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9436283) q[2];
sx q[2];
rz(-1.6683656) q[2];
sx q[2];
rz(0.53331214) q[2];
rz(-2.7126281) q[3];
sx q[3];
rz(-0.52754378) q[3];
sx q[3];
rz(2.0146446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0241942) q[0];
sx q[0];
rz(-2.0485931) q[0];
sx q[0];
rz(2.5999516) q[0];
rz(-2.533124) q[1];
sx q[1];
rz(-2.922373) q[1];
sx q[1];
rz(2.0419962) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2521508) q[0];
sx q[0];
rz(-2.5398206) q[0];
sx q[0];
rz(2.2579231) q[0];
rz(2.9526688) q[2];
sx q[2];
rz(-2.0940229) q[2];
sx q[2];
rz(2.1649233) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3179143) q[1];
sx q[1];
rz(-0.70033973) q[1];
sx q[1];
rz(-0.27842303) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8058365) q[3];
sx q[3];
rz(-0.81695518) q[3];
sx q[3];
rz(0.18248617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.6301443) q[2];
sx q[2];
rz(-1.5254598) q[2];
sx q[2];
rz(-2.8586094) q[2];
rz(2.4354368) q[3];
sx q[3];
rz(-1.599584) q[3];
sx q[3];
rz(-2.506536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.054984897) q[0];
sx q[0];
rz(-0.98751706) q[0];
sx q[0];
rz(-1.5299861) q[0];
rz(-2.4967172) q[1];
sx q[1];
rz(-2.6480643) q[1];
sx q[1];
rz(2.9842916) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10282117) q[0];
sx q[0];
rz(-2.4024995) q[0];
sx q[0];
rz(3.0022013) q[0];
rz(2.3072725) q[2];
sx q[2];
rz(-1.7322707) q[2];
sx q[2];
rz(-1.8404567) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.09307043) q[1];
sx q[1];
rz(-1.7659566) q[1];
sx q[1];
rz(0.94961571) q[1];
rz(0.62742426) q[3];
sx q[3];
rz(-1.809123) q[3];
sx q[3];
rz(1.081092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1823696) q[2];
sx q[2];
rz(-2.5490641) q[2];
sx q[2];
rz(-0.83089337) q[2];
rz(-3.0977141) q[3];
sx q[3];
rz(-1.2336122) q[3];
sx q[3];
rz(0.20251814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6458994) q[0];
sx q[0];
rz(-2.2024787) q[0];
sx q[0];
rz(0.010852531) q[0];
rz(-0.27663484) q[1];
sx q[1];
rz(-0.77181548) q[1];
sx q[1];
rz(0.29327926) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5921708) q[0];
sx q[0];
rz(-2.5991419) q[0];
sx q[0];
rz(-3.1300934) q[0];
rz(0.16565928) q[2];
sx q[2];
rz(-0.93369166) q[2];
sx q[2];
rz(0.13656244) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5756702) q[1];
sx q[1];
rz(-0.88469632) q[1];
sx q[1];
rz(-2.5887606) q[1];
x q[2];
rz(0.46383143) q[3];
sx q[3];
rz(-1.2004735) q[3];
sx q[3];
rz(2.2390389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.04348065) q[2];
sx q[2];
rz(-3.0467693) q[2];
sx q[2];
rz(3.0528255) q[2];
rz(-2.8914715) q[3];
sx q[3];
rz(-2.0177896) q[3];
sx q[3];
rz(-2.5626101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0042689) q[0];
sx q[0];
rz(-0.22432888) q[0];
sx q[0];
rz(-1.0639169) q[0];
rz(-0.44257277) q[1];
sx q[1];
rz(-1.7174218) q[1];
sx q[1];
rz(-1.3508505) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3063072) q[0];
sx q[0];
rz(-1.6500236) q[0];
sx q[0];
rz(-1.4863192) q[0];
x q[1];
rz(-0.66019085) q[2];
sx q[2];
rz(-1.635951) q[2];
sx q[2];
rz(2.8934663) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0296214) q[1];
sx q[1];
rz(-2.0051149) q[1];
sx q[1];
rz(-0.017945826) q[1];
x q[2];
rz(1.3488995) q[3];
sx q[3];
rz(-1.4699234) q[3];
sx q[3];
rz(-1.5125121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0315447) q[2];
sx q[2];
rz(-1.8507379) q[2];
sx q[2];
rz(0.44719493) q[2];
rz(-1.7556919) q[3];
sx q[3];
rz(-2.8409676) q[3];
sx q[3];
rz(-0.51469222) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31496012) q[0];
sx q[0];
rz(-1.5216014) q[0];
sx q[0];
rz(-0.78053027) q[0];
rz(2.7138117) q[1];
sx q[1];
rz(-0.3573187) q[1];
sx q[1];
rz(-2.9719877) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1113838) q[0];
sx q[0];
rz(-1.2134524) q[0];
sx q[0];
rz(-1.0650728) q[0];
x q[1];
rz(0.054762997) q[2];
sx q[2];
rz(-0.4224531) q[2];
sx q[2];
rz(-1.7516608) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6360213) q[1];
sx q[1];
rz(-2.1850056) q[1];
sx q[1];
rz(0.057227055) q[1];
x q[2];
rz(2.9986936) q[3];
sx q[3];
rz(-1.8892965) q[3];
sx q[3];
rz(0.83574502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2910989) q[2];
sx q[2];
rz(-1.1824181) q[2];
sx q[2];
rz(-0.69236857) q[2];
rz(-0.46323562) q[3];
sx q[3];
rz(-1.4866456) q[3];
sx q[3];
rz(1.108981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2904084) q[0];
sx q[0];
rz(-1.6155227) q[0];
sx q[0];
rz(0.68646705) q[0];
rz(-0.51207536) q[1];
sx q[1];
rz(-2.8072186) q[1];
sx q[1];
rz(1.0288815) q[1];
rz(-1.2693263) q[2];
sx q[2];
rz(-1.3133247) q[2];
sx q[2];
rz(1.7016344) q[2];
rz(-1.061822) q[3];
sx q[3];
rz(-2.90425) q[3];
sx q[3];
rz(-1.4958285) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
