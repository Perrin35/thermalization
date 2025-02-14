OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.65741444) q[0];
sx q[0];
rz(1.2793469) q[0];
sx q[0];
rz(10.068324) q[0];
rz(0.14856385) q[1];
sx q[1];
rz(-2.4951976) q[1];
sx q[1];
rz(2.5656011) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86344415) q[0];
sx q[0];
rz(-1.2605901) q[0];
sx q[0];
rz(-1.006406) q[0];
rz(-0.17259105) q[2];
sx q[2];
rz(-1.8629889) q[2];
sx q[2];
rz(-1.2248696) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6580087) q[1];
sx q[1];
rz(-1.600126) q[1];
sx q[1];
rz(1.4714249) q[1];
rz(-pi) q[2];
rz(-0.20078076) q[3];
sx q[3];
rz(-0.83856486) q[3];
sx q[3];
rz(1.4522566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6286455) q[2];
sx q[2];
rz(-0.60428667) q[2];
sx q[2];
rz(2.4224572) q[2];
rz(2.5288845) q[3];
sx q[3];
rz(-2.6477974) q[3];
sx q[3];
rz(-1.6814211) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7270108) q[0];
sx q[0];
rz(-2.5717323) q[0];
sx q[0];
rz(1.3563096) q[0];
rz(-0.59056979) q[1];
sx q[1];
rz(-1.5547724) q[1];
sx q[1];
rz(-2.2812567) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38274327) q[0];
sx q[0];
rz(-2.4573987) q[0];
sx q[0];
rz(-1.280715) q[0];
rz(-pi) q[1];
rz(-1.3957295) q[2];
sx q[2];
rz(-1.0222728) q[2];
sx q[2];
rz(1.107094) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2812432) q[1];
sx q[1];
rz(-0.45595995) q[1];
sx q[1];
rz(-0.79001536) q[1];
rz(1.4707056) q[3];
sx q[3];
rz(-2.2948569) q[3];
sx q[3];
rz(-2.049618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7343973) q[2];
sx q[2];
rz(-2.037214) q[2];
sx q[2];
rz(-0.26367903) q[2];
rz(0.10144357) q[3];
sx q[3];
rz(-2.0439309) q[3];
sx q[3];
rz(-1.6775848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9145255) q[0];
sx q[0];
rz(-2.3093746) q[0];
sx q[0];
rz(-1.2365923) q[0];
rz(0.49691686) q[1];
sx q[1];
rz(-0.91941112) q[1];
sx q[1];
rz(-0.20981728) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.041417) q[0];
sx q[0];
rz(-3.0094973) q[0];
sx q[0];
rz(-2.9939674) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0210002) q[2];
sx q[2];
rz(-1.4471902) q[2];
sx q[2];
rz(0.093401366) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4821149) q[1];
sx q[1];
rz(-2.6365197) q[1];
sx q[1];
rz(-1.6467846) q[1];
x q[2];
rz(-2.9852226) q[3];
sx q[3];
rz(-2.2982882) q[3];
sx q[3];
rz(-2.7263977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0158374) q[2];
sx q[2];
rz(-2.0070984) q[2];
sx q[2];
rz(-2.8718359) q[2];
rz(-0.45675993) q[3];
sx q[3];
rz(-0.41958198) q[3];
sx q[3];
rz(0.30341283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.03511) q[0];
sx q[0];
rz(-2.7484317) q[0];
sx q[0];
rz(0.092057236) q[0];
rz(-2.5606142) q[1];
sx q[1];
rz(-2.5216504) q[1];
sx q[1];
rz(-0.43749014) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5212183) q[0];
sx q[0];
rz(-0.58914864) q[0];
sx q[0];
rz(-1.8626088) q[0];
rz(-pi) q[1];
rz(1.0863414) q[2];
sx q[2];
rz(-0.24366442) q[2];
sx q[2];
rz(-0.59020868) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.60595861) q[1];
sx q[1];
rz(-1.8211963) q[1];
sx q[1];
rz(1.1634924) q[1];
rz(-pi) q[2];
rz(1.9482671) q[3];
sx q[3];
rz(-0.74006343) q[3];
sx q[3];
rz(-2.7365794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.94365326) q[2];
sx q[2];
rz(-1.3760171) q[2];
sx q[2];
rz(0.1680689) q[2];
rz(-0.69593143) q[3];
sx q[3];
rz(-1.9441425) q[3];
sx q[3];
rz(1.3958989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54774061) q[0];
sx q[0];
rz(-0.20741367) q[0];
sx q[0];
rz(0.70163027) q[0];
rz(0.54948366) q[1];
sx q[1];
rz(-1.0797078) q[1];
sx q[1];
rz(0.93542498) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12865484) q[0];
sx q[0];
rz(-1.8194981) q[0];
sx q[0];
rz(-2.0123345) q[0];
rz(-pi) q[1];
rz(3.0183724) q[2];
sx q[2];
rz(-1.5180523) q[2];
sx q[2];
rz(-3.0739917) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4758318) q[1];
sx q[1];
rz(-2.40959) q[1];
sx q[1];
rz(0.93722645) q[1];
rz(-pi) q[2];
rz(-0.071604972) q[3];
sx q[3];
rz(-1.5061028) q[3];
sx q[3];
rz(2.477579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.51297411) q[2];
sx q[2];
rz(-1.3209891) q[2];
sx q[2];
rz(-2.4690907) q[2];
rz(2.3682112) q[3];
sx q[3];
rz(-2.0956109) q[3];
sx q[3];
rz(2.8092303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32909876) q[0];
sx q[0];
rz(-0.90730539) q[0];
sx q[0];
rz(1.4688274) q[0];
rz(0.86917296) q[1];
sx q[1];
rz(-2.4425127) q[1];
sx q[1];
rz(-2.8156143) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7237968) q[0];
sx q[0];
rz(-1.4362596) q[0];
sx q[0];
rz(-0.30067443) q[0];
x q[1];
rz(-13/(2*pi)) q[2];
sx q[2];
rz(-2.0788134) q[2];
sx q[2];
rz(-1.9624868) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.068858154) q[1];
sx q[1];
rz(-0.78887296) q[1];
sx q[1];
rz(-2.4383208) q[1];
x q[2];
rz(1.7160837) q[3];
sx q[3];
rz(-0.43244468) q[3];
sx q[3];
rz(-0.43859453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.088923067) q[2];
sx q[2];
rz(-1.5648877) q[2];
sx q[2];
rz(-0.79724533) q[2];
rz(0.097660216) q[3];
sx q[3];
rz(-2.3982513) q[3];
sx q[3];
rz(1.2231479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2753817) q[0];
sx q[0];
rz(-0.99358639) q[0];
sx q[0];
rz(-0.98268253) q[0];
rz(-1.4879701) q[1];
sx q[1];
rz(-2.5762312) q[1];
sx q[1];
rz(2.8178094) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97456565) q[0];
sx q[0];
rz(-1.2858316) q[0];
sx q[0];
rz(1.6581013) q[0];
rz(-pi) q[1];
rz(-2.8558015) q[2];
sx q[2];
rz(-2.1634844) q[2];
sx q[2];
rz(2.1140849) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6205755) q[1];
sx q[1];
rz(-1.4572826) q[1];
sx q[1];
rz(-0.6599636) q[1];
rz(-pi) q[2];
rz(-1.2305831) q[3];
sx q[3];
rz(-1.0939301) q[3];
sx q[3];
rz(0.25097706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.71161485) q[2];
sx q[2];
rz(-2.6275676) q[2];
sx q[2];
rz(-0.035710486) q[2];
rz(2.4014373) q[3];
sx q[3];
rz(-1.6868846) q[3];
sx q[3];
rz(1.9945701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78298727) q[0];
sx q[0];
rz(-2.7645223) q[0];
sx q[0];
rz(-2.9407035) q[0];
rz(-2.5783077) q[1];
sx q[1];
rz(-1.3689684) q[1];
sx q[1];
rz(-2.7422781) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4056643) q[0];
sx q[0];
rz(-0.20819323) q[0];
sx q[0];
rz(-2.3407779) q[0];
rz(-pi) q[1];
rz(1.0901997) q[2];
sx q[2];
rz(-1.5342964) q[2];
sx q[2];
rz(-0.69611204) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6714564) q[1];
sx q[1];
rz(-2.5920007) q[1];
sx q[1];
rz(2.2983645) q[1];
rz(-pi) q[2];
rz(-0.42291989) q[3];
sx q[3];
rz(-0.81133553) q[3];
sx q[3];
rz(2.3911311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1839972) q[2];
sx q[2];
rz(-0.60773578) q[2];
sx q[2];
rz(0.11873928) q[2];
rz(-0.27058288) q[3];
sx q[3];
rz(-2.1515473) q[3];
sx q[3];
rz(-2.9873007) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0305233) q[0];
sx q[0];
rz(-1.1708165) q[0];
sx q[0];
rz(-0.4554553) q[0];
rz(0.2332553) q[1];
sx q[1];
rz(-2.3979954) q[1];
sx q[1];
rz(-0.0049678405) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7457897) q[0];
sx q[0];
rz(-0.088287778) q[0];
sx q[0];
rz(-2.7448065) q[0];
x q[1];
rz(-1.2521148) q[2];
sx q[2];
rz(-1.9405085) q[2];
sx q[2];
rz(1.8869417) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.50434166) q[1];
sx q[1];
rz(-1.1105065) q[1];
sx q[1];
rz(-2.1408726) q[1];
rz(-0.016214215) q[3];
sx q[3];
rz(-1.1450279) q[3];
sx q[3];
rz(0.92416422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3576971) q[2];
sx q[2];
rz(-1.8507439) q[2];
sx q[2];
rz(2.2515187) q[2];
rz(-2.2642073) q[3];
sx q[3];
rz(-0.93534094) q[3];
sx q[3];
rz(2.4963511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18168618) q[0];
sx q[0];
rz(-3.063995) q[0];
sx q[0];
rz(2.087387) q[0];
rz(2.8554845) q[1];
sx q[1];
rz(-1.0461297) q[1];
sx q[1];
rz(-1.8792763) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4842997) q[0];
sx q[0];
rz(-1.2861006) q[0];
sx q[0];
rz(-2.4980252) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9364519) q[2];
sx q[2];
rz(-1.3134505) q[2];
sx q[2];
rz(2.8068451) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9156277) q[1];
sx q[1];
rz(-2.0061992) q[1];
sx q[1];
rz(-1.3691203) q[1];
rz(-1.7055776) q[3];
sx q[3];
rz(-0.55897903) q[3];
sx q[3];
rz(-2.9121484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1288278) q[2];
sx q[2];
rz(-2.3598857) q[2];
sx q[2];
rz(0.26415602) q[2];
rz(0.24010298) q[3];
sx q[3];
rz(-2.0930347) q[3];
sx q[3];
rz(0.31606328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24513292) q[0];
sx q[0];
rz(-0.75440732) q[0];
sx q[0];
rz(-0.9040133) q[0];
rz(2.8527507) q[1];
sx q[1];
rz(-1.387351) q[1];
sx q[1];
rz(3.0659061) q[1];
rz(-1.8175387) q[2];
sx q[2];
rz(-1.8810234) q[2];
sx q[2];
rz(-0.16253139) q[2];
rz(1.6116517) q[3];
sx q[3];
rz(-1.4993555) q[3];
sx q[3];
rz(-2.7460668) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
