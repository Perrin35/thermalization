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
rz(-1.8622458) q[0];
sx q[0];
rz(-0.6435464) q[0];
rz(0.14856385) q[1];
sx q[1];
rz(3.7879877) q[1];
sx q[1];
rz(11.990379) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9849562) q[0];
sx q[0];
rz(-0.63581563) q[0];
sx q[0];
rz(1.030907) q[0];
x q[1];
rz(-1.2744489) q[2];
sx q[2];
rz(-1.7360032) q[2];
sx q[2];
rz(-2.8458386) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9427355) q[1];
sx q[1];
rz(-0.1035957) q[1];
sx q[1];
rz(-1.8583244) q[1];
rz(1.3524782) q[3];
sx q[3];
rz(-2.3872833) q[3];
sx q[3];
rz(-1.393817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6286455) q[2];
sx q[2];
rz(-0.60428667) q[2];
sx q[2];
rz(-0.71913546) q[2];
rz(2.5288845) q[3];
sx q[3];
rz(-0.49379525) q[3];
sx q[3];
rz(1.6814211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
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
rz(-2.5510229) q[1];
sx q[1];
rz(-1.5547724) q[1];
sx q[1];
rz(-0.86033598) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96073396) q[0];
sx q[0];
rz(-1.7525808) q[0];
sx q[0];
rz(-2.2341827) q[0];
x q[1];
rz(0.55540076) q[2];
sx q[2];
rz(-1.7199708) q[2];
sx q[2];
rz(-0.55567105) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.115827) q[1];
sx q[1];
rz(-1.2526667) q[1];
sx q[1];
rz(-0.33237388) q[1];
rz(-3.0290696) q[3];
sx q[3];
rz(-0.72970245) q[3];
sx q[3];
rz(-0.94151783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7343973) q[2];
sx q[2];
rz(-1.1043786) q[2];
sx q[2];
rz(2.8779136) q[2];
rz(-0.10144357) q[3];
sx q[3];
rz(-1.0976617) q[3];
sx q[3];
rz(-1.6775848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2270671) q[0];
sx q[0];
rz(-0.83221808) q[0];
sx q[0];
rz(1.9050003) q[0];
rz(0.49691686) q[1];
sx q[1];
rz(-0.91941112) q[1];
sx q[1];
rz(2.9317754) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38302177) q[0];
sx q[0];
rz(-1.5514217) q[0];
sx q[0];
rz(3.0109177) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0210002) q[2];
sx q[2];
rz(-1.6944024) q[2];
sx q[2];
rz(-0.093401366) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1637437) q[1];
sx q[1];
rz(-1.6075378) q[1];
sx q[1];
rz(-1.0669462) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7439458) q[3];
sx q[3];
rz(-0.741091) q[3];
sx q[3];
rz(-0.64797621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0158374) q[2];
sx q[2];
rz(-2.0070984) q[2];
sx q[2];
rz(-2.8718359) q[2];
rz(-2.6848327) q[3];
sx q[3];
rz(-0.41958198) q[3];
sx q[3];
rz(2.8381798) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1064827) q[0];
sx q[0];
rz(-2.7484317) q[0];
sx q[0];
rz(3.0495354) q[0];
rz(0.58097845) q[1];
sx q[1];
rz(-0.61994225) q[1];
sx q[1];
rz(0.43749014) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62037435) q[0];
sx q[0];
rz(-2.552444) q[0];
sx q[0];
rz(1.8626088) q[0];
rz(-1.0863414) q[2];
sx q[2];
rz(-0.24366442) q[2];
sx q[2];
rz(0.59020868) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.60595861) q[1];
sx q[1];
rz(-1.3203964) q[1];
sx q[1];
rz(-1.1634924) q[1];
rz(-pi) q[2];
rz(-1.1933255) q[3];
sx q[3];
rz(-2.4015292) q[3];
sx q[3];
rz(-0.40501324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.94365326) q[2];
sx q[2];
rz(-1.3760171) q[2];
sx q[2];
rz(0.1680689) q[2];
rz(0.69593143) q[3];
sx q[3];
rz(-1.9441425) q[3];
sx q[3];
rz(1.7456938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54774061) q[0];
sx q[0];
rz(-0.20741367) q[0];
sx q[0];
rz(-0.70163027) q[0];
rz(-0.54948366) q[1];
sx q[1];
rz(-2.0618849) q[1];
sx q[1];
rz(0.93542498) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12865484) q[0];
sx q[0];
rz(-1.8194981) q[0];
sx q[0];
rz(-1.1292581) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5176501) q[2];
sx q[2];
rz(-1.4477483) q[2];
sx q[2];
rz(1.4966663) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4051263) q[1];
sx q[1];
rz(-1.163981) q[1];
sx q[1];
rz(-0.94405014) q[1];
rz(-0.071604972) q[3];
sx q[3];
rz(-1.6354899) q[3];
sx q[3];
rz(-2.477579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.51297411) q[2];
sx q[2];
rz(-1.3209891) q[2];
sx q[2];
rz(-2.4690907) q[2];
rz(-2.3682112) q[3];
sx q[3];
rz(-1.0459817) q[3];
sx q[3];
rz(2.8092303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8124939) q[0];
sx q[0];
rz(-2.2342873) q[0];
sx q[0];
rz(1.4688274) q[0];
rz(0.86917296) q[1];
sx q[1];
rz(-2.4425127) q[1];
sx q[1];
rz(-2.8156143) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9470254) q[0];
sx q[0];
rz(-1.86867) q[0];
sx q[0];
rz(-1.711571) q[0];
x q[1];
rz(-2.5766854) q[2];
sx q[2];
rz(-2.0015026) q[2];
sx q[2];
rz(-0.6503833) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1010998) q[1];
sx q[1];
rz(-1.0940576) q[1];
sx q[1];
rz(0.65495305) q[1];
rz(-pi) q[2];
x q[2];
rz(0.066727344) q[3];
sx q[3];
rz(-1.9983833) q[3];
sx q[3];
rz(-0.59837435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0526696) q[2];
sx q[2];
rz(-1.5648877) q[2];
sx q[2];
rz(2.3443473) q[2];
rz(0.097660216) q[3];
sx q[3];
rz(-0.74334136) q[3];
sx q[3];
rz(1.9184448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2753817) q[0];
sx q[0];
rz(-2.1480063) q[0];
sx q[0];
rz(-2.1589101) q[0];
rz(1.6536225) q[1];
sx q[1];
rz(-2.5762312) q[1];
sx q[1];
rz(2.8178094) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57162962) q[0];
sx q[0];
rz(-1.6545719) q[0];
sx q[0];
rz(-0.28599583) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.95877842) q[2];
sx q[2];
rz(-1.806815) q[2];
sx q[2];
rz(-0.70597202) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1947136) q[1];
sx q[1];
rz(-0.66821287) q[1];
sx q[1];
rz(2.9577423) q[1];
rz(-2.6402506) q[3];
sx q[3];
rz(-1.2698113) q[3];
sx q[3];
rz(-1.6607064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.71161485) q[2];
sx q[2];
rz(-2.6275676) q[2];
sx q[2];
rz(3.1058822) q[2];
rz(-2.4014373) q[3];
sx q[3];
rz(-1.454708) q[3];
sx q[3];
rz(1.9945701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3586054) q[0];
sx q[0];
rz(-0.37707034) q[0];
sx q[0];
rz(-0.20088917) q[0];
rz(-2.5783077) q[1];
sx q[1];
rz(-1.3689684) q[1];
sx q[1];
rz(0.39931452) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7359283) q[0];
sx q[0];
rz(-2.9333994) q[0];
sx q[0];
rz(0.8008147) q[0];
rz(-pi) q[1];
rz(-1.6496193) q[2];
sx q[2];
rz(-0.48187253) q[2];
sx q[2];
rz(-2.3367867) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.47013624) q[1];
sx q[1];
rz(-0.54959198) q[1];
sx q[1];
rz(-2.2983645) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7186728) q[3];
sx q[3];
rz(-2.3302571) q[3];
sx q[3];
rz(2.3911311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1839972) q[2];
sx q[2];
rz(-2.5338569) q[2];
sx q[2];
rz(0.11873928) q[2];
rz(-0.27058288) q[3];
sx q[3];
rz(-0.99004531) q[3];
sx q[3];
rz(-0.15429193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1110693) q[0];
sx q[0];
rz(-1.9707762) q[0];
sx q[0];
rz(-0.4554553) q[0];
rz(-2.9083374) q[1];
sx q[1];
rz(-0.74359727) q[1];
sx q[1];
rz(0.0049678405) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.395803) q[0];
sx q[0];
rz(-3.0533049) q[0];
sx q[0];
rz(-2.7448065) q[0];
x q[1];
rz(-1.2521148) q[2];
sx q[2];
rz(-1.2010842) q[2];
sx q[2];
rz(1.2546509) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.50434166) q[1];
sx q[1];
rz(-2.0310861) q[1];
sx q[1];
rz(-1.00072) q[1];
rz(-pi) q[2];
rz(-0.016214215) q[3];
sx q[3];
rz(-1.9965648) q[3];
sx q[3];
rz(-0.92416422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7838955) q[2];
sx q[2];
rz(-1.8507439) q[2];
sx q[2];
rz(-0.89007393) q[2];
rz(2.2642073) q[3];
sx q[3];
rz(-0.93534094) q[3];
sx q[3];
rz(-2.4963511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18168618) q[0];
sx q[0];
rz(-3.063995) q[0];
sx q[0];
rz(-2.087387) q[0];
rz(-0.28610817) q[1];
sx q[1];
rz(-1.0461297) q[1];
sx q[1];
rz(-1.8792763) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12114502) q[0];
sx q[0];
rz(-2.1845) q[0];
sx q[0];
rz(1.2200939) q[0];
x q[1];
rz(-1.3081864) q[2];
sx q[2];
rz(-1.3725026) q[2];
sx q[2];
rz(-1.8526371) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8827829) q[1];
sx q[1];
rz(-1.38816) q[1];
sx q[1];
rz(0.44328494) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1256874) q[3];
sx q[3];
rz(-1.6421179) q[3];
sx q[3];
rz(-1.9147022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.012764843) q[2];
sx q[2];
rz(-2.3598857) q[2];
sx q[2];
rz(0.26415602) q[2];
rz(2.9014897) q[3];
sx q[3];
rz(-2.0930347) q[3];
sx q[3];
rz(-0.31606328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(0.24513292) q[0];
sx q[0];
rz(-0.75440732) q[0];
sx q[0];
rz(-0.9040133) q[0];
rz(-2.8527507) q[1];
sx q[1];
rz(-1.7542417) q[1];
sx q[1];
rz(-0.075686553) q[1];
rz(-2.8223128) q[2];
sx q[2];
rz(-1.8055331) q[2];
sx q[2];
rz(-1.8100678) q[2];
rz(-3.0700923) q[3];
sx q[3];
rz(-1.6115474) q[3];
sx q[3];
rz(-1.1781884) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
