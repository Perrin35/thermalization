OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6501939) q[0];
sx q[0];
rz(-2.8770652) q[0];
sx q[0];
rz(-2.7471623) q[0];
rz(-3.1354304) q[1];
sx q[1];
rz(-2.8013464) q[1];
sx q[1];
rz(-1.9415829) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.278468) q[0];
sx q[0];
rz(-1.5764563) q[0];
sx q[0];
rz(-1.6186884) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0179022) q[2];
sx q[2];
rz(-1.2667709) q[2];
sx q[2];
rz(2.6927039) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8438699) q[1];
sx q[1];
rz(-0.6134609) q[1];
sx q[1];
rz(-0.86617275) q[1];
rz(-1.6148189) q[3];
sx q[3];
rz(-1.9485954) q[3];
sx q[3];
rz(1.6954741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8346943) q[2];
sx q[2];
rz(-1.4935741) q[2];
sx q[2];
rz(-2.8519894) q[2];
rz(-2.2662207) q[3];
sx q[3];
rz(-2.1353728) q[3];
sx q[3];
rz(-0.059710596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61525476) q[0];
sx q[0];
rz(-0.86741388) q[0];
sx q[0];
rz(2.7976024) q[0];
rz(3.0572609) q[1];
sx q[1];
rz(-2.4721959) q[1];
sx q[1];
rz(-1.7864236) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7536613) q[0];
sx q[0];
rz(-2.5270215) q[0];
sx q[0];
rz(-1.4580926) q[0];
x q[1];
rz(-0.16263527) q[2];
sx q[2];
rz(-1.6911104) q[2];
sx q[2];
rz(1.4974809) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4210216) q[1];
sx q[1];
rz(-0.45383006) q[1];
sx q[1];
rz(-1.4245016) q[1];
rz(-pi) q[2];
rz(1.9676898) q[3];
sx q[3];
rz(-0.19860425) q[3];
sx q[3];
rz(-1.393115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.29558674) q[2];
sx q[2];
rz(-2.3082374) q[2];
sx q[2];
rz(-2.6039092) q[2];
rz(2.6387571) q[3];
sx q[3];
rz(-0.42481315) q[3];
sx q[3];
rz(-2.1311549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4780592) q[0];
sx q[0];
rz(-0.53776598) q[0];
sx q[0];
rz(-0.64087254) q[0];
rz(2.3928941) q[1];
sx q[1];
rz(-2.1332032) q[1];
sx q[1];
rz(-2.0764988) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0707866) q[0];
sx q[0];
rz(-0.35937989) q[0];
sx q[0];
rz(-1.5887567) q[0];
x q[1];
rz(0.89945729) q[2];
sx q[2];
rz(-2.3231635) q[2];
sx q[2];
rz(-0.98086548) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.47485891) q[1];
sx q[1];
rz(-1.534728) q[1];
sx q[1];
rz(0.95400793) q[1];
rz(-pi) q[2];
rz(0.84824003) q[3];
sx q[3];
rz(-1.3812997) q[3];
sx q[3];
rz(-3.0998067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.76434) q[2];
sx q[2];
rz(-1.179402) q[2];
sx q[2];
rz(0.71737814) q[2];
rz(0.68850368) q[3];
sx q[3];
rz(-2.513956) q[3];
sx q[3];
rz(0.29754105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3142969) q[0];
sx q[0];
rz(-1.941444) q[0];
sx q[0];
rz(0.24969077) q[0];
rz(-2.1266134) q[1];
sx q[1];
rz(-2.8453638) q[1];
sx q[1];
rz(-3.1304741) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7751559) q[0];
sx q[0];
rz(-1.1661134) q[0];
sx q[0];
rz(-2.0648271) q[0];
rz(-pi) q[1];
rz(-2.501802) q[2];
sx q[2];
rz(-2.1961229) q[2];
sx q[2];
rz(-1.7196136) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.74777714) q[1];
sx q[1];
rz(-1.2443466) q[1];
sx q[1];
rz(-2.087489) q[1];
x q[2];
rz(2.8233077) q[3];
sx q[3];
rz(-1.7537698) q[3];
sx q[3];
rz(-2.0103679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.49542385) q[2];
sx q[2];
rz(-1.255722) q[2];
sx q[2];
rz(-2.8584976) q[2];
rz(-0.66343534) q[3];
sx q[3];
rz(-0.59316558) q[3];
sx q[3];
rz(2.2535113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11113142) q[0];
sx q[0];
rz(-0.24065329) q[0];
sx q[0];
rz(0.33183137) q[0];
rz(-0.49452531) q[1];
sx q[1];
rz(-1.8437513) q[1];
sx q[1];
rz(1.3269075) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3834284) q[0];
sx q[0];
rz(-2.313662) q[0];
sx q[0];
rz(1.3616256) q[0];
x q[1];
rz(-2.5577776) q[2];
sx q[2];
rz(-2.2540255) q[2];
sx q[2];
rz(-0.77606397) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8618968) q[1];
sx q[1];
rz(-1.8021291) q[1];
sx q[1];
rz(-0.33613236) q[1];
rz(-2.1760686) q[3];
sx q[3];
rz(-11*pi/13) q[3];
sx q[3];
rz(1.0799288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.999324) q[2];
sx q[2];
rz(-0.48196718) q[2];
sx q[2];
rz(-1.8959321) q[2];
rz(1.2549531) q[3];
sx q[3];
rz(-0.84147036) q[3];
sx q[3];
rz(-2.3033223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6923043) q[0];
sx q[0];
rz(-0.0031539991) q[0];
sx q[0];
rz(0.6814878) q[0];
rz(-2.9340414) q[1];
sx q[1];
rz(-0.46336585) q[1];
sx q[1];
rz(1.1157657) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.241248) q[0];
sx q[0];
rz(-1.9755409) q[0];
sx q[0];
rz(0.83165283) q[0];
rz(-pi) q[1];
x q[1];
rz(0.90066465) q[2];
sx q[2];
rz(-1.9158944) q[2];
sx q[2];
rz(-0.64054856) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7129732) q[1];
sx q[1];
rz(-1.5884807) q[1];
sx q[1];
rz(-2.5531205) q[1];
rz(-pi) q[2];
rz(1.3917771) q[3];
sx q[3];
rz(-0.90913032) q[3];
sx q[3];
rz(-0.20565198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.21268022) q[2];
sx q[2];
rz(-1.8447515) q[2];
sx q[2];
rz(2.7872655) q[2];
rz(-2.8220693) q[3];
sx q[3];
rz(-1.1525681) q[3];
sx q[3];
rz(2.7697146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1324683) q[0];
sx q[0];
rz(-2.9102944) q[0];
sx q[0];
rz(2.4672467) q[0];
rz(1.1122423) q[1];
sx q[1];
rz(-0.66450417) q[1];
sx q[1];
rz(2.5792714) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6816662) q[0];
sx q[0];
rz(-2.473218) q[0];
sx q[0];
rz(1.5468803) q[0];
rz(-pi) q[1];
rz(-0.33195915) q[2];
sx q[2];
rz(-1.4753642) q[2];
sx q[2];
rz(-1.3704259) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.79531419) q[1];
sx q[1];
rz(-0.242713) q[1];
sx q[1];
rz(-0.64530356) q[1];
rz(-pi) q[2];
rz(3.0543442) q[3];
sx q[3];
rz(-2.1656519) q[3];
sx q[3];
rz(-1.113254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0397296) q[2];
sx q[2];
rz(-0.9938643) q[2];
sx q[2];
rz(-0.34004655) q[2];
rz(0.21751054) q[3];
sx q[3];
rz(-2.1931931) q[3];
sx q[3];
rz(-2.8229009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6927004) q[0];
sx q[0];
rz(-0.046148766) q[0];
sx q[0];
rz(0.39644077) q[0];
rz(0.13892826) q[1];
sx q[1];
rz(-2.6815092) q[1];
sx q[1];
rz(1.6202392) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6103219) q[0];
sx q[0];
rz(-1.5025508) q[0];
sx q[0];
rz(0.1971498) q[0];
rz(2.7395758) q[2];
sx q[2];
rz(-2.0094299) q[2];
sx q[2];
rz(1.5356262) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8926881) q[1];
sx q[1];
rz(-2.3612594) q[1];
sx q[1];
rz(0.08463879) q[1];
rz(-pi) q[2];
rz(-1.7117386) q[3];
sx q[3];
rz(-0.68416506) q[3];
sx q[3];
rz(-1.9977026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5788995) q[2];
sx q[2];
rz(-1.0572628) q[2];
sx q[2];
rz(2.8472624) q[2];
rz(-2.0108022) q[3];
sx q[3];
rz(-1.7629938) q[3];
sx q[3];
rz(0.99564266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6482553) q[0];
sx q[0];
rz(-0.91518891) q[0];
sx q[0];
rz(2.7822568) q[0];
rz(2.1954779) q[1];
sx q[1];
rz(-0.39603907) q[1];
sx q[1];
rz(-2.8709581) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7662738) q[0];
sx q[0];
rz(-3.0941071) q[0];
sx q[0];
rz(0.71592285) q[0];
rz(-pi) q[1];
rz(0.46220772) q[2];
sx q[2];
rz(-2.2773829) q[2];
sx q[2];
rz(-0.30009899) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.088085727) q[1];
sx q[1];
rz(-0.063789531) q[1];
sx q[1];
rz(-2.3699058) q[1];
x q[2];
rz(-2.1569096) q[3];
sx q[3];
rz(-0.95041785) q[3];
sx q[3];
rz(1.8970722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3140807) q[2];
sx q[2];
rz(-0.22958799) q[2];
sx q[2];
rz(-0.45544004) q[2];
rz(-2.3296302) q[3];
sx q[3];
rz(-2.2596695) q[3];
sx q[3];
rz(0.5493831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(-2.6311326) q[0];
sx q[0];
rz(-1.5083418) q[0];
sx q[0];
rz(2.4023138) q[0];
rz(2.9108858) q[1];
sx q[1];
rz(-2.6735327) q[1];
sx q[1];
rz(0.49490067) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4821645) q[0];
sx q[0];
rz(-0.48766252) q[0];
sx q[0];
rz(-1.6517261) q[0];
rz(2.8137389) q[2];
sx q[2];
rz(-1.577539) q[2];
sx q[2];
rz(1.1525796) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1315688) q[1];
sx q[1];
rz(-1.6820757) q[1];
sx q[1];
rz(2.2113423) q[1];
rz(-pi) q[2];
rz(0.87402113) q[3];
sx q[3];
rz(-1.6439983) q[3];
sx q[3];
rz(0.30833581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.13359244) q[2];
sx q[2];
rz(-1.091489) q[2];
sx q[2];
rz(-2.8137394) q[2];
rz(-2.949529) q[3];
sx q[3];
rz(-2.8938507) q[3];
sx q[3];
rz(-1.0333992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0192169) q[0];
sx q[0];
rz(-1.6288971) q[0];
sx q[0];
rz(-1.6678641) q[0];
rz(-2.3090251) q[1];
sx q[1];
rz(-1.7201798) q[1];
sx q[1];
rz(-1.7725772) q[1];
rz(-2.1140425) q[2];
sx q[2];
rz(-1.3941358) q[2];
sx q[2];
rz(0.91640581) q[2];
rz(1.3151863) q[3];
sx q[3];
rz(-1.8288463) q[3];
sx q[3];
rz(2.0127206) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
