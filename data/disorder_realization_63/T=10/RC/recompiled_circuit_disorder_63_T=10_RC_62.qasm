OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.9665943) q[0];
sx q[0];
rz(-2.7881665) q[0];
sx q[0];
rz(2.0768291) q[0];
rz(-2.3454173) q[1];
sx q[1];
rz(-1.2086955) q[1];
sx q[1];
rz(2.605521) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.85815) q[0];
sx q[0];
rz(-1.3387696) q[0];
sx q[0];
rz(-2.0767077) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4921161) q[2];
sx q[2];
rz(-2.4298551) q[2];
sx q[2];
rz(-2.5770238) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.88692611) q[1];
sx q[1];
rz(-2.1517422) q[1];
sx q[1];
rz(1.3455774) q[1];
x q[2];
rz(0.53341289) q[3];
sx q[3];
rz(-0.80848137) q[3];
sx q[3];
rz(-1.4097139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4094231) q[2];
sx q[2];
rz(-2.7004341) q[2];
sx q[2];
rz(2.1146963) q[2];
rz(-0.25201592) q[3];
sx q[3];
rz(-1.1427294) q[3];
sx q[3];
rz(2.2759329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.058763) q[0];
sx q[0];
rz(-2.7936462) q[0];
sx q[0];
rz(-1.3902364) q[0];
rz(2.305796) q[1];
sx q[1];
rz(-0.73671571) q[1];
sx q[1];
rz(0.70835152) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3268711) q[0];
sx q[0];
rz(-0.46778361) q[0];
sx q[0];
rz(-2.4844556) q[0];
x q[1];
rz(-2.1721241) q[2];
sx q[2];
rz(-1.6919486) q[2];
sx q[2];
rz(-3.0960992) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0185512) q[1];
sx q[1];
rz(-2.0611144) q[1];
sx q[1];
rz(-2.3398188) q[1];
rz(-pi) q[2];
rz(0.80531081) q[3];
sx q[3];
rz(-0.66262965) q[3];
sx q[3];
rz(0.90028541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0248802) q[2];
sx q[2];
rz(-2.34989) q[2];
sx q[2];
rz(2.8498245) q[2];
rz(-3.0388888) q[3];
sx q[3];
rz(-1.7385959) q[3];
sx q[3];
rz(-1.5244938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3290688) q[0];
sx q[0];
rz(-3.0561495) q[0];
sx q[0];
rz(0.30971757) q[0];
rz(1.6614871) q[1];
sx q[1];
rz(-1.8146351) q[1];
sx q[1];
rz(-0.57166878) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11274352) q[0];
sx q[0];
rz(-1.0380121) q[0];
sx q[0];
rz(1.4677731) q[0];
x q[1];
rz(-2.0057726) q[2];
sx q[2];
rz(-1.8990714) q[2];
sx q[2];
rz(-0.091718397) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6211105) q[1];
sx q[1];
rz(-0.74319327) q[1];
sx q[1];
rz(-0.1964257) q[1];
rz(0.27922697) q[3];
sx q[3];
rz(-2.673827) q[3];
sx q[3];
rz(1.2665018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9812575) q[2];
sx q[2];
rz(-0.91032878) q[2];
sx q[2];
rz(0.70880115) q[2];
rz(-2.839084) q[3];
sx q[3];
rz(-1.6995647) q[3];
sx q[3];
rz(-1.1423053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4335094) q[0];
sx q[0];
rz(-2.0286562) q[0];
sx q[0];
rz(-2.3420912) q[0];
rz(-0.049830534) q[1];
sx q[1];
rz(-2.2466876) q[1];
sx q[1];
rz(-2.9611011) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4806992) q[0];
sx q[0];
rz(-1.1850712) q[0];
sx q[0];
rz(1.6688523) q[0];
x q[1];
rz(-1.7454342) q[2];
sx q[2];
rz(-1.9348113) q[2];
sx q[2];
rz(0.82496914) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4304639) q[1];
sx q[1];
rz(-2.8006878) q[1];
sx q[1];
rz(-2.7845963) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0666215) q[3];
sx q[3];
rz(-1.5600481) q[3];
sx q[3];
rz(1.4533952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0522456) q[2];
sx q[2];
rz(-1.1648488) q[2];
sx q[2];
rz(1.5578516) q[2];
rz(-1.7662988) q[3];
sx q[3];
rz(-1.3170653) q[3];
sx q[3];
rz(-2.2089675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76064008) q[0];
sx q[0];
rz(-2.8044658) q[0];
sx q[0];
rz(-0.98651648) q[0];
rz(-1.1622693) q[1];
sx q[1];
rz(-1.2170075) q[1];
sx q[1];
rz(-0.25156897) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7895296) q[0];
sx q[0];
rz(-1.7460199) q[0];
sx q[0];
rz(-0.87900889) q[0];
rz(0.17390522) q[2];
sx q[2];
rz(-0.46602962) q[2];
sx q[2];
rz(-0.58338651) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.561589) q[1];
sx q[1];
rz(-1.5981234) q[1];
sx q[1];
rz(0.1954397) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.86815636) q[3];
sx q[3];
rz(-1.9725102) q[3];
sx q[3];
rz(0.86740869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.95785561) q[2];
sx q[2];
rz(-0.541406) q[2];
sx q[2];
rz(1.0305369) q[2];
rz(1.41097) q[3];
sx q[3];
rz(-2.1822699) q[3];
sx q[3];
rz(2.5103536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75491607) q[0];
sx q[0];
rz(-2.6658391) q[0];
sx q[0];
rz(-0.85012287) q[0];
rz(1.1823581) q[1];
sx q[1];
rz(-1.2344924) q[1];
sx q[1];
rz(2.7485671) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.128199) q[0];
sx q[0];
rz(-1.734874) q[0];
sx q[0];
rz(1.8830995) q[0];
rz(-pi) q[1];
rz(-2.9057301) q[2];
sx q[2];
rz(-0.84025192) q[2];
sx q[2];
rz(2.0116531) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1054354) q[1];
sx q[1];
rz(-1.5762377) q[1];
sx q[1];
rz(-0.91962199) q[1];
rz(1.2971446) q[3];
sx q[3];
rz(-1.8575462) q[3];
sx q[3];
rz(-1.2424038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7810016) q[2];
sx q[2];
rz(-2.0157308) q[2];
sx q[2];
rz(0.97314107) q[2];
rz(-0.94758236) q[3];
sx q[3];
rz(-2.0411453) q[3];
sx q[3];
rz(-1.1943641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4457552) q[0];
sx q[0];
rz(-0.42036244) q[0];
sx q[0];
rz(1.5234891) q[0];
rz(-0.37480005) q[1];
sx q[1];
rz(-1.8811036) q[1];
sx q[1];
rz(-3.135625) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80124827) q[0];
sx q[0];
rz(-0.72754117) q[0];
sx q[0];
rz(1.1447385) q[0];
x q[1];
rz(-2.0231831) q[2];
sx q[2];
rz(-0.87086073) q[2];
sx q[2];
rz(1.6391022) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.70049268) q[1];
sx q[1];
rz(-1.0171434) q[1];
sx q[1];
rz(3.0667552) q[1];
rz(-pi) q[2];
rz(-3.0217516) q[3];
sx q[3];
rz(-0.37444515) q[3];
sx q[3];
rz(3.0344809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5696047) q[2];
sx q[2];
rz(-2.5964952) q[2];
sx q[2];
rz(-2.9515284) q[2];
rz(2.7251785) q[3];
sx q[3];
rz(-0.9581241) q[3];
sx q[3];
rz(-0.73474187) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1104601) q[0];
sx q[0];
rz(-2.0955595) q[0];
sx q[0];
rz(0.53623143) q[0];
rz(-0.93332943) q[1];
sx q[1];
rz(-1.5737165) q[1];
sx q[1];
rz(2.1933864) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3535239) q[0];
sx q[0];
rz(-2.3234962) q[0];
sx q[0];
rz(-1.3226932) q[0];
rz(2.2657329) q[2];
sx q[2];
rz(-0.63820733) q[2];
sx q[2];
rz(-1.9583595) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1369689) q[1];
sx q[1];
rz(-1.7925295) q[1];
sx q[1];
rz(-1.5066513) q[1];
rz(-0.20235297) q[3];
sx q[3];
rz(-0.36365299) q[3];
sx q[3];
rz(2.3285151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5006717) q[2];
sx q[2];
rz(-1.5189974) q[2];
sx q[2];
rz(0.24027696) q[2];
rz(2.7622973) q[3];
sx q[3];
rz(-1.0255739) q[3];
sx q[3];
rz(1.7933638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.76294476) q[0];
sx q[0];
rz(-0.33319107) q[0];
sx q[0];
rz(-2.8906524) q[0];
rz(0.25302408) q[1];
sx q[1];
rz(-1.3809985) q[1];
sx q[1];
rz(0.35266638) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34516016) q[0];
sx q[0];
rz(-0.16819084) q[0];
sx q[0];
rz(-2.3229984) q[0];
x q[1];
rz(-0.30013957) q[2];
sx q[2];
rz(-1.7087414) q[2];
sx q[2];
rz(2.2696242) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7749412) q[1];
sx q[1];
rz(-1.0171486) q[1];
sx q[1];
rz(1.304438) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0274067) q[3];
sx q[3];
rz(-0.81334844) q[3];
sx q[3];
rz(0.020997626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.76715604) q[2];
sx q[2];
rz(-2.0220951) q[2];
sx q[2];
rz(-0.68332589) q[2];
rz(-1.4871037) q[3];
sx q[3];
rz(-2.7023102) q[3];
sx q[3];
rz(1.1504415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7651354) q[0];
sx q[0];
rz(-2.6177804) q[0];
sx q[0];
rz(1.8413683) q[0];
rz(0.70872778) q[1];
sx q[1];
rz(-2.6649902) q[1];
sx q[1];
rz(-2.2050819) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4721601) q[0];
sx q[0];
rz(-2.5494808) q[0];
sx q[0];
rz(2.6955312) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2859225) q[2];
sx q[2];
rz(-1.7459918) q[2];
sx q[2];
rz(-2.3223562) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9156956) q[1];
sx q[1];
rz(-2.2749593) q[1];
sx q[1];
rz(-2.8812863) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8285633) q[3];
sx q[3];
rz(-2.1306681) q[3];
sx q[3];
rz(-2.5525639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5986754) q[2];
sx q[2];
rz(-2.8024555) q[2];
sx q[2];
rz(0.014952095) q[2];
rz(-0.22710083) q[3];
sx q[3];
rz(-2.1588219) q[3];
sx q[3];
rz(1.8792413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1508355) q[0];
sx q[0];
rz(-2.3322454) q[0];
sx q[0];
rz(-1.1582751) q[0];
rz(1.5630209) q[1];
sx q[1];
rz(-2.38588) q[1];
sx q[1];
rz(-0.26185782) q[1];
rz(-1.0583393) q[2];
sx q[2];
rz(-2.8556311) q[2];
sx q[2];
rz(-0.43559504) q[2];
rz(-1.3068009) q[3];
sx q[3];
rz(-2.0859857) q[3];
sx q[3];
rz(1.2386238) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
