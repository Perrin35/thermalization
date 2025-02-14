OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3736149) q[0];
sx q[0];
rz(-0.64426214) q[0];
sx q[0];
rz(-0.86384073) q[0];
rz(-1.5683501) q[1];
sx q[1];
rz(-1.22998) q[1];
sx q[1];
rz(-0.65043515) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5395376) q[0];
sx q[0];
rz(-1.3663174) q[0];
sx q[0];
rz(0.50907593) q[0];
rz(-pi) q[1];
rz(3.0502599) q[2];
sx q[2];
rz(-0.31158456) q[2];
sx q[2];
rz(0.1908737) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.12604683) q[1];
sx q[1];
rz(-2.0728557) q[1];
sx q[1];
rz(2.253336) q[1];
rz(-pi) q[2];
rz(-0.69919805) q[3];
sx q[3];
rz(-0.16575925) q[3];
sx q[3];
rz(2.0729271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9170561) q[2];
sx q[2];
rz(-1.7660331) q[2];
sx q[2];
rz(1.3517693) q[2];
rz(-2.3228862) q[3];
sx q[3];
rz(-1.4592146) q[3];
sx q[3];
rz(1.5514577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8390035) q[0];
sx q[0];
rz(-2.4636457) q[0];
sx q[0];
rz(-0.57574058) q[0];
rz(-0.88306824) q[1];
sx q[1];
rz(-1.6022316) q[1];
sx q[1];
rz(1.8227122) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42717375) q[0];
sx q[0];
rz(-1.321081) q[0];
sx q[0];
rz(2.3840133) q[0];
rz(-pi) q[1];
rz(0.78470604) q[2];
sx q[2];
rz(-1.5883351) q[2];
sx q[2];
rz(0.49388805) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.79440642) q[1];
sx q[1];
rz(-0.81450317) q[1];
sx q[1];
rz(1.0760572) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8208597) q[3];
sx q[3];
rz(-1.7400636) q[3];
sx q[3];
rz(1.8603659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.64112249) q[2];
sx q[2];
rz(-2.1239855) q[2];
sx q[2];
rz(-0.23915954) q[2];
rz(0.33254361) q[3];
sx q[3];
rz(-1.4556689) q[3];
sx q[3];
rz(1.0913764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7964433) q[0];
sx q[0];
rz(-0.74757663) q[0];
sx q[0];
rz(2.8969966) q[0];
rz(0.39464513) q[1];
sx q[1];
rz(-2.5375745) q[1];
sx q[1];
rz(-1.2044725) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3217425) q[0];
sx q[0];
rz(-1.5816147) q[0];
sx q[0];
rz(-3.130581) q[0];
x q[1];
rz(2.7543805) q[2];
sx q[2];
rz(-2.0287598) q[2];
sx q[2];
rz(-0.94113038) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.50379163) q[1];
sx q[1];
rz(-1.0560913) q[1];
sx q[1];
rz(-0.33333244) q[1];
rz(0.40562907) q[3];
sx q[3];
rz(-1.105471) q[3];
sx q[3];
rz(2.3387208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.34387732) q[2];
sx q[2];
rz(-2.158973) q[2];
sx q[2];
rz(-2.9003411) q[2];
rz(1.3681715) q[3];
sx q[3];
rz(-1.3894812) q[3];
sx q[3];
rz(2.1696137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57206804) q[0];
sx q[0];
rz(-1.854874) q[0];
sx q[0];
rz(2.0401814) q[0];
rz(-1.6481579) q[1];
sx q[1];
rz(-2.8547574) q[1];
sx q[1];
rz(2.1410087) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7582486) q[0];
sx q[0];
rz(-0.42294082) q[0];
sx q[0];
rz(-1.7420041) q[0];
rz(1.5865636) q[2];
sx q[2];
rz(-0.55634275) q[2];
sx q[2];
rz(0.76130262) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8297208) q[1];
sx q[1];
rz(-1.545136) q[1];
sx q[1];
rz(3.0405634) q[1];
rz(-pi) q[2];
rz(-1.1492997) q[3];
sx q[3];
rz(-1.7022082) q[3];
sx q[3];
rz(-0.67883713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3558041) q[2];
sx q[2];
rz(-1.641909) q[2];
sx q[2];
rz(-0.22264063) q[2];
rz(1.4786134) q[3];
sx q[3];
rz(-0.40545774) q[3];
sx q[3];
rz(-1.2610029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3123689) q[0];
sx q[0];
rz(-2.0033422) q[0];
sx q[0];
rz(-2.9421222) q[0];
rz(-1.5169187) q[1];
sx q[1];
rz(-1.3360887) q[1];
sx q[1];
rz(0.46474251) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9896444) q[0];
sx q[0];
rz(-1.2825263) q[0];
sx q[0];
rz(-3.0040859) q[0];
rz(-pi) q[1];
rz(-0.22146341) q[2];
sx q[2];
rz(-1.5128947) q[2];
sx q[2];
rz(-1.0137978) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6256959) q[1];
sx q[1];
rz(-2.7346238) q[1];
sx q[1];
rz(0.76677983) q[1];
rz(-0.0031849234) q[3];
sx q[3];
rz(-2.4441458) q[3];
sx q[3];
rz(1.7980597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0988079) q[2];
sx q[2];
rz(-1.3036737) q[2];
sx q[2];
rz(-2.888077) q[2];
rz(2.2809196) q[3];
sx q[3];
rz(-0.52152514) q[3];
sx q[3];
rz(2.0025939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0694224) q[0];
sx q[0];
rz(-1.5266029) q[0];
sx q[0];
rz(1.7813064) q[0];
rz(2.5379429) q[1];
sx q[1];
rz(-0.79704469) q[1];
sx q[1];
rz(-2.6545677) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85350689) q[0];
sx q[0];
rz(-0.65847337) q[0];
sx q[0];
rz(3.0229324) q[0];
x q[1];
rz(1.4757882) q[2];
sx q[2];
rz(-0.83432799) q[2];
sx q[2];
rz(0.67259865) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4364061) q[1];
sx q[1];
rz(-1.7366221) q[1];
sx q[1];
rz(-1.2064183) q[1];
rz(-pi) q[2];
rz(1.9726874) q[3];
sx q[3];
rz(-1.783566) q[3];
sx q[3];
rz(-2.5032942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6765678) q[2];
sx q[2];
rz(-1.0176696) q[2];
sx q[2];
rz(2.6505995) q[2];
rz(-2.6077014) q[3];
sx q[3];
rz(-0.54755727) q[3];
sx q[3];
rz(3.0965613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.295306) q[0];
sx q[0];
rz(-2.5745109) q[0];
sx q[0];
rz(-1.2714161) q[0];
rz(-2.2170587) q[1];
sx q[1];
rz(-1.9969767) q[1];
sx q[1];
rz(-2.1163993) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7780925) q[0];
sx q[0];
rz(-0.94787593) q[0];
sx q[0];
rz(1.3062551) q[0];
x q[1];
rz(1.2751638) q[2];
sx q[2];
rz(-0.64648333) q[2];
sx q[2];
rz(1.6103086) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5793663) q[1];
sx q[1];
rz(-2.0831897) q[1];
sx q[1];
rz(2.775869) q[1];
rz(-pi) q[2];
rz(-0.7067755) q[3];
sx q[3];
rz(-1.3825582) q[3];
sx q[3];
rz(-2.1911603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9292494) q[2];
sx q[2];
rz(-0.29444567) q[2];
sx q[2];
rz(-0.92246169) q[2];
rz(0.38659066) q[3];
sx q[3];
rz(-1.8542733) q[3];
sx q[3];
rz(-1.6283901) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5492822) q[0];
sx q[0];
rz(-3.1070502) q[0];
sx q[0];
rz(2.4336245) q[0];
rz(-2.6225846) q[1];
sx q[1];
rz(-0.52878562) q[1];
sx q[1];
rz(-2.4755898) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7836664) q[0];
sx q[0];
rz(-2.5979498) q[0];
sx q[0];
rz(-1.7500135) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.11614271) q[2];
sx q[2];
rz(-1.2525038) q[2];
sx q[2];
rz(-1.7081941) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4909926) q[1];
sx q[1];
rz(-1.5751936) q[1];
sx q[1];
rz(0.66908361) q[1];
rz(-pi) q[2];
x q[2];
rz(0.95548463) q[3];
sx q[3];
rz(-1.59366) q[3];
sx q[3];
rz(-0.61256204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6718601) q[2];
sx q[2];
rz(-1.9817151) q[2];
sx q[2];
rz(1.0603909) q[2];
rz(-1.4276069) q[3];
sx q[3];
rz(-1.5771882) q[3];
sx q[3];
rz(-0.76712999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1308744) q[0];
sx q[0];
rz(-0.81881443) q[0];
sx q[0];
rz(-0.70097104) q[0];
rz(-2.2531807) q[1];
sx q[1];
rz(-2.7263434) q[1];
sx q[1];
rz(-1.7078687) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5007223) q[0];
sx q[0];
rz(-1.6789915) q[0];
sx q[0];
rz(-2.6263424) q[0];
rz(-1.5979684) q[2];
sx q[2];
rz(-1.9765696) q[2];
sx q[2];
rz(2.9494065) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.82821694) q[1];
sx q[1];
rz(-3.113548) q[1];
sx q[1];
rz(-0.7181479) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4096595) q[3];
sx q[3];
rz(-2.3022392) q[3];
sx q[3];
rz(-1.5459154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1514757) q[2];
sx q[2];
rz(-1.9141804) q[2];
sx q[2];
rz(-1.4768451) q[2];
rz(0.19674033) q[3];
sx q[3];
rz(-1.1367831) q[3];
sx q[3];
rz(2.2018532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6075341) q[0];
sx q[0];
rz(-1.5628096) q[0];
sx q[0];
rz(2.0205355) q[0];
rz(2.6634482) q[1];
sx q[1];
rz(-2.4604764) q[1];
sx q[1];
rz(-0.71472439) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9909458) q[0];
sx q[0];
rz(-2.0523221) q[0];
sx q[0];
rz(1.7981547) q[0];
x q[1];
rz(2.1529229) q[2];
sx q[2];
rz(-2.5091362) q[2];
sx q[2];
rz(2.0159231) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.27036941) q[1];
sx q[1];
rz(-1.1308693) q[1];
sx q[1];
rz(-1.0031071) q[1];
rz(-pi) q[2];
x q[2];
rz(0.86831324) q[3];
sx q[3];
rz(-0.83783093) q[3];
sx q[3];
rz(-0.7759076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5152682) q[2];
sx q[2];
rz(-1.3268665) q[2];
sx q[2];
rz(1.2088306) q[2];
rz(-1.442046) q[3];
sx q[3];
rz(-0.88651005) q[3];
sx q[3];
rz(-0.065571688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47990738) q[0];
sx q[0];
rz(-1.3085145) q[0];
sx q[0];
rz(3.0611962) q[0];
rz(2.5936364) q[1];
sx q[1];
rz(-0.68956551) q[1];
sx q[1];
rz(-1.7908304) q[1];
rz(-2.0399848) q[2];
sx q[2];
rz(-1.3528578) q[2];
sx q[2];
rz(2.8982671) q[2];
rz(-0.73765909) q[3];
sx q[3];
rz(-2.548703) q[3];
sx q[3];
rz(-0.92787837) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
