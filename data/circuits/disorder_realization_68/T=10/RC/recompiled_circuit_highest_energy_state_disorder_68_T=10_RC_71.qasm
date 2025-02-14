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
rz(0.23213586) q[0];
sx q[0];
rz(2.870626) q[0];
sx q[0];
rz(10.597064) q[0];
rz(1.5341893) q[1];
sx q[1];
rz(4.6309957) q[1];
sx q[1];
rz(9.9729436) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48972691) q[0];
sx q[0];
rz(-1.4604367) q[0];
sx q[0];
rz(1.2953899) q[0];
rz(-pi) q[1];
x q[1];
rz(0.51082533) q[2];
sx q[2];
rz(-1.0721954) q[2];
sx q[2];
rz(-0.076534903) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1480268) q[1];
sx q[1];
rz(-1.2669592) q[1];
sx q[1];
rz(-2.6152454) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0788623) q[3];
sx q[3];
rz(-2.0196805) q[3];
sx q[3];
rz(-2.1470438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.24070172) q[2];
sx q[2];
rz(-1.9343932) q[2];
sx q[2];
rz(1.594225) q[2];
rz(-2.2089925) q[3];
sx q[3];
rz(-2.2280732) q[3];
sx q[3];
rz(-2.3368321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43657434) q[0];
sx q[0];
rz(-2.5387634) q[0];
sx q[0];
rz(-0.21827179) q[0];
rz(-1.6328579) q[1];
sx q[1];
rz(-2.2684596) q[1];
sx q[1];
rz(2.0192718) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6607165) q[0];
sx q[0];
rz(-0.97597741) q[0];
sx q[0];
rz(0.52158458) q[0];
rz(-pi) q[1];
rz(0.74064765) q[2];
sx q[2];
rz(-1.7193931) q[2];
sx q[2];
rz(-2.1987178) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6186445) q[1];
sx q[1];
rz(-1.0241593) q[1];
sx q[1];
rz(1.0376105) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2551941) q[3];
sx q[3];
rz(-1.8132993) q[3];
sx q[3];
rz(2.8959003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0108769) q[2];
sx q[2];
rz(-1.9515832) q[2];
sx q[2];
rz(2.7575764) q[2];
rz(2.5859517) q[3];
sx q[3];
rz(-0.58755392) q[3];
sx q[3];
rz(2.2491992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4246849) q[0];
sx q[0];
rz(-2.6838344) q[0];
sx q[0];
rz(-2.6440788) q[0];
rz(2.6566907) q[1];
sx q[1];
rz(-1.5496016) q[1];
sx q[1];
rz(0.44152322) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8988441) q[0];
sx q[0];
rz(-1.4024807) q[0];
sx q[0];
rz(2.2344913) q[0];
rz(1.9633358) q[2];
sx q[2];
rz(-0.52670331) q[2];
sx q[2];
rz(-1.8369499) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.42558266) q[1];
sx q[1];
rz(-2.0507318) q[1];
sx q[1];
rz(0.77253567) q[1];
rz(-pi) q[2];
x q[2];
rz(0.86526342) q[3];
sx q[3];
rz(-0.40832106) q[3];
sx q[3];
rz(2.7029519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8170844) q[2];
sx q[2];
rz(-2.0980947) q[2];
sx q[2];
rz(-2.1868152) q[2];
rz(3.1214516) q[3];
sx q[3];
rz(-2.3432799) q[3];
sx q[3];
rz(0.32625833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37498736) q[0];
sx q[0];
rz(-0.52244455) q[0];
sx q[0];
rz(-0.0047542714) q[0];
rz(-2.7004554) q[1];
sx q[1];
rz(-3.0405567) q[1];
sx q[1];
rz(1.4099247) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4213181) q[0];
sx q[0];
rz(-1.4393115) q[0];
sx q[0];
rz(-1.907503) q[0];
rz(-pi) q[1];
rz(2.8660315) q[2];
sx q[2];
rz(-0.86121923) q[2];
sx q[2];
rz(2.1103275) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.04206229) q[1];
sx q[1];
rz(-0.5501318) q[1];
sx q[1];
rz(2.0127556) q[1];
rz(1.7594654) q[3];
sx q[3];
rz(-1.5975286) q[3];
sx q[3];
rz(-2.5875497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.29177523) q[2];
sx q[2];
rz(-1.8287649) q[2];
sx q[2];
rz(2.4449091) q[2];
rz(1.6126532) q[3];
sx q[3];
rz(-0.63642514) q[3];
sx q[3];
rz(0.67970413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1729537) q[0];
sx q[0];
rz(-0.30690673) q[0];
sx q[0];
rz(0.90721834) q[0];
rz(0.13547678) q[1];
sx q[1];
rz(-1.2716581) q[1];
sx q[1];
rz(-2.4438593) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1013779) q[0];
sx q[0];
rz(-1.4411949) q[0];
sx q[0];
rz(1.965353) q[0];
rz(1.5595857) q[2];
sx q[2];
rz(-0.19489842) q[2];
sx q[2];
rz(2.8417808) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.88416022) q[1];
sx q[1];
rz(-1.8990277) q[1];
sx q[1];
rz(1.3882257) q[1];
rz(-pi) q[2];
rz(1.5248564) q[3];
sx q[3];
rz(-1.5410863) q[3];
sx q[3];
rz(1.8159831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.1873143) q[2];
sx q[2];
rz(-0.66437393) q[2];
sx q[2];
rz(0.3453671) q[2];
rz(2.6464388) q[3];
sx q[3];
rz(-0.59585714) q[3];
sx q[3];
rz(-0.12721795) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0536163) q[0];
sx q[0];
rz(-1.6642267) q[0];
sx q[0];
rz(1.9867058) q[0];
rz(2.3467973) q[1];
sx q[1];
rz(-1.844901) q[1];
sx q[1];
rz(2.6577139) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41394407) q[0];
sx q[0];
rz(-0.92287529) q[0];
sx q[0];
rz(-2.6242058) q[0];
x q[1];
rz(1.0740499) q[2];
sx q[2];
rz(-1.4306465) q[2];
sx q[2];
rz(-1.6983918) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.7495854) q[1];
sx q[1];
rz(-1.3806496) q[1];
sx q[1];
rz(-2.0431678) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.57829469) q[3];
sx q[3];
rz(-1.2242172) q[3];
sx q[3];
rz(0.23988304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.38857073) q[2];
sx q[2];
rz(-1.2349671) q[2];
sx q[2];
rz(2.0335601) q[2];
rz(1.5293416) q[3];
sx q[3];
rz(-0.11950167) q[3];
sx q[3];
rz(2.748238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.954708) q[0];
sx q[0];
rz(-2.055838) q[0];
sx q[0];
rz(0.67436522) q[0];
rz(2.2652594) q[1];
sx q[1];
rz(-0.74462157) q[1];
sx q[1];
rz(-0.032329917) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6276777) q[0];
sx q[0];
rz(-1.6925294) q[0];
sx q[0];
rz(-1.4343778) q[0];
rz(-pi) q[1];
rz(-1.7252847) q[2];
sx q[2];
rz(-1.7927899) q[2];
sx q[2];
rz(0.40419182) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9738071) q[1];
sx q[1];
rz(-1.2330489) q[1];
sx q[1];
rz(-2.5925497) q[1];
rz(-0.97801925) q[3];
sx q[3];
rz(-1.7306574) q[3];
sx q[3];
rz(2.7474952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.18675599) q[2];
sx q[2];
rz(-2.7080471) q[2];
sx q[2];
rz(2.1046861) q[2];
rz(-1.9131276) q[3];
sx q[3];
rz(-2.2322673) q[3];
sx q[3];
rz(2.9508446) q[3];
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
rz(-2.4036338) q[0];
sx q[0];
rz(-0.027733138) q[0];
sx q[0];
rz(-2.8926358) q[0];
rz(0.20949334) q[1];
sx q[1];
rz(-1.341235) q[1];
sx q[1];
rz(0.6627717) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4512206) q[0];
sx q[0];
rz(-1.5826384) q[0];
sx q[0];
rz(-1.5607587) q[0];
rz(-2.3784654) q[2];
sx q[2];
rz(-1.9312806) q[2];
sx q[2];
rz(0.57661001) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6697093) q[1];
sx q[1];
rz(-0.85195573) q[1];
sx q[1];
rz(-2.813176) q[1];
rz(-pi) q[2];
rz(1.2609188) q[3];
sx q[3];
rz(-2.204748) q[3];
sx q[3];
rz(0.11989633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2304307) q[2];
sx q[2];
rz(-1.570236) q[2];
sx q[2];
rz(-2.7479808) q[2];
rz(-0.39129928) q[3];
sx q[3];
rz(-2.6063882) q[3];
sx q[3];
rz(-0.75214255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0126208) q[0];
sx q[0];
rz(-2.9550645) q[0];
sx q[0];
rz(-2.5616126) q[0];
rz(-0.36803666) q[1];
sx q[1];
rz(-0.8152222) q[1];
sx q[1];
rz(-0.11485242) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84282473) q[0];
sx q[0];
rz(-1.4359763) q[0];
sx q[0];
rz(-2.299286) q[0];
rz(-2.9909953) q[2];
sx q[2];
rz(-1.1335177) q[2];
sx q[2];
rz(2.8018746) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1424574) q[1];
sx q[1];
rz(-1.8977381) q[1];
sx q[1];
rz(1.8349951) q[1];
rz(-pi) q[2];
rz(-0.87646342) q[3];
sx q[3];
rz(-1.8871347) q[3];
sx q[3];
rz(2.9522459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7241235) q[2];
sx q[2];
rz(-0.33895156) q[2];
sx q[2];
rz(2.4453956) q[2];
rz(-2.0255069) q[3];
sx q[3];
rz(-0.79056549) q[3];
sx q[3];
rz(2.4000786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29257947) q[0];
sx q[0];
rz(-0.96939033) q[0];
sx q[0];
rz(-0.27266362) q[0];
rz(0.69372454) q[1];
sx q[1];
rz(-1.1169746) q[1];
sx q[1];
rz(2.5781217) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54105896) q[0];
sx q[0];
rz(-0.1002914) q[0];
sx q[0];
rz(1.1931399) q[0];
x q[1];
rz(2.5921037) q[2];
sx q[2];
rz(-0.27590431) q[2];
sx q[2];
rz(0.23978182) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.60937) q[1];
sx q[1];
rz(-0.32222363) q[1];
sx q[1];
rz(1.8300309) q[1];
x q[2];
rz(0.93592398) q[3];
sx q[3];
rz(-2.1800123) q[3];
sx q[3];
rz(0.09906957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1571265) q[2];
sx q[2];
rz(-2.0968585) q[2];
sx q[2];
rz(0.54894471) q[2];
rz(2.1460311) q[3];
sx q[3];
rz(-0.82058161) q[3];
sx q[3];
rz(-2.5115749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41848771) q[0];
sx q[0];
rz(-0.99011078) q[0];
sx q[0];
rz(-0.44575442) q[0];
rz(0.86722974) q[1];
sx q[1];
rz(-2.0384616) q[1];
sx q[1];
rz(-1.9834317) q[1];
rz(0.42648496) q[2];
sx q[2];
rz(-2.6250962) q[2];
sx q[2];
rz(0.442183) q[2];
rz(1.4250725) q[3];
sx q[3];
rz(-2.1139664) q[3];
sx q[3];
rz(-2.3900087) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
