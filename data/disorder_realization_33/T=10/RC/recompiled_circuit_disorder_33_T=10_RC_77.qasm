OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6361976) q[0];
sx q[0];
rz(6.0072748) q[0];
sx q[0];
rz(10.732565) q[0];
rz(1.1360599) q[1];
sx q[1];
rz(-0.93568957) q[1];
sx q[1];
rz(1.5703262) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9309064) q[0];
sx q[0];
rz(-1.9084198) q[0];
sx q[0];
rz(2.7772285) q[0];
rz(-pi) q[1];
rz(2.3762796) q[2];
sx q[2];
rz(-2.1047449) q[2];
sx q[2];
rz(3.090976) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.624144) q[1];
sx q[1];
rz(-1.9074719) q[1];
sx q[1];
rz(1.3144073) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0516112) q[3];
sx q[3];
rz(-2.7365723) q[3];
sx q[3];
rz(0.68457505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2661665) q[2];
sx q[2];
rz(-2.8484919) q[2];
sx q[2];
rz(1.1323294) q[2];
rz(1.6752361) q[3];
sx q[3];
rz(-1.8050067) q[3];
sx q[3];
rz(-1.0124538) q[3];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19673008) q[0];
sx q[0];
rz(-0.20962993) q[0];
sx q[0];
rz(0.18584132) q[0];
rz(0.56022412) q[1];
sx q[1];
rz(-1.8461684) q[1];
sx q[1];
rz(-2.9247608) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8502055) q[0];
sx q[0];
rz(-2.4225525) q[0];
sx q[0];
rz(-2.015381) q[0];
rz(-pi) q[1];
x q[1];
rz(0.47302834) q[2];
sx q[2];
rz(-2.0748667) q[2];
sx q[2];
rz(-2.8216528) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.20570457) q[1];
sx q[1];
rz(-1.7824031) q[1];
sx q[1];
rz(0.88797027) q[1];
rz(2.2651477) q[3];
sx q[3];
rz(-1.0507686) q[3];
sx q[3];
rz(-2.2976573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.310114) q[2];
sx q[2];
rz(-0.82565132) q[2];
sx q[2];
rz(-1.8537834) q[2];
rz(2.3790322) q[3];
sx q[3];
rz(-1.1688787) q[3];
sx q[3];
rz(2.8365703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4644311) q[0];
sx q[0];
rz(-0.34496775) q[0];
sx q[0];
rz(2.537354) q[0];
rz(-1.8151981) q[1];
sx q[1];
rz(-1.7809968) q[1];
sx q[1];
rz(2.2089829) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5843825) q[0];
sx q[0];
rz(-0.32214468) q[0];
sx q[0];
rz(-1.4003217) q[0];
rz(0.18205299) q[2];
sx q[2];
rz(-1.3472054) q[2];
sx q[2];
rz(-1.6857266) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.697726) q[1];
sx q[1];
rz(-1.4523456) q[1];
sx q[1];
rz(-0.5603793) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.014939) q[3];
sx q[3];
rz(-1.1851289) q[3];
sx q[3];
rz(-0.25394299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9937667) q[2];
sx q[2];
rz(-2.0596762) q[2];
sx q[2];
rz(2.0489342) q[2];
rz(-0.5422194) q[3];
sx q[3];
rz(-2.0565624) q[3];
sx q[3];
rz(0.96737635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7820691) q[0];
sx q[0];
rz(-3.0451267) q[0];
sx q[0];
rz(-2.6413667) q[0];
rz(2.3362828) q[1];
sx q[1];
rz(-1.9814682) q[1];
sx q[1];
rz(-1.4979699) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7786176) q[0];
sx q[0];
rz(-2.551429) q[0];
sx q[0];
rz(-1.0233364) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0849886) q[2];
sx q[2];
rz(-1.3805693) q[2];
sx q[2];
rz(-2.9375926) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5768891) q[1];
sx q[1];
rz(-0.35208811) q[1];
sx q[1];
rz(-2.4114154) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0235396) q[3];
sx q[3];
rz(-1.9948043) q[3];
sx q[3];
rz(-0.31183576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3952289) q[2];
sx q[2];
rz(-2.5791898) q[2];
sx q[2];
rz(-2.4397819) q[2];
rz(-2.3102405) q[3];
sx q[3];
rz(-2.1777007) q[3];
sx q[3];
rz(2.519616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(0.9005301) q[0];
sx q[0];
rz(-2.5456972) q[0];
sx q[0];
rz(-2.3262614) q[0];
rz(-1.6197846) q[1];
sx q[1];
rz(-0.83414572) q[1];
sx q[1];
rz(1.048208) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74131504) q[0];
sx q[0];
rz(-1.321723) q[0];
sx q[0];
rz(1.2988017) q[0];
x q[1];
rz(-1.5585209) q[2];
sx q[2];
rz(-0.91931146) q[2];
sx q[2];
rz(-0.94142454) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9427467) q[1];
sx q[1];
rz(-1.2578576) q[1];
sx q[1];
rz(-1.8206157) q[1];
x q[2];
rz(-2.1861595) q[3];
sx q[3];
rz(-1.2293929) q[3];
sx q[3];
rz(-2.0379025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6158225) q[2];
sx q[2];
rz(-2.5746391) q[2];
sx q[2];
rz(-2.0416416) q[2];
rz(2.3163017) q[3];
sx q[3];
rz(-2.0992978) q[3];
sx q[3];
rz(-2.2560789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0734171) q[0];
sx q[0];
rz(-0.59403479) q[0];
sx q[0];
rz(-0.90240479) q[0];
rz(-2.1249318) q[1];
sx q[1];
rz(-1.0598176) q[1];
sx q[1];
rz(3.0117603) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6274525) q[0];
sx q[0];
rz(-0.99533004) q[0];
sx q[0];
rz(2.0155725) q[0];
rz(1.0340704) q[2];
sx q[2];
rz(-0.64646361) q[2];
sx q[2];
rz(1.549364) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1631158) q[1];
sx q[1];
rz(-0.69677959) q[1];
sx q[1];
rz(2.5549868) q[1];
rz(-pi) q[2];
rz(-0.31452175) q[3];
sx q[3];
rz(-2.5701227) q[3];
sx q[3];
rz(2.275327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8292024) q[2];
sx q[2];
rz(-2.1924993) q[2];
sx q[2];
rz(-2.9373346) q[2];
rz(1.9355109) q[3];
sx q[3];
rz(-1.6198502) q[3];
sx q[3];
rz(2.9061785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4181353) q[0];
sx q[0];
rz(-1.8122939) q[0];
sx q[0];
rz(-1.4468505) q[0];
rz(-1.8824668) q[1];
sx q[1];
rz(-0.99021688) q[1];
sx q[1];
rz(0.68626219) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2080363) q[0];
sx q[0];
rz(-0.66970034) q[0];
sx q[0];
rz(2.3963388) q[0];
rz(-pi) q[1];
rz(-2.5574066) q[2];
sx q[2];
rz(-0.91911941) q[2];
sx q[2];
rz(0.78648957) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2346748) q[1];
sx q[1];
rz(-1.385681) q[1];
sx q[1];
rz(3.0659552) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7191914) q[3];
sx q[3];
rz(-1.2761315) q[3];
sx q[3];
rz(-0.15448031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.69616047) q[2];
sx q[2];
rz(-1.7636718) q[2];
sx q[2];
rz(-3.1398204) q[2];
rz(-2.5799675) q[3];
sx q[3];
rz(-0.91149819) q[3];
sx q[3];
rz(-1.5047489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5381662) q[0];
sx q[0];
rz(-0.68646938) q[0];
sx q[0];
rz(1.4461393) q[0];
rz(-0.7810477) q[1];
sx q[1];
rz(-1.8361517) q[1];
sx q[1];
rz(1.5015645) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9672464) q[0];
sx q[0];
rz(-1.5513199) q[0];
sx q[0];
rz(-2.0635701) q[0];
x q[1];
rz(1.1649706) q[2];
sx q[2];
rz(-0.067194447) q[2];
sx q[2];
rz(0.42025987) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0280684) q[1];
sx q[1];
rz(-1.4816195) q[1];
sx q[1];
rz(-2.813617) q[1];
x q[2];
rz(-0.73236671) q[3];
sx q[3];
rz(-2.3673956) q[3];
sx q[3];
rz(-2.458651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7897196) q[2];
sx q[2];
rz(-1.3616273) q[2];
sx q[2];
rz(1.8224576) q[2];
rz(-1.9296648) q[3];
sx q[3];
rz(-1.2865678) q[3];
sx q[3];
rz(-2.8222728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.33655745) q[0];
sx q[0];
rz(-0.55258495) q[0];
sx q[0];
rz(-1.2040899) q[0];
rz(2.7583292) q[1];
sx q[1];
rz(-0.52572322) q[1];
sx q[1];
rz(-0.35167545) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8512745) q[0];
sx q[0];
rz(-1.1997249) q[0];
sx q[0];
rz(2.0593658) q[0];
rz(-pi) q[1];
rz(2.4497689) q[2];
sx q[2];
rz(-1.3335388) q[2];
sx q[2];
rz(2.0123864) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0723567) q[1];
sx q[1];
rz(-1.623739) q[1];
sx q[1];
rz(-0.80979053) q[1];
rz(-pi) q[2];
rz(-3.0319801) q[3];
sx q[3];
rz(-1.622756) q[3];
sx q[3];
rz(-2.6250641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.3433156) q[2];
sx q[2];
rz(-2.0337992) q[2];
sx q[2];
rz(1.2822255) q[2];
rz(1.4964237) q[3];
sx q[3];
rz(-1.6069501) q[3];
sx q[3];
rz(-2.0857247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6431817) q[0];
sx q[0];
rz(-1.8739941) q[0];
sx q[0];
rz(-2.9472651) q[0];
rz(-1.0378029) q[1];
sx q[1];
rz(-0.56832814) q[1];
sx q[1];
rz(-2.1077572) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94668418) q[0];
sx q[0];
rz(-0.91741981) q[0];
sx q[0];
rz(-2.9772467) q[0];
rz(-pi) q[1];
rz(-2.6331484) q[2];
sx q[2];
rz(-0.93548453) q[2];
sx q[2];
rz(2.9490162) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2227877) q[1];
sx q[1];
rz(-0.58090392) q[1];
sx q[1];
rz(1.7564299) q[1];
rz(-pi) q[2];
x q[2];
rz(0.25063534) q[3];
sx q[3];
rz(-2.4126629) q[3];
sx q[3];
rz(-2.7540516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0620492) q[2];
sx q[2];
rz(-0.94576183) q[2];
sx q[2];
rz(-2.5058084) q[2];
rz(0.27030269) q[3];
sx q[3];
rz(-2.342194) q[3];
sx q[3];
rz(1.5283782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4476267) q[0];
sx q[0];
rz(-1.8287369) q[0];
sx q[0];
rz(1.0736314) q[0];
rz(-1.7059965) q[1];
sx q[1];
rz(-1.5789079) q[1];
sx q[1];
rz(0.78067738) q[1];
rz(-3.1096935) q[2];
sx q[2];
rz(-0.96822856) q[2];
sx q[2];
rz(-0.4005489) q[2];
rz(1.3676436) q[3];
sx q[3];
rz(-2.805134) q[3];
sx q[3];
rz(0.72611879) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
