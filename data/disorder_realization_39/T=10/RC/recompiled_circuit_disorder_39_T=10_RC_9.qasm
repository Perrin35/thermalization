OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.0467779) q[0];
sx q[0];
rz(-1.0682286) q[0];
sx q[0];
rz(2.6775223) q[0];
rz(1.9595454) q[1];
sx q[1];
rz(6.2160677) q[1];
sx q[1];
rz(10.709229) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0425134) q[0];
sx q[0];
rz(-2.6455542) q[0];
sx q[0];
rz(-2.2201559) q[0];
rz(-pi) q[1];
rz(-1.60381) q[2];
sx q[2];
rz(-1.0812949) q[2];
sx q[2];
rz(-0.93544338) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7317036) q[1];
sx q[1];
rz(-1.1176425) q[1];
sx q[1];
rz(-1.975592) q[1];
x q[2];
rz(3.0331217) q[3];
sx q[3];
rz(-0.15158187) q[3];
sx q[3];
rz(0.82419318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7444732) q[2];
sx q[2];
rz(-2.2019272) q[2];
sx q[2];
rz(2.822067) q[2];
rz(0.5784936) q[3];
sx q[3];
rz(-2.6632023) q[3];
sx q[3];
rz(0.67392504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4085061) q[0];
sx q[0];
rz(-1.5717614) q[0];
sx q[0];
rz(-0.52655667) q[0];
rz(-0.56354848) q[1];
sx q[1];
rz(-1.5310409) q[1];
sx q[1];
rz(-0.79663509) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5986886) q[0];
sx q[0];
rz(-1.2447378) q[0];
sx q[0];
rz(-0.90755264) q[0];
rz(0.30226207) q[2];
sx q[2];
rz(-0.59603359) q[2];
sx q[2];
rz(-2.9994534) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7044201) q[1];
sx q[1];
rz(-2.17802) q[1];
sx q[1];
rz(1.5551268) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1809686) q[3];
sx q[3];
rz(-1.1517797) q[3];
sx q[3];
rz(-0.39715365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.18156302) q[2];
sx q[2];
rz(-1.7627565) q[2];
sx q[2];
rz(2.3550418) q[2];
rz(-2.6484047) q[3];
sx q[3];
rz(-1.9599873) q[3];
sx q[3];
rz(0.44737679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2717473) q[0];
sx q[0];
rz(-2.2646876) q[0];
sx q[0];
rz(1.440381) q[0];
rz(-0.72021833) q[1];
sx q[1];
rz(-2.5455988) q[1];
sx q[1];
rz(0.70297855) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5876578) q[0];
sx q[0];
rz(-0.30447391) q[0];
sx q[0];
rz(-1.1712043) q[0];
rz(-pi) q[1];
rz(-2.5419652) q[2];
sx q[2];
rz(-2.5587974) q[2];
sx q[2];
rz(1.8011013) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8382593) q[1];
sx q[1];
rz(-1.5240182) q[1];
sx q[1];
rz(0.10951885) q[1];
rz(-0.94743518) q[3];
sx q[3];
rz(-1.7121592) q[3];
sx q[3];
rz(-2.3879104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8748223) q[2];
sx q[2];
rz(-1.4845994) q[2];
sx q[2];
rz(0.91153574) q[2];
rz(-0.91397816) q[3];
sx q[3];
rz(-1.711288) q[3];
sx q[3];
rz(-2.3142557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5660969) q[0];
sx q[0];
rz(-0.63593447) q[0];
sx q[0];
rz(0.77600586) q[0];
rz(-1.2674747) q[1];
sx q[1];
rz(-1.1088561) q[1];
sx q[1];
rz(-2.5783096) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.027772) q[0];
sx q[0];
rz(-1.2978683) q[0];
sx q[0];
rz(0.91822894) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2345384) q[2];
sx q[2];
rz(-0.90869892) q[2];
sx q[2];
rz(-1.7788356) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5967484) q[1];
sx q[1];
rz(-1.2990284) q[1];
sx q[1];
rz(1.5775058) q[1];
x q[2];
rz(2.8068845) q[3];
sx q[3];
rz(-2.3957806) q[3];
sx q[3];
rz(-2.0832182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6528066) q[2];
sx q[2];
rz(-0.93549171) q[2];
sx q[2];
rz(2.9768067) q[2];
rz(-0.22848836) q[3];
sx q[3];
rz(-2.861664) q[3];
sx q[3];
rz(-0.94648186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27424681) q[0];
sx q[0];
rz(-0.86425257) q[0];
sx q[0];
rz(-2.2891323) q[0];
rz(0.35119855) q[1];
sx q[1];
rz(-0.57792592) q[1];
sx q[1];
rz(1.16211) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0162109) q[0];
sx q[0];
rz(-2.1124766) q[0];
sx q[0];
rz(1.1075695) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.39761333) q[2];
sx q[2];
rz(-1.0828472) q[2];
sx q[2];
rz(0.60148009) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4967242) q[1];
sx q[1];
rz(-0.86360303) q[1];
sx q[1];
rz(-1.0134539) q[1];
rz(-pi) q[2];
rz(-1.5557628) q[3];
sx q[3];
rz(-2.2923277) q[3];
sx q[3];
rz(0.3258457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6953485) q[2];
sx q[2];
rz(-2.7480795) q[2];
sx q[2];
rz(1.2472786) q[2];
rz(3.128483) q[3];
sx q[3];
rz(-1.7316827) q[3];
sx q[3];
rz(0.22918992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2728249) q[0];
sx q[0];
rz(-1.8936963) q[0];
sx q[0];
rz(0.75063467) q[0];
rz(-2.0320832) q[1];
sx q[1];
rz(-0.37056071) q[1];
sx q[1];
rz(-1.8355339) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7795777) q[0];
sx q[0];
rz(-1.0581731) q[0];
sx q[0];
rz(-2.6984452) q[0];
x q[1];
rz(-0.38988955) q[2];
sx q[2];
rz(-0.49553686) q[2];
sx q[2];
rz(1.748566) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8522779) q[1];
sx q[1];
rz(-1.5783974) q[1];
sx q[1];
rz(0.72035933) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9206198) q[3];
sx q[3];
rz(-0.60022012) q[3];
sx q[3];
rz(2.5252987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3874454) q[2];
sx q[2];
rz(-2.1345317) q[2];
sx q[2];
rz(-0.60069096) q[2];
rz(1.0026275) q[3];
sx q[3];
rz(-0.53835821) q[3];
sx q[3];
rz(-0.35282648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7464741) q[0];
sx q[0];
rz(-1.4557319) q[0];
sx q[0];
rz(1.0466928) q[0];
rz(-1.612161) q[1];
sx q[1];
rz(-1.7144831) q[1];
sx q[1];
rz(-2.7244862) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1244303) q[0];
sx q[0];
rz(-0.22148795) q[0];
sx q[0];
rz(1.682196) q[0];
x q[1];
rz(-2.8344526) q[2];
sx q[2];
rz(-1.2899613) q[2];
sx q[2];
rz(-1.839523) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6542146) q[1];
sx q[1];
rz(-1.1114792) q[1];
sx q[1];
rz(-2.3801801) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8700637) q[3];
sx q[3];
rz(-2.2021658) q[3];
sx q[3];
rz(-2.805998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.0059011857) q[2];
sx q[2];
rz(-2.6689172) q[2];
sx q[2];
rz(-1.7741514) q[2];
rz(-2.588429) q[3];
sx q[3];
rz(-1.7382712) q[3];
sx q[3];
rz(0.52545351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32507867) q[0];
sx q[0];
rz(-1.8186318) q[0];
sx q[0];
rz(1.1897855) q[0];
rz(1.7143543) q[1];
sx q[1];
rz(-2.863527) q[1];
sx q[1];
rz(-0.11238012) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1988797) q[0];
sx q[0];
rz(-0.50857022) q[0];
sx q[0];
rz(1.2099427) q[0];
x q[1];
rz(2.8078812) q[2];
sx q[2];
rz(-1.6534272) q[2];
sx q[2];
rz(-0.73252788) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.97265128) q[1];
sx q[1];
rz(-0.87247889) q[1];
sx q[1];
rz(0.45291839) q[1];
x q[2];
rz(1.3120679) q[3];
sx q[3];
rz(-1.6885763) q[3];
sx q[3];
rz(-0.044101276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0671976) q[2];
sx q[2];
rz(-1.9125568) q[2];
sx q[2];
rz(1.3486264) q[2];
rz(-1.9366692) q[3];
sx q[3];
rz(-1.7347074) q[3];
sx q[3];
rz(-2.9437734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39712054) q[0];
sx q[0];
rz(-1.4646126) q[0];
sx q[0];
rz(-0.034974139) q[0];
rz(-2.2947625) q[1];
sx q[1];
rz(-1.7085107) q[1];
sx q[1];
rz(0.91167489) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1393136) q[0];
sx q[0];
rz(-2.7873435) q[0];
sx q[0];
rz(2.2646963) q[0];
rz(2.9334925) q[2];
sx q[2];
rz(-2.5884183) q[2];
sx q[2];
rz(-0.91536872) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.72494353) q[1];
sx q[1];
rz(-1.1133725) q[1];
sx q[1];
rz(-1.3436505) q[1];
rz(-pi) q[2];
rz(0.53578844) q[3];
sx q[3];
rz(-1.1402604) q[3];
sx q[3];
rz(0.99539869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.70790616) q[2];
sx q[2];
rz(-2.6277379) q[2];
sx q[2];
rz(-0.24924499) q[2];
rz(0.76672673) q[3];
sx q[3];
rz(-1.8079575) q[3];
sx q[3];
rz(-2.7419817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3578167) q[0];
sx q[0];
rz(-1.0483085) q[0];
sx q[0];
rz(-0.37208474) q[0];
rz(2.5601939) q[1];
sx q[1];
rz(-0.77722469) q[1];
sx q[1];
rz(-1.4153597) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5370731) q[0];
sx q[0];
rz(-1.8329289) q[0];
sx q[0];
rz(-1.3634691) q[0];
x q[1];
rz(-2.084311) q[2];
sx q[2];
rz(-1.0143177) q[2];
sx q[2];
rz(1.8884115) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2409754) q[1];
sx q[1];
rz(-1.4957349) q[1];
sx q[1];
rz(-1.6402628) q[1];
rz(-pi) q[2];
x q[2];
rz(0.98202242) q[3];
sx q[3];
rz(-2.5128799) q[3];
sx q[3];
rz(0.94129896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.32594484) q[2];
sx q[2];
rz(-0.14020136) q[2];
sx q[2];
rz(-0.20467219) q[2];
rz(-1.4137911) q[3];
sx q[3];
rz(-1.9982343) q[3];
sx q[3];
rz(1.0958825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(0.50080147) q[0];
sx q[0];
rz(-1.4688107) q[0];
sx q[0];
rz(-2.2039913) q[0];
rz(1.5851371) q[1];
sx q[1];
rz(-1.7592447) q[1];
sx q[1];
rz(-1.6978282) q[1];
rz(-1.3265058) q[2];
sx q[2];
rz(-0.46637022) q[2];
sx q[2];
rz(-2.7844219) q[2];
rz(1.6987883) q[3];
sx q[3];
rz(-0.57590579) q[3];
sx q[3];
rz(-0.87344575) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];