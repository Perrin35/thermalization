OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.0948148) q[0];
sx q[0];
rz(-2.0733641) q[0];
sx q[0];
rz(0.46407035) q[0];
rz(-1.1820473) q[1];
sx q[1];
rz(-3.074475) q[1];
sx q[1];
rz(-1.2844515) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0425134) q[0];
sx q[0];
rz(-0.49603841) q[0];
sx q[0];
rz(2.2201559) q[0];
rz(-pi) q[1];
x q[1];
rz(0.061878248) q[2];
sx q[2];
rz(-0.49052325) q[2];
sx q[2];
rz(2.2762736) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1661108) q[1];
sx q[1];
rz(-1.9327285) q[1];
sx q[1];
rz(-0.48717498) q[1];
x q[2];
rz(-1.5873317) q[3];
sx q[3];
rz(-1.7214805) q[3];
sx q[3];
rz(-2.2076803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.39711943) q[2];
sx q[2];
rz(-0.93966547) q[2];
sx q[2];
rz(-2.822067) q[2];
rz(-2.5630991) q[3];
sx q[3];
rz(-2.6632023) q[3];
sx q[3];
rz(0.67392504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7330866) q[0];
sx q[0];
rz(-1.5717614) q[0];
sx q[0];
rz(-2.615036) q[0];
rz(-2.5780442) q[1];
sx q[1];
rz(-1.6105517) q[1];
sx q[1];
rz(-0.79663509) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54290402) q[0];
sx q[0];
rz(-1.2447378) q[0];
sx q[0];
rz(-2.23404) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.57467069) q[2];
sx q[2];
rz(-1.4029014) q[2];
sx q[2];
rz(-1.4603953) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7044201) q[1];
sx q[1];
rz(-2.17802) q[1];
sx q[1];
rz(1.5551268) q[1];
rz(-pi) q[2];
rz(-2.1809686) q[3];
sx q[3];
rz(-1.1517797) q[3];
sx q[3];
rz(-2.744439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9600296) q[2];
sx q[2];
rz(-1.3788362) q[2];
sx q[2];
rz(2.3550418) q[2];
rz(2.6484047) q[3];
sx q[3];
rz(-1.9599873) q[3];
sx q[3];
rz(2.6942159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86984533) q[0];
sx q[0];
rz(-2.2646876) q[0];
sx q[0];
rz(1.7012117) q[0];
rz(-0.72021833) q[1];
sx q[1];
rz(-2.5455988) q[1];
sx q[1];
rz(-2.4386141) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5876578) q[0];
sx q[0];
rz(-2.8371187) q[0];
sx q[0];
rz(1.9703883) q[0];
x q[1];
rz(2.5419652) q[2];
sx q[2];
rz(-0.58279524) q[2];
sx q[2];
rz(1.8011013) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2726047) q[1];
sx q[1];
rz(-1.4613978) q[1];
sx q[1];
rz(-1.5237367) q[1];
rz(1.3316783) q[3];
sx q[3];
rz(-2.5044887) q[3];
sx q[3];
rz(-2.1309731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8748223) q[2];
sx q[2];
rz(-1.4845994) q[2];
sx q[2];
rz(-0.91153574) q[2];
rz(0.91397816) q[3];
sx q[3];
rz(-1.4303047) q[3];
sx q[3];
rz(0.82733697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5660969) q[0];
sx q[0];
rz(-2.5056582) q[0];
sx q[0];
rz(-0.77600586) q[0];
rz(-1.874118) q[1];
sx q[1];
rz(-1.1088561) q[1];
sx q[1];
rz(-0.56328303) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.027772) q[0];
sx q[0];
rz(-1.2978683) q[0];
sx q[0];
rz(-2.2233637) q[0];
x q[1];
rz(1.9070542) q[2];
sx q[2];
rz(-2.2328937) q[2];
sx q[2];
rz(1.3627571) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.56983419) q[1];
sx q[1];
rz(-0.27184871) q[1];
sx q[1];
rz(0.024072577) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4241583) q[3];
sx q[3];
rz(-1.7955901) q[3];
sx q[3];
rz(2.3790529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6528066) q[2];
sx q[2];
rz(-0.93549171) q[2];
sx q[2];
rz(-0.164786) q[2];
rz(-0.22848836) q[3];
sx q[3];
rz(-2.861664) q[3];
sx q[3];
rz(-0.94648186) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8673458) q[0];
sx q[0];
rz(-2.2773401) q[0];
sx q[0];
rz(-2.2891323) q[0];
rz(-2.7903941) q[1];
sx q[1];
rz(-0.57792592) q[1];
sx q[1];
rz(-1.9794827) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0162109) q[0];
sx q[0];
rz(-1.0291161) q[0];
sx q[0];
rz(2.0340232) q[0];
rz(-pi) q[1];
x q[1];
rz(0.94050546) q[2];
sx q[2];
rz(-2.5224707) q[2];
sx q[2];
rz(1.3319912) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.1197549) q[1];
sx q[1];
rz(-0.86958414) q[1];
sx q[1];
rz(-0.55418684) q[1];
x q[2];
rz(-3.1245072) q[3];
sx q[3];
rz(-2.4199329) q[3];
sx q[3];
rz(-0.34860308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.44624415) q[2];
sx q[2];
rz(-2.7480795) q[2];
sx q[2];
rz(-1.8943141) q[2];
rz(-3.128483) q[3];
sx q[3];
rz(-1.7316827) q[3];
sx q[3];
rz(2.9124027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2728249) q[0];
sx q[0];
rz(-1.8936963) q[0];
sx q[0];
rz(2.390958) q[0];
rz(-1.1095095) q[1];
sx q[1];
rz(-2.7710319) q[1];
sx q[1];
rz(1.3060588) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.437498) q[0];
sx q[0];
rz(-1.1878345) q[0];
sx q[0];
rz(-1.0136481) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6779804) q[2];
sx q[2];
rz(-1.3890651) q[2];
sx q[2];
rz(-2.6169427) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.27481025) q[1];
sx q[1];
rz(-0.85046235) q[1];
sx q[1];
rz(1.5809098) q[1];
x q[2];
rz(2.9206198) q[3];
sx q[3];
rz(-2.5413725) q[3];
sx q[3];
rz(0.61629399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3874454) q[2];
sx q[2];
rz(-1.0070609) q[2];
sx q[2];
rz(-0.60069096) q[2];
rz(-1.0026275) q[3];
sx q[3];
rz(-0.53835821) q[3];
sx q[3];
rz(-2.7887662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
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
rz(0.41710645) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0171623) q[0];
sx q[0];
rz(-2.9201047) q[0];
sx q[0];
rz(-1.682196) q[0];
x q[1];
rz(-1.276937) q[2];
sx q[2];
rz(-1.8655348) q[2];
sx q[2];
rz(2.7851832) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6542146) q[1];
sx q[1];
rz(-2.0301135) q[1];
sx q[1];
rz(2.3801801) q[1];
rz(1.2715289) q[3];
sx q[3];
rz(-0.93942681) q[3];
sx q[3];
rz(-2.805998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.0059011857) q[2];
sx q[2];
rz(-2.6689172) q[2];
sx q[2];
rz(1.7741514) q[2];
rz(-2.588429) q[3];
sx q[3];
rz(-1.7382712) q[3];
sx q[3];
rz(0.52545351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
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
rz(-2.816514) q[0];
sx q[0];
rz(-1.8186318) q[0];
sx q[0];
rz(-1.9518071) q[0];
rz(1.4272383) q[1];
sx q[1];
rz(-0.27806565) q[1];
sx q[1];
rz(3.0292125) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6067137) q[0];
sx q[0];
rz(-1.0977854) q[0];
sx q[0];
rz(0.19434778) q[0];
x q[1];
rz(0.24765315) q[2];
sx q[2];
rz(-0.34341771) q[2];
sx q[2];
rz(-2.0695956) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1689414) q[1];
sx q[1];
rz(-0.87247889) q[1];
sx q[1];
rz(-2.6886743) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1376082) q[3];
sx q[3];
rz(-2.8578651) q[3];
sx q[3];
rz(-1.9445436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0743951) q[2];
sx q[2];
rz(-1.9125568) q[2];
sx q[2];
rz(1.7929662) q[2];
rz(1.9366692) q[3];
sx q[3];
rz(-1.7347074) q[3];
sx q[3];
rz(2.9437734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(2.7444721) q[0];
sx q[0];
rz(-1.4646126) q[0];
sx q[0];
rz(-3.1066185) q[0];
rz(-0.84683013) q[1];
sx q[1];
rz(-1.433082) q[1];
sx q[1];
rz(0.91167489) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9104886) q[0];
sx q[0];
rz(-1.3470874) q[0];
sx q[0];
rz(1.847812) q[0];
rz(-pi) q[1];
x q[1];
rz(0.54347221) q[2];
sx q[2];
rz(-1.4620355) q[2];
sx q[2];
rz(0.47765884) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1940143) q[1];
sx q[1];
rz(-1.36735) q[1];
sx q[1];
rz(-2.6737763) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0802286) q[3];
sx q[3];
rz(-1.0883696) q[3];
sx q[3];
rz(0.8182984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4336865) q[2];
sx q[2];
rz(-0.5138548) q[2];
sx q[2];
rz(-2.8923477) q[2];
rz(2.3748659) q[3];
sx q[3];
rz(-1.3336351) q[3];
sx q[3];
rz(0.39961091) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3578167) q[0];
sx q[0];
rz(-2.0932842) q[0];
sx q[0];
rz(0.37208474) q[0];
rz(0.58139873) q[1];
sx q[1];
rz(-2.364368) q[1];
sx q[1];
rz(1.7262329) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5370731) q[0];
sx q[0];
rz(-1.3086638) q[0];
sx q[0];
rz(-1.3634691) q[0];
rz(-0.62016983) q[2];
sx q[2];
rz(-2.0010741) q[2];
sx q[2];
rz(-2.534453) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.806554) q[1];
sx q[1];
rz(-1.5015258) q[1];
sx q[1];
rz(-0.075242234) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1595702) q[3];
sx q[3];
rz(-0.62871274) q[3];
sx q[3];
rz(0.94129896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8156478) q[2];
sx q[2];
rz(-0.14020136) q[2];
sx q[2];
rz(2.9369205) q[2];
rz(1.7278016) q[3];
sx q[3];
rz(-1.1433583) q[3];
sx q[3];
rz(2.0457101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50080147) q[0];
sx q[0];
rz(-1.6727819) q[0];
sx q[0];
rz(0.93760136) q[0];
rz(-1.5851371) q[1];
sx q[1];
rz(-1.382348) q[1];
sx q[1];
rz(1.4437645) q[1];
rz(3.0204308) q[2];
sx q[2];
rz(-1.1193174) q[2];
sx q[2];
rz(0.08502273) q[2];
rz(-0.99863573) q[3];
sx q[3];
rz(-1.5012267) q[3];
sx q[3];
rz(-2.5517626) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];