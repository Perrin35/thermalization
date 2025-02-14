OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.98915339) q[0];
sx q[0];
rz(-1.5664772) q[0];
sx q[0];
rz(2.0164665) q[0];
rz(-2.1195124) q[1];
sx q[1];
rz(-2.4740969) q[1];
sx q[1];
rz(-0.8134841) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73972244) q[0];
sx q[0];
rz(-2.7474294) q[0];
sx q[0];
rz(-1.9456844) q[0];
rz(1.9314693) q[2];
sx q[2];
rz(-1.1481592) q[2];
sx q[2];
rz(-0.3061184) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1853987) q[1];
sx q[1];
rz(-0.46727249) q[1];
sx q[1];
rz(-2.8059792) q[1];
rz(-0.090510719) q[3];
sx q[3];
rz(-1.7337243) q[3];
sx q[3];
rz(2.3855057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.4890613) q[2];
sx q[2];
rz(-1.6901313) q[2];
sx q[2];
rz(0.92864621) q[2];
rz(-1.5993902) q[3];
sx q[3];
rz(-1.3336811) q[3];
sx q[3];
rz(-1.1716051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20392513) q[0];
sx q[0];
rz(-1.7610022) q[0];
sx q[0];
rz(-3.0145338) q[0];
rz(-2.1584885) q[1];
sx q[1];
rz(-1.3652722) q[1];
sx q[1];
rz(0.7712706) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8255804) q[0];
sx q[0];
rz(-2.2716652) q[0];
sx q[0];
rz(-2.6722145) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7423058) q[2];
sx q[2];
rz(-2.7993188) q[2];
sx q[2];
rz(0.31708395) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2576372) q[1];
sx q[1];
rz(-1.2246545) q[1];
sx q[1];
rz(0.37948541) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3985474) q[3];
sx q[3];
rz(-2.4186169) q[3];
sx q[3];
rz(-1.6188177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1473006) q[2];
sx q[2];
rz(-1.2039801) q[2];
sx q[2];
rz(-3.1331983) q[2];
rz(-0.66347915) q[3];
sx q[3];
rz(-1.2365664) q[3];
sx q[3];
rz(0.26507637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0405149) q[0];
sx q[0];
rz(-2.3150257) q[0];
sx q[0];
rz(-2.7000632) q[0];
rz(-2.1532374) q[1];
sx q[1];
rz(-1.1343196) q[1];
sx q[1];
rz(0.13557869) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39865935) q[0];
sx q[0];
rz(-3.1236095) q[0];
sx q[0];
rz(2.6989486) q[0];
x q[1];
rz(1.1889815) q[2];
sx q[2];
rz(-2.2509607) q[2];
sx q[2];
rz(-1.9651246) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6171268) q[1];
sx q[1];
rz(-1.2167261) q[1];
sx q[1];
rz(-1.0431784) q[1];
rz(-pi) q[2];
rz(1.7884004) q[3];
sx q[3];
rz(-0.74624589) q[3];
sx q[3];
rz(-1.7640132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2077937) q[2];
sx q[2];
rz(-1.4321045) q[2];
sx q[2];
rz(0.03820339) q[2];
rz(-0.52538747) q[3];
sx q[3];
rz(-0.62756687) q[3];
sx q[3];
rz(-0.6853404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66048375) q[0];
sx q[0];
rz(-2.3479192) q[0];
sx q[0];
rz(-2.5307181) q[0];
rz(-1.8065709) q[1];
sx q[1];
rz(-1.3742615) q[1];
sx q[1];
rz(0.23922051) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51620519) q[0];
sx q[0];
rz(-0.62407485) q[0];
sx q[0];
rz(1.1136057) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6692713) q[2];
sx q[2];
rz(-1.9237674) q[2];
sx q[2];
rz(-1.2833088) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.94915056) q[1];
sx q[1];
rz(-2.572963) q[1];
sx q[1];
rz(-2.5123764) q[1];
x q[2];
rz(1.3705105) q[3];
sx q[3];
rz(-2.0242175) q[3];
sx q[3];
rz(-3.0016921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7188344) q[2];
sx q[2];
rz(-1.6966635) q[2];
sx q[2];
rz(-0.0017496721) q[2];
rz(0.29378978) q[3];
sx q[3];
rz(-1.1894476) q[3];
sx q[3];
rz(-2.8747115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
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
rz(2.9002429) q[0];
sx q[0];
rz(-1.7647864) q[0];
sx q[0];
rz(-0.84306651) q[0];
rz(-0.33755606) q[1];
sx q[1];
rz(-1.2755716) q[1];
sx q[1];
rz(1.5379803) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.085097236) q[0];
sx q[0];
rz(-0.81696327) q[0];
sx q[0];
rz(-0.2635862) q[0];
x q[1];
rz(2.847371) q[2];
sx q[2];
rz(-0.18202848) q[2];
sx q[2];
rz(-2.1862669) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.1142769) q[1];
sx q[1];
rz(-1.7661347) q[1];
sx q[1];
rz(0.26135011) q[1];
rz(-pi) q[2];
rz(1.6513651) q[3];
sx q[3];
rz(-1.779883) q[3];
sx q[3];
rz(-0.63864708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.43188492) q[2];
sx q[2];
rz(-2.0543435) q[2];
sx q[2];
rz(-1.0464) q[2];
rz(-2.8228068) q[3];
sx q[3];
rz(-2.4304978) q[3];
sx q[3];
rz(-2.5939202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.043561291) q[0];
sx q[0];
rz(-1.8836319) q[0];
sx q[0];
rz(0.64796722) q[0];
rz(-1.3994392) q[1];
sx q[1];
rz(-1.49767) q[1];
sx q[1];
rz(-1.3290149) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9861222) q[0];
sx q[0];
rz(-2.8948445) q[0];
sx q[0];
rz(2.0339436) q[0];
rz(-1.7850661) q[2];
sx q[2];
rz(-1.8545215) q[2];
sx q[2];
rz(-1.2905215) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9860898) q[1];
sx q[1];
rz(-2.2515319) q[1];
sx q[1];
rz(0.40386856) q[1];
rz(-pi) q[2];
rz(0.59698128) q[3];
sx q[3];
rz(-2.5425445) q[3];
sx q[3];
rz(-2.0839276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7602188) q[2];
sx q[2];
rz(-1.2083283) q[2];
sx q[2];
rz(-1.5927429) q[2];
rz(-3.0717487) q[3];
sx q[3];
rz(-1.2589688) q[3];
sx q[3];
rz(0.86863345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0233362) q[0];
sx q[0];
rz(-1.9255487) q[0];
sx q[0];
rz(-2.1997531) q[0];
rz(-0.16009227) q[1];
sx q[1];
rz(-1.4804877) q[1];
sx q[1];
rz(-2.9494185) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1051467) q[0];
sx q[0];
rz(-0.22686401) q[0];
sx q[0];
rz(1.4099717) q[0];
rz(-pi) q[1];
x q[1];
rz(0.040080796) q[2];
sx q[2];
rz(-2.0922086) q[2];
sx q[2];
rz(1.5837976) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1035489) q[1];
sx q[1];
rz(-2.2749341) q[1];
sx q[1];
rz(2.0226993) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5302079) q[3];
sx q[3];
rz(-0.49792624) q[3];
sx q[3];
rz(0.17526173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.75752246) q[2];
sx q[2];
rz(-1.2382058) q[2];
sx q[2];
rz(2.7080217) q[2];
rz(-1.5571669) q[3];
sx q[3];
rz(-1.4945364) q[3];
sx q[3];
rz(-0.43237329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6381391) q[0];
sx q[0];
rz(-1.3766377) q[0];
sx q[0];
rz(1.7973416) q[0];
rz(1.9380219) q[1];
sx q[1];
rz(-1.1886339) q[1];
sx q[1];
rz(-1.3628091) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8840795) q[0];
sx q[0];
rz(-1.6005033) q[0];
sx q[0];
rz(-0.021935181) q[0];
rz(-1.7596471) q[2];
sx q[2];
rz(-0.7536234) q[2];
sx q[2];
rz(-0.058503956) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4688022) q[1];
sx q[1];
rz(-1.6957449) q[1];
sx q[1];
rz(0.23016696) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8239488) q[3];
sx q[3];
rz(-2.2322725) q[3];
sx q[3];
rz(-2.3931488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5726996) q[2];
sx q[2];
rz(-2.3685444) q[2];
sx q[2];
rz(-2.1577238) q[2];
rz(-0.24108663) q[3];
sx q[3];
rz(-0.9328931) q[3];
sx q[3];
rz(-2.4510395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47138658) q[0];
sx q[0];
rz(-0.79600483) q[0];
sx q[0];
rz(-0.27467003) q[0];
rz(1.7999016) q[1];
sx q[1];
rz(-2.5119753) q[1];
sx q[1];
rz(1.0460269) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24580851) q[0];
sx q[0];
rz(-1.2465451) q[0];
sx q[0];
rz(0.62600531) q[0];
x q[1];
rz(-1.3730896) q[2];
sx q[2];
rz(-1.4914163) q[2];
sx q[2];
rz(2.6153836) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.14635) q[1];
sx q[1];
rz(-2.1533794) q[1];
sx q[1];
rz(-1.6487153) q[1];
rz(1.0677797) q[3];
sx q[3];
rz(-1.5340337) q[3];
sx q[3];
rz(3.0779148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.86299738) q[2];
sx q[2];
rz(-1.7509165) q[2];
sx q[2];
rz(2.7421303) q[2];
rz(-2.5755889) q[3];
sx q[3];
rz(-2.1750906) q[3];
sx q[3];
rz(2.32617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28867662) q[0];
sx q[0];
rz(-1.4194019) q[0];
sx q[0];
rz(0.5823108) q[0];
rz(2.9583926) q[1];
sx q[1];
rz(-0.87545005) q[1];
sx q[1];
rz(-0.80642548) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69003478) q[0];
sx q[0];
rz(-1.4831044) q[0];
sx q[0];
rz(1.466485) q[0];
rz(-pi) q[1];
rz(0.33558947) q[2];
sx q[2];
rz(-1.8238471) q[2];
sx q[2];
rz(1.0831246) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.90102531) q[1];
sx q[1];
rz(-1.3677003) q[1];
sx q[1];
rz(-2.3939449) q[1];
x q[2];
rz(-0.34727283) q[3];
sx q[3];
rz(-1.4662305) q[3];
sx q[3];
rz(2.4109651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.23239423) q[2];
sx q[2];
rz(-2.4269203) q[2];
sx q[2];
rz(-2.3363028) q[2];
rz(0.14690873) q[3];
sx q[3];
rz(-2.3555136) q[3];
sx q[3];
rz(-0.73827353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2522226) q[0];
sx q[0];
rz(-1.7563553) q[0];
sx q[0];
rz(1.8808543) q[0];
rz(1.5421142) q[1];
sx q[1];
rz(-1.5442994) q[1];
sx q[1];
rz(-1.50179) q[1];
rz(-2.2940214) q[2];
sx q[2];
rz(-1.6568274) q[2];
sx q[2];
rz(-0.33100707) q[2];
rz(-0.64503786) q[3];
sx q[3];
rz(-1.8492263) q[3];
sx q[3];
rz(3.1254569) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
