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
rz(0.30492914) q[0];
sx q[0];
rz(-0.066216901) q[0];
sx q[0];
rz(-2.9856292) q[0];
rz(-1.9744385) q[1];
sx q[1];
rz(-0.65216291) q[1];
sx q[1];
rz(-2.6551533) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4733361) q[0];
sx q[0];
rz(-1.0205246) q[0];
sx q[0];
rz(-2.2815198) q[0];
rz(-pi) q[1];
x q[1];
rz(0.42066853) q[2];
sx q[2];
rz(-1.0887869) q[2];
sx q[2];
rz(1.2139621) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.26775441) q[1];
sx q[1];
rz(-2.7044786) q[1];
sx q[1];
rz(2.3972191) q[1];
rz(-pi) q[2];
rz(0.21297314) q[3];
sx q[3];
rz(-0.83847441) q[3];
sx q[3];
rz(0.71668437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8638986) q[2];
sx q[2];
rz(-2.4319477) q[2];
sx q[2];
rz(-2.7126183) q[2];
rz(-3.0947558) q[3];
sx q[3];
rz(-2.7497079) q[3];
sx q[3];
rz(-2.528791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.2752537) q[0];
sx q[0];
rz(-0.99228042) q[0];
sx q[0];
rz(-0.66542768) q[0];
rz(2.0003419) q[1];
sx q[1];
rz(-1.3803866) q[1];
sx q[1];
rz(0.64249396) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35303799) q[0];
sx q[0];
rz(-0.73563535) q[0];
sx q[0];
rz(-3.0658196) q[0];
rz(2.4464954) q[2];
sx q[2];
rz(-1.8029658) q[2];
sx q[2];
rz(-2.1921981) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.68243945) q[1];
sx q[1];
rz(-0.87284589) q[1];
sx q[1];
rz(-1.5345011) q[1];
rz(-pi) q[2];
x q[2];
rz(2.965224) q[3];
sx q[3];
rz(-2.1049728) q[3];
sx q[3];
rz(1.4458223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3147754) q[2];
sx q[2];
rz(-2.3544957) q[2];
sx q[2];
rz(-2.5891506) q[2];
rz(1.6636482) q[3];
sx q[3];
rz(-1.2976357) q[3];
sx q[3];
rz(-0.67599952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8058572) q[0];
sx q[0];
rz(-0.72000802) q[0];
sx q[0];
rz(-2.3190401) q[0];
rz(-2.3303253) q[1];
sx q[1];
rz(-0.30458105) q[1];
sx q[1];
rz(1.4623581) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84483355) q[0];
sx q[0];
rz(-0.61723304) q[0];
sx q[0];
rz(0.69302859) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.64531) q[2];
sx q[2];
rz(-1.9828116) q[2];
sx q[2];
rz(-2.3262784) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2671859) q[1];
sx q[1];
rz(-1.864434) q[1];
sx q[1];
rz(1.0974357) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8331746) q[3];
sx q[3];
rz(-1.8659411) q[3];
sx q[3];
rz(-2.4663188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2751969) q[2];
sx q[2];
rz(-2.3123645) q[2];
sx q[2];
rz(-0.59494507) q[2];
rz(0.40469894) q[3];
sx q[3];
rz(-1.1070822) q[3];
sx q[3];
rz(1.8428724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73612708) q[0];
sx q[0];
rz(-1.2768224) q[0];
sx q[0];
rz(-2.8986616) q[0];
rz(-2.2566707) q[1];
sx q[1];
rz(-2.4156069) q[1];
sx q[1];
rz(3.1387709) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1412239) q[0];
sx q[0];
rz(-1.9544811) q[0];
sx q[0];
rz(0.82392086) q[0];
x q[1];
rz(-2.1921333) q[2];
sx q[2];
rz(-2.0120125) q[2];
sx q[2];
rz(-0.37829933) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2828212) q[1];
sx q[1];
rz(-2.438218) q[1];
sx q[1];
rz(-1.4618631) q[1];
rz(-pi) q[2];
rz(-2.1636064) q[3];
sx q[3];
rz(-2.1007256) q[3];
sx q[3];
rz(-2.0265614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.60488492) q[2];
sx q[2];
rz(-1.9900091) q[2];
sx q[2];
rz(-2.4082157) q[2];
rz(2.4865161) q[3];
sx q[3];
rz(-2.8337182) q[3];
sx q[3];
rz(2.5797599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4578399) q[0];
sx q[0];
rz(-2.8003052) q[0];
sx q[0];
rz(2.999268) q[0];
rz(-1.3617474) q[1];
sx q[1];
rz(-2.7889377) q[1];
sx q[1];
rz(-0.49837643) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.179006) q[0];
sx q[0];
rz(-1.5091584) q[0];
sx q[0];
rz(-2.2844102) q[0];
x q[1];
rz(1.519573) q[2];
sx q[2];
rz(-2.3290259) q[2];
sx q[2];
rz(1.0513154) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0432669) q[1];
sx q[1];
rz(-1.2778751) q[1];
sx q[1];
rz(0.75324159) q[1];
rz(-pi) q[2];
x q[2];
rz(0.77985127) q[3];
sx q[3];
rz(-1.3289467) q[3];
sx q[3];
rz(-0.72834258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0087697) q[2];
sx q[2];
rz(-1.2558197) q[2];
sx q[2];
rz(-2.447017) q[2];
rz(0.83235598) q[3];
sx q[3];
rz(-1.1135626) q[3];
sx q[3];
rz(-0.25025234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-2.71448) q[0];
sx q[0];
rz(-2.6321593) q[0];
sx q[0];
rz(0.66202778) q[0];
rz(1.9561249) q[1];
sx q[1];
rz(-1.0292425) q[1];
sx q[1];
rz(-1.9756636) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93454516) q[0];
sx q[0];
rz(-0.70431346) q[0];
sx q[0];
rz(0.31426104) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.60816216) q[2];
sx q[2];
rz(-1.4443304) q[2];
sx q[2];
rz(-1.6236931) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5422029) q[1];
sx q[1];
rz(-1.521061) q[1];
sx q[1];
rz(-1.6433783) q[1];
rz(-2.3656769) q[3];
sx q[3];
rz(-2.3563623) q[3];
sx q[3];
rz(-1.6104346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7621496) q[2];
sx q[2];
rz(-1.2078441) q[2];
sx q[2];
rz(2.4318648) q[2];
rz(2.659667) q[3];
sx q[3];
rz(-0.47526264) q[3];
sx q[3];
rz(0.019439241) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67126453) q[0];
sx q[0];
rz(-2.1605991) q[0];
sx q[0];
rz(1.5671267) q[0];
rz(-1.4944685) q[1];
sx q[1];
rz(-0.44691214) q[1];
sx q[1];
rz(2.2692197) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4851661) q[0];
sx q[0];
rz(-1.4012294) q[0];
sx q[0];
rz(1.2575968) q[0];
x q[1];
rz(-1.7581697) q[2];
sx q[2];
rz(-2.2843666) q[2];
sx q[2];
rz(-3.0564412) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.10441594) q[1];
sx q[1];
rz(-3.0311916) q[1];
sx q[1];
rz(1.2913741) q[1];
rz(3.116901) q[3];
sx q[3];
rz(-0.91626142) q[3];
sx q[3];
rz(-0.86658044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3288154) q[2];
sx q[2];
rz(-0.8605364) q[2];
sx q[2];
rz(1.1278197) q[2];
rz(-2.7382216) q[3];
sx q[3];
rz(-1.4453459) q[3];
sx q[3];
rz(0.71322125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23680747) q[0];
sx q[0];
rz(-1.1431563) q[0];
sx q[0];
rz(1.2971725) q[0];
rz(-0.5212658) q[1];
sx q[1];
rz(-1.2626941) q[1];
sx q[1];
rz(-0.062779471) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2817665) q[0];
sx q[0];
rz(-1.6370356) q[0];
sx q[0];
rz(-1.1270004) q[0];
rz(1.3387483) q[2];
sx q[2];
rz(-1.9221483) q[2];
sx q[2];
rz(-2.6356489) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.52555841) q[1];
sx q[1];
rz(-1.3241321) q[1];
sx q[1];
rz(2.9933369) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.217123) q[3];
sx q[3];
rz(-1.1271584) q[3];
sx q[3];
rz(2.8868589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2433743) q[2];
sx q[2];
rz(-2.2369907) q[2];
sx q[2];
rz(0.96837366) q[2];
rz(0.51356703) q[3];
sx q[3];
rz(-1.3526724) q[3];
sx q[3];
rz(-0.33094049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53720713) q[0];
sx q[0];
rz(-0.66050118) q[0];
sx q[0];
rz(2.3576417) q[0];
rz(-1.9369269) q[1];
sx q[1];
rz(-2.7684863) q[1];
sx q[1];
rz(1.3752259) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88343745) q[0];
sx q[0];
rz(-2.7984507) q[0];
sx q[0];
rz(-0.68455066) q[0];
rz(-pi) q[1];
rz(1.398574) q[2];
sx q[2];
rz(-0.15744124) q[2];
sx q[2];
rz(0.96913183) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7837401) q[1];
sx q[1];
rz(-2.2688818) q[1];
sx q[1];
rz(-0.74905209) q[1];
x q[2];
rz(-2.9638644) q[3];
sx q[3];
rz(-1.2627708) q[3];
sx q[3];
rz(2.7975688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.84460008) q[2];
sx q[2];
rz(-2.832909) q[2];
sx q[2];
rz(-2.6958579) q[2];
rz(0.60278696) q[3];
sx q[3];
rz(-1.2647537) q[3];
sx q[3];
rz(2.2980237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26850253) q[0];
sx q[0];
rz(-0.52274811) q[0];
sx q[0];
rz(-2.6173746) q[0];
rz(1.9408608) q[1];
sx q[1];
rz(-1.8914765) q[1];
sx q[1];
rz(1.9573617) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1969057) q[0];
sx q[0];
rz(-1.2254929) q[0];
sx q[0];
rz(1.7717591) q[0];
rz(-pi) q[1];
rz(3.0734407) q[2];
sx q[2];
rz(-1.8859204) q[2];
sx q[2];
rz(1.5308876) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5752856) q[1];
sx q[1];
rz(-3.008209) q[1];
sx q[1];
rz(-3.0493983) q[1];
rz(2.5126825) q[3];
sx q[3];
rz(-1.5428318) q[3];
sx q[3];
rz(2.4816251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9253917) q[2];
sx q[2];
rz(-0.4500469) q[2];
sx q[2];
rz(1.0518543) q[2];
rz(2.9205186) q[3];
sx q[3];
rz(-1.8408006) q[3];
sx q[3];
rz(-2.2033447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5939519) q[0];
sx q[0];
rz(-1.6407536) q[0];
sx q[0];
rz(-3.092691) q[0];
rz(-0.37638695) q[1];
sx q[1];
rz(-2.2793437) q[1];
sx q[1];
rz(1.9752165) q[1];
rz(-0.57009956) q[2];
sx q[2];
rz(-1.5975633) q[2];
sx q[2];
rz(-0.057889197) q[2];
rz(2.4276499) q[3];
sx q[3];
rz(-1.9458013) q[3];
sx q[3];
rz(-1.8251606) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
