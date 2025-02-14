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
rz(-0.90801348) q[2];
sx q[2];
rz(-2.5129299) q[2];
sx q[2];
rz(-2.6952621) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9984879) q[1];
sx q[1];
rz(-1.8616901) q[1];
sx q[1];
rz(2.8105633) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8016889) q[3];
sx q[3];
rz(-2.3844815) q[3];
sx q[3];
rz(-2.7377303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.277694) q[2];
sx q[2];
rz(-0.70964491) q[2];
sx q[2];
rz(-2.7126183) q[2];
rz(0.046836827) q[3];
sx q[3];
rz(-0.39188477) q[3];
sx q[3];
rz(-0.61280167) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
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
rz(-1.7612061) q[1];
sx q[1];
rz(-0.64249396) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35303799) q[0];
sx q[0];
rz(-2.4059573) q[0];
sx q[0];
rz(0.075773043) q[0];
rz(-2.4464954) q[2];
sx q[2];
rz(-1.8029658) q[2];
sx q[2];
rz(2.1921981) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2765668) q[1];
sx q[1];
rz(-1.5429909) q[1];
sx q[1];
rz(2.4433178) q[1];
rz(2.965224) q[3];
sx q[3];
rz(-1.0366199) q[3];
sx q[3];
rz(-1.4458223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3147754) q[2];
sx q[2];
rz(-0.78709698) q[2];
sx q[2];
rz(-2.5891506) q[2];
rz(-1.6636482) q[3];
sx q[3];
rz(-1.843957) q[3];
sx q[3];
rz(2.4655931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33573547) q[0];
sx q[0];
rz(-2.4215846) q[0];
sx q[0];
rz(0.82255256) q[0];
rz(-2.3303253) q[1];
sx q[1];
rz(-0.30458105) q[1];
sx q[1];
rz(-1.6792345) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8203634) q[0];
sx q[0];
rz(-1.9495533) q[0];
sx q[0];
rz(2.6418153) q[0];
rz(-pi) q[1];
rz(1.4962827) q[2];
sx q[2];
rz(-1.1587811) q[2];
sx q[2];
rz(2.3262784) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9595527) q[1];
sx q[1];
rz(-2.5904852) q[1];
sx q[1];
rz(0.98513794) q[1];
rz(-pi) q[2];
rz(-2.3554166) q[3];
sx q[3];
rz(-2.7179931) q[3];
sx q[3];
rz(-2.9860403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.86639577) q[2];
sx q[2];
rz(-0.82922816) q[2];
sx q[2];
rz(-0.59494507) q[2];
rz(2.7368937) q[3];
sx q[3];
rz(-2.0345104) q[3];
sx q[3];
rz(-1.2987202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73612708) q[0];
sx q[0];
rz(-1.2768224) q[0];
sx q[0];
rz(2.8986616) q[0];
rz(2.2566707) q[1];
sx q[1];
rz(-0.72598571) q[1];
sx q[1];
rz(3.1387709) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76317518) q[0];
sx q[0];
rz(-2.2522914) q[0];
sx q[0];
rz(-0.50294089) q[0];
x q[1];
rz(-2.6153938) q[2];
sx q[2];
rz(-2.1251273) q[2];
sx q[2];
rz(-2.2458011) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1403919) q[1];
sx q[1];
rz(-0.87243783) q[1];
sx q[1];
rz(0.091940885) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3798928) q[3];
sx q[3];
rz(-2.3683067) q[3];
sx q[3];
rz(-0.1879711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5367077) q[2];
sx q[2];
rz(-1.1515836) q[2];
sx q[2];
rz(2.4082157) q[2];
rz(-2.4865161) q[3];
sx q[3];
rz(-2.8337182) q[3];
sx q[3];
rz(0.56183279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68375278) q[0];
sx q[0];
rz(-0.34128749) q[0];
sx q[0];
rz(-0.14232464) q[0];
rz(1.3617474) q[1];
sx q[1];
rz(-2.7889377) q[1];
sx q[1];
rz(-2.6432162) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67923421) q[0];
sx q[0];
rz(-0.7158044) q[0];
sx q[0];
rz(-1.4767893) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3827078) q[2];
sx q[2];
rz(-1.6079796) q[2];
sx q[2];
rz(-2.5868724) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.26311647) q[1];
sx q[1];
rz(-0.85678393) q[1];
sx q[1];
rz(-1.962838) q[1];
rz(-pi) q[2];
rz(0.77985127) q[3];
sx q[3];
rz(-1.3289467) q[3];
sx q[3];
rz(-0.72834258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.132823) q[2];
sx q[2];
rz(-1.2558197) q[2];
sx q[2];
rz(-2.447017) q[2];
rz(-2.3092367) q[3];
sx q[3];
rz(-1.1135626) q[3];
sx q[3];
rz(2.8913403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
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
rz(0.4271127) q[0];
sx q[0];
rz(-2.6321593) q[0];
sx q[0];
rz(0.66202778) q[0];
rz(-1.9561249) q[1];
sx q[1];
rz(-2.1123501) q[1];
sx q[1];
rz(1.1659291) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93454516) q[0];
sx q[0];
rz(-0.70431346) q[0];
sx q[0];
rz(0.31426104) q[0];
rz(-0.60816216) q[2];
sx q[2];
rz(-1.6972622) q[2];
sx q[2];
rz(1.6236931) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1665713) q[1];
sx q[1];
rz(-1.6432883) q[1];
sx q[1];
rz(-0.049866421) q[1];
rz(-0.61975584) q[3];
sx q[3];
rz(-1.0527851) q[3];
sx q[3];
rz(-2.5745846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.379443) q[2];
sx q[2];
rz(-1.9337485) q[2];
sx q[2];
rz(-0.70972788) q[2];
rz(2.659667) q[3];
sx q[3];
rz(-2.66633) q[3];
sx q[3];
rz(-0.019439241) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4703281) q[0];
sx q[0];
rz(-0.98099357) q[0];
sx q[0];
rz(-1.5671267) q[0];
rz(1.6471242) q[1];
sx q[1];
rz(-0.44691214) q[1];
sx q[1];
rz(-0.87237298) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4851661) q[0];
sx q[0];
rz(-1.7403633) q[0];
sx q[0];
rz(1.2575968) q[0];
rz(1.383423) q[2];
sx q[2];
rz(-0.85722605) q[2];
sx q[2];
rz(3.0564412) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7561314) q[1];
sx q[1];
rz(-1.4646936) q[1];
sx q[1];
rz(3.1110292) q[1];
x q[2];
rz(0.91611417) q[3];
sx q[3];
rz(-1.5512084) q[3];
sx q[3];
rz(2.4223428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3288154) q[2];
sx q[2];
rz(-0.8605364) q[2];
sx q[2];
rz(-1.1278197) q[2];
rz(2.7382216) q[3];
sx q[3];
rz(-1.4453459) q[3];
sx q[3];
rz(-0.71322125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9047852) q[0];
sx q[0];
rz(-1.9984364) q[0];
sx q[0];
rz(-1.2971725) q[0];
rz(0.5212658) q[1];
sx q[1];
rz(-1.2626941) q[1];
sx q[1];
rz(0.062779471) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8598261) q[0];
sx q[0];
rz(-1.504557) q[0];
sx q[0];
rz(1.1270004) q[0];
x q[1];
rz(0.56030484) q[2];
sx q[2];
rz(-0.41839278) q[2];
sx q[2];
rz(0.095730893) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0599036) q[1];
sx q[1];
rz(-1.7145331) q[1];
sx q[1];
rz(-1.3215076) q[1];
x q[2];
rz(-2.6046329) q[3];
sx q[3];
rz(-0.99565047) q[3];
sx q[3];
rz(1.6292269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2433743) q[2];
sx q[2];
rz(-0.90460193) q[2];
sx q[2];
rz(-2.173219) q[2];
rz(-2.6280256) q[3];
sx q[3];
rz(-1.3526724) q[3];
sx q[3];
rz(-0.33094049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53720713) q[0];
sx q[0];
rz(-2.4810915) q[0];
sx q[0];
rz(-0.78395098) q[0];
rz(1.2046658) q[1];
sx q[1];
rz(-0.37310633) q[1];
sx q[1];
rz(-1.3752259) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7989144) q[0];
sx q[0];
rz(-1.7851789) q[0];
sx q[0];
rz(0.27002295) q[0];
rz(-pi) q[1];
rz(-3.1143931) q[2];
sx q[2];
rz(-1.7258894) q[2];
sx q[2];
rz(1.9981249) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3578525) q[1];
sx q[1];
rz(-0.87271089) q[1];
sx q[1];
rz(0.74905209) q[1];
rz(-1.0635942) q[3];
sx q[3];
rz(-0.35420277) q[3];
sx q[3];
rz(0.87888792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.84460008) q[2];
sx q[2];
rz(-0.30868369) q[2];
sx q[2];
rz(-2.6958579) q[2];
rz(-0.60278696) q[3];
sx q[3];
rz(-1.2647537) q[3];
sx q[3];
rz(-2.2980237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8730901) q[0];
sx q[0];
rz(-0.52274811) q[0];
sx q[0];
rz(0.52421808) q[0];
rz(-1.2007319) q[1];
sx q[1];
rz(-1.8914765) q[1];
sx q[1];
rz(-1.1842309) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55726526) q[0];
sx q[0];
rz(-1.7597489) q[0];
sx q[0];
rz(0.35182955) q[0];
rz(1.7767364) q[2];
sx q[2];
rz(-0.32216926) q[2];
sx q[2];
rz(1.8274771) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5752856) q[1];
sx q[1];
rz(-3.008209) q[1];
sx q[1];
rz(3.0493983) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0940787) q[3];
sx q[3];
rz(-0.62944747) q[3];
sx q[3];
rz(2.269182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.21620096) q[2];
sx q[2];
rz(-0.4500469) q[2];
sx q[2];
rz(2.0897384) q[2];
rz(-2.9205186) q[3];
sx q[3];
rz(-1.3007921) q[3];
sx q[3];
rz(-2.2033447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.5939519) q[0];
sx q[0];
rz(-1.5008391) q[0];
sx q[0];
rz(0.048901625) q[0];
rz(2.7652057) q[1];
sx q[1];
rz(-2.2793437) q[1];
sx q[1];
rz(1.9752165) q[1];
rz(-1.5390039) q[2];
sx q[2];
rz(-2.1406662) q[2];
sx q[2];
rz(-1.6458423) q[2];
rz(-2.6003414) q[3];
sx q[3];
rz(-0.79081906) q[3];
sx q[3];
rz(-2.9959903) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
