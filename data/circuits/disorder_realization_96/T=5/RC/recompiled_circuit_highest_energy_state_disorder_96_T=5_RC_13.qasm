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
rz(0.86007287) q[0];
x q[1];
rz(-0.90801348) q[2];
sx q[2];
rz(-2.5129299) q[2];
sx q[2];
rz(-2.6952621) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9984879) q[1];
sx q[1];
rz(-1.2799026) q[1];
sx q[1];
rz(-2.8105633) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.82709978) q[3];
sx q[3];
rz(-1.4129644) q[3];
sx q[3];
rz(-2.1438897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.277694) q[2];
sx q[2];
rz(-0.70964491) q[2];
sx q[2];
rz(-0.42897439) q[2];
rz(-0.046836827) q[3];
sx q[3];
rz(-2.7497079) q[3];
sx q[3];
rz(-0.61280167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.866339) q[0];
sx q[0];
rz(-2.1493122) q[0];
sx q[0];
rz(0.66542768) q[0];
rz(-2.0003419) q[1];
sx q[1];
rz(-1.7612061) q[1];
sx q[1];
rz(-2.4990987) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25099558) q[0];
sx q[0];
rz(-0.83775508) q[0];
sx q[0];
rz(1.6392073) q[0];
rz(-1.8694473) q[2];
sx q[2];
rz(-2.2437216) q[2];
sx q[2];
rz(-2.330614) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4027109) q[1];
sx q[1];
rz(-2.4428574) q[1];
sx q[1];
rz(-0.043234392) q[1];
x q[2];
rz(-1.2824545) q[3];
sx q[3];
rz(-0.55984646) q[3];
sx q[3];
rz(2.0324872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.82681727) q[2];
sx q[2];
rz(-0.78709698) q[2];
sx q[2];
rz(-0.55244201) q[2];
rz(-1.6636482) q[3];
sx q[3];
rz(-1.843957) q[3];
sx q[3];
rz(-0.67599952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8058572) q[0];
sx q[0];
rz(-2.4215846) q[0];
sx q[0];
rz(0.82255256) q[0];
rz(0.81126732) q[1];
sx q[1];
rz(-2.8370116) q[1];
sx q[1];
rz(1.6792345) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84483355) q[0];
sx q[0];
rz(-2.5243596) q[0];
sx q[0];
rz(0.69302859) q[0];
rz(-pi) q[1];
x q[1];
rz(0.16872318) q[2];
sx q[2];
rz(-0.41831917) q[2];
sx q[2];
rz(-0.63101286) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9595527) q[1];
sx q[1];
rz(-2.5904852) q[1];
sx q[1];
rz(0.98513794) q[1];
rz(-pi) q[2];
rz(2.3554166) q[3];
sx q[3];
rz(-2.7179931) q[3];
sx q[3];
rz(2.9860403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2751969) q[2];
sx q[2];
rz(-0.82922816) q[2];
sx q[2];
rz(-2.5466476) q[2];
rz(0.40469894) q[3];
sx q[3];
rz(-1.1070822) q[3];
sx q[3];
rz(-1.2987202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4054656) q[0];
sx q[0];
rz(-1.8647702) q[0];
sx q[0];
rz(2.8986616) q[0];
rz(2.2566707) q[1];
sx q[1];
rz(-2.4156069) q[1];
sx q[1];
rz(-3.1387709) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0962727) q[0];
sx q[0];
rz(-2.3191873) q[0];
sx q[0];
rz(-2.1069645) q[0];
rz(-pi) q[1];
rz(-2.6153938) q[2];
sx q[2];
rz(-2.1251273) q[2];
sx q[2];
rz(0.89579158) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5127848) q[1];
sx q[1];
rz(-1.50042) q[1];
sx q[1];
rz(-2.2712399) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.76169984) q[3];
sx q[3];
rz(-2.3683067) q[3];
sx q[3];
rz(2.9536216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.60488492) q[2];
sx q[2];
rz(-1.1515836) q[2];
sx q[2];
rz(2.4082157) q[2];
rz(-2.4865161) q[3];
sx q[3];
rz(-2.8337182) q[3];
sx q[3];
rz(-2.5797599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4578399) q[0];
sx q[0];
rz(-2.8003052) q[0];
sx q[0];
rz(-0.14232464) q[0];
rz(1.3617474) q[1];
sx q[1];
rz(-2.7889377) q[1];
sx q[1];
rz(0.49837643) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4623584) q[0];
sx q[0];
rz(-0.7158044) q[0];
sx q[0];
rz(1.6648034) q[0];
rz(-pi) q[1];
rz(-0.054008897) q[2];
sx q[2];
rz(-2.3819792) q[2];
sx q[2];
rz(-0.97689607) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0432669) q[1];
sx q[1];
rz(-1.2778751) q[1];
sx q[1];
rz(-0.75324159) q[1];
rz(-0.77985127) q[3];
sx q[3];
rz(-1.3289467) q[3];
sx q[3];
rz(-2.4132501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.132823) q[2];
sx q[2];
rz(-1.2558197) q[2];
sx q[2];
rz(-0.69457561) q[2];
rz(2.3092367) q[3];
sx q[3];
rz(-1.1135626) q[3];
sx q[3];
rz(-2.8913403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.71448) q[0];
sx q[0];
rz(-2.6321593) q[0];
sx q[0];
rz(-0.66202778) q[0];
rz(-1.9561249) q[1];
sx q[1];
rz(-1.0292425) q[1];
sx q[1];
rz(-1.1659291) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53138779) q[0];
sx q[0];
rz(-0.90739668) q[0];
sx q[0];
rz(-1.8276455) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4170961) q[2];
sx q[2];
rz(-0.96818334) q[2];
sx q[2];
rz(3.1069047) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.59938972) q[1];
sx q[1];
rz(-1.521061) q[1];
sx q[1];
rz(-1.6433783) q[1];
rz(-pi) q[2];
rz(-2.1816129) q[3];
sx q[3];
rz(-2.0997542) q[3];
sx q[3];
rz(-0.66416603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7621496) q[2];
sx q[2];
rz(-1.9337485) q[2];
sx q[2];
rz(2.4318648) q[2];
rz(-0.48192561) q[3];
sx q[3];
rz(-0.47526264) q[3];
sx q[3];
rz(0.019439241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-2.6946805) q[1];
sx q[1];
rz(-2.2692197) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.746856) q[0];
sx q[0];
rz(-2.7867705) q[0];
sx q[0];
rz(1.0635934) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.21193223) q[2];
sx q[2];
rz(-0.733558) q[2];
sx q[2];
rz(-0.19679815) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7561314) q[1];
sx q[1];
rz(-1.676899) q[1];
sx q[1];
rz(-0.030563449) q[1];
rz(-pi) q[2];
rz(-0.024691651) q[3];
sx q[3];
rz(-0.91626142) q[3];
sx q[3];
rz(-0.86658044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.81277728) q[2];
sx q[2];
rz(-0.8605364) q[2];
sx q[2];
rz(-2.013773) q[2];
rz(0.40337107) q[3];
sx q[3];
rz(-1.4453459) q[3];
sx q[3];
rz(-2.4283714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9047852) q[0];
sx q[0];
rz(-1.1431563) q[0];
sx q[0];
rz(1.2971725) q[0];
rz(-2.6203268) q[1];
sx q[1];
rz(-1.8788985) q[1];
sx q[1];
rz(3.0788132) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.990898) q[0];
sx q[0];
rz(-2.6932062) q[0];
sx q[0];
rz(-1.4175116) q[0];
x q[1];
rz(-0.56030484) q[2];
sx q[2];
rz(-0.41839278) q[2];
sx q[2];
rz(-0.095730893) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.023397327) q[1];
sx q[1];
rz(-0.28701008) q[1];
sx q[1];
rz(1.0402518) q[1];
rz(2.217123) q[3];
sx q[3];
rz(-2.0144343) q[3];
sx q[3];
rz(-0.2547338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2433743) q[2];
sx q[2];
rz(-0.90460193) q[2];
sx q[2];
rz(0.96837366) q[2];
rz(-2.6280256) q[3];
sx q[3];
rz(-1.7889203) q[3];
sx q[3];
rz(0.33094049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6043855) q[0];
sx q[0];
rz(-2.4810915) q[0];
sx q[0];
rz(-0.78395098) q[0];
rz(-1.2046658) q[1];
sx q[1];
rz(-2.7684863) q[1];
sx q[1];
rz(1.7663667) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9722906) q[0];
sx q[0];
rz(-1.307104) q[0];
sx q[0];
rz(1.792981) q[0];
rz(1.398574) q[2];
sx q[2];
rz(-0.15744124) q[2];
sx q[2];
rz(-2.1724608) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7837401) q[1];
sx q[1];
rz(-2.2688818) q[1];
sx q[1];
rz(-2.3925406) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8834388) q[3];
sx q[3];
rz(-1.7400768) q[3];
sx q[3];
rz(-1.1723684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2969926) q[2];
sx q[2];
rz(-0.30868369) q[2];
sx q[2];
rz(0.44573477) q[2];
rz(-0.60278696) q[3];
sx q[3];
rz(-1.2647537) q[3];
sx q[3];
rz(-2.2980237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26850253) q[0];
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
rz(2.1969057) q[0];
sx q[0];
rz(-1.9160998) q[0];
sx q[0];
rz(-1.7717591) q[0];
rz(-1.3648562) q[2];
sx q[2];
rz(-0.32216926) q[2];
sx q[2];
rz(-1.3141156) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5752856) q[1];
sx q[1];
rz(-3.008209) q[1];
sx q[1];
rz(-0.092194362) q[1];
x q[2];
rz(1.6053725) q[3];
sx q[3];
rz(-2.1994221) q[3];
sx q[3];
rz(2.2104267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.21620096) q[2];
sx q[2];
rz(-0.4500469) q[2];
sx q[2];
rz(-1.0518543) q[2];
rz(0.22107407) q[3];
sx q[3];
rz(-1.8408006) q[3];
sx q[3];
rz(-0.93824798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5476407) q[0];
sx q[0];
rz(-1.5008391) q[0];
sx q[0];
rz(0.048901625) q[0];
rz(2.7652057) q[1];
sx q[1];
rz(-2.2793437) q[1];
sx q[1];
rz(1.9752165) q[1];
rz(1.5390039) q[2];
sx q[2];
rz(-1.0009264) q[2];
sx q[2];
rz(1.4957503) q[2];
rz(-2.4276499) q[3];
sx q[3];
rz(-1.1957914) q[3];
sx q[3];
rz(1.316432) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
