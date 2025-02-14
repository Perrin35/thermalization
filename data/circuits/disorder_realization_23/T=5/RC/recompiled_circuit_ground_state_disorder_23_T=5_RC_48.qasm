OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.57103676) q[0];
sx q[0];
rz(-2.4693627) q[0];
sx q[0];
rz(1.6663405) q[0];
rz(-2.1885459) q[1];
sx q[1];
rz(-0.16521984) q[1];
sx q[1];
rz(1.8728949) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27086285) q[0];
sx q[0];
rz(-2.6755973) q[0];
sx q[0];
rz(-2.6129338) q[0];
rz(-pi) q[1];
rz(-2.181594) q[2];
sx q[2];
rz(-2.3899601) q[2];
sx q[2];
rz(-0.014873504) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4339612) q[1];
sx q[1];
rz(-2.5565256) q[1];
sx q[1];
rz(-1.332704) q[1];
rz(-pi) q[2];
rz(1.3399063) q[3];
sx q[3];
rz(-0.82726631) q[3];
sx q[3];
rz(-0.32318599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1573726) q[2];
sx q[2];
rz(-0.34175384) q[2];
sx q[2];
rz(-0.79322195) q[2];
rz(-0.67304099) q[3];
sx q[3];
rz(-2.0561736) q[3];
sx q[3];
rz(-1.2696772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15758812) q[0];
sx q[0];
rz(-2.3389811) q[0];
sx q[0];
rz(0.60930914) q[0];
rz(-2.9318103) q[1];
sx q[1];
rz(-0.79133004) q[1];
sx q[1];
rz(0.9598859) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4835371) q[0];
sx q[0];
rz(-1.5262881) q[0];
sx q[0];
rz(-3.0625212) q[0];
rz(3.0716607) q[2];
sx q[2];
rz(-2.0581823) q[2];
sx q[2];
rz(-0.63475459) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.4757735) q[1];
sx q[1];
rz(-0.83516652) q[1];
sx q[1];
rz(-3.0569265) q[1];
rz(1.0632564) q[3];
sx q[3];
rz(-1.1843269) q[3];
sx q[3];
rz(-1.8563351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3557055) q[2];
sx q[2];
rz(-0.83928078) q[2];
sx q[2];
rz(2.6500224) q[2];
rz(-3.135904) q[3];
sx q[3];
rz(-2.0523043) q[3];
sx q[3];
rz(2.6277241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(3.046535) q[0];
sx q[0];
rz(-0.53545606) q[0];
sx q[0];
rz(0.74263483) q[0];
rz(-0.62478089) q[1];
sx q[1];
rz(-2.2014328) q[1];
sx q[1];
rz(-1.0661485) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4542592) q[0];
sx q[0];
rz(-1.3479184) q[0];
sx q[0];
rz(2.9965626) q[0];
x q[1];
rz(-0.19045388) q[2];
sx q[2];
rz(-1.2553658) q[2];
sx q[2];
rz(-0.73376943) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6934476) q[1];
sx q[1];
rz(-1.168201) q[1];
sx q[1];
rz(-3.0047396) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9227299) q[3];
sx q[3];
rz(-1.0893679) q[3];
sx q[3];
rz(0.66314135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8891958) q[2];
sx q[2];
rz(-2.9434581) q[2];
sx q[2];
rz(-2.5890217) q[2];
rz(0.50734723) q[3];
sx q[3];
rz(-1.0891958) q[3];
sx q[3];
rz(2.5721917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4813389) q[0];
sx q[0];
rz(-2.5642671) q[0];
sx q[0];
rz(-1.2465771) q[0];
rz(-1.735911) q[1];
sx q[1];
rz(-2.3697) q[1];
sx q[1];
rz(2.4235922) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8193977) q[0];
sx q[0];
rz(-2.2912663) q[0];
sx q[0];
rz(0.87089355) q[0];
rz(-pi) q[1];
rz(1.202232) q[2];
sx q[2];
rz(-2.5466726) q[2];
sx q[2];
rz(1.0234317) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.64252428) q[1];
sx q[1];
rz(-2.32076) q[1];
sx q[1];
rz(1.7700511) q[1];
rz(-pi) q[2];
rz(2.4627732) q[3];
sx q[3];
rz(-0.20256895) q[3];
sx q[3];
rz(2.5819419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8743073) q[2];
sx q[2];
rz(-2.3882046) q[2];
sx q[2];
rz(2.4370952) q[2];
rz(2.7590175) q[3];
sx q[3];
rz(-1.4380598) q[3];
sx q[3];
rz(-0.85324919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.48069561) q[0];
sx q[0];
rz(-0.040204164) q[0];
sx q[0];
rz(2.8327292) q[0];
rz(0.94912306) q[1];
sx q[1];
rz(-0.32052761) q[1];
sx q[1];
rz(2.5905051) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8723485) q[0];
sx q[0];
rz(-1.9595032) q[0];
sx q[0];
rz(1.7354911) q[0];
x q[1];
rz(-1.6597719) q[2];
sx q[2];
rz(-0.65976652) q[2];
sx q[2];
rz(-2.8967186) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.84173735) q[1];
sx q[1];
rz(-1.7256712) q[1];
sx q[1];
rz(-2.9938404) q[1];
rz(-2.9948009) q[3];
sx q[3];
rz(-0.78189584) q[3];
sx q[3];
rz(1.1679389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.71517313) q[2];
sx q[2];
rz(-2.3931563) q[2];
sx q[2];
rz(0.61881649) q[2];
rz(2.5185781) q[3];
sx q[3];
rz(-2.2996733) q[3];
sx q[3];
rz(2.714341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1022559) q[0];
sx q[0];
rz(-0.66348851) q[0];
sx q[0];
rz(-2.6797507) q[0];
rz(-2.6448008) q[1];
sx q[1];
rz(-1.8180314) q[1];
sx q[1];
rz(-1.7740446) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38765466) q[0];
sx q[0];
rz(-2.4552422) q[0];
sx q[0];
rz(1.4453056) q[0];
rz(-pi) q[1];
rz(-1.073764) q[2];
sx q[2];
rz(-1.4364527) q[2];
sx q[2];
rz(1.4984992) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.16713472) q[1];
sx q[1];
rz(-0.30213812) q[1];
sx q[1];
rz(-1.664851) q[1];
rz(-pi) q[2];
rz(-2.1928093) q[3];
sx q[3];
rz(-1.439073) q[3];
sx q[3];
rz(-1.6831116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0355012) q[2];
sx q[2];
rz(-2.0095299) q[2];
sx q[2];
rz(0.79037017) q[2];
rz(2.5370989) q[3];
sx q[3];
rz(-1.0249745) q[3];
sx q[3];
rz(-2.7214971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4496434) q[0];
sx q[0];
rz(-2.255891) q[0];
sx q[0];
rz(2.7272136) q[0];
rz(2.5212133) q[1];
sx q[1];
rz(-1.9199771) q[1];
sx q[1];
rz(-1.5588123) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4844783) q[0];
sx q[0];
rz(-0.14095356) q[0];
sx q[0];
rz(0.64224215) q[0];
x q[1];
rz(-1.5887268) q[2];
sx q[2];
rz(-1.6758783) q[2];
sx q[2];
rz(-1.1370657) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1190289) q[1];
sx q[1];
rz(-0.058534082) q[1];
sx q[1];
rz(-1.4083393) q[1];
rz(-pi) q[2];
rz(-1.2101733) q[3];
sx q[3];
rz(-2.7322331) q[3];
sx q[3];
rz(-1.6334074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9862426) q[2];
sx q[2];
rz(-0.68779951) q[2];
sx q[2];
rz(0.17779329) q[2];
rz(-2.6627461) q[3];
sx q[3];
rz(-1.1136473) q[3];
sx q[3];
rz(-3.0732885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9947522) q[0];
sx q[0];
rz(-2.1629592) q[0];
sx q[0];
rz(-0.11257182) q[0];
rz(0.62537891) q[1];
sx q[1];
rz(-2.0140779) q[1];
sx q[1];
rz(-0.88923997) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0028232) q[0];
sx q[0];
rz(-0.089931503) q[0];
sx q[0];
rz(2.0683209) q[0];
rz(-pi) q[1];
rz(1.6499405) q[2];
sx q[2];
rz(-1.835485) q[2];
sx q[2];
rz(1.9078329) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9611284) q[1];
sx q[1];
rz(-0.58149177) q[1];
sx q[1];
rz(1.1322458) q[1];
rz(1.806385) q[3];
sx q[3];
rz(-1.5054678) q[3];
sx q[3];
rz(0.8222841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2832977) q[2];
sx q[2];
rz(-0.12792835) q[2];
sx q[2];
rz(-2.8683786) q[2];
rz(-2.9122399) q[3];
sx q[3];
rz(-1.554824) q[3];
sx q[3];
rz(1.2685512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(-0.44723085) q[0];
sx q[0];
rz(-0.87089592) q[0];
sx q[0];
rz(0.91621512) q[0];
rz(-2.9592196) q[1];
sx q[1];
rz(-2.9353751) q[1];
sx q[1];
rz(-1.1590385) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.050738036) q[0];
sx q[0];
rz(-1.3373475) q[0];
sx q[0];
rz(-1.3755193) q[0];
rz(-pi) q[1];
x q[1];
rz(0.095511212) q[2];
sx q[2];
rz(-1.5198802) q[2];
sx q[2];
rz(1.3595622) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.43913821) q[1];
sx q[1];
rz(-1.411644) q[1];
sx q[1];
rz(1.8887146) q[1];
rz(-pi) q[2];
x q[2];
rz(0.93709903) q[3];
sx q[3];
rz(-0.94985147) q[3];
sx q[3];
rz(1.6976732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3850022) q[2];
sx q[2];
rz(-2.341541) q[2];
sx q[2];
rz(2.8263367) q[2];
rz(-2.5473525) q[3];
sx q[3];
rz(-0.88163328) q[3];
sx q[3];
rz(-0.63410223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2600128) q[0];
sx q[0];
rz(-2.4729112) q[0];
sx q[0];
rz(-0.15923937) q[0];
rz(1.0881933) q[1];
sx q[1];
rz(-0.91251487) q[1];
sx q[1];
rz(2.870627) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9054026) q[0];
sx q[0];
rz(-1.8216678) q[0];
sx q[0];
rz(1.0948927) q[0];
rz(-0.12853821) q[2];
sx q[2];
rz(-2.4646799) q[2];
sx q[2];
rz(2.235379) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3828363) q[1];
sx q[1];
rz(-2.3985632) q[1];
sx q[1];
rz(-0.50937517) q[1];
rz(0.58862092) q[3];
sx q[3];
rz(-0.56097066) q[3];
sx q[3];
rz(0.64631337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3403885) q[2];
sx q[2];
rz(-0.77216721) q[2];
sx q[2];
rz(0.32269746) q[2];
rz(0.05803756) q[3];
sx q[3];
rz(-2.3400584) q[3];
sx q[3];
rz(0.29840741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67644453) q[0];
sx q[0];
rz(-1.5100751) q[0];
sx q[0];
rz(-0.81612192) q[0];
rz(2.8836518) q[1];
sx q[1];
rz(-1.1358658) q[1];
sx q[1];
rz(1.6075016) q[1];
rz(-0.28889523) q[2];
sx q[2];
rz(-1.8452273) q[2];
sx q[2];
rz(0.065963521) q[2];
rz(-2.0097575) q[3];
sx q[3];
rz(-0.96228941) q[3];
sx q[3];
rz(-1.8214772) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
