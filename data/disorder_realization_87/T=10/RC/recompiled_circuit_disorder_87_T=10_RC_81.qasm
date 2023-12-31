OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9053469) q[0];
sx q[0];
rz(-0.72609225) q[0];
sx q[0];
rz(-0.2015764) q[0];
rz(0.4959313) q[1];
sx q[1];
rz(-0.5402686) q[1];
sx q[1];
rz(-0.93710605) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0721604) q[0];
sx q[0];
rz(-1.9919112) q[0];
sx q[0];
rz(2.518597) q[0];
x q[1];
rz(0.65242282) q[2];
sx q[2];
rz(-2.4180275) q[2];
sx q[2];
rz(-3.0384118) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1569251) q[1];
sx q[1];
rz(-0.19913864) q[1];
sx q[1];
rz(1.1652633) q[1];
x q[2];
rz(1.3017544) q[3];
sx q[3];
rz(-1.4770513) q[3];
sx q[3];
rz(-2.4274488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.15930882) q[2];
sx q[2];
rz(-0.26580492) q[2];
sx q[2];
rz(-3.0554331) q[2];
rz(-2.384095) q[3];
sx q[3];
rz(-0.75452724) q[3];
sx q[3];
rz(1.0936201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5646097) q[0];
sx q[0];
rz(-1.6596721) q[0];
sx q[0];
rz(1.3264867) q[0];
rz(-1.2558698) q[1];
sx q[1];
rz(-1.5763667) q[1];
sx q[1];
rz(-2.870141) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87646987) q[0];
sx q[0];
rz(-2.070283) q[0];
sx q[0];
rz(0.45544099) q[0];
rz(1.8797305) q[2];
sx q[2];
rz(-1.075282) q[2];
sx q[2];
rz(-0.77906424) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2084864) q[1];
sx q[1];
rz(-0.83175627) q[1];
sx q[1];
rz(1.0528802) q[1];
rz(-pi) q[2];
rz(1.9430964) q[3];
sx q[3];
rz(-2.3447403) q[3];
sx q[3];
rz(-0.70355319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0945956) q[2];
sx q[2];
rz(-1.7589898) q[2];
sx q[2];
rz(2.5644152) q[2];
rz(-0.92352891) q[3];
sx q[3];
rz(-2.2369592) q[3];
sx q[3];
rz(1.4097759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24401027) q[0];
sx q[0];
rz(-2.0799461) q[0];
sx q[0];
rz(-0.14257167) q[0];
rz(1.3525195) q[1];
sx q[1];
rz(-1.0357772) q[1];
sx q[1];
rz(-0.20908633) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3269743) q[0];
sx q[0];
rz(-2.1332392) q[0];
sx q[0];
rz(-2.4787865) q[0];
rz(2.694391) q[2];
sx q[2];
rz(-0.7913835) q[2];
sx q[2];
rz(0.42391047) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0774539) q[1];
sx q[1];
rz(-1.9874548) q[1];
sx q[1];
rz(1.0799079) q[1];
rz(-pi) q[2];
rz(-0.61061065) q[3];
sx q[3];
rz(-0.36452499) q[3];
sx q[3];
rz(-3.1162457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3645939) q[2];
sx q[2];
rz(-0.89407388) q[2];
sx q[2];
rz(-2.7374632) q[2];
rz(-1.8557619) q[3];
sx q[3];
rz(-1.1288246) q[3];
sx q[3];
rz(0.16734853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7053112) q[0];
sx q[0];
rz(-3.0349338) q[0];
sx q[0];
rz(2.1787815) q[0];
rz(0.46936938) q[1];
sx q[1];
rz(-0.58987394) q[1];
sx q[1];
rz(-0.00096360047) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92361802) q[0];
sx q[0];
rz(-2.1305363) q[0];
sx q[0];
rz(-2.2773507) q[0];
rz(1.5713111) q[2];
sx q[2];
rz(-1.6948023) q[2];
sx q[2];
rz(2.4027783) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5597792) q[1];
sx q[1];
rz(-2.1372791) q[1];
sx q[1];
rz(-0.15484667) q[1];
rz(-pi) q[2];
rz(-1.6603052) q[3];
sx q[3];
rz(-1.0032017) q[3];
sx q[3];
rz(1.2168509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0659539) q[2];
sx q[2];
rz(-1.6299738) q[2];
sx q[2];
rz(0.11165079) q[2];
rz(2.3305437) q[3];
sx q[3];
rz(-2.6871197) q[3];
sx q[3];
rz(-3.1276935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(1.951293) q[0];
sx q[0];
rz(-2.2409029) q[0];
sx q[0];
rz(0.35650373) q[0];
rz(0.50645343) q[1];
sx q[1];
rz(-1.7849779) q[1];
sx q[1];
rz(2.8809663) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8286752) q[0];
sx q[0];
rz(-0.14794359) q[0];
sx q[0];
rz(-0.28767985) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8904301) q[2];
sx q[2];
rz(-1.3248982) q[2];
sx q[2];
rz(2.9159301) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.063559859) q[1];
sx q[1];
rz(-1.5884001) q[1];
sx q[1];
rz(2.5672008) q[1];
rz(-pi) q[2];
rz(-2.3652606) q[3];
sx q[3];
rz(-2.9388802) q[3];
sx q[3];
rz(2.4622038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9115209) q[2];
sx q[2];
rz(-1.6610961) q[2];
sx q[2];
rz(-2.9122706) q[2];
rz(-0.54245943) q[3];
sx q[3];
rz(-0.31033236) q[3];
sx q[3];
rz(1.0361766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77263537) q[0];
sx q[0];
rz(-0.90894037) q[0];
sx q[0];
rz(1.4916346) q[0];
rz(-1.0391327) q[1];
sx q[1];
rz(-1.8455448) q[1];
sx q[1];
rz(-1.414149) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4947858) q[0];
sx q[0];
rz(-1.2959058) q[0];
sx q[0];
rz(1.6199714) q[0];
rz(-pi) q[1];
rz(-0.51775198) q[2];
sx q[2];
rz(-1.3377829) q[2];
sx q[2];
rz(-2.0875967) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0318839) q[1];
sx q[1];
rz(-2.0307396) q[1];
sx q[1];
rz(2.1214532) q[1];
rz(-0.075501637) q[3];
sx q[3];
rz(-1.7864979) q[3];
sx q[3];
rz(-2.2942033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.002939) q[2];
sx q[2];
rz(-2.5230375) q[2];
sx q[2];
rz(-3.1138528) q[2];
rz(2.6489143) q[3];
sx q[3];
rz(-1.490482) q[3];
sx q[3];
rz(-1.3940575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7572927) q[0];
sx q[0];
rz(-1.5526271) q[0];
sx q[0];
rz(0.12167715) q[0];
rz(1.9901468) q[1];
sx q[1];
rz(-2.6897488) q[1];
sx q[1];
rz(-0.40245232) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2936195) q[0];
sx q[0];
rz(-1.8403887) q[0];
sx q[0];
rz(-2.2718391) q[0];
rz(-pi) q[1];
rz(1.5624814) q[2];
sx q[2];
rz(-1.6168211) q[2];
sx q[2];
rz(-0.41551057) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.35343364) q[1];
sx q[1];
rz(-2.806059) q[1];
sx q[1];
rz(-2.5787756) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6472858) q[3];
sx q[3];
rz(-2.7923931) q[3];
sx q[3];
rz(-0.66222144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.014331269) q[2];
sx q[2];
rz(-1.971259) q[2];
sx q[2];
rz(-0.86223117) q[2];
rz(-0.47752738) q[3];
sx q[3];
rz(-1.1815485) q[3];
sx q[3];
rz(0.91807085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35571337) q[0];
sx q[0];
rz(-0.15545758) q[0];
sx q[0];
rz(-2.4086337) q[0];
rz(0.14006242) q[1];
sx q[1];
rz(-0.99761325) q[1];
sx q[1];
rz(1.0345116) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93787545) q[0];
sx q[0];
rz(-2.8563742) q[0];
sx q[0];
rz(-2.2144149) q[0];
rz(-pi) q[1];
rz(-0.10266281) q[2];
sx q[2];
rz(-1.6409988) q[2];
sx q[2];
rz(1.7662802) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4701177) q[1];
sx q[1];
rz(-2.5744994) q[1];
sx q[1];
rz(-1.6398029) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4626059) q[3];
sx q[3];
rz(-1.9868317) q[3];
sx q[3];
rz(0.070852208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.90116477) q[2];
sx q[2];
rz(-1.1952091) q[2];
sx q[2];
rz(-2.365716) q[2];
rz(-0.72426978) q[3];
sx q[3];
rz(-0.80562076) q[3];
sx q[3];
rz(1.4340713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4416606) q[0];
sx q[0];
rz(-0.76403809) q[0];
sx q[0];
rz(-0.016816703) q[0];
rz(-0.018521221) q[1];
sx q[1];
rz(-1.8073945) q[1];
sx q[1];
rz(-2.3628078) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9176365) q[0];
sx q[0];
rz(-2.7178239) q[0];
sx q[0];
rz(-0.64026041) q[0];
x q[1];
rz(2.8674556) q[2];
sx q[2];
rz(-2.0275653) q[2];
sx q[2];
rz(1.0711087) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.76965895) q[1];
sx q[1];
rz(-1.2049335) q[1];
sx q[1];
rz(-0.2497754) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3498464) q[3];
sx q[3];
rz(-1.6337992) q[3];
sx q[3];
rz(-0.27730478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4497711) q[2];
sx q[2];
rz(-1.8372953) q[2];
sx q[2];
rz(-2.4251535) q[2];
rz(1.4572432) q[3];
sx q[3];
rz(-1.4102035) q[3];
sx q[3];
rz(-1.289207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.107782) q[0];
sx q[0];
rz(-2.6066715) q[0];
sx q[0];
rz(0.15144908) q[0];
rz(-2.9653213) q[1];
sx q[1];
rz(-1.9468032) q[1];
sx q[1];
rz(-0.72296468) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82397126) q[0];
sx q[0];
rz(-0.184632) q[0];
sx q[0];
rz(-0.95438192) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9480013) q[2];
sx q[2];
rz(-1.7282439) q[2];
sx q[2];
rz(-2.8709656) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1180229) q[1];
sx q[1];
rz(-1.6260864) q[1];
sx q[1];
rz(1.0667332) q[1];
rz(2.3530657) q[3];
sx q[3];
rz(-1.2796113) q[3];
sx q[3];
rz(-2.9332719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.67939776) q[2];
sx q[2];
rz(-1.2827736) q[2];
sx q[2];
rz(-2.6160713) q[2];
rz(2.8578791) q[3];
sx q[3];
rz(-1.107639) q[3];
sx q[3];
rz(2.7450558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5491966) q[0];
sx q[0];
rz(-2.017673) q[0];
sx q[0];
rz(2.8155433) q[0];
rz(-1.0992959) q[1];
sx q[1];
rz(-2.9922843) q[1];
sx q[1];
rz(-0.86984787) q[1];
rz(0.99156689) q[2];
sx q[2];
rz(-2.3835923) q[2];
sx q[2];
rz(0.53072416) q[2];
rz(-1.3463734) q[3];
sx q[3];
rz(-1.2524458) q[3];
sx q[3];
rz(2.8129775) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
