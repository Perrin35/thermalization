OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5705559) q[0];
sx q[0];
rz(-0.67222995) q[0];
sx q[0];
rz(-1.6663405) q[0];
rz(-2.1885459) q[1];
sx q[1];
rz(-0.16521984) q[1];
sx q[1];
rz(-1.2686977) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7808963) q[0];
sx q[0];
rz(-1.3421881) q[0];
sx q[0];
rz(-2.7318888) q[0];
rz(2.6495291) q[2];
sx q[2];
rz(-2.1644219) q[2];
sx q[2];
rz(0.77897859) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.42434537) q[1];
sx q[1];
rz(-1.0043036) q[1];
sx q[1];
rz(2.9866108) q[1];
x q[2];
rz(1.8016863) q[3];
sx q[3];
rz(-0.82726631) q[3];
sx q[3];
rz(-2.8184067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1573726) q[2];
sx q[2];
rz(-0.34175384) q[2];
sx q[2];
rz(0.79322195) q[2];
rz(-2.4685517) q[3];
sx q[3];
rz(-1.0854191) q[3];
sx q[3];
rz(-1.2696772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15758812) q[0];
sx q[0];
rz(-0.80261153) q[0];
sx q[0];
rz(0.60930914) q[0];
rz(-2.9318103) q[1];
sx q[1];
rz(-2.3502626) q[1];
sx q[1];
rz(-0.9598859) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6580556) q[0];
sx q[0];
rz(-1.6153045) q[0];
sx q[0];
rz(-3.0625212) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4397214) q[2];
sx q[2];
rz(-0.49197754) q[2];
sx q[2];
rz(-2.358369) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6658191) q[1];
sx q[1];
rz(-0.83516652) q[1];
sx q[1];
rz(-0.084666157) q[1];
rz(-pi) q[2];
rz(-2.0783362) q[3];
sx q[3];
rz(-1.1843269) q[3];
sx q[3];
rz(-1.8563351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3557055) q[2];
sx q[2];
rz(-0.83928078) q[2];
sx q[2];
rz(0.4915702) q[2];
rz(0.0056886557) q[3];
sx q[3];
rz(-2.0523043) q[3];
sx q[3];
rz(-0.51386851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.046535) q[0];
sx q[0];
rz(-0.53545606) q[0];
sx q[0];
rz(-2.3989578) q[0];
rz(2.5168118) q[1];
sx q[1];
rz(-0.94015986) q[1];
sx q[1];
rz(1.0661485) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4542592) q[0];
sx q[0];
rz(-1.3479184) q[0];
sx q[0];
rz(0.14503004) q[0];
rz(2.0964618) q[2];
sx q[2];
rz(-0.36681766) q[2];
sx q[2];
rz(-0.17772533) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6934476) q[1];
sx q[1];
rz(-1.9733917) q[1];
sx q[1];
rz(-0.13685302) q[1];
rz(-pi) q[2];
rz(-0.2188628) q[3];
sx q[3];
rz(-1.0893679) q[3];
sx q[3];
rz(2.4784513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2523969) q[2];
sx q[2];
rz(-0.19813457) q[2];
sx q[2];
rz(2.5890217) q[2];
rz(2.6342454) q[3];
sx q[3];
rz(-1.0891958) q[3];
sx q[3];
rz(0.56940091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4813389) q[0];
sx q[0];
rz(-2.5642671) q[0];
sx q[0];
rz(1.2465771) q[0];
rz(1.4056816) q[1];
sx q[1];
rz(-0.77189267) q[1];
sx q[1];
rz(-2.4235922) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8193977) q[0];
sx q[0];
rz(-0.85032636) q[0];
sx q[0];
rz(2.2706991) q[0];
rz(-pi) q[1];
rz(1.9393607) q[2];
sx q[2];
rz(-0.59492009) q[2];
sx q[2];
rz(1.0234317) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.93054616) q[1];
sx q[1];
rz(-0.77096838) q[1];
sx q[1];
rz(-0.20937415) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.67881949) q[3];
sx q[3];
rz(-2.9390237) q[3];
sx q[3];
rz(0.55965078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8743073) q[2];
sx q[2];
rz(-2.3882046) q[2];
sx q[2];
rz(2.4370952) q[2];
rz(0.38257515) q[3];
sx q[3];
rz(-1.7035328) q[3];
sx q[3];
rz(-0.85324919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.660897) q[0];
sx q[0];
rz(-0.040204164) q[0];
sx q[0];
rz(0.30886343) q[0];
rz(-2.1924696) q[1];
sx q[1];
rz(-0.32052761) q[1];
sx q[1];
rz(2.5905051) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9029459) q[0];
sx q[0];
rz(-1.723105) q[0];
sx q[0];
rz(0.39350827) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6597719) q[2];
sx q[2];
rz(-2.4818261) q[2];
sx q[2];
rz(-2.8967186) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6092564) q[1];
sx q[1];
rz(-2.9279531) q[1];
sx q[1];
rz(0.81476252) q[1];
rz(-pi) q[2];
rz(2.3651028) q[3];
sx q[3];
rz(-1.6740419) q[3];
sx q[3];
rz(2.6342027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4264195) q[2];
sx q[2];
rz(-2.3931563) q[2];
sx q[2];
rz(0.61881649) q[2];
rz(2.5185781) q[3];
sx q[3];
rz(-0.84191936) q[3];
sx q[3];
rz(-2.714341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.039336786) q[0];
sx q[0];
rz(-2.4781041) q[0];
sx q[0];
rz(-2.6797507) q[0];
rz(0.49679187) q[1];
sx q[1];
rz(-1.3235612) q[1];
sx q[1];
rz(1.7740446) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.753938) q[0];
sx q[0];
rz(-0.68635041) q[0];
sx q[0];
rz(1.4453056) q[0];
x q[1];
rz(-0.15256792) q[2];
sx q[2];
rz(-2.0629473) q[2];
sx q[2];
rz(-2.9967665) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9744579) q[1];
sx q[1];
rz(-0.30213812) q[1];
sx q[1];
rz(1.664851) q[1];
rz(-pi) q[2];
x q[2];
rz(0.9487834) q[3];
sx q[3];
rz(-1.7025196) q[3];
sx q[3];
rz(-1.458481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.10609145) q[2];
sx q[2];
rz(-1.1320628) q[2];
sx q[2];
rz(-2.3512225) q[2];
rz(0.60449374) q[3];
sx q[3];
rz(-1.0249745) q[3];
sx q[3];
rz(-0.42009556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6919493) q[0];
sx q[0];
rz(-2.255891) q[0];
sx q[0];
rz(-2.7272136) q[0];
rz(-0.62037933) q[1];
sx q[1];
rz(-1.2216156) q[1];
sx q[1];
rz(1.5588123) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.010084933) q[0];
sx q[0];
rz(-1.6835308) q[0];
sx q[0];
rz(-1.4860064) q[0];
rz(-0.10509872) q[2];
sx q[2];
rz(-1.5529648) q[2];
sx q[2];
rz(-2.7097429) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0225637) q[1];
sx q[1];
rz(-0.058534082) q[1];
sx q[1];
rz(1.7332533) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2101733) q[3];
sx q[3];
rz(-0.40935959) q[3];
sx q[3];
rz(-1.5081852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9862426) q[2];
sx q[2];
rz(-0.68779951) q[2];
sx q[2];
rz(0.17779329) q[2];
rz(-0.47884652) q[3];
sx q[3];
rz(-2.0279453) q[3];
sx q[3];
rz(0.068304121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14684045) q[0];
sx q[0];
rz(-2.1629592) q[0];
sx q[0];
rz(-0.11257182) q[0];
rz(0.62537891) q[1];
sx q[1];
rz(-2.0140779) q[1];
sx q[1];
rz(2.2523527) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5035985) q[0];
sx q[0];
rz(-1.6498008) q[0];
sx q[0];
rz(-3.0985831) q[0];
rz(-2.8577789) q[2];
sx q[2];
rz(-0.2760016) q[2];
sx q[2];
rz(2.2021879) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.016805033) q[1];
sx q[1];
rz(-1.3353925) q[1];
sx q[1];
rz(1.0339869) q[1];
x q[2];
rz(-1.806385) q[3];
sx q[3];
rz(-1.6361248) q[3];
sx q[3];
rz(-2.3193086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2832977) q[2];
sx q[2];
rz(-3.0136643) q[2];
sx q[2];
rz(0.27321401) q[2];
rz(0.22935271) q[3];
sx q[3];
rz(-1.5867686) q[3];
sx q[3];
rz(-1.2685512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6943618) q[0];
sx q[0];
rz(-0.87089592) q[0];
sx q[0];
rz(0.91621512) q[0];
rz(-0.18237309) q[1];
sx q[1];
rz(-0.20621754) q[1];
sx q[1];
rz(1.9825541) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6672598) q[0];
sx q[0];
rz(-1.3808819) q[0];
sx q[0];
rz(2.9037881) q[0];
rz(-pi) q[1];
rz(2.6508337) q[2];
sx q[2];
rz(-3.033394) q[2];
sx q[2];
rz(0.27709093) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4587443) q[1];
sx q[1];
rz(-2.7872751) q[1];
sx q[1];
rz(2.0451727) q[1];
rz(-pi) q[2];
x q[2];
rz(0.69144122) q[3];
sx q[3];
rz(-0.85606282) q[3];
sx q[3];
rz(-0.5428398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
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
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88157982) q[0];
sx q[0];
rz(-2.4729112) q[0];
sx q[0];
rz(-0.15923937) q[0];
rz(1.0881933) q[1];
sx q[1];
rz(-0.91251487) q[1];
sx q[1];
rz(2.870627) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20734678) q[0];
sx q[0];
rz(-2.0306315) q[0];
sx q[0];
rz(-2.8608972) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.12853821) q[2];
sx q[2];
rz(-0.67691278) q[2];
sx q[2];
rz(0.9062137) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.20226711) q[1];
sx q[1];
rz(-1.9069873) q[1];
sx q[1];
rz(-0.675981) q[1];
rz(-pi) q[2];
rz(0.48153523) q[3];
sx q[3];
rz(-1.2709444) q[3];
sx q[3];
rz(1.7026342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8012041) q[2];
sx q[2];
rz(-2.3694254) q[2];
sx q[2];
rz(0.32269746) q[2];
rz(-3.0835551) q[3];
sx q[3];
rz(-0.80153424) q[3];
sx q[3];
rz(2.8431852) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4651481) q[0];
sx q[0];
rz(-1.6315176) q[0];
sx q[0];
rz(2.3254707) q[0];
rz(0.25794087) q[1];
sx q[1];
rz(-2.0057269) q[1];
sx q[1];
rz(-1.534091) q[1];
rz(0.77946812) q[2];
sx q[2];
rz(-0.39579724) q[2];
sx q[2];
rz(-0.76553065) q[2];
rz(0.54775379) q[3];
sx q[3];
rz(-0.73368254) q[3];
sx q[3];
rz(2.0077326) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
