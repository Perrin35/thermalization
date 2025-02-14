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
rz(-0.7402339) q[0];
sx q[0];
rz(4.8010173) q[0];
sx q[0];
rz(12.231449) q[0];
rz(0.51796335) q[1];
sx q[1];
rz(5.2809102) q[1];
sx q[1];
rz(10.032293) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54017055) q[0];
sx q[0];
rz(-0.54781944) q[0];
sx q[0];
rz(1.1146077) q[0];
rz(-pi) q[1];
rz(-1.0614228) q[2];
sx q[2];
rz(-1.4685838) q[2];
sx q[2];
rz(-1.4776023) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.40558896) q[1];
sx q[1];
rz(-2.0479255) q[1];
sx q[1];
rz(-1.3911584) q[1];
x q[2];
rz(-0.88413864) q[3];
sx q[3];
rz(-0.55301412) q[3];
sx q[3];
rz(-0.64698863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.5241549) q[2];
sx q[2];
rz(-1.2158771) q[2];
sx q[2];
rz(2.4784135) q[2];
rz(-0.080862008) q[3];
sx q[3];
rz(-0.20564779) q[3];
sx q[3];
rz(-1.9614356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8993768) q[0];
sx q[0];
rz(-1.7689393) q[0];
sx q[0];
rz(2.2826165) q[0];
rz(1.8513177) q[1];
sx q[1];
rz(-1.4651508) q[1];
sx q[1];
rz(1.4069517) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1024723) q[0];
sx q[0];
rz(-0.539383) q[0];
sx q[0];
rz(-1.1155737) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2386049) q[2];
sx q[2];
rz(-1.536035) q[2];
sx q[2];
rz(-1.3369651) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.42900742) q[1];
sx q[1];
rz(-2.6204094) q[1];
sx q[1];
rz(1.3888478) q[1];
rz(0.59515335) q[3];
sx q[3];
rz(-0.66616026) q[3];
sx q[3];
rz(2.7056138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0591639) q[2];
sx q[2];
rz(-0.51247207) q[2];
sx q[2];
rz(-1.814369) q[2];
rz(-1.0446769) q[3];
sx q[3];
rz(-2.4278214) q[3];
sx q[3];
rz(0.41675848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5752207) q[0];
sx q[0];
rz(-2.2307668) q[0];
sx q[0];
rz(0.4253934) q[0];
rz(1.7644024) q[1];
sx q[1];
rz(-1.6433989) q[1];
sx q[1];
rz(-1.707071) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1578428) q[0];
sx q[0];
rz(-0.93261496) q[0];
sx q[0];
rz(1.7179836) q[0];
rz(-pi) q[1];
rz(2.4332831) q[2];
sx q[2];
rz(-1.5454486) q[2];
sx q[2];
rz(-0.011653221) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.19615281) q[1];
sx q[1];
rz(-1.9827794) q[1];
sx q[1];
rz(-2.4293071) q[1];
rz(-1.9886964) q[3];
sx q[3];
rz(-0.66158453) q[3];
sx q[3];
rz(1.0583056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7003358) q[2];
sx q[2];
rz(-1.2563027) q[2];
sx q[2];
rz(1.8505992) q[2];
rz(1.357632) q[3];
sx q[3];
rz(-1.5057526) q[3];
sx q[3];
rz(-1.4972081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
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
rz(0.62035471) q[0];
sx q[0];
rz(-2.1339895) q[0];
sx q[0];
rz(0.77019101) q[0];
rz(-2.1417446) q[1];
sx q[1];
rz(-2.5412173) q[1];
sx q[1];
rz(-1.8005449) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5504549) q[0];
sx q[0];
rz(-1.7831752) q[0];
sx q[0];
rz(2.5881833) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2571428) q[2];
sx q[2];
rz(-1.2053066) q[2];
sx q[2];
rz(-0.7892424) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.25136687) q[1];
sx q[1];
rz(-2.736675) q[1];
sx q[1];
rz(-0.7271073) q[1];
rz(-pi) q[2];
rz(-1.8561268) q[3];
sx q[3];
rz(-1.4754681) q[3];
sx q[3];
rz(0.25432977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3406713) q[2];
sx q[2];
rz(-1.069671) q[2];
sx q[2];
rz(-2.3731903) q[2];
rz(-3.0411804) q[3];
sx q[3];
rz(-1.5407591) q[3];
sx q[3];
rz(1.2787308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1860344) q[0];
sx q[0];
rz(-2.8362507) q[0];
sx q[0];
rz(2.9146063) q[0];
rz(-1.3817894) q[1];
sx q[1];
rz(-2.558936) q[1];
sx q[1];
rz(-1.7452128) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3811548) q[0];
sx q[0];
rz(-1.210307) q[0];
sx q[0];
rz(-1.9442149) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.41047217) q[2];
sx q[2];
rz(-1.7232401) q[2];
sx q[2];
rz(-1.6688011) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.760031) q[1];
sx q[1];
rz(-0.9388105) q[1];
sx q[1];
rz(1.6595592) q[1];
rz(-1.6038069) q[3];
sx q[3];
rz(-2.3950658) q[3];
sx q[3];
rz(1.3804264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9153626) q[2];
sx q[2];
rz(-1.9910944) q[2];
sx q[2];
rz(0.076233141) q[2];
rz(-1.7539615) q[3];
sx q[3];
rz(-2.1662655) q[3];
sx q[3];
rz(1.4001728) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48080322) q[0];
sx q[0];
rz(-1.5810409) q[0];
sx q[0];
rz(-1.3355108) q[0];
rz(-2.238359) q[1];
sx q[1];
rz(-1.3390373) q[1];
sx q[1];
rz(-0.18641557) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6680697) q[0];
sx q[0];
rz(-1.1374363) q[0];
sx q[0];
rz(2.1587203) q[0];
rz(-0.30419402) q[2];
sx q[2];
rz(-0.59202164) q[2];
sx q[2];
rz(-1.5409868) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.820436) q[1];
sx q[1];
rz(-1.3296491) q[1];
sx q[1];
rz(2.2187869) q[1];
rz(-pi) q[2];
rz(-2.647688) q[3];
sx q[3];
rz(-1.1260238) q[3];
sx q[3];
rz(2.1879823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.67941252) q[2];
sx q[2];
rz(-1.1932411) q[2];
sx q[2];
rz(-1.3235486) q[2];
rz(1.9994252) q[3];
sx q[3];
rz(-0.78502941) q[3];
sx q[3];
rz(0.89046684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6776176) q[0];
sx q[0];
rz(-0.1828201) q[0];
sx q[0];
rz(-0.63419813) q[0];
rz(1.0448666) q[1];
sx q[1];
rz(-0.8909145) q[1];
sx q[1];
rz(0.96010906) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38518181) q[0];
sx q[0];
rz(-2.7419081) q[0];
sx q[0];
rz(1.5002285) q[0];
rz(-pi) q[1];
rz(-1.5788001) q[2];
sx q[2];
rz(-1.8718613) q[2];
sx q[2];
rz(2.0640399) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.81183103) q[1];
sx q[1];
rz(-1.0454181) q[1];
sx q[1];
rz(0.01339042) q[1];
rz(0.76310254) q[3];
sx q[3];
rz(-1.0013442) q[3];
sx q[3];
rz(-2.4214937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1977957) q[2];
sx q[2];
rz(-0.65239492) q[2];
sx q[2];
rz(1.5459527) q[2];
rz(1.4173077) q[3];
sx q[3];
rz(-2.3167819) q[3];
sx q[3];
rz(-1.6212757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6087795) q[0];
sx q[0];
rz(-2.9488035) q[0];
sx q[0];
rz(0.97484318) q[0];
rz(-0.10803647) q[1];
sx q[1];
rz(-1.255475) q[1];
sx q[1];
rz(1.9727762) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2377396) q[0];
sx q[0];
rz(-1.0971591) q[0];
sx q[0];
rz(2.6698378) q[0];
rz(-pi) q[1];
rz(-1.4270876e-05) q[2];
sx q[2];
rz(-0.76185267) q[2];
sx q[2];
rz(0.95214168) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7119693) q[1];
sx q[1];
rz(-1.483146) q[1];
sx q[1];
rz(-1.3369249) q[1];
rz(-pi) q[2];
x q[2];
rz(0.036691477) q[3];
sx q[3];
rz(-0.49202575) q[3];
sx q[3];
rz(-1.8560611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.10890266) q[2];
sx q[2];
rz(-0.87778512) q[2];
sx q[2];
rz(2.4857793) q[2];
rz(3.0374895) q[3];
sx q[3];
rz(-1.3452353) q[3];
sx q[3];
rz(2.9534705) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8769237) q[0];
sx q[0];
rz(-1.8380565) q[0];
sx q[0];
rz(-1.6498097) q[0];
rz(-1.0393633) q[1];
sx q[1];
rz(-2.3014258) q[1];
sx q[1];
rz(2.6731491) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.428498) q[0];
sx q[0];
rz(-1.5665496) q[0];
sx q[0];
rz(0.18761401) q[0];
rz(-pi) q[1];
rz(-2.5889187) q[2];
sx q[2];
rz(-0.50853339) q[2];
sx q[2];
rz(-2.7486211) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0026928) q[1];
sx q[1];
rz(-1.7322092) q[1];
sx q[1];
rz(-2.3494865) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3571068) q[3];
sx q[3];
rz(-0.38803369) q[3];
sx q[3];
rz(1.3875543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.046772) q[2];
sx q[2];
rz(-1.5609317) q[2];
sx q[2];
rz(2.2136733) q[2];
rz(-1.8038484) q[3];
sx q[3];
rz(-2.6313621) q[3];
sx q[3];
rz(1.6524338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.054166404) q[0];
sx q[0];
rz(-1.0324284) q[0];
sx q[0];
rz(1.077865) q[0];
rz(2.7742591) q[1];
sx q[1];
rz(-1.3307738) q[1];
sx q[1];
rz(1.8574538) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36791641) q[0];
sx q[0];
rz(-1.8225637) q[0];
sx q[0];
rz(2.0579703) q[0];
x q[1];
rz(-1.7214891) q[2];
sx q[2];
rz(-0.65905276) q[2];
sx q[2];
rz(2.6113752) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.974677) q[1];
sx q[1];
rz(-1.484718) q[1];
sx q[1];
rz(-2.9831216) q[1];
rz(-pi) q[2];
rz(-2.908091) q[3];
sx q[3];
rz(-0.74849183) q[3];
sx q[3];
rz(-0.12602636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1401356) q[2];
sx q[2];
rz(-0.45150253) q[2];
sx q[2];
rz(-0.4293116) q[2];
rz(0.15520994) q[3];
sx q[3];
rz(-2.8838172) q[3];
sx q[3];
rz(1.8163053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8975288) q[0];
sx q[0];
rz(-1.5915992) q[0];
sx q[0];
rz(-1.5912548) q[0];
rz(-0.86391972) q[1];
sx q[1];
rz(-0.37352957) q[1];
sx q[1];
rz(-1.4600798) q[1];
rz(0.39045329) q[2];
sx q[2];
rz(-2.5627372) q[2];
sx q[2];
rz(2.550203) q[2];
rz(0.94384296) q[3];
sx q[3];
rz(-1.8494434) q[3];
sx q[3];
rz(-0.2260126) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
