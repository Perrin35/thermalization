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
rz(0.024700392) q[0];
sx q[0];
rz(8.165897) q[0];
sx q[0];
rz(11.293647) q[0];
rz(-4.89115) q[1];
sx q[1];
rz(1.7311544) q[1];
sx q[1];
rz(16.601736) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0860505) q[0];
sx q[0];
rz(-0.95261935) q[0];
sx q[0];
rz(2.2937141) q[0];
x q[1];
rz(1.3145163) q[2];
sx q[2];
rz(-1.8602784) q[2];
sx q[2];
rz(0.50965259) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.3013348) q[1];
sx q[1];
rz(-2.4221759) q[1];
sx q[1];
rz(0.69566984) q[1];
rz(-pi) q[2];
rz(-0.36811604) q[3];
sx q[3];
rz(-1.5541346) q[3];
sx q[3];
rz(-2.1429495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7749403) q[2];
sx q[2];
rz(-1.6795936) q[2];
sx q[2];
rz(-2.6095663) q[2];
rz(-2.4074647) q[3];
sx q[3];
rz(-2.8668154) q[3];
sx q[3];
rz(2.0348569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4183913) q[0];
sx q[0];
rz(-0.98576236) q[0];
sx q[0];
rz(3.0611839) q[0];
rz(0.53781992) q[1];
sx q[1];
rz(-2.0729013) q[1];
sx q[1];
rz(-2.417876) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2576669) q[0];
sx q[0];
rz(-0.91086713) q[0];
sx q[0];
rz(-1.0038478) q[0];
rz(-pi) q[1];
rz(-2.2369464) q[2];
sx q[2];
rz(-1.1084741) q[2];
sx q[2];
rz(-1.2733627) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.98772939) q[1];
sx q[1];
rz(-1.5692668) q[1];
sx q[1];
rz(0.9823003) q[1];
rz(-pi) q[2];
rz(-0.81775093) q[3];
sx q[3];
rz(-1.1766542) q[3];
sx q[3];
rz(-0.5241636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.85372743) q[2];
sx q[2];
rz(-2.008805) q[2];
sx q[2];
rz(0.90616027) q[2];
rz(-0.3624889) q[3];
sx q[3];
rz(-1.5972842) q[3];
sx q[3];
rz(-0.046886142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.872181) q[0];
sx q[0];
rz(-1.325664) q[0];
sx q[0];
rz(1.0915225) q[0];
rz(-0.12256924) q[1];
sx q[1];
rz(-1.4361607) q[1];
sx q[1];
rz(1.3465808) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65660921) q[0];
sx q[0];
rz(-2.7843028) q[0];
sx q[0];
rz(1.6553418) q[0];
x q[1];
rz(1.8982688) q[2];
sx q[2];
rz(-2.2041614) q[2];
sx q[2];
rz(-2.4302389) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1344188) q[1];
sx q[1];
rz(-1.6225623) q[1];
sx q[1];
rz(1.7758796) q[1];
rz(-pi) q[2];
rz(-1.4580549) q[3];
sx q[3];
rz(-2.097297) q[3];
sx q[3];
rz(2.2248588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4284105) q[2];
sx q[2];
rz(-1.1298263) q[2];
sx q[2];
rz(0.23207363) q[2];
rz(1.170916) q[3];
sx q[3];
rz(-2.5210095) q[3];
sx q[3];
rz(0.28461972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21122268) q[0];
sx q[0];
rz(-2.3486597) q[0];
sx q[0];
rz(-1.5027745) q[0];
rz(-0.59208313) q[1];
sx q[1];
rz(-1.1555669) q[1];
sx q[1];
rz(-2.2519462) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29763313) q[0];
sx q[0];
rz(-1.0464764) q[0];
sx q[0];
rz(0.53589937) q[0];
rz(0.2438898) q[2];
sx q[2];
rz(-2.0497344) q[2];
sx q[2];
rz(-2.6289491) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.46737503) q[1];
sx q[1];
rz(-1.2548037) q[1];
sx q[1];
rz(-2.9520079) q[1];
rz(1.9806421) q[3];
sx q[3];
rz(-1.92627) q[3];
sx q[3];
rz(0.4996757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.64638102) q[2];
sx q[2];
rz(-1.4633598) q[2];
sx q[2];
rz(0.67592534) q[2];
rz(1.4365139) q[3];
sx q[3];
rz(-2.8016475) q[3];
sx q[3];
rz(2.9743312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9306358) q[0];
sx q[0];
rz(-0.3322596) q[0];
sx q[0];
rz(-0.19530547) q[0];
rz(0.12061128) q[1];
sx q[1];
rz(-0.21574012) q[1];
sx q[1];
rz(-2.9511071) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2094177) q[0];
sx q[0];
rz(-3.0324192) q[0];
sx q[0];
rz(-2.5072844) q[0];
rz(0.67660013) q[2];
sx q[2];
rz(-1.6381421) q[2];
sx q[2];
rz(1.2852033) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1200285) q[1];
sx q[1];
rz(-3.0282674) q[1];
sx q[1];
rz(0.47685195) q[1];
rz(-pi) q[2];
rz(1.1498012) q[3];
sx q[3];
rz(-0.99830571) q[3];
sx q[3];
rz(1.5508625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6327989) q[2];
sx q[2];
rz(-1.6570647) q[2];
sx q[2];
rz(1.18139) q[2];
rz(-1.7975636) q[3];
sx q[3];
rz(-0.7414147) q[3];
sx q[3];
rz(-0.77176362) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8089499) q[0];
sx q[0];
rz(-1.7985666) q[0];
sx q[0];
rz(-2.0020265) q[0];
rz(2.471916) q[1];
sx q[1];
rz(-2.2256336) q[1];
sx q[1];
rz(1.12961) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.35992) q[0];
sx q[0];
rz(-1.8501213) q[0];
sx q[0];
rz(-1.2646227) q[0];
rz(1.3048346) q[2];
sx q[2];
rz(-0.87087357) q[2];
sx q[2];
rz(2.093932) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.46907963) q[1];
sx q[1];
rz(-1.7653078) q[1];
sx q[1];
rz(-0.37176337) q[1];
rz(-pi) q[2];
rz(1.4127168) q[3];
sx q[3];
rz(-1.3922833) q[3];
sx q[3];
rz(1.8399389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.42664042) q[2];
sx q[2];
rz(-1.7317829) q[2];
sx q[2];
rz(-2.7362774) q[2];
rz(-2.3146368) q[3];
sx q[3];
rz(-0.6260286) q[3];
sx q[3];
rz(2.2820182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7760794) q[0];
sx q[0];
rz(-3.1190393) q[0];
sx q[0];
rz(-0.93589163) q[0];
rz(-2.2731958) q[1];
sx q[1];
rz(-0.52803841) q[1];
sx q[1];
rz(1.7220928) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0763797) q[0];
sx q[0];
rz(-1.2166436) q[0];
sx q[0];
rz(2.5962418) q[0];
x q[1];
rz(2.9090857) q[2];
sx q[2];
rz(-1.075313) q[2];
sx q[2];
rz(-1.2459618) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.30314919) q[1];
sx q[1];
rz(-1.1080964) q[1];
sx q[1];
rz(0.82973231) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9299632) q[3];
sx q[3];
rz(-0.0018001477) q[3];
sx q[3];
rz(1.7626007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.3203656) q[2];
sx q[2];
rz(-1.7419523) q[2];
sx q[2];
rz(-2.121675) q[2];
rz(0.50926456) q[3];
sx q[3];
rz(-1.441322) q[3];
sx q[3];
rz(0.078484623) q[3];
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
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34198636) q[0];
sx q[0];
rz(-1.6260363) q[0];
sx q[0];
rz(1.9827783) q[0];
rz(0.47053567) q[1];
sx q[1];
rz(-1.7262986) q[1];
sx q[1];
rz(-1.4656969) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9500831) q[0];
sx q[0];
rz(-1.6217163) q[0];
sx q[0];
rz(-1.5945192) q[0];
x q[1];
rz(0.41185015) q[2];
sx q[2];
rz(-2.4076216) q[2];
sx q[2];
rz(-3.0628772) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.640721) q[1];
sx q[1];
rz(-1.6476008) q[1];
sx q[1];
rz(0.66052628) q[1];
rz(-pi) q[2];
rz(-1.8578148) q[3];
sx q[3];
rz(-0.94285184) q[3];
sx q[3];
rz(0.43691853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9965685) q[2];
sx q[2];
rz(-1.327876) q[2];
sx q[2];
rz(-0.81929755) q[2];
rz(-0.79646349) q[3];
sx q[3];
rz(-2.7345246) q[3];
sx q[3];
rz(1.6527294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2040937) q[0];
sx q[0];
rz(-3.0452073) q[0];
sx q[0];
rz(2.3199484) q[0];
rz(0.67283982) q[1];
sx q[1];
rz(-0.57603374) q[1];
sx q[1];
rz(-2.0988665) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0053528) q[0];
sx q[0];
rz(-1.3658821) q[0];
sx q[0];
rz(1.4968275) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.94093948) q[2];
sx q[2];
rz(-0.90412882) q[2];
sx q[2];
rz(0.43064865) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1112177) q[1];
sx q[1];
rz(-2.1407323) q[1];
sx q[1];
rz(2.6968677) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6058873) q[3];
sx q[3];
rz(-0.2920449) q[3];
sx q[3];
rz(0.10242187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3599856) q[2];
sx q[2];
rz(-1.8693417) q[2];
sx q[2];
rz(-2.2992415) q[2];
rz(0.53884566) q[3];
sx q[3];
rz(-1.4875965) q[3];
sx q[3];
rz(0.52411383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6738324) q[0];
sx q[0];
rz(-2.9957275) q[0];
sx q[0];
rz(-1.9061331) q[0];
rz(2.1927059) q[1];
sx q[1];
rz(-0.83273879) q[1];
sx q[1];
rz(-1.6611151) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47280374) q[0];
sx q[0];
rz(-1.848671) q[0];
sx q[0];
rz(0.86937277) q[0];
rz(-pi) q[1];
rz(-2.5790501) q[2];
sx q[2];
rz(-2.1058561) q[2];
sx q[2];
rz(1.3356653) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.93714156) q[1];
sx q[1];
rz(-1.6455263) q[1];
sx q[1];
rz(-0.16507574) q[1];
rz(-2.5148095) q[3];
sx q[3];
rz(-2.7434182) q[3];
sx q[3];
rz(-2.6836723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.71198717) q[2];
sx q[2];
rz(-1.7179855) q[2];
sx q[2];
rz(-2.831366) q[2];
rz(-0.45983908) q[3];
sx q[3];
rz(-0.55791563) q[3];
sx q[3];
rz(-1.3742113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9643758) q[0];
sx q[0];
rz(-1.1342659) q[0];
sx q[0];
rz(-2.076617) q[0];
rz(-2.2433157) q[1];
sx q[1];
rz(-2.5986462) q[1];
sx q[1];
rz(2.6959261) q[1];
rz(1.0679792) q[2];
sx q[2];
rz(-0.72952727) q[2];
sx q[2];
rz(-2.9222957) q[2];
rz(-2.0528005) q[3];
sx q[3];
rz(-2.4203151) q[3];
sx q[3];
rz(-2.388777) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
