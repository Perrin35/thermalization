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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0451242) q[0];
sx q[0];
rz(-0.91334263) q[0];
sx q[0];
rz(2.3923001) q[0];
rz(-1.3145163) q[2];
sx q[2];
rz(-1.2813142) q[2];
sx q[2];
rz(-2.6319401) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8302248) q[1];
sx q[1];
rz(-1.1347949) q[1];
sx q[1];
rz(0.59200759) q[1];
x q[2];
rz(2.7734766) q[3];
sx q[3];
rz(-1.5541346) q[3];
sx q[3];
rz(0.9986432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7749403) q[2];
sx q[2];
rz(-1.6795936) q[2];
sx q[2];
rz(2.6095663) q[2];
rz(0.73412791) q[3];
sx q[3];
rz(-0.2747772) q[3];
sx q[3];
rz(1.1067357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.4183913) q[0];
sx q[0];
rz(-2.1558303) q[0];
sx q[0];
rz(0.080408737) q[0];
rz(0.53781992) q[1];
sx q[1];
rz(-2.0729013) q[1];
sx q[1];
rz(-2.417876) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0619556) q[0];
sx q[0];
rz(-2.300206) q[0];
sx q[0];
rz(2.5361912) q[0];
x q[1];
rz(-0.89214708) q[2];
sx q[2];
rz(-2.3513459) q[2];
sx q[2];
rz(-0.21871601) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1538633) q[1];
sx q[1];
rz(-1.5692668) q[1];
sx q[1];
rz(2.1592924) q[1];
rz(-pi) q[2];
rz(-0.51809727) q[3];
sx q[3];
rz(-2.2541917) q[3];
sx q[3];
rz(-0.70113126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2878652) q[2];
sx q[2];
rz(-2.008805) q[2];
sx q[2];
rz(2.2354324) q[2];
rz(0.3624889) q[3];
sx q[3];
rz(-1.5972842) q[3];
sx q[3];
rz(0.046886142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.872181) q[0];
sx q[0];
rz(-1.8159287) q[0];
sx q[0];
rz(-1.0915225) q[0];
rz(-3.0190234) q[1];
sx q[1];
rz(-1.7054319) q[1];
sx q[1];
rz(-1.7950119) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65660921) q[0];
sx q[0];
rz(-0.35728982) q[0];
sx q[0];
rz(-1.4862509) q[0];
x q[1];
rz(0.41287072) q[2];
sx q[2];
rz(-2.4390568) q[2];
sx q[2];
rz(1.9091878) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.57438499) q[1];
sx q[1];
rz(-1.775601) q[1];
sx q[1];
rz(0.052871953) q[1];
rz(0.52927203) q[3];
sx q[3];
rz(-1.6682169) q[3];
sx q[3];
rz(2.4306963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.71318212) q[2];
sx q[2];
rz(-1.1298263) q[2];
sx q[2];
rz(0.23207363) q[2];
rz(1.9706767) q[3];
sx q[3];
rz(-2.5210095) q[3];
sx q[3];
rz(-0.28461972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21122268) q[0];
sx q[0];
rz(-2.3486597) q[0];
sx q[0];
rz(1.5027745) q[0];
rz(-0.59208313) q[1];
sx q[1];
rz(-1.9860257) q[1];
sx q[1];
rz(2.2519462) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9735721) q[0];
sx q[0];
rz(-2.4102927) q[0];
sx q[0];
rz(-0.84748737) q[0];
rz(-pi) q[1];
rz(-1.1355033) q[2];
sx q[2];
rz(-2.6084628) q[2];
sx q[2];
rz(-1.0077623) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1210416) q[1];
sx q[1];
rz(-0.36686037) q[1];
sx q[1];
rz(1.0479142) q[1];
x q[2];
rz(0.82075228) q[3];
sx q[3];
rz(-2.6057557) q[3];
sx q[3];
rz(-1.395063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.64638102) q[2];
sx q[2];
rz(-1.4633598) q[2];
sx q[2];
rz(2.4656673) q[2];
rz(-1.4365139) q[3];
sx q[3];
rz(-2.8016475) q[3];
sx q[3];
rz(0.16726141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-1.9306358) q[0];
sx q[0];
rz(-2.8093331) q[0];
sx q[0];
rz(-2.9462872) q[0];
rz(3.0209814) q[1];
sx q[1];
rz(-0.21574012) q[1];
sx q[1];
rz(-0.19048555) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.932175) q[0];
sx q[0];
rz(-3.0324192) q[0];
sx q[0];
rz(0.63430826) q[0];
rz(-1.6570857) q[2];
sx q[2];
rz(-0.89601529) q[2];
sx q[2];
rz(0.23160411) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.45791679) q[1];
sx q[1];
rz(-1.4701587) q[1];
sx q[1];
rz(1.6229873) q[1];
rz(-pi) q[2];
rz(1.9917914) q[3];
sx q[3];
rz(-0.99830571) q[3];
sx q[3];
rz(-1.5508625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6327989) q[2];
sx q[2];
rz(-1.6570647) q[2];
sx q[2];
rz(-1.9602027) q[2];
rz(-1.3440291) q[3];
sx q[3];
rz(-2.400178) q[3];
sx q[3];
rz(2.369829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.8089499) q[0];
sx q[0];
rz(-1.3430261) q[0];
sx q[0];
rz(-1.1395662) q[0];
rz(0.66967669) q[1];
sx q[1];
rz(-0.91595903) q[1];
sx q[1];
rz(-2.0119827) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2655339) q[0];
sx q[0];
rz(-1.8647412) q[0];
sx q[0];
rz(-2.8493899) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8367581) q[2];
sx q[2];
rz(-2.2707191) q[2];
sx q[2];
rz(2.093932) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.46907963) q[1];
sx q[1];
rz(-1.3762849) q[1];
sx q[1];
rz(-0.37176337) q[1];
rz(-0.7173921) q[3];
sx q[3];
rz(-2.903707) q[3];
sx q[3];
rz(-0.57008509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.42664042) q[2];
sx q[2];
rz(-1.7317829) q[2];
sx q[2];
rz(-2.7362774) q[2];
rz(-2.3146368) q[3];
sx q[3];
rz(-2.5155641) q[3];
sx q[3];
rz(-2.2820182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3655132) q[0];
sx q[0];
rz(-3.1190393) q[0];
sx q[0];
rz(0.93589163) q[0];
rz(2.2731958) q[1];
sx q[1];
rz(-2.6135542) q[1];
sx q[1];
rz(-1.4194999) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.065212981) q[0];
sx q[0];
rz(-1.924949) q[0];
sx q[0];
rz(-0.54535086) q[0];
rz(1.1677891) q[2];
sx q[2];
rz(-2.5984077) q[2];
sx q[2];
rz(-0.78389558) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3277927) q[1];
sx q[1];
rz(-0.84974242) q[1];
sx q[1];
rz(-0.93443296) q[1];
rz(-1.5711745) q[3];
sx q[3];
rz(-1.5725563) q[3];
sx q[3];
rz(1.5906217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.821227) q[2];
sx q[2];
rz(-1.3996404) q[2];
sx q[2];
rz(-1.0199176) q[2];
rz(-0.50926456) q[3];
sx q[3];
rz(-1.7002707) q[3];
sx q[3];
rz(0.078484623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34198636) q[0];
sx q[0];
rz(-1.5155563) q[0];
sx q[0];
rz(-1.9827783) q[0];
rz(0.47053567) q[1];
sx q[1];
rz(-1.4152941) q[1];
sx q[1];
rz(1.4656969) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5138595) q[0];
sx q[0];
rz(-0.056170551) q[0];
sx q[0];
rz(-0.43561952) q[0];
x q[1];
rz(-0.41185015) q[2];
sx q[2];
rz(-0.73397103) q[2];
sx q[2];
rz(0.07871544) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.50087167) q[1];
sx q[1];
rz(-1.6476008) q[1];
sx q[1];
rz(-0.66052628) q[1];
rz(0.64792525) q[3];
sx q[3];
rz(-1.8019391) q[3];
sx q[3];
rz(-2.1794139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9965685) q[2];
sx q[2];
rz(-1.327876) q[2];
sx q[2];
rz(0.81929755) q[2];
rz(-2.3451292) q[3];
sx q[3];
rz(-0.4070681) q[3];
sx q[3];
rz(1.6527294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9374989) q[0];
sx q[0];
rz(-3.0452073) q[0];
sx q[0];
rz(-2.3199484) q[0];
rz(-2.4687528) q[1];
sx q[1];
rz(-2.5655589) q[1];
sx q[1];
rz(2.0988665) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0053528) q[0];
sx q[0];
rz(-1.7757105) q[0];
sx q[0];
rz(1.6447652) q[0];
x q[1];
rz(-0.64260245) q[2];
sx q[2];
rz(-2.2589141) q[2];
sx q[2];
rz(-2.7049899) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.71127001) q[1];
sx q[1];
rz(-1.9414328) q[1];
sx q[1];
rz(2.1881585) q[1];
rz(0.25298499) q[3];
sx q[3];
rz(-1.4232985) q[3];
sx q[3];
rz(1.1563039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3599856) q[2];
sx q[2];
rz(-1.2722509) q[2];
sx q[2];
rz(2.2992415) q[2];
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
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.46776029) q[0];
sx q[0];
rz(-2.9957275) q[0];
sx q[0];
rz(-1.2354596) q[0];
rz(0.94888672) q[1];
sx q[1];
rz(-2.3088539) q[1];
sx q[1];
rz(1.4804776) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4119856) q[0];
sx q[0];
rz(-0.74568891) q[0];
sx q[0];
rz(1.1545769) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.5625426) q[2];
sx q[2];
rz(-1.0357366) q[2];
sx q[2];
rz(-1.8059274) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.64609194) q[1];
sx q[1];
rz(-1.4061855) q[1];
sx q[1];
rz(-1.4950404) q[1];
rz(-2.5148095) q[3];
sx q[3];
rz(-2.7434182) q[3];
sx q[3];
rz(-2.6836723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4296055) q[2];
sx q[2];
rz(-1.4236071) q[2];
sx q[2];
rz(-0.31022662) q[2];
rz(-0.45983908) q[3];
sx q[3];
rz(-0.55791563) q[3];
sx q[3];
rz(1.7673813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9643758) q[0];
sx q[0];
rz(-2.0073267) q[0];
sx q[0];
rz(1.0649756) q[0];
rz(0.89827697) q[1];
sx q[1];
rz(-2.5986462) q[1];
sx q[1];
rz(2.6959261) q[1];
rz(-1.0679792) q[2];
sx q[2];
rz(-2.4120654) q[2];
sx q[2];
rz(0.21929693) q[2];
rz(1.0887922) q[3];
sx q[3];
rz(-2.4203151) q[3];
sx q[3];
rz(-2.388777) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
