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
rz(2.0848701) q[0];
sx q[0];
rz(-1.8615847) q[0];
sx q[0];
rz(-2.4155389) q[0];
rz(0.59745204) q[1];
sx q[1];
rz(3.6570972) q[1];
sx q[1];
rz(7.3154156) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9361794) q[0];
sx q[0];
rz(-1.6717311) q[0];
sx q[0];
rz(-2.4839782) q[0];
rz(-2.0886253) q[2];
sx q[2];
rz(-1.6047239) q[2];
sx q[2];
rz(-2.9133493) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.91019708) q[1];
sx q[1];
rz(-0.33581844) q[1];
sx q[1];
rz(0.37792716) q[1];
rz(-2.0809523) q[3];
sx q[3];
rz(-0.83612305) q[3];
sx q[3];
rz(1.5165129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7008179) q[2];
sx q[2];
rz(-0.60443193) q[2];
sx q[2];
rz(0.087892858) q[2];
rz(-1.6739738) q[3];
sx q[3];
rz(-1.847495) q[3];
sx q[3];
rz(0.99883336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0198332) q[0];
sx q[0];
rz(-1.3061981) q[0];
sx q[0];
rz(1.649296) q[0];
rz(-2.3536033) q[1];
sx q[1];
rz(-1.9580656) q[1];
sx q[1];
rz(-2.1764887) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0168646) q[0];
sx q[0];
rz(-1.6627321) q[0];
sx q[0];
rz(-0.73081907) q[0];
rz(1.4282706) q[2];
sx q[2];
rz(-2.3025142) q[2];
sx q[2];
rz(-2.1676262) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1652884) q[1];
sx q[1];
rz(-1.7563662) q[1];
sx q[1];
rz(-1.1052422) q[1];
x q[2];
rz(2.6718465) q[3];
sx q[3];
rz(-1.489991) q[3];
sx q[3];
rz(-1.4572136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1284156) q[2];
sx q[2];
rz(-0.33856496) q[2];
sx q[2];
rz(-1.817912) q[2];
rz(-0.88661083) q[3];
sx q[3];
rz(-1.1611791) q[3];
sx q[3];
rz(2.9429341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.897268) q[0];
sx q[0];
rz(-2.1710945) q[0];
sx q[0];
rz(-2.4460728) q[0];
rz(2.3059402) q[1];
sx q[1];
rz(-1.7213768) q[1];
sx q[1];
rz(-0.64170352) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8678829) q[0];
sx q[0];
rz(-0.84016227) q[0];
sx q[0];
rz(-2.9466936) q[0];
rz(-pi) q[1];
rz(0.19410816) q[2];
sx q[2];
rz(-2.420937) q[2];
sx q[2];
rz(-2.1168762) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1759626) q[1];
sx q[1];
rz(-1.2568736) q[1];
sx q[1];
rz(2.0759275) q[1];
rz(-0.68797994) q[3];
sx q[3];
rz(-1.1065346) q[3];
sx q[3];
rz(-1.9339428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9782763) q[2];
sx q[2];
rz(-1.1288319) q[2];
sx q[2];
rz(1.3345435) q[2];
rz(1.4000019) q[3];
sx q[3];
rz(-1.497523) q[3];
sx q[3];
rz(2.002772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.9982933) q[0];
sx q[0];
rz(-2.2704953) q[0];
sx q[0];
rz(-0.77348462) q[0];
rz(0.086291226) q[1];
sx q[1];
rz(-2.6373865) q[1];
sx q[1];
rz(-2.6533244) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.047866743) q[0];
sx q[0];
rz(-0.4673839) q[0];
sx q[0];
rz(2.9899238) q[0];
rz(-pi) q[1];
rz(-2.1147229) q[2];
sx q[2];
rz(-2.1928117) q[2];
sx q[2];
rz(2.029947) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3662373) q[1];
sx q[1];
rz(-1.3455433) q[1];
sx q[1];
rz(-2.3083592) q[1];
x q[2];
rz(0.78924307) q[3];
sx q[3];
rz(-2.5960138) q[3];
sx q[3];
rz(0.18462791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.46828541) q[2];
sx q[2];
rz(-1.1136592) q[2];
sx q[2];
rz(-2.9735273) q[2];
rz(0.42260653) q[3];
sx q[3];
rz(-2.1674619) q[3];
sx q[3];
rz(-2.9882123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76822686) q[0];
sx q[0];
rz(-2.234937) q[0];
sx q[0];
rz(-1.9510829) q[0];
rz(2.6137784) q[1];
sx q[1];
rz(-2.0595136) q[1];
sx q[1];
rz(-3.1324918) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4531698) q[0];
sx q[0];
rz(-2.1707782) q[0];
sx q[0];
rz(-2.1196516) q[0];
rz(-pi) q[1];
rz(1.509769) q[2];
sx q[2];
rz(-1.6313071) q[2];
sx q[2];
rz(0.64693816) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.121351) q[1];
sx q[1];
rz(-0.96881908) q[1];
sx q[1];
rz(0.51687981) q[1];
rz(-pi) q[2];
rz(-0.9420932) q[3];
sx q[3];
rz(-1.9607984) q[3];
sx q[3];
rz(-0.77446924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3147543) q[2];
sx q[2];
rz(-2.645731) q[2];
sx q[2];
rz(0.66519386) q[2];
rz(1.0264171) q[3];
sx q[3];
rz(-1.0746936) q[3];
sx q[3];
rz(0.78072602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6219532) q[0];
sx q[0];
rz(-2.6060947) q[0];
sx q[0];
rz(2.7435379) q[0];
rz(1.4705426) q[1];
sx q[1];
rz(-2.2649951) q[1];
sx q[1];
rz(0.7801396) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.018293457) q[0];
sx q[0];
rz(-1.9255877) q[0];
sx q[0];
rz(1.0650207) q[0];
x q[1];
rz(-2.8988437) q[2];
sx q[2];
rz(-1.374259) q[2];
sx q[2];
rz(2.8145973) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4698339) q[1];
sx q[1];
rz(-0.90171725) q[1];
sx q[1];
rz(-1.9356739) q[1];
x q[2];
rz(-2.8242802) q[3];
sx q[3];
rz(-0.91083357) q[3];
sx q[3];
rz(2.8803006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.857343) q[2];
sx q[2];
rz(-1.1025068) q[2];
sx q[2];
rz(-0.16172376) q[2];
rz(-3.0297847) q[3];
sx q[3];
rz(-0.72607741) q[3];
sx q[3];
rz(2.4017754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.9991456) q[0];
sx q[0];
rz(-1.4302) q[0];
sx q[0];
rz(2.4430742) q[0];
rz(-1.0922208) q[1];
sx q[1];
rz(-0.75463808) q[1];
sx q[1];
rz(-3.0879367) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2630477) q[0];
sx q[0];
rz(-1.8302655) q[0];
sx q[0];
rz(0.83570133) q[0];
rz(-0.38044615) q[2];
sx q[2];
rz(-1.9815677) q[2];
sx q[2];
rz(0.10027567) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0567017) q[1];
sx q[1];
rz(-2.6680601) q[1];
sx q[1];
rz(-0.11736203) q[1];
rz(-pi) q[2];
rz(-0.92412432) q[3];
sx q[3];
rz(-1.4507626) q[3];
sx q[3];
rz(1.2965352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0003537) q[2];
sx q[2];
rz(-2.443479) q[2];
sx q[2];
rz(0.21256438) q[2];
rz(-2.8138748) q[3];
sx q[3];
rz(-0.98723427) q[3];
sx q[3];
rz(-1.3535961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.245529) q[0];
sx q[0];
rz(-2.0241757) q[0];
sx q[0];
rz(-2.8118706) q[0];
rz(2.3700736) q[1];
sx q[1];
rz(-2.3604269) q[1];
sx q[1];
rz(0.82569295) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86733183) q[0];
sx q[0];
rz(-1.6578339) q[0];
sx q[0];
rz(2.4768193) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4554899) q[2];
sx q[2];
rz(-0.59944154) q[2];
sx q[2];
rz(-1.2931461) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.030071478) q[1];
sx q[1];
rz(-2.3276405) q[1];
sx q[1];
rz(1.4210644) q[1];
rz(2.8844374) q[3];
sx q[3];
rz(-1.3793334) q[3];
sx q[3];
rz(1.934727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.19994911) q[2];
sx q[2];
rz(-1.4120833) q[2];
sx q[2];
rz(1.6597718) q[2];
rz(-3.0485349) q[3];
sx q[3];
rz(-1.2991644) q[3];
sx q[3];
rz(-0.53946462) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.900443) q[0];
sx q[0];
rz(-2.7261782) q[0];
sx q[0];
rz(2.7442617) q[0];
rz(0.98995248) q[1];
sx q[1];
rz(-1.3497458) q[1];
sx q[1];
rz(-2.1053402) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5118714) q[0];
sx q[0];
rz(-0.40490926) q[0];
sx q[0];
rz(0.029442336) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7099047) q[2];
sx q[2];
rz(-1.221162) q[2];
sx q[2];
rz(-1.0316499) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.73309988) q[1];
sx q[1];
rz(-1.6599791) q[1];
sx q[1];
rz(-0.21987889) q[1];
x q[2];
rz(2.0898706) q[3];
sx q[3];
rz(-1.241893) q[3];
sx q[3];
rz(2.2954706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1194666) q[2];
sx q[2];
rz(-2.4194722) q[2];
sx q[2];
rz(-1.3163346) q[2];
rz(-0.25257603) q[3];
sx q[3];
rz(-1.794869) q[3];
sx q[3];
rz(0.36309567) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0040358) q[0];
sx q[0];
rz(-2.3513849) q[0];
sx q[0];
rz(0.23605119) q[0];
rz(0.099418489) q[1];
sx q[1];
rz(-0.75629083) q[1];
sx q[1];
rz(-1.8831467) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6981068) q[0];
sx q[0];
rz(-2.5763566) q[0];
sx q[0];
rz(-0.14552571) q[0];
x q[1];
rz(-2.4781371) q[2];
sx q[2];
rz(-1.1642312) q[2];
sx q[2];
rz(2.6321509) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.26791456) q[1];
sx q[1];
rz(-0.89338747) q[1];
sx q[1];
rz(1.4235086) q[1];
rz(0.80983617) q[3];
sx q[3];
rz(-2.5812529) q[3];
sx q[3];
rz(-0.59691012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9749757) q[2];
sx q[2];
rz(-0.48305837) q[2];
sx q[2];
rz(2.100259) q[2];
rz(0.58639041) q[3];
sx q[3];
rz(-2.6643463) q[3];
sx q[3];
rz(2.3311116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7753684) q[0];
sx q[0];
rz(-1.0564221) q[0];
sx q[0];
rz(-1.139241) q[0];
rz(2.1495023) q[1];
sx q[1];
rz(-1.5403668) q[1];
sx q[1];
rz(-1.5425727) q[1];
rz(-2.8236961) q[2];
sx q[2];
rz(-1.3513202) q[2];
sx q[2];
rz(2.8693595) q[2];
rz(0.6057779) q[3];
sx q[3];
rz(-0.68690261) q[3];
sx q[3];
rz(1.5920832) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
