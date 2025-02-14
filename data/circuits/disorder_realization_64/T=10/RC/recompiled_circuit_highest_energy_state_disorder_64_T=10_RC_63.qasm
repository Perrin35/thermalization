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
rz(4.4216006) q[0];
sx q[0];
rz(13.292424) q[0];
rz(-2.5441406) q[1];
sx q[1];
rz(-0.51550454) q[1];
sx q[1];
rz(2.1093624) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55695483) q[0];
sx q[0];
rz(-0.9171066) q[0];
sx q[0];
rz(-1.6980706) q[0];
rz(-1.502336) q[2];
sx q[2];
rz(-0.51883829) q[2];
sx q[2];
rz(1.8585132) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1223678) q[1];
sx q[1];
rz(-1.4488954) q[1];
sx q[1];
rz(2.8278964) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.80268245) q[3];
sx q[3];
rz(-1.2000053) q[3];
sx q[3];
rz(-0.30457531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.44077474) q[2];
sx q[2];
rz(-0.60443193) q[2];
sx q[2];
rz(0.087892858) q[2];
rz(1.6739738) q[3];
sx q[3];
rz(-1.847495) q[3];
sx q[3];
rz(2.1427593) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1217594) q[0];
sx q[0];
rz(-1.8353945) q[0];
sx q[0];
rz(1.649296) q[0];
rz(-0.78798931) q[1];
sx q[1];
rz(-1.183527) q[1];
sx q[1];
rz(0.96510395) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3639587) q[0];
sx q[0];
rz(-0.84375644) q[0];
sx q[0];
rz(-1.6939837) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4048058) q[2];
sx q[2];
rz(-1.6766785) q[2];
sx q[2];
rz(-2.4491765) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3841159) q[1];
sx q[1];
rz(-0.49864492) q[1];
sx q[1];
rz(1.1747141) q[1];
rz(-1.6613668) q[3];
sx q[3];
rz(-1.1027059) q[3];
sx q[3];
rz(-0.15453574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.013177055) q[2];
sx q[2];
rz(-2.8030277) q[2];
sx q[2];
rz(1.3236807) q[2];
rz(0.88661083) q[3];
sx q[3];
rz(-1.9804136) q[3];
sx q[3];
rz(2.9429341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.897268) q[0];
sx q[0];
rz(-2.1710945) q[0];
sx q[0];
rz(-2.4460728) q[0];
rz(0.8356525) q[1];
sx q[1];
rz(-1.7213768) q[1];
sx q[1];
rz(0.64170352) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27370975) q[0];
sx q[0];
rz(-0.84016227) q[0];
sx q[0];
rz(2.9466936) q[0];
rz(-pi) q[1];
rz(1.738606) q[2];
sx q[2];
rz(-0.8665167) q[2];
sx q[2];
rz(-1.86098) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5776331) q[1];
sx q[1];
rz(-2.0490987) q[1];
sx q[1];
rz(2.7863414) q[1];
rz(-0.99567497) q[3];
sx q[3];
rz(-0.96700562) q[3];
sx q[3];
rz(3.1310981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9782763) q[2];
sx q[2];
rz(-2.0127608) q[2];
sx q[2];
rz(1.8070492) q[2];
rz(-1.7415907) q[3];
sx q[3];
rz(-1.6440697) q[3];
sx q[3];
rz(-2.002772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9982933) q[0];
sx q[0];
rz(-0.87109733) q[0];
sx q[0];
rz(0.77348462) q[0];
rz(-0.086291226) q[1];
sx q[1];
rz(-2.6373865) q[1];
sx q[1];
rz(2.6533244) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9241656) q[0];
sx q[0];
rz(-1.1091976) q[0];
sx q[0];
rz(1.6469) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1147229) q[2];
sx q[2];
rz(-0.94878093) q[2];
sx q[2];
rz(-1.1116456) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.5544071) q[1];
sx q[1];
rz(-2.3766365) q[1];
sx q[1];
rz(1.8991963) q[1];
x q[2];
rz(-2.7375324) q[3];
sx q[3];
rz(-1.9480138) q[3];
sx q[3];
rz(-0.67507832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.46828541) q[2];
sx q[2];
rz(-2.0279334) q[2];
sx q[2];
rz(-0.16806531) q[2];
rz(-2.7189861) q[3];
sx q[3];
rz(-0.97413078) q[3];
sx q[3];
rz(2.9882123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76822686) q[0];
sx q[0];
rz(-0.90665561) q[0];
sx q[0];
rz(-1.9510829) q[0];
rz(-0.52781421) q[1];
sx q[1];
rz(-2.0595136) q[1];
sx q[1];
rz(-3.1324918) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4531698) q[0];
sx q[0];
rz(-0.97081447) q[0];
sx q[0];
rz(-1.021941) q[0];
rz(-3.0809693) q[2];
sx q[2];
rz(-1.5098808) q[2];
sx q[2];
rz(0.9275533) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3745388) q[1];
sx q[1];
rz(-0.77195814) q[1];
sx q[1];
rz(-0.94725398) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1808257) q[3];
sx q[3];
rz(-2.4159263) q[3];
sx q[3];
rz(0.31455597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3147543) q[2];
sx q[2];
rz(-2.645731) q[2];
sx q[2];
rz(-2.4763988) q[2];
rz(-2.1151755) q[3];
sx q[3];
rz(-1.0746936) q[3];
sx q[3];
rz(-2.3608666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(1.5196395) q[0];
sx q[0];
rz(-2.6060947) q[0];
sx q[0];
rz(2.7435379) q[0];
rz(1.4705426) q[1];
sx q[1];
rz(-0.87659756) q[1];
sx q[1];
rz(2.3614531) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0288643) q[0];
sx q[0];
rz(-2.5327589) q[0];
sx q[0];
rz(-0.91797773) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7731116) q[2];
sx q[2];
rz(-1.8087808) q[2];
sx q[2];
rz(-1.9461103) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.91785039) q[1];
sx q[1];
rz(-0.7484352) q[1];
sx q[1];
rz(2.7176932) q[1];
rz(-1.9530746) q[3];
sx q[3];
rz(-2.4197289) q[3];
sx q[3];
rz(-2.3885299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.857343) q[2];
sx q[2];
rz(-1.1025068) q[2];
sx q[2];
rz(-0.16172376) q[2];
rz(3.0297847) q[3];
sx q[3];
rz(-0.72607741) q[3];
sx q[3];
rz(-2.4017754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9991456) q[0];
sx q[0];
rz(-1.4302) q[0];
sx q[0];
rz(2.4430742) q[0];
rz(1.0922208) q[1];
sx q[1];
rz(-0.75463808) q[1];
sx q[1];
rz(3.0879367) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41574902) q[0];
sx q[0];
rz(-0.77139716) q[0];
sx q[0];
rz(-1.9476858) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.8648134) q[2];
sx q[2];
rz(-2.5891782) q[2];
sx q[2];
rz(-0.68589003) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0567017) q[1];
sx q[1];
rz(-2.6680601) q[1];
sx q[1];
rz(-0.11736203) q[1];
x q[2];
rz(2.9916006) q[3];
sx q[3];
rz(-0.92954554) q[3];
sx q[3];
rz(2.7771726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1412389) q[2];
sx q[2];
rz(-2.443479) q[2];
sx q[2];
rz(2.9290283) q[2];
rz(-2.8138748) q[3];
sx q[3];
rz(-2.1543584) q[3];
sx q[3];
rz(1.3535961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.245529) q[0];
sx q[0];
rz(-1.117417) q[0];
sx q[0];
rz(0.32972202) q[0];
rz(2.3700736) q[1];
sx q[1];
rz(-0.78116575) q[1];
sx q[1];
rz(-0.82569295) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86733183) q[0];
sx q[0];
rz(-1.4837588) q[0];
sx q[0];
rz(-0.66477338) q[0];
rz(-0.078455047) q[2];
sx q[2];
rz(-0.9758853) q[2];
sx q[2];
rz(-1.9878146) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9553884) q[1];
sx q[1];
rz(-2.3729747) q[1];
sx q[1];
rz(-0.15665084) q[1];
rz(-pi) q[2];
rz(2.8844374) q[3];
sx q[3];
rz(-1.3793334) q[3];
sx q[3];
rz(-1.2068656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.19994911) q[2];
sx q[2];
rz(-1.4120833) q[2];
sx q[2];
rz(1.4818209) q[2];
rz(0.09305772) q[3];
sx q[3];
rz(-1.8424282) q[3];
sx q[3];
rz(-2.602128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2411497) q[0];
sx q[0];
rz(-0.41541442) q[0];
sx q[0];
rz(-2.7442617) q[0];
rz(2.1516402) q[1];
sx q[1];
rz(-1.3497458) q[1];
sx q[1];
rz(-1.0362524) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91401228) q[0];
sx q[0];
rz(-1.5591994) q[0];
sx q[0];
rz(-2.7368403) q[0];
rz(-pi) q[1];
rz(2.7781899) q[2];
sx q[2];
rz(-2.7663603) q[2];
sx q[2];
rz(2.4979532) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.73309988) q[1];
sx q[1];
rz(-1.6599791) q[1];
sx q[1];
rz(2.9217138) q[1];
x q[2];
rz(-0.96817044) q[3];
sx q[3];
rz(-2.535248) q[3];
sx q[3];
rz(1.9023638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.1194666) q[2];
sx q[2];
rz(-2.4194722) q[2];
sx q[2];
rz(1.3163346) q[2];
rz(0.25257603) q[3];
sx q[3];
rz(-1.794869) q[3];
sx q[3];
rz(-0.36309567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0040358) q[0];
sx q[0];
rz(-2.3513849) q[0];
sx q[0];
rz(2.9055415) q[0];
rz(3.0421742) q[1];
sx q[1];
rz(-0.75629083) q[1];
sx q[1];
rz(1.8831467) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1374219) q[0];
sx q[0];
rz(-1.4930471) q[0];
sx q[0];
rz(0.56044062) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5314517) q[2];
sx q[2];
rz(-0.76172667) q[2];
sx q[2];
rz(-1.6118647) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8736781) q[1];
sx q[1];
rz(-0.89338747) q[1];
sx q[1];
rz(-1.4235086) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.40835898) q[3];
sx q[3];
rz(-1.9658739) q[3];
sx q[3];
rz(-1.7008821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1666169) q[2];
sx q[2];
rz(-0.48305837) q[2];
sx q[2];
rz(-1.0413337) q[2];
rz(-0.58639041) q[3];
sx q[3];
rz(-0.47724637) q[3];
sx q[3];
rz(-0.81048107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7753684) q[0];
sx q[0];
rz(-1.0564221) q[0];
sx q[0];
rz(-1.139241) q[0];
rz(0.99209039) q[1];
sx q[1];
rz(-1.6012259) q[1];
sx q[1];
rz(1.59902) q[1];
rz(1.3401399) q[2];
sx q[2];
rz(-1.2607831) q[2];
sx q[2];
rz(1.2270437) q[2];
rz(2.0076892) q[3];
sx q[3];
rz(-2.1191759) q[3];
sx q[3];
rz(0.86158781) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
