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
rz(-1.0322303) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55695483) q[0];
sx q[0];
rz(-0.9171066) q[0];
sx q[0];
rz(-1.6980706) q[0];
rz(-pi) q[1];
rz(0.039041877) q[2];
sx q[2];
rz(-2.0882975) q[2];
sx q[2];
rz(1.7797178) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.5121442) q[1];
sx q[1];
rz(-1.2595065) q[1];
sx q[1];
rz(-1.4427079) q[1];
x q[2];
rz(2.0809523) q[3];
sx q[3];
rz(-2.3054696) q[3];
sx q[3];
rz(1.5165129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.44077474) q[2];
sx q[2];
rz(-0.60443193) q[2];
sx q[2];
rz(-0.087892858) q[2];
rz(1.4676189) q[3];
sx q[3];
rz(-1.847495) q[3];
sx q[3];
rz(-2.1427593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1217594) q[0];
sx q[0];
rz(-1.8353945) q[0];
sx q[0];
rz(-1.649296) q[0];
rz(-0.78798931) q[1];
sx q[1];
rz(-1.9580656) q[1];
sx q[1];
rz(-0.96510395) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0168646) q[0];
sx q[0];
rz(-1.6627321) q[0];
sx q[0];
rz(2.4107736) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9847174) q[2];
sx q[2];
rz(-0.74293908) q[2];
sx q[2];
rz(2.3791831) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9763042) q[1];
sx q[1];
rz(-1.3852264) q[1];
sx q[1];
rz(2.0363505) q[1];
rz(-pi) q[2];
x q[2];
rz(0.177029) q[3];
sx q[3];
rz(-2.6654589) q[3];
sx q[3];
rz(3.0974922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.013177055) q[2];
sx q[2];
rz(-2.8030277) q[2];
sx q[2];
rz(1.3236807) q[2];
rz(2.2549818) q[3];
sx q[3];
rz(-1.9804136) q[3];
sx q[3];
rz(0.19865856) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.897268) q[0];
sx q[0];
rz(-0.97049814) q[0];
sx q[0];
rz(-0.6955198) q[0];
rz(-2.3059402) q[1];
sx q[1];
rz(-1.4202159) q[1];
sx q[1];
rz(2.4998891) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7135237) q[0];
sx q[0];
rz(-1.4260573) q[0];
sx q[0];
rz(2.3109396) q[0];
rz(1.4029866) q[2];
sx q[2];
rz(-0.8665167) q[2];
sx q[2];
rz(-1.2806127) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5776331) q[1];
sx q[1];
rz(-1.0924939) q[1];
sx q[1];
rz(-2.7863414) q[1];
rz(-pi) q[2];
x q[2];
rz(0.99567497) q[3];
sx q[3];
rz(-0.96700562) q[3];
sx q[3];
rz(-3.1310981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9782763) q[2];
sx q[2];
rz(-1.1288319) q[2];
sx q[2];
rz(-1.3345435) q[2];
rz(1.7415907) q[3];
sx q[3];
rz(-1.6440697) q[3];
sx q[3];
rz(-1.1388206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9982933) q[0];
sx q[0];
rz(-2.2704953) q[0];
sx q[0];
rz(0.77348462) q[0];
rz(3.0553014) q[1];
sx q[1];
rz(-2.6373865) q[1];
sx q[1];
rz(2.6533244) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9241656) q[0];
sx q[0];
rz(-2.0323951) q[0];
sx q[0];
rz(1.4946926) q[0];
rz(-1.0268698) q[2];
sx q[2];
rz(-2.1928117) q[2];
sx q[2];
rz(-2.029947) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.5544071) q[1];
sx q[1];
rz(-2.3766365) q[1];
sx q[1];
rz(1.8991963) q[1];
x q[2];
rz(1.1639488) q[3];
sx q[3];
rz(-1.9449807) q[3];
sx q[3];
rz(-2.0896623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.46828541) q[2];
sx q[2];
rz(-2.0279334) q[2];
sx q[2];
rz(-2.9735273) q[2];
rz(-2.7189861) q[3];
sx q[3];
rz(-0.97413078) q[3];
sx q[3];
rz(2.9882123) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3733658) q[0];
sx q[0];
rz(-2.234937) q[0];
sx q[0];
rz(1.1905097) q[0];
rz(-2.6137784) q[1];
sx q[1];
rz(-1.0820791) q[1];
sx q[1];
rz(-3.1324918) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9267488) q[0];
sx q[0];
rz(-1.1256477) q[0];
sx q[0];
rz(2.4657004) q[0];
rz(2.3528655) q[2];
sx q[2];
rz(-0.085914748) q[2];
sx q[2];
rz(-0.14363657) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8619662) q[1];
sx q[1];
rz(-1.1512966) q[1];
sx q[1];
rz(0.90200938) q[1];
rz(-2.1808257) q[3];
sx q[3];
rz(-0.72566635) q[3];
sx q[3];
rz(-2.8270367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8268383) q[2];
sx q[2];
rz(-0.49586168) q[2];
sx q[2];
rz(-0.66519386) q[2];
rz(-2.1151755) q[3];
sx q[3];
rz(-1.0746936) q[3];
sx q[3];
rz(-2.3608666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.5196395) q[0];
sx q[0];
rz(-0.53549796) q[0];
sx q[0];
rz(0.39805472) q[0];
rz(-1.4705426) q[1];
sx q[1];
rz(-2.2649951) q[1];
sx q[1];
rz(2.3614531) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3624293) q[0];
sx q[0];
rz(-2.0423959) q[0];
sx q[0];
rz(-2.7410024) q[0];
rz(1.3684811) q[2];
sx q[2];
rz(-1.8087808) q[2];
sx q[2];
rz(-1.1954824) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6717588) q[1];
sx q[1];
rz(-2.2398754) q[1];
sx q[1];
rz(-1.9356739) q[1];
rz(1.1885181) q[3];
sx q[3];
rz(-0.72186379) q[3];
sx q[3];
rz(-0.75306276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.857343) q[2];
sx q[2];
rz(-2.0390859) q[2];
sx q[2];
rz(-2.9798689) q[2];
rz(-3.0297847) q[3];
sx q[3];
rz(-2.4155152) q[3];
sx q[3];
rz(-2.4017754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1424471) q[0];
sx q[0];
rz(-1.4302) q[0];
sx q[0];
rz(-0.69851843) q[0];
rz(1.0922208) q[1];
sx q[1];
rz(-0.75463808) q[1];
sx q[1];
rz(3.0879367) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7258436) q[0];
sx q[0];
rz(-0.77139716) q[0];
sx q[0];
rz(1.1939069) q[0];
x q[1];
rz(-1.1321819) q[2];
sx q[2];
rz(-1.2234067) q[2];
sx q[2];
rz(-1.5127104) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7602348) q[1];
sx q[1];
rz(-1.6242199) q[1];
sx q[1];
rz(-2.670856) q[1];
rz(0.14999203) q[3];
sx q[3];
rz(-2.2120471) q[3];
sx q[3];
rz(2.7771726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1412389) q[2];
sx q[2];
rz(-2.443479) q[2];
sx q[2];
rz(-0.21256438) q[2];
rz(0.32771787) q[3];
sx q[3];
rz(-2.1543584) q[3];
sx q[3];
rz(1.3535961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89606365) q[0];
sx q[0];
rz(-1.117417) q[0];
sx q[0];
rz(-2.8118706) q[0];
rz(0.7715191) q[1];
sx q[1];
rz(-2.3604269) q[1];
sx q[1];
rz(-0.82569295) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3701009) q[0];
sx q[0];
rz(-0.90898642) q[0];
sx q[0];
rz(-1.460381) q[0];
rz(-pi) q[1];
rz(1.6861028) q[2];
sx q[2];
rz(-0.59944154) q[2];
sx q[2];
rz(-1.8484465) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.18620423) q[1];
sx q[1];
rz(-2.3729747) q[1];
sx q[1];
rz(-0.15665084) q[1];
rz(1.7686033) q[3];
sx q[3];
rz(-1.8231492) q[3];
sx q[3];
rz(-0.41393241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.19994911) q[2];
sx q[2];
rz(-1.4120833) q[2];
sx q[2];
rz(1.6597718) q[2];
rz(3.0485349) q[3];
sx q[3];
rz(-1.2991644) q[3];
sx q[3];
rz(-2.602128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-2.1516402) q[1];
sx q[1];
rz(-1.7918469) q[1];
sx q[1];
rz(-1.0362524) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2275804) q[0];
sx q[0];
rz(-1.5591994) q[0];
sx q[0];
rz(2.7368403) q[0];
rz(-pi) q[1];
x q[1];
rz(0.36340275) q[2];
sx q[2];
rz(-2.7663603) q[2];
sx q[2];
rz(-2.4979532) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9246722) q[1];
sx q[1];
rz(-2.9045871) q[1];
sx q[1];
rz(2.7525178) q[1];
rz(-pi) q[2];
x q[2];
rz(0.37452368) q[3];
sx q[3];
rz(-1.0820884) q[3];
sx q[3];
rz(-0.54218369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1194666) q[2];
sx q[2];
rz(-2.4194722) q[2];
sx q[2];
rz(-1.3163346) q[2];
rz(-2.8890166) q[3];
sx q[3];
rz(-1.3467237) q[3];
sx q[3];
rz(0.36309567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13755688) q[0];
sx q[0];
rz(-2.3513849) q[0];
sx q[0];
rz(0.23605119) q[0];
rz(-3.0421742) q[1];
sx q[1];
rz(-2.3853018) q[1];
sx q[1];
rz(1.8831467) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44348587) q[0];
sx q[0];
rz(-2.5763566) q[0];
sx q[0];
rz(0.14552571) q[0];
rz(-pi) q[1];
rz(-0.61014097) q[2];
sx q[2];
rz(-2.379866) q[2];
sx q[2];
rz(-1.6118647) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.035485766) q[1];
sx q[1];
rz(-0.69075876) q[1];
sx q[1];
rz(0.18045998) q[1];
x q[2];
rz(-1.997273) q[3];
sx q[3];
rz(-1.1955373) q[3];
sx q[3];
rz(-0.29510185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1666169) q[2];
sx q[2];
rz(-2.6585343) q[2];
sx q[2];
rz(-2.100259) q[2];
rz(-0.58639041) q[3];
sx q[3];
rz(-2.6643463) q[3];
sx q[3];
rz(-2.3311116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36622421) q[0];
sx q[0];
rz(-2.0851705) q[0];
sx q[0];
rz(2.0023517) q[0];
rz(-0.99209039) q[1];
sx q[1];
rz(-1.5403668) q[1];
sx q[1];
rz(-1.5425727) q[1];
rz(0.31789657) q[2];
sx q[2];
rz(-1.3513202) q[2];
sx q[2];
rz(2.8693595) q[2];
rz(1.1339034) q[3];
sx q[3];
rz(-1.0224167) q[3];
sx q[3];
rz(-2.2800048) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
