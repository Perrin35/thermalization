OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.57269078) q[0];
sx q[0];
rz(-2.1532018) q[0];
sx q[0];
rz(0.44901499) q[0];
rz(-0.16945893) q[1];
sx q[1];
rz(-3.0100477) q[1];
sx q[1];
rz(-2.0102672) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22085359) q[0];
sx q[0];
rz(-1.8628696) q[0];
sx q[0];
rz(2.1983485) q[0];
rz(-pi) q[1];
rz(-2.4822818) q[2];
sx q[2];
rz(-0.74987312) q[2];
sx q[2];
rz(-1.4091968) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3851812) q[1];
sx q[1];
rz(-2.1527421) q[1];
sx q[1];
rz(0.20828282) q[1];
x q[2];
rz(-1.851397) q[3];
sx q[3];
rz(-0.56125703) q[3];
sx q[3];
rz(0.11229501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1863056) q[2];
sx q[2];
rz(-0.39262843) q[2];
sx q[2];
rz(0.10360959) q[2];
rz(2.7747532) q[3];
sx q[3];
rz(-1.6678526) q[3];
sx q[3];
rz(1.8240671) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1845301) q[0];
sx q[0];
rz(-0.4758895) q[0];
sx q[0];
rz(2.6153508) q[0];
rz(1.8980252) q[1];
sx q[1];
rz(-1.4105816) q[1];
sx q[1];
rz(-0.19613656) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1271106) q[0];
sx q[0];
rz(-2.0271993) q[0];
sx q[0];
rz(3.0672202) q[0];
rz(-1.4043442) q[2];
sx q[2];
rz(-0.86814082) q[2];
sx q[2];
rz(0.26674092) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.89333234) q[1];
sx q[1];
rz(-2.6702017) q[1];
sx q[1];
rz(1.8880647) q[1];
rz(-pi) q[2];
rz(-2.4200581) q[3];
sx q[3];
rz(-0.01577687) q[3];
sx q[3];
rz(1.1392913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5940932) q[2];
sx q[2];
rz(-0.27442351) q[2];
sx q[2];
rz(0.37626949) q[2];
rz(-0.88092342) q[3];
sx q[3];
rz(-1.8062402) q[3];
sx q[3];
rz(2.3295565) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6905717) q[0];
sx q[0];
rz(-2.8132827) q[0];
sx q[0];
rz(2.1300533) q[0];
rz(-0.66863376) q[1];
sx q[1];
rz(-1.1625544) q[1];
sx q[1];
rz(-2.5943601) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1062463) q[0];
sx q[0];
rz(-1.5515455) q[0];
sx q[0];
rz(0.20142844) q[0];
x q[1];
rz(1.9220566) q[2];
sx q[2];
rz(-2.329956) q[2];
sx q[2];
rz(1.7691624) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8142952) q[1];
sx q[1];
rz(-1.3298099) q[1];
sx q[1];
rz(2.2404033) q[1];
rz(-2.3272334) q[3];
sx q[3];
rz(-1.7118771) q[3];
sx q[3];
rz(-2.1997978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.68613595) q[2];
sx q[2];
rz(-0.28527173) q[2];
sx q[2];
rz(1.6395462) q[2];
rz(-0.57602588) q[3];
sx q[3];
rz(-1.4833769) q[3];
sx q[3];
rz(-0.36147931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6308924) q[0];
sx q[0];
rz(-0.80131131) q[0];
sx q[0];
rz(-2.8019688) q[0];
rz(2.6009808) q[1];
sx q[1];
rz(-2.4410591) q[1];
sx q[1];
rz(2.2115754) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1620063) q[0];
sx q[0];
rz(-1.7135059) q[0];
sx q[0];
rz(-3.0703785) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5901718) q[2];
sx q[2];
rz(-0.65931554) q[2];
sx q[2];
rz(2.2267411) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2423289) q[1];
sx q[1];
rz(-1.5771958) q[1];
sx q[1];
rz(2.0501627) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.633938) q[3];
sx q[3];
rz(-2.1498907) q[3];
sx q[3];
rz(-2.7217025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.12863079) q[2];
sx q[2];
rz(-1.4350812) q[2];
sx q[2];
rz(1.5738515) q[2];
rz(0.74357998) q[3];
sx q[3];
rz(-0.79135528) q[3];
sx q[3];
rz(-0.96405205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6745233) q[0];
sx q[0];
rz(-0.85551298) q[0];
sx q[0];
rz(-3.0849482) q[0];
rz(-1.4757587) q[1];
sx q[1];
rz(-2.00878) q[1];
sx q[1];
rz(1.3585565) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3751793) q[0];
sx q[0];
rz(-2.4046899) q[0];
sx q[0];
rz(0.76216682) q[0];
rz(-0.26014056) q[2];
sx q[2];
rz(-0.67408757) q[2];
sx q[2];
rz(-2.0345794) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.92734071) q[1];
sx q[1];
rz(-2.3433629) q[1];
sx q[1];
rz(-2.4386474) q[1];
x q[2];
rz(-1.6886466) q[3];
sx q[3];
rz(-1.0123583) q[3];
sx q[3];
rz(-1.8359566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8654827) q[2];
sx q[2];
rz(-1.7296187) q[2];
sx q[2];
rz(-1.126368) q[2];
rz(-1.3364835) q[3];
sx q[3];
rz(-2.0650605) q[3];
sx q[3];
rz(-1.0895464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4514076) q[0];
sx q[0];
rz(-1.0253588) q[0];
sx q[0];
rz(-3.0080646) q[0];
rz(2.1639157) q[1];
sx q[1];
rz(-2.1656499) q[1];
sx q[1];
rz(2.8547063) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3505976) q[0];
sx q[0];
rz(-0.7028326) q[0];
sx q[0];
rz(-1.4969456) q[0];
rz(-1.8989766) q[2];
sx q[2];
rz(-1.6115453) q[2];
sx q[2];
rz(-2.2475257) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9911954) q[1];
sx q[1];
rz(-2.5270473) q[1];
sx q[1];
rz(-1.9646364) q[1];
rz(-0.37128504) q[3];
sx q[3];
rz(-0.6273191) q[3];
sx q[3];
rz(0.016022779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4213244) q[2];
sx q[2];
rz(-1.4807533) q[2];
sx q[2];
rz(1.486091) q[2];
rz(2.7739575) q[3];
sx q[3];
rz(-1.0272107) q[3];
sx q[3];
rz(1.0953085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7169645) q[0];
sx q[0];
rz(-2.3724738) q[0];
sx q[0];
rz(-2.6677483) q[0];
rz(2.4773856) q[1];
sx q[1];
rz(-1.217548) q[1];
sx q[1];
rz(-1.7582105) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5055588) q[0];
sx q[0];
rz(-2.2657388) q[0];
sx q[0];
rz(0.64993422) q[0];
x q[1];
rz(1.6331312) q[2];
sx q[2];
rz(-0.41878715) q[2];
sx q[2];
rz(0.58277786) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3480839) q[1];
sx q[1];
rz(-1.3474819) q[1];
sx q[1];
rz(2.9053754) q[1];
x q[2];
rz(0.9035191) q[3];
sx q[3];
rz(-2.857803) q[3];
sx q[3];
rz(-1.9061778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.15709269) q[2];
sx q[2];
rz(-1.7326771) q[2];
sx q[2];
rz(-2.816693) q[2];
rz(2.7609008) q[3];
sx q[3];
rz(-0.84563962) q[3];
sx q[3];
rz(-1.0977753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0710881) q[0];
sx q[0];
rz(-0.66362137) q[0];
sx q[0];
rz(-0.93884236) q[0];
rz(-0.70010575) q[1];
sx q[1];
rz(-2.5036) q[1];
sx q[1];
rz(0.28900388) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5225713) q[0];
sx q[0];
rz(-2.8998525) q[0];
sx q[0];
rz(-1.8999343) q[0];
x q[1];
rz(-2.1949057) q[2];
sx q[2];
rz(-0.74926976) q[2];
sx q[2];
rz(2.5890337) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.49456143) q[1];
sx q[1];
rz(-2.6115578) q[1];
sx q[1];
rz(-1.3840335) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4988171) q[3];
sx q[3];
rz(-1.9508617) q[3];
sx q[3];
rz(-0.38004181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2603904) q[2];
sx q[2];
rz(-1.6198879) q[2];
sx q[2];
rz(0.41268665) q[2];
rz(2.2212501) q[3];
sx q[3];
rz(-2.6350382) q[3];
sx q[3];
rz(-1.2296366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5145787) q[0];
sx q[0];
rz(-0.98015061) q[0];
sx q[0];
rz(0.29943109) q[0];
rz(2.0191655) q[1];
sx q[1];
rz(-1.2010682) q[1];
sx q[1];
rz(-2.1103512) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5297197) q[0];
sx q[0];
rz(-2.7752989) q[0];
sx q[0];
rz(-1.6008911) q[0];
rz(-0.38487969) q[2];
sx q[2];
rz(-0.49321929) q[2];
sx q[2];
rz(2.9259591) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.37027994) q[1];
sx q[1];
rz(-2.2590933) q[1];
sx q[1];
rz(-2.7629333) q[1];
x q[2];
rz(-2.5599285) q[3];
sx q[3];
rz(-2.1278893) q[3];
sx q[3];
rz(-0.34639369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9791744) q[2];
sx q[2];
rz(-1.7981497) q[2];
sx q[2];
rz(-0.24270414) q[2];
rz(1.3001214) q[3];
sx q[3];
rz(-2.0827115) q[3];
sx q[3];
rz(0.21568957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8278787) q[0];
sx q[0];
rz(-0.21610459) q[0];
sx q[0];
rz(-1.4208273) q[0];
rz(-2.8315262) q[1];
sx q[1];
rz(-2.2646751) q[1];
sx q[1];
rz(-1.5975331) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2313854) q[0];
sx q[0];
rz(-0.69150309) q[0];
sx q[0];
rz(-0.67759902) q[0];
rz(-pi) q[1];
rz(1.9126911) q[2];
sx q[2];
rz(-1.6544621) q[2];
sx q[2];
rz(-0.25632206) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.76031715) q[1];
sx q[1];
rz(-0.3418498) q[1];
sx q[1];
rz(-2.443497) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9821854) q[3];
sx q[3];
rz(-1.4324354) q[3];
sx q[3];
rz(-2.0970999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1342423) q[2];
sx q[2];
rz(-1.4126567) q[2];
sx q[2];
rz(1.3664112) q[2];
rz(-0.8775231) q[3];
sx q[3];
rz(-1.6759422) q[3];
sx q[3];
rz(0.88959488) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9919745) q[0];
sx q[0];
rz(-1.5852954) q[0];
sx q[0];
rz(-1.4854767) q[0];
rz(0.86434518) q[1];
sx q[1];
rz(-1.0135916) q[1];
sx q[1];
rz(2.7085173) q[1];
rz(-0.50044671) q[2];
sx q[2];
rz(-1.3896349) q[2];
sx q[2];
rz(-1.096772) q[2];
rz(-1.6771562) q[3];
sx q[3];
rz(-0.48135664) q[3];
sx q[3];
rz(2.4873747) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
