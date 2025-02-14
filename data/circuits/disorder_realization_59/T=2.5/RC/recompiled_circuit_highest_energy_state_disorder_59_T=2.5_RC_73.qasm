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
rz(-0.72379273) q[0];
sx q[0];
rz(4.2552636) q[0];
sx q[0];
rz(10.389513) q[0];
rz(-0.073702987) q[1];
sx q[1];
rz(-0.6414203) q[1];
sx q[1];
rz(-1.4944271) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7429356) q[0];
sx q[0];
rz(-1.1857067) q[0];
sx q[0];
rz(2.4214817) q[0];
x q[1];
rz(-0.10434465) q[2];
sx q[2];
rz(-1.2547224) q[2];
sx q[2];
rz(2.2869273) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9656417) q[1];
sx q[1];
rz(-2.0276105) q[1];
sx q[1];
rz(0.58689249) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0953856) q[3];
sx q[3];
rz(-1.0418834) q[3];
sx q[3];
rz(-2.5387704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.73704314) q[2];
sx q[2];
rz(-1.5207542) q[2];
sx q[2];
rz(0.21657319) q[2];
rz(-0.51566044) q[3];
sx q[3];
rz(-1.0786846) q[3];
sx q[3];
rz(3.0576341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11230042) q[0];
sx q[0];
rz(-3.0700505) q[0];
sx q[0];
rz(1.4805502) q[0];
rz(-2.5484565) q[1];
sx q[1];
rz(-0.99855223) q[1];
sx q[1];
rz(0.28233972) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38587785) q[0];
sx q[0];
rz(-0.1308724) q[0];
sx q[0];
rz(-0.22481052) q[0];
rz(-pi) q[1];
x q[1];
rz(3.00782) q[2];
sx q[2];
rz(-1.3743128) q[2];
sx q[2];
rz(-1.6272194) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8806539) q[1];
sx q[1];
rz(-2.1092514) q[1];
sx q[1];
rz(1.7314584) q[1];
x q[2];
rz(0.37455885) q[3];
sx q[3];
rz(-1.1653596) q[3];
sx q[3];
rz(0.64664155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4234408) q[2];
sx q[2];
rz(-0.56090063) q[2];
sx q[2];
rz(2.1924428) q[2];
rz(-0.078350457) q[3];
sx q[3];
rz(-1.6716985) q[3];
sx q[3];
rz(-1.2300434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(-2.3345728) q[0];
sx q[0];
rz(-0.0581352) q[0];
sx q[0];
rz(2.9491501) q[0];
rz(-2.1678534) q[1];
sx q[1];
rz(-1.0304008) q[1];
sx q[1];
rz(-1.2724426) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5010492) q[0];
sx q[0];
rz(-1.328381) q[0];
sx q[0];
rz(0.69037504) q[0];
rz(1.3272304) q[2];
sx q[2];
rz(-2.3911016) q[2];
sx q[2];
rz(1.0433973) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1509708) q[1];
sx q[1];
rz(-0.34268296) q[1];
sx q[1];
rz(-2.3938378) q[1];
rz(1.4951823) q[3];
sx q[3];
rz(-2.6845884) q[3];
sx q[3];
rz(1.5357032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5729313) q[2];
sx q[2];
rz(-2.2243786) q[2];
sx q[2];
rz(2.5407963) q[2];
rz(0.83493799) q[3];
sx q[3];
rz(-0.91885126) q[3];
sx q[3];
rz(0.44343534) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3804974) q[0];
sx q[0];
rz(-2.0090071) q[0];
sx q[0];
rz(-1.4196716) q[0];
rz(-2.8555866) q[1];
sx q[1];
rz(-2.0068469) q[1];
sx q[1];
rz(2.6223415) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4409281) q[0];
sx q[0];
rz(-1.6608149) q[0];
sx q[0];
rz(0.16914455) q[0];
rz(-pi) q[1];
rz(-0.51651603) q[2];
sx q[2];
rz(-0.20843796) q[2];
sx q[2];
rz(1.0112273) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1468526) q[1];
sx q[1];
rz(-2.6255872) q[1];
sx q[1];
rz(-0.13928646) q[1];
rz(-1.6029699) q[3];
sx q[3];
rz(-1.0250499) q[3];
sx q[3];
rz(1.4226899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.1661735) q[2];
sx q[2];
rz(-0.68048802) q[2];
sx q[2];
rz(-0.55230459) q[2];
rz(2.35899) q[3];
sx q[3];
rz(-1.3026594) q[3];
sx q[3];
rz(2.9569614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(0.85823378) q[0];
sx q[0];
rz(-0.62223804) q[0];
sx q[0];
rz(1.8587814) q[0];
rz(1.2239617) q[1];
sx q[1];
rz(-1.6061648) q[1];
sx q[1];
rz(2.431869) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8741651) q[0];
sx q[0];
rz(-1.5418959) q[0];
sx q[0];
rz(-1.5927107) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.100509) q[2];
sx q[2];
rz(-2.6729995) q[2];
sx q[2];
rz(3.1236609) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.97204527) q[1];
sx q[1];
rz(-1.2575713) q[1];
sx q[1];
rz(0.07446988) q[1];
rz(-1.9842667) q[3];
sx q[3];
rz(-1.7352823) q[3];
sx q[3];
rz(2.2604347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2964786) q[2];
sx q[2];
rz(-0.60135403) q[2];
sx q[2];
rz(-0.55310407) q[2];
rz(-2.3002355) q[3];
sx q[3];
rz(-2.0935121) q[3];
sx q[3];
rz(-1.1355737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-0.76310054) q[0];
sx q[0];
rz(-1.0569514) q[0];
sx q[0];
rz(0.77504778) q[0];
rz(1.4729602) q[1];
sx q[1];
rz(-2.5170363) q[1];
sx q[1];
rz(2.3042309) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96356364) q[0];
sx q[0];
rz(-2.508381) q[0];
sx q[0];
rz(-1.7342315) q[0];
x q[1];
rz(2.6858575) q[2];
sx q[2];
rz(-1.4066469) q[2];
sx q[2];
rz(1.8388621) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.130116) q[1];
sx q[1];
rz(-1.1728047) q[1];
sx q[1];
rz(0.46020646) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9894137) q[3];
sx q[3];
rz(-1.8004724) q[3];
sx q[3];
rz(-0.21462378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2396635) q[2];
sx q[2];
rz(-1.9749125) q[2];
sx q[2];
rz(-2.3789294) q[2];
rz(2.934382) q[3];
sx q[3];
rz(-1.8800294) q[3];
sx q[3];
rz(-2.5337849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36829456) q[0];
sx q[0];
rz(-0.93450707) q[0];
sx q[0];
rz(-1.6873129) q[0];
rz(1.2794718) q[1];
sx q[1];
rz(-0.8431294) q[1];
sx q[1];
rz(1.4131193) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2103677) q[0];
sx q[0];
rz(-0.43950167) q[0];
sx q[0];
rz(2.6511433) q[0];
rz(-pi) q[1];
rz(2.3033875) q[2];
sx q[2];
rz(-2.0186485) q[2];
sx q[2];
rz(0.034016646) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.6225902) q[1];
sx q[1];
rz(-2.0068002) q[1];
sx q[1];
rz(-1.2335797) q[1];
rz(-pi) q[2];
x q[2];
rz(0.90274324) q[3];
sx q[3];
rz(-0.53865005) q[3];
sx q[3];
rz(-0.99942943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1450682) q[2];
sx q[2];
rz(-0.86864305) q[2];
sx q[2];
rz(-0.16150148) q[2];
rz(-0.7343556) q[3];
sx q[3];
rz(-0.86638325) q[3];
sx q[3];
rz(1.8374779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9228009) q[0];
sx q[0];
rz(-0.53023338) q[0];
sx q[0];
rz(-0.23736048) q[0];
rz(-2.3573719) q[1];
sx q[1];
rz(-1.3619245) q[1];
sx q[1];
rz(1.4201737) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.500017) q[0];
sx q[0];
rz(-0.58860129) q[0];
sx q[0];
rz(0.55843784) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2591822) q[2];
sx q[2];
rz(-1.3974481) q[2];
sx q[2];
rz(-1.2899931) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3846459) q[1];
sx q[1];
rz(-1.8638041) q[1];
sx q[1];
rz(-0.56045597) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3444207) q[3];
sx q[3];
rz(-2.4690921) q[3];
sx q[3];
rz(-0.029880015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.64804849) q[2];
sx q[2];
rz(-1.8914787) q[2];
sx q[2];
rz(-1.2250712) q[2];
rz(-1.3263634) q[3];
sx q[3];
rz(-2.1882961) q[3];
sx q[3];
rz(-1.1939322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61106435) q[0];
sx q[0];
rz(-3.0623797) q[0];
sx q[0];
rz(-2.5031669) q[0];
rz(3.0910659) q[1];
sx q[1];
rz(-1.9332998) q[1];
sx q[1];
rz(-0.60727492) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80562185) q[0];
sx q[0];
rz(-2.6329649) q[0];
sx q[0];
rz(2.2679099) q[0];
x q[1];
rz(-0.061640306) q[2];
sx q[2];
rz(-2.2547005) q[2];
sx q[2];
rz(-2.7126171) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6029507) q[1];
sx q[1];
rz(-1.2455229) q[1];
sx q[1];
rz(3.0673749) q[1];
rz(-pi) q[2];
rz(-1.7631986) q[3];
sx q[3];
rz(-1.3891313) q[3];
sx q[3];
rz(0.84027973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4813469) q[2];
sx q[2];
rz(-1.1911743) q[2];
sx q[2];
rz(0.78835362) q[2];
rz(-2.3039019) q[3];
sx q[3];
rz(-0.82101429) q[3];
sx q[3];
rz(1.2825509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4072708) q[0];
sx q[0];
rz(-0.73198524) q[0];
sx q[0];
rz(-0.5087854) q[0];
rz(-1.58163) q[1];
sx q[1];
rz(-1.0204693) q[1];
sx q[1];
rz(0.4090974) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74241766) q[0];
sx q[0];
rz(-1.3361201) q[0];
sx q[0];
rz(0.15813903) q[0];
x q[1];
rz(0.81014388) q[2];
sx q[2];
rz(-0.64254566) q[2];
sx q[2];
rz(2.0755529) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0171788) q[1];
sx q[1];
rz(-0.22089566) q[1];
sx q[1];
rz(3.0934889) q[1];
rz(-pi) q[2];
rz(0.38842885) q[3];
sx q[3];
rz(-1.9216537) q[3];
sx q[3];
rz(2.4728123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6009377) q[2];
sx q[2];
rz(-1.6161852) q[2];
sx q[2];
rz(1.8589004) q[2];
rz(0.013785275) q[3];
sx q[3];
rz(-2.0163592) q[3];
sx q[3];
rz(-1.3330601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8448821) q[0];
sx q[0];
rz(-1.409335) q[0];
sx q[0];
rz(0.62406337) q[0];
rz(-0.86112549) q[1];
sx q[1];
rz(-2.594941) q[1];
sx q[1];
rz(-3.0088967) q[1];
rz(-2.843905) q[2];
sx q[2];
rz(-0.33391446) q[2];
sx q[2];
rz(-3.1261974) q[2];
rz(0.63152159) q[3];
sx q[3];
rz(-2.1248264) q[3];
sx q[3];
rz(-0.2223224) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
