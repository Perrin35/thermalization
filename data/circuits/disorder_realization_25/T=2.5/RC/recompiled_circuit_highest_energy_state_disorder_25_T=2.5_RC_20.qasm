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
rz(-1.1005041) q[0];
sx q[0];
rz(-2.4325844) q[0];
sx q[0];
rz(-0.79431835) q[0];
rz(-2.2494443) q[1];
sx q[1];
rz(-1.4798857) q[1];
sx q[1];
rz(-2.7977112) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2457857) q[0];
sx q[0];
rz(-0.26660472) q[0];
sx q[0];
rz(-0.33651276) q[0];
rz(-pi) q[1];
rz(0.8240118) q[2];
sx q[2];
rz(-2.4617947) q[2];
sx q[2];
rz(-0.714091) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.77293198) q[1];
sx q[1];
rz(-1.6999869) q[1];
sx q[1];
rz(-0.15993273) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0214693) q[3];
sx q[3];
rz(-1.7645482) q[3];
sx q[3];
rz(-1.4258949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.5188562) q[2];
sx q[2];
rz(-1.7665665) q[2];
sx q[2];
rz(0.95109099) q[2];
rz(2.2226492) q[3];
sx q[3];
rz(-0.94075847) q[3];
sx q[3];
rz(2.9520891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49239355) q[0];
sx q[0];
rz(-2.352582) q[0];
sx q[0];
rz(1.5845818) q[0];
rz(0.7400662) q[1];
sx q[1];
rz(-2.2941755) q[1];
sx q[1];
rz(1.7795631) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.073512162) q[0];
sx q[0];
rz(-2.425504) q[0];
sx q[0];
rz(-2.3329122) q[0];
rz(1.884489) q[2];
sx q[2];
rz(-2.7195632) q[2];
sx q[2];
rz(1.817746) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3950065) q[1];
sx q[1];
rz(-0.94330245) q[1];
sx q[1];
rz(2.021241) q[1];
rz(2.2828987) q[3];
sx q[3];
rz(-1.5017548) q[3];
sx q[3];
rz(-0.25180975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.43807277) q[2];
sx q[2];
rz(-1.1601996) q[2];
sx q[2];
rz(2.4959219) q[2];
rz(3.1072295) q[3];
sx q[3];
rz(-2.2785701) q[3];
sx q[3];
rz(0.063095108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7647758) q[0];
sx q[0];
rz(-2.9314628) q[0];
sx q[0];
rz(0.78395504) q[0];
rz(0.3166554) q[1];
sx q[1];
rz(-1.6461992) q[1];
sx q[1];
rz(2.1276316) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4120042) q[0];
sx q[0];
rz(-2.9656898) q[0];
sx q[0];
rz(0.094502016) q[0];
rz(2.9152206) q[2];
sx q[2];
rz(-1.1034411) q[2];
sx q[2];
rz(-2.4025092) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8617312) q[1];
sx q[1];
rz(-0.35898924) q[1];
sx q[1];
rz(2.3718114) q[1];
rz(-pi) q[2];
rz(2.2714642) q[3];
sx q[3];
rz(-2.4149272) q[3];
sx q[3];
rz(0.18584968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.97518146) q[2];
sx q[2];
rz(-1.3813436) q[2];
sx q[2];
rz(-0.4057917) q[2];
rz(0.25519145) q[3];
sx q[3];
rz(-1.2602592) q[3];
sx q[3];
rz(-2.3658128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.880421) q[0];
sx q[0];
rz(-0.71735993) q[0];
sx q[0];
rz(0.98947155) q[0];
rz(-0.62670341) q[1];
sx q[1];
rz(-0.53755212) q[1];
sx q[1];
rz(-2.7129043) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8256879) q[0];
sx q[0];
rz(-1.6839714) q[0];
sx q[0];
rz(-0.76430385) q[0];
rz(-pi) q[1];
rz(1.7103042) q[2];
sx q[2];
rz(-1.9006471) q[2];
sx q[2];
rz(-0.2784066) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.76213065) q[1];
sx q[1];
rz(-2.3281347) q[1];
sx q[1];
rz(0.7239119) q[1];
x q[2];
rz(-1.8460817) q[3];
sx q[3];
rz(-0.87486736) q[3];
sx q[3];
rz(2.8462209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.89877146) q[2];
sx q[2];
rz(-3.0148744) q[2];
sx q[2];
rz(0.80647331) q[2];
rz(1.9857231) q[3];
sx q[3];
rz(-2.3838796) q[3];
sx q[3];
rz(2.9466467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4194346) q[0];
sx q[0];
rz(-0.86287713) q[0];
sx q[0];
rz(-2.8353598) q[0];
rz(-2.560794) q[1];
sx q[1];
rz(-0.88169801) q[1];
sx q[1];
rz(2.4961848) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5991678) q[0];
sx q[0];
rz(-1.4087447) q[0];
sx q[0];
rz(-0.57181825) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.077351582) q[2];
sx q[2];
rz(-2.0708551) q[2];
sx q[2];
rz(-1.3429221) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4821808) q[1];
sx q[1];
rz(-2.0320091) q[1];
sx q[1];
rz(0.97395514) q[1];
rz(1.4939088) q[3];
sx q[3];
rz(-2.9103177) q[3];
sx q[3];
rz(0.62406711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.2885651) q[2];
sx q[2];
rz(-1.7662851) q[2];
sx q[2];
rz(2.1993401) q[2];
rz(-0.45091584) q[3];
sx q[3];
rz(-2.1555566) q[3];
sx q[3];
rz(-2.8773785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0577724) q[0];
sx q[0];
rz(-1.271893) q[0];
sx q[0];
rz(0.94495946) q[0];
rz(-2.4522929) q[1];
sx q[1];
rz(-2.0336475) q[1];
sx q[1];
rz(1.1030997) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5985918) q[0];
sx q[0];
rz(-0.75816407) q[0];
sx q[0];
rz(-0.36587327) q[0];
rz(-1.9486261) q[2];
sx q[2];
rz(-2.6891064) q[2];
sx q[2];
rz(0.19271344) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0660432) q[1];
sx q[1];
rz(-2.4166921) q[1];
sx q[1];
rz(-2.3182475) q[1];
rz(0.50838105) q[3];
sx q[3];
rz(-1.0219921) q[3];
sx q[3];
rz(2.3387952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5276706) q[2];
sx q[2];
rz(-0.066378243) q[2];
sx q[2];
rz(-1.684368) q[2];
rz(2.1704604) q[3];
sx q[3];
rz(-1.6077653) q[3];
sx q[3];
rz(-0.11626135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74294746) q[0];
sx q[0];
rz(-1.255641) q[0];
sx q[0];
rz(0.1980814) q[0];
rz(2.6598341) q[1];
sx q[1];
rz(-1.2629197) q[1];
sx q[1];
rz(-1.6162965) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2577615) q[0];
sx q[0];
rz(-3.1320268) q[0];
sx q[0];
rz(2.7394858) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5781949) q[2];
sx q[2];
rz(-0.45646898) q[2];
sx q[2];
rz(-2.9204863) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.35514606) q[1];
sx q[1];
rz(-2.0679367) q[1];
sx q[1];
rz(1.0721779) q[1];
rz(-pi) q[2];
rz(1.7500739) q[3];
sx q[3];
rz(-2.6411169) q[3];
sx q[3];
rz(-1.2795606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.6306448) q[2];
sx q[2];
rz(-1.3628565) q[2];
sx q[2];
rz(0.93713588) q[2];
rz(-0.42738327) q[3];
sx q[3];
rz(-1.00777) q[3];
sx q[3];
rz(1.5905323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9464924) q[0];
sx q[0];
rz(-1.4668523) q[0];
sx q[0];
rz(1.3910158) q[0];
rz(-1.4637949) q[1];
sx q[1];
rz(-2.548806) q[1];
sx q[1];
rz(-0.52267271) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7792845) q[0];
sx q[0];
rz(-1.2433508) q[0];
sx q[0];
rz(0.069652005) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2874882) q[2];
sx q[2];
rz(-2.6649545) q[2];
sx q[2];
rz(-2.625537) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9190977) q[1];
sx q[1];
rz(-1.7593707) q[1];
sx q[1];
rz(-1.7401845) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0460879) q[3];
sx q[3];
rz(-1.7874092) q[3];
sx q[3];
rz(-0.10817402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3576144) q[2];
sx q[2];
rz(-1.7366333) q[2];
sx q[2];
rz(-3.0493375) q[2];
rz(-2.9914894) q[3];
sx q[3];
rz(-1.152252) q[3];
sx q[3];
rz(-0.86103719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0263696) q[0];
sx q[0];
rz(-2.8404591) q[0];
sx q[0];
rz(1.2979771) q[0];
rz(-2.2932529) q[1];
sx q[1];
rz(-1.4804761) q[1];
sx q[1];
rz(-0.27378219) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9605676) q[0];
sx q[0];
rz(-1.3342764) q[0];
sx q[0];
rz(-1.5351415) q[0];
x q[1];
rz(2.8086814) q[2];
sx q[2];
rz(-0.71117461) q[2];
sx q[2];
rz(0.43363562) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6322475) q[1];
sx q[1];
rz(-1.6654105) q[1];
sx q[1];
rz(-2.3791887) q[1];
rz(-pi) q[2];
x q[2];
rz(0.43183434) q[3];
sx q[3];
rz(-0.46647989) q[3];
sx q[3];
rz(-1.3386809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0691284) q[2];
sx q[2];
rz(-2.0764565) q[2];
sx q[2];
rz(-0.96159846) q[2];
rz(1.0570863) q[3];
sx q[3];
rz(-2.0545394) q[3];
sx q[3];
rz(0.45007733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3707651) q[0];
sx q[0];
rz(-1.2130883) q[0];
sx q[0];
rz(-0.83592498) q[0];
rz(-1.7302552) q[1];
sx q[1];
rz(-2.5524499) q[1];
sx q[1];
rz(0.020847598) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.037558) q[0];
sx q[0];
rz(-1.7171894) q[0];
sx q[0];
rz(0.030265042) q[0];
rz(-pi) q[1];
rz(-2.8785275) q[2];
sx q[2];
rz(-2.6537958) q[2];
sx q[2];
rz(-0.2476902) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4642568) q[1];
sx q[1];
rz(-1.8337103) q[1];
sx q[1];
rz(-2.8347375) q[1];
x q[2];
rz(-1.3278058) q[3];
sx q[3];
rz(-1.3919481) q[3];
sx q[3];
rz(0.91142765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.47936049) q[2];
sx q[2];
rz(-1.6785097) q[2];
sx q[2];
rz(1.4987017) q[2];
rz(-2.6920476) q[3];
sx q[3];
rz(-2.7142664) q[3];
sx q[3];
rz(0.27165616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0679006) q[0];
sx q[0];
rz(-1.3561159) q[0];
sx q[0];
rz(2.7664716) q[0];
rz(0.66315229) q[1];
sx q[1];
rz(-1.8795492) q[1];
sx q[1];
rz(2.1660027) q[1];
rz(-3.01928) q[2];
sx q[2];
rz(-1.3795102) q[2];
sx q[2];
rz(-1.1176197) q[2];
rz(-1.258044) q[3];
sx q[3];
rz(-2.0087852) q[3];
sx q[3];
rz(1.164418) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
