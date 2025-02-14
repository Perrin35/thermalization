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
rz(2.0410886) q[0];
sx q[0];
rz(-0.70900822) q[0];
sx q[0];
rz(0.79431835) q[0];
rz(0.89214832) q[1];
sx q[1];
rz(-1.661707) q[1];
sx q[1];
rz(2.7977112) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5479255) q[0];
sx q[0];
rz(-1.3194783) q[0];
sx q[0];
rz(-1.4808606) q[0];
rz(2.1062102) q[2];
sx q[2];
rz(-1.1296002) q[2];
sx q[2];
rz(0.23278415) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.12394373) q[1];
sx q[1];
rz(-0.20524612) q[1];
sx q[1];
rz(2.4573047) q[1];
rz(-pi) q[2];
rz(-1.7659193) q[3];
sx q[3];
rz(-1.6886615) q[3];
sx q[3];
rz(3.0199277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6227365) q[2];
sx q[2];
rz(-1.7665665) q[2];
sx q[2];
rz(-2.1905017) q[2];
rz(-2.2226492) q[3];
sx q[3];
rz(-2.2008342) q[3];
sx q[3];
rz(-0.18950352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49239355) q[0];
sx q[0];
rz(-2.352582) q[0];
sx q[0];
rz(-1.5570109) q[0];
rz(-2.4015265) q[1];
sx q[1];
rz(-2.2941755) q[1];
sx q[1];
rz(1.7795631) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3131378) q[0];
sx q[0];
rz(-1.0759996) q[0];
sx q[0];
rz(2.6005756) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1671503) q[2];
sx q[2];
rz(-1.4440618) q[2];
sx q[2];
rz(-0.53467804) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7465861) q[1];
sx q[1];
rz(-2.1982902) q[1];
sx q[1];
rz(1.1203517) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2828987) q[3];
sx q[3];
rz(-1.5017548) q[3];
sx q[3];
rz(2.8897829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7035199) q[2];
sx q[2];
rz(-1.9813931) q[2];
sx q[2];
rz(0.6456708) q[2];
rz(-0.034363184) q[3];
sx q[3];
rz(-2.2785701) q[3];
sx q[3];
rz(0.063095108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
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
rz(-2.8249373) q[1];
sx q[1];
rz(-1.6461992) q[1];
sx q[1];
rz(-1.0139611) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82556258) q[0];
sx q[0];
rz(-1.7459062) q[0];
sx q[0];
rz(1.5540261) q[0];
rz(-2.9152206) q[2];
sx q[2];
rz(-2.0381515) q[2];
sx q[2];
rz(0.73908347) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.55401504) q[1];
sx q[1];
rz(-1.8178175) q[1];
sx q[1];
rz(-2.8783911) q[1];
rz(-pi) q[2];
rz(-0.87012847) q[3];
sx q[3];
rz(-0.72666541) q[3];
sx q[3];
rz(2.955743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.97518146) q[2];
sx q[2];
rz(-1.3813436) q[2];
sx q[2];
rz(0.4057917) q[2];
rz(-0.25519145) q[3];
sx q[3];
rz(-1.8813335) q[3];
sx q[3];
rz(0.7757799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.880421) q[0];
sx q[0];
rz(-2.4242327) q[0];
sx q[0];
rz(-2.1521211) q[0];
rz(-0.62670341) q[1];
sx q[1];
rz(-2.6040405) q[1];
sx q[1];
rz(2.7129043) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1470448) q[0];
sx q[0];
rz(-0.8126077) q[0];
sx q[0];
rz(-1.414624) q[0];
rz(-pi) q[1];
rz(0.33285411) q[2];
sx q[2];
rz(-1.4388545) q[2];
sx q[2];
rz(-1.8037519) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.379462) q[1];
sx q[1];
rz(-2.3281347) q[1];
sx q[1];
rz(2.4176808) q[1];
rz(-pi) q[2];
rz(1.8460817) q[3];
sx q[3];
rz(-0.87486736) q[3];
sx q[3];
rz(0.29537173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.89877146) q[2];
sx q[2];
rz(-0.12671825) q[2];
sx q[2];
rz(2.3351193) q[2];
rz(-1.9857231) q[3];
sx q[3];
rz(-2.3838796) q[3];
sx q[3];
rz(-2.9466467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(2.4194346) q[0];
sx q[0];
rz(-0.86287713) q[0];
sx q[0];
rz(0.3062329) q[0];
rz(-2.560794) q[1];
sx q[1];
rz(-2.2598946) q[1];
sx q[1];
rz(0.64540783) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54242486) q[0];
sx q[0];
rz(-1.732848) q[0];
sx q[0];
rz(2.5697744) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7112947) q[2];
sx q[2];
rz(-2.6360883) q[2];
sx q[2];
rz(1.5031832) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6594118) q[1];
sx q[1];
rz(-1.1095835) q[1];
sx q[1];
rz(-0.97395514) q[1];
rz(0.018086334) q[3];
sx q[3];
rz(-1.340217) q[3];
sx q[3];
rz(-0.54508506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.2885651) q[2];
sx q[2];
rz(-1.7662851) q[2];
sx q[2];
rz(-0.94225252) q[2];
rz(-2.6906768) q[3];
sx q[3];
rz(-0.98603606) q[3];
sx q[3];
rz(0.26421419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0577724) q[0];
sx q[0];
rz(-1.271893) q[0];
sx q[0];
rz(-0.94495946) q[0];
rz(2.4522929) q[1];
sx q[1];
rz(-2.0336475) q[1];
sx q[1];
rz(2.0384929) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5430008) q[0];
sx q[0];
rz(-0.75816407) q[0];
sx q[0];
rz(2.7757194) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9641453) q[2];
sx q[2];
rz(-1.1523231) q[2];
sx q[2];
rz(0.22289101) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0755495) q[1];
sx q[1];
rz(-2.4166921) q[1];
sx q[1];
rz(-0.82334519) q[1];
rz(-pi) q[2];
rz(0.89845851) q[3];
sx q[3];
rz(-0.72988011) q[3];
sx q[3];
rz(1.5204423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5276706) q[2];
sx q[2];
rz(-3.0752144) q[2];
sx q[2];
rz(-1.684368) q[2];
rz(2.1704604) q[3];
sx q[3];
rz(-1.5338273) q[3];
sx q[3];
rz(-3.0253313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
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
rz(-2.9435112) q[0];
rz(0.48175851) q[1];
sx q[1];
rz(-1.878673) q[1];
sx q[1];
rz(-1.6162965) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2859545) q[0];
sx q[0];
rz(-1.5619935) q[0];
sx q[0];
rz(-1.5745401) q[0];
rz(-pi) q[1];
rz(1.1143382) q[2];
sx q[2];
rz(-1.5740574) q[2];
sx q[2];
rz(-1.3563311) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1800148) q[1];
sx q[1];
rz(-1.1369929) q[1];
sx q[1];
rz(-2.5881744) q[1];
x q[2];
rz(-2.0644998) q[3];
sx q[3];
rz(-1.6564661) q[3];
sx q[3];
rz(-0.44891294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.6306448) q[2];
sx q[2];
rz(-1.7787361) q[2];
sx q[2];
rz(-0.93713588) q[2];
rz(0.42738327) q[3];
sx q[3];
rz(-2.1338227) q[3];
sx q[3];
rz(1.5905323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9464924) q[0];
sx q[0];
rz(-1.6747403) q[0];
sx q[0];
rz(-1.7505769) q[0];
rz(1.6777978) q[1];
sx q[1];
rz(-0.59278667) q[1];
sx q[1];
rz(0.52267271) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2309225) q[0];
sx q[0];
rz(-1.5048507) q[0];
sx q[0];
rz(1.8989819) q[0];
rz(-pi) q[1];
rz(-1.2874882) q[2];
sx q[2];
rz(-2.6649545) q[2];
sx q[2];
rz(0.51605564) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1795038) q[1];
sx q[1];
rz(-0.25280646) q[1];
sx q[1];
rz(-0.72353717) q[1];
rz(-pi) q[2];
rz(-0.24899616) q[3];
sx q[3];
rz(-1.0595624) q[3];
sx q[3];
rz(1.5552023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.78397822) q[2];
sx q[2];
rz(-1.7366333) q[2];
sx q[2];
rz(0.092255175) q[2];
rz(-0.15010321) q[3];
sx q[3];
rz(-1.152252) q[3];
sx q[3];
rz(-2.2805555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0263696) q[0];
sx q[0];
rz(-0.30113354) q[0];
sx q[0];
rz(1.8436155) q[0];
rz(0.84833974) q[1];
sx q[1];
rz(-1.4804761) q[1];
sx q[1];
rz(2.8678105) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39812931) q[0];
sx q[0];
rz(-1.6054581) q[0];
sx q[0];
rz(-0.23666478) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8086814) q[2];
sx q[2];
rz(-0.71117461) q[2];
sx q[2];
rz(0.43363562) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1514329) q[1];
sx q[1];
rz(-2.3289375) q[1];
sx q[1];
rz(-1.7012738) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7785092) q[3];
sx q[3];
rz(-1.150032) q[3];
sx q[3];
rz(-1.3265691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0691284) q[2];
sx q[2];
rz(-2.0764565) q[2];
sx q[2];
rz(-2.1799942) q[2];
rz(-1.0570863) q[3];
sx q[3];
rz(-1.0870533) q[3];
sx q[3];
rz(-2.6915153) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7708275) q[0];
sx q[0];
rz(-1.2130883) q[0];
sx q[0];
rz(-2.3056677) q[0];
rz(-1.4113374) q[1];
sx q[1];
rz(-0.58914271) q[1];
sx q[1];
rz(-3.1207451) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3086714) q[0];
sx q[0];
rz(-0.14946689) q[0];
sx q[0];
rz(1.368379) q[0];
x q[1];
rz(1.4336939) q[2];
sx q[2];
rz(-1.1011656) q[2];
sx q[2];
rz(-0.048222311) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1659662) q[1];
sx q[1];
rz(-1.2748112) q[1];
sx q[1];
rz(1.8459666) q[1];
x q[2];
rz(0.18414149) q[3];
sx q[3];
rz(-1.3317579) q[3];
sx q[3];
rz(-2.4381541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6622322) q[2];
sx q[2];
rz(-1.4630829) q[2];
sx q[2];
rz(-1.4987017) q[2];
rz(-0.44954506) q[3];
sx q[3];
rz(-2.7142664) q[3];
sx q[3];
rz(-0.27165616) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0736921) q[0];
sx q[0];
rz(-1.7854767) q[0];
sx q[0];
rz(-0.37512107) q[0];
rz(2.4784404) q[1];
sx q[1];
rz(-1.2620435) q[1];
sx q[1];
rz(-0.97558998) q[1];
rz(-2.1330053) q[2];
sx q[2];
rz(-0.22664438) q[2];
sx q[2];
rz(1.4500153) q[2];
rz(0.58126311) q[3];
sx q[3];
rz(-0.5322335) q[3];
sx q[3];
rz(1.8158326) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
