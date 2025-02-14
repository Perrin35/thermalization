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
rz(0.99474466) q[0];
sx q[0];
rz(-2.5543307) q[0];
sx q[0];
rz(0.62556148) q[0];
rz(-2.5211531) q[1];
sx q[1];
rz(-0.99045366) q[1];
sx q[1];
rz(2.3576696) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6990358) q[0];
sx q[0];
rz(-2.7074506) q[0];
sx q[0];
rz(0.0084369466) q[0];
rz(-pi) q[1];
rz(0.13909362) q[2];
sx q[2];
rz(-1.4882097) q[2];
sx q[2];
rz(2.6001584) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3202438) q[1];
sx q[1];
rz(-2.26198) q[1];
sx q[1];
rz(2.3752179) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.089139912) q[3];
sx q[3];
rz(-2.3618792) q[3];
sx q[3];
rz(2.4453795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1938532) q[2];
sx q[2];
rz(-2.096602) q[2];
sx q[2];
rz(-0.6673153) q[2];
rz(-2.8269178) q[3];
sx q[3];
rz(-0.59942013) q[3];
sx q[3];
rz(-2.2144894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2941403) q[0];
sx q[0];
rz(-2.2951234) q[0];
sx q[0];
rz(-2.8417929) q[0];
rz(-2.1369797) q[1];
sx q[1];
rz(-2.424898) q[1];
sx q[1];
rz(-2.2915548) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3687821) q[0];
sx q[0];
rz(-1.4748992) q[0];
sx q[0];
rz(-0.96536388) q[0];
x q[1];
rz(0.17053594) q[2];
sx q[2];
rz(-1.6493622) q[2];
sx q[2];
rz(0.54135347) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1239224) q[1];
sx q[1];
rz(-2.7233989) q[1];
sx q[1];
rz(-2.1060418) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5570693) q[3];
sx q[3];
rz(-1.9753595) q[3];
sx q[3];
rz(-1.8006397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5726418) q[2];
sx q[2];
rz(-0.53223842) q[2];
sx q[2];
rz(0.66298318) q[2];
rz(2.3891383) q[3];
sx q[3];
rz(-1.821725) q[3];
sx q[3];
rz(-0.19645709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74333423) q[0];
sx q[0];
rz(-0.4011811) q[0];
sx q[0];
rz(0.67984003) q[0];
rz(2.4196449) q[1];
sx q[1];
rz(-0.55689055) q[1];
sx q[1];
rz(-0.78071761) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3975383) q[0];
sx q[0];
rz(-2.7408532) q[0];
sx q[0];
rz(1.5512054) q[0];
rz(-pi) q[1];
rz(-1.5268097) q[2];
sx q[2];
rz(-3.0553748) q[2];
sx q[2];
rz(-1.1607064) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.91619766) q[1];
sx q[1];
rz(-2.19529) q[1];
sx q[1];
rz(-2.0382814) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9227681) q[3];
sx q[3];
rz(-0.84753643) q[3];
sx q[3];
rz(1.7909602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.22689936) q[2];
sx q[2];
rz(-1.9136027) q[2];
sx q[2];
rz(0.66960382) q[2];
rz(2.2412444) q[3];
sx q[3];
rz(-2.1293631) q[3];
sx q[3];
rz(0.77696925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(2.8825619) q[0];
sx q[0];
rz(-0.43864033) q[0];
sx q[0];
rz(2.2341109) q[0];
rz(-0.59440815) q[1];
sx q[1];
rz(-2.2120357) q[1];
sx q[1];
rz(2.0118735) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8639899) q[0];
sx q[0];
rz(-2.4765795) q[0];
sx q[0];
rz(-1.2590842) q[0];
x q[1];
rz(0.14950651) q[2];
sx q[2];
rz(-1.0440011) q[2];
sx q[2];
rz(0.12509987) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4817244) q[1];
sx q[1];
rz(-2.4818083) q[1];
sx q[1];
rz(-2.6901812) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.80813415) q[3];
sx q[3];
rz(-2.7459347) q[3];
sx q[3];
rz(1.3612703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.10924673) q[2];
sx q[2];
rz(-1.3999908) q[2];
sx q[2];
rz(0.67231154) q[2];
rz(0.0191056) q[3];
sx q[3];
rz(-0.34135434) q[3];
sx q[3];
rz(-3.0066971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30447176) q[0];
sx q[0];
rz(-0.79455513) q[0];
sx q[0];
rz(2.9735612) q[0];
rz(-1.2270323) q[1];
sx q[1];
rz(-1.9214168) q[1];
sx q[1];
rz(-0.29774818) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3963602) q[0];
sx q[0];
rz(-2.0848635) q[0];
sx q[0];
rz(2.5395509) q[0];
rz(2.3565294) q[2];
sx q[2];
rz(-0.74399555) q[2];
sx q[2];
rz(-0.97563484) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6444466) q[1];
sx q[1];
rz(-2.9342817) q[1];
sx q[1];
rz(0.12122341) q[1];
x q[2];
rz(1.7295804) q[3];
sx q[3];
rz(-1.5701236) q[3];
sx q[3];
rz(-2.1070045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1905404) q[2];
sx q[2];
rz(-2.0420044) q[2];
sx q[2];
rz(2.5229048) q[2];
rz(-0.80858532) q[3];
sx q[3];
rz(-2.6638668) q[3];
sx q[3];
rz(1.3396858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65281463) q[0];
sx q[0];
rz(-1.6365016) q[0];
sx q[0];
rz(-2.6763647) q[0];
rz(2.6564927) q[1];
sx q[1];
rz(-0.41050375) q[1];
sx q[1];
rz(0.088265158) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6104077) q[0];
sx q[0];
rz(-1.0779128) q[0];
sx q[0];
rz(0.51778527) q[0];
rz(-pi) q[1];
rz(-1.5563929) q[2];
sx q[2];
rz(-2.4126137) q[2];
sx q[2];
rz(2.5165503) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7317071) q[1];
sx q[1];
rz(-2.3165858) q[1];
sx q[1];
rz(-1.7959926) q[1];
rz(0.28782423) q[3];
sx q[3];
rz(-1.4860286) q[3];
sx q[3];
rz(-2.1957977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.39468592) q[2];
sx q[2];
rz(-0.8553108) q[2];
sx q[2];
rz(2.6439903) q[2];
rz(-2.6615214) q[3];
sx q[3];
rz(-0.92104715) q[3];
sx q[3];
rz(0.068643071) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8467872) q[0];
sx q[0];
rz(-0.46600309) q[0];
sx q[0];
rz(1.4303327) q[0];
rz(1.3149186) q[1];
sx q[1];
rz(-2.7480835) q[1];
sx q[1];
rz(1.5550782) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67051149) q[0];
sx q[0];
rz(-1.1789315) q[0];
sx q[0];
rz(-3.0125822) q[0];
rz(-pi) q[1];
x q[1];
rz(0.58760037) q[2];
sx q[2];
rz(-0.94600224) q[2];
sx q[2];
rz(1.8163565) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.79826971) q[1];
sx q[1];
rz(-2.1499691) q[1];
sx q[1];
rz(-3.0208241) q[1];
rz(-pi) q[2];
rz(-1.4659381) q[3];
sx q[3];
rz(-0.63408454) q[3];
sx q[3];
rz(2.255267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2985208) q[2];
sx q[2];
rz(-0.71121794) q[2];
sx q[2];
rz(0.92740518) q[2];
rz(-1.983042) q[3];
sx q[3];
rz(-0.68803334) q[3];
sx q[3];
rz(-0.096175171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4226828) q[0];
sx q[0];
rz(-0.89253187) q[0];
sx q[0];
rz(0.476015) q[0];
rz(2.1503275) q[1];
sx q[1];
rz(-2.3920993) q[1];
sx q[1];
rz(-3.057726) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.335673) q[0];
sx q[0];
rz(-1.085338) q[0];
sx q[0];
rz(-0.40747633) q[0];
x q[1];
rz(-2.4886683) q[2];
sx q[2];
rz(-2.5819375) q[2];
sx q[2];
rz(-1.6728354) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2460275) q[1];
sx q[1];
rz(-1.9894587) q[1];
sx q[1];
rz(0.093105969) q[1];
x q[2];
rz(0.08633498) q[3];
sx q[3];
rz(-1.9342967) q[3];
sx q[3];
rz(-2.8192161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.91741651) q[2];
sx q[2];
rz(-2.8755964) q[2];
sx q[2];
rz(-2.4931397) q[2];
rz(-2.2367541) q[3];
sx q[3];
rz(-1.4584533) q[3];
sx q[3];
rz(0.31074935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34173486) q[0];
sx q[0];
rz(-2.4315727) q[0];
sx q[0];
rz(-2.572686) q[0];
rz(-2.3078602) q[1];
sx q[1];
rz(-1.8582452) q[1];
sx q[1];
rz(-2.1368829) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2434655) q[0];
sx q[0];
rz(-1.3782189) q[0];
sx q[0];
rz(0.63611998) q[0];
rz(-pi) q[1];
rz(-2.9113681) q[2];
sx q[2];
rz(-2.0966999) q[2];
sx q[2];
rz(-2.3221225) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8160287) q[1];
sx q[1];
rz(-0.70239151) q[1];
sx q[1];
rz(1.6827453) q[1];
rz(-pi) q[2];
rz(-0.58752693) q[3];
sx q[3];
rz(-1.6391132) q[3];
sx q[3];
rz(-0.63623601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7052762) q[2];
sx q[2];
rz(-2.269561) q[2];
sx q[2];
rz(-0.3479859) q[2];
rz(2.3157388) q[3];
sx q[3];
rz(-2.9233942) q[3];
sx q[3];
rz(2.5364449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0572877) q[0];
sx q[0];
rz(-1.0658406) q[0];
sx q[0];
rz(0.2070981) q[0];
rz(-0.53889489) q[1];
sx q[1];
rz(-0.62108827) q[1];
sx q[1];
rz(0.60212392) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46901921) q[0];
sx q[0];
rz(-1.0946602) q[0];
sx q[0];
rz(-1.0350521) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5670577) q[2];
sx q[2];
rz(-0.67107397) q[2];
sx q[2];
rz(0.35071638) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.20032) q[1];
sx q[1];
rz(-1.6459904) q[1];
sx q[1];
rz(-1.1635416) q[1];
rz(-pi) q[2];
rz(1.8982417) q[3];
sx q[3];
rz(-2.2909938) q[3];
sx q[3];
rz(2.7011288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.78486097) q[2];
sx q[2];
rz(-0.95180231) q[2];
sx q[2];
rz(1.3447364) q[2];
rz(-0.52062672) q[3];
sx q[3];
rz(-0.27025637) q[3];
sx q[3];
rz(2.3394913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6945334) q[0];
sx q[0];
rz(-0.72493989) q[0];
sx q[0];
rz(-1.3865393) q[0];
rz(-0.78320349) q[1];
sx q[1];
rz(-1.7192817) q[1];
sx q[1];
rz(-1.5789938) q[1];
rz(0.65880792) q[2];
sx q[2];
rz(-0.42429069) q[2];
sx q[2];
rz(-1.9784603) q[2];
rz(2.9014723) q[3];
sx q[3];
rz(-1.5594395) q[3];
sx q[3];
rz(-1.6394564) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
