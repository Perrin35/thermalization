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
rz(-2.146848) q[0];
sx q[0];
rz(-0.58726197) q[0];
sx q[0];
rz(-0.62556148) q[0];
rz(0.62043959) q[1];
sx q[1];
rz(-2.151139) q[1];
sx q[1];
rz(0.78392309) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44255689) q[0];
sx q[0];
rz(-0.43414206) q[0];
sx q[0];
rz(3.1331557) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.002499) q[2];
sx q[2];
rz(-1.6533829) q[2];
sx q[2];
rz(-2.6001584) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.82134889) q[1];
sx q[1];
rz(-0.87961266) q[1];
sx q[1];
rz(-0.76637474) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.089139912) q[3];
sx q[3];
rz(-2.3618792) q[3];
sx q[3];
rz(-0.69621315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.94773942) q[2];
sx q[2];
rz(-1.0449907) q[2];
sx q[2];
rz(2.4742773) q[2];
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
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2941403) q[0];
sx q[0];
rz(-0.84646928) q[0];
sx q[0];
rz(-0.2997998) q[0];
rz(-1.0046129) q[1];
sx q[1];
rz(-2.424898) q[1];
sx q[1];
rz(2.2915548) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93545234) q[0];
sx q[0];
rz(-0.61203921) q[0];
sx q[0];
rz(-1.403359) q[0];
rz(-pi) q[1];
rz(0.17053594) q[2];
sx q[2];
rz(-1.4922304) q[2];
sx q[2];
rz(2.6002392) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5483173) q[1];
sx q[1];
rz(-1.2139582) q[1];
sx q[1];
rz(-2.9186881) q[1];
rz(-pi) q[2];
rz(-3.1095457) q[3];
sx q[3];
rz(-0.40478313) q[3];
sx q[3];
rz(-1.7657775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5726418) q[2];
sx q[2];
rz(-0.53223842) q[2];
sx q[2];
rz(-2.4786095) q[2];
rz(2.3891383) q[3];
sx q[3];
rz(-1.821725) q[3];
sx q[3];
rz(-0.19645709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74333423) q[0];
sx q[0];
rz(-2.7404116) q[0];
sx q[0];
rz(-2.4617526) q[0];
rz(-0.72194779) q[1];
sx q[1];
rz(-0.55689055) q[1];
sx q[1];
rz(2.360875) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3975383) q[0];
sx q[0];
rz(-2.7408532) q[0];
sx q[0];
rz(1.5512054) q[0];
x q[1];
rz(-1.6569312) q[2];
sx q[2];
rz(-1.5670098) q[2];
sx q[2];
rz(0.45391339) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6284077) q[1];
sx q[1];
rz(-2.3806913) q[1];
sx q[1];
rz(2.582798) q[1];
rz(-1.3296919) q[3];
sx q[3];
rz(-2.3917197) q[3];
sx q[3];
rz(-2.1151224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9146933) q[2];
sx q[2];
rz(-1.2279899) q[2];
sx q[2];
rz(-0.66960382) q[2];
rz(-2.2412444) q[3];
sx q[3];
rz(-2.1293631) q[3];
sx q[3];
rz(-0.77696925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8825619) q[0];
sx q[0];
rz(-2.7029523) q[0];
sx q[0];
rz(-0.90748179) q[0];
rz(2.5471845) q[1];
sx q[1];
rz(-2.2120357) q[1];
sx q[1];
rz(-1.1297191) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0967207) q[0];
sx q[0];
rz(-1.3803998) q[0];
sx q[0];
rz(-2.2119766) q[0];
rz(-pi) q[1];
x q[1];
rz(0.14950651) q[2];
sx q[2];
rz(-1.0440011) q[2];
sx q[2];
rz(0.12509987) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0320624) q[1];
sx q[1];
rz(-0.9865762) q[1];
sx q[1];
rz(-1.24448) q[1];
rz(-pi) q[2];
rz(-1.2775189) q[3];
sx q[3];
rz(-1.8403111) q[3];
sx q[3];
rz(-2.6285009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.10924673) q[2];
sx q[2];
rz(-1.7416019) q[2];
sx q[2];
rz(-2.4692811) q[2];
rz(-3.1224871) q[3];
sx q[3];
rz(-2.8002383) q[3];
sx q[3];
rz(-0.13489558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30447176) q[0];
sx q[0];
rz(-2.3470375) q[0];
sx q[0];
rz(2.9735612) q[0];
rz(-1.2270323) q[1];
sx q[1];
rz(-1.2201759) q[1];
sx q[1];
rz(-2.8438445) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7452324) q[0];
sx q[0];
rz(-1.0567291) q[0];
sx q[0];
rz(2.5395509) q[0];
rz(-pi) q[1];
rz(-2.1476168) q[2];
sx q[2];
rz(-2.0703531) q[2];
sx q[2];
rz(3.1021038) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6444466) q[1];
sx q[1];
rz(-2.9342817) q[1];
sx q[1];
rz(-0.12122341) q[1];
x q[2];
rz(-3.1409114) q[3];
sx q[3];
rz(-1.7295803) q[3];
sx q[3];
rz(2.6052767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1905404) q[2];
sx q[2];
rz(-2.0420044) q[2];
sx q[2];
rz(2.5229048) q[2];
rz(-0.80858532) q[3];
sx q[3];
rz(-0.4777258) q[3];
sx q[3];
rz(1.8019069) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65281463) q[0];
sx q[0];
rz(-1.5050911) q[0];
sx q[0];
rz(-0.46522796) q[0];
rz(0.48509994) q[1];
sx q[1];
rz(-2.7310889) q[1];
sx q[1];
rz(-3.0533275) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6104077) q[0];
sx q[0];
rz(-1.0779128) q[0];
sx q[0];
rz(-0.51778527) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5851998) q[2];
sx q[2];
rz(-0.72897899) q[2];
sx q[2];
rz(0.62504238) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0572964) q[1];
sx q[1];
rz(-2.3688593) q[1];
sx q[1];
rz(0.23717662) q[1];
rz(-pi) q[2];
rz(-0.29084716) q[3];
sx q[3];
rz(-2.8418782) q[3];
sx q[3];
rz(-0.3464454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7469067) q[2];
sx q[2];
rz(-2.2862819) q[2];
sx q[2];
rz(-0.49760231) q[2];
rz(0.48007128) q[3];
sx q[3];
rz(-2.2205455) q[3];
sx q[3];
rz(3.0729496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8467872) q[0];
sx q[0];
rz(-2.6755896) q[0];
sx q[0];
rz(-1.4303327) q[0];
rz(1.826674) q[1];
sx q[1];
rz(-0.39350915) q[1];
sx q[1];
rz(1.5550782) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85077943) q[0];
sx q[0];
rz(-1.4516136) q[0];
sx q[0];
rz(-1.1759773) q[0];
rz(-2.2847963) q[2];
sx q[2];
rz(-1.1044377) q[2];
sx q[2];
rz(-0.12596065) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.83885289) q[1];
sx q[1];
rz(-1.4697971) q[1];
sx q[1];
rz(-2.1533226) q[1];
rz(-2.2022543) q[3];
sx q[3];
rz(-1.6328446) q[3];
sx q[3];
rz(0.59988672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.84307182) q[2];
sx q[2];
rz(-0.71121794) q[2];
sx q[2];
rz(0.92740518) q[2];
rz(-1.983042) q[3];
sx q[3];
rz(-0.68803334) q[3];
sx q[3];
rz(3.0454175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.4226828) q[0];
sx q[0];
rz(-2.2490608) q[0];
sx q[0];
rz(2.6655777) q[0];
rz(0.99126518) q[1];
sx q[1];
rz(-0.74949336) q[1];
sx q[1];
rz(0.083866619) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.335673) q[0];
sx q[0];
rz(-1.085338) q[0];
sx q[0];
rz(0.40747633) q[0];
rz(-pi) q[1];
rz(0.4617347) q[2];
sx q[2];
rz(-1.2424011) q[2];
sx q[2];
rz(-0.4730313) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1213346) q[1];
sx q[1];
rz(-2.7133006) q[1];
sx q[1];
rz(-1.7767724) q[1];
rz(-pi) q[2];
rz(1.2060542) q[3];
sx q[3];
rz(-1.6514773) q[3];
sx q[3];
rz(1.9239359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2241761) q[2];
sx q[2];
rz(-2.8755964) q[2];
sx q[2];
rz(0.64845294) q[2];
rz(-2.2367541) q[3];
sx q[3];
rz(-1.4584533) q[3];
sx q[3];
rz(0.31074935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34173486) q[0];
sx q[0];
rz(-0.71001995) q[0];
sx q[0];
rz(-2.572686) q[0];
rz(-2.3078602) q[1];
sx q[1];
rz(-1.8582452) q[1];
sx q[1];
rz(-2.1368829) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2434655) q[0];
sx q[0];
rz(-1.7633738) q[0];
sx q[0];
rz(-0.63611998) q[0];
rz(-2.1083852) q[2];
sx q[2];
rz(-1.3721319) q[2];
sx q[2];
rz(0.63420701) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6697996) q[1];
sx q[1];
rz(-0.87369117) q[1];
sx q[1];
rz(-3.0473188) q[1];
rz(-pi) q[2];
rz(-1.6528204) q[3];
sx q[3];
rz(-0.98482212) q[3];
sx q[3];
rz(-2.2524633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7052762) q[2];
sx q[2];
rz(-0.87203163) q[2];
sx q[2];
rz(-0.3479859) q[2];
rz(2.3157388) q[3];
sx q[3];
rz(-0.21819849) q[3];
sx q[3];
rz(0.60514778) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0572877) q[0];
sx q[0];
rz(-2.075752) q[0];
sx q[0];
rz(2.9344946) q[0];
rz(-0.53889489) q[1];
sx q[1];
rz(-0.62108827) q[1];
sx q[1];
rz(0.60212392) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83610632) q[0];
sx q[0];
rz(-2.0417111) q[0];
sx q[0];
rz(2.6013646) q[0];
rz(-pi) q[1];
rz(0.57453491) q[2];
sx q[2];
rz(-2.4705187) q[2];
sx q[2];
rz(2.7908763) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.94127266) q[1];
sx q[1];
rz(-1.6459904) q[1];
sx q[1];
rz(-1.1635416) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.243351) q[3];
sx q[3];
rz(-0.85059887) q[3];
sx q[3];
rz(-2.7011288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3567317) q[2];
sx q[2];
rz(-0.95180231) q[2];
sx q[2];
rz(1.7968563) q[2];
rz(-2.6209659) q[3];
sx q[3];
rz(-0.27025637) q[3];
sx q[3];
rz(0.80210137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44705924) q[0];
sx q[0];
rz(-0.72493989) q[0];
sx q[0];
rz(-1.3865393) q[0];
rz(2.3583892) q[1];
sx q[1];
rz(-1.7192817) q[1];
sx q[1];
rz(-1.5789938) q[1];
rz(0.34306768) q[2];
sx q[2];
rz(-1.3160327) q[2];
sx q[2];
rz(2.119488) q[2];
rz(-0.24012031) q[3];
sx q[3];
rz(-1.5594395) q[3];
sx q[3];
rz(-1.6394564) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
