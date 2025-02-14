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
rz(4.1320463) q[1];
sx q[1];
rz(10.208701) q[1];
sx q[2];
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
rz(-1.4874081) q[2];
sx q[2];
rz(-1.4321799) q[2];
sx q[2];
rz(2.1237789) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8417518) q[1];
sx q[1];
rz(-1.0070485) q[1];
sx q[1];
rz(2.4251516) q[1];
rz(-pi) q[2];
x q[2];
rz(0.089139912) q[3];
sx q[3];
rz(-2.3618792) q[3];
sx q[3];
rz(0.69621315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.94773942) q[2];
sx q[2];
rz(-1.0449907) q[2];
sx q[2];
rz(0.6673153) q[2];
rz(2.8269178) q[3];
sx q[3];
rz(-2.5421725) q[3];
sx q[3];
rz(0.92710322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
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
rz(-0.71669465) q[1];
sx q[1];
rz(2.2915548) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93545234) q[0];
sx q[0];
rz(-0.61203921) q[0];
sx q[0];
rz(-1.7382337) q[0];
rz(2.9710567) q[2];
sx q[2];
rz(-1.4922304) q[2];
sx q[2];
rz(0.54135347) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1239224) q[1];
sx q[1];
rz(-2.7233989) q[1];
sx q[1];
rz(-2.1060418) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.032046958) q[3];
sx q[3];
rz(-0.40478313) q[3];
sx q[3];
rz(-1.3758152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5726418) q[2];
sx q[2];
rz(-0.53223842) q[2];
sx q[2];
rz(2.4786095) q[2];
rz(0.7524544) q[3];
sx q[3];
rz(-1.3198676) q[3];
sx q[3];
rz(-0.19645709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3982584) q[0];
sx q[0];
rz(-0.4011811) q[0];
sx q[0];
rz(0.67984003) q[0];
rz(0.72194779) q[1];
sx q[1];
rz(-0.55689055) q[1];
sx q[1];
rz(0.78071761) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80870285) q[0];
sx q[0];
rz(-1.5631544) q[0];
sx q[0];
rz(-1.9714668) q[0];
rz(-1.6147829) q[2];
sx q[2];
rz(-3.0553748) q[2];
sx q[2];
rz(1.1607064) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.91619766) q[1];
sx q[1];
rz(-2.19529) q[1];
sx q[1];
rz(1.1033113) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9227681) q[3];
sx q[3];
rz(-2.2940562) q[3];
sx q[3];
rz(-1.3506324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9146933) q[2];
sx q[2];
rz(-1.2279899) q[2];
sx q[2];
rz(-2.4719888) q[2];
rz(2.2412444) q[3];
sx q[3];
rz(-1.0122296) q[3];
sx q[3];
rz(2.3646234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8825619) q[0];
sx q[0];
rz(-2.7029523) q[0];
sx q[0];
rz(2.2341109) q[0];
rz(-2.5471845) q[1];
sx q[1];
rz(-2.2120357) q[1];
sx q[1];
rz(1.1297191) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8639899) q[0];
sx q[0];
rz(-2.4765795) q[0];
sx q[0];
rz(-1.8825085) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8215034) q[2];
sx q[2];
rz(-2.5959229) q[2];
sx q[2];
rz(-0.16599338) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.65986827) q[1];
sx q[1];
rz(-2.4818083) q[1];
sx q[1];
rz(2.6901812) q[1];
rz(1.2775189) q[3];
sx q[3];
rz(-1.3012816) q[3];
sx q[3];
rz(0.5130918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.10924673) q[2];
sx q[2];
rz(-1.3999908) q[2];
sx q[2];
rz(2.4692811) q[2];
rz(0.0191056) q[3];
sx q[3];
rz(-2.8002383) q[3];
sx q[3];
rz(-0.13489558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8371209) q[0];
sx q[0];
rz(-2.3470375) q[0];
sx q[0];
rz(-2.9735612) q[0];
rz(1.9145603) q[1];
sx q[1];
rz(-1.2201759) q[1];
sx q[1];
rz(0.29774818) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79556161) q[0];
sx q[0];
rz(-2.3712284) q[0];
sx q[0];
rz(2.3576231) q[0];
x q[1];
rz(-0.99397583) q[2];
sx q[2];
rz(-2.0703531) q[2];
sx q[2];
rz(-3.1021038) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7682957) q[1];
sx q[1];
rz(-1.776564) q[1];
sx q[1];
rz(-1.5453669) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.00068126531) q[3];
sx q[3];
rz(-1.4120123) q[3];
sx q[3];
rz(-0.53631594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1905404) q[2];
sx q[2];
rz(-2.0420044) q[2];
sx q[2];
rz(-2.5229048) q[2];
rz(-0.80858532) q[3];
sx q[3];
rz(-2.6638668) q[3];
sx q[3];
rz(1.3396858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.488778) q[0];
sx q[0];
rz(-1.5050911) q[0];
sx q[0];
rz(-2.6763647) q[0];
rz(-2.6564927) q[1];
sx q[1];
rz(-2.7310889) q[1];
sx q[1];
rz(-3.0533275) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8387091) q[0];
sx q[0];
rz(-2.0219936) q[0];
sx q[0];
rz(-1.0171577) q[0];
rz(-pi) q[1];
rz(-1.5563929) q[2];
sx q[2];
rz(-2.4126137) q[2];
sx q[2];
rz(-0.62504238) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0572964) q[1];
sx q[1];
rz(-2.3688593) q[1];
sx q[1];
rz(0.23717662) q[1];
rz(-pi) q[2];
rz(-2.8507455) q[3];
sx q[3];
rz(-2.8418782) q[3];
sx q[3];
rz(-2.7951473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7469067) q[2];
sx q[2];
rz(-2.2862819) q[2];
sx q[2];
rz(-0.49760231) q[2];
rz(2.6615214) q[3];
sx q[3];
rz(-2.2205455) q[3];
sx q[3];
rz(0.068643071) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8467872) q[0];
sx q[0];
rz(-2.6755896) q[0];
sx q[0];
rz(1.71126) q[0];
rz(1.3149186) q[1];
sx q[1];
rz(-0.39350915) q[1];
sx q[1];
rz(-1.5550782) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99796999) q[0];
sx q[0];
rz(-0.41150948) q[0];
sx q[0];
rz(-1.8726148) q[0];
rz(-0.58760037) q[2];
sx q[2];
rz(-0.94600224) q[2];
sx q[2];
rz(-1.8163565) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3433229) q[1];
sx q[1];
rz(-2.1499691) q[1];
sx q[1];
rz(3.0208241) q[1];
rz(-pi) q[2];
rz(0.076818941) q[3];
sx q[3];
rz(-2.2008476) q[3];
sx q[3];
rz(-1.0162285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.84307182) q[2];
sx q[2];
rz(-2.4303747) q[2];
sx q[2];
rz(-2.2141875) q[2];
rz(1.1585506) q[3];
sx q[3];
rz(-2.4535593) q[3];
sx q[3];
rz(0.096175171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4226828) q[0];
sx q[0];
rz(-2.2490608) q[0];
sx q[0];
rz(-0.476015) q[0];
rz(-0.99126518) q[1];
sx q[1];
rz(-0.74949336) q[1];
sx q[1];
rz(-0.083866619) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0821486) q[0];
sx q[0];
rz(-2.5184439) q[0];
sx q[0];
rz(-2.2150458) q[0];
rz(-pi) q[1];
rz(-0.4617347) q[2];
sx q[2];
rz(-1.8991915) q[2];
sx q[2];
rz(-0.4730313) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2460275) q[1];
sx q[1];
rz(-1.152134) q[1];
sx q[1];
rz(-3.0484867) q[1];
x q[2];
rz(0.08633498) q[3];
sx q[3];
rz(-1.207296) q[3];
sx q[3];
rz(2.8192161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.91741651) q[2];
sx q[2];
rz(-0.26599628) q[2];
sx q[2];
rz(2.4931397) q[2];
rz(-0.90483856) q[3];
sx q[3];
rz(-1.6831393) q[3];
sx q[3];
rz(-2.8308433) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34173486) q[0];
sx q[0];
rz(-0.71001995) q[0];
sx q[0];
rz(2.572686) q[0];
rz(-2.3078602) q[1];
sx q[1];
rz(-1.2833475) q[1];
sx q[1];
rz(2.1368829) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6093401) q[0];
sx q[0];
rz(-2.1933317) q[0];
sx q[0];
rz(-1.3329766) q[0];
rz(1.1962074) q[2];
sx q[2];
rz(-0.56972144) q[2];
sx q[2];
rz(-1.2562871) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9819599) q[1];
sx q[1];
rz(-1.4985604) q[1];
sx q[1];
rz(-0.87149974) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6528204) q[3];
sx q[3];
rz(-2.1567705) q[3];
sx q[3];
rz(0.88912933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7052762) q[2];
sx q[2];
rz(-0.87203163) q[2];
sx q[2];
rz(0.3479859) q[2];
rz(2.3157388) q[3];
sx q[3];
rz(-2.9233942) q[3];
sx q[3];
rz(2.5364449) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0843049) q[0];
sx q[0];
rz(-2.075752) q[0];
sx q[0];
rz(-2.9344946) q[0];
rz(2.6026978) q[1];
sx q[1];
rz(-2.5205044) q[1];
sx q[1];
rz(2.5394687) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83610632) q[0];
sx q[0];
rz(-1.0998816) q[0];
sx q[0];
rz(0.54022809) q[0];
rz(-pi) q[1];
rz(0.58790284) q[2];
sx q[2];
rz(-1.9155115) q[2];
sx q[2];
rz(1.6893139) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6844897) q[1];
sx q[1];
rz(-0.41375638) q[1];
sx q[1];
rz(1.7587508) q[1];
rz(-1.8982417) q[3];
sx q[3];
rz(-0.85059887) q[3];
sx q[3];
rz(2.7011288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.78486097) q[2];
sx q[2];
rz(-2.1897903) q[2];
sx q[2];
rz(1.3447364) q[2];
rz(2.6209659) q[3];
sx q[3];
rz(-0.27025637) q[3];
sx q[3];
rz(2.3394913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6945334) q[0];
sx q[0];
rz(-0.72493989) q[0];
sx q[0];
rz(-1.3865393) q[0];
rz(0.78320349) q[1];
sx q[1];
rz(-1.422311) q[1];
sx q[1];
rz(1.5625988) q[1];
rz(-2.4827847) q[2];
sx q[2];
rz(-0.42429069) q[2];
sx q[2];
rz(-1.9784603) q[2];
rz(-2.9014723) q[3];
sx q[3];
rz(-1.5821531) q[3];
sx q[3];
rz(1.5021363) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
