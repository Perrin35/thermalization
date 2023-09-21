OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.1922167) q[0];
sx q[0];
rz(-2.0944216) q[0];
sx q[0];
rz(-0.068724364) q[0];
rz(-1.3955431) q[1];
sx q[1];
rz(-1.5323324) q[1];
sx q[1];
rz(-1.9332164) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8954593) q[0];
sx q[0];
rz(-1.6172505) q[0];
sx q[0];
rz(-1.6669271) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0460998) q[2];
sx q[2];
rz(-0.12636939) q[2];
sx q[2];
rz(-0.40590826) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5282643) q[1];
sx q[1];
rz(-1.4600888) q[1];
sx q[1];
rz(-0.61943357) q[1];
x q[2];
rz(-0.52238676) q[3];
sx q[3];
rz(-1.1484227) q[3];
sx q[3];
rz(1.1952323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8157114) q[2];
sx q[2];
rz(-1.4006961) q[2];
sx q[2];
rz(-1.4665843) q[2];
rz(-2.4438434) q[3];
sx q[3];
rz(-2.0402699) q[3];
sx q[3];
rz(-2.3944323) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0242457) q[0];
sx q[0];
rz(-2.0139366) q[0];
sx q[0];
rz(-1.9673989) q[0];
rz(-2.970447) q[1];
sx q[1];
rz(-2.0967963) q[1];
sx q[1];
rz(0.29719621) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3127808) q[0];
sx q[0];
rz(-3.045407) q[0];
sx q[0];
rz(2.0975153) q[0];
x q[1];
rz(1.0969639) q[2];
sx q[2];
rz(-2.2261438) q[2];
sx q[2];
rz(-1.4807793) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9651523) q[1];
sx q[1];
rz(-1.0607751) q[1];
sx q[1];
rz(0.39217197) q[1];
rz(1.9678715) q[3];
sx q[3];
rz(-0.84871549) q[3];
sx q[3];
rz(1.811036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.16195665) q[2];
sx q[2];
rz(-1.9571783) q[2];
sx q[2];
rz(0.61398181) q[2];
rz(-0.87614122) q[3];
sx q[3];
rz(-2.4738779) q[3];
sx q[3];
rz(2.5045625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0456332) q[0];
sx q[0];
rz(-1.8563844) q[0];
sx q[0];
rz(2.8339548) q[0];
rz(0.74854198) q[1];
sx q[1];
rz(-2.8100439) q[1];
sx q[1];
rz(0.83980733) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1727985) q[0];
sx q[0];
rz(-1.286924) q[0];
sx q[0];
rz(-2.7127405) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4762127) q[2];
sx q[2];
rz(-1.2328096) q[2];
sx q[2];
rz(-0.72159492) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7786583) q[1];
sx q[1];
rz(-1.1204473) q[1];
sx q[1];
rz(1.910701) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3358467) q[3];
sx q[3];
rz(-0.64628212) q[3];
sx q[3];
rz(-0.39138734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2741189) q[2];
sx q[2];
rz(-1.3511191) q[2];
sx q[2];
rz(1.8236558) q[2];
rz(1.2157724) q[3];
sx q[3];
rz(-0.35651818) q[3];
sx q[3];
rz(1.5095476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76386219) q[0];
sx q[0];
rz(-1.3803991) q[0];
sx q[0];
rz(0.44556251) q[0];
rz(2.5207649) q[1];
sx q[1];
rz(-1.2460243) q[1];
sx q[1];
rz(0.96558085) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25528204) q[0];
sx q[0];
rz(-0.83034407) q[0];
sx q[0];
rz(1.1925973) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1104286) q[2];
sx q[2];
rz(-0.84261299) q[2];
sx q[2];
rz(-2.503501) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.671215) q[1];
sx q[1];
rz(-0.2556076) q[1];
sx q[1];
rz(1.4443936) q[1];
x q[2];
rz(-0.86430092) q[3];
sx q[3];
rz(-2.9923277) q[3];
sx q[3];
rz(0.60993689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.821637) q[2];
sx q[2];
rz(-0.76239061) q[2];
sx q[2];
rz(0.92932534) q[2];
rz(-2.4980513) q[3];
sx q[3];
rz(-1.0361592) q[3];
sx q[3];
rz(-1.003456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
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
rz(-1.4716361) q[0];
sx q[0];
rz(-1.5820553) q[0];
sx q[0];
rz(1.2840282) q[0];
rz(0.28981003) q[1];
sx q[1];
rz(-0.73957864) q[1];
sx q[1];
rz(1.0481542) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50131932) q[0];
sx q[0];
rz(-0.87448705) q[0];
sx q[0];
rz(-2.2107844) q[0];
rz(-pi) q[1];
rz(-0.94190188) q[2];
sx q[2];
rz(-2.2976029) q[2];
sx q[2];
rz(-2.294159) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.54742868) q[1];
sx q[1];
rz(-0.33572008) q[1];
sx q[1];
rz(2.8572542) q[1];
rz(-pi) q[2];
rz(-2.2506511) q[3];
sx q[3];
rz(-1.1853301) q[3];
sx q[3];
rz(0.95213529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2709048) q[2];
sx q[2];
rz(-2.1257766) q[2];
sx q[2];
rz(0.90083814) q[2];
rz(-2.0488996) q[3];
sx q[3];
rz(-2.1381502) q[3];
sx q[3];
rz(-1.9074915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.1822405) q[0];
sx q[0];
rz(-1.2077967) q[0];
sx q[0];
rz(0.59610468) q[0];
rz(-1.4959363) q[1];
sx q[1];
rz(-2.2438965) q[1];
sx q[1];
rz(-1.8966819) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0743474) q[0];
sx q[0];
rz(-1.0676358) q[0];
sx q[0];
rz(2.7748681) q[0];
x q[1];
rz(-1.8876569) q[2];
sx q[2];
rz(-2.0358026) q[2];
sx q[2];
rz(0.76380619) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.53828) q[1];
sx q[1];
rz(-2.0671751) q[1];
sx q[1];
rz(-0.75116317) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.75002589) q[3];
sx q[3];
rz(-2.8083028) q[3];
sx q[3];
rz(-0.94543524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.231455) q[2];
sx q[2];
rz(-2.4427876) q[2];
sx q[2];
rz(-2.9166252) q[2];
rz(-3.0531626) q[3];
sx q[3];
rz(-1.4513353) q[3];
sx q[3];
rz(2.6627873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.46463075) q[0];
sx q[0];
rz(-1.2372274) q[0];
sx q[0];
rz(2.8616469) q[0];
rz(-1.6784558) q[1];
sx q[1];
rz(-1.8755553) q[1];
sx q[1];
rz(0.25269145) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0137018) q[0];
sx q[0];
rz(-1.863592) q[0];
sx q[0];
rz(-0.049833628) q[0];
rz(-pi) q[1];
rz(-0.018888868) q[2];
sx q[2];
rz(-0.9364555) q[2];
sx q[2];
rz(-2.6332476) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.686324) q[1];
sx q[1];
rz(-0.9372006) q[1];
sx q[1];
rz(-2.2151161) q[1];
x q[2];
rz(2.3971862) q[3];
sx q[3];
rz(-2.230847) q[3];
sx q[3];
rz(0.50124121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8213886) q[2];
sx q[2];
rz(-0.89367047) q[2];
sx q[2];
rz(-1.5318711) q[2];
rz(1.1931217) q[3];
sx q[3];
rz(-2.2333998) q[3];
sx q[3];
rz(1.0866722) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10483345) q[0];
sx q[0];
rz(-1.3830673) q[0];
sx q[0];
rz(-1.6554792) q[0];
rz(0.2688109) q[1];
sx q[1];
rz(-1.1299645) q[1];
sx q[1];
rz(0.2789467) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69142197) q[0];
sx q[0];
rz(-0.84934635) q[0];
sx q[0];
rz(-2.2539317) q[0];
x q[1];
rz(-1.0686764) q[2];
sx q[2];
rz(-1.4197822) q[2];
sx q[2];
rz(0.47889027) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.61356269) q[1];
sx q[1];
rz(-2.2655305) q[1];
sx q[1];
rz(1.0747521) q[1];
rz(-pi) q[2];
rz(0.99854462) q[3];
sx q[3];
rz(-1.4456985) q[3];
sx q[3];
rz(0.78748122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3387317) q[2];
sx q[2];
rz(-1.677745) q[2];
sx q[2];
rz(1.5489244) q[2];
rz(-2.0907949) q[3];
sx q[3];
rz(-1.2069586) q[3];
sx q[3];
rz(-1.1999493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46362296) q[0];
sx q[0];
rz(-0.93358731) q[0];
sx q[0];
rz(-1.8883702) q[0];
rz(2.5190917) q[1];
sx q[1];
rz(-1.6815192) q[1];
sx q[1];
rz(1.1463096) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2749467) q[0];
sx q[0];
rz(-1.6410315) q[0];
sx q[0];
rz(0.74830351) q[0];
rz(-pi) q[1];
rz(0.88845466) q[2];
sx q[2];
rz(-2.1092215) q[2];
sx q[2];
rz(-0.74935645) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6282564) q[1];
sx q[1];
rz(-0.69111052) q[1];
sx q[1];
rz(-1.6860794) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1504437) q[3];
sx q[3];
rz(-2.162809) q[3];
sx q[3];
rz(2.1136485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.15381947) q[2];
sx q[2];
rz(-2.1058857) q[2];
sx q[2];
rz(-1.0158319) q[2];
rz(2.231797) q[3];
sx q[3];
rz(-2.4715021) q[3];
sx q[3];
rz(2.773496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1290454) q[0];
sx q[0];
rz(-2.5279901) q[0];
sx q[0];
rz(-3.1273499) q[0];
rz(2.2968538) q[1];
sx q[1];
rz(-0.95294398) q[1];
sx q[1];
rz(1.7600118) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89109126) q[0];
sx q[0];
rz(-1.4282465) q[0];
sx q[0];
rz(3.13091) q[0];
rz(-pi) q[1];
rz(1.472677) q[2];
sx q[2];
rz(-0.58110229) q[2];
sx q[2];
rz(-0.89400089) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6784918) q[1];
sx q[1];
rz(-1.4271724) q[1];
sx q[1];
rz(-2.1302057) q[1];
rz(1.7900449) q[3];
sx q[3];
rz(-2.0937113) q[3];
sx q[3];
rz(-0.067283665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0649197) q[2];
sx q[2];
rz(-1.2023456) q[2];
sx q[2];
rz(2.1515965) q[2];
rz(-0.89407095) q[3];
sx q[3];
rz(-2.6529513) q[3];
sx q[3];
rz(0.22542424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
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
rz(2.0902696) q[0];
sx q[0];
rz(-0.6482424) q[0];
sx q[0];
rz(2.0363664) q[0];
rz(-1.3399711) q[1];
sx q[1];
rz(-2.5201288) q[1];
sx q[1];
rz(-2.7531243) q[1];
rz(-0.52836616) q[2];
sx q[2];
rz(-2.3522204) q[2];
sx q[2];
rz(1.5632202) q[2];
rz(-0.43511012) q[3];
sx q[3];
rz(-0.995244) q[3];
sx q[3];
rz(1.4959195) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
