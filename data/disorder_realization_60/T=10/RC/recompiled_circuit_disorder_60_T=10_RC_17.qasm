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
rz(3.0728683) q[0];
rz(-1.3955431) q[1];
sx q[1];
rz(-1.5323324) q[1];
sx q[1];
rz(-1.9332164) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24613334) q[0];
sx q[0];
rz(-1.6172505) q[0];
sx q[0];
rz(1.4746656) q[0];
rz(3.0460998) q[2];
sx q[2];
rz(-0.12636939) q[2];
sx q[2];
rz(-0.40590826) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2627416) q[1];
sx q[1];
rz(-0.95572119) q[1];
sx q[1];
rz(-1.7064852) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0923907) q[3];
sx q[3];
rz(-1.0983101) q[3];
sx q[3];
rz(2.9977968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3258813) q[2];
sx q[2];
rz(-1.4006961) q[2];
sx q[2];
rz(-1.4665843) q[2];
rz(0.69774929) q[3];
sx q[3];
rz(-2.0402699) q[3];
sx q[3];
rz(-2.3944323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.117347) q[0];
sx q[0];
rz(-2.0139366) q[0];
sx q[0];
rz(-1.9673989) q[0];
rz(0.17114561) q[1];
sx q[1];
rz(-1.0447964) q[1];
sx q[1];
rz(-0.29719621) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78404616) q[0];
sx q[0];
rz(-1.6539126) q[0];
sx q[0];
rz(0.048464171) q[0];
x q[1];
rz(-0.53571312) q[2];
sx q[2];
rz(-0.78768724) q[2];
sx q[2];
rz(0.96131575) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9651523) q[1];
sx q[1];
rz(-2.0808176) q[1];
sx q[1];
rz(-2.7494207) q[1];
rz(-pi) q[2];
rz(-2.3791749) q[3];
sx q[3];
rz(-1.276351) q[3];
sx q[3];
rz(0.03014119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.16195665) q[2];
sx q[2];
rz(-1.1844144) q[2];
sx q[2];
rz(-2.5276108) q[2];
rz(-2.2654514) q[3];
sx q[3];
rz(-0.66771475) q[3];
sx q[3];
rz(-0.63703018) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0456332) q[0];
sx q[0];
rz(-1.2852083) q[0];
sx q[0];
rz(2.8339548) q[0];
rz(-0.74854198) q[1];
sx q[1];
rz(-2.8100439) q[1];
sx q[1];
rz(-0.83980733) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9899983) q[0];
sx q[0];
rz(-0.50938207) q[0];
sx q[0];
rz(2.529782) q[0];
rz(-pi) q[1];
x q[1];
rz(1.66538) q[2];
sx q[2];
rz(-1.2328096) q[2];
sx q[2];
rz(-0.72159492) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0452022) q[1];
sx q[1];
rz(-2.5844816) q[1];
sx q[1];
rz(-0.60369173) q[1];
rz(1.3358467) q[3];
sx q[3];
rz(-0.64628212) q[3];
sx q[3];
rz(-0.39138734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2741189) q[2];
sx q[2];
rz(-1.3511191) q[2];
sx q[2];
rz(-1.8236558) q[2];
rz(-1.9258202) q[3];
sx q[3];
rz(-0.35651818) q[3];
sx q[3];
rz(-1.6320451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76386219) q[0];
sx q[0];
rz(-1.3803991) q[0];
sx q[0];
rz(-2.6960301) q[0];
rz(-0.62082779) q[1];
sx q[1];
rz(-1.8955684) q[1];
sx q[1];
rz(-0.96558085) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8646116) q[0];
sx q[0];
rz(-0.81482139) q[0];
sx q[0];
rz(-0.38397249) q[0];
rz(-0.84237174) q[2];
sx q[2];
rz(-1.5475376) q[2];
sx q[2];
rz(-2.1881441) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4703776) q[1];
sx q[1];
rz(-0.2556076) q[1];
sx q[1];
rz(1.6971991) q[1];
rz(-pi) q[2];
x q[2];
rz(0.097316381) q[3];
sx q[3];
rz(-1.4574377) q[3];
sx q[3];
rz(-1.8196343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3199557) q[2];
sx q[2];
rz(-0.76239061) q[2];
sx q[2];
rz(-0.92932534) q[2];
rz(-2.4980513) q[3];
sx q[3];
rz(-2.1054335) q[3];
sx q[3];
rz(1.003456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6699566) q[0];
sx q[0];
rz(-1.5820553) q[0];
sx q[0];
rz(1.2840282) q[0];
rz(2.8517826) q[1];
sx q[1];
rz(-0.73957864) q[1];
sx q[1];
rz(-1.0481542) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6402733) q[0];
sx q[0];
rz(-0.87448705) q[0];
sx q[0];
rz(2.2107844) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1996908) q[2];
sx q[2];
rz(-2.2976029) q[2];
sx q[2];
rz(-2.294159) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2939799) q[1];
sx q[1];
rz(-1.8925397) q[1];
sx q[1];
rz(1.4732248) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.9976451) q[3];
sx q[3];
rz(-2.3754658) q[3];
sx q[3];
rz(-1.0539953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8706878) q[2];
sx q[2];
rz(-1.015816) q[2];
sx q[2];
rz(-2.2407545) q[2];
rz(-2.0488996) q[3];
sx q[3];
rz(-2.1381502) q[3];
sx q[3];
rz(-1.9074915) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1822405) q[0];
sx q[0];
rz(-1.933796) q[0];
sx q[0];
rz(0.59610468) q[0];
rz(-1.6456564) q[1];
sx q[1];
rz(-0.89769617) q[1];
sx q[1];
rz(1.2449107) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4017063) q[0];
sx q[0];
rz(-0.6131999) q[0];
sx q[0];
rz(-0.99341157) q[0];
x q[1];
rz(-0.4857829) q[2];
sx q[2];
rz(-1.8530288) q[2];
sx q[2];
rz(-2.1886052) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5925496) q[1];
sx q[1];
rz(-0.92714308) q[1];
sx q[1];
rz(2.2085269) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.2480898) q[3];
sx q[3];
rz(-1.3458985) q[3];
sx q[3];
rz(0.096506462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9101377) q[2];
sx q[2];
rz(-2.4427876) q[2];
sx q[2];
rz(0.22496741) q[2];
rz(-0.088430017) q[3];
sx q[3];
rz(-1.4513353) q[3];
sx q[3];
rz(-2.6627873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6769619) q[0];
sx q[0];
rz(-1.2372274) q[0];
sx q[0];
rz(2.8616469) q[0];
rz(1.4631368) q[1];
sx q[1];
rz(-1.2660374) q[1];
sx q[1];
rz(-0.25269145) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8425927) q[0];
sx q[0];
rz(-2.8447066) q[0];
sx q[0];
rz(-1.7345558) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2052223) q[2];
sx q[2];
rz(-1.5555824) q[2];
sx q[2];
rz(-2.0903367) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6934769) q[1];
sx q[1];
rz(-0.87065334) q[1];
sx q[1];
rz(-0.6853939) q[1];
rz(-pi) q[2];
rz(0.75848363) q[3];
sx q[3];
rz(-2.1355724) q[3];
sx q[3];
rz(-1.5578711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8213886) q[2];
sx q[2];
rz(-0.89367047) q[2];
sx q[2];
rz(-1.6097216) q[2];
rz(1.948471) q[3];
sx q[3];
rz(-0.90819287) q[3];
sx q[3];
rz(-2.0549205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10483345) q[0];
sx q[0];
rz(-1.7585254) q[0];
sx q[0];
rz(-1.6554792) q[0];
rz(0.2688109) q[1];
sx q[1];
rz(-2.0116282) q[1];
sx q[1];
rz(2.862646) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9439518) q[0];
sx q[0];
rz(-0.949172) q[0];
sx q[0];
rz(-0.62244121) q[0];
x q[1];
rz(1.8770304) q[2];
sx q[2];
rz(-0.5224723) q[2];
sx q[2];
rz(2.317121) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.61356269) q[1];
sx q[1];
rz(-0.87606214) q[1];
sx q[1];
rz(-2.0668405) q[1];
x q[2];
rz(-2.143048) q[3];
sx q[3];
rz(-1.4456985) q[3];
sx q[3];
rz(-2.3541114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.802861) q[2];
sx q[2];
rz(-1.4638476) q[2];
sx q[2];
rz(1.5489244) q[2];
rz(1.0507978) q[3];
sx q[3];
rz(-1.934634) q[3];
sx q[3];
rz(1.1999493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6779697) q[0];
sx q[0];
rz(-0.93358731) q[0];
sx q[0];
rz(-1.2532225) q[0];
rz(-0.62250096) q[1];
sx q[1];
rz(-1.6815192) q[1];
sx q[1];
rz(1.1463096) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86664591) q[0];
sx q[0];
rz(-1.6410315) q[0];
sx q[0];
rz(-0.74830351) q[0];
rz(-pi) q[1];
x q[1];
rz(0.65593221) q[2];
sx q[2];
rz(-0.99870517) q[2];
sx q[2];
rz(2.7149372) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5133363) q[1];
sx q[1];
rz(-2.4504821) q[1];
sx q[1];
rz(-1.4555132) q[1];
rz(1.9911489) q[3];
sx q[3];
rz(-0.97878362) q[3];
sx q[3];
rz(-1.0279442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.15381947) q[2];
sx q[2];
rz(-1.035707) q[2];
sx q[2];
rz(1.0158319) q[2];
rz(-0.9097957) q[3];
sx q[3];
rz(-2.4715021) q[3];
sx q[3];
rz(-0.36809665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0125473) q[0];
sx q[0];
rz(-0.61360252) q[0];
sx q[0];
rz(3.1273499) q[0];
rz(0.8447389) q[1];
sx q[1];
rz(-0.95294398) q[1];
sx q[1];
rz(1.3815809) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81603564) q[0];
sx q[0];
rz(-2.9986458) q[0];
sx q[0];
rz(-1.6450892) q[0];
rz(-pi) q[1];
rz(-0.99190418) q[2];
sx q[2];
rz(-1.6245981) q[2];
sx q[2];
rz(0.75888854) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.018316293) q[1];
sx q[1];
rz(-2.1237719) q[1];
sx q[1];
rz(-2.9725914) q[1];
x q[2];
rz(-1.7900449) q[3];
sx q[3];
rz(-1.0478813) q[3];
sx q[3];
rz(-0.067283665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0649197) q[2];
sx q[2];
rz(-1.2023456) q[2];
sx q[2];
rz(0.98999611) q[2];
rz(-0.89407095) q[3];
sx q[3];
rz(-2.6529513) q[3];
sx q[3];
rz(-2.9161684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(1.8016215) q[1];
sx q[1];
rz(-2.5201288) q[1];
sx q[1];
rz(-2.7531243) q[1];
rz(0.52836616) q[2];
sx q[2];
rz(-0.7893723) q[2];
sx q[2];
rz(-1.5783725) q[2];
rz(2.1918478) q[3];
sx q[3];
rz(-1.2093778) q[3];
sx q[3];
rz(-2.9686684) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
