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
rz(1.7460495) q[1];
sx q[1];
rz(-1.6092602) q[1];
sx q[1];
rz(-1.2083763) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3291408) q[0];
sx q[0];
rz(-1.666823) q[0];
sx q[0];
rz(-0.046669331) q[0];
rz(0.12579972) q[2];
sx q[2];
rz(-1.5828136) q[2];
sx q[2];
rz(-1.8819686) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0304058) q[1];
sx q[1];
rz(-0.62796794) q[1];
sx q[1];
rz(2.9524132) q[1];
rz(-pi) q[2];
rz(-0.73322202) q[3];
sx q[3];
rz(-0.65922046) q[3];
sx q[3];
rz(-0.99430195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8157114) q[2];
sx q[2];
rz(-1.4006961) q[2];
sx q[2];
rz(-1.6750083) q[2];
rz(-2.4438434) q[3];
sx q[3];
rz(-1.1013228) q[3];
sx q[3];
rz(2.3944323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.117347) q[0];
sx q[0];
rz(-2.0139366) q[0];
sx q[0];
rz(-1.9673989) q[0];
rz(2.970447) q[1];
sx q[1];
rz(-2.0967963) q[1];
sx q[1];
rz(-0.29719621) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3127808) q[0];
sx q[0];
rz(-0.096185616) q[0];
sx q[0];
rz(-1.0440774) q[0];
rz(0.53571312) q[2];
sx q[2];
rz(-2.3539054) q[2];
sx q[2];
rz(-2.1802769) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9651523) q[1];
sx q[1];
rz(-2.0808176) q[1];
sx q[1];
rz(2.7494207) q[1];
x q[2];
rz(-2.3791749) q[3];
sx q[3];
rz(-1.8652417) q[3];
sx q[3];
rz(-0.03014119) q[3];
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
rz(-0.87614122) q[3];
sx q[3];
rz(-0.66771475) q[3];
sx q[3];
rz(0.63703018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0456332) q[0];
sx q[0];
rz(-1.2852083) q[0];
sx q[0];
rz(-0.30763787) q[0];
rz(-2.3930507) q[1];
sx q[1];
rz(-2.8100439) q[1];
sx q[1];
rz(-2.3017853) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9899983) q[0];
sx q[0];
rz(-0.50938207) q[0];
sx q[0];
rz(-2.529782) q[0];
x q[1];
rz(2.8022021) q[2];
sx q[2];
rz(-1.4815785) q[2];
sx q[2];
rz(0.81775507) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0963904) q[1];
sx q[1];
rz(-2.5844816) q[1];
sx q[1];
rz(0.60369173) q[1];
rz(0.93785502) q[3];
sx q[3];
rz(-1.7114534) q[3];
sx q[3];
rz(1.7733639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2741189) q[2];
sx q[2];
rz(-1.3511191) q[2];
sx q[2];
rz(1.8236558) q[2];
rz(1.9258202) q[3];
sx q[3];
rz(-2.7850745) q[3];
sx q[3];
rz(-1.6320451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.3777305) q[0];
sx q[0];
rz(-1.7611935) q[0];
sx q[0];
rz(-0.44556251) q[0];
rz(2.5207649) q[1];
sx q[1];
rz(-1.8955684) q[1];
sx q[1];
rz(-0.96558085) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27698101) q[0];
sx q[0];
rz(-2.3267713) q[0];
sx q[0];
rz(-0.38397249) q[0];
x q[1];
rz(0.031164073) q[2];
sx q[2];
rz(-0.84261299) q[2];
sx q[2];
rz(0.63809168) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6009778) q[1];
sx q[1];
rz(-1.317273) q[1];
sx q[1];
rz(3.1086604) q[1];
rz(-1.4569034) q[3];
sx q[3];
rz(-1.4741065) q[3];
sx q[3];
rz(-2.8817122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.821637) q[2];
sx q[2];
rz(-0.76239061) q[2];
sx q[2];
rz(0.92932534) q[2];
rz(-0.6435414) q[3];
sx q[3];
rz(-1.0361592) q[3];
sx q[3];
rz(-2.1381366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4716361) q[0];
sx q[0];
rz(-1.5820553) q[0];
sx q[0];
rz(1.8575645) q[0];
rz(-0.28981003) q[1];
sx q[1];
rz(-0.73957864) q[1];
sx q[1];
rz(-1.0481542) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7832344) q[0];
sx q[0];
rz(-0.90792197) q[0];
sx q[0];
rz(0.62028424) q[0];
rz(-pi) q[1];
rz(-0.94190188) q[2];
sx q[2];
rz(-0.84398979) q[2];
sx q[2];
rz(-0.84743365) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.84761274) q[1];
sx q[1];
rz(-1.249053) q[1];
sx q[1];
rz(1.6683679) q[1];
rz(-pi) q[2];
rz(-2.2506511) q[3];
sx q[3];
rz(-1.1853301) q[3];
sx q[3];
rz(0.95213529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8706878) q[2];
sx q[2];
rz(-2.1257766) q[2];
sx q[2];
rz(2.2407545) q[2];
rz(1.0926931) q[3];
sx q[3];
rz(-1.0034424) q[3];
sx q[3];
rz(1.9074915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1822405) q[0];
sx q[0];
rz(-1.2077967) q[0];
sx q[0];
rz(0.59610468) q[0];
rz(1.4959363) q[1];
sx q[1];
rz(-2.2438965) q[1];
sx q[1];
rz(-1.2449107) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7398864) q[0];
sx q[0];
rz(-2.5283928) q[0];
sx q[0];
rz(-0.99341157) q[0];
rz(-pi) q[1];
rz(0.55576022) q[2];
sx q[2];
rz(-0.55609497) q[2];
sx q[2];
rz(-0.13242002) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.53828) q[1];
sx q[1];
rz(-2.0671751) q[1];
sx q[1];
rz(-2.3904295) q[1];
rz(-pi) q[2];
rz(2.3915668) q[3];
sx q[3];
rz(-2.8083028) q[3];
sx q[3];
rz(-0.94543524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9101377) q[2];
sx q[2];
rz(-0.69880501) q[2];
sx q[2];
rz(-0.22496741) q[2];
rz(0.088430017) q[3];
sx q[3];
rz(-1.6902573) q[3];
sx q[3];
rz(-2.6627873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6769619) q[0];
sx q[0];
rz(-1.2372274) q[0];
sx q[0];
rz(-2.8616469) q[0];
rz(1.6784558) q[1];
sx q[1];
rz(-1.2660374) q[1];
sx q[1];
rz(-2.8889012) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1278909) q[0];
sx q[0];
rz(-1.2780006) q[0];
sx q[0];
rz(-0.049833628) q[0];
rz(-pi) q[1];
rz(1.5451317) q[2];
sx q[2];
rz(-0.63458323) q[2];
sx q[2];
rz(2.6013825) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6934769) q[1];
sx q[1];
rz(-2.2709393) q[1];
sx q[1];
rz(0.6853939) q[1];
rz(-pi) q[2];
x q[2];
rz(0.85315506) q[3];
sx q[3];
rz(-0.95082885) q[3];
sx q[3];
rz(2.6593047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.32020405) q[2];
sx q[2];
rz(-0.89367047) q[2];
sx q[2];
rz(-1.6097216) q[2];
rz(-1.948471) q[3];
sx q[3];
rz(-2.2333998) q[3];
sx q[3];
rz(1.0866722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10483345) q[0];
sx q[0];
rz(-1.3830673) q[0];
sx q[0];
rz(1.4861134) q[0];
rz(-2.8727818) q[1];
sx q[1];
rz(-2.0116282) q[1];
sx q[1];
rz(-0.2789467) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3725961) q[0];
sx q[0];
rz(-1.0770174) q[0];
sx q[0];
rz(2.2934224) q[0];
rz(-1.2645623) q[2];
sx q[2];
rz(-2.6191204) q[2];
sx q[2];
rz(-2.317121) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5178691) q[1];
sx q[1];
rz(-1.9451127) q[1];
sx q[1];
rz(2.3831297) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.14848498) q[3];
sx q[3];
rz(-1.0035702) q[3];
sx q[3];
rz(-2.4384769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.802861) q[2];
sx q[2];
rz(-1.4638476) q[2];
sx q[2];
rz(-1.5926682) q[2];
rz(-1.0507978) q[3];
sx q[3];
rz(-1.2069586) q[3];
sx q[3];
rz(-1.9416434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6779697) q[0];
sx q[0];
rz(-2.2080053) q[0];
sx q[0];
rz(1.2532225) q[0];
rz(2.5190917) q[1];
sx q[1];
rz(-1.4600735) q[1];
sx q[1];
rz(1.995283) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3619986) q[0];
sx q[0];
rz(-0.75095526) q[0];
sx q[0];
rz(-3.0385618) q[0];
rz(-2.253138) q[2];
sx q[2];
rz(-2.1092215) q[2];
sx q[2];
rz(-0.74935645) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6282564) q[1];
sx q[1];
rz(-2.4504821) q[1];
sx q[1];
rz(1.4555132) q[1];
rz(1.1504437) q[3];
sx q[3];
rz(-2.162809) q[3];
sx q[3];
rz(2.1136485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9877732) q[2];
sx q[2];
rz(-1.035707) q[2];
sx q[2];
rz(2.1257607) q[2];
rz(2.231797) q[3];
sx q[3];
rz(-0.67009059) q[3];
sx q[3];
rz(-2.773496) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1290454) q[0];
sx q[0];
rz(-0.61360252) q[0];
sx q[0];
rz(-3.1273499) q[0];
rz(2.2968538) q[1];
sx q[1];
rz(-0.95294398) q[1];
sx q[1];
rz(1.7600118) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.325557) q[0];
sx q[0];
rz(-2.9986458) q[0];
sx q[0];
rz(-1.6450892) q[0];
rz(-2.1496885) q[2];
sx q[2];
rz(-1.5169946) q[2];
sx q[2];
rz(0.75888854) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8091619) q[1];
sx q[1];
rz(-0.5756439) q[1];
sx q[1];
rz(1.3047421) q[1];
x q[2];
rz(2.6081309) q[3];
sx q[3];
rz(-1.3812314) q[3];
sx q[3];
rz(-1.61434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0766729) q[2];
sx q[2];
rz(-1.939247) q[2];
sx q[2];
rz(0.98999611) q[2];
rz(0.89407095) q[3];
sx q[3];
rz(-0.48864135) q[3];
sx q[3];
rz(-2.9161684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
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
rz(2.6132265) q[2];
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