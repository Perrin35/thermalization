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
rz(1.7460495) q[1];
sx q[1];
rz(4.6739251) q[1];
sx q[1];
rz(8.2164017) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87586227) q[0];
sx q[0];
rz(-3.0348572) q[0];
sx q[0];
rz(2.021832) q[0];
rz(-3/(10*pi)) q[2];
sx q[2];
rz(-3.0152233) q[2];
sx q[2];
rz(0.40590826) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2627416) q[1];
sx q[1];
rz(-0.95572119) q[1];
sx q[1];
rz(-1.7064852) q[1];
rz(0.73322202) q[3];
sx q[3];
rz(-0.65922046) q[3];
sx q[3];
rz(-2.1472907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3258813) q[2];
sx q[2];
rz(-1.4006961) q[2];
sx q[2];
rz(-1.6750083) q[2];
rz(-0.69774929) q[3];
sx q[3];
rz(-2.0402699) q[3];
sx q[3];
rz(-0.74716032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0242457) q[0];
sx q[0];
rz(-1.1276561) q[0];
sx q[0];
rz(-1.9673989) q[0];
rz(2.970447) q[1];
sx q[1];
rz(-2.0967963) q[1];
sx q[1];
rz(2.8443964) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78404616) q[0];
sx q[0];
rz(-1.4876801) q[0];
sx q[0];
rz(-0.048464171) q[0];
rz(0.71248033) q[2];
sx q[2];
rz(-1.2006294) q[2];
sx q[2];
rz(-2.9287101) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9651523) q[1];
sx q[1];
rz(-1.0607751) q[1];
sx q[1];
rz(-0.39217197) q[1];
rz(-pi) q[2];
x q[2];
rz(0.76241775) q[3];
sx q[3];
rz(-1.276351) q[3];
sx q[3];
rz(-3.1114515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.16195665) q[2];
sx q[2];
rz(-1.1844144) q[2];
sx q[2];
rz(0.61398181) q[2];
rz(0.87614122) q[3];
sx q[3];
rz(-0.66771475) q[3];
sx q[3];
rz(-0.63703018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0959594) q[0];
sx q[0];
rz(-1.2852083) q[0];
sx q[0];
rz(-0.30763787) q[0];
rz(-0.74854198) q[1];
sx q[1];
rz(-2.8100439) q[1];
sx q[1];
rz(2.3017853) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15159431) q[0];
sx q[0];
rz(-0.50938207) q[0];
sx q[0];
rz(-0.61181061) q[0];
x q[1];
rz(0.26250458) q[2];
sx q[2];
rz(-2.7911107) q[2];
sx q[2];
rz(2.1413435) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0864598) q[1];
sx q[1];
rz(-1.8756525) q[1];
sx q[1];
rz(2.6677368) q[1];
rz(-pi) q[2];
x q[2];
rz(0.17383667) q[3];
sx q[3];
rz(-0.94508119) q[3];
sx q[3];
rz(3.0415149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
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
rz(-0.62082779) q[1];
sx q[1];
rz(-1.8955684) q[1];
sx q[1];
rz(-0.96558085) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25528204) q[0];
sx q[0];
rz(-0.83034407) q[0];
sx q[0];
rz(-1.9489954) q[0];
rz(0.031164073) q[2];
sx q[2];
rz(-0.84261299) q[2];
sx q[2];
rz(-2.503501) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6009778) q[1];
sx q[1];
rz(-1.8243196) q[1];
sx q[1];
rz(3.1086604) q[1];
rz(1.4569034) q[3];
sx q[3];
rz(-1.4741065) q[3];
sx q[3];
rz(2.8817122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3199557) q[2];
sx q[2];
rz(-2.379202) q[2];
sx q[2];
rz(0.92932534) q[2];
rz(-2.4980513) q[3];
sx q[3];
rz(-1.0361592) q[3];
sx q[3];
rz(2.1381366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6699566) q[0];
sx q[0];
rz(-1.5595373) q[0];
sx q[0];
rz(-1.8575645) q[0];
rz(-0.28981003) q[1];
sx q[1];
rz(-2.402014) q[1];
sx q[1];
rz(1.0481542) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5149887) q[0];
sx q[0];
rz(-1.0948613) q[0];
sx q[0];
rz(0.80608741) q[0];
rz(-pi) q[1];
x q[1];
rz(0.94190188) q[2];
sx q[2];
rz(-0.84398979) q[2];
sx q[2];
rz(0.84743365) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.594164) q[1];
sx q[1];
rz(-2.8058726) q[1];
sx q[1];
rz(-0.28433849) q[1];
rz(-pi) q[2];
rz(0.9976451) q[3];
sx q[3];
rz(-0.76612681) q[3];
sx q[3];
rz(-1.0539953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2709048) q[2];
sx q[2];
rz(-1.015816) q[2];
sx q[2];
rz(-0.90083814) q[2];
rz(1.0926931) q[3];
sx q[3];
rz(-2.1381502) q[3];
sx q[3];
rz(1.2341011) q[3];
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
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1822405) q[0];
sx q[0];
rz(-1.2077967) q[0];
sx q[0];
rz(-2.545488) q[0];
rz(1.4959363) q[1];
sx q[1];
rz(-0.89769617) q[1];
sx q[1];
rz(1.2449107) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7398864) q[0];
sx q[0];
rz(-2.5283928) q[0];
sx q[0];
rz(-2.1481811) q[0];
rz(-pi) q[1];
rz(-0.4857829) q[2];
sx q[2];
rz(-1.2885639) q[2];
sx q[2];
rz(2.1886052) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4391172) q[1];
sx q[1];
rz(-2.2687952) q[1];
sx q[1];
rz(-2.4707787) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8935029) q[3];
sx q[3];
rz(-1.7956942) q[3];
sx q[3];
rz(0.096506462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9101377) q[2];
sx q[2];
rz(-2.4427876) q[2];
sx q[2];
rz(-0.22496741) q[2];
rz(-0.088430017) q[3];
sx q[3];
rz(-1.4513353) q[3];
sx q[3];
rz(0.47880539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6769619) q[0];
sx q[0];
rz(-1.2372274) q[0];
sx q[0];
rz(-0.27994573) q[0];
rz(1.4631368) q[1];
sx q[1];
rz(-1.8755553) q[1];
sx q[1];
rz(-2.8889012) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42851105) q[0];
sx q[0];
rz(-1.5230852) q[0];
sx q[0];
rz(-1.8639355) q[0];
rz(0.018888868) q[2];
sx q[2];
rz(-0.9364555) q[2];
sx q[2];
rz(-0.50834507) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.44811571) q[1];
sx q[1];
rz(-0.87065334) q[1];
sx q[1];
rz(-0.6853939) q[1];
rz(-pi) q[2];
rz(-2.3971862) q[3];
sx q[3];
rz(-0.91074569) q[3];
sx q[3];
rz(0.50124121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8213886) q[2];
sx q[2];
rz(-0.89367047) q[2];
sx q[2];
rz(-1.5318711) q[2];
rz(1.1931217) q[3];
sx q[3];
rz(-0.90819287) q[3];
sx q[3];
rz(2.0549205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0367592) q[0];
sx q[0];
rz(-1.3830673) q[0];
sx q[0];
rz(-1.6554792) q[0];
rz(-0.2688109) q[1];
sx q[1];
rz(-2.0116282) q[1];
sx q[1];
rz(0.2789467) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3725961) q[0];
sx q[0];
rz(-1.0770174) q[0];
sx q[0];
rz(-2.2934224) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9697044) q[2];
sx q[2];
rz(-2.066678) q[2];
sx q[2];
rz(-1.1743197) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.62372359) q[1];
sx q[1];
rz(-1.1964799) q[1];
sx q[1];
rz(2.3831297) q[1];
x q[2];
rz(-0.14848498) q[3];
sx q[3];
rz(-1.0035702) q[3];
sx q[3];
rz(0.70311577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.802861) q[2];
sx q[2];
rz(-1.4638476) q[2];
sx q[2];
rz(1.5926682) q[2];
rz(2.0907949) q[3];
sx q[3];
rz(-1.934634) q[3];
sx q[3];
rz(-1.1999493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6779697) q[0];
sx q[0];
rz(-0.93358731) q[0];
sx q[0];
rz(1.8883702) q[0];
rz(-0.62250096) q[1];
sx q[1];
rz(-1.6815192) q[1];
sx q[1];
rz(1.1463096) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.779594) q[0];
sx q[0];
rz(-0.75095526) q[0];
sx q[0];
rz(-3.0385618) q[0];
x q[1];
rz(0.88845466) q[2];
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
sx q[0];
rz(pi/2) q[0];
rz(1.6282564) q[1];
sx q[1];
rz(-2.4504821) q[1];
sx q[1];
rz(1.4555132) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5961612) q[3];
sx q[3];
rz(-2.430393) q[3];
sx q[3];
rz(1.4382854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.15381947) q[2];
sx q[2];
rz(-2.1058857) q[2];
sx q[2];
rz(2.1257607) q[2];
rz(-0.9097957) q[3];
sx q[3];
rz(-2.4715021) q[3];
sx q[3];
rz(-0.36809665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0125473) q[0];
sx q[0];
rz(-0.61360252) q[0];
sx q[0];
rz(-3.1273499) q[0];
rz(-0.8447389) q[1];
sx q[1];
rz(-2.1886487) q[1];
sx q[1];
rz(-1.7600118) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.325557) q[0];
sx q[0];
rz(-0.14294681) q[0];
sx q[0];
rz(-1.6450892) q[0];
rz(2.1496885) q[2];
sx q[2];
rz(-1.5169946) q[2];
sx q[2];
rz(2.3827041) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4631008) q[1];
sx q[1];
rz(-1.7144202) q[1];
sx q[1];
rz(-1.0113869) q[1];
rz(-pi) q[2];
rz(-0.36079447) q[3];
sx q[3];
rz(-2.5785355) q[3];
sx q[3];
rz(0.35239708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0649197) q[2];
sx q[2];
rz(-1.2023456) q[2];
sx q[2];
rz(-2.1515965) q[2];
rz(-0.89407095) q[3];
sx q[3];
rz(-2.6529513) q[3];
sx q[3];
rz(-2.9161684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0902696) q[0];
sx q[0];
rz(-0.6482424) q[0];
sx q[0];
rz(2.0363664) q[0];
rz(1.3399711) q[1];
sx q[1];
rz(-0.62146386) q[1];
sx q[1];
rz(0.38846831) q[1];
rz(2.0409394) q[2];
sx q[2];
rz(-2.2307776) q[2];
sx q[2];
rz(-2.270436) q[2];
rz(2.7064825) q[3];
sx q[3];
rz(-0.995244) q[3];
sx q[3];
rz(1.4959195) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
