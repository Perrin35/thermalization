OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.94937593) q[0];
sx q[0];
rz(-1.047171) q[0];
sx q[0];
rz(0.068724364) q[0];
rz(-1.3955431) q[1];
sx q[1];
rz(-1.5323324) q[1];
sx q[1];
rz(-1.9332164) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24613334) q[0];
sx q[0];
rz(-1.6172505) q[0];
sx q[0];
rz(-1.4746656) q[0];
x q[1];
rz(-3.0157929) q[2];
sx q[2];
rz(-1.558779) q[2];
sx q[2];
rz(-1.2596241) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2627416) q[1];
sx q[1];
rz(-2.1858715) q[1];
sx q[1];
rz(-1.4351074) q[1];
x q[2];
rz(-0.73322202) q[3];
sx q[3];
rz(-2.4823722) q[3];
sx q[3];
rz(0.99430195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3258813) q[2];
sx q[2];
rz(-1.7408966) q[2];
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
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.117347) q[0];
sx q[0];
rz(-2.0139366) q[0];
sx q[0];
rz(-1.9673989) q[0];
rz(-2.970447) q[1];
sx q[1];
rz(-1.0447964) q[1];
sx q[1];
rz(2.8443964) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78272351) q[0];
sx q[0];
rz(-1.5224996) q[0];
sx q[0];
rz(1.6540098) q[0];
rz(-pi) q[1];
rz(0.53571312) q[2];
sx q[2];
rz(-0.78768724) q[2];
sx q[2];
rz(2.1802769) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.593593) q[1];
sx q[1];
rz(-1.2307234) q[1];
sx q[1];
rz(2.1151357) q[1];
rz(-pi) q[2];
rz(2.3791749) q[3];
sx q[3];
rz(-1.8652417) q[3];
sx q[3];
rz(0.03014119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.979636) q[2];
sx q[2];
rz(-1.1844144) q[2];
sx q[2];
rz(-0.61398181) q[2];
rz(0.87614122) q[3];
sx q[3];
rz(-0.66771475) q[3];
sx q[3];
rz(-0.63703018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0959594) q[0];
sx q[0];
rz(-1.2852083) q[0];
sx q[0];
rz(-2.8339548) q[0];
rz(-0.74854198) q[1];
sx q[1];
rz(-0.33154878) q[1];
sx q[1];
rz(-2.3017853) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15159431) q[0];
sx q[0];
rz(-0.50938207) q[0];
sx q[0];
rz(-0.61181061) q[0];
rz(2.8790881) q[2];
sx q[2];
rz(-0.35048198) q[2];
sx q[2];
rz(2.1413435) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3629344) q[1];
sx q[1];
rz(-2.0211453) q[1];
sx q[1];
rz(-1.2308916) q[1];
rz(2.2037376) q[3];
sx q[3];
rz(-1.4301392) q[3];
sx q[3];
rz(-1.3682287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.86747375) q[2];
sx q[2];
rz(-1.7904736) q[2];
sx q[2];
rz(-1.8236558) q[2];
rz(-1.2157724) q[3];
sx q[3];
rz(-0.35651818) q[3];
sx q[3];
rz(-1.5095476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-0.76386219) q[0];
sx q[0];
rz(-1.7611935) q[0];
sx q[0];
rz(-0.44556251) q[0];
rz(0.62082779) q[1];
sx q[1];
rz(-1.8955684) q[1];
sx q[1];
rz(-2.1760118) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8863106) q[0];
sx q[0];
rz(-0.83034407) q[0];
sx q[0];
rz(1.9489954) q[0];
x q[1];
rz(0.031164073) q[2];
sx q[2];
rz(-0.84261299) q[2];
sx q[2];
rz(0.63809168) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.1196739) q[1];
sx q[1];
rz(-1.5389171) q[1];
sx q[1];
rz(1.3171413) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4569034) q[3];
sx q[3];
rz(-1.4741065) q[3];
sx q[3];
rz(-0.25988042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.821637) q[2];
sx q[2];
rz(-2.379202) q[2];
sx q[2];
rz(-0.92932534) q[2];
rz(0.6435414) q[3];
sx q[3];
rz(-1.0361592) q[3];
sx q[3];
rz(2.1381366) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6699566) q[0];
sx q[0];
rz(-1.5595373) q[0];
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
rz(0.50131932) q[0];
sx q[0];
rz(-0.87448705) q[0];
sx q[0];
rz(2.2107844) q[0];
rz(-pi) q[1];
rz(-0.94190188) q[2];
sx q[2];
rz(-2.2976029) q[2];
sx q[2];
rz(0.84743365) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3874665) q[1];
sx q[1];
rz(-1.6633463) q[1];
sx q[1];
rz(2.8184163) q[1];
rz(0.89094152) q[3];
sx q[3];
rz(-1.1853301) q[3];
sx q[3];
rz(0.95213529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8706878) q[2];
sx q[2];
rz(-1.015816) q[2];
sx q[2];
rz(-0.90083814) q[2];
rz(-2.0488996) q[3];
sx q[3];
rz(-2.1381502) q[3];
sx q[3];
rz(1.2341011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95935217) q[0];
sx q[0];
rz(-1.933796) q[0];
sx q[0];
rz(-2.545488) q[0];
rz(-1.6456564) q[1];
sx q[1];
rz(-0.89769617) q[1];
sx q[1];
rz(1.2449107) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0743474) q[0];
sx q[0];
rz(-1.0676358) q[0];
sx q[0];
rz(2.7748681) q[0];
rz(1.2539358) q[2];
sx q[2];
rz(-1.1057901) q[2];
sx q[2];
rz(2.3777865) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7024755) q[1];
sx q[1];
rz(-2.2687952) q[1];
sx q[1];
rz(-0.67081397) q[1];
x q[2];
rz(2.8935029) q[3];
sx q[3];
rz(-1.3458985) q[3];
sx q[3];
rz(0.096506462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.231455) q[2];
sx q[2];
rz(-2.4427876) q[2];
sx q[2];
rz(0.22496741) q[2];
rz(-0.088430017) q[3];
sx q[3];
rz(-1.6902573) q[3];
sx q[3];
rz(2.6627873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46463075) q[0];
sx q[0];
rz(-1.2372274) q[0];
sx q[0];
rz(-0.27994573) q[0];
rz(1.6784558) q[1];
sx q[1];
rz(-1.8755553) q[1];
sx q[1];
rz(2.8889012) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0137018) q[0];
sx q[0];
rz(-1.863592) q[0];
sx q[0];
rz(-3.091759) q[0];
rz(-pi) q[1];
rz(1.5964609) q[2];
sx q[2];
rz(-0.63458323) q[2];
sx q[2];
rz(0.54021013) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5340427) q[1];
sx q[1];
rz(-2.0760963) q[1];
sx q[1];
rz(0.74313785) q[1];
rz(0.75848363) q[3];
sx q[3];
rz(-2.1355724) q[3];
sx q[3];
rz(1.5837216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8213886) q[2];
sx q[2];
rz(-2.2479222) q[2];
sx q[2];
rz(1.5318711) q[2];
rz(-1.1931217) q[3];
sx q[3];
rz(-0.90819287) q[3];
sx q[3];
rz(-2.0549205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0367592) q[0];
sx q[0];
rz(-1.3830673) q[0];
sx q[0];
rz(1.6554792) q[0];
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
rz(1.7689965) q[0];
sx q[0];
rz(-2.0645752) q[0];
sx q[0];
rz(-2.2934224) q[0];
rz(-pi) q[1];
rz(2.9697044) q[2];
sx q[2];
rz(-2.066678) q[2];
sx q[2];
rz(-1.9672729) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5178691) q[1];
sx q[1];
rz(-1.9451127) q[1];
sx q[1];
rz(0.75846292) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7989743) q[3];
sx q[3];
rz(-2.5573213) q[3];
sx q[3];
rz(-2.1669471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3387317) q[2];
sx q[2];
rz(-1.677745) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6779697) q[0];
sx q[0];
rz(-2.2080053) q[0];
sx q[0];
rz(-1.8883702) q[0];
rz(0.62250096) q[1];
sx q[1];
rz(-1.6815192) q[1];
sx q[1];
rz(1.995283) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86664591) q[0];
sx q[0];
rz(-1.5005611) q[0];
sx q[0];
rz(-2.3932891) q[0];
rz(-pi) q[1];
rz(2.4856604) q[2];
sx q[2];
rz(-0.99870517) q[2];
sx q[2];
rz(-2.7149372) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.1100626) q[1];
sx q[1];
rz(-1.6441802) q[1];
sx q[1];
rz(0.88295464) q[1];
rz(-0.54543145) q[3];
sx q[3];
rz(-2.430393) q[3];
sx q[3];
rz(-1.4382854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.15381947) q[2];
sx q[2];
rz(-1.035707) q[2];
sx q[2];
rz(-2.1257607) q[2];
rz(0.9097957) q[3];
sx q[3];
rz(-0.67009059) q[3];
sx q[3];
rz(2.773496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1290454) q[0];
sx q[0];
rz(-2.5279901) q[0];
sx q[0];
rz(-3.1273499) q[0];
rz(-2.2968538) q[1];
sx q[1];
rz(-2.1886487) q[1];
sx q[1];
rz(1.7600118) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.325557) q[0];
sx q[0];
rz(-0.14294681) q[0];
sx q[0];
rz(-1.4965034) q[0];
x q[1];
rz(1.472677) q[2];
sx q[2];
rz(-0.58110229) q[2];
sx q[2];
rz(2.2475918) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4631008) q[1];
sx q[1];
rz(-1.4271724) q[1];
sx q[1];
rz(-1.0113869) q[1];
rz(-pi) q[2];
rz(0.53346177) q[3];
sx q[3];
rz(-1.7603612) q[3];
sx q[3];
rz(-1.61434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0649197) q[2];
sx q[2];
rz(-1.939247) q[2];
sx q[2];
rz(2.1515965) q[2];
rz(-2.2475217) q[3];
sx q[3];
rz(-0.48864135) q[3];
sx q[3];
rz(0.22542424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
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
rz(1.1006533) q[2];
sx q[2];
rz(-0.91081506) q[2];
sx q[2];
rz(0.87115661) q[2];
rz(-0.99466956) q[3];
sx q[3];
rz(-0.70636402) q[3];
sx q[3];
rz(-0.93887226) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
