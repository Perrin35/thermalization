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
rz(5.2360143) q[0];
sx q[0];
rz(9.4935023) q[0];
rz(1.7460495) q[1];
sx q[1];
rz(4.6739251) q[1];
sx q[1];
rz(8.2164017) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2657304) q[0];
sx q[0];
rz(-3.0348572) q[0];
sx q[0];
rz(1.1197607) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0157929) q[2];
sx q[2];
rz(-1.558779) q[2];
sx q[2];
rz(1.8819686) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.61332834) q[1];
sx q[1];
rz(-1.4600888) q[1];
sx q[1];
rz(2.5221591) q[1];
x q[2];
rz(2.6192059) q[3];
sx q[3];
rz(-1.99317) q[3];
sx q[3];
rz(1.9463604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3258813) q[2];
sx q[2];
rz(-1.4006961) q[2];
sx q[2];
rz(-1.4665843) q[2];
rz(-0.69774929) q[3];
sx q[3];
rz(-2.0402699) q[3];
sx q[3];
rz(-0.74716032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0242457) q[0];
sx q[0];
rz(-1.1276561) q[0];
sx q[0];
rz(-1.1741937) q[0];
rz(2.970447) q[1];
sx q[1];
rz(-1.0447964) q[1];
sx q[1];
rz(-2.8443964) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78272351) q[0];
sx q[0];
rz(-1.6190931) q[0];
sx q[0];
rz(1.4875828) q[0];
rz(-pi) q[1];
rz(2.6058795) q[2];
sx q[2];
rz(-2.3539054) q[2];
sx q[2];
rz(2.1802769) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9651523) q[1];
sx q[1];
rz(-1.0607751) q[1];
sx q[1];
rz(0.39217197) q[1];
rz(-pi) q[2];
rz(1.9678715) q[3];
sx q[3];
rz(-2.2928772) q[3];
sx q[3];
rz(-1.811036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.16195665) q[2];
sx q[2];
rz(-1.1844144) q[2];
sx q[2];
rz(-0.61398181) q[2];
rz(-2.2654514) q[3];
sx q[3];
rz(-0.66771475) q[3];
sx q[3];
rz(2.5045625) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0456332) q[0];
sx q[0];
rz(-1.2852083) q[0];
sx q[0];
rz(-2.8339548) q[0];
rz(2.3930507) q[1];
sx q[1];
rz(-0.33154878) q[1];
sx q[1];
rz(0.83980733) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6162286) q[0];
sx q[0];
rz(-1.1601686) q[0];
sx q[0];
rz(1.8812268) q[0];
x q[1];
rz(2.8790881) q[2];
sx q[2];
rz(-2.7911107) q[2];
sx q[2];
rz(-2.1413435) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0452022) q[1];
sx q[1];
rz(-2.5844816) q[1];
sx q[1];
rz(-0.60369173) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.17383667) q[3];
sx q[3];
rz(-2.1965115) q[3];
sx q[3];
rz(-0.10007773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.86747375) q[2];
sx q[2];
rz(-1.3511191) q[2];
sx q[2];
rz(1.3179368) q[2];
rz(-1.9258202) q[3];
sx q[3];
rz(-2.7850745) q[3];
sx q[3];
rz(1.6320451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-2.1760118) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8646116) q[0];
sx q[0];
rz(-0.81482139) q[0];
sx q[0];
rz(-2.7576202) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1104286) q[2];
sx q[2];
rz(-0.84261299) q[2];
sx q[2];
rz(0.63809168) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.021918745) q[1];
sx q[1];
rz(-1.6026755) q[1];
sx q[1];
rz(-1.3171413) q[1];
rz(-1.4569034) q[3];
sx q[3];
rz(-1.6674862) q[3];
sx q[3];
rz(-0.25988042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3199557) q[2];
sx q[2];
rz(-0.76239061) q[2];
sx q[2];
rz(2.2122673) q[2];
rz(0.6435414) q[3];
sx q[3];
rz(-1.0361592) q[3];
sx q[3];
rz(-1.003456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4716361) q[0];
sx q[0];
rz(-1.5595373) q[0];
sx q[0];
rz(1.8575645) q[0];
rz(-2.8517826) q[1];
sx q[1];
rz(-2.402014) q[1];
sx q[1];
rz(-1.0481542) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3583583) q[0];
sx q[0];
rz(-0.90792197) q[0];
sx q[0];
rz(-2.5213084) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1996908) q[2];
sx q[2];
rz(-0.84398979) q[2];
sx q[2];
rz(0.84743365) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.594164) q[1];
sx q[1];
rz(-0.33572008) q[1];
sx q[1];
rz(-0.28433849) q[1];
x q[2];
rz(-0.89094152) q[3];
sx q[3];
rz(-1.1853301) q[3];
sx q[3];
rz(-0.95213529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8706878) q[2];
sx q[2];
rz(-1.015816) q[2];
sx q[2];
rz(2.2407545) q[2];
rz(2.0488996) q[3];
sx q[3];
rz(-2.1381502) q[3];
sx q[3];
rz(1.9074915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.95935217) q[0];
sx q[0];
rz(-1.933796) q[0];
sx q[0];
rz(2.545488) q[0];
rz(1.4959363) q[1];
sx q[1];
rz(-0.89769617) q[1];
sx q[1];
rz(1.2449107) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0672452) q[0];
sx q[0];
rz(-2.0739569) q[0];
sx q[0];
rz(-0.36672451) q[0];
x q[1];
rz(2.6558098) q[2];
sx q[2];
rz(-1.8530288) q[2];
sx q[2];
rz(0.9529875) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4391172) q[1];
sx q[1];
rz(-0.87279746) q[1];
sx q[1];
rz(-0.67081397) q[1];
x q[2];
rz(0.75002589) q[3];
sx q[3];
rz(-2.8083028) q[3];
sx q[3];
rz(-2.1961574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9101377) q[2];
sx q[2];
rz(-2.4427876) q[2];
sx q[2];
rz(-2.9166252) q[2];
rz(3.0531626) q[3];
sx q[3];
rz(-1.4513353) q[3];
sx q[3];
rz(0.47880539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6769619) q[0];
sx q[0];
rz(-1.9043652) q[0];
sx q[0];
rz(2.8616469) q[0];
rz(-1.4631368) q[1];
sx q[1];
rz(-1.8755553) q[1];
sx q[1];
rz(2.8889012) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7130816) q[0];
sx q[0];
rz(-1.5230852) q[0];
sx q[0];
rz(1.8639355) q[0];
rz(-pi) q[1];
rz(0.018888868) q[2];
sx q[2];
rz(-0.9364555) q[2];
sx q[2];
rz(-0.50834507) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.60755) q[1];
sx q[1];
rz(-2.0760963) q[1];
sx q[1];
rz(-0.74313785) q[1];
rz(-pi) q[2];
x q[2];
rz(2.383109) q[3];
sx q[3];
rz(-1.0060203) q[3];
sx q[3];
rz(1.5837216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.32020405) q[2];
sx q[2];
rz(-2.2479222) q[2];
sx q[2];
rz(1.5318711) q[2];
rz(1.948471) q[3];
sx q[3];
rz(-2.2333998) q[3];
sx q[3];
rz(2.0549205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0367592) q[0];
sx q[0];
rz(-1.7585254) q[0];
sx q[0];
rz(1.4861134) q[0];
rz(0.2688109) q[1];
sx q[1];
rz(-1.1299645) q[1];
sx q[1];
rz(0.2789467) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3725961) q[0];
sx q[0];
rz(-1.0770174) q[0];
sx q[0];
rz(-2.2934224) q[0];
rz(-pi) q[1];
rz(-0.17188822) q[2];
sx q[2];
rz(-1.0749146) q[2];
sx q[2];
rz(1.9672729) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.61356269) q[1];
sx q[1];
rz(-2.2655305) q[1];
sx q[1];
rz(2.0668405) q[1];
rz(-pi) q[2];
x q[2];
rz(0.14848498) q[3];
sx q[3];
rz(-2.1380224) q[3];
sx q[3];
rz(-2.4384769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3387317) q[2];
sx q[2];
rz(-1.677745) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46362296) q[0];
sx q[0];
rz(-2.2080053) q[0];
sx q[0];
rz(1.8883702) q[0];
rz(-2.5190917) q[1];
sx q[1];
rz(-1.4600735) q[1];
sx q[1];
rz(1.1463096) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5025053) q[0];
sx q[0];
rz(-2.3168132) q[0];
sx q[0];
rz(1.6665002) q[0];
x q[1];
rz(-2.3290645) q[2];
sx q[2];
rz(-2.3000237) q[2];
sx q[2];
rz(-1.7571882) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6282564) q[1];
sx q[1];
rz(-2.4504821) q[1];
sx q[1];
rz(-1.6860794) q[1];
rz(2.5067234) q[3];
sx q[3];
rz(-1.9162617) q[3];
sx q[3];
rz(-0.29840252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.15381947) q[2];
sx q[2];
rz(-2.1058857) q[2];
sx q[2];
rz(2.1257607) q[2];
rz(0.9097957) q[3];
sx q[3];
rz(-0.67009059) q[3];
sx q[3];
rz(-0.36809665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(1.1290454) q[0];
sx q[0];
rz(-2.5279901) q[0];
sx q[0];
rz(-0.01424271) q[0];
rz(0.8447389) q[1];
sx q[1];
rz(-0.95294398) q[1];
sx q[1];
rz(1.3815809) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4603699) q[0];
sx q[0];
rz(-1.5813706) q[0];
sx q[0];
rz(1.4282385) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6689156) q[2];
sx q[2];
rz(-0.58110229) q[2];
sx q[2];
rz(-0.89400089) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6784918) q[1];
sx q[1];
rz(-1.4271724) q[1];
sx q[1];
rz(-2.1302057) q[1];
rz(-pi) q[2];
rz(1.7900449) q[3];
sx q[3];
rz(-2.0937113) q[3];
sx q[3];
rz(3.074309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0766729) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0902696) q[0];
sx q[0];
rz(-0.6482424) q[0];
sx q[0];
rz(2.0363664) q[0];
rz(-1.8016215) q[1];
sx q[1];
rz(-0.62146386) q[1];
sx q[1];
rz(0.38846831) q[1];
rz(1.1006533) q[2];
sx q[2];
rz(-0.91081506) q[2];
sx q[2];
rz(0.87115661) q[2];
rz(-2.7064825) q[3];
sx q[3];
rz(-2.1463487) q[3];
sx q[3];
rz(-1.6456732) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];