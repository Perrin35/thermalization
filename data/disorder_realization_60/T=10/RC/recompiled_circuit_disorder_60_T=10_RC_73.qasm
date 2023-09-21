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
rz(1.7460495) q[1];
sx q[1];
rz(-1.6092602) q[1];
sx q[1];
rz(-1.2083763) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3291408) q[0];
sx q[0];
rz(-1.4747696) q[0];
sx q[0];
rz(-0.046669331) q[0];
x q[1];
rz(3.0460998) q[2];
sx q[2];
rz(-0.12636939) q[2];
sx q[2];
rz(-0.40590826) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1111869) q[1];
sx q[1];
rz(-0.62796794) q[1];
sx q[1];
rz(0.18917947) q[1];
rz(-0.52238676) q[3];
sx q[3];
rz(-1.99317) q[3];
sx q[3];
rz(-1.1952323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3258813) q[2];
sx q[2];
rz(-1.7408966) q[2];
sx q[2];
rz(1.4665843) q[2];
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
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0242457) q[0];
sx q[0];
rz(-1.1276561) q[0];
sx q[0];
rz(1.1741937) q[0];
rz(0.17114561) q[1];
sx q[1];
rz(-2.0967963) q[1];
sx q[1];
rz(-2.8443964) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78272351) q[0];
sx q[0];
rz(-1.6190931) q[0];
sx q[0];
rz(1.6540098) q[0];
rz(-pi) q[1];
rz(-0.53571312) q[2];
sx q[2];
rz(-2.3539054) q[2];
sx q[2];
rz(-0.96131575) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.17644037) q[1];
sx q[1];
rz(-2.0808176) q[1];
sx q[1];
rz(2.7494207) q[1];
x q[2];
rz(1.1737212) q[3];
sx q[3];
rz(-2.2928772) q[3];
sx q[3];
rz(-1.3305566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.16195665) q[2];
sx q[2];
rz(-1.1844144) q[2];
sx q[2];
rz(-0.61398181) q[2];
rz(-0.87614122) q[3];
sx q[3];
rz(-0.66771475) q[3];
sx q[3];
rz(-2.5045625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0456332) q[0];
sx q[0];
rz(-1.8563844) q[0];
sx q[0];
rz(0.30763787) q[0];
rz(-0.74854198) q[1];
sx q[1];
rz(-0.33154878) q[1];
sx q[1];
rz(-2.3017853) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15159431) q[0];
sx q[0];
rz(-0.50938207) q[0];
sx q[0];
rz(2.529782) q[0];
x q[1];
rz(0.26250458) q[2];
sx q[2];
rz(-0.35048198) q[2];
sx q[2];
rz(-2.1413435) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0963904) q[1];
sx q[1];
rz(-2.5844816) q[1];
sx q[1];
rz(-0.60369173) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8057459) q[3];
sx q[3];
rz(-0.64628212) q[3];
sx q[3];
rz(0.39138734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2741189) q[2];
sx q[2];
rz(-1.3511191) q[2];
sx q[2];
rz(-1.3179368) q[2];
rz(1.9258202) q[3];
sx q[3];
rz(-2.7850745) q[3];
sx q[3];
rz(-1.6320451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3777305) q[0];
sx q[0];
rz(-1.7611935) q[0];
sx q[0];
rz(2.6960301) q[0];
rz(0.62082779) q[1];
sx q[1];
rz(-1.2460243) q[1];
sx q[1];
rz(2.1760118) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25528204) q[0];
sx q[0];
rz(-2.3112486) q[0];
sx q[0];
rz(1.9489954) q[0];
rz(-0.031164073) q[2];
sx q[2];
rz(-2.2989797) q[2];
sx q[2];
rz(0.63809168) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.671215) q[1];
sx q[1];
rz(-0.2556076) q[1];
sx q[1];
rz(-1.6971991) q[1];
x q[2];
rz(-0.86430092) q[3];
sx q[3];
rz(-0.14926499) q[3];
sx q[3];
rz(-0.60993689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3199557) q[2];
sx q[2];
rz(-0.76239061) q[2];
sx q[2];
rz(2.2122673) q[2];
rz(-0.6435414) q[3];
sx q[3];
rz(-1.0361592) q[3];
sx q[3];
rz(-2.1381366) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
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
rz(-2.402014) q[1];
sx q[1];
rz(-1.0481542) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50131932) q[0];
sx q[0];
rz(-2.2671056) q[0];
sx q[0];
rz(0.93080824) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1996908) q[2];
sx q[2];
rz(-0.84398979) q[2];
sx q[2];
rz(0.84743365) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.75412616) q[1];
sx q[1];
rz(-1.6633463) q[1];
sx q[1];
rz(-2.8184163) q[1];
rz(-pi) q[2];
rz(0.89094152) q[3];
sx q[3];
rz(-1.9562625) q[3];
sx q[3];
rz(-0.95213529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8706878) q[2];
sx q[2];
rz(-2.1257766) q[2];
sx q[2];
rz(-0.90083814) q[2];
rz(2.0488996) q[3];
sx q[3];
rz(-1.0034424) q[3];
sx q[3];
rz(-1.9074915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95935217) q[0];
sx q[0];
rz(-1.933796) q[0];
sx q[0];
rz(-0.59610468) q[0];
rz(-1.4959363) q[1];
sx q[1];
rz(-0.89769617) q[1];
sx q[1];
rz(-1.2449107) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7398864) q[0];
sx q[0];
rz(-0.6131999) q[0];
sx q[0];
rz(-0.99341157) q[0];
x q[1];
rz(1.2539358) q[2];
sx q[2];
rz(-1.1057901) q[2];
sx q[2];
rz(-0.76380619) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.60331261) q[1];
sx q[1];
rz(-2.0671751) q[1];
sx q[1];
rz(2.3904295) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3390433) q[3];
sx q[3];
rz(-1.3290805) q[3];
sx q[3];
rz(-1.7237323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9101377) q[2];
sx q[2];
rz(-0.69880501) q[2];
sx q[2];
rz(0.22496741) q[2];
rz(0.088430017) q[3];
sx q[3];
rz(-1.4513353) q[3];
sx q[3];
rz(-0.47880539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46463075) q[0];
sx q[0];
rz(-1.9043652) q[0];
sx q[0];
rz(-2.8616469) q[0];
rz(-1.4631368) q[1];
sx q[1];
rz(-1.8755553) q[1];
sx q[1];
rz(-0.25269145) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1278909) q[0];
sx q[0];
rz(-1.2780006) q[0];
sx q[0];
rz(-3.091759) q[0];
rz(-pi) q[1];
rz(-1.5964609) q[2];
sx q[2];
rz(-2.5070094) q[2];
sx q[2];
rz(0.54021013) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.686324) q[1];
sx q[1];
rz(-0.9372006) q[1];
sx q[1];
rz(-0.92647657) q[1];
x q[2];
rz(-0.85315506) q[3];
sx q[3];
rz(-0.95082885) q[3];
sx q[3];
rz(-2.6593047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8213886) q[2];
sx q[2];
rz(-2.2479222) q[2];
sx q[2];
rz(1.6097216) q[2];
rz(1.948471) q[3];
sx q[3];
rz(-2.2333998) q[3];
sx q[3];
rz(2.0549205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0367592) q[0];
sx q[0];
rz(-1.3830673) q[0];
sx q[0];
rz(1.6554792) q[0];
rz(-0.2688109) q[1];
sx q[1];
rz(-1.1299645) q[1];
sx q[1];
rz(-0.2789467) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3725961) q[0];
sx q[0];
rz(-2.0645752) q[0];
sx q[0];
rz(-0.84817024) q[0];
rz(-2.9697044) q[2];
sx q[2];
rz(-2.066678) q[2];
sx q[2];
rz(-1.1743197) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.61356269) q[1];
sx q[1];
rz(-2.2655305) q[1];
sx q[1];
rz(-2.0668405) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.143048) q[3];
sx q[3];
rz(-1.6958941) q[3];
sx q[3];
rz(-0.78748122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.802861) q[2];
sx q[2];
rz(-1.677745) q[2];
sx q[2];
rz(-1.5926682) q[2];
rz(-1.0507978) q[3];
sx q[3];
rz(-1.934634) q[3];
sx q[3];
rz(1.9416434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46362296) q[0];
sx q[0];
rz(-0.93358731) q[0];
sx q[0];
rz(1.2532225) q[0];
rz(-2.5190917) q[1];
sx q[1];
rz(-1.4600735) q[1];
sx q[1];
rz(1.1463096) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.779594) q[0];
sx q[0];
rz(-2.3906374) q[0];
sx q[0];
rz(0.10303084) q[0];
x q[1];
rz(2.4856604) q[2];
sx q[2];
rz(-0.99870517) q[2];
sx q[2];
rz(-2.7149372) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6282564) q[1];
sx q[1];
rz(-0.69111052) q[1];
sx q[1];
rz(-1.4555132) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9911489) q[3];
sx q[3];
rz(-0.97878362) q[3];
sx q[3];
rz(2.1136485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.15381947) q[2];
sx q[2];
rz(-2.1058857) q[2];
sx q[2];
rz(2.1257607) q[2];
rz(-0.9097957) q[3];
sx q[3];
rz(-2.4715021) q[3];
sx q[3];
rz(2.773496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
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
rz(-0.8447389) q[1];
sx q[1];
rz(-2.1886487) q[1];
sx q[1];
rz(-1.7600118) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2505014) q[0];
sx q[0];
rz(-1.7133461) q[0];
sx q[0];
rz(-3.13091) q[0];
rz(0.99190418) q[2];
sx q[2];
rz(-1.5169946) q[2];
sx q[2];
rz(0.75888854) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6784918) q[1];
sx q[1];
rz(-1.4271724) q[1];
sx q[1];
rz(-1.0113869) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7900449) q[3];
sx q[3];
rz(-1.0478813) q[3];
sx q[3];
rz(-0.067283665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0649197) q[2];
sx q[2];
rz(-1.939247) q[2];
sx q[2];
rz(-0.98999611) q[2];
rz(2.2475217) q[3];
sx q[3];
rz(-2.6529513) q[3];
sx q[3];
rz(0.22542424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-1.0513231) q[0];
sx q[0];
rz(-0.6482424) q[0];
sx q[0];
rz(2.0363664) q[0];
rz(-1.8016215) q[1];
sx q[1];
rz(-0.62146386) q[1];
sx q[1];
rz(0.38846831) q[1];
rz(2.6132265) q[2];
sx q[2];
rz(-2.3522204) q[2];
sx q[2];
rz(1.5632202) q[2];
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
