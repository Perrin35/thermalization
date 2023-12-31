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
rz(-1.3955431) q[1];
sx q[1];
rz(-1.5323324) q[1];
sx q[1];
rz(1.2083763) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3291408) q[0];
sx q[0];
rz(-1.4747696) q[0];
sx q[0];
rz(-3.0949233) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0460998) q[2];
sx q[2];
rz(-3.0152233) q[2];
sx q[2];
rz(0.40590826) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0304058) q[1];
sx q[1];
rz(-2.5136247) q[1];
sx q[1];
rz(-0.18917947) q[1];
x q[2];
rz(-2.4083706) q[3];
sx q[3];
rz(-2.4823722) q[3];
sx q[3];
rz(-0.99430195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3258813) q[2];
sx q[2];
rz(-1.4006961) q[2];
sx q[2];
rz(1.4665843) q[2];
rz(-2.4438434) q[3];
sx q[3];
rz(-1.1013228) q[3];
sx q[3];
rz(-0.74716032) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0242457) q[0];
sx q[0];
rz(-1.1276561) q[0];
sx q[0];
rz(-1.1741937) q[0];
rz(-2.970447) q[1];
sx q[1];
rz(-1.0447964) q[1];
sx q[1];
rz(-0.29719621) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3588691) q[0];
sx q[0];
rz(-1.6190931) q[0];
sx q[0];
rz(-1.4875828) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6058795) q[2];
sx q[2];
rz(-0.78768724) q[2];
sx q[2];
rz(-2.1802769) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.17644037) q[1];
sx q[1];
rz(-1.0607751) q[1];
sx q[1];
rz(-0.39217197) q[1];
rz(-pi) q[2];
rz(2.3791749) q[3];
sx q[3];
rz(-1.8652417) q[3];
sx q[3];
rz(0.03014119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.979636) q[2];
sx q[2];
rz(-1.1844144) q[2];
sx q[2];
rz(-2.5276108) q[2];
rz(0.87614122) q[3];
sx q[3];
rz(-2.4738779) q[3];
sx q[3];
rz(-2.5045625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0456332) q[0];
sx q[0];
rz(-1.8563844) q[0];
sx q[0];
rz(-0.30763787) q[0];
rz(-2.3930507) q[1];
sx q[1];
rz(-0.33154878) q[1];
sx q[1];
rz(-0.83980733) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15159431) q[0];
sx q[0];
rz(-2.6322106) q[0];
sx q[0];
rz(-2.529782) q[0];
rz(-pi) q[1];
rz(1.4762127) q[2];
sx q[2];
rz(-1.2328096) q[2];
sx q[2];
rz(-2.4199977) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0452022) q[1];
sx q[1];
rz(-0.55711105) q[1];
sx q[1];
rz(-2.5379009) q[1];
rz(-pi) q[2];
rz(-1.8057459) q[3];
sx q[3];
rz(-0.64628212) q[3];
sx q[3];
rz(2.7502053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2741189) q[2];
sx q[2];
rz(-1.7904736) q[2];
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
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3777305) q[0];
sx q[0];
rz(-1.3803991) q[0];
sx q[0];
rz(2.6960301) q[0];
rz(2.5207649) q[1];
sx q[1];
rz(-1.8955684) q[1];
sx q[1];
rz(-0.96558085) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25528204) q[0];
sx q[0];
rz(-0.83034407) q[0];
sx q[0];
rz(-1.9489954) q[0];
x q[1];
rz(-3.1104286) q[2];
sx q[2];
rz(-2.2989797) q[2];
sx q[2];
rz(-0.63809168) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5406148) q[1];
sx q[1];
rz(-1.317273) q[1];
sx q[1];
rz(-3.1086604) q[1];
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
rz(-1.821637) q[2];
sx q[2];
rz(-2.379202) q[2];
sx q[2];
rz(-2.2122673) q[2];
rz(0.6435414) q[3];
sx q[3];
rz(-1.0361592) q[3];
sx q[3];
rz(2.1381366) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4716361) q[0];
sx q[0];
rz(-1.5820553) q[0];
sx q[0];
rz(-1.8575645) q[0];
rz(2.8517826) q[1];
sx q[1];
rz(-2.402014) q[1];
sx q[1];
rz(-2.0934385) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50131932) q[0];
sx q[0];
rz(-2.2671056) q[0];
sx q[0];
rz(2.2107844) q[0];
rz(2.5571312) q[2];
sx q[2];
rz(-2.2197154) q[2];
sx q[2];
rz(1.6780168) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3874665) q[1];
sx q[1];
rz(-1.4782463) q[1];
sx q[1];
rz(2.8184163) q[1];
rz(2.1439476) q[3];
sx q[3];
rz(-2.3754658) q[3];
sx q[3];
rz(-1.0539953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8706878) q[2];
sx q[2];
rz(-2.1257766) q[2];
sx q[2];
rz(-2.2407545) q[2];
rz(1.0926931) q[3];
sx q[3];
rz(-1.0034424) q[3];
sx q[3];
rz(-1.2341011) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1822405) q[0];
sx q[0];
rz(-1.2077967) q[0];
sx q[0];
rz(0.59610468) q[0];
rz(1.6456564) q[1];
sx q[1];
rz(-2.2438965) q[1];
sx q[1];
rz(1.2449107) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67958528) q[0];
sx q[0];
rz(-1.2512659) q[0];
sx q[0];
rz(-1.0380448) q[0];
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
rz(-2.5925496) q[1];
sx q[1];
rz(-0.92714308) q[1];
sx q[1];
rz(-2.2085269) q[1];
rz(-pi) q[2];
x q[2];
rz(0.75002589) q[3];
sx q[3];
rz(-2.8083028) q[3];
sx q[3];
rz(0.94543524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9101377) q[2];
sx q[2];
rz(-2.4427876) q[2];
sx q[2];
rz(0.22496741) q[2];
rz(3.0531626) q[3];
sx q[3];
rz(-1.6902573) q[3];
sx q[3];
rz(-0.47880539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6769619) q[0];
sx q[0];
rz(-1.2372274) q[0];
sx q[0];
rz(0.27994573) q[0];
rz(-1.6784558) q[1];
sx q[1];
rz(-1.8755553) q[1];
sx q[1];
rz(0.25269145) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1278909) q[0];
sx q[0];
rz(-1.863592) q[0];
sx q[0];
rz(0.049833628) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5451317) q[2];
sx q[2];
rz(-0.63458323) q[2];
sx q[2];
rz(0.54021013) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5340427) q[1];
sx q[1];
rz(-1.0654963) q[1];
sx q[1];
rz(2.3984548) q[1];
rz(0.74440646) q[3];
sx q[3];
rz(-0.91074569) q[3];
sx q[3];
rz(-2.6403514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8213886) q[2];
sx q[2];
rz(-2.2479222) q[2];
sx q[2];
rz(-1.5318711) q[2];
rz(-1.948471) q[3];
sx q[3];
rz(-0.90819287) q[3];
sx q[3];
rz(-1.0866722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10483345) q[0];
sx q[0];
rz(-1.3830673) q[0];
sx q[0];
rz(1.4861134) q[0];
rz(2.8727818) q[1];
sx q[1];
rz(-2.0116282) q[1];
sx q[1];
rz(0.2789467) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19764087) q[0];
sx q[0];
rz(-0.949172) q[0];
sx q[0];
rz(0.62244121) q[0];
x q[1];
rz(2.0729162) q[2];
sx q[2];
rz(-1.7218105) q[2];
sx q[2];
rz(2.6627024) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.61356269) q[1];
sx q[1];
rz(-0.87606214) q[1];
sx q[1];
rz(1.0747521) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7989743) q[3];
sx q[3];
rz(-0.5842714) q[3];
sx q[3];
rz(-2.1669471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3387317) q[2];
sx q[2];
rz(-1.677745) q[2];
sx q[2];
rz(1.5926682) q[2];
rz(-1.0507978) q[3];
sx q[3];
rz(-1.934634) q[3];
sx q[3];
rz(-1.1999493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46362296) q[0];
sx q[0];
rz(-2.2080053) q[0];
sx q[0];
rz(-1.8883702) q[0];
rz(0.62250096) q[1];
sx q[1];
rz(-1.6815192) q[1];
sx q[1];
rz(1.995283) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.779594) q[0];
sx q[0];
rz(-0.75095526) q[0];
sx q[0];
rz(0.10303084) q[0];
x q[1];
rz(2.4856604) q[2];
sx q[2];
rz(-0.99870517) q[2];
sx q[2];
rz(-2.7149372) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.031530023) q[1];
sx q[1];
rz(-1.4974125) q[1];
sx q[1];
rz(-2.258638) q[1];
rz(-0.54543145) q[3];
sx q[3];
rz(-0.71119961) q[3];
sx q[3];
rz(-1.7033073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9877732) q[2];
sx q[2];
rz(-2.1058857) q[2];
sx q[2];
rz(1.0158319) q[2];
rz(2.231797) q[3];
sx q[3];
rz(-0.67009059) q[3];
sx q[3];
rz(0.36809665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0125473) q[0];
sx q[0];
rz(-2.5279901) q[0];
sx q[0];
rz(-3.1273499) q[0];
rz(2.2968538) q[1];
sx q[1];
rz(-2.1886487) q[1];
sx q[1];
rz(1.3815809) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81603564) q[0];
sx q[0];
rz(-0.14294681) q[0];
sx q[0];
rz(1.4965034) q[0];
x q[1];
rz(-2.1496885) q[2];
sx q[2];
rz(-1.5169946) q[2];
sx q[2];
rz(0.75888854) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.018316293) q[1];
sx q[1];
rz(-1.0178207) q[1];
sx q[1];
rz(0.1690013) q[1];
x q[2];
rz(1.7900449) q[3];
sx q[3];
rz(-1.0478813) q[3];
sx q[3];
rz(0.067283665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0649197) q[2];
sx q[2];
rz(-1.2023456) q[2];
sx q[2];
rz(-0.98999611) q[2];
rz(-2.2475217) q[3];
sx q[3];
rz(-2.6529513) q[3];
sx q[3];
rz(2.9161684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0513231) q[0];
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
rz(-2.1469231) q[3];
sx q[3];
rz(-2.4352286) q[3];
sx q[3];
rz(2.2027204) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
