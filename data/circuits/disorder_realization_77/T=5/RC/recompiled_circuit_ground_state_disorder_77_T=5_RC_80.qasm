OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.90280688) q[0];
sx q[0];
rz(-1.1385318) q[0];
sx q[0];
rz(1.0983374) q[0];
rz(1.1605473) q[1];
sx q[1];
rz(2.0055973) q[1];
sx q[1];
rz(10.027167) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9928474) q[0];
sx q[0];
rz(-1.1365599) q[0];
sx q[0];
rz(0.57172267) q[0];
rz(-1.6260398) q[2];
sx q[2];
rz(-0.85795244) q[2];
sx q[2];
rz(-1.6875658) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.25665313) q[1];
sx q[1];
rz(-0.49044427) q[1];
sx q[1];
rz(-1.3262733) q[1];
rz(-pi) q[2];
rz(-0.037558719) q[3];
sx q[3];
rz(-1.9989487) q[3];
sx q[3];
rz(1.960609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.41523734) q[2];
sx q[2];
rz(-2.0550315) q[2];
sx q[2];
rz(1.6437257) q[2];
rz(-3.030153) q[3];
sx q[3];
rz(-1.2846416) q[3];
sx q[3];
rz(2.1845412) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9621256) q[0];
sx q[0];
rz(-0.92573708) q[0];
sx q[0];
rz(0.16394462) q[0];
rz(-2.1666849) q[1];
sx q[1];
rz(-2.4539852) q[1];
sx q[1];
rz(-1.6418537) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83368084) q[0];
sx q[0];
rz(-1.320086) q[0];
sx q[0];
rz(-0.83969076) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.80625466) q[2];
sx q[2];
rz(-2.0948176) q[2];
sx q[2];
rz(-1.1576705) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.837599) q[1];
sx q[1];
rz(-1.493206) q[1];
sx q[1];
rz(2.4291888) q[1];
rz(-pi) q[2];
rz(2.1466931) q[3];
sx q[3];
rz(-2.7585976) q[3];
sx q[3];
rz(-2.8378311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6770596) q[2];
sx q[2];
rz(-2.634282) q[2];
sx q[2];
rz(-0.23059174) q[2];
rz(0.683189) q[3];
sx q[3];
rz(-0.69470996) q[3];
sx q[3];
rz(0.75146875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36000073) q[0];
sx q[0];
rz(-1.3104023) q[0];
sx q[0];
rz(-0.42136425) q[0];
rz(1.6257809) q[1];
sx q[1];
rz(-0.49585626) q[1];
sx q[1];
rz(2.1741507) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28855265) q[0];
sx q[0];
rz(-1.3167736) q[0];
sx q[0];
rz(1.2328531) q[0];
x q[1];
rz(-1.9052299) q[2];
sx q[2];
rz(-2.6504271) q[2];
sx q[2];
rz(-0.14310357) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6104504) q[1];
sx q[1];
rz(-1.4262661) q[1];
sx q[1];
rz(1.8462314) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.5021867) q[3];
sx q[3];
rz(-0.62231502) q[3];
sx q[3];
rz(2.6414504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4030054) q[2];
sx q[2];
rz(-1.0120729) q[2];
sx q[2];
rz(-0.15815132) q[2];
rz(-2.3824597) q[3];
sx q[3];
rz(-1.4594376) q[3];
sx q[3];
rz(-2.7706743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.89449969) q[0];
sx q[0];
rz(-2.5573754) q[0];
sx q[0];
rz(-1.1868813) q[0];
rz(-0.12313708) q[1];
sx q[1];
rz(-1.5703399) q[1];
sx q[1];
rz(1.311696) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7191221) q[0];
sx q[0];
rz(-1.718169) q[0];
sx q[0];
rz(-3.0585994) q[0];
x q[1];
rz(-3.0391215) q[2];
sx q[2];
rz(-1.6471286) q[2];
sx q[2];
rz(1.1427405) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0801436) q[1];
sx q[1];
rz(-1.9654566) q[1];
sx q[1];
rz(0.45545642) q[1];
rz(1.5748686) q[3];
sx q[3];
rz(-1.5713672) q[3];
sx q[3];
rz(-0.79889311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9436403) q[2];
sx q[2];
rz(-1.8215048) q[2];
sx q[2];
rz(-0.8521592) q[2];
rz(-0.85396829) q[3];
sx q[3];
rz(-1.9210509) q[3];
sx q[3];
rz(1.0796237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-2.3581486) q[0];
sx q[0];
rz(-1.7572948) q[0];
sx q[0];
rz(-3.070991) q[0];
rz(1.0380896) q[1];
sx q[1];
rz(-1.9990653) q[1];
sx q[1];
rz(-0.80931726) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22915086) q[0];
sx q[0];
rz(-2.1437867) q[0];
sx q[0];
rz(1.6567162) q[0];
rz(2.4368183) q[2];
sx q[2];
rz(-0.88181978) q[2];
sx q[2];
rz(-2.5344684) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0496724) q[1];
sx q[1];
rz(-0.61245239) q[1];
sx q[1];
rz(-2.6785707) q[1];
rz(-1.3594987) q[3];
sx q[3];
rz(-0.70970067) q[3];
sx q[3];
rz(3.078408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.91850963) q[2];
sx q[2];
rz(-0.56879908) q[2];
sx q[2];
rz(-1.1893547) q[2];
rz(0.14063028) q[3];
sx q[3];
rz(-0.46969241) q[3];
sx q[3];
rz(2.752059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4024432) q[0];
sx q[0];
rz(-2.2279255) q[0];
sx q[0];
rz(-1.5341349) q[0];
rz(0.88273478) q[1];
sx q[1];
rz(-0.84461707) q[1];
sx q[1];
rz(-0.61558634) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7587335) q[0];
sx q[0];
rz(-2.737597) q[0];
sx q[0];
rz(-0.23601471) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.83258038) q[2];
sx q[2];
rz(-2.2924149) q[2];
sx q[2];
rz(0.20038062) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2193377) q[1];
sx q[1];
rz(-0.23592792) q[1];
sx q[1];
rz(0.94249819) q[1];
rz(1.5546404) q[3];
sx q[3];
rz(-0.62045853) q[3];
sx q[3];
rz(-0.66974332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9756056) q[2];
sx q[2];
rz(-0.96419445) q[2];
sx q[2];
rz(-1.9624422) q[2];
rz(-0.21373448) q[3];
sx q[3];
rz(-1.7515747) q[3];
sx q[3];
rz(-2.8960622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24446503) q[0];
sx q[0];
rz(-1.1470969) q[0];
sx q[0];
rz(0.044064673) q[0];
rz(-2.7524718) q[1];
sx q[1];
rz(-2.6521284) q[1];
sx q[1];
rz(0.74180952) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76352966) q[0];
sx q[0];
rz(-0.46398315) q[0];
sx q[0];
rz(-0.65391175) q[0];
x q[1];
rz(2.7905383) q[2];
sx q[2];
rz(-0.85509713) q[2];
sx q[2];
rz(2.9911016) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.68315893) q[1];
sx q[1];
rz(-1.6583879) q[1];
sx q[1];
rz(0.32282543) q[1];
rz(1.7735305) q[3];
sx q[3];
rz(-2.2623767) q[3];
sx q[3];
rz(1.4940624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.47702181) q[2];
sx q[2];
rz(-1.7844113) q[2];
sx q[2];
rz(2.6396497) q[2];
rz(1.4255514) q[3];
sx q[3];
rz(-2.612096) q[3];
sx q[3];
rz(0.79290032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75941706) q[0];
sx q[0];
rz(-1.4074396) q[0];
sx q[0];
rz(2.7446246) q[0];
rz(-1.5396897) q[1];
sx q[1];
rz(-2.868728) q[1];
sx q[1];
rz(2.4698965) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26123372) q[0];
sx q[0];
rz(-1.6971998) q[0];
sx q[0];
rz(-2.0194299) q[0];
rz(-pi) q[1];
rz(-2.0136859) q[2];
sx q[2];
rz(-1.9077565) q[2];
sx q[2];
rz(-2.2786587) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4038329) q[1];
sx q[1];
rz(-1.8813526) q[1];
sx q[1];
rz(-2.4871189) q[1];
rz(-0.83064744) q[3];
sx q[3];
rz(-1.1201356) q[3];
sx q[3];
rz(-0.78307952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4063065) q[2];
sx q[2];
rz(-2.6267509) q[2];
sx q[2];
rz(-1.7740645) q[2];
rz(-0.96274084) q[3];
sx q[3];
rz(-1.5509501) q[3];
sx q[3];
rz(2.4875557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3108567) q[0];
sx q[0];
rz(-1.5850569) q[0];
sx q[0];
rz(2.8209525) q[0];
rz(0.93797183) q[1];
sx q[1];
rz(-1.9357977) q[1];
sx q[1];
rz(1.4042312) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7119634) q[0];
sx q[0];
rz(-2.9949463) q[0];
sx q[0];
rz(2.3578927) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5455568) q[2];
sx q[2];
rz(-1.847731) q[2];
sx q[2];
rz(2.0933088) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2102846) q[1];
sx q[1];
rz(-0.043593229) q[1];
sx q[1];
rz(2.0491854) q[1];
rz(-0.76227607) q[3];
sx q[3];
rz(-1.395985) q[3];
sx q[3];
rz(-2.0559514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2001026) q[2];
sx q[2];
rz(-0.76277554) q[2];
sx q[2];
rz(2.4007559) q[2];
rz(-1.0624933) q[3];
sx q[3];
rz(-1.1081568) q[3];
sx q[3];
rz(-1.9443996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2049388) q[0];
sx q[0];
rz(-2.7118201) q[0];
sx q[0];
rz(-2.3858261) q[0];
rz(1.3142122) q[1];
sx q[1];
rz(-1.5145489) q[1];
sx q[1];
rz(-2.4155713) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7814815) q[0];
sx q[0];
rz(-1.5364293) q[0];
sx q[0];
rz(-0.78256677) q[0];
rz(-pi) q[1];
rz(3.0376932) q[2];
sx q[2];
rz(-2.2533721) q[2];
sx q[2];
rz(2.7524591) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0997114) q[1];
sx q[1];
rz(-1.9181632) q[1];
sx q[1];
rz(-0.16176407) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0199225) q[3];
sx q[3];
rz(-1.4360997) q[3];
sx q[3];
rz(-1.0659892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3296335) q[2];
sx q[2];
rz(-0.52874955) q[2];
sx q[2];
rz(2.192396) q[2];
rz(2.9694929) q[3];
sx q[3];
rz(-2.1642919) q[3];
sx q[3];
rz(-2.7026091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
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
rz(3.1127472) q[0];
sx q[0];
rz(-1.5486568) q[0];
sx q[0];
rz(-1.5966709) q[0];
rz(-1.9817837) q[1];
sx q[1];
rz(-1.8198967) q[1];
sx q[1];
rz(1.5130704) q[1];
rz(-1.415435) q[2];
sx q[2];
rz(-1.0280357) q[2];
sx q[2];
rz(0.24161777) q[2];
rz(0.070318182) q[3];
sx q[3];
rz(-2.0587772) q[3];
sx q[3];
rz(-1.8785431) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
