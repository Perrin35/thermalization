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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24613334) q[0];
sx q[0];
rz(-1.5243422) q[0];
sx q[0];
rz(1.4746656) q[0];
rz(-pi) q[1];
rz(-1.5586833) q[2];
sx q[2];
rz(-1.4450057) q[2];
sx q[2];
rz(2.8319401) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1111869) q[1];
sx q[1];
rz(-0.62796794) q[1];
sx q[1];
rz(0.18917947) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0923907) q[3];
sx q[3];
rz(-2.0432825) q[3];
sx q[3];
rz(-0.14379584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8157114) q[2];
sx q[2];
rz(-1.7408966) q[2];
sx q[2];
rz(1.4665843) q[2];
rz(0.69774929) q[3];
sx q[3];
rz(-2.0402699) q[3];
sx q[3];
rz(-2.3944323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.117347) q[0];
sx q[0];
rz(-2.0139366) q[0];
sx q[0];
rz(-1.1741937) q[0];
rz(0.17114561) q[1];
sx q[1];
rz(-2.0967963) q[1];
sx q[1];
rz(0.29719621) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8288119) q[0];
sx q[0];
rz(-3.045407) q[0];
sx q[0];
rz(-2.0975153) q[0];
rz(-pi) q[1];
rz(-2.4291123) q[2];
sx q[2];
rz(-1.9409632) q[2];
sx q[2];
rz(2.9287101) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9651523) q[1];
sx q[1];
rz(-1.0607751) q[1];
sx q[1];
rz(-0.39217197) q[1];
rz(-0.41373613) q[3];
sx q[3];
rz(-0.80647751) q[3];
sx q[3];
rz(1.2456576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.979636) q[2];
sx q[2];
rz(-1.9571783) q[2];
sx q[2];
rz(0.61398181) q[2];
rz(2.2654514) q[3];
sx q[3];
rz(-0.66771475) q[3];
sx q[3];
rz(-2.5045625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
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
rz(-2.0456332) q[0];
sx q[0];
rz(-1.2852083) q[0];
sx q[0];
rz(-0.30763787) q[0];
rz(2.3930507) q[1];
sx q[1];
rz(-2.8100439) q[1];
sx q[1];
rz(2.3017853) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52536406) q[0];
sx q[0];
rz(-1.9814241) q[0];
sx q[0];
rz(1.2603659) q[0];
rz(-pi) q[1];
x q[1];
rz(0.26250458) q[2];
sx q[2];
rz(-0.35048198) q[2];
sx q[2];
rz(1.0002491) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.055132853) q[1];
sx q[1];
rz(-1.2659402) q[1];
sx q[1];
rz(-0.47385584) q[1];
rz(-pi) q[2];
rz(-2.967756) q[3];
sx q[3];
rz(-0.94508119) q[3];
sx q[3];
rz(-0.10007773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2741189) q[2];
sx q[2];
rz(-1.3511191) q[2];
sx q[2];
rz(1.8236558) q[2];
rz(-1.9258202) q[3];
sx q[3];
rz(-0.35651818) q[3];
sx q[3];
rz(1.5095476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
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
rz(-2.6960301) q[0];
rz(0.62082779) q[1];
sx q[1];
rz(-1.2460243) q[1];
sx q[1];
rz(2.1760118) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25528204) q[0];
sx q[0];
rz(-0.83034407) q[0];
sx q[0];
rz(1.9489954) q[0];
rz(-pi) q[1];
rz(2.2992209) q[2];
sx q[2];
rz(-1.5475376) q[2];
sx q[2];
rz(-2.1881441) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.021918745) q[1];
sx q[1];
rz(-1.5389171) q[1];
sx q[1];
rz(-1.3171413) q[1];
x q[2];
rz(-0.86430092) q[3];
sx q[3];
rz(-0.14926499) q[3];
sx q[3];
rz(2.5316558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.821637) q[2];
sx q[2];
rz(-0.76239061) q[2];
sx q[2];
rz(-2.2122673) q[2];
rz(2.4980513) q[3];
sx q[3];
rz(-1.0361592) q[3];
sx q[3];
rz(1.003456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4716361) q[0];
sx q[0];
rz(-1.5820553) q[0];
sx q[0];
rz(1.2840282) q[0];
rz(-0.28981003) q[1];
sx q[1];
rz(-0.73957864) q[1];
sx q[1];
rz(-1.0481542) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7832344) q[0];
sx q[0];
rz(-2.2336707) q[0];
sx q[0];
rz(0.62028424) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.83277793) q[2];
sx q[2];
rz(-1.1156429) q[2];
sx q[2];
rz(0.27311329) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.84761274) q[1];
sx q[1];
rz(-1.8925397) q[1];
sx q[1];
rz(1.6683679) q[1];
rz(-pi) q[2];
x q[2];
rz(0.9976451) q[3];
sx q[3];
rz(-0.76612681) q[3];
sx q[3];
rz(-1.0539953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8706878) q[2];
sx q[2];
rz(-2.1257766) q[2];
sx q[2];
rz(0.90083814) q[2];
rz(1.0926931) q[3];
sx q[3];
rz(-2.1381502) q[3];
sx q[3];
rz(1.2341011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95935217) q[0];
sx q[0];
rz(-1.933796) q[0];
sx q[0];
rz(-0.59610468) q[0];
rz(-1.6456564) q[1];
sx q[1];
rz(-2.2438965) q[1];
sx q[1];
rz(-1.2449107) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0743474) q[0];
sx q[0];
rz(-1.0676358) q[0];
sx q[0];
rz(-0.36672451) q[0];
rz(0.4857829) q[2];
sx q[2];
rz(-1.2885639) q[2];
sx q[2];
rz(0.9529875) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7024755) q[1];
sx q[1];
rz(-0.87279746) q[1];
sx q[1];
rz(2.4707787) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3390433) q[3];
sx q[3];
rz(-1.8125121) q[3];
sx q[3];
rz(-1.7237323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.231455) q[2];
sx q[2];
rz(-0.69880501) q[2];
sx q[2];
rz(0.22496741) q[2];
rz(-3.0531626) q[3];
sx q[3];
rz(-1.6902573) q[3];
sx q[3];
rz(0.47880539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46463075) q[0];
sx q[0];
rz(-1.2372274) q[0];
sx q[0];
rz(-2.8616469) q[0];
rz(1.4631368) q[1];
sx q[1];
rz(-1.2660374) q[1];
sx q[1];
rz(2.8889012) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1278909) q[0];
sx q[0];
rz(-1.863592) q[0];
sx q[0];
rz(3.091759) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1227038) q[2];
sx q[2];
rz(-0.9364555) q[2];
sx q[2];
rz(2.6332476) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.44811571) q[1];
sx q[1];
rz(-2.2709393) q[1];
sx q[1];
rz(-2.4561988) q[1];
rz(-pi) q[2];
x q[2];
rz(0.75848363) q[3];
sx q[3];
rz(-1.0060203) q[3];
sx q[3];
rz(-1.5837216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8213886) q[2];
sx q[2];
rz(-2.2479222) q[2];
sx q[2];
rz(1.6097216) q[2];
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
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0367592) q[0];
sx q[0];
rz(-1.7585254) q[0];
sx q[0];
rz(-1.4861134) q[0];
rz(2.8727818) q[1];
sx q[1];
rz(-2.0116282) q[1];
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
rz(1.2645623) q[2];
sx q[2];
rz(-0.5224723) q[2];
sx q[2];
rz(-2.317121) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8262417) q[1];
sx q[1];
rz(-2.3126174) q[1];
sx q[1];
rz(0.51893236) q[1];
rz(-1.7989743) q[3];
sx q[3];
rz(-0.5842714) q[3];
sx q[3];
rz(2.1669471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.802861) q[2];
sx q[2];
rz(-1.677745) q[2];
sx q[2];
rz(1.5926682) q[2];
rz(-2.0907949) q[3];
sx q[3];
rz(-1.934634) q[3];
sx q[3];
rz(-1.9416434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46362296) q[0];
sx q[0];
rz(-0.93358731) q[0];
sx q[0];
rz(1.8883702) q[0];
rz(0.62250096) q[1];
sx q[1];
rz(-1.6815192) q[1];
sx q[1];
rz(-1.1463096) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86664591) q[0];
sx q[0];
rz(-1.6410315) q[0];
sx q[0];
rz(0.74830351) q[0];
rz(-pi) q[1];
x q[1];
rz(0.88845466) q[2];
sx q[2];
rz(-2.1092215) q[2];
sx q[2];
rz(2.3922362) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.1100626) q[1];
sx q[1];
rz(-1.6441802) q[1];
sx q[1];
rz(-0.88295464) q[1];
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
rz(pi/2) q[1];
sx q[1];
rz(-2.9877732) q[2];
sx q[2];
rz(-2.1058857) q[2];
sx q[2];
rz(-1.0158319) q[2];
rz(2.231797) q[3];
sx q[3];
rz(-0.67009059) q[3];
sx q[3];
rz(0.36809665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.0125473) q[0];
sx q[0];
rz(-0.61360252) q[0];
sx q[0];
rz(-3.1273499) q[0];
rz(-2.2968538) q[1];
sx q[1];
rz(-0.95294398) q[1];
sx q[1];
rz(-1.7600118) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.325557) q[0];
sx q[0];
rz(-0.14294681) q[0];
sx q[0];
rz(-1.6450892) q[0];
rz(-pi) q[1];
rz(0.064247473) q[2];
sx q[2];
rz(-2.148743) q[2];
sx q[2];
rz(0.77677514) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.33243079) q[1];
sx q[1];
rz(-2.5659487) q[1];
sx q[1];
rz(1.8368506) q[1];
rz(-pi) q[2];
rz(-0.53346177) q[3];
sx q[3];
rz(-1.7603612) q[3];
sx q[3];
rz(-1.5272527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(1.0513231) q[0];
sx q[0];
rz(-2.4933503) q[0];
sx q[0];
rz(-1.1052263) q[0];
rz(1.3399711) q[1];
sx q[1];
rz(-0.62146386) q[1];
sx q[1];
rz(0.38846831) q[1];
rz(2.4253035) q[2];
sx q[2];
rz(-1.9367957) q[2];
sx q[2];
rz(-0.39762485) q[2];
rz(0.94974489) q[3];
sx q[3];
rz(-1.9322149) q[3];
sx q[3];
rz(0.17292427) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
