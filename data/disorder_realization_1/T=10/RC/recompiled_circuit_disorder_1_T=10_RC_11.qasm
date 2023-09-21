OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.82436615) q[0];
sx q[0];
rz(-1.1146201) q[0];
sx q[0];
rz(-0.00014076509) q[0];
rz(1.3340985) q[1];
sx q[1];
rz(-2.1773832) q[1];
sx q[1];
rz(1.1934086) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0573187) q[0];
sx q[0];
rz(-2.8017375) q[0];
sx q[0];
rz(2.1291332) q[0];
rz(-pi) q[1];
rz(0.5483746) q[2];
sx q[2];
rz(-1.3142685) q[2];
sx q[2];
rz(-2.246849) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5692917) q[1];
sx q[1];
rz(-2.3094059) q[1];
sx q[1];
rz(0.65170793) q[1];
rz(-pi) q[2];
rz(0.10981202) q[3];
sx q[3];
rz(-1.3545274) q[3];
sx q[3];
rz(0.11085489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6821735) q[2];
sx q[2];
rz(-0.023962263) q[2];
sx q[2];
rz(1.2288644) q[2];
rz(-1.4131644) q[3];
sx q[3];
rz(-1.1011522) q[3];
sx q[3];
rz(1.6536973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5380149) q[0];
sx q[0];
rz(-1.6390272) q[0];
sx q[0];
rz(1.0128101) q[0];
rz(-3.1139328) q[1];
sx q[1];
rz(-0.67359567) q[1];
sx q[1];
rz(-1.123463) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7682122) q[0];
sx q[0];
rz(-1.1102311) q[0];
sx q[0];
rz(3.0773613) q[0];
x q[1];
rz(0.77983071) q[2];
sx q[2];
rz(-1.5260328) q[2];
sx q[2];
rz(1.1080527) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9065735) q[1];
sx q[1];
rz(-2.0779013) q[1];
sx q[1];
rz(-0.97096918) q[1];
rz(-pi) q[2];
rz(0.68617679) q[3];
sx q[3];
rz(-2.3585329) q[3];
sx q[3];
rz(0.99883294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3479487) q[2];
sx q[2];
rz(-2.0517893) q[2];
sx q[2];
rz(2.222555) q[2];
rz(0.67409003) q[3];
sx q[3];
rz(-0.6522817) q[3];
sx q[3];
rz(1.6154217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8640901) q[0];
sx q[0];
rz(-2.9798177) q[0];
sx q[0];
rz(-1.2751689) q[0];
rz(-2.4480942) q[1];
sx q[1];
rz(-1.8854515) q[1];
sx q[1];
rz(-1.1330053) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1904859) q[0];
sx q[0];
rz(-1.5740984) q[0];
sx q[0];
rz(-1.4285054) q[0];
rz(-pi) q[1];
rz(-2.6373367) q[2];
sx q[2];
rz(-2.2849053) q[2];
sx q[2];
rz(-2.2222663) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9085711) q[1];
sx q[1];
rz(-2.1070707) q[1];
sx q[1];
rz(-2.8051393) q[1];
x q[2];
rz(2.1902309) q[3];
sx q[3];
rz(-2.1054483) q[3];
sx q[3];
rz(1.9922647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.2514078) q[2];
sx q[2];
rz(-0.79139411) q[2];
sx q[2];
rz(-1.8481002) q[2];
rz(-3.1022762) q[3];
sx q[3];
rz(-1.2189564) q[3];
sx q[3];
rz(-1.8815276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2599729) q[0];
sx q[0];
rz(-3.0631174) q[0];
sx q[0];
rz(-1.1608634) q[0];
rz(-0.89598957) q[1];
sx q[1];
rz(-1.7005824) q[1];
sx q[1];
rz(-3.0060351) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3121658) q[0];
sx q[0];
rz(-0.78071785) q[0];
sx q[0];
rz(-2.6902945) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9549471) q[2];
sx q[2];
rz(-1.640056) q[2];
sx q[2];
rz(0.70089507) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.23952661) q[1];
sx q[1];
rz(-0.24589989) q[1];
sx q[1];
rz(-0.86218254) q[1];
rz(-pi) q[2];
rz(2.0849864) q[3];
sx q[3];
rz(-2.8286472) q[3];
sx q[3];
rz(1.3514951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9049412) q[2];
sx q[2];
rz(-2.1950978) q[2];
sx q[2];
rz(-0.87990749) q[2];
rz(3.0974292) q[3];
sx q[3];
rz(-1.5019838) q[3];
sx q[3];
rz(2.8529609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1039466) q[0];
sx q[0];
rz(-0.3750616) q[0];
sx q[0];
rz(1.0132382) q[0];
rz(3.0918616) q[1];
sx q[1];
rz(-2.2278992) q[1];
sx q[1];
rz(1.0838881) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5686544) q[0];
sx q[0];
rz(-1.3875811) q[0];
sx q[0];
rz(1.818548) q[0];
x q[1];
rz(-2.3124144) q[2];
sx q[2];
rz(-0.36311705) q[2];
sx q[2];
rz(-1.6104289) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0368082) q[1];
sx q[1];
rz(-1.0462865) q[1];
sx q[1];
rz(0.12496897) q[1];
rz(-pi) q[2];
rz(3.0524592) q[3];
sx q[3];
rz(-2.6260178) q[3];
sx q[3];
rz(0.46407947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.23285398) q[2];
sx q[2];
rz(-2.8149657) q[2];
sx q[2];
rz(-2.8971635) q[2];
rz(-2.7092253) q[3];
sx q[3];
rz(-1.3997388) q[3];
sx q[3];
rz(-0.50306815) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8571092) q[0];
sx q[0];
rz(-1.720022) q[0];
sx q[0];
rz(3.0474512) q[0];
rz(-0.17177467) q[1];
sx q[1];
rz(-2.005902) q[1];
sx q[1];
rz(0.89541268) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58127922) q[0];
sx q[0];
rz(-0.34438294) q[0];
sx q[0];
rz(3.0292105) q[0];
rz(-1.8108098) q[2];
sx q[2];
rz(-0.95108205) q[2];
sx q[2];
rz(-1.0735219) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6678644) q[1];
sx q[1];
rz(-1.6147991) q[1];
sx q[1];
rz(1.4021137) q[1];
rz(-1.496109) q[3];
sx q[3];
rz(-1.5731305) q[3];
sx q[3];
rz(0.96935779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.133693) q[2];
sx q[2];
rz(-2.7329625) q[2];
sx q[2];
rz(2.3383979) q[2];
rz(-1.9512272) q[3];
sx q[3];
rz(-1.2322216) q[3];
sx q[3];
rz(2.7289594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.068709277) q[0];
sx q[0];
rz(-2.9769653) q[0];
sx q[0];
rz(0.51914006) q[0];
rz(-0.58147645) q[1];
sx q[1];
rz(-2.0362208) q[1];
sx q[1];
rz(-1.2566459) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.089766895) q[0];
sx q[0];
rz(-1.8221812) q[0];
sx q[0];
rz(-3.0999523) q[0];
rz(-pi) q[1];
rz(-0.2101375) q[2];
sx q[2];
rz(-2.0268831) q[2];
sx q[2];
rz(1.8535341) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9659121) q[1];
sx q[1];
rz(-1.5064081) q[1];
sx q[1];
rz(2.8596911) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.7092181) q[3];
sx q[3];
rz(-1.3243444) q[3];
sx q[3];
rz(-2.646288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1533623) q[2];
sx q[2];
rz(-1.0299269) q[2];
sx q[2];
rz(-1.777565) q[2];
rz(0.91056943) q[3];
sx q[3];
rz(-1.986859) q[3];
sx q[3];
rz(1.5301269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0664739) q[0];
sx q[0];
rz(-2.5771038) q[0];
sx q[0];
rz(2.8334154) q[0];
rz(-0.072487436) q[1];
sx q[1];
rz(-2.1283573) q[1];
sx q[1];
rz(2.7546308) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21803074) q[0];
sx q[0];
rz(-2.1523928) q[0];
sx q[0];
rz(-2.8704355) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.87256356) q[2];
sx q[2];
rz(-1.9930895) q[2];
sx q[2];
rz(0.97908212) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3124233) q[1];
sx q[1];
rz(-0.90065354) q[1];
sx q[1];
rz(1.6664684) q[1];
rz(-pi) q[2];
rz(-2.7510795) q[3];
sx q[3];
rz(-0.53005855) q[3];
sx q[3];
rz(-2.183225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.62362921) q[2];
sx q[2];
rz(-1.7771746) q[2];
sx q[2];
rz(-2.7015838) q[2];
rz(0.7157588) q[3];
sx q[3];
rz(-1.7093168) q[3];
sx q[3];
rz(-1.0796775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0004262) q[0];
sx q[0];
rz(-0.74581242) q[0];
sx q[0];
rz(2.0429042) q[0];
rz(-0.72775841) q[1];
sx q[1];
rz(-0.37574238) q[1];
sx q[1];
rz(-3.0922906) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54287275) q[0];
sx q[0];
rz(-2.1763986) q[0];
sx q[0];
rz(2.4332895) q[0];
rz(0.98408913) q[2];
sx q[2];
rz(-2.3684635) q[2];
sx q[2];
rz(0.099260515) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.8763435) q[1];
sx q[1];
rz(-2.0298404) q[1];
sx q[1];
rz(-2.6190119) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0827761) q[3];
sx q[3];
rz(-0.89142311) q[3];
sx q[3];
rz(-2.8454012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0631642) q[2];
sx q[2];
rz(-0.59946632) q[2];
sx q[2];
rz(0.72193974) q[2];
rz(2.1980964) q[3];
sx q[3];
rz(-2.3908581) q[3];
sx q[3];
rz(-2.8872484) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9938875) q[0];
sx q[0];
rz(-1.1575971) q[0];
sx q[0];
rz(1.0797427) q[0];
rz(-1.059277) q[1];
sx q[1];
rz(-2.9187027) q[1];
sx q[1];
rz(-1.4019029) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3106874) q[0];
sx q[0];
rz(-1.3584104) q[0];
sx q[0];
rz(-0.12916818) q[0];
x q[1];
rz(-1.5435113) q[2];
sx q[2];
rz(-1.3840884) q[2];
sx q[2];
rz(-2.7207431) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8733752) q[1];
sx q[1];
rz(-2.2443319) q[1];
sx q[1];
rz(2.7482277) q[1];
x q[2];
rz(2.510342) q[3];
sx q[3];
rz(-1.4025941) q[3];
sx q[3];
rz(1.9025308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6587276) q[2];
sx q[2];
rz(-1.7752703) q[2];
sx q[2];
rz(-1.6213017) q[2];
rz(-2.5907717) q[3];
sx q[3];
rz(-0.80533022) q[3];
sx q[3];
rz(2.4441392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14810066) q[0];
sx q[0];
rz(-1.8363331) q[0];
sx q[0];
rz(1.6114417) q[0];
rz(0.91611721) q[1];
sx q[1];
rz(-0.59090186) q[1];
sx q[1];
rz(-0.59060243) q[1];
rz(2.4089087) q[2];
sx q[2];
rz(-0.97326836) q[2];
sx q[2];
rz(1.2748981) q[2];
rz(-3.1237596) q[3];
sx q[3];
rz(-1.0537536) q[3];
sx q[3];
rz(-0.23585933) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
