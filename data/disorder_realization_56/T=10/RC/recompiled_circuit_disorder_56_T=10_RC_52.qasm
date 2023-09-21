OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7528485) q[0];
sx q[0];
rz(-0.53628439) q[0];
sx q[0];
rz(2.1938238) q[0];
rz(1.8127958) q[1];
sx q[1];
rz(-1.2674018) q[1];
sx q[1];
rz(-1.0277494) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42067895) q[0];
sx q[0];
rz(-1.5479711) q[0];
sx q[0];
rz(-0.37585293) q[0];
rz(1.1711636) q[2];
sx q[2];
rz(-0.52414775) q[2];
sx q[2];
rz(1.5607967) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.23119588) q[1];
sx q[1];
rz(-0.80363552) q[1];
sx q[1];
rz(-2.5956144) q[1];
rz(-0.62698934) q[3];
sx q[3];
rz(-1.1999745) q[3];
sx q[3];
rz(-2.2771319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0119005) q[2];
sx q[2];
rz(-1.4346069) q[2];
sx q[2];
rz(-2.091308) q[2];
rz(-2.0283279) q[3];
sx q[3];
rz(-2.2498825) q[3];
sx q[3];
rz(3.0734857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0691836) q[0];
sx q[0];
rz(-1.2453112) q[0];
sx q[0];
rz(-0.29775277) q[0];
rz(-0.61966664) q[1];
sx q[1];
rz(-2.1344118) q[1];
sx q[1];
rz(1.108095) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0034804) q[0];
sx q[0];
rz(-1.3618042) q[0];
sx q[0];
rz(-0.94603993) q[0];
rz(-0.42995288) q[2];
sx q[2];
rz(-1.0064831) q[2];
sx q[2];
rz(2.058409) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.88156466) q[1];
sx q[1];
rz(-1.9238872) q[1];
sx q[1];
rz(-0.41207037) q[1];
rz(2.8262029) q[3];
sx q[3];
rz(-0.72761256) q[3];
sx q[3];
rz(1.7244347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.46253282) q[2];
sx q[2];
rz(-2.1647537) q[2];
sx q[2];
rz(-0.97529808) q[2];
rz(0.9179999) q[3];
sx q[3];
rz(-1.8564329) q[3];
sx q[3];
rz(2.8454034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96238962) q[0];
sx q[0];
rz(-0.85150349) q[0];
sx q[0];
rz(0.54291022) q[0];
rz(0.88223282) q[1];
sx q[1];
rz(-2.0062607) q[1];
sx q[1];
rz(-0.96484819) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1436413) q[0];
sx q[0];
rz(-2.7086341) q[0];
sx q[0];
rz(-2.4233682) q[0];
rz(-pi) q[1];
rz(-2.8446571) q[2];
sx q[2];
rz(-1.6116217) q[2];
sx q[2];
rz(-0.21619913) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.45005349) q[1];
sx q[1];
rz(-1.69194) q[1];
sx q[1];
rz(3.0768081) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1539677) q[3];
sx q[3];
rz(-1.3891451) q[3];
sx q[3];
rz(-2.3895398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0470011) q[2];
sx q[2];
rz(-0.61085218) q[2];
sx q[2];
rz(1.1331406) q[2];
rz(2.9099693) q[3];
sx q[3];
rz(-1.2730205) q[3];
sx q[3];
rz(0.75727063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27947458) q[0];
sx q[0];
rz(-3.1311488) q[0];
sx q[0];
rz(1.7650771) q[0];
rz(-0.51849413) q[1];
sx q[1];
rz(-1.2644178) q[1];
sx q[1];
rz(0.24212295) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3926485) q[0];
sx q[0];
rz(-2.2376275) q[0];
sx q[0];
rz(-1.7783661) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1130586) q[2];
sx q[2];
rz(-0.42978537) q[2];
sx q[2];
rz(-0.98791771) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8230069) q[1];
sx q[1];
rz(-1.4011523) q[1];
sx q[1];
rz(-2.6881048) q[1];
rz(-pi) q[2];
rz(-1.4297156) q[3];
sx q[3];
rz(-1.1564848) q[3];
sx q[3];
rz(-0.29841081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.018192856) q[2];
sx q[2];
rz(-0.96863666) q[2];
sx q[2];
rz(3.0467395) q[2];
rz(1.3421966) q[3];
sx q[3];
rz(-1.3972524) q[3];
sx q[3];
rz(2.7024787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3577394) q[0];
sx q[0];
rz(-1.8308324) q[0];
sx q[0];
rz(-2.7807996) q[0];
rz(-1.7533253) q[1];
sx q[1];
rz(-1.3307945) q[1];
sx q[1];
rz(-1.1345908) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15123385) q[0];
sx q[0];
rz(-0.84999527) q[0];
sx q[0];
rz(-2.4609341) q[0];
x q[1];
rz(-1.1410494) q[2];
sx q[2];
rz(-1.8998713) q[2];
sx q[2];
rz(-2.3171901) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5135009) q[1];
sx q[1];
rz(-0.47422945) q[1];
sx q[1];
rz(-2.7952816) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7048022) q[3];
sx q[3];
rz(-1.646864) q[3];
sx q[3];
rz(1.2341183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0052884) q[2];
sx q[2];
rz(-1.7207928) q[2];
sx q[2];
rz(-0.57265442) q[2];
rz(2.2128361) q[3];
sx q[3];
rz(-0.52162617) q[3];
sx q[3];
rz(-1.1675534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3271493) q[0];
sx q[0];
rz(-1.0961908) q[0];
sx q[0];
rz(-1.8922528) q[0];
rz(1.918474) q[1];
sx q[1];
rz(-1.616281) q[1];
sx q[1];
rz(1.1522326) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1359065) q[0];
sx q[0];
rz(-1.583477) q[0];
sx q[0];
rz(-1.4655051) q[0];
rz(-0.50243369) q[2];
sx q[2];
rz(-1.2049434) q[2];
sx q[2];
rz(-2.1937214) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2684979) q[1];
sx q[1];
rz(-2.9151306) q[1];
sx q[1];
rz(-0.48392673) q[1];
x q[2];
rz(2.4640718) q[3];
sx q[3];
rz(-1.0761677) q[3];
sx q[3];
rz(0.90466162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6727009) q[2];
sx q[2];
rz(-1.2142618) q[2];
sx q[2];
rz(0.99096283) q[2];
rz(-0.64783603) q[3];
sx q[3];
rz(-2.1760553) q[3];
sx q[3];
rz(0.34753862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8836442) q[0];
sx q[0];
rz(-0.22709665) q[0];
sx q[0];
rz(-3.0793072) q[0];
rz(-0.1858055) q[1];
sx q[1];
rz(-1.4567016) q[1];
sx q[1];
rz(-0.3947765) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.018935238) q[0];
sx q[0];
rz(-1.1766953) q[0];
sx q[0];
rz(1.3144486) q[0];
rz(-pi) q[1];
x q[1];
rz(0.6997509) q[2];
sx q[2];
rz(-2.664898) q[2];
sx q[2];
rz(-1.9112019) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.77183206) q[1];
sx q[1];
rz(-2.6409915) q[1];
sx q[1];
rz(0.53497772) q[1];
rz(0.15345927) q[3];
sx q[3];
rz(-1.5342661) q[3];
sx q[3];
rz(3.1094482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.64849598) q[2];
sx q[2];
rz(-2.811537) q[2];
sx q[2];
rz(0.27302343) q[2];
rz(1.8388883) q[3];
sx q[3];
rz(-1.8282993) q[3];
sx q[3];
rz(-2.8295512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(2.7438695) q[0];
sx q[0];
rz(-2.0070772) q[0];
sx q[0];
rz(-1.8564818) q[0];
rz(1.6400736) q[1];
sx q[1];
rz(-1.39095) q[1];
sx q[1];
rz(-1.3407019) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9467981) q[0];
sx q[0];
rz(-1.1800982) q[0];
sx q[0];
rz(2.7892053) q[0];
rz(-pi) q[1];
x q[1];
rz(0.060872002) q[2];
sx q[2];
rz(-0.69673046) q[2];
sx q[2];
rz(0.75887647) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9400085) q[1];
sx q[1];
rz(-1.9218947) q[1];
sx q[1];
rz(1.1854978) q[1];
rz(-pi) q[2];
rz(1.1287442) q[3];
sx q[3];
rz(-1.8574517) q[3];
sx q[3];
rz(3.0451881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0354707) q[2];
sx q[2];
rz(-2.3351228) q[2];
sx q[2];
rz(2.1179874) q[2];
rz(-0.18493955) q[3];
sx q[3];
rz(-2.7513294) q[3];
sx q[3];
rz(-0.24188724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0446562) q[0];
sx q[0];
rz(-2.1450295) q[0];
sx q[0];
rz(1.5203083) q[0];
rz(2.8114491) q[1];
sx q[1];
rz(-1.9338927) q[1];
sx q[1];
rz(-0.83713371) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22568208) q[0];
sx q[0];
rz(-1.3689694) q[0];
sx q[0];
rz(0.68737824) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5651921) q[2];
sx q[2];
rz(-1.3199558) q[2];
sx q[2];
rz(2.5121411) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3763807) q[1];
sx q[1];
rz(-1.7227255) q[1];
sx q[1];
rz(0.10581776) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8505441) q[3];
sx q[3];
rz(-2.0488727) q[3];
sx q[3];
rz(0.92791286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5403486) q[2];
sx q[2];
rz(-2.0566172) q[2];
sx q[2];
rz(2.9619651) q[2];
rz(0.99572292) q[3];
sx q[3];
rz(-1.2487753) q[3];
sx q[3];
rz(-1.3109591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3778465) q[0];
sx q[0];
rz(-2.7959931) q[0];
sx q[0];
rz(-2.0843704) q[0];
rz(-0.10748848) q[1];
sx q[1];
rz(-1.2534393) q[1];
sx q[1];
rz(2.1616139) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36469034) q[0];
sx q[0];
rz(-1.5199465) q[0];
sx q[0];
rz(-1.4450253) q[0];
rz(-pi) q[1];
rz(0.76358958) q[2];
sx q[2];
rz(-0.66765235) q[2];
sx q[2];
rz(-1.180336) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.53303888) q[1];
sx q[1];
rz(-1.8095008) q[1];
sx q[1];
rz(-2.0134316) q[1];
x q[2];
rz(-0.1006871) q[3];
sx q[3];
rz(-1.2486613) q[3];
sx q[3];
rz(1.5299357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8367299) q[2];
sx q[2];
rz(-1.9367846) q[2];
sx q[2];
rz(-2.1137962) q[2];
rz(-1.3868388) q[3];
sx q[3];
rz(-1.8324865) q[3];
sx q[3];
rz(-0.28361472) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86823157) q[0];
sx q[0];
rz(-1.0366476) q[0];
sx q[0];
rz(-1.114053) q[0];
rz(2.3095619) q[1];
sx q[1];
rz(-2.6770626) q[1];
sx q[1];
rz(-2.4774036) q[1];
rz(2.9639347) q[2];
sx q[2];
rz(-1.0189609) q[2];
sx q[2];
rz(-0.56448274) q[2];
rz(1.9527312) q[3];
sx q[3];
rz(-2.0406796) q[3];
sx q[3];
rz(-1.1036967) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];