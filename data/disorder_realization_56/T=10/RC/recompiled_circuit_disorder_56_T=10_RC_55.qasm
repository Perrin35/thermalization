OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.38874415) q[0];
sx q[0];
rz(-2.6053083) q[0];
sx q[0];
rz(0.94776881) q[0];
rz(1.8127958) q[1];
sx q[1];
rz(-1.2674018) q[1];
sx q[1];
rz(-1.0277494) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2078903) q[0];
sx q[0];
rz(-0.37651248) q[0];
sx q[0];
rz(-3.0794789) q[0];
x q[1];
rz(2.920354) q[2];
sx q[2];
rz(-1.0916296) q[2];
sx q[2];
rz(1.1269119) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4029018) q[1];
sx q[1];
rz(-1.1876904) q[1];
sx q[1];
rz(0.72523592) q[1];
x q[2];
rz(-0.62698934) q[3];
sx q[3];
rz(-1.9416182) q[3];
sx q[3];
rz(-0.8644608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1296922) q[2];
sx q[2];
rz(-1.4346069) q[2];
sx q[2];
rz(2.091308) q[2];
rz(2.0283279) q[3];
sx q[3];
rz(-2.2498825) q[3];
sx q[3];
rz(-3.0734857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.072409078) q[0];
sx q[0];
rz(-1.8962815) q[0];
sx q[0];
rz(0.29775277) q[0];
rz(-0.61966664) q[1];
sx q[1];
rz(-2.1344118) q[1];
sx q[1];
rz(-2.0334977) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28716921) q[0];
sx q[0];
rz(-0.65432917) q[0];
sx q[0];
rz(1.2229344) q[0];
x q[1];
rz(-0.42995288) q[2];
sx q[2];
rz(-1.0064831) q[2];
sx q[2];
rz(-1.0831837) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.019921692) q[1];
sx q[1];
rz(-0.53598511) q[1];
sx q[1];
rz(2.3977445) q[1];
rz(-pi) q[2];
x q[2];
rz(0.70257367) q[3];
sx q[3];
rz(-1.3630023) q[3];
sx q[3];
rz(-0.39263615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6790598) q[2];
sx q[2];
rz(-2.1647537) q[2];
sx q[2];
rz(-0.97529808) q[2];
rz(2.2235928) q[3];
sx q[3];
rz(-1.2851597) q[3];
sx q[3];
rz(-0.29618922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.179203) q[0];
sx q[0];
rz(-0.85150349) q[0];
sx q[0];
rz(0.54291022) q[0];
rz(0.88223282) q[1];
sx q[1];
rz(-2.0062607) q[1];
sx q[1];
rz(-0.96484819) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0438) q[0];
sx q[0];
rz(-1.291073) q[0];
sx q[0];
rz(2.8066737) q[0];
rz(1.5281048) q[2];
sx q[2];
rz(-1.8674769) q[2];
sx q[2];
rz(1.3670849) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1988586) q[1];
sx q[1];
rz(-0.13730362) q[1];
sx q[1];
rz(1.0820461) q[1];
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
rz(-pi) q[1];
rz(-3.0470011) q[2];
sx q[2];
rz(-2.5307405) q[2];
sx q[2];
rz(2.0084521) q[2];
rz(-0.23162332) q[3];
sx q[3];
rz(-1.2730205) q[3];
sx q[3];
rz(-2.384322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27947458) q[0];
sx q[0];
rz(-3.1311488) q[0];
sx q[0];
rz(1.7650771) q[0];
rz(-2.6230985) q[1];
sx q[1];
rz(-1.2644178) q[1];
sx q[1];
rz(2.8994697) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69232363) q[0];
sx q[0];
rz(-1.4081435) q[0];
sx q[0];
rz(0.67740324) q[0];
x q[1];
rz(-3.1130586) q[2];
sx q[2];
rz(-2.7118073) q[2];
sx q[2];
rz(-2.1536749) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0603234) q[1];
sx q[1];
rz(-0.4821061) q[1];
sx q[1];
rz(0.37270765) q[1];
rz(-pi) q[2];
rz(0.30947134) q[3];
sx q[3];
rz(-0.43635338) q[3];
sx q[3];
rz(-2.5040124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.018192856) q[2];
sx q[2];
rz(-0.96863666) q[2];
sx q[2];
rz(3.0467395) q[2];
rz(-1.3421966) q[3];
sx q[3];
rz(-1.3972524) q[3];
sx q[3];
rz(-2.7024787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3577394) q[0];
sx q[0];
rz(-1.3107603) q[0];
sx q[0];
rz(-2.7807996) q[0];
rz(-1.7533253) q[1];
sx q[1];
rz(-1.3307945) q[1];
sx q[1];
rz(2.0070019) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15123385) q[0];
sx q[0];
rz(-2.2915974) q[0];
sx q[0];
rz(-0.6806586) q[0];
rz(0.88419948) q[2];
sx q[2];
rz(-2.6066385) q[2];
sx q[2];
rz(-1.7810437) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6280917) q[1];
sx q[1];
rz(-2.6673632) q[1];
sx q[1];
rz(-0.34631108) q[1];
x q[2];
rz(-2.9633425) q[3];
sx q[3];
rz(-2.6986487) q[3];
sx q[3];
rz(2.6435542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0052884) q[2];
sx q[2];
rz(-1.4207999) q[2];
sx q[2];
rz(-0.57265442) q[2];
rz(2.2128361) q[3];
sx q[3];
rz(-2.6199665) q[3];
sx q[3];
rz(-1.9740392) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3271493) q[0];
sx q[0];
rz(-2.0454018) q[0];
sx q[0];
rz(-1.8922528) q[0];
rz(-1.918474) q[1];
sx q[1];
rz(-1.616281) q[1];
sx q[1];
rz(1.9893601) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5637701) q[0];
sx q[0];
rz(-1.676079) q[0];
sx q[0];
rz(0.012751243) q[0];
rz(-0.67201519) q[2];
sx q[2];
rz(-0.61215559) q[2];
sx q[2];
rz(1.9415346) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7631543) q[1];
sx q[1];
rz(-1.3707146) q[1];
sx q[1];
rz(-1.4640019) q[1];
rz(-pi) q[2];
rz(-2.4640718) q[3];
sx q[3];
rz(-2.065425) q[3];
sx q[3];
rz(-2.236931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.46889177) q[2];
sx q[2];
rz(-1.2142618) q[2];
sx q[2];
rz(2.1506298) q[2];
rz(2.4937566) q[3];
sx q[3];
rz(-0.9655374) q[3];
sx q[3];
rz(-0.34753862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25794849) q[0];
sx q[0];
rz(-0.22709665) q[0];
sx q[0];
rz(-0.062285475) q[0];
rz(2.9557872) q[1];
sx q[1];
rz(-1.6848911) q[1];
sx q[1];
rz(0.3947765) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.018935238) q[0];
sx q[0];
rz(-1.1766953) q[0];
sx q[0];
rz(1.3144486) q[0];
rz(-pi) q[1];
rz(-0.37624069) q[2];
sx q[2];
rz(-1.8707841) q[2];
sx q[2];
rz(2.8397727) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2784087) q[1];
sx q[1];
rz(-1.3235958) q[1];
sx q[1];
rz(0.43988887) q[1];
rz(0.23468252) q[3];
sx q[3];
rz(-0.15771401) q[3];
sx q[3];
rz(1.7705256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.64849598) q[2];
sx q[2];
rz(-0.33005565) q[2];
sx q[2];
rz(-2.8685692) q[2];
rz(-1.3027044) q[3];
sx q[3];
rz(-1.8282993) q[3];
sx q[3];
rz(0.31204143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7438695) q[0];
sx q[0];
rz(-2.0070772) q[0];
sx q[0];
rz(1.8564818) q[0];
rz(-1.5015191) q[1];
sx q[1];
rz(-1.7506426) q[1];
sx q[1];
rz(1.3407019) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9467981) q[0];
sx q[0];
rz(-1.9614944) q[0];
sx q[0];
rz(0.35238738) q[0];
rz(3.0807207) q[2];
sx q[2];
rz(-2.4448622) q[2];
sx q[2];
rz(0.75887647) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.33412877) q[1];
sx q[1];
rz(-2.6263116) q[1];
sx q[1];
rz(-2.3433102) q[1];
x q[2];
rz(2.1741381) q[3];
sx q[3];
rz(-0.52166044) q[3];
sx q[3];
rz(2.012901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1061219) q[2];
sx q[2];
rz(-2.3351228) q[2];
sx q[2];
rz(-2.1179874) q[2];
rz(-0.18493955) q[3];
sx q[3];
rz(-2.7513294) q[3];
sx q[3];
rz(-0.24188724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0446562) q[0];
sx q[0];
rz(-0.99656314) q[0];
sx q[0];
rz(-1.6212844) q[0];
rz(-2.8114491) q[1];
sx q[1];
rz(-1.2076999) q[1];
sx q[1];
rz(-0.83713371) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9159106) q[0];
sx q[0];
rz(-1.7726232) q[0];
sx q[0];
rz(2.4542144) q[0];
rz(-pi) q[1];
rz(1.2741954) q[2];
sx q[2];
rz(-1.0146078) q[2];
sx q[2];
rz(2.0402758) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7652119) q[1];
sx q[1];
rz(-1.7227255) q[1];
sx q[1];
rz(0.10581776) q[1];
x q[2];
rz(-1.0749531) q[3];
sx q[3];
rz(-1.8284203) q[3];
sx q[3];
rz(-2.6356634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.60124406) q[2];
sx q[2];
rz(-2.0566172) q[2];
sx q[2];
rz(2.9619651) q[2];
rz(2.1458697) q[3];
sx q[3];
rz(-1.2487753) q[3];
sx q[3];
rz(1.3109591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76374617) q[0];
sx q[0];
rz(-0.34559956) q[0];
sx q[0];
rz(2.0843704) q[0];
rz(-0.10748848) q[1];
sx q[1];
rz(-1.2534393) q[1];
sx q[1];
rz(2.1616139) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9419132) q[0];
sx q[0];
rz(-1.6964039) q[0];
sx q[0];
rz(3.0903387) q[0];
rz(-2.0699632) q[2];
sx q[2];
rz(-1.1071148) q[2];
sx q[2];
rz(1.0774563) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.92614782) q[1];
sx q[1];
rz(-2.0000334) q[1];
sx q[1];
rz(0.26305671) q[1];
rz(-1.8633217) q[3];
sx q[3];
rz(-0.33698002) q[3];
sx q[3];
rz(-1.8388336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3048627) q[2];
sx q[2];
rz(-1.9367846) q[2];
sx q[2];
rz(1.0277964) q[2];
rz(-1.7547539) q[3];
sx q[3];
rz(-1.3091062) q[3];
sx q[3];
rz(-0.28361472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2733611) q[0];
sx q[0];
rz(-2.1049451) q[0];
sx q[0];
rz(2.0275397) q[0];
rz(-0.83203075) q[1];
sx q[1];
rz(-2.6770626) q[1];
sx q[1];
rz(-2.4774036) q[1];
rz(-2.1297395) q[2];
sx q[2];
rz(-1.7218628) q[2];
sx q[2];
rz(-2.0414258) q[2];
rz(2.5084393) q[3];
sx q[3];
rz(-2.5452151) q[3];
sx q[3];
rz(2.7635318) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
