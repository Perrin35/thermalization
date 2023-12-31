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
rz(2.1138432) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1411113) q[0];
sx q[0];
rz(-1.9465465) q[0];
sx q[0];
rz(-1.5953338) q[0];
rz(-0.22123863) q[2];
sx q[2];
rz(-2.049963) q[2];
sx q[2];
rz(2.0146807) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9103968) q[1];
sx q[1];
rz(-0.80363552) q[1];
sx q[1];
rz(0.5459783) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5146033) q[3];
sx q[3];
rz(-1.9416182) q[3];
sx q[3];
rz(-0.8644608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0119005) q[2];
sx q[2];
rz(-1.4346069) q[2];
sx q[2];
rz(1.0502846) q[2];
rz(2.0283279) q[3];
sx q[3];
rz(-2.2498825) q[3];
sx q[3];
rz(0.068107001) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0691836) q[0];
sx q[0];
rz(-1.8962815) q[0];
sx q[0];
rz(2.8438399) q[0];
rz(0.61966664) q[1];
sx q[1];
rz(-1.0071808) q[1];
sx q[1];
rz(1.108095) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8544234) q[0];
sx q[0];
rz(-0.65432917) q[0];
sx q[0];
rz(1.9186583) q[0];
rz(0.98845311) q[2];
sx q[2];
rz(-2.4465912) q[2];
sx q[2];
rz(2.7671791) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.83924343) q[1];
sx q[1];
rz(-1.1855372) q[1];
sx q[1];
rz(-1.9531996) q[1];
rz(-pi) q[2];
x q[2];
rz(0.31538972) q[3];
sx q[3];
rz(-0.72761256) q[3];
sx q[3];
rz(1.417158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.46253282) q[2];
sx q[2];
rz(-2.1647537) q[2];
sx q[2];
rz(-0.97529808) q[2];
rz(2.2235928) q[3];
sx q[3];
rz(-1.8564329) q[3];
sx q[3];
rz(-2.8454034) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.179203) q[0];
sx q[0];
rz(-2.2900892) q[0];
sx q[0];
rz(2.5986824) q[0];
rz(2.2593598) q[1];
sx q[1];
rz(-2.0062607) q[1];
sx q[1];
rz(0.96484819) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99795139) q[0];
sx q[0];
rz(-2.7086341) q[0];
sx q[0];
rz(-2.4233682) q[0];
rz(-pi) q[1];
rz(2.8446571) q[2];
sx q[2];
rz(-1.6116217) q[2];
sx q[2];
rz(-2.9253935) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.112903) q[1];
sx q[1];
rz(-1.5064872) q[1];
sx q[1];
rz(1.4494004) q[1];
rz(-pi) q[2];
x q[2];
rz(0.19823719) q[3];
sx q[3];
rz(-1.1612411) q[3];
sx q[3];
rz(2.4026681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.09459153) q[2];
sx q[2];
rz(-0.61085218) q[2];
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
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27947458) q[0];
sx q[0];
rz(-0.010443895) q[0];
sx q[0];
rz(-1.3765155) q[0];
rz(2.6230985) q[1];
sx q[1];
rz(-1.8771749) q[1];
sx q[1];
rz(-0.24212295) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74894416) q[0];
sx q[0];
rz(-2.2376275) q[0];
sx q[0];
rz(-1.3632266) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7119615) q[2];
sx q[2];
rz(-1.5826844) q[2];
sx q[2];
rz(2.5846543) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.081269216) q[1];
sx q[1];
rz(-0.4821061) q[1];
sx q[1];
rz(-2.768885) q[1];
rz(-pi) q[2];
rz(0.30947134) q[3];
sx q[3];
rz(-2.7052393) q[3];
sx q[3];
rz(2.5040124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.1233998) q[2];
sx q[2];
rz(-0.96863666) q[2];
sx q[2];
rz(-3.0467395) q[2];
rz(1.3421966) q[3];
sx q[3];
rz(-1.7443402) q[3];
sx q[3];
rz(0.43911394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3577394) q[0];
sx q[0];
rz(-1.3107603) q[0];
sx q[0];
rz(2.7807996) q[0];
rz(1.7533253) q[1];
sx q[1];
rz(-1.8107982) q[1];
sx q[1];
rz(2.0070019) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.40588) q[0];
sx q[0];
rz(-2.1942733) q[0];
sx q[0];
rz(0.9491802) q[0];
rz(-pi) q[1];
rz(0.88419948) q[2];
sx q[2];
rz(-0.53495416) q[2];
sx q[2];
rz(-1.3605489) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5135009) q[1];
sx q[1];
rz(-0.47422945) q[1];
sx q[1];
rz(-0.34631108) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7048022) q[3];
sx q[3];
rz(-1.646864) q[3];
sx q[3];
rz(1.9074744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1363042) q[2];
sx q[2];
rz(-1.4207999) q[2];
sx q[2];
rz(2.5689382) q[2];
rz(-2.2128361) q[3];
sx q[3];
rz(-0.52162617) q[3];
sx q[3];
rz(1.1675534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8144433) q[0];
sx q[0];
rz(-2.0454018) q[0];
sx q[0];
rz(-1.8922528) q[0];
rz(-1.918474) q[1];
sx q[1];
rz(-1.5253116) q[1];
sx q[1];
rz(-1.9893601) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5778225) q[0];
sx q[0];
rz(-1.676079) q[0];
sx q[0];
rz(3.1288414) q[0];
rz(-pi) q[1];
rz(2.4695775) q[2];
sx q[2];
rz(-0.61215559) q[2];
sx q[2];
rz(1.9415346) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2684979) q[1];
sx q[1];
rz(-2.9151306) q[1];
sx q[1];
rz(2.6576659) q[1];
x q[2];
rz(-2.1762987) q[3];
sx q[3];
rz(-0.98635841) q[3];
sx q[3];
rz(0.30129978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.46889177) q[2];
sx q[2];
rz(-1.2142618) q[2];
sx q[2];
rz(2.1506298) q[2];
rz(-2.4937566) q[3];
sx q[3];
rz(-2.1760553) q[3];
sx q[3];
rz(-0.34753862) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8836442) q[0];
sx q[0];
rz(-2.914496) q[0];
sx q[0];
rz(0.062285475) q[0];
rz(0.1858055) q[1];
sx q[1];
rz(-1.4567016) q[1];
sx q[1];
rz(0.3947765) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4515526) q[0];
sx q[0];
rz(-1.3344904) q[0];
sx q[0];
rz(2.7355746) q[0];
rz(-pi) q[1];
rz(-1.8918745) q[2];
sx q[2];
rz(-1.9294538) q[2];
sx q[2];
rz(-1.1527588) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.77183206) q[1];
sx q[1];
rz(-0.50060111) q[1];
sx q[1];
rz(-0.53497772) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5338321) q[3];
sx q[3];
rz(-1.7241524) q[3];
sx q[3];
rz(1.5330029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4930967) q[2];
sx q[2];
rz(-0.33005565) q[2];
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
x q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7438695) q[0];
sx q[0];
rz(-2.0070772) q[0];
sx q[0];
rz(1.2851108) q[0];
rz(-1.5015191) q[1];
sx q[1];
rz(-1.7506426) q[1];
sx q[1];
rz(1.3407019) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1947945) q[0];
sx q[0];
rz(-1.9614944) q[0];
sx q[0];
rz(2.7892053) q[0];
rz(1.6216535) q[2];
sx q[2];
rz(-0.87561456) q[2];
sx q[2];
rz(2.4620172) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2015842) q[1];
sx q[1];
rz(-1.9218947) q[1];
sx q[1];
rz(-1.1854978) q[1];
rz(-pi) q[2];
rz(-0.96745456) q[3];
sx q[3];
rz(-0.52166044) q[3];
sx q[3];
rz(-1.1286917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0354707) q[2];
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
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0969365) q[0];
sx q[0];
rz(-0.99656314) q[0];
sx q[0];
rz(-1.6212844) q[0];
rz(0.3301436) q[1];
sx q[1];
rz(-1.2076999) q[1];
sx q[1];
rz(-0.83713371) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1820113) q[0];
sx q[0];
rz(-2.2416229) q[0];
sx q[0];
rz(-1.3120033) q[0];
rz(-0.43948549) q[2];
sx q[2];
rz(-2.5186933) q[2];
sx q[2];
rz(2.5650131) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7652119) q[1];
sx q[1];
rz(-1.4188671) q[1];
sx q[1];
rz(-0.10581776) q[1];
rz(-pi) q[2];
rz(0.29104851) q[3];
sx q[3];
rz(-1.0927199) q[3];
sx q[3];
rz(-0.92791286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5403486) q[2];
sx q[2];
rz(-2.0566172) q[2];
sx q[2];
rz(-2.9619651) q[2];
rz(2.1458697) q[3];
sx q[3];
rz(-1.8928173) q[3];
sx q[3];
rz(-1.3109591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3778465) q[0];
sx q[0];
rz(-2.7959931) q[0];
sx q[0];
rz(2.0843704) q[0];
rz(3.0341042) q[1];
sx q[1];
rz(-1.2534393) q[1];
sx q[1];
rz(2.1616139) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36469034) q[0];
sx q[0];
rz(-1.5199465) q[0];
sx q[0];
rz(-1.6965673) q[0];
rz(-pi) q[1];
x q[1];
rz(0.51771848) q[2];
sx q[2];
rz(-2.0132408) q[2];
sx q[2];
rz(0.25416086) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.500463) q[1];
sx q[1];
rz(-0.49912057) q[1];
sx q[1];
rz(-1.0541381) q[1];
rz(-pi) q[2];
rz(-1.2471334) q[3];
sx q[3];
rz(-1.4753046) q[3];
sx q[3];
rz(3.132706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8367299) q[2];
sx q[2];
rz(-1.2048081) q[2];
sx q[2];
rz(-1.0277964) q[2];
rz(1.7547539) q[3];
sx q[3];
rz(-1.8324865) q[3];
sx q[3];
rz(2.8579779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86823157) q[0];
sx q[0];
rz(-2.1049451) q[0];
sx q[0];
rz(2.0275397) q[0];
rz(-2.3095619) q[1];
sx q[1];
rz(-0.46453005) q[1];
sx q[1];
rz(0.66418905) q[1];
rz(-2.9639347) q[2];
sx q[2];
rz(-2.1226317) q[2];
sx q[2];
rz(2.5771099) q[2];
rz(2.6408623) q[3];
sx q[3];
rz(-1.9095608) q[3];
sx q[3];
rz(0.6469971) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
