OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.5503791) q[0];
sx q[0];
rz(3.1449218) q[0];
sx q[0];
rz(9.7700906) q[0];
rz(-1.3357063) q[1];
sx q[1];
rz(-0.3392646) q[1];
sx q[1];
rz(-2.8621434) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.011443519) q[0];
sx q[0];
rz(-0.37405095) q[0];
sx q[0];
rz(-2.1610297) q[0];
rz(-pi) q[1];
rz(0.34502132) q[2];
sx q[2];
rz(-0.81232386) q[2];
sx q[2];
rz(1.617384) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.63424078) q[1];
sx q[1];
rz(-1.5574641) q[1];
sx q[1];
rz(1.7608545) q[1];
rz(-3.0885987) q[3];
sx q[3];
rz(-0.64265673) q[3];
sx q[3];
rz(1.3621804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7001069) q[2];
sx q[2];
rz(-1.6690648) q[2];
sx q[2];
rz(-0.95735615) q[2];
rz(-0.29933128) q[3];
sx q[3];
rz(-0.39756164) q[3];
sx q[3];
rz(2.7295952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
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
rz(1.1608202) q[0];
sx q[0];
rz(-0.94962025) q[0];
sx q[0];
rz(-0.41369307) q[0];
rz(1.3445688) q[1];
sx q[1];
rz(-0.78318703) q[1];
sx q[1];
rz(-0.63562524) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7660852) q[0];
sx q[0];
rz(-1.5264395) q[0];
sx q[0];
rz(1.7031329) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.046819709) q[2];
sx q[2];
rz(-1.4390107) q[2];
sx q[2];
rz(1.552396) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.41308584) q[1];
sx q[1];
rz(-1.8131078) q[1];
sx q[1];
rz(2.8788024) q[1];
rz(2.2222744) q[3];
sx q[3];
rz(-2.161536) q[3];
sx q[3];
rz(-2.7880993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.7522493) q[2];
sx q[2];
rz(-1.2717335) q[2];
sx q[2];
rz(0.52865571) q[2];
rz(-1.901249) q[3];
sx q[3];
rz(-0.35959187) q[3];
sx q[3];
rz(-0.24578978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5813331) q[0];
sx q[0];
rz(-1.154705) q[0];
sx q[0];
rz(0.2581968) q[0];
rz(1.5064346) q[1];
sx q[1];
rz(-0.55748993) q[1];
sx q[1];
rz(-2.7071276) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70179825) q[0];
sx q[0];
rz(-2.5520303) q[0];
sx q[0];
rz(2.2586285) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.065504727) q[2];
sx q[2];
rz(-2.0993877) q[2];
sx q[2];
rz(-0.62578177) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.30063054) q[1];
sx q[1];
rz(-2.0218439) q[1];
sx q[1];
rz(1.8039963) q[1];
rz(-pi) q[2];
rz(-2.4934019) q[3];
sx q[3];
rz(-2.18581) q[3];
sx q[3];
rz(-2.0730413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7210641) q[2];
sx q[2];
rz(-1.3238182) q[2];
sx q[2];
rz(-0.055796441) q[2];
rz(-0.93786401) q[3];
sx q[3];
rz(-2.8580229) q[3];
sx q[3];
rz(-0.24308932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0980804) q[0];
sx q[0];
rz(-0.94399095) q[0];
sx q[0];
rz(-0.14053024) q[0];
rz(-0.17164104) q[1];
sx q[1];
rz(-1.834603) q[1];
sx q[1];
rz(0.26352873) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4123654) q[0];
sx q[0];
rz(-2.1365039) q[0];
sx q[0];
rz(1.4501146) q[0];
x q[1];
rz(-0.61435917) q[2];
sx q[2];
rz(-0.76459568) q[2];
sx q[2];
rz(1.489153) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.22562379) q[1];
sx q[1];
rz(-1.2265424) q[1];
sx q[1];
rz(0.72361372) q[1];
rz(1.2372381) q[3];
sx q[3];
rz(-1.8309438) q[3];
sx q[3];
rz(-0.0098269193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0984829) q[2];
sx q[2];
rz(-2.7973599) q[2];
sx q[2];
rz(-1.8919224) q[2];
rz(1.9469056) q[3];
sx q[3];
rz(-2.1134816) q[3];
sx q[3];
rz(0.67888129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2018305) q[0];
sx q[0];
rz(-1.8562466) q[0];
sx q[0];
rz(-3.1125267) q[0];
rz(1.7395696) q[1];
sx q[1];
rz(-2.7840835) q[1];
sx q[1];
rz(-2.8129541) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5111361) q[0];
sx q[0];
rz(-1.6726613) q[0];
sx q[0];
rz(1.0826375) q[0];
rz(-pi) q[1];
rz(-2.5673037) q[2];
sx q[2];
rz(-2.8238378) q[2];
sx q[2];
rz(-0.53900915) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6824274) q[1];
sx q[1];
rz(-0.43773663) q[1];
sx q[1];
rz(-0.15037219) q[1];
x q[2];
rz(-0.7469437) q[3];
sx q[3];
rz(-1.1960104) q[3];
sx q[3];
rz(2.4657616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.11792004) q[2];
sx q[2];
rz(-0.4549883) q[2];
sx q[2];
rz(2.9809791) q[2];
rz(-1.9953856) q[3];
sx q[3];
rz(-1.8278154) q[3];
sx q[3];
rz(2.5207991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49204957) q[0];
sx q[0];
rz(-0.8594802) q[0];
sx q[0];
rz(0.78654003) q[0];
rz(0.37711626) q[1];
sx q[1];
rz(-2.2943594) q[1];
sx q[1];
rz(0.4424817) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2690951) q[0];
sx q[0];
rz(-1.2051393) q[0];
sx q[0];
rz(0.38537607) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6857576) q[2];
sx q[2];
rz(-2.5172533) q[2];
sx q[2];
rz(0.81941831) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.22873951) q[1];
sx q[1];
rz(-0.40816669) q[1];
sx q[1];
rz(-2.2250882) q[1];
x q[2];
rz(-3.0480012) q[3];
sx q[3];
rz(-1.2564661) q[3];
sx q[3];
rz(-0.67925727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.68142146) q[2];
sx q[2];
rz(-0.57839102) q[2];
sx q[2];
rz(0.9712514) q[2];
rz(2.4462637) q[3];
sx q[3];
rz(-2.9054502) q[3];
sx q[3];
rz(-2.7142081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5044395) q[0];
sx q[0];
rz(-1.210286) q[0];
sx q[0];
rz(2.1321645) q[0];
rz(-0.55074739) q[1];
sx q[1];
rz(-1.8661205) q[1];
sx q[1];
rz(-1.9783463) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7834085) q[0];
sx q[0];
rz(-1.9005214) q[0];
sx q[0];
rz(-0.14915906) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8909056) q[2];
sx q[2];
rz(-2.2633268) q[2];
sx q[2];
rz(-0.78679774) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.35112652) q[1];
sx q[1];
rz(-2.8664221) q[1];
sx q[1];
rz(-0.64063425) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.30927741) q[3];
sx q[3];
rz(-0.89865548) q[3];
sx q[3];
rz(2.6291763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.358868) q[2];
sx q[2];
rz(-1.1512558) q[2];
sx q[2];
rz(-0.7981832) q[2];
rz(-0.77945566) q[3];
sx q[3];
rz(-2.6051086) q[3];
sx q[3];
rz(-1.1727758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.1786132) q[0];
sx q[0];
rz(-2.6586752) q[0];
sx q[0];
rz(3.0187507) q[0];
rz(-0.12610647) q[1];
sx q[1];
rz(-1.5051196) q[1];
sx q[1];
rz(1.2164446) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13623304) q[0];
sx q[0];
rz(-1.1598806) q[0];
sx q[0];
rz(0.67816011) q[0];
x q[1];
rz(-1.4037651) q[2];
sx q[2];
rz(-2.0207496) q[2];
sx q[2];
rz(3.0917167) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6718037) q[1];
sx q[1];
rz(-2.7046013) q[1];
sx q[1];
rz(2.6595518) q[1];
rz(-pi) q[2];
x q[2];
rz(0.59206483) q[3];
sx q[3];
rz(-1.7267623) q[3];
sx q[3];
rz(2.0162984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8806261) q[2];
sx q[2];
rz(-2.5239021) q[2];
sx q[2];
rz(0.043126062) q[2];
rz(-2.96636) q[3];
sx q[3];
rz(-2.2641116) q[3];
sx q[3];
rz(1.5568679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50424987) q[0];
sx q[0];
rz(-3.0872587) q[0];
sx q[0];
rz(-0.20877008) q[0];
rz(1.4978706) q[1];
sx q[1];
rz(-1.6521963) q[1];
sx q[1];
rz(2.0170905) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18692423) q[0];
sx q[0];
rz(-0.19321975) q[0];
sx q[0];
rz(0.78878553) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4268609) q[2];
sx q[2];
rz(-1.3681612) q[2];
sx q[2];
rz(-1.6172025) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.323656) q[1];
sx q[1];
rz(-1.4374104) q[1];
sx q[1];
rz(-0.0011841983) q[1];
rz(-pi) q[2];
rz(-0.48891588) q[3];
sx q[3];
rz(-0.93547869) q[3];
sx q[3];
rz(-2.2152701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1407397) q[2];
sx q[2];
rz(-2.4177987) q[2];
sx q[2];
rz(-0.17803426) q[2];
rz(1.5464276) q[3];
sx q[3];
rz(-1.3083357) q[3];
sx q[3];
rz(-2.6509638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35266018) q[0];
sx q[0];
rz(-1.9910318) q[0];
sx q[0];
rz(-1.0349405) q[0];
rz(2.3433698) q[1];
sx q[1];
rz(-0.98033506) q[1];
sx q[1];
rz(2.9916874) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7184188) q[0];
sx q[0];
rz(-0.75773865) q[0];
sx q[0];
rz(-0.54990479) q[0];
rz(1.5871928) q[2];
sx q[2];
rz(-1.2842442) q[2];
sx q[2];
rz(0.80114844) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0059709) q[1];
sx q[1];
rz(-1.0711728) q[1];
sx q[1];
rz(2.6702704) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8372739) q[3];
sx q[3];
rz(-2.439194) q[3];
sx q[3];
rz(1.9866634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7297111) q[2];
sx q[2];
rz(-1.231266) q[2];
sx q[2];
rz(2.7330772) q[2];
rz(-0.92489964) q[3];
sx q[3];
rz(-0.81245208) q[3];
sx q[3];
rz(-2.6598721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9020486) q[0];
sx q[0];
rz(-1.5477381) q[0];
sx q[0];
rz(1.5292194) q[0];
rz(1.7655903) q[1];
sx q[1];
rz(-1.1767495) q[1];
sx q[1];
rz(-1.8935988) q[1];
rz(-2.3464936) q[2];
sx q[2];
rz(-1.3394525) q[2];
sx q[2];
rz(0.48223334) q[2];
rz(1.3281214) q[3];
sx q[3];
rz(-2.1767053) q[3];
sx q[3];
rz(-2.3613031) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
