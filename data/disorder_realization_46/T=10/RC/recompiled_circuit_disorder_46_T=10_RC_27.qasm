OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.5912136) q[0];
sx q[0];
rz(-0.0033291078) q[0];
sx q[0];
rz(2.79628) q[0];
rz(1.8058864) q[1];
sx q[1];
rz(-2.8023281) q[1];
sx q[1];
rz(-0.27944922) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1301491) q[0];
sx q[0];
rz(-2.7675417) q[0];
sx q[0];
rz(-0.98056294) q[0];
rz(-pi) q[1];
rz(0.34502132) q[2];
sx q[2];
rz(-0.81232386) q[2];
sx q[2];
rz(1.617384) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.86736673) q[1];
sx q[1];
rz(-2.9510731) q[1];
sx q[1];
rz(1.6412559) q[1];
rz(-pi) q[2];
rz(0.052993943) q[3];
sx q[3];
rz(-0.64265673) q[3];
sx q[3];
rz(1.3621804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.44148579) q[2];
sx q[2];
rz(-1.4725279) q[2];
sx q[2];
rz(0.95735615) q[2];
rz(2.8422614) q[3];
sx q[3];
rz(-2.744031) q[3];
sx q[3];
rz(-2.7295952) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9807724) q[0];
sx q[0];
rz(-0.94962025) q[0];
sx q[0];
rz(2.7278996) q[0];
rz(-1.7970239) q[1];
sx q[1];
rz(-0.78318703) q[1];
sx q[1];
rz(-0.63562524) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5168415) q[0];
sx q[0];
rz(-0.13953129) q[0];
sx q[0];
rz(-1.2463039) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9102155) q[2];
sx q[2];
rz(-3.0017827) q[2];
sx q[2];
rz(-1.2466873) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0932514) q[1];
sx q[1];
rz(-1.8257358) q[1];
sx q[1];
rz(-1.3202207) q[1];
rz(-pi) q[2];
rz(-0.73511519) q[3];
sx q[3];
rz(-0.8494091) q[3];
sx q[3];
rz(1.848156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.7522493) q[2];
sx q[2];
rz(-1.2717335) q[2];
sx q[2];
rz(2.6129369) q[2];
rz(1.901249) q[3];
sx q[3];
rz(-2.7820008) q[3];
sx q[3];
rz(2.8958029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56025958) q[0];
sx q[0];
rz(-1.154705) q[0];
sx q[0];
rz(-0.2581968) q[0];
rz(1.6351581) q[1];
sx q[1];
rz(-2.5841027) q[1];
sx q[1];
rz(-2.7071276) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4397944) q[0];
sx q[0];
rz(-0.5895624) q[0];
sx q[0];
rz(2.2586285) q[0];
x q[1];
rz(-0.065504727) q[2];
sx q[2];
rz(-1.0422049) q[2];
sx q[2];
rz(-2.5158109) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3420521) q[1];
sx q[1];
rz(-2.6375348) q[1];
sx q[1];
rz(-2.6964158) q[1];
rz(0.64819077) q[3];
sx q[3];
rz(-2.18581) q[3];
sx q[3];
rz(-2.0730413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.42052856) q[2];
sx q[2];
rz(-1.8177744) q[2];
sx q[2];
rz(-0.055796441) q[2];
rz(2.2037286) q[3];
sx q[3];
rz(-2.8580229) q[3];
sx q[3];
rz(2.8985033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0980804) q[0];
sx q[0];
rz(-0.94399095) q[0];
sx q[0];
rz(-3.0010624) q[0];
rz(-2.9699516) q[1];
sx q[1];
rz(-1.834603) q[1];
sx q[1];
rz(-0.26352873) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7292273) q[0];
sx q[0];
rz(-1.0050887) q[0];
sx q[0];
rz(-1.4501146) q[0];
rz(-0.61435917) q[2];
sx q[2];
rz(-0.76459568) q[2];
sx q[2];
rz(-1.6524397) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.22562379) q[1];
sx q[1];
rz(-1.9150503) q[1];
sx q[1];
rz(-0.72361372) q[1];
rz(-pi) q[2];
rz(0.88818355) q[3];
sx q[3];
rz(-2.7215951) q[3];
sx q[3];
rz(-0.94204599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0984829) q[2];
sx q[2];
rz(-0.3442328) q[2];
sx q[2];
rz(1.8919224) q[2];
rz(1.9469056) q[3];
sx q[3];
rz(-2.1134816) q[3];
sx q[3];
rz(-2.4627114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.93976218) q[0];
sx q[0];
rz(-1.285346) q[0];
sx q[0];
rz(-0.029065954) q[0];
rz(-1.4020231) q[1];
sx q[1];
rz(-2.7840835) q[1];
sx q[1];
rz(0.32863858) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11360725) q[0];
sx q[0];
rz(-2.0562045) q[0];
sx q[0];
rz(3.0263682) q[0];
rz(-pi) q[1];
rz(2.5673037) q[2];
sx q[2];
rz(-2.8238378) q[2];
sx q[2];
rz(0.53900915) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6824274) q[1];
sx q[1];
rz(-2.703856) q[1];
sx q[1];
rz(2.9912205) q[1];
rz(-pi) q[2];
rz(2.0629115) q[3];
sx q[3];
rz(-2.2552367) q[3];
sx q[3];
rz(-1.2217611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0236726) q[2];
sx q[2];
rz(-0.4549883) q[2];
sx q[2];
rz(0.16061352) q[2];
rz(1.9953856) q[3];
sx q[3];
rz(-1.8278154) q[3];
sx q[3];
rz(0.62079352) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49204957) q[0];
sx q[0];
rz(-2.2821125) q[0];
sx q[0];
rz(2.3550526) q[0];
rz(0.37711626) q[1];
sx q[1];
rz(-0.84723324) q[1];
sx q[1];
rz(2.699111) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15765685) q[0];
sx q[0];
rz(-1.2120976) q[0];
sx q[0];
rz(1.1789807) q[0];
x q[1];
rz(2.6857576) q[2];
sx q[2];
rz(-0.6243394) q[2];
sx q[2];
rz(-0.81941831) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.22873951) q[1];
sx q[1];
rz(-2.733426) q[1];
sx q[1];
rz(2.2250882) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0480012) q[3];
sx q[3];
rz(-1.8851265) q[3];
sx q[3];
rz(2.4623354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4601712) q[2];
sx q[2];
rz(-0.57839102) q[2];
sx q[2];
rz(2.1703413) q[2];
rz(2.4462637) q[3];
sx q[3];
rz(-2.9054502) q[3];
sx q[3];
rz(-2.7142081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-2.5044395) q[0];
sx q[0];
rz(-1.210286) q[0];
sx q[0];
rz(2.1321645) q[0];
rz(-2.5908453) q[1];
sx q[1];
rz(-1.2754722) q[1];
sx q[1];
rz(-1.9783463) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.163994) q[0];
sx q[0];
rz(-1.7118651) q[0];
sx q[0];
rz(-1.2376357) q[0];
rz(-pi) q[1];
rz(1.2802358) q[2];
sx q[2];
rz(-0.72939789) q[2];
sx q[2];
rz(-1.1682208) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7904661) q[1];
sx q[1];
rz(-0.27517056) q[1];
sx q[1];
rz(2.5009584) q[1];
rz(-1.2054687) q[3];
sx q[3];
rz(-2.4118773) q[3];
sx q[3];
rz(0.038289379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.78272468) q[2];
sx q[2];
rz(-1.9903368) q[2];
sx q[2];
rz(-2.3434095) q[2];
rz(0.77945566) q[3];
sx q[3];
rz(-0.53648406) q[3];
sx q[3];
rz(-1.1727758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9629795) q[0];
sx q[0];
rz(-2.6586752) q[0];
sx q[0];
rz(0.12284199) q[0];
rz(-3.0154862) q[1];
sx q[1];
rz(-1.5051196) q[1];
sx q[1];
rz(-1.2164446) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0053596) q[0];
sx q[0];
rz(-1.981712) q[0];
sx q[0];
rz(2.4634325) q[0];
rz(-pi) q[1];
x q[1];
rz(0.45546592) q[2];
sx q[2];
rz(-1.4205237) q[2];
sx q[2];
rz(-1.5474743) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0876906) q[1];
sx q[1];
rz(-1.9551827) q[1];
sx q[1];
rz(1.3575421) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.27463953) q[3];
sx q[3];
rz(-2.5317149) q[3];
sx q[3];
rz(-0.21851893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8806261) q[2];
sx q[2];
rz(-2.5239021) q[2];
sx q[2];
rz(-0.043126062) q[2];
rz(2.96636) q[3];
sx q[3];
rz(-0.8774811) q[3];
sx q[3];
rz(1.5568679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6373428) q[0];
sx q[0];
rz(-3.0872587) q[0];
sx q[0];
rz(-2.9328226) q[0];
rz(1.4978706) q[1];
sx q[1];
rz(-1.6521963) q[1];
sx q[1];
rz(2.0170905) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9546684) q[0];
sx q[0];
rz(-0.19321975) q[0];
sx q[0];
rz(2.3528071) q[0];
x q[1];
rz(0.71473177) q[2];
sx q[2];
rz(-1.7734314) q[2];
sx q[2];
rz(1.5243901) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.323656) q[1];
sx q[1];
rz(-1.7041823) q[1];
sx q[1];
rz(0.0011841983) q[1];
x q[2];
rz(-0.48891588) q[3];
sx q[3];
rz(-2.206114) q[3];
sx q[3];
rz(-0.92632252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0008529) q[2];
sx q[2];
rz(-0.72379392) q[2];
sx q[2];
rz(-2.9635584) q[2];
rz(1.595165) q[3];
sx q[3];
rz(-1.833257) q[3];
sx q[3];
rz(0.49062887) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7889325) q[0];
sx q[0];
rz(-1.9910318) q[0];
sx q[0];
rz(-2.1066522) q[0];
rz(0.79822284) q[1];
sx q[1];
rz(-2.1612576) q[1];
sx q[1];
rz(-0.14990526) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7222544) q[0];
sx q[0];
rz(-2.1968578) q[0];
sx q[0];
rz(-1.1115848) q[0];
x q[1];
rz(-0.28658861) q[2];
sx q[2];
rz(-1.5865241) q[2];
sx q[2];
rz(-2.36731) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.18969892) q[1];
sx q[1];
rz(-2.4687597) q[1];
sx q[1];
rz(0.876902) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8191765) q[3];
sx q[3];
rz(-0.90667778) q[3];
sx q[3];
rz(0.76457232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7297111) q[2];
sx q[2];
rz(-1.9103266) q[2];
sx q[2];
rz(-2.7330772) q[2];
rz(-2.216693) q[3];
sx q[3];
rz(-2.3291406) q[3];
sx q[3];
rz(-2.6598721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2395441) q[0];
sx q[0];
rz(-1.5477381) q[0];
sx q[0];
rz(1.5292194) q[0];
rz(1.3760024) q[1];
sx q[1];
rz(-1.9648432) q[1];
sx q[1];
rz(1.2479938) q[1];
rz(-2.3464936) q[2];
sx q[2];
rz(-1.3394525) q[2];
sx q[2];
rz(0.48223334) q[2];
rz(-0.61990191) q[3];
sx q[3];
rz(-1.3719659) q[3];
sx q[3];
rz(-0.93056783) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
