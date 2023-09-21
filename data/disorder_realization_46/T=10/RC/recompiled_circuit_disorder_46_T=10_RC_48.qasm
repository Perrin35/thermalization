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
rz(1.8058864) q[1];
sx q[1];
rz(-2.8023281) q[1];
sx q[1];
rz(-0.27944922) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61235185) q[0];
sx q[0];
rz(-1.2623598) q[0];
sx q[0];
rz(-0.2150857) q[0];
x q[1];
rz(-2.7965713) q[2];
sx q[2];
rz(-0.81232386) q[2];
sx q[2];
rz(1.617384) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.86736673) q[1];
sx q[1];
rz(-2.9510731) q[1];
sx q[1];
rz(1.5003367) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6104326) q[3];
sx q[3];
rz(-0.92919022) q[3];
sx q[3];
rz(1.2960145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.44148579) q[2];
sx q[2];
rz(-1.4725279) q[2];
sx q[2];
rz(-0.95735615) q[2];
rz(2.8422614) q[3];
sx q[3];
rz(-2.744031) q[3];
sx q[3];
rz(0.41199747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1608202) q[0];
sx q[0];
rz(-2.1919724) q[0];
sx q[0];
rz(-2.7278996) q[0];
rz(1.7970239) q[1];
sx q[1];
rz(-2.3584056) q[1];
sx q[1];
rz(-0.63562524) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6247511) q[0];
sx q[0];
rz(-0.13953129) q[0];
sx q[0];
rz(1.2463039) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4388678) q[2];
sx q[2];
rz(-1.6172098) q[2];
sx q[2];
rz(0.012243587) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0483413) q[1];
sx q[1];
rz(-1.3158568) q[1];
sx q[1];
rz(-1.3202207) q[1];
rz(-pi) q[2];
rz(-2.440968) q[3];
sx q[3];
rz(-2.0985588) q[3];
sx q[3];
rz(-2.3259195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.7522493) q[2];
sx q[2];
rz(-1.2717335) q[2];
sx q[2];
rz(0.52865571) q[2];
rz(1.901249) q[3];
sx q[3];
rz(-2.7820008) q[3];
sx q[3];
rz(2.8958029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56025958) q[0];
sx q[0];
rz(-1.9868877) q[0];
sx q[0];
rz(-0.2581968) q[0];
rz(1.6351581) q[1];
sx q[1];
rz(-2.5841027) q[1];
sx q[1];
rz(-2.7071276) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8718078) q[0];
sx q[0];
rz(-1.2100394) q[0];
sx q[0];
rz(-1.0937793) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1003232) q[2];
sx q[2];
rz(-1.5142421) q[2];
sx q[2];
rz(2.2296485) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3420521) q[1];
sx q[1];
rz(-2.6375348) q[1];
sx q[1];
rz(-2.6964158) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2779949) q[3];
sx q[3];
rz(-0.86175418) q[3];
sx q[3];
rz(-1.9883224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.42052856) q[2];
sx q[2];
rz(-1.8177744) q[2];
sx q[2];
rz(-0.055796441) q[2];
rz(0.93786401) q[3];
sx q[3];
rz(-2.8580229) q[3];
sx q[3];
rz(0.24308932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0980804) q[0];
sx q[0];
rz(-2.1976017) q[0];
sx q[0];
rz(-3.0010624) q[0];
rz(2.9699516) q[1];
sx q[1];
rz(-1.3069897) q[1];
sx q[1];
rz(2.8780639) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7292273) q[0];
sx q[0];
rz(-2.1365039) q[0];
sx q[0];
rz(-1.6914781) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0758923) q[2];
sx q[2];
rz(-2.1720338) q[2];
sx q[2];
rz(0.87841735) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9159689) q[1];
sx q[1];
rz(-1.2265424) q[1];
sx q[1];
rz(0.72361372) q[1];
rz(-pi) q[2];
rz(-1.2372381) q[3];
sx q[3];
rz(-1.3106489) q[3];
sx q[3];
rz(3.1317657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0984829) q[2];
sx q[2];
rz(-0.3442328) q[2];
sx q[2];
rz(1.2496703) q[2];
rz(1.9469056) q[3];
sx q[3];
rz(-2.1134816) q[3];
sx q[3];
rz(-2.4627114) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2018305) q[0];
sx q[0];
rz(-1.285346) q[0];
sx q[0];
rz(0.029065954) q[0];
rz(-1.7395696) q[1];
sx q[1];
rz(-2.7840835) q[1];
sx q[1];
rz(2.8129541) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5111361) q[0];
sx q[0];
rz(-1.4689313) q[0];
sx q[0];
rz(-2.0589552) q[0];
x q[1];
rz(-0.57428898) q[2];
sx q[2];
rz(-0.31775489) q[2];
sx q[2];
rz(2.6025835) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8481816) q[1];
sx q[1];
rz(-2.003258) q[1];
sx q[1];
rz(-1.6407938) q[1];
rz(-pi) q[2];
rz(-0.7469437) q[3];
sx q[3];
rz(-1.9455823) q[3];
sx q[3];
rz(-2.4657616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.11792004) q[2];
sx q[2];
rz(-2.6866044) q[2];
sx q[2];
rz(-2.9809791) q[2];
rz(-1.1462071) q[3];
sx q[3];
rz(-1.3137772) q[3];
sx q[3];
rz(-0.62079352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-0.49204957) q[0];
sx q[0];
rz(-0.8594802) q[0];
sx q[0];
rz(-2.3550526) q[0];
rz(-0.37711626) q[1];
sx q[1];
rz(-2.2943594) q[1];
sx q[1];
rz(2.699111) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9839358) q[0];
sx q[0];
rz(-1.2120976) q[0];
sx q[0];
rz(-1.1789807) q[0];
rz(-pi) q[1];
rz(-1.2636678) q[2];
sx q[2];
rz(-2.1233635) q[2];
sx q[2];
rz(0.27586684) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.46736273) q[1];
sx q[1];
rz(-1.2503887) q[1];
sx q[1];
rz(-2.8842389) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0480012) q[3];
sx q[3];
rz(-1.8851265) q[3];
sx q[3];
rz(-0.67925727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4601712) q[2];
sx q[2];
rz(-0.57839102) q[2];
sx q[2];
rz(0.9712514) q[2];
rz(-0.69532895) q[3];
sx q[3];
rz(-0.23614241) q[3];
sx q[3];
rz(2.7142081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(0.63715315) q[0];
sx q[0];
rz(-1.210286) q[0];
sx q[0];
rz(2.1321645) q[0];
rz(2.5908453) q[1];
sx q[1];
rz(-1.8661205) q[1];
sx q[1];
rz(-1.9783463) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3488719) q[0];
sx q[0];
rz(-0.36076818) q[0];
sx q[0];
rz(1.1611206) q[0];
x q[1];
rz(1.8613569) q[2];
sx q[2];
rz(-2.4121948) q[2];
sx q[2];
rz(1.9733719) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7904661) q[1];
sx q[1];
rz(-2.8664221) q[1];
sx q[1];
rz(2.5009584) q[1];
rz(-pi) q[2];
rz(0.87485119) q[3];
sx q[3];
rz(-1.8112744) q[3];
sx q[3];
rz(1.2547753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.358868) q[2];
sx q[2];
rz(-1.1512558) q[2];
sx q[2];
rz(-2.3434095) q[2];
rz(-0.77945566) q[3];
sx q[3];
rz(-0.53648406) q[3];
sx q[3];
rz(-1.9688169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1786132) q[0];
sx q[0];
rz(-2.6586752) q[0];
sx q[0];
rz(0.12284199) q[0];
rz(-0.12610647) q[1];
sx q[1];
rz(-1.6364731) q[1];
sx q[1];
rz(1.925148) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2467263) q[0];
sx q[0];
rz(-2.3658532) q[0];
sx q[0];
rz(2.5345483) q[0];
rz(1.7378275) q[2];
sx q[2];
rz(-2.0207496) q[2];
sx q[2];
rz(-0.049875967) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6718037) q[1];
sx q[1];
rz(-0.43699139) q[1];
sx q[1];
rz(0.48204084) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.59206483) q[3];
sx q[3];
rz(-1.7267623) q[3];
sx q[3];
rz(1.1252943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8806261) q[2];
sx q[2];
rz(-0.61769056) q[2];
sx q[2];
rz(-3.0984666) q[2];
rz(-0.17523266) q[3];
sx q[3];
rz(-2.2641116) q[3];
sx q[3];
rz(-1.5568679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50424987) q[0];
sx q[0];
rz(-3.0872587) q[0];
sx q[0];
rz(-0.20877008) q[0];
rz(1.6437221) q[1];
sx q[1];
rz(-1.6521963) q[1];
sx q[1];
rz(1.1245022) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61125206) q[0];
sx q[0];
rz(-1.4350622) q[0];
sx q[0];
rz(-1.432857) q[0];
rz(0.71473177) q[2];
sx q[2];
rz(-1.3681612) q[2];
sx q[2];
rz(-1.5243901) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.323656) q[1];
sx q[1];
rz(-1.4374104) q[1];
sx q[1];
rz(-3.1404085) q[1];
rz(-pi) q[2];
rz(-0.87499683) q[3];
sx q[3];
rz(-1.9584624) q[3];
sx q[3];
rz(-2.8029203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1407397) q[2];
sx q[2];
rz(-0.72379392) q[2];
sx q[2];
rz(0.17803426) q[2];
rz(-1.5464276) q[3];
sx q[3];
rz(-1.833257) q[3];
sx q[3];
rz(0.49062887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35266018) q[0];
sx q[0];
rz(-1.9910318) q[0];
sx q[0];
rz(1.0349405) q[0];
rz(-2.3433698) q[1];
sx q[1];
rz(-0.98033506) q[1];
sx q[1];
rz(-2.9916874) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4231739) q[0];
sx q[0];
rz(-0.75773865) q[0];
sx q[0];
rz(-0.54990479) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.055585102) q[2];
sx q[2];
rz(-0.28700799) q[2];
sx q[2];
rz(0.74319786) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3372911) q[1];
sx q[1];
rz(-1.1608487) q[1];
sx q[1];
rz(2.1204342) q[1];
x q[2];
rz(-0.67930119) q[3];
sx q[3];
rz(-1.7656109) q[3];
sx q[3];
rz(-0.65115813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4118816) q[2];
sx q[2];
rz(-1.9103266) q[2];
sx q[2];
rz(2.7330772) q[2];
rz(-2.216693) q[3];
sx q[3];
rz(-2.3291406) q[3];
sx q[3];
rz(0.48172054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2395441) q[0];
sx q[0];
rz(-1.5938546) q[0];
sx q[0];
rz(-1.6123733) q[0];
rz(-1.7655903) q[1];
sx q[1];
rz(-1.9648432) q[1];
sx q[1];
rz(1.2479938) q[1];
rz(-0.79509905) q[2];
sx q[2];
rz(-1.8021402) q[2];
sx q[2];
rz(-2.6593593) q[2];
rz(2.5216907) q[3];
sx q[3];
rz(-1.3719659) q[3];
sx q[3];
rz(-0.93056783) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
