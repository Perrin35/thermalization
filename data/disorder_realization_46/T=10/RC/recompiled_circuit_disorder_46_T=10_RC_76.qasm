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
rz(-3.1382635) q[0];
sx q[0];
rz(0.34531265) q[0];
rz(-1.3357063) q[1];
sx q[1];
rz(-0.3392646) q[1];
sx q[1];
rz(0.27944922) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5292408) q[0];
sx q[0];
rz(-1.8792329) q[0];
sx q[0];
rz(2.926507) q[0];
rz(0.7819671) q[2];
sx q[2];
rz(-1.8188393) q[2];
sx q[2];
rz(0.19575191) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2024723) q[1];
sx q[1];
rz(-1.7608374) q[1];
sx q[1];
rz(3.128016) q[1];
rz(-pi) q[2];
x q[2];
rz(0.052993943) q[3];
sx q[3];
rz(-0.64265673) q[3];
sx q[3];
rz(-1.7794123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7001069) q[2];
sx q[2];
rz(-1.6690648) q[2];
sx q[2];
rz(2.1842365) q[2];
rz(2.8422614) q[3];
sx q[3];
rz(-0.39756164) q[3];
sx q[3];
rz(-0.41199747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9807724) q[0];
sx q[0];
rz(-0.94962025) q[0];
sx q[0];
rz(2.7278996) q[0];
rz(1.7970239) q[1];
sx q[1];
rz(-0.78318703) q[1];
sx q[1];
rz(-2.5059674) q[1];
rz(-pi) q[2];
x q[2];
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
rz(1.7027249) q[2];
sx q[2];
rz(-1.6172098) q[2];
sx q[2];
rz(-3.1293491) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7285068) q[1];
sx q[1];
rz(-1.8131078) q[1];
sx q[1];
rz(-0.26279022) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2222744) q[3];
sx q[3];
rz(-0.98005664) q[3];
sx q[3];
rz(0.35349333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.7522493) q[2];
sx q[2];
rz(-1.2717335) q[2];
sx q[2];
rz(-2.6129369) q[2];
rz(1.2403437) q[3];
sx q[3];
rz(-2.7820008) q[3];
sx q[3];
rz(0.24578978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5813331) q[0];
sx q[0];
rz(-1.9868877) q[0];
sx q[0];
rz(-2.8833959) q[0];
rz(-1.5064346) q[1];
sx q[1];
rz(-2.5841027) q[1];
sx q[1];
rz(0.43446508) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8718078) q[0];
sx q[0];
rz(-1.9315533) q[0];
sx q[0];
rz(2.0478134) q[0];
rz(-pi) q[1];
rz(-3.0760879) q[2];
sx q[2];
rz(-1.0422049) q[2];
sx q[2];
rz(2.5158109) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.30063054) q[1];
sx q[1];
rz(-1.1197487) q[1];
sx q[1];
rz(1.8039963) q[1];
x q[2];
rz(0.86359777) q[3];
sx q[3];
rz(-0.86175418) q[3];
sx q[3];
rz(-1.9883224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7210641) q[2];
sx q[2];
rz(-1.8177744) q[2];
sx q[2];
rz(-3.0857962) q[2];
rz(-2.2037286) q[3];
sx q[3];
rz(-0.28356975) q[3];
sx q[3];
rz(-0.24308932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.043512251) q[0];
sx q[0];
rz(-2.1976017) q[0];
sx q[0];
rz(-0.14053024) q[0];
rz(2.9699516) q[1];
sx q[1];
rz(-1.834603) q[1];
sx q[1];
rz(-2.8780639) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22334252) q[0];
sx q[0];
rz(-1.468987) q[0];
sx q[0];
rz(2.5725767) q[0];
rz(1.0657004) q[2];
sx q[2];
rz(-2.1720338) q[2];
sx q[2];
rz(-0.87841735) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0553953) q[1];
sx q[1];
rz(-0.89790422) q[1];
sx q[1];
rz(-1.124568) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.88818355) q[3];
sx q[3];
rz(-0.41999751) q[3];
sx q[3];
rz(-0.94204599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.043109743) q[2];
sx q[2];
rz(-2.7973599) q[2];
sx q[2];
rz(-1.2496703) q[2];
rz(1.9469056) q[3];
sx q[3];
rz(-2.1134816) q[3];
sx q[3];
rz(-2.4627114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2018305) q[0];
sx q[0];
rz(-1.8562466) q[0];
sx q[0];
rz(-0.029065954) q[0];
rz(-1.7395696) q[1];
sx q[1];
rz(-2.7840835) q[1];
sx q[1];
rz(2.8129541) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0279854) q[0];
sx q[0];
rz(-1.0853882) q[0];
sx q[0];
rz(-0.11522449) q[0];
rz(-pi) q[1];
rz(1.7476014) q[2];
sx q[2];
rz(-1.8362152) q[2];
sx q[2];
rz(3.0826498) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8935834) q[1];
sx q[1];
rz(-1.6343405) q[1];
sx q[1];
rz(-2.7081972) q[1];
rz(-pi) q[2];
rz(1.0786812) q[3];
sx q[3];
rz(-2.2552367) q[3];
sx q[3];
rz(1.2217611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0236726) q[2];
sx q[2];
rz(-0.4549883) q[2];
sx q[2];
rz(-0.16061352) q[2];
rz(-1.9953856) q[3];
sx q[3];
rz(-1.3137772) q[3];
sx q[3];
rz(0.62079352) q[3];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6495431) q[0];
sx q[0];
rz(-0.8594802) q[0];
sx q[0];
rz(2.3550526) q[0];
rz(0.37711626) q[1];
sx q[1];
rz(-2.2943594) q[1];
sx q[1];
rz(-2.699111) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0241644) q[0];
sx q[0];
rz(-2.6167343) q[0];
sx q[0];
rz(0.79458046) q[0];
x q[1];
rz(2.5673893) q[2];
sx q[2];
rz(-1.8310556) q[2];
sx q[2];
rz(-1.1299709) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6742299) q[1];
sx q[1];
rz(-1.2503887) q[1];
sx q[1];
rz(0.25735374) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8864185) q[3];
sx q[3];
rz(-1.481803) q[3];
sx q[3];
rz(0.86252585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.68142146) q[2];
sx q[2];
rz(-2.5632016) q[2];
sx q[2];
rz(-0.9712514) q[2];
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
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5044395) q[0];
sx q[0];
rz(-1.210286) q[0];
sx q[0];
rz(-2.1321645) q[0];
rz(2.5908453) q[1];
sx q[1];
rz(-1.8661205) q[1];
sx q[1];
rz(-1.9783463) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7834085) q[0];
sx q[0];
rz(-1.9005214) q[0];
sx q[0];
rz(-2.9924336) q[0];
rz(2.2789754) q[2];
sx q[2];
rz(-1.3786945) q[2];
sx q[2];
rz(2.5196599) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5442859) q[1];
sx q[1];
rz(-1.7339216) q[1];
sx q[1];
rz(2.9189928) q[1];
rz(1.2054687) q[3];
sx q[3];
rz(-2.4118773) q[3];
sx q[3];
rz(3.1033033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.358868) q[2];
sx q[2];
rz(-1.9903368) q[2];
sx q[2];
rz(0.7981832) q[2];
rz(2.362137) q[3];
sx q[3];
rz(-2.6051086) q[3];
sx q[3];
rz(-1.1727758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9629795) q[0];
sx q[0];
rz(-0.48291746) q[0];
sx q[0];
rz(3.0187507) q[0];
rz(0.12610647) q[1];
sx q[1];
rz(-1.5051196) q[1];
sx q[1];
rz(-1.2164446) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0053596) q[0];
sx q[0];
rz(-1.981712) q[0];
sx q[0];
rz(-2.4634325) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8100796) q[2];
sx q[2];
rz(-2.6636332) q[2];
sx q[2];
rz(2.8216459) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.46978894) q[1];
sx q[1];
rz(-2.7046013) q[1];
sx q[1];
rz(-0.48204084) q[1];
rz(0.27463953) q[3];
sx q[3];
rz(-0.60987771) q[3];
sx q[3];
rz(-0.21851893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8806261) q[2];
sx q[2];
rz(-0.61769056) q[2];
sx q[2];
rz(-3.0984666) q[2];
rz(-0.17523266) q[3];
sx q[3];
rz(-2.2641116) q[3];
sx q[3];
rz(1.5847248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50424987) q[0];
sx q[0];
rz(-3.0872587) q[0];
sx q[0];
rz(-2.9328226) q[0];
rz(1.4978706) q[1];
sx q[1];
rz(-1.6521963) q[1];
sx q[1];
rz(-1.1245022) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18692423) q[0];
sx q[0];
rz(-2.9483729) q[0];
sx q[0];
rz(-2.3528071) q[0];
rz(-pi) q[1];
rz(2.4268609) q[2];
sx q[2];
rz(-1.7734314) q[2];
sx q[2];
rz(-1.5243901) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.323656) q[1];
sx q[1];
rz(-1.7041823) q[1];
sx q[1];
rz(-3.1404085) q[1];
rz(-pi) q[2];
rz(2.6526768) q[3];
sx q[3];
rz(-0.93547869) q[3];
sx q[3];
rz(0.92632252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1407397) q[2];
sx q[2];
rz(-0.72379392) q[2];
sx q[2];
rz(-0.17803426) q[2];
rz(-1.595165) q[3];
sx q[3];
rz(-1.3083357) q[3];
sx q[3];
rz(0.49062887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7889325) q[0];
sx q[0];
rz(-1.9910318) q[0];
sx q[0];
rz(-1.0349405) q[0];
rz(0.79822284) q[1];
sx q[1];
rz(-2.1612576) q[1];
sx q[1];
rz(-0.14990526) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7184188) q[0];
sx q[0];
rz(-0.75773865) q[0];
sx q[0];
rz(2.5916879) q[0];
rz(-pi) q[1];
rz(-0.28658861) q[2];
sx q[2];
rz(-1.5550685) q[2];
sx q[2];
rz(-0.77428267) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3372911) q[1];
sx q[1];
rz(-1.1608487) q[1];
sx q[1];
rz(-1.0211584) q[1];
rz(-pi) q[2];
x q[2];
rz(0.30431872) q[3];
sx q[3];
rz(-0.70239866) q[3];
sx q[3];
rz(-1.9866634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7297111) q[2];
sx q[2];
rz(-1.231266) q[2];
sx q[2];
rz(0.40851545) q[2];
rz(-2.216693) q[3];
sx q[3];
rz(-2.3291406) q[3];
sx q[3];
rz(-2.6598721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.9020486) q[0];
sx q[0];
rz(-1.5477381) q[0];
sx q[0];
rz(1.5292194) q[0];
rz(-1.3760024) q[1];
sx q[1];
rz(-1.1767495) q[1];
sx q[1];
rz(-1.8935988) q[1];
rz(0.31870141) q[2];
sx q[2];
rz(-0.8209043) q[2];
sx q[2];
rz(2.2742297) q[2];
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
