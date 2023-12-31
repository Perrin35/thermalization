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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61235185) q[0];
sx q[0];
rz(-1.2623598) q[0];
sx q[0];
rz(0.2150857) q[0];
rz(-pi) q[1];
rz(2.3596256) q[2];
sx q[2];
rz(-1.8188393) q[2];
sx q[2];
rz(-0.19575191) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.86736673) q[1];
sx q[1];
rz(-2.9510731) q[1];
sx q[1];
rz(1.5003367) q[1];
rz(-pi) q[2];
x q[2];
rz(0.052993943) q[3];
sx q[3];
rz(-2.4989359) q[3];
sx q[3];
rz(1.7794123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.44148579) q[2];
sx q[2];
rz(-1.4725279) q[2];
sx q[2];
rz(0.95735615) q[2];
rz(-2.8422614) q[3];
sx q[3];
rz(-2.744031) q[3];
sx q[3];
rz(-0.41199747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1608202) q[0];
sx q[0];
rz(-2.1919724) q[0];
sx q[0];
rz(0.41369307) q[0];
rz(1.3445688) q[1];
sx q[1];
rz(-2.3584056) q[1];
sx q[1];
rz(-2.5059674) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18938633) q[0];
sx q[0];
rz(-1.4385907) q[0];
sx q[0];
rz(3.0968451) q[0];
x q[1];
rz(3.0947729) q[2];
sx q[2];
rz(-1.702582) q[2];
sx q[2];
rz(-1.552396) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8860652) q[1];
sx q[1];
rz(-2.7860502) q[1];
sx q[1];
rz(-0.76053263) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4064775) q[3];
sx q[3];
rz(-0.8494091) q[3];
sx q[3];
rz(-1.848156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.7522493) q[2];
sx q[2];
rz(-1.8698591) q[2];
sx q[2];
rz(2.6129369) q[2];
rz(-1.901249) q[3];
sx q[3];
rz(-2.7820008) q[3];
sx q[3];
rz(0.24578978) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56025958) q[0];
sx q[0];
rz(-1.154705) q[0];
sx q[0];
rz(0.2581968) q[0];
rz(-1.5064346) q[1];
sx q[1];
rz(-0.55748993) q[1];
sx q[1];
rz(2.7071276) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70179825) q[0];
sx q[0];
rz(-0.5895624) q[0];
sx q[0];
rz(-0.88296417) q[0];
rz(-pi) q[1];
rz(-1.0412695) q[2];
sx q[2];
rz(-1.5142421) q[2];
sx q[2];
rz(2.2296485) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.79954051) q[1];
sx q[1];
rz(-0.5040579) q[1];
sx q[1];
rz(-0.44517681) q[1];
x q[2];
rz(-0.86359777) q[3];
sx q[3];
rz(-0.86175418) q[3];
sx q[3];
rz(1.9883224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7210641) q[2];
sx q[2];
rz(-1.3238182) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(3.0980804) q[0];
sx q[0];
rz(-0.94399095) q[0];
sx q[0];
rz(-0.14053024) q[0];
rz(0.17164104) q[1];
sx q[1];
rz(-1.834603) q[1];
sx q[1];
rz(2.8780639) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4123654) q[0];
sx q[0];
rz(-1.0050887) q[0];
sx q[0];
rz(1.4501146) q[0];
rz(-0.61435917) q[2];
sx q[2];
rz(-2.376997) q[2];
sx q[2];
rz(-1.489153) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0553953) q[1];
sx q[1];
rz(-0.89790422) q[1];
sx q[1];
rz(-1.124568) q[1];
x q[2];
rz(2.2534091) q[3];
sx q[3];
rz(-0.41999751) q[3];
sx q[3];
rz(2.1995467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0984829) q[2];
sx q[2];
rz(-2.7973599) q[2];
sx q[2];
rz(-1.8919224) q[2];
rz(1.9469056) q[3];
sx q[3];
rz(-1.0281111) q[3];
sx q[3];
rz(2.4627114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2018305) q[0];
sx q[0];
rz(-1.8562466) q[0];
sx q[0];
rz(0.029065954) q[0];
rz(1.7395696) q[1];
sx q[1];
rz(-0.35750917) q[1];
sx q[1];
rz(2.8129541) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11360725) q[0];
sx q[0];
rz(-2.0562045) q[0];
sx q[0];
rz(-0.11522449) q[0];
rz(-1.7476014) q[2];
sx q[2];
rz(-1.3053775) q[2];
sx q[2];
rz(-0.058942827) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.45916522) q[1];
sx q[1];
rz(-2.703856) q[1];
sx q[1];
rz(-2.9912205) q[1];
rz(-pi) q[2];
x q[2];
rz(0.7469437) q[3];
sx q[3];
rz(-1.9455823) q[3];
sx q[3];
rz(2.4657616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0236726) q[2];
sx q[2];
rz(-0.4549883) q[2];
sx q[2];
rz(-2.9809791) q[2];
rz(-1.9953856) q[3];
sx q[3];
rz(-1.3137772) q[3];
sx q[3];
rz(0.62079352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-0.84723324) q[1];
sx q[1];
rz(2.699111) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1174283) q[0];
sx q[0];
rz(-2.6167343) q[0];
sx q[0];
rz(2.3470122) q[0];
rz(-pi) q[1];
rz(-0.45583506) q[2];
sx q[2];
rz(-2.5172533) q[2];
sx q[2];
rz(0.81941831) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6742299) q[1];
sx q[1];
rz(-1.891204) q[1];
sx q[1];
rz(-0.25735374) q[1];
rz(3.0480012) q[3];
sx q[3];
rz(-1.8851265) q[3];
sx q[3];
rz(-0.67925727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4601712) q[2];
sx q[2];
rz(-2.5632016) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
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
rz(-1.2754722) q[1];
sx q[1];
rz(1.9783463) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7834085) q[0];
sx q[0];
rz(-1.2410713) q[0];
sx q[0];
rz(-2.9924336) q[0];
rz(-1.8613569) q[2];
sx q[2];
rz(-0.72939789) q[2];
sx q[2];
rz(-1.1682208) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.59730676) q[1];
sx q[1];
rz(-1.407671) q[1];
sx q[1];
rz(0.2225999) q[1];
rz(0.87485119) q[3];
sx q[3];
rz(-1.8112744) q[3];
sx q[3];
rz(-1.8868173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.358868) q[2];
sx q[2];
rz(-1.1512558) q[2];
sx q[2];
rz(0.7981832) q[2];
rz(-0.77945566) q[3];
sx q[3];
rz(-2.6051086) q[3];
sx q[3];
rz(1.9688169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1786132) q[0];
sx q[0];
rz(-0.48291746) q[0];
sx q[0];
rz(-3.0187507) q[0];
rz(3.0154862) q[1];
sx q[1];
rz(-1.6364731) q[1];
sx q[1];
rz(-1.2164446) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8948664) q[0];
sx q[0];
rz(-0.77573949) q[0];
sx q[0];
rz(-2.5345483) q[0];
rz(-1.7378275) q[2];
sx q[2];
rz(-2.0207496) q[2];
sx q[2];
rz(-3.0917167) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.053902102) q[1];
sx q[1];
rz(-1.1864099) q[1];
sx q[1];
rz(-1.3575421) q[1];
rz(-2.5495278) q[3];
sx q[3];
rz(-1.4148303) q[3];
sx q[3];
rz(1.1252943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8806261) q[2];
sx q[2];
rz(-0.61769056) q[2];
sx q[2];
rz(-3.0984666) q[2];
rz(2.96636) q[3];
sx q[3];
rz(-2.2641116) q[3];
sx q[3];
rz(1.5847248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.6373428) q[0];
sx q[0];
rz(-3.0872587) q[0];
sx q[0];
rz(2.9328226) q[0];
rz(1.6437221) q[1];
sx q[1];
rz(-1.6521963) q[1];
sx q[1];
rz(-2.0170905) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18692423) q[0];
sx q[0];
rz(-0.19321975) q[0];
sx q[0];
rz(2.3528071) q[0];
rz(0.30377109) q[2];
sx q[2];
rz(-0.73799947) q[2];
sx q[2];
rz(-2.8673025) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.75301718) q[1];
sx q[1];
rz(-1.57197) q[1];
sx q[1];
rz(-1.4374103) q[1];
rz(-pi) q[2];
rz(-2.1379925) q[3];
sx q[3];
rz(-2.3611259) q[3];
sx q[3];
rz(1.6571852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0008529) q[2];
sx q[2];
rz(-2.4177987) q[2];
sx q[2];
rz(0.17803426) q[2];
rz(-1.595165) q[3];
sx q[3];
rz(-1.3083357) q[3];
sx q[3];
rz(0.49062887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
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
rz(1.0349405) q[0];
rz(-2.3433698) q[1];
sx q[1];
rz(-0.98033506) q[1];
sx q[1];
rz(-2.9916874) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4193383) q[0];
sx q[0];
rz(-0.94473487) q[0];
sx q[0];
rz(2.0300079) q[0];
rz(-pi) q[1];
rz(0.28658861) q[2];
sx q[2];
rz(-1.5865241) q[2];
sx q[2];
rz(2.36731) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9518937) q[1];
sx q[1];
rz(-0.67283291) q[1];
sx q[1];
rz(-2.2646907) q[1];
rz(-pi) q[2];
rz(-1.8191765) q[3];
sx q[3];
rz(-2.2349149) q[3];
sx q[3];
rz(2.3770203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7297111) q[2];
sx q[2];
rz(-1.231266) q[2];
sx q[2];
rz(-2.7330772) q[2];
rz(-0.92489964) q[3];
sx q[3];
rz(-0.81245208) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
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
rz(-1.7655903) q[1];
sx q[1];
rz(-1.9648432) q[1];
sx q[1];
rz(1.2479938) q[1];
rz(2.3464936) q[2];
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
