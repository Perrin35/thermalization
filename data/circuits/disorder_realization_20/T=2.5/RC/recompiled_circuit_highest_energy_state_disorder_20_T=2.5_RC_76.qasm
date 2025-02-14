OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.72467726) q[0];
sx q[0];
rz(4.4530498) q[0];
sx q[0];
rz(11.761576) q[0];
rz(-2.5465487) q[1];
sx q[1];
rz(-0.19785985) q[1];
sx q[1];
rz(-1.2208389) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0752995) q[0];
sx q[0];
rz(-1.6907215) q[0];
sx q[0];
rz(-1.6060702) q[0];
rz(1.7096504) q[2];
sx q[2];
rz(-1.1119016) q[2];
sx q[2];
rz(-0.51990055) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2455622) q[1];
sx q[1];
rz(-2.3136825) q[1];
sx q[1];
rz(1.6309392) q[1];
rz(-2.9555129) q[3];
sx q[3];
rz(-0.64129378) q[3];
sx q[3];
rz(-2.5945714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4294942) q[2];
sx q[2];
rz(-2.6338989) q[2];
sx q[2];
rz(1.3362159) q[2];
rz(1.3974238) q[3];
sx q[3];
rz(-1.6349399) q[3];
sx q[3];
rz(1.5558745) q[3];
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
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76032388) q[0];
sx q[0];
rz(-1.5070494) q[0];
sx q[0];
rz(0.67833483) q[0];
rz(0.61596576) q[1];
sx q[1];
rz(-1.0266285) q[1];
sx q[1];
rz(-2.9982627) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45637886) q[0];
sx q[0];
rz(-2.7806578) q[0];
sx q[0];
rz(1.1120615) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7023588) q[2];
sx q[2];
rz(-0.72267294) q[2];
sx q[2];
rz(-2.6844048) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.7166087) q[1];
sx q[1];
rz(-2.0758481) q[1];
sx q[1];
rz(-1.6249955) q[1];
rz(-pi) q[2];
rz(-0.42558221) q[3];
sx q[3];
rz(-1.2493361) q[3];
sx q[3];
rz(-1.5150573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6501288) q[2];
sx q[2];
rz(-0.77892196) q[2];
sx q[2];
rz(-1.2196563) q[2];
rz(2.8236112) q[3];
sx q[3];
rz(-0.82428437) q[3];
sx q[3];
rz(-0.73941755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7431444) q[0];
sx q[0];
rz(-2.2205181) q[0];
sx q[0];
rz(0.63968023) q[0];
rz(-1.8094874) q[1];
sx q[1];
rz(-2.2001241) q[1];
sx q[1];
rz(-0.21523062) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2202323) q[0];
sx q[0];
rz(-1.7841993) q[0];
sx q[0];
rz(2.6409297) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5528312) q[2];
sx q[2];
rz(-0.28436545) q[2];
sx q[2];
rz(-0.21321061) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.135268) q[1];
sx q[1];
rz(-2.746114) q[1];
sx q[1];
rz(1.5436103) q[1];
rz(2.027719) q[3];
sx q[3];
rz(-1.6789241) q[3];
sx q[3];
rz(-2.7939452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3053863) q[2];
sx q[2];
rz(-2.6658194) q[2];
sx q[2];
rz(0.22030182) q[2];
rz(-0.78271714) q[3];
sx q[3];
rz(-1.3681151) q[3];
sx q[3];
rz(0.14557423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43133217) q[0];
sx q[0];
rz(-0.98589698) q[0];
sx q[0];
rz(1.8248935) q[0];
rz(0.45007625) q[1];
sx q[1];
rz(-1.4349667) q[1];
sx q[1];
rz(-0.11134527) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8985626) q[0];
sx q[0];
rz(-2.4353409) q[0];
sx q[0];
rz(0.39470048) q[0];
rz(-1.3100805) q[2];
sx q[2];
rz(-2.8088154) q[2];
sx q[2];
rz(2.8944601) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7584472) q[1];
sx q[1];
rz(-2.1229593) q[1];
sx q[1];
rz(-0.32917413) q[1];
rz(2.0599635) q[3];
sx q[3];
rz(-2.0738154) q[3];
sx q[3];
rz(-2.4693927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4075809) q[2];
sx q[2];
rz(-2.0110726) q[2];
sx q[2];
rz(-0.13047516) q[2];
rz(0.16034165) q[3];
sx q[3];
rz(-0.66157833) q[3];
sx q[3];
rz(0.17016889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57123667) q[0];
sx q[0];
rz(-0.2061051) q[0];
sx q[0];
rz(0.11827949) q[0];
rz(2.2453399) q[1];
sx q[1];
rz(-1.4786485) q[1];
sx q[1];
rz(0.55312696) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.129707) q[0];
sx q[0];
rz(-1.0965523) q[0];
sx q[0];
rz(-1.9408731) q[0];
rz(-2.4055355) q[2];
sx q[2];
rz(-1.9070574) q[2];
sx q[2];
rz(-2.2081809) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0748079) q[1];
sx q[1];
rz(-0.60943595) q[1];
sx q[1];
rz(1.273073) q[1];
rz(-pi) q[2];
rz(-0.40377577) q[3];
sx q[3];
rz(-2.194768) q[3];
sx q[3];
rz(1.7221019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7913738) q[2];
sx q[2];
rz(-0.72177902) q[2];
sx q[2];
rz(-1.0132033) q[2];
rz(2.8241217) q[3];
sx q[3];
rz(-1.6976796) q[3];
sx q[3];
rz(0.95399323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0233699) q[0];
sx q[0];
rz(-2.2582) q[0];
sx q[0];
rz(-0.50514847) q[0];
rz(0.39816868) q[1];
sx q[1];
rz(-1.265637) q[1];
sx q[1];
rz(-1.3291976) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5564541) q[0];
sx q[0];
rz(-2.6449361) q[0];
sx q[0];
rz(1.2206379) q[0];
rz(1.2114491) q[2];
sx q[2];
rz(-1.0654176) q[2];
sx q[2];
rz(2.0063248) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8201645) q[1];
sx q[1];
rz(-2.1315931) q[1];
sx q[1];
rz(3.012077) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8425214) q[3];
sx q[3];
rz(-1.8193805) q[3];
sx q[3];
rz(0.18928537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.087611) q[2];
sx q[2];
rz(-1.8113281) q[2];
sx q[2];
rz(0.83016738) q[2];
rz(-1.6603445) q[3];
sx q[3];
rz(-1.0426499) q[3];
sx q[3];
rz(1.164485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2723349) q[0];
sx q[0];
rz(-2.2255958) q[0];
sx q[0];
rz(1.2087615) q[0];
rz(-2.8973978) q[1];
sx q[1];
rz(-1.9701651) q[1];
sx q[1];
rz(-0.29921439) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1080782) q[0];
sx q[0];
rz(-0.21768269) q[0];
sx q[0];
rz(1.449145) q[0];
rz(2.1019548) q[2];
sx q[2];
rz(-2.1130145) q[2];
sx q[2];
rz(1.416666) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.17822325) q[1];
sx q[1];
rz(-2.2223964) q[1];
sx q[1];
rz(3.0804407) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.70280178) q[3];
sx q[3];
rz(-1.8918103) q[3];
sx q[3];
rz(-1.4218083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9485335) q[2];
sx q[2];
rz(-1.4024573) q[2];
sx q[2];
rz(-2.194727) q[2];
rz(-1.3777422) q[3];
sx q[3];
rz(-2.6008714) q[3];
sx q[3];
rz(-0.64259678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68719012) q[0];
sx q[0];
rz(-0.43455046) q[0];
sx q[0];
rz(-2.6026671) q[0];
rz(2.9680805) q[1];
sx q[1];
rz(-0.75650802) q[1];
sx q[1];
rz(0.38421806) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.980775) q[0];
sx q[0];
rz(-1.4133102) q[0];
sx q[0];
rz(2.6927267) q[0];
rz(2.539345) q[2];
sx q[2];
rz(-1.7302697) q[2];
sx q[2];
rz(0.98724706) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9335971) q[1];
sx q[1];
rz(-1.5446481) q[1];
sx q[1];
rz(-1.5118062) q[1];
x q[2];
rz(-2.124765) q[3];
sx q[3];
rz(-1.5477383) q[3];
sx q[3];
rz(0.35820828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.82181278) q[2];
sx q[2];
rz(-0.76924789) q[2];
sx q[2];
rz(-0.5963076) q[2];
rz(-2.1042306) q[3];
sx q[3];
rz(-2.3682902) q[3];
sx q[3];
rz(-1.1843225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3768815) q[0];
sx q[0];
rz(-2.3887971) q[0];
sx q[0];
rz(-2.9994614) q[0];
rz(-2.0243952) q[1];
sx q[1];
rz(-1.486472) q[1];
sx q[1];
rz(-1.9140917) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.759533) q[0];
sx q[0];
rz(-0.56175023) q[0];
sx q[0];
rz(0.53680935) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4332716) q[2];
sx q[2];
rz(-1.5548717) q[2];
sx q[2];
rz(-0.28057306) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.043316226) q[1];
sx q[1];
rz(-2.195561) q[1];
sx q[1];
rz(-0.079778133) q[1];
rz(-pi) q[2];
rz(2.7015721) q[3];
sx q[3];
rz(-2.2624863) q[3];
sx q[3];
rz(-1.3662076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0907937) q[2];
sx q[2];
rz(-2.8992081) q[2];
sx q[2];
rz(2.3432815) q[2];
rz(-1.5618886) q[3];
sx q[3];
rz(-1.3825994) q[3];
sx q[3];
rz(-0.17670512) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0853737) q[0];
sx q[0];
rz(-2.4985785) q[0];
sx q[0];
rz(0.0052848919) q[0];
rz(0.082848631) q[1];
sx q[1];
rz(-1.1033892) q[1];
sx q[1];
rz(0.26758912) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2076715) q[0];
sx q[0];
rz(-0.90822847) q[0];
sx q[0];
rz(1.9703869) q[0];
x q[1];
rz(-2.204699) q[2];
sx q[2];
rz(-0.46910252) q[2];
sx q[2];
rz(-1.3718951) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5983683) q[1];
sx q[1];
rz(-2.7562047) q[1];
sx q[1];
rz(2.1699127) q[1];
rz(-pi) q[2];
rz(-1.0571207) q[3];
sx q[3];
rz(-1.9909715) q[3];
sx q[3];
rz(-0.021180245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.573632) q[2];
sx q[2];
rz(-2.6601514) q[2];
sx q[2];
rz(-1.2130223) q[2];
rz(1.2447478) q[3];
sx q[3];
rz(-0.48143482) q[3];
sx q[3];
rz(-1.9511694) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5382814) q[0];
sx q[0];
rz(-0.9700226) q[0];
sx q[0];
rz(0.097401311) q[0];
rz(3.1311323) q[1];
sx q[1];
rz(-0.86581007) q[1];
sx q[1];
rz(2.5444358) q[1];
rz(0.14870208) q[2];
sx q[2];
rz(-0.54051334) q[2];
sx q[2];
rz(1.9220026) q[2];
rz(2.1599471) q[3];
sx q[3];
rz(-1.4245778) q[3];
sx q[3];
rz(2.1733976) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
