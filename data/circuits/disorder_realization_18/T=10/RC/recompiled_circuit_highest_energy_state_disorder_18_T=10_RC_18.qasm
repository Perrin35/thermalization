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
rz(-1.350116) q[0];
sx q[0];
rz(-2.465829) q[0];
sx q[0];
rz(0.0078553353) q[0];
rz(0.76771843) q[1];
sx q[1];
rz(1.6176728) q[1];
sx q[1];
rz(9.6674506) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4123101) q[0];
sx q[0];
rz(-0.72107102) q[0];
sx q[0];
rz(-1.3789873) q[0];
rz(-pi) q[1];
rz(1.2062293) q[2];
sx q[2];
rz(-2.0040214) q[2];
sx q[2];
rz(-1.792576) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3625177) q[1];
sx q[1];
rz(-1.9240161) q[1];
sx q[1];
rz(1.5959969) q[1];
rz(-pi) q[2];
rz(-1.8170361) q[3];
sx q[3];
rz(-0.9199577) q[3];
sx q[3];
rz(-0.090373978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0365888) q[2];
sx q[2];
rz(-1.236311) q[2];
sx q[2];
rz(-0.71199065) q[2];
rz(2.8365734) q[3];
sx q[3];
rz(-0.22235338) q[3];
sx q[3];
rz(-0.045507889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63118339) q[0];
sx q[0];
rz(-1.864186) q[0];
sx q[0];
rz(-1.7741868) q[0];
rz(0.039904682) q[1];
sx q[1];
rz(-0.80723643) q[1];
sx q[1];
rz(2.5423999) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73030534) q[0];
sx q[0];
rz(-1.4559901) q[0];
sx q[0];
rz(0.47700096) q[0];
rz(-2.3223898) q[2];
sx q[2];
rz(-1.1408198) q[2];
sx q[2];
rz(-2.5366304) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2476462) q[1];
sx q[1];
rz(-2.1987913) q[1];
sx q[1];
rz(1.1537854) q[1];
rz(2.4559095) q[3];
sx q[3];
rz(-1.7910379) q[3];
sx q[3];
rz(1.1789279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0590608) q[2];
sx q[2];
rz(-2.579687) q[2];
sx q[2];
rz(-1.3085636) q[2];
rz(0.023905309) q[3];
sx q[3];
rz(-1.5289565) q[3];
sx q[3];
rz(1.5628975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6858653) q[0];
sx q[0];
rz(-2.4754334) q[0];
sx q[0];
rz(0.84079963) q[0];
rz(2.1127545) q[1];
sx q[1];
rz(-2.7327171) q[1];
sx q[1];
rz(-1.6811446) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.235229) q[0];
sx q[0];
rz(-2.0301986) q[0];
sx q[0];
rz(1.4109341) q[0];
x q[1];
rz(0.65030789) q[2];
sx q[2];
rz(-2.569649) q[2];
sx q[2];
rz(0.24947333) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.19059316) q[1];
sx q[1];
rz(-1.237932) q[1];
sx q[1];
rz(-0.69575633) q[1];
rz(-pi) q[2];
rz(1.5140686) q[3];
sx q[3];
rz(-2.1773844) q[3];
sx q[3];
rz(1.1067672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8847522) q[2];
sx q[2];
rz(-1.931087) q[2];
sx q[2];
rz(0.15667008) q[2];
rz(-1.6414075) q[3];
sx q[3];
rz(-1.898396) q[3];
sx q[3];
rz(2.959804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0750065) q[0];
sx q[0];
rz(-1.8445419) q[0];
sx q[0];
rz(-0.080667607) q[0];
rz(-2.0196041) q[1];
sx q[1];
rz(-0.64067084) q[1];
sx q[1];
rz(1.5544308) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98469668) q[0];
sx q[0];
rz(-1.6365572) q[0];
sx q[0];
rz(0.82729152) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.3360668) q[2];
sx q[2];
rz(-0.12305752) q[2];
sx q[2];
rz(2.6390136) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.746628) q[1];
sx q[1];
rz(-1.0735128) q[1];
sx q[1];
rz(-2.4056466) q[1];
rz(-pi) q[2];
rz(-0.60685632) q[3];
sx q[3];
rz(-1.8491866) q[3];
sx q[3];
rz(-2.3902219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.88177219) q[2];
sx q[2];
rz(-2.0505003) q[2];
sx q[2];
rz(2.1176977) q[2];
rz(-2.3648868) q[3];
sx q[3];
rz(-1.9909765) q[3];
sx q[3];
rz(2.1385433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40989947) q[0];
sx q[0];
rz(-0.85407805) q[0];
sx q[0];
rz(-0.62752974) q[0];
rz(-0.69951406) q[1];
sx q[1];
rz(-2.1414089) q[1];
sx q[1];
rz(-0.90739179) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32226792) q[0];
sx q[0];
rz(-1.6108247) q[0];
sx q[0];
rz(2.3884474) q[0];
rz(-2.444958) q[2];
sx q[2];
rz(-0.24605909) q[2];
sx q[2];
rz(1.5738894) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.34037922) q[1];
sx q[1];
rz(-2.0836909) q[1];
sx q[1];
rz(2.7199634) q[1];
rz(0.058919546) q[3];
sx q[3];
rz(-1.0353966) q[3];
sx q[3];
rz(2.8297092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0044535) q[2];
sx q[2];
rz(-1.8567825) q[2];
sx q[2];
rz(-2.2966906) q[2];
rz(0.6984624) q[3];
sx q[3];
rz(-0.74028492) q[3];
sx q[3];
rz(2.3499878) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38248211) q[0];
sx q[0];
rz(-1.4610721) q[0];
sx q[0];
rz(-1.2680898) q[0];
rz(-1.6146487) q[1];
sx q[1];
rz(-2.0497649) q[1];
sx q[1];
rz(2.7472034) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1295107) q[0];
sx q[0];
rz(-0.92716588) q[0];
sx q[0];
rz(0.056931007) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4066448) q[2];
sx q[2];
rz(-1.5089261) q[2];
sx q[2];
rz(2.8501373) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0002546) q[1];
sx q[1];
rz(-0.91857498) q[1];
sx q[1];
rz(-1.5144996) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9378618) q[3];
sx q[3];
rz(-2.7739848) q[3];
sx q[3];
rz(2.550761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7093198) q[2];
sx q[2];
rz(-3.0797112) q[2];
sx q[2];
rz(-0.56149948) q[2];
rz(1.1310486) q[3];
sx q[3];
rz(-0.80315042) q[3];
sx q[3];
rz(-0.48857442) q[3];
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
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5973709) q[0];
sx q[0];
rz(-0.18246305) q[0];
sx q[0];
rz(0.23319787) q[0];
rz(3.0274262) q[1];
sx q[1];
rz(-2.3515067) q[1];
sx q[1];
rz(2.9679969) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1602992) q[0];
sx q[0];
rz(-0.7795142) q[0];
sx q[0];
rz(1.0407838) q[0];
x q[1];
rz(-2.0162705) q[2];
sx q[2];
rz(-1.2230754) q[2];
sx q[2];
rz(1.7537774) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1039338) q[1];
sx q[1];
rz(-1.3500431) q[1];
sx q[1];
rz(0.16640618) q[1];
rz(-pi) q[2];
rz(1.5585445) q[3];
sx q[3];
rz(-1.3821332) q[3];
sx q[3];
rz(2.0607299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0357828) q[2];
sx q[2];
rz(-2.1828987) q[2];
sx q[2];
rz(-2.2868273) q[2];
rz(-0.0828951) q[3];
sx q[3];
rz(-1.9708743) q[3];
sx q[3];
rz(-0.054072592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4891124) q[0];
sx q[0];
rz(-1.2615477) q[0];
sx q[0];
rz(-1.9626807) q[0];
rz(1.0014125) q[1];
sx q[1];
rz(-1.8779571) q[1];
sx q[1];
rz(-2.9875535) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1106873) q[0];
sx q[0];
rz(-0.89546916) q[0];
sx q[0];
rz(1.9934325) q[0];
x q[1];
rz(0.87100864) q[2];
sx q[2];
rz(-2.1453875) q[2];
sx q[2];
rz(0.99829295) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2630059) q[1];
sx q[1];
rz(-1.1060426) q[1];
sx q[1];
rz(2.8057363) q[1];
x q[2];
rz(-2.4064111) q[3];
sx q[3];
rz(-0.38215853) q[3];
sx q[3];
rz(-0.89565403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4607294) q[2];
sx q[2];
rz(-2.0622084) q[2];
sx q[2];
rz(-0.66780773) q[2];
rz(1.7364511) q[3];
sx q[3];
rz(-1.3677771) q[3];
sx q[3];
rz(2.5051129) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3619096) q[0];
sx q[0];
rz(-1.4754262) q[0];
sx q[0];
rz(-0.59378004) q[0];
rz(-1.8608015) q[1];
sx q[1];
rz(-2.180876) q[1];
sx q[1];
rz(1.8437754) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.058859874) q[0];
sx q[0];
rz(-1.642859) q[0];
sx q[0];
rz(2.5565992) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.4164575) q[2];
sx q[2];
rz(-0.67085941) q[2];
sx q[2];
rz(-0.79029467) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1480963) q[1];
sx q[1];
rz(-1.2132753) q[1];
sx q[1];
rz(-1.7054249) q[1];
rz(-1.7421977) q[3];
sx q[3];
rz(-0.77091445) q[3];
sx q[3];
rz(-0.41219974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7456776) q[2];
sx q[2];
rz(-3.119645) q[2];
sx q[2];
rz(1.5094666) q[2];
rz(-0.11219003) q[3];
sx q[3];
rz(-1.0271007) q[3];
sx q[3];
rz(-1.3794911) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6701732) q[0];
sx q[0];
rz(-1.0059953) q[0];
sx q[0];
rz(0.26563409) q[0];
rz(-0.98948014) q[1];
sx q[1];
rz(-1.2629291) q[1];
sx q[1];
rz(-2.8094453) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6393747) q[0];
sx q[0];
rz(-1.9042339) q[0];
sx q[0];
rz(3.1375454) q[0];
x q[1];
rz(2.1493456) q[2];
sx q[2];
rz(-1.8162811) q[2];
sx q[2];
rz(-0.13260376) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.89236081) q[1];
sx q[1];
rz(-2.7977825) q[1];
sx q[1];
rz(-1.7147786) q[1];
rz(-pi) q[2];
rz(-1.321855) q[3];
sx q[3];
rz(-2.948108) q[3];
sx q[3];
rz(-1.8850808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0418479) q[2];
sx q[2];
rz(-2.0691278) q[2];
sx q[2];
rz(0.15038807) q[2];
rz(0.8485052) q[3];
sx q[3];
rz(-0.34981194) q[3];
sx q[3];
rz(-1.8457886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083241845) q[0];
sx q[0];
rz(-1.8235089) q[0];
sx q[0];
rz(1.9705809) q[0];
rz(1.8104443) q[1];
sx q[1];
rz(-2.34453) q[1];
sx q[1];
rz(-1.1042368) q[1];
rz(0.83126478) q[2];
sx q[2];
rz(-1.5223885) q[2];
sx q[2];
rz(1.6123966) q[2];
rz(-2.4348197) q[3];
sx q[3];
rz(-1.5512244) q[3];
sx q[3];
rz(2.649818) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
