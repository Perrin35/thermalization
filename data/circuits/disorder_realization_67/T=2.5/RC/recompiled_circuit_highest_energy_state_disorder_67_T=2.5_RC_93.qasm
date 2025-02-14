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
rz(1.2713852) q[0];
sx q[0];
rz(-0.013590824) q[0];
sx q[0];
rz(3.115227) q[0];
rz(0.68323505) q[1];
sx q[1];
rz(4.9024138) q[1];
sx q[1];
rz(9.496357) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.723322) q[0];
sx q[0];
rz(-1.5958565) q[0];
sx q[0];
rz(1.5610831) q[0];
rz(-pi) q[1];
rz(-1.3355005) q[2];
sx q[2];
rz(-1.9517731) q[2];
sx q[2];
rz(1.345696) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.84715473) q[1];
sx q[1];
rz(-1.5658453) q[1];
sx q[1];
rz(-1.5895542) q[1];
rz(0.79238331) q[3];
sx q[3];
rz(-0.43607831) q[3];
sx q[3];
rz(2.53873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2804395) q[2];
sx q[2];
rz(-0.70715487) q[2];
sx q[2];
rz(0.46671483) q[2];
rz(0.47433445) q[3];
sx q[3];
rz(-3.1204087) q[3];
sx q[3];
rz(-0.02136136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83653432) q[0];
sx q[0];
rz(-0.49764043) q[0];
sx q[0];
rz(3.1217788) q[0];
rz(1.548798) q[1];
sx q[1];
rz(-0.21983799) q[1];
sx q[1];
rz(-1.4733431) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23035717) q[0];
sx q[0];
rz(-0.94382554) q[0];
sx q[0];
rz(0.66642739) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2884533) q[2];
sx q[2];
rz(-0.38333508) q[2];
sx q[2];
rz(2.9789973) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.87621385) q[1];
sx q[1];
rz(-0.20890954) q[1];
sx q[1];
rz(1.9494328) q[1];
rz(-pi) q[2];
x q[2];
rz(0.83062828) q[3];
sx q[3];
rz(-2.1601956) q[3];
sx q[3];
rz(-2.7089861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.3026498) q[2];
sx q[2];
rz(-2.5431716) q[2];
sx q[2];
rz(-1.8518651) q[2];
rz(1.9047811) q[3];
sx q[3];
rz(-2.808282) q[3];
sx q[3];
rz(2.439177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97238338) q[0];
sx q[0];
rz(-1.9820259) q[0];
sx q[0];
rz(1.5706536) q[0];
rz(1.4717357) q[1];
sx q[1];
rz(-1.6153299) q[1];
sx q[1];
rz(-2.6872046) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1412068) q[0];
sx q[0];
rz(-2.1659746) q[0];
sx q[0];
rz(2.9808874) q[0];
rz(-pi) q[1];
rz(1.7128031) q[2];
sx q[2];
rz(-1.596984) q[2];
sx q[2];
rz(2.533503) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.66167605) q[1];
sx q[1];
rz(-3.0212086) q[1];
sx q[1];
rz(1.7604339) q[1];
x q[2];
rz(-2.4942355) q[3];
sx q[3];
rz(-2.1212) q[3];
sx q[3];
rz(3.0542706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.74139524) q[2];
sx q[2];
rz(-3.1252842) q[2];
sx q[2];
rz(0.34747094) q[2];
rz(-2.6853284) q[3];
sx q[3];
rz(-0.01472344) q[3];
sx q[3];
rz(2.1485476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4149813) q[0];
sx q[0];
rz(-1.240629) q[0];
sx q[0];
rz(1.7353143) q[0];
rz(-2.6982488) q[1];
sx q[1];
rz(-1.025238) q[1];
sx q[1];
rz(1.5700856) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9260028) q[0];
sx q[0];
rz(-2.2099751) q[0];
sx q[0];
rz(-1.502468) q[0];
rz(-pi) q[1];
rz(2.3894044) q[2];
sx q[2];
rz(-0.094399422) q[2];
sx q[2];
rz(-2.029325) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0896596) q[1];
sx q[1];
rz(-1.5610236) q[1];
sx q[1];
rz(-1.8490514) q[1];
rz(-1.0138233) q[3];
sx q[3];
rz(-2.644745) q[3];
sx q[3];
rz(2.5574977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7683679) q[2];
sx q[2];
rz(-2.7478605) q[2];
sx q[2];
rz(3.0623398) q[2];
rz(-1.1894038) q[3];
sx q[3];
rz(-1.7704084) q[3];
sx q[3];
rz(-1.6325379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(0.70334148) q[0];
sx q[0];
rz(-0.51333135) q[0];
sx q[0];
rz(-0.80605036) q[0];
rz(-2.3015859) q[1];
sx q[1];
rz(-3.1286616) q[1];
sx q[1];
rz(-2.3654225) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0426725) q[0];
sx q[0];
rz(-1.4399043) q[0];
sx q[0];
rz(2.5517273) q[0];
rz(-3.0018535) q[2];
sx q[2];
rz(-0.010541803) q[2];
sx q[2];
rz(3.000562) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8766206) q[1];
sx q[1];
rz(-1.435346) q[1];
sx q[1];
rz(1.5587864) q[1];
rz(-pi) q[2];
x q[2];
rz(0.25960323) q[3];
sx q[3];
rz(-1.6791376) q[3];
sx q[3];
rz(-1.3110127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1157896) q[2];
sx q[2];
rz(-1.5805406) q[2];
sx q[2];
rz(2.4157794) q[2];
rz(-0.18802655) q[3];
sx q[3];
rz(-0.060421061) q[3];
sx q[3];
rz(-0.70992011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1793154) q[0];
sx q[0];
rz(-0.55914068) q[0];
sx q[0];
rz(2.6002) q[0];
rz(0.18556449) q[1];
sx q[1];
rz(-1.5927529) q[1];
sx q[1];
rz(-3.013179) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5914766) q[0];
sx q[0];
rz(-1.8719561) q[0];
sx q[0];
rz(2.7635283) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1368581) q[2];
sx q[2];
rz(-1.4494697) q[2];
sx q[2];
rz(-1.5743299) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7306108) q[1];
sx q[1];
rz(-1.1514947) q[1];
sx q[1];
rz(-1.6602519) q[1];
rz(-pi) q[2];
rz(1.7592247) q[3];
sx q[3];
rz(-1.5604094) q[3];
sx q[3];
rz(-2.5760108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.386117) q[2];
sx q[2];
rz(-0.057689276) q[2];
sx q[2];
rz(2.3003787) q[2];
rz(-0.20127131) q[3];
sx q[3];
rz(-1.5320675) q[3];
sx q[3];
rz(0.25963983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5801308) q[0];
sx q[0];
rz(-0.78829563) q[0];
sx q[0];
rz(1.5607675) q[0];
rz(-2.4687817) q[1];
sx q[1];
rz(-1.4038059) q[1];
sx q[1];
rz(-0.028884551) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2530841) q[0];
sx q[0];
rz(-2.7612503) q[0];
sx q[0];
rz(-2.4230291) q[0];
x q[1];
rz(1.3703652) q[2];
sx q[2];
rz(-1.1005963) q[2];
sx q[2];
rz(1.7261795) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1950732) q[1];
sx q[1];
rz(-1.9819489) q[1];
sx q[1];
rz(2.8821936) q[1];
x q[2];
rz(-0.77856346) q[3];
sx q[3];
rz(-1.9664008) q[3];
sx q[3];
rz(-2.9958519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2334571) q[2];
sx q[2];
rz(-0.59671777) q[2];
sx q[2];
rz(0.92823589) q[2];
rz(-0.2975896) q[3];
sx q[3];
rz(-0.15407763) q[3];
sx q[3];
rz(2.2969864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.069227844) q[0];
sx q[0];
rz(-2.897825) q[0];
sx q[0];
rz(0.099076554) q[0];
rz(0.98723269) q[1];
sx q[1];
rz(-1.8376553) q[1];
sx q[1];
rz(-2.6027021) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9469556) q[0];
sx q[0];
rz(-1.4879003) q[0];
sx q[0];
rz(0.006587365) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6201934) q[2];
sx q[2];
rz(-1.5954752) q[2];
sx q[2];
rz(0.53804735) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3618631) q[1];
sx q[1];
rz(-1.2113844) q[1];
sx q[1];
rz(-0.56198175) q[1];
rz(-pi) q[2];
rz(1.8935558) q[3];
sx q[3];
rz(-0.071167067) q[3];
sx q[3];
rz(0.038480345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.99523669) q[2];
sx q[2];
rz(-0.015492798) q[2];
sx q[2];
rz(-0.30919477) q[2];
rz(2.501798) q[3];
sx q[3];
rz(-3.1412509) q[3];
sx q[3];
rz(0.063808002) q[3];
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
rz(-2.8362506) q[0];
sx q[0];
rz(-2.5449365) q[0];
sx q[0];
rz(3.0869361) q[0];
rz(-1.1326185) q[1];
sx q[1];
rz(-1.9839958) q[1];
sx q[1];
rz(-1.2861015) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4187546) q[0];
sx q[0];
rz(-1.6143454) q[0];
sx q[0];
rz(-1.7489793) q[0];
x q[1];
rz(-2.340592) q[2];
sx q[2];
rz(-0.05861662) q[2];
sx q[2];
rz(2.336077) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.4234679) q[1];
sx q[1];
rz(-1.3455952) q[1];
sx q[1];
rz(1.6632715) q[1];
rz(0.07225424) q[3];
sx q[3];
rz(-2.0569306) q[3];
sx q[3];
rz(2.2421601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9201811) q[2];
sx q[2];
rz(-2.5765918) q[2];
sx q[2];
rz(-1.1037702) q[2];
rz(1.6086027) q[3];
sx q[3];
rz(-3.0975603) q[3];
sx q[3];
rz(-0.65582961) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1405545) q[0];
sx q[0];
rz(-0.1796722) q[0];
sx q[0];
rz(-0.0044862577) q[0];
rz(-1.5890315) q[1];
sx q[1];
rz(-1.4493425) q[1];
sx q[1];
rz(0.057131279) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82027869) q[0];
sx q[0];
rz(-1.4749267) q[0];
sx q[0];
rz(-1.3906327) q[0];
rz(1.42679) q[2];
sx q[2];
rz(-1.7911573) q[2];
sx q[2];
rz(1.2978467) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2559214) q[1];
sx q[1];
rz(-0.90270611) q[1];
sx q[1];
rz(2.9873965) q[1];
x q[2];
rz(0.83052333) q[3];
sx q[3];
rz(-1.4581513) q[3];
sx q[3];
rz(-0.5121246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.39157465) q[2];
sx q[2];
rz(-3.1142758) q[2];
sx q[2];
rz(-2.2726783) q[2];
rz(-0.9817552) q[3];
sx q[3];
rz(-3.1120286) q[3];
sx q[3];
rz(-2.6386007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0960196) q[0];
sx q[0];
rz(-1.6572784) q[0];
sx q[0];
rz(-1.4838765) q[0];
rz(-0.45687301) q[1];
sx q[1];
rz(-0.15468205) q[1];
sx q[1];
rz(-0.044943132) q[1];
rz(-0.1931242) q[2];
sx q[2];
rz(-2.3681691) q[2];
sx q[2];
rz(-2.9409627) q[2];
rz(0.18263541) q[3];
sx q[3];
rz(-0.67125139) q[3];
sx q[3];
rz(-2.8677979) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
