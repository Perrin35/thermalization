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
rz(3.1280018) q[0];
sx q[0];
rz(9.4511436) q[0];
rz(-2.4583576) q[1];
sx q[1];
rz(-1.7608211) q[1];
sx q[1];
rz(3.0700136) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3535093) q[0];
sx q[0];
rz(-3.1147163) q[0];
sx q[0];
rz(-2.7719017) q[0];
rz(0.52710345) q[2];
sx q[2];
rz(-2.6968535) q[2];
sx q[2];
rz(-1.9183733) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.72373448) q[1];
sx q[1];
rz(-1.5520387) q[1];
sx q[1];
rz(-3.1366407) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.79238331) q[3];
sx q[3];
rz(-2.7055143) q[3];
sx q[3];
rz(2.53873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8611531) q[2];
sx q[2];
rz(-0.70715487) q[2];
sx q[2];
rz(0.46671483) q[2];
rz(-2.6672582) q[3];
sx q[3];
rz(-0.021183906) q[3];
sx q[3];
rz(0.02136136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3050583) q[0];
sx q[0];
rz(-2.6439522) q[0];
sx q[0];
rz(3.1217788) q[0];
rz(1.5927947) q[1];
sx q[1];
rz(-0.21983799) q[1];
sx q[1];
rz(1.4733431) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1602185) q[0];
sx q[0];
rz(-2.2607973) q[0];
sx q[0];
rz(0.86440683) q[0];
rz(1.8657487) q[2];
sx q[2];
rz(-1.3222857) q[2];
sx q[2];
rz(-2.4139443) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.32343601) q[1];
sx q[1];
rz(-1.4940573) q[1];
sx q[1];
rz(1.7652926) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3519312) q[3];
sx q[3];
rz(-2.2316311) q[3];
sx q[3];
rz(2.5500962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.3026498) q[2];
sx q[2];
rz(-2.5431716) q[2];
sx q[2];
rz(1.2897276) q[2];
rz(1.9047811) q[3];
sx q[3];
rz(-0.33331063) q[3];
sx q[3];
rz(-2.439177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97238338) q[0];
sx q[0];
rz(-1.9820259) q[0];
sx q[0];
rz(-1.5709391) q[0];
rz(-1.4717357) q[1];
sx q[1];
rz(-1.5262628) q[1];
sx q[1];
rz(0.45438802) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7189302) q[0];
sx q[0];
rz(-0.61394982) q[0];
sx q[0];
rz(1.3387247) q[0];
rz(-pi) q[1];
rz(1.387792) q[2];
sx q[2];
rz(-2.9972074) q[2];
sx q[2];
rz(2.3600277) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4799166) q[1];
sx q[1];
rz(-0.12038409) q[1];
sx q[1];
rz(-1.3811587) q[1];
x q[2];
rz(2.347499) q[3];
sx q[3];
rz(-0.82320582) q[3];
sx q[3];
rz(-2.2633207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4001974) q[2];
sx q[2];
rz(-0.016308451) q[2];
sx q[2];
rz(-2.7941217) q[2];
rz(0.45626429) q[3];
sx q[3];
rz(-0.01472344) q[3];
sx q[3];
rz(-0.99304503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4149813) q[0];
sx q[0];
rz(-1.9009637) q[0];
sx q[0];
rz(1.4062784) q[0];
rz(-0.44334385) q[1];
sx q[1];
rz(-2.1163546) q[1];
sx q[1];
rz(1.5700856) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8271874) q[0];
sx q[0];
rz(-1.5159722) q[0];
sx q[0];
rz(2.501295) q[0];
rz(3.072567) q[2];
sx q[2];
rz(-1.6352425) q[2];
sx q[2];
rz(2.850159) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6199377) q[1];
sx q[1];
rz(-1.2925549) q[1];
sx q[1];
rz(-0.010163608) q[1];
rz(-pi) q[2];
rz(0.27914234) q[3];
sx q[3];
rz(-1.1542392) q[3];
sx q[3];
rz(0.032241658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7683679) q[2];
sx q[2];
rz(-2.7478605) q[2];
sx q[2];
rz(-3.0623398) q[2];
rz(-1.1894038) q[3];
sx q[3];
rz(-1.7704084) q[3];
sx q[3];
rz(-1.6325379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70334148) q[0];
sx q[0];
rz(-2.6282613) q[0];
sx q[0];
rz(-2.3355423) q[0];
rz(-2.3015859) q[1];
sx q[1];
rz(-0.012931074) q[1];
sx q[1];
rz(2.3654225) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7568593) q[0];
sx q[0];
rz(-0.98664588) q[0];
sx q[0];
rz(1.7279051) q[0];
rz(-pi) q[1];
rz(1.5722647) q[2];
sx q[2];
rz(-1.5603573) q[2];
sx q[2];
rz(-2.8608152) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8766206) q[1];
sx q[1];
rz(-1.435346) q[1];
sx q[1];
rz(1.5587864) q[1];
rz(-pi) q[2];
rz(0.40078409) q[3];
sx q[3];
rz(-0.28083193) q[3];
sx q[3];
rz(0.1268445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.025803056) q[2];
sx q[2];
rz(-1.5805406) q[2];
sx q[2];
rz(-2.4157794) q[2];
rz(-2.9535661) q[3];
sx q[3];
rz(-3.0811716) q[3];
sx q[3];
rz(2.4316725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1793154) q[0];
sx q[0];
rz(-2.582452) q[0];
sx q[0];
rz(2.6002) q[0];
rz(-0.18556449) q[1];
sx q[1];
rz(-1.5927529) q[1];
sx q[1];
rz(-0.12841368) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0036412) q[0];
sx q[0];
rz(-1.9310474) q[0];
sx q[0];
rz(-1.8933348) q[0];
rz(-pi) q[1];
rz(1.5319848) q[2];
sx q[2];
rz(-0.12141849) q[2];
sx q[2];
rz(-1.613429) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0182788) q[1];
sx q[1];
rz(-1.4891081) q[1];
sx q[1];
rz(0.42079349) q[1];
rz(1.382368) q[3];
sx q[3];
rz(-1.5604094) q[3];
sx q[3];
rz(2.5760108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7554756) q[2];
sx q[2];
rz(-3.0839034) q[2];
sx q[2];
rz(2.3003787) q[2];
rz(-2.9403213) q[3];
sx q[3];
rz(-1.5320675) q[3];
sx q[3];
rz(2.8819528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5801308) q[0];
sx q[0];
rz(-2.353297) q[0];
sx q[0];
rz(1.5607675) q[0];
rz(-2.4687817) q[1];
sx q[1];
rz(-1.7377868) q[1];
sx q[1];
rz(-3.1127081) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3643234) q[0];
sx q[0];
rz(-1.3239081) q[0];
sx q[0];
rz(0.29233513) q[0];
x q[1];
rz(-2.6631764) q[2];
sx q[2];
rz(-1.3923651) q[2];
sx q[2];
rz(-2.894424) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.72994443) q[1];
sx q[1];
rz(-1.8081417) q[1];
sx q[1];
rz(-1.9945595) q[1];
x q[2];
rz(2.1013124) q[3];
sx q[3];
rz(-0.86581745) q[3];
sx q[3];
rz(-1.353273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9081356) q[2];
sx q[2];
rz(-0.59671777) q[2];
sx q[2];
rz(-0.92823589) q[2];
rz(-2.8440031) q[3];
sx q[3];
rz(-0.15407763) q[3];
sx q[3];
rz(-2.2969864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0723648) q[0];
sx q[0];
rz(-0.2437676) q[0];
sx q[0];
rz(3.0425161) q[0];
rz(0.98723269) q[1];
sx q[1];
rz(-1.3039373) q[1];
sx q[1];
rz(2.6027021) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27402747) q[0];
sx q[0];
rz(-3.0584359) q[0];
sx q[0];
rz(-1.4916791) q[0];
rz(-pi) q[1];
rz(-3.1168836) q[2];
sx q[2];
rz(-1.5214143) q[2];
sx q[2];
rz(-2.1076237) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7797295) q[1];
sx q[1];
rz(-1.9302082) q[1];
sx q[1];
rz(-0.56198175) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5032926) q[3];
sx q[3];
rz(-1.5933523) q[3];
sx q[3];
rz(1.2872788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.146356) q[2];
sx q[2];
rz(-3.1260999) q[2];
sx q[2];
rz(0.30919477) q[2];
rz(2.501798) q[3];
sx q[3];
rz(-3.1412509) q[3];
sx q[3];
rz(0.063808002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
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
rz(-1.1575969) q[1];
sx q[1];
rz(1.2861015) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.085145935) q[0];
sx q[0];
rz(-2.9582199) q[0];
sx q[0];
rz(-1.3297179) q[0];
rz(-pi) q[1];
rz(-0.80100061) q[2];
sx q[2];
rz(-3.082976) q[2];
sx q[2];
rz(2.336077) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.029833492) q[1];
sx q[1];
rz(-0.24315029) q[1];
sx q[1];
rz(-0.38317005) q[1];
rz(-pi) q[2];
rz(1.083582) q[3];
sx q[3];
rz(-1.5069252) q[3];
sx q[3];
rz(2.4364249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9201811) q[2];
sx q[2];
rz(-2.5765918) q[2];
sx q[2];
rz(1.1037702) q[2];
rz(-1.6086027) q[3];
sx q[3];
rz(-3.0975603) q[3];
sx q[3];
rz(0.65582961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.00103818) q[0];
sx q[0];
rz(-0.1796722) q[0];
sx q[0];
rz(0.0044862577) q[0];
rz(-1.5525612) q[1];
sx q[1];
rz(-1.4493425) q[1];
sx q[1];
rz(3.0844614) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9071974) q[0];
sx q[0];
rz(-0.20383862) q[0];
sx q[0];
rz(2.06334) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9190019) q[2];
sx q[2];
rz(-1.711297) q[2];
sx q[2];
rz(2.8369571) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0100952) q[1];
sx q[1];
rz(-2.4586077) q[1];
sx q[1];
rz(-1.7630153) q[1];
rz(-pi) q[2];
rz(1.404625) q[3];
sx q[3];
rz(-0.7471841) q[3];
sx q[3];
rz(-2.2053444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.750018) q[2];
sx q[2];
rz(-3.1142758) q[2];
sx q[2];
rz(2.2726783) q[2];
rz(2.1598375) q[3];
sx q[3];
rz(-0.029564094) q[3];
sx q[3];
rz(-0.50299197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.045573087) q[0];
sx q[0];
rz(-1.4843142) q[0];
sx q[0];
rz(1.6577161) q[0];
rz(0.45687301) q[1];
sx q[1];
rz(-2.9869106) q[1];
sx q[1];
rz(3.0966495) q[1];
rz(1.7560319) q[2];
sx q[2];
rz(-2.3262263) q[2];
sx q[2];
rz(-0.066166886) q[2];
rz(-0.18263541) q[3];
sx q[3];
rz(-2.4703413) q[3];
sx q[3];
rz(0.27379476) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
