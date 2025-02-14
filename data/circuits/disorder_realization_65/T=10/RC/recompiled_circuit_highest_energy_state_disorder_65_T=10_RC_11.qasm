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
rz(0.49864545) q[0];
sx q[0];
rz(3.6874229) q[0];
sx q[0];
rz(9.3064718) q[0];
rz(0.91953295) q[1];
sx q[1];
rz(4.4448648) q[1];
sx q[1];
rz(8.2104609) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81033731) q[0];
sx q[0];
rz(-2.6542568) q[0];
sx q[0];
rz(-1.1054071) q[0];
rz(-pi) q[1];
x q[1];
rz(0.9779529) q[2];
sx q[2];
rz(-2.1101769) q[2];
sx q[2];
rz(0.025503615) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2976297) q[1];
sx q[1];
rz(-1.4342035) q[1];
sx q[1];
rz(0.94874391) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.494058) q[3];
sx q[3];
rz(-2.8676448) q[3];
sx q[3];
rz(-0.63989598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.168657) q[2];
sx q[2];
rz(-1.1309705) q[2];
sx q[2];
rz(-2.0449779) q[2];
rz(-2.9003669) q[3];
sx q[3];
rz(-1.3696407) q[3];
sx q[3];
rz(-0.84478861) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.089652561) q[0];
sx q[0];
rz(-1.0900499) q[0];
sx q[0];
rz(0.3845149) q[0];
rz(0.3081201) q[1];
sx q[1];
rz(-0.48079753) q[1];
sx q[1];
rz(-0.62517977) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.168361) q[0];
sx q[0];
rz(-0.28795469) q[0];
sx q[0];
rz(0.11174996) q[0];
rz(1.4433631) q[2];
sx q[2];
rz(-0.857876) q[2];
sx q[2];
rz(2.8486203) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3846085) q[1];
sx q[1];
rz(-0.9238657) q[1];
sx q[1];
rz(1.4642843) q[1];
rz(-2.2143638) q[3];
sx q[3];
rz(-0.72690287) q[3];
sx q[3];
rz(0.34360269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1326617) q[2];
sx q[2];
rz(-0.93175685) q[2];
sx q[2];
rz(2.8399732) q[2];
rz(-2.610142) q[3];
sx q[3];
rz(-2.2558432) q[3];
sx q[3];
rz(-1.1227932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6279491) q[0];
sx q[0];
rz(-1.0114089) q[0];
sx q[0];
rz(-1.9568141) q[0];
rz(-0.14492497) q[1];
sx q[1];
rz(-1.2624319) q[1];
sx q[1];
rz(-1.5473993) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74565174) q[0];
sx q[0];
rz(-0.81123039) q[0];
sx q[0];
rz(-0.70180362) q[0];
rz(-pi) q[1];
rz(0.1833488) q[2];
sx q[2];
rz(-1.9446704) q[2];
sx q[2];
rz(-3.0484859) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5963187) q[1];
sx q[1];
rz(-0.68877586) q[1];
sx q[1];
rz(0.94249814) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2482485) q[3];
sx q[3];
rz(-1.9574021) q[3];
sx q[3];
rz(0.69535386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.80921119) q[2];
sx q[2];
rz(-0.96497649) q[2];
sx q[2];
rz(1.8154209) q[2];
rz(-0.91340804) q[3];
sx q[3];
rz(-1.923442) q[3];
sx q[3];
rz(-2.6083045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7186385) q[0];
sx q[0];
rz(-0.47684968) q[0];
sx q[0];
rz(1.780321) q[0];
rz(2.4286229) q[1];
sx q[1];
rz(-2.0013335) q[1];
sx q[1];
rz(-2.4580809) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8806802) q[0];
sx q[0];
rz(-0.99755462) q[0];
sx q[0];
rz(-3.1295958) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7943125) q[2];
sx q[2];
rz(-1.3663251) q[2];
sx q[2];
rz(2.7140537) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5826368) q[1];
sx q[1];
rz(-0.92933944) q[1];
sx q[1];
rz(-2.2333686) q[1];
x q[2];
rz(1.470742) q[3];
sx q[3];
rz(-1.3276427) q[3];
sx q[3];
rz(0.2167165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6354562) q[2];
sx q[2];
rz(-2.7244302) q[2];
sx q[2];
rz(0.48640856) q[2];
rz(0.5168612) q[3];
sx q[3];
rz(-1.1812527) q[3];
sx q[3];
rz(-1.0335056) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.085623398) q[0];
sx q[0];
rz(-1.509868) q[0];
sx q[0];
rz(-1.7686718) q[0];
rz(1.7347451) q[1];
sx q[1];
rz(-0.56766784) q[1];
sx q[1];
rz(0.47703201) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1005362) q[0];
sx q[0];
rz(-1.1541799) q[0];
sx q[0];
rz(-2.9265397) q[0];
rz(-pi) q[1];
rz(2.1978756) q[2];
sx q[2];
rz(-2.5260128) q[2];
sx q[2];
rz(-2.1114388) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7030695) q[1];
sx q[1];
rz(-1.0756936) q[1];
sx q[1];
rz(0.34214051) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.72735765) q[3];
sx q[3];
rz(-0.92669707) q[3];
sx q[3];
rz(-1.3535045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.10739747) q[2];
sx q[2];
rz(-0.98661462) q[2];
sx q[2];
rz(2.6743215) q[2];
rz(0.16921903) q[3];
sx q[3];
rz(-1.4110112) q[3];
sx q[3];
rz(2.8362078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53448236) q[0];
sx q[0];
rz(-2.2619673) q[0];
sx q[0];
rz(-0.18950263) q[0];
rz(-2.8684008) q[1];
sx q[1];
rz(-2.0676282) q[1];
sx q[1];
rz(-2.3409519) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1976182) q[0];
sx q[0];
rz(-1.5276143) q[0];
sx q[0];
rz(-1.8042685) q[0];
rz(-1.4802981) q[2];
sx q[2];
rz(-1.3234089) q[2];
sx q[2];
rz(0.43898459) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.26594083) q[1];
sx q[1];
rz(-1.7663301) q[1];
sx q[1];
rz(0.13130782) q[1];
rz(-pi) q[2];
rz(2.6012035) q[3];
sx q[3];
rz(-2.0107993) q[3];
sx q[3];
rz(-2.7203181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2022986) q[2];
sx q[2];
rz(-2.9677128) q[2];
sx q[2];
rz(-0.536971) q[2];
rz(1.8592853) q[3];
sx q[3];
rz(-1.8351646) q[3];
sx q[3];
rz(-1.6622701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.031161664) q[0];
sx q[0];
rz(-2.3261676) q[0];
sx q[0];
rz(-1.331331) q[0];
rz(-1.4415461) q[1];
sx q[1];
rz(-1.6621637) q[1];
sx q[1];
rz(1.2225245) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0095897) q[0];
sx q[0];
rz(-1.1545657) q[0];
sx q[0];
rz(0.25968174) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6454846) q[2];
sx q[2];
rz(-2.1213581) q[2];
sx q[2];
rz(0.84361831) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4536087) q[1];
sx q[1];
rz(-0.88064945) q[1];
sx q[1];
rz(0.28120561) q[1];
x q[2];
rz(1.3243213) q[3];
sx q[3];
rz(-0.45650864) q[3];
sx q[3];
rz(-2.498561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.20087251) q[2];
sx q[2];
rz(-0.49786374) q[2];
sx q[2];
rz(2.8087924) q[2];
rz(-2.8424272) q[3];
sx q[3];
rz(-1.2406415) q[3];
sx q[3];
rz(-3.1201194) q[3];
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
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0126295) q[0];
sx q[0];
rz(-2.3733932) q[0];
sx q[0];
rz(2.8666038) q[0];
rz(-3.0601652) q[1];
sx q[1];
rz(-0.80137253) q[1];
sx q[1];
rz(1.3224695) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8997661) q[0];
sx q[0];
rz(-1.5112919) q[0];
sx q[0];
rz(1.8023947) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.36825387) q[2];
sx q[2];
rz(-2.2241631) q[2];
sx q[2];
rz(-2.7453932) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9857707) q[1];
sx q[1];
rz(-2.3233674) q[1];
sx q[1];
rz(-2.6236412) q[1];
rz(-pi) q[2];
rz(-0.72663088) q[3];
sx q[3];
rz(-1.0612349) q[3];
sx q[3];
rz(-1.0973615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.72121173) q[2];
sx q[2];
rz(-3.0215441) q[2];
sx q[2];
rz(-1.4018641) q[2];
rz(1.2841355) q[3];
sx q[3];
rz(-1.1893585) q[3];
sx q[3];
rz(3.106626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.5597124) q[0];
sx q[0];
rz(-1.548883) q[0];
sx q[0];
rz(-0.48932073) q[0];
rz(-3.1210461) q[1];
sx q[1];
rz(-1.2362044) q[1];
sx q[1];
rz(2.4403341) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4160181) q[0];
sx q[0];
rz(-1.8113935) q[0];
sx q[0];
rz(-2.1558236) q[0];
rz(-0.047766165) q[2];
sx q[2];
rz(-2.2738308) q[2];
sx q[2];
rz(-0.060279708) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8119231) q[1];
sx q[1];
rz(-0.81336248) q[1];
sx q[1];
rz(1.8977036) q[1];
x q[2];
rz(2.4550426) q[3];
sx q[3];
rz(-1.9420164) q[3];
sx q[3];
rz(-2.8438501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7265085) q[2];
sx q[2];
rz(-0.61843151) q[2];
sx q[2];
rz(0.48039594) q[2];
rz(1.4091617) q[3];
sx q[3];
rz(-1.3536072) q[3];
sx q[3];
rz(2.8044146) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82520634) q[0];
sx q[0];
rz(-1.9012863) q[0];
sx q[0];
rz(2.7301042) q[0];
rz(-1.8972137) q[1];
sx q[1];
rz(-2.0020516) q[1];
sx q[1];
rz(0.32434514) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.024687) q[0];
sx q[0];
rz(-1.9431264) q[0];
sx q[0];
rz(3.0850379) q[0];
rz(-pi) q[1];
rz(2.7997843) q[2];
sx q[2];
rz(-1.9035619) q[2];
sx q[2];
rz(1.8012799) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.577475) q[1];
sx q[1];
rz(-2.1337957) q[1];
sx q[1];
rz(3.0433118) q[1];
x q[2];
rz(-0.10735546) q[3];
sx q[3];
rz(-0.12190652) q[3];
sx q[3];
rz(-2.0481244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1856498) q[2];
sx q[2];
rz(-2.0273384) q[2];
sx q[2];
rz(-3.0246217) q[2];
rz(-2.1864435) q[3];
sx q[3];
rz(-2.9626466) q[3];
sx q[3];
rz(0.54168934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2627926) q[0];
sx q[0];
rz(-2.265082) q[0];
sx q[0];
rz(-2.2646917) q[0];
rz(-2.0211438) q[1];
sx q[1];
rz(-0.81605492) q[1];
sx q[1];
rz(1.1647404) q[1];
rz(-1.0638857) q[2];
sx q[2];
rz(-1.756466) q[2];
sx q[2];
rz(-1.9274439) q[2];
rz(-1.6616115) q[3];
sx q[3];
rz(-2.869893) q[3];
sx q[3];
rz(0.13520959) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
