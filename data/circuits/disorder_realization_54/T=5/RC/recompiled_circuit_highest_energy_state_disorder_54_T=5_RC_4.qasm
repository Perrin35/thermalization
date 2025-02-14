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
rz(-2.4009268) q[0];
sx q[0];
rz(-1.8494777) q[0];
sx q[0];
rz(-2.3722755) q[0];
rz(-1.5462592) q[1];
sx q[1];
rz(-0.21472628) q[1];
sx q[1];
rz(-2.5200342) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3174048) q[0];
sx q[0];
rz(-0.82302588) q[0];
sx q[0];
rz(-0.70931859) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5311422) q[2];
sx q[2];
rz(-0.41799212) q[2];
sx q[2];
rz(-2.5681873) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4159704) q[1];
sx q[1];
rz(-0.86639222) q[1];
sx q[1];
rz(-2.7654057) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0075843) q[3];
sx q[3];
rz(-1.4913627) q[3];
sx q[3];
rz(1.1955737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8605211) q[2];
sx q[2];
rz(-0.29025429) q[2];
sx q[2];
rz(0.93112913) q[2];
rz(0.18533254) q[3];
sx q[3];
rz(-2.1087746) q[3];
sx q[3];
rz(-1.7271656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(2.0463878) q[0];
sx q[0];
rz(-0.33400184) q[0];
sx q[0];
rz(0.97292501) q[0];
rz(0.75353777) q[1];
sx q[1];
rz(-0.85644186) q[1];
sx q[1];
rz(-0.10261745) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0282183) q[0];
sx q[0];
rz(-1.4732142) q[0];
sx q[0];
rz(-0.085208864) q[0];
rz(0.13394103) q[2];
sx q[2];
rz(-1.6689577) q[2];
sx q[2];
rz(2.5636473) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.18801461) q[1];
sx q[1];
rz(-0.17718592) q[1];
sx q[1];
rz(1.0832975) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2351949) q[3];
sx q[3];
rz(-1.4398267) q[3];
sx q[3];
rz(1.5168911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7332581) q[2];
sx q[2];
rz(-1.7034986) q[2];
sx q[2];
rz(-0.56780887) q[2];
rz(-2.7693977) q[3];
sx q[3];
rz(-1.8122383) q[3];
sx q[3];
rz(-1.6224434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55783015) q[0];
sx q[0];
rz(-0.74302858) q[0];
sx q[0];
rz(1.8999735) q[0];
rz(0.81635967) q[1];
sx q[1];
rz(-2.7528449) q[1];
sx q[1];
rz(-2.1519318) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1802487) q[0];
sx q[0];
rz(-1.3575866) q[0];
sx q[0];
rz(0.38589392) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0895592) q[2];
sx q[2];
rz(-2.9567869) q[2];
sx q[2];
rz(-1.081274) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6523) q[1];
sx q[1];
rz(-1.1649141) q[1];
sx q[1];
rz(-1.3764705) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9450467) q[3];
sx q[3];
rz(-1.47577) q[3];
sx q[3];
rz(-0.72685164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.376754) q[2];
sx q[2];
rz(-1.2412485) q[2];
sx q[2];
rz(-1.494701) q[2];
rz(-2.3946848) q[3];
sx q[3];
rz(-2.5266095) q[3];
sx q[3];
rz(-0.99228215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84585369) q[0];
sx q[0];
rz(-2.5339412) q[0];
sx q[0];
rz(2.1920152) q[0];
rz(1.1664561) q[1];
sx q[1];
rz(-1.8828853) q[1];
sx q[1];
rz(0.20225254) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1631781) q[0];
sx q[0];
rz(-1.6014153) q[0];
sx q[0];
rz(-3.091673) q[0];
rz(-0.20258383) q[2];
sx q[2];
rz(-2.4420945) q[2];
sx q[2];
rz(1.0765178) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7864322) q[1];
sx q[1];
rz(-1.0626241) q[1];
sx q[1];
rz(2.2340005) q[1];
x q[2];
rz(0.99503354) q[3];
sx q[3];
rz(-2.0821794) q[3];
sx q[3];
rz(-0.56727876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4919081) q[2];
sx q[2];
rz(-1.7077571) q[2];
sx q[2];
rz(1.3362308) q[2];
rz(0.46184552) q[3];
sx q[3];
rz(-1.0577842) q[3];
sx q[3];
rz(2.113078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.860054) q[0];
sx q[0];
rz(-0.10801948) q[0];
sx q[0];
rz(-2.7463013) q[0];
rz(-2.1832502) q[1];
sx q[1];
rz(-1.7610565) q[1];
sx q[1];
rz(-0.39452943) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0392691) q[0];
sx q[0];
rz(-1.8365055) q[0];
sx q[0];
rz(-0.018331176) q[0];
x q[1];
rz(-1.8498603) q[2];
sx q[2];
rz(-1.9968281) q[2];
sx q[2];
rz(-0.59471496) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3580618) q[1];
sx q[1];
rz(-2.3435842) q[1];
sx q[1];
rz(-1.111387) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7845527) q[3];
sx q[3];
rz(-0.96352623) q[3];
sx q[3];
rz(-0.34363817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.40256527) q[2];
sx q[2];
rz(-1.6141011) q[2];
sx q[2];
rz(-1.019545) q[2];
rz(1.839365) q[3];
sx q[3];
rz(-3.0890833) q[3];
sx q[3];
rz(0.70560613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.5604945) q[0];
sx q[0];
rz(-2.5769233) q[0];
sx q[0];
rz(1.0127006) q[0];
rz(-0.46214354) q[1];
sx q[1];
rz(-1.7299165) q[1];
sx q[1];
rz(3.0949458) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.93067) q[0];
sx q[0];
rz(-1.4477237) q[0];
sx q[0];
rz(0.010753429) q[0];
rz(-pi) q[1];
rz(1.2744802) q[2];
sx q[2];
rz(-1.8159397) q[2];
sx q[2];
rz(2.5398538) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.44429) q[1];
sx q[1];
rz(-1.238916) q[1];
sx q[1];
rz(2.1730971) q[1];
x q[2];
rz(3.1174692) q[3];
sx q[3];
rz(-1.1937448) q[3];
sx q[3];
rz(1.1386865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3095653) q[2];
sx q[2];
rz(-3.053061) q[2];
sx q[2];
rz(2.1868165) q[2];
rz(-0.66983062) q[3];
sx q[3];
rz(-1.9570743) q[3];
sx q[3];
rz(-0.36439103) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1466115) q[0];
sx q[0];
rz(-0.28047383) q[0];
sx q[0];
rz(-1.1501508) q[0];
rz(3.0448044) q[1];
sx q[1];
rz(-1.140118) q[1];
sx q[1];
rz(1.6614301) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1065001) q[0];
sx q[0];
rz(-2.273002) q[0];
sx q[0];
rz(2.0766792) q[0];
rz(-pi) q[1];
rz(-0.63235967) q[2];
sx q[2];
rz(-0.39146921) q[2];
sx q[2];
rz(-0.3492569) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0310049) q[1];
sx q[1];
rz(-1.2801941) q[1];
sx q[1];
rz(0.1433934) q[1];
rz(1.1365165) q[3];
sx q[3];
rz(-1.4308813) q[3];
sx q[3];
rz(-0.29014465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4074771) q[2];
sx q[2];
rz(-0.87782562) q[2];
sx q[2];
rz(-2.5131098) q[2];
rz(2.8001522) q[3];
sx q[3];
rz(-0.98613685) q[3];
sx q[3];
rz(0.83610523) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8516561) q[0];
sx q[0];
rz(-2.4735232) q[0];
sx q[0];
rz(-3.0776899) q[0];
rz(0.54197657) q[1];
sx q[1];
rz(-1.130645) q[1];
sx q[1];
rz(2.7625387) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.097103216) q[0];
sx q[0];
rz(-1.3987204) q[0];
sx q[0];
rz(1.8015566) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7070243) q[2];
sx q[2];
rz(-0.76644015) q[2];
sx q[2];
rz(-2.5640783) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2868834) q[1];
sx q[1];
rz(-1.5552551) q[1];
sx q[1];
rz(-2.5980224) q[1];
rz(-pi) q[2];
rz(2.1695215) q[3];
sx q[3];
rz(-2.4643371) q[3];
sx q[3];
rz(-2.2680091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1698251) q[2];
sx q[2];
rz(-1.1173893) q[2];
sx q[2];
rz(-2.8525412) q[2];
rz(1.0158094) q[3];
sx q[3];
rz(-2.2063875) q[3];
sx q[3];
rz(-2.9142006) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0357901) q[0];
sx q[0];
rz(-2.7547024) q[0];
sx q[0];
rz(-0.33045688) q[0];
rz(0.48078787) q[1];
sx q[1];
rz(-1.5590706) q[1];
sx q[1];
rz(-2.6343583) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15132771) q[0];
sx q[0];
rz(-2.2033407) q[0];
sx q[0];
rz(-0.91077478) q[0];
x q[1];
rz(2.4749996) q[2];
sx q[2];
rz(-0.5616411) q[2];
sx q[2];
rz(-3.0009342) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1363234) q[1];
sx q[1];
rz(-1.1526967) q[1];
sx q[1];
rz(-2.9015002) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5413569) q[3];
sx q[3];
rz(-1.5215989) q[3];
sx q[3];
rz(0.95358301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.5759739) q[2];
sx q[2];
rz(-2.2810292) q[2];
sx q[2];
rz(-2.7536075) q[2];
rz(3.0787789) q[3];
sx q[3];
rz(-2.4020436) q[3];
sx q[3];
rz(-1.1886103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0906618) q[0];
sx q[0];
rz(-0.020314038) q[0];
sx q[0];
rz(3.0860197) q[0];
rz(-2.9203501) q[1];
sx q[1];
rz(-2.1317) q[1];
sx q[1];
rz(-1.6411068) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7941754) q[0];
sx q[0];
rz(-2.5567434) q[0];
sx q[0];
rz(1.0561159) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5211721) q[2];
sx q[2];
rz(-1.7073586) q[2];
sx q[2];
rz(0.23601664) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6506429) q[1];
sx q[1];
rz(-1.888928) q[1];
sx q[1];
rz(1.9716119) q[1];
rz(-2.5440823) q[3];
sx q[3];
rz(-2.198285) q[3];
sx q[3];
rz(-1.6399872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1859493) q[2];
sx q[2];
rz(-0.187749) q[2];
sx q[2];
rz(-1.1090247) q[2];
rz(3.1254613) q[3];
sx q[3];
rz(-1.4343836) q[3];
sx q[3];
rz(-0.46512887) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7598509) q[0];
sx q[0];
rz(-1.009059) q[0];
sx q[0];
rz(-0.41643634) q[0];
rz(-2.3460559) q[1];
sx q[1];
rz(-1.3858613) q[1];
sx q[1];
rz(-0.081079986) q[1];
rz(1.9151081) q[2];
sx q[2];
rz(-1.7034265) q[2];
sx q[2];
rz(-1.8272057) q[2];
rz(-1.0352739) q[3];
sx q[3];
rz(-2.2761619) q[3];
sx q[3];
rz(-0.63009562) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
