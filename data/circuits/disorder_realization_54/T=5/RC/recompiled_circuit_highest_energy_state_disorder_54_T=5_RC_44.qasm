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
rz(0.74066585) q[0];
sx q[0];
rz(1.8494777) q[0];
sx q[0];
rz(10.194095) q[0];
rz(-1.5462592) q[1];
sx q[1];
rz(-0.21472628) q[1];
sx q[1];
rz(-2.5200342) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7249994) q[0];
sx q[0];
rz(-2.160797) q[0];
sx q[0];
rz(-2.1830465) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1530959) q[2];
sx q[2];
rz(-1.5547032) q[2];
sx q[2];
rz(0.96114767) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.5947508) q[1];
sx q[1];
rz(-1.2870645) q[1];
sx q[1];
rz(2.3111514) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0075843) q[3];
sx q[3];
rz(-1.6502299) q[3];
sx q[3];
rz(1.1955737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2810716) q[2];
sx q[2];
rz(-0.29025429) q[2];
sx q[2];
rz(2.2104635) q[2];
rz(-2.9562601) q[3];
sx q[3];
rz(-2.1087746) q[3];
sx q[3];
rz(1.4144271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.0463878) q[0];
sx q[0];
rz(-0.33400184) q[0];
sx q[0];
rz(2.1686676) q[0];
rz(2.3880549) q[1];
sx q[1];
rz(-0.85644186) q[1];
sx q[1];
rz(0.10261745) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1133744) q[0];
sx q[0];
rz(-1.6683785) q[0];
sx q[0];
rz(-3.0563838) q[0];
rz(-pi) q[1];
rz(-3.0076516) q[2];
sx q[2];
rz(-1.472635) q[2];
sx q[2];
rz(0.5779454) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.18801461) q[1];
sx q[1];
rz(-0.17718592) q[1];
sx q[1];
rz(-2.0582951) q[1];
rz(-0.13861175) q[3];
sx q[3];
rz(-1.2381805) q[3];
sx q[3];
rz(-3.1332071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.40833452) q[2];
sx q[2];
rz(-1.7034986) q[2];
sx q[2];
rz(0.56780887) q[2];
rz(-0.37219498) q[3];
sx q[3];
rz(-1.3293543) q[3];
sx q[3];
rz(1.5191493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55783015) q[0];
sx q[0];
rz(-0.74302858) q[0];
sx q[0];
rz(1.8999735) q[0];
rz(2.325233) q[1];
sx q[1];
rz(-2.7528449) q[1];
sx q[1];
rz(2.1519318) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12927221) q[0];
sx q[0];
rz(-2.7033157) q[0];
sx q[0];
rz(-2.6196036) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0520334) q[2];
sx q[2];
rz(-2.9567869) q[2];
sx q[2];
rz(2.0603186) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.1376393) q[1];
sx q[1];
rz(-1.7491566) q[1];
sx q[1];
rz(-0.41282005) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9450467) q[3];
sx q[3];
rz(-1.47577) q[3];
sx q[3];
rz(-2.414741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7648387) q[2];
sx q[2];
rz(-1.2412485) q[2];
sx q[2];
rz(1.494701) q[2];
rz(-0.74690789) q[3];
sx q[3];
rz(-2.5266095) q[3];
sx q[3];
rz(-2.1493105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.295739) q[0];
sx q[0];
rz(-0.60765147) q[0];
sx q[0];
rz(-0.94957748) q[0];
rz(-1.1664561) q[1];
sx q[1];
rz(-1.8828853) q[1];
sx q[1];
rz(2.9393401) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14213053) q[0];
sx q[0];
rz(-3.0830374) q[0];
sx q[0];
rz(2.5910795) q[0];
x q[1];
rz(-1.4030898) q[2];
sx q[2];
rz(-0.88838345) q[2];
sx q[2];
rz(-0.81425999) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7864322) q[1];
sx q[1];
rz(-2.0789685) q[1];
sx q[1];
rz(2.2340005) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.77029404) q[3];
sx q[3];
rz(-2.3912734) q[3];
sx q[3];
rz(1.4920774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.64968455) q[2];
sx q[2];
rz(-1.4338355) q[2];
sx q[2];
rz(1.8053619) q[2];
rz(-0.46184552) q[3];
sx q[3];
rz(-2.0838085) q[3];
sx q[3];
rz(2.113078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.860054) q[0];
sx q[0];
rz(-0.10801948) q[0];
sx q[0];
rz(0.39529133) q[0];
rz(2.1832502) q[1];
sx q[1];
rz(-1.3805362) q[1];
sx q[1];
rz(2.7470632) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4636587) q[0];
sx q[0];
rz(-1.5531085) q[0];
sx q[0];
rz(1.3050446) q[0];
rz(2.5960693) q[2];
sx q[2];
rz(-2.6370271) q[2];
sx q[2];
rz(-1.9406174) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9745512) q[1];
sx q[1];
rz(-2.2675505) q[1];
sx q[1];
rz(-0.42679326) q[1];
rz(2.5234918) q[3];
sx q[3];
rz(-1.7458946) q[3];
sx q[3];
rz(-1.7911946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.40256527) q[2];
sx q[2];
rz(-1.5274916) q[2];
sx q[2];
rz(-2.1220477) q[2];
rz(-1.3022276) q[3];
sx q[3];
rz(-3.0890833) q[3];
sx q[3];
rz(-2.4359865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5810982) q[0];
sx q[0];
rz(-0.56466931) q[0];
sx q[0];
rz(-2.1288921) q[0];
rz(0.46214354) q[1];
sx q[1];
rz(-1.7299165) q[1];
sx q[1];
rz(0.046646811) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.93067) q[0];
sx q[0];
rz(-1.6938689) q[0];
sx q[0];
rz(-3.1308392) q[0];
rz(-pi) q[1];
rz(2.279206) q[2];
sx q[2];
rz(-2.7593331) q[2];
sx q[2];
rz(1.5008937) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0477284) q[1];
sx q[1];
rz(-1.0055491) q[1];
sx q[1];
rz(-0.39611343) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1936452) q[3];
sx q[3];
rz(-1.5483678) q[3];
sx q[3];
rz(-2.7183661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8320273) q[2];
sx q[2];
rz(-0.088531606) q[2];
sx q[2];
rz(-2.1868165) q[2];
rz(0.66983062) q[3];
sx q[3];
rz(-1.9570743) q[3];
sx q[3];
rz(-2.7772016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-1.1466115) q[0];
sx q[0];
rz(-2.8611188) q[0];
sx q[0];
rz(1.9914419) q[0];
rz(-0.096788302) q[1];
sx q[1];
rz(-1.140118) q[1];
sx q[1];
rz(-1.4801625) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8154248) q[0];
sx q[0];
rz(-0.83954869) q[0];
sx q[0];
rz(0.52013784) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3315013) q[2];
sx q[2];
rz(-1.2579489) q[2];
sx q[2];
rz(-2.8205736) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7227308) q[1];
sx q[1];
rz(-1.7081385) q[1];
sx q[1];
rz(1.2773499) q[1];
rz(0.15401675) q[3];
sx q[3];
rz(-1.1410442) q[3];
sx q[3];
rz(1.3452443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.73411554) q[2];
sx q[2];
rz(-0.87782562) q[2];
sx q[2];
rz(-2.5131098) q[2];
rz(2.8001522) q[3];
sx q[3];
rz(-2.1554558) q[3];
sx q[3];
rz(-0.83610523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8516561) q[0];
sx q[0];
rz(-2.4735232) q[0];
sx q[0];
rz(-0.063902721) q[0];
rz(-0.54197657) q[1];
sx q[1];
rz(-2.0109476) q[1];
sx q[1];
rz(-0.37905395) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0444894) q[0];
sx q[0];
rz(-1.3987204) q[0];
sx q[0];
rz(-1.8015566) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3325868) q[2];
sx q[2];
rz(-1.6651285) q[2];
sx q[2];
rz(1.091711) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8547092) q[1];
sx q[1];
rz(-1.5552551) q[1];
sx q[1];
rz(2.5980224) q[1];
rz(-pi) q[2];
rz(2.1695215) q[3];
sx q[3];
rz(-2.4643371) q[3];
sx q[3];
rz(-2.2680091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9717676) q[2];
sx q[2];
rz(-2.0242033) q[2];
sx q[2];
rz(-2.8525412) q[2];
rz(1.0158094) q[3];
sx q[3];
rz(-0.93520516) q[3];
sx q[3];
rz(2.9142006) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10580258) q[0];
sx q[0];
rz(-0.38689026) q[0];
sx q[0];
rz(0.33045688) q[0];
rz(2.6608048) q[1];
sx q[1];
rz(-1.5590706) q[1];
sx q[1];
rz(2.6343583) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15132771) q[0];
sx q[0];
rz(-2.2033407) q[0];
sx q[0];
rz(2.2308179) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1997517) q[2];
sx q[2];
rz(-1.1389274) q[2];
sx q[2];
rz(0.88954207) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5937336) q[1];
sx q[1];
rz(-0.47858176) q[1];
sx q[1];
rz(2.0622158) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5112002) q[3];
sx q[3];
rz(-0.97138849) q[3];
sx q[3];
rz(-0.58356482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.5759739) q[2];
sx q[2];
rz(-2.2810292) q[2];
sx q[2];
rz(0.38798517) q[2];
rz(-0.062813736) q[3];
sx q[3];
rz(-2.4020436) q[3];
sx q[3];
rz(-1.1886103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0906618) q[0];
sx q[0];
rz(-3.1212786) q[0];
sx q[0];
rz(0.055572979) q[0];
rz(-2.9203501) q[1];
sx q[1];
rz(-1.0098927) q[1];
sx q[1];
rz(1.6411068) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.198198) q[0];
sx q[0];
rz(-2.072081) q[0];
sx q[0];
rz(0.31507612) q[0];
x q[1];
rz(-1.6204206) q[2];
sx q[2];
rz(-1.7073586) q[2];
sx q[2];
rz(-0.23601664) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.55602395) q[1];
sx q[1];
rz(-2.635286) q[1];
sx q[1];
rz(2.271818) q[1];
rz(-pi) q[2];
rz(-0.9110578) q[3];
sx q[3];
rz(-0.83759901) q[3];
sx q[3];
rz(-0.78105951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9556433) q[2];
sx q[2];
rz(-0.187749) q[2];
sx q[2];
rz(-2.032568) q[2];
rz(-0.016131314) q[3];
sx q[3];
rz(-1.707209) q[3];
sx q[3];
rz(0.46512887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3817417) q[0];
sx q[0];
rz(-1.009059) q[0];
sx q[0];
rz(-0.41643634) q[0];
rz(2.3460559) q[1];
sx q[1];
rz(-1.7557314) q[1];
sx q[1];
rz(3.0605127) q[1];
rz(0.14079413) q[2];
sx q[2];
rz(-1.2296321) q[2];
sx q[2];
rz(2.932569) q[2];
rz(1.0352739) q[3];
sx q[3];
rz(-0.86543075) q[3];
sx q[3];
rz(2.511497) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
