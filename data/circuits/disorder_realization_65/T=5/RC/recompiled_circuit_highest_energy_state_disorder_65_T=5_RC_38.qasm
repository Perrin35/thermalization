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
rz(-1.3574358) q[0];
sx q[0];
rz(-2.5166002) q[0];
sx q[0];
rz(-1.5224737) q[0];
rz(1.1433262) q[1];
sx q[1];
rz(-2.5050617) q[1];
sx q[1];
rz(-1.358939) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9516884) q[0];
sx q[0];
rz(-2.1225327) q[0];
sx q[0];
rz(0.34073982) q[0];
rz(-pi) q[1];
rz(-2.236249) q[2];
sx q[2];
rz(-0.45768379) q[2];
sx q[2];
rz(-0.25027088) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2083019) q[1];
sx q[1];
rz(-1.33889) q[1];
sx q[1];
rz(1.8342706) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.10349689) q[3];
sx q[3];
rz(-2.0465133) q[3];
sx q[3];
rz(2.5535943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6978567) q[2];
sx q[2];
rz(-0.82252684) q[2];
sx q[2];
rz(-2.6653384) q[2];
rz(3.0016628) q[3];
sx q[3];
rz(-1.5014476) q[3];
sx q[3];
rz(-0.073277624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3530465) q[0];
sx q[0];
rz(-1.5209501) q[0];
sx q[0];
rz(0.36343685) q[0];
rz(-0.67047554) q[1];
sx q[1];
rz(-1.8164219) q[1];
sx q[1];
rz(-1.0796116) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99936679) q[0];
sx q[0];
rz(-1.8086664) q[0];
sx q[0];
rz(-1.0269985) q[0];
x q[1];
rz(3.011376) q[2];
sx q[2];
rz(-1.0417582) q[2];
sx q[2];
rz(-2.0646273) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9448589) q[1];
sx q[1];
rz(-0.24538876) q[1];
sx q[1];
rz(-1.3887482) q[1];
x q[2];
rz(0.3877181) q[3];
sx q[3];
rz(-2.5250375) q[3];
sx q[3];
rz(-1.8247557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4662027) q[2];
sx q[2];
rz(-2.4066996) q[2];
sx q[2];
rz(-2.7030763) q[2];
rz(2.130326) q[3];
sx q[3];
rz(-0.68532419) q[3];
sx q[3];
rz(1.6605759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48280516) q[0];
sx q[0];
rz(-0.52849448) q[0];
sx q[0];
rz(-2.3129789) q[0];
rz(-1.8939182) q[1];
sx q[1];
rz(-1.9745461) q[1];
sx q[1];
rz(2.8819328) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8369401) q[0];
sx q[0];
rz(-0.90422179) q[0];
sx q[0];
rz(2.7949692) q[0];
rz(1.2010048) q[2];
sx q[2];
rz(-1.2134873) q[2];
sx q[2];
rz(1.0543038) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6344172) q[1];
sx q[1];
rz(-0.89752642) q[1];
sx q[1];
rz(2.5310764) q[1];
x q[2];
rz(1.9454221) q[3];
sx q[3];
rz(-2.4536479) q[3];
sx q[3];
rz(1.3796395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9344249) q[2];
sx q[2];
rz(-1.5639037) q[2];
sx q[2];
rz(-0.820532) q[2];
rz(-0.82646838) q[3];
sx q[3];
rz(-0.99244899) q[3];
sx q[3];
rz(-1.3124527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25927037) q[0];
sx q[0];
rz(-0.81040183) q[0];
sx q[0];
rz(-2.209254) q[0];
rz(-1.3177634) q[1];
sx q[1];
rz(-1.1612786) q[1];
sx q[1];
rz(0.12277776) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49157351) q[0];
sx q[0];
rz(-0.54560018) q[0];
sx q[0];
rz(-0.16262098) q[0];
rz(-pi) q[1];
rz(-1.2447692) q[2];
sx q[2];
rz(-0.78067987) q[2];
sx q[2];
rz(2.8512466) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9726561) q[1];
sx q[1];
rz(-2.3984809) q[1];
sx q[1];
rz(-1.572322) q[1];
rz(-pi) q[2];
rz(0.56468954) q[3];
sx q[3];
rz(-1.3211234) q[3];
sx q[3];
rz(-0.76631977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.54470283) q[2];
sx q[2];
rz(-2.2201316) q[2];
sx q[2];
rz(0.60453647) q[2];
rz(2.4321411) q[3];
sx q[3];
rz(-0.76990288) q[3];
sx q[3];
rz(2.3179222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17968793) q[0];
sx q[0];
rz(-3.0178495) q[0];
sx q[0];
rz(1.7748348) q[0];
rz(-2.1429515) q[1];
sx q[1];
rz(-0.81592453) q[1];
sx q[1];
rz(-2.6008115) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3635938) q[0];
sx q[0];
rz(-0.25083298) q[0];
sx q[0];
rz(0.59441113) q[0];
rz(2.3983058) q[2];
sx q[2];
rz(-2.1319445) q[2];
sx q[2];
rz(2.4744792) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.65366983) q[1];
sx q[1];
rz(-1.2752295) q[1];
sx q[1];
rz(1.3699156) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.56638797) q[3];
sx q[3];
rz(-1.4175804) q[3];
sx q[3];
rz(-0.18463102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.46405408) q[2];
sx q[2];
rz(-0.81391922) q[2];
sx q[2];
rz(0.82826725) q[2];
rz(-1.6895435) q[3];
sx q[3];
rz(-1.4938846) q[3];
sx q[3];
rz(-0.90788466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0238817) q[0];
sx q[0];
rz(-1.1957059) q[0];
sx q[0];
rz(1.9697795) q[0];
rz(-0.68012971) q[1];
sx q[1];
rz(-1.2260194) q[1];
sx q[1];
rz(0.5562869) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6345917) q[0];
sx q[0];
rz(-1.9456353) q[0];
sx q[0];
rz(-0.27347538) q[0];
x q[1];
rz(-1.4242709) q[2];
sx q[2];
rz(-1.5735448) q[2];
sx q[2];
rz(2.3607852) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.57842931) q[1];
sx q[1];
rz(-0.99908864) q[1];
sx q[1];
rz(1.6469547) q[1];
rz(-pi) q[2];
rz(-0.61304355) q[3];
sx q[3];
rz(-2.4073232) q[3];
sx q[3];
rz(-3.120852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3351626) q[2];
sx q[2];
rz(-1.282629) q[2];
sx q[2];
rz(0.20509091) q[2];
rz(2.3388376) q[3];
sx q[3];
rz(-2.0358678) q[3];
sx q[3];
rz(-2.5870489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71451181) q[0];
sx q[0];
rz(-1.0733805) q[0];
sx q[0];
rz(0.22739534) q[0];
rz(-2.3511476) q[1];
sx q[1];
rz(-0.77722725) q[1];
sx q[1];
rz(-0.8367742) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28115434) q[0];
sx q[0];
rz(-1.3520762) q[0];
sx q[0];
rz(0.52978306) q[0];
rz(1.927194) q[2];
sx q[2];
rz(-2.1292392) q[2];
sx q[2];
rz(-1.3791549) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0322276) q[1];
sx q[1];
rz(-1.9739474) q[1];
sx q[1];
rz(1.228986) q[1];
rz(-pi) q[2];
rz(2.3804139) q[3];
sx q[3];
rz(-1.2952779) q[3];
sx q[3];
rz(-0.31106424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5922015) q[2];
sx q[2];
rz(-1.6624007) q[2];
sx q[2];
rz(-0.57360348) q[2];
rz(-2.8163689) q[3];
sx q[3];
rz(-2.2461788) q[3];
sx q[3];
rz(0.4044683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0322872) q[0];
sx q[0];
rz(-0.17806299) q[0];
sx q[0];
rz(-2.4626515) q[0];
rz(3.0067054) q[1];
sx q[1];
rz(-1.0604246) q[1];
sx q[1];
rz(-1.7178242) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2436062) q[0];
sx q[0];
rz(-0.73315128) q[0];
sx q[0];
rz(-0.85444684) q[0];
rz(0.94332327) q[2];
sx q[2];
rz(-2.3296142) q[2];
sx q[2];
rz(-1.8705778) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0404741) q[1];
sx q[1];
rz(-2.3777739) q[1];
sx q[1];
rz(-1.5274672) q[1];
x q[2];
rz(-2.0349488) q[3];
sx q[3];
rz(-1.8309793) q[3];
sx q[3];
rz(-1.2342208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5941102) q[2];
sx q[2];
rz(-0.91071931) q[2];
sx q[2];
rz(0.2317079) q[2];
rz(-2.1314651) q[3];
sx q[3];
rz(-1.2332656) q[3];
sx q[3];
rz(-2.2505984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6758839) q[0];
sx q[0];
rz(-0.75513419) q[0];
sx q[0];
rz(-0.51374197) q[0];
rz(1.7087917) q[1];
sx q[1];
rz(-1.2756196) q[1];
sx q[1];
rz(-1.6339711) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20376539) q[0];
sx q[0];
rz(-1.7047593) q[0];
sx q[0];
rz(1.521827) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8880226) q[2];
sx q[2];
rz(-1.7268983) q[2];
sx q[2];
rz(1.9954322) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6027939) q[1];
sx q[1];
rz(-1.4414296) q[1];
sx q[1];
rz(1.4461229) q[1];
x q[2];
rz(-1.6674897) q[3];
sx q[3];
rz(-2.8887199) q[3];
sx q[3];
rz(-2.0677572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.12019176) q[2];
sx q[2];
rz(-0.84453619) q[2];
sx q[2];
rz(2.8821442) q[2];
rz(1.9868959) q[3];
sx q[3];
rz(-2.3754933) q[3];
sx q[3];
rz(1.8416789) q[3];
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
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65086377) q[0];
sx q[0];
rz(-1.6614953) q[0];
sx q[0];
rz(-1.65253) q[0];
rz(0.64000714) q[1];
sx q[1];
rz(-1.0970486) q[1];
sx q[1];
rz(2.1515813) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31371597) q[0];
sx q[0];
rz(-1.4826688) q[0];
sx q[0];
rz(0.14493305) q[0];
x q[1];
rz(2.0406538) q[2];
sx q[2];
rz(-2.2820916) q[2];
sx q[2];
rz(-0.92352042) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1460625) q[1];
sx q[1];
rz(-1.7319253) q[1];
sx q[1];
rz(-1.1481768) q[1];
rz(-pi) q[2];
rz(-1.3481989) q[3];
sx q[3];
rz(-1.5951459) q[3];
sx q[3];
rz(-2.3241732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9330357) q[2];
sx q[2];
rz(-1.2988337) q[2];
sx q[2];
rz(2.9162858) q[2];
rz(-0.92132583) q[3];
sx q[3];
rz(-2.7511629) q[3];
sx q[3];
rz(0.050617378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0113572) q[0];
sx q[0];
rz(-0.36778944) q[0];
sx q[0];
rz(0.854048) q[0];
rz(-1.8232952) q[1];
sx q[1];
rz(-1.5250991) q[1];
sx q[1];
rz(-0.93223882) q[1];
rz(-1.3110326) q[2];
sx q[2];
rz(-1.4767892) q[2];
sx q[2];
rz(-1.3108419) q[2];
rz(2.1193567) q[3];
sx q[3];
rz(-1.4132186) q[3];
sx q[3];
rz(1.3506387) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
