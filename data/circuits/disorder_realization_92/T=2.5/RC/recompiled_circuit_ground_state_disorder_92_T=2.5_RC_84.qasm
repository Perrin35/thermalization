OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.2607245) q[0];
sx q[0];
rz(-0.27624929) q[0];
sx q[0];
rz(-2.2978388) q[0];
rz(-2.1387956) q[1];
sx q[1];
rz(-0.83477867) q[1];
sx q[1];
rz(-0.53258449) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7363889) q[0];
sx q[0];
rz(-1.4215934) q[0];
sx q[0];
rz(0.91513855) q[0];
rz(2.631408) q[2];
sx q[2];
rz(-1.6932339) q[2];
sx q[2];
rz(-2.7363026) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3644173) q[1];
sx q[1];
rz(-0.88349408) q[1];
sx q[1];
rz(2.9455094) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3477676) q[3];
sx q[3];
rz(-0.50758368) q[3];
sx q[3];
rz(0.48963293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8325995) q[2];
sx q[2];
rz(-1.7367881) q[2];
sx q[2];
rz(3.0187606) q[2];
rz(-3.099856) q[3];
sx q[3];
rz(-1.583497) q[3];
sx q[3];
rz(0.48833716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9915344) q[0];
sx q[0];
rz(-0.3287065) q[0];
sx q[0];
rz(-1.1154037) q[0];
rz(-0.54488048) q[1];
sx q[1];
rz(-1.7439525) q[1];
sx q[1];
rz(2.5741408) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79115907) q[0];
sx q[0];
rz(-1.6164403) q[0];
sx q[0];
rz(-1.537039) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5390776) q[2];
sx q[2];
rz(-2.3467772) q[2];
sx q[2];
rz(-1.0358126) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9055674) q[1];
sx q[1];
rz(-1.4705338) q[1];
sx q[1];
rz(-0.57159337) q[1];
rz(-pi) q[2];
rz(1.2699158) q[3];
sx q[3];
rz(-1.1071651) q[3];
sx q[3];
rz(-0.6248354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6313717) q[2];
sx q[2];
rz(-3.0192182) q[2];
sx q[2];
rz(-0.1046293) q[2];
rz(-2.6369324) q[3];
sx q[3];
rz(-1.3480836) q[3];
sx q[3];
rz(-0.73205718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1402682) q[0];
sx q[0];
rz(-2.9502385) q[0];
sx q[0];
rz(1.485317) q[0];
rz(0.56400076) q[1];
sx q[1];
rz(-1.0537078) q[1];
sx q[1];
rz(-0.82411134) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38496537) q[0];
sx q[0];
rz(-2.5838296) q[0];
sx q[0];
rz(1.8855479) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1698157) q[2];
sx q[2];
rz(-2.2702771) q[2];
sx q[2];
rz(-2.5339507) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.66131683) q[1];
sx q[1];
rz(-1.7719898) q[1];
sx q[1];
rz(0.43432216) q[1];
rz(-3.0715354) q[3];
sx q[3];
rz(-1.3590711) q[3];
sx q[3];
rz(2.5222561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.092502681) q[2];
sx q[2];
rz(-2.8503032) q[2];
sx q[2];
rz(-1.5352486) q[2];
rz(-2.2384079) q[3];
sx q[3];
rz(-1.817037) q[3];
sx q[3];
rz(-2.7170392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9400738) q[0];
sx q[0];
rz(-0.98648447) q[0];
sx q[0];
rz(0.41341138) q[0];
rz(-0.19551936) q[1];
sx q[1];
rz(-1.0370516) q[1];
sx q[1];
rz(-1.4541385) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66178807) q[0];
sx q[0];
rz(-2.256752) q[0];
sx q[0];
rz(0.16459008) q[0];
rz(0.35572534) q[2];
sx q[2];
rz(-1.3918096) q[2];
sx q[2];
rz(0.21790522) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5315872) q[1];
sx q[1];
rz(-1.0949228) q[1];
sx q[1];
rz(0.31756084) q[1];
rz(-pi) q[2];
rz(-2.4471483) q[3];
sx q[3];
rz(-1.4141091) q[3];
sx q[3];
rz(0.82686916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2426408) q[2];
sx q[2];
rz(-0.61598888) q[2];
sx q[2];
rz(1.8757437) q[2];
rz(-0.85593456) q[3];
sx q[3];
rz(-1.7549691) q[3];
sx q[3];
rz(-1.2148733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4796546) q[0];
sx q[0];
rz(-2.9582773) q[0];
sx q[0];
rz(1.4187752) q[0];
rz(1.4310369) q[1];
sx q[1];
rz(-1.5242256) q[1];
sx q[1];
rz(0.72351825) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.60155) q[0];
sx q[0];
rz(-1.7474084) q[0];
sx q[0];
rz(2.5536182) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3143598) q[2];
sx q[2];
rz(-2.1783354) q[2];
sx q[2];
rz(0.22996685) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4505087) q[1];
sx q[1];
rz(-1.5732068) q[1];
sx q[1];
rz(0.86402135) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5635023) q[3];
sx q[3];
rz(-1.9230584) q[3];
sx q[3];
rz(-1.9520813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7001223) q[2];
sx q[2];
rz(-2.1258326) q[2];
sx q[2];
rz(2.7777242) q[2];
rz(-1.8338592) q[3];
sx q[3];
rz(-1.4738844) q[3];
sx q[3];
rz(1.7893121) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7135007) q[0];
sx q[0];
rz(-2.2245753) q[0];
sx q[0];
rz(-2.2416903) q[0];
rz(1.046754) q[1];
sx q[1];
rz(-2.7622107) q[1];
sx q[1];
rz(2.5894763) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22495843) q[0];
sx q[0];
rz(-2.8021268) q[0];
sx q[0];
rz(0.084966226) q[0];
x q[1];
rz(-2.1119166) q[2];
sx q[2];
rz(-0.94391247) q[2];
sx q[2];
rz(-2.1205001) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0830517) q[1];
sx q[1];
rz(-1.6612435) q[1];
sx q[1];
rz(1.3988711) q[1];
rz(-pi) q[2];
rz(1.6523727) q[3];
sx q[3];
rz(-2.1595567) q[3];
sx q[3];
rz(2.0643864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.8206574) q[2];
sx q[2];
rz(-2.4652017) q[2];
sx q[2];
rz(0.19860849) q[2];
rz(-1.1292388) q[3];
sx q[3];
rz(-1.4720474) q[3];
sx q[3];
rz(-2.3614597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8010913) q[0];
sx q[0];
rz(-1.8393562) q[0];
sx q[0];
rz(-2.4275725) q[0];
rz(1.9188312) q[1];
sx q[1];
rz(-2.265265) q[1];
sx q[1];
rz(-0.27539918) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.794316) q[0];
sx q[0];
rz(-1.9660608) q[0];
sx q[0];
rz(-2.9121132) q[0];
rz(2.8267639) q[2];
sx q[2];
rz(-2.1595975) q[2];
sx q[2];
rz(-2.3662629) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.022897331) q[1];
sx q[1];
rz(-1.1015019) q[1];
sx q[1];
rz(2.1170627) q[1];
x q[2];
rz(0.72526284) q[3];
sx q[3];
rz(-1.5607087) q[3];
sx q[3];
rz(0.46848224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8347281) q[2];
sx q[2];
rz(-1.6708259) q[2];
sx q[2];
rz(-1.4677706) q[2];
rz(1.5214527) q[3];
sx q[3];
rz(-0.68632564) q[3];
sx q[3];
rz(-0.93792382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15243212) q[0];
sx q[0];
rz(-2.2280362) q[0];
sx q[0];
rz(2.1754225) q[0];
rz(0.78978157) q[1];
sx q[1];
rz(-0.25895324) q[1];
sx q[1];
rz(0.97377473) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0278661) q[0];
sx q[0];
rz(-1.4513512) q[0];
sx q[0];
rz(-1.3238504) q[0];
rz(-pi) q[1];
rz(-0.42736407) q[2];
sx q[2];
rz(-1.9081429) q[2];
sx q[2];
rz(0.34572476) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3605347) q[1];
sx q[1];
rz(-1.2670749) q[1];
sx q[1];
rz(0.93704929) q[1];
rz(-pi) q[2];
rz(-1.6867181) q[3];
sx q[3];
rz(-2.6629857) q[3];
sx q[3];
rz(-0.49540621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8345653) q[2];
sx q[2];
rz(-1.4296738) q[2];
sx q[2];
rz(-2.5359421) q[2];
rz(1.0904795) q[3];
sx q[3];
rz(-2.8993789) q[3];
sx q[3];
rz(0.098701326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7379446) q[0];
sx q[0];
rz(-0.048796766) q[0];
sx q[0];
rz(-0.57156372) q[0];
rz(-0.38772186) q[1];
sx q[1];
rz(-2.3523836) q[1];
sx q[1];
rz(-2.2472084) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8799202) q[0];
sx q[0];
rz(-0.74194569) q[0];
sx q[0];
rz(2.3319753) q[0];
rz(-3.0249075) q[2];
sx q[2];
rz(-1.0067847) q[2];
sx q[2];
rz(-2.2224768) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9891178) q[1];
sx q[1];
rz(-0.47429171) q[1];
sx q[1];
rz(-2.1211339) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4924516) q[3];
sx q[3];
rz(-2.0781029) q[3];
sx q[3];
rz(-2.0149734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9992708) q[2];
sx q[2];
rz(-2.2597376) q[2];
sx q[2];
rz(-0.87232653) q[2];
rz(0.074660389) q[3];
sx q[3];
rz(-1.229769) q[3];
sx q[3];
rz(-2.5653896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1416624) q[0];
sx q[0];
rz(-2.7003728) q[0];
sx q[0];
rz(0.73750752) q[0];
rz(-2.4420338) q[1];
sx q[1];
rz(-2.12205) q[1];
sx q[1];
rz(0.84651822) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.072089) q[0];
sx q[0];
rz(-2.3583625) q[0];
sx q[0];
rz(-2.5602719) q[0];
rz(-1.4581095) q[2];
sx q[2];
rz(-0.42851617) q[2];
sx q[2];
rz(-0.90135114) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6168931) q[1];
sx q[1];
rz(-1.9977101) q[1];
sx q[1];
rz(2.9906143) q[1];
x q[2];
rz(-0.41937866) q[3];
sx q[3];
rz(-0.43361317) q[3];
sx q[3];
rz(-0.24672844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.49125853) q[2];
sx q[2];
rz(-1.9776055) q[2];
sx q[2];
rz(2.7872046) q[2];
rz(-2.3645511) q[3];
sx q[3];
rz(-0.6207501) q[3];
sx q[3];
rz(1.0907772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.95018321) q[0];
sx q[0];
rz(-1.2744899) q[0];
sx q[0];
rz(0.50344678) q[0];
rz(1.7783816) q[1];
sx q[1];
rz(-0.53047219) q[1];
sx q[1];
rz(-0.06906876) q[1];
rz(0.10748482) q[2];
sx q[2];
rz(-0.86089118) q[2];
sx q[2];
rz(-0.64313342) q[2];
rz(-1.1152399) q[3];
sx q[3];
rz(-1.2610049) q[3];
sx q[3];
rz(-0.42381248) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
