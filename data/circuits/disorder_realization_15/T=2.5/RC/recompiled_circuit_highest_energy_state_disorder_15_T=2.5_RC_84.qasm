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
rz(-2.5315142) q[0];
sx q[0];
rz(-0.22992034) q[0];
sx q[0];
rz(-2.4218986) q[0];
rz(2.5201058) q[1];
sx q[1];
rz(-1.7401594) q[1];
sx q[1];
rz(-3.016234) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2456477) q[0];
sx q[0];
rz(-0.60310793) q[0];
sx q[0];
rz(-2.3853517) q[0];
rz(-2.8530082) q[2];
sx q[2];
rz(-1.6615531) q[2];
sx q[2];
rz(1.8708269) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4957089) q[1];
sx q[1];
rz(-1.570374) q[1];
sx q[1];
rz(1.5721815) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.8475581) q[3];
sx q[3];
rz(-1.2985909) q[3];
sx q[3];
rz(2.6295626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7980935) q[2];
sx q[2];
rz(-2.7332879) q[2];
sx q[2];
rz(2.2957392) q[2];
rz(-0.79126233) q[3];
sx q[3];
rz(-3.1281804) q[3];
sx q[3];
rz(-3.0748034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7283519) q[0];
sx q[0];
rz(-2.6744196) q[0];
sx q[0];
rz(3.0979284) q[0];
rz(-1.5664258) q[1];
sx q[1];
rz(-1.3730201) q[1];
sx q[1];
rz(1.498819) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7455755) q[0];
sx q[0];
rz(-1.5748576) q[0];
sx q[0];
rz(-2.1745357) q[0];
rz(-pi) q[1];
rz(1.5506707) q[2];
sx q[2];
rz(-2.5701036) q[2];
sx q[2];
rz(-3.1352941) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.0017569204) q[1];
sx q[1];
rz(-3.0593605) q[1];
sx q[1];
rz(2.7539192) q[1];
x q[2];
rz(-1.2263857) q[3];
sx q[3];
rz(-1.5395172) q[3];
sx q[3];
rz(2.7117604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0955536) q[2];
sx q[2];
rz(-0.150103) q[2];
sx q[2];
rz(2.6138439) q[2];
rz(2.2857417) q[3];
sx q[3];
rz(-3.1400883) q[3];
sx q[3];
rz(1.9201479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7432778) q[0];
sx q[0];
rz(-0.97284955) q[0];
sx q[0];
rz(1.995218) q[0];
rz(-1.758681) q[1];
sx q[1];
rz(-2.8490366) q[1];
sx q[1];
rz(-3.0376099) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87097528) q[0];
sx q[0];
rz(-1.7954602) q[0];
sx q[0];
rz(-1.490834) q[0];
x q[1];
rz(-1.8431333) q[2];
sx q[2];
rz(-1.587279) q[2];
sx q[2];
rz(-2.4695244) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5755363) q[1];
sx q[1];
rz(-2.3588099) q[1];
sx q[1];
rz(-2.8498883) q[1];
rz(0.35211925) q[3];
sx q[3];
rz(-2.8330292) q[3];
sx q[3];
rz(1.4328014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.18342239) q[2];
sx q[2];
rz(-0.0068181097) q[2];
sx q[2];
rz(2.5799694) q[2];
rz(3.0567452) q[3];
sx q[3];
rz(-3.1361129) q[3];
sx q[3];
rz(-0.00076278846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4168004) q[0];
sx q[0];
rz(-0.264072) q[0];
sx q[0];
rz(-0.52077878) q[0];
rz(-0.15781038) q[1];
sx q[1];
rz(-2.4746555) q[1];
sx q[1];
rz(-3.0657943) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28215274) q[0];
sx q[0];
rz(-1.4569868) q[0];
sx q[0];
rz(2.5591617) q[0];
rz(-1.5694322) q[2];
sx q[2];
rz(-1.5701541) q[2];
sx q[2];
rz(-0.13880348) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0204649) q[1];
sx q[1];
rz(-1.9473828) q[1];
sx q[1];
rz(0.15213206) q[1];
rz(-0.10785477) q[3];
sx q[3];
rz(-1.1198992) q[3];
sx q[3];
rz(2.0863365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6303404) q[2];
sx q[2];
rz(-0.016009089) q[2];
sx q[2];
rz(1.7208257) q[2];
rz(-0.00682791) q[3];
sx q[3];
rz(-0.029191645) q[3];
sx q[3];
rz(1.6543057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0066978) q[0];
sx q[0];
rz(-0.20773523) q[0];
sx q[0];
rz(-2.8790706) q[0];
rz(0.94995704) q[1];
sx q[1];
rz(-0.078153178) q[1];
sx q[1];
rz(-0.21771678) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1658459) q[0];
sx q[0];
rz(-1.4474156) q[0];
sx q[0];
rz(-1.1749772) q[0];
rz(-pi) q[1];
rz(1.4870013) q[2];
sx q[2];
rz(-1.4838654) q[2];
sx q[2];
rz(-3.0695335) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7568183) q[1];
sx q[1];
rz(-2.9692602) q[1];
sx q[1];
rz(-1.5832572) q[1];
rz(-1.6895377) q[3];
sx q[3];
rz(-2.6104527) q[3];
sx q[3];
rz(2.0599928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.21506423) q[2];
sx q[2];
rz(-3.1319517) q[2];
sx q[2];
rz(-0.24791524) q[2];
rz(-1.3748112) q[3];
sx q[3];
rz(-0.050516613) q[3];
sx q[3];
rz(1.4352528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0404102) q[0];
sx q[0];
rz(-1.269416) q[0];
sx q[0];
rz(2.5293479) q[0];
rz(-2.9712037) q[1];
sx q[1];
rz(-0.082516106) q[1];
sx q[1];
rz(-1.6498529) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1338324) q[0];
sx q[0];
rz(-2.2585224) q[0];
sx q[0];
rz(0.69761116) q[0];
rz(-pi) q[1];
rz(3.1270624) q[2];
sx q[2];
rz(-1.5666012) q[2];
sx q[2];
rz(2.1452877) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1759169) q[1];
sx q[1];
rz(-1.8634444) q[1];
sx q[1];
rz(-1.6403557) q[1];
x q[2];
rz(1.2475523) q[3];
sx q[3];
rz(-1.3969143) q[3];
sx q[3];
rz(-1.6798879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.49543574) q[2];
sx q[2];
rz(-0.010224552) q[2];
sx q[2];
rz(2.9025027) q[2];
rz(-0.80860364) q[3];
sx q[3];
rz(-3.1294398) q[3];
sx q[3];
rz(-1.808572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7072153) q[0];
sx q[0];
rz(-0.093952976) q[0];
sx q[0];
rz(-2.3648426) q[0];
rz(0.010919318) q[1];
sx q[1];
rz(-0.24591406) q[1];
sx q[1];
rz(1.6687261) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68395972) q[0];
sx q[0];
rz(-1.3575956) q[0];
sx q[0];
rz(-0.85473608) q[0];
x q[1];
rz(3.1180834) q[2];
sx q[2];
rz(-1.5644367) q[2];
sx q[2];
rz(-1.3964749) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5653021) q[1];
sx q[1];
rz(-1.6716752) q[1];
sx q[1];
rz(-2.3957344) q[1];
rz(-pi) q[2];
rz(1.3192572) q[3];
sx q[3];
rz(-0.29485475) q[3];
sx q[3];
rz(-2.7406462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1482859) q[2];
sx q[2];
rz(-3.1243656) q[2];
sx q[2];
rz(-0.41603184) q[2];
rz(-2.3038583) q[3];
sx q[3];
rz(-0.13844027) q[3];
sx q[3];
rz(-1.4778888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1763022) q[0];
sx q[0];
rz(-0.039881341) q[0];
sx q[0];
rz(-2.1862929) q[0];
rz(-1.2166066) q[1];
sx q[1];
rz(-2.80426) q[1];
sx q[1];
rz(2.7359656) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55014706) q[0];
sx q[0];
rz(-2.2713775) q[0];
sx q[0];
rz(-1.1484543) q[0];
rz(1.4816059) q[2];
sx q[2];
rz(-1.729166) q[2];
sx q[2];
rz(-1.0430481) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4279528) q[1];
sx q[1];
rz(-1.7312839) q[1];
sx q[1];
rz(-1.6482216) q[1];
x q[2];
rz(-2.9285223) q[3];
sx q[3];
rz(-1.4858507) q[3];
sx q[3];
rz(-1.5379048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8568628) q[2];
sx q[2];
rz(-0.036981985) q[2];
sx q[2];
rz(2.3178318) q[2];
rz(3.01037) q[3];
sx q[3];
rz(-2.8516912) q[3];
sx q[3];
rz(-0.32740617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4284978) q[0];
sx q[0];
rz(-2.9976124) q[0];
sx q[0];
rz(-2.4289828) q[0];
rz(2.5932942) q[1];
sx q[1];
rz(-2.8916841) q[1];
sx q[1];
rz(2.9329494) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31748235) q[0];
sx q[0];
rz(-1.1718996) q[0];
sx q[0];
rz(2.0092416) q[0];
rz(-pi) q[1];
rz(1.0717355) q[2];
sx q[2];
rz(-0.03943561) q[2];
sx q[2];
rz(-2.5372504) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.123173) q[1];
sx q[1];
rz(-0.55572617) q[1];
sx q[1];
rz(1.5798142) q[1];
rz(-2.3712158) q[3];
sx q[3];
rz(-2.340014) q[3];
sx q[3];
rz(2.1397391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6930406) q[2];
sx q[2];
rz(-0.081144944) q[2];
sx q[2];
rz(0.21716675) q[2];
rz(-2.7740357) q[3];
sx q[3];
rz(-0.03438545) q[3];
sx q[3];
rz(-2.1448081) q[3];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9656068) q[0];
sx q[0];
rz(-0.088986926) q[0];
sx q[0];
rz(-2.9678645) q[0];
rz(1.732775) q[1];
sx q[1];
rz(-1.6830187) q[1];
sx q[1];
rz(-1.6857612) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.78409) q[0];
sx q[0];
rz(-1.0561868) q[0];
sx q[0];
rz(1.3823439) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.13952026) q[2];
sx q[2];
rz(-1.7762842) q[2];
sx q[2];
rz(0.71256283) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.10911726) q[1];
sx q[1];
rz(-0.87504234) q[1];
sx q[1];
rz(-1.7233707) q[1];
rz(1.1697269) q[3];
sx q[3];
rz(-0.26905832) q[3];
sx q[3];
rz(-0.45462409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9659861) q[2];
sx q[2];
rz(-0.49314988) q[2];
sx q[2];
rz(1.7612339) q[2];
rz(0.68584758) q[3];
sx q[3];
rz(-3.139747) q[3];
sx q[3];
rz(2.4628911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7420409) q[0];
sx q[0];
rz(-2.4492332) q[0];
sx q[0];
rz(-1.5204182) q[0];
rz(1.5581268) q[1];
sx q[1];
rz(-1.5065267) q[1];
sx q[1];
rz(-2.9361257) q[1];
rz(1.4420527) q[2];
sx q[2];
rz(-3.0151571) q[2];
sx q[2];
rz(-2.7966316) q[2];
rz(-1.2711502) q[3];
sx q[3];
rz(-1.4574403) q[3];
sx q[3];
rz(1.311655) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
