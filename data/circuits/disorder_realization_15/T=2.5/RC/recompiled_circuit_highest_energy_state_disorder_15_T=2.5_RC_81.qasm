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
rz(0.71969405) q[0];
rz(2.5201058) q[1];
sx q[1];
rz(-1.7401594) q[1];
sx q[1];
rz(-3.016234) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1272362) q[0];
sx q[0];
rz(-1.171021) q[0];
sx q[0];
rz(0.46443224) q[0];
rz(-0.28858443) q[2];
sx q[2];
rz(-1.4800396) q[2];
sx q[2];
rz(-1.2707658) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5125864) q[1];
sx q[1];
rz(-0.0014481469) q[1];
sx q[1];
rz(1.8667029) q[1];
rz(-0.35644021) q[3];
sx q[3];
rz(-2.2619777) q[3];
sx q[3];
rz(-0.82572848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3434992) q[2];
sx q[2];
rz(-0.40830475) q[2];
sx q[2];
rz(-2.2957392) q[2];
rz(-2.3503303) q[3];
sx q[3];
rz(-3.1281804) q[3];
sx q[3];
rz(-0.066789269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4132408) q[0];
sx q[0];
rz(-2.6744196) q[0];
sx q[0];
rz(3.0979284) q[0];
rz(1.5664258) q[1];
sx q[1];
rz(-1.3730201) q[1];
sx q[1];
rz(-1.498819) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18066809) q[0];
sx q[0];
rz(-2.5378413) q[0];
sx q[0];
rz(-1.5636428) q[0];
x q[1];
rz(2.1421932) q[2];
sx q[2];
rz(-1.5816814) q[2];
sx q[2];
rz(1.581426) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.39061645) q[1];
sx q[1];
rz(-1.6469139) q[1];
sx q[1];
rz(1.5396495) q[1];
rz(-pi) q[2];
rz(1.6632027) q[3];
sx q[3];
rz(-2.795821) q[3];
sx q[3];
rz(-1.2279296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0955536) q[2];
sx q[2];
rz(-0.150103) q[2];
sx q[2];
rz(-2.6138439) q[2];
rz(-2.2857417) q[3];
sx q[3];
rz(-0.0015043613) q[3];
sx q[3];
rz(-1.2214448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39831487) q[0];
sx q[0];
rz(-2.1687431) q[0];
sx q[0];
rz(-1.1463746) q[0];
rz(1.3829117) q[1];
sx q[1];
rz(-2.8490366) q[1];
sx q[1];
rz(0.10398277) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71767112) q[0];
sx q[0];
rz(-1.648745) q[0];
sx q[0];
rz(2.9162327) q[0];
rz(-pi) q[1];
rz(-1.632003) q[2];
sx q[2];
rz(-0.27282295) q[2];
sx q[2];
rz(0.83977985) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5755363) q[1];
sx q[1];
rz(-2.3588099) q[1];
sx q[1];
rz(-2.8498883) q[1];
rz(-pi) q[2];
rz(1.6802877) q[3];
sx q[3];
rz(-1.8598622) q[3];
sx q[3];
rz(1.3407202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.18342239) q[2];
sx q[2];
rz(-3.1347745) q[2];
sx q[2];
rz(-2.5799694) q[2];
rz(-3.0567452) q[3];
sx q[3];
rz(-0.0054797879) q[3];
sx q[3];
rz(-0.00076278846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7247923) q[0];
sx q[0];
rz(-2.8775207) q[0];
sx q[0];
rz(2.6208139) q[0];
rz(2.9837823) q[1];
sx q[1];
rz(-2.4746555) q[1];
sx q[1];
rz(0.075798362) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8594399) q[0];
sx q[0];
rz(-1.6846058) q[0];
sx q[0];
rz(0.58243097) q[0];
x q[1];
rz(1.5694322) q[2];
sx q[2];
rz(-1.5701541) q[2];
sx q[2];
rz(0.13880348) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0204649) q[1];
sx q[1];
rz(-1.1942099) q[1];
sx q[1];
rz(-0.15213206) q[1];
rz(1.1176093) q[3];
sx q[3];
rz(-1.6678358) q[3];
sx q[3];
rz(0.56268855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6303404) q[2];
sx q[2];
rz(-3.1255836) q[2];
sx q[2];
rz(1.7208257) q[2];
rz(3.1347647) q[3];
sx q[3];
rz(-3.112401) q[3];
sx q[3];
rz(-1.6543057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1348949) q[0];
sx q[0];
rz(-2.9338574) q[0];
sx q[0];
rz(2.8790706) q[0];
rz(0.94995704) q[1];
sx q[1];
rz(-3.0634395) q[1];
sx q[1];
rz(0.21771678) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1658459) q[0];
sx q[0];
rz(-1.4474156) q[0];
sx q[0];
rz(1.9666155) q[0];
rz(-pi) q[1];
rz(1.6545914) q[2];
sx q[2];
rz(-1.6577273) q[2];
sx q[2];
rz(0.072059137) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.1982983) q[1];
sx q[1];
rz(-1.5729331) q[1];
sx q[1];
rz(1.7431156) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0721139) q[3];
sx q[3];
rz(-1.0437878) q[3];
sx q[3];
rz(-1.219092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9265284) q[2];
sx q[2];
rz(-0.0096409163) q[2];
sx q[2];
rz(2.8936774) q[2];
rz(-1.3748112) q[3];
sx q[3];
rz(-0.050516613) q[3];
sx q[3];
rz(1.4352528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1011825) q[0];
sx q[0];
rz(-1.269416) q[0];
sx q[0];
rz(-2.5293479) q[0];
rz(0.170389) q[1];
sx q[1];
rz(-0.082516106) q[1];
sx q[1];
rz(1.4917397) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9303974) q[0];
sx q[0];
rz(-2.204551) q[0];
sx q[0];
rz(2.2344345) q[0];
rz(1.5666008) q[2];
sx q[2];
rz(-1.5562662) q[2];
sx q[2];
rz(2.5670403) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.93894847) q[1];
sx q[1];
rz(-0.30057014) q[1];
sx q[1];
rz(-2.9148797) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2475523) q[3];
sx q[3];
rz(-1.7446784) q[3];
sx q[3];
rz(1.4617048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6461569) q[2];
sx q[2];
rz(-3.1313681) q[2];
sx q[2];
rz(-0.23908991) q[2];
rz(-0.80860364) q[3];
sx q[3];
rz(-0.012152823) q[3];
sx q[3];
rz(1.808572) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4343774) q[0];
sx q[0];
rz(-0.093952976) q[0];
sx q[0];
rz(-0.77675003) q[0];
rz(-0.010919318) q[1];
sx q[1];
rz(-0.24591406) q[1];
sx q[1];
rz(-1.6687261) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70476711) q[0];
sx q[0];
rz(-2.2673207) q[0];
sx q[0];
rz(-2.8621307) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1180834) q[2];
sx q[2];
rz(-1.5644367) q[2];
sx q[2];
rz(1.7451177) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.90172988) q[1];
sx q[1];
rz(-0.82962501) q[1];
sx q[1];
rz(1.4338486) q[1];
x q[2];
rz(-0.075447791) q[3];
sx q[3];
rz(-1.2854854) q[3];
sx q[3];
rz(2.4782654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9933068) q[2];
sx q[2];
rz(-0.017227087) q[2];
sx q[2];
rz(0.41603184) q[2];
rz(-0.83773437) q[3];
sx q[3];
rz(-0.13844027) q[3];
sx q[3];
rz(-1.6637038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96529043) q[0];
sx q[0];
rz(-0.039881341) q[0];
sx q[0];
rz(-2.1862929) q[0];
rz(-1.924986) q[1];
sx q[1];
rz(-0.33733264) q[1];
sx q[1];
rz(-0.40562707) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5914456) q[0];
sx q[0];
rz(-0.87021512) q[0];
sx q[0];
rz(-1.9931384) q[0];
rz(-pi) q[1];
x q[1];
rz(0.15899105) q[2];
sx q[2];
rz(-1.6588677) q[2];
sx q[2];
rz(-2.627947) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1656087) q[1];
sx q[1];
rz(-2.9635495) q[1];
sx q[1];
rz(2.6958532) q[1];
rz(-pi) q[2];
rz(1.4838951) q[3];
sx q[3];
rz(-1.358506) q[3];
sx q[3];
rz(0.014537285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8568628) q[2];
sx q[2];
rz(-3.1046107) q[2];
sx q[2];
rz(2.3178318) q[2];
rz(3.01037) q[3];
sx q[3];
rz(-0.28990144) q[3];
sx q[3];
rz(0.32740617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4284978) q[0];
sx q[0];
rz(-2.9976124) q[0];
sx q[0];
rz(0.71260989) q[0];
rz(2.5932942) q[1];
sx q[1];
rz(-2.8916841) q[1];
sx q[1];
rz(-0.20864329) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0684194) q[0];
sx q[0];
rz(-1.9727339) q[0];
sx q[0];
rz(-2.7059024) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6054262) q[2];
sx q[2];
rz(-1.5519262) q[2];
sx q[2];
rz(-0.46771995) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.1125579) q[1];
sx q[1];
rz(-1.0150954) q[1];
sx q[1];
rz(-3.1359926) q[1];
x q[2];
rz(0.63788148) q[3];
sx q[3];
rz(-1.0468093) q[3];
sx q[3];
rz(1.978911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6930406) q[2];
sx q[2];
rz(-0.081144944) q[2];
sx q[2];
rz(0.21716675) q[2];
rz(-0.36755696) q[3];
sx q[3];
rz(-0.03438545) q[3];
sx q[3];
rz(-0.99678451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9656068) q[0];
sx q[0];
rz(-3.0526057) q[0];
sx q[0];
rz(-2.9678645) q[0];
rz(1.732775) q[1];
sx q[1];
rz(-1.6830187) q[1];
sx q[1];
rz(-1.6857612) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98783606) q[0];
sx q[0];
rz(-2.5964886) q[0];
sx q[0];
rz(-0.31994064) q[0];
rz(-3.0020724) q[2];
sx q[2];
rz(-1.7762842) q[2];
sx q[2];
rz(-0.71256283) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5816702) q[1];
sx q[1];
rz(-1.6877203) q[1];
sx q[1];
rz(2.440084) q[1];
rz(-pi) q[2];
rz(-3.0343541) q[3];
sx q[3];
rz(-1.8180366) q[3];
sx q[3];
rz(0.040217248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9659861) q[2];
sx q[2];
rz(-0.49314988) q[2];
sx q[2];
rz(1.7612339) q[2];
rz(-0.68584758) q[3];
sx q[3];
rz(-3.139747) q[3];
sx q[3];
rz(0.67870158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3995517) q[0];
sx q[0];
rz(-2.4492332) q[0];
sx q[0];
rz(-1.5204182) q[0];
rz(-1.5581268) q[1];
sx q[1];
rz(-1.635066) q[1];
sx q[1];
rz(0.20546694) q[1];
rz(-1.6995399) q[2];
sx q[2];
rz(-3.0151571) q[2];
sx q[2];
rz(-2.7966316) q[2];
rz(-1.9388922) q[3];
sx q[3];
rz(-2.8218269) q[3];
sx q[3];
rz(0.091824986) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
