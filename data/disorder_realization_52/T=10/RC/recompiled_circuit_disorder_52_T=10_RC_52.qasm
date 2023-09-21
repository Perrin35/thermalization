OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.62587005) q[0];
sx q[0];
rz(-2.5928901) q[0];
sx q[0];
rz(0.8843511) q[0];
rz(-1.7110775) q[1];
sx q[1];
rz(-0.95354748) q[1];
sx q[1];
rz(1.6391099) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26412548) q[0];
sx q[0];
rz(-1.929951) q[0];
sx q[0];
rz(2.9496664) q[0];
x q[1];
rz(0.56514481) q[2];
sx q[2];
rz(-1.5752715) q[2];
sx q[2];
rz(0.43484136) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1391746) q[1];
sx q[1];
rz(-2.8867509) q[1];
sx q[1];
rz(-2.8789218) q[1];
rz(0.86125918) q[3];
sx q[3];
rz(-1.2951295) q[3];
sx q[3];
rz(-1.1112801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3551336) q[2];
sx q[2];
rz(-0.81374514) q[2];
sx q[2];
rz(-2.4856429) q[2];
rz(1.9338699) q[3];
sx q[3];
rz(-1.175712) q[3];
sx q[3];
rz(2.1470127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(1.8939963) q[0];
sx q[0];
rz(-0.42034724) q[0];
sx q[0];
rz(-2.7080652) q[0];
rz(2.9128089) q[1];
sx q[1];
rz(-2.7119535) q[1];
sx q[1];
rz(-3.1343592) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1996961) q[0];
sx q[0];
rz(-1.4775839) q[0];
sx q[0];
rz(-1.1774363) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2556633) q[2];
sx q[2];
rz(-1.2909781) q[2];
sx q[2];
rz(1.8880106) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0011598) q[1];
sx q[1];
rz(-1.7908485) q[1];
sx q[1];
rz(-0.37186719) q[1];
rz(2.238027) q[3];
sx q[3];
rz(-1.1840608) q[3];
sx q[3];
rz(0.29153338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6146415) q[2];
sx q[2];
rz(-2.3336637) q[2];
sx q[2];
rz(2.4439404) q[2];
rz(3.0200322) q[3];
sx q[3];
rz(-1.9024885) q[3];
sx q[3];
rz(0.30383032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6737297) q[0];
sx q[0];
rz(-2.2255852) q[0];
sx q[0];
rz(-1.3695705) q[0];
rz(-1.2415775) q[1];
sx q[1];
rz(-1.7280271) q[1];
sx q[1];
rz(-0.26161584) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81944377) q[0];
sx q[0];
rz(-1.561164) q[0];
sx q[0];
rz(1.5887512) q[0];
rz(-pi) q[1];
rz(-1.4357655) q[2];
sx q[2];
rz(-0.74880744) q[2];
sx q[2];
rz(0.43040648) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2178104) q[1];
sx q[1];
rz(-2.3465956) q[1];
sx q[1];
rz(-1.7407106) q[1];
x q[2];
rz(3.0858585) q[3];
sx q[3];
rz(-0.51364964) q[3];
sx q[3];
rz(-1.0518215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6802784) q[2];
sx q[2];
rz(-1.6363982) q[2];
sx q[2];
rz(0.20351163) q[2];
rz(-0.92173785) q[3];
sx q[3];
rz(-1.2676055) q[3];
sx q[3];
rz(2.8620499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1638284) q[0];
sx q[0];
rz(-1.5777359) q[0];
sx q[0];
rz(1.5699566) q[0];
rz(-1.0034026) q[1];
sx q[1];
rz(-1.3137716) q[1];
sx q[1];
rz(1.2483695) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49444775) q[0];
sx q[0];
rz(-1.6361423) q[0];
sx q[0];
rz(-1.4739743) q[0];
rz(-pi) q[1];
rz(1.7522881) q[2];
sx q[2];
rz(-1.776812) q[2];
sx q[2];
rz(-2.4455621) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.26989386) q[1];
sx q[1];
rz(-1.4596241) q[1];
sx q[1];
rz(-2.486869) q[1];
x q[2];
rz(-2.8473833) q[3];
sx q[3];
rz(-0.74511408) q[3];
sx q[3];
rz(-1.7780768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1127597) q[2];
sx q[2];
rz(-1.8303822) q[2];
sx q[2];
rz(-2.3045585) q[2];
rz(1.2083496) q[3];
sx q[3];
rz(-1.874606) q[3];
sx q[3];
rz(2.3560431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1355302) q[0];
sx q[0];
rz(-2.0879789) q[0];
sx q[0];
rz(-0.77520448) q[0];
rz(2.7397621) q[1];
sx q[1];
rz(-2.1907175) q[1];
sx q[1];
rz(-2.2391589) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9426614) q[0];
sx q[0];
rz(-0.64093243) q[0];
sx q[0];
rz(-0.076365691) q[0];
rz(-pi) q[1];
rz(2.6685153) q[2];
sx q[2];
rz(-2.6396751) q[2];
sx q[2];
rz(2.9830473) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9296226) q[1];
sx q[1];
rz(-2.3096482) q[1];
sx q[1];
rz(-2.3421939) q[1];
rz(-3.1387572) q[3];
sx q[3];
rz(-1.8875202) q[3];
sx q[3];
rz(-1.0032652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.19501413) q[2];
sx q[2];
rz(-2.6533551) q[2];
sx q[2];
rz(1.9449332) q[2];
rz(1.442391) q[3];
sx q[3];
rz(-1.2714352) q[3];
sx q[3];
rz(-1.6825914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3301795) q[0];
sx q[0];
rz(-0.17689642) q[0];
sx q[0];
rz(-0.50317558) q[0];
rz(1.4563837) q[1];
sx q[1];
rz(-1.0738942) q[1];
sx q[1];
rz(-0.20176372) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1267032) q[0];
sx q[0];
rz(-3.0198583) q[0];
sx q[0];
rz(0.99862167) q[0];
rz(-pi) q[1];
rz(-0.47735881) q[2];
sx q[2];
rz(-1.1672033) q[2];
sx q[2];
rz(-1.8237643) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4469874) q[1];
sx q[1];
rz(-1.1519377) q[1];
sx q[1];
rz(0.5254196) q[1];
rz(-pi) q[2];
rz(-1.1062578) q[3];
sx q[3];
rz(-2.2482276) q[3];
sx q[3];
rz(-2.5708452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9266944) q[2];
sx q[2];
rz(-1.639067) q[2];
sx q[2];
rz(-2.8743437) q[2];
rz(2.3184508) q[3];
sx q[3];
rz(-3.0624793) q[3];
sx q[3];
rz(0.87583035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8478407) q[0];
sx q[0];
rz(-0.1593312) q[0];
sx q[0];
rz(-3.0840432) q[0];
rz(1.6607704) q[1];
sx q[1];
rz(-1.7819504) q[1];
sx q[1];
rz(-0.94271359) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38521117) q[0];
sx q[0];
rz(-0.74059534) q[0];
sx q[0];
rz(-0.60473196) q[0];
rz(-pi) q[1];
rz(2.9808688) q[2];
sx q[2];
rz(-1.3458999) q[2];
sx q[2];
rz(2.7149534) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.47989935) q[1];
sx q[1];
rz(-2.6091895) q[1];
sx q[1];
rz(-0.9049306) q[1];
x q[2];
rz(1.1293291) q[3];
sx q[3];
rz(-2.5848476) q[3];
sx q[3];
rz(2.7260775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1192347) q[2];
sx q[2];
rz(-2.0998349) q[2];
sx q[2];
rz(-0.38267246) q[2];
rz(2.102397) q[3];
sx q[3];
rz(-1.8363876) q[3];
sx q[3];
rz(-0.97810811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4246178) q[0];
sx q[0];
rz(-0.096352339) q[0];
sx q[0];
rz(-0.27012816) q[0];
rz(-0.62942901) q[1];
sx q[1];
rz(-2.4286178) q[1];
sx q[1];
rz(0.28392917) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2487508) q[0];
sx q[0];
rz(-2.1584956) q[0];
sx q[0];
rz(-0.52389223) q[0];
x q[1];
rz(-2.5416652) q[2];
sx q[2];
rz(-1.6919961) q[2];
sx q[2];
rz(-1.3231414) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.200951) q[1];
sx q[1];
rz(-2.0307699) q[1];
sx q[1];
rz(0.56401395) q[1];
rz(-pi) q[2];
rz(1.5408526) q[3];
sx q[3];
rz(-2.2303914) q[3];
sx q[3];
rz(-2.3823882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.55390629) q[2];
sx q[2];
rz(-2.9710785) q[2];
sx q[2];
rz(-1.2109057) q[2];
rz(-2.8816913) q[3];
sx q[3];
rz(-2.5172958) q[3];
sx q[3];
rz(0.62817812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.07638409) q[0];
sx q[0];
rz(-0.56088352) q[0];
sx q[0];
rz(2.912345) q[0];
rz(-0.30300888) q[1];
sx q[1];
rz(-1.3906994) q[1];
sx q[1];
rz(-1.4607666) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0871353) q[0];
sx q[0];
rz(-1.5532877) q[0];
sx q[0];
rz(-0.28355916) q[0];
rz(-pi) q[1];
x q[1];
rz(0.37947189) q[2];
sx q[2];
rz(-1.4037637) q[2];
sx q[2];
rz(0.52969474) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8742121) q[1];
sx q[1];
rz(-2.0955718) q[1];
sx q[1];
rz(-1.5966148) q[1];
x q[2];
rz(-1.1838412) q[3];
sx q[3];
rz(-2.317252) q[3];
sx q[3];
rz(1.6798306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.71172697) q[2];
sx q[2];
rz(-1.2167598) q[2];
sx q[2];
rz(-2.4460068) q[2];
rz(0.43186489) q[3];
sx q[3];
rz(-2.6769107) q[3];
sx q[3];
rz(-0.7152043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.573695) q[0];
sx q[0];
rz(-1.2117813) q[0];
sx q[0];
rz(0.33690548) q[0];
rz(0.20740549) q[1];
sx q[1];
rz(-1.0284871) q[1];
sx q[1];
rz(-2.7609603) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1487409) q[0];
sx q[0];
rz(-1.3694351) q[0];
sx q[0];
rz(-3.0701748) q[0];
rz(-pi) q[1];
rz(-0.98999087) q[2];
sx q[2];
rz(-0.91772807) q[2];
sx q[2];
rz(0.68819118) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.759728) q[1];
sx q[1];
rz(-1.3761531) q[1];
sx q[1];
rz(2.0774283) q[1];
rz(-pi) q[2];
rz(-0.59488876) q[3];
sx q[3];
rz(-0.66685646) q[3];
sx q[3];
rz(-2.8951333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.0011065817) q[2];
sx q[2];
rz(-1.1662741) q[2];
sx q[2];
rz(-2.005119) q[2];
rz(-0.040955695) q[3];
sx q[3];
rz(-2.332873) q[3];
sx q[3];
rz(1.3142746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.8052335) q[0];
sx q[0];
rz(-0.90538607) q[0];
sx q[0];
rz(-0.3120099) q[0];
rz(1.0271172) q[1];
sx q[1];
rz(-1.849091) q[1];
sx q[1];
rz(-1.0277933) q[1];
rz(1.3409875) q[2];
sx q[2];
rz(-1.3290214) q[2];
sx q[2];
rz(2.7711836) q[2];
rz(-2.9410578) q[3];
sx q[3];
rz(-2.0024828) q[3];
sx q[3];
rz(1.9867292) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];