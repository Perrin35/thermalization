OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6686749) q[0];
sx q[0];
rz(-0.023107419) q[0];
sx q[0];
rz(0.90149108) q[0];
rz(-1.787552) q[1];
sx q[1];
rz(-1.6156337) q[1];
sx q[1];
rz(-1.807133) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62518277) q[0];
sx q[0];
rz(-1.7889272) q[0];
sx q[0];
rz(1.6263917) q[0];
rz(-pi) q[1];
x q[1];
rz(2.790144) q[2];
sx q[2];
rz(-0.16021591) q[2];
sx q[2];
rz(1.7376228) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.035369594) q[1];
sx q[1];
rz(-1.5179885) q[1];
sx q[1];
rz(1.04671) q[1];
x q[2];
rz(0.18906148) q[3];
sx q[3];
rz(-0.24654085) q[3];
sx q[3];
rz(0.93545675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2287075) q[2];
sx q[2];
rz(-3.0444453) q[2];
sx q[2];
rz(0.65931064) q[2];
rz(-0.77624503) q[3];
sx q[3];
rz(-3.1228437) q[3];
sx q[3];
rz(-0.68500486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9126251) q[0];
sx q[0];
rz(-1.2094867) q[0];
sx q[0];
rz(1.1932766) q[0];
rz(0.044366447) q[1];
sx q[1];
rz(-0.012244789) q[1];
sx q[1];
rz(0.22656974) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3831155) q[0];
sx q[0];
rz(-1.5673914) q[0];
sx q[0];
rz(-0.0016511516) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.042828538) q[2];
sx q[2];
rz(-0.37938243) q[2];
sx q[2];
rz(-1.6216175) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4227155) q[1];
sx q[1];
rz(-1.1402292) q[1];
sx q[1];
rz(-1.6490472) q[1];
rz(-2.2983589) q[3];
sx q[3];
rz(-1.4012464) q[3];
sx q[3];
rz(0.42593004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.442753) q[2];
sx q[2];
rz(-1.5719465) q[2];
sx q[2];
rz(-1.5296096) q[2];
rz(2.2085341) q[3];
sx q[3];
rz(-1.6589087) q[3];
sx q[3];
rz(-0.23811594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17732349) q[0];
sx q[0];
rz(-3.1260335) q[0];
sx q[0];
rz(2.9413057) q[0];
rz(-0.00042032584) q[1];
sx q[1];
rz(-2.2047408) q[1];
sx q[1];
rz(0.013484152) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6967675) q[0];
sx q[0];
rz(-2.7119141) q[0];
sx q[0];
rz(1.7264771) q[0];
rz(-pi) q[1];
x q[1];
rz(0.066670074) q[2];
sx q[2];
rz(-1.5278897) q[2];
sx q[2];
rz(1.5868452) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2993907) q[1];
sx q[1];
rz(-2.9889571) q[1];
sx q[1];
rz(-2.045689) q[1];
rz(-pi) q[2];
rz(-0.30463574) q[3];
sx q[3];
rz(-0.17028642) q[3];
sx q[3];
rz(2.0326322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9424332) q[2];
sx q[2];
rz(-1.5853256) q[2];
sx q[2];
rz(-1.5996492) q[2];
rz(1.0017627) q[3];
sx q[3];
rz(-2.9250513) q[3];
sx q[3];
rz(0.17222968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7672985) q[0];
sx q[0];
rz(-0.16574398) q[0];
sx q[0];
rz(-2.7941008) q[0];
rz(-0.63684288) q[1];
sx q[1];
rz(-3.1359735) q[1];
sx q[1];
rz(1.9510795) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74032846) q[0];
sx q[0];
rz(-1.6303827) q[0];
sx q[0];
rz(-1.4387095) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1338155) q[2];
sx q[2];
rz(-1.6396804) q[2];
sx q[2];
rz(-1.5635671) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1381629) q[1];
sx q[1];
rz(-2.3210879) q[1];
sx q[1];
rz(1.8292887) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.52996323) q[3];
sx q[3];
rz(-0.41614446) q[3];
sx q[3];
rz(-2.346739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5570598) q[2];
sx q[2];
rz(-0.0396885) q[2];
sx q[2];
rz(-1.7812799) q[2];
rz(1.4761866) q[3];
sx q[3];
rz(-1.5740266) q[3];
sx q[3];
rz(-0.40569693) q[3];
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
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0535468) q[0];
sx q[0];
rz(-0.73943728) q[0];
sx q[0];
rz(0.0060225688) q[0];
rz(1.4206403) q[1];
sx q[1];
rz(-0.050844897) q[1];
sx q[1];
rz(0.086070148) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6688441) q[0];
sx q[0];
rz(-1.5540344) q[0];
sx q[0];
rz(-3.0344323) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.055069607) q[2];
sx q[2];
rz(-1.5968284) q[2];
sx q[2];
rz(-0.30773417) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3427375) q[1];
sx q[1];
rz(-1.6621843) q[1];
sx q[1];
rz(-1.6080086) q[1];
rz(-pi) q[2];
rz(-3.0962493) q[3];
sx q[3];
rz(-1.5102949) q[3];
sx q[3];
rz(0.44024548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.033279557) q[2];
sx q[2];
rz(-2.4187708) q[2];
sx q[2];
rz(1.7576199) q[2];
rz(-0.39517394) q[3];
sx q[3];
rz(-3.0925909) q[3];
sx q[3];
rz(-1.1672195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60642099) q[0];
sx q[0];
rz(-2.9223154) q[0];
sx q[0];
rz(-1.0004591) q[0];
rz(2.4011627) q[1];
sx q[1];
rz(-0.43559566) q[1];
sx q[1];
rz(2.7620517) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1943037) q[0];
sx q[0];
rz(-1.6662046) q[0];
sx q[0];
rz(-1.545115) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8431115) q[2];
sx q[2];
rz(-1.1109931) q[2];
sx q[2];
rz(-2.1260171) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4011456) q[1];
sx q[1];
rz(-1.5974495) q[1];
sx q[1];
rz(-1.0668287) q[1];
x q[2];
rz(-1.1218171) q[3];
sx q[3];
rz(-1.4662577) q[3];
sx q[3];
rz(3.0825101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.1257816) q[2];
sx q[2];
rz(-2.8534079) q[2];
sx q[2];
rz(-3.06456) q[2];
rz(-3.1214664) q[3];
sx q[3];
rz(-3.088981) q[3];
sx q[3];
rz(-2.3441815) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.035148419) q[0];
sx q[0];
rz(-0.012520944) q[0];
sx q[0];
rz(1.7092108) q[0];
rz(-2.8445981) q[1];
sx q[1];
rz(-0.15428267) q[1];
sx q[1];
rz(-3.0800379) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2765574) q[0];
sx q[0];
rz(-1.1662959) q[0];
sx q[0];
rz(-2.1767031) q[0];
x q[1];
rz(-0.02587895) q[2];
sx q[2];
rz(-2.9270083) q[2];
sx q[2];
rz(-1.9865004) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6541432) q[1];
sx q[1];
rz(-1.2795078) q[1];
sx q[1];
rz(1.6421516) q[1];
rz(-pi) q[2];
rz(0.43470862) q[3];
sx q[3];
rz(-1.604933) q[3];
sx q[3];
rz(-0.0061090547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.84294549) q[2];
sx q[2];
rz(-0.25235287) q[2];
sx q[2];
rz(-1.418815) q[2];
rz(1.336054) q[3];
sx q[3];
rz(-0.027848363) q[3];
sx q[3];
rz(1.7347887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1143188) q[0];
sx q[0];
rz(-2.424746) q[0];
sx q[0];
rz(-1.5007716) q[0];
rz(1.1227135) q[1];
sx q[1];
rz(-2.8148459) q[1];
sx q[1];
rz(1.8480802) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8282765) q[0];
sx q[0];
rz(-1.8819081) q[0];
sx q[0];
rz(0.88514741) q[0];
rz(-1.229153) q[2];
sx q[2];
rz(-0.90589303) q[2];
sx q[2];
rz(-1.5835012) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1652917) q[1];
sx q[1];
rz(-1.2982076) q[1];
sx q[1];
rz(1.7458245) q[1];
rz(-pi) q[2];
rz(-2.294306) q[3];
sx q[3];
rz(-2.1186817) q[3];
sx q[3];
rz(0.63036608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5844172) q[2];
sx q[2];
rz(-0.36089218) q[2];
sx q[2];
rz(-1.7822251) q[2];
rz(-2.6946097) q[3];
sx q[3];
rz(-0.042526571) q[3];
sx q[3];
rz(0.16714787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8166703) q[0];
sx q[0];
rz(-1.8055547) q[0];
sx q[0];
rz(-1.100612) q[0];
rz(2.0657516) q[1];
sx q[1];
rz(-2.4918719) q[1];
sx q[1];
rz(-0.67695391) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8146801) q[0];
sx q[0];
rz(-1.956358) q[0];
sx q[0];
rz(-1.7837693) q[0];
x q[1];
rz(0.1157427) q[2];
sx q[2];
rz(-1.6925768) q[2];
sx q[2];
rz(-2.1809721) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7266287) q[1];
sx q[1];
rz(-1.5708956) q[1];
sx q[1];
rz(-1.5703771) q[1];
rz(1.7221801) q[3];
sx q[3];
rz(-2.0213691) q[3];
sx q[3];
rz(-2.7344658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0000275) q[2];
sx q[2];
rz(-0.0024777369) q[2];
sx q[2];
rz(-0.72762093) q[2];
rz(1.8992807) q[3];
sx q[3];
rz(-3.1051903) q[3];
sx q[3];
rz(-1.9696994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.2321371) q[0];
sx q[0];
rz(-2.1196892) q[0];
sx q[0];
rz(0.68956462) q[0];
rz(1.6652971) q[1];
sx q[1];
rz(-0.25054014) q[1];
sx q[1];
rz(2.9857059) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4094226) q[0];
sx q[0];
rz(-0.66930938) q[0];
sx q[0];
rz(-0.21440345) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7047809) q[2];
sx q[2];
rz(-2.8653318) q[2];
sx q[2];
rz(-2.7154684) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2624596) q[1];
sx q[1];
rz(-3.1371318) q[1];
sx q[1];
rz(0.35426472) q[1];
x q[2];
rz(1.4149551) q[3];
sx q[3];
rz(-1.4771802) q[3];
sx q[3];
rz(-2.7839212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9547687) q[2];
sx q[2];
rz(-0.12207741) q[2];
sx q[2];
rz(-2.1464777) q[2];
rz(-2.9156445) q[3];
sx q[3];
rz(-0.048361691) q[3];
sx q[3];
rz(2.3154955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36289287) q[0];
sx q[0];
rz(-0.97465546) q[0];
sx q[0];
rz(1.4018651) q[0];
rz(1.4317935) q[1];
sx q[1];
rz(-1.2964389) q[1];
sx q[1];
rz(-2.5242205) q[1];
rz(-1.7079034) q[2];
sx q[2];
rz(-2.2513486) q[2];
sx q[2];
rz(0.90604102) q[2];
rz(-1.9598087) q[3];
sx q[3];
rz(-1.4445732) q[3];
sx q[3];
rz(-2.5635535) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
