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
rz(0.61007845) q[0];
sx q[0];
rz(-2.9116723) q[0];
sx q[0];
rz(2.4218986) q[0];
rz(2.5201058) q[1];
sx q[1];
rz(1.4014333) q[1];
sx q[1];
rz(9.2994193) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1272362) q[0];
sx q[0];
rz(-1.9705716) q[0];
sx q[0];
rz(-0.46443224) q[0];
rz(-0.30949952) q[2];
sx q[2];
rz(-0.30213854) q[2];
sx q[2];
rz(-3.1378656) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.6290063) q[1];
sx q[1];
rz(-0.0014481469) q[1];
sx q[1];
rz(1.2748898) q[1];
rz(2.2940346) q[3];
sx q[3];
rz(-1.8430018) q[3];
sx q[3];
rz(0.51203007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3434992) q[2];
sx q[2];
rz(-0.40830475) q[2];
sx q[2];
rz(0.84585345) q[2];
rz(-2.3503303) q[3];
sx q[3];
rz(-0.013412272) q[3];
sx q[3];
rz(0.066789269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4132408) q[0];
sx q[0];
rz(-2.6744196) q[0];
sx q[0];
rz(-0.043664232) q[0];
rz(-1.5664258) q[1];
sx q[1];
rz(-1.3730201) q[1];
sx q[1];
rz(1.498819) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3960172) q[0];
sx q[0];
rz(-1.5748576) q[0];
sx q[0];
rz(0.96705695) q[0];
x q[1];
rz(1.5506707) q[2];
sx q[2];
rz(-0.57148904) q[2];
sx q[2];
rz(3.1352941) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1825492) q[1];
sx q[1];
rz(-1.5397397) q[1];
sx q[1];
rz(0.076154321) q[1];
rz(-pi) q[2];
rz(1.915207) q[3];
sx q[3];
rz(-1.6020755) q[3];
sx q[3];
rz(0.42983228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0955536) q[2];
sx q[2];
rz(-0.150103) q[2];
sx q[2];
rz(0.52774876) q[2];
rz(0.85585099) q[3];
sx q[3];
rz(-3.1400883) q[3];
sx q[3];
rz(-1.9201479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7432778) q[0];
sx q[0];
rz(-2.1687431) q[0];
sx q[0];
rz(1.995218) q[0];
rz(1.3829117) q[1];
sx q[1];
rz(-0.29255602) q[1];
sx q[1];
rz(-0.10398277) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4239215) q[0];
sx q[0];
rz(-1.4928476) q[0];
sx q[0];
rz(-0.22535997) q[0];
rz(0.017113233) q[2];
sx q[2];
rz(-1.2984973) q[2];
sx q[2];
rz(-2.2382617) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.56605634) q[1];
sx q[1];
rz(-2.3588099) q[1];
sx q[1];
rz(-0.29170431) q[1];
x q[2];
rz(-1.461305) q[3];
sx q[3];
rz(-1.2817304) q[3];
sx q[3];
rz(-1.3407202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9581703) q[2];
sx q[2];
rz(-3.1347745) q[2];
sx q[2];
rz(0.56162322) q[2];
rz(3.0567452) q[3];
sx q[3];
rz(-0.0054797879) q[3];
sx q[3];
rz(-3.1408299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4168004) q[0];
sx q[0];
rz(-0.264072) q[0];
sx q[0];
rz(-0.52077878) q[0];
rz(0.15781038) q[1];
sx q[1];
rz(-0.66693711) q[1];
sx q[1];
rz(0.075798362) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4593845) q[0];
sx q[0];
rz(-0.59218107) q[0];
sx q[0];
rz(2.9367052) q[0];
x q[1];
rz(2.0108156) q[2];
sx q[2];
rz(-0.001507757) q[2];
sx q[2];
rz(-1.269581) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0204649) q[1];
sx q[1];
rz(-1.1942099) q[1];
sx q[1];
rz(-2.9894606) q[1];
x q[2];
rz(0.10785477) q[3];
sx q[3];
rz(-2.0216935) q[3];
sx q[3];
rz(-1.0552561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.51125222) q[2];
sx q[2];
rz(-3.1255836) q[2];
sx q[2];
rz(1.7208257) q[2];
rz(-0.00682791) q[3];
sx q[3];
rz(-3.112401) q[3];
sx q[3];
rz(-1.6543057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1348949) q[0];
sx q[0];
rz(-2.9338574) q[0];
sx q[0];
rz(-2.8790706) q[0];
rz(-0.94995704) q[1];
sx q[1];
rz(-3.0634395) q[1];
sx q[1];
rz(-0.21771678) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1658459) q[0];
sx q[0];
rz(-1.6941771) q[0];
sx q[0];
rz(1.1749772) q[0];
x q[1];
rz(0.08723549) q[2];
sx q[2];
rz(-1.4873184) q[2];
sx q[2];
rz(-1.491445) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7694665) q[1];
sx q[1];
rz(-1.3984774) q[1];
sx q[1];
rz(-3.1394238) q[1];
rz(-1.0427371) q[3];
sx q[3];
rz(-1.6308356) q[3];
sx q[3];
rz(2.754902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.21506423) q[2];
sx q[2];
rz(-3.1319517) q[2];
sx q[2];
rz(0.24791524) q[2];
rz(-1.7667814) q[3];
sx q[3];
rz(-3.091076) q[3];
sx q[3];
rz(-1.7063399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0404102) q[0];
sx q[0];
rz(-1.8721767) q[0];
sx q[0];
rz(-0.61224473) q[0];
rz(-2.9712037) q[1];
sx q[1];
rz(-3.0590765) q[1];
sx q[1];
rz(-1.4917397) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0077602) q[0];
sx q[0];
rz(-2.2585224) q[0];
sx q[0];
rz(2.4439815) q[0];
x q[1];
rz(0.014530226) q[2];
sx q[2];
rz(-1.5749914) q[2];
sx q[2];
rz(2.1452877) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1759169) q[1];
sx q[1];
rz(-1.8634444) q[1];
sx q[1];
rz(-1.6403557) q[1];
rz(-2.0759367) q[3];
sx q[3];
rz(-0.36559421) q[3];
sx q[3];
rz(2.5557809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7072153) q[0];
sx q[0];
rz(-0.093952976) q[0];
sx q[0];
rz(2.3648426) q[0];
rz(0.010919318) q[1];
sx q[1];
rz(-0.24591406) q[1];
sx q[1];
rz(1.6687261) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1253788) q[0];
sx q[0];
rz(-0.74170602) q[0];
sx q[0];
rz(-1.8893695) q[0];
rz(1.5771578) q[2];
sx q[2];
rz(-1.5943051) q[2];
sx q[2];
rz(2.9671217) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.57629055) q[1];
sx q[1];
rz(-1.6716752) q[1];
sx q[1];
rz(-2.3957344) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3192572) q[3];
sx q[3];
rz(-2.8467379) q[3];
sx q[3];
rz(0.40094646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96529043) q[0];
sx q[0];
rz(-3.1017113) q[0];
sx q[0];
rz(-2.1862929) q[0];
rz(-1.924986) q[1];
sx q[1];
rz(-2.80426) q[1];
sx q[1];
rz(-2.7359656) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3026311) q[0];
sx q[0];
rz(-1.8895188) q[0];
sx q[0];
rz(-2.3953954) q[0];
x q[1];
rz(2.9826016) q[2];
sx q[2];
rz(-1.6588677) q[2];
sx q[2];
rz(-0.51364567) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1656087) q[1];
sx q[1];
rz(-0.17804314) q[1];
sx q[1];
rz(-2.6958532) q[1];
x q[2];
rz(0.38281103) q[3];
sx q[3];
rz(-2.9124526) q[3];
sx q[3];
rz(-0.40661302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2847298) q[2];
sx q[2];
rz(-3.1046107) q[2];
sx q[2];
rz(0.82376087) q[2];
rz(-3.01037) q[3];
sx q[3];
rz(-0.28990144) q[3];
sx q[3];
rz(2.8141865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4284978) q[0];
sx q[0];
rz(-2.9976124) q[0];
sx q[0];
rz(-0.71260989) q[0];
rz(2.5932942) q[1];
sx q[1];
rz(-0.24990853) q[1];
sx q[1];
rz(0.20864329) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31748235) q[0];
sx q[0];
rz(-1.1718996) q[0];
sx q[0];
rz(-1.1323511) q[0];
rz(-2.0698572) q[2];
sx q[2];
rz(-3.102157) q[2];
sx q[2];
rz(-0.60434228) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5447158) q[1];
sx q[1];
rz(-1.5660389) q[1];
sx q[1];
rz(2.1265043) q[1];
rz(-pi) q[2];
rz(-0.94722469) q[3];
sx q[3];
rz(-2.1125018) q[3];
sx q[3];
rz(0.052963363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.448552) q[2];
sx q[2];
rz(-3.0604477) q[2];
sx q[2];
rz(-0.21716675) q[2];
rz(2.7740357) q[3];
sx q[3];
rz(-0.03438545) q[3];
sx q[3];
rz(2.1448081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9656068) q[0];
sx q[0];
rz(-3.0526057) q[0];
sx q[0];
rz(-0.17372818) q[0];
rz(-1.4088176) q[1];
sx q[1];
rz(-1.4585739) q[1];
sx q[1];
rz(-1.4558314) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30688863) q[0];
sx q[0];
rz(-1.406989) q[0];
sx q[0];
rz(0.52229806) q[0];
rz(-0.98241229) q[2];
sx q[2];
rz(-2.8937701) q[2];
sx q[2];
rz(-3.0318225) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.12633623) q[1];
sx q[1];
rz(-0.70955196) q[1];
sx q[1];
rz(-2.9615598) q[1];
rz(3.0343541) q[3];
sx q[3];
rz(-1.3235561) q[3];
sx q[3];
rz(-3.1013754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.17560656) q[2];
sx q[2];
rz(-2.6484428) q[2];
sx q[2];
rz(1.7612339) q[2];
rz(-0.68584758) q[3];
sx q[3];
rz(-0.0018456056) q[3];
sx q[3];
rz(2.4628911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-3.1252742) q[2];
sx q[2];
rz(-1.6961799) q[2];
sx q[2];
rz(0.21519306) q[2];
rz(1.9388922) q[3];
sx q[3];
rz(-0.31976578) q[3];
sx q[3];
rz(-3.0497677) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
