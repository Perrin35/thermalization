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
rz(1.6506305) q[0];
sx q[0];
rz(-1.4181674) q[0];
sx q[0];
rz(11.081063) q[0];
rz(-2.0909042) q[1];
sx q[1];
rz(4.8839999) q[1];
sx q[1];
rz(13.304741) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90621072) q[0];
sx q[0];
rz(-0.28774187) q[0];
sx q[0];
rz(0.55858992) q[0];
rz(-pi) q[1];
rz(1.9638871) q[2];
sx q[2];
rz(-1.7121268) q[2];
sx q[2];
rz(1.3648975) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3787427) q[1];
sx q[1];
rz(-2.0113011) q[1];
sx q[1];
rz(1.4164914) q[1];
rz(-1.2655696) q[3];
sx q[3];
rz(-2.1681227) q[3];
sx q[3];
rz(1.0553774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.38306132) q[2];
sx q[2];
rz(-2.8585275) q[2];
sx q[2];
rz(2.7281249) q[2];
rz(-2.5773898) q[3];
sx q[3];
rz(-1.6744303) q[3];
sx q[3];
rz(-1.9570501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41459945) q[0];
sx q[0];
rz(-1.0597543) q[0];
sx q[0];
rz(-3.0376036) q[0];
rz(0.14101401) q[1];
sx q[1];
rz(-0.68599373) q[1];
sx q[1];
rz(0.41735059) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.418103) q[0];
sx q[0];
rz(-1.8716629) q[0];
sx q[0];
rz(1.0194433) q[0];
x q[1];
rz(0.012243791) q[2];
sx q[2];
rz(-1.2070388) q[2];
sx q[2];
rz(0.084244339) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2552143) q[1];
sx q[1];
rz(-0.30579771) q[1];
sx q[1];
rz(-1.7096108) q[1];
x q[2];
rz(0.82666918) q[3];
sx q[3];
rz(-2.429768) q[3];
sx q[3];
rz(0.28361187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.523681) q[2];
sx q[2];
rz(-1.0289501) q[2];
sx q[2];
rz(-2.9376612) q[2];
rz(-1.5150874) q[3];
sx q[3];
rz(-0.47658673) q[3];
sx q[3];
rz(-1.6319857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5390891) q[0];
sx q[0];
rz(-1.6747549) q[0];
sx q[0];
rz(-2.7432826) q[0];
rz(1.0771982) q[1];
sx q[1];
rz(-0.64882433) q[1];
sx q[1];
rz(0.55346742) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9170415) q[0];
sx q[0];
rz(-1.5603491) q[0];
sx q[0];
rz(-2.4071669) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8936549) q[2];
sx q[2];
rz(-0.72681475) q[2];
sx q[2];
rz(-3.0023129) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0981507) q[1];
sx q[1];
rz(-2.4484608) q[1];
sx q[1];
rz(-2.5767703) q[1];
rz(2.3970277) q[3];
sx q[3];
rz(-0.25601632) q[3];
sx q[3];
rz(2.721173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0674949) q[2];
sx q[2];
rz(-2.7390538) q[2];
sx q[2];
rz(1.925776) q[2];
rz(2.1408234) q[3];
sx q[3];
rz(-1.3832904) q[3];
sx q[3];
rz(1.8892939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2860586) q[0];
sx q[0];
rz(-1.0581886) q[0];
sx q[0];
rz(-1.704294) q[0];
rz(-1.8057436) q[1];
sx q[1];
rz(-2.5156486) q[1];
sx q[1];
rz(-2.6864973) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2360596) q[0];
sx q[0];
rz(-0.62232557) q[0];
sx q[0];
rz(-3.0041191) q[0];
rz(-pi) q[1];
rz(-1.0735157) q[2];
sx q[2];
rz(-2.7233363) q[2];
sx q[2];
rz(-2.4474395) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0245106) q[1];
sx q[1];
rz(-0.41943892) q[1];
sx q[1];
rz(0.6371577) q[1];
rz(0.42174815) q[3];
sx q[3];
rz(-0.22796002) q[3];
sx q[3];
rz(0.82928073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8670292) q[2];
sx q[2];
rz(-0.44963351) q[2];
sx q[2];
rz(-0.8320128) q[2];
rz(0.653382) q[3];
sx q[3];
rz(-2.2218406) q[3];
sx q[3];
rz(-2.466989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6298544) q[0];
sx q[0];
rz(-1.3146223) q[0];
sx q[0];
rz(-0.6889371) q[0];
rz(2.9298933) q[1];
sx q[1];
rz(-1.7023106) q[1];
sx q[1];
rz(2.6874218) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9849761) q[0];
sx q[0];
rz(-2.1097221) q[0];
sx q[0];
rz(-0.52103736) q[0];
rz(-0.37068292) q[2];
sx q[2];
rz(-2.1739568) q[2];
sx q[2];
rz(-0.93450817) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8555357) q[1];
sx q[1];
rz(-2.2814676) q[1];
sx q[1];
rz(-2.2362806) q[1];
rz(-pi) q[2];
x q[2];
rz(0.16726475) q[3];
sx q[3];
rz(-1.8297046) q[3];
sx q[3];
rz(2.6791999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6058558) q[2];
sx q[2];
rz(-1.7822632) q[2];
sx q[2];
rz(0.5109171) q[2];
rz(-1.9262675) q[3];
sx q[3];
rz(-1.5769703) q[3];
sx q[3];
rz(0.63466614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0257492) q[0];
sx q[0];
rz(-2.8307493) q[0];
sx q[0];
rz(2.6281443) q[0];
rz(2.2250941) q[1];
sx q[1];
rz(-1.1767574) q[1];
sx q[1];
rz(1.7880012) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5149263) q[0];
sx q[0];
rz(-1.1276704) q[0];
sx q[0];
rz(-1.4791545) q[0];
rz(-0.84557477) q[2];
sx q[2];
rz(-0.52547272) q[2];
sx q[2];
rz(-2.2487244) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.533355) q[1];
sx q[1];
rz(-1.5441193) q[1];
sx q[1];
rz(0.4877301) q[1];
rz(-pi) q[2];
rz(0.065033536) q[3];
sx q[3];
rz(-1.4819996) q[3];
sx q[3];
rz(-0.2296344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6650271) q[2];
sx q[2];
rz(-2.5883784) q[2];
sx q[2];
rz(3.108007) q[2];
rz(0.50470662) q[3];
sx q[3];
rz(-1.4387771) q[3];
sx q[3];
rz(0.75468841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1228834) q[0];
sx q[0];
rz(-1.9174734) q[0];
sx q[0];
rz(1.1705742) q[0];
rz(-2.9508044) q[1];
sx q[1];
rz(-0.46563322) q[1];
sx q[1];
rz(0.64291397) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87827728) q[0];
sx q[0];
rz(-1.8366792) q[0];
sx q[0];
rz(-3.0356867) q[0];
rz(-pi) q[1];
x q[1];
rz(0.91251164) q[2];
sx q[2];
rz(-2.7325919) q[2];
sx q[2];
rz(0.3106948) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2819124) q[1];
sx q[1];
rz(-1.873906) q[1];
sx q[1];
rz(1.8213085) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6899818) q[3];
sx q[3];
rz(-1.9728807) q[3];
sx q[3];
rz(-1.6329671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1428895) q[2];
sx q[2];
rz(-0.10991749) q[2];
sx q[2];
rz(1.9996803) q[2];
rz(0.029684639) q[3];
sx q[3];
rz(-1.6073062) q[3];
sx q[3];
rz(2.3619385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90758816) q[0];
sx q[0];
rz(-1.1827844) q[0];
sx q[0];
rz(-1.8335861) q[0];
rz(2.0246778) q[1];
sx q[1];
rz(-1.6447379) q[1];
sx q[1];
rz(-2.9294779) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4329047) q[0];
sx q[0];
rz(-1.8861265) q[0];
sx q[0];
rz(-0.84724119) q[0];
x q[1];
rz(0.30435698) q[2];
sx q[2];
rz(-1.7946912) q[2];
sx q[2];
rz(-0.18098132) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.37112889) q[1];
sx q[1];
rz(-2.0088638) q[1];
sx q[1];
rz(-0.99645946) q[1];
rz(-pi) q[2];
rz(2.4808933) q[3];
sx q[3];
rz(-0.96008077) q[3];
sx q[3];
rz(-2.2211521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9461296) q[2];
sx q[2];
rz(-1.7886432) q[2];
sx q[2];
rz(-1.2809666) q[2];
rz(-1.403275) q[3];
sx q[3];
rz(-0.91126982) q[3];
sx q[3];
rz(0.72428552) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69287777) q[0];
sx q[0];
rz(-2.4249478) q[0];
sx q[0];
rz(0.88248673) q[0];
rz(2.8485883) q[1];
sx q[1];
rz(-2.2182783) q[1];
sx q[1];
rz(-2.015347) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4868797) q[0];
sx q[0];
rz(-0.77201399) q[0];
sx q[0];
rz(-2.604542) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.0056929767) q[2];
sx q[2];
rz(-2.5416592) q[2];
sx q[2];
rz(-1.1545187) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3224015) q[1];
sx q[1];
rz(-1.5171851) q[1];
sx q[1];
rz(-2.3399669) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2542488) q[3];
sx q[3];
rz(-1.1831565) q[3];
sx q[3];
rz(-0.46771184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5249411) q[2];
sx q[2];
rz(-0.61351073) q[2];
sx q[2];
rz(2.9743527) q[2];
rz(-0.031115726) q[3];
sx q[3];
rz(-1.9626706) q[3];
sx q[3];
rz(-0.64605609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6328218) q[0];
sx q[0];
rz(-1.333586) q[0];
sx q[0];
rz(-0.81047812) q[0];
rz(-1.3088016) q[1];
sx q[1];
rz(-0.93593132) q[1];
sx q[1];
rz(-3.0442309) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3980078) q[0];
sx q[0];
rz(-2.3438489) q[0];
sx q[0];
rz(-1.8037075) q[0];
rz(-pi) q[1];
rz(-1.6426769) q[2];
sx q[2];
rz(-0.38205636) q[2];
sx q[2];
rz(0.17554131) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3607935) q[1];
sx q[1];
rz(-2.153984) q[1];
sx q[1];
rz(-1.6287644) q[1];
x q[2];
rz(-0.11730365) q[3];
sx q[3];
rz(-1.3051093) q[3];
sx q[3];
rz(0.22646389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.29791609) q[2];
sx q[2];
rz(-2.127779) q[2];
sx q[2];
rz(1.8664912) q[2];
rz(-2.7013333) q[3];
sx q[3];
rz(-0.85846725) q[3];
sx q[3];
rz(-0.27537235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.020582747) q[0];
sx q[0];
rz(-2.5204211) q[0];
sx q[0];
rz(-0.68516635) q[0];
rz(-2.5195925) q[1];
sx q[1];
rz(-1.8432462) q[1];
sx q[1];
rz(-1.8274399) q[1];
rz(-1.6820108) q[2];
sx q[2];
rz(-1.6485452) q[2];
sx q[2];
rz(-2.9948276) q[2];
rz(-0.22777423) q[3];
sx q[3];
rz(-2.4854599) q[3];
sx q[3];
rz(1.5067185) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
