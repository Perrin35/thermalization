OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3044843) q[0];
sx q[0];
rz(-1.6882856) q[0];
sx q[0];
rz(2.8300571) q[0];
rz(2.7040634) q[1];
sx q[1];
rz(4.9649927) q[1];
sx q[1];
rz(8.8658219) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11249609) q[0];
sx q[0];
rz(-1.4315499) q[0];
sx q[0];
rz(-1.990156) q[0];
rz(-pi) q[1];
x q[1];
rz(0.55732255) q[2];
sx q[2];
rz(-1.5601336) q[2];
sx q[2];
rz(-0.23920857) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.47145876) q[1];
sx q[1];
rz(-2.11073) q[1];
sx q[1];
rz(-0.88175168) q[1];
x q[2];
rz(2.7028014) q[3];
sx q[3];
rz(-1.8306797) q[3];
sx q[3];
rz(-0.70513844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7493593) q[2];
sx q[2];
rz(-1.8584676) q[2];
sx q[2];
rz(-2.5048845) q[2];
rz(-0.84896815) q[3];
sx q[3];
rz(-0.62148062) q[3];
sx q[3];
rz(2.9076715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7648776) q[0];
sx q[0];
rz(-2.8945518) q[0];
sx q[0];
rz(-0.15287457) q[0];
rz(-0.75694594) q[1];
sx q[1];
rz(-1.5870973) q[1];
sx q[1];
rz(-2.1551932) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61030412) q[0];
sx q[0];
rz(-0.072040759) q[0];
sx q[0];
rz(0.59285592) q[0];
rz(-pi) q[1];
rz(1.0436922) q[2];
sx q[2];
rz(-0.88000789) q[2];
sx q[2];
rz(-1.3259128) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0806597) q[1];
sx q[1];
rz(-0.45127171) q[1];
sx q[1];
rz(-2.5732451) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.13568474) q[3];
sx q[3];
rz(-1.8968582) q[3];
sx q[3];
rz(0.56860926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5622921) q[2];
sx q[2];
rz(-1.219517) q[2];
sx q[2];
rz(-2.3584649) q[2];
rz(0.018571818) q[3];
sx q[3];
rz(-1.5037856) q[3];
sx q[3];
rz(-2.7338681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.0531533) q[0];
sx q[0];
rz(-2.8494371) q[0];
sx q[0];
rz(-0.96167481) q[0];
rz(0.36034521) q[1];
sx q[1];
rz(-1.1018437) q[1];
sx q[1];
rz(3.0128984) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0592244) q[0];
sx q[0];
rz(-1.5229862) q[0];
sx q[0];
rz(-1.4882908) q[0];
x q[1];
rz(-0.80758904) q[2];
sx q[2];
rz(-1.979504) q[2];
sx q[2];
rz(1.7847716) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.62540185) q[1];
sx q[1];
rz(-1.1266202) q[1];
sx q[1];
rz(-1.2605915) q[1];
rz(-1.4314753) q[3];
sx q[3];
rz(-2.3686667) q[3];
sx q[3];
rz(-1.7804002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.06015691) q[2];
sx q[2];
rz(-1.9443941) q[2];
sx q[2];
rz(1.241768) q[2];
rz(2.5545819) q[3];
sx q[3];
rz(-0.95654064) q[3];
sx q[3];
rz(2.164042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.509165) q[0];
sx q[0];
rz(-2.2594663) q[0];
sx q[0];
rz(-1.084491) q[0];
rz(1.4831316) q[1];
sx q[1];
rz(-2.5741534) q[1];
sx q[1];
rz(-3.049057) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4400892) q[0];
sx q[0];
rz(-1.8532231) q[0];
sx q[0];
rz(-0.28042067) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1763457) q[2];
sx q[2];
rz(-2.5778228) q[2];
sx q[2];
rz(-3.0639067) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.42379984) q[1];
sx q[1];
rz(-1.51263) q[1];
sx q[1];
rz(-2.7983626) q[1];
rz(-2.9805698) q[3];
sx q[3];
rz(-0.92217126) q[3];
sx q[3];
rz(0.30740689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.83355054) q[2];
sx q[2];
rz(-1.4414859) q[2];
sx q[2];
rz(2.8095424) q[2];
rz(2.0856693) q[3];
sx q[3];
rz(-0.27763405) q[3];
sx q[3];
rz(0.61029303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32245359) q[0];
sx q[0];
rz(-1.8503014) q[0];
sx q[0];
rz(-0.21155587) q[0];
rz(1.8353204) q[1];
sx q[1];
rz(-1.2439367) q[1];
sx q[1];
rz(-2.4938915) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.081213148) q[0];
sx q[0];
rz(-0.19506422) q[0];
sx q[0];
rz(2.5293406) q[0];
rz(0.34824246) q[2];
sx q[2];
rz(-1.4255382) q[2];
sx q[2];
rz(0.20656221) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.065121) q[1];
sx q[1];
rz(-1.3741125) q[1];
sx q[1];
rz(-2.0842488) q[1];
rz(3.1224498) q[3];
sx q[3];
rz(-1.585841) q[3];
sx q[3];
rz(2.6263833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5806879) q[2];
sx q[2];
rz(-0.40955341) q[2];
sx q[2];
rz(-0.69331759) q[2];
rz(0.66926113) q[3];
sx q[3];
rz(-1.7047434) q[3];
sx q[3];
rz(3.1366689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3174021) q[0];
sx q[0];
rz(-0.22739246) q[0];
sx q[0];
rz(1.2325226) q[0];
rz(1.0725853) q[1];
sx q[1];
rz(-1.0718081) q[1];
sx q[1];
rz(-2.9673064) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6264544) q[0];
sx q[0];
rz(-0.8695375) q[0];
sx q[0];
rz(-2.7324972) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.8308438) q[2];
sx q[2];
rz(-0.85692642) q[2];
sx q[2];
rz(-1.4884399) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.41387687) q[1];
sx q[1];
rz(-1.9503647) q[1];
sx q[1];
rz(-2.3351401) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2234736) q[3];
sx q[3];
rz(-2.8109549) q[3];
sx q[3];
rz(-0.91352458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3198513) q[2];
sx q[2];
rz(-2.3345626) q[2];
sx q[2];
rz(0.19763395) q[2];
rz(2.8526784) q[3];
sx q[3];
rz(-2.2646326) q[3];
sx q[3];
rz(1.7355841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0777247) q[0];
sx q[0];
rz(-0.48982319) q[0];
sx q[0];
rz(-2.9329964) q[0];
rz(2.1754307) q[1];
sx q[1];
rz(-1.9816793) q[1];
sx q[1];
rz(1.5055515) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.063232139) q[0];
sx q[0];
rz(-2.1713543) q[0];
sx q[0];
rz(2.2230704) q[0];
rz(-pi) q[1];
rz(2.7111972) q[2];
sx q[2];
rz(-1.5376523) q[2];
sx q[2];
rz(2.7660649) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.70430763) q[1];
sx q[1];
rz(-2.0211126) q[1];
sx q[1];
rz(0.98547658) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.3397293) q[3];
sx q[3];
rz(-0.30258426) q[3];
sx q[3];
rz(-1.8031977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2839526) q[2];
sx q[2];
rz(-2.3985034) q[2];
sx q[2];
rz(-0.097578438) q[2];
rz(1.7476667) q[3];
sx q[3];
rz(-1.3503617) q[3];
sx q[3];
rz(-0.71684366) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39032787) q[0];
sx q[0];
rz(-1.3289691) q[0];
sx q[0];
rz(-2.6275997) q[0];
rz(-0.12318525) q[1];
sx q[1];
rz(-0.24736483) q[1];
sx q[1];
rz(2.2095912) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4829464) q[0];
sx q[0];
rz(-2.6771149) q[0];
sx q[0];
rz(-1.2645725) q[0];
rz(2.2064662) q[2];
sx q[2];
rz(-1.3526275) q[2];
sx q[2];
rz(2.2114434) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8646647) q[1];
sx q[1];
rz(-1.0229467) q[1];
sx q[1];
rz(-0.95363708) q[1];
x q[2];
rz(-1.7979513) q[3];
sx q[3];
rz(-0.72941581) q[3];
sx q[3];
rz(-1.2136572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0294068) q[2];
sx q[2];
rz(-1.1108578) q[2];
sx q[2];
rz(-0.4294447) q[2];
rz(-1.2094234) q[3];
sx q[3];
rz(-0.37241396) q[3];
sx q[3];
rz(0.69303524) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3747303) q[0];
sx q[0];
rz(-1.6748036) q[0];
sx q[0];
rz(-1.3388348) q[0];
rz(-0.70612899) q[1];
sx q[1];
rz(-1.2495722) q[1];
sx q[1];
rz(0.95058092) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50981748) q[0];
sx q[0];
rz(-1.063856) q[0];
sx q[0];
rz(-0.99960534) q[0];
rz(-pi) q[1];
rz(1.739903) q[2];
sx q[2];
rz(-2.7204614) q[2];
sx q[2];
rz(3.054045) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3705759) q[1];
sx q[1];
rz(-2.5713213) q[1];
sx q[1];
rz(-1.9221406) q[1];
rz(-0.25456984) q[3];
sx q[3];
rz(-1.1708784) q[3];
sx q[3];
rz(-3.1018156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4799698) q[2];
sx q[2];
rz(-1.6500436) q[2];
sx q[2];
rz(-2.6573112) q[2];
rz(-0.92710036) q[3];
sx q[3];
rz(-1.3092594) q[3];
sx q[3];
rz(-0.62121975) q[3];
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
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9973688) q[0];
sx q[0];
rz(-3.0529418) q[0];
sx q[0];
rz(-0.22928672) q[0];
rz(-2.7067822) q[1];
sx q[1];
rz(-1.2229342) q[1];
sx q[1];
rz(-0.71892175) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15764788) q[0];
sx q[0];
rz(-2.5074208) q[0];
sx q[0];
rz(1.7303403) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.11833338) q[2];
sx q[2];
rz(-2.1552857) q[2];
sx q[2];
rz(3.0362533) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9467266) q[1];
sx q[1];
rz(-2.7416639) q[1];
sx q[1];
rz(2.7547794) q[1];
rz(-pi) q[2];
rz(2.6370254) q[3];
sx q[3];
rz(-2.1310398) q[3];
sx q[3];
rz(-1.954078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3832613) q[2];
sx q[2];
rz(-1.9394082) q[2];
sx q[2];
rz(0.74679217) q[2];
rz(-0.87219277) q[3];
sx q[3];
rz(-2.3132497) q[3];
sx q[3];
rz(0.94223589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9185716) q[0];
sx q[0];
rz(-1.7445607) q[0];
sx q[0];
rz(1.8013409) q[0];
rz(-2.7643798) q[1];
sx q[1];
rz(-1.6709534) q[1];
sx q[1];
rz(0.51660641) q[1];
rz(-1.8833075) q[2];
sx q[2];
rz(-2.0682813) q[2];
sx q[2];
rz(-2.993519) q[2];
rz(1.7198346) q[3];
sx q[3];
rz(-2.2825713) q[3];
sx q[3];
rz(1.5406516) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
