OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.83710837) q[0];
sx q[0];
rz(-1.4533071) q[0];
sx q[0];
rz(0.31153554) q[0];
rz(2.7040634) q[1];
sx q[1];
rz(4.9649927) q[1];
sx q[1];
rz(8.8658219) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1564724) q[0];
sx q[0];
rz(-2.7010242) q[0];
sx q[0];
rz(1.23929) q[0];
x q[1];
rz(-1.5582325) q[2];
sx q[2];
rz(-2.1280834) q[2];
sx q[2];
rz(1.3249427) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.47145876) q[1];
sx q[1];
rz(-2.11073) q[1];
sx q[1];
rz(-0.88175168) q[1];
rz(-1.2851089) q[3];
sx q[3];
rz(-1.9938855) q[3];
sx q[3];
rz(0.74564122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7493593) q[2];
sx q[2];
rz(-1.2831251) q[2];
sx q[2];
rz(0.63670811) q[2];
rz(-0.84896815) q[3];
sx q[3];
rz(-2.520112) q[3];
sx q[3];
rz(0.23392114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7648776) q[0];
sx q[0];
rz(-2.8945518) q[0];
sx q[0];
rz(-2.9887181) q[0];
rz(-0.75694594) q[1];
sx q[1];
rz(-1.5544954) q[1];
sx q[1];
rz(-0.98639948) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5312885) q[0];
sx q[0];
rz(-3.0695519) q[0];
sx q[0];
rz(-0.59285592) q[0];
rz(-pi) q[1];
x q[1];
rz(2.37843) q[2];
sx q[2];
rz(-1.1726511) q[2];
sx q[2];
rz(2.5415908) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.03147787) q[1];
sx q[1];
rz(-1.8077393) q[1];
sx q[1];
rz(0.38775878) q[1];
rz(-1.8996703) q[3];
sx q[3];
rz(-1.4423014) q[3];
sx q[3];
rz(-1.0458898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5793005) q[2];
sx q[2];
rz(-1.9220756) q[2];
sx q[2];
rz(-0.78312773) q[2];
rz(-3.1230208) q[3];
sx q[3];
rz(-1.5037856) q[3];
sx q[3];
rz(-2.7338681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0884393) q[0];
sx q[0];
rz(-0.29215559) q[0];
sx q[0];
rz(-0.96167481) q[0];
rz(2.7812474) q[1];
sx q[1];
rz(-2.0397489) q[1];
sx q[1];
rz(3.0128984) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0592244) q[0];
sx q[0];
rz(-1.6186065) q[0];
sx q[0];
rz(1.4882908) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3340036) q[2];
sx q[2];
rz(-1.1620887) q[2];
sx q[2];
rz(-1.7847716) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.62540185) q[1];
sx q[1];
rz(-1.1266202) q[1];
sx q[1];
rz(1.2605915) q[1];
rz(2.3388561) q[3];
sx q[3];
rz(-1.667913) q[3];
sx q[3];
rz(-0.10955284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0814357) q[2];
sx q[2];
rz(-1.1971985) q[2];
sx q[2];
rz(1.241768) q[2];
rz(-0.5870108) q[3];
sx q[3];
rz(-0.95654064) q[3];
sx q[3];
rz(2.164042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.509165) q[0];
sx q[0];
rz(-2.2594663) q[0];
sx q[0];
rz(-1.084491) q[0];
rz(1.658461) q[1];
sx q[1];
rz(-0.56743923) q[1];
sx q[1];
rz(-3.049057) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94938577) q[0];
sx q[0];
rz(-1.3017676) q[0];
sx q[0];
rz(1.2775248) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1763457) q[2];
sx q[2];
rz(-2.5778228) q[2];
sx q[2];
rz(3.0639067) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1558518) q[1];
sx q[1];
rz(-0.34793138) q[1];
sx q[1];
rz(2.9702529) q[1];
x q[2];
rz(2.9805698) q[3];
sx q[3];
rz(-0.92217126) q[3];
sx q[3];
rz(2.8341858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3080421) q[2];
sx q[2];
rz(-1.4414859) q[2];
sx q[2];
rz(0.33205024) q[2];
rz(-2.0856693) q[3];
sx q[3];
rz(-2.8639586) q[3];
sx q[3];
rz(0.61029303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8191391) q[0];
sx q[0];
rz(-1.8503014) q[0];
sx q[0];
rz(-0.21155587) q[0];
rz(-1.3062723) q[1];
sx q[1];
rz(-1.2439367) q[1];
sx q[1];
rz(0.64770118) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0928597) q[0];
sx q[0];
rz(-1.4591685) q[0];
sx q[0];
rz(2.9812921) q[0];
rz(-pi) q[1];
rz(1.7251882) q[2];
sx q[2];
rz(-1.9152181) q[2];
sx q[2];
rz(-1.4167348) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9808637) q[1];
sx q[1];
rz(-2.594922) q[1];
sx q[1];
rz(-1.1854118) q[1];
x q[2];
rz(-3.1224498) q[3];
sx q[3];
rz(-1.585841) q[3];
sx q[3];
rz(-2.6263833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5806879) q[2];
sx q[2];
rz(-2.7320392) q[2];
sx q[2];
rz(0.69331759) q[2];
rz(-0.66926113) q[3];
sx q[3];
rz(-1.4368493) q[3];
sx q[3];
rz(-0.0049237331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3174021) q[0];
sx q[0];
rz(-0.22739246) q[0];
sx q[0];
rz(-1.2325226) q[0];
rz(2.0690074) q[1];
sx q[1];
rz(-1.0718081) q[1];
sx q[1];
rz(-0.17428621) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6264544) q[0];
sx q[0];
rz(-2.2720552) q[0];
sx q[0];
rz(-0.40909543) q[0];
rz(-pi) q[1];
rz(2.3107489) q[2];
sx q[2];
rz(-0.85692642) q[2];
sx q[2];
rz(-1.4884399) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.78813997) q[1];
sx q[1];
rz(-2.3056245) q[1];
sx q[1];
rz(1.047903) q[1];
rz(-pi) q[2];
rz(-1.2234736) q[3];
sx q[3];
rz(-0.33063774) q[3];
sx q[3];
rz(0.91352458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3198513) q[2];
sx q[2];
rz(-0.80703002) q[2];
sx q[2];
rz(2.9439587) q[2];
rz(0.28891426) q[3];
sx q[3];
rz(-0.87696004) q[3];
sx q[3];
rz(1.7355841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(-1.063868) q[0];
sx q[0];
rz(-2.6517695) q[0];
sx q[0];
rz(0.20859627) q[0];
rz(-2.1754307) q[1];
sx q[1];
rz(-1.9816793) q[1];
sx q[1];
rz(1.6360412) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9150328) q[0];
sx q[0];
rz(-1.0462927) q[0];
sx q[0];
rz(2.4302308) q[0];
x q[1];
rz(-0.43039544) q[2];
sx q[2];
rz(-1.6039404) q[2];
sx q[2];
rz(0.37552777) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.437285) q[1];
sx q[1];
rz(-1.1204801) q[1];
sx q[1];
rz(-2.1561161) q[1];
x q[2];
rz(-2.8553477) q[3];
sx q[3];
rz(-1.4713333) q[3];
sx q[3];
rz(3.0486097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.85764) q[2];
sx q[2];
rz(-2.3985034) q[2];
sx q[2];
rz(-0.097578438) q[2];
rz(-1.7476667) q[3];
sx q[3];
rz(-1.3503617) q[3];
sx q[3];
rz(0.71684366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39032787) q[0];
sx q[0];
rz(-1.3289691) q[0];
sx q[0];
rz(2.6275997) q[0];
rz(-3.0184074) q[1];
sx q[1];
rz(-2.8942278) q[1];
sx q[1];
rz(2.2095912) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31873736) q[0];
sx q[0];
rz(-2.0120976) q[0];
sx q[0];
rz(0.1499099) q[0];
x q[1];
rz(-2.8727505) q[2];
sx q[2];
rz(-0.95250597) q[2];
sx q[2];
rz(2.6593069) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8646647) q[1];
sx q[1];
rz(-1.0229467) q[1];
sx q[1];
rz(-0.95363708) q[1];
rz(1.3436414) q[3];
sx q[3];
rz(-0.72941581) q[3];
sx q[3];
rz(1.9279355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0294068) q[2];
sx q[2];
rz(-1.1108578) q[2];
sx q[2];
rz(2.712148) q[2];
rz(1.9321692) q[3];
sx q[3];
rz(-0.37241396) q[3];
sx q[3];
rz(0.69303524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76686239) q[0];
sx q[0];
rz(-1.466789) q[0];
sx q[0];
rz(1.3388348) q[0];
rz(0.70612899) q[1];
sx q[1];
rz(-1.8920205) q[1];
sx q[1];
rz(0.95058092) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50981748) q[0];
sx q[0];
rz(-2.0777367) q[0];
sx q[0];
rz(-0.99960534) q[0];
x q[1];
rz(3.0663475) q[2];
sx q[2];
rz(-1.9855472) q[2];
sx q[2];
rz(-3.0441949) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3705759) q[1];
sx q[1];
rz(-0.5702714) q[1];
sx q[1];
rz(1.2194521) q[1];
rz(-pi) q[2];
x q[2];
rz(0.25456984) q[3];
sx q[3];
rz(-1.9707142) q[3];
sx q[3];
rz(0.039777048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6616228) q[2];
sx q[2];
rz(-1.6500436) q[2];
sx q[2];
rz(2.6573112) q[2];
rz(2.2144923) q[3];
sx q[3];
rz(-1.8323332) q[3];
sx q[3];
rz(-2.5203729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.1442239) q[0];
sx q[0];
rz(-3.0529418) q[0];
sx q[0];
rz(2.9123059) q[0];
rz(0.43481049) q[1];
sx q[1];
rz(-1.2229342) q[1];
sx q[1];
rz(2.4226709) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15764788) q[0];
sx q[0];
rz(-0.63417182) q[0];
sx q[0];
rz(1.4112524) q[0];
rz(-pi) q[1];
x q[1];
rz(0.98307307) q[2];
sx q[2];
rz(-1.6694153) q[2];
sx q[2];
rz(-1.3999511) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.016991888) q[1];
sx q[1];
rz(-1.4233839) q[1];
sx q[1];
rz(2.7684545) q[1];
rz(-pi) q[2];
x q[2];
rz(0.9142011) q[3];
sx q[3];
rz(-0.73540348) q[3];
sx q[3];
rz(0.38287336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3832613) q[2];
sx q[2];
rz(-1.2021844) q[2];
sx q[2];
rz(-0.74679217) q[2];
rz(-0.87219277) q[3];
sx q[3];
rz(-2.3132497) q[3];
sx q[3];
rz(-2.1993568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(1.223021) q[0];
sx q[0];
rz(-1.7445607) q[0];
sx q[0];
rz(1.8013409) q[0];
rz(-2.7643798) q[1];
sx q[1];
rz(-1.6709534) q[1];
sx q[1];
rz(0.51660641) q[1];
rz(1.8833075) q[2];
sx q[2];
rz(-1.0733114) q[2];
sx q[2];
rz(0.14807362) q[2];
rz(-0.17046455) q[3];
sx q[3];
rz(-2.4170605) q[3];
sx q[3];
rz(-1.3749882) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
