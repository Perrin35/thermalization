OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.74270785) q[0];
sx q[0];
rz(-0.78200114) q[0];
sx q[0];
rz(1.8703823) q[0];
rz(0.27702364) q[1];
sx q[1];
rz(-0.47200534) q[1];
sx q[1];
rz(-3.1402052) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9936375) q[0];
sx q[0];
rz(-1.4076828) q[0];
sx q[0];
rz(1.1975343) q[0];
rz(-pi) q[1];
rz(1.5288562) q[2];
sx q[2];
rz(-1.1239927) q[2];
sx q[2];
rz(-2.1697793) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6359771) q[1];
sx q[1];
rz(-0.50230366) q[1];
sx q[1];
rz(-2.99519) q[1];
rz(-pi) q[2];
rz(-2.8858521) q[3];
sx q[3];
rz(-1.2107953) q[3];
sx q[3];
rz(2.6933302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9871621) q[2];
sx q[2];
rz(-2.5240832) q[2];
sx q[2];
rz(2.3922065) q[2];
rz(1.0162214) q[3];
sx q[3];
rz(-1.1775492) q[3];
sx q[3];
rz(2.7367676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(1.4352903) q[0];
sx q[0];
rz(-2.3162233) q[0];
sx q[0];
rz(-2.170927) q[0];
rz(2.1043815) q[1];
sx q[1];
rz(-1.4379921) q[1];
sx q[1];
rz(2.326139) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52112752) q[0];
sx q[0];
rz(-1.5107811) q[0];
sx q[0];
rz(-2.492766) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8380727) q[2];
sx q[2];
rz(-0.75280658) q[2];
sx q[2];
rz(-0.34740651) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.38824575) q[1];
sx q[1];
rz(-1.8716295) q[1];
sx q[1];
rz(1.5335598) q[1];
rz(2.4957982) q[3];
sx q[3];
rz(-1.7495973) q[3];
sx q[3];
rz(-1.7088695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4619535) q[2];
sx q[2];
rz(-1.5660428) q[2];
sx q[2];
rz(2.5088076) q[2];
rz(1.1535545) q[3];
sx q[3];
rz(-0.76806918) q[3];
sx q[3];
rz(2.8320584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84045029) q[0];
sx q[0];
rz(-1.9294894) q[0];
sx q[0];
rz(-2.2667623) q[0];
rz(-1.3300928) q[1];
sx q[1];
rz(-1.4346088) q[1];
sx q[1];
rz(2.1420746) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32854983) q[0];
sx q[0];
rz(-1.8694436) q[0];
sx q[0];
rz(-0.15655984) q[0];
rz(-pi) q[1];
rz(2.3109762) q[2];
sx q[2];
rz(-2.3306371) q[2];
sx q[2];
rz(-2.9142771) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.98993436) q[1];
sx q[1];
rz(-2.5807568) q[1];
sx q[1];
rz(2.6480617) q[1];
rz(-pi) q[2];
rz(-3.0619377) q[3];
sx q[3];
rz(-2.0783391) q[3];
sx q[3];
rz(1.5910651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.53753608) q[2];
sx q[2];
rz(-2.2194922) q[2];
sx q[2];
rz(-0.58004722) q[2];
rz(0.81702685) q[3];
sx q[3];
rz(-1.3823119) q[3];
sx q[3];
rz(-1.9918611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76628768) q[0];
sx q[0];
rz(-1.5910609) q[0];
sx q[0];
rz(-0.91039175) q[0];
rz(-0.45122775) q[1];
sx q[1];
rz(-1.5952361) q[1];
sx q[1];
rz(-2.8667563) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1823605) q[0];
sx q[0];
rz(-1.649547) q[0];
sx q[0];
rz(-1.4819281) q[0];
rz(-3.1399973) q[2];
sx q[2];
rz(-1.110382) q[2];
sx q[2];
rz(0.38412016) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.53510016) q[1];
sx q[1];
rz(-2.1537158) q[1];
sx q[1];
rz(-0.54775723) q[1];
rz(-pi) q[2];
rz(0.87423012) q[3];
sx q[3];
rz(-1.4733553) q[3];
sx q[3];
rz(-1.8271354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3466907) q[2];
sx q[2];
rz(-2.0236423) q[2];
sx q[2];
rz(1.6332731) q[2];
rz(-1.9968962) q[3];
sx q[3];
rz(-0.73994023) q[3];
sx q[3];
rz(2.9798853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65790025) q[0];
sx q[0];
rz(-1.2150486) q[0];
sx q[0];
rz(2.143798) q[0];
rz(0.18355852) q[1];
sx q[1];
rz(-1.4869556) q[1];
sx q[1];
rz(1.516974) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.330634) q[0];
sx q[0];
rz(-1.9348382) q[0];
sx q[0];
rz(-0.90322687) q[0];
rz(-pi) q[1];
rz(0.60913182) q[2];
sx q[2];
rz(-2.5172148) q[2];
sx q[2];
rz(2.6513211) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5114054) q[1];
sx q[1];
rz(-1.1668219) q[1];
sx q[1];
rz(-1.1955111) q[1];
rz(-pi) q[2];
rz(-0.59668221) q[3];
sx q[3];
rz(-2.0697069) q[3];
sx q[3];
rz(0.98017207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8020442) q[2];
sx q[2];
rz(-0.99207726) q[2];
sx q[2];
rz(2.8175763) q[2];
rz(-1.8185395) q[3];
sx q[3];
rz(-2.3855305) q[3];
sx q[3];
rz(1.6103305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3762387) q[0];
sx q[0];
rz(-2.1026251) q[0];
sx q[0];
rz(1.8776241) q[0];
rz(2.2309247) q[1];
sx q[1];
rz(-1.940454) q[1];
sx q[1];
rz(-2.8009159) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1893038) q[0];
sx q[0];
rz(-1.6366819) q[0];
sx q[0];
rz(0.031866372) q[0];
rz(0.75603007) q[2];
sx q[2];
rz(-1.7871734) q[2];
sx q[2];
rz(-1.1083958) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1813982) q[1];
sx q[1];
rz(-1.8006693) q[1];
sx q[1];
rz(2.6541436) q[1];
x q[2];
rz(-0.42231456) q[3];
sx q[3];
rz(-1.1186244) q[3];
sx q[3];
rz(-0.57002588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.95057758) q[2];
sx q[2];
rz(-2.5577929) q[2];
sx q[2];
rz(-0.77159709) q[2];
rz(-0.54780444) q[3];
sx q[3];
rz(-0.9698202) q[3];
sx q[3];
rz(0.56345338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9724378) q[0];
sx q[0];
rz(-1.6945524) q[0];
sx q[0];
rz(2.7959438) q[0];
rz(3.0787643) q[1];
sx q[1];
rz(-0.47880104) q[1];
sx q[1];
rz(-0.46494928) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60663215) q[0];
sx q[0];
rz(-0.50919845) q[0];
sx q[0];
rz(1.9726994) q[0];
rz(2.2642235) q[2];
sx q[2];
rz(-2.4420028) q[2];
sx q[2];
rz(0.96044651) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9629434) q[1];
sx q[1];
rz(-1.4649676) q[1];
sx q[1];
rz(-2.9220198) q[1];
rz(-pi) q[2];
rz(-1.0602337) q[3];
sx q[3];
rz(-0.052882346) q[3];
sx q[3];
rz(-2.3898861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0039625) q[2];
sx q[2];
rz(-1.5045065) q[2];
sx q[2];
rz(-2.8239992) q[2];
rz(2.5701304) q[3];
sx q[3];
rz(-2.1025889) q[3];
sx q[3];
rz(2.8542744) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6222318) q[0];
sx q[0];
rz(-1.8472291) q[0];
sx q[0];
rz(2.8572594) q[0];
rz(-2.590086) q[1];
sx q[1];
rz(-0.14177828) q[1];
sx q[1];
rz(-3.0632339) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0067622234) q[0];
sx q[0];
rz(-1.5015331) q[0];
sx q[0];
rz(2.3368821) q[0];
rz(-pi) q[1];
rz(-2.9044754) q[2];
sx q[2];
rz(-1.3318921) q[2];
sx q[2];
rz(2.4799926) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6137177) q[1];
sx q[1];
rz(-2.4637239) q[1];
sx q[1];
rz(-1.2916958) q[1];
rz(1.7275229) q[3];
sx q[3];
rz(-1.6046451) q[3];
sx q[3];
rz(2.1050997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4006965) q[2];
sx q[2];
rz(-2.233278) q[2];
sx q[2];
rz(0.25137869) q[2];
rz(-2.5583983) q[3];
sx q[3];
rz(-2.0299032) q[3];
sx q[3];
rz(-1.73197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-0.79779977) q[0];
sx q[0];
rz(-0.082158953) q[0];
sx q[0];
rz(3.0902241) q[0];
rz(-2.2180166) q[1];
sx q[1];
rz(-2.4802465) q[1];
sx q[1];
rz(0.87402469) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6459991) q[0];
sx q[0];
rz(-0.80276239) q[0];
sx q[0];
rz(-1.5587224) q[0];
rz(-2.2888695) q[2];
sx q[2];
rz(-2.6211779) q[2];
sx q[2];
rz(-0.66327099) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.28305001) q[1];
sx q[1];
rz(-2.0701323) q[1];
sx q[1];
rz(-3.0602171) q[1];
rz(-pi) q[2];
rz(-1.9130575) q[3];
sx q[3];
rz(-0.98981333) q[3];
sx q[3];
rz(2.5653605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.41436568) q[2];
sx q[2];
rz(-0.7545158) q[2];
sx q[2];
rz(-0.61974636) q[2];
rz(-1.184458) q[3];
sx q[3];
rz(-1.8871566) q[3];
sx q[3];
rz(-1.363389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0062362) q[0];
sx q[0];
rz(-2.0993435) q[0];
sx q[0];
rz(0.7243048) q[0];
rz(-0.18877098) q[1];
sx q[1];
rz(-2.962208) q[1];
sx q[1];
rz(1.1788517) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8783979) q[0];
sx q[0];
rz(-0.34391719) q[0];
sx q[0];
rz(-0.92349903) q[0];
x q[1];
rz(-0.63209052) q[2];
sx q[2];
rz(-2.9846016) q[2];
sx q[2];
rz(2.0096411) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5621592) q[1];
sx q[1];
rz(-2.4907618) q[1];
sx q[1];
rz(1.8821554) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.52272777) q[3];
sx q[3];
rz(-1.3739112) q[3];
sx q[3];
rz(0.47375351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0835691) q[2];
sx q[2];
rz(-1.0408137) q[2];
sx q[2];
rz(-2.2422092) q[2];
rz(-2.2670238) q[3];
sx q[3];
rz(-0.42566291) q[3];
sx q[3];
rz(-1.6806867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5205004) q[0];
sx q[0];
rz(-2.7324471) q[0];
sx q[0];
rz(0.24656217) q[0];
rz(0.75795603) q[1];
sx q[1];
rz(-1.6592204) q[1];
sx q[1];
rz(1.6827676) q[1];
rz(-1.1307217) q[2];
sx q[2];
rz(-1.6025087) q[2];
sx q[2];
rz(-0.58379731) q[2];
rz(0.48537985) q[3];
sx q[3];
rz(-0.30967064) q[3];
sx q[3];
rz(-2.526024) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
