OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7528485) q[0];
sx q[0];
rz(-0.53628439) q[0];
sx q[0];
rz(-0.94776881) q[0];
rz(-1.3287969) q[1];
sx q[1];
rz(-1.8741908) q[1];
sx q[1];
rz(1.0277494) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0004814) q[0];
sx q[0];
rz(-1.1950462) q[0];
sx q[0];
rz(1.5462589) q[0];
rz(-pi) q[1];
rz(2.920354) q[2];
sx q[2];
rz(-2.049963) q[2];
sx q[2];
rz(2.0146807) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4029018) q[1];
sx q[1];
rz(-1.1876904) q[1];
sx q[1];
rz(-0.72523592) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5563452) q[3];
sx q[3];
rz(-0.71551502) q[3];
sx q[3];
rz(-1.1701208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1296922) q[2];
sx q[2];
rz(-1.7069858) q[2];
sx q[2];
rz(-2.091308) q[2];
rz(-2.0283279) q[3];
sx q[3];
rz(-2.2498825) q[3];
sx q[3];
rz(-0.068107001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0691836) q[0];
sx q[0];
rz(-1.8962815) q[0];
sx q[0];
rz(-2.8438399) q[0];
rz(2.521926) q[1];
sx q[1];
rz(-2.1344118) q[1];
sx q[1];
rz(-2.0334977) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28716921) q[0];
sx q[0];
rz(-2.4872635) q[0];
sx q[0];
rz(-1.2229344) q[0];
rz(-2.7116398) q[2];
sx q[2];
rz(-2.1351095) q[2];
sx q[2];
rz(-1.0831837) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.260028) q[1];
sx q[1];
rz(-1.2177055) q[1];
sx q[1];
rz(0.41207037) q[1];
x q[2];
rz(1.8403345) q[3];
sx q[3];
rz(-2.2552935) q[3];
sx q[3];
rz(1.0052296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6790598) q[2];
sx q[2];
rz(-0.97683895) q[2];
sx q[2];
rz(0.97529808) q[2];
rz(-0.9179999) q[3];
sx q[3];
rz(-1.2851597) q[3];
sx q[3];
rz(-0.29618922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.179203) q[0];
sx q[0];
rz(-2.2900892) q[0];
sx q[0];
rz(-2.5986824) q[0];
rz(0.88223282) q[1];
sx q[1];
rz(-1.135332) q[1];
sx q[1];
rz(0.96484819) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.097792625) q[0];
sx q[0];
rz(-1.291073) q[0];
sx q[0];
rz(-2.8066737) q[0];
rz(-pi) q[1];
rz(2.8446571) q[2];
sx q[2];
rz(-1.529971) q[2];
sx q[2];
rz(2.9253935) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.94273401) q[1];
sx q[1];
rz(-3.004289) q[1];
sx q[1];
rz(1.0820461) q[1];
x q[2];
rz(-0.19823719) q[3];
sx q[3];
rz(-1.1612411) q[3];
sx q[3];
rz(-2.4026681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0470011) q[2];
sx q[2];
rz(-0.61085218) q[2];
sx q[2];
rz(2.0084521) q[2];
rz(0.23162332) q[3];
sx q[3];
rz(-1.8685721) q[3];
sx q[3];
rz(0.75727063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8621181) q[0];
sx q[0];
rz(-3.1311488) q[0];
sx q[0];
rz(1.7650771) q[0];
rz(0.51849413) q[1];
sx q[1];
rz(-1.2644178) q[1];
sx q[1];
rz(-0.24212295) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74894416) q[0];
sx q[0];
rz(-2.2376275) q[0];
sx q[0];
rz(-1.3632266) q[0];
rz(-3.1130586) q[2];
sx q[2];
rz(-0.42978537) q[2];
sx q[2];
rz(-0.98791771) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8230069) q[1];
sx q[1];
rz(-1.4011523) q[1];
sx q[1];
rz(0.45348788) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7235892) q[3];
sx q[3];
rz(-1.4417218) q[3];
sx q[3];
rz(-1.2152745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.1233998) q[2];
sx q[2];
rz(-2.172956) q[2];
sx q[2];
rz(0.094853178) q[2];
rz(1.799396) q[3];
sx q[3];
rz(-1.3972524) q[3];
sx q[3];
rz(0.43911394) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3577394) q[0];
sx q[0];
rz(-1.8308324) q[0];
sx q[0];
rz(-2.7807996) q[0];
rz(1.7533253) q[1];
sx q[1];
rz(-1.3307945) q[1];
sx q[1];
rz(-2.0070019) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9903588) q[0];
sx q[0];
rz(-0.84999527) q[0];
sx q[0];
rz(-0.6806586) q[0];
rz(-pi) q[1];
rz(0.88419948) q[2];
sx q[2];
rz(-0.53495416) q[2];
sx q[2];
rz(1.7810437) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7736588) q[1];
sx q[1];
rz(-1.415167) q[1];
sx q[1];
rz(-2.6917798) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7048022) q[3];
sx q[3];
rz(-1.646864) q[3];
sx q[3];
rz(1.9074744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1363042) q[2];
sx q[2];
rz(-1.7207928) q[2];
sx q[2];
rz(-2.5689382) q[2];
rz(0.92875656) q[3];
sx q[3];
rz(-0.52162617) q[3];
sx q[3];
rz(-1.9740392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3271493) q[0];
sx q[0];
rz(-1.0961908) q[0];
sx q[0];
rz(1.2493398) q[0];
rz(1.2231187) q[1];
sx q[1];
rz(-1.5253116) q[1];
sx q[1];
rz(-1.9893601) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6845247) q[0];
sx q[0];
rz(-3.0355434) q[0];
sx q[0];
rz(-1.4507136) q[0];
rz(-pi) q[1];
rz(0.67201519) q[2];
sx q[2];
rz(-0.61215559) q[2];
sx q[2];
rz(1.200058) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7631543) q[1];
sx q[1];
rz(-1.7708781) q[1];
sx q[1];
rz(1.4640019) q[1];
rz(-0.67752083) q[3];
sx q[3];
rz(-2.065425) q[3];
sx q[3];
rz(-0.90466162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.46889177) q[2];
sx q[2];
rz(-1.2142618) q[2];
sx q[2];
rz(2.1506298) q[2];
rz(-2.4937566) q[3];
sx q[3];
rz(-2.1760553) q[3];
sx q[3];
rz(-0.34753862) q[3];
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
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25794849) q[0];
sx q[0];
rz(-2.914496) q[0];
sx q[0];
rz(-3.0793072) q[0];
rz(-2.9557872) q[1];
sx q[1];
rz(-1.4567016) q[1];
sx q[1];
rz(-2.7468162) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1226574) q[0];
sx q[0];
rz(-1.9648974) q[0];
sx q[0];
rz(-1.8271441) q[0];
rz(2.765352) q[2];
sx q[2];
rz(-1.8707841) q[2];
sx q[2];
rz(2.8397727) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2784087) q[1];
sx q[1];
rz(-1.3235958) q[1];
sx q[1];
rz(0.43988887) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.15345927) q[3];
sx q[3];
rz(-1.5342661) q[3];
sx q[3];
rz(-3.1094482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.64849598) q[2];
sx q[2];
rz(-2.811537) q[2];
sx q[2];
rz(-0.27302343) q[2];
rz(-1.8388883) q[3];
sx q[3];
rz(-1.3132934) q[3];
sx q[3];
rz(0.31204143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7438695) q[0];
sx q[0];
rz(-2.0070772) q[0];
sx q[0];
rz(1.8564818) q[0];
rz(-1.5015191) q[1];
sx q[1];
rz(-1.7506426) q[1];
sx q[1];
rz(1.3407019) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9467981) q[0];
sx q[0];
rz(-1.9614944) q[0];
sx q[0];
rz(0.35238738) q[0];
rz(-pi) q[1];
rz(-2.4457744) q[2];
sx q[2];
rz(-1.5317481) q[2];
sx q[2];
rz(-0.85862904) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2015842) q[1];
sx q[1];
rz(-1.9218947) q[1];
sx q[1];
rz(-1.9560948) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1287442) q[3];
sx q[3];
rz(-1.2841409) q[3];
sx q[3];
rz(0.096404508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0354707) q[2];
sx q[2];
rz(-0.80646986) q[2];
sx q[2];
rz(1.0236053) q[2];
rz(0.18493955) q[3];
sx q[3];
rz(-0.39026323) q[3];
sx q[3];
rz(2.8997054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0969365) q[0];
sx q[0];
rz(-0.99656314) q[0];
sx q[0];
rz(1.6212844) q[0];
rz(-0.3301436) q[1];
sx q[1];
rz(-1.2076999) q[1];
sx q[1];
rz(-2.3044589) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5569816) q[0];
sx q[0];
rz(-0.71174445) q[0];
sx q[0];
rz(0.31194375) q[0];
rz(-1.8673973) q[2];
sx q[2];
rz(-2.1269848) q[2];
sx q[2];
rz(-2.0402758) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.21048927) q[1];
sx q[1];
rz(-1.466202) q[1];
sx q[1];
rz(-1.4180257) q[1];
rz(-pi) q[2];
rz(-2.0666396) q[3];
sx q[3];
rz(-1.3131724) q[3];
sx q[3];
rz(0.50592929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.60124406) q[2];
sx q[2];
rz(-1.0849755) q[2];
sx q[2];
rz(-2.9619651) q[2];
rz(2.1458697) q[3];
sx q[3];
rz(-1.8928173) q[3];
sx q[3];
rz(-1.3109591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76374617) q[0];
sx q[0];
rz(-0.34559956) q[0];
sx q[0];
rz(-2.0843704) q[0];
rz(-0.10748848) q[1];
sx q[1];
rz(-1.2534393) q[1];
sx q[1];
rz(-0.9799788) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7769023) q[0];
sx q[0];
rz(-1.5199465) q[0];
sx q[0];
rz(-1.4450253) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6238742) q[2];
sx q[2];
rz(-1.1283518) q[2];
sx q[2];
rz(0.25416086) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.92614782) q[1];
sx q[1];
rz(-1.1415592) q[1];
sx q[1];
rz(2.8785359) q[1];
rz(1.8944593) q[3];
sx q[3];
rz(-1.666288) q[3];
sx q[3];
rz(0.0088866339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3048627) q[2];
sx q[2];
rz(-1.9367846) q[2];
sx q[2];
rz(-1.0277964) q[2];
rz(1.3868388) q[3];
sx q[3];
rz(-1.8324865) q[3];
sx q[3];
rz(-2.8579779) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2733611) q[0];
sx q[0];
rz(-2.1049451) q[0];
sx q[0];
rz(2.0275397) q[0];
rz(-2.3095619) q[1];
sx q[1];
rz(-0.46453005) q[1];
sx q[1];
rz(0.66418905) q[1];
rz(2.1297395) q[2];
sx q[2];
rz(-1.4197299) q[2];
sx q[2];
rz(1.1001669) q[2];
rz(-2.5084393) q[3];
sx q[3];
rz(-0.59637759) q[3];
sx q[3];
rz(-0.37806088) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
