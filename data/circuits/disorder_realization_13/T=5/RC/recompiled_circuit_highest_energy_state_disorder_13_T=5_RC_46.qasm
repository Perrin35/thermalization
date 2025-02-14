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
rz(-2.9450671) q[0];
sx q[0];
rz(-2.3152469) q[0];
sx q[0];
rz(2.2235121) q[0];
rz(-0.11153587) q[1];
sx q[1];
rz(-1.3476975) q[1];
sx q[1];
rz(-1.5741875) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8188326) q[0];
sx q[0];
rz(-0.010060223) q[0];
sx q[0];
rz(-1.9181817) q[0];
x q[1];
rz(1.1996881) q[2];
sx q[2];
rz(-2.1016309) q[2];
sx q[2];
rz(-2.4715854) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7899129) q[1];
sx q[1];
rz(-1.4654298) q[1];
sx q[1];
rz(-1.8743452) q[1];
rz(2.6484916) q[3];
sx q[3];
rz(-1.6776909) q[3];
sx q[3];
rz(-1.3809539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.427318) q[2];
sx q[2];
rz(-1.6552507) q[2];
sx q[2];
rz(1.2662668) q[2];
rz(2.0990939) q[3];
sx q[3];
rz(-1.5407341) q[3];
sx q[3];
rz(1.7460167) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0365486) q[0];
sx q[0];
rz(-2.5140913) q[0];
sx q[0];
rz(0.096916048) q[0];
rz(-2.3333683) q[1];
sx q[1];
rz(-2.7327635) q[1];
sx q[1];
rz(1.2867297) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5135749) q[0];
sx q[0];
rz(-2.6903408) q[0];
sx q[0];
rz(-1.4892764) q[0];
rz(-pi) q[1];
rz(-0.64435763) q[2];
sx q[2];
rz(-2.3080809) q[2];
sx q[2];
rz(-0.65284158) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.60139537) q[1];
sx q[1];
rz(-2.2608065) q[1];
sx q[1];
rz(-2.6371535) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3318491) q[3];
sx q[3];
rz(-2.6698723) q[3];
sx q[3];
rz(2.7641486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5010995) q[2];
sx q[2];
rz(-2.2576136) q[2];
sx q[2];
rz(2.7596149) q[2];
rz(0.52465087) q[3];
sx q[3];
rz(-1.4068406) q[3];
sx q[3];
rz(-0.14373556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3749738) q[0];
sx q[0];
rz(-1.5856278) q[0];
sx q[0];
rz(0.68614706) q[0];
rz(-2.0740017) q[1];
sx q[1];
rz(-0.99718863) q[1];
sx q[1];
rz(3.0608665) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9322745) q[0];
sx q[0];
rz(-0.22695146) q[0];
sx q[0];
rz(-1.7448241) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3600134) q[2];
sx q[2];
rz(-0.8197166) q[2];
sx q[2];
rz(2.9224861) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0856493) q[1];
sx q[1];
rz(-2.9102444) q[1];
sx q[1];
rz(3.0401708) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0423546) q[3];
sx q[3];
rz(-0.84065372) q[3];
sx q[3];
rz(1.366758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2637691) q[2];
sx q[2];
rz(-2.4407385) q[2];
sx q[2];
rz(2.7027255) q[2];
rz(1.8435439) q[3];
sx q[3];
rz(-0.38709199) q[3];
sx q[3];
rz(-0.41518655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91442672) q[0];
sx q[0];
rz(-2.5974847) q[0];
sx q[0];
rz(-2.4650204) q[0];
rz(0.1380955) q[1];
sx q[1];
rz(-2.0510249) q[1];
sx q[1];
rz(-2.9409883) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5745068) q[0];
sx q[0];
rz(-2.4100523) q[0];
sx q[0];
rz(1.0163496) q[0];
rz(0.62941636) q[2];
sx q[2];
rz(-1.417629) q[2];
sx q[2];
rz(-1.5063733) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.84465161) q[1];
sx q[1];
rz(-1.703023) q[1];
sx q[1];
rz(0.6118212) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4273321) q[3];
sx q[3];
rz(-1.0142361) q[3];
sx q[3];
rz(1.0866764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5033919) q[2];
sx q[2];
rz(-1.8663422) q[2];
sx q[2];
rz(1.4208043) q[2];
rz(-0.11451379) q[3];
sx q[3];
rz(-1.4736466) q[3];
sx q[3];
rz(0.99944559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81922174) q[0];
sx q[0];
rz(-2.350816) q[0];
sx q[0];
rz(-2.1694699) q[0];
rz(0.32360336) q[1];
sx q[1];
rz(-1.4515896) q[1];
sx q[1];
rz(1.9591029) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13045584) q[0];
sx q[0];
rz(-0.19936518) q[0];
sx q[0];
rz(1.7749271) q[0];
x q[1];
rz(-3.1321822) q[2];
sx q[2];
rz(-2.5713671) q[2];
sx q[2];
rz(-1.648099) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.8684291) q[1];
sx q[1];
rz(-1.4123962) q[1];
sx q[1];
rz(-0.077977254) q[1];
rz(2.5706815) q[3];
sx q[3];
rz(-1.1215253) q[3];
sx q[3];
rz(0.8415287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9779196) q[2];
sx q[2];
rz(-1.369643) q[2];
sx q[2];
rz(-0.53691205) q[2];
rz(0.72530693) q[3];
sx q[3];
rz(-0.0497497) q[3];
sx q[3];
rz(0.79881001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86730114) q[0];
sx q[0];
rz(-2.0766356) q[0];
sx q[0];
rz(1.6987479) q[0];
rz(2.9310215) q[1];
sx q[1];
rz(-1.4434394) q[1];
sx q[1];
rz(-2.494536) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0083425) q[0];
sx q[0];
rz(-1.9444591) q[0];
sx q[0];
rz(1.8461807) q[0];
x q[1];
rz(-1.4616632) q[2];
sx q[2];
rz(-1.0169378) q[2];
sx q[2];
rz(0.46968974) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3270204) q[1];
sx q[1];
rz(-1.9144192) q[1];
sx q[1];
rz(2.6966641) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8043205) q[3];
sx q[3];
rz(-1.8303855) q[3];
sx q[3];
rz(1.0675499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0183384) q[2];
sx q[2];
rz(-2.047057) q[2];
sx q[2];
rz(-1.9617762) q[2];
rz(-1.8566462) q[3];
sx q[3];
rz(-0.40950567) q[3];
sx q[3];
rz(0.063260945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0354075) q[0];
sx q[0];
rz(-1.4466865) q[0];
sx q[0];
rz(-2.2440198) q[0];
rz(1.3342185) q[1];
sx q[1];
rz(-1.370627) q[1];
sx q[1];
rz(2.1468377) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75324995) q[0];
sx q[0];
rz(-0.80034798) q[0];
sx q[0];
rz(2.0092756) q[0];
rz(-pi) q[1];
rz(2.686196) q[2];
sx q[2];
rz(-1.0173305) q[2];
sx q[2];
rz(1.139251) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.82349087) q[1];
sx q[1];
rz(-1.8629322) q[1];
sx q[1];
rz(1.1000164) q[1];
rz(-pi) q[2];
rz(-3.1174421) q[3];
sx q[3];
rz(-1.2901879) q[3];
sx q[3];
rz(0.58089549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.083954088) q[2];
sx q[2];
rz(-1.8123764) q[2];
sx q[2];
rz(0.29339054) q[2];
rz(0.053226274) q[3];
sx q[3];
rz(-2.2727727) q[3];
sx q[3];
rz(1.7085541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5695802) q[0];
sx q[0];
rz(-1.7520289) q[0];
sx q[0];
rz(-2.8386175) q[0];
rz(-1.6527893) q[1];
sx q[1];
rz(-1.3158512) q[1];
sx q[1];
rz(2.4724919) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6953872) q[0];
sx q[0];
rz(-1.2406557) q[0];
sx q[0];
rz(-0.47254011) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1630493) q[2];
sx q[2];
rz(-2.7251232) q[2];
sx q[2];
rz(1.5451252) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4035621) q[1];
sx q[1];
rz(-2.4564017) q[1];
sx q[1];
rz(0.86802796) q[1];
rz(-pi) q[2];
rz(2.6458098) q[3];
sx q[3];
rz(-2.6675329) q[3];
sx q[3];
rz(2.1366675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1423433) q[2];
sx q[2];
rz(-1.3331022) q[2];
sx q[2];
rz(-0.86714253) q[2];
rz(2.0182746) q[3];
sx q[3];
rz(-0.85223782) q[3];
sx q[3];
rz(-1.999202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-1.2808696) q[0];
sx q[0];
rz(-0.46361247) q[0];
sx q[0];
rz(3.0238357) q[0];
rz(-0.87751687) q[1];
sx q[1];
rz(-2.120647) q[1];
sx q[1];
rz(-2.1868475) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60459899) q[0];
sx q[0];
rz(-2.1736896) q[0];
sx q[0];
rz(-2.0379809) q[0];
rz(-pi) q[1];
rz(-2.0659201) q[2];
sx q[2];
rz(-1.9946163) q[2];
sx q[2];
rz(2.0002805) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.89867184) q[1];
sx q[1];
rz(-1.7864711) q[1];
sx q[1];
rz(0.87529239) q[1];
rz(-pi) q[2];
rz(-2.4525663) q[3];
sx q[3];
rz(-2.195916) q[3];
sx q[3];
rz(2.2415015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.8942326) q[2];
sx q[2];
rz(-0.33984137) q[2];
sx q[2];
rz(1.4840688) q[2];
rz(1.5225211) q[3];
sx q[3];
rz(-1.3831474) q[3];
sx q[3];
rz(1.9671666) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5574944) q[0];
sx q[0];
rz(-1.431594) q[0];
sx q[0];
rz(-0.75310055) q[0];
rz(1.9901216) q[1];
sx q[1];
rz(-2.617372) q[1];
sx q[1];
rz(1.5392736) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4351411) q[0];
sx q[0];
rz(-2.0192084) q[0];
sx q[0];
rz(-2.9382692) q[0];
rz(-pi) q[1];
rz(-0.87290092) q[2];
sx q[2];
rz(-1.97556) q[2];
sx q[2];
rz(-1.0412316) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0311097) q[1];
sx q[1];
rz(-1.3785161) q[1];
sx q[1];
rz(2.7483205) q[1];
rz(-pi) q[2];
rz(1.3063823) q[3];
sx q[3];
rz(-0.9328649) q[3];
sx q[3];
rz(2.8273945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3573542) q[2];
sx q[2];
rz(-0.18111649) q[2];
sx q[2];
rz(0.46585807) q[2];
rz(1.0957796) q[3];
sx q[3];
rz(-0.992479) q[3];
sx q[3];
rz(2.5994658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0453542) q[0];
sx q[0];
rz(-1.5753373) q[0];
sx q[0];
rz(-1.5684431) q[0];
rz(2.0737598) q[1];
sx q[1];
rz(-2.7012431) q[1];
sx q[1];
rz(2.3213097) q[1];
rz(-1.8041415) q[2];
sx q[2];
rz(-1.8276311) q[2];
sx q[2];
rz(1.7414321) q[2];
rz(2.0910083) q[3];
sx q[3];
rz(-1.0903842) q[3];
sx q[3];
rz(0.36538418) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
