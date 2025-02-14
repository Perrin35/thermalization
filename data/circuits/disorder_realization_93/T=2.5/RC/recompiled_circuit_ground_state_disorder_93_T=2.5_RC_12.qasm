OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.20004162) q[0];
sx q[0];
rz(3.3257421) q[0];
sx q[0];
rz(10.781547) q[0];
rz(-1.6954724) q[1];
sx q[1];
rz(-1.5741916) q[1];
sx q[1];
rz(2.8653436) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8987792) q[0];
sx q[0];
rz(-1.5540637) q[0];
sx q[0];
rz(-1.3748598) q[0];
rz(1.5612565) q[2];
sx q[2];
rz(-1.6838639) q[2];
sx q[2];
rz(0.58361562) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3018926) q[1];
sx q[1];
rz(-3.0954315) q[1];
sx q[1];
rz(-2.7916174) q[1];
x q[2];
rz(1.5050234) q[3];
sx q[3];
rz(-1.2322174) q[3];
sx q[3];
rz(-2.8985646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.19286284) q[2];
sx q[2];
rz(-2.7728045) q[2];
sx q[2];
rz(-0.31664872) q[2];
rz(-0.51490319) q[3];
sx q[3];
rz(-0.0080990214) q[3];
sx q[3];
rz(2.3308339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1279476) q[0];
sx q[0];
rz(-2.7982496) q[0];
sx q[0];
rz(-3.1105594) q[0];
rz(-1.7015999) q[1];
sx q[1];
rz(-0.80039918) q[1];
sx q[1];
rz(1.4061141) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8588327) q[0];
sx q[0];
rz(-1.7525391) q[0];
sx q[0];
rz(1.9980998) q[0];
rz(0.038766301) q[2];
sx q[2];
rz(-1.7749466) q[2];
sx q[2];
rz(-0.37224712) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5819491) q[1];
sx q[1];
rz(-0.20155263) q[1];
sx q[1];
rz(-0.45362107) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.77649103) q[3];
sx q[3];
rz(-1.6776976) q[3];
sx q[3];
rz(2.2872146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.42355737) q[2];
sx q[2];
rz(-2.9325298) q[2];
sx q[2];
rz(-2.5162856) q[2];
rz(1.5710255) q[3];
sx q[3];
rz(-0.26206854) q[3];
sx q[3];
rz(3.1077793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5686947) q[0];
sx q[0];
rz(-0.53302065) q[0];
sx q[0];
rz(2.0030588) q[0];
rz(-0.89433995) q[1];
sx q[1];
rz(-3.1411451) q[1];
sx q[1];
rz(1.0630039) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3819696) q[0];
sx q[0];
rz(-2.2142525) q[0];
sx q[0];
rz(-0.30610419) q[0];
rz(-0.74805061) q[2];
sx q[2];
rz(-1.4292293) q[2];
sx q[2];
rz(0.54439616) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.56617006) q[1];
sx q[1];
rz(-2.992792) q[1];
sx q[1];
rz(-2.0729561) q[1];
x q[2];
rz(-2.9195905) q[3];
sx q[3];
rz(-1.7736122) q[3];
sx q[3];
rz(-2.0012282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.87463093) q[2];
sx q[2];
rz(-0.78445542) q[2];
sx q[2];
rz(-1.1355012) q[2];
rz(1.9973495) q[3];
sx q[3];
rz(-3.0502697) q[3];
sx q[3];
rz(1.7448366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8046232) q[0];
sx q[0];
rz(-2.7380044) q[0];
sx q[0];
rz(-1.1176156) q[0];
rz(0.090967372) q[1];
sx q[1];
rz(-3.1035943) q[1];
sx q[1];
rz(-1.5100347) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6081597) q[0];
sx q[0];
rz(-1.3638487) q[0];
sx q[0];
rz(2.1536835) q[0];
x q[1];
rz(-2.9363459) q[2];
sx q[2];
rz(-1.4637498) q[2];
sx q[2];
rz(-2.7595124) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7638844) q[1];
sx q[1];
rz(-0.33033366) q[1];
sx q[1];
rz(-1.1829472) q[1];
rz(-pi) q[2];
x q[2];
rz(0.61057874) q[3];
sx q[3];
rz(-2.0269666) q[3];
sx q[3];
rz(1.3464286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9522004) q[2];
sx q[2];
rz(-3.0244833) q[2];
sx q[2];
rz(-1.2516775) q[2];
rz(2.7513454) q[3];
sx q[3];
rz(-0.016642112) q[3];
sx q[3];
rz(0.19002953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2802884) q[0];
sx q[0];
rz(-2.5087293) q[0];
sx q[0];
rz(1.2907775) q[0];
rz(-0.0064918594) q[1];
sx q[1];
rz(-2.8752654) q[1];
sx q[1];
rz(-1.806462) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2119009) q[0];
sx q[0];
rz(-1.9121355) q[0];
sx q[0];
rz(-1.4262188) q[0];
x q[1];
rz(-0.13258055) q[2];
sx q[2];
rz(-2.6040316) q[2];
sx q[2];
rz(-2.3258296) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.74400988) q[1];
sx q[1];
rz(-2.8967794) q[1];
sx q[1];
rz(2.6086065) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9560821) q[3];
sx q[3];
rz(-1.9855567) q[3];
sx q[3];
rz(-1.1832373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.409965) q[2];
sx q[2];
rz(-0.16606398) q[2];
sx q[2];
rz(-0.49327332) q[2];
rz(-1.0924245) q[3];
sx q[3];
rz(-0.0041996669) q[3];
sx q[3];
rz(-2.5815729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6888371) q[0];
sx q[0];
rz(-1.1863363) q[0];
sx q[0];
rz(-1.3415927) q[0];
rz(1.5780559) q[1];
sx q[1];
rz(-0.23461239) q[1];
sx q[1];
rz(-0.87306195) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36045433) q[0];
sx q[0];
rz(-2.6037881) q[0];
sx q[0];
rz(0.27307684) q[0];
x q[1];
rz(0.43967562) q[2];
sx q[2];
rz(-1.5182677) q[2];
sx q[2];
rz(2.479468) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7028919) q[1];
sx q[1];
rz(-0.14335614) q[1];
sx q[1];
rz(-2.1208255) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8382243) q[3];
sx q[3];
rz(-1.7828327) q[3];
sx q[3];
rz(0.39006453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0688429) q[2];
sx q[2];
rz(-1.9708956) q[2];
sx q[2];
rz(1.076131) q[2];
rz(-2.7592646) q[3];
sx q[3];
rz(-3.101442) q[3];
sx q[3];
rz(2.373915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1165987) q[0];
sx q[0];
rz(-3.0285663) q[0];
sx q[0];
rz(-0.22079994) q[0];
rz(2.2258672) q[1];
sx q[1];
rz(-0.16685367) q[1];
sx q[1];
rz(-1.6545666) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90900786) q[0];
sx q[0];
rz(-0.56031432) q[0];
sx q[0];
rz(-0.73580297) q[0];
rz(-2.0528232) q[2];
sx q[2];
rz(-2.7332167) q[2];
sx q[2];
rz(-2.0118304) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2579563) q[1];
sx q[1];
rz(-2.2410434) q[1];
sx q[1];
rz(0.61780091) q[1];
x q[2];
rz(1.7127895) q[3];
sx q[3];
rz(-1.4513223) q[3];
sx q[3];
rz(-0.17585299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9080092) q[2];
sx q[2];
rz(-0.073069409) q[2];
sx q[2];
rz(-1.118534) q[2];
rz(1.3633049) q[3];
sx q[3];
rz(-3.0647291) q[3];
sx q[3];
rz(1.3479056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8175279) q[0];
sx q[0];
rz(-2.9370566) q[0];
sx q[0];
rz(-1.8472141) q[0];
rz(-1.4140465) q[1];
sx q[1];
rz(-3.0906446) q[1];
sx q[1];
rz(-0.22305138) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18666727) q[0];
sx q[0];
rz(-2.3340805) q[0];
sx q[0];
rz(-0.6874716) q[0];
rz(-pi) q[1];
rz(0.18977087) q[2];
sx q[2];
rz(-1.4354354) q[2];
sx q[2];
rz(-0.0046829872) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1522882) q[1];
sx q[1];
rz(-0.76123255) q[1];
sx q[1];
rz(-2.1851995) q[1];
rz(-0.87568385) q[3];
sx q[3];
rz(-1.6611665) q[3];
sx q[3];
rz(-0.68708693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.53905067) q[2];
sx q[2];
rz(-1.9571303) q[2];
sx q[2];
rz(-0.43262532) q[2];
rz(1.0299579) q[3];
sx q[3];
rz(-0.10051388) q[3];
sx q[3];
rz(0.71737927) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75281993) q[0];
sx q[0];
rz(-2.6322375) q[0];
sx q[0];
rz(2.9260522) q[0];
rz(-1.132025) q[1];
sx q[1];
rz(-0.80765453) q[1];
sx q[1];
rz(1.6708803) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4819538) q[0];
sx q[0];
rz(-0.99876548) q[0];
sx q[0];
rz(-0.37449775) q[0];
rz(0.37781654) q[2];
sx q[2];
rz(-1.1616128) q[2];
sx q[2];
rz(2.5175187) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.52016693) q[1];
sx q[1];
rz(-1.5158719) q[1];
sx q[1];
rz(-1.6009694) q[1];
rz(-pi) q[2];
rz(0.4827667) q[3];
sx q[3];
rz(-1.996187) q[3];
sx q[3];
rz(-1.1541884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1620764) q[2];
sx q[2];
rz(-3.1023878) q[2];
sx q[2];
rz(-0.51852477) q[2];
rz(-2.3833852) q[3];
sx q[3];
rz(-2.3617305) q[3];
sx q[3];
rz(-1.4832179) q[3];
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
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6766847) q[0];
sx q[0];
rz(-1.2228891) q[0];
sx q[0];
rz(1.2567047) q[0];
rz(2.8769809) q[1];
sx q[1];
rz(-3.1393317) q[1];
sx q[1];
rz(1.9432) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1408069) q[0];
sx q[0];
rz(-0.0065751271) q[0];
sx q[0];
rz(3.0964609) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2368868) q[2];
sx q[2];
rz(-1.5057202) q[2];
sx q[2];
rz(-1.8462311) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2923515) q[1];
sx q[1];
rz(-1.3363943) q[1];
sx q[1];
rz(1.5567907) q[1];
x q[2];
rz(1.7400473) q[3];
sx q[3];
rz(-1.608132) q[3];
sx q[3];
rz(2.7558655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4635224) q[2];
sx q[2];
rz(-3.1375725) q[2];
sx q[2];
rz(0.013997495) q[2];
rz(0.30063453) q[3];
sx q[3];
rz(-3.1379353) q[3];
sx q[3];
rz(1.5699056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32776253) q[0];
sx q[0];
rz(-1.7131085) q[0];
sx q[0];
rz(-1.499927) q[0];
rz(3.1060524) q[1];
sx q[1];
rz(-2.7062369) q[1];
sx q[1];
rz(-2.8027262) q[1];
rz(-1.8076996) q[2];
sx q[2];
rz(-2.2541004) q[2];
sx q[2];
rz(-2.83082) q[2];
rz(0.20094674) q[3];
sx q[3];
rz(-2.3280716) q[3];
sx q[3];
rz(0.18542326) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
