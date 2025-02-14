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
rz(-2.9574432) q[0];
sx q[0];
rz(1.3567691) q[0];
rz(1.4461203) q[1];
sx q[1];
rz(-1.5674011) q[1];
sx q[1];
rz(-2.8653436) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8169308) q[0];
sx q[0];
rz(-1.766705) q[0];
sx q[0];
rz(-0.017058998) q[0];
rz(-pi) q[1];
rz(-0.11307271) q[2];
sx q[2];
rz(-1.5802752) q[2];
sx q[2];
rz(0.98825708) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6522112) q[1];
sx q[1];
rz(-1.5274352) q[1];
sx q[1];
rz(-1.5549589) q[1];
rz(-pi) q[2];
rz(-1.5050234) q[3];
sx q[3];
rz(-1.2322174) q[3];
sx q[3];
rz(-0.24302806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9487298) q[2];
sx q[2];
rz(-2.7728045) q[2];
sx q[2];
rz(0.31664872) q[2];
rz(0.51490319) q[3];
sx q[3];
rz(-3.1334936) q[3];
sx q[3];
rz(2.3308339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0136451) q[0];
sx q[0];
rz(-2.7982496) q[0];
sx q[0];
rz(-3.1105594) q[0];
rz(-1.4399928) q[1];
sx q[1];
rz(-0.80039918) q[1];
sx q[1];
rz(-1.4061141) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2314082) q[0];
sx q[0];
rz(-2.6794464) q[0];
sx q[0];
rz(-1.9881835) q[0];
rz(-pi) q[1];
rz(-0.038766301) q[2];
sx q[2];
rz(-1.7749466) q[2];
sx q[2];
rz(-2.7693455) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1202218) q[1];
sx q[1];
rz(-1.7517266) q[1];
sx q[1];
rz(1.4814936) q[1];
rz(1.7201028) q[3];
sx q[3];
rz(-2.341695) q[3];
sx q[3];
rz(0.61198583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.42355737) q[2];
sx q[2];
rz(-0.2090629) q[2];
sx q[2];
rz(-0.62530708) q[2];
rz(1.5710255) q[3];
sx q[3];
rz(-2.8795241) q[3];
sx q[3];
rz(0.03381332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5728979) q[0];
sx q[0];
rz(-0.53302065) q[0];
sx q[0];
rz(2.0030588) q[0];
rz(-2.2472527) q[1];
sx q[1];
rz(-0.00044755512) q[1];
sx q[1];
rz(-2.0785887) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.244449) q[0];
sx q[0];
rz(-2.4384866) q[0];
sx q[0];
rz(1.9528815) q[0];
rz(-pi) q[1];
x q[1];
rz(0.20653609) q[2];
sx q[2];
rz(-0.75876941) q[2];
sx q[2];
rz(-1.9643003) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.059308278) q[1];
sx q[1];
rz(-1.7011146) q[1];
sx q[1];
rz(-3.0695634) q[1];
x q[2];
rz(0.7512666) q[3];
sx q[3];
rz(-0.29956529) q[3];
sx q[3];
rz(1.159016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2669617) q[2];
sx q[2];
rz(-0.78445542) q[2];
sx q[2];
rz(-1.1355012) q[2];
rz(1.1442432) q[3];
sx q[3];
rz(-3.0502697) q[3];
sx q[3];
rz(-1.7448366) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8046232) q[0];
sx q[0];
rz(-2.7380044) q[0];
sx q[0];
rz(-2.023977) q[0];
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
rz(2.8768432) q[0];
sx q[0];
rz(-0.6144982) q[0];
sx q[0];
rz(1.2063998) q[0];
rz(0.20524673) q[2];
sx q[2];
rz(-1.4637498) q[2];
sx q[2];
rz(-2.7595124) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3561894) q[1];
sx q[1];
rz(-1.265824) q[1];
sx q[1];
rz(-3.0126291) q[1];
rz(-2.1104061) q[3];
sx q[3];
rz(-1.0301208) q[3];
sx q[3];
rz(-0.074652076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9522004) q[2];
sx q[2];
rz(-0.11710937) q[2];
sx q[2];
rz(-1.8899151) q[2];
rz(0.39024726) q[3];
sx q[3];
rz(-0.016642112) q[3];
sx q[3];
rz(2.9515631) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8613043) q[0];
sx q[0];
rz(-0.6328634) q[0];
sx q[0];
rz(1.8508152) q[0];
rz(-0.0064918594) q[1];
sx q[1];
rz(-2.8752654) q[1];
sx q[1];
rz(-1.806462) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3399451) q[0];
sx q[0];
rz(-2.772004) q[0];
sx q[0];
rz(0.38529372) q[0];
rz(-pi) q[1];
rz(3.0090121) q[2];
sx q[2];
rz(-0.53756102) q[2];
sx q[2];
rz(2.3258296) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2903394) q[1];
sx q[1];
rz(-1.3604934) q[1];
sx q[1];
rz(1.444534) q[1];
rz(-pi) q[2];
x q[2];
rz(0.18551056) q[3];
sx q[3];
rz(-1.9855567) q[3];
sx q[3];
rz(1.9583553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.73162762) q[2];
sx q[2];
rz(-0.16606398) q[2];
sx q[2];
rz(0.49327332) q[2];
rz(-1.0924245) q[3];
sx q[3];
rz(-3.137393) q[3];
sx q[3];
rz(-0.56001979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6888371) q[0];
sx q[0];
rz(-1.9552564) q[0];
sx q[0];
rz(1.8) q[0];
rz(-1.5780559) q[1];
sx q[1];
rz(-0.23461239) q[1];
sx q[1];
rz(0.87306195) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7811383) q[0];
sx q[0];
rz(-2.6037881) q[0];
sx q[0];
rz(2.8685158) q[0];
rz(-2.701917) q[2];
sx q[2];
rz(-1.5182677) q[2];
sx q[2];
rz(2.479468) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6775408) q[1];
sx q[1];
rz(-1.4960491) q[1];
sx q[1];
rz(-1.6932373) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8382243) q[3];
sx q[3];
rz(-1.7828327) q[3];
sx q[3];
rz(2.7515281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0688429) q[2];
sx q[2];
rz(-1.9708956) q[2];
sx q[2];
rz(-1.076131) q[2];
rz(2.7592646) q[3];
sx q[3];
rz(-3.101442) q[3];
sx q[3];
rz(0.7676777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1165987) q[0];
sx q[0];
rz(-3.0285663) q[0];
sx q[0];
rz(-0.22079994) q[0];
rz(0.91572541) q[1];
sx q[1];
rz(-0.16685367) q[1];
sx q[1];
rz(1.6545666) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2325848) q[0];
sx q[0];
rz(-0.56031432) q[0];
sx q[0];
rz(0.73580297) q[0];
rz(-pi) q[1];
rz(1.9369097) q[2];
sx q[2];
rz(-1.7559474) q[2];
sx q[2];
rz(-2.2529035) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.2714752) q[1];
sx q[1];
rz(-2.041973) q[1];
sx q[1];
rz(0.79939009) q[1];
x q[2];
rz(2.2742918) q[3];
sx q[3];
rz(-0.18531042) q[3];
sx q[3];
rz(-1.0516604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9080092) q[2];
sx q[2];
rz(-0.073069409) q[2];
sx q[2];
rz(1.118534) q[2];
rz(-1.7782878) q[3];
sx q[3];
rz(-0.076863591) q[3];
sx q[3];
rz(1.793687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3240647) q[0];
sx q[0];
rz(-0.20453608) q[0];
sx q[0];
rz(1.8472141) q[0];
rz(1.4140465) q[1];
sx q[1];
rz(-3.0906446) q[1];
sx q[1];
rz(0.22305138) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9004001) q[0];
sx q[0];
rz(-2.0471388) q[0];
sx q[0];
rz(2.4621055) q[0];
rz(0.62533929) q[2];
sx q[2];
rz(-2.9089667) q[2];
sx q[2];
rz(-2.1785506) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1522882) q[1];
sx q[1];
rz(-0.76123255) q[1];
sx q[1];
rz(-2.1851995) q[1];
x q[2];
rz(-3.024142) q[3];
sx q[3];
rz(-0.87908213) q[3];
sx q[3];
rz(-2.1827616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.53905067) q[2];
sx q[2];
rz(-1.9571303) q[2];
sx q[2];
rz(2.7089673) q[2];
rz(1.0299579) q[3];
sx q[3];
rz(-0.10051388) q[3];
sx q[3];
rz(-2.4242134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75281993) q[0];
sx q[0];
rz(-0.50935519) q[0];
sx q[0];
rz(0.21554047) q[0];
rz(-2.0095677) q[1];
sx q[1];
rz(-0.80765453) q[1];
sx q[1];
rz(1.4707123) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4819538) q[0];
sx q[0];
rz(-2.1428272) q[0];
sx q[0];
rz(-2.7670949) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1342513) q[2];
sx q[2];
rz(-1.9160532) q[2];
sx q[2];
rz(-0.79011142) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.017478213) q[1];
sx q[1];
rz(-0.062659293) q[1];
sx q[1];
rz(2.6397328) q[1];
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
rz(0.97951621) q[2];
sx q[2];
rz(-3.1023878) q[2];
sx q[2];
rz(-2.6230679) q[2];
rz(-2.3833852) q[3];
sx q[3];
rz(-2.3617305) q[3];
sx q[3];
rz(1.6583748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6766847) q[0];
sx q[0];
rz(-1.2228891) q[0];
sx q[0];
rz(-1.884888) q[0];
rz(2.8769809) q[1];
sx q[1];
rz(-3.1393317) q[1];
sx q[1];
rz(-1.1983926) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0007858) q[0];
sx q[0];
rz(-0.0065751271) q[0];
sx q[0];
rz(3.0964609) q[0];
x q[1];
rz(1.7670777) q[2];
sx q[2];
rz(-0.33995865) q[2];
sx q[2];
rz(-2.6808471) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9094734) q[1];
sx q[1];
rz(-0.23481242) q[1];
sx q[1];
rz(-0.058583184) q[1];
x q[2];
rz(-1.7890187) q[3];
sx q[3];
rz(-2.9683111) q[3];
sx q[3];
rz(-0.97001433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6780702) q[2];
sx q[2];
rz(-0.0040201298) q[2];
sx q[2];
rz(3.1275952) q[2];
rz(-0.30063453) q[3];
sx q[3];
rz(-3.1379353) q[3];
sx q[3];
rz(1.5716871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8138301) q[0];
sx q[0];
rz(-1.4284842) q[0];
sx q[0];
rz(1.6416657) q[0];
rz(-0.035540237) q[1];
sx q[1];
rz(-2.7062369) q[1];
sx q[1];
rz(-2.8027262) q[1];
rz(2.4443807) q[2];
sx q[2];
rz(-1.3877735) q[2];
sx q[2];
rz(-1.4112926) q[2];
rz(-0.80336848) q[3];
sx q[3];
rz(-1.4252335) q[3];
sx q[3];
rz(-1.2463481) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
