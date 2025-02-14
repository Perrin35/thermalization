OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.941551) q[0];
sx q[0];
rz(-0.18414944) q[0];
sx q[0];
rz(1.7848236) q[0];
rz(1.4461203) q[1];
sx q[1];
rz(-1.5674011) q[1];
sx q[1];
rz(0.27624908) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8169308) q[0];
sx q[0];
rz(-1.766705) q[0];
sx q[0];
rz(-3.1245337) q[0];
rz(-pi) q[1];
rz(3.0285199) q[2];
sx q[2];
rz(-1.5802752) q[2];
sx q[2];
rz(0.98825708) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3018926) q[1];
sx q[1];
rz(-3.0954315) q[1];
sx q[1];
rz(0.34997527) q[1];
rz(1.5050234) q[3];
sx q[3];
rz(-1.9093752) q[3];
sx q[3];
rz(2.8985646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.19286284) q[2];
sx q[2];
rz(-2.7728045) q[2];
sx q[2];
rz(-2.8249439) q[2];
rz(0.51490319) q[3];
sx q[3];
rz(-0.0080990214) q[3];
sx q[3];
rz(-2.3308339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0136451) q[0];
sx q[0];
rz(-0.34334308) q[0];
sx q[0];
rz(-0.031033255) q[0];
rz(1.4399928) q[1];
sx q[1];
rz(-0.80039918) q[1];
sx q[1];
rz(1.4061141) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8588327) q[0];
sx q[0];
rz(-1.7525391) q[0];
sx q[0];
rz(1.1434929) q[0];
x q[1];
rz(-0.038766301) q[2];
sx q[2];
rz(-1.7749466) q[2];
sx q[2];
rz(-2.7693455) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5819491) q[1];
sx q[1];
rz(-2.94004) q[1];
sx q[1];
rz(-2.6879716) q[1];
x q[2];
rz(-2.9896432) q[3];
sx q[3];
rz(-0.78228509) q[3];
sx q[3];
rz(-2.3169828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.42355737) q[2];
sx q[2];
rz(-0.2090629) q[2];
sx q[2];
rz(0.62530708) q[2];
rz(-1.5705671) q[3];
sx q[3];
rz(-0.26206854) q[3];
sx q[3];
rz(3.1077793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5728979) q[0];
sx q[0];
rz(-2.608572) q[0];
sx q[0];
rz(2.0030588) q[0];
rz(2.2472527) q[1];
sx q[1];
rz(-0.00044755512) q[1];
sx q[1];
rz(2.0785887) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5178058) q[0];
sx q[0];
rz(-1.3273125) q[0];
sx q[0];
rz(2.2372449) q[0];
rz(-pi) q[1];
rz(2.9350566) q[2];
sx q[2];
rz(-0.75876941) q[2];
sx q[2];
rz(1.9643003) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5021119) q[1];
sx q[1];
rz(-1.6422137) q[1];
sx q[1];
rz(-1.4401431) q[1];
rz(-pi) q[2];
rz(2.9195905) q[3];
sx q[3];
rz(-1.3679805) q[3];
sx q[3];
rz(1.1403644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2669617) q[2];
sx q[2];
rz(-2.3571372) q[2];
sx q[2];
rz(2.0060914) q[2];
rz(1.9973495) q[3];
sx q[3];
rz(-0.091322986) q[3];
sx q[3];
rz(-1.7448366) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3369694) q[0];
sx q[0];
rz(-0.40358821) q[0];
sx q[0];
rz(1.1176156) q[0];
rz(0.090967372) q[1];
sx q[1];
rz(-3.1035943) q[1];
sx q[1];
rz(-1.5100347) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.172013) q[0];
sx q[0];
rz(-1.0019127) q[0];
sx q[0];
rz(-2.8952231) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6801198) q[2];
sx q[2];
rz(-1.7748516) q[2];
sx q[2];
rz(1.2109546) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3561894) q[1];
sx q[1];
rz(-1.8757687) q[1];
sx q[1];
rz(3.0126291) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0311866) q[3];
sx q[3];
rz(-1.0301208) q[3];
sx q[3];
rz(-3.0669406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1893922) q[2];
sx q[2];
rz(-3.0244833) q[2];
sx q[2];
rz(-1.8899151) q[2];
rz(-0.39024726) q[3];
sx q[3];
rz(-0.016642112) q[3];
sx q[3];
rz(0.19002953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2802884) q[0];
sx q[0];
rz(-2.5087293) q[0];
sx q[0];
rz(1.8508152) q[0];
rz(-3.1351008) q[1];
sx q[1];
rz(-0.26632729) q[1];
sx q[1];
rz(1.3351306) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3399451) q[0];
sx q[0];
rz(-0.3695887) q[0];
sx q[0];
rz(-0.38529372) q[0];
rz(-pi) q[1];
rz(-1.4921564) q[2];
sx q[2];
rz(-2.1031339) q[2];
sx q[2];
rz(0.66173205) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.74400988) q[1];
sx q[1];
rz(-2.8967794) q[1];
sx q[1];
rz(0.53298612) q[1];
rz(1.1496159) q[3];
sx q[3];
rz(-1.4011746) q[3];
sx q[3];
rz(2.6785525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.73162762) q[2];
sx q[2];
rz(-2.9755287) q[2];
sx q[2];
rz(0.49327332) q[2];
rz(-2.0491681) q[3];
sx q[3];
rz(-3.137393) q[3];
sx q[3];
rz(-2.5815729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6888371) q[0];
sx q[0];
rz(-1.9552564) q[0];
sx q[0];
rz(-1.8) q[0];
rz(-1.5780559) q[1];
sx q[1];
rz(-2.9069803) q[1];
sx q[1];
rz(2.2685307) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36045433) q[0];
sx q[0];
rz(-2.6037881) q[0];
sx q[0];
rz(-0.27307684) q[0];
x q[1];
rz(3.0186924) q[2];
sx q[2];
rz(-2.6989938) q[2];
sx q[2];
rz(-2.1217608) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.11593379) q[1];
sx q[1];
rz(-1.6928937) q[1];
sx q[1];
rz(3.0662838) q[1];
x q[2];
rz(-0.21960658) q[3];
sx q[3];
rz(-1.3094995) q[3];
sx q[3];
rz(1.1231339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0727498) q[2];
sx q[2];
rz(-1.170697) q[2];
sx q[2];
rz(-2.0654616) q[2];
rz(-0.38232803) q[3];
sx q[3];
rz(-3.101442) q[3];
sx q[3];
rz(0.7676777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.024993984) q[0];
sx q[0];
rz(-0.11302639) q[0];
sx q[0];
rz(-2.9207927) q[0];
rz(-2.2258672) q[1];
sx q[1];
rz(-2.974739) q[1];
sx q[1];
rz(1.487026) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0512569) q[0];
sx q[0];
rz(-1.1658586) q[0];
sx q[0];
rz(1.1722408) q[0];
rz(-pi) q[1];
rz(-1.0887695) q[2];
sx q[2];
rz(-0.408376) q[2];
sx q[2];
rz(1.1297623) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.88363632) q[1];
sx q[1];
rz(-0.90054926) q[1];
sx q[1];
rz(-0.61780091) q[1];
x q[2];
rz(2.2742918) q[3];
sx q[3];
rz(-0.18531042) q[3];
sx q[3];
rz(2.0899322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9080092) q[2];
sx q[2];
rz(-0.073069409) q[2];
sx q[2];
rz(2.0230587) q[2];
rz(-1.7782878) q[3];
sx q[3];
rz(-0.076863591) q[3];
sx q[3];
rz(-1.3479056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68434277) q[0];
sx q[0];
rz(-0.97829223) q[0];
sx q[0];
rz(0.98513435) q[0];
rz(-0.62533929) q[2];
sx q[2];
rz(-2.9089667) q[2];
sx q[2];
rz(-0.96304202) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0878764) q[1];
sx q[1];
rz(-1.9797592) q[1];
sx q[1];
rz(2.2323204) q[1];
rz(-pi) q[2];
rz(-1.7113481) q[3];
sx q[3];
rz(-0.69999123) q[3];
sx q[3];
rz(-2.3656781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.53905067) q[2];
sx q[2];
rz(-1.9571303) q[2];
sx q[2];
rz(2.7089673) q[2];
rz(1.0299579) q[3];
sx q[3];
rz(-3.0410788) q[3];
sx q[3];
rz(-0.71737927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75281993) q[0];
sx q[0];
rz(-2.6322375) q[0];
sx q[0];
rz(-0.21554047) q[0];
rz(2.0095677) q[1];
sx q[1];
rz(-0.80765453) q[1];
sx q[1];
rz(-1.4707123) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12080321) q[0];
sx q[0];
rz(-1.258158) q[0];
sx q[0];
rz(-0.96560574) q[0];
rz(2.2756654) q[2];
sx q[2];
rz(-0.54958639) q[2];
sx q[2];
rz(-1.7333502) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.017478213) q[1];
sx q[1];
rz(-0.062659293) q[1];
sx q[1];
rz(2.6397328) q[1];
x q[2];
rz(-2.3683794) q[3];
sx q[3];
rz(-2.5094731) q[3];
sx q[3];
rz(-0.25037294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1620764) q[2];
sx q[2];
rz(-0.039204892) q[2];
sx q[2];
rz(0.51852477) q[2];
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
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6766847) q[0];
sx q[0];
rz(-1.2228891) q[0];
sx q[0];
rz(1.884888) q[0];
rz(2.8769809) q[1];
sx q[1];
rz(-0.0022609641) q[1];
sx q[1];
rz(-1.9432) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6664724) q[0];
sx q[0];
rz(-1.571093) q[0];
sx q[0];
rz(-3.1350242) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.374515) q[2];
sx q[2];
rz(-0.33995865) q[2];
sx q[2];
rz(-2.6808471) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.28169802) q[1];
sx q[1];
rz(-1.5571737) q[1];
sx q[1];
rz(0.23442421) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7400473) q[3];
sx q[3];
rz(-1.5334606) q[3];
sx q[3];
rz(-0.38572714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4635224) q[2];
sx q[2];
rz(-0.0040201298) q[2];
sx q[2];
rz(-3.1275952) q[2];
rz(-2.8409581) q[3];
sx q[3];
rz(-0.0036573452) q[3];
sx q[3];
rz(-1.5699056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-2.8138301) q[0];
sx q[0];
rz(-1.4284842) q[0];
sx q[0];
rz(1.6416657) q[0];
rz(0.035540237) q[1];
sx q[1];
rz(-0.4353558) q[1];
sx q[1];
rz(0.33886649) q[1];
rz(-0.69721194) q[2];
sx q[2];
rz(-1.3877735) q[2];
sx q[2];
rz(-1.4112926) q[2];
rz(-2.3382242) q[3];
sx q[3];
rz(-1.7163591) q[3];
sx q[3];
rz(1.8952445) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
