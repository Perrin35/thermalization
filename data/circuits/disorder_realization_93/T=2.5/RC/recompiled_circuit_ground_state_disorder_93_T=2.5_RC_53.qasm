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
rz(-1.3567691) q[0];
rz(-1.6954724) q[1];
sx q[1];
rz(-1.5741916) q[1];
sx q[1];
rz(2.8653436) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8987792) q[0];
sx q[0];
rz(-1.5540637) q[0];
sx q[0];
rz(-1.7667328) q[0];
rz(-0.11307271) q[2];
sx q[2];
rz(-1.5613174) q[2];
sx q[2];
rz(-0.98825708) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3018926) q[1];
sx q[1];
rz(-0.046161126) q[1];
sx q[1];
rz(-0.34997527) q[1];
rz(-1.6365693) q[3];
sx q[3];
rz(-1.9093752) q[3];
sx q[3];
rz(2.8985646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.19286284) q[2];
sx q[2];
rz(-2.7728045) q[2];
sx q[2];
rz(2.8249439) q[2];
rz(-2.6266895) q[3];
sx q[3];
rz(-0.0080990214) q[3];
sx q[3];
rz(0.81075877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1279476) q[0];
sx q[0];
rz(-0.34334308) q[0];
sx q[0];
rz(3.1105594) q[0];
rz(-1.7015999) q[1];
sx q[1];
rz(-0.80039918) q[1];
sx q[1];
rz(1.4061141) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3701551) q[0];
sx q[0];
rz(-1.9906128) q[0];
sx q[0];
rz(-2.9423454) q[0];
rz(-pi) q[1];
rz(0.038766301) q[2];
sx q[2];
rz(-1.3666461) q[2];
sx q[2];
rz(0.37224712) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0213709) q[1];
sx q[1];
rz(-1.3898661) q[1];
sx q[1];
rz(-1.660099) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3651016) q[3];
sx q[3];
rz(-1.463895) q[3];
sx q[3];
rz(-2.2872146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7180353) q[2];
sx q[2];
rz(-2.9325298) q[2];
sx q[2];
rz(-2.5162856) q[2];
rz(-1.5710255) q[3];
sx q[3];
rz(-0.26206854) q[3];
sx q[3];
rz(0.03381332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5728979) q[0];
sx q[0];
rz(-2.608572) q[0];
sx q[0];
rz(1.1385338) q[0];
rz(-0.89433995) q[1];
sx q[1];
rz(-3.1411451) q[1];
sx q[1];
rz(1.0630039) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8971436) q[0];
sx q[0];
rz(-0.70310601) q[0];
sx q[0];
rz(1.9528815) q[0];
x q[1];
rz(1.7628305) q[2];
sx q[2];
rz(-2.3096032) q[2];
sx q[2];
rz(-2.2453824) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.059308278) q[1];
sx q[1];
rz(-1.7011146) q[1];
sx q[1];
rz(-3.0695634) q[1];
x q[2];
rz(-1.3630217) q[3];
sx q[3];
rz(-1.3534184) q[3];
sx q[3];
rz(-0.3849963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.87463093) q[2];
sx q[2];
rz(-0.78445542) q[2];
sx q[2];
rz(-2.0060914) q[2];
rz(-1.1442432) q[3];
sx q[3];
rz(-0.091322986) q[3];
sx q[3];
rz(1.3967561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8046232) q[0];
sx q[0];
rz(-0.40358821) q[0];
sx q[0];
rz(1.1176156) q[0];
rz(-3.0506253) q[1];
sx q[1];
rz(-3.1035943) q[1];
sx q[1];
rz(-1.5100347) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5334329) q[0];
sx q[0];
rz(-1.3638487) q[0];
sx q[0];
rz(-2.1536835) q[0];
x q[1];
rz(-0.20524673) q[2];
sx q[2];
rz(-1.4637498) q[2];
sx q[2];
rz(2.7595124) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7638844) q[1];
sx q[1];
rz(-2.811259) q[1];
sx q[1];
rz(-1.1829472) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0311866) q[3];
sx q[3];
rz(-2.1114719) q[3];
sx q[3];
rz(-0.074652076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9522004) q[2];
sx q[2];
rz(-0.11710937) q[2];
sx q[2];
rz(-1.8899151) q[2];
rz(0.39024726) q[3];
sx q[3];
rz(-3.1249505) q[3];
sx q[3];
rz(-2.9515631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8613043) q[0];
sx q[0];
rz(-2.5087293) q[0];
sx q[0];
rz(1.8508152) q[0];
rz(0.0064918594) q[1];
sx q[1];
rz(-0.26632729) q[1];
sx q[1];
rz(1.3351306) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80164755) q[0];
sx q[0];
rz(-2.772004) q[0];
sx q[0];
rz(0.38529372) q[0];
x q[1];
rz(-3.0090121) q[2];
sx q[2];
rz(-2.6040316) q[2];
sx q[2];
rz(2.3258296) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2903394) q[1];
sx q[1];
rz(-1.7810993) q[1];
sx q[1];
rz(-1.6970587) q[1];
rz(-0.18551056) q[3];
sx q[3];
rz(-1.1560359) q[3];
sx q[3];
rz(1.9583553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.73162762) q[2];
sx q[2];
rz(-2.9755287) q[2];
sx q[2];
rz(-2.6483193) q[2];
rz(-1.0924245) q[3];
sx q[3];
rz(-0.0041996669) q[3];
sx q[3];
rz(-2.5815729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6888371) q[0];
sx q[0];
rz(-1.1863363) q[0];
sx q[0];
rz(1.8) q[0];
rz(1.5635368) q[1];
sx q[1];
rz(-2.9069803) q[1];
sx q[1];
rz(2.2685307) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4463947) q[0];
sx q[0];
rz(-1.4322011) q[0];
sx q[0];
rz(2.6202484) q[0];
x q[1];
rz(1.5127584) q[2];
sx q[2];
rz(-1.1317694) q[2];
sx q[2];
rz(0.88397938) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.4387008) q[1];
sx q[1];
rz(-0.14335614) q[1];
sx q[1];
rz(1.0207671) q[1];
x q[2];
rz(-0.88718702) q[3];
sx q[3];
rz(-0.33970133) q[3];
sx q[3];
rz(-1.3058794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0688429) q[2];
sx q[2];
rz(-1.170697) q[2];
sx q[2];
rz(1.076131) q[2];
rz(2.7592646) q[3];
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
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1165987) q[0];
sx q[0];
rz(-0.11302639) q[0];
sx q[0];
rz(0.22079994) q[0];
rz(-0.91572541) q[1];
sx q[1];
rz(-2.974739) q[1];
sx q[1];
rz(-1.487026) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90900786) q[0];
sx q[0];
rz(-0.56031432) q[0];
sx q[0];
rz(-0.73580297) q[0];
x q[1];
rz(-2.9436297) q[2];
sx q[2];
rz(-1.930365) q[2];
sx q[2];
rz(0.61164226) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4055192) q[1];
sx q[1];
rz(-2.2637475) q[1];
sx q[1];
rz(-2.2018717) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2742918) q[3];
sx q[3];
rz(-0.18531042) q[3];
sx q[3];
rz(-1.0516604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.23358341) q[2];
sx q[2];
rz(-3.0685232) q[2];
sx q[2];
rz(1.118534) q[2];
rz(-1.7782878) q[3];
sx q[3];
rz(-3.0647291) q[3];
sx q[3];
rz(-1.793687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3240647) q[0];
sx q[0];
rz(-2.9370566) q[0];
sx q[0];
rz(1.8472141) q[0];
rz(-1.4140465) q[1];
sx q[1];
rz(-0.050948016) q[1];
sx q[1];
rz(0.22305138) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9004001) q[0];
sx q[0];
rz(-1.0944538) q[0];
sx q[0];
rz(-2.4621055) q[0];
x q[1];
rz(0.62533929) q[2];
sx q[2];
rz(-2.9089667) q[2];
sx q[2];
rz(-2.1785506) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9893045) q[1];
sx q[1];
rz(-0.76123255) q[1];
sx q[1];
rz(2.1851995) q[1];
rz(-pi) q[2];
rz(-0.87568385) q[3];
sx q[3];
rz(-1.6611665) q[3];
sx q[3];
rz(2.4545057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.53905067) q[2];
sx q[2];
rz(-1.9571303) q[2];
sx q[2];
rz(-2.7089673) q[2];
rz(1.0299579) q[3];
sx q[3];
rz(-0.10051388) q[3];
sx q[3];
rz(-2.4242134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75281993) q[0];
sx q[0];
rz(-0.50935519) q[0];
sx q[0];
rz(-2.9260522) q[0];
rz(-1.132025) q[1];
sx q[1];
rz(-2.3339381) q[1];
sx q[1];
rz(-1.6708803) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0316353) q[0];
sx q[0];
rz(-0.67207974) q[0];
sx q[0];
rz(2.0874778) q[0];
rz(-2.0073414) q[2];
sx q[2];
rz(-1.9160532) q[2];
sx q[2];
rz(-0.79011142) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.1241144) q[1];
sx q[1];
rz(-0.062659293) q[1];
sx q[1];
rz(-2.6397328) q[1];
rz(-pi) q[2];
rz(-2.3683794) q[3];
sx q[3];
rz(-2.5094731) q[3];
sx q[3];
rz(2.8912197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1620764) q[2];
sx q[2];
rz(-0.039204892) q[2];
sx q[2];
rz(2.6230679) q[2];
rz(-2.3833852) q[3];
sx q[3];
rz(-2.3617305) q[3];
sx q[3];
rz(-1.4832179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46490797) q[0];
sx q[0];
rz(-1.9187036) q[0];
sx q[0];
rz(1.2567047) q[0];
rz(-0.26461178) q[1];
sx q[1];
rz(-3.1393317) q[1];
sx q[1];
rz(1.9432) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4751202) q[0];
sx q[0];
rz(-1.571093) q[0];
sx q[0];
rz(-3.1350242) q[0];
rz(-1.2368868) q[2];
sx q[2];
rz(-1.6358725) q[2];
sx q[2];
rz(-1.8462311) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9094734) q[1];
sx q[1];
rz(-2.9067802) q[1];
sx q[1];
rz(-3.0830095) q[1];
rz(-pi) q[2];
rz(1.7400473) q[3];
sx q[3];
rz(-1.608132) q[3];
sx q[3];
rz(2.7558655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4635224) q[2];
sx q[2];
rz(-0.0040201298) q[2];
sx q[2];
rz(-3.1275952) q[2];
rz(0.30063453) q[3];
sx q[3];
rz(-3.1379353) q[3];
sx q[3];
rz(1.5699056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(0.28066228) q[2];
sx q[2];
rz(-2.4246695) q[2];
sx q[2];
rz(-0.05447745) q[2];
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
