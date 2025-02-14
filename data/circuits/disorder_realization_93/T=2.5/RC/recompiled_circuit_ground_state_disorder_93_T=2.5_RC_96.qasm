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
rz(1.4461203) q[1];
sx q[1];
rz(-1.5674011) q[1];
sx q[1];
rz(-2.8653436) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32466185) q[0];
sx q[0];
rz(-1.766705) q[0];
sx q[0];
rz(-0.017058998) q[0];
x q[1];
rz(-3.0285199) q[2];
sx q[2];
rz(-1.5613174) q[2];
sx q[2];
rz(-2.1533356) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3018926) q[1];
sx q[1];
rz(-0.046161126) q[1];
sx q[1];
rz(0.34997527) q[1];
x q[2];
rz(0.18452355) q[3];
sx q[3];
rz(-0.34466668) q[3];
sx q[3];
rz(-2.7027948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.19286284) q[2];
sx q[2];
rz(-2.7728045) q[2];
sx q[2];
rz(0.31664872) q[2];
rz(2.6266895) q[3];
sx q[3];
rz(-3.1334936) q[3];
sx q[3];
rz(-2.3308339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0136451) q[0];
sx q[0];
rz(-0.34334308) q[0];
sx q[0];
rz(0.031033255) q[0];
rz(-1.4399928) q[1];
sx q[1];
rz(-2.3411935) q[1];
sx q[1];
rz(1.4061141) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91018447) q[0];
sx q[0];
rz(-2.6794464) q[0];
sx q[0];
rz(1.9881835) q[0];
x q[1];
rz(-3.1028264) q[2];
sx q[2];
rz(-1.3666461) q[2];
sx q[2];
rz(0.37224712) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1202218) q[1];
sx q[1];
rz(-1.3898661) q[1];
sx q[1];
rz(-1.4814936) q[1];
rz(-pi) q[2];
rz(-2.3651016) q[3];
sx q[3];
rz(-1.6776976) q[3];
sx q[3];
rz(0.85437802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7180353) q[2];
sx q[2];
rz(-0.2090629) q[2];
sx q[2];
rz(-2.5162856) q[2];
rz(1.5705671) q[3];
sx q[3];
rz(-2.8795241) q[3];
sx q[3];
rz(-0.03381332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5728979) q[0];
sx q[0];
rz(-0.53302065) q[0];
sx q[0];
rz(2.0030588) q[0];
rz(2.2472527) q[1];
sx q[1];
rz(-0.00044755512) q[1];
sx q[1];
rz(2.0785887) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5178058) q[0];
sx q[0];
rz(-1.8142801) q[0];
sx q[0];
rz(0.90434771) q[0];
rz(-pi) q[1];
rz(1.7628305) q[2];
sx q[2];
rz(-0.83198943) q[2];
sx q[2];
rz(-0.89621025) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.059308278) q[1];
sx q[1];
rz(-1.7011146) q[1];
sx q[1];
rz(-0.072029217) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.390326) q[3];
sx q[3];
rz(-0.29956529) q[3];
sx q[3];
rz(-1.9825767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.87463093) q[2];
sx q[2];
rz(-2.3571372) q[2];
sx q[2];
rz(2.0060914) q[2];
rz(1.1442432) q[3];
sx q[3];
rz(-3.0502697) q[3];
sx q[3];
rz(-1.7448366) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3369694) q[0];
sx q[0];
rz(-0.40358821) q[0];
sx q[0];
rz(-1.1176156) q[0];
rz(0.090967372) q[1];
sx q[1];
rz(-0.03799835) q[1];
sx q[1];
rz(1.5100347) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26474947) q[0];
sx q[0];
rz(-0.6144982) q[0];
sx q[0];
rz(-1.9351929) q[0];
rz(-pi) q[1];
rz(-0.20524673) q[2];
sx q[2];
rz(-1.6778429) q[2];
sx q[2];
rz(0.38208026) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7638844) q[1];
sx q[1];
rz(-0.33033366) q[1];
sx q[1];
rz(1.9586455) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.7078739) q[3];
sx q[3];
rz(-0.74433577) q[3];
sx q[3];
rz(-0.78619781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9522004) q[2];
sx q[2];
rz(-3.0244833) q[2];
sx q[2];
rz(1.8899151) q[2];
rz(-0.39024726) q[3];
sx q[3];
rz(-0.016642112) q[3];
sx q[3];
rz(0.19002953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2802884) q[0];
sx q[0];
rz(-0.6328634) q[0];
sx q[0];
rz(1.2907775) q[0];
rz(3.1351008) q[1];
sx q[1];
rz(-2.8752654) q[1];
sx q[1];
rz(1.3351306) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4075942) q[0];
sx q[0];
rz(-1.4346135) q[0];
sx q[0];
rz(-2.7969317) q[0];
x q[1];
rz(-0.13258055) q[2];
sx q[2];
rz(-2.6040316) q[2];
sx q[2];
rz(0.81576306) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8512533) q[1];
sx q[1];
rz(-1.7810993) q[1];
sx q[1];
rz(1.6970587) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9919768) q[3];
sx q[3];
rz(-1.7404181) q[3];
sx q[3];
rz(2.6785525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.409965) q[2];
sx q[2];
rz(-0.16606398) q[2];
sx q[2];
rz(0.49327332) q[2];
rz(-2.0491681) q[3];
sx q[3];
rz(-3.137393) q[3];
sx q[3];
rz(0.56001979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6888371) q[0];
sx q[0];
rz(-1.9552564) q[0];
sx q[0];
rz(1.8) q[0];
rz(-1.5635368) q[1];
sx q[1];
rz(-0.23461239) q[1];
sx q[1];
rz(-0.87306195) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7811383) q[0];
sx q[0];
rz(-0.53780452) q[0];
sx q[0];
rz(2.8685158) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0186924) q[2];
sx q[2];
rz(-0.44259884) q[2];
sx q[2];
rz(2.1217608) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4640518) q[1];
sx q[1];
rz(-1.4960491) q[1];
sx q[1];
rz(1.6932373) q[1];
x q[2];
rz(-1.3033684) q[3];
sx q[3];
rz(-1.7828327) q[3];
sx q[3];
rz(-2.7515281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0688429) q[2];
sx q[2];
rz(-1.9708956) q[2];
sx q[2];
rz(2.0654616) q[2];
rz(-2.7592646) q[3];
sx q[3];
rz(-3.101442) q[3];
sx q[3];
rz(-0.7676777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.024993984) q[0];
sx q[0];
rz(-3.0285663) q[0];
sx q[0];
rz(2.9207927) q[0];
rz(2.2258672) q[1];
sx q[1];
rz(-0.16685367) q[1];
sx q[1];
rz(1.487026) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3160639) q[0];
sx q[0];
rz(-1.9355312) q[0];
sx q[0];
rz(0.43532128) q[0];
x q[1];
rz(1.0887695) q[2];
sx q[2];
rz(-2.7332167) q[2];
sx q[2];
rz(-2.0118304) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2579563) q[1];
sx q[1];
rz(-0.90054926) q[1];
sx q[1];
rz(-2.5237917) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2742918) q[3];
sx q[3];
rz(-2.9562822) q[3];
sx q[3];
rz(-2.0899322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.23358341) q[2];
sx q[2];
rz(-3.0685232) q[2];
sx q[2];
rz(-1.118534) q[2];
rz(-1.3633049) q[3];
sx q[3];
rz(-3.0647291) q[3];
sx q[3];
rz(-1.3479056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(1.8175279) q[0];
sx q[0];
rz(-2.9370566) q[0];
sx q[0];
rz(-1.8472141) q[0];
rz(-1.7275461) q[1];
sx q[1];
rz(-3.0906446) q[1];
sx q[1];
rz(0.22305138) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2411925) q[0];
sx q[0];
rz(-2.0471388) q[0];
sx q[0];
rz(-2.4621055) q[0];
rz(-1.4329918) q[2];
sx q[2];
rz(-1.3827822) q[2];
sx q[2];
rz(-1.540198) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2168426) q[1];
sx q[1];
rz(-0.97198689) q[1];
sx q[1];
rz(0.50227614) q[1];
x q[2];
rz(-0.87568385) q[3];
sx q[3];
rz(-1.4804262) q[3];
sx q[3];
rz(-2.4545057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.53905067) q[2];
sx q[2];
rz(-1.1844623) q[2];
sx q[2];
rz(-0.43262532) q[2];
rz(-1.0299579) q[3];
sx q[3];
rz(-0.10051388) q[3];
sx q[3];
rz(-0.71737927) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
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
rz(1.6708803) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0316353) q[0];
sx q[0];
rz(-2.4695129) q[0];
sx q[0];
rz(1.0541148) q[0];
rz(2.7637761) q[2];
sx q[2];
rz(-1.1616128) q[2];
sx q[2];
rz(-2.5175187) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1241144) q[1];
sx q[1];
rz(-0.062659293) q[1];
sx q[1];
rz(2.6397328) q[1];
rz(-pi) q[2];
rz(0.4827667) q[3];
sx q[3];
rz(-1.1454057) q[3];
sx q[3];
rz(-1.9874043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1620764) q[2];
sx q[2];
rz(-3.1023878) q[2];
sx q[2];
rz(2.6230679) q[2];
rz(2.3833852) q[3];
sx q[3];
rz(-2.3617305) q[3];
sx q[3];
rz(1.4832179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46490797) q[0];
sx q[0];
rz(-1.2228891) q[0];
sx q[0];
rz(1.2567047) q[0];
rz(-0.26461178) q[1];
sx q[1];
rz(-3.1393317) q[1];
sx q[1];
rz(1.9432) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0459185) q[0];
sx q[0];
rz(-1.5642279) q[0];
sx q[0];
rz(1.571093) q[0];
rz(1.2368868) q[2];
sx q[2];
rz(-1.5057202) q[2];
sx q[2];
rz(1.2953616) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2923515) q[1];
sx q[1];
rz(-1.3363943) q[1];
sx q[1];
rz(1.584802) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4015454) q[3];
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
rz(1.6780702) q[2];
sx q[2];
rz(-3.1375725) q[2];
sx q[2];
rz(0.013997495) q[2];
rz(0.30063453) q[3];
sx q[3];
rz(-0.0036573452) q[3];
sx q[3];
rz(1.5716871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32776253) q[0];
sx q[0];
rz(-1.4284842) q[0];
sx q[0];
rz(1.6416657) q[0];
rz(-3.1060524) q[1];
sx q[1];
rz(-0.4353558) q[1];
sx q[1];
rz(0.33886649) q[1];
rz(-1.333893) q[2];
sx q[2];
rz(-0.8874923) q[2];
sx q[2];
rz(0.31077261) q[2];
rz(-0.20094674) q[3];
sx q[3];
rz(-0.81352109) q[3];
sx q[3];
rz(-2.9561694) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
