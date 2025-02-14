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
rz(-0.79688537) q[0];
sx q[0];
rz(4.8719811) q[0];
sx q[0];
rz(13.48579) q[0];
rz(1.3031651) q[1];
sx q[1];
rz(-0.04920955) q[1];
sx q[1];
rz(-0.71370178) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8827191) q[0];
sx q[0];
rz(-1.8451515) q[0];
sx q[0];
rz(0.15299847) q[0];
rz(-1.6993148) q[2];
sx q[2];
rz(-1.5084195) q[2];
sx q[2];
rz(-0.85278748) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9992396) q[1];
sx q[1];
rz(-2.0217784) q[1];
sx q[1];
rz(1.5692479) q[1];
rz(-3.1072989) q[3];
sx q[3];
rz(-1.099713) q[3];
sx q[3];
rz(-0.33609875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0539703) q[2];
sx q[2];
rz(-2.7555608) q[2];
sx q[2];
rz(-0.36960441) q[2];
rz(-1.9965648) q[3];
sx q[3];
rz(-1.6686882) q[3];
sx q[3];
rz(-2.7692774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.10577781) q[0];
sx q[0];
rz(-1.6749629) q[0];
sx q[0];
rz(2.5685487) q[0];
rz(-1.0967968) q[1];
sx q[1];
rz(-1.5410475) q[1];
sx q[1];
rz(-0.47164741) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4147316) q[0];
sx q[0];
rz(-1.4607753) q[0];
sx q[0];
rz(2.1741406) q[0];
x q[1];
rz(2.422633) q[2];
sx q[2];
rz(-2.3078354) q[2];
sx q[2];
rz(1.7938839) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9976366) q[1];
sx q[1];
rz(-0.9846217) q[1];
sx q[1];
rz(-2.2464941) q[1];
rz(-pi) q[2];
x q[2];
rz(0.92949683) q[3];
sx q[3];
rz(-2.5166582) q[3];
sx q[3];
rz(-0.24111023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4196709) q[2];
sx q[2];
rz(-2.2025043) q[2];
sx q[2];
rz(-1.088885) q[2];
rz(2.0770843) q[3];
sx q[3];
rz(-1.1176611) q[3];
sx q[3];
rz(-2.815912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0786781) q[0];
sx q[0];
rz(-2.9525472) q[0];
sx q[0];
rz(-3.1245226) q[0];
rz(2.4315289) q[1];
sx q[1];
rz(-2.3080669) q[1];
sx q[1];
rz(0.56627083) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0294993) q[0];
sx q[0];
rz(-1.5819966) q[0];
sx q[0];
rz(-2.311934) q[0];
x q[1];
rz(-2.5336419) q[2];
sx q[2];
rz(-0.81038108) q[2];
sx q[2];
rz(-1.9160401) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.10668392) q[1];
sx q[1];
rz(-1.7572173) q[1];
sx q[1];
rz(1.4843462) q[1];
rz(-pi) q[2];
rz(-0.2416824) q[3];
sx q[3];
rz(-0.2727601) q[3];
sx q[3];
rz(2.5916416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.337734) q[2];
sx q[2];
rz(-2.2457819) q[2];
sx q[2];
rz(-0.0090948661) q[2];
rz(-3.005262) q[3];
sx q[3];
rz(-2.3970042) q[3];
sx q[3];
rz(-1.2135308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.424054) q[0];
sx q[0];
rz(-0.67988765) q[0];
sx q[0];
rz(2.9754382) q[0];
rz(-2.0280139) q[1];
sx q[1];
rz(-0.49003092) q[1];
sx q[1];
rz(-2.9794433) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55570468) q[0];
sx q[0];
rz(-1.3909719) q[0];
sx q[0];
rz(3.0469045) q[0];
x q[1];
rz(-1.5954564) q[2];
sx q[2];
rz(-1.0219821) q[2];
sx q[2];
rz(-1.2545214) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.051013324) q[1];
sx q[1];
rz(-0.83003488) q[1];
sx q[1];
rz(-0.49511893) q[1];
rz(-pi) q[2];
rz(-2.477263) q[3];
sx q[3];
rz(-1.2264612) q[3];
sx q[3];
rz(-0.20928247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.98264155) q[2];
sx q[2];
rz(-1.0794159) q[2];
sx q[2];
rz(0.52784935) q[2];
rz(0.85150254) q[3];
sx q[3];
rz(-0.32367555) q[3];
sx q[3];
rz(-2.1984656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8714137) q[0];
sx q[0];
rz(-1.542792) q[0];
sx q[0];
rz(2.340509) q[0];
rz(-0.78394765) q[1];
sx q[1];
rz(-2.3990217) q[1];
sx q[1];
rz(-0.98091006) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1005412) q[0];
sx q[0];
rz(-1.9773736) q[0];
sx q[0];
rz(-1.340614) q[0];
rz(2.7655914) q[2];
sx q[2];
rz(-1.094295) q[2];
sx q[2];
rz(-0.27793542) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.66554932) q[1];
sx q[1];
rz(-0.11075039) q[1];
sx q[1];
rz(1.0057463) q[1];
rz(-pi) q[2];
rz(1.1889691) q[3];
sx q[3];
rz(-2.5562048) q[3];
sx q[3];
rz(2.8849734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0826147) q[2];
sx q[2];
rz(-1.5951944) q[2];
sx q[2];
rz(2.2583466) q[2];
rz(0.20032459) q[3];
sx q[3];
rz(-0.89503461) q[3];
sx q[3];
rz(-1.0909572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9969295) q[0];
sx q[0];
rz(-1.0522333) q[0];
sx q[0];
rz(-0.99217478) q[0];
rz(0.80157533) q[1];
sx q[1];
rz(-1.293332) q[1];
sx q[1];
rz(1.7209524) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.84467) q[0];
sx q[0];
rz(-1.8886744) q[0];
sx q[0];
rz(0.65638377) q[0];
rz(-1.0384485) q[2];
sx q[2];
rz(-1.6611929) q[2];
sx q[2];
rz(0.20648512) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0040505) q[1];
sx q[1];
rz(-0.95064771) q[1];
sx q[1];
rz(-1.9778848) q[1];
x q[2];
rz(-0.037461683) q[3];
sx q[3];
rz(-1.473663) q[3];
sx q[3];
rz(-0.85696062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8617323) q[2];
sx q[2];
rz(-1.7143152) q[2];
sx q[2];
rz(2.6045065) q[2];
rz(3.1308657) q[3];
sx q[3];
rz(-2.3948632) q[3];
sx q[3];
rz(0.85095325) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2456197) q[0];
sx q[0];
rz(-2.8697822) q[0];
sx q[0];
rz(0.33988345) q[0];
rz(-1.3019568) q[1];
sx q[1];
rz(-0.94568959) q[1];
sx q[1];
rz(-1.1713015) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3816649) q[0];
sx q[0];
rz(-1.6070843) q[0];
sx q[0];
rz(2.8445303) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7843855) q[2];
sx q[2];
rz(-1.7860951) q[2];
sx q[2];
rz(2.4395296) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4170467) q[1];
sx q[1];
rz(-1.5929211) q[1];
sx q[1];
rz(1.538365) q[1];
x q[2];
rz(0.62452353) q[3];
sx q[3];
rz(-1.4281264) q[3];
sx q[3];
rz(2.4017911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.12477144) q[2];
sx q[2];
rz(-1.7243959) q[2];
sx q[2];
rz(2.4701414) q[2];
rz(-0.5365544) q[3];
sx q[3];
rz(-0.96009976) q[3];
sx q[3];
rz(-1.051739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8272098) q[0];
sx q[0];
rz(-1.0541414) q[0];
sx q[0];
rz(-2.2453454) q[0];
rz(0.25262901) q[1];
sx q[1];
rz(-0.8600421) q[1];
sx q[1];
rz(1.0505229) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80879656) q[0];
sx q[0];
rz(-1.044863) q[0];
sx q[0];
rz(-1.512708) q[0];
rz(-pi) q[1];
rz(1.4398378) q[2];
sx q[2];
rz(-0.90990657) q[2];
sx q[2];
rz(0.63643989) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.76681644) q[1];
sx q[1];
rz(-1.4207867) q[1];
sx q[1];
rz(-1.1039724) q[1];
x q[2];
rz(-1.518431) q[3];
sx q[3];
rz(-0.81365055) q[3];
sx q[3];
rz(-2.6770075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.62319055) q[2];
sx q[2];
rz(-2.9727327) q[2];
sx q[2];
rz(1.5625578) q[2];
rz(-1.7744428) q[3];
sx q[3];
rz(-1.046215) q[3];
sx q[3];
rz(-0.76213837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52687454) q[0];
sx q[0];
rz(-1.7521097) q[0];
sx q[0];
rz(-0.38994625) q[0];
rz(-2.3155616) q[1];
sx q[1];
rz(-0.69458687) q[1];
sx q[1];
rz(1.0521851) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56985939) q[0];
sx q[0];
rz(-0.89825199) q[0];
sx q[0];
rz(2.74617) q[0];
x q[1];
rz(2.7446943) q[2];
sx q[2];
rz(-1.4087311) q[2];
sx q[2];
rz(-1.6610749) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3996399) q[1];
sx q[1];
rz(-1.9583079) q[1];
sx q[1];
rz(2.8273929) q[1];
rz(0.94433208) q[3];
sx q[3];
rz(-1.9667224) q[3];
sx q[3];
rz(-2.348071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0803926) q[2];
sx q[2];
rz(-1.281851) q[2];
sx q[2];
rz(-2.6123987) q[2];
rz(-0.75096327) q[3];
sx q[3];
rz(-2.1801528) q[3];
sx q[3];
rz(-0.19415893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26218364) q[0];
sx q[0];
rz(-0.59469596) q[0];
sx q[0];
rz(-1.9770812) q[0];
rz(-1.2871845) q[1];
sx q[1];
rz(-0.6856122) q[1];
sx q[1];
rz(2.6374292) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39162779) q[0];
sx q[0];
rz(-2.7996783) q[0];
sx q[0];
rz(2.2750399) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5515366) q[2];
sx q[2];
rz(-1.1941205) q[2];
sx q[2];
rz(-2.6835359) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2626942) q[1];
sx q[1];
rz(-2.8126908) q[1];
sx q[1];
rz(0.94543381) q[1];
rz(-pi) q[2];
x q[2];
rz(0.21728094) q[3];
sx q[3];
rz(-2.2437527) q[3];
sx q[3];
rz(0.5289883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.21504271) q[2];
sx q[2];
rz(-1.5861009) q[2];
sx q[2];
rz(-3.1070993) q[2];
rz(-1.5283594) q[3];
sx q[3];
rz(-2.4210763) q[3];
sx q[3];
rz(0.31149402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.16019776) q[0];
sx q[0];
rz(-1.5934516) q[0];
sx q[0];
rz(-0.39504575) q[0];
rz(2.1389217) q[1];
sx q[1];
rz(-1.6802588) q[1];
sx q[1];
rz(-0.37793876) q[1];
rz(-2.7591697) q[2];
sx q[2];
rz(-0.48880807) q[2];
sx q[2];
rz(0.53570408) q[2];
rz(0.32714365) q[3];
sx q[3];
rz(-1.5113304) q[3];
sx q[3];
rz(-2.3979901) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
