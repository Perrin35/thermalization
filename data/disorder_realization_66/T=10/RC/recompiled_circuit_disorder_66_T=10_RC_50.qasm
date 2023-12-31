OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.80823094) q[0];
sx q[0];
rz(-2.1699177) q[0];
sx q[0];
rz(1.4810286) q[0];
rz(-0.12227585) q[1];
sx q[1];
rz(-0.08637698) q[1];
sx q[1];
rz(0.02286214) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1296596) q[0];
sx q[0];
rz(-1.4225905) q[0];
sx q[0];
rz(-0.11797842) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.68732287) q[2];
sx q[2];
rz(-2.5097373) q[2];
sx q[2];
rz(0.99422115) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.248977) q[1];
sx q[1];
rz(-1.9141344) q[1];
sx q[1];
rz(-2.9825319) q[1];
rz(-0.95211897) q[3];
sx q[3];
rz(-1.3631571) q[3];
sx q[3];
rz(-1.729897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7227398) q[2];
sx q[2];
rz(-0.8272233) q[2];
sx q[2];
rz(-2.375405) q[2];
rz(2.9700759) q[3];
sx q[3];
rz(-2.4092716) q[3];
sx q[3];
rz(2.5550301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8121174) q[0];
sx q[0];
rz(-0.43918878) q[0];
sx q[0];
rz(-3.0549333) q[0];
rz(0.86241972) q[1];
sx q[1];
rz(-2.5984867) q[1];
sx q[1];
rz(-3.0564953) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21762411) q[0];
sx q[0];
rz(-2.1930165) q[0];
sx q[0];
rz(-1.5909965) q[0];
rz(-pi) q[1];
rz(0.30479635) q[2];
sx q[2];
rz(-0.99064231) q[2];
sx q[2];
rz(3.0457029) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.47463402) q[1];
sx q[1];
rz(-0.14161319) q[1];
sx q[1];
rz(2.2255564) q[1];
rz(-pi) q[2];
rz(0.99382932) q[3];
sx q[3];
rz(-2.5964063) q[3];
sx q[3];
rz(-2.546437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8791085) q[2];
sx q[2];
rz(-1.5510677) q[2];
sx q[2];
rz(-1.6764199) q[2];
rz(-1.9654467) q[3];
sx q[3];
rz(-0.90548235) q[3];
sx q[3];
rz(-0.24648497) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6572606) q[0];
sx q[0];
rz(-3.0069139) q[0];
sx q[0];
rz(-2.235967) q[0];
rz(-0.33272818) q[1];
sx q[1];
rz(-0.85488027) q[1];
sx q[1];
rz(1.3844301) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8850088) q[0];
sx q[0];
rz(-2.6589767) q[0];
sx q[0];
rz(-2.9628739) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2134174) q[2];
sx q[2];
rz(-2.536142) q[2];
sx q[2];
rz(1.6124992) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.82921925) q[1];
sx q[1];
rz(-2.3260498) q[1];
sx q[1];
rz(1.6014598) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.609385) q[3];
sx q[3];
rz(-2.0720707) q[3];
sx q[3];
rz(-1.7621079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.24421346) q[2];
sx q[2];
rz(-2.801557) q[2];
sx q[2];
rz(1.8133694) q[2];
rz(-1.3978847) q[3];
sx q[3];
rz(-0.54026794) q[3];
sx q[3];
rz(0.11169294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
rz(-0.06552799) q[0];
sx q[0];
rz(-0.17834839) q[0];
sx q[0];
rz(0.67101014) q[0];
rz(-1.7975851) q[1];
sx q[1];
rz(-1.9799415) q[1];
sx q[1];
rz(-2.9342594) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57732338) q[0];
sx q[0];
rz(-3.0992357) q[0];
sx q[0];
rz(1.7839412) q[0];
rz(1.7447128) q[2];
sx q[2];
rz(-1.6007649) q[2];
sx q[2];
rz(-2.7001691) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.748121) q[1];
sx q[1];
rz(-1.3337743) q[1];
sx q[1];
rz(-0.15740983) q[1];
x q[2];
rz(-1.5563577) q[3];
sx q[3];
rz(-1.0591905) q[3];
sx q[3];
rz(2.5479941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8308782) q[2];
sx q[2];
rz(-1.0881311) q[2];
sx q[2];
rz(-0.80835289) q[2];
rz(2.3338142) q[3];
sx q[3];
rz(-2.419796) q[3];
sx q[3];
rz(-2.9479153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5508674) q[0];
sx q[0];
rz(-0.71722513) q[0];
sx q[0];
rz(0.89548683) q[0];
rz(2.8764309) q[1];
sx q[1];
rz(-2.3064955) q[1];
sx q[1];
rz(2.6079544) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95414872) q[0];
sx q[0];
rz(-1.8549518) q[0];
sx q[0];
rz(0.28676333) q[0];
x q[1];
rz(-1.7970423) q[2];
sx q[2];
rz(-2.3257253) q[2];
sx q[2];
rz(-2.0362542) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0698318) q[1];
sx q[1];
rz(-1.7904209) q[1];
sx q[1];
rz(1.4211618) q[1];
x q[2];
rz(-0.94138937) q[3];
sx q[3];
rz(-1.7377059) q[3];
sx q[3];
rz(-0.0080136673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2395893) q[2];
sx q[2];
rz(-1.5633554) q[2];
sx q[2];
rz(-0.56838244) q[2];
rz(-1.9249632) q[3];
sx q[3];
rz(-0.47074461) q[3];
sx q[3];
rz(0.34354982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50826532) q[0];
sx q[0];
rz(-0.87843043) q[0];
sx q[0];
rz(-1.1761965) q[0];
rz(1.5856702) q[1];
sx q[1];
rz(-2.1891179) q[1];
sx q[1];
rz(1.0046545) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6001119) q[0];
sx q[0];
rz(-1.2280557) q[0];
sx q[0];
rz(-0.39370763) q[0];
rz(2.7447694) q[2];
sx q[2];
rz(-1.4962215) q[2];
sx q[2];
rz(0.65149388) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.55014729) q[1];
sx q[1];
rz(-2.1131556) q[1];
sx q[1];
rz(-1.28246) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5980914) q[3];
sx q[3];
rz(-2.7276162) q[3];
sx q[3];
rz(-0.15347813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4040318) q[2];
sx q[2];
rz(-0.99016756) q[2];
sx q[2];
rz(2.693434) q[2];
rz(-0.29799497) q[3];
sx q[3];
rz(-2.4500676) q[3];
sx q[3];
rz(2.7929849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77666831) q[0];
sx q[0];
rz(-1.8402599) q[0];
sx q[0];
rz(-2.4834494) q[0];
rz(1.8572042) q[1];
sx q[1];
rz(-2.6874266) q[1];
sx q[1];
rz(2.4694494) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37049609) q[0];
sx q[0];
rz(-1.578195) q[0];
sx q[0];
rz(-1.7367944) q[0];
rz(1.1291814) q[2];
sx q[2];
rz(-1.9600944) q[2];
sx q[2];
rz(-2.0814975) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.044683177) q[1];
sx q[1];
rz(-1.9367366) q[1];
sx q[1];
rz(0.52176042) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7147195) q[3];
sx q[3];
rz(-2.2669499) q[3];
sx q[3];
rz(1.5211943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.60823524) q[2];
sx q[2];
rz(-0.63295263) q[2];
sx q[2];
rz(2.7427924) q[2];
rz(2.3170025) q[3];
sx q[3];
rz(-1.7493068) q[3];
sx q[3];
rz(-0.59857541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(0.44052112) q[0];
sx q[0];
rz(-2.7109787) q[0];
sx q[0];
rz(-2.9833802) q[0];
rz(-1.0868866) q[1];
sx q[1];
rz(-2.465076) q[1];
sx q[1];
rz(-2.3349082) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0772484) q[0];
sx q[0];
rz(-1.6883388) q[0];
sx q[0];
rz(2.006152) q[0];
rz(-2.6629692) q[2];
sx q[2];
rz(-2.9106986) q[2];
sx q[2];
rz(0.17639562) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1572691) q[1];
sx q[1];
rz(-1.2829797) q[1];
sx q[1];
rz(-0.70043522) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2699548) q[3];
sx q[3];
rz(-1.5402208) q[3];
sx q[3];
rz(1.9105063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5640101) q[2];
sx q[2];
rz(-1.492604) q[2];
sx q[2];
rz(-0.99964833) q[2];
rz(0.10351652) q[3];
sx q[3];
rz(-0.13705702) q[3];
sx q[3];
rz(2.0172393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(3.0209811) q[0];
sx q[0];
rz(-2.3636901) q[0];
sx q[0];
rz(0.65752423) q[0];
rz(0.28655562) q[1];
sx q[1];
rz(-2.2189238) q[1];
sx q[1];
rz(-2.8009169) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7332471) q[0];
sx q[0];
rz(-0.78662965) q[0];
sx q[0];
rz(1.1128845) q[0];
rz(-0.22206837) q[2];
sx q[2];
rz(-1.0169612) q[2];
sx q[2];
rz(2.1332707) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3714911) q[1];
sx q[1];
rz(-1.0189459) q[1];
sx q[1];
rz(2.9279207) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.70298985) q[3];
sx q[3];
rz(-0.25745108) q[3];
sx q[3];
rz(-1.5233745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3866773) q[2];
sx q[2];
rz(-1.9415386) q[2];
sx q[2];
rz(0.49794751) q[2];
rz(0.30161101) q[3];
sx q[3];
rz(-2.7409654) q[3];
sx q[3];
rz(0.51914674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72934812) q[0];
sx q[0];
rz(-3.0817741) q[0];
sx q[0];
rz(2.8630032) q[0];
rz(2.5623698) q[1];
sx q[1];
rz(-0.93943739) q[1];
sx q[1];
rz(-3.0648807) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68966507) q[0];
sx q[0];
rz(-1.8849012) q[0];
sx q[0];
rz(-0.12977022) q[0];
rz(-2.9701091) q[2];
sx q[2];
rz(-1.0043317) q[2];
sx q[2];
rz(2.7330074) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1288651) q[1];
sx q[1];
rz(-2.7755133) q[1];
sx q[1];
rz(-3.0831343) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4823227) q[3];
sx q[3];
rz(-1.7469329) q[3];
sx q[3];
rz(-2.4000771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.77438337) q[2];
sx q[2];
rz(-0.70720208) q[2];
sx q[2];
rz(-0.45551604) q[2];
rz(2.7101743) q[3];
sx q[3];
rz(-3.016267) q[3];
sx q[3];
rz(-0.063522696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9615622) q[0];
sx q[0];
rz(-1.4008235) q[0];
sx q[0];
rz(0.93749198) q[0];
rz(-0.58615276) q[1];
sx q[1];
rz(-1.6393607) q[1];
sx q[1];
rz(1.6929109) q[1];
rz(1.0473245) q[2];
sx q[2];
rz(-0.53146711) q[2];
sx q[2];
rz(-0.55234595) q[2];
rz(-2.7145731) q[3];
sx q[3];
rz(-1.9058766) q[3];
sx q[3];
rz(1.8591892) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
