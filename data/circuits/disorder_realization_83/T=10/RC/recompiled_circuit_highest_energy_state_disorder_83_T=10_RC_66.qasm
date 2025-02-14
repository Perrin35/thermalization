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
rz(-0.42521617) q[0];
sx q[0];
rz(4.4758237) q[0];
sx q[0];
rz(9.0444179) q[0];
rz(1.7827787) q[1];
sx q[1];
rz(-0.18667297) q[1];
sx q[1];
rz(-2.2169854) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7036669) q[0];
sx q[0];
rz(-0.80090085) q[0];
sx q[0];
rz(-2.6573703) q[0];
rz(-2.9745462) q[2];
sx q[2];
rz(-1.0606597) q[2];
sx q[2];
rz(1.8710006) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8421665) q[1];
sx q[1];
rz(-1.5314845) q[1];
sx q[1];
rz(-2.2996603) q[1];
x q[2];
rz(-1.5798453) q[3];
sx q[3];
rz(-1.442872) q[3];
sx q[3];
rz(1.8198899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0164464) q[2];
sx q[2];
rz(-0.95637286) q[2];
sx q[2];
rz(-2.1535786) q[2];
rz(2.9347349) q[3];
sx q[3];
rz(-1.7346953) q[3];
sx q[3];
rz(-2.701581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4619231) q[0];
sx q[0];
rz(-1.4689057) q[0];
sx q[0];
rz(2.8894506) q[0];
rz(-0.02154669) q[1];
sx q[1];
rz(-0.69308678) q[1];
sx q[1];
rz(-2.2861939) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48777521) q[0];
sx q[0];
rz(-3.1145373) q[0];
sx q[0];
rz(1.6271126) q[0];
rz(-0.62521817) q[2];
sx q[2];
rz(-2.4897235) q[2];
sx q[2];
rz(2.9232581) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2829078) q[1];
sx q[1];
rz(-2.5577684) q[1];
sx q[1];
rz(-1.9782542) q[1];
rz(0.68287373) q[3];
sx q[3];
rz(-1.8337909) q[3];
sx q[3];
rz(-2.9337286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1713193) q[2];
sx q[2];
rz(-3.1372941) q[2];
sx q[2];
rz(-2.8805736) q[2];
rz(3.1018992) q[3];
sx q[3];
rz(-1.7429765) q[3];
sx q[3];
rz(-0.54350129) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4259341) q[0];
sx q[0];
rz(-1.4733227) q[0];
sx q[0];
rz(-0.69450992) q[0];
rz(2.014324) q[1];
sx q[1];
rz(-2.290461) q[1];
sx q[1];
rz(0.99004254) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58082579) q[0];
sx q[0];
rz(-0.35835727) q[0];
sx q[0];
rz(2.906515) q[0];
rz(-pi) q[1];
rz(-0.98857356) q[2];
sx q[2];
rz(-1.4043839) q[2];
sx q[2];
rz(-2.5022282) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7496365) q[1];
sx q[1];
rz(-2.5534494) q[1];
sx q[1];
rz(0.68444499) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6022791) q[3];
sx q[3];
rz(-0.94180543) q[3];
sx q[3];
rz(-1.1719538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8433044) q[2];
sx q[2];
rz(-0.97323209) q[2];
sx q[2];
rz(-0.49015552) q[2];
rz(-2.7847024) q[3];
sx q[3];
rz(-0.39339742) q[3];
sx q[3];
rz(1.8753768) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.056034293) q[0];
sx q[0];
rz(-1.9558676) q[0];
sx q[0];
rz(0.84247843) q[0];
rz(2.339824) q[1];
sx q[1];
rz(-0.27920488) q[1];
sx q[1];
rz(0.78559771) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88670896) q[0];
sx q[0];
rz(-2.2629316) q[0];
sx q[0];
rz(1.9002302) q[0];
rz(-pi) q[1];
rz(-0.17274348) q[2];
sx q[2];
rz(-2.4170791) q[2];
sx q[2];
rz(2.8959664) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3093058) q[1];
sx q[1];
rz(-2.5046299) q[1];
sx q[1];
rz(-0.5354521) q[1];
x q[2];
rz(0.20468851) q[3];
sx q[3];
rz(-2.31143) q[3];
sx q[3];
rz(1.1799174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0164612) q[2];
sx q[2];
rz(-2.3931914) q[2];
sx q[2];
rz(-1.008519) q[2];
rz(-0.85159167) q[3];
sx q[3];
rz(-2.3840756) q[3];
sx q[3];
rz(2.0612702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9652902) q[0];
sx q[0];
rz(-2.1602614) q[0];
sx q[0];
rz(-1.0705795) q[0];
rz(-2.3640682) q[1];
sx q[1];
rz(-2.2309525) q[1];
sx q[1];
rz(-2.1308965) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.010653822) q[0];
sx q[0];
rz(-1.9294327) q[0];
sx q[0];
rz(2.6333195) q[0];
rz(-pi) q[1];
rz(1.3959051) q[2];
sx q[2];
rz(-1.5064459) q[2];
sx q[2];
rz(-1.5140057) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4324716) q[1];
sx q[1];
rz(-0.84408954) q[1];
sx q[1];
rz(-0.62950397) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.15701954) q[3];
sx q[3];
rz(-2.5738705) q[3];
sx q[3];
rz(-2.1322676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8396478) q[2];
sx q[2];
rz(-0.10493111) q[2];
sx q[2];
rz(-1.5860175) q[2];
rz(-2.1258866) q[3];
sx q[3];
rz(-1.1414707) q[3];
sx q[3];
rz(1.4122081) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3509336) q[0];
sx q[0];
rz(-0.18336329) q[0];
sx q[0];
rz(2.3405128) q[0];
rz(0.13310295) q[1];
sx q[1];
rz(-1.8201273) q[1];
sx q[1];
rz(2.7395693) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9511356) q[0];
sx q[0];
rz(-1.9044838) q[0];
sx q[0];
rz(-0.78898375) q[0];
rz(-pi) q[1];
rz(1.5779675) q[2];
sx q[2];
rz(-0.29080393) q[2];
sx q[2];
rz(-2.4343642) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1014598) q[1];
sx q[1];
rz(-1.8454843) q[1];
sx q[1];
rz(1.0708315) q[1];
rz(-pi) q[2];
rz(-2.4834255) q[3];
sx q[3];
rz(-1.4716545) q[3];
sx q[3];
rz(2.1226573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5037527) q[2];
sx q[2];
rz(-1.4053586) q[2];
sx q[2];
rz(0.77581882) q[2];
rz(1.6860298) q[3];
sx q[3];
rz(-1.9676696) q[3];
sx q[3];
rz(-2.1337401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0966454) q[0];
sx q[0];
rz(-0.89770397) q[0];
sx q[0];
rz(2.4626379) q[0];
rz(1.8567765) q[1];
sx q[1];
rz(-0.7130475) q[1];
sx q[1];
rz(-1.3791893) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5011936) q[0];
sx q[0];
rz(-1.3974579) q[0];
sx q[0];
rz(-2.3480183) q[0];
rz(-pi) q[1];
rz(3.0581362) q[2];
sx q[2];
rz(-1.8493358) q[2];
sx q[2];
rz(2.4065774) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.13869446) q[1];
sx q[1];
rz(-2.0636673) q[1];
sx q[1];
rz(-1.0652131) q[1];
rz(-pi) q[2];
rz(0.055886397) q[3];
sx q[3];
rz(-1.7706141) q[3];
sx q[3];
rz(-0.36147396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3369559) q[2];
sx q[2];
rz(-2.3865484) q[2];
sx q[2];
rz(-1.2389368) q[2];
rz(1.3321446) q[3];
sx q[3];
rz(-2.381031) q[3];
sx q[3];
rz(2.8968887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4181344) q[0];
sx q[0];
rz(-1.9754388) q[0];
sx q[0];
rz(0.13846692) q[0];
rz(-0.18374099) q[1];
sx q[1];
rz(-2.467149) q[1];
sx q[1];
rz(2.8750681) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0040695) q[0];
sx q[0];
rz(-2.7752987) q[0];
sx q[0];
rz(1.1031325) q[0];
x q[1];
rz(2.1378072) q[2];
sx q[2];
rz(-1.3409753) q[2];
sx q[2];
rz(0.15574317) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3764972) q[1];
sx q[1];
rz(-1.9899211) q[1];
sx q[1];
rz(0.061582743) q[1];
x q[2];
rz(1.5795535) q[3];
sx q[3];
rz(-2.2976795) q[3];
sx q[3];
rz(-1.7576287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0515685) q[2];
sx q[2];
rz(-1.8847621) q[2];
sx q[2];
rz(-2.9070692) q[2];
rz(-1.6656434) q[3];
sx q[3];
rz(-1.539295) q[3];
sx q[3];
rz(-0.75638151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(2.804857) q[0];
sx q[0];
rz(-2.2754301) q[0];
sx q[0];
rz(-0.79972237) q[0];
rz(2.5526478) q[1];
sx q[1];
rz(-2.2942693) q[1];
sx q[1];
rz(-0.2074997) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8556535) q[0];
sx q[0];
rz(-0.41623273) q[0];
sx q[0];
rz(-0.32815423) q[0];
x q[1];
rz(-2.2830711) q[2];
sx q[2];
rz(-1.3337161) q[2];
sx q[2];
rz(-1.3905987) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.60002414) q[1];
sx q[1];
rz(-2.3316521) q[1];
sx q[1];
rz(-1.0066973) q[1];
rz(2.0405681) q[3];
sx q[3];
rz(-1.54514) q[3];
sx q[3];
rz(-0.86650833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3346682) q[2];
sx q[2];
rz(-0.88403264) q[2];
sx q[2];
rz(2.7743288) q[2];
rz(-1.4372829) q[3];
sx q[3];
rz(-1.9673012) q[3];
sx q[3];
rz(-1.1737163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87840286) q[0];
sx q[0];
rz(-2.1387565) q[0];
sx q[0];
rz(2.3314085) q[0];
rz(-2.0267678) q[1];
sx q[1];
rz(-2.7572542) q[1];
sx q[1];
rz(-0.55327639) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5779752) q[0];
sx q[0];
rz(-1.6262691) q[0];
sx q[0];
rz(2.2091921) q[0];
x q[1];
rz(-0.61364321) q[2];
sx q[2];
rz(-2.2479703) q[2];
sx q[2];
rz(2.8751862) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7825494) q[1];
sx q[1];
rz(-1.7573253) q[1];
sx q[1];
rz(1.2839509) q[1];
x q[2];
rz(-0.35240473) q[3];
sx q[3];
rz(-0.81768546) q[3];
sx q[3];
rz(-1.529983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1222003) q[2];
sx q[2];
rz(-2.4330752) q[2];
sx q[2];
rz(-0.23966399) q[2];
rz(1.8137118) q[3];
sx q[3];
rz(-2.0769104) q[3];
sx q[3];
rz(2.6740668) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.079530579) q[0];
sx q[0];
rz(-1.3539599) q[0];
sx q[0];
rz(0.39394105) q[0];
rz(-2.486034) q[1];
sx q[1];
rz(-1.1759023) q[1];
sx q[1];
rz(-2.1942153) q[1];
rz(0.66329931) q[2];
sx q[2];
rz(-0.51021432) q[2];
sx q[2];
rz(-2.9205657) q[2];
rz(2.7648457) q[3];
sx q[3];
rz(-2.527311) q[3];
sx q[3];
rz(-0.75025076) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
