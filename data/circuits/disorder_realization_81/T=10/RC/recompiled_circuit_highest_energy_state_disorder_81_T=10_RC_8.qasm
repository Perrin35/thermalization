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
rz(0.089224815) q[0];
sx q[0];
rz(-0.74639809) q[0];
sx q[0];
rz(-2.449692) q[0];
rz(0.46299419) q[1];
sx q[1];
rz(-1.6166592) q[1];
sx q[1];
rz(-1.0301627) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93744103) q[0];
sx q[0];
rz(-0.61249706) q[0];
sx q[0];
rz(-2.7322311) q[0];
rz(-pi) q[1];
rz(1.9857731) q[2];
sx q[2];
rz(-0.43104592) q[2];
sx q[2];
rz(-2.8243714) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9249484) q[1];
sx q[1];
rz(-2.7949667) q[1];
sx q[1];
rz(-1.9199172) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6923156) q[3];
sx q[3];
rz(-0.88551023) q[3];
sx q[3];
rz(-0.44997643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2618711) q[2];
sx q[2];
rz(-2.3679831) q[2];
sx q[2];
rz(-0.36641463) q[2];
rz(-2.053818) q[3];
sx q[3];
rz(-2.352114) q[3];
sx q[3];
rz(-2.4521949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3159897) q[0];
sx q[0];
rz(-0.082254224) q[0];
sx q[0];
rz(2.2218623) q[0];
rz(2.3986744) q[1];
sx q[1];
rz(-1.2580322) q[1];
sx q[1];
rz(2.0222576) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18765981) q[0];
sx q[0];
rz(-1.4260248) q[0];
sx q[0];
rz(0.06186084) q[0];
rz(1.7529923) q[2];
sx q[2];
rz(-0.86002738) q[2];
sx q[2];
rz(1.7882535) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.169848) q[1];
sx q[1];
rz(-1.3339284) q[1];
sx q[1];
rz(0.74044098) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1503937) q[3];
sx q[3];
rz(-1.1288191) q[3];
sx q[3];
rz(-0.22721618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4982345) q[2];
sx q[2];
rz(-1.1946119) q[2];
sx q[2];
rz(3.0894792) q[2];
rz(-1.107996) q[3];
sx q[3];
rz(-0.85927695) q[3];
sx q[3];
rz(-0.064854709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68010083) q[0];
sx q[0];
rz(-1.9100186) q[0];
sx q[0];
rz(1.289806) q[0];
rz(0.75311226) q[1];
sx q[1];
rz(-1.6878004) q[1];
sx q[1];
rz(2.6444816) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6996993) q[0];
sx q[0];
rz(-1.9632247) q[0];
sx q[0];
rz(-2.5470069) q[0];
rz(2.2154664) q[2];
sx q[2];
rz(-2.1417649) q[2];
sx q[2];
rz(-2.1059011) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.048084413) q[1];
sx q[1];
rz(-1.1968379) q[1];
sx q[1];
rz(-1.5300586) q[1];
x q[2];
rz(0.65264355) q[3];
sx q[3];
rz(-1.4941998) q[3];
sx q[3];
rz(-2.9216539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6734753) q[2];
sx q[2];
rz(-0.27286369) q[2];
sx q[2];
rz(2.4172778) q[2];
rz(-2.916548) q[3];
sx q[3];
rz(-1.2757755) q[3];
sx q[3];
rz(0.79506522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.986079) q[0];
sx q[0];
rz(-0.88548311) q[0];
sx q[0];
rz(-0.28050637) q[0];
rz(-1.4836503) q[1];
sx q[1];
rz(-1.9397441) q[1];
sx q[1];
rz(2.6703506) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12898937) q[0];
sx q[0];
rz(-1.5935197) q[0];
sx q[0];
rz(1.9381592) q[0];
rz(-pi) q[1];
rz(1.4578826) q[2];
sx q[2];
rz(-1.6364003) q[2];
sx q[2];
rz(-2.0173617) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.9556826) q[1];
sx q[1];
rz(-1.5985399) q[1];
sx q[1];
rz(-1.8415336) q[1];
rz(-0.97906487) q[3];
sx q[3];
rz(-2.2869898) q[3];
sx q[3];
rz(1.5997052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9435141) q[2];
sx q[2];
rz(-2.5835865) q[2];
sx q[2];
rz(-1.0151218) q[2];
rz(-0.73457581) q[3];
sx q[3];
rz(-1.3646804) q[3];
sx q[3];
rz(-0.53810292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15767844) q[0];
sx q[0];
rz(-1.9240802) q[0];
sx q[0];
rz(2.1814046) q[0];
rz(0.045711191) q[1];
sx q[1];
rz(-2.263133) q[1];
sx q[1];
rz(0.033230573) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8321182) q[0];
sx q[0];
rz(-0.72369472) q[0];
sx q[0];
rz(0.81731082) q[0];
rz(-pi) q[1];
rz(1.6258036) q[2];
sx q[2];
rz(-0.26794456) q[2];
sx q[2];
rz(2.1017113) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7864125) q[1];
sx q[1];
rz(-0.35874736) q[1];
sx q[1];
rz(-1.9500109) q[1];
rz(-pi) q[2];
rz(0.73992336) q[3];
sx q[3];
rz(-2.5725274) q[3];
sx q[3];
rz(-0.25404938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1632605) q[2];
sx q[2];
rz(-0.99908081) q[2];
sx q[2];
rz(0.31326374) q[2];
rz(2.7836109) q[3];
sx q[3];
rz(-2.101818) q[3];
sx q[3];
rz(0.09859214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98704308) q[0];
sx q[0];
rz(-1.6313666) q[0];
sx q[0];
rz(2.6556067) q[0];
rz(-1.9588574) q[1];
sx q[1];
rz(-2.4458838) q[1];
sx q[1];
rz(2.8114496) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3257623) q[0];
sx q[0];
rz(-1.6005524) q[0];
sx q[0];
rz(1.8682006) q[0];
x q[1];
rz(0.80578226) q[2];
sx q[2];
rz(-1.0892765) q[2];
sx q[2];
rz(1.0891506) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3482326) q[1];
sx q[1];
rz(-0.72472445) q[1];
sx q[1];
rz(2.3824005) q[1];
rz(-2.7948912) q[3];
sx q[3];
rz(-1.3606405) q[3];
sx q[3];
rz(-0.29276785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5064454) q[2];
sx q[2];
rz(-0.61656419) q[2];
sx q[2];
rz(0.35663024) q[2];
rz(2.5470274) q[3];
sx q[3];
rz(-1.4831355) q[3];
sx q[3];
rz(2.4026292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4643788) q[0];
sx q[0];
rz(-0.76736275) q[0];
sx q[0];
rz(-2.5546524) q[0];
rz(-2.7691973) q[1];
sx q[1];
rz(-1.3601235) q[1];
sx q[1];
rz(-0.24758235) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1598561) q[0];
sx q[0];
rz(-1.1061108) q[0];
sx q[0];
rz(-3.1064731) q[0];
rz(1.814194) q[2];
sx q[2];
rz(-2.1888424) q[2];
sx q[2];
rz(2.7400561) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7786488) q[1];
sx q[1];
rz(-1.2374068) q[1];
sx q[1];
rz(-0.89493805) q[1];
x q[2];
rz(1.8195175) q[3];
sx q[3];
rz(-1.910672) q[3];
sx q[3];
rz(-2.9074977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8783012) q[2];
sx q[2];
rz(-1.2211439) q[2];
sx q[2];
rz(0.2090052) q[2];
rz(-0.59669295) q[3];
sx q[3];
rz(-0.34691072) q[3];
sx q[3];
rz(-1.6056812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7335032) q[0];
sx q[0];
rz(-0.96036378) q[0];
sx q[0];
rz(0.39304131) q[0];
rz(0.6706925) q[1];
sx q[1];
rz(-1.1436983) q[1];
sx q[1];
rz(-1.447698) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5825621) q[0];
sx q[0];
rz(-0.5120638) q[0];
sx q[0];
rz(-1.8265695) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1426439) q[2];
sx q[2];
rz(-2.482126) q[2];
sx q[2];
rz(-3.1127549) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.32998521) q[1];
sx q[1];
rz(-1.0701655) q[1];
sx q[1];
rz(-0.58793427) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0441229) q[3];
sx q[3];
rz(-0.83126634) q[3];
sx q[3];
rz(2.748718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4286917) q[2];
sx q[2];
rz(-0.8608326) q[2];
sx q[2];
rz(3.0248896) q[2];
rz(-2.1278925) q[3];
sx q[3];
rz(-1.9096392) q[3];
sx q[3];
rz(-1.8359756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3882465) q[0];
sx q[0];
rz(-1.6755063) q[0];
sx q[0];
rz(1.4353132) q[0];
rz(0.995579) q[1];
sx q[1];
rz(-1.3187871) q[1];
sx q[1];
rz(-0.54042712) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30975809) q[0];
sx q[0];
rz(-1.09344) q[0];
sx q[0];
rz(-1.274513) q[0];
rz(-pi) q[1];
rz(-0.58806555) q[2];
sx q[2];
rz(-2.1508475) q[2];
sx q[2];
rz(-1.5813476) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3167773) q[1];
sx q[1];
rz(-1.1939923) q[1];
sx q[1];
rz(2.3945859) q[1];
rz(-pi) q[2];
rz(-2.5164817) q[3];
sx q[3];
rz(-1.2338936) q[3];
sx q[3];
rz(0.57527104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.871375) q[2];
sx q[2];
rz(-1.8367218) q[2];
sx q[2];
rz(-2.3649141) q[2];
rz(-3.1089879) q[3];
sx q[3];
rz(-1.6083345) q[3];
sx q[3];
rz(-1.5672654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33083522) q[0];
sx q[0];
rz(-0.36892712) q[0];
sx q[0];
rz(-0.43357968) q[0];
rz(2.4724204) q[1];
sx q[1];
rz(-1.1829665) q[1];
sx q[1];
rz(-0.42632595) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2992914) q[0];
sx q[0];
rz(-2.8482245) q[0];
sx q[0];
rz(2.3175453) q[0];
rz(2.3947515) q[2];
sx q[2];
rz(-1.3051118) q[2];
sx q[2];
rz(-2.8248252) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5475216) q[1];
sx q[1];
rz(-1.2202377) q[1];
sx q[1];
rz(2.622753) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0570404) q[3];
sx q[3];
rz(-1.0412239) q[3];
sx q[3];
rz(-2.2467006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3907889) q[2];
sx q[2];
rz(-2.8533253) q[2];
sx q[2];
rz(2.0354347) q[2];
rz(3.028051) q[3];
sx q[3];
rz(-2.2083486) q[3];
sx q[3];
rz(0.043702628) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9857367) q[0];
sx q[0];
rz(-2.2423797) q[0];
sx q[0];
rz(2.5806497) q[0];
rz(2.7583495) q[1];
sx q[1];
rz(-2.0189197) q[1];
sx q[1];
rz(-0.22519208) q[1];
rz(-0.46950317) q[2];
sx q[2];
rz(-2.532685) q[2];
sx q[2];
rz(-1.572123) q[2];
rz(-2.0416971) q[3];
sx q[3];
rz(-0.80977898) q[3];
sx q[3];
rz(1.5017189) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
