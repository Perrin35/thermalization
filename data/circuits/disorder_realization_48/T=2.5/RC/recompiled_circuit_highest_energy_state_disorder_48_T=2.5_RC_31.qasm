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
rz(-0.35204044) q[0];
sx q[0];
rz(-0.8249324) q[0];
sx q[0];
rz(-2.6068249) q[0];
rz(-2.1575902) q[1];
sx q[1];
rz(-0.50552955) q[1];
sx q[1];
rz(-1.2619789) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13684863) q[0];
sx q[0];
rz(-1.3082349) q[0];
sx q[0];
rz(0.62341086) q[0];
rz(0.87902576) q[2];
sx q[2];
rz(-0.75415416) q[2];
sx q[2];
rz(-1.1471495) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.76319219) q[1];
sx q[1];
rz(-2.1175368) q[1];
sx q[1];
rz(-1.631447) q[1];
rz(0.38930157) q[3];
sx q[3];
rz(-1.3817781) q[3];
sx q[3];
rz(2.1019328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.24903211) q[2];
sx q[2];
rz(-2.6292215) q[2];
sx q[2];
rz(1.4830291) q[2];
rz(0.28462166) q[3];
sx q[3];
rz(-2.5183545) q[3];
sx q[3];
rz(2.0774138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8677218) q[0];
sx q[0];
rz(-0.72672788) q[0];
sx q[0];
rz(3.100585) q[0];
rz(-1.9339336) q[1];
sx q[1];
rz(-0.227808) q[1];
sx q[1];
rz(-1.6809195) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8475094) q[0];
sx q[0];
rz(-2.9929711) q[0];
sx q[0];
rz(-2.9411208) q[0];
x q[1];
rz(0.76717214) q[2];
sx q[2];
rz(-1.9843272) q[2];
sx q[2];
rz(-2.3656379) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.25819689) q[1];
sx q[1];
rz(-0.25942311) q[1];
sx q[1];
rz(1.7218462) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8051265) q[3];
sx q[3];
rz(-2.4357492) q[3];
sx q[3];
rz(0.88283759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.66551208) q[2];
sx q[2];
rz(-1.4613232) q[2];
sx q[2];
rz(-1.7512789) q[2];
rz(0.0811854) q[3];
sx q[3];
rz(-0.66892162) q[3];
sx q[3];
rz(-1.4359052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
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
rz(1.7471033) q[0];
sx q[0];
rz(-3.1293226) q[0];
sx q[0];
rz(0.61007208) q[0];
rz(1.3129781) q[1];
sx q[1];
rz(-2.4459631) q[1];
sx q[1];
rz(-2.5909766) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2064821) q[0];
sx q[0];
rz(-1.9294104) q[0];
sx q[0];
rz(2.4367843) q[0];
rz(2.8189893) q[2];
sx q[2];
rz(-1.8292973) q[2];
sx q[2];
rz(2.5532364) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0455835) q[1];
sx q[1];
rz(-2.7852045) q[1];
sx q[1];
rz(2.136904) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.14282823) q[3];
sx q[3];
rz(-2.3108644) q[3];
sx q[3];
rz(-0.21634858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.47784558) q[2];
sx q[2];
rz(-2.9691594) q[2];
sx q[2];
rz(2.0752068) q[2];
rz(1.1586698) q[3];
sx q[3];
rz(-1.3830769) q[3];
sx q[3];
rz(0.042044736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8313507) q[0];
sx q[0];
rz(-0.61315918) q[0];
sx q[0];
rz(0.9170652) q[0];
rz(-2.4861368) q[1];
sx q[1];
rz(-0.90924811) q[1];
sx q[1];
rz(-2.7520032) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.675466) q[0];
sx q[0];
rz(-1.6524757) q[0];
sx q[0];
rz(-1.7577175) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2245313) q[2];
sx q[2];
rz(-0.79981632) q[2];
sx q[2];
rz(-0.66674495) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1234189) q[1];
sx q[1];
rz(-1.5416073) q[1];
sx q[1];
rz(1.1694714) q[1];
x q[2];
rz(-0.54549952) q[3];
sx q[3];
rz(-1.0441458) q[3];
sx q[3];
rz(2.5712476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.83170825) q[2];
sx q[2];
rz(-2.1990621) q[2];
sx q[2];
rz(-1.7892828) q[2];
rz(2.3934707) q[3];
sx q[3];
rz(-1.8375405) q[3];
sx q[3];
rz(-1.5266533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30369502) q[0];
sx q[0];
rz(-2.5717773) q[0];
sx q[0];
rz(1.3414398) q[0];
rz(-2.1603284) q[1];
sx q[1];
rz(-1.860268) q[1];
sx q[1];
rz(0.70995465) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9729346) q[0];
sx q[0];
rz(-0.32635411) q[0];
sx q[0];
rz(-0.059684233) q[0];
rz(-pi) q[1];
rz(-1.9290438) q[2];
sx q[2];
rz(-2.5351371) q[2];
sx q[2];
rz(-2.0403634) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.1011184) q[1];
sx q[1];
rz(-1.7832563) q[1];
sx q[1];
rz(-2.967415) q[1];
rz(1.4758797) q[3];
sx q[3];
rz(-2.1664985) q[3];
sx q[3];
rz(1.9323521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.95825163) q[2];
sx q[2];
rz(-2.3812582) q[2];
sx q[2];
rz(2.23488) q[2];
rz(-2.9592311) q[3];
sx q[3];
rz(-1.7014818) q[3];
sx q[3];
rz(0.85842925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7957423) q[0];
sx q[0];
rz(-2.1868571) q[0];
sx q[0];
rz(0.20172754) q[0];
rz(-1.3564823) q[1];
sx q[1];
rz(-2.4701665) q[1];
sx q[1];
rz(-1.1511525) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4326103) q[0];
sx q[0];
rz(-1.994619) q[0];
sx q[0];
rz(2.9124385) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.62080748) q[2];
sx q[2];
rz(-2.4352031) q[2];
sx q[2];
rz(-0.28677973) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5804739) q[1];
sx q[1];
rz(-2.4137133) q[1];
sx q[1];
rz(-2.5125458) q[1];
rz(0.51064193) q[3];
sx q[3];
rz(-1.770442) q[3];
sx q[3];
rz(-0.24195237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.98728937) q[2];
sx q[2];
rz(-2.5550948) q[2];
sx q[2];
rz(-0.34745535) q[2];
rz(-2.3197428) q[3];
sx q[3];
rz(-1.3653267) q[3];
sx q[3];
rz(1.6177026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3407985) q[0];
sx q[0];
rz(-2.4947385) q[0];
sx q[0];
rz(-1.1232173) q[0];
rz(2.7128291) q[1];
sx q[1];
rz(-2.2518497) q[1];
sx q[1];
rz(2.9926328) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70388795) q[0];
sx q[0];
rz(-1.0539248) q[0];
sx q[0];
rz(-0.025512841) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.5908509) q[2];
sx q[2];
rz(-0.35195165) q[2];
sx q[2];
rz(-2.1972434) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.80854446) q[1];
sx q[1];
rz(-0.77189779) q[1];
sx q[1];
rz(2.9289362) q[1];
rz(-pi) q[2];
rz(3.1066057) q[3];
sx q[3];
rz(-1.1652428) q[3];
sx q[3];
rz(-1.7008894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1296156) q[2];
sx q[2];
rz(-0.35795438) q[2];
sx q[2];
rz(0.32148662) q[2];
rz(2.0251515) q[3];
sx q[3];
rz(-2.2827086) q[3];
sx q[3];
rz(-1.4139676) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1534934) q[0];
sx q[0];
rz(-1.7970002) q[0];
sx q[0];
rz(-0.022911428) q[0];
rz(1.5921536) q[1];
sx q[1];
rz(-1.0433082) q[1];
sx q[1];
rz(-1.2124088) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7717944) q[0];
sx q[0];
rz(-1.3098469) q[0];
sx q[0];
rz(-1.0124932) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3018478) q[2];
sx q[2];
rz(-1.296412) q[2];
sx q[2];
rz(-1.7626761) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.089253292) q[1];
sx q[1];
rz(-1.2683378) q[1];
sx q[1];
rz(0.41771981) q[1];
rz(-pi) q[2];
rz(3.0453835) q[3];
sx q[3];
rz(-2.2819977) q[3];
sx q[3];
rz(2.8195153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2108078) q[2];
sx q[2];
rz(-2.942694) q[2];
sx q[2];
rz(0.99736324) q[2];
rz(1.2584244) q[3];
sx q[3];
rz(-1.9505898) q[3];
sx q[3];
rz(-1.2303801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1579943) q[0];
sx q[0];
rz(-0.65584922) q[0];
sx q[0];
rz(-2.6672145) q[0];
rz(0.55167088) q[1];
sx q[1];
rz(-0.76849476) q[1];
sx q[1];
rz(-0.65753585) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6303355) q[0];
sx q[0];
rz(-1.610582) q[0];
sx q[0];
rz(1.5944832) q[0];
rz(-1.9507061) q[2];
sx q[2];
rz(-0.43336855) q[2];
sx q[2];
rz(0.72335183) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.20831693) q[1];
sx q[1];
rz(-2.7728656) q[1];
sx q[1];
rz(1.6160203) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7153431) q[3];
sx q[3];
rz(-1.4667505) q[3];
sx q[3];
rz(-0.34510976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.23003301) q[2];
sx q[2];
rz(-2.7342789) q[2];
sx q[2];
rz(-1.3113021) q[2];
rz(2.4274872) q[3];
sx q[3];
rz(-1.8592535) q[3];
sx q[3];
rz(2.5717521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0062200935) q[0];
sx q[0];
rz(-2.1923809) q[0];
sx q[0];
rz(-2.5332992) q[0];
rz(0.81781203) q[1];
sx q[1];
rz(-1.9655971) q[1];
sx q[1];
rz(-0.83293319) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3522763) q[0];
sx q[0];
rz(-2.7171405) q[0];
sx q[0];
rz(-3.0158774) q[0];
rz(-1.8274177) q[2];
sx q[2];
rz(-1.875281) q[2];
sx q[2];
rz(1.2198795) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2257833) q[1];
sx q[1];
rz(-1.7034917) q[1];
sx q[1];
rz(2.5578294) q[1];
rz(-pi) q[2];
rz(0.046499297) q[3];
sx q[3];
rz(-2.5480418) q[3];
sx q[3];
rz(1.4607185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1828764) q[2];
sx q[2];
rz(-2.7887838) q[2];
sx q[2];
rz(2.555441) q[2];
rz(2.2385249) q[3];
sx q[3];
rz(-0.62246263) q[3];
sx q[3];
rz(-3.0557475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9027973) q[0];
sx q[0];
rz(-1.1302523) q[0];
sx q[0];
rz(-0.94373066) q[0];
rz(-2.0441652) q[1];
sx q[1];
rz(-0.25249093) q[1];
sx q[1];
rz(2.9846334) q[1];
rz(1.1949933) q[2];
sx q[2];
rz(-2.7766418) q[2];
sx q[2];
rz(3.1212213) q[2];
rz(-0.4914183) q[3];
sx q[3];
rz(-0.93486356) q[3];
sx q[3];
rz(-0.60318666) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
