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
rz(-1.6662958) q[0];
sx q[0];
rz(-1.8721606) q[0];
sx q[0];
rz(2.4694634) q[0];
rz(0.44828662) q[1];
sx q[1];
rz(-1.5223794) q[1];
sx q[1];
rz(0.7575922) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8286228) q[0];
sx q[0];
rz(-1.4951473) q[0];
sx q[0];
rz(-1.4609758) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.937061) q[2];
sx q[2];
rz(-1.0984813) q[2];
sx q[2];
rz(-0.83845316) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4861826) q[1];
sx q[1];
rz(-1.3453801) q[1];
sx q[1];
rz(2.7214526) q[1];
rz(-pi) q[2];
rz(2.4663062) q[3];
sx q[3];
rz(-0.8199586) q[3];
sx q[3];
rz(-0.91778558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0932833) q[2];
sx q[2];
rz(-1.3292686) q[2];
sx q[2];
rz(1.7993571) q[2];
rz(-1.7736769) q[3];
sx q[3];
rz(-2.048309) q[3];
sx q[3];
rz(1.5526519) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20464483) q[0];
sx q[0];
rz(-1.0138252) q[0];
sx q[0];
rz(2.192705) q[0];
rz(-3.0401547) q[1];
sx q[1];
rz(-2.0815492) q[1];
sx q[1];
rz(-0.92996517) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0895871) q[0];
sx q[0];
rz(-1.6723335) q[0];
sx q[0];
rz(0.2137645) q[0];
rz(0.013596046) q[2];
sx q[2];
rz(-2.6363723) q[2];
sx q[2];
rz(1.9590953) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8023493) q[1];
sx q[1];
rz(-1.4096469) q[1];
sx q[1];
rz(1.4646962) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.61701507) q[3];
sx q[3];
rz(-2.2301144) q[3];
sx q[3];
rz(-0.78711817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.13964222) q[2];
sx q[2];
rz(-1.6763687) q[2];
sx q[2];
rz(0.742221) q[2];
rz(-0.46418515) q[3];
sx q[3];
rz(-1.3451385) q[3];
sx q[3];
rz(1.5884885) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4699698) q[0];
sx q[0];
rz(-0.18593423) q[0];
sx q[0];
rz(-1.6999014) q[0];
rz(2.9761159) q[1];
sx q[1];
rz(-0.73935699) q[1];
sx q[1];
rz(-2.2742719) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.353426) q[0];
sx q[0];
rz(-1.5702308) q[0];
sx q[0];
rz(-3.1415523) q[0];
x q[1];
rz(-1.1050842) q[2];
sx q[2];
rz(-2.5869479) q[2];
sx q[2];
rz(1.9618193) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1731889) q[1];
sx q[1];
rz(-0.43668742) q[1];
sx q[1];
rz(-0.71370947) q[1];
rz(-2.1380566) q[3];
sx q[3];
rz(-2.2972882) q[3];
sx q[3];
rz(-0.15641016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1778339) q[2];
sx q[2];
rz(-1.1744262) q[2];
sx q[2];
rz(0.7589232) q[2];
rz(-2.3376076) q[3];
sx q[3];
rz(-2.1920429) q[3];
sx q[3];
rz(-0.72833958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8300962) q[0];
sx q[0];
rz(-1.8599956) q[0];
sx q[0];
rz(2.6575644) q[0];
rz(1.5361891) q[1];
sx q[1];
rz(-1.7264629) q[1];
sx q[1];
rz(-1.4253634) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7361476) q[0];
sx q[0];
rz(-2.510294) q[0];
sx q[0];
rz(-0.59595705) q[0];
rz(-0.62230939) q[2];
sx q[2];
rz(-0.37823411) q[2];
sx q[2];
rz(2.6985199) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0565383) q[1];
sx q[1];
rz(-0.94496545) q[1];
sx q[1];
rz(0.27295785) q[1];
x q[2];
rz(0.45141545) q[3];
sx q[3];
rz(-1.2299996) q[3];
sx q[3];
rz(-1.6926852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.75590762) q[2];
sx q[2];
rz(-2.0734831) q[2];
sx q[2];
rz(0.29328406) q[2];
rz(-0.048132345) q[3];
sx q[3];
rz(-1.8354514) q[3];
sx q[3];
rz(-0.3046681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83679477) q[0];
sx q[0];
rz(-1.4076819) q[0];
sx q[0];
rz(1.1037214) q[0];
rz(-2.2301105) q[1];
sx q[1];
rz(-1.5244923) q[1];
sx q[1];
rz(0.10890659) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1128588) q[0];
sx q[0];
rz(-2.2153531) q[0];
sx q[0];
rz(2.1875406) q[0];
rz(-pi) q[1];
rz(-2.8418903) q[2];
sx q[2];
rz(-2.5863918) q[2];
sx q[2];
rz(-0.19367684) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2786633) q[1];
sx q[1];
rz(-0.68611523) q[1];
sx q[1];
rz(-1.3115694) q[1];
x q[2];
rz(-0.44556983) q[3];
sx q[3];
rz(-2.516149) q[3];
sx q[3];
rz(-1.7991964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4917422) q[2];
sx q[2];
rz(-1.3567341) q[2];
sx q[2];
rz(-2.8835127) q[2];
rz(0.38170013) q[3];
sx q[3];
rz(-2.2038867) q[3];
sx q[3];
rz(-2.5842353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3968762) q[0];
sx q[0];
rz(-0.98137403) q[0];
sx q[0];
rz(-1.9816403) q[0];
rz(2.5634735) q[1];
sx q[1];
rz(-1.6696397) q[1];
sx q[1];
rz(2.8181308) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8628629) q[0];
sx q[0];
rz(-2.5040309) q[0];
sx q[0];
rz(2.3414073) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6450364) q[2];
sx q[2];
rz(-1.7944031) q[2];
sx q[2];
rz(-1.8211435) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0486794) q[1];
sx q[1];
rz(-1.3447176) q[1];
sx q[1];
rz(0.58633713) q[1];
x q[2];
rz(0.46053912) q[3];
sx q[3];
rz(-1.8611188) q[3];
sx q[3];
rz(-3.0825305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9385927) q[2];
sx q[2];
rz(-2.3679569) q[2];
sx q[2];
rz(-0.8482376) q[2];
rz(1.2674468) q[3];
sx q[3];
rz(-1.7402612) q[3];
sx q[3];
rz(-2.447824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11442014) q[0];
sx q[0];
rz(-2.4878451) q[0];
sx q[0];
rz(0.87345901) q[0];
rz(-1.9050441) q[1];
sx q[1];
rz(-2.1658587) q[1];
sx q[1];
rz(-3.0580318) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77503428) q[0];
sx q[0];
rz(-1.8904443) q[0];
sx q[0];
rz(-1.4421706) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1764063) q[2];
sx q[2];
rz(-2.5766234) q[2];
sx q[2];
rz(0.51630563) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.096683) q[1];
sx q[1];
rz(-0.88777855) q[1];
sx q[1];
rz(-2.6234496) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.6912937) q[3];
sx q[3];
rz(-2.6861827) q[3];
sx q[3];
rz(-0.56870715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7053335) q[2];
sx q[2];
rz(-1.6062364) q[2];
sx q[2];
rz(-1.7748888) q[2];
rz(0.27240917) q[3];
sx q[3];
rz(-1.0206157) q[3];
sx q[3];
rz(-0.75850707) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2347539) q[0];
sx q[0];
rz(-0.82643569) q[0];
sx q[0];
rz(-0.93836623) q[0];
rz(2.3387108) q[1];
sx q[1];
rz(-2.1736841) q[1];
sx q[1];
rz(-2.3209007) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7308189) q[0];
sx q[0];
rz(-1.2026498) q[0];
sx q[0];
rz(-0.14342043) q[0];
rz(0.16751473) q[2];
sx q[2];
rz(-1.2021087) q[2];
sx q[2];
rz(2.9978254) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.332224) q[1];
sx q[1];
rz(-1.6767772) q[1];
sx q[1];
rz(-0.94632228) q[1];
rz(-pi) q[2];
rz(1.1701297) q[3];
sx q[3];
rz(-1.7224215) q[3];
sx q[3];
rz(2.0057297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.20251033) q[2];
sx q[2];
rz(-0.15190092) q[2];
sx q[2];
rz(2.0738156) q[2];
rz(2.0104525) q[3];
sx q[3];
rz(-0.83429566) q[3];
sx q[3];
rz(-0.080032674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.985567) q[0];
sx q[0];
rz(-1.9859059) q[0];
sx q[0];
rz(2.735403) q[0];
rz(1.73229) q[1];
sx q[1];
rz(-1.0180165) q[1];
sx q[1];
rz(-1.0850151) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.800348) q[0];
sx q[0];
rz(-0.64697504) q[0];
sx q[0];
rz(-0.2917618) q[0];
rz(-pi) q[1];
rz(-1.5753393) q[2];
sx q[2];
rz(-0.96049236) q[2];
sx q[2];
rz(1.6729599) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.20445261) q[1];
sx q[1];
rz(-1.8984183) q[1];
sx q[1];
rz(2.3069068) q[1];
rz(2.5625251) q[3];
sx q[3];
rz(-2.0430312) q[3];
sx q[3];
rz(-3.0621665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5857508) q[2];
sx q[2];
rz(-2.5297574) q[2];
sx q[2];
rz(2.2108868) q[2];
rz(-1.7454923) q[3];
sx q[3];
rz(-1.3627108) q[3];
sx q[3];
rz(-3.0220368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68468204) q[0];
sx q[0];
rz(-2.1645808) q[0];
sx q[0];
rz(1.3071625) q[0];
rz(-2.7556509) q[1];
sx q[1];
rz(-1.7470876) q[1];
sx q[1];
rz(-1.0345667) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2725814) q[0];
sx q[0];
rz(-2.9391461) q[0];
sx q[0];
rz(1.0342717) q[0];
rz(-1.7770281) q[2];
sx q[2];
rz(-1.1185371) q[2];
sx q[2];
rz(0.79712501) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.93232357) q[1];
sx q[1];
rz(-1.1205562) q[1];
sx q[1];
rz(-2.2850014) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7004556) q[3];
sx q[3];
rz(-2.377452) q[3];
sx q[3];
rz(1.0636927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4474386) q[2];
sx q[2];
rz(-2.4093781) q[2];
sx q[2];
rz(2.7745957) q[2];
rz(-2.0224109) q[3];
sx q[3];
rz(-2.7435591) q[3];
sx q[3];
rz(-1.5252349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3601396) q[0];
sx q[0];
rz(-1.9522788) q[0];
sx q[0];
rz(2.0576394) q[0];
rz(-3.0837334) q[1];
sx q[1];
rz(-1.2670988) q[1];
sx q[1];
rz(1.7115464) q[1];
rz(1.6690785) q[2];
sx q[2];
rz(-1.3986893) q[2];
sx q[2];
rz(-1.324011) q[2];
rz(-0.69969768) q[3];
sx q[3];
rz(-1.2917923) q[3];
sx q[3];
rz(2.7421799) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
