OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3818504) q[0];
sx q[0];
rz(-0.83431017) q[0];
sx q[0];
rz(-0.068339737) q[0];
rz(2.4812658) q[1];
sx q[1];
rz(3.9897444) q[1];
sx q[1];
rz(6.2453649) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8605182) q[0];
sx q[0];
rz(-1.7639419) q[0];
sx q[0];
rz(0.079695745) q[0];
rz(-1.7224563) q[2];
sx q[2];
rz(-0.66705396) q[2];
sx q[2];
rz(2.5487713) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6112069) q[1];
sx q[1];
rz(-1.5198277) q[1];
sx q[1];
rz(2.3947122) q[1];
x q[2];
rz(1.171265) q[3];
sx q[3];
rz(-0.37326187) q[3];
sx q[3];
rz(1.319862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7906856) q[2];
sx q[2];
rz(-2.6280792) q[2];
sx q[2];
rz(-1.2834056) q[2];
rz(0.12617271) q[3];
sx q[3];
rz(-1.4161371) q[3];
sx q[3];
rz(0.090601966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.693817) q[0];
sx q[0];
rz(-0.87647804) q[0];
sx q[0];
rz(0.59666657) q[0];
rz(-1.5555351) q[1];
sx q[1];
rz(-1.3417473) q[1];
sx q[1];
rz(1.7780875) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1024433) q[0];
sx q[0];
rz(-1.1457232) q[0];
sx q[0];
rz(-2.5545679) q[0];
rz(0.081625799) q[2];
sx q[2];
rz(-0.55742369) q[2];
sx q[2];
rz(-0.73478414) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5433568) q[1];
sx q[1];
rz(-0.9423061) q[1];
sx q[1];
rz(1.223982) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.76916839) q[3];
sx q[3];
rz(-0.90220074) q[3];
sx q[3];
rz(1.234496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3229225) q[2];
sx q[2];
rz(-2.1408036) q[2];
sx q[2];
rz(2.4070516) q[2];
rz(2.5143886) q[3];
sx q[3];
rz(-1.6278798) q[3];
sx q[3];
rz(-0.74497765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7221786) q[0];
sx q[0];
rz(-2.3086771) q[0];
sx q[0];
rz(2.869379) q[0];
rz(-0.84725562) q[1];
sx q[1];
rz(-1.3191185) q[1];
sx q[1];
rz(0.31262696) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49293672) q[0];
sx q[0];
rz(-1.3515633) q[0];
sx q[0];
rz(1.5159025) q[0];
x q[1];
rz(-1.1001415) q[2];
sx q[2];
rz(-1.0725642) q[2];
sx q[2];
rz(-0.91980308) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.15309139) q[1];
sx q[1];
rz(-1.3107436) q[1];
sx q[1];
rz(3.0696763) q[1];
rz(-0.39194312) q[3];
sx q[3];
rz(-0.80229811) q[3];
sx q[3];
rz(-2.9585569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.12038885) q[2];
sx q[2];
rz(-0.52336064) q[2];
sx q[2];
rz(-2.3266501) q[2];
rz(-0.96873823) q[3];
sx q[3];
rz(-2.0961943) q[3];
sx q[3];
rz(-2.708784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3880436) q[0];
sx q[0];
rz(-2.9409445) q[0];
sx q[0];
rz(2.5182305) q[0];
rz(0.81758824) q[1];
sx q[1];
rz(-0.31660429) q[1];
sx q[1];
rz(3.0923016) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.277963) q[0];
sx q[0];
rz(-2.3290714) q[0];
sx q[0];
rz(-0.066594007) q[0];
rz(2.2823015) q[2];
sx q[2];
rz(-1.8154732) q[2];
sx q[2];
rz(1.7161075) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.1218159) q[1];
sx q[1];
rz(-0.41408086) q[1];
sx q[1];
rz(2.2112234) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.44932511) q[3];
sx q[3];
rz(-2.0530564) q[3];
sx q[3];
rz(-1.1376017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4108882) q[2];
sx q[2];
rz(-1.5521908) q[2];
sx q[2];
rz(-2.1566186) q[2];
rz(-2.2385521) q[3];
sx q[3];
rz(-1.1629546) q[3];
sx q[3];
rz(-2.8682017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34489283) q[0];
sx q[0];
rz(-1.6051689) q[0];
sx q[0];
rz(-3.0539736) q[0];
rz(0.15631974) q[1];
sx q[1];
rz(-2.5904398) q[1];
sx q[1];
rz(0.87096754) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6773416) q[0];
sx q[0];
rz(-0.28872492) q[0];
sx q[0];
rz(0.10084734) q[0];
rz(-pi) q[1];
rz(-0.15337495) q[2];
sx q[2];
rz(-1.9605637) q[2];
sx q[2];
rz(-0.49488059) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4082143) q[1];
sx q[1];
rz(-0.92792643) q[1];
sx q[1];
rz(-1.2182359) q[1];
rz(-2.9669168) q[3];
sx q[3];
rz(-0.81411618) q[3];
sx q[3];
rz(-2.0811618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0481723) q[2];
sx q[2];
rz(-2.293736) q[2];
sx q[2];
rz(-0.35287738) q[2];
rz(-2.2757163) q[3];
sx q[3];
rz(-2.4398949) q[3];
sx q[3];
rz(-0.29233366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6784994) q[0];
sx q[0];
rz(-2.0149639) q[0];
sx q[0];
rz(3.034806) q[0];
rz(1.1865901) q[1];
sx q[1];
rz(-2.1148966) q[1];
sx q[1];
rz(-1.0245163) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2811919) q[0];
sx q[0];
rz(-1.8587347) q[0];
sx q[0];
rz(-0.47173758) q[0];
rz(-2.0265987) q[2];
sx q[2];
rz(-2.1264646) q[2];
sx q[2];
rz(-1.9089886) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5441106) q[1];
sx q[1];
rz(-2.1253617) q[1];
sx q[1];
rz(-2.5764562) q[1];
x q[2];
rz(0.49196524) q[3];
sx q[3];
rz(-1.2428189) q[3];
sx q[3];
rz(-2.5043951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7496877) q[2];
sx q[2];
rz(-1.2140032) q[2];
sx q[2];
rz(2.270703) q[2];
rz(1.8093367) q[3];
sx q[3];
rz(-2.199316) q[3];
sx q[3];
rz(-1.4107305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1720599) q[0];
sx q[0];
rz(-1.4071858) q[0];
sx q[0];
rz(-2.5906738) q[0];
rz(2.6761966) q[1];
sx q[1];
rz(-2.4702563) q[1];
sx q[1];
rz(-0.23682061) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2315002) q[0];
sx q[0];
rz(-1.5331368) q[0];
sx q[0];
rz(1.670027) q[0];
rz(-pi) q[1];
rz(2.4004164) q[2];
sx q[2];
rz(-0.20222649) q[2];
sx q[2];
rz(1.9066332) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9501818) q[1];
sx q[1];
rz(-1.6457874) q[1];
sx q[1];
rz(-1.300315) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2348433) q[3];
sx q[3];
rz(-0.7630322) q[3];
sx q[3];
rz(1.8900161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6992496) q[2];
sx q[2];
rz(-2.4599059) q[2];
sx q[2];
rz(0.44011763) q[2];
rz(1.0951428) q[3];
sx q[3];
rz(-1.6710072) q[3];
sx q[3];
rz(-1.7949036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31297627) q[0];
sx q[0];
rz(-1.799311) q[0];
sx q[0];
rz(-3.112088) q[0];
rz(2.1942031) q[1];
sx q[1];
rz(-2.3209929) q[1];
sx q[1];
rz(3.0227919) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1664365) q[0];
sx q[0];
rz(-1.6140037) q[0];
sx q[0];
rz(-2.8309566) q[0];
rz(1.0256836) q[2];
sx q[2];
rz(-1.9695749) q[2];
sx q[2];
rz(2.039197) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.58662215) q[1];
sx q[1];
rz(-1.8645617) q[1];
sx q[1];
rz(-2.657386) q[1];
rz(-pi) q[2];
rz(0.16468594) q[3];
sx q[3];
rz(-2.0593004) q[3];
sx q[3];
rz(0.96795852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3999195) q[2];
sx q[2];
rz(-0.46367773) q[2];
sx q[2];
rz(-1.5734394) q[2];
rz(-1.9741156) q[3];
sx q[3];
rz(-1.4195331) q[3];
sx q[3];
rz(1.2742111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90299273) q[0];
sx q[0];
rz(-2.3165343) q[0];
sx q[0];
rz(2.9883244) q[0];
rz(-2.080147) q[1];
sx q[1];
rz(-0.57466424) q[1];
sx q[1];
rz(-1.9877888) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.823258) q[0];
sx q[0];
rz(-1.9360006) q[0];
sx q[0];
rz(1.2814786) q[0];
rz(-pi) q[1];
rz(-0.54406464) q[2];
sx q[2];
rz(-0.96368921) q[2];
sx q[2];
rz(-2.0008759) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.02281636) q[1];
sx q[1];
rz(-0.21636848) q[1];
sx q[1];
rz(-1.8323891) q[1];
rz(-pi) q[2];
rz(2.2581402) q[3];
sx q[3];
rz(-1.6646107) q[3];
sx q[3];
rz(-1.0132307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.29685059) q[2];
sx q[2];
rz(-1.9694318) q[2];
sx q[2];
rz(-1.5273757) q[2];
rz(-0.99689233) q[3];
sx q[3];
rz(-1.6485873) q[3];
sx q[3];
rz(-2.2629288) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8940354) q[0];
sx q[0];
rz(-2.0443125) q[0];
sx q[0];
rz(2.9515008) q[0];
rz(-0.67063531) q[1];
sx q[1];
rz(-1.9452483) q[1];
sx q[1];
rz(-1.6533096) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3538441) q[0];
sx q[0];
rz(-1.4453381) q[0];
sx q[0];
rz(1.5397443) q[0];
rz(-pi) q[1];
rz(0.52942099) q[2];
sx q[2];
rz(-0.501907) q[2];
sx q[2];
rz(1.5159964) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.30257672) q[1];
sx q[1];
rz(-2.1864428) q[1];
sx q[1];
rz(-2.2018196) q[1];
x q[2];
rz(0.87741239) q[3];
sx q[3];
rz(-1.1023695) q[3];
sx q[3];
rz(1.7978158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.0037447475) q[2];
sx q[2];
rz(-0.90128428) q[2];
sx q[2];
rz(-2.2424662) q[2];
rz(-0.51504618) q[3];
sx q[3];
rz(-2.6656272) q[3];
sx q[3];
rz(-0.17764828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0777733) q[0];
sx q[0];
rz(-1.3615006) q[0];
sx q[0];
rz(2.5831945) q[0];
rz(2.8521815) q[1];
sx q[1];
rz(-2.2402973) q[1];
sx q[1];
rz(-1.4351861) q[1];
rz(2.8171956) q[2];
sx q[2];
rz(-2.4051721) q[2];
sx q[2];
rz(0.13798513) q[2];
rz(-1.9990986) q[3];
sx q[3];
rz(-0.29853587) q[3];
sx q[3];
rz(1.5034911) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
