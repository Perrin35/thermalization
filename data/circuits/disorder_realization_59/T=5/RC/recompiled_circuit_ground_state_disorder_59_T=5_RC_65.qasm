OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.6049603) q[0];
sx q[0];
rz(3.3024238) q[0];
sx q[0];
rz(9.816058) q[0];
rz(3.0985576) q[1];
sx q[1];
rz(-1.6208836) q[1];
sx q[1];
rz(-2.8032141) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3489549) q[0];
sx q[0];
rz(-0.39039877) q[0];
sx q[0];
rz(-0.81116311) q[0];
rz(1.7853043) q[2];
sx q[2];
rz(-1.5212998) q[2];
sx q[2];
rz(1.6550198) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1115626) q[1];
sx q[1];
rz(-1.3805559) q[1];
sx q[1];
rz(-2.2493786) q[1];
rz(-pi) q[2];
rz(2.8998221) q[3];
sx q[3];
rz(-0.42754506) q[3];
sx q[3];
rz(1.1113785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.36838621) q[2];
sx q[2];
rz(-2.1619004) q[2];
sx q[2];
rz(-2.0835908) q[2];
rz(3.0392785) q[3];
sx q[3];
rz(-1.6175858) q[3];
sx q[3];
rz(1.7412831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31829396) q[0];
sx q[0];
rz(-2.3866391) q[0];
sx q[0];
rz(-2.9276983) q[0];
rz(1.3338044) q[1];
sx q[1];
rz(-0.50024453) q[1];
sx q[1];
rz(0.28876367) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5014415) q[0];
sx q[0];
rz(-2.1955865) q[0];
sx q[0];
rz(-2.9286317) q[0];
x q[1];
rz(0.83728055) q[2];
sx q[2];
rz(-2.0990666) q[2];
sx q[2];
rz(-2.3579896) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.97280661) q[1];
sx q[1];
rz(-1.791535) q[1];
sx q[1];
rz(-0.03038637) q[1];
rz(-pi) q[2];
rz(1.2330194) q[3];
sx q[3];
rz(-1.1572517) q[3];
sx q[3];
rz(0.097584399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.6194965) q[2];
sx q[2];
rz(-0.82021004) q[2];
sx q[2];
rz(1.8572469) q[2];
rz(-1.8340825) q[3];
sx q[3];
rz(-1.8239572) q[3];
sx q[3];
rz(-2.4431084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4035325) q[0];
sx q[0];
rz(-1.0665251) q[0];
sx q[0];
rz(2.2962978) q[0];
rz(-0.90557939) q[1];
sx q[1];
rz(-2.5366492) q[1];
sx q[1];
rz(1.9639429) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3628927) q[0];
sx q[0];
rz(-1.4757627) q[0];
sx q[0];
rz(1.5508077) q[0];
x q[1];
rz(1.0217754) q[2];
sx q[2];
rz(-2.0284101) q[2];
sx q[2];
rz(2.3624376) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3556087) q[1];
sx q[1];
rz(-2.5289383) q[1];
sx q[1];
rz(2.775928) q[1];
rz(-2.2159141) q[3];
sx q[3];
rz(-1.4429907) q[3];
sx q[3];
rz(-2.5296581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5251069) q[2];
sx q[2];
rz(-2.1973886) q[2];
sx q[2];
rz(-3.0109829) q[2];
rz(-1.3579926) q[3];
sx q[3];
rz(-1.0087174) q[3];
sx q[3];
rz(1.85359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4865049) q[0];
sx q[0];
rz(-2.8772652) q[0];
sx q[0];
rz(2.7274912) q[0];
rz(2.2762903) q[1];
sx q[1];
rz(-1.2369786) q[1];
sx q[1];
rz(2.5383331) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7287528) q[0];
sx q[0];
rz(-1.9590634) q[0];
sx q[0];
rz(-0.30465841) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9500721) q[2];
sx q[2];
rz(-0.3580123) q[2];
sx q[2];
rz(-0.70068089) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.016876) q[1];
sx q[1];
rz(-1.6137329) q[1];
sx q[1];
rz(0.34833585) q[1];
rz(-pi) q[2];
x q[2];
rz(0.11988414) q[3];
sx q[3];
rz(-1.2601687) q[3];
sx q[3];
rz(-1.436123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6174751) q[2];
sx q[2];
rz(-1.5323428) q[2];
sx q[2];
rz(-2.4821846) q[2];
rz(-0.13255969) q[3];
sx q[3];
rz(-2.5949251) q[3];
sx q[3];
rz(-0.70422712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34683126) q[0];
sx q[0];
rz(-0.95585388) q[0];
sx q[0];
rz(-0.26926789) q[0];
rz(1.6412093) q[1];
sx q[1];
rz(-1.2790479) q[1];
sx q[1];
rz(1.0795116) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37845072) q[0];
sx q[0];
rz(-2.2104163) q[0];
sx q[0];
rz(-1.3467953) q[0];
rz(-pi) q[1];
rz(0.66840489) q[2];
sx q[2];
rz(-0.89819095) q[2];
sx q[2];
rz(2.2318537) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.36240807) q[1];
sx q[1];
rz(-2.0498865) q[1];
sx q[1];
rz(-0.37491231) q[1];
x q[2];
rz(0.57860878) q[3];
sx q[3];
rz(-0.37453412) q[3];
sx q[3];
rz(0.43363627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.72011224) q[2];
sx q[2];
rz(-2.5613027) q[2];
sx q[2];
rz(-0.87023467) q[2];
rz(-1.3358215) q[3];
sx q[3];
rz(-1.8060874) q[3];
sx q[3];
rz(0.67572063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1338335) q[0];
sx q[0];
rz(-0.02820153) q[0];
sx q[0];
rz(2.5303685) q[0];
rz(2.4080343) q[1];
sx q[1];
rz(-1.6707784) q[1];
sx q[1];
rz(-2.7582817) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54827842) q[0];
sx q[0];
rz(-0.76799521) q[0];
sx q[0];
rz(2.7357714) q[0];
rz(1.352241) q[2];
sx q[2];
rz(-2.3862184) q[2];
sx q[2];
rz(-2.4807376) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.53364175) q[1];
sx q[1];
rz(-1.4429394) q[1];
sx q[1];
rz(-0.92725384) q[1];
rz(-pi) q[2];
rz(1.1556968) q[3];
sx q[3];
rz(-2.2533247) q[3];
sx q[3];
rz(1.1235588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.40206259) q[2];
sx q[2];
rz(-0.5548839) q[2];
sx q[2];
rz(-2.128672) q[2];
rz(2.5439475) q[3];
sx q[3];
rz(-1.7326071) q[3];
sx q[3];
rz(0.92596084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0830072) q[0];
sx q[0];
rz(-0.42917955) q[0];
sx q[0];
rz(0.035932628) q[0];
rz(1.9195456) q[1];
sx q[1];
rz(-2.4655894) q[1];
sx q[1];
rz(2.8935208) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67001024) q[0];
sx q[0];
rz(-1.8267794) q[0];
sx q[0];
rz(2.94245) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8808458) q[2];
sx q[2];
rz(-2.6364821) q[2];
sx q[2];
rz(1.3253761) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9315459) q[1];
sx q[1];
rz(-1.4074874) q[1];
sx q[1];
rz(1.1985267) q[1];
x q[2];
rz(0.024240243) q[3];
sx q[3];
rz(-3.0366237) q[3];
sx q[3];
rz(-1.6258383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0566569) q[2];
sx q[2];
rz(-2.0818905) q[2];
sx q[2];
rz(-2.6094931) q[2];
rz(2.687124) q[3];
sx q[3];
rz(-1.5957417) q[3];
sx q[3];
rz(1.2940297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71659511) q[0];
sx q[0];
rz(-1.0497365) q[0];
sx q[0];
rz(0.24018921) q[0];
rz(-2.5364618) q[1];
sx q[1];
rz(-1.3477707) q[1];
sx q[1];
rz(-0.64632195) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8511348) q[0];
sx q[0];
rz(-2.2735828) q[0];
sx q[0];
rz(-2.7443307) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.209561) q[2];
sx q[2];
rz(-1.5363599) q[2];
sx q[2];
rz(-3.051935) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5390215) q[1];
sx q[1];
rz(-2.638431) q[1];
sx q[1];
rz(2.7064067) q[1];
rz(-0.21625285) q[3];
sx q[3];
rz(-0.92666221) q[3];
sx q[3];
rz(0.82051051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1120844) q[2];
sx q[2];
rz(-1.6479475) q[2];
sx q[2];
rz(-3.0599111) q[2];
rz(-2.2937842) q[3];
sx q[3];
rz(-0.40926465) q[3];
sx q[3];
rz(-0.22296396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34143701) q[0];
sx q[0];
rz(-0.84469047) q[0];
sx q[0];
rz(-2.2499114) q[0];
rz(1.910123) q[1];
sx q[1];
rz(-2.6530118) q[1];
sx q[1];
rz(-0.60985342) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9187136) q[0];
sx q[0];
rz(-1.0787258) q[0];
sx q[0];
rz(2.6124873) q[0];
rz(-pi) q[1];
rz(-3.1270486) q[2];
sx q[2];
rz(-0.78476465) q[2];
sx q[2];
rz(2.5958928) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5130723) q[1];
sx q[1];
rz(-1.2781464) q[1];
sx q[1];
rz(0.77192941) q[1];
rz(-pi) q[2];
rz(1.9159007) q[3];
sx q[3];
rz(-1.5550348) q[3];
sx q[3];
rz(-1.0712766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2764728) q[2];
sx q[2];
rz(-0.96729326) q[2];
sx q[2];
rz(-2.8709732) q[2];
rz(-2.596415) q[3];
sx q[3];
rz(-1.3278278) q[3];
sx q[3];
rz(0.24817185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4973064) q[0];
sx q[0];
rz(-1.3077284) q[0];
sx q[0];
rz(2.9551031) q[0];
rz(-1.2093774) q[1];
sx q[1];
rz(-1.3554363) q[1];
sx q[1];
rz(3.1219416) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6015897) q[0];
sx q[0];
rz(-2.181291) q[0];
sx q[0];
rz(-2.6363346) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8397464) q[2];
sx q[2];
rz(-1.8230613) q[2];
sx q[2];
rz(2.4076574) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.84148592) q[1];
sx q[1];
rz(-1.4177313) q[1];
sx q[1];
rz(0.35565214) q[1];
rz(-2.4830706) q[3];
sx q[3];
rz(-0.88044518) q[3];
sx q[3];
rz(0.041680574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.43202117) q[2];
sx q[2];
rz(-0.95546466) q[2];
sx q[2];
rz(-1.552399) q[2];
rz(0.10913695) q[3];
sx q[3];
rz(-1.2723294) q[3];
sx q[3];
rz(-2.8652338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5149665) q[0];
sx q[0];
rz(-1.4780541) q[0];
sx q[0];
rz(1.9274101) q[0];
rz(-1.7264438) q[1];
sx q[1];
rz(-2.1198004) q[1];
sx q[1];
rz(-1.459495) q[1];
rz(-1.8564275) q[2];
sx q[2];
rz(-2.0002596) q[2];
sx q[2];
rz(3.132435) q[2];
rz(-1.6676025) q[3];
sx q[3];
rz(-1.9514241) q[3];
sx q[3];
rz(-1.5497691) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
