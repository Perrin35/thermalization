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
rz(-2.8615992) q[0];
sx q[0];
rz(-2.1092829) q[0];
sx q[0];
rz(-1.9814459) q[0];
rz(-0.87814826) q[1];
sx q[1];
rz(-2.703009) q[1];
sx q[1];
rz(2.261472) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0297247) q[0];
sx q[0];
rz(-2.1639185) q[0];
sx q[0];
rz(0.2082227) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3841259) q[2];
sx q[2];
rz(-0.40928939) q[2];
sx q[2];
rz(3.0011506) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5118106) q[1];
sx q[1];
rz(-1.4079878) q[1];
sx q[1];
rz(2.1753499) q[1];
x q[2];
rz(1.6900469) q[3];
sx q[3];
rz(-1.1195983) q[3];
sx q[3];
rz(0.50332068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0264414) q[2];
sx q[2];
rz(-0.5618962) q[2];
sx q[2];
rz(2.0217516) q[2];
rz(-2.1378873) q[3];
sx q[3];
rz(-0.34990889) q[3];
sx q[3];
rz(-2.2899535) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.411946) q[0];
sx q[0];
rz(-2.248652) q[0];
sx q[0];
rz(-2.6935284) q[0];
rz(-1.0753151) q[1];
sx q[1];
rz(-2.0649464) q[1];
sx q[1];
rz(0.040464673) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69175289) q[0];
sx q[0];
rz(-1.7621855) q[0];
sx q[0];
rz(2.2474849) q[0];
rz(-3.05182) q[2];
sx q[2];
rz(-1.1589638) q[2];
sx q[2];
rz(1.8805945) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4869747) q[1];
sx q[1];
rz(-2.8697213) q[1];
sx q[1];
rz(-0.15403037) q[1];
rz(-0.31261698) q[3];
sx q[3];
rz(-1.5789869) q[3];
sx q[3];
rz(0.20144486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.86541092) q[2];
sx q[2];
rz(-2.2409596) q[2];
sx q[2];
rz(0.5033699) q[2];
rz(-1.6940176) q[3];
sx q[3];
rz(-1.2332799) q[3];
sx q[3];
rz(0.90731049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(2.4564948) q[0];
sx q[0];
rz(-2.1879897) q[0];
sx q[0];
rz(2.8431235) q[0];
rz(1.586277) q[1];
sx q[1];
rz(-1.7319873) q[1];
sx q[1];
rz(1.635199) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8776309) q[0];
sx q[0];
rz(-1.4957424) q[0];
sx q[0];
rz(2.8427426) q[0];
x q[1];
rz(-1.9048018) q[2];
sx q[2];
rz(-0.71454218) q[2];
sx q[2];
rz(1.8264019) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8019234) q[1];
sx q[1];
rz(-2.443772) q[1];
sx q[1];
rz(-0.20734792) q[1];
rz(-pi) q[2];
rz(1.8465145) q[3];
sx q[3];
rz(-0.19681588) q[3];
sx q[3];
rz(-2.4829602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7895268) q[2];
sx q[2];
rz(-0.76849285) q[2];
sx q[2];
rz(-2.3699769) q[2];
rz(-1.7302892) q[3];
sx q[3];
rz(-1.227042) q[3];
sx q[3];
rz(2.3119149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70374933) q[0];
sx q[0];
rz(-2.4898536) q[0];
sx q[0];
rz(-1.8072476) q[0];
rz(1.5127381) q[1];
sx q[1];
rz(-1.520243) q[1];
sx q[1];
rz(-0.13660647) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.893394) q[0];
sx q[0];
rz(-1.8780439) q[0];
sx q[0];
rz(2.224438) q[0];
rz(-pi) q[1];
rz(-3.0380546) q[2];
sx q[2];
rz(-1.8788726) q[2];
sx q[2];
rz(3.0335226) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1699511) q[1];
sx q[1];
rz(-1.882269) q[1];
sx q[1];
rz(2.4258652) q[1];
rz(1.9196687) q[3];
sx q[3];
rz(-1.8851981) q[3];
sx q[3];
rz(1.507198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.85285774) q[2];
sx q[2];
rz(-1.4790269) q[2];
sx q[2];
rz(2.9020201) q[2];
rz(-2.1824956) q[3];
sx q[3];
rz(-1.164091) q[3];
sx q[3];
rz(-1.0378999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94944823) q[0];
sx q[0];
rz(-1.0553772) q[0];
sx q[0];
rz(1.5437833) q[0];
rz(-1.1100618) q[1];
sx q[1];
rz(-1.9381356) q[1];
sx q[1];
rz(-1.69453) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8444203) q[0];
sx q[0];
rz(-1.8691952) q[0];
sx q[0];
rz(1.5633388) q[0];
rz(-pi) q[1];
rz(-1.9138811) q[2];
sx q[2];
rz(-1.6091847) q[2];
sx q[2];
rz(-2.6643201) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8548818) q[1];
sx q[1];
rz(-2.4457702) q[1];
sx q[1];
rz(-2.8665101) q[1];
rz(-pi) q[2];
rz(-1.9915699) q[3];
sx q[3];
rz(-2.07486) q[3];
sx q[3];
rz(2.1114388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1453104) q[2];
sx q[2];
rz(-2.0029533) q[2];
sx q[2];
rz(0.81926695) q[2];
rz(-0.94775689) q[3];
sx q[3];
rz(-2.3534333) q[3];
sx q[3];
rz(1.88571) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6325697) q[0];
sx q[0];
rz(-3.0036354) q[0];
sx q[0];
rz(-2.5914958) q[0];
rz(-2.8944648) q[1];
sx q[1];
rz(-1.988966) q[1];
sx q[1];
rz(2.1955042) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9086015) q[0];
sx q[0];
rz(-2.6938022) q[0];
sx q[0];
rz(0.93771387) q[0];
rz(-1.7839512) q[2];
sx q[2];
rz(-1.6317131) q[2];
sx q[2];
rz(-2.3062381) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4934686) q[1];
sx q[1];
rz(-1.6065242) q[1];
sx q[1];
rz(-0.46983998) q[1];
rz(-pi) q[2];
x q[2];
rz(0.40762679) q[3];
sx q[3];
rz(-0.25926155) q[3];
sx q[3];
rz(1.892499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3863824) q[2];
sx q[2];
rz(-2.5911665) q[2];
sx q[2];
rz(-0.68106252) q[2];
rz(2.3303473) q[3];
sx q[3];
rz(-1.0439876) q[3];
sx q[3];
rz(-2.5913141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.022245971) q[0];
sx q[0];
rz(-2.5211054) q[0];
sx q[0];
rz(-0.74561179) q[0];
rz(1.9781808) q[1];
sx q[1];
rz(-2.5201576) q[1];
sx q[1];
rz(2.9643639) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1354373) q[0];
sx q[0];
rz(-2.7072142) q[0];
sx q[0];
rz(0.88157369) q[0];
x q[1];
rz(2.178995) q[2];
sx q[2];
rz(-3.0311435) q[2];
sx q[2];
rz(2.3640574) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9504576) q[1];
sx q[1];
rz(-2.0759575) q[1];
sx q[1];
rz(2.7209366) q[1];
rz(-0.050595119) q[3];
sx q[3];
rz(-2.4462163) q[3];
sx q[3];
rz(-0.68948244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5916198) q[2];
sx q[2];
rz(-2.1592996) q[2];
sx q[2];
rz(1.4493235) q[2];
rz(-2.5189404) q[3];
sx q[3];
rz(-2.6144274) q[3];
sx q[3];
rz(1.8221633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1249238) q[0];
sx q[0];
rz(-1.387384) q[0];
sx q[0];
rz(-1.5119875) q[0];
rz(-2.2177057) q[1];
sx q[1];
rz(-1.4199363) q[1];
sx q[1];
rz(-0.61663827) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9074207) q[0];
sx q[0];
rz(-1.7230095) q[0];
sx q[0];
rz(-2.175075) q[0];
rz(1.6144912) q[2];
sx q[2];
rz(-2.0382883) q[2];
sx q[2];
rz(-1.0861008) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5376327) q[1];
sx q[1];
rz(-0.62509552) q[1];
sx q[1];
rz(-2.9949466) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.736896) q[3];
sx q[3];
rz(-2.5202978) q[3];
sx q[3];
rz(-2.9018847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3672678) q[2];
sx q[2];
rz(-0.35949817) q[2];
sx q[2];
rz(-0.16505879) q[2];
rz(0.74602357) q[3];
sx q[3];
rz(-1.5187998) q[3];
sx q[3];
rz(3.0987958) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5999254) q[0];
sx q[0];
rz(-2.8494819) q[0];
sx q[0];
rz(0.41411972) q[0];
rz(-2.7060624) q[1];
sx q[1];
rz(-2.6256517) q[1];
sx q[1];
rz(-1.863265) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.965344) q[0];
sx q[0];
rz(-2.8581736) q[0];
sx q[0];
rz(2.118628) q[0];
rz(-pi) q[1];
rz(1.7444853) q[2];
sx q[2];
rz(-1.2635482) q[2];
sx q[2];
rz(-1.1905313) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8929249) q[1];
sx q[1];
rz(-1.5589412) q[1];
sx q[1];
rz(-0.067889386) q[1];
rz(-2.0986365) q[3];
sx q[3];
rz(-1.5164638) q[3];
sx q[3];
rz(-2.899037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4182959) q[2];
sx q[2];
rz(-1.7593242) q[2];
sx q[2];
rz(-0.74271512) q[2];
rz(0.79139477) q[3];
sx q[3];
rz(-0.68358889) q[3];
sx q[3];
rz(1.4092103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.047569711) q[0];
sx q[0];
rz(-0.3293193) q[0];
sx q[0];
rz(2.2299715) q[0];
rz(-0.98512638) q[1];
sx q[1];
rz(-2.7148425) q[1];
sx q[1];
rz(-1.8034579) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3365173) q[0];
sx q[0];
rz(-2.0455845) q[0];
sx q[0];
rz(-2.7874448) q[0];
rz(-pi) q[1];
rz(0.88078946) q[2];
sx q[2];
rz(-2.6652209) q[2];
sx q[2];
rz(-0.42467372) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.12305189) q[1];
sx q[1];
rz(-1.9705074) q[1];
sx q[1];
rz(-2.2247061) q[1];
rz(-pi) q[2];
x q[2];
rz(0.84566859) q[3];
sx q[3];
rz(-2.659043) q[3];
sx q[3];
rz(-0.023511271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.35475874) q[2];
sx q[2];
rz(-1.4749196) q[2];
sx q[2];
rz(1.1653398) q[2];
rz(-1.745584) q[3];
sx q[3];
rz(-1.952012) q[3];
sx q[3];
rz(1.5918119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5042481) q[0];
sx q[0];
rz(-1.5865542) q[0];
sx q[0];
rz(-2.4028548) q[0];
rz(-0.63738102) q[1];
sx q[1];
rz(-1.5370054) q[1];
sx q[1];
rz(-0.76269033) q[1];
rz(-1.1470096) q[2];
sx q[2];
rz(-1.7041361) q[2];
sx q[2];
rz(2.3455483) q[2];
rz(2.021234) q[3];
sx q[3];
rz(-1.0584581) q[3];
sx q[3];
rz(-1.2367482) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
