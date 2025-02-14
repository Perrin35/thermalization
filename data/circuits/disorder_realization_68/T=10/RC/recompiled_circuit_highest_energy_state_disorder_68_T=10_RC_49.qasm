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
rz(0.23213586) q[0];
sx q[0];
rz(2.870626) q[0];
sx q[0];
rz(10.597064) q[0];
rz(1.5341893) q[1];
sx q[1];
rz(-1.6521896) q[1];
sx q[1];
rz(0.54816562) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0294103) q[0];
sx q[0];
rz(-1.8444841) q[0];
sx q[0];
rz(3.0269483) q[0];
x q[1];
rz(1.0127958) q[2];
sx q[2];
rz(-1.1269953) q[2];
sx q[2];
rz(-1.9091801) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.89828789) q[1];
sx q[1];
rz(-0.60051781) q[1];
sx q[1];
rz(0.55796786) q[1];
x q[2];
rz(0.062730363) q[3];
sx q[3];
rz(-1.1219121) q[3];
sx q[3];
rz(-2.1470438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9008909) q[2];
sx q[2];
rz(-1.9343932) q[2];
sx q[2];
rz(-1.5473676) q[2];
rz(0.93260014) q[3];
sx q[3];
rz(-0.91351944) q[3];
sx q[3];
rz(-0.80476052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7050183) q[0];
sx q[0];
rz(-0.60282928) q[0];
sx q[0];
rz(0.21827179) q[0];
rz(1.5087347) q[1];
sx q[1];
rz(-2.2684596) q[1];
sx q[1];
rz(2.0192718) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86271226) q[0];
sx q[0];
rz(-0.76970184) q[0];
sx q[0];
rz(-0.93602009) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7709184) q[2];
sx q[2];
rz(-2.3014143) q[2];
sx q[2];
rz(0.49337988) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3716723) q[1];
sx q[1];
rz(-0.74410106) q[1];
sx q[1];
rz(-2.4456853) q[1];
rz(0.89790475) q[3];
sx q[3];
rz(-2.7460685) q[3];
sx q[3];
rz(-0.69098587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0108769) q[2];
sx q[2];
rz(-1.1900095) q[2];
sx q[2];
rz(-0.38401628) q[2];
rz(-2.5859517) q[3];
sx q[3];
rz(-0.58755392) q[3];
sx q[3];
rz(-2.2491992) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4246849) q[0];
sx q[0];
rz(-2.6838344) q[0];
sx q[0];
rz(-2.6440788) q[0];
rz(0.48490191) q[1];
sx q[1];
rz(-1.5496016) q[1];
sx q[1];
rz(-0.44152322) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2427485) q[0];
sx q[0];
rz(-1.4024807) q[0];
sx q[0];
rz(0.90710137) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1782568) q[2];
sx q[2];
rz(-0.52670331) q[2];
sx q[2];
rz(1.3046427) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.42558266) q[1];
sx q[1];
rz(-2.0507318) q[1];
sx q[1];
rz(0.77253567) q[1];
rz(-2.2763292) q[3];
sx q[3];
rz(-2.7332716) q[3];
sx q[3];
rz(-2.7029519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8170844) q[2];
sx q[2];
rz(-2.0980947) q[2];
sx q[2];
rz(0.95477742) q[2];
rz(0.020141007) q[3];
sx q[3];
rz(-2.3432799) q[3];
sx q[3];
rz(-0.32625833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7666053) q[0];
sx q[0];
rz(-2.6191481) q[0];
sx q[0];
rz(0.0047542714) q[0];
rz(-2.7004554) q[1];
sx q[1];
rz(-0.10103592) q[1];
sx q[1];
rz(-1.4099247) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4213181) q[0];
sx q[0];
rz(-1.4393115) q[0];
sx q[0];
rz(1.2340896) q[0];
rz(2.2994322) q[2];
sx q[2];
rz(-1.7787063) q[2];
sx q[2];
rz(2.4198857) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6769242) q[1];
sx q[1];
rz(-2.0629971) q[1];
sx q[1];
rz(-2.8850624) q[1];
rz(-pi) q[2];
rz(3.1143777) q[3];
sx q[3];
rz(-1.7593972) q[3];
sx q[3];
rz(2.1299429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8498174) q[2];
sx q[2];
rz(-1.8287649) q[2];
sx q[2];
rz(0.69668359) q[2];
rz(-1.5289395) q[3];
sx q[3];
rz(-0.63642514) q[3];
sx q[3];
rz(-2.4618885) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9686389) q[0];
sx q[0];
rz(-0.30690673) q[0];
sx q[0];
rz(0.90721834) q[0];
rz(-3.0061159) q[1];
sx q[1];
rz(-1.2716581) q[1];
sx q[1];
rz(0.69773331) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22963069) q[0];
sx q[0];
rz(-0.41423395) q[0];
sx q[0];
rz(1.8976865) q[0];
x q[1];
rz(-1.582007) q[2];
sx q[2];
rz(-0.19489842) q[2];
sx q[2];
rz(-0.29981183) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.88416022) q[1];
sx q[1];
rz(-1.8990277) q[1];
sx q[1];
rz(1.753367) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1118513) q[3];
sx q[3];
rz(-1.5248767) q[3];
sx q[3];
rz(2.8950402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.1873143) q[2];
sx q[2];
rz(-0.66437393) q[2];
sx q[2];
rz(0.3453671) q[2];
rz(0.49515381) q[3];
sx q[3];
rz(-2.5457355) q[3];
sx q[3];
rz(-0.12721795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0879764) q[0];
sx q[0];
rz(-1.6642267) q[0];
sx q[0];
rz(-1.9867058) q[0];
rz(-0.79479533) q[1];
sx q[1];
rz(-1.2966917) q[1];
sx q[1];
rz(-2.6577139) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4876955) q[0];
sx q[0];
rz(-1.9761818) q[0];
sx q[0];
rz(-0.8542819) q[0];
x q[1];
rz(-2.0675428) q[2];
sx q[2];
rz(-1.4306465) q[2];
sx q[2];
rz(-1.6983918) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.7495854) q[1];
sx q[1];
rz(-1.3806496) q[1];
sx q[1];
rz(-1.0984248) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.563298) q[3];
sx q[3];
rz(-1.9173754) q[3];
sx q[3];
rz(-2.9017096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.38857073) q[2];
sx q[2];
rz(-1.9066255) q[2];
sx q[2];
rz(1.1080326) q[2];
rz(-1.612251) q[3];
sx q[3];
rz(-0.11950167) q[3];
sx q[3];
rz(-0.39335462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18688467) q[0];
sx q[0];
rz(-1.0857546) q[0];
sx q[0];
rz(-2.4672274) q[0];
rz(-2.2652594) q[1];
sx q[1];
rz(-2.3969711) q[1];
sx q[1];
rz(3.1092627) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5139149) q[0];
sx q[0];
rz(-1.4490632) q[0];
sx q[0];
rz(-1.7072149) q[0];
rz(-pi) q[1];
rz(1.416308) q[2];
sx q[2];
rz(-1.7927899) q[2];
sx q[2];
rz(0.40419182) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2029973) q[1];
sx q[1];
rz(-1.055966) q[1];
sx q[1];
rz(1.1802303) q[1];
x q[2];
rz(2.9495839) q[3];
sx q[3];
rz(-2.1550094) q[3];
sx q[3];
rz(1.2835128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.18675599) q[2];
sx q[2];
rz(-0.43354553) q[2];
sx q[2];
rz(-1.0369066) q[2];
rz(1.2284651) q[3];
sx q[3];
rz(-2.2322673) q[3];
sx q[3];
rz(2.9508446) q[3];
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
rz(-2.4036338) q[0];
sx q[0];
rz(-3.1138595) q[0];
sx q[0];
rz(-0.24895689) q[0];
rz(-2.9320993) q[1];
sx q[1];
rz(-1.8003576) q[1];
sx q[1];
rz(2.478821) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4512206) q[0];
sx q[0];
rz(-1.5589542) q[0];
sx q[0];
rz(1.5607587) q[0];
x q[1];
rz(-0.49928645) q[2];
sx q[2];
rz(-2.3134276) q[2];
sx q[2];
rz(-1.7940831) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4718834) q[1];
sx q[1];
rz(-0.85195573) q[1];
sx q[1];
rz(-0.32841668) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2609188) q[3];
sx q[3];
rz(-2.204748) q[3];
sx q[3];
rz(3.0216963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2304307) q[2];
sx q[2];
rz(-1.5713567) q[2];
sx q[2];
rz(-2.7479808) q[2];
rz(0.39129928) q[3];
sx q[3];
rz(-0.53520447) q[3];
sx q[3];
rz(2.3894501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0126208) q[0];
sx q[0];
rz(-2.9550645) q[0];
sx q[0];
rz(-2.5616126) q[0];
rz(-0.36803666) q[1];
sx q[1];
rz(-0.8152222) q[1];
sx q[1];
rz(-0.11485242) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57844394) q[0];
sx q[0];
rz(-2.4029893) q[0];
sx q[0];
rz(-1.369801) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9909953) q[2];
sx q[2];
rz(-1.1335177) q[2];
sx q[2];
rz(0.33971805) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.65832135) q[1];
sx q[1];
rz(-1.8206925) q[1];
sx q[1];
rz(2.8037594) q[1];
x q[2];
rz(2.7389218) q[3];
sx q[3];
rz(-2.2244903) q[3];
sx q[3];
rz(1.634906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4174691) q[2];
sx q[2];
rz(-2.8026411) q[2];
sx q[2];
rz(-2.4453956) q[2];
rz(-2.0255069) q[3];
sx q[3];
rz(-2.3510272) q[3];
sx q[3];
rz(-2.4000786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29257947) q[0];
sx q[0];
rz(-0.96939033) q[0];
sx q[0];
rz(-2.868929) q[0];
rz(2.4478681) q[1];
sx q[1];
rz(-2.0246181) q[1];
sx q[1];
rz(-0.56347096) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7359223) q[0];
sx q[0];
rz(-1.5338681) q[0];
sx q[0];
rz(-1.4775299) q[0];
x q[1];
rz(2.9046802) q[2];
sx q[2];
rz(-1.713551) q[2];
sx q[2];
rz(-2.3430489) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.28501355) q[1];
sx q[1];
rz(-1.4895298) q[1];
sx q[1];
rz(1.2586428) q[1];
rz(-0.70448168) q[3];
sx q[3];
rz(-2.292013) q[3];
sx q[3];
rz(0.81126708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1571265) q[2];
sx q[2];
rz(-1.0447341) q[2];
sx q[2];
rz(-2.5926479) q[2];
rz(0.99556154) q[3];
sx q[3];
rz(-0.82058161) q[3];
sx q[3];
rz(2.5115749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7231049) q[0];
sx q[0];
rz(-0.99011078) q[0];
sx q[0];
rz(-0.44575442) q[0];
rz(-2.2743629) q[1];
sx q[1];
rz(-2.0384616) q[1];
sx q[1];
rz(-1.9834317) q[1];
rz(-1.340048) q[2];
sx q[2];
rz(-1.1044772) q[2];
sx q[2];
rz(0.92369631) q[2];
rz(1.7165202) q[3];
sx q[3];
rz(-1.0276262) q[3];
sx q[3];
rz(0.75158391) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
