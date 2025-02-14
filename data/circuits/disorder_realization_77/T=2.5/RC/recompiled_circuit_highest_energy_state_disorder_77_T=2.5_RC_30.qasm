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
rz(-0.4975118) q[0];
sx q[0];
rz(4.4805718) q[0];
sx q[0];
rz(6.5988402) q[0];
rz(1.1701801) q[1];
sx q[1];
rz(2.7538731) q[1];
sx q[1];
rz(8.8207689) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16306388) q[0];
sx q[0];
rz(-1.2685568) q[0];
sx q[0];
rz(0.053044293) q[0];
rz(2.6409615) q[2];
sx q[2];
rz(-2.1785469) q[2];
sx q[2];
rz(-1.0333824) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.042272557) q[1];
sx q[1];
rz(-1.0612773) q[1];
sx q[1];
rz(-0.38739844) q[1];
rz(-pi) q[2];
rz(1.2420869) q[3];
sx q[3];
rz(-2.8162873) q[3];
sx q[3];
rz(-1.6360375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.151256) q[2];
sx q[2];
rz(-1.9467111) q[2];
sx q[2];
rz(-1.3585496) q[2];
rz(-1.5938866) q[3];
sx q[3];
rz(-2.2113776) q[3];
sx q[3];
rz(-1.5096629) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5556521) q[0];
sx q[0];
rz(-2.5037615) q[0];
sx q[0];
rz(-1.1974539) q[0];
rz(-2.6400631) q[1];
sx q[1];
rz(-1.5634368) q[1];
sx q[1];
rz(-1.4362358) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89051188) q[0];
sx q[0];
rz(-2.5174826) q[0];
sx q[0];
rz(-1.1680383) q[0];
rz(-pi) q[1];
rz(2.4386232) q[2];
sx q[2];
rz(-1.4948934) q[2];
sx q[2];
rz(-0.10026201) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9402071) q[1];
sx q[1];
rz(-1.4680982) q[1];
sx q[1];
rz(0.19725712) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6858176) q[3];
sx q[3];
rz(-1.3823798) q[3];
sx q[3];
rz(0.51793232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2531835) q[2];
sx q[2];
rz(-1.9056355) q[2];
sx q[2];
rz(-3.1186228) q[2];
rz(0.90211558) q[3];
sx q[3];
rz(-1.2227367) q[3];
sx q[3];
rz(0.42660108) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4334634) q[0];
sx q[0];
rz(-1.7925649) q[0];
sx q[0];
rz(-0.9915114) q[0];
rz(2.1082711) q[1];
sx q[1];
rz(-0.28305498) q[1];
sx q[1];
rz(1.2352157) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9571494) q[0];
sx q[0];
rz(-2.1547124) q[0];
sx q[0];
rz(-1.8690442) q[0];
rz(-pi) q[1];
rz(2.5318145) q[2];
sx q[2];
rz(-1.9710566) q[2];
sx q[2];
rz(-0.61851172) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8950303) q[1];
sx q[1];
rz(-0.8324648) q[1];
sx q[1];
rz(-2.8042996) q[1];
x q[2];
rz(-2.3091812) q[3];
sx q[3];
rz(-2.4042272) q[3];
sx q[3];
rz(-0.6399782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.027355) q[2];
sx q[2];
rz(-0.75216746) q[2];
sx q[2];
rz(2.1503964) q[2];
rz(1.2973971) q[3];
sx q[3];
rz(-1.2605366) q[3];
sx q[3];
rz(-1.7245002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7436413) q[0];
sx q[0];
rz(-0.93248168) q[0];
sx q[0];
rz(2.3241296) q[0];
rz(0.3262597) q[1];
sx q[1];
rz(-1.1810818) q[1];
sx q[1];
rz(-2.356333) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1251251) q[0];
sx q[0];
rz(-2.7332691) q[0];
sx q[0];
rz(-0.84618469) q[0];
rz(-2.4982959) q[2];
sx q[2];
rz(-1.729082) q[2];
sx q[2];
rz(0.23223755) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7689159) q[1];
sx q[1];
rz(-2.0119889) q[1];
sx q[1];
rz(2.6720474) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6503851) q[3];
sx q[3];
rz(-1.2168988) q[3];
sx q[3];
rz(0.45512629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0205445) q[2];
sx q[2];
rz(-0.64457568) q[2];
sx q[2];
rz(-2.5784967) q[2];
rz(-0.8935039) q[3];
sx q[3];
rz(-1.1734791) q[3];
sx q[3];
rz(-1.9068498) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9017482) q[0];
sx q[0];
rz(-0.32308602) q[0];
sx q[0];
rz(-1.595994) q[0];
rz(1.0218989) q[1];
sx q[1];
rz(-2.0183759) q[1];
sx q[1];
rz(-2.9603069) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5820532) q[0];
sx q[0];
rz(-1.3783558) q[0];
sx q[0];
rz(2.0700702) q[0];
x q[1];
rz(-2.6659545) q[2];
sx q[2];
rz(-1.6132832) q[2];
sx q[2];
rz(-1.795639) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5904894) q[1];
sx q[1];
rz(-0.22399513) q[1];
sx q[1];
rz(1.444054) q[1];
x q[2];
rz(2.8195404) q[3];
sx q[3];
rz(-0.46370927) q[3];
sx q[3];
rz(2.2820594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1229317) q[2];
sx q[2];
rz(-0.67492008) q[2];
sx q[2];
rz(-1.2459416) q[2];
rz(0.59257007) q[3];
sx q[3];
rz(-1.4453245) q[3];
sx q[3];
rz(0.84111253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5331921) q[0];
sx q[0];
rz(-2.1493981) q[0];
sx q[0];
rz(-2.8327508) q[0];
rz(2.8187075) q[1];
sx q[1];
rz(-2.7643118) q[1];
sx q[1];
rz(3.0016532) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2353781) q[0];
sx q[0];
rz(-1.8188634) q[0];
sx q[0];
rz(0.95309044) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1148648) q[2];
sx q[2];
rz(-1.4856292) q[2];
sx q[2];
rz(2.2625201) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4637101) q[1];
sx q[1];
rz(-1.2695128) q[1];
sx q[1];
rz(0.82656411) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2298584) q[3];
sx q[3];
rz(-1.3305802) q[3];
sx q[3];
rz(1.5347337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.66190019) q[2];
sx q[2];
rz(-2.4883344) q[2];
sx q[2];
rz(2.8996186) q[2];
rz(2.3954929) q[3];
sx q[3];
rz(-1.7753764) q[3];
sx q[3];
rz(-1.7297176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.022843) q[0];
sx q[0];
rz(-1.043909) q[0];
sx q[0];
rz(-1.8705077) q[0];
rz(-0.56308833) q[1];
sx q[1];
rz(-0.92910281) q[1];
sx q[1];
rz(1.7005327) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8898687) q[0];
sx q[0];
rz(-1.8419203) q[0];
sx q[0];
rz(2.1878178) q[0];
x q[1];
rz(1.7092429) q[2];
sx q[2];
rz(-1.1927774) q[2];
sx q[2];
rz(-2.0496617) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.064441) q[1];
sx q[1];
rz(-1.7930601) q[1];
sx q[1];
rz(-0.42776107) q[1];
x q[2];
rz(0.7251803) q[3];
sx q[3];
rz(-1.4075052) q[3];
sx q[3];
rz(-2.0189328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.21961221) q[2];
sx q[2];
rz(-1.2954804) q[2];
sx q[2];
rz(1.8249576) q[2];
rz(2.5804139) q[3];
sx q[3];
rz(-2.1771274) q[3];
sx q[3];
rz(-1.6130028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0372666) q[0];
sx q[0];
rz(-2.9236561) q[0];
sx q[0];
rz(0.88687801) q[0];
rz(0.2001702) q[1];
sx q[1];
rz(-1.472241) q[1];
sx q[1];
rz(1.0221457) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5258785) q[0];
sx q[0];
rz(-1.4396884) q[0];
sx q[0];
rz(0.025630533) q[0];
rz(-pi) q[1];
rz(2.4046477) q[2];
sx q[2];
rz(-0.96074694) q[2];
sx q[2];
rz(1.7231154) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8238942) q[1];
sx q[1];
rz(-2.3168457) q[1];
sx q[1];
rz(0.0024248799) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6273055) q[3];
sx q[3];
rz(-0.64046851) q[3];
sx q[3];
rz(-2.2286704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.48978051) q[2];
sx q[2];
rz(-0.45280364) q[2];
sx q[2];
rz(-2.32302) q[2];
rz(2.1583648) q[3];
sx q[3];
rz(-0.95663095) q[3];
sx q[3];
rz(0.53808588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.035456903) q[0];
sx q[0];
rz(-1.1470733) q[0];
sx q[0];
rz(-1.5552833) q[0];
rz(-0.69681329) q[1];
sx q[1];
rz(-0.15334829) q[1];
sx q[1];
rz(-2.3568025) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0138088) q[0];
sx q[0];
rz(-1.133636) q[0];
sx q[0];
rz(-2.5378835) q[0];
rz(-pi) q[1];
rz(-0.92918877) q[2];
sx q[2];
rz(-2.8904473) q[2];
sx q[2];
rz(-0.26640688) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.49336772) q[1];
sx q[1];
rz(-2.7052212) q[1];
sx q[1];
rz(-3.0619951) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1510714) q[3];
sx q[3];
rz(-1.0015206) q[3];
sx q[3];
rz(1.1476393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.67184225) q[2];
sx q[2];
rz(-1.8639114) q[2];
sx q[2];
rz(2.3308241) q[2];
rz(-1.9226711) q[3];
sx q[3];
rz(-2.6007077) q[3];
sx q[3];
rz(-1.0651275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3588381) q[0];
sx q[0];
rz(-2.8420119) q[0];
sx q[0];
rz(-1.4240356) q[0];
rz(2.3639823) q[1];
sx q[1];
rz(-0.97336665) q[1];
sx q[1];
rz(-2.3861859) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1574673) q[0];
sx q[0];
rz(-1.3608772) q[0];
sx q[0];
rz(-0.0046878417) q[0];
x q[1];
rz(3.099775) q[2];
sx q[2];
rz(-2.7448069) q[2];
sx q[2];
rz(2.7355742) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3302119) q[1];
sx q[1];
rz(-0.63618681) q[1];
sx q[1];
rz(-2.5274171) q[1];
rz(1.8978664) q[3];
sx q[3];
rz(-2.2199515) q[3];
sx q[3];
rz(2.7130733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5648254) q[2];
sx q[2];
rz(-2.8456523) q[2];
sx q[2];
rz(0.15288615) q[2];
rz(1.2711924) q[3];
sx q[3];
rz(-1.252424) q[3];
sx q[3];
rz(-2.474474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85970238) q[0];
sx q[0];
rz(-0.49260456) q[0];
sx q[0];
rz(2.3717666) q[0];
rz(1.9181171) q[1];
sx q[1];
rz(-1.2212831) q[1];
sx q[1];
rz(2.4881359) q[1];
rz(2.5443947) q[2];
sx q[2];
rz(-1.2033249) q[2];
sx q[2];
rz(-2.0133599) q[2];
rz(-2.7889403) q[3];
sx q[3];
rz(-1.1897539) q[3];
sx q[3];
rz(-2.6170058) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
