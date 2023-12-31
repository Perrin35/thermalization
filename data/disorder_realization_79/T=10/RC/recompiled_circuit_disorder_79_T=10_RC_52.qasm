OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.22566158) q[0];
sx q[0];
rz(-2.2731279) q[0];
sx q[0];
rz(-2.948569) q[0];
rz(1.141619) q[1];
sx q[1];
rz(-0.42998278) q[1];
sx q[1];
rz(2.4584682) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8773168) q[0];
sx q[0];
rz(-1.5125934) q[0];
sx q[0];
rz(0.067406128) q[0];
x q[1];
rz(-2.9002951) q[2];
sx q[2];
rz(-1.9120875) q[2];
sx q[2];
rz(-1.2933921) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.623917) q[1];
sx q[1];
rz(-0.78849925) q[1];
sx q[1];
rz(2.4967525) q[1];
x q[2];
rz(-2.1360374) q[3];
sx q[3];
rz(-1.2473277) q[3];
sx q[3];
rz(-2.2068057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1258939) q[2];
sx q[2];
rz(-1.761972) q[2];
sx q[2];
rz(-2.0430298) q[2];
rz(2.0627608) q[3];
sx q[3];
rz(-0.94509411) q[3];
sx q[3];
rz(-1.1014972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1612448) q[0];
sx q[0];
rz(-0.315936) q[0];
sx q[0];
rz(2.9336477) q[0];
rz(-2.5646599) q[1];
sx q[1];
rz(-0.88795841) q[1];
sx q[1];
rz(1.6764486) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4629048) q[0];
sx q[0];
rz(-1.573274) q[0];
sx q[0];
rz(0.72420995) q[0];
x q[1];
rz(-2.0589774) q[2];
sx q[2];
rz(-1.7980051) q[2];
sx q[2];
rz(0.28085923) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.84491731) q[1];
sx q[1];
rz(-2.1293318) q[1];
sx q[1];
rz(-2.1293473) q[1];
x q[2];
rz(0.37104718) q[3];
sx q[3];
rz(-0.78073946) q[3];
sx q[3];
rz(2.3649529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1229822) q[2];
sx q[2];
rz(-0.98615042) q[2];
sx q[2];
rz(-3.0318276) q[2];
rz(2.5189853) q[3];
sx q[3];
rz(-0.37125769) q[3];
sx q[3];
rz(1.4573147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1784172) q[0];
sx q[0];
rz(-0.85997471) q[0];
sx q[0];
rz(-0.51613581) q[0];
rz(-0.57488817) q[1];
sx q[1];
rz(-2.2153885) q[1];
sx q[1];
rz(2.3410472) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8974509) q[0];
sx q[0];
rz(-1.4276917) q[0];
sx q[0];
rz(-1.1600526) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.66652253) q[2];
sx q[2];
rz(-2.1608673) q[2];
sx q[2];
rz(0.18596622) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2785981) q[1];
sx q[1];
rz(-1.9837712) q[1];
sx q[1];
rz(-1.5763271) q[1];
x q[2];
rz(1.21739) q[3];
sx q[3];
rz(-2.4089775) q[3];
sx q[3];
rz(-2.4917847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4642554) q[2];
sx q[2];
rz(-2.8223473) q[2];
sx q[2];
rz(1.3595954) q[2];
rz(-2.1740186) q[3];
sx q[3];
rz(-1.8680957) q[3];
sx q[3];
rz(-1.7165855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1490705) q[0];
sx q[0];
rz(-1.8795805) q[0];
sx q[0];
rz(2.676679) q[0];
rz(-0.34856302) q[1];
sx q[1];
rz(-0.26270738) q[1];
sx q[1];
rz(-1.0850614) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1151428) q[0];
sx q[0];
rz(-1.0466252) q[0];
sx q[0];
rz(1.3036222) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.24809804) q[2];
sx q[2];
rz(-1.9271701) q[2];
sx q[2];
rz(2.6756289) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3241868) q[1];
sx q[1];
rz(-1.7122014) q[1];
sx q[1];
rz(-1.6721339) q[1];
x q[2];
rz(-0.55235858) q[3];
sx q[3];
rz(-2.4529152) q[3];
sx q[3];
rz(0.83113447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.00502914) q[2];
sx q[2];
rz(-2.0932784) q[2];
sx q[2];
rz(-2.4604649) q[2];
rz(2.629771) q[3];
sx q[3];
rz(-0.32326439) q[3];
sx q[3];
rz(0.26369035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27914771) q[0];
sx q[0];
rz(-2.6036766) q[0];
sx q[0];
rz(-1.7329247) q[0];
rz(0.43235835) q[1];
sx q[1];
rz(-2.2996348) q[1];
sx q[1];
rz(2.1599105) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78445804) q[0];
sx q[0];
rz(-1.1103837) q[0];
sx q[0];
rz(-0.61607342) q[0];
rz(1.328674) q[2];
sx q[2];
rz(-1.7244581) q[2];
sx q[2];
rz(0.15649934) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2770734) q[1];
sx q[1];
rz(-2.652087) q[1];
sx q[1];
rz(0.31788748) q[1];
x q[2];
rz(0.97335191) q[3];
sx q[3];
rz(-1.5095599) q[3];
sx q[3];
rz(0.51613584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0354054) q[2];
sx q[2];
rz(-1.4865439) q[2];
sx q[2];
rz(-2.6521818) q[2];
rz(-1.0148467) q[3];
sx q[3];
rz(-0.5265407) q[3];
sx q[3];
rz(-1.3283407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6569825) q[0];
sx q[0];
rz(-2.9841612) q[0];
sx q[0];
rz(2.7691675) q[0];
rz(-1.8107481) q[1];
sx q[1];
rz(-1.0354038) q[1];
sx q[1];
rz(2.9763124) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80711354) q[0];
sx q[0];
rz(-1.665859) q[0];
sx q[0];
rz(-1.6001742) q[0];
rz(-2.4539102) q[2];
sx q[2];
rz(-1.624794) q[2];
sx q[2];
rz(2.8962367) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2416934) q[1];
sx q[1];
rz(-2.4447933) q[1];
sx q[1];
rz(0.65710575) q[1];
x q[2];
rz(-1.1092471) q[3];
sx q[3];
rz(-2.647532) q[3];
sx q[3];
rz(2.9500614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.91339397) q[2];
sx q[2];
rz(-1.0281576) q[2];
sx q[2];
rz(-1.4432663) q[2];
rz(2.7741487) q[3];
sx q[3];
rz(-1.2747217) q[3];
sx q[3];
rz(-2.929556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(1.3570324) q[0];
sx q[0];
rz(-1.0914047) q[0];
sx q[0];
rz(-0.34926397) q[0];
rz(0.7473942) q[1];
sx q[1];
rz(-0.29577297) q[1];
sx q[1];
rz(0.73648891) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1897141) q[0];
sx q[0];
rz(-1.6197546) q[0];
sx q[0];
rz(-0.017107054) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2320802) q[2];
sx q[2];
rz(-1.806353) q[2];
sx q[2];
rz(-0.92600694) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0349717) q[1];
sx q[1];
rz(-1.7502516) q[1];
sx q[1];
rz(-0.76645318) q[1];
x q[2];
rz(-1.1398846) q[3];
sx q[3];
rz(-2.0322554) q[3];
sx q[3];
rz(2.0737089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0050469) q[2];
sx q[2];
rz(-1.9133291) q[2];
sx q[2];
rz(-2.6100256) q[2];
rz(0.68938869) q[3];
sx q[3];
rz(-1.4644943) q[3];
sx q[3];
rz(1.2954856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.078995973) q[0];
sx q[0];
rz(-1.2051219) q[0];
sx q[0];
rz(-1.2980365) q[0];
rz(-2.334306) q[1];
sx q[1];
rz(-1.1785945) q[1];
sx q[1];
rz(-0.92179006) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4066276) q[0];
sx q[0];
rz(-1.4359183) q[0];
sx q[0];
rz(-2.448003) q[0];
rz(-pi) q[1];
rz(-1.2577031) q[2];
sx q[2];
rz(-2.3790092) q[2];
sx q[2];
rz(-1.8956172) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.39107716) q[1];
sx q[1];
rz(-2.8120496) q[1];
sx q[1];
rz(2.4604172) q[1];
rz(-pi) q[2];
rz(-2.2790518) q[3];
sx q[3];
rz(-1.4338014) q[3];
sx q[3];
rz(0.18602895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.90199295) q[2];
sx q[2];
rz(-2.5805876) q[2];
sx q[2];
rz(-1.9699338) q[2];
rz(1.7840067) q[3];
sx q[3];
rz(-1.4566908) q[3];
sx q[3];
rz(-1.6931036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8885324) q[0];
sx q[0];
rz(-1.1760412) q[0];
sx q[0];
rz(1.6660447) q[0];
rz(-1.6015923) q[1];
sx q[1];
rz(-1.4614636) q[1];
sx q[1];
rz(-2.4618861) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25526991) q[0];
sx q[0];
rz(-0.5605883) q[0];
sx q[0];
rz(-0.38829304) q[0];
rz(-pi) q[1];
rz(-1.3557415) q[2];
sx q[2];
rz(-2.4348223) q[2];
sx q[2];
rz(0.099345318) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3956086) q[1];
sx q[1];
rz(-1.9095699) q[1];
sx q[1];
rz(0.70396522) q[1];
rz(-pi) q[2];
rz(0.71984843) q[3];
sx q[3];
rz(-0.82092972) q[3];
sx q[3];
rz(1.5035226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1264964) q[2];
sx q[2];
rz(-1.6588147) q[2];
sx q[2];
rz(-2.2311907) q[2];
rz(-0.67534584) q[3];
sx q[3];
rz(-2.2048435) q[3];
sx q[3];
rz(2.2909686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(3.0891721) q[0];
sx q[0];
rz(-1.9071254) q[0];
sx q[0];
rz(-2.4269379) q[0];
rz(-2.4275298) q[1];
sx q[1];
rz(-0.95497447) q[1];
sx q[1];
rz(-1.1766599) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3424073) q[0];
sx q[0];
rz(-0.17051324) q[0];
sx q[0];
rz(1.8605581) q[0];
x q[1];
rz(2.5942957) q[2];
sx q[2];
rz(-1.0101057) q[2];
sx q[2];
rz(-1.0123569) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1706108) q[1];
sx q[1];
rz(-0.5628399) q[1];
sx q[1];
rz(-1.4303722) q[1];
x q[2];
rz(0.51104607) q[3];
sx q[3];
rz(-1.4941477) q[3];
sx q[3];
rz(-0.49728909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8979793) q[2];
sx q[2];
rz(-1.672013) q[2];
sx q[2];
rz(-2.0142377) q[2];
rz(2.7838498) q[3];
sx q[3];
rz(-2.2704411) q[3];
sx q[3];
rz(2.0991142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42416278) q[0];
sx q[0];
rz(-1.2107727) q[0];
sx q[0];
rz(-2.6834224) q[0];
rz(0.39623109) q[1];
sx q[1];
rz(-0.025066499) q[1];
sx q[1];
rz(0.33096663) q[1];
rz(1.8489807) q[2];
sx q[2];
rz(-0.85360151) q[2];
sx q[2];
rz(-3.1384946) q[2];
rz(1.5626004) q[3];
sx q[3];
rz(-1.2200439) q[3];
sx q[3];
rz(-0.72190819) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
