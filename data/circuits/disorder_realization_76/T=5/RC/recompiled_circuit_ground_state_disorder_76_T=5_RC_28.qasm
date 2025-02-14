OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6948833) q[0];
sx q[0];
rz(2.0599685) q[0];
sx q[0];
rz(10.164227) q[0];
rz(2.3820355) q[1];
sx q[1];
rz(-1.3146725) q[1];
sx q[1];
rz(-1.7159599) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8872966) q[0];
sx q[0];
rz(-1.0440517) q[0];
sx q[0];
rz(-0.33381427) q[0];
rz(0.47503586) q[2];
sx q[2];
rz(-2.5940707) q[2];
sx q[2];
rz(0.56625596) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.718218) q[1];
sx q[1];
rz(-1.6483232) q[1];
sx q[1];
rz(1.7192057) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5208779) q[3];
sx q[3];
rz(-0.71943362) q[3];
sx q[3];
rz(2.0572544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.90929675) q[2];
sx q[2];
rz(-0.84315073) q[2];
sx q[2];
rz(-2.9006531) q[2];
rz(-0.029189261) q[3];
sx q[3];
rz(-1.3386644) q[3];
sx q[3];
rz(-1.9624814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9772684) q[0];
sx q[0];
rz(-1.9376396) q[0];
sx q[0];
rz(-0.48467317) q[0];
rz(-2.4618705) q[1];
sx q[1];
rz(-1.8599963) q[1];
sx q[1];
rz(-1.9814804) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80473111) q[0];
sx q[0];
rz(-1.8747678) q[0];
sx q[0];
rz(1.5901523) q[0];
rz(-pi) q[1];
rz(-0.71573513) q[2];
sx q[2];
rz(-1.01075) q[2];
sx q[2];
rz(2.3079688) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1801123) q[1];
sx q[1];
rz(-0.65458502) q[1];
sx q[1];
rz(-2.727319) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2413361) q[3];
sx q[3];
rz(-1.0715535) q[3];
sx q[3];
rz(1.3105621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.29740563) q[2];
sx q[2];
rz(-0.30218267) q[2];
sx q[2];
rz(-2.6837132) q[2];
rz(2.0186021) q[3];
sx q[3];
rz(-0.99959683) q[3];
sx q[3];
rz(2.5751953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8726525) q[0];
sx q[0];
rz(-2.432423) q[0];
sx q[0];
rz(-0.031524468) q[0];
rz(2.8541376) q[1];
sx q[1];
rz(-0.87702409) q[1];
sx q[1];
rz(1.215975) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77338868) q[0];
sx q[0];
rz(-0.88951123) q[0];
sx q[0];
rz(0.70685203) q[0];
rz(-0.38162614) q[2];
sx q[2];
rz(-2.6256621) q[2];
sx q[2];
rz(-0.53781539) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.85092227) q[1];
sx q[1];
rz(-1.8462204) q[1];
sx q[1];
rz(1.8551793) q[1];
rz(-pi) q[2];
rz(0.9312882) q[3];
sx q[3];
rz(-2.1289325) q[3];
sx q[3];
rz(2.1344413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.03881255) q[2];
sx q[2];
rz(-1.5422042) q[2];
sx q[2];
rz(2.3056324) q[2];
rz(1.7838259) q[3];
sx q[3];
rz(-0.72503763) q[3];
sx q[3];
rz(-0.74762216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8589856) q[0];
sx q[0];
rz(-1.405412) q[0];
sx q[0];
rz(3.0614241) q[0];
rz(-1.0824341) q[1];
sx q[1];
rz(-0.23324649) q[1];
sx q[1];
rz(-1.3287883) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4445368) q[0];
sx q[0];
rz(-1.2460684) q[0];
sx q[0];
rz(0.53931196) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3983509) q[2];
sx q[2];
rz(-1.2524464) q[2];
sx q[2];
rz(-2.6686252) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.13212407) q[1];
sx q[1];
rz(-1.4259725) q[1];
sx q[1];
rz(2.0028466) q[1];
rz(-pi) q[2];
rz(2.8058047) q[3];
sx q[3];
rz(-1.3779252) q[3];
sx q[3];
rz(1.646281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7666011) q[2];
sx q[2];
rz(-0.98321715) q[2];
sx q[2];
rz(-0.90744606) q[2];
rz(2.4713016) q[3];
sx q[3];
rz(-0.82796103) q[3];
sx q[3];
rz(1.7214187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2403253) q[0];
sx q[0];
rz(-2.2690052) q[0];
sx q[0];
rz(0.43148828) q[0];
rz(2.0101428) q[1];
sx q[1];
rz(-1.1791469) q[1];
sx q[1];
rz(-2.7010837) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90135306) q[0];
sx q[0];
rz(-2.6640737) q[0];
sx q[0];
rz(-1.3130472) q[0];
rz(-pi) q[1];
rz(-1.155987) q[2];
sx q[2];
rz(-1.7602663) q[2];
sx q[2];
rz(-2.1304325) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2825644) q[1];
sx q[1];
rz(-0.33278123) q[1];
sx q[1];
rz(-0.12092332) q[1];
x q[2];
rz(-2.308485) q[3];
sx q[3];
rz(-0.78189497) q[3];
sx q[3];
rz(1.9946919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0195007) q[2];
sx q[2];
rz(-2.2152405) q[2];
sx q[2];
rz(-0.41395536) q[2];
rz(-1.3881989) q[3];
sx q[3];
rz(-1.5940758) q[3];
sx q[3];
rz(-0.12510124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-2.6269161) q[0];
sx q[0];
rz(-1.7358945) q[0];
sx q[0];
rz(0.93389121) q[0];
rz(-2.4712708) q[1];
sx q[1];
rz(-1.8972081) q[1];
sx q[1];
rz(-1.8062887) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9066228) q[0];
sx q[0];
rz(-0.9714533) q[0];
sx q[0];
rz(-3.0204828) q[0];
rz(-pi) q[1];
rz(-1.5960632) q[2];
sx q[2];
rz(-1.7284365) q[2];
sx q[2];
rz(-2.9767286) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.28999968) q[1];
sx q[1];
rz(-1.4581175) q[1];
sx q[1];
rz(0.022701724) q[1];
rz(-0.63702668) q[3];
sx q[3];
rz(-1.5638132) q[3];
sx q[3];
rz(1.03656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.89207092) q[2];
sx q[2];
rz(-2.8688909) q[2];
sx q[2];
rz(1.8708694) q[2];
rz(1.7719841) q[3];
sx q[3];
rz(-1.9000051) q[3];
sx q[3];
rz(2.6197267) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18913604) q[0];
sx q[0];
rz(-2.554775) q[0];
sx q[0];
rz(-2.9879046) q[0];
rz(-2.9391089) q[1];
sx q[1];
rz(-1.9053562) q[1];
sx q[1];
rz(-2.1930146) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5178793) q[0];
sx q[0];
rz(-2.1800632) q[0];
sx q[0];
rz(1.3246956) q[0];
rz(-1.1467298) q[2];
sx q[2];
rz(-0.68409195) q[2];
sx q[2];
rz(0.37307326) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.13742971) q[1];
sx q[1];
rz(-1.3118801) q[1];
sx q[1];
rz(-2.9015885) q[1];
rz(-pi) q[2];
rz(3.0074869) q[3];
sx q[3];
rz(-0.92536345) q[3];
sx q[3];
rz(0.2303309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.004868) q[2];
sx q[2];
rz(-1.0978881) q[2];
sx q[2];
rz(3.1401805) q[2];
rz(0.062156113) q[3];
sx q[3];
rz(-1.7319738) q[3];
sx q[3];
rz(-0.96127659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.38178) q[0];
sx q[0];
rz(-0.75730046) q[0];
sx q[0];
rz(1.9267474) q[0];
rz(2.3836721) q[1];
sx q[1];
rz(-0.52752033) q[1];
sx q[1];
rz(2.5835999) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76282952) q[0];
sx q[0];
rz(-0.42233322) q[0];
sx q[0];
rz(-0.93393737) q[0];
rz(2.5900389) q[2];
sx q[2];
rz(-1.6948912) q[2];
sx q[2];
rz(-2.6089904) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.8415738) q[1];
sx q[1];
rz(-1.2928559) q[1];
sx q[1];
rz(2.071408) q[1];
rz(-1.2936866) q[3];
sx q[3];
rz(-1.596611) q[3];
sx q[3];
rz(1.9896979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5743635) q[2];
sx q[2];
rz(-1.1939253) q[2];
sx q[2];
rz(-0.26926678) q[2];
rz(-1.9783798) q[3];
sx q[3];
rz(-2.5389157) q[3];
sx q[3];
rz(-2.9009624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6729386) q[0];
sx q[0];
rz(-1.0373632) q[0];
sx q[0];
rz(1.0733676) q[0];
rz(-2.0138373) q[1];
sx q[1];
rz(-2.914371) q[1];
sx q[1];
rz(-0.55666298) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8506691) q[0];
sx q[0];
rz(-1.2525096) q[0];
sx q[0];
rz(-2.9922263) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5684909) q[2];
sx q[2];
rz(-1.0759883) q[2];
sx q[2];
rz(-0.061740969) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4924188) q[1];
sx q[1];
rz(-0.88090501) q[1];
sx q[1];
rz(2.0418732) q[1];
rz(-0.39770522) q[3];
sx q[3];
rz(-2.8880062) q[3];
sx q[3];
rz(1.8613147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2399981) q[2];
sx q[2];
rz(-2.0549213) q[2];
sx q[2];
rz(-1.1617917) q[2];
rz(0.53317541) q[3];
sx q[3];
rz(-0.52459255) q[3];
sx q[3];
rz(-2.9840792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6530782) q[0];
sx q[0];
rz(-0.97452679) q[0];
sx q[0];
rz(0.59481204) q[0];
rz(0.57890233) q[1];
sx q[1];
rz(-1.4756823) q[1];
sx q[1];
rz(-2.4822809) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80583094) q[0];
sx q[0];
rz(-1.7044401) q[0];
sx q[0];
rz(2.5138096) q[0];
x q[1];
rz(2.5393007) q[2];
sx q[2];
rz(-2.4458439) q[2];
sx q[2];
rz(0.93590036) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6281471) q[1];
sx q[1];
rz(-1.7127348) q[1];
sx q[1];
rz(-1.8177086) q[1];
rz(-pi) q[2];
rz(-0.80307428) q[3];
sx q[3];
rz(-2.7702906) q[3];
sx q[3];
rz(2.7367531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.351563) q[2];
sx q[2];
rz(-1.943925) q[2];
sx q[2];
rz(-3.0885546) q[2];
rz(-1.3430345) q[3];
sx q[3];
rz(-1.2767295) q[3];
sx q[3];
rz(3.067335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2863083) q[0];
sx q[0];
rz(-1.7779779) q[0];
sx q[0];
rz(1.5386982) q[0];
rz(-3.042649) q[1];
sx q[1];
rz(-1.3700486) q[1];
sx q[1];
rz(0.055421967) q[1];
rz(1.5794051) q[2];
sx q[2];
rz(-1.9762403) q[2];
sx q[2];
rz(2.5019675) q[2];
rz(-2.4002038) q[3];
sx q[3];
rz(-1.6869873) q[3];
sx q[3];
rz(1.6968873) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
