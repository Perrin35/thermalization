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
rz(-0.79987502) q[0];
sx q[0];
rz(-0.81698155) q[0];
sx q[0];
rz(2.6843827) q[0];
rz(0.41181052) q[1];
sx q[1];
rz(-1.5195941) q[1];
sx q[1];
rz(-2.5817459) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7871817) q[0];
sx q[0];
rz(-2.5776064) q[0];
sx q[0];
rz(2.5111879) q[0];
rz(-2.0737846) q[2];
sx q[2];
rz(-1.3185698) q[2];
sx q[2];
rz(-2.5258738) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3774041) q[1];
sx q[1];
rz(-1.2930096) q[1];
sx q[1];
rz(-0.51731108) q[1];
rz(-pi) q[2];
rz(-2.9212037) q[3];
sx q[3];
rz(-1.7145836) q[3];
sx q[3];
rz(-1.4291562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7655699) q[2];
sx q[2];
rz(-2.1799808) q[2];
sx q[2];
rz(-1.8454856) q[2];
rz(2.5206595) q[3];
sx q[3];
rz(-2.4758078) q[3];
sx q[3];
rz(-0.92500979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45944443) q[0];
sx q[0];
rz(-0.58732533) q[0];
sx q[0];
rz(2.4272954) q[0];
rz(0.722305) q[1];
sx q[1];
rz(-1.0853094) q[1];
sx q[1];
rz(1.3246271) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2091529) q[0];
sx q[0];
rz(-0.55190933) q[0];
sx q[0];
rz(2.1789684) q[0];
x q[1];
rz(-1.6257951) q[2];
sx q[2];
rz(-0.83874133) q[2];
sx q[2];
rz(0.74541237) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.63571841) q[1];
sx q[1];
rz(-1.4117452) q[1];
sx q[1];
rz(-2.768874) q[1];
x q[2];
rz(1.7198661) q[3];
sx q[3];
rz(-2.4507634) q[3];
sx q[3];
rz(2.8579503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.1076727) q[2];
sx q[2];
rz(-1.718797) q[2];
sx q[2];
rz(-0.99204341) q[2];
rz(2.3658559) q[3];
sx q[3];
rz(-0.78191596) q[3];
sx q[3];
rz(-2.0591263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42763212) q[0];
sx q[0];
rz(-1.3466703) q[0];
sx q[0];
rz(-2.2295075) q[0];
rz(2.7569547) q[1];
sx q[1];
rz(-2.1802528) q[1];
sx q[1];
rz(-0.49547637) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81264436) q[0];
sx q[0];
rz(-1.5248796) q[0];
sx q[0];
rz(1.4981235) q[0];
rz(-0.74929535) q[2];
sx q[2];
rz(-1.7118508) q[2];
sx q[2];
rz(-0.8666477) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.8532995) q[1];
sx q[1];
rz(-1.6213196) q[1];
sx q[1];
rz(0.72447296) q[1];
rz(-pi) q[2];
rz(-3.1285237) q[3];
sx q[3];
rz(-2.1926697) q[3];
sx q[3];
rz(2.6789846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.89874011) q[2];
sx q[2];
rz(-1.1246559) q[2];
sx q[2];
rz(2.7070572) q[2];
rz(-0.09856002) q[3];
sx q[3];
rz(-1.7193272) q[3];
sx q[3];
rz(-2.3677473) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74715215) q[0];
sx q[0];
rz(-2.9065865) q[0];
sx q[0];
rz(-2.8727942) q[0];
rz(-2.0647743) q[1];
sx q[1];
rz(-1.7348758) q[1];
sx q[1];
rz(1.9487618) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4556325) q[0];
sx q[0];
rz(-1.7436899) q[0];
sx q[0];
rz(2.8894823) q[0];
rz(-pi) q[1];
rz(-0.13612408) q[2];
sx q[2];
rz(-1.7002215) q[2];
sx q[2];
rz(-2.4905175) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.079568) q[1];
sx q[1];
rz(-2.8753198) q[1];
sx q[1];
rz(-0.63527271) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9029023) q[3];
sx q[3];
rz(-2.4260739) q[3];
sx q[3];
rz(-1.3022002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.25527915) q[2];
sx q[2];
rz(-2.2463319) q[2];
sx q[2];
rz(-2.8423584) q[2];
rz(-2.8507161) q[3];
sx q[3];
rz(-1.2521005) q[3];
sx q[3];
rz(-0.65557426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84151477) q[0];
sx q[0];
rz(-2.1617007) q[0];
sx q[0];
rz(0.078068659) q[0];
rz(2.1878751) q[1];
sx q[1];
rz(-2.6225312) q[1];
sx q[1];
rz(1.5740707) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8918607) q[0];
sx q[0];
rz(-1.390269) q[0];
sx q[0];
rz(-2.2523227) q[0];
rz(-pi) q[1];
rz(-0.40549739) q[2];
sx q[2];
rz(-0.73842305) q[2];
sx q[2];
rz(0.15106311) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5213274) q[1];
sx q[1];
rz(-2.5286739) q[1];
sx q[1];
rz(1.6271724) q[1];
x q[2];
rz(-1.5007581) q[3];
sx q[3];
rz(-2.2626855) q[3];
sx q[3];
rz(-0.50229154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8629525) q[2];
sx q[2];
rz(-2.7095257) q[2];
sx q[2];
rz(-2.8734109) q[2];
rz(-2.6523318) q[3];
sx q[3];
rz(-2.0475976) q[3];
sx q[3];
rz(1.4485654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0972524) q[0];
sx q[0];
rz(-3.1179929) q[0];
sx q[0];
rz(-2.6771255) q[0];
rz(-2.6181472) q[1];
sx q[1];
rz(-2.372066) q[1];
sx q[1];
rz(0.73062599) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36408851) q[0];
sx q[0];
rz(-0.83151649) q[0];
sx q[0];
rz(-0.19130188) q[0];
x q[1];
rz(2.4695005) q[2];
sx q[2];
rz(-1.8037829) q[2];
sx q[2];
rz(-1.9547878) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0981104) q[1];
sx q[1];
rz(-0.93096369) q[1];
sx q[1];
rz(-1.0366862) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8898289) q[3];
sx q[3];
rz(-1.6986966) q[3];
sx q[3];
rz(2.3826016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0940493) q[2];
sx q[2];
rz(-1.0796248) q[2];
sx q[2];
rz(3.003982) q[2];
rz(-1.1329457) q[3];
sx q[3];
rz(-1.9652941) q[3];
sx q[3];
rz(1.8784116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9889744) q[0];
sx q[0];
rz(-1.512383) q[0];
sx q[0];
rz(3.1076987) q[0];
rz(-2.7821817) q[1];
sx q[1];
rz(-1.5243328) q[1];
sx q[1];
rz(0.83438897) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2330025) q[0];
sx q[0];
rz(-0.49783537) q[0];
sx q[0];
rz(-2.4497238) q[0];
rz(-pi) q[1];
rz(2.0584201) q[2];
sx q[2];
rz(-2.1994123) q[2];
sx q[2];
rz(2.4982029) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0800054) q[1];
sx q[1];
rz(-1.3179895) q[1];
sx q[1];
rz(3.04313) q[1];
x q[2];
rz(2.4700048) q[3];
sx q[3];
rz(-1.62) q[3];
sx q[3];
rz(-0.68873304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.413201) q[2];
sx q[2];
rz(-0.30343702) q[2];
sx q[2];
rz(-2.0835853) q[2];
rz(-2.6672065) q[3];
sx q[3];
rz(-0.79550231) q[3];
sx q[3];
rz(-1.0559731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
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
rz(-1.0742842) q[0];
sx q[0];
rz(-1.5165167) q[0];
sx q[0];
rz(1.3577331) q[0];
rz(-2.7108497) q[1];
sx q[1];
rz(-1.4746702) q[1];
sx q[1];
rz(2.6745083) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8170094) q[0];
sx q[0];
rz(-1.5276434) q[0];
sx q[0];
rz(3.140121) q[0];
x q[1];
rz(0.70286669) q[2];
sx q[2];
rz(-1.7184259) q[2];
sx q[2];
rz(0.63181782) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.43710762) q[1];
sx q[1];
rz(-0.52981716) q[1];
sx q[1];
rz(3.1049314) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6428493) q[3];
sx q[3];
rz(-0.63910881) q[3];
sx q[3];
rz(-0.022834384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0507386) q[2];
sx q[2];
rz(-0.41768062) q[2];
sx q[2];
rz(-2.0609071) q[2];
rz(-2.4596227) q[3];
sx q[3];
rz(-0.8050279) q[3];
sx q[3];
rz(-2.4431156) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4526116) q[0];
sx q[0];
rz(-2.6230951) q[0];
sx q[0];
rz(-2.3338351) q[0];
rz(-1.0001596) q[1];
sx q[1];
rz(-0.88614416) q[1];
sx q[1];
rz(2.3449786) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33045443) q[0];
sx q[0];
rz(-1.66798) q[0];
sx q[0];
rz(-3.099346) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3299204) q[2];
sx q[2];
rz(-0.95930525) q[2];
sx q[2];
rz(0.73630737) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1353242) q[1];
sx q[1];
rz(-2.709143) q[1];
sx q[1];
rz(-1.2342374) q[1];
x q[2];
rz(1.8851938) q[3];
sx q[3];
rz(-1.3290429) q[3];
sx q[3];
rz(-0.35255656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.54958582) q[2];
sx q[2];
rz(-1.0909811) q[2];
sx q[2];
rz(0.31527147) q[2];
rz(-1.5677876) q[3];
sx q[3];
rz(-1.9939634) q[3];
sx q[3];
rz(-1.1228336) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10298097) q[0];
sx q[0];
rz(-0.6655612) q[0];
sx q[0];
rz(-2.3200206) q[0];
rz(-2.6577677) q[1];
sx q[1];
rz(-1.7211823) q[1];
sx q[1];
rz(2.9169567) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.064817682) q[0];
sx q[0];
rz(-2.7833287) q[0];
sx q[0];
rz(1.129877) q[0];
rz(-2.5098652) q[2];
sx q[2];
rz(-1.016482) q[2];
sx q[2];
rz(0.71507031) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0269751) q[1];
sx q[1];
rz(-1.3624422) q[1];
sx q[1];
rz(2.2592553) q[1];
x q[2];
rz(1.7849237) q[3];
sx q[3];
rz(-1.6928634) q[3];
sx q[3];
rz(1.7904953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3839174) q[2];
sx q[2];
rz(-0.11666798) q[2];
sx q[2];
rz(3.0465916) q[2];
rz(1.654024) q[3];
sx q[3];
rz(-0.58949685) q[3];
sx q[3];
rz(2.3768363) q[3];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2471531) q[0];
sx q[0];
rz(-0.40774397) q[0];
sx q[0];
rz(2.8751873) q[0];
rz(-2.5447625) q[1];
sx q[1];
rz(-1.6886371) q[1];
sx q[1];
rz(1.4773038) q[1];
rz(-1.1372139) q[2];
sx q[2];
rz(-2.4162393) q[2];
sx q[2];
rz(-2.2888714) q[2];
rz(0.0048051759) q[3];
sx q[3];
rz(-0.37552308) q[3];
sx q[3];
rz(-2.7442591) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
