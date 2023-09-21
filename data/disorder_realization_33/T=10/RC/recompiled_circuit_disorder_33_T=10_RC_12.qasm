OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5053951) q[0];
sx q[0];
rz(-2.8656821) q[0];
sx q[0];
rz(1.8338058) q[0];
rz(1.1360599) q[1];
sx q[1];
rz(-0.93568957) q[1];
sx q[1];
rz(-1.5712665) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6558134) q[0];
sx q[0];
rz(-1.2278779) q[0];
sx q[0];
rz(-1.2113843) q[0];
rz(-pi) q[1];
rz(-0.76531305) q[2];
sx q[2];
rz(-2.1047449) q[2];
sx q[2];
rz(-0.050616654) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.95333245) q[1];
sx q[1];
rz(-2.7213875) q[1];
sx q[1];
rz(0.62700595) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0516112) q[3];
sx q[3];
rz(-2.7365723) q[3];
sx q[3];
rz(-0.68457505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2661665) q[2];
sx q[2];
rz(-0.29310075) q[2];
sx q[2];
rz(1.1323294) q[2];
rz(-1.6752361) q[3];
sx q[3];
rz(-1.8050067) q[3];
sx q[3];
rz(-2.1291389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9448626) q[0];
sx q[0];
rz(-0.20962993) q[0];
sx q[0];
rz(-2.9557513) q[0];
rz(0.56022412) q[1];
sx q[1];
rz(-1.8461684) q[1];
sx q[1];
rz(0.21683189) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2730392) q[0];
sx q[0];
rz(-0.93398636) q[0];
sx q[0];
rz(2.7815232) q[0];
rz(-pi) q[1];
rz(0.88044135) q[2];
sx q[2];
rz(-0.67697064) q[2];
sx q[2];
rz(-1.134269) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.617802) q[1];
sx q[1];
rz(-2.4317867) q[1];
sx q[1];
rz(-1.8989423) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.300755) q[3];
sx q[3];
rz(-0.84078046) q[3];
sx q[3];
rz(-2.9527612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.310114) q[2];
sx q[2];
rz(-0.82565132) q[2];
sx q[2];
rz(1.8537834) q[2];
rz(0.76256049) q[3];
sx q[3];
rz(-1.972714) q[3];
sx q[3];
rz(-0.30502239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
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
rz(1.4644311) q[0];
sx q[0];
rz(-0.34496775) q[0];
sx q[0];
rz(2.537354) q[0];
rz(1.3263946) q[1];
sx q[1];
rz(-1.3605958) q[1];
sx q[1];
rz(0.93260971) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37768294) q[0];
sx q[0];
rz(-1.8881067) q[0];
sx q[0];
rz(-0.056563932) q[0];
x q[1];
rz(0.18205299) q[2];
sx q[2];
rz(-1.7943873) q[2];
sx q[2];
rz(1.6857266) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0886791) q[1];
sx q[1];
rz(-2.1267849) q[1];
sx q[1];
rz(1.710379) q[1];
rz(-pi) q[2];
rz(2.695735) q[3];
sx q[3];
rz(-1.0599531) q[3];
sx q[3];
rz(1.0872935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.147826) q[2];
sx q[2];
rz(-2.0596762) q[2];
sx q[2];
rz(-1.0926584) q[2];
rz(2.5993733) q[3];
sx q[3];
rz(-1.0850302) q[3];
sx q[3];
rz(-0.96737635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7820691) q[0];
sx q[0];
rz(-3.0451267) q[0];
sx q[0];
rz(2.6413667) q[0];
rz(2.3362828) q[1];
sx q[1];
rz(-1.9814682) q[1];
sx q[1];
rz(-1.4979699) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7786176) q[0];
sx q[0];
rz(-2.551429) q[0];
sx q[0];
rz(-2.1182563) q[0];
x q[1];
rz(0.056604071) q[2];
sx q[2];
rz(-1.7610234) q[2];
sx q[2];
rz(0.20400001) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.56470358) q[1];
sx q[1];
rz(-2.7895045) q[1];
sx q[1];
rz(-0.7301773) q[1];
rz(0.76969947) q[3];
sx q[3];
rz(-2.5315428) q[3];
sx q[3];
rz(1.9609914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3952289) q[2];
sx q[2];
rz(-2.5791898) q[2];
sx q[2];
rz(-0.70181075) q[2];
rz(0.83135215) q[3];
sx q[3];
rz(-2.1777007) q[3];
sx q[3];
rz(-0.62197661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-0.9005301) q[0];
sx q[0];
rz(-2.5456972) q[0];
sx q[0];
rz(2.3262614) q[0];
rz(-1.5218081) q[1];
sx q[1];
rz(-2.3074469) q[1];
sx q[1];
rz(1.048208) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3807555) q[0];
sx q[0];
rz(-1.3073982) q[0];
sx q[0];
rz(2.8834228) q[0];
x q[1];
rz(0.65152119) q[2];
sx q[2];
rz(-1.5610352) q[2];
sx q[2];
rz(-0.62192813) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.293562) q[1];
sx q[1];
rz(-1.8082431) q[1];
sx q[1];
rz(0.32229396) q[1];
rz(-pi) q[2];
rz(2.1225554) q[3];
sx q[3];
rz(-0.69283797) q[3];
sx q[3];
rz(-2.2321731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.52577019) q[2];
sx q[2];
rz(-0.56695357) q[2];
sx q[2];
rz(1.099951) q[2];
rz(-2.3163017) q[3];
sx q[3];
rz(-1.0422948) q[3];
sx q[3];
rz(0.88551372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0734171) q[0];
sx q[0];
rz(-0.59403479) q[0];
sx q[0];
rz(-2.2391879) q[0];
rz(1.0166608) q[1];
sx q[1];
rz(-2.0817751) q[1];
sx q[1];
rz(-3.0117603) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8311365) q[0];
sx q[0];
rz(-1.2015011) q[0];
sx q[0];
rz(-0.62311689) q[0];
rz(-pi) q[1];
rz(-2.1461357) q[2];
sx q[2];
rz(-1.8838922) q[2];
sx q[2];
rz(2.7196333) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.063720238) q[1];
sx q[1];
rz(-1.2076326) q[1];
sx q[1];
rz(-0.6086463) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5927605) q[3];
sx q[3];
rz(-1.7389113) q[3];
sx q[3];
rz(-2.7041534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8292024) q[2];
sx q[2];
rz(-0.94909334) q[2];
sx q[2];
rz(-0.20425805) q[2];
rz(-1.2060818) q[3];
sx q[3];
rz(-1.5217425) q[3];
sx q[3];
rz(0.23541418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4181353) q[0];
sx q[0];
rz(-1.8122939) q[0];
sx q[0];
rz(1.6947421) q[0];
rz(1.8824668) q[1];
sx q[1];
rz(-0.99021688) q[1];
sx q[1];
rz(2.4553305) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9335564) q[0];
sx q[0];
rz(-2.4718923) q[0];
sx q[0];
rz(-0.74525381) q[0];
rz(-pi) q[1];
rz(0.58418602) q[2];
sx q[2];
rz(-0.91911941) q[2];
sx q[2];
rz(-2.3551031) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4916617) q[1];
sx q[1];
rz(-1.6451391) q[1];
sx q[1];
rz(-1.3851623) q[1];
x q[2];
rz(0.63728441) q[3];
sx q[3];
rz(-0.50989671) q[3];
sx q[3];
rz(-0.84264681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.69616047) q[2];
sx q[2];
rz(-1.7636718) q[2];
sx q[2];
rz(3.1398204) q[2];
rz(-0.56162515) q[3];
sx q[3];
rz(-2.2300945) q[3];
sx q[3];
rz(1.6368438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5381662) q[0];
sx q[0];
rz(-0.68646938) q[0];
sx q[0];
rz(-1.4461393) q[0];
rz(0.7810477) q[1];
sx q[1];
rz(-1.3054409) q[1];
sx q[1];
rz(1.5015645) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43270375) q[0];
sx q[0];
rz(-2.6484657) q[0];
sx q[0];
rz(-1.6119484) q[0];
rz(-pi) q[1];
rz(-3.115032) q[2];
sx q[2];
rz(-1.6325258) q[2];
sx q[2];
rz(-0.82690566) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1135243) q[1];
sx q[1];
rz(-1.4816195) q[1];
sx q[1];
rz(0.32797565) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5128965) q[3];
sx q[3];
rz(-1.0843715) q[3];
sx q[3];
rz(-1.6823671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7897196) q[2];
sx q[2];
rz(-1.7799653) q[2];
sx q[2];
rz(1.3191351) q[2];
rz(1.9296648) q[3];
sx q[3];
rz(-1.2865678) q[3];
sx q[3];
rz(2.8222728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8050352) q[0];
sx q[0];
rz(-0.55258495) q[0];
sx q[0];
rz(1.9375027) q[0];
rz(0.38326344) q[1];
sx q[1];
rz(-2.6158694) q[1];
sx q[1];
rz(-0.35167545) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4708913) q[0];
sx q[0];
rz(-1.1180709) q[0];
sx q[0];
rz(-0.41505138) q[0];
rz(-pi) q[1];
x q[1];
rz(0.36231626) q[2];
sx q[2];
rz(-10*pi/13) q[2];
sx q[2];
rz(-2.97646) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0723567) q[1];
sx q[1];
rz(-1.5178536) q[1];
sx q[1];
rz(2.3318021) q[1];
rz(-1.5185235) q[3];
sx q[3];
rz(-1.6802603) q[3];
sx q[3];
rz(-2.0930406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.3433156) q[2];
sx q[2];
rz(-2.0337992) q[2];
sx q[2];
rz(-1.2822255) q[2];
rz(1.6451689) q[3];
sx q[3];
rz(-1.5346425) q[3];
sx q[3];
rz(-2.0857247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4984109) q[0];
sx q[0];
rz(-1.2675985) q[0];
sx q[0];
rz(2.9472651) q[0];
rz(-2.1037897) q[1];
sx q[1];
rz(-0.56832814) q[1];
sx q[1];
rz(-1.0338354) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94668418) q[0];
sx q[0];
rz(-0.91741981) q[0];
sx q[0];
rz(0.16434591) q[0];
x q[1];
rz(2.1543703) q[2];
sx q[2];
rz(-2.3505031) q[2];
sx q[2];
rz(2.1949878) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1397814) q[1];
sx q[1];
rz(-2.1404631) q[1];
sx q[1];
rz(0.12057481) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4284027) q[3];
sx q[3];
rz(-1.7367559) q[3];
sx q[3];
rz(2.1470269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0620492) q[2];
sx q[2];
rz(-2.1958308) q[2];
sx q[2];
rz(-0.6357843) q[2];
rz(-2.87129) q[3];
sx q[3];
rz(-0.79939866) q[3];
sx q[3];
rz(-1.5283782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4476267) q[0];
sx q[0];
rz(-1.8287369) q[0];
sx q[0];
rz(1.0736314) q[0];
rz(-1.7059965) q[1];
sx q[1];
rz(-1.5789079) q[1];
sx q[1];
rz(0.78067738) q[1];
rz(-1.6171261) q[2];
sx q[2];
rz(-0.6033069) q[2];
sx q[2];
rz(2.6848007) q[2];
rz(1.900832) q[3];
sx q[3];
rz(-1.5041372) q[3];
sx q[3];
rz(2.1048673) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
