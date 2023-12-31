OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5151514) q[0];
sx q[0];
rz(-0.03349537) q[0];
sx q[0];
rz(1.7749696) q[0];
rz(-1.051149) q[1];
sx q[1];
rz(-1.4895952) q[1];
sx q[1];
rz(1.1319914) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53349797) q[0];
sx q[0];
rz(-2.29288) q[0];
sx q[0];
rz(2.6736834) q[0];
rz(-2.9801324) q[2];
sx q[2];
rz(-1.7345288) q[2];
sx q[2];
rz(1.8979567) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1177897) q[1];
sx q[1];
rz(-2.1315247) q[1];
sx q[1];
rz(0.48052131) q[1];
rz(-pi) q[2];
rz(-0.55478401) q[3];
sx q[3];
rz(-1.3265508) q[3];
sx q[3];
rz(-2.3158405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.91036096) q[2];
sx q[2];
rz(-1.8138764) q[2];
sx q[2];
rz(1.988391) q[2];
rz(-0.48405805) q[3];
sx q[3];
rz(-2.3279133) q[3];
sx q[3];
rz(-0.4593862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3020878) q[0];
sx q[0];
rz(-0.33058259) q[0];
sx q[0];
rz(2.6385345) q[0];
rz(-1.5548276) q[1];
sx q[1];
rz(-0.70650548) q[1];
sx q[1];
rz(-2.9876626) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.670823) q[0];
sx q[0];
rz(-1.5851067) q[0];
sx q[0];
rz(-2.0736681) q[0];
x q[1];
rz(0.92127992) q[2];
sx q[2];
rz(-1.3037455) q[2];
sx q[2];
rz(-2.9608179) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9736204) q[1];
sx q[1];
rz(-0.64860839) q[1];
sx q[1];
rz(0.82778511) q[1];
rz(-0.0014988621) q[3];
sx q[3];
rz(-2.631819) q[3];
sx q[3];
rz(1.9382167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2543891) q[2];
sx q[2];
rz(-0.35787359) q[2];
sx q[2];
rz(-1.3228234) q[2];
rz(1.6555188) q[3];
sx q[3];
rz(-1.6008987) q[3];
sx q[3];
rz(2.4310908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1948497) q[0];
sx q[0];
rz(-2.0236334) q[0];
sx q[0];
rz(-0.90993607) q[0];
rz(-0.77727708) q[1];
sx q[1];
rz(-2.3060019) q[1];
sx q[1];
rz(2.1562703) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9604208) q[0];
sx q[0];
rz(-0.67502484) q[0];
sx q[0];
rz(1.8280562) q[0];
rz(-pi) q[1];
rz(-2.2276332) q[2];
sx q[2];
rz(-1.8310391) q[2];
sx q[2];
rz(2.7514806) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6775551) q[1];
sx q[1];
rz(-2.8863393) q[1];
sx q[1];
rz(-2.9222701) q[1];
rz(-pi) q[2];
rz(0.76936929) q[3];
sx q[3];
rz(-0.83988512) q[3];
sx q[3];
rz(1.4589256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.862792) q[2];
sx q[2];
rz(-1.1547487) q[2];
sx q[2];
rz(-1.8939691) q[2];
rz(-2.6990081) q[3];
sx q[3];
rz(-1.3675888) q[3];
sx q[3];
rz(-0.32143337) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9001532) q[0];
sx q[0];
rz(-1.0192008) q[0];
sx q[0];
rz(1.1234269) q[0];
rz(2.5627047) q[1];
sx q[1];
rz(-1.6826948) q[1];
sx q[1];
rz(-1.3935864) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14560315) q[0];
sx q[0];
rz(-1.1568501) q[0];
sx q[0];
rz(2.4664509) q[0];
rz(-pi) q[1];
rz(-1.7384999) q[2];
sx q[2];
rz(-2.0003013) q[2];
sx q[2];
rz(1.475856) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2323944) q[1];
sx q[1];
rz(-1.1189788) q[1];
sx q[1];
rz(0.74814817) q[1];
x q[2];
rz(-1.3136775) q[3];
sx q[3];
rz(-2.3895279) q[3];
sx q[3];
rz(2.5079692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0397772) q[2];
sx q[2];
rz(-1.5856051) q[2];
sx q[2];
rz(-0.0021136443) q[2];
rz(2.5801616) q[3];
sx q[3];
rz(-2.1241472) q[3];
sx q[3];
rz(0.58825618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.026022) q[0];
sx q[0];
rz(-1.1239115) q[0];
sx q[0];
rz(1.9157238) q[0];
rz(1.4670124) q[1];
sx q[1];
rz(-1.8672698) q[1];
sx q[1];
rz(1.3668758) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5850692) q[0];
sx q[0];
rz(-1.9559304) q[0];
sx q[0];
rz(1.3655846) q[0];
x q[1];
rz(-0.74722228) q[2];
sx q[2];
rz(-1.5929654) q[2];
sx q[2];
rz(-1.5378466) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3102476) q[1];
sx q[1];
rz(-2.275895) q[1];
sx q[1];
rz(0.10465937) q[1];
x q[2];
rz(2.6392691) q[3];
sx q[3];
rz(-1.886743) q[3];
sx q[3];
rz(-1.5918921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.86429578) q[2];
sx q[2];
rz(-2.0531451) q[2];
sx q[2];
rz(-2.7362291) q[2];
rz(-0.46164414) q[3];
sx q[3];
rz(-2.3159537) q[3];
sx q[3];
rz(1.5464787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49155238) q[0];
sx q[0];
rz(-1.7164282) q[0];
sx q[0];
rz(-0.50338411) q[0];
rz(0.21884306) q[1];
sx q[1];
rz(-1.2501406) q[1];
sx q[1];
rz(0.65178451) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8021278) q[0];
sx q[0];
rz(-1.8101242) q[0];
sx q[0];
rz(1.7665187) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6794372) q[2];
sx q[2];
rz(-0.75192736) q[2];
sx q[2];
rz(1.3736563) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.030901959) q[1];
sx q[1];
rz(-1.7052403) q[1];
sx q[1];
rz(0.052255222) q[1];
rz(-pi) q[2];
rz(2.8586646) q[3];
sx q[3];
rz(-1.471721) q[3];
sx q[3];
rz(2.2968452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6289604) q[2];
sx q[2];
rz(-1.7444538) q[2];
sx q[2];
rz(-2.7499278) q[2];
rz(-2.9351249) q[3];
sx q[3];
rz(-0.73409096) q[3];
sx q[3];
rz(0.68968836) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3702635) q[0];
sx q[0];
rz(-1.6917133) q[0];
sx q[0];
rz(-0.96965924) q[0];
rz(0.58352739) q[1];
sx q[1];
rz(-1.1279761) q[1];
sx q[1];
rz(-0.13024174) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1286344) q[0];
sx q[0];
rz(-2.7854714) q[0];
sx q[0];
rz(2.998259) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.66871754) q[2];
sx q[2];
rz(-2.0520743) q[2];
sx q[2];
rz(2.1858093) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.11634593) q[1];
sx q[1];
rz(-1.7882008) q[1];
sx q[1];
rz(2.6245963) q[1];
x q[2];
rz(2.7534915) q[3];
sx q[3];
rz(-2.1183876) q[3];
sx q[3];
rz(-0.91528085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6234201) q[2];
sx q[2];
rz(-1.5815846) q[2];
sx q[2];
rz(1.0858034) q[2];
rz(-3.0893677) q[3];
sx q[3];
rz(-1.4445392) q[3];
sx q[3];
rz(1.0857371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(2.4645585) q[0];
sx q[0];
rz(-2.1573986) q[0];
sx q[0];
rz(1.1664671) q[0];
rz(2.4160066) q[1];
sx q[1];
rz(-1.8478994) q[1];
sx q[1];
rz(-2.3988147) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9966492) q[0];
sx q[0];
rz(-1.3784694) q[0];
sx q[0];
rz(-0.74914519) q[0];
x q[1];
rz(-0.6767512) q[2];
sx q[2];
rz(-1.9988212) q[2];
sx q[2];
rz(-1.8523491) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.038866) q[1];
sx q[1];
rz(-2.3042149) q[1];
sx q[1];
rz(2.1329761) q[1];
rz(-pi) q[2];
rz(0.86117427) q[3];
sx q[3];
rz(-2.3434601) q[3];
sx q[3];
rz(2.8482311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5143738) q[2];
sx q[2];
rz(-1.5292995) q[2];
sx q[2];
rz(-0.26088866) q[2];
rz(2.0339113) q[3];
sx q[3];
rz(-1.2872144) q[3];
sx q[3];
rz(1.0816983) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7850007) q[0];
sx q[0];
rz(-0.78105015) q[0];
sx q[0];
rz(-2.2055431) q[0];
rz(0.014135663) q[1];
sx q[1];
rz(-1.8363258) q[1];
sx q[1];
rz(0.65151185) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4441898) q[0];
sx q[0];
rz(-1.5664926) q[0];
sx q[0];
rz(-1.5944907) q[0];
rz(2.198095) q[2];
sx q[2];
rz(-2.128696) q[2];
sx q[2];
rz(0.055659143) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1956605) q[1];
sx q[1];
rz(-1.4504499) q[1];
sx q[1];
rz(-1.2502928) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2048079) q[3];
sx q[3];
rz(-1.9848739) q[3];
sx q[3];
rz(2.8556292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8462048) q[2];
sx q[2];
rz(-0.69869763) q[2];
sx q[2];
rz(2.6182168) q[2];
rz(-0.30803099) q[3];
sx q[3];
rz(-0.57294661) q[3];
sx q[3];
rz(1.3380922) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.734252) q[0];
sx q[0];
rz(-0.14695209) q[0];
sx q[0];
rz(1.7379606) q[0];
rz(2.667528) q[1];
sx q[1];
rz(-1.4539376) q[1];
sx q[1];
rz(1.7636991) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5869851) q[0];
sx q[0];
rz(-2.4025612) q[0];
sx q[0];
rz(2.2686195) q[0];
x q[1];
rz(0.79499598) q[2];
sx q[2];
rz(-2.1254052) q[2];
sx q[2];
rz(0.57658741) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.024152012) q[1];
sx q[1];
rz(-1.1785058) q[1];
sx q[1];
rz(-1.3082318) q[1];
rz(-1.4078478) q[3];
sx q[3];
rz(-2.1221707) q[3];
sx q[3];
rz(2.5577302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3906117) q[2];
sx q[2];
rz(-1.6335952) q[2];
sx q[2];
rz(-1.5314468) q[2];
rz(-1.6837998) q[3];
sx q[3];
rz(-1.0554353) q[3];
sx q[3];
rz(-0.84993258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(2.4703341) q[0];
sx q[0];
rz(-0.27161921) q[0];
sx q[0];
rz(-2.0441396) q[0];
rz(2.509027) q[1];
sx q[1];
rz(-2.0805151) q[1];
sx q[1];
rz(0.17593304) q[1];
rz(-2.73545) q[2];
sx q[2];
rz(-0.71229013) q[2];
sx q[2];
rz(-2.9980414) q[2];
rz(-2.2652763) q[3];
sx q[3];
rz(-2.1764168) q[3];
sx q[3];
rz(-0.62302667) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
