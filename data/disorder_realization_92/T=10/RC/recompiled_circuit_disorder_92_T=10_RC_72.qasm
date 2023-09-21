OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.6264412) q[0];
sx q[0];
rz(3.175088) q[0];
sx q[0];
rz(7.6498084) q[0];
rz(-1.051149) q[1];
sx q[1];
rz(-1.4895952) q[1];
sx q[1];
rz(1.1319914) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11925919) q[0];
sx q[0];
rz(-0.83689892) q[0];
sx q[0];
rz(-1.0975305) q[0];
rz(-pi) q[1];
x q[1];
rz(0.79905431) q[2];
sx q[2];
rz(-2.91215) q[2];
sx q[2];
rz(-0.45861751) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1177897) q[1];
sx q[1];
rz(-1.0100679) q[1];
sx q[1];
rz(-2.6610713) q[1];
x q[2];
rz(2.6996783) q[3];
sx q[3];
rz(-0.60097296) q[3];
sx q[3];
rz(-1.1170944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.91036096) q[2];
sx q[2];
rz(-1.8138764) q[2];
sx q[2];
rz(-1.988391) q[2];
rz(2.6575346) q[3];
sx q[3];
rz(-0.81367937) q[3];
sx q[3];
rz(-2.6822065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83950481) q[0];
sx q[0];
rz(-0.33058259) q[0];
sx q[0];
rz(2.6385345) q[0];
rz(1.5548276) q[1];
sx q[1];
rz(-0.70650548) q[1];
sx q[1];
rz(-0.15393004) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0921558) q[0];
sx q[0];
rz(-1.0679809) q[0];
sx q[0];
rz(-3.1252607) q[0];
x q[1];
rz(1.9956079) q[2];
sx q[2];
rz(-2.4467391) q[2];
sx q[2];
rz(1.4171464) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.82984) q[1];
sx q[1];
rz(-1.1097739) q[1];
sx q[1];
rz(-0.47383576) q[1];
rz(-1.5699584) q[3];
sx q[3];
rz(-2.0805693) q[3];
sx q[3];
rz(-1.2050932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2543891) q[2];
sx q[2];
rz(-2.7837191) q[2];
sx q[2];
rz(-1.8187693) q[2];
rz(-1.4860738) q[3];
sx q[3];
rz(-1.6008987) q[3];
sx q[3];
rz(2.4310908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94674295) q[0];
sx q[0];
rz(-2.0236334) q[0];
sx q[0];
rz(-0.90993607) q[0];
rz(-2.3643156) q[1];
sx q[1];
rz(-0.83559075) q[1];
sx q[1];
rz(2.1562703) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9545427) q[0];
sx q[0];
rz(-1.7304725) q[0];
sx q[0];
rz(-2.2295582) q[0];
x q[1];
rz(-2.2276332) q[2];
sx q[2];
rz(-1.8310391) q[2];
sx q[2];
rz(2.7514806) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.31919033) q[1];
sx q[1];
rz(-1.625758) q[1];
sx q[1];
rz(2.8922006) q[1];
rz(-2.2306973) q[3];
sx q[3];
rz(-2.1351372) q[3];
sx q[3];
rz(2.6499555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.862792) q[2];
sx q[2];
rz(-1.1547487) q[2];
sx q[2];
rz(1.8939691) q[2];
rz(2.6990081) q[3];
sx q[3];
rz(-1.3675888) q[3];
sx q[3];
rz(0.32143337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(1.2414395) q[0];
sx q[0];
rz(-1.0192008) q[0];
sx q[0];
rz(1.1234269) q[0];
rz(-2.5627047) q[1];
sx q[1];
rz(-1.4588979) q[1];
sx q[1];
rz(-1.3935864) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8907341) q[0];
sx q[0];
rz(-0.77461857) q[0];
sx q[0];
rz(-0.61268341) q[0];
x q[1];
rz(-1.4030928) q[2];
sx q[2];
rz(-2.0003013) q[2];
sx q[2];
rz(1.6657366) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2323944) q[1];
sx q[1];
rz(-1.1189788) q[1];
sx q[1];
rz(-2.3934445) q[1];
x q[2];
rz(1.3136775) q[3];
sx q[3];
rz(-0.75206471) q[3];
sx q[3];
rz(-0.63362345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0397772) q[2];
sx q[2];
rz(-1.5559876) q[2];
sx q[2];
rz(-0.0021136443) q[2];
rz(2.5801616) q[3];
sx q[3];
rz(-2.1241472) q[3];
sx q[3];
rz(-2.5533365) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.026022) q[0];
sx q[0];
rz(-2.0176812) q[0];
sx q[0];
rz(-1.9157238) q[0];
rz(1.6745802) q[1];
sx q[1];
rz(-1.8672698) q[1];
sx q[1];
rz(1.7747169) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0910049) q[0];
sx q[0];
rz(-0.43397167) q[0];
sx q[0];
rz(2.6758053) q[0];
x q[1];
rz(-1.6010124) q[2];
sx q[2];
rz(-0.82380166) q[2];
sx q[2];
rz(-3.1291762) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.99202427) q[1];
sx q[1];
rz(-0.7115041) q[1];
sx q[1];
rz(-1.6929388) q[1];
x q[2];
rz(-1.2138052) q[3];
sx q[3];
rz(-1.095466) q[3];
sx q[3];
rz(-2.9514422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2772969) q[2];
sx q[2];
rz(-1.0884476) q[2];
sx q[2];
rz(0.40536353) q[2];
rz(0.46164414) q[3];
sx q[3];
rz(-2.3159537) q[3];
sx q[3];
rz(-1.5464787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49155238) q[0];
sx q[0];
rz(-1.4251645) q[0];
sx q[0];
rz(2.6382085) q[0];
rz(2.9227496) q[1];
sx q[1];
rz(-1.8914521) q[1];
sx q[1];
rz(-2.4898081) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1843684) q[0];
sx q[0];
rz(-1.7608709) q[0];
sx q[0];
rz(-0.24380542) q[0];
x q[1];
rz(1.6794372) q[2];
sx q[2];
rz(-0.75192736) q[2];
sx q[2];
rz(1.3736563) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5946878) q[1];
sx q[1];
rz(-1.6225796) q[1];
sx q[1];
rz(1.4361708) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8586646) q[3];
sx q[3];
rz(-1.6698716) q[3];
sx q[3];
rz(-2.2968452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.51263222) q[2];
sx q[2];
rz(-1.7444538) q[2];
sx q[2];
rz(0.39166489) q[2];
rz(-0.20646778) q[3];
sx q[3];
rz(-2.4075017) q[3];
sx q[3];
rz(0.68968836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3702635) q[0];
sx q[0];
rz(-1.4498793) q[0];
sx q[0];
rz(-0.96965924) q[0];
rz(-0.58352739) q[1];
sx q[1];
rz(-2.0136166) q[1];
sx q[1];
rz(-0.13024174) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9758494) q[0];
sx q[0];
rz(-1.9231057) q[0];
sx q[0];
rz(-1.5177112) q[0];
rz(0.70003216) q[2];
sx q[2];
rz(-2.3398952) q[2];
sx q[2];
rz(-1.1449555) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.11634593) q[1];
sx q[1];
rz(-1.3533918) q[1];
sx q[1];
rz(-2.6245963) q[1];
x q[2];
rz(2.7534915) q[3];
sx q[3];
rz(-1.023205) q[3];
sx q[3];
rz(-2.2263118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.51817259) q[2];
sx q[2];
rz(-1.5815846) q[2];
sx q[2];
rz(-2.0557892) q[2];
rz(3.0893677) q[3];
sx q[3];
rz(-1.6970535) q[3];
sx q[3];
rz(1.0857371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4645585) q[0];
sx q[0];
rz(-0.98419404) q[0];
sx q[0];
rz(1.9751256) q[0];
rz(0.72558609) q[1];
sx q[1];
rz(-1.2936932) q[1];
sx q[1];
rz(-2.3988147) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9966492) q[0];
sx q[0];
rz(-1.3784694) q[0];
sx q[0];
rz(-2.3924475) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5119945) q[2];
sx q[2];
rz(-2.3592735) q[2];
sx q[2];
rz(2.9462189) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2837977) q[1];
sx q[1];
rz(-0.89110121) q[1];
sx q[1];
rz(2.6074175) q[1];
rz(-pi) q[2];
rz(2.2320896) q[3];
sx q[3];
rz(-1.0854183) q[3];
sx q[3];
rz(-2.4042326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.62721884) q[2];
sx q[2];
rz(-1.6122931) q[2];
sx q[2];
rz(-0.26088866) q[2];
rz(-1.1076814) q[3];
sx q[3];
rz(-1.2872144) q[3];
sx q[3];
rz(1.0816983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35659197) q[0];
sx q[0];
rz(-2.3605425) q[0];
sx q[0];
rz(0.93604952) q[0];
rz(0.014135663) q[1];
sx q[1];
rz(-1.8363258) q[1];
sx q[1];
rz(0.65151185) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4441898) q[0];
sx q[0];
rz(-1.5664926) q[0];
sx q[0];
rz(-1.5944907) q[0];
x q[1];
rz(2.3867943) q[2];
sx q[2];
rz(-0.81364606) q[2];
sx q[2];
rz(2.1458643) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1134909) q[1];
sx q[1];
rz(-0.3416225) q[1];
sx q[1];
rz(-1.2042868) q[1];
rz(-1.9367847) q[3];
sx q[3];
rz(-1.9848739) q[3];
sx q[3];
rz(-2.8556292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8462048) q[2];
sx q[2];
rz(-0.69869763) q[2];
sx q[2];
rz(0.52337581) q[2];
rz(0.30803099) q[3];
sx q[3];
rz(-2.568646) q[3];
sx q[3];
rz(-1.8035005) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.734252) q[0];
sx q[0];
rz(-2.9946406) q[0];
sx q[0];
rz(1.7379606) q[0];
rz(2.667528) q[1];
sx q[1];
rz(-1.6876551) q[1];
sx q[1];
rz(-1.7636991) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5705469) q[0];
sx q[0];
rz(-2.018398) q[0];
sx q[0];
rz(0.96121995) q[0];
rz(-pi) q[1];
rz(-0.7147185) q[2];
sx q[2];
rz(-2.208459) q[2];
sx q[2];
rz(1.4710466) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6973411) q[1];
sx q[1];
rz(-1.8129983) q[1];
sx q[1];
rz(-0.40476207) q[1];
rz(-2.58425) q[3];
sx q[3];
rz(-1.7094269) q[3];
sx q[3];
rz(2.0687452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3906117) q[2];
sx q[2];
rz(-1.5079974) q[2];
sx q[2];
rz(-1.6101458) q[2];
rz(1.6837998) q[3];
sx q[3];
rz(-2.0861574) q[3];
sx q[3];
rz(-0.84993258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67125852) q[0];
sx q[0];
rz(-0.27161921) q[0];
sx q[0];
rz(-2.0441396) q[0];
rz(2.509027) q[1];
sx q[1];
rz(-2.0805151) q[1];
sx q[1];
rz(0.17593304) q[1];
rz(-1.2420281) q[2];
sx q[2];
rz(-2.2148) q[2];
sx q[2];
rz(2.7684341) q[2];
rz(0.73344161) q[3];
sx q[3];
rz(-2.1248795) q[3];
sx q[3];
rz(-2.636573) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
