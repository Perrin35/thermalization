OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9893875) q[0];
sx q[0];
rz(-2.9289991) q[0];
sx q[0];
rz(2.4074182) q[0];
rz(-1.9250159) q[1];
sx q[1];
rz(2.7956378) q[1];
sx q[1];
rz(12.591127) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35015955) q[0];
sx q[0];
rz(-1.6140249) q[0];
sx q[0];
rz(-3.079209) q[0];
rz(-pi) q[1];
rz(1.8046298) q[2];
sx q[2];
rz(-1.6811835) q[2];
sx q[2];
rz(1.8821007) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.36944446) q[1];
sx q[1];
rz(-1.2719063) q[1];
sx q[1];
rz(1.5159831) q[1];
rz(1.3583899) q[3];
sx q[3];
rz(-1.6254566) q[3];
sx q[3];
rz(-0.89059356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2572702) q[2];
sx q[2];
rz(-1.4929644) q[2];
sx q[2];
rz(-0.30882588) q[2];
rz(-0.8257927) q[3];
sx q[3];
rz(-2.1335996) q[3];
sx q[3];
rz(-2.3776313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9708213) q[0];
sx q[0];
rz(-1.229137) q[0];
sx q[0];
rz(-1.0276851) q[0];
rz(0.81614256) q[1];
sx q[1];
rz(-1.1183389) q[1];
sx q[1];
rz(2.3291086) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1493133) q[0];
sx q[0];
rz(-1.9198787) q[0];
sx q[0];
rz(-2.9724253) q[0];
rz(-pi) q[1];
rz(-0.37336739) q[2];
sx q[2];
rz(-1.5840216) q[2];
sx q[2];
rz(-2.2227299) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8390439) q[1];
sx q[1];
rz(-1.4329646) q[1];
sx q[1];
rz(2.0853569) q[1];
rz(2.7926718) q[3];
sx q[3];
rz(-1.9066208) q[3];
sx q[3];
rz(0.15062697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5673148) q[2];
sx q[2];
rz(-1.4894166) q[2];
sx q[2];
rz(-0.91892773) q[2];
rz(2.610176) q[3];
sx q[3];
rz(-1.0454949) q[3];
sx q[3];
rz(-0.04118583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70319217) q[0];
sx q[0];
rz(-2.0719318) q[0];
sx q[0];
rz(-1.137314) q[0];
rz(-1.7510341) q[1];
sx q[1];
rz(-2.5847692) q[1];
sx q[1];
rz(0.73296076) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82405381) q[0];
sx q[0];
rz(-1.8630872) q[0];
sx q[0];
rz(-2.7089416) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.82436136) q[2];
sx q[2];
rz(-2.3869053) q[2];
sx q[2];
rz(3.1122308) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2021091) q[1];
sx q[1];
rz(-1.1014448) q[1];
sx q[1];
rz(-2.6497627) q[1];
x q[2];
rz(1.5571345) q[3];
sx q[3];
rz(-2.6928296) q[3];
sx q[3];
rz(-2.6828533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.0040032337) q[2];
sx q[2];
rz(-1.7294451) q[2];
sx q[2];
rz(-1.0388733) q[2];
rz(2.1393356) q[3];
sx q[3];
rz(-1.301731) q[3];
sx q[3];
rz(0.29016289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20659474) q[0];
sx q[0];
rz(-2.0552141) q[0];
sx q[0];
rz(2.2344053) q[0];
rz(0.15779237) q[1];
sx q[1];
rz(-1.0148427) q[1];
sx q[1];
rz(-1.2812322) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3003101) q[0];
sx q[0];
rz(-2.1774946) q[0];
sx q[0];
rz(-2.2187895) q[0];
x q[1];
rz(0.39878504) q[2];
sx q[2];
rz(-2.246309) q[2];
sx q[2];
rz(-0.28093064) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5990178) q[1];
sx q[1];
rz(-0.78831023) q[1];
sx q[1];
rz(-2.5224204) q[1];
rz(-2.2415555) q[3];
sx q[3];
rz(-2.0931912) q[3];
sx q[3];
rz(-0.99985048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8944051) q[2];
sx q[2];
rz(-0.93623585) q[2];
sx q[2];
rz(-1.2897162) q[2];
rz(1.3442518) q[3];
sx q[3];
rz(-1.684609) q[3];
sx q[3];
rz(-2.4414506) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46071389) q[0];
sx q[0];
rz(-2.6432156) q[0];
sx q[0];
rz(1.0182678) q[0];
rz(1.1030039) q[1];
sx q[1];
rz(-0.7119199) q[1];
sx q[1];
rz(-1.8011372) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0809652) q[0];
sx q[0];
rz(-0.61765352) q[0];
sx q[0];
rz(0.65184848) q[0];
rz(-1.6685772) q[2];
sx q[2];
rz(-1.1927529) q[2];
sx q[2];
rz(-0.60286544) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3892438) q[1];
sx q[1];
rz(-1.443814) q[1];
sx q[1];
rz(1.701322) q[1];
rz(-pi) q[2];
rz(-3.0839447) q[3];
sx q[3];
rz(-1.0312005) q[3];
sx q[3];
rz(3.0209013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5420142) q[2];
sx q[2];
rz(-0.88853637) q[2];
sx q[2];
rz(2.1395903) q[2];
rz(0.69139785) q[3];
sx q[3];
rz(-1.9611497) q[3];
sx q[3];
rz(-1.6208167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5387251) q[0];
sx q[0];
rz(-2.0724917) q[0];
sx q[0];
rz(-2.5163203) q[0];
rz(-1.193115) q[1];
sx q[1];
rz(-2.4517877) q[1];
sx q[1];
rz(-0.56732059) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35463542) q[0];
sx q[0];
rz(-2.3428095) q[0];
sx q[0];
rz(1.3243933) q[0];
rz(-pi) q[1];
rz(-2.3151933) q[2];
sx q[2];
rz(-1.5763088) q[2];
sx q[2];
rz(-0.369584) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9822944) q[1];
sx q[1];
rz(-1.5448227) q[1];
sx q[1];
rz(2.1893255) q[1];
rz(-pi) q[2];
rz(0.17622275) q[3];
sx q[3];
rz(-1.7473975) q[3];
sx q[3];
rz(0.11208243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.83711964) q[2];
sx q[2];
rz(-1.8581055) q[2];
sx q[2];
rz(-1.5817969) q[2];
rz(-0.61257735) q[3];
sx q[3];
rz(-1.9809096) q[3];
sx q[3];
rz(-0.15577236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0462069) q[0];
sx q[0];
rz(-1.950773) q[0];
sx q[0];
rz(2.0147391) q[0];
rz(1.2696179) q[1];
sx q[1];
rz(-1.9536628) q[1];
sx q[1];
rz(0.52106214) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2349486) q[0];
sx q[0];
rz(-2.368133) q[0];
sx q[0];
rz(-2.0818772) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6667188) q[2];
sx q[2];
rz(-0.33039364) q[2];
sx q[2];
rz(1.0164193) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4849629) q[1];
sx q[1];
rz(-2.1466675) q[1];
sx q[1];
rz(2.5099975) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.54102202) q[3];
sx q[3];
rz(-1.8086026) q[3];
sx q[3];
rz(2.1062807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1178939) q[2];
sx q[2];
rz(-1.3641027) q[2];
sx q[2];
rz(-1.2269616) q[2];
rz(1.6795233) q[3];
sx q[3];
rz(-1.6242124) q[3];
sx q[3];
rz(2.7000361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-3.0284477) q[0];
sx q[0];
rz(-2.1118836) q[0];
sx q[0];
rz(2.5944769) q[0];
rz(-3.0534577) q[1];
sx q[1];
rz(-1.3476177) q[1];
sx q[1];
rz(0.42253447) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2937909) q[0];
sx q[0];
rz(-1.592684) q[0];
sx q[0];
rz(1.2897911) q[0];
x q[1];
rz(-0.59178517) q[2];
sx q[2];
rz(-0.6066423) q[2];
sx q[2];
rz(-3.0559412) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.40228981) q[1];
sx q[1];
rz(-1.6268432) q[1];
sx q[1];
rz(1.6399553) q[1];
rz(-1.7541774) q[3];
sx q[3];
rz(-1.4327345) q[3];
sx q[3];
rz(2.3473397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9162468) q[2];
sx q[2];
rz(-1.6371181) q[2];
sx q[2];
rz(-0.28727356) q[2];
rz(0.54287994) q[3];
sx q[3];
rz(-0.77016872) q[3];
sx q[3];
rz(-3.0834901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3995689) q[0];
sx q[0];
rz(-2.4574807) q[0];
sx q[0];
rz(-2.3642819) q[0];
rz(0.9264535) q[1];
sx q[1];
rz(-1.1421721) q[1];
sx q[1];
rz(1.7466338) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05714942) q[0];
sx q[0];
rz(-1.2532338) q[0];
sx q[0];
rz(1.9801514) q[0];
rz(-pi) q[1];
rz(-1.3205166) q[2];
sx q[2];
rz(-2.1553073) q[2];
sx q[2];
rz(0.78148851) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3923126) q[1];
sx q[1];
rz(-0.5764851) q[1];
sx q[1];
rz(-2.2063401) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.93392196) q[3];
sx q[3];
rz(-1.0524629) q[3];
sx q[3];
rz(0.76494382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7159783) q[2];
sx q[2];
rz(-1.7691111) q[2];
sx q[2];
rz(2.051029) q[2];
rz(-2.1770832) q[3];
sx q[3];
rz(-2.1301853) q[3];
sx q[3];
rz(-1.2151659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57508093) q[0];
sx q[0];
rz(-2.8031741) q[0];
sx q[0];
rz(-1.6315208) q[0];
rz(-3.1072726) q[1];
sx q[1];
rz(-1.8020554) q[1];
sx q[1];
rz(-2.7095749) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62160174) q[0];
sx q[0];
rz(-1.5120602) q[0];
sx q[0];
rz(-2.1099173) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9218742) q[2];
sx q[2];
rz(-1.6475999) q[2];
sx q[2];
rz(1.6258282) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0610785) q[1];
sx q[1];
rz(-0.17560683) q[1];
sx q[1];
rz(-0.59934111) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.26507399) q[3];
sx q[3];
rz(-1.005583) q[3];
sx q[3];
rz(-3.0967747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4711275) q[2];
sx q[2];
rz(-1.9894783) q[2];
sx q[2];
rz(-2.067789) q[2];
rz(-2.9527169) q[3];
sx q[3];
rz(-0.60868588) q[3];
sx q[3];
rz(-0.3824189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.043924532) q[0];
sx q[0];
rz(-2.8207939) q[0];
sx q[0];
rz(2.4378142) q[0];
rz(1.7882998) q[1];
sx q[1];
rz(-2.4663993) q[1];
sx q[1];
rz(-0.83723062) q[1];
rz(-0.7424217) q[2];
sx q[2];
rz(-1.0992285) q[2];
sx q[2];
rz(0.64115094) q[2];
rz(-1.6394284) q[3];
sx q[3];
rz(-1.1123688) q[3];
sx q[3];
rz(0.49751626) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
