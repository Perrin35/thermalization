OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6586128) q[0];
sx q[0];
rz(-0.40402544) q[0];
sx q[0];
rz(0.39028302) q[0];
rz(3.1080988) q[1];
sx q[1];
rz(-0.53116763) q[1];
sx q[1];
rz(-0.20294987) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34863198) q[0];
sx q[0];
rz(-1.9237776) q[0];
sx q[0];
rz(-0.99162942) q[0];
rz(-pi) q[1];
x q[1];
rz(3.13596) q[2];
sx q[2];
rz(-1.6597103) q[2];
sx q[2];
rz(-2.2064457) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2704248) q[1];
sx q[1];
rz(-2.1121889) q[1];
sx q[1];
rz(1.7893214) q[1];
x q[2];
rz(-1.6018576) q[3];
sx q[3];
rz(-0.72664562) q[3];
sx q[3];
rz(2.2769312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1927294) q[2];
sx q[2];
rz(-0.84550965) q[2];
sx q[2];
rz(-1.753099) q[2];
rz(-0.071831547) q[3];
sx q[3];
rz(-2.6303232) q[3];
sx q[3];
rz(-0.82740074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
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
rz(0.097505957) q[0];
sx q[0];
rz(-0.16479099) q[0];
sx q[0];
rz(2.6440115) q[0];
rz(1.3961821) q[1];
sx q[1];
rz(-2.1131682) q[1];
sx q[1];
rz(-2.5713249) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0222063) q[0];
sx q[0];
rz(-2.019068) q[0];
sx q[0];
rz(1.4239076) q[0];
rz(-pi) q[1];
rz(0.20469529) q[2];
sx q[2];
rz(-1.7768915) q[2];
sx q[2];
rz(0.88155277) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.1062938) q[1];
sx q[1];
rz(-0.72162823) q[1];
sx q[1];
rz(0.053573805) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2482613) q[3];
sx q[3];
rz(-0.81557298) q[3];
sx q[3];
rz(1.7120685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1121062) q[2];
sx q[2];
rz(-1.4542397) q[2];
sx q[2];
rz(0.71375978) q[2];
rz(1.7589689) q[3];
sx q[3];
rz(-2.6191923) q[3];
sx q[3];
rz(0.4471603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56938982) q[0];
sx q[0];
rz(-0.97219205) q[0];
sx q[0];
rz(-2.8045281) q[0];
rz(2.6481248) q[1];
sx q[1];
rz(-0.69030535) q[1];
sx q[1];
rz(-2.2148671) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42905607) q[0];
sx q[0];
rz(-1.829603) q[0];
sx q[0];
rz(1.2919798) q[0];
rz(2.3078467) q[2];
sx q[2];
rz(-1.0930702) q[2];
sx q[2];
rz(1.9557949) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.265343) q[1];
sx q[1];
rz(-0.45062765) q[1];
sx q[1];
rz(1.9540811) q[1];
rz(-pi) q[2];
rz(0.89348578) q[3];
sx q[3];
rz(-1.15128) q[3];
sx q[3];
rz(-0.25279395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9530764) q[2];
sx q[2];
rz(-0.86543721) q[2];
sx q[2];
rz(-0.68149978) q[2];
rz(1.513688) q[3];
sx q[3];
rz(-1.8047921) q[3];
sx q[3];
rz(0.17076913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3156768) q[0];
sx q[0];
rz(-0.1033533) q[0];
sx q[0];
rz(-2.090825) q[0];
rz(0.61327618) q[1];
sx q[1];
rz(-2.3502217) q[1];
sx q[1];
rz(1.223986) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67658778) q[0];
sx q[0];
rz(-0.95116827) q[0];
sx q[0];
rz(1.1628435) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2649369) q[2];
sx q[2];
rz(-0.75738615) q[2];
sx q[2];
rz(0.11981431) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5526841) q[1];
sx q[1];
rz(-0.97637227) q[1];
sx q[1];
rz(-1.5859222) q[1];
rz(-pi) q[2];
rz(2.5494416) q[3];
sx q[3];
rz(-1.4363406) q[3];
sx q[3];
rz(2.8176702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0497389) q[2];
sx q[2];
rz(-0.78511304) q[2];
sx q[2];
rz(0.66748691) q[2];
rz(-1.5751754) q[3];
sx q[3];
rz(-2.9508041) q[3];
sx q[3];
rz(3.0070087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5176373) q[0];
sx q[0];
rz(-1.0409545) q[0];
sx q[0];
rz(2.5868296) q[0];
rz(1.293921) q[1];
sx q[1];
rz(-0.41176739) q[1];
sx q[1];
rz(0.88465869) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8644476) q[0];
sx q[0];
rz(-0.53098035) q[0];
sx q[0];
rz(1.619521) q[0];
rz(0.53699643) q[2];
sx q[2];
rz(-0.77672138) q[2];
sx q[2];
rz(-3.1077488) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.64120871) q[1];
sx q[1];
rz(-0.51189089) q[1];
sx q[1];
rz(2.0140532) q[1];
rz(-pi) q[2];
rz(-2.5949535) q[3];
sx q[3];
rz(-2.6363576) q[3];
sx q[3];
rz(0.11696493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1137997) q[2];
sx q[2];
rz(-2.5774559) q[2];
sx q[2];
rz(-2.0417058) q[2];
rz(0.22282985) q[3];
sx q[3];
rz(-1.9375216) q[3];
sx q[3];
rz(0.13404624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1564388) q[0];
sx q[0];
rz(-0.30010656) q[0];
sx q[0];
rz(1.9578178) q[0];
rz(-1.8249594) q[1];
sx q[1];
rz(-2.1778409) q[1];
sx q[1];
rz(-2.1060941) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.079662948) q[0];
sx q[0];
rz(-1.890914) q[0];
sx q[0];
rz(2.0826262) q[0];
rz(-pi) q[1];
rz(2.7155128) q[2];
sx q[2];
rz(-2.1392864) q[2];
sx q[2];
rz(-0.3884494) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.246418) q[1];
sx q[1];
rz(-1.3988462) q[1];
sx q[1];
rz(2.9587967) q[1];
rz(0.075087021) q[3];
sx q[3];
rz(-2.3737566) q[3];
sx q[3];
rz(1.7164643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3731132) q[2];
sx q[2];
rz(-2.7805685) q[2];
sx q[2];
rz(2.7601472) q[2];
rz(1.2162195) q[3];
sx q[3];
rz(-0.65665025) q[3];
sx q[3];
rz(0.60429627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65378791) q[0];
sx q[0];
rz(-0.31823802) q[0];
sx q[0];
rz(2.8741264) q[0];
rz(-1.6194612) q[1];
sx q[1];
rz(-0.61264241) q[1];
sx q[1];
rz(2.6681275) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8657314) q[0];
sx q[0];
rz(-1.9124219) q[0];
sx q[0];
rz(-2.7462358) q[0];
x q[1];
rz(2.6619745) q[2];
sx q[2];
rz(-1.4066182) q[2];
sx q[2];
rz(-0.93725433) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5468041) q[1];
sx q[1];
rz(-1.4546397) q[1];
sx q[1];
rz(2.919477) q[1];
x q[2];
rz(-1.7340937) q[3];
sx q[3];
rz(-1.2015523) q[3];
sx q[3];
rz(0.45344093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9240616) q[2];
sx q[2];
rz(-1.4950098) q[2];
sx q[2];
rz(1.3727429) q[2];
rz(-0.30630201) q[3];
sx q[3];
rz(-2.114571) q[3];
sx q[3];
rz(-2.8021804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31996763) q[0];
sx q[0];
rz(-0.29682934) q[0];
sx q[0];
rz(2.6676275) q[0];
rz(-0.78556806) q[1];
sx q[1];
rz(-2.5626917) q[1];
sx q[1];
rz(-0.89326352) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2300917) q[0];
sx q[0];
rz(-2.8727838) q[0];
sx q[0];
rz(-1.5066766) q[0];
rz(-0.39126663) q[2];
sx q[2];
rz(-2.3377468) q[2];
sx q[2];
rz(2.2101457) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.66679044) q[1];
sx q[1];
rz(-1.7795085) q[1];
sx q[1];
rz(1.6041715) q[1];
x q[2];
rz(1.8056554) q[3];
sx q[3];
rz(-2.7357172) q[3];
sx q[3];
rz(-2.0879951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6713509) q[2];
sx q[2];
rz(-1.6157776) q[2];
sx q[2];
rz(-0.5980171) q[2];
rz(-0.20448576) q[3];
sx q[3];
rz(-3.0308767) q[3];
sx q[3];
rz(-0.71271768) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0996899) q[0];
sx q[0];
rz(-2.1429017) q[0];
sx q[0];
rz(-2.4424851) q[0];
rz(-0.39012575) q[1];
sx q[1];
rz(-2.4603619) q[1];
sx q[1];
rz(-0.9763388) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5792865) q[0];
sx q[0];
rz(-1.7111943) q[0];
sx q[0];
rz(1.3112351) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.43760145) q[2];
sx q[2];
rz(-1.8695306) q[2];
sx q[2];
rz(1.4288115) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7568552) q[1];
sx q[1];
rz(-2.6275674) q[1];
sx q[1];
rz(-1.4207178) q[1];
rz(0.41918079) q[3];
sx q[3];
rz(-2.5435735) q[3];
sx q[3];
rz(2.647915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1583027) q[2];
sx q[2];
rz(-2.172564) q[2];
sx q[2];
rz(0.14652531) q[2];
rz(2.8807785) q[3];
sx q[3];
rz(-2.0050037) q[3];
sx q[3];
rz(0.28723106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2448267) q[0];
sx q[0];
rz(-0.34938669) q[0];
sx q[0];
rz(-2.8283327) q[0];
rz(2.3406155) q[1];
sx q[1];
rz(-1.6311092) q[1];
sx q[1];
rz(-2.7371178) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6666732) q[0];
sx q[0];
rz(-1.4316779) q[0];
sx q[0];
rz(0.1864726) q[0];
rz(-pi) q[1];
rz(-0.54623917) q[2];
sx q[2];
rz(-1.4010385) q[2];
sx q[2];
rz(1.7129829) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8449895) q[1];
sx q[1];
rz(-2.2063886) q[1];
sx q[1];
rz(-2.2172026) q[1];
x q[2];
rz(-1.9076212) q[3];
sx q[3];
rz(-1.4416896) q[3];
sx q[3];
rz(0.58671236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.45058382) q[2];
sx q[2];
rz(-0.49387026) q[2];
sx q[2];
rz(-0.79130006) q[2];
rz(-2.6386236) q[3];
sx q[3];
rz(-1.0478323) q[3];
sx q[3];
rz(0.80016971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7293411) q[0];
sx q[0];
rz(-1.3164192) q[0];
sx q[0];
rz(-1.3715716) q[0];
rz(0.62660632) q[1];
sx q[1];
rz(-1.659844) q[1];
sx q[1];
rz(-1.0214092) q[1];
rz(-0.37143636) q[2];
sx q[2];
rz(-1.4196266) q[2];
sx q[2];
rz(-3.1279903) q[2];
rz(-2.9081591) q[3];
sx q[3];
rz(-1.558565) q[3];
sx q[3];
rz(-0.19812921) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
