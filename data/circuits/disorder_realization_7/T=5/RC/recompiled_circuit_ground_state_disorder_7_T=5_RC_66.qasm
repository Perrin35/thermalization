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
rz(-2.7513096) q[0];
rz(-0.033493869) q[1];
sx q[1];
rz(3.6727603) q[1];
sx q[1];
rz(9.6277278) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6970878) q[0];
sx q[0];
rz(-1.0314419) q[0];
sx q[0];
rz(-0.41467337) q[0];
rz(-pi) q[1];
rz(1.5076981) q[2];
sx q[2];
rz(-0.089091688) q[2];
sx q[2];
rz(-2.2697946) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.81357274) q[1];
sx q[1];
rz(-1.3839233) q[1];
sx q[1];
rz(0.55208167) q[1];
rz(-pi) q[2];
rz(1.539735) q[3];
sx q[3];
rz(-2.414947) q[3];
sx q[3];
rz(-2.2769312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.94886327) q[2];
sx q[2];
rz(-0.84550965) q[2];
sx q[2];
rz(1.3884937) q[2];
rz(0.071831547) q[3];
sx q[3];
rz(-2.6303232) q[3];
sx q[3];
rz(-2.3141919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0440867) q[0];
sx q[0];
rz(-0.16479099) q[0];
sx q[0];
rz(2.6440115) q[0];
rz(-1.7454106) q[1];
sx q[1];
rz(-2.1131682) q[1];
sx q[1];
rz(0.57026774) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48455301) q[0];
sx q[0];
rz(-1.4385106) q[0];
sx q[0];
rz(-0.45251493) q[0];
rz(-1.7811586) q[2];
sx q[2];
rz(-1.3704925) q[2];
sx q[2];
rz(0.64678538) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.03601053) q[1];
sx q[1];
rz(-2.2911628) q[1];
sx q[1];
rz(1.5237113) q[1];
rz(2.2621798) q[3];
sx q[3];
rz(-2.0447404) q[3];
sx q[3];
rz(-2.4963801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1121062) q[2];
sx q[2];
rz(-1.4542397) q[2];
sx q[2];
rz(-2.4278329) q[2];
rz(-1.3826238) q[3];
sx q[3];
rz(-2.6191923) q[3];
sx q[3];
rz(0.4471603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56938982) q[0];
sx q[0];
rz(-2.1694006) q[0];
sx q[0];
rz(2.8045281) q[0];
rz(-2.6481248) q[1];
sx q[1];
rz(-2.4512873) q[1];
sx q[1];
rz(-2.2148671) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42905607) q[0];
sx q[0];
rz(-1.3119896) q[0];
sx q[0];
rz(-1.8496129) q[0];
rz(-pi) q[1];
rz(0.83374597) q[2];
sx q[2];
rz(-1.0930702) q[2];
sx q[2];
rz(1.1857978) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.34636075) q[1];
sx q[1];
rz(-1.7344001) q[1];
sx q[1];
rz(-1.1490046) q[1];
rz(-pi) q[2];
rz(-0.89348578) q[3];
sx q[3];
rz(-1.9903127) q[3];
sx q[3];
rz(2.8887987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1885163) q[2];
sx q[2];
rz(-0.86543721) q[2];
sx q[2];
rz(0.68149978) q[2];
rz(1.513688) q[3];
sx q[3];
rz(-1.8047921) q[3];
sx q[3];
rz(0.17076913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82591581) q[0];
sx q[0];
rz(-3.0382394) q[0];
sx q[0];
rz(1.0507677) q[0];
rz(-0.61327618) q[1];
sx q[1];
rz(-0.79137099) q[1];
sx q[1];
rz(1.223986) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67658778) q[0];
sx q[0];
rz(-0.95116827) q[0];
sx q[0];
rz(1.1628435) q[0];
rz(-0.94237071) q[2];
sx q[2];
rz(-2.0258459) q[2];
sx q[2];
rz(-1.1466743) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5526841) q[1];
sx q[1];
rz(-0.97637227) q[1];
sx q[1];
rz(1.5556704) q[1];
x q[2];
rz(-2.9038186) q[3];
sx q[3];
rz(-0.60543767) q[3];
sx q[3];
rz(-1.0501705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0918538) q[2];
sx q[2];
rz(-0.78511304) q[2];
sx q[2];
rz(-0.66748691) q[2];
rz(1.5751754) q[3];
sx q[3];
rz(-2.9508041) q[3];
sx q[3];
rz(-3.0070087) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62395537) q[0];
sx q[0];
rz(-1.0409545) q[0];
sx q[0];
rz(0.55476302) q[0];
rz(-1.8476716) q[1];
sx q[1];
rz(-2.7298253) q[1];
sx q[1];
rz(-0.88465869) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.805917) q[0];
sx q[0];
rz(-1.5461304) q[0];
sx q[0];
rz(-1.0403344) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4403202) q[2];
sx q[2];
rz(-1.2040569) q[2];
sx q[2];
rz(-1.203095) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.64120871) q[1];
sx q[1];
rz(-0.51189089) q[1];
sx q[1];
rz(1.1275395) q[1];
x q[2];
rz(-0.54663916) q[3];
sx q[3];
rz(-0.50523509) q[3];
sx q[3];
rz(-3.0246277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.027792949) q[2];
sx q[2];
rz(-2.5774559) q[2];
sx q[2];
rz(-2.0417058) q[2];
rz(-2.9187628) q[3];
sx q[3];
rz(-1.9375216) q[3];
sx q[3];
rz(-3.0075464) q[3];
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
rz(1.9851538) q[0];
sx q[0];
rz(-2.8414861) q[0];
sx q[0];
rz(-1.1837748) q[0];
rz(1.8249594) q[1];
sx q[1];
rz(-0.96375179) q[1];
sx q[1];
rz(1.0354985) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4754921) q[0];
sx q[0];
rz(-1.0872835) q[0];
sx q[0];
rz(0.36336459) q[0];
rz(-pi) q[1];
rz(-2.7155128) q[2];
sx q[2];
rz(-2.1392864) q[2];
sx q[2];
rz(0.3884494) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.70724132) q[1];
sx q[1];
rz(-1.3907258) q[1];
sx q[1];
rz(-1.7456013) q[1];
rz(-2.3751665) q[3];
sx q[3];
rz(-1.5186678) q[3];
sx q[3];
rz(2.9418569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3731132) q[2];
sx q[2];
rz(-0.36102411) q[2];
sx q[2];
rz(-2.7601472) q[2];
rz(1.9253731) q[3];
sx q[3];
rz(-2.4849424) q[3];
sx q[3];
rz(0.60429627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4878047) q[0];
sx q[0];
rz(-2.8233546) q[0];
sx q[0];
rz(0.26746622) q[0];
rz(1.6194612) q[1];
sx q[1];
rz(-2.5289502) q[1];
sx q[1];
rz(2.6681275) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7077443) q[0];
sx q[0];
rz(-1.9421541) q[0];
sx q[0];
rz(-1.2030364) q[0];
rz(2.6619745) q[2];
sx q[2];
rz(-1.4066182) q[2];
sx q[2];
rz(-0.93725433) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5468041) q[1];
sx q[1];
rz(-1.686953) q[1];
sx q[1];
rz(-0.22211566) q[1];
x q[2];
rz(-1.4074989) q[3];
sx q[3];
rz(-1.9400404) q[3];
sx q[3];
rz(-2.6881517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9240616) q[2];
sx q[2];
rz(-1.4950098) q[2];
sx q[2];
rz(1.3727429) q[2];
rz(-0.30630201) q[3];
sx q[3];
rz(-2.114571) q[3];
sx q[3];
rz(0.33941227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(-0.31996763) q[0];
sx q[0];
rz(-0.29682934) q[0];
sx q[0];
rz(0.4739652) q[0];
rz(-0.78556806) q[1];
sx q[1];
rz(-2.5626917) q[1];
sx q[1];
rz(2.2483291) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59747229) q[0];
sx q[0];
rz(-1.553778) q[0];
sx q[0];
rz(-1.3025137) q[0];
rz(0.76456531) q[2];
sx q[2];
rz(-1.8489601) q[2];
sx q[2];
rz(2.7810627) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.50702679) q[1];
sx q[1];
rz(-0.21132547) q[1];
sx q[1];
rz(2.9853249) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8056554) q[3];
sx q[3];
rz(-2.7357172) q[3];
sx q[3];
rz(-1.0535976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.47024176) q[2];
sx q[2];
rz(-1.525815) q[2];
sx q[2];
rz(-2.5435756) q[2];
rz(-2.9371069) q[3];
sx q[3];
rz(-3.0308767) q[3];
sx q[3];
rz(-2.428875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.04190271) q[0];
sx q[0];
rz(-0.99869096) q[0];
sx q[0];
rz(0.69910753) q[0];
rz(-2.7514669) q[1];
sx q[1];
rz(-0.68123078) q[1];
sx q[1];
rz(2.1652538) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6481144) q[0];
sx q[0];
rz(-2.8472487) q[0];
sx q[0];
rz(1.0674547) q[0];
rz(-pi) q[1];
rz(1.2430698) q[2];
sx q[2];
rz(-1.9877745) q[2];
sx q[2];
rz(-2.8627739) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.38473746) q[1];
sx q[1];
rz(-0.51402521) q[1];
sx q[1];
rz(1.4207178) q[1];
rz(1.8412718) q[3];
sx q[3];
rz(-2.1109443) q[3];
sx q[3];
rz(3.1407874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1583027) q[2];
sx q[2];
rz(-2.172564) q[2];
sx q[2];
rz(-2.9950673) q[2];
rz(0.26081416) q[3];
sx q[3];
rz(-2.0050037) q[3];
sx q[3];
rz(-0.28723106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89676595) q[0];
sx q[0];
rz(-0.34938669) q[0];
sx q[0];
rz(2.8283327) q[0];
rz(2.3406155) q[1];
sx q[1];
rz(-1.5104834) q[1];
sx q[1];
rz(-0.40447485) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4749194) q[0];
sx q[0];
rz(-1.4316779) q[0];
sx q[0];
rz(-2.9551201) q[0];
rz(0.31871291) q[2];
sx q[2];
rz(-2.57215) q[2];
sx q[2];
rz(-0.12886831) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.39242649) q[1];
sx q[1];
rz(-2.268384) q[1];
sx q[1];
rz(-0.68470271) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1964818) q[3];
sx q[3];
rz(-0.35983837) q[3];
sx q[3];
rz(-1.3364178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.45058382) q[2];
sx q[2];
rz(-2.6477224) q[2];
sx q[2];
rz(0.79130006) q[2];
rz(-2.6386236) q[3];
sx q[3];
rz(-2.0937604) q[3];
sx q[3];
rz(-0.80016971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7293411) q[0];
sx q[0];
rz(-1.3164192) q[0];
sx q[0];
rz(-1.3715716) q[0];
rz(2.5149863) q[1];
sx q[1];
rz(-1.4817487) q[1];
sx q[1];
rz(2.1201835) q[1];
rz(0.39737293) q[2];
sx q[2];
rz(-2.7418991) q[2];
sx q[2];
rz(-1.9261123) q[2];
rz(0.23343352) q[3];
sx q[3];
rz(-1.558565) q[3];
sx q[3];
rz(-0.19812921) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
