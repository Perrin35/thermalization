OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.4829798) q[0];
sx q[0];
rz(3.5456181) q[0];
sx q[0];
rz(9.0344949) q[0];
rz(3.1080988) q[1];
sx q[1];
rz(-0.53116763) q[1];
sx q[1];
rz(-0.20294987) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7929607) q[0];
sx q[0];
rz(-1.2178151) q[0];
sx q[0];
rz(-0.99162942) q[0];
rz(-pi) q[1];
rz(-1.481881) q[2];
sx q[2];
rz(-1.5764067) q[2];
sx q[2];
rz(2.5054431) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2704248) q[1];
sx q[1];
rz(-2.1121889) q[1];
sx q[1];
rz(1.3522712) q[1];
rz(-pi) q[2];
rz(0.84439028) q[3];
sx q[3];
rz(-1.5501621) q[3];
sx q[3];
rz(0.72935361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1927294) q[2];
sx q[2];
rz(-2.296083) q[2];
sx q[2];
rz(1.3884937) q[2];
rz(-0.071831547) q[3];
sx q[3];
rz(-0.51126945) q[3];
sx q[3];
rz(-2.3141919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0440867) q[0];
sx q[0];
rz(-0.16479099) q[0];
sx q[0];
rz(2.6440115) q[0];
rz(1.7454106) q[1];
sx q[1];
rz(-2.1131682) q[1];
sx q[1];
rz(-0.57026774) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.790417) q[0];
sx q[0];
rz(-0.47016682) q[0];
sx q[0];
rz(-2.8461661) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.79945081) q[2];
sx q[2];
rz(-2.852147) q[2];
sx q[2];
rz(1.4674526) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.035298881) q[1];
sx q[1];
rz(-0.72162823) q[1];
sx q[1];
rz(3.0880188) q[1];
x q[2];
rz(2.2482613) q[3];
sx q[3];
rz(-2.3260197) q[3];
sx q[3];
rz(1.7120685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0294864) q[2];
sx q[2];
rz(-1.4542397) q[2];
sx q[2];
rz(2.4278329) q[2];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5722028) q[0];
sx q[0];
rz(-2.1694006) q[0];
sx q[0];
rz(-2.8045281) q[0];
rz(-2.6481248) q[1];
sx q[1];
rz(-2.4512873) q[1];
sx q[1];
rz(-2.2148671) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7125366) q[0];
sx q[0];
rz(-1.829603) q[0];
sx q[0];
rz(-1.8496129) q[0];
rz(-pi) q[1];
rz(-0.83374597) q[2];
sx q[2];
rz(-2.0485224) q[2];
sx q[2];
rz(-1.9557949) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2973916) q[1];
sx q[1];
rz(-1.1549885) q[1];
sx q[1];
rz(-2.9625921) q[1];
rz(-pi) q[2];
rz(-2.2481069) q[3];
sx q[3];
rz(-1.9903127) q[3];
sx q[3];
rz(-2.8887987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9530764) q[2];
sx q[2];
rz(-2.2761554) q[2];
sx q[2];
rz(0.68149978) q[2];
rz(-1.513688) q[3];
sx q[3];
rz(-1.3368006) q[3];
sx q[3];
rz(0.17076913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82591581) q[0];
sx q[0];
rz(-3.0382394) q[0];
sx q[0];
rz(-2.090825) q[0];
rz(2.5283165) q[1];
sx q[1];
rz(-2.3502217) q[1];
sx q[1];
rz(1.9176066) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.036788615) q[0];
sx q[0];
rz(-2.4147644) q[0];
sx q[0];
rz(0.50755551) q[0];
rz(-pi) q[1];
rz(-0.54398016) q[2];
sx q[2];
rz(-2.1270985) q[2];
sx q[2];
rz(-2.408319) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.026583662) q[1];
sx q[1];
rz(-1.5582651) q[1];
sx q[1];
rz(-0.59447713) q[1];
x q[2];
rz(0.59215109) q[3];
sx q[3];
rz(-1.705252) q[3];
sx q[3];
rz(-0.32392247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0918538) q[2];
sx q[2];
rz(-2.3564796) q[2];
sx q[2];
rz(-0.66748691) q[2];
rz(-1.5751754) q[3];
sx q[3];
rz(-2.9508041) q[3];
sx q[3];
rz(-0.13458399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62395537) q[0];
sx q[0];
rz(-2.1006382) q[0];
sx q[0];
rz(0.55476302) q[0];
rz(-1.8476716) q[1];
sx q[1];
rz(-0.41176739) q[1];
sx q[1];
rz(0.88465869) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8644476) q[0];
sx q[0];
rz(-0.53098035) q[0];
sx q[0];
rz(1.5220716) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6045962) q[2];
sx q[2];
rz(-0.77672138) q[2];
sx q[2];
rz(-0.033843856) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.64120871) q[1];
sx q[1];
rz(-0.51189089) q[1];
sx q[1];
rz(2.0140532) q[1];
x q[2];
rz(-2.7001722) q[3];
sx q[3];
rz(-1.316464) q[3];
sx q[3];
rz(-2.1770432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.1137997) q[2];
sx q[2];
rz(-2.5774559) q[2];
sx q[2];
rz(-1.0998868) q[2];
rz(2.9187628) q[3];
sx q[3];
rz(-1.204071) q[3];
sx q[3];
rz(0.13404624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9851538) q[0];
sx q[0];
rz(-2.8414861) q[0];
sx q[0];
rz(1.1837748) q[0];
rz(-1.8249594) q[1];
sx q[1];
rz(-2.1778409) q[1];
sx q[1];
rz(1.0354985) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4754921) q[0];
sx q[0];
rz(-1.0872835) q[0];
sx q[0];
rz(2.7782281) q[0];
x q[1];
rz(2.7155128) q[2];
sx q[2];
rz(-1.0023062) q[2];
sx q[2];
rz(0.3884494) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.07044) q[1];
sx q[1];
rz(-0.25030085) q[1];
sx q[1];
rz(2.3790199) q[1];
x q[2];
rz(-1.4984958) q[3];
sx q[3];
rz(-2.3359155) q[3];
sx q[3];
rz(-1.8206545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.7684795) q[2];
sx q[2];
rz(-0.36102411) q[2];
sx q[2];
rz(2.7601472) q[2];
rz(-1.2162195) q[3];
sx q[3];
rz(-0.65665025) q[3];
sx q[3];
rz(-0.60429627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.65378791) q[0];
sx q[0];
rz(-0.31823802) q[0];
sx q[0];
rz(-0.26746622) q[0];
rz(1.5221315) q[1];
sx q[1];
rz(-0.61264241) q[1];
sx q[1];
rz(2.6681275) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43384837) q[0];
sx q[0];
rz(-1.9421541) q[0];
sx q[0];
rz(-1.9385563) q[0];
x q[1];
rz(-2.6619745) q[2];
sx q[2];
rz(-1.7349744) q[2];
sx q[2];
rz(-0.93725433) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.50187868) q[1];
sx q[1];
rz(-2.8913829) q[1];
sx q[1];
rz(0.48709695) q[1];
rz(-1.7340937) q[3];
sx q[3];
rz(-1.2015523) q[3];
sx q[3];
rz(-2.6881517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.21753103) q[2];
sx q[2];
rz(-1.4950098) q[2];
sx q[2];
rz(-1.3727429) q[2];
rz(-0.30630201) q[3];
sx q[3];
rz(-1.0270216) q[3];
sx q[3];
rz(-0.33941227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.821625) q[0];
sx q[0];
rz(-0.29682934) q[0];
sx q[0];
rz(0.4739652) q[0];
rz(2.3560246) q[1];
sx q[1];
rz(-0.57890099) q[1];
sx q[1];
rz(-2.2483291) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2300917) q[0];
sx q[0];
rz(-2.8727838) q[0];
sx q[0];
rz(1.6349161) q[0];
rz(-pi) q[1];
rz(2.3770273) q[2];
sx q[2];
rz(-1.2926326) q[2];
sx q[2];
rz(-0.36052997) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.230669) q[1];
sx q[1];
rz(-1.5381457) q[1];
sx q[1];
rz(0.20882512) q[1];
rz(-1.8056554) q[3];
sx q[3];
rz(-0.40587546) q[3];
sx q[3];
rz(1.0535976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6713509) q[2];
sx q[2];
rz(-1.6157776) q[2];
sx q[2];
rz(2.5435756) q[2];
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
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.04190271) q[0];
sx q[0];
rz(-2.1429017) q[0];
sx q[0];
rz(0.69910753) q[0];
rz(0.39012575) q[1];
sx q[1];
rz(-0.68123078) q[1];
sx q[1];
rz(2.1652538) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5623062) q[0];
sx q[0];
rz(-1.7111943) q[0];
sx q[0];
rz(-1.8303575) q[0];
rz(-0.43760145) q[2];
sx q[2];
rz(-1.8695306) q[2];
sx q[2];
rz(-1.7127812) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0551377) q[1];
sx q[1];
rz(-1.4972151) q[1];
sx q[1];
rz(-1.0615968) q[1];
rz(-2.5849708) q[3];
sx q[3];
rz(-1.8020013) q[3];
sx q[3];
rz(-1.7116261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1583027) q[2];
sx q[2];
rz(-0.96902865) q[2];
sx q[2];
rz(-0.14652531) q[2];
rz(-0.26081416) q[3];
sx q[3];
rz(-2.0050037) q[3];
sx q[3];
rz(0.28723106) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2448267) q[0];
sx q[0];
rz(-0.34938669) q[0];
sx q[0];
rz(-2.8283327) q[0];
rz(0.80097711) q[1];
sx q[1];
rz(-1.5104834) q[1];
sx q[1];
rz(0.40447485) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72971463) q[0];
sx q[0];
rz(-2.9094271) q[0];
sx q[0];
rz(-0.64686717) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7687665) q[2];
sx q[2];
rz(-2.1083197) q[2];
sx q[2];
rz(-0.24453577) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7491662) q[1];
sx q[1];
rz(-0.87320864) q[1];
sx q[1];
rz(-0.68470271) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9076212) q[3];
sx q[3];
rz(-1.699903) q[3];
sx q[3];
rz(-0.58671236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6910088) q[2];
sx q[2];
rz(-2.6477224) q[2];
sx q[2];
rz(0.79130006) q[2];
rz(-0.50296909) q[3];
sx q[3];
rz(-1.0478323) q[3];
sx q[3];
rz(-0.80016971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41225152) q[0];
sx q[0];
rz(-1.8251735) q[0];
sx q[0];
rz(1.770021) q[0];
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
rz(1.5833686) q[3];
sx q[3];
rz(-1.8042121) q[3];
sx q[3];
rz(1.369759) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
