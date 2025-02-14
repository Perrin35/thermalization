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
rz(-2.610425) q[1];
sx q[1];
rz(0.20294987) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73589486) q[0];
sx q[0];
rz(-0.66758388) q[0];
sx q[0];
rz(-0.9783469) q[0];
x q[1];
rz(-0.0056326515) q[2];
sx q[2];
rz(-1.6597103) q[2];
sx q[2];
rz(-2.2064457) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3280199) q[1];
sx q[1];
rz(-1.3839233) q[1];
sx q[1];
rz(-0.55208167) q[1];
rz(-pi) q[2];
rz(3.1139939) q[3];
sx q[3];
rz(-0.8445794) q[3];
sx q[3];
rz(0.82311326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1927294) q[2];
sx q[2];
rz(-2.296083) q[2];
sx q[2];
rz(-1.3884937) q[2];
rz(0.071831547) q[3];
sx q[3];
rz(-0.51126945) q[3];
sx q[3];
rz(2.3141919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.097505957) q[0];
sx q[0];
rz(-0.16479099) q[0];
sx q[0];
rz(-0.49758115) q[0];
rz(1.7454106) q[1];
sx q[1];
rz(-1.0284245) q[1];
sx q[1];
rz(-2.5713249) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3511757) q[0];
sx q[0];
rz(-2.6714258) q[0];
sx q[0];
rz(-2.8461661) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.360434) q[2];
sx q[2];
rz(-1.7711002) q[2];
sx q[2];
rz(-2.4948073) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.1062938) q[1];
sx q[1];
rz(-0.72162823) q[1];
sx q[1];
rz(-3.0880188) q[1];
x q[2];
rz(0.58742828) q[3];
sx q[3];
rz(-0.96754388) q[3];
sx q[3];
rz(2.5771844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0294864) q[2];
sx q[2];
rz(-1.4542397) q[2];
sx q[2];
rz(2.4278329) q[2];
rz(1.7589689) q[3];
sx q[3];
rz(-0.52240038) q[3];
sx q[3];
rz(-0.4471603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5722028) q[0];
sx q[0];
rz(-0.97219205) q[0];
sx q[0];
rz(0.33706459) q[0];
rz(-0.49346787) q[1];
sx q[1];
rz(-2.4512873) q[1];
sx q[1];
rz(-0.92672551) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0729869) q[0];
sx q[0];
rz(-1.8400917) q[0];
sx q[0];
rz(-0.26871839) q[0];
x q[1];
rz(-2.2271635) q[2];
sx q[2];
rz(-0.85322748) q[2];
sx q[2];
rz(2.2877483) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.87624967) q[1];
sx q[1];
rz(-0.45062765) q[1];
sx q[1];
rz(1.1875115) q[1];
rz(-0.89348578) q[3];
sx q[3];
rz(-1.15128) q[3];
sx q[3];
rz(-2.8887987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9530764) q[2];
sx q[2];
rz(-0.86543721) q[2];
sx q[2];
rz(2.4600929) q[2];
rz(-1.513688) q[3];
sx q[3];
rz(-1.3368006) q[3];
sx q[3];
rz(-2.9708235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82591581) q[0];
sx q[0];
rz(-3.0382394) q[0];
sx q[0];
rz(1.0507677) q[0];
rz(2.5283165) q[1];
sx q[1];
rz(-0.79137099) q[1];
sx q[1];
rz(-1.9176066) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0014718) q[0];
sx q[0];
rz(-1.8996692) q[0];
sx q[0];
rz(-0.66063459) q[0];
rz(-pi) q[1];
rz(2.1992219) q[2];
sx q[2];
rz(-2.0258459) q[2];
sx q[2];
rz(1.9949184) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5889086) q[1];
sx q[1];
rz(-2.1652204) q[1];
sx q[1];
rz(1.5859222) q[1];
x q[2];
rz(1.4091893) q[3];
sx q[3];
rz(-2.1568885) q[3];
sx q[3];
rz(-1.8047892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0497389) q[2];
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
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62395537) q[0];
sx q[0];
rz(-2.1006382) q[0];
sx q[0];
rz(2.5868296) q[0];
rz(-1.293921) q[1];
sx q[1];
rz(-2.7298253) q[1];
sx q[1];
rz(0.88465869) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22065565) q[0];
sx q[0];
rz(-2.1010799) q[0];
sx q[0];
rz(-3.1129989) q[0];
rz(-2.4403202) q[2];
sx q[2];
rz(-1.2040569) q[2];
sx q[2];
rz(1.203095) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.5371478) q[1];
sx q[1];
rz(-1.7824518) q[1];
sx q[1];
rz(1.1011293) q[1];
x q[2];
rz(-1.2908261) q[3];
sx q[3];
rz(-1.9970702) q[3];
sx q[3];
rz(2.4170162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.027792949) q[2];
sx q[2];
rz(-2.5774559) q[2];
sx q[2];
rz(-2.0417058) q[2];
rz(-0.22282985) q[3];
sx q[3];
rz(-1.204071) q[3];
sx q[3];
rz(-3.0075464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9851538) q[0];
sx q[0];
rz(-2.8414861) q[0];
sx q[0];
rz(1.9578178) q[0];
rz(1.8249594) q[1];
sx q[1];
rz(-0.96375179) q[1];
sx q[1];
rz(-2.1060941) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4754921) q[0];
sx q[0];
rz(-2.0543092) q[0];
sx q[0];
rz(2.7782281) q[0];
x q[1];
rz(2.7155128) q[2];
sx q[2];
rz(-2.1392864) q[2];
sx q[2];
rz(-0.3884494) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.70724132) q[1];
sx q[1];
rz(-1.7508669) q[1];
sx q[1];
rz(1.7456013) q[1];
x q[2];
rz(-3.0665056) q[3];
sx q[3];
rz(-2.3737566) q[3];
sx q[3];
rz(-1.4251284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.7684795) q[2];
sx q[2];
rz(-0.36102411) q[2];
sx q[2];
rz(2.7601472) q[2];
rz(1.2162195) q[3];
sx q[3];
rz(-0.65665025) q[3];
sx q[3];
rz(-2.5372964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4878047) q[0];
sx q[0];
rz(-0.31823802) q[0];
sx q[0];
rz(2.8741264) q[0];
rz(1.5221315) q[1];
sx q[1];
rz(-2.5289502) q[1];
sx q[1];
rz(0.47346514) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2758613) q[0];
sx q[0];
rz(-1.9124219) q[0];
sx q[0];
rz(-0.39535687) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6619745) q[2];
sx q[2];
rz(-1.4066182) q[2];
sx q[2];
rz(2.2043383) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.50187868) q[1];
sx q[1];
rz(-2.8913829) q[1];
sx q[1];
rz(-0.48709695) q[1];
rz(-2.7438873) q[3];
sx q[3];
rz(-2.7393712) q[3];
sx q[3];
rz(-0.88170748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.21753103) q[2];
sx q[2];
rz(-1.4950098) q[2];
sx q[2];
rz(-1.7688497) q[2];
rz(0.30630201) q[3];
sx q[3];
rz(-1.0270216) q[3];
sx q[3];
rz(-2.8021804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31996763) q[0];
sx q[0];
rz(-0.29682934) q[0];
sx q[0];
rz(-0.4739652) q[0];
rz(-0.78556806) q[1];
sx q[1];
rz(-2.5626917) q[1];
sx q[1];
rz(2.2483291) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91150099) q[0];
sx q[0];
rz(-2.8727838) q[0];
sx q[0];
rz(-1.5066766) q[0];
x q[1];
rz(2.750326) q[2];
sx q[2];
rz(-0.80384582) q[2];
sx q[2];
rz(0.93144691) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6345659) q[1];
sx q[1];
rz(-2.9302672) q[1];
sx q[1];
rz(-0.15626772) q[1];
rz(-pi) q[2];
x q[2];
rz(0.099670815) q[3];
sx q[3];
rz(-1.1766889) q[3];
sx q[3];
rz(1.3083713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.47024176) q[2];
sx q[2];
rz(-1.6157776) q[2];
sx q[2];
rz(0.5980171) q[2];
rz(2.9371069) q[3];
sx q[3];
rz(-3.0308767) q[3];
sx q[3];
rz(2.428875) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.04190271) q[0];
sx q[0];
rz(-2.1429017) q[0];
sx q[0];
rz(-2.4424851) q[0];
rz(-0.39012575) q[1];
sx q[1];
rz(-2.4603619) q[1];
sx q[1];
rz(2.1652538) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1129393) q[0];
sx q[0];
rz(-1.3138471) q[0];
sx q[0];
rz(2.9963958) q[0];
x q[1];
rz(-0.43760145) q[2];
sx q[2];
rz(-1.8695306) q[2];
sx q[2];
rz(-1.7127812) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7568552) q[1];
sx q[1];
rz(-2.6275674) q[1];
sx q[1];
rz(-1.4207178) q[1];
x q[2];
rz(-2.5849708) q[3];
sx q[3];
rz(-1.3395914) q[3];
sx q[3];
rz(-1.4299666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.98328996) q[2];
sx q[2];
rz(-2.172564) q[2];
sx q[2];
rz(-2.9950673) q[2];
rz(-2.8807785) q[3];
sx q[3];
rz(-1.1365889) q[3];
sx q[3];
rz(-2.8543616) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89676595) q[0];
sx q[0];
rz(-2.792206) q[0];
sx q[0];
rz(-0.31325999) q[0];
rz(0.80097711) q[1];
sx q[1];
rz(-1.5104834) q[1];
sx q[1];
rz(0.40447485) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.069720751) q[0];
sx q[0];
rz(-1.3861462) q[0];
sx q[0];
rz(-1.429256) q[0];
x q[1];
rz(-2.8228797) q[2];
sx q[2];
rz(-2.57215) q[2];
sx q[2];
rz(3.0127243) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8449895) q[1];
sx q[1];
rz(-0.93520404) q[1];
sx q[1];
rz(0.92439009) q[1];
rz(-pi) q[2];
rz(-1.9076212) q[3];
sx q[3];
rz(-1.4416896) q[3];
sx q[3];
rz(-2.5548803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6910088) q[2];
sx q[2];
rz(-2.6477224) q[2];
sx q[2];
rz(-0.79130006) q[2];
rz(0.50296909) q[3];
sx q[3];
rz(-2.0937604) q[3];
sx q[3];
rz(2.3414229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7293411) q[0];
sx q[0];
rz(-1.3164192) q[0];
sx q[0];
rz(-1.3715716) q[0];
rz(-2.5149863) q[1];
sx q[1];
rz(-1.659844) q[1];
sx q[1];
rz(-1.0214092) q[1];
rz(-2.7442197) q[2];
sx q[2];
rz(-2.7418991) q[2];
sx q[2];
rz(-1.9261123) q[2];
rz(-1.5833686) q[3];
sx q[3];
rz(-1.3373806) q[3];
sx q[3];
rz(-1.7718337) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
