OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.26602715) q[0];
sx q[0];
rz(-0.53524435) q[0];
sx q[0];
rz(-2.3875561) q[0];
rz(-2.3513878) q[1];
sx q[1];
rz(-1.9145929) q[1];
sx q[1];
rz(-1.9807281) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1140808) q[0];
sx q[0];
rz(-1.0923166) q[0];
sx q[0];
rz(1.6041243) q[0];
rz(-pi) q[1];
rz(2.0487469) q[2];
sx q[2];
rz(-1.7149601) q[2];
sx q[2];
rz(0.62010566) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7272211) q[1];
sx q[1];
rz(-1.1611658) q[1];
sx q[1];
rz(-1.0135256) q[1];
rz(-2.394862) q[3];
sx q[3];
rz(-0.24198469) q[3];
sx q[3];
rz(0.63499588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6291818) q[2];
sx q[2];
rz(-1.8449731) q[2];
sx q[2];
rz(-0.66649246) q[2];
rz(0.50981057) q[3];
sx q[3];
rz(-1.1504268) q[3];
sx q[3];
rz(1.8681017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3614685) q[0];
sx q[0];
rz(-1.7391917) q[0];
sx q[0];
rz(2.0289039) q[0];
rz(0.15377046) q[1];
sx q[1];
rz(-2.230481) q[1];
sx q[1];
rz(-1.3382834) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7614377) q[0];
sx q[0];
rz(-2.2614711) q[0];
sx q[0];
rz(1.3209016) q[0];
rz(0.83805214) q[2];
sx q[2];
rz(-1.8997955) q[2];
sx q[2];
rz(-2.7555562) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3160623) q[1];
sx q[1];
rz(-1.1989374) q[1];
sx q[1];
rz(-1.4447681) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9017341) q[3];
sx q[3];
rz(-0.60088241) q[3];
sx q[3];
rz(-2.2440653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0587557) q[2];
sx q[2];
rz(-0.78410167) q[2];
sx q[2];
rz(1.9632957) q[2];
rz(0.96238771) q[3];
sx q[3];
rz(-2.048384) q[3];
sx q[3];
rz(-0.66550955) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1798379) q[0];
sx q[0];
rz(-2.6340155) q[0];
sx q[0];
rz(-1.5270365) q[0];
rz(0.64287341) q[1];
sx q[1];
rz(-1.0765272) q[1];
sx q[1];
rz(-2.8082074) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3231455) q[0];
sx q[0];
rz(-0.72431699) q[0];
sx q[0];
rz(1.2562654) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1836428) q[2];
sx q[2];
rz(-1.6946812) q[2];
sx q[2];
rz(-1.3882335) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.15409878) q[1];
sx q[1];
rz(-1.0990267) q[1];
sx q[1];
rz(0.091847329) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3473347) q[3];
sx q[3];
rz(-2.1922605) q[3];
sx q[3];
rz(-2.7623451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0201515) q[2];
sx q[2];
rz(-1.8276428) q[2];
sx q[2];
rz(-0.80231673) q[2];
rz(0.24117593) q[3];
sx q[3];
rz(-0.69883385) q[3];
sx q[3];
rz(0.10087092) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0872831) q[0];
sx q[0];
rz(-1.5232975) q[0];
sx q[0];
rz(-0.56418443) q[0];
rz(0.57812771) q[1];
sx q[1];
rz(-1.4923948) q[1];
sx q[1];
rz(-2.633458) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6619686) q[0];
sx q[0];
rz(-1.9920252) q[0];
sx q[0];
rz(2.7022916) q[0];
rz(-pi) q[1];
rz(2.5035985) q[2];
sx q[2];
rz(-1.5533414) q[2];
sx q[2];
rz(-0.88776112) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9307738) q[1];
sx q[1];
rz(-2.3465119) q[1];
sx q[1];
rz(-2.7706183) q[1];
x q[2];
rz(-2.9594801) q[3];
sx q[3];
rz(-1.84387) q[3];
sx q[3];
rz(1.9338716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6886787) q[2];
sx q[2];
rz(-0.33334175) q[2];
sx q[2];
rz(-0.63344669) q[2];
rz(-2.5417035) q[3];
sx q[3];
rz(-1.9918631) q[3];
sx q[3];
rz(-1.6413123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7552345) q[0];
sx q[0];
rz(-0.30325493) q[0];
sx q[0];
rz(0.19609837) q[0];
rz(1.261699) q[1];
sx q[1];
rz(-2.321107) q[1];
sx q[1];
rz(1.0713779) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3065942) q[0];
sx q[0];
rz(-2.1822565) q[0];
sx q[0];
rz(2.6131265) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.41197889) q[2];
sx q[2];
rz(-2.2569071) q[2];
sx q[2];
rz(-2.0737322) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6043678) q[1];
sx q[1];
rz(-2.62694) q[1];
sx q[1];
rz(-2.1746219) q[1];
x q[2];
rz(2.0282182) q[3];
sx q[3];
rz(-1.1856106) q[3];
sx q[3];
rz(0.53264632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.034996899) q[2];
sx q[2];
rz(-1.804616) q[2];
sx q[2];
rz(-2.1595947) q[2];
rz(-2.9563831) q[3];
sx q[3];
rz(-0.84398142) q[3];
sx q[3];
rz(1.265032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-1.6035325) q[0];
sx q[0];
rz(-2.2221727) q[0];
sx q[0];
rz(2.61125) q[0];
rz(1.2999339) q[1];
sx q[1];
rz(-1.8118186) q[1];
sx q[1];
rz(-2.9249654) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.075561698) q[0];
sx q[0];
rz(-2.3530585) q[0];
sx q[0];
rz(-1.9476452) q[0];
rz(-1.2712619) q[2];
sx q[2];
rz(-0.96200633) q[2];
sx q[2];
rz(-1.8744206) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4649012) q[1];
sx q[1];
rz(-1.5790107) q[1];
sx q[1];
rz(1.2354922) q[1];
rz(-pi) q[2];
x q[2];
rz(0.90029193) q[3];
sx q[3];
rz(-1.4196463) q[3];
sx q[3];
rz(-2.8063262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4914322) q[2];
sx q[2];
rz(-1.6314793) q[2];
sx q[2];
rz(2.1177297) q[2];
rz(-0.8762382) q[3];
sx q[3];
rz(-0.70228464) q[3];
sx q[3];
rz(1.8700301) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4890471) q[0];
sx q[0];
rz(-1.1627731) q[0];
sx q[0];
rz(0.52893692) q[0];
rz(1.6128929) q[1];
sx q[1];
rz(-1.9493608) q[1];
sx q[1];
rz(-2.0524009) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74635909) q[0];
sx q[0];
rz(-2.8440209) q[0];
sx q[0];
rz(-1.5901106) q[0];
rz(-pi) q[1];
x q[1];
rz(0.94605298) q[2];
sx q[2];
rz(-2.3240528) q[2];
sx q[2];
rz(2.3600876) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.91428077) q[1];
sx q[1];
rz(-3.1078048) q[1];
sx q[1];
rz(-2.950331) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0508556) q[3];
sx q[3];
rz(-1.6999348) q[3];
sx q[3];
rz(-1.4108301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0163394) q[2];
sx q[2];
rz(-2.3374127) q[2];
sx q[2];
rz(1.9160697) q[2];
rz(1.6493753) q[3];
sx q[3];
rz(-1.4368613) q[3];
sx q[3];
rz(0.15587458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4847223) q[0];
sx q[0];
rz(-0.06198922) q[0];
sx q[0];
rz(-0.86762506) q[0];
rz(3.0743657) q[1];
sx q[1];
rz(-2.1104689) q[1];
sx q[1];
rz(0.19518383) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0941489) q[0];
sx q[0];
rz(-0.35956811) q[0];
sx q[0];
rz(2.5387499) q[0];
x q[1];
rz(-0.40114258) q[2];
sx q[2];
rz(-1.5578798) q[2];
sx q[2];
rz(-1.0912947) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6703549) q[1];
sx q[1];
rz(-2.0002623) q[1];
sx q[1];
rz(3.0958789) q[1];
rz(-pi) q[2];
rz(1.0953426) q[3];
sx q[3];
rz(-1.2519149) q[3];
sx q[3];
rz(1.1974602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8937257) q[2];
sx q[2];
rz(-0.21394955) q[2];
sx q[2];
rz(-0.90551886) q[2];
rz(1.1577822) q[3];
sx q[3];
rz(-1.9730622) q[3];
sx q[3];
rz(0.15914966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4147707) q[0];
sx q[0];
rz(-2.045571) q[0];
sx q[0];
rz(1.0409521) q[0];
rz(3.0629311) q[1];
sx q[1];
rz(-2.9610596) q[1];
sx q[1];
rz(-2.7862766) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61475635) q[0];
sx q[0];
rz(-0.38072452) q[0];
sx q[0];
rz(1.8371131) q[0];
rz(-1.9955194) q[2];
sx q[2];
rz(-0.5364843) q[2];
sx q[2];
rz(1.9076965) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.73341093) q[1];
sx q[1];
rz(-1.8942041) q[1];
sx q[1];
rz(-3.0820993) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5411685) q[3];
sx q[3];
rz(-1.416559) q[3];
sx q[3];
rz(-1.0004117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.67655247) q[2];
sx q[2];
rz(-2.5341454) q[2];
sx q[2];
rz(-1.7379649) q[2];
rz(3.1353531) q[3];
sx q[3];
rz(-1.3429567) q[3];
sx q[3];
rz(1.6041554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86826098) q[0];
sx q[0];
rz(-1.9737759) q[0];
sx q[0];
rz(1.8433174) q[0];
rz(2.6955993) q[1];
sx q[1];
rz(-1.1152209) q[1];
sx q[1];
rz(-2.840852) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3107662) q[0];
sx q[0];
rz(-0.58861536) q[0];
sx q[0];
rz(-1.1525843) q[0];
rz(0.61873318) q[2];
sx q[2];
rz(-0.64090568) q[2];
sx q[2];
rz(0.59782019) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2749394) q[1];
sx q[1];
rz(-1.5544974) q[1];
sx q[1];
rz(0.29124041) q[1];
rz(2.2311287) q[3];
sx q[3];
rz(-1.2372036) q[3];
sx q[3];
rz(3.0433082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.62844244) q[2];
sx q[2];
rz(-1.0978038) q[2];
sx q[2];
rz(-1.1395617) q[2];
rz(1.3509753) q[3];
sx q[3];
rz(-2.5438178) q[3];
sx q[3];
rz(-1.9679507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.448485) q[0];
sx q[0];
rz(-1.9036475) q[0];
sx q[0];
rz(0.87686476) q[0];
rz(1.7383472) q[1];
sx q[1];
rz(-1.2259903) q[1];
sx q[1];
rz(-1.2798053) q[1];
rz(3.0856037) q[2];
sx q[2];
rz(-1.5724788) q[2];
sx q[2];
rz(2.4168617) q[2];
rz(0.26294796) q[3];
sx q[3];
rz(-0.96649747) q[3];
sx q[3];
rz(-1.4154712) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
