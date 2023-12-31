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
rz(0.75403655) q[0];
rz(-5.4929805) q[1];
sx q[1];
rz(5.0561855) q[1];
sx q[1];
rz(8.2639134) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1140808) q[0];
sx q[0];
rz(-1.0923166) q[0];
sx q[0];
rz(1.5374684) q[0];
rz(-pi) q[1];
rz(-0.16205807) q[2];
sx q[2];
rz(-2.0433807) q[2];
sx q[2];
rz(-0.87640793) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4143715) q[1];
sx q[1];
rz(-1.9804269) q[1];
sx q[1];
rz(-1.0135256) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9623904) q[3];
sx q[3];
rz(-1.407302) q[3];
sx q[3];
rz(-1.4737827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.51241088) q[2];
sx q[2];
rz(-1.8449731) q[2];
sx q[2];
rz(2.4751002) q[2];
rz(2.6317821) q[3];
sx q[3];
rz(-1.9911659) q[3];
sx q[3];
rz(-1.2734909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3614685) q[0];
sx q[0];
rz(-1.7391917) q[0];
sx q[0];
rz(-1.1126888) q[0];
rz(2.9878222) q[1];
sx q[1];
rz(-2.230481) q[1];
sx q[1];
rz(-1.8033093) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99909821) q[0];
sx q[0];
rz(-2.4141443) q[0];
sx q[0];
rz(2.8508458) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0427225) q[2];
sx q[2];
rz(-0.79052351) q[2];
sx q[2];
rz(2.3015442) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6515793) q[1];
sx q[1];
rz(-2.7499008) q[1];
sx q[1];
rz(-2.8298122) q[1];
x q[2];
rz(1.409378) q[3];
sx q[3];
rz(-0.98940778) q[3];
sx q[3];
rz(0.60928173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0587557) q[2];
sx q[2];
rz(-2.357491) q[2];
sx q[2];
rz(-1.9632957) q[2];
rz(0.96238771) q[3];
sx q[3];
rz(-1.0932086) q[3];
sx q[3];
rz(-2.4760831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1798379) q[0];
sx q[0];
rz(-2.6340155) q[0];
sx q[0];
rz(-1.6145561) q[0];
rz(2.4987192) q[1];
sx q[1];
rz(-2.0650654) q[1];
sx q[1];
rz(-2.8082074) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51334914) q[0];
sx q[0];
rz(-1.7772563) q[0];
sx q[0];
rz(0.87135656) q[0];
rz(-pi) q[1];
rz(-1.9579499) q[2];
sx q[2];
rz(-1.6946812) q[2];
sx q[2];
rz(1.7533592) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6830605) q[1];
sx q[1];
rz(-1.6525869) q[1];
sx q[1];
rz(1.0973147) q[1];
rz(-1.3473347) q[3];
sx q[3];
rz(-2.1922605) q[3];
sx q[3];
rz(-2.7623451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0201515) q[2];
sx q[2];
rz(-1.8276428) q[2];
sx q[2];
rz(-2.3392759) q[2];
rz(-2.9004167) q[3];
sx q[3];
rz(-0.69883385) q[3];
sx q[3];
rz(-3.0407217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.0543095) q[0];
sx q[0];
rz(-1.6182951) q[0];
sx q[0];
rz(-2.5774082) q[0];
rz(0.57812771) q[1];
sx q[1];
rz(-1.6491978) q[1];
sx q[1];
rz(-0.50813466) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.479624) q[0];
sx q[0];
rz(-1.9920252) q[0];
sx q[0];
rz(2.7022916) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1122909) q[2];
sx q[2];
rz(-2.5033931) q[2];
sx q[2];
rz(-0.65949856) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8453464) q[1];
sx q[1];
rz(-2.2987662) q[1];
sx q[1];
rz(-1.216757) q[1];
x q[2];
rz(-0.99677892) q[3];
sx q[3];
rz(-0.32696163) q[3];
sx q[3];
rz(-0.60861482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4529139) q[2];
sx q[2];
rz(-2.8082509) q[2];
sx q[2];
rz(-0.63344669) q[2];
rz(0.59988919) q[3];
sx q[3];
rz(-1.1497295) q[3];
sx q[3];
rz(-1.5002804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3863581) q[0];
sx q[0];
rz(-2.8383377) q[0];
sx q[0];
rz(-0.19609837) q[0];
rz(-1.261699) q[1];
sx q[1];
rz(-0.82048565) q[1];
sx q[1];
rz(-2.0702147) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6288554) q[0];
sx q[0];
rz(-0.78539408) q[0];
sx q[0];
rz(2.1942755) q[0];
rz(-2.3000556) q[2];
sx q[2];
rz(-1.255799) q[2];
sx q[2];
rz(-0.77301651) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6043678) q[1];
sx q[1];
rz(-2.62694) q[1];
sx q[1];
rz(-0.96697076) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.82810546) q[3];
sx q[3];
rz(-2.5525186) q[3];
sx q[3];
rz(0.38618726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.1065958) q[2];
sx q[2];
rz(-1.804616) q[2];
sx q[2];
rz(-0.98199797) q[2];
rz(2.9563831) q[3];
sx q[3];
rz(-0.84398142) q[3];
sx q[3];
rz(-1.265032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6035325) q[0];
sx q[0];
rz(-0.91941994) q[0];
sx q[0];
rz(-0.53034267) q[0];
rz(1.8416587) q[1];
sx q[1];
rz(-1.329774) q[1];
sx q[1];
rz(0.21662724) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3743065) q[0];
sx q[0];
rz(-1.3067129) q[0];
sx q[0];
rz(2.3229775) q[0];
x q[1];
rz(-2.5112721) q[2];
sx q[2];
rz(-1.3263055) q[2];
sx q[2];
rz(-2.663161) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0121289) q[1];
sx q[1];
rz(-0.33540091) q[1];
sx q[1];
rz(-1.5458376) q[1];
x q[2];
rz(1.8111749) q[3];
sx q[3];
rz(-0.68475311) q[3];
sx q[3];
rz(-1.7184337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4914322) q[2];
sx q[2];
rz(-1.6314793) q[2];
sx q[2];
rz(2.1177297) q[2];
rz(0.8762382) q[3];
sx q[3];
rz(-0.70228464) q[3];
sx q[3];
rz(-1.8700301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4890471) q[0];
sx q[0];
rz(-1.1627731) q[0];
sx q[0];
rz(2.6126557) q[0];
rz(1.6128929) q[1];
sx q[1];
rz(-1.1922319) q[1];
sx q[1];
rz(-1.0891917) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74635909) q[0];
sx q[0];
rz(-2.8440209) q[0];
sx q[0];
rz(-1.5901106) q[0];
rz(-pi) q[1];
rz(-2.2839374) q[2];
sx q[2];
rz(-2.0115888) q[2];
sx q[2];
rz(0.33106523) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.293922) q[1];
sx q[1];
rz(-1.5772181) q[1];
sx q[1];
rz(0.033172219) q[1];
rz(-pi) q[2];
rz(1.7004622) q[3];
sx q[3];
rz(-1.4808169) q[3];
sx q[3];
rz(-2.9699096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.12525325) q[2];
sx q[2];
rz(-0.80417997) q[2];
sx q[2];
rz(1.9160697) q[2];
rz(-1.4922173) q[3];
sx q[3];
rz(-1.4368613) q[3];
sx q[3];
rz(-2.9857181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4847223) q[0];
sx q[0];
rz(-3.0796034) q[0];
sx q[0];
rz(2.2739676) q[0];
rz(0.067226974) q[1];
sx q[1];
rz(-1.0311238) q[1];
sx q[1];
rz(0.19518383) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0458826) q[0];
sx q[0];
rz(-1.3699431) q[0];
sx q[0];
rz(-0.30028371) q[0];
rz(1.5567661) q[2];
sx q[2];
rz(-1.1696891) q[2];
sx q[2];
rz(2.6566128) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6703549) q[1];
sx q[1];
rz(-1.1413304) q[1];
sx q[1];
rz(0.045713748) q[1];
rz(2.786039) q[3];
sx q[3];
rz(-1.1211294) q[3];
sx q[3];
rz(-2.928283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.247867) q[2];
sx q[2];
rz(-2.9276431) q[2];
sx q[2];
rz(0.90551886) q[2];
rz(-1.1577822) q[3];
sx q[3];
rz(-1.9730622) q[3];
sx q[3];
rz(2.982443) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4147707) q[0];
sx q[0];
rz(-2.045571) q[0];
sx q[0];
rz(1.0409521) q[0];
rz(-3.0629311) q[1];
sx q[1];
rz(-0.18053308) q[1];
sx q[1];
rz(-2.7862766) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.937505) q[0];
sx q[0];
rz(-1.6687487) q[0];
sx q[0];
rz(-1.9393001) q[0];
x q[1];
rz(1.9955194) q[2];
sx q[2];
rz(-0.5364843) q[2];
sx q[2];
rz(-1.9076965) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.81845835) q[1];
sx q[1];
rz(-1.6272021) q[1];
sx q[1];
rz(1.2468546) q[1];
x q[2];
rz(-0.26852946) q[3];
sx q[3];
rz(-2.5240457) q[3];
sx q[3];
rz(2.7919046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.67655247) q[2];
sx q[2];
rz(-0.6074473) q[2];
sx q[2];
rz(1.7379649) q[2];
rz(-3.1353531) q[3];
sx q[3];
rz(-1.3429567) q[3];
sx q[3];
rz(-1.6041554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86826098) q[0];
sx q[0];
rz(-1.1678168) q[0];
sx q[0];
rz(1.8433174) q[0];
rz(0.4459933) q[1];
sx q[1];
rz(-2.0263717) q[1];
sx q[1];
rz(0.30074063) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8200127) q[0];
sx q[0];
rz(-1.0386779) q[0];
sx q[0];
rz(2.8768455) q[0];
rz(0.61873318) q[2];
sx q[2];
rz(-0.64090568) q[2];
sx q[2];
rz(-2.5437725) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4423351) q[1];
sx q[1];
rz(-1.2795957) q[1];
sx q[1];
rz(-1.5878116) q[1];
rz(-pi) q[2];
rz(-1.0565287) q[3];
sx q[3];
rz(-2.4132055) q[3];
sx q[3];
rz(1.2700833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.62844244) q[2];
sx q[2];
rz(-1.0978038) q[2];
sx q[2];
rz(1.1395617) q[2];
rz(-1.3509753) q[3];
sx q[3];
rz(-2.5438178) q[3];
sx q[3];
rz(-1.1736419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.448485) q[0];
sx q[0];
rz(-1.9036475) q[0];
sx q[0];
rz(0.87686476) q[0];
rz(1.4032455) q[1];
sx q[1];
rz(-1.9156024) q[1];
sx q[1];
rz(1.8617873) q[1];
rz(1.5724814) q[2];
sx q[2];
rz(-1.6267852) q[2];
sx q[2];
rz(0.84597107) q[2];
rz(-2.1915477) q[3];
sx q[3];
rz(-1.3552356) q[3];
sx q[3];
rz(-3.1380359) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
