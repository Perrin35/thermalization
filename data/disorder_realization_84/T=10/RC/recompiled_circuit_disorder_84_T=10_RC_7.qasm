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
rz(0.79020483) q[1];
sx q[1];
rz(-1.2269998) q[1];
sx q[1];
rz(-1.1608646) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0275118) q[0];
sx q[0];
rz(-2.049276) q[0];
sx q[0];
rz(-1.5374684) q[0];
rz(2.0487469) q[2];
sx q[2];
rz(-1.4266326) q[2];
sx q[2];
rz(-0.62010566) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7418993) q[1];
sx q[1];
rz(-2.0772935) q[1];
sx q[1];
rz(0.472881) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4046895) q[3];
sx q[3];
rz(-1.7475834) q[3];
sx q[3];
rz(3.0151031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6291818) q[2];
sx q[2];
rz(-1.8449731) q[2];
sx q[2];
rz(2.4751002) q[2];
rz(-2.6317821) q[3];
sx q[3];
rz(-1.1504268) q[3];
sx q[3];
rz(1.8681017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78012413) q[0];
sx q[0];
rz(-1.402401) q[0];
sx q[0];
rz(-2.0289039) q[0];
rz(2.9878222) q[1];
sx q[1];
rz(-0.91111168) q[1];
sx q[1];
rz(-1.3382834) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99909821) q[0];
sx q[0];
rz(-0.72744838) q[0];
sx q[0];
rz(0.29074685) q[0];
rz(0.83805214) q[2];
sx q[2];
rz(-1.2417972) q[2];
sx q[2];
rz(-0.38603644) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.2087304) q[1];
sx q[1];
rz(-1.6881697) q[1];
sx q[1];
rz(-0.37456234) q[1];
x q[2];
rz(-1.409378) q[3];
sx q[3];
rz(-2.1521849) q[3];
sx q[3];
rz(0.60928173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0587557) q[2];
sx q[2];
rz(-0.78410167) q[2];
sx q[2];
rz(1.178297) q[2];
rz(0.96238771) q[3];
sx q[3];
rz(-2.048384) q[3];
sx q[3];
rz(-0.66550955) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1798379) q[0];
sx q[0];
rz(-0.50757718) q[0];
sx q[0];
rz(-1.5270365) q[0];
rz(-0.64287341) q[1];
sx q[1];
rz(-2.0650654) q[1];
sx q[1];
rz(0.33338526) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9133559) q[0];
sx q[0];
rz(-0.88909273) q[0];
sx q[0];
rz(-2.8744254) q[0];
rz(1.1836428) q[2];
sx q[2];
rz(-1.4469115) q[2];
sx q[2];
rz(-1.7533592) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4585321) q[1];
sx q[1];
rz(-1.4890057) q[1];
sx q[1];
rz(1.0973147) q[1];
rz(-pi) q[2];
rz(1.3473347) q[3];
sx q[3];
rz(-0.94933214) q[3];
sx q[3];
rz(-2.7623451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0201515) q[2];
sx q[2];
rz(-1.3139498) q[2];
sx q[2];
rz(0.80231673) q[2];
rz(0.24117593) q[3];
sx q[3];
rz(-0.69883385) q[3];
sx q[3];
rz(0.10087092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0543095) q[0];
sx q[0];
rz(-1.5232975) q[0];
sx q[0];
rz(0.56418443) q[0];
rz(-0.57812771) q[1];
sx q[1];
rz(-1.4923948) q[1];
sx q[1];
rz(-0.50813466) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3343737) q[0];
sx q[0];
rz(-0.59893543) q[0];
sx q[0];
rz(0.81143023) q[0];
rz(-pi) q[1];
rz(-0.029301734) q[2];
sx q[2];
rz(-2.5033931) q[2];
sx q[2];
rz(-0.65949856) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5156887) q[1];
sx q[1];
rz(-1.3090033) q[1];
sx q[1];
rz(-0.7598676) q[1];
x q[2];
rz(-2.1448137) q[3];
sx q[3];
rz(-2.814631) q[3];
sx q[3];
rz(-0.60861482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4529139) q[2];
sx q[2];
rz(-0.33334175) q[2];
sx q[2];
rz(-2.508146) q[2];
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
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3863581) q[0];
sx q[0];
rz(-2.8383377) q[0];
sx q[0];
rz(0.19609837) q[0];
rz(1.261699) q[1];
sx q[1];
rz(-0.82048565) q[1];
sx q[1];
rz(-1.0713779) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3065942) q[0];
sx q[0];
rz(-2.1822565) q[0];
sx q[0];
rz(0.5284662) q[0];
x q[1];
rz(-2.0256261) q[2];
sx q[2];
rz(-0.78274667) q[2];
sx q[2];
rz(-2.6775529) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.86712671) q[1];
sx q[1];
rz(-1.1536088) q[1];
sx q[1];
rz(-2.8309114) q[1];
x q[2];
rz(-0.82810546) q[3];
sx q[3];
rz(-2.5525186) q[3];
sx q[3];
rz(-2.7554054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1065958) q[2];
sx q[2];
rz(-1.804616) q[2];
sx q[2];
rz(-0.98199797) q[2];
rz(-0.18520959) q[3];
sx q[3];
rz(-0.84398142) q[3];
sx q[3];
rz(-1.265032) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5380602) q[0];
sx q[0];
rz(-0.91941994) q[0];
sx q[0];
rz(-0.53034267) q[0];
rz(1.2999339) q[1];
sx q[1];
rz(-1.329774) q[1];
sx q[1];
rz(2.9249654) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7672862) q[0];
sx q[0];
rz(-1.3067129) q[0];
sx q[0];
rz(0.81861511) q[0];
rz(-pi) q[1];
x q[1];
rz(0.40041133) q[2];
sx q[2];
rz(-2.4715804) q[2];
sx q[2];
rz(-2.3695721) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0121289) q[1];
sx q[1];
rz(-2.8061917) q[1];
sx q[1];
rz(1.5458376) q[1];
rz(1.8111749) q[3];
sx q[3];
rz(-2.4568395) q[3];
sx q[3];
rz(1.7184337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.65016046) q[2];
sx q[2];
rz(-1.6314793) q[2];
sx q[2];
rz(1.023863) q[2];
rz(-0.8762382) q[3];
sx q[3];
rz(-2.439308) q[3];
sx q[3];
rz(-1.8700301) q[3];
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
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4890471) q[0];
sx q[0];
rz(-1.9788195) q[0];
sx q[0];
rz(2.6126557) q[0];
rz(1.6128929) q[1];
sx q[1];
rz(-1.9493608) q[1];
sx q[1];
rz(1.0891917) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.335621) q[0];
sx q[0];
rz(-1.5764589) q[0];
sx q[0];
rz(1.8683158) q[0];
rz(-pi) q[1];
rz(-0.5577001) q[2];
sx q[2];
rz(-0.93765646) q[2];
sx q[2];
rz(1.5483088) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4186801) q[1];
sx q[1];
rz(-1.5376248) q[1];
sx q[1];
rz(-1.5772217) q[1];
x q[2];
rz(-1.7004622) q[3];
sx q[3];
rz(-1.6607758) q[3];
sx q[3];
rz(-2.9699096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0163394) q[2];
sx q[2];
rz(-0.80417997) q[2];
sx q[2];
rz(1.2255229) q[2];
rz(1.6493753) q[3];
sx q[3];
rz(-1.4368613) q[3];
sx q[3];
rz(-2.9857181) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4847223) q[0];
sx q[0];
rz(-3.0796034) q[0];
sx q[0];
rz(0.86762506) q[0];
rz(0.067226974) q[1];
sx q[1];
rz(-1.0311238) q[1];
sx q[1];
rz(0.19518383) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.09571) q[0];
sx q[0];
rz(-1.3699431) q[0];
sx q[0];
rz(-2.8413089) q[0];
x q[1];
rz(1.5567661) q[2];
sx q[2];
rz(-1.9719035) q[2];
sx q[2];
rz(0.48497981) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.36181297) q[1];
sx q[1];
rz(-0.43174141) q[1];
sx q[1];
rz(-1.6702495) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.04625) q[3];
sx q[3];
rz(-1.2519149) q[3];
sx q[3];
rz(-1.9441324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.247867) q[2];
sx q[2];
rz(-0.21394955) q[2];
sx q[2];
rz(-2.2360738) q[2];
rz(1.1577822) q[3];
sx q[3];
rz(-1.1685305) q[3];
sx q[3];
rz(-0.15914966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4147707) q[0];
sx q[0];
rz(-2.045571) q[0];
sx q[0];
rz(-2.1006405) q[0];
rz(3.0629311) q[1];
sx q[1];
rz(-2.9610596) q[1];
sx q[1];
rz(-2.7862766) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.937505) q[0];
sx q[0];
rz(-1.4728439) q[0];
sx q[0];
rz(1.9393001) q[0];
rz(-pi) q[1];
rz(2.0673429) q[2];
sx q[2];
rz(-1.7830007) q[2];
sx q[2];
rz(-2.4339536) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2229059) q[1];
sx q[1];
rz(-2.8129473) q[1];
sx q[1];
rz(-1.7463643) q[1];
x q[2];
rz(-2.8730632) q[3];
sx q[3];
rz(-0.61754698) q[3];
sx q[3];
rz(-0.34968801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.67655247) q[2];
sx q[2];
rz(-0.6074473) q[2];
sx q[2];
rz(1.4036277) q[2];
rz(-3.1353531) q[3];
sx q[3];
rz(-1.798636) q[3];
sx q[3];
rz(-1.5374373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2733317) q[0];
sx q[0];
rz(-1.9737759) q[0];
sx q[0];
rz(1.2982752) q[0];
rz(-2.6955993) q[1];
sx q[1];
rz(-2.0263717) q[1];
sx q[1];
rz(-2.840852) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38590955) q[0];
sx q[0];
rz(-1.798238) q[0];
sx q[0];
rz(2.118551) q[0];
rz(-0.54603521) q[2];
sx q[2];
rz(-1.9249501) q[2];
sx q[2];
rz(1.6500157) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.69925752) q[1];
sx q[1];
rz(-1.8619969) q[1];
sx q[1];
rz(1.5878116) q[1];
x q[2];
rz(-2.7281076) q[3];
sx q[3];
rz(-2.1889912) q[3];
sx q[3];
rz(1.2234883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.62844244) q[2];
sx q[2];
rz(-1.0978038) q[2];
sx q[2];
rz(-2.002031) q[2];
rz(-1.7906174) q[3];
sx q[3];
rz(-0.59777483) q[3];
sx q[3];
rz(-1.1736419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
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
rz(-1.5724814) q[2];
sx q[2];
rz(-1.5148074) q[2];
sx q[2];
rz(-2.2956216) q[2];
rz(-1.2107522) q[3];
sx q[3];
rz(-2.4891709) q[3];
sx q[3];
rz(-1.8579033) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
