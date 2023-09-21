OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.9010889) q[0];
sx q[0];
rz(-0.17740372) q[0];
sx q[0];
rz(1.1343962) q[0];
rz(-1.9534684) q[1];
sx q[1];
rz(-1.0367353) q[1];
sx q[1];
rz(0.66361767) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7952607) q[0];
sx q[0];
rz(-1.5471317) q[0];
sx q[0];
rz(0.46121009) q[0];
x q[1];
rz(-2.7825836) q[2];
sx q[2];
rz(-2.3308672) q[2];
sx q[2];
rz(-2.5100978) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8780958) q[1];
sx q[1];
rz(-2.7090008) q[1];
sx q[1];
rz(2.2357975) q[1];
rz(0.10856467) q[3];
sx q[3];
rz(-1.3285471) q[3];
sx q[3];
rz(-3.0755175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2628281) q[2];
sx q[2];
rz(-2.6972013) q[2];
sx q[2];
rz(3.0905241) q[2];
rz(-2.5845394) q[3];
sx q[3];
rz(-0.80018187) q[3];
sx q[3];
rz(1.5548271) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5550845) q[0];
sx q[0];
rz(-2.3095135) q[0];
sx q[0];
rz(-2.5449975) q[0];
rz(2.3157628) q[1];
sx q[1];
rz(-1.4412216) q[1];
sx q[1];
rz(1.9155496) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37218371) q[0];
sx q[0];
rz(-0.49159494) q[0];
sx q[0];
rz(-2.7093637) q[0];
rz(2.366757) q[2];
sx q[2];
rz(-2.2171387) q[2];
sx q[2];
rz(1.6006084) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1184428) q[1];
sx q[1];
rz(-0.21986248) q[1];
sx q[1];
rz(-1.5688194) q[1];
rz(-pi) q[2];
rz(-2.1434104) q[3];
sx q[3];
rz(-1.5787573) q[3];
sx q[3];
rz(-1.5919459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.79919672) q[2];
sx q[2];
rz(-1.9689955) q[2];
sx q[2];
rz(-2.8105695) q[2];
rz(0.80667574) q[3];
sx q[3];
rz(-1.2253864) q[3];
sx q[3];
rz(1.4276918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1598635) q[0];
sx q[0];
rz(-1.9112497) q[0];
sx q[0];
rz(-1.249041) q[0];
rz(-3.0535835) q[1];
sx q[1];
rz(-2.0188589) q[1];
sx q[1];
rz(-2.1121315) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2217076) q[0];
sx q[0];
rz(-1.789933) q[0];
sx q[0];
rz(0.94421454) q[0];
rz(2.7715893) q[2];
sx q[2];
rz(-2.4577933) q[2];
sx q[2];
rz(1.5497269) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0440812) q[1];
sx q[1];
rz(-1.71002) q[1];
sx q[1];
rz(0.15686762) q[1];
rz(-pi) q[2];
rz(-2.401628) q[3];
sx q[3];
rz(-1.660534) q[3];
sx q[3];
rz(1.2147853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.133698) q[2];
sx q[2];
rz(-1.7241314) q[2];
sx q[2];
rz(2.6339445) q[2];
rz(1.3890022) q[3];
sx q[3];
rz(-0.32998431) q[3];
sx q[3];
rz(1.9784137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65748173) q[0];
sx q[0];
rz(-0.27814516) q[0];
sx q[0];
rz(-1.5456276) q[0];
rz(-2.0987299) q[1];
sx q[1];
rz(-1.1735801) q[1];
sx q[1];
rz(1.625659) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.078331) q[0];
sx q[0];
rz(-2.445979) q[0];
sx q[0];
rz(0.22380933) q[0];
x q[1];
rz(-1.707294) q[2];
sx q[2];
rz(-2.578306) q[2];
sx q[2];
rz(-2.0493281) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3486466) q[1];
sx q[1];
rz(-2.209084) q[1];
sx q[1];
rz(-1.9116886) q[1];
x q[2];
rz(2.5528615) q[3];
sx q[3];
rz(-1.7541459) q[3];
sx q[3];
rz(-0.4962894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.76379124) q[2];
sx q[2];
rz(-0.8447454) q[2];
sx q[2];
rz(1.2949004) q[2];
rz(-0.14136782) q[3];
sx q[3];
rz(-0.54261345) q[3];
sx q[3];
rz(0.98658371) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7871053) q[0];
sx q[0];
rz(-1.1802477) q[0];
sx q[0];
rz(-1.6217344) q[0];
rz(2.5095818) q[1];
sx q[1];
rz(-2.4193587) q[1];
sx q[1];
rz(-2.246726) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8653523) q[0];
sx q[0];
rz(-1.1603174) q[0];
sx q[0];
rz(-0.68666896) q[0];
rz(-pi) q[1];
rz(1.0084444) q[2];
sx q[2];
rz(-2.214606) q[2];
sx q[2];
rz(1.2959727) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6579073) q[1];
sx q[1];
rz(-0.92805082) q[1];
sx q[1];
rz(-2.583858) q[1];
rz(-2.8551293) q[3];
sx q[3];
rz(-1.5634007) q[3];
sx q[3];
rz(1.3163819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.52508369) q[2];
sx q[2];
rz(-2.775165) q[2];
sx q[2];
rz(1.1425225) q[2];
rz(-2.3948495) q[3];
sx q[3];
rz(-1.2544422) q[3];
sx q[3];
rz(-2.0660627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58364761) q[0];
sx q[0];
rz(-2.8201411) q[0];
sx q[0];
rz(1.3775795) q[0];
rz(0.47239834) q[1];
sx q[1];
rz(-0.51858416) q[1];
sx q[1];
rz(-0.46498743) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0156292) q[0];
sx q[0];
rz(-2.0032126) q[0];
sx q[0];
rz(-2.7727491) q[0];
x q[1];
rz(-1.497252) q[2];
sx q[2];
rz(-1.2300756) q[2];
sx q[2];
rz(-1.3982915) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2496693) q[1];
sx q[1];
rz(-1.6185456) q[1];
sx q[1];
rz(-2.2062917) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.432514) q[3];
sx q[3];
rz(-1.4104098) q[3];
sx q[3];
rz(1.3080314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.49089367) q[2];
sx q[2];
rz(-1.2630254) q[2];
sx q[2];
rz(-2.6965551) q[2];
rz(0.93368357) q[3];
sx q[3];
rz(-1.4327587) q[3];
sx q[3];
rz(-0.26708189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12748195) q[0];
sx q[0];
rz(-1.575379) q[0];
sx q[0];
rz(1.7215464) q[0];
rz(-0.02380112) q[1];
sx q[1];
rz(-0.61444608) q[1];
sx q[1];
rz(0.15596095) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10931817) q[0];
sx q[0];
rz(-0.30800691) q[0];
sx q[0];
rz(-2.5449635) q[0];
rz(2.4620352) q[2];
sx q[2];
rz(-2.5153) q[2];
sx q[2];
rz(-1.2223513) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5739463) q[1];
sx q[1];
rz(-1.5168377) q[1];
sx q[1];
rz(-1.3957363) q[1];
x q[2];
rz(-0.74219269) q[3];
sx q[3];
rz(-0.67897292) q[3];
sx q[3];
rz(1.2364482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0722787) q[2];
sx q[2];
rz(-2.5645655) q[2];
sx q[2];
rz(1.0726661) q[2];
rz(0.32564751) q[3];
sx q[3];
rz(-1.1762534) q[3];
sx q[3];
rz(-1.6252888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1919365) q[0];
sx q[0];
rz(-2.7392116) q[0];
sx q[0];
rz(0.33777133) q[0];
rz(-2.0514964) q[1];
sx q[1];
rz(-2.1665159) q[1];
sx q[1];
rz(-0.24857323) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6876858) q[0];
sx q[0];
rz(-2.0810063) q[0];
sx q[0];
rz(1.8574255) q[0];
rz(-pi) q[1];
rz(0.48034251) q[2];
sx q[2];
rz(-2.0311653) q[2];
sx q[2];
rz(-2.4127221) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2182902) q[1];
sx q[1];
rz(-2.3308838) q[1];
sx q[1];
rz(-0.87358012) q[1];
rz(1.1484654) q[3];
sx q[3];
rz(-2.9557807) q[3];
sx q[3];
rz(1.2687792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2934072) q[2];
sx q[2];
rz(-2.6082787) q[2];
sx q[2];
rz(-3.0548813) q[2];
rz(-0.48197204) q[3];
sx q[3];
rz(-2.635699) q[3];
sx q[3];
rz(-1.1530676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6329353) q[0];
sx q[0];
rz(-1.0077227) q[0];
sx q[0];
rz(-0.28276643) q[0];
rz(0.70156082) q[1];
sx q[1];
rz(-2.3200254) q[1];
sx q[1];
rz(1.3185906) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9160999) q[0];
sx q[0];
rz(-2.6927104) q[0];
sx q[0];
rz(-1.7569957) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.73762383) q[2];
sx q[2];
rz(-2.3361623) q[2];
sx q[2];
rz(1.1951624) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1773771) q[1];
sx q[1];
rz(-2.4417158) q[1];
sx q[1];
rz(1.3436951) q[1];
x q[2];
rz(2.7301844) q[3];
sx q[3];
rz(-0.79380006) q[3];
sx q[3];
rz(1.0994764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7302154) q[2];
sx q[2];
rz(-0.29835478) q[2];
sx q[2];
rz(-2.7424157) q[2];
rz(2.2579851) q[3];
sx q[3];
rz(-1.6688321) q[3];
sx q[3];
rz(1.2214899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76319641) q[0];
sx q[0];
rz(-2.3662687) q[0];
sx q[0];
rz(1.5426853) q[0];
rz(-2.0762766) q[1];
sx q[1];
rz(-0.97384802) q[1];
sx q[1];
rz(-1.8803966) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1051837) q[0];
sx q[0];
rz(-1.1622218) q[0];
sx q[0];
rz(1.6856134) q[0];
rz(1.5485498) q[2];
sx q[2];
rz(-2.1689479) q[2];
sx q[2];
rz(-1.0215789) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.8775455) q[1];
sx q[1];
rz(-1.0730768) q[1];
sx q[1];
rz(1.643814) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.68393771) q[3];
sx q[3];
rz(-1.6947019) q[3];
sx q[3];
rz(-2.370196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5836872) q[2];
sx q[2];
rz(-0.62290278) q[2];
sx q[2];
rz(-0.56979257) q[2];
rz(-1.9231046) q[3];
sx q[3];
rz(-0.64703882) q[3];
sx q[3];
rz(-2.5789554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31310836) q[0];
sx q[0];
rz(-2.2611571) q[0];
sx q[0];
rz(-1.8631998) q[0];
rz(0.60824153) q[1];
sx q[1];
rz(-2.6651762) q[1];
sx q[1];
rz(-2.6574635) q[1];
rz(2.1424978) q[2];
sx q[2];
rz(-1.5487164) q[2];
sx q[2];
rz(-1.6185417) q[2];
rz(-2.5504997) q[3];
sx q[3];
rz(-1.686284) q[3];
sx q[3];
rz(-1.4808663) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
