OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.1228696) q[0];
sx q[0];
rz(-2.2824204) q[0];
sx q[0];
rz(0.81508842) q[0];
rz(1.9563142) q[1];
sx q[1];
rz(4.5524608) q[1];
sx q[1];
rz(8.3571385) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0420494) q[0];
sx q[0];
rz(-1.1394412) q[0];
sx q[0];
rz(-3.1154446) q[0];
rz(-pi) q[1];
x q[1];
rz(0.52402988) q[2];
sx q[2];
rz(-1.9240148) q[2];
sx q[2];
rz(-0.86862446) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.63032615) q[1];
sx q[1];
rz(-2.8095803) q[1];
sx q[1];
rz(1.270833) q[1];
x q[2];
rz(2.2566363) q[3];
sx q[3];
rz(-1.7772632) q[3];
sx q[3];
rz(-1.5625169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.69433576) q[2];
sx q[2];
rz(-0.7434291) q[2];
sx q[2];
rz(-1.2747964) q[2];
rz(-2.6796807) q[3];
sx q[3];
rz(-0.67449823) q[3];
sx q[3];
rz(-2.0089202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79134113) q[0];
sx q[0];
rz(-0.30838648) q[0];
sx q[0];
rz(1.2530918) q[0];
rz(-2.99627) q[1];
sx q[1];
rz(-1.7456313) q[1];
sx q[1];
rz(-1.0911509) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4900794) q[0];
sx q[0];
rz(-1.7848178) q[0];
sx q[0];
rz(-1.1351311) q[0];
rz(2.4895489) q[2];
sx q[2];
rz(-2.5863159) q[2];
sx q[2];
rz(-2.9342143) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5538283) q[1];
sx q[1];
rz(-1.3979646) q[1];
sx q[1];
rz(0.21765222) q[1];
rz(-0.66295538) q[3];
sx q[3];
rz(-2.5409449) q[3];
sx q[3];
rz(-1.8523703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.48065177) q[2];
sx q[2];
rz(-1.7512243) q[2];
sx q[2];
rz(2.6118028) q[2];
rz(-2.3482813) q[3];
sx q[3];
rz(-1.6042234) q[3];
sx q[3];
rz(-2.5180499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48524258) q[0];
sx q[0];
rz(-2.1915477) q[0];
sx q[0];
rz(2.6128838) q[0];
rz(2.5953925) q[1];
sx q[1];
rz(-2.1827953) q[1];
sx q[1];
rz(2.8012457) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2961194) q[0];
sx q[0];
rz(-1.8438135) q[0];
sx q[0];
rz(-0.34015981) q[0];
rz(-pi) q[1];
rz(-2.9067847) q[2];
sx q[2];
rz(-1.4291414) q[2];
sx q[2];
rz(1.0740785) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7527377) q[1];
sx q[1];
rz(-1.5740663) q[1];
sx q[1];
rz(-3.0470303) q[1];
rz(-pi) q[2];
rz(0.82967088) q[3];
sx q[3];
rz(-1.8541012) q[3];
sx q[3];
rz(2.0633351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9693552) q[2];
sx q[2];
rz(-2.7293971) q[2];
sx q[2];
rz(-2.7117512) q[2];
rz(2.3853081) q[3];
sx q[3];
rz(-2.9679306) q[3];
sx q[3];
rz(-2.3156796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38465685) q[0];
sx q[0];
rz(-1.5058368) q[0];
sx q[0];
rz(-1.4935619) q[0];
rz(-2.3736296) q[1];
sx q[1];
rz(-2.6710644) q[1];
sx q[1];
rz(-2.5951662) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38619216) q[0];
sx q[0];
rz(-1.6419171) q[0];
sx q[0];
rz(3.0767308) q[0];
rz(-2.3713263) q[2];
sx q[2];
rz(-0.67521836) q[2];
sx q[2];
rz(-0.16083052) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5560234) q[1];
sx q[1];
rz(-1.6680191) q[1];
sx q[1];
rz(2.9929964) q[1];
x q[2];
rz(-0.91937842) q[3];
sx q[3];
rz(-0.99777824) q[3];
sx q[3];
rz(1.023055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4297318) q[2];
sx q[2];
rz(-1.6542566) q[2];
sx q[2];
rz(0.30409733) q[2];
rz(1.3652623) q[3];
sx q[3];
rz(-1.2208166) q[3];
sx q[3];
rz(1.0767267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6292608) q[0];
sx q[0];
rz(-2.6973695) q[0];
sx q[0];
rz(0.00057922676) q[0];
rz(-2.0594788) q[1];
sx q[1];
rz(-0.52434701) q[1];
sx q[1];
rz(-2.2023315) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2596596) q[0];
sx q[0];
rz(-1.1355917) q[0];
sx q[0];
rz(-2.4324904) q[0];
rz(-pi) q[1];
x q[1];
rz(0.3498431) q[2];
sx q[2];
rz(-0.98825422) q[2];
sx q[2];
rz(-2.2900641) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4537332) q[1];
sx q[1];
rz(-2.5304768) q[1];
sx q[1];
rz(1.9495717) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1904508) q[3];
sx q[3];
rz(-2.7627146) q[3];
sx q[3];
rz(-0.71228107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.038736343) q[2];
sx q[2];
rz(-1.880371) q[2];
sx q[2];
rz(-0.2612513) q[2];
rz(0.86483613) q[3];
sx q[3];
rz(-1.7295001) q[3];
sx q[3];
rz(0.29204667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84457266) q[0];
sx q[0];
rz(-1.7140056) q[0];
sx q[0];
rz(2.936506) q[0];
rz(-2.0948441) q[1];
sx q[1];
rz(-1.3958684) q[1];
sx q[1];
rz(-2.7632025) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16941026) q[0];
sx q[0];
rz(-0.43146389) q[0];
sx q[0];
rz(-1.1110825) q[0];
rz(-pi) q[1];
rz(2.8257337) q[2];
sx q[2];
rz(-2.6489885) q[2];
sx q[2];
rz(1.619018) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.3832964) q[1];
sx q[1];
rz(-1.2219011) q[1];
sx q[1];
rz(-2.1339244) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6343253) q[3];
sx q[3];
rz(-2.3823822) q[3];
sx q[3];
rz(-1.678148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.64688524) q[2];
sx q[2];
rz(-1.1581706) q[2];
sx q[2];
rz(1.0542487) q[2];
rz(-0.65822893) q[3];
sx q[3];
rz(-2.4963278) q[3];
sx q[3];
rz(0.48404199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9996416) q[0];
sx q[0];
rz(-1.208409) q[0];
sx q[0];
rz(0.64055881) q[0];
rz(-0.41796747) q[1];
sx q[1];
rz(-1.2203981) q[1];
sx q[1];
rz(-0.80498615) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98950878) q[0];
sx q[0];
rz(-1.594127) q[0];
sx q[0];
rz(-0.70744608) q[0];
rz(-pi) q[1];
rz(0.4067602) q[2];
sx q[2];
rz(-0.8443588) q[2];
sx q[2];
rz(-0.11876362) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0143483) q[1];
sx q[1];
rz(-1.5906723) q[1];
sx q[1];
rz(1.5991421) q[1];
x q[2];
rz(2.715519) q[3];
sx q[3];
rz(-0.63562993) q[3];
sx q[3];
rz(-0.15492188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.54887041) q[2];
sx q[2];
rz(-0.9950811) q[2];
sx q[2];
rz(-0.76118809) q[2];
rz(2.393764) q[3];
sx q[3];
rz(-1.8239832) q[3];
sx q[3];
rz(-0.67659155) q[3];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51650301) q[0];
sx q[0];
rz(-2.857132) q[0];
sx q[0];
rz(1.4768584) q[0];
rz(2.6853216) q[1];
sx q[1];
rz(-1.7218593) q[1];
sx q[1];
rz(0.87108535) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81646279) q[0];
sx q[0];
rz(-1.5832381) q[0];
sx q[0];
rz(-2.5274656) q[0];
rz(-pi) q[1];
rz(1.4159059) q[2];
sx q[2];
rz(-1.911507) q[2];
sx q[2];
rz(-2.129385) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5567315) q[1];
sx q[1];
rz(-1.24354) q[1];
sx q[1];
rz(-1.3393988) q[1];
rz(3.0814287) q[3];
sx q[3];
rz(-1.9964661) q[3];
sx q[3];
rz(-0.7983467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2386834) q[2];
sx q[2];
rz(-2.6148655) q[2];
sx q[2];
rz(-2.5229559) q[2];
rz(3.1213308) q[3];
sx q[3];
rz(-2.198115) q[3];
sx q[3];
rz(2.3616135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.063754931) q[0];
sx q[0];
rz(-1.5851333) q[0];
sx q[0];
rz(-0.42386398) q[0];
rz(2.9546812) q[1];
sx q[1];
rz(-0.82004768) q[1];
sx q[1];
rz(-1.4580457) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8992726) q[0];
sx q[0];
rz(-1.5849216) q[0];
sx q[0];
rz(-0.11798162) q[0];
rz(-pi) q[1];
rz(-2.0293268) q[2];
sx q[2];
rz(-2.5800309) q[2];
sx q[2];
rz(2.4379345) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2413689) q[1];
sx q[1];
rz(-1.7094862) q[1];
sx q[1];
rz(1.5189511) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.35590509) q[3];
sx q[3];
rz(-1.6545466) q[3];
sx q[3];
rz(0.29415302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.75108782) q[2];
sx q[2];
rz(-2.0743399) q[2];
sx q[2];
rz(2.0474153) q[2];
rz(-1.68082) q[3];
sx q[3];
rz(-1.7655617) q[3];
sx q[3];
rz(-1.8028397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8193034) q[0];
sx q[0];
rz(-1.3811454) q[0];
sx q[0];
rz(-1.5018916) q[0];
rz(0.27885258) q[1];
sx q[1];
rz(-2.1809705) q[1];
sx q[1];
rz(-2.1174812) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0251966) q[0];
sx q[0];
rz(-1.0828583) q[0];
sx q[0];
rz(1.4665718) q[0];
rz(-pi) q[1];
rz(1.845593) q[2];
sx q[2];
rz(-1.6195903) q[2];
sx q[2];
rz(-2.9279207) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1361724) q[1];
sx q[1];
rz(-1.3498303) q[1];
sx q[1];
rz(1.3653838) q[1];
rz(-pi) q[2];
rz(1.9387705) q[3];
sx q[3];
rz(-2.362613) q[3];
sx q[3];
rz(0.72768962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9378822) q[2];
sx q[2];
rz(-1.2869765) q[2];
sx q[2];
rz(0.53696519) q[2];
rz(1.9133866) q[3];
sx q[3];
rz(-0.95913404) q[3];
sx q[3];
rz(-0.29135191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0252329) q[0];
sx q[0];
rz(-1.3675084) q[0];
sx q[0];
rz(1.0687923) q[0];
rz(2.4346726) q[1];
sx q[1];
rz(-1.128935) q[1];
sx q[1];
rz(3.1405906) q[1];
rz(-0.42264414) q[2];
sx q[2];
rz(-0.44036897) q[2];
sx q[2];
rz(-1.8457495) q[2];
rz(-0.14866004) q[3];
sx q[3];
rz(-1.3547263) q[3];
sx q[3];
rz(2.3480036) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
