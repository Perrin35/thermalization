OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg meas[4];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.27999347448349) q[0];
sx q[0];
rz(2.10928288300569) q[0];
sx q[0];
rz(10.5849246740262) q[0];
rz(-0.878148257732391) q[1];
sx q[1];
rz(3.58017635543878) q[1];
sx q[1];
rz(11.68624994754) q[1];
cx q[1],q[0];
rz(3.25346064567566) q[0];
sx q[0];
rz(2.16391852696473) q[0];
sx q[0];
rz(8.06220433711215) q[0];
rz(3.3282630443573) q[2];
sx q[2];
rz(2.73230326374108) q[2];
sx q[2];
rz(9.28433591722652) q[2];
cx q[2],q[1];
rz(1.97012007236481) q[1];
sx q[1];
rz(2.51815411646897) q[1];
sx q[1];
rz(9.70610669850513) q[1];
rz(4.47162055969238) q[3];
sx q[3];
rz(3.60723313887651) q[3];
sx q[3];
rz(9.18962516485854) q[3];
cx q[3],q[2];
rz(4.25674390792847) q[2];
sx q[2];
rz(3.70348885853822) q[2];
sx q[2];
rz(7.40302631854221) q[2];
rz(-1.00370538234711) q[3];
sx q[3];
rz(3.49150154192979) q[3];
sx q[3];
rz(7.13482449053928) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.729646623134613) q[0];
sx q[0];
rz(4.03453335364396) q[0];
sx q[0];
rz(9.87284216880008) q[0];
rz(1.075315117836) q[1];
sx q[1];
rz(5.20653906662995) q[1];
sx q[1];
rz(9.46524263396069) q[1];
cx q[1],q[0];
rz(2.41416358947754) q[0];
sx q[0];
rz(4.05028125842149) q[0];
sx q[0];
rz(10.7519865989606) q[0];
rz(1.66056895256042) q[2];
sx q[2];
rz(4.30055645306642) q[2];
sx q[2];
rz(10.685776090614) q[2];
cx q[2],q[1];
rz(0.654617965221405) q[1];
sx q[1];
rz(3.41346398194367) q[1];
sx q[1];
rz(11.1496046543042) q[1];
rz(-1.54417026042938) q[3];
sx q[3];
rz(3.45431339939172) q[3];
sx q[3];
rz(10.7687941551129) q[3];
cx q[3],q[2];
rz(-0.865410923957825) q[2];
sx q[2];
rz(4.04222575028474) q[2];
sx q[2];
rz(13.069740509979) q[2];
rz(1.6940176486969) q[3];
sx q[3];
rz(5.04990545113618) q[3];
sx q[3];
rz(8.51746747492954) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.685097873210907) q[0];
sx q[0];
rz(5.32958236535127) q[0];
sx q[0];
rz(9.12630877494022) q[0];
rz(1.58627700805664) q[1];
sx q[1];
rz(4.55119803746278) q[1];
sx q[1];
rz(11.0599769115369) q[1];
cx q[1],q[0];
rz(2.85785460472107) q[0];
sx q[0];
rz(1.86877921422059) q[0];
sx q[0];
rz(9.34625703691646) q[0];
rz(-2.45506262779236) q[2];
sx q[2];
rz(7.63747659524018) q[2];
sx q[2];
rz(8.91281358002826) q[2];
cx q[2],q[1];
rz(1.75064444541931) q[1];
sx q[1];
rz(1.43812886078889) q[1];
sx q[1];
rz(8.54116258620425) q[1];
rz(0.275718152523041) q[3];
sx q[3];
rz(3.33840853174264) q[3];
sx q[3];
rz(11.907738184921) q[3];
cx q[3],q[2];
rz(5.93111944198608) q[2];
sx q[2];
rz(5.51469245751435) q[2];
sx q[2];
rz(10.1963937640111) q[2];
rz(1.73028922080994) q[3];
sx q[3];
rz(4.36863461335237) q[3];
sx q[3];
rz(8.59510020016834) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.43784332275391) q[0];
sx q[0];
rz(2.48985359271104) q[0];
sx q[0];
rz(11.2320255994718) q[0];
rz(-1.62885451316833) q[1];
sx q[1];
rz(4.76294234593446) q[1];
sx q[1];
rz(9.28817149101897) q[1];
cx q[1],q[0];
rz(-0.301349818706512) q[0];
sx q[0];
rz(5.57063523133332) q[0];
sx q[0];
rz(9.90571109055682) q[0];
rz(-0.30963122844696) q[2];
sx q[2];
rz(7.75533452828462) q[2];
sx q[2];
rz(11.0721466302793) q[2];
cx q[2],q[1];
rz(0.260306686162949) q[1];
sx q[1];
rz(3.91102388699586) q[1];
sx q[1];
rz(11.4517187833707) q[1];
rz(-1.903932929039) q[3];
sx q[3];
rz(4.38129058678681) q[3];
sx q[3];
rz(9.47319784983202) q[3];
cx q[3],q[2];
rz(2.28873491287231) q[2];
sx q[2];
rz(4.62061956723268) q[2];
sx q[2];
rz(9.1852054208438) q[2];
rz(0.95909708738327) q[3];
sx q[3];
rz(5.11909428437287) q[3];
sx q[3];
rz(8.38687810897037) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.949448227882385) q[0];
sx q[0];
rz(5.22780814965303) q[0];
sx q[0];
rz(10.9685612678449) q[0];
rz(1.1100617647171) q[1];
sx q[1];
rz(5.07972827752168) q[1];
sx q[1];
rz(10.8718405723493) q[1];
cx q[1],q[0];
rz(0.271810680627823) q[0];
sx q[0];
rz(2.84310341079766) q[0];
sx q[0];
rz(11.0198149442594) q[0];
rz(0.343084782361984) q[2];
sx q[2];
rz(1.6091846545511) q[2];
sx q[2];
rz(8.94750536083385) q[2];
cx q[2],q[1];
rz(1.28671085834503) q[1];
sx q[1];
rz(2.44577023585374) q[1];
sx q[1];
rz(10.7204917430799) q[1];
rz(-2.20820569992065) q[3];
sx q[3];
rz(3.78634444077546) q[3];
sx q[3];
rz(11.2022102832715) q[3];
cx q[3],q[2];
rz(0.996282279491425) q[2];
sx q[2];
rz(5.14454594452912) q[2];
sx q[2];
rz(11.7471036672513) q[2];
rz(0.947756886482239) q[3];
sx q[3];
rz(2.35343334277207) q[3];
sx q[3];
rz(7.53906795977756) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.63256967067719) q[0];
sx q[0];
rz(3.2795498838001) q[0];
sx q[0];
rz(9.97487477063342) q[0];
rz(0.247127860784531) q[1];
sx q[1];
rz(4.2942192872339) q[1];
sx q[1];
rz(11.620282149307) q[1];
cx q[1],q[0];
rz(-0.450296342372894) q[0];
sx q[0];
rz(4.35581293900544) q[0];
sx q[0];
rz(8.13086280821964) q[0];
rz(-1.63312005996704) q[2];
sx q[2];
rz(1.78354969819123) q[2];
sx q[2];
rz(10.1733957886617) q[2];
cx q[2],q[1];
rz(0.94080775976181) q[1];
sx q[1];
rz(5.18190494378144) q[1];
sx q[1];
rz(9.46484338342353) q[1];
rz(-1.97842311859131) q[3];
sx q[3];
rz(3.40085420210893) q[3];
sx q[3];
rz(11.3172769308011) q[3];
cx q[3],q[2];
rz(0.755210220813751) q[2];
sx q[2];
rz(5.73275914986665) q[2];
sx q[2];
rz(10.1058404803197) q[2];
rz(2.33034729957581) q[3];
sx q[3];
rz(2.097605021792) q[3];
sx q[3];
rz(12.0160920381467) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(3.16383862495422) q[0];
sx q[0];
rz(2.52110538085038) q[0];
sx q[0];
rz(10.1703897476117) q[0];
rz(1.97818076610565) q[1];
sx q[1];
rz(2.52015760739381) q[1];
sx q[1];
rz(9.24754922687217) q[1];
cx q[1],q[0];
rz(-1.13543725013733) q[0];
sx q[0];
rz(3.57597112854058) q[0];
sx q[0];
rz(10.114000594608) q[0];
rz(0.0907635763287544) q[2];
sx q[2];
rz(4.77541199524934) q[2];
sx q[2];
rz(11.1677739381711) q[2];
cx q[2],q[1];
rz(-1.19113504886627) q[1];
sx q[1];
rz(5.21755019028718) q[1];
sx q[1];
rz(11.4162302970807) q[1];
rz(0.04217429459095) q[3];
sx q[3];
rz(2.26510533888871) q[3];
sx q[3];
rz(8.66943947076007) q[3];
cx q[3],q[2];
rz(3.69156551361084) q[2];
sx q[2];
rz(0.982293041544505) q[2];
sx q[2];
rz(4.59091613291904) q[2];
rz(-2.51894044876099) q[3];
sx q[3];
rz(3.66875788767869) q[3];
sx q[3];
rz(11.2469413041989) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(3.15826153755188) q[0];
sx q[0];
rz(1.3873839696222) q[0];
sx q[0];
rz(11.054383134834) q[0];
rz(-2.21770572662354) q[1];
sx q[1];
rz(1.72165635426576) q[1];
sx q[1];
rz(13.1830088853757) q[1];
cx q[1],q[0];
rz(1.44093263149261) q[0];
sx q[0];
rz(2.1671138723665) q[0];
sx q[0];
rz(10.8112787961881) q[0];
rz(-1.65710997581482) q[2];
sx q[2];
rz(3.61097181041772) q[2];
sx q[2];
rz(11.5769886732022) q[2];
cx q[2],q[1];
rz(-0.603959977626801) q[1];
sx q[1];
rz(5.65808978875215) q[1];
sx q[1];
rz(7.70733556746646) q[1];
rz(0.988790273666382) q[3];
sx q[3];
rz(4.4811419566446) q[3];
sx q[3];
rz(10.4207116722982) q[3];
cx q[3],q[2];
rz(-1.77432489395142) q[2];
sx q[2];
rz(3.50109082658822) q[2];
sx q[2];
rz(15.873022055618) q[2];
rz(2.39556908607483) q[3];
sx q[3];
rz(4.76438549359376) q[3];
sx q[3];
rz(9.46757481469914) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.541667222976685) q[0];
sx q[0];
rz(3.4337034543329) q[0];
sx q[0];
rz(9.01065824031039) q[0];
rz(-0.435530215501785) q[1];
sx q[1];
rz(2.62565168936784) q[1];
sx q[1];
rz(8.14645025729343) q[1];
cx q[1],q[0];
rz(1.17624866962433) q[0];
sx q[0];
rz(2.85817355115945) q[0];
sx q[0];
rz(6.83101699351474) q[0];
rz(0.173689022660255) q[2];
sx q[2];
rz(-1.26354822318023) q[2];
sx q[2];
rz(8.23424670695468) q[2];
cx q[2],q[1];
rz(3.46291518211365) q[1];
sx q[1];
rz(1.50291171868379) q[1];
sx q[1];
rz(9.41289550381854) q[1];
rz(-0.107563078403473) q[3];
sx q[3];
rz(3.67195835907991) q[3];
sx q[3];
rz(7.70430383681461) q[3];
cx q[3],q[2];
rz(0.723296761512756) q[2];
sx q[2];
rz(4.52386108239228) q[2];
sx q[2];
rz(10.1674930810849) q[2];
rz(-0.791394770145416) q[3];
sx q[3];
rz(5.59959641297395) q[3];
sx q[3];
rz(8.01556763648196) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.0475697107613087) q[0];
sx q[0];
rz(3.47091195185716) q[0];
sx q[0];
rz(10.336399114124) q[0];
rz(0.985126376152039) q[1];
sx q[1];
rz(3.56834277709062) q[1];
sx q[1];
rz(11.2282358169477) q[1];
cx q[1],q[0];
rz(2.01659226417542) q[0];
sx q[0];
rz(5.69906249840791) q[0];
sx q[0];
rz(10.0183115959088) q[0];
rz(5.02974796295166) q[2];
sx q[2];
rz(4.35090902646119) q[2];
sx q[2];
rz(8.25163004397556) q[2];
cx q[2],q[1];
rz(0.123051889240742) q[1];
sx q[1];
rz(4.31267789204652) q[1];
sx q[1];
rz(10.0786877036016) q[1];
rz(0.373644143342972) q[3];
sx q[3];
rz(1.2579529603296) q[3];
sx q[3];
rz(10.3065403461377) q[3];
cx q[3],q[2];
rz(0.354758739471436) q[2];
sx q[2];
rz(4.6165122111612) q[2];
sx q[2];
rz(13.7317104101102) q[2];
rz(1.3960086107254) q[3];
sx q[3];
rz(1.95201197464997) q[3];
sx q[3];
rz(7.87499723433658) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.0521978251636028) q[0];
sx q[0];
rz(5.45101371605928) q[0];
sx q[0];
rz(9.44609055704578) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(0.965688049793243) q[1];
sx q[1];
rz(0.808651359873362) q[1];
sx q[1];
rz(9.47150011210843) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(-1.85960257053375) q[2];
sx q[2];
rz(5.49940720398957) q[2];
sx q[2];
rz(12.7558197736661) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(-1.6839154958725) q[3];
sx q[3];
rz(4.1184048970514) q[3];
sx q[3];
rz(13.8096427678983) q[3];
rz(pi/2) q[3];
sx q[3];
rz(3*pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> meas[0];
measure q[1] -> meas[1];
measure q[2] -> meas[2];
measure q[3] -> meas[3];
