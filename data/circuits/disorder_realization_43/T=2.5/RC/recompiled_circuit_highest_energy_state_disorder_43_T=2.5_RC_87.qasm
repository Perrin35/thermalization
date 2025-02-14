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
rz(-3.0946331) q[0];
sx q[0];
rz(-1.5927915) q[0];
sx q[0];
rz(-1.3943075) q[0];
rz(2.6810763) q[1];
sx q[1];
rz(-1.2682275) q[1];
sx q[1];
rz(0.4761129) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5232485) q[0];
sx q[0];
rz(-1.9501513) q[0];
sx q[0];
rz(0.98486395) q[0];
x q[1];
rz(0.85501315) q[2];
sx q[2];
rz(-2.7913793) q[2];
sx q[2];
rz(2.8976655) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.81836593) q[1];
sx q[1];
rz(-1.4225679) q[1];
sx q[1];
rz(1.5532975) q[1];
rz(-pi) q[2];
rz(1.634811) q[3];
sx q[3];
rz(-1.598017) q[3];
sx q[3];
rz(-0.79827362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8146879) q[2];
sx q[2];
rz(-1.9196332) q[2];
sx q[2];
rz(-1.8587221) q[2];
rz(3.0021216) q[3];
sx q[3];
rz(-1.8255511) q[3];
sx q[3];
rz(0.68525806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71444702) q[0];
sx q[0];
rz(-2.7834263) q[0];
sx q[0];
rz(2.8156679) q[0];
rz(-0.70949316) q[1];
sx q[1];
rz(-2.8804417) q[1];
sx q[1];
rz(-0.66169468) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70764953) q[0];
sx q[0];
rz(-0.012059472) q[0];
sx q[0];
rz(-0.90977408) q[0];
x q[1];
rz(-0.9307894) q[2];
sx q[2];
rz(-2.402596) q[2];
sx q[2];
rz(1.3698824) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0256309) q[1];
sx q[1];
rz(-0.44646663) q[1];
sx q[1];
rz(2.129351) q[1];
rz(-pi) q[2];
rz(0.98425166) q[3];
sx q[3];
rz(-0.97347608) q[3];
sx q[3];
rz(-0.95109361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0207396) q[2];
sx q[2];
rz(-1.8031305) q[2];
sx q[2];
rz(-2.8815114) q[2];
rz(0.86380473) q[3];
sx q[3];
rz(-1.6983906) q[3];
sx q[3];
rz(-1.1687733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23673713) q[0];
sx q[0];
rz(-0.78176347) q[0];
sx q[0];
rz(-0.28808638) q[0];
rz(-1.2068564) q[1];
sx q[1];
rz(-1.1129881) q[1];
sx q[1];
rz(-2.3634214) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.208565) q[0];
sx q[0];
rz(-1.3817245) q[0];
sx q[0];
rz(-2.2651432) q[0];
rz(-pi) q[1];
rz(1.6323998) q[2];
sx q[2];
rz(-2.3280316) q[2];
sx q[2];
rz(0.98682994) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.84904387) q[1];
sx q[1];
rz(-2.0660225) q[1];
sx q[1];
rz(-1.0775113) q[1];
rz(-pi) q[2];
rz(0.92375375) q[3];
sx q[3];
rz(-1.2525932) q[3];
sx q[3];
rz(1.3830087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.244016) q[2];
sx q[2];
rz(-2.0147169) q[2];
sx q[2];
rz(0.24813949) q[2];
rz(-2.1085619) q[3];
sx q[3];
rz(-1.4890198) q[3];
sx q[3];
rz(-1.8831133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7018062) q[0];
sx q[0];
rz(-2.2070364) q[0];
sx q[0];
rz(-2.5820861) q[0];
rz(1.4065546) q[1];
sx q[1];
rz(-2.1969257) q[1];
sx q[1];
rz(1.14538) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1702207) q[0];
sx q[0];
rz(-1.6165501) q[0];
sx q[0];
rz(-2.6849823) q[0];
x q[1];
rz(-2.4539646) q[2];
sx q[2];
rz(-2.3634644) q[2];
sx q[2];
rz(-1.8541921) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8662474) q[1];
sx q[1];
rz(-1.8653578) q[1];
sx q[1];
rz(-0.7742851) q[1];
rz(1.000993) q[3];
sx q[3];
rz(-1.5581597) q[3];
sx q[3];
rz(-2.3848052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.8765325) q[2];
sx q[2];
rz(-1.0141076) q[2];
sx q[2];
rz(1.8939135) q[2];
rz(1.0809336) q[3];
sx q[3];
rz(-1.388988) q[3];
sx q[3];
rz(1.2601674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0718229) q[0];
sx q[0];
rz(-0.48957303) q[0];
sx q[0];
rz(0.74482942) q[0];
rz(-1.3818332) q[1];
sx q[1];
rz(-1.1188743) q[1];
sx q[1];
rz(-0.098085731) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9811889) q[0];
sx q[0];
rz(-1.2645928) q[0];
sx q[0];
rz(2.7723377) q[0];
rz(-pi) q[1];
rz(0.44386835) q[2];
sx q[2];
rz(-1.265268) q[2];
sx q[2];
rz(-1.533184) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.67691117) q[1];
sx q[1];
rz(-1.5524852) q[1];
sx q[1];
rz(0.69377884) q[1];
rz(1.7616055) q[3];
sx q[3];
rz(-2.6529783) q[3];
sx q[3];
rz(0.055475108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.44094008) q[2];
sx q[2];
rz(-0.38267371) q[2];
sx q[2];
rz(1.0901701) q[2];
rz(-2.8160461) q[3];
sx q[3];
rz(-1.5109589) q[3];
sx q[3];
rz(-0.24698273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22661041) q[0];
sx q[0];
rz(-2.60422) q[0];
sx q[0];
rz(-1.2363303) q[0];
rz(0.63713282) q[1];
sx q[1];
rz(-2.5814711) q[1];
sx q[1];
rz(0.93322745) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97605825) q[0];
sx q[0];
rz(-1.2090431) q[0];
sx q[0];
rz(-2.1657337) q[0];
x q[1];
rz(-1.1514329) q[2];
sx q[2];
rz(-1.3977111) q[2];
sx q[2];
rz(-1.8669548) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5781587) q[1];
sx q[1];
rz(-1.8616042) q[1];
sx q[1];
rz(-1.6573424) q[1];
rz(-pi) q[2];
rz(2.0863462) q[3];
sx q[3];
rz(-1.7988867) q[3];
sx q[3];
rz(2.4683964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8895662) q[2];
sx q[2];
rz(-1.7539975) q[2];
sx q[2];
rz(-1.265556) q[2];
rz(-1.7172074) q[3];
sx q[3];
rz(-2.6129183) q[3];
sx q[3];
rz(2.2831447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3739361) q[0];
sx q[0];
rz(-2.7601384) q[0];
sx q[0];
rz(-2.9071627) q[0];
rz(-0.43342057) q[1];
sx q[1];
rz(-1.0973009) q[1];
sx q[1];
rz(-2.0030599) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94110452) q[0];
sx q[0];
rz(-0.88747665) q[0];
sx q[0];
rz(0.79534689) q[0];
x q[1];
rz(3.1267371) q[2];
sx q[2];
rz(-2.6508923) q[2];
sx q[2];
rz(2.4717836) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5229021) q[1];
sx q[1];
rz(-1.5444718) q[1];
sx q[1];
rz(2.0731889) q[1];
x q[2];
rz(-0.68692211) q[3];
sx q[3];
rz(-1.7217155) q[3];
sx q[3];
rz(0.17450813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.66702691) q[2];
sx q[2];
rz(-2.462025) q[2];
sx q[2];
rz(0.95728528) q[2];
rz(-2.3840617) q[3];
sx q[3];
rz(-1.4742955) q[3];
sx q[3];
rz(0.1434513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2738508) q[0];
sx q[0];
rz(-1.1026646) q[0];
sx q[0];
rz(-2.962501) q[0];
rz(2.9415019) q[1];
sx q[1];
rz(-2.4724019) q[1];
sx q[1];
rz(-0.77654138) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0514798) q[0];
sx q[0];
rz(-0.5585554) q[0];
sx q[0];
rz(1.9165975) q[0];
rz(-pi) q[1];
rz(0.29751038) q[2];
sx q[2];
rz(-2.1665467) q[2];
sx q[2];
rz(1.7864625) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.071826) q[1];
sx q[1];
rz(-0.36301314) q[1];
sx q[1];
rz(-2.9337204) q[1];
x q[2];
rz(-1.4158049) q[3];
sx q[3];
rz(-2.8202882) q[3];
sx q[3];
rz(2.6329272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.17579235) q[2];
sx q[2];
rz(-2.0696023) q[2];
sx q[2];
rz(-1.779186) q[2];
rz(-1.3029178) q[3];
sx q[3];
rz(-1.6630273) q[3];
sx q[3];
rz(1.038507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
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
rz(0.68504828) q[0];
sx q[0];
rz(-1.1948723) q[0];
sx q[0];
rz(-0.1142647) q[0];
rz(1.9396797) q[1];
sx q[1];
rz(-2.6004531) q[1];
sx q[1];
rz(3.0466383) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1916434) q[0];
sx q[0];
rz(-1.3943293) q[0];
sx q[0];
rz(-2.1032189) q[0];
x q[1];
rz(2.0906779) q[2];
sx q[2];
rz(-0.78453817) q[2];
sx q[2];
rz(0.21518165) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.488738) q[1];
sx q[1];
rz(-1.6269557) q[1];
sx q[1];
rz(-0.55027669) q[1];
rz(-pi) q[2];
rz(0.90387266) q[3];
sx q[3];
rz(-0.61792497) q[3];
sx q[3];
rz(2.7300837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8408884) q[2];
sx q[2];
rz(-1.7965094) q[2];
sx q[2];
rz(2.4697206) q[2];
rz(0.68615174) q[3];
sx q[3];
rz(-2.4785564) q[3];
sx q[3];
rz(-1.2175951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1662064) q[0];
sx q[0];
rz(-0.19151846) q[0];
sx q[0];
rz(-1.6084877) q[0];
rz(-2.8168822) q[1];
sx q[1];
rz(-1.8638116) q[1];
sx q[1];
rz(-0.3130354) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7164118) q[0];
sx q[0];
rz(-1.6508174) q[0];
sx q[0];
rz(-2.1437313) q[0];
rz(-pi) q[1];
rz(3.0803142) q[2];
sx q[2];
rz(-1.3512502) q[2];
sx q[2];
rz(0.82876182) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.5165007) q[1];
sx q[1];
rz(-0.91462574) q[1];
sx q[1];
rz(2.28338) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5418986) q[3];
sx q[3];
rz(-2.4238677) q[3];
sx q[3];
rz(-2.5120171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2031871) q[2];
sx q[2];
rz(-1.9878191) q[2];
sx q[2];
rz(-1.3333092) q[2];
rz(2.5238254) q[3];
sx q[3];
rz(-2.1963162) q[3];
sx q[3];
rz(1.1355404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7137322) q[0];
sx q[0];
rz(-1.328631) q[0];
sx q[0];
rz(-2.791639) q[0];
rz(-0.89823828) q[1];
sx q[1];
rz(-1.0844834) q[1];
sx q[1];
rz(-0.35379298) q[1];
rz(-2.4785715) q[2];
sx q[2];
rz(-0.86227476) q[2];
sx q[2];
rz(-3.0515565) q[2];
rz(-1.5946424) q[3];
sx q[3];
rz(-2.3004354) q[3];
sx q[3];
rz(-2.7247747) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
