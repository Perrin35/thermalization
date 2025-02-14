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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61834419) q[0];
sx q[0];
rz(-1.1914413) q[0];
sx q[0];
rz(-0.98486395) q[0];
rz(-0.23525518) q[2];
sx q[2];
rz(-1.8326743) q[2];
sx q[2];
rz(-2.1507598) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3232267) q[1];
sx q[1];
rz(-1.7190248) q[1];
sx q[1];
rz(-1.5532975) q[1];
rz(-pi) q[2];
rz(1.5067817) q[3];
sx q[3];
rz(-1.5435757) q[3];
sx q[3];
rz(-0.79827362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.32690471) q[2];
sx q[2];
rz(-1.9196332) q[2];
sx q[2];
rz(-1.2828705) q[2];
rz(3.0021216) q[3];
sx q[3];
rz(-1.3160416) q[3];
sx q[3];
rz(-0.68525806) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4271456) q[0];
sx q[0];
rz(-0.35816631) q[0];
sx q[0];
rz(-0.32592475) q[0];
rz(-2.4320995) q[1];
sx q[1];
rz(-2.8804417) q[1];
sx q[1];
rz(0.66169468) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.046592043) q[0];
sx q[0];
rz(-1.5612771) q[0];
sx q[0];
rz(0.0074038238) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.49837947) q[2];
sx q[2];
rz(-1.0000668) q[2];
sx q[2];
rz(0.98266593) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0831415) q[1];
sx q[1];
rz(-1.8016691) q[1];
sx q[1];
rz(1.1851694) q[1];
rz(-2.4567408) q[3];
sx q[3];
rz(-2.0461444) q[3];
sx q[3];
rz(-0.97739391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1208531) q[2];
sx q[2];
rz(-1.8031305) q[2];
sx q[2];
rz(-2.8815114) q[2];
rz(2.2777879) q[3];
sx q[3];
rz(-1.443202) q[3];
sx q[3];
rz(1.9728194) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23673713) q[0];
sx q[0];
rz(-2.3598292) q[0];
sx q[0];
rz(-0.28808638) q[0];
rz(-1.2068564) q[1];
sx q[1];
rz(-1.1129881) q[1];
sx q[1];
rz(-2.3634214) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2818277) q[0];
sx q[0];
rz(-2.4261103) q[0];
sx q[0];
rz(-1.2802109) q[0];
rz(-pi) q[1];
rz(0.065041754) q[2];
sx q[2];
rz(-0.75924004) q[2];
sx q[2];
rz(-2.0652079) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.4716123) q[1];
sx q[1];
rz(-2.0005352) q[1];
sx q[1];
rz(2.5915036) q[1];
rz(-pi) q[2];
rz(0.92375375) q[3];
sx q[3];
rz(-1.8889995) q[3];
sx q[3];
rz(1.7585839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8975767) q[2];
sx q[2];
rz(-2.0147169) q[2];
sx q[2];
rz(-0.24813949) q[2];
rz(-1.0330307) q[3];
sx q[3];
rz(-1.6525729) q[3];
sx q[3];
rz(1.2584794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4397864) q[0];
sx q[0];
rz(-2.2070364) q[0];
sx q[0];
rz(2.5820861) q[0];
rz(1.4065546) q[1];
sx q[1];
rz(-0.94466698) q[1];
sx q[1];
rz(1.9962126) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9713719) q[0];
sx q[0];
rz(-1.6165501) q[0];
sx q[0];
rz(0.45661033) q[0];
rz(-pi) q[1];
rz(2.4907095) q[2];
sx q[2];
rz(-2.0325629) q[2];
sx q[2];
rz(-0.24591638) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.0068906764) q[1];
sx q[1];
rz(-2.3242053) q[1];
sx q[1];
rz(0.40938739) q[1];
x q[2];
rz(1.5473751) q[3];
sx q[3];
rz(-2.5716647) q[3];
sx q[3];
rz(-2.3078634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2650602) q[2];
sx q[2];
rz(-1.0141076) q[2];
sx q[2];
rz(1.2476791) q[2];
rz(-1.0809336) q[3];
sx q[3];
rz(-1.7526046) q[3];
sx q[3];
rz(-1.8814253) q[3];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0718229) q[0];
sx q[0];
rz(-2.6520196) q[0];
sx q[0];
rz(0.74482942) q[0];
rz(-1.7597594) q[1];
sx q[1];
rz(-2.0227183) q[1];
sx q[1];
rz(-0.098085731) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25135219) q[0];
sx q[0];
rz(-2.6664263) q[0];
sx q[0];
rz(0.71936468) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9068001) q[2];
sx q[2];
rz(-1.9927551) q[2];
sx q[2];
rz(0.10445751) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4646815) q[1];
sx q[1];
rz(-1.5891074) q[1];
sx q[1];
rz(-0.69377884) q[1];
rz(-pi) q[2];
rz(1.0897343) q[3];
sx q[3];
rz(-1.4816545) q[3];
sx q[3];
rz(-1.3463904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.44094008) q[2];
sx q[2];
rz(-2.7589189) q[2];
sx q[2];
rz(-2.0514226) q[2];
rz(0.32554659) q[3];
sx q[3];
rz(-1.6306337) q[3];
sx q[3];
rz(-2.8946099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22661041) q[0];
sx q[0];
rz(-0.53737265) q[0];
sx q[0];
rz(-1.9052624) q[0];
rz(2.5044598) q[1];
sx q[1];
rz(-2.5814711) q[1];
sx q[1];
rz(-0.93322745) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35966152) q[0];
sx q[0];
rz(-1.0190411) q[0];
sx q[0];
rz(-2.7130068) q[0];
rz(1.1514329) q[2];
sx q[2];
rz(-1.7438816) q[2];
sx q[2];
rz(-1.8669548) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2843202) q[1];
sx q[1];
rz(-0.30306268) q[1];
sx q[1];
rz(2.8604355) q[1];
rz(2.0108443) q[3];
sx q[3];
rz(-2.5820316) q[3];
sx q[3];
rz(-2.6236629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8895662) q[2];
sx q[2];
rz(-1.7539975) q[2];
sx q[2];
rz(-1.265556) q[2];
rz(-1.7172074) q[3];
sx q[3];
rz(-0.5286743) q[3];
sx q[3];
rz(0.85844794) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3739361) q[0];
sx q[0];
rz(-0.38145426) q[0];
sx q[0];
rz(-0.23442991) q[0];
rz(2.7081721) q[1];
sx q[1];
rz(-2.0442918) q[1];
sx q[1];
rz(-1.1385328) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2018801) q[0];
sx q[0];
rz(-2.1576886) q[0];
sx q[0];
rz(-0.71016772) q[0];
rz(3.1267371) q[2];
sx q[2];
rz(-2.6508923) q[2];
sx q[2];
rz(-0.66980904) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2039472) q[1];
sx q[1];
rz(-2.0729985) q[1];
sx q[1];
rz(0.03003386) q[1];
x q[2];
rz(-1.3765923) q[3];
sx q[3];
rz(-2.248431) q[3];
sx q[3];
rz(-1.5189796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4745657) q[2];
sx q[2];
rz(-2.462025) q[2];
sx q[2];
rz(-2.1843074) q[2];
rz(0.75753093) q[3];
sx q[3];
rz(-1.4742955) q[3];
sx q[3];
rz(-2.9981414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8677419) q[0];
sx q[0];
rz(-2.0389281) q[0];
sx q[0];
rz(-0.17909166) q[0];
rz(-2.9415019) q[1];
sx q[1];
rz(-0.6691907) q[1];
sx q[1];
rz(2.3650513) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9574253) q[0];
sx q[0];
rz(-1.7514075) q[0];
sx q[0];
rz(2.1023048) q[0];
x q[1];
rz(0.29751038) q[2];
sx q[2];
rz(-2.1665467) q[2];
sx q[2];
rz(-1.3551301) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4458865) q[1];
sx q[1];
rz(-1.6441455) q[1];
sx q[1];
rz(-2.785745) q[1];
x q[2];
rz(1.8885047) q[3];
sx q[3];
rz(-1.5220257) q[3];
sx q[3];
rz(-1.2093076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9658003) q[2];
sx q[2];
rz(-1.0719904) q[2];
sx q[2];
rz(1.3624066) q[2];
rz(-1.8386748) q[3];
sx q[3];
rz(-1.6630273) q[3];
sx q[3];
rz(2.1030857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68504828) q[0];
sx q[0];
rz(-1.1948723) q[0];
sx q[0];
rz(3.0273279) q[0];
rz(1.201913) q[1];
sx q[1];
rz(-2.6004531) q[1];
sx q[1];
rz(-3.0466383) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1916434) q[0];
sx q[0];
rz(-1.3943293) q[0];
sx q[0];
rz(-1.0383738) q[0];
rz(0.46038119) q[2];
sx q[2];
rz(-2.2308439) q[2];
sx q[2];
rz(-0.4649064) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6528546) q[1];
sx q[1];
rz(-1.514637) q[1];
sx q[1];
rz(-2.591316) q[1];
rz(-pi) q[2];
x q[2];
rz(0.41422959) q[3];
sx q[3];
rz(-2.043402) q[3];
sx q[3];
rz(-0.35652682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.3007043) q[2];
sx q[2];
rz(-1.3450832) q[2];
sx q[2];
rz(-2.4697206) q[2];
rz(-0.68615174) q[3];
sx q[3];
rz(-0.66303623) q[3];
sx q[3];
rz(1.9239976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9753863) q[0];
sx q[0];
rz(-2.9500742) q[0];
sx q[0];
rz(1.6084877) q[0];
rz(-2.8168822) q[1];
sx q[1];
rz(-1.8638116) q[1];
sx q[1];
rz(2.8285573) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1192577) q[0];
sx q[0];
rz(-2.5637163) q[0];
sx q[0];
rz(1.7176601) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3508505) q[2];
sx q[2];
rz(-1.6306021) q[2];
sx q[2];
rz(-2.4129197) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.43949142) q[1];
sx q[1];
rz(-0.92787039) q[1];
sx q[1];
rz(2.4376252) q[1];
rz(-2.5171017) q[3];
sx q[3];
rz(-1.9510799) q[3];
sx q[3];
rz(2.6759669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2031871) q[2];
sx q[2];
rz(-1.9878191) q[2];
sx q[2];
rz(1.8082834) q[2];
rz(-0.61776727) q[3];
sx q[3];
rz(-2.1963162) q[3];
sx q[3];
rz(-2.0060523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4278605) q[0];
sx q[0];
rz(-1.328631) q[0];
sx q[0];
rz(-2.791639) q[0];
rz(-2.2433544) q[1];
sx q[1];
rz(-2.0571092) q[1];
sx q[1];
rz(2.7877997) q[1];
rz(-0.94793574) q[2];
sx q[2];
rz(-2.2123631) q[2];
sx q[2];
rz(2.3555059) q[2];
rz(-2.4118123) q[3];
sx q[3];
rz(-1.5885708) q[3];
sx q[3];
rz(1.9717168) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
