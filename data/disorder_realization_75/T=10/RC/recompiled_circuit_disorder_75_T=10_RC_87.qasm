OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5443213) q[0];
sx q[0];
rz(3.5452329) q[0];
sx q[0];
rz(9.7950254) q[0];
rz(0.20180841) q[1];
sx q[1];
rz(1.8887853) q[1];
sx q[1];
rz(10.843756) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0050632) q[0];
sx q[0];
rz(-2.8423474) q[0];
sx q[0];
rz(-2.5314999) q[0];
x q[1];
rz(2.6909157) q[2];
sx q[2];
rz(-0.55826954) q[2];
sx q[2];
rz(2.6435341) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7367243) q[1];
sx q[1];
rz(-1.0213335) q[1];
sx q[1];
rz(0.0615555) q[1];
x q[2];
rz(1.8889514) q[3];
sx q[3];
rz(-2.3080024) q[3];
sx q[3];
rz(2.6300501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7449164) q[2];
sx q[2];
rz(-1.4531206) q[2];
sx q[2];
rz(2.7837226) q[2];
rz(-2.9499124) q[3];
sx q[3];
rz(-2.7087757) q[3];
sx q[3];
rz(2.5884957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0033922694) q[0];
sx q[0];
rz(-1.9906185) q[0];
sx q[0];
rz(3.0733118) q[0];
rz(2.1349019) q[1];
sx q[1];
rz(-0.12193646) q[1];
sx q[1];
rz(0.65111792) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90322813) q[0];
sx q[0];
rz(-0.67175409) q[0];
sx q[0];
rz(0.3268468) q[0];
x q[1];
rz(0.041243677) q[2];
sx q[2];
rz(-2.6154499) q[2];
sx q[2];
rz(2.0522842) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.880289) q[1];
sx q[1];
rz(-2.0118879) q[1];
sx q[1];
rz(-2.1217568) q[1];
x q[2];
rz(2.7495456) q[3];
sx q[3];
rz(-0.52467504) q[3];
sx q[3];
rz(2.3384561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6887001) q[2];
sx q[2];
rz(-2.9850027) q[2];
sx q[2];
rz(2.2382656) q[2];
rz(2.3870758) q[3];
sx q[3];
rz(-1.6468331) q[3];
sx q[3];
rz(-2.8675458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84905255) q[0];
sx q[0];
rz(-0.64340574) q[0];
sx q[0];
rz(2.6480411) q[0];
rz(1.4986528) q[1];
sx q[1];
rz(-2.7354) q[1];
sx q[1];
rz(-2.1123871) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.670186) q[0];
sx q[0];
rz(-2.0051415) q[0];
sx q[0];
rz(1.9451408) q[0];
rz(-pi) q[1];
rz(1.1495972) q[2];
sx q[2];
rz(-1.2427949) q[2];
sx q[2];
rz(1.0331819) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5268847) q[1];
sx q[1];
rz(-1.3612559) q[1];
sx q[1];
rz(-2.8309612) q[1];
rz(-pi) q[2];
rz(-3.0164099) q[3];
sx q[3];
rz(-1.0297965) q[3];
sx q[3];
rz(-0.5806877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9323953) q[2];
sx q[2];
rz(-2.5019427) q[2];
sx q[2];
rz(-1.2970682) q[2];
rz(0.57859892) q[3];
sx q[3];
rz(-1.2207299) q[3];
sx q[3];
rz(-0.67563081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2066752) q[0];
sx q[0];
rz(-2.4949555) q[0];
sx q[0];
rz(-0.40859616) q[0];
rz(-1.4273377) q[1];
sx q[1];
rz(-1.6427549) q[1];
sx q[1];
rz(-0.88159195) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0094385) q[0];
sx q[0];
rz(-0.75320019) q[0];
sx q[0];
rz(2.8892345) q[0];
x q[1];
rz(-2.4741715) q[2];
sx q[2];
rz(-2.1758658) q[2];
sx q[2];
rz(2.8990926) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9095163) q[1];
sx q[1];
rz(-1.0788003) q[1];
sx q[1];
rz(2.3823649) q[1];
rz(2.7241957) q[3];
sx q[3];
rz(-1.5098803) q[3];
sx q[3];
rz(0.062373769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0980229) q[2];
sx q[2];
rz(-1.4901525) q[2];
sx q[2];
rz(-2.0289452) q[2];
rz(0.55551314) q[3];
sx q[3];
rz(-1.7881309) q[3];
sx q[3];
rz(1.6120733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(2.9920138) q[0];
sx q[0];
rz(-1.1554138) q[0];
sx q[0];
rz(-1.7768815) q[0];
rz(-2.8240906) q[1];
sx q[1];
rz(-2.1798539) q[1];
sx q[1];
rz(3.024335) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.011615959) q[0];
sx q[0];
rz(-2.5419309) q[0];
sx q[0];
rz(0.2325124) q[0];
rz(1.3567032) q[2];
sx q[2];
rz(-0.97852409) q[2];
sx q[2];
rz(1.5562197) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3634503) q[1];
sx q[1];
rz(-2.8906129) q[1];
sx q[1];
rz(-0.66448934) q[1];
x q[2];
rz(1.6860784) q[3];
sx q[3];
rz(-0.59024631) q[3];
sx q[3];
rz(-0.78199996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.34510288) q[2];
sx q[2];
rz(-1.779665) q[2];
sx q[2];
rz(1.7618746) q[2];
rz(-1.951925) q[3];
sx q[3];
rz(-2.9822571) q[3];
sx q[3];
rz(3.0670847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18774408) q[0];
sx q[0];
rz(-1.501361) q[0];
sx q[0];
rz(-0.70621079) q[0];
rz(2.0300991) q[1];
sx q[1];
rz(-2.5128384) q[1];
sx q[1];
rz(-0.10791735) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4554493) q[0];
sx q[0];
rz(-2.3685072) q[0];
sx q[0];
rz(-2.7129052) q[0];
rz(-pi) q[1];
rz(-1.5594257) q[2];
sx q[2];
rz(-1.6108542) q[2];
sx q[2];
rz(2.1795321) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.94090998) q[1];
sx q[1];
rz(-2.4001277) q[1];
sx q[1];
rz(-0.42592589) q[1];
rz(-0.48072731) q[3];
sx q[3];
rz(-0.93989621) q[3];
sx q[3];
rz(2.7709099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8191021) q[2];
sx q[2];
rz(-2.1671447) q[2];
sx q[2];
rz(2.0325913) q[2];
rz(-1.2935982) q[3];
sx q[3];
rz(-1.3556017) q[3];
sx q[3];
rz(3.0343645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7727707) q[0];
sx q[0];
rz(-1.6596376) q[0];
sx q[0];
rz(0.014904508) q[0];
rz(0.42124721) q[1];
sx q[1];
rz(-2.0887471) q[1];
sx q[1];
rz(0.79963911) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22916238) q[0];
sx q[0];
rz(-1.8254571) q[0];
sx q[0];
rz(2.4486662) q[0];
x q[1];
rz(-2.4104426) q[2];
sx q[2];
rz(-0.41229782) q[2];
sx q[2];
rz(-1.9129802) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0886503) q[1];
sx q[1];
rz(-0.155641) q[1];
sx q[1];
rz(2.774653) q[1];
rz(1.0054672) q[3];
sx q[3];
rz(-1.6871095) q[3];
sx q[3];
rz(3.1073991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5853167) q[2];
sx q[2];
rz(-0.63627807) q[2];
sx q[2];
rz(0.9220534) q[2];
rz(1.8317892) q[3];
sx q[3];
rz(-1.9366879) q[3];
sx q[3];
rz(-2.183765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9668982) q[0];
sx q[0];
rz(-2.0860724) q[0];
sx q[0];
rz(-0.2640557) q[0];
rz(-1.4008201) q[1];
sx q[1];
rz(-1.6903279) q[1];
sx q[1];
rz(-0.31731269) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89462751) q[0];
sx q[0];
rz(-1.4980226) q[0];
sx q[0];
rz(1.9586246) q[0];
x q[1];
rz(1.9900471) q[2];
sx q[2];
rz(-2.7630685) q[2];
sx q[2];
rz(-2.7839157) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6416157) q[1];
sx q[1];
rz(-0.14650211) q[1];
sx q[1];
rz(-2.4661857) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4044912) q[3];
sx q[3];
rz(-0.93207031) q[3];
sx q[3];
rz(-2.2685662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7060966) q[2];
sx q[2];
rz(-2.2075682) q[2];
sx q[2];
rz(1.502011) q[2];
rz(2.8912985) q[3];
sx q[3];
rz(-1.7264265) q[3];
sx q[3];
rz(1.1423473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2450927) q[0];
sx q[0];
rz(-2.8332906) q[0];
sx q[0];
rz(1.4244351) q[0];
rz(-0.4793438) q[1];
sx q[1];
rz(-1.4762676) q[1];
sx q[1];
rz(3.0260578) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45142052) q[0];
sx q[0];
rz(-2.2884986) q[0];
sx q[0];
rz(0.79057981) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.71474448) q[2];
sx q[2];
rz(-1.359316) q[2];
sx q[2];
rz(-1.5541058) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2781196) q[1];
sx q[1];
rz(-1.9895456) q[1];
sx q[1];
rz(0.31660415) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9980738) q[3];
sx q[3];
rz(-0.70129881) q[3];
sx q[3];
rz(1.4876175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4775548) q[2];
sx q[2];
rz(-1.1096191) q[2];
sx q[2];
rz(1.4257365) q[2];
rz(1.7410295) q[3];
sx q[3];
rz(-1.0481342) q[3];
sx q[3];
rz(0.28276309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7857159) q[0];
sx q[0];
rz(-2.4079005) q[0];
sx q[0];
rz(-0.26915959) q[0];
rz(-2.0458938) q[1];
sx q[1];
rz(-2.2311189) q[1];
sx q[1];
rz(1.3795308) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1901967) q[0];
sx q[0];
rz(-1.405389) q[0];
sx q[0];
rz(3.0347546) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1534259) q[2];
sx q[2];
rz(-1.4977314) q[2];
sx q[2];
rz(-0.086581143) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3456612) q[1];
sx q[1];
rz(-0.17594166) q[1];
sx q[1];
rz(1.3962586) q[1];
rz(-2.3926211) q[3];
sx q[3];
rz(-2.7016771) q[3];
sx q[3];
rz(-2.9797152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0925838) q[2];
sx q[2];
rz(-0.2139341) q[2];
sx q[2];
rz(-1.6258378) q[2];
rz(1.1670636) q[3];
sx q[3];
rz(-1.5853106) q[3];
sx q[3];
rz(1.0313755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9243069) q[0];
sx q[0];
rz(-1.3047682) q[0];
sx q[0];
rz(-2.7346942) q[0];
rz(2.6869607) q[1];
sx q[1];
rz(-1.1063207) q[1];
sx q[1];
rz(2.8938821) q[1];
rz(-2.968593) q[2];
sx q[2];
rz(-1.1038008) q[2];
sx q[2];
rz(-1.7660869) q[2];
rz(-1.7975939) q[3];
sx q[3];
rz(-0.26567017) q[3];
sx q[3];
rz(-2.9220823) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
