OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5972714) q[0];
sx q[0];
rz(-0.40364021) q[0];
sx q[0];
rz(-0.37024745) q[0];
rz(0.20180841) q[1];
sx q[1];
rz(-1.2528074) q[1];
sx q[1];
rz(1.7226146) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1184074) q[0];
sx q[0];
rz(-1.4010795) q[0];
sx q[0];
rz(2.8939308) q[0];
x q[1];
rz(-2.6909157) q[2];
sx q[2];
rz(-0.55826954) q[2];
sx q[2];
rz(0.4980586) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.19810361) q[1];
sx q[1];
rz(-1.6232821) q[1];
sx q[1];
rz(-1.0204888) q[1];
x q[2];
rz(-2.3787093) q[3];
sx q[3];
rz(-1.3370822) q[3];
sx q[3];
rz(1.2771311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7449164) q[2];
sx q[2];
rz(-1.6884721) q[2];
sx q[2];
rz(2.7837226) q[2];
rz(-2.9499124) q[3];
sx q[3];
rz(-2.7087757) q[3];
sx q[3];
rz(-0.55309692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0033922694) q[0];
sx q[0];
rz(-1.1509742) q[0];
sx q[0];
rz(0.068280846) q[0];
rz(-2.1349019) q[1];
sx q[1];
rz(-0.12193646) q[1];
sx q[1];
rz(-0.65111792) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7333974) q[0];
sx q[0];
rz(-1.7719643) q[0];
sx q[0];
rz(2.4961619) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5468555) q[2];
sx q[2];
rz(-2.0964453) q[2];
sx q[2];
rz(1.0416232) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.56602851) q[1];
sx q[1];
rz(-2.063942) q[1];
sx q[1];
rz(-2.6355987) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7495456) q[3];
sx q[3];
rz(-2.6169176) q[3];
sx q[3];
rz(2.3384561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.45289257) q[2];
sx q[2];
rz(-0.15658997) q[2];
sx q[2];
rz(2.2382656) q[2];
rz(-0.7545169) q[3];
sx q[3];
rz(-1.4947596) q[3];
sx q[3];
rz(-0.27404684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84905255) q[0];
sx q[0];
rz(-2.4981869) q[0];
sx q[0];
rz(-2.6480411) q[0];
rz(1.4986528) q[1];
sx q[1];
rz(-2.7354) q[1];
sx q[1];
rz(1.0292056) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8619336) q[0];
sx q[0];
rz(-2.5760981) q[0];
sx q[0];
rz(-2.4740567) q[0];
rz(-pi) q[1];
rz(-1.9919954) q[2];
sx q[2];
rz(-1.2427949) q[2];
sx q[2];
rz(1.0331819) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.610565) q[1];
sx q[1];
rz(-2.768801) q[1];
sx q[1];
rz(-0.60786604) q[1];
rz(-1.7756996) q[3];
sx q[3];
rz(-0.5538867) q[3];
sx q[3];
rz(-0.82034558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2091973) q[2];
sx q[2];
rz(-0.63964996) q[2];
sx q[2];
rz(1.2970682) q[2];
rz(0.57859892) q[3];
sx q[3];
rz(-1.9208627) q[3];
sx q[3];
rz(0.67563081) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9349174) q[0];
sx q[0];
rz(-0.6466372) q[0];
sx q[0];
rz(0.40859616) q[0];
rz(-1.4273377) q[1];
sx q[1];
rz(-1.6427549) q[1];
sx q[1];
rz(-0.88159195) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13215412) q[0];
sx q[0];
rz(-0.75320019) q[0];
sx q[0];
rz(-2.8892345) q[0];
x q[1];
rz(-2.300823) q[2];
sx q[2];
rz(-0.86849125) q[2];
sx q[2];
rz(2.4384987) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.2320764) q[1];
sx q[1];
rz(-1.0788003) q[1];
sx q[1];
rz(-2.3823649) q[1];
x q[2];
rz(-2.7241957) q[3];
sx q[3];
rz(-1.6317123) q[3];
sx q[3];
rz(0.062373769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0435698) q[2];
sx q[2];
rz(-1.4901525) q[2];
sx q[2];
rz(1.1126474) q[2];
rz(-2.5860795) q[3];
sx q[3];
rz(-1.3534618) q[3];
sx q[3];
rz(1.5295193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1495789) q[0];
sx q[0];
rz(-1.1554138) q[0];
sx q[0];
rz(1.7768815) q[0];
rz(-0.31750202) q[1];
sx q[1];
rz(-0.96173871) q[1];
sx q[1];
rz(-0.11725765) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3661341) q[0];
sx q[0];
rz(-1.4403847) q[0];
sx q[0];
rz(-2.5545757) q[0];
rz(-pi) q[1];
x q[1];
rz(0.60301493) q[2];
sx q[2];
rz(-1.3935967) q[2];
sx q[2];
rz(0.10620968) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.77814233) q[1];
sx q[1];
rz(-2.8906129) q[1];
sx q[1];
rz(-2.4771033) q[1];
rz(1.6860784) q[3];
sx q[3];
rz(-2.5513463) q[3];
sx q[3];
rz(0.78199996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7964898) q[2];
sx q[2];
rz(-1.779665) q[2];
sx q[2];
rz(1.7618746) q[2];
rz(-1.1896677) q[3];
sx q[3];
rz(-0.15933557) q[3];
sx q[3];
rz(-0.074507944) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18774408) q[0];
sx q[0];
rz(-1.501361) q[0];
sx q[0];
rz(0.70621079) q[0];
rz(-1.1114936) q[1];
sx q[1];
rz(-0.62875426) q[1];
sx q[1];
rz(0.10791735) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88718016) q[0];
sx q[0];
rz(-2.2590056) q[0];
sx q[0];
rz(1.1855017) q[0];
rz(1.5594257) q[2];
sx q[2];
rz(-1.6108542) q[2];
sx q[2];
rz(0.96206059) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.38938746) q[1];
sx q[1];
rz(-0.90837332) q[1];
sx q[1];
rz(-1.9325158) q[1];
rz(-pi) q[2];
rz(2.1351486) q[3];
sx q[3];
rz(-2.3688149) q[3];
sx q[3];
rz(-0.3533065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.32249054) q[2];
sx q[2];
rz(-2.1671447) q[2];
sx q[2];
rz(-1.1090013) q[2];
rz(1.2935982) q[3];
sx q[3];
rz(-1.3556017) q[3];
sx q[3];
rz(-3.0343645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.368822) q[0];
sx q[0];
rz(-1.6596376) q[0];
sx q[0];
rz(-0.014904508) q[0];
rz(0.42124721) q[1];
sx q[1];
rz(-2.0887471) q[1];
sx q[1];
rz(-2.3419535) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22916238) q[0];
sx q[0];
rz(-1.8254571) q[0];
sx q[0];
rz(-2.4486662) q[0];
rz(-pi) q[1];
x q[1];
rz(1.854935) q[2];
sx q[2];
rz(-1.8737027) q[2];
sx q[2];
rz(0.45380935) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0886503) q[1];
sx q[1];
rz(-2.9859516) q[1];
sx q[1];
rz(0.36693962) q[1];
rz(0.13749595) q[3];
sx q[3];
rz(-1.0097479) q[3];
sx q[3];
rz(-1.4631127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.556276) q[2];
sx q[2];
rz(-0.63627807) q[2];
sx q[2];
rz(-0.9220534) q[2];
rz(1.3098035) q[3];
sx q[3];
rz(-1.9366879) q[3];
sx q[3];
rz(2.183765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17469445) q[0];
sx q[0];
rz(-1.0555203) q[0];
sx q[0];
rz(-2.877537) q[0];
rz(-1.7407725) q[1];
sx q[1];
rz(-1.6903279) q[1];
sx q[1];
rz(-2.82428) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6415629) q[0];
sx q[0];
rz(-2.7473358) q[0];
sx q[0];
rz(1.3803598) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9900471) q[2];
sx q[2];
rz(-2.7630685) q[2];
sx q[2];
rz(0.35767698) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6416157) q[1];
sx q[1];
rz(-2.9950905) q[1];
sx q[1];
rz(-2.4661857) q[1];
x q[2];
rz(0.21934261) q[3];
sx q[3];
rz(-2.4845124) q[3];
sx q[3];
rz(0.59857644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.43549609) q[2];
sx q[2];
rz(-0.93402445) q[2];
sx q[2];
rz(-1.502011) q[2];
rz(2.8912985) q[3];
sx q[3];
rz(-1.4151662) q[3];
sx q[3];
rz(-1.1423473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2450927) q[0];
sx q[0];
rz(-2.8332906) q[0];
sx q[0];
rz(-1.7171575) q[0];
rz(-0.4793438) q[1];
sx q[1];
rz(-1.4762676) q[1];
sx q[1];
rz(-0.11553484) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53287017) q[0];
sx q[0];
rz(-2.1358129) q[0];
sx q[0];
rz(-2.4633521) q[0];
x q[1];
rz(1.2938415) q[2];
sx q[2];
rz(-0.87522725) q[2];
sx q[2];
rz(2.9447174) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2781196) q[1];
sx q[1];
rz(-1.9895456) q[1];
sx q[1];
rz(2.8249885) q[1];
x q[2];
rz(-0.1435189) q[3];
sx q[3];
rz(-0.70129881) q[3];
sx q[3];
rz(-1.4876175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.66403786) q[2];
sx q[2];
rz(-2.0319735) q[2];
sx q[2];
rz(1.7158562) q[2];
rz(1.4005631) q[3];
sx q[3];
rz(-2.0934584) q[3];
sx q[3];
rz(-2.8588296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35587674) q[0];
sx q[0];
rz(-0.73369217) q[0];
sx q[0];
rz(-2.8724331) q[0];
rz(2.0458938) q[1];
sx q[1];
rz(-2.2311189) q[1];
sx q[1];
rz(-1.3795308) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7675161) q[0];
sx q[0];
rz(-0.19664581) q[0];
sx q[0];
rz(-1.0023414) q[0];
rz(-pi) q[1];
x q[1];
rz(0.079897957) q[2];
sx q[2];
rz(-1.1546087) q[2];
sx q[2];
rz(-1.4518567) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3456612) q[1];
sx q[1];
rz(-0.17594166) q[1];
sx q[1];
rz(-1.7453341) q[1];
rz(-pi) q[2];
x q[2];
rz(0.74897154) q[3];
sx q[3];
rz(-0.43991551) q[3];
sx q[3];
rz(-0.16187748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0925838) q[2];
sx q[2];
rz(-0.2139341) q[2];
sx q[2];
rz(-1.5157549) q[2];
rz(1.9745291) q[3];
sx q[3];
rz(-1.5853106) q[3];
sx q[3];
rz(2.1102171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2172858) q[0];
sx q[0];
rz(-1.3047682) q[0];
sx q[0];
rz(-2.7346942) q[0];
rz(-0.45463195) q[1];
sx q[1];
rz(-1.1063207) q[1];
sx q[1];
rz(2.8938821) q[1];
rz(-2.968593) q[2];
sx q[2];
rz(-1.1038008) q[2];
sx q[2];
rz(-1.7660869) q[2];
rz(1.8299673) q[3];
sx q[3];
rz(-1.6298686) q[3];
sx q[3];
rz(2.009404) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
