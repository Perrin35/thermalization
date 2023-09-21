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
rz(-2.7379524) q[0];
sx q[0];
rz(0.37024745) q[0];
rz(0.20180841) q[1];
sx q[1];
rz(-1.2528074) q[1];
sx q[1];
rz(-1.4189781) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0050632) q[0];
sx q[0];
rz(-2.8423474) q[0];
sx q[0];
rz(-0.61009272) q[0];
rz(-pi) q[1];
rz(0.51214829) q[2];
sx q[2];
rz(-1.3379659) q[2];
sx q[2];
rz(1.6793959) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4048684) q[1];
sx q[1];
rz(-1.0213335) q[1];
sx q[1];
rz(3.0800372) q[1];
rz(-2.3787093) q[3];
sx q[3];
rz(-1.8045104) q[3];
sx q[3];
rz(1.8644615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3966763) q[2];
sx q[2];
rz(-1.6884721) q[2];
sx q[2];
rz(2.7837226) q[2];
rz(-0.19168028) q[3];
sx q[3];
rz(-2.7087757) q[3];
sx q[3];
rz(-2.5884957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40819528) q[0];
sx q[0];
rz(-1.3696284) q[0];
sx q[0];
rz(-2.4961619) q[0];
rz(-pi) q[1];
rz(1.5468555) q[2];
sx q[2];
rz(-2.0964453) q[2];
sx q[2];
rz(-1.0416232) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.56602851) q[1];
sx q[1];
rz(-2.063942) q[1];
sx q[1];
rz(-2.6355987) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.39204709) q[3];
sx q[3];
rz(-2.6169176) q[3];
sx q[3];
rz(0.80313659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.45289257) q[2];
sx q[2];
rz(-2.9850027) q[2];
sx q[2];
rz(-2.2382656) q[2];
rz(-2.3870758) q[3];
sx q[3];
rz(-1.6468331) q[3];
sx q[3];
rz(2.8675458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84905255) q[0];
sx q[0];
rz(-2.4981869) q[0];
sx q[0];
rz(2.6480411) q[0];
rz(-1.6429398) q[1];
sx q[1];
rz(-2.7354) q[1];
sx q[1];
rz(1.0292056) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27965901) q[0];
sx q[0];
rz(-0.56549457) q[0];
sx q[0];
rz(-2.4740567) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2649293) q[2];
sx q[2];
rz(-0.52774094) q[2];
sx q[2];
rz(-1.9805679) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2521816) q[1];
sx q[1];
rz(-1.8744138) q[1];
sx q[1];
rz(1.7905411) q[1];
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
sx q[1];
rz(pi/2) q[1];
rz(1.9323953) q[2];
sx q[2];
rz(-0.63964996) q[2];
sx q[2];
rz(-1.2970682) q[2];
rz(-2.5629937) q[3];
sx q[3];
rz(-1.9208627) q[3];
sx q[3];
rz(0.67563081) q[3];
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
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9349174) q[0];
sx q[0];
rz(-2.4949555) q[0];
sx q[0];
rz(2.7329965) q[0];
rz(1.4273377) q[1];
sx q[1];
rz(-1.4988377) q[1];
sx q[1];
rz(2.2600007) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20760575) q[0];
sx q[0];
rz(-2.2946977) q[0];
sx q[0];
rz(1.3408322) q[0];
x q[1];
rz(-0.84882952) q[2];
sx q[2];
rz(-1.0366882) q[2];
sx q[2];
rz(1.7497077) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9095163) q[1];
sx q[1];
rz(-2.0627923) q[1];
sx q[1];
rz(-2.3823649) q[1];
rz(1.5041755) q[3];
sx q[3];
rz(-1.1542218) q[3];
sx q[3];
rz(1.4814324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0435698) q[2];
sx q[2];
rz(-1.6514401) q[2];
sx q[2];
rz(1.1126474) q[2];
rz(-0.55551314) q[3];
sx q[3];
rz(-1.7881309) q[3];
sx q[3];
rz(-1.6120733) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9920138) q[0];
sx q[0];
rz(-1.1554138) q[0];
sx q[0];
rz(1.3647112) q[0];
rz(2.8240906) q[1];
sx q[1];
rz(-2.1798539) q[1];
sx q[1];
rz(-3.024335) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.011615959) q[0];
sx q[0];
rz(-2.5419309) q[0];
sx q[0];
rz(-0.2325124) q[0];
x q[1];
rz(0.30585441) q[2];
sx q[2];
rz(-2.5161985) q[2];
sx q[2];
rz(1.2139699) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.14355625) q[1];
sx q[1];
rz(-1.4170425) q[1];
sx q[1];
rz(2.9424332) q[1];
rz(-pi) q[2];
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
rz(0.34510288) q[2];
sx q[2];
rz(-1.779665) q[2];
sx q[2];
rz(1.3797181) q[2];
rz(-1.1896677) q[3];
sx q[3];
rz(-0.15933557) q[3];
sx q[3];
rz(-0.074507944) q[3];
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
rz(-0.18774408) q[0];
sx q[0];
rz(-1.501361) q[0];
sx q[0];
rz(-2.4353819) q[0];
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
rz(0.43150049) q[0];
sx q[0];
rz(-1.8653231) q[0];
sx q[0];
rz(-2.4157903) q[0];
rz(-pi) q[1];
rz(1.582167) q[2];
sx q[2];
rz(-1.5307384) q[2];
sx q[2];
rz(0.96206059) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1888215) q[1];
sx q[1];
rz(-1.2880039) q[1];
sx q[1];
rz(-2.4464843) q[1];
rz(2.6608653) q[3];
sx q[3];
rz(-2.2016964) q[3];
sx q[3];
rz(0.37068278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8191021) q[2];
sx q[2];
rz(-0.97444797) q[2];
sx q[2];
rz(1.1090013) q[2];
rz(1.2935982) q[3];
sx q[3];
rz(-1.785991) q[3];
sx q[3];
rz(-0.10722815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.368822) q[0];
sx q[0];
rz(-1.4819551) q[0];
sx q[0];
rz(0.014904508) q[0];
rz(2.7203454) q[1];
sx q[1];
rz(-2.0887471) q[1];
sx q[1];
rz(2.3419535) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9124303) q[0];
sx q[0];
rz(-1.3161356) q[0];
sx q[0];
rz(2.4486662) q[0];
x q[1];
rz(0.31475474) q[2];
sx q[2];
rz(-1.2999279) q[2];
sx q[2];
rz(-1.0300919) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0529424) q[1];
sx q[1];
rz(-2.9859516) q[1];
sx q[1];
rz(-0.36693962) q[1];
rz(-pi) q[2];
rz(1.0054672) q[3];
sx q[3];
rz(-1.4544832) q[3];
sx q[3];
rz(0.034193589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.556276) q[2];
sx q[2];
rz(-2.5053146) q[2];
sx q[2];
rz(-2.2195393) q[2];
rz(-1.8317892) q[3];
sx q[3];
rz(-1.2049048) q[3];
sx q[3];
rz(0.95782763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(0.2640557) q[0];
rz(-1.4008201) q[1];
sx q[1];
rz(-1.4512647) q[1];
sx q[1];
rz(-2.82428) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89462751) q[0];
sx q[0];
rz(-1.4980226) q[0];
sx q[0];
rz(1.9586246) q[0];
rz(1.9900471) q[2];
sx q[2];
rz(-0.37852415) q[2];
sx q[2];
rz(-0.35767698) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.81930868) q[1];
sx q[1];
rz(-1.456619) q[1];
sx q[1];
rz(1.6627922) q[1];
rz(-2.4962037) q[3];
sx q[3];
rz(-1.7040952) q[3];
sx q[3];
rz(-0.79750878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7060966) q[2];
sx q[2];
rz(-0.93402445) q[2];
sx q[2];
rz(-1.502011) q[2];
rz(-0.25029415) q[3];
sx q[3];
rz(-1.7264265) q[3];
sx q[3];
rz(1.1423473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2450927) q[0];
sx q[0];
rz(-0.30830202) q[0];
sx q[0];
rz(1.7171575) q[0];
rz(0.4793438) q[1];
sx q[1];
rz(-1.665325) q[1];
sx q[1];
rz(-0.11553484) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53287017) q[0];
sx q[0];
rz(-1.0057797) q[0];
sx q[0];
rz(-2.4633521) q[0];
x q[1];
rz(1.8477511) q[2];
sx q[2];
rz(-2.2663654) q[2];
sx q[2];
rz(2.9447174) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.301831) q[1];
sx q[1];
rz(-1.2823732) q[1];
sx q[1];
rz(-2.0088058) q[1];
rz(-pi) q[2];
rz(0.1435189) q[3];
sx q[3];
rz(-0.70129881) q[3];
sx q[3];
rz(-1.6539751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.66403786) q[2];
sx q[2];
rz(-1.1096191) q[2];
sx q[2];
rz(1.4257365) q[2];
rz(-1.4005631) q[3];
sx q[3];
rz(-2.0934584) q[3];
sx q[3];
rz(2.8588296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35587674) q[0];
sx q[0];
rz(-2.4079005) q[0];
sx q[0];
rz(-2.8724331) q[0];
rz(-2.0458938) q[1];
sx q[1];
rz(-0.91047374) q[1];
sx q[1];
rz(1.7620618) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9513959) q[0];
sx q[0];
rz(-1.7362036) q[0];
sx q[0];
rz(0.10683807) q[0];
rz(1.9881667) q[2];
sx q[2];
rz(-1.4977314) q[2];
sx q[2];
rz(3.0550115) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3456612) q[1];
sx q[1];
rz(-0.17594166) q[1];
sx q[1];
rz(1.3962586) q[1];
rz(-pi) q[2];
rz(2.8096301) q[3];
sx q[3];
rz(-1.276607) q[3];
sx q[3];
rz(-0.709579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.049008869) q[2];
sx q[2];
rz(-2.9276586) q[2];
sx q[2];
rz(-1.5157549) q[2];
rz(-1.9745291) q[3];
sx q[3];
rz(-1.5853106) q[3];
sx q[3];
rz(-2.1102171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
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
rz(-1.8998014) q[2];
sx q[2];
rz(-0.49578373) q[2];
sx q[2];
rz(1.0052581) q[2];
rz(1.7975939) q[3];
sx q[3];
rz(-2.8759225) q[3];
sx q[3];
rz(0.21951036) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
