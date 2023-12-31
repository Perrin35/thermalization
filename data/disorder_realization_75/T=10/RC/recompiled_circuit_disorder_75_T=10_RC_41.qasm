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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0050632) q[0];
sx q[0];
rz(-2.8423474) q[0];
sx q[0];
rz(2.5314999) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3051891) q[2];
sx q[2];
rz(-2.0678389) q[2];
sx q[2];
rz(3.1211987) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4048684) q[1];
sx q[1];
rz(-2.1202592) q[1];
sx q[1];
rz(3.0800372) q[1];
rz(-pi) q[2];
rz(-1.2526413) q[3];
sx q[3];
rz(-2.3080024) q[3];
sx q[3];
rz(2.6300501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3966763) q[2];
sx q[2];
rz(-1.6884721) q[2];
sx q[2];
rz(0.35787004) q[2];
rz(0.19168028) q[3];
sx q[3];
rz(-0.43281698) q[3];
sx q[3];
rz(-2.5884957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1382004) q[0];
sx q[0];
rz(-1.9906185) q[0];
sx q[0];
rz(3.0733118) q[0];
rz(2.1349019) q[1];
sx q[1];
rz(-3.0196562) q[1];
sx q[1];
rz(-0.65111792) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3119445) q[0];
sx q[0];
rz(-0.94046578) q[0];
sx q[0];
rz(-1.3208566) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6158193) q[2];
sx q[2];
rz(-1.5500881) q[2];
sx q[2];
rz(-2.6244342) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5755641) q[1];
sx q[1];
rz(-2.063942) q[1];
sx q[1];
rz(-0.50599392) q[1];
rz(-1.7884368) q[3];
sx q[3];
rz(-1.0895035) q[3];
sx q[3];
rz(-2.7841115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6887001) q[2];
sx q[2];
rz(-2.9850027) q[2];
sx q[2];
rz(-0.90332705) q[2];
rz(-2.3870758) q[3];
sx q[3];
rz(-1.4947596) q[3];
sx q[3];
rz(0.27404684) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2925401) q[0];
sx q[0];
rz(-0.64340574) q[0];
sx q[0];
rz(-2.6480411) q[0];
rz(-1.6429398) q[1];
sx q[1];
rz(-0.40619266) q[1];
sx q[1];
rz(-1.0292056) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.670186) q[0];
sx q[0];
rz(-1.1364511) q[0];
sx q[0];
rz(-1.9451408) q[0];
x q[1];
rz(0.87666338) q[2];
sx q[2];
rz(-0.52774094) q[2];
sx q[2];
rz(-1.1610247) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5268847) q[1];
sx q[1];
rz(-1.7803368) q[1];
sx q[1];
rz(-0.31063147) q[1];
rz(3.0164099) q[3];
sx q[3];
rz(-1.0297965) q[3];
sx q[3];
rz(-2.560905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2091973) q[2];
sx q[2];
rz(-2.5019427) q[2];
sx q[2];
rz(-1.8445245) q[2];
rz(2.5629937) q[3];
sx q[3];
rz(-1.2207299) q[3];
sx q[3];
rz(-2.4659618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9349174) q[0];
sx q[0];
rz(-2.4949555) q[0];
sx q[0];
rz(0.40859616) q[0];
rz(1.4273377) q[1];
sx q[1];
rz(-1.4988377) q[1];
sx q[1];
rz(-0.88159195) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13215412) q[0];
sx q[0];
rz(-2.3883925) q[0];
sx q[0];
rz(-0.25235812) q[0];
x q[1];
rz(0.84882952) q[2];
sx q[2];
rz(-2.1049044) q[2];
sx q[2];
rz(1.7497077) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3410014) q[1];
sx q[1];
rz(-0.87716555) q[1];
sx q[1];
rz(2.4800406) q[1];
x q[2];
rz(-2.9922585) q[3];
sx q[3];
rz(-0.42156005) q[3];
sx q[3];
rz(-1.4967407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0435698) q[2];
sx q[2];
rz(-1.4901525) q[2];
sx q[2];
rz(1.1126474) q[2];
rz(2.5860795) q[3];
sx q[3];
rz(-1.3534618) q[3];
sx q[3];
rz(1.6120733) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1495789) q[0];
sx q[0];
rz(-1.1554138) q[0];
sx q[0];
rz(1.7768815) q[0];
rz(-2.8240906) q[1];
sx q[1];
rz(-2.1798539) q[1];
sx q[1];
rz(-0.11725765) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2909572) q[0];
sx q[0];
rz(-0.98941776) q[0];
sx q[0];
rz(-1.7270389) q[0];
rz(-pi) q[1];
rz(-0.30585441) q[2];
sx q[2];
rz(-0.62539414) q[2];
sx q[2];
rz(1.2139699) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.458141) q[1];
sx q[1];
rz(-1.7675753) q[1];
sx q[1];
rz(1.4139921) q[1];
rz(0.98362555) q[3];
sx q[3];
rz(-1.5067325) q[3];
sx q[3];
rz(-0.6928882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7964898) q[2];
sx q[2];
rz(-1.779665) q[2];
sx q[2];
rz(1.3797181) q[2];
rz(-1.1896677) q[3];
sx q[3];
rz(-2.9822571) q[3];
sx q[3];
rz(0.074507944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4554493) q[0];
sx q[0];
rz(-2.3685072) q[0];
sx q[0];
rz(-2.7129052) q[0];
rz(-1.5594257) q[2];
sx q[2];
rz(-1.5307384) q[2];
sx q[2];
rz(0.96206059) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1888215) q[1];
sx q[1];
rz(-1.8535887) q[1];
sx q[1];
rz(-0.69510831) q[1];
rz(-pi) q[2];
rz(2.2599254) q[3];
sx q[3];
rz(-1.1881184) q[3];
sx q[3];
rz(-1.6430287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.32249054) q[2];
sx q[2];
rz(-0.97444797) q[2];
sx q[2];
rz(2.0325913) q[2];
rz(-1.2935982) q[3];
sx q[3];
rz(-1.785991) q[3];
sx q[3];
rz(-3.0343645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7727707) q[0];
sx q[0];
rz(-1.6596376) q[0];
sx q[0];
rz(-0.014904508) q[0];
rz(-0.42124721) q[1];
sx q[1];
rz(-2.0887471) q[1];
sx q[1];
rz(-0.79963911) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0945471) q[0];
sx q[0];
rz(-0.73091113) q[0];
sx q[0];
rz(0.38696179) q[0];
rz(-pi) q[1];
rz(1.2866576) q[2];
sx q[2];
rz(-1.26789) q[2];
sx q[2];
rz(0.45380935) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7176197) q[1];
sx q[1];
rz(-1.4255925) q[1];
sx q[1];
rz(1.6270301) q[1];
rz(-pi) q[2];
rz(-3.0040967) q[3];
sx q[3];
rz(-2.1318448) q[3];
sx q[3];
rz(1.4631127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5853167) q[2];
sx q[2];
rz(-0.63627807) q[2];
sx q[2];
rz(0.9220534) q[2];
rz(-1.3098035) q[3];
sx q[3];
rz(-1.2049048) q[3];
sx q[3];
rz(2.183765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9668982) q[0];
sx q[0];
rz(-2.0860724) q[0];
sx q[0];
rz(2.877537) q[0];
rz(1.4008201) q[1];
sx q[1];
rz(-1.4512647) q[1];
sx q[1];
rz(2.82428) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2469651) q[0];
sx q[0];
rz(-1.4980226) q[0];
sx q[0];
rz(-1.1829681) q[0];
rz(1.1515456) q[2];
sx q[2];
rz(-0.37852415) q[2];
sx q[2];
rz(0.35767698) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.322284) q[1];
sx q[1];
rz(-1.6849736) q[1];
sx q[1];
rz(1.4788005) q[1];
rz(-pi) q[2];
rz(-0.64538892) q[3];
sx q[3];
rz(-1.4374975) q[3];
sx q[3];
rz(2.3440839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7060966) q[2];
sx q[2];
rz(-2.2075682) q[2];
sx q[2];
rz(-1.502011) q[2];
rz(2.8912985) q[3];
sx q[3];
rz(-1.4151662) q[3];
sx q[3];
rz(1.9992453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
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
rz(-2.6622488) q[1];
sx q[1];
rz(-1.665325) q[1];
sx q[1];
rz(3.0260578) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.696366) q[0];
sx q[0];
rz(-2.1292902) q[0];
sx q[0];
rz(0.88748705) q[0];
rz(-0.71474448) q[2];
sx q[2];
rz(-1.7822767) q[2];
sx q[2];
rz(-1.5874869) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.83976165) q[1];
sx q[1];
rz(-1.2823732) q[1];
sx q[1];
rz(-1.1327869) q[1];
rz(-1.450591) q[3];
sx q[3];
rz(-2.2634441) q[3];
sx q[3];
rz(-1.8409178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.66403786) q[2];
sx q[2];
rz(-2.0319735) q[2];
sx q[2];
rz(-1.4257365) q[2];
rz(1.7410295) q[3];
sx q[3];
rz(-2.0934584) q[3];
sx q[3];
rz(-0.28276309) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7857159) q[0];
sx q[0];
rz(-0.73369217) q[0];
sx q[0];
rz(2.8724331) q[0];
rz(-1.0956988) q[1];
sx q[1];
rz(-2.2311189) q[1];
sx q[1];
rz(-1.3795308) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3740765) q[0];
sx q[0];
rz(-2.9449468) q[0];
sx q[0];
rz(2.1392512) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3921521) q[2];
sx q[2];
rz(-2.7182455) q[2];
sx q[2];
rz(1.4942126) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5386242) q[1];
sx q[1];
rz(-1.5403962) q[1];
sx q[1];
rz(1.3974742) q[1];
rz(-1.2606603) q[3];
sx q[3];
rz(-1.2536067) q[3];
sx q[3];
rz(-2.1807501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0925838) q[2];
sx q[2];
rz(-2.9276586) q[2];
sx q[2];
rz(-1.6258378) q[2];
rz(-1.9745291) q[3];
sx q[3];
rz(-1.556282) q[3];
sx q[3];
rz(-1.0313755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2172858) q[0];
sx q[0];
rz(-1.3047682) q[0];
sx q[0];
rz(-2.7346942) q[0];
rz(2.6869607) q[1];
sx q[1];
rz(-1.1063207) q[1];
sx q[1];
rz(2.8938821) q[1];
rz(2.968593) q[2];
sx q[2];
rz(-2.0377918) q[2];
sx q[2];
rz(1.3755058) q[2];
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
