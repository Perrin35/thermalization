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
rz(0.87822479) q[0];
sx q[0];
rz(3.8642519) q[0];
sx q[0];
rz(10.470358) q[0];
rz(0.0027520952) q[1];
sx q[1];
rz(-1.4581008) q[1];
sx q[1];
rz(-0.72088617) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65675019) q[0];
sx q[0];
rz(-0.67615254) q[0];
sx q[0];
rz(-1.7184896) q[0];
x q[1];
rz(-1.7621668) q[2];
sx q[2];
rz(-1.6211133) q[2];
sx q[2];
rz(2.2309395) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.406479) q[1];
sx q[1];
rz(-1.8476474) q[1];
sx q[1];
rz(0.39445725) q[1];
rz(-1.4414532) q[3];
sx q[3];
rz(-2.3096519) q[3];
sx q[3];
rz(1.3373614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2846994) q[2];
sx q[2];
rz(-1.2130986) q[2];
sx q[2];
rz(-2.6846679) q[2];
rz(0.65447718) q[3];
sx q[3];
rz(-2.579253) q[3];
sx q[3];
rz(0.42403179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65861312) q[0];
sx q[0];
rz(-2.787866) q[0];
sx q[0];
rz(0.60802996) q[0];
rz(1.1990625) q[1];
sx q[1];
rz(-2.6741195) q[1];
sx q[1];
rz(2.3466568) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6632961) q[0];
sx q[0];
rz(-1.2353871) q[0];
sx q[0];
rz(-1.8778503) q[0];
x q[1];
rz(3.0367467) q[2];
sx q[2];
rz(-2.8730132) q[2];
sx q[2];
rz(-1.9263445) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7424189) q[1];
sx q[1];
rz(-1.8493506) q[1];
sx q[1];
rz(-2.2225077) q[1];
x q[2];
rz(0.62165798) q[3];
sx q[3];
rz(-1.835653) q[3];
sx q[3];
rz(1.6408629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.46637154) q[2];
sx q[2];
rz(-1.0543062) q[2];
sx q[2];
rz(-1.6611453) q[2];
rz(-2.772707) q[3];
sx q[3];
rz(-1.2494913) q[3];
sx q[3];
rz(-2.4108385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22496255) q[0];
sx q[0];
rz(-1.0039622) q[0];
sx q[0];
rz(-2.4555901) q[0];
rz(-1.2432159) q[1];
sx q[1];
rz(-0.22626433) q[1];
sx q[1];
rz(-2.539198) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7334977) q[0];
sx q[0];
rz(-2.0614834) q[0];
sx q[0];
rz(0.48510929) q[0];
rz(2.2445065) q[2];
sx q[2];
rz(-3.0897475) q[2];
sx q[2];
rz(-2.12814) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2418306) q[1];
sx q[1];
rz(-1.4736403) q[1];
sx q[1];
rz(-1.001557) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1021578) q[3];
sx q[3];
rz(-1.8640362) q[3];
sx q[3];
rz(1.6654429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.38597044) q[2];
sx q[2];
rz(-1.8381511) q[2];
sx q[2];
rz(-1.8899567) q[2];
rz(-2.5899467) q[3];
sx q[3];
rz(-2.9988204) q[3];
sx q[3];
rz(0.84757203) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30152339) q[0];
sx q[0];
rz(-1.2395549) q[0];
sx q[0];
rz(0.10313343) q[0];
rz(0.81941191) q[1];
sx q[1];
rz(-0.90140072) q[1];
sx q[1];
rz(-0.55319667) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90097809) q[0];
sx q[0];
rz(-1.2831956) q[0];
sx q[0];
rz(1.1905627) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0248839) q[2];
sx q[2];
rz(-1.8736412) q[2];
sx q[2];
rz(-0.63797039) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4955289) q[1];
sx q[1];
rz(-0.49139443) q[1];
sx q[1];
rz(-1.8629406) q[1];
rz(-pi) q[2];
rz(-1.1392713) q[3];
sx q[3];
rz(-1.478704) q[3];
sx q[3];
rz(1.0669277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.72054046) q[2];
sx q[2];
rz(-2.4002176) q[2];
sx q[2];
rz(0.53483024) q[2];
rz(-2.363291) q[3];
sx q[3];
rz(-2.4141267) q[3];
sx q[3];
rz(2.9261869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68742037) q[0];
sx q[0];
rz(-1.146831) q[0];
sx q[0];
rz(-2.4046894) q[0];
rz(-3.0829644) q[1];
sx q[1];
rz(-0.77719378) q[1];
sx q[1];
rz(0.58707213) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.042087) q[0];
sx q[0];
rz(-1.5387968) q[0];
sx q[0];
rz(0.93670292) q[0];
rz(1.8002073) q[2];
sx q[2];
rz(-1.5093813) q[2];
sx q[2];
rz(0.0030850911) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.88564155) q[1];
sx q[1];
rz(-1.5382861) q[1];
sx q[1];
rz(-1.2684275) q[1];
rz(1.8892509) q[3];
sx q[3];
rz(-1.8618004) q[3];
sx q[3];
rz(0.97614563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.29207841) q[2];
sx q[2];
rz(-2.2167315) q[2];
sx q[2];
rz(-0.64797956) q[2];
rz(2.5255711) q[3];
sx q[3];
rz(-1.5024622) q[3];
sx q[3];
rz(-0.047671635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7306526) q[0];
sx q[0];
rz(-0.76496449) q[0];
sx q[0];
rz(-2.1867645) q[0];
rz(-3.094589) q[1];
sx q[1];
rz(-2.4689597) q[1];
sx q[1];
rz(2.1254553) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8076626) q[0];
sx q[0];
rz(-0.74650812) q[0];
sx q[0];
rz(-2.4115415) q[0];
rz(2.1121641) q[2];
sx q[2];
rz(-2.2185791) q[2];
sx q[2];
rz(1.364653) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.25364843) q[1];
sx q[1];
rz(-1.5852815) q[1];
sx q[1];
rz(0.50555857) q[1];
rz(-pi) q[2];
rz(-1.5136857) q[3];
sx q[3];
rz(-2.1942003) q[3];
sx q[3];
rz(0.36950612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5434692) q[2];
sx q[2];
rz(-1.2657961) q[2];
sx q[2];
rz(-0.050696105) q[2];
rz(0.051004574) q[3];
sx q[3];
rz(-1.3011322) q[3];
sx q[3];
rz(-0.72371975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8766668) q[0];
sx q[0];
rz(-0.51487881) q[0];
sx q[0];
rz(1.4798973) q[0];
rz(-1.121608) q[1];
sx q[1];
rz(-1.879296) q[1];
sx q[1];
rz(-2.6720537) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2377937) q[0];
sx q[0];
rz(-1.617476) q[0];
sx q[0];
rz(-1.6782128) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9856244) q[2];
sx q[2];
rz(-1.1824338) q[2];
sx q[2];
rz(-2.591382) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.96315876) q[1];
sx q[1];
rz(-1.7443774) q[1];
sx q[1];
rz(1.0976237) q[1];
x q[2];
rz(-0.64823635) q[3];
sx q[3];
rz(-2.2508374) q[3];
sx q[3];
rz(1.4095962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.40662128) q[2];
sx q[2];
rz(-1.5512286) q[2];
sx q[2];
rz(-2.0242019) q[2];
rz(-2.0495074) q[3];
sx q[3];
rz(-1.4039682) q[3];
sx q[3];
rz(-2.7058069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1209054) q[0];
sx q[0];
rz(-1.1130604) q[0];
sx q[0];
rz(-3.0764965) q[0];
rz(2.8096325) q[1];
sx q[1];
rz(-1.8354225) q[1];
sx q[1];
rz(-1.2750221) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1972361) q[0];
sx q[0];
rz(-1.5612537) q[0];
sx q[0];
rz(-1.3413603) q[0];
rz(2.0135856) q[2];
sx q[2];
rz(-1.4472359) q[2];
sx q[2];
rz(2.3905001) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.79679497) q[1];
sx q[1];
rz(-1.5196504) q[1];
sx q[1];
rz(-0.68294345) q[1];
rz(-pi) q[2];
rz(0.80294369) q[3];
sx q[3];
rz(-1.6979473) q[3];
sx q[3];
rz(-1.91636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.28978213) q[2];
sx q[2];
rz(-1.6921356) q[2];
sx q[2];
rz(0.97274485) q[2];
rz(-0.70074493) q[3];
sx q[3];
rz(-0.86207977) q[3];
sx q[3];
rz(-0.8249445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7026611) q[0];
sx q[0];
rz(-0.40533608) q[0];
sx q[0];
rz(0.40913707) q[0];
rz(1.1459972) q[1];
sx q[1];
rz(-1.0640249) q[1];
sx q[1];
rz(1.2629868) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0008903) q[0];
sx q[0];
rz(-1.5626199) q[0];
sx q[0];
rz(0.00012881669) q[0];
rz(-0.97102286) q[2];
sx q[2];
rz(-2.6749938) q[2];
sx q[2];
rz(-2.3393037) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.990155) q[1];
sx q[1];
rz(-1.7977771) q[1];
sx q[1];
rz(-0.87363665) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.13554445) q[3];
sx q[3];
rz(-2.5672847) q[3];
sx q[3];
rz(1.7611461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.24141773) q[2];
sx q[2];
rz(-1.9200385) q[2];
sx q[2];
rz(-2.9202666) q[2];
rz(-2.6545702) q[3];
sx q[3];
rz(-0.86921391) q[3];
sx q[3];
rz(-1.937449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0814334) q[0];
sx q[0];
rz(-2.8964323) q[0];
sx q[0];
rz(-1.2231476) q[0];
rz(-3.0606015) q[1];
sx q[1];
rz(-1.5206189) q[1];
sx q[1];
rz(-1.6193259) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8318849) q[0];
sx q[0];
rz(-1.4060347) q[0];
sx q[0];
rz(2.805183) q[0];
x q[1];
rz(0.72725216) q[2];
sx q[2];
rz(-1.1071812) q[2];
sx q[2];
rz(1.1905244) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4829068) q[1];
sx q[1];
rz(-1.7665837) q[1];
sx q[1];
rz(1.8808603) q[1];
x q[2];
rz(1.7459938) q[3];
sx q[3];
rz(-2.0591551) q[3];
sx q[3];
rz(-1.7217896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.12485518) q[2];
sx q[2];
rz(-1.3913466) q[2];
sx q[2];
rz(2.4511231) q[2];
rz(2.8705583) q[3];
sx q[3];
rz(-0.46509585) q[3];
sx q[3];
rz(-1.5976228) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4616213) q[0];
sx q[0];
rz(-2.2053056) q[0];
sx q[0];
rz(2.9619138) q[0];
rz(-2.8996254) q[1];
sx q[1];
rz(-2.2503743) q[1];
sx q[1];
rz(1.888884) q[1];
rz(1.3623357) q[2];
sx q[2];
rz(-1.9849384) q[2];
sx q[2];
rz(-2.7445856) q[2];
rz(-2.252752) q[3];
sx q[3];
rz(-1.8163637) q[3];
sx q[3];
rz(1.502542) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
