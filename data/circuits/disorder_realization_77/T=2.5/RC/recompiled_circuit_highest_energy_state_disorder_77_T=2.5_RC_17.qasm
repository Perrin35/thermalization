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
rz(2.6440808) q[0];
sx q[0];
rz(-1.3389791) q[0];
sx q[0];
rz(2.8259377) q[0];
rz(1.1701801) q[1];
sx q[1];
rz(-0.38771954) q[1];
sx q[1];
rz(0.60400909) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3919298) q[0];
sx q[0];
rz(-1.5201585) q[0];
sx q[0];
rz(-1.873436) q[0];
rz(-pi) q[1];
x q[1];
rz(0.9004132) q[2];
sx q[2];
rz(-1.9758103) q[2];
sx q[2];
rz(-0.23460282) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0993201) q[1];
sx q[1];
rz(-1.0612773) q[1];
sx q[1];
rz(-2.7541942) q[1];
rz(0.1084566) q[3];
sx q[3];
rz(-1.2634957) q[3];
sx q[3];
rz(-1.1600329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9903367) q[2];
sx q[2];
rz(-1.9467111) q[2];
sx q[2];
rz(1.7830431) q[2];
rz(1.5938866) q[3];
sx q[3];
rz(-2.2113776) q[3];
sx q[3];
rz(-1.6319298) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5556521) q[0];
sx q[0];
rz(-0.63783115) q[0];
sx q[0];
rz(1.9441388) q[0];
rz(2.6400631) q[1];
sx q[1];
rz(-1.5781559) q[1];
sx q[1];
rz(-1.4362358) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2510808) q[0];
sx q[0];
rz(-0.62411004) q[0];
sx q[0];
rz(1.1680383) q[0];
rz(-0.11709638) q[2];
sx q[2];
rz(-2.4352322) q[2];
sx q[2];
rz(-1.381284) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.20138559) q[1];
sx q[1];
rz(-1.4680982) q[1];
sx q[1];
rz(-0.19725712) q[1];
x q[2];
rz(2.7327939) q[3];
sx q[3];
rz(-0.49064562) q[3];
sx q[3];
rz(-1.4178432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2531835) q[2];
sx q[2];
rz(-1.2359572) q[2];
sx q[2];
rz(-3.1186228) q[2];
rz(-2.2394771) q[3];
sx q[3];
rz(-1.918856) q[3];
sx q[3];
rz(-0.42660108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.4334634) q[0];
sx q[0];
rz(-1.7925649) q[0];
sx q[0];
rz(2.1500812) q[0];
rz(-1.0333215) q[1];
sx q[1];
rz(-0.28305498) q[1];
sx q[1];
rz(1.2352157) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69312364) q[0];
sx q[0];
rz(-0.64769324) q[0];
sx q[0];
rz(2.7231587) q[0];
x q[1];
rz(2.5318145) q[2];
sx q[2];
rz(-1.1705361) q[2];
sx q[2];
rz(-2.5230809) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.72690645) q[1];
sx q[1];
rz(-0.79830805) q[1];
sx q[1];
rz(1.2220135) q[1];
x q[2];
rz(-0.97930538) q[3];
sx q[3];
rz(-2.0404173) q[3];
sx q[3];
rz(-0.33794935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.027355) q[2];
sx q[2];
rz(-2.3894252) q[2];
sx q[2];
rz(-2.1503964) q[2];
rz(1.8441955) q[3];
sx q[3];
rz(-1.8810561) q[3];
sx q[3];
rz(1.4170925) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39795136) q[0];
sx q[0];
rz(-0.93248168) q[0];
sx q[0];
rz(2.3241296) q[0];
rz(-0.3262597) q[1];
sx q[1];
rz(-1.1810818) q[1];
sx q[1];
rz(2.356333) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0164675) q[0];
sx q[0];
rz(-2.7332691) q[0];
sx q[0];
rz(0.84618469) q[0];
rz(-pi) q[1];
rz(-2.4982959) q[2];
sx q[2];
rz(-1.729082) q[2];
sx q[2];
rz(0.23223755) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8977114) q[1];
sx q[1];
rz(-2.5088025) q[1];
sx q[1];
rz(-2.3348319) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6503851) q[3];
sx q[3];
rz(-1.9246939) q[3];
sx q[3];
rz(-2.6864664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0205445) q[2];
sx q[2];
rz(-0.64457568) q[2];
sx q[2];
rz(-2.5784967) q[2];
rz(-0.8935039) q[3];
sx q[3];
rz(-1.1734791) q[3];
sx q[3];
rz(1.2347429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9017482) q[0];
sx q[0];
rz(-2.8185066) q[0];
sx q[0];
rz(-1.5455986) q[0];
rz(-2.1196938) q[1];
sx q[1];
rz(-1.1232168) q[1];
sx q[1];
rz(2.9603069) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55953944) q[0];
sx q[0];
rz(-1.3783558) q[0];
sx q[0];
rz(-1.0715225) q[0];
rz(0.47563817) q[2];
sx q[2];
rz(-1.5283094) q[2];
sx q[2];
rz(-1.3459537) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0376772) q[1];
sx q[1];
rz(-1.5988776) q[1];
sx q[1];
rz(1.3485391) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.44293176) q[3];
sx q[3];
rz(-1.4287523) q[3];
sx q[3];
rz(1.0012817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.01866092) q[2];
sx q[2];
rz(-2.4666726) q[2];
sx q[2];
rz(1.8956511) q[2];
rz(-2.5490226) q[3];
sx q[3];
rz(-1.6962681) q[3];
sx q[3];
rz(2.3004801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5331921) q[0];
sx q[0];
rz(-2.1493981) q[0];
sx q[0];
rz(2.8327508) q[0];
rz(-2.8187075) q[1];
sx q[1];
rz(-0.37728089) q[1];
sx q[1];
rz(3.0016532) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16271834) q[0];
sx q[0];
rz(-0.97467438) q[0];
sx q[0];
rz(0.30124248) q[0];
x q[1];
rz(1.3792737) q[2];
sx q[2];
rz(-0.46326783) q[2];
sx q[2];
rz(-2.6216216) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9819239) q[1];
sx q[1];
rz(-2.2743011) q[1];
sx q[1];
rz(-0.39969201) q[1];
rz(-pi) q[2];
rz(1.1902963) q[3];
sx q[3];
rz(-2.4463013) q[3];
sx q[3];
rz(-2.8075346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.66190019) q[2];
sx q[2];
rz(-2.4883344) q[2];
sx q[2];
rz(-0.24197401) q[2];
rz(0.74609977) q[3];
sx q[3];
rz(-1.3662162) q[3];
sx q[3];
rz(-1.7297176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11874966) q[0];
sx q[0];
rz(-1.043909) q[0];
sx q[0];
rz(-1.2710849) q[0];
rz(-2.5785043) q[1];
sx q[1];
rz(-2.2124898) q[1];
sx q[1];
rz(-1.44106) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25172397) q[0];
sx q[0];
rz(-1.8419203) q[0];
sx q[0];
rz(-2.1878178) q[0];
rz(2.7602642) q[2];
sx q[2];
rz(-1.6994119) q[2];
sx q[2];
rz(0.53024697) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5477834) q[1];
sx q[1];
rz(-1.9873706) q[1];
sx q[1];
rz(1.814247) q[1];
rz(0.7251803) q[3];
sx q[3];
rz(-1.7340875) q[3];
sx q[3];
rz(2.0189328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9219804) q[2];
sx q[2];
rz(-1.8461123) q[2];
sx q[2];
rz(-1.316635) q[2];
rz(-2.5804139) q[3];
sx q[3];
rz(-2.1771274) q[3];
sx q[3];
rz(1.6130028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0372666) q[0];
sx q[0];
rz(-2.9236561) q[0];
sx q[0];
rz(-0.88687801) q[0];
rz(2.9414224) q[1];
sx q[1];
rz(-1.6693516) q[1];
sx q[1];
rz(1.0221457) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0933233) q[0];
sx q[0];
rz(-1.5962068) q[0];
sx q[0];
rz(1.4396458) q[0];
rz(-pi) q[1];
rz(-0.736945) q[2];
sx q[2];
rz(-0.96074694) q[2];
sx q[2];
rz(-1.4184773) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2547438) q[1];
sx q[1];
rz(-1.5690156) q[1];
sx q[1];
rz(-0.82474553) q[1];
rz(1.2194013) q[3];
sx q[3];
rz(-1.0236275) q[3];
sx q[3];
rz(1.6148293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6518121) q[2];
sx q[2];
rz(-2.688789) q[2];
sx q[2];
rz(-0.81857267) q[2];
rz(0.98322785) q[3];
sx q[3];
rz(-0.95663095) q[3];
sx q[3];
rz(2.6035068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.035456903) q[0];
sx q[0];
rz(-1.9945194) q[0];
sx q[0];
rz(1.5863093) q[0];
rz(-2.4447794) q[1];
sx q[1];
rz(-0.15334829) q[1];
sx q[1];
rz(2.3568025) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1479758) q[0];
sx q[0];
rz(-2.412556) q[0];
sx q[0];
rz(0.6887129) q[0];
x q[1];
rz(0.92918877) q[2];
sx q[2];
rz(-0.25114533) q[2];
sx q[2];
rz(-0.26640688) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1363298) q[1];
sx q[1];
rz(-1.6044093) q[1];
sx q[1];
rz(-2.7064347) q[1];
rz(-0.56702964) q[3];
sx q[3];
rz(-2.4483557) q[3];
sx q[3];
rz(2.6854613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4697504) q[2];
sx q[2];
rz(-1.2776813) q[2];
sx q[2];
rz(0.81076852) q[2];
rz(1.9226711) q[3];
sx q[3];
rz(-0.54088497) q[3];
sx q[3];
rz(-1.0651275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3588381) q[0];
sx q[0];
rz(-2.8420119) q[0];
sx q[0];
rz(1.4240356) q[0];
rz(-2.3639823) q[1];
sx q[1];
rz(-0.97336665) q[1];
sx q[1];
rz(-0.75540677) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1799602) q[0];
sx q[0];
rz(-0.20997071) q[0];
sx q[0];
rz(-1.5487973) q[0];
rz(-pi) q[1];
x q[1];
rz(0.39647409) q[2];
sx q[2];
rz(-1.586953) q[2];
sx q[2];
rz(-1.1262058) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.81138071) q[1];
sx q[1];
rz(-2.5054058) q[1];
sx q[1];
rz(2.5274171) q[1];
rz(-pi) q[2];
rz(-1.8978664) q[3];
sx q[3];
rz(-0.9216412) q[3];
sx q[3];
rz(-0.4285194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5648254) q[2];
sx q[2];
rz(-0.29594031) q[2];
sx q[2];
rz(0.15288615) q[2];
rz(-1.8704002) q[3];
sx q[3];
rz(-1.8891687) q[3];
sx q[3];
rz(-0.66711867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2818903) q[0];
sx q[0];
rz(-0.49260456) q[0];
sx q[0];
rz(2.3717666) q[0];
rz(-1.9181171) q[1];
sx q[1];
rz(-1.9203095) q[1];
sx q[1];
rz(-0.65345678) q[1];
rz(-1.1350994) q[2];
sx q[2];
rz(-1.0183327) q[2];
sx q[2];
rz(-0.20295126) q[2];
rz(0.85930227) q[3];
sx q[3];
rz(-0.51325428) q[3];
sx q[3];
rz(-0.2556066) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
