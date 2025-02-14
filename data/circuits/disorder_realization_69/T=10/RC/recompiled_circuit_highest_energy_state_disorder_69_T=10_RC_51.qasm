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
rz(0.82290736) q[0];
sx q[0];
rz(-0.35879254) q[0];
sx q[0];
rz(0.86831492) q[0];
rz(1.7262285) q[1];
sx q[1];
rz(4.2754824) q[1];
sx q[1];
rz(10.418539) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70282763) q[0];
sx q[0];
rz(-1.1701705) q[0];
sx q[0];
rz(-1.6933269) q[0];
rz(-pi) q[1];
rz(2.4325718) q[2];
sx q[2];
rz(-0.62554255) q[2];
sx q[2];
rz(0.65777212) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6958625) q[1];
sx q[1];
rz(-0.63803405) q[1];
sx q[1];
rz(2.3444434) q[1];
rz(-pi) q[2];
rz(2.8959729) q[3];
sx q[3];
rz(-1.1455451) q[3];
sx q[3];
rz(0.58477816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.22780861) q[2];
sx q[2];
rz(-1.3334393) q[2];
sx q[2];
rz(-2.559973) q[2];
rz(0.73451129) q[3];
sx q[3];
rz(-1.651265) q[3];
sx q[3];
rz(-0.41828004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2468579) q[0];
sx q[0];
rz(-1.3792091) q[0];
sx q[0];
rz(2.6439164) q[0];
rz(1.0408164) q[1];
sx q[1];
rz(-0.40319315) q[1];
sx q[1];
rz(1.9042447) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5631249) q[0];
sx q[0];
rz(-1.5284561) q[0];
sx q[0];
rz(2.4189831) q[0];
rz(-pi) q[1];
rz(-1.1038647) q[2];
sx q[2];
rz(-0.40019372) q[2];
sx q[2];
rz(-0.23506842) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.633017) q[1];
sx q[1];
rz(-0.45387065) q[1];
sx q[1];
rz(-1.8967486) q[1];
rz(-pi) q[2];
rz(-0.26555221) q[3];
sx q[3];
rz(-1.7803811) q[3];
sx q[3];
rz(0.095106212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.70025468) q[2];
sx q[2];
rz(-0.73107084) q[2];
sx q[2];
rz(2.8311484) q[2];
rz(1.1566409) q[3];
sx q[3];
rz(-2.0168596) q[3];
sx q[3];
rz(1.1667075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1336841) q[0];
sx q[0];
rz(-2.5569361) q[0];
sx q[0];
rz(0.34580082) q[0];
rz(2.8969104) q[1];
sx q[1];
rz(-0.9318277) q[1];
sx q[1];
rz(-1.7049047) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0416243) q[0];
sx q[0];
rz(-1.9910553) q[0];
sx q[0];
rz(0.34456518) q[0];
rz(-pi) q[1];
x q[1];
rz(2.56836) q[2];
sx q[2];
rz(-2.4893508) q[2];
sx q[2];
rz(2.531372) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.50592283) q[1];
sx q[1];
rz(-1.782592) q[1];
sx q[1];
rz(0.79896547) q[1];
rz(-pi) q[2];
rz(0.47759836) q[3];
sx q[3];
rz(-1.9061889) q[3];
sx q[3];
rz(2.9085577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0077385) q[2];
sx q[2];
rz(-1.1275007) q[2];
sx q[2];
rz(2.5109042) q[2];
rz(0.33356365) q[3];
sx q[3];
rz(-2.1400698) q[3];
sx q[3];
rz(-2.1400616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-3.0665322) q[0];
sx q[0];
rz(-2.1933031) q[0];
sx q[0];
rz(1.4045658) q[0];
rz(0.36918494) q[1];
sx q[1];
rz(-1.4063947) q[1];
sx q[1];
rz(-0.085748347) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2959901) q[0];
sx q[0];
rz(-1.1628502) q[0];
sx q[0];
rz(-1.8970117) q[0];
rz(-pi) q[1];
rz(1.9783114) q[2];
sx q[2];
rz(-1.9980944) q[2];
sx q[2];
rz(-0.80941647) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.71158108) q[1];
sx q[1];
rz(-2.0329355) q[1];
sx q[1];
rz(0.26600809) q[1];
x q[2];
rz(-2.4749807) q[3];
sx q[3];
rz(-0.48891196) q[3];
sx q[3];
rz(2.0036445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.839445) q[2];
sx q[2];
rz(-1.9044694) q[2];
sx q[2];
rz(-0.5298003) q[2];
rz(-0.038330404) q[3];
sx q[3];
rz(-2.4119792) q[3];
sx q[3];
rz(-1.5995601) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2235276) q[0];
sx q[0];
rz(-0.46850884) q[0];
sx q[0];
rz(-1.9388306) q[0];
rz(2.862152) q[1];
sx q[1];
rz(-1.025082) q[1];
sx q[1];
rz(-2.3675945) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27085486) q[0];
sx q[0];
rz(-1.8841198) q[0];
sx q[0];
rz(2.3148651) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4882795) q[2];
sx q[2];
rz(-1.3386209) q[2];
sx q[2];
rz(-1.9407995) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7420089) q[1];
sx q[1];
rz(-3.0019433) q[1];
sx q[1];
rz(1.6634455) q[1];
rz(-pi) q[2];
rz(-0.0027917248) q[3];
sx q[3];
rz(-1.0353902) q[3];
sx q[3];
rz(-1.1100779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1191795) q[2];
sx q[2];
rz(-0.1410307) q[2];
sx q[2];
rz(-0.5640344) q[2];
rz(0.65308475) q[3];
sx q[3];
rz(-1.1354732) q[3];
sx q[3];
rz(-2.691332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7804724) q[0];
sx q[0];
rz(-0.55279624) q[0];
sx q[0];
rz(2.4984388) q[0];
rz(1.9505352) q[1];
sx q[1];
rz(-1.6845614) q[1];
sx q[1];
rz(-0.65470421) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7793286) q[0];
sx q[0];
rz(-1.6342499) q[0];
sx q[0];
rz(0.16880798) q[0];
rz(1.6088609) q[2];
sx q[2];
rz(-0.88721472) q[2];
sx q[2];
rz(0.75929196) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7186942) q[1];
sx q[1];
rz(-2.5942583) q[1];
sx q[1];
rz(1.4446299) q[1];
x q[2];
rz(0.07043802) q[3];
sx q[3];
rz(-1.1813191) q[3];
sx q[3];
rz(2.4088483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5074629) q[2];
sx q[2];
rz(-2.7644988) q[2];
sx q[2];
rz(1.6667574) q[2];
rz(-2.6569488) q[3];
sx q[3];
rz(-2.1438997) q[3];
sx q[3];
rz(1.8528329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6396879) q[0];
sx q[0];
rz(-0.26308331) q[0];
sx q[0];
rz(-0.68341533) q[0];
rz(-3.0112093) q[1];
sx q[1];
rz(-1.5905292) q[1];
sx q[1];
rz(3.108976) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0158247) q[0];
sx q[0];
rz(-0.59249632) q[0];
sx q[0];
rz(1.2529208) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4559385) q[2];
sx q[2];
rz(-1.7259806) q[2];
sx q[2];
rz(-0.48587605) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.72717818) q[1];
sx q[1];
rz(-1.4041931) q[1];
sx q[1];
rz(1.8397306) q[1];
rz(-2.272761) q[3];
sx q[3];
rz(-1.371939) q[3];
sx q[3];
rz(1.5031888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.42328295) q[2];
sx q[2];
rz(-2.0330567) q[2];
sx q[2];
rz(-0.56524593) q[2];
rz(1.7806753) q[3];
sx q[3];
rz(-2.8748685) q[3];
sx q[3];
rz(-0.81418532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77124202) q[0];
sx q[0];
rz(-3.0841565) q[0];
sx q[0];
rz(0.29956079) q[0];
rz(-1.6948505) q[1];
sx q[1];
rz(-2.2694777) q[1];
sx q[1];
rz(-0.95132336) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0438536) q[0];
sx q[0];
rz(-1.7758177) q[0];
sx q[0];
rz(3.1084188) q[0];
x q[1];
rz(-1.6137684) q[2];
sx q[2];
rz(-2.9298721) q[2];
sx q[2];
rz(0.73390244) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7129977) q[1];
sx q[1];
rz(-2.7168437) q[1];
sx q[1];
rz(0.049475706) q[1];
rz(0.087835066) q[3];
sx q[3];
rz(-1.8316275) q[3];
sx q[3];
rz(1.5276599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.99913725) q[2];
sx q[2];
rz(-2.3951525) q[2];
sx q[2];
rz(-2.8738521) q[2];
rz(-2.1382051) q[3];
sx q[3];
rz(-2.181874) q[3];
sx q[3];
rz(-2.3333534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35850152) q[0];
sx q[0];
rz(-2.6434904) q[0];
sx q[0];
rz(2.1844693) q[0];
rz(-0.3793017) q[1];
sx q[1];
rz(-2.7298268) q[1];
sx q[1];
rz(2.9764825) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50277621) q[0];
sx q[0];
rz(-2.2383225) q[0];
sx q[0];
rz(-2.0124042) q[0];
x q[1];
rz(-1.6696641) q[2];
sx q[2];
rz(-0.47438388) q[2];
sx q[2];
rz(2.6883467) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.86805008) q[1];
sx q[1];
rz(-2.1434577) q[1];
sx q[1];
rz(2.9258201) q[1];
rz(-pi) q[2];
rz(1.2312789) q[3];
sx q[3];
rz(-0.63106189) q[3];
sx q[3];
rz(1.5691882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2740606) q[2];
sx q[2];
rz(-1.2313077) q[2];
sx q[2];
rz(0.77862281) q[2];
rz(2.8042931) q[3];
sx q[3];
rz(-1.0236579) q[3];
sx q[3];
rz(-1.5844828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2839277) q[0];
sx q[0];
rz(-1.090467) q[0];
sx q[0];
rz(-2.0528059) q[0];
rz(2.7133443) q[1];
sx q[1];
rz(-2.1183522) q[1];
sx q[1];
rz(2.4931989) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89598362) q[0];
sx q[0];
rz(-1.147384) q[0];
sx q[0];
rz(1.8603252) q[0];
rz(1.4444541) q[2];
sx q[2];
rz(-0.34046945) q[2];
sx q[2];
rz(0.44035092) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9251483) q[1];
sx q[1];
rz(-0.92324644) q[1];
sx q[1];
rz(1.6236033) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5779453) q[3];
sx q[3];
rz(-1.4741338) q[3];
sx q[3];
rz(-3.0954297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.26583656) q[2];
sx q[2];
rz(-2.0529604) q[2];
sx q[2];
rz(-2.9618373) q[2];
rz(-1.1579375) q[3];
sx q[3];
rz(-1.6854743) q[3];
sx q[3];
rz(1.2402844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6801878) q[0];
sx q[0];
rz(-2.9869933) q[0];
sx q[0];
rz(-2.3416478) q[0];
rz(1.1663306) q[1];
sx q[1];
rz(-2.1569398) q[1];
sx q[1];
rz(0.91469761) q[1];
rz(0.88506076) q[2];
sx q[2];
rz(-0.37715465) q[2];
sx q[2];
rz(-2.5436795) q[2];
rz(-1.3553741) q[3];
sx q[3];
rz(-1.2385291) q[3];
sx q[3];
rz(1.877906) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
