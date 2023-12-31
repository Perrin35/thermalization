OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.0467779) q[0];
sx q[0];
rz(-1.0682286) q[0];
sx q[0];
rz(-0.46407035) q[0];
rz(1.9595454) q[1];
sx q[1];
rz(6.2160677) q[1];
sx q[1];
rz(10.709229) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5285437) q[0];
sx q[0];
rz(-1.9595946) q[0];
sx q[0];
rz(-2.8253428) q[0];
rz(-pi) q[1];
rz(3.0797144) q[2];
sx q[2];
rz(-0.49052325) q[2];
sx q[2];
rz(0.86531901) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9570934) q[1];
sx q[1];
rz(-2.543499) q[1];
sx q[1];
rz(-2.4615272) q[1];
rz(-pi) q[2];
rz(-0.10847096) q[3];
sx q[3];
rz(-2.9900108) q[3];
sx q[3];
rz(-0.82419318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7444732) q[2];
sx q[2];
rz(-0.93966547) q[2];
sx q[2];
rz(2.822067) q[2];
rz(2.5630991) q[3];
sx q[3];
rz(-2.6632023) q[3];
sx q[3];
rz(-0.67392504) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4085061) q[0];
sx q[0];
rz(-1.5698313) q[0];
sx q[0];
rz(2.615036) q[0];
rz(2.5780442) q[1];
sx q[1];
rz(-1.5310409) q[1];
sx q[1];
rz(2.3449576) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5986886) q[0];
sx q[0];
rz(-1.2447378) q[0];
sx q[0];
rz(-0.90755264) q[0];
rz(-pi) q[1];
x q[1];
rz(0.30226207) q[2];
sx q[2];
rz(-0.59603359) q[2];
sx q[2];
rz(0.1421393) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9990275) q[1];
sx q[1];
rz(-1.5579281) q[1];
sx q[1];
rz(0.60728118) q[1];
rz(-pi) q[2];
rz(-2.1809686) q[3];
sx q[3];
rz(-1.989813) q[3];
sx q[3];
rz(-0.39715365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.18156302) q[2];
sx q[2];
rz(-1.7627565) q[2];
sx q[2];
rz(-2.3550418) q[2];
rz(-2.6484047) q[3];
sx q[3];
rz(-1.9599873) q[3];
sx q[3];
rz(-2.6942159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2717473) q[0];
sx q[0];
rz(-2.2646876) q[0];
sx q[0];
rz(-1.7012117) q[0];
rz(2.4213743) q[1];
sx q[1];
rz(-2.5455988) q[1];
sx q[1];
rz(0.70297855) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7754606) q[0];
sx q[0];
rz(-1.4538987) q[0];
sx q[0];
rz(1.8525821) q[0];
x q[1];
rz(0.49836068) q[2];
sx q[2];
rz(-1.2549855) q[2];
sx q[2];
rz(-0.74893803) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.30333334) q[1];
sx q[1];
rz(-1.6175744) q[1];
sx q[1];
rz(-0.10951885) q[1];
rz(0.17351563) q[3];
sx q[3];
rz(-0.95458889) q[3];
sx q[3];
rz(-0.71615744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8748223) q[2];
sx q[2];
rz(-1.4845994) q[2];
sx q[2];
rz(2.2300569) q[2];
rz(2.2276145) q[3];
sx q[3];
rz(-1.4303047) q[3];
sx q[3];
rz(-0.82733697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5660969) q[0];
sx q[0];
rz(-0.63593447) q[0];
sx q[0];
rz(0.77600586) q[0];
rz(-1.2674747) q[1];
sx q[1];
rz(-2.0327366) q[1];
sx q[1];
rz(2.5783096) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3455428) q[0];
sx q[0];
rz(-2.4420218) q[0];
sx q[0];
rz(1.1388586) q[0];
x q[1];
rz(1.2345384) q[2];
sx q[2];
rz(-0.90869892) q[2];
sx q[2];
rz(1.3627571) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5717585) q[1];
sx q[1];
rz(-2.8697439) q[1];
sx q[1];
rz(-0.024072577) q[1];
rz(1.2761649) q[3];
sx q[3];
rz(-0.8751103) q[3];
sx q[3];
rz(-2.5254315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.48878601) q[2];
sx q[2];
rz(-0.93549171) q[2];
sx q[2];
rz(-2.9768067) q[2];
rz(2.9131043) q[3];
sx q[3];
rz(-0.27992862) q[3];
sx q[3];
rz(0.94648186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27424681) q[0];
sx q[0];
rz(-2.2773401) q[0];
sx q[0];
rz(-2.2891323) q[0];
rz(0.35119855) q[1];
sx q[1];
rz(-0.57792592) q[1];
sx q[1];
rz(1.16211) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4441372) q[0];
sx q[0];
rz(-1.1778957) q[0];
sx q[0];
rz(-0.59209728) q[0];
rz(-0.94050546) q[2];
sx q[2];
rz(-2.5224707) q[2];
sx q[2];
rz(1.8096015) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.64486849) q[1];
sx q[1];
rz(-2.2779896) q[1];
sx q[1];
rz(1.0134539) q[1];
rz(-3.1245072) q[3];
sx q[3];
rz(-0.72165976) q[3];
sx q[3];
rz(0.34860308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6953485) q[2];
sx q[2];
rz(-2.7480795) q[2];
sx q[2];
rz(1.8943141) q[2];
rz(3.128483) q[3];
sx q[3];
rz(-1.40991) q[3];
sx q[3];
rz(2.9124027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2728249) q[0];
sx q[0];
rz(-1.2478963) q[0];
sx q[0];
rz(2.390958) q[0];
rz(1.1095095) q[1];
sx q[1];
rz(-0.37056071) q[1];
sx q[1];
rz(-1.8355339) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7040946) q[0];
sx q[0];
rz(-1.1878345) q[0];
sx q[0];
rz(-1.0136481) q[0];
rz(-1.7734217) q[2];
sx q[2];
rz(-1.1154004) q[2];
sx q[2];
rz(-2.1855598) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8514511) q[1];
sx q[1];
rz(-0.72039225) q[1];
sx q[1];
rz(-0.011522567) q[1];
rz(-pi) q[2];
rz(-0.2209729) q[3];
sx q[3];
rz(-0.60022012) q[3];
sx q[3];
rz(-0.61629399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3874454) q[2];
sx q[2];
rz(-1.0070609) q[2];
sx q[2];
rz(-2.5409017) q[2];
rz(-2.1389652) q[3];
sx q[3];
rz(-2.6032344) q[3];
sx q[3];
rz(0.35282648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3951185) q[0];
sx q[0];
rz(-1.6858608) q[0];
sx q[0];
rz(2.0948998) q[0];
rz(1.612161) q[1];
sx q[1];
rz(-1.7144831) q[1];
sx q[1];
rz(-0.41710645) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0171623) q[0];
sx q[0];
rz(-0.22148795) q[0];
sx q[0];
rz(1.682196) q[0];
x q[1];
rz(0.76191683) q[2];
sx q[2];
rz(-2.7284405) q[2];
sx q[2];
rz(-0.44943902) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6542146) q[1];
sx q[1];
rz(-2.0301135) q[1];
sx q[1];
rz(2.3801801) q[1];
rz(-pi) q[2];
rz(-1.8700637) q[3];
sx q[3];
rz(-0.93942681) q[3];
sx q[3];
rz(0.33559468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.0059011857) q[2];
sx q[2];
rz(-0.4726755) q[2];
sx q[2];
rz(1.7741514) q[2];
rz(-2.588429) q[3];
sx q[3];
rz(-1.4033214) q[3];
sx q[3];
rz(2.6161391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.816514) q[0];
sx q[0];
rz(-1.8186318) q[0];
sx q[0];
rz(1.1897855) q[0];
rz(1.7143543) q[1];
sx q[1];
rz(-2.863527) q[1];
sx q[1];
rz(3.0292125) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.534879) q[0];
sx q[0];
rz(-1.0977854) q[0];
sx q[0];
rz(0.19434778) q[0];
rz(-1.6582279) q[2];
sx q[2];
rz(-1.9033252) q[2];
sx q[2];
rz(2.33193) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1689414) q[1];
sx q[1];
rz(-0.87247889) q[1];
sx q[1];
rz(-2.6886743) q[1];
rz(0.12179575) q[3];
sx q[3];
rz(-1.313901) q[3];
sx q[3];
rz(1.6459873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0671976) q[2];
sx q[2];
rz(-1.2290359) q[2];
sx q[2];
rz(-1.7929662) q[2];
rz(1.9366692) q[3];
sx q[3];
rz(-1.4068853) q[3];
sx q[3];
rz(0.19781923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7444721) q[0];
sx q[0];
rz(-1.4646126) q[0];
sx q[0];
rz(0.034974139) q[0];
rz(-0.84683013) q[1];
sx q[1];
rz(-1.433082) q[1];
sx q[1];
rz(-2.2299178) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23110403) q[0];
sx q[0];
rz(-1.3470874) q[0];
sx q[0];
rz(1.847812) q[0];
rz(-pi) q[1];
rz(-1.4439092) q[2];
sx q[2];
rz(-1.03089) q[2];
sx q[2];
rz(1.1586231) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1940143) q[1];
sx q[1];
rz(-1.36735) q[1];
sx q[1];
rz(0.46781637) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.53578844) q[3];
sx q[3];
rz(-2.0013323) q[3];
sx q[3];
rz(-2.146194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4336865) q[2];
sx q[2];
rz(-2.6277379) q[2];
sx q[2];
rz(2.8923477) q[2];
rz(-2.3748659) q[3];
sx q[3];
rz(-1.8079575) q[3];
sx q[3];
rz(-2.7419817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7837759) q[0];
sx q[0];
rz(-1.0483085) q[0];
sx q[0];
rz(-2.7695079) q[0];
rz(-2.5601939) q[1];
sx q[1];
rz(-0.77722469) q[1];
sx q[1];
rz(1.4153597) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6045195) q[0];
sx q[0];
rz(-1.3086638) q[0];
sx q[0];
rz(-1.7781236) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4731589) q[2];
sx q[2];
rz(-0.738315) q[2];
sx q[2];
rz(2.7065606) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.806554) q[1];
sx q[1];
rz(-1.5015258) q[1];
sx q[1];
rz(-0.075242234) q[1];
x q[2];
rz(-1.0269208) q[3];
sx q[3];
rz(-1.9034991) q[3];
sx q[3];
rz(1.1247016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8156478) q[2];
sx q[2];
rz(-0.14020136) q[2];
sx q[2];
rz(0.20467219) q[2];
rz(1.4137911) q[3];
sx q[3];
rz(-1.9982343) q[3];
sx q[3];
rz(2.0457101) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6407912) q[0];
sx q[0];
rz(-1.4688107) q[0];
sx q[0];
rz(-2.2039913) q[0];
rz(-1.5564556) q[1];
sx q[1];
rz(-1.7592447) q[1];
sx q[1];
rz(-1.6978282) q[1];
rz(2.0251705) q[2];
sx q[2];
rz(-1.4618256) q[2];
sx q[2];
rz(1.7088919) q[2];
rz(1.6987883) q[3];
sx q[3];
rz(-0.57590579) q[3];
sx q[3];
rz(-0.87344575) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
