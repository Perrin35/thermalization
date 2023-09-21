OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5151514) q[0];
sx q[0];
rz(-0.03349537) q[0];
sx q[0];
rz(-1.3666231) q[0];
rz(2.0904436) q[1];
sx q[1];
rz(-1.6519974) q[1];
sx q[1];
rz(-1.1319914) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53349797) q[0];
sx q[0];
rz(-2.29288) q[0];
sx q[0];
rz(0.46790926) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7366473) q[2];
sx q[2];
rz(-1.4115141) q[2];
sx q[2];
rz(-2.7878891) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1177897) q[1];
sx q[1];
rz(-2.1315247) q[1];
sx q[1];
rz(-2.6610713) q[1];
rz(-pi) q[2];
rz(-0.55478401) q[3];
sx q[3];
rz(-1.3265508) q[3];
sx q[3];
rz(0.82575219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2312317) q[2];
sx q[2];
rz(-1.8138764) q[2];
sx q[2];
rz(1.988391) q[2];
rz(-0.48405805) q[3];
sx q[3];
rz(-2.3279133) q[3];
sx q[3];
rz(-0.4593862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3020878) q[0];
sx q[0];
rz(-2.8110101) q[0];
sx q[0];
rz(-2.6385345) q[0];
rz(-1.5548276) q[1];
sx q[1];
rz(-2.4350872) q[1];
sx q[1];
rz(-0.15393004) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.670823) q[0];
sx q[0];
rz(-1.5851067) q[0];
sx q[0];
rz(-2.0736681) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.33090654) q[2];
sx q[2];
rz(-2.1936596) q[2];
sx q[2];
rz(1.9493584) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1068374) q[1];
sx q[1];
rz(-1.991786) q[1];
sx q[1];
rz(-2.0799328) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5699584) q[3];
sx q[3];
rz(-1.0610233) q[3];
sx q[3];
rz(1.2050932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2543891) q[2];
sx q[2];
rz(-2.7837191) q[2];
sx q[2];
rz(1.3228234) q[2];
rz(-1.6555188) q[3];
sx q[3];
rz(-1.6008987) q[3];
sx q[3];
rz(-2.4310908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94674295) q[0];
sx q[0];
rz(-2.0236334) q[0];
sx q[0];
rz(0.90993607) q[0];
rz(2.3643156) q[1];
sx q[1];
rz(-2.3060019) q[1];
sx q[1];
rz(-0.98532239) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9604208) q[0];
sx q[0];
rz(-0.67502484) q[0];
sx q[0];
rz(1.3135364) q[0];
rz(-pi) q[1];
x q[1];
rz(0.3243685) q[2];
sx q[2];
rz(-2.2019221) q[2];
sx q[2];
rz(0.98482519) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2376155) q[1];
sx q[1];
rz(-1.3217889) q[1];
sx q[1];
rz(-1.6275089) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4661029) q[3];
sx q[3];
rz(-2.1152861) q[3];
sx q[3];
rz(-2.4558223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.862792) q[2];
sx q[2];
rz(-1.1547487) q[2];
sx q[2];
rz(1.2476236) q[2];
rz(-0.4425846) q[3];
sx q[3];
rz(-1.7740039) q[3];
sx q[3];
rz(2.8201593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9001532) q[0];
sx q[0];
rz(-2.1223919) q[0];
sx q[0];
rz(-2.0181657) q[0];
rz(-0.57888794) q[1];
sx q[1];
rz(-1.4588979) q[1];
sx q[1];
rz(-1.7480063) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1136365) q[0];
sx q[0];
rz(-0.96158577) q[0];
sx q[0];
rz(-2.0834126) q[0];
x q[1];
rz(-1.7384999) q[2];
sx q[2];
rz(-1.1412914) q[2];
sx q[2];
rz(1.6657366) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.77809282) q[1];
sx q[1];
rz(-2.2908195) q[1];
sx q[1];
rz(-2.5219445) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3062069) q[3];
sx q[3];
rz(-1.3961892) q[3];
sx q[3];
rz(1.1268827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0397772) q[2];
sx q[2];
rz(-1.5559876) q[2];
sx q[2];
rz(3.139479) q[2];
rz(2.5801616) q[3];
sx q[3];
rz(-2.1241472) q[3];
sx q[3];
rz(0.58825618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.026022) q[0];
sx q[0];
rz(-2.0176812) q[0];
sx q[0];
rz(1.9157238) q[0];
rz(-1.4670124) q[1];
sx q[1];
rz(-1.2743228) q[1];
sx q[1];
rz(1.3668758) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0923094) q[0];
sx q[0];
rz(-1.7607848) q[0];
sx q[0];
rz(-2.7490194) q[0];
rz(-pi) q[1];
rz(-0.032614313) q[2];
sx q[2];
rz(-0.74748749) q[2];
sx q[2];
rz(-0.05687296) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.83134507) q[1];
sx q[1];
rz(-2.275895) q[1];
sx q[1];
rz(-3.0369333) q[1];
rz(-pi) q[2];
rz(-1.2138052) q[3];
sx q[3];
rz(-1.095466) q[3];
sx q[3];
rz(0.19015042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2772969) q[2];
sx q[2];
rz(-2.0531451) q[2];
sx q[2];
rz(-0.40536353) q[2];
rz(-2.6799485) q[3];
sx q[3];
rz(-0.82563892) q[3];
sx q[3];
rz(-1.595114) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6500403) q[0];
sx q[0];
rz(-1.7164282) q[0];
sx q[0];
rz(-0.50338411) q[0];
rz(0.21884306) q[1];
sx q[1];
rz(-1.2501406) q[1];
sx q[1];
rz(0.65178451) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8021278) q[0];
sx q[0];
rz(-1.8101242) q[0];
sx q[0];
rz(1.7665187) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6794372) q[2];
sx q[2];
rz(-0.75192736) q[2];
sx q[2];
rz(1.3736563) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.3411322) q[1];
sx q[1];
rz(-0.144185) q[1];
sx q[1];
rz(1.9393117) q[1];
rz(-1.4676474) q[3];
sx q[3];
rz(-1.2892937) q[3];
sx q[3];
rz(0.75479773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6289604) q[2];
sx q[2];
rz(-1.3971389) q[2];
sx q[2];
rz(2.7499278) q[2];
rz(0.20646778) q[3];
sx q[3];
rz(-2.4075017) q[3];
sx q[3];
rz(-0.68968836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7713292) q[0];
sx q[0];
rz(-1.6917133) q[0];
sx q[0];
rz(0.96965924) q[0];
rz(0.58352739) q[1];
sx q[1];
rz(-2.0136166) q[1];
sx q[1];
rz(-3.0113509) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42338615) q[0];
sx q[0];
rz(-1.5209746) q[0];
sx q[0];
rz(0.35276619) q[0];
rz(-2.4728751) q[2];
sx q[2];
rz(-1.0895184) q[2];
sx q[2];
rz(-0.95578335) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0252467) q[1];
sx q[1];
rz(-1.7882008) q[1];
sx q[1];
rz(-0.51699637) q[1];
rz(-pi) q[2];
rz(2.1533265) q[3];
sx q[3];
rz(-1.8997972) q[3];
sx q[3];
rz(-2.695801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6234201) q[2];
sx q[2];
rz(-1.5600081) q[2];
sx q[2];
rz(1.0858034) q[2];
rz(0.052224934) q[3];
sx q[3];
rz(-1.6970535) q[3];
sx q[3];
rz(-1.0857371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6770342) q[0];
sx q[0];
rz(-2.1573986) q[0];
sx q[0];
rz(-1.1664671) q[0];
rz(-0.72558609) q[1];
sx q[1];
rz(-1.8478994) q[1];
sx q[1];
rz(0.74277791) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22314534) q[0];
sx q[0];
rz(-0.76876516) q[0];
sx q[0];
rz(0.2785152) q[0];
rz(-pi) q[1];
rz(-2.5119945) q[2];
sx q[2];
rz(-0.78231914) q[2];
sx q[2];
rz(0.19537374) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2837977) q[1];
sx q[1];
rz(-0.89110121) q[1];
sx q[1];
rz(-0.53417511) q[1];
rz(-2.2320896) q[3];
sx q[3];
rz(-1.0854183) q[3];
sx q[3];
rz(-0.73736008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5143738) q[2];
sx q[2];
rz(-1.5292995) q[2];
sx q[2];
rz(-2.880704) q[2];
rz(-2.0339113) q[3];
sx q[3];
rz(-1.2872144) q[3];
sx q[3];
rz(-1.0816983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35659197) q[0];
sx q[0];
rz(-0.78105015) q[0];
sx q[0];
rz(-0.93604952) q[0];
rz(0.014135663) q[1];
sx q[1];
rz(-1.3052669) q[1];
sx q[1];
rz(2.4900808) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12670853) q[0];
sx q[0];
rz(-1.5471022) q[0];
sx q[0];
rz(-3.1372877) q[0];
rz(-pi) q[1];
x q[1];
rz(0.65666171) q[2];
sx q[2];
rz(-2.0920394) q[2];
sx q[2];
rz(-1.9929287) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.41496719) q[1];
sx q[1];
rz(-1.8888998) q[1];
sx q[1];
rz(3.0148562) q[1];
x q[2];
rz(-1.9367847) q[3];
sx q[3];
rz(-1.9848739) q[3];
sx q[3];
rz(0.28596349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8462048) q[2];
sx q[2];
rz(-0.69869763) q[2];
sx q[2];
rz(-0.52337581) q[2];
rz(0.30803099) q[3];
sx q[3];
rz(-2.568646) q[3];
sx q[3];
rz(1.3380922) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.734252) q[0];
sx q[0];
rz(-0.14695209) q[0];
sx q[0];
rz(1.7379606) q[0];
rz(2.667528) q[1];
sx q[1];
rz(-1.6876551) q[1];
sx q[1];
rz(1.3778936) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5869851) q[0];
sx q[0];
rz(-0.73903144) q[0];
sx q[0];
rz(-0.87297312) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3465967) q[2];
sx q[2];
rz(-1.0161875) q[2];
sx q[2];
rz(-2.5650052) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.024152012) q[1];
sx q[1];
rz(-1.9630868) q[1];
sx q[1];
rz(1.8333608) q[1];
rz(-pi) q[2];
rz(1.7337448) q[3];
sx q[3];
rz(-1.019422) q[3];
sx q[3];
rz(-2.5577302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.750981) q[2];
sx q[2];
rz(-1.5079974) q[2];
sx q[2];
rz(-1.5314468) q[2];
rz(-1.4577929) q[3];
sx q[3];
rz(-2.0861574) q[3];
sx q[3];
rz(-0.84993258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4703341) q[0];
sx q[0];
rz(-2.8699734) q[0];
sx q[0];
rz(1.097453) q[0];
rz(0.63256565) q[1];
sx q[1];
rz(-1.0610776) q[1];
sx q[1];
rz(-2.9656596) q[1];
rz(1.2420281) q[2];
sx q[2];
rz(-0.9267926) q[2];
sx q[2];
rz(-0.37315858) q[2];
rz(-2.2652763) q[3];
sx q[3];
rz(-2.1764168) q[3];
sx q[3];
rz(-0.62302667) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
