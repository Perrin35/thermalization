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
rz(-0.80225575) q[0];
sx q[0];
rz(-1.7576317) q[0];
sx q[0];
rz(1.2686165) q[0];
rz(0.4624548) q[1];
sx q[1];
rz(-3.0825244) q[1];
sx q[1];
rz(1.4055835) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48412886) q[0];
sx q[0];
rz(-1.1384369) q[0];
sx q[0];
rz(-1.0493438) q[0];
rz(-pi) q[1];
rz(1.8008158) q[2];
sx q[2];
rz(-0.91311306) q[2];
sx q[2];
rz(0.44154007) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.51912145) q[1];
sx q[1];
rz(-0.85825413) q[1];
sx q[1];
rz(-2.126241) q[1];
x q[2];
rz(-2.9146092) q[3];
sx q[3];
rz(-2.2344518) q[3];
sx q[3];
rz(0.047078156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1543701) q[2];
sx q[2];
rz(-1.1673678) q[2];
sx q[2];
rz(2.3517189) q[2];
rz(3.1176944) q[3];
sx q[3];
rz(-1.3393211) q[3];
sx q[3];
rz(-0.53643119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8928878) q[0];
sx q[0];
rz(-1.4905812) q[0];
sx q[0];
rz(-1.8875246) q[0];
rz(-1.4887384) q[1];
sx q[1];
rz(-2.2658927) q[1];
sx q[1];
rz(-2.658433) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2452908) q[0];
sx q[0];
rz(-1.5190796) q[0];
sx q[0];
rz(1.381676) q[0];
rz(-pi) q[1];
rz(0.73504059) q[2];
sx q[2];
rz(-1.9695821) q[2];
sx q[2];
rz(2.6533108) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.70139685) q[1];
sx q[1];
rz(-2.6876175) q[1];
sx q[1];
rz(-2.5099436) q[1];
x q[2];
rz(0.29045136) q[3];
sx q[3];
rz(-1.402352) q[3];
sx q[3];
rz(-0.12663933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2529605) q[2];
sx q[2];
rz(-1.6855719) q[2];
sx q[2];
rz(1.6271094) q[2];
rz(0.0557946) q[3];
sx q[3];
rz(-1.5972219) q[3];
sx q[3];
rz(-1.4896711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2273939) q[0];
sx q[0];
rz(-0.44454235) q[0];
sx q[0];
rz(-0.12810531) q[0];
rz(-2.5654492) q[1];
sx q[1];
rz(-0.73155254) q[1];
sx q[1];
rz(-0.76549706) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2554995) q[0];
sx q[0];
rz(-2.1247201) q[0];
sx q[0];
rz(-1.9785247) q[0];
rz(-2.4548495) q[2];
sx q[2];
rz(-1.2740861) q[2];
sx q[2];
rz(2.3627757) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.004359) q[1];
sx q[1];
rz(-2.4437014) q[1];
sx q[1];
rz(-0.76142074) q[1];
rz(-pi) q[2];
rz(-1.5452905) q[3];
sx q[3];
rz(-1.7141986) q[3];
sx q[3];
rz(-2.0633782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3196044) q[2];
sx q[2];
rz(-2.2687948) q[2];
sx q[2];
rz(-0.86307159) q[2];
rz(0.34267628) q[3];
sx q[3];
rz(-0.42049146) q[3];
sx q[3];
rz(-1.4884523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2494025) q[0];
sx q[0];
rz(-0.96977314) q[0];
sx q[0];
rz(2.6595907) q[0];
rz(0.30872289) q[1];
sx q[1];
rz(-0.32544193) q[1];
sx q[1];
rz(2.5293005) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68744874) q[0];
sx q[0];
rz(-2.4061235) q[0];
sx q[0];
rz(0.81540458) q[0];
rz(-pi) q[1];
x q[1];
rz(0.014043645) q[2];
sx q[2];
rz(-1.6545452) q[2];
sx q[2];
rz(-1.3824826) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4548851) q[1];
sx q[1];
rz(-2.5338182) q[1];
sx q[1];
rz(1.6493504) q[1];
rz(-3.0975625) q[3];
sx q[3];
rz(-0.14942871) q[3];
sx q[3];
rz(1.7946515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4838532) q[2];
sx q[2];
rz(-2.190399) q[2];
sx q[2];
rz(1.8335906) q[2];
rz(-0.43404239) q[3];
sx q[3];
rz(-1.7900034) q[3];
sx q[3];
rz(-1.2691809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16172116) q[0];
sx q[0];
rz(-1.5890108) q[0];
sx q[0];
rz(0.27444926) q[0];
rz(1.5212003) q[1];
sx q[1];
rz(-1.043964) q[1];
sx q[1];
rz(-1.257198) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1311943) q[0];
sx q[0];
rz(-1.3290414) q[0];
sx q[0];
rz(-0.25126033) q[0];
rz(-2.1499356) q[2];
sx q[2];
rz(-2.5265387) q[2];
sx q[2];
rz(-1.312254) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3538227) q[1];
sx q[1];
rz(-2.1102043) q[1];
sx q[1];
rz(2.8899075) q[1];
x q[2];
rz(-1.9433504) q[3];
sx q[3];
rz(-0.74001678) q[3];
sx q[3];
rz(-1.7905025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7460798) q[2];
sx q[2];
rz(-1.505625) q[2];
sx q[2];
rz(0.5557605) q[2];
rz(2.3505576) q[3];
sx q[3];
rz(-1.0020703) q[3];
sx q[3];
rz(-1.5033495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24285862) q[0];
sx q[0];
rz(-1.4396311) q[0];
sx q[0];
rz(-0.71204251) q[0];
rz(0.60246077) q[1];
sx q[1];
rz(-2.7280877) q[1];
sx q[1];
rz(0.079040225) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9475381) q[0];
sx q[0];
rz(-1.506516) q[0];
sx q[0];
rz(-3.0966731) q[0];
rz(-pi) q[1];
rz(-3.0661568) q[2];
sx q[2];
rz(-1.4860538) q[2];
sx q[2];
rz(2.2817979) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6106633) q[1];
sx q[1];
rz(-2.3641415) q[1];
sx q[1];
rz(2.8270097) q[1];
rz(-pi) q[2];
rz(-2.8435264) q[3];
sx q[3];
rz(-1.6692999) q[3];
sx q[3];
rz(-0.57386604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1003803) q[2];
sx q[2];
rz(-0.93116394) q[2];
sx q[2];
rz(-2.1616914) q[2];
rz(-1.6117217) q[3];
sx q[3];
rz(-1.2414705) q[3];
sx q[3];
rz(-2.8821168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.74705446) q[0];
sx q[0];
rz(-1.2815481) q[0];
sx q[0];
rz(0.10093149) q[0];
rz(-2.586567) q[1];
sx q[1];
rz(-2.2075768) q[1];
sx q[1];
rz(-1.2947882) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40833631) q[0];
sx q[0];
rz(-2.2332158) q[0];
sx q[0];
rz(0.29710476) q[0];
rz(1.5739977) q[2];
sx q[2];
rz(-0.29433196) q[2];
sx q[2];
rz(-0.41121361) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.50404149) q[1];
sx q[1];
rz(-2.6890272) q[1];
sx q[1];
rz(2.8082147) q[1];
rz(-pi) q[2];
rz(-1.0247308) q[3];
sx q[3];
rz(-0.52123681) q[3];
sx q[3];
rz(-1.2574399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9635705) q[2];
sx q[2];
rz(-0.87248412) q[2];
sx q[2];
rz(-2.8998609) q[2];
rz(-0.47197765) q[3];
sx q[3];
rz(-0.21502544) q[3];
sx q[3];
rz(-1.5836466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80411512) q[0];
sx q[0];
rz(-2.1601456) q[0];
sx q[0];
rz(-0.61019439) q[0];
rz(-2.923851) q[1];
sx q[1];
rz(-2.1861031) q[1];
sx q[1];
rz(0.82966963) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.037704) q[0];
sx q[0];
rz(-1.0173326) q[0];
sx q[0];
rz(-0.82413701) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.83135624) q[2];
sx q[2];
rz(-2.6717253) q[2];
sx q[2];
rz(0.63127764) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.024914537) q[1];
sx q[1];
rz(-2.6905511) q[1];
sx q[1];
rz(-0.43629676) q[1];
rz(-pi) q[2];
rz(1.4368966) q[3];
sx q[3];
rz(-2.4256676) q[3];
sx q[3];
rz(-1.9980824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3024451) q[2];
sx q[2];
rz(-1.382901) q[2];
sx q[2];
rz(2.006532) q[2];
rz(0.48842397) q[3];
sx q[3];
rz(-1.6596183) q[3];
sx q[3];
rz(-1.9058913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9957073) q[0];
sx q[0];
rz(-1.1111525) q[0];
sx q[0];
rz(-2.1671894) q[0];
rz(-2.3459332) q[1];
sx q[1];
rz(-0.37745044) q[1];
sx q[1];
rz(-0.61765751) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5159588) q[0];
sx q[0];
rz(-3.0581116) q[0];
sx q[0];
rz(-0.030103695) q[0];
x q[1];
rz(1.1716003) q[2];
sx q[2];
rz(-2.0621057) q[2];
sx q[2];
rz(1.0688865) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5292735) q[1];
sx q[1];
rz(-0.63748432) q[1];
sx q[1];
rz(0.46302621) q[1];
rz(-pi) q[2];
rz(-2.1626804) q[3];
sx q[3];
rz(-2.5911461) q[3];
sx q[3];
rz(-1.5959671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.016102942) q[2];
sx q[2];
rz(-1.2697271) q[2];
sx q[2];
rz(1.3886836) q[2];
rz(1.3056508) q[3];
sx q[3];
rz(-0.7235705) q[3];
sx q[3];
rz(1.8587662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0131123) q[0];
sx q[0];
rz(-0.73902577) q[0];
sx q[0];
rz(-2.541743) q[0];
rz(-2.5685617) q[1];
sx q[1];
rz(-2.532798) q[1];
sx q[1];
rz(-1.0240239) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1896418) q[0];
sx q[0];
rz(-1.1770403) q[0];
sx q[0];
rz(1.1313361) q[0];
rz(-2.4889718) q[2];
sx q[2];
rz(-1.5372477) q[2];
sx q[2];
rz(1.6217569) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3249859) q[1];
sx q[1];
rz(-1.6546975) q[1];
sx q[1];
rz(2.744426) q[1];
x q[2];
rz(-1.4198331) q[3];
sx q[3];
rz(-1.9486041) q[3];
sx q[3];
rz(2.092416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.92274252) q[2];
sx q[2];
rz(-1.3877733) q[2];
sx q[2];
rz(-2.3663523) q[2];
rz(1.5357337) q[3];
sx q[3];
rz(-1.9550867) q[3];
sx q[3];
rz(-1.8839914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6105462) q[0];
sx q[0];
rz(-0.95120593) q[0];
sx q[0];
rz(-1.9134941) q[0];
rz(1.7465406) q[1];
sx q[1];
rz(-1.4735305) q[1];
sx q[1];
rz(2.7458618) q[1];
rz(1.878123) q[2];
sx q[2];
rz(-2.0057445) q[2];
sx q[2];
rz(2.4705171) q[2];
rz(3.1093194) q[3];
sx q[3];
rz(-2.3393031) q[3];
sx q[3];
rz(2.5325251) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
