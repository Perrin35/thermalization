OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.87941909) q[0];
sx q[0];
rz(4.8882422) q[0];
sx q[0];
rz(12.56765) q[0];
rz(-1.6969504) q[1];
sx q[1];
rz(-2.0386219) q[1];
sx q[1];
rz(-2.3666518) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89844184) q[0];
sx q[0];
rz(-1.713322) q[0];
sx q[0];
rz(1.5205163) q[0];
rz(-2.2968729) q[2];
sx q[2];
rz(-1.2436927) q[2];
sx q[2];
rz(-0.5288045) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.70667627) q[1];
sx q[1];
rz(-1.7129363) q[1];
sx q[1];
rz(-3.1013156) q[1];
x q[2];
rz(-0.28489057) q[3];
sx q[3];
rz(-1.7171515) q[3];
sx q[3];
rz(-0.26998664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0960192) q[2];
sx q[2];
rz(-1.1117659) q[2];
sx q[2];
rz(1.1958896) q[2];
rz(-1.9879509) q[3];
sx q[3];
rz(-0.78912815) q[3];
sx q[3];
rz(-1.7529863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4202704) q[0];
sx q[0];
rz(-0.35636154) q[0];
sx q[0];
rz(-1.8700245) q[0];
rz(1.0999854) q[1];
sx q[1];
rz(-2.0309235) q[1];
sx q[1];
rz(1.3756479) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.401706) q[0];
sx q[0];
rz(-1.6580083) q[0];
sx q[0];
rz(-2.7313822) q[0];
rz(0.045250821) q[2];
sx q[2];
rz(-1.3327193) q[2];
sx q[2];
rz(-2.0766052) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.83982044) q[1];
sx q[1];
rz(-0.41026792) q[1];
sx q[1];
rz(-2.6778585) q[1];
rz(-1.0248915) q[3];
sx q[3];
rz(-2.5149269) q[3];
sx q[3];
rz(-2.6729667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9282844) q[2];
sx q[2];
rz(-1.3457315) q[2];
sx q[2];
rz(2.6518872) q[2];
rz(1.1335763) q[3];
sx q[3];
rz(-2.9462892) q[3];
sx q[3];
rz(1.0666696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5415444) q[0];
sx q[0];
rz(-1.8627889) q[0];
sx q[0];
rz(2.9009853) q[0];
rz(2.799017) q[1];
sx q[1];
rz(-0.97476417) q[1];
sx q[1];
rz(1.2352357) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2152527) q[0];
sx q[0];
rz(-1.0010166) q[0];
sx q[0];
rz(2.0195872) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8232934) q[2];
sx q[2];
rz(-1.2090948) q[2];
sx q[2];
rz(-2.8351438) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0864799) q[1];
sx q[1];
rz(-1.3172473) q[1];
sx q[1];
rz(-2.8716645) q[1];
x q[2];
rz(-1.0838305) q[3];
sx q[3];
rz(-2.1420797) q[3];
sx q[3];
rz(-0.42698241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3774595) q[2];
sx q[2];
rz(-1.6230134) q[2];
sx q[2];
rz(-1.4366478) q[2];
rz(1.4012339) q[3];
sx q[3];
rz(-1.876372) q[3];
sx q[3];
rz(0.60825545) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.065598) q[0];
sx q[0];
rz(-1.5290715) q[0];
sx q[0];
rz(2.3377989) q[0];
rz(2.1919788) q[1];
sx q[1];
rz(-1.6751553) q[1];
sx q[1];
rz(-0.11985699) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6896497) q[0];
sx q[0];
rz(-1.6065238) q[0];
sx q[0];
rz(-2.1769051) q[0];
rz(-0.36074952) q[2];
sx q[2];
rz(-1.4348239) q[2];
sx q[2];
rz(3.0038358) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4704628) q[1];
sx q[1];
rz(-1.5531335) q[1];
sx q[1];
rz(1.3385685) q[1];
x q[2];
rz(-2.5424764) q[3];
sx q[3];
rz(-1.1154419) q[3];
sx q[3];
rz(1.276254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5084761) q[2];
sx q[2];
rz(-1.8957596) q[2];
sx q[2];
rz(3.0692549) q[2];
rz(2.7667601) q[3];
sx q[3];
rz(-0.65632498) q[3];
sx q[3];
rz(1.1981296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6937834) q[0];
sx q[0];
rz(-2.5700975) q[0];
sx q[0];
rz(0.48686349) q[0];
rz(-0.72987366) q[1];
sx q[1];
rz(-0.90881538) q[1];
sx q[1];
rz(1.240085) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86868984) q[0];
sx q[0];
rz(-1.5839974) q[0];
sx q[0];
rz(-2.9024283) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.24057062) q[2];
sx q[2];
rz(-1.3153937) q[2];
sx q[2];
rz(2.0672928) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.33927321) q[1];
sx q[1];
rz(-2.5335651) q[1];
sx q[1];
rz(0.1165216) q[1];
rz(-pi) q[2];
rz(2.6579882) q[3];
sx q[3];
rz(-2.1198453) q[3];
sx q[3];
rz(-2.1476114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1029677) q[2];
sx q[2];
rz(-0.90927783) q[2];
sx q[2];
rz(0.70303482) q[2];
rz(-1.4098343) q[3];
sx q[3];
rz(-1.792428) q[3];
sx q[3];
rz(-0.15771244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1773961) q[0];
sx q[0];
rz(-1.8494158) q[0];
sx q[0];
rz(0.99037209) q[0];
rz(3.0888427) q[1];
sx q[1];
rz(-2.2149448) q[1];
sx q[1];
rz(1.4809158) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4617417) q[0];
sx q[0];
rz(-1.4977507) q[0];
sx q[0];
rz(1.9528271) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0159857) q[2];
sx q[2];
rz(-1.5700969) q[2];
sx q[2];
rz(2.6631151) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.78696886) q[1];
sx q[1];
rz(-2.0957392) q[1];
sx q[1];
rz(-0.24137361) q[1];
rz(-pi) q[2];
rz(-1.6310286) q[3];
sx q[3];
rz(-2.4619953) q[3];
sx q[3];
rz(-2.4356902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.0084373077) q[2];
sx q[2];
rz(-1.6494273) q[2];
sx q[2];
rz(-0.56224242) q[2];
rz(-2.0810614) q[3];
sx q[3];
rz(-0.74354592) q[3];
sx q[3];
rz(-2.8806768) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9437207) q[0];
sx q[0];
rz(-1.7892388) q[0];
sx q[0];
rz(-1.4690171) q[0];
rz(-1.0143657) q[1];
sx q[1];
rz(-1.0275774) q[1];
sx q[1];
rz(1.3247103) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4949188) q[0];
sx q[0];
rz(-3.0122628) q[0];
sx q[0];
rz(-0.17636756) q[0];
x q[1];
rz(0.84682805) q[2];
sx q[2];
rz(-2.4412051) q[2];
sx q[2];
rz(0.91947739) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.054338) q[1];
sx q[1];
rz(-1.9458276) q[1];
sx q[1];
rz(0.43941984) q[1];
x q[2];
rz(-2.0539829) q[3];
sx q[3];
rz(-2.4251068) q[3];
sx q[3];
rz(2.0487752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5801195) q[2];
sx q[2];
rz(-1.7529528) q[2];
sx q[2];
rz(2.6980147) q[2];
rz(-0.94868547) q[3];
sx q[3];
rz(-1.1921927) q[3];
sx q[3];
rz(-2.366812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7396486) q[0];
sx q[0];
rz(-2.1304603) q[0];
sx q[0];
rz(-0.46052128) q[0];
rz(-3.0415688) q[1];
sx q[1];
rz(-0.99383751) q[1];
sx q[1];
rz(-1.8519648) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3562718) q[0];
sx q[0];
rz(-0.7542146) q[0];
sx q[0];
rz(2.450599) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4256267) q[2];
sx q[2];
rz(-0.95574524) q[2];
sx q[2];
rz(1.8383319) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.76028484) q[1];
sx q[1];
rz(-0.21796255) q[1];
sx q[1];
rz(2.1078307) q[1];
rz(-pi) q[2];
rz(2.9309978) q[3];
sx q[3];
rz(-0.67376332) q[3];
sx q[3];
rz(-1.3099561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.148968) q[2];
sx q[2];
rz(-2.830539) q[2];
sx q[2];
rz(-1.0827433) q[2];
rz(-0.096171245) q[3];
sx q[3];
rz(-1.440719) q[3];
sx q[3];
rz(1.910803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37725317) q[0];
sx q[0];
rz(-0.23582533) q[0];
sx q[0];
rz(2.8073231) q[0];
rz(-1.9175247) q[1];
sx q[1];
rz(-1.5609488) q[1];
sx q[1];
rz(-0.28265488) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5196913) q[0];
sx q[0];
rz(-2.0265409) q[0];
sx q[0];
rz(-0.68463188) q[0];
rz(-pi) q[1];
rz(-0.33090584) q[2];
sx q[2];
rz(-1.7754284) q[2];
sx q[2];
rz(-2.3757039) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4947195) q[1];
sx q[1];
rz(-0.44534007) q[1];
sx q[1];
rz(2.7237707) q[1];
rz(-pi) q[2];
x q[2];
rz(0.88528021) q[3];
sx q[3];
rz(-1.6618528) q[3];
sx q[3];
rz(-2.1212999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9294372) q[2];
sx q[2];
rz(-1.9781457) q[2];
sx q[2];
rz(-0.7412509) q[2];
rz(0.50179982) q[3];
sx q[3];
rz(-1.3391677) q[3];
sx q[3];
rz(-0.58399502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.042645) q[0];
sx q[0];
rz(-2.2462923) q[0];
sx q[0];
rz(2.6328971) q[0];
rz(-0.11518654) q[1];
sx q[1];
rz(-0.70786628) q[1];
sx q[1];
rz(-0.68181109) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9671229) q[0];
sx q[0];
rz(-0.57708626) q[0];
sx q[0];
rz(-1.8814927) q[0];
rz(2.4585312) q[2];
sx q[2];
rz(-0.94711727) q[2];
sx q[2];
rz(-1.9649486) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0360003) q[1];
sx q[1];
rz(-0.75165527) q[1];
sx q[1];
rz(2.4113301) q[1];
x q[2];
rz(-2.3144249) q[3];
sx q[3];
rz(-0.70561545) q[3];
sx q[3];
rz(-0.08882113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.29601413) q[2];
sx q[2];
rz(-1.5073551) q[2];
sx q[2];
rz(-0.34109035) q[2];
rz(-1.0836481) q[3];
sx q[3];
rz(-0.79569474) q[3];
sx q[3];
rz(-2.9705689) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9733799) q[0];
sx q[0];
rz(-1.4773049) q[0];
sx q[0];
rz(2.1422577) q[0];
rz(-0.60733168) q[1];
sx q[1];
rz(-2.2317531) q[1];
sx q[1];
rz(-2.9324525) q[1];
rz(1.0495937) q[2];
sx q[2];
rz(-0.97432077) q[2];
sx q[2];
rz(-1.5885098) q[2];
rz(2.7568983) q[3];
sx q[3];
rz(-1.7582498) q[3];
sx q[3];
rz(2.869217) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
