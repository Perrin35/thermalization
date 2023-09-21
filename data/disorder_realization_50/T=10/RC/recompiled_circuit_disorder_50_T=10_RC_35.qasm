OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2621736) q[0];
sx q[0];
rz(-1.7466495) q[0];
sx q[0];
rz(3.1403132) q[0];
rz(-1.6969504) q[1];
sx q[1];
rz(4.2445634) q[1];
sx q[1];
rz(7.0581262) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2389195) q[0];
sx q[0];
rz(-0.15107778) q[0];
sx q[0];
rz(-2.8047049) q[0];
x q[1];
rz(-1.0983659) q[2];
sx q[2];
rz(-0.78394267) q[2];
sx q[2];
rz(2.4468165) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.85841132) q[1];
sx q[1];
rz(-1.5309257) q[1];
sx q[1];
rz(1.7130501) q[1];
rz(-pi) q[2];
rz(1.7232056) q[3];
sx q[3];
rz(-1.2890352) q[3];
sx q[3];
rz(1.8834653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0455735) q[2];
sx q[2];
rz(-1.1117659) q[2];
sx q[2];
rz(-1.1958896) q[2];
rz(1.9879509) q[3];
sx q[3];
rz(-2.3524645) q[3];
sx q[3];
rz(1.3886064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7213223) q[0];
sx q[0];
rz(-2.7852311) q[0];
sx q[0];
rz(-1.2715682) q[0];
rz(-2.0416073) q[1];
sx q[1];
rz(-1.1106691) q[1];
sx q[1];
rz(1.7659448) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028582024) q[0];
sx q[0];
rz(-2.7227289) q[0];
sx q[0];
rz(-2.9257665) q[0];
rz(1.3865115) q[2];
sx q[2];
rz(-2.8993336) q[2];
sx q[2];
rz(-0.87528961) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.34054204) q[1];
sx q[1];
rz(-1.935563) q[1];
sx q[1];
rz(1.7629452) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7820285) q[3];
sx q[3];
rz(-1.0457977) q[3];
sx q[3];
rz(1.1121225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.21330825) q[2];
sx q[2];
rz(-1.3457315) q[2];
sx q[2];
rz(-0.48970547) q[2];
rz(-2.0080163) q[3];
sx q[3];
rz(-0.19530345) q[3];
sx q[3];
rz(2.074923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5415444) q[0];
sx q[0];
rz(-1.8627889) q[0];
sx q[0];
rz(0.24060732) q[0];
rz(2.799017) q[1];
sx q[1];
rz(-2.1668285) q[1];
sx q[1];
rz(1.906357) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19757195) q[0];
sx q[0];
rz(-0.70957843) q[0];
sx q[0];
rz(2.5463085) q[0];
x q[1];
rz(1.8232934) q[2];
sx q[2];
rz(-1.2090948) q[2];
sx q[2];
rz(-2.8351438) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0864799) q[1];
sx q[1];
rz(-1.8243454) q[1];
sx q[1];
rz(-2.8716645) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.62883212) q[3];
sx q[3];
rz(-1.1662081) q[3];
sx q[3];
rz(-1.42266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3774595) q[2];
sx q[2];
rz(-1.6230134) q[2];
sx q[2];
rz(-1.7049449) q[2];
rz(-1.4012339) q[3];
sx q[3];
rz(-1.876372) q[3];
sx q[3];
rz(2.5333372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.075994611) q[0];
sx q[0];
rz(-1.5290715) q[0];
sx q[0];
rz(2.3377989) q[0];
rz(-0.94961387) q[1];
sx q[1];
rz(-1.4664374) q[1];
sx q[1];
rz(-3.0217357) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45194295) q[0];
sx q[0];
rz(-1.6065238) q[0];
sx q[0];
rz(2.1769051) q[0];
x q[1];
rz(-0.36074952) q[2];
sx q[2];
rz(-1.7067688) q[2];
sx q[2];
rz(-3.0038358) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4704628) q[1];
sx q[1];
rz(-1.5531335) q[1];
sx q[1];
rz(1.8030241) q[1];
rz(-0.5991163) q[3];
sx q[3];
rz(-1.1154419) q[3];
sx q[3];
rz(-1.276254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.63311657) q[2];
sx q[2];
rz(-1.2458331) q[2];
sx q[2];
rz(-0.072337739) q[2];
rz(-2.7667601) q[3];
sx q[3];
rz(-2.4852677) q[3];
sx q[3];
rz(1.1981296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4478093) q[0];
sx q[0];
rz(-0.57149514) q[0];
sx q[0];
rz(-2.6547292) q[0];
rz(0.72987366) q[1];
sx q[1];
rz(-0.90881538) q[1];
sx q[1];
rz(1.9015076) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70532521) q[0];
sx q[0];
rz(-1.3316532) q[0];
sx q[0];
rz(-1.5843841) q[0];
x q[1];
rz(2.3104722) q[2];
sx q[2];
rz(-0.34905012) q[2];
sx q[2];
rz(2.8380053) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6606635) q[1];
sx q[1];
rz(-2.1741121) q[1];
sx q[1];
rz(-1.4900581) q[1];
rz(-pi) q[2];
rz(-2.6579882) q[3];
sx q[3];
rz(-2.1198453) q[3];
sx q[3];
rz(2.1476114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1029677) q[2];
sx q[2];
rz(-2.2323148) q[2];
sx q[2];
rz(0.70303482) q[2];
rz(-1.7317584) q[3];
sx q[3];
rz(-1.3491646) q[3];
sx q[3];
rz(-0.15771244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96419656) q[0];
sx q[0];
rz(-1.2921768) q[0];
sx q[0];
rz(2.1512206) q[0];
rz(3.0888427) q[1];
sx q[1];
rz(-2.2149448) q[1];
sx q[1];
rz(1.4809158) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92025837) q[0];
sx q[0];
rz(-1.9517559) q[0];
sx q[0];
rz(-0.078698054) q[0];
x q[1];
rz(-1.5724206) q[2];
sx q[2];
rz(-0.44518984) q[2];
sx q[2];
rz(1.0908529) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.78696886) q[1];
sx q[1];
rz(-2.0957392) q[1];
sx q[1];
rz(-0.24137361) q[1];
x q[2];
rz(-3.092993) q[3];
sx q[3];
rz(-0.89266333) q[3];
sx q[3];
rz(2.5130659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.1331553) q[2];
sx q[2];
rz(-1.6494273) q[2];
sx q[2];
rz(-0.56224242) q[2];
rz(2.0810614) q[3];
sx q[3];
rz(-2.3980467) q[3];
sx q[3];
rz(0.26091584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19787191) q[0];
sx q[0];
rz(-1.7892388) q[0];
sx q[0];
rz(-1.6725756) q[0];
rz(2.127227) q[1];
sx q[1];
rz(-1.0275774) q[1];
sx q[1];
rz(-1.8168824) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64667386) q[0];
sx q[0];
rz(-3.0122628) q[0];
sx q[0];
rz(0.17636756) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1340738) q[2];
sx q[2];
rz(-1.129732) q[2];
sx q[2];
rz(-1.2457459) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3130256) q[1];
sx q[1];
rz(-1.9777858) q[1];
sx q[1];
rz(1.160497) q[1];
rz(-pi) q[2];
rz(-0.91388254) q[3];
sx q[3];
rz(-1.2607288) q[3];
sx q[3];
rz(0.1012181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5801195) q[2];
sx q[2];
rz(-1.3886398) q[2];
sx q[2];
rz(2.6980147) q[2];
rz(-0.94868547) q[3];
sx q[3];
rz(-1.9493999) q[3];
sx q[3];
rz(-0.77478066) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.401944) q[0];
sx q[0];
rz(-1.0111324) q[0];
sx q[0];
rz(0.46052128) q[0];
rz(-0.1000239) q[1];
sx q[1];
rz(-0.99383751) q[1];
sx q[1];
rz(1.8519648) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.063110654) q[0];
sx q[0];
rz(-1.0149628) q[0];
sx q[0];
rz(-1.0312992) q[0];
x q[1];
rz(2.4256267) q[2];
sx q[2];
rz(-0.95574524) q[2];
sx q[2];
rz(-1.3032608) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3370812) q[1];
sx q[1];
rz(-1.6816499) q[1];
sx q[1];
rz(-1.3827419) q[1];
rz(1.7361705) q[3];
sx q[3];
rz(-2.227042) q[3];
sx q[3];
rz(-2.0986433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.9926247) q[2];
sx q[2];
rz(-2.830539) q[2];
sx q[2];
rz(1.0827433) q[2];
rz(0.096171245) q[3];
sx q[3];
rz(-1.7008737) q[3];
sx q[3];
rz(-1.2307897) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37725317) q[0];
sx q[0];
rz(-0.23582533) q[0];
sx q[0];
rz(2.8073231) q[0];
rz(1.224068) q[1];
sx q[1];
rz(-1.5806438) q[1];
sx q[1];
rz(0.28265488) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6982272) q[0];
sx q[0];
rz(-0.8015612) q[0];
sx q[0];
rz(-0.65936868) q[0];
rz(2.5731509) q[2];
sx q[2];
rz(-0.38707765) q[2];
sx q[2];
rz(-2.8708411) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5987451) q[1];
sx q[1];
rz(-1.3951021) q[1];
sx q[1];
rz(-0.41136841) q[1];
rz(-pi) q[2];
rz(-1.4275527) q[3];
sx q[3];
rz(-2.4510265) q[3];
sx q[3];
rz(-2.7018202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9294372) q[2];
sx q[2];
rz(-1.1634469) q[2];
sx q[2];
rz(2.4003417) q[2];
rz(2.6397928) q[3];
sx q[3];
rz(-1.3391677) q[3];
sx q[3];
rz(-2.5575976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(1.0989477) q[0];
sx q[0];
rz(-0.89530033) q[0];
sx q[0];
rz(2.6328971) q[0];
rz(0.11518654) q[1];
sx q[1];
rz(-2.4337264) q[1];
sx q[1];
rz(-0.68181109) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54036056) q[0];
sx q[0];
rz(-2.1170179) q[0];
sx q[0];
rz(2.9451314) q[0];
x q[1];
rz(-2.4585312) q[2];
sx q[2];
rz(-0.94711727) q[2];
sx q[2];
rz(-1.1766441) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9224285) q[1];
sx q[1];
rz(-2.104496) q[1];
sx q[1];
rz(-2.1283172) q[1];
rz(-pi) q[2];
rz(-2.1308594) q[3];
sx q[3];
rz(-2.0252953) q[3];
sx q[3];
rz(1.0487995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.29601413) q[2];
sx q[2];
rz(-1.5073551) q[2];
sx q[2];
rz(0.34109035) q[2];
rz(-2.0579445) q[3];
sx q[3];
rz(-0.79569474) q[3];
sx q[3];
rz(2.9705689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9733799) q[0];
sx q[0];
rz(-1.6642878) q[0];
sx q[0];
rz(-0.99933495) q[0];
rz(-2.534261) q[1];
sx q[1];
rz(-0.9098396) q[1];
sx q[1];
rz(0.20914016) q[1];
rz(-0.63275679) q[2];
sx q[2];
rz(-2.3709595) q[2];
sx q[2];
rz(-0.79216935) q[2];
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
