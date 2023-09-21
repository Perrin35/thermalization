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
rz(1.4446422) q[1];
sx q[1];
rz(-1.1029707) q[1];
sx q[1];
rz(-0.77494088) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2431508) q[0];
sx q[0];
rz(-1.713322) q[0];
sx q[0];
rz(-1.6210763) q[0];
rz(-pi) q[1];
rz(-0.84471976) q[2];
sx q[2];
rz(-1.2436927) q[2];
sx q[2];
rz(-2.6127882) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.85841132) q[1];
sx q[1];
rz(-1.6106669) q[1];
sx q[1];
rz(-1.7130501) q[1];
rz(-pi) q[2];
rz(-1.418387) q[3];
sx q[3];
rz(-1.2890352) q[3];
sx q[3];
rz(1.8834653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0960192) q[2];
sx q[2];
rz(-2.0298268) q[2];
sx q[2];
rz(-1.9457031) q[2];
rz(1.9879509) q[3];
sx q[3];
rz(-2.3524645) q[3];
sx q[3];
rz(-1.7529863) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7213223) q[0];
sx q[0];
rz(-2.7852311) q[0];
sx q[0];
rz(-1.8700245) q[0];
rz(-1.0999854) q[1];
sx q[1];
rz(-1.1106691) q[1];
sx q[1];
rz(1.3756479) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028582024) q[0];
sx q[0];
rz(-2.7227289) q[0];
sx q[0];
rz(-2.9257665) q[0];
rz(1.7550811) q[2];
sx q[2];
rz(-2.8993336) q[2];
sx q[2];
rz(0.87528961) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.83982044) q[1];
sx q[1];
rz(-0.41026792) q[1];
sx q[1];
rz(0.46373414) q[1];
x q[2];
rz(0.35956412) q[3];
sx q[3];
rz(-1.0457977) q[3];
sx q[3];
rz(-2.0294702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9282844) q[2];
sx q[2];
rz(-1.3457315) q[2];
sx q[2];
rz(0.48970547) q[2];
rz(-2.0080163) q[3];
sx q[3];
rz(-2.9462892) q[3];
sx q[3];
rz(-2.074923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
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
rz(0.60004822) q[0];
sx q[0];
rz(-1.2788037) q[0];
sx q[0];
rz(0.24060732) q[0];
rz(-2.799017) q[1];
sx q[1];
rz(-0.97476417) q[1];
sx q[1];
rz(1.906357) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89861682) q[0];
sx q[0];
rz(-1.1968062) q[0];
sx q[0];
rz(-0.6181194) q[0];
rz(-pi) q[1];
rz(-2.5580102) q[2];
sx q[2];
rz(-2.7036813) q[2];
sx q[2];
rz(2.8180518) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5879844) q[1];
sx q[1];
rz(-1.3097035) q[1];
sx q[1];
rz(-1.8334465) q[1];
rz(-pi) q[2];
rz(-0.62883212) q[3];
sx q[3];
rz(-1.9753846) q[3];
sx q[3];
rz(1.42266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3774595) q[2];
sx q[2];
rz(-1.5185792) q[2];
sx q[2];
rz(1.7049449) q[2];
rz(-1.7403587) q[3];
sx q[3];
rz(-1.2652206) q[3];
sx q[3];
rz(2.5333372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
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
rz(-0.80379379) q[0];
rz(0.94961387) q[1];
sx q[1];
rz(-1.4664374) q[1];
sx q[1];
rz(3.0217357) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0474931) q[0];
sx q[0];
rz(-2.1764628) q[0];
sx q[0];
rz(0.043461965) q[0];
x q[1];
rz(0.36977936) q[2];
sx q[2];
rz(-2.7571207) q[2];
sx q[2];
rz(1.0880926) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2377492) q[1];
sx q[1];
rz(-1.3386054) q[1];
sx q[1];
rz(3.1234427) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4265392) q[3];
sx q[3];
rz(-2.4063769) q[3];
sx q[3];
rz(-0.27763593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5084761) q[2];
sx q[2];
rz(-1.8957596) q[2];
sx q[2];
rz(0.072337739) q[2];
rz(-2.7667601) q[3];
sx q[3];
rz(-2.4852677) q[3];
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
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6937834) q[0];
sx q[0];
rz(-0.57149514) q[0];
sx q[0];
rz(-2.6547292) q[0];
rz(-0.72987366) q[1];
sx q[1];
rz(-0.90881538) q[1];
sx q[1];
rz(-1.9015076) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2729028) q[0];
sx q[0];
rz(-1.5575952) q[0];
sx q[0];
rz(-0.23916434) q[0];
rz(0.83112049) q[2];
sx q[2];
rz(-2.7925425) q[2];
sx q[2];
rz(-0.30358735) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.33927321) q[1];
sx q[1];
rz(-0.60802751) q[1];
sx q[1];
rz(0.1165216) q[1];
x q[2];
rz(0.48360444) q[3];
sx q[3];
rz(-1.0217474) q[3];
sx q[3];
rz(-2.1476114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1029677) q[2];
sx q[2];
rz(-2.2323148) q[2];
sx q[2];
rz(0.70303482) q[2];
rz(-1.4098343) q[3];
sx q[3];
rz(-1.3491646) q[3];
sx q[3];
rz(0.15771244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96419656) q[0];
sx q[0];
rz(-1.8494158) q[0];
sx q[0];
rz(0.99037209) q[0];
rz(-3.0888427) q[1];
sx q[1];
rz(-2.2149448) q[1];
sx q[1];
rz(1.6606768) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2213343) q[0];
sx q[0];
rz(-1.9517559) q[0];
sx q[0];
rz(-0.078698054) q[0];
x q[1];
rz(3.1408177) q[2];
sx q[2];
rz(-2.0159855) q[2];
sx q[2];
rz(1.0926525) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3546238) q[1];
sx q[1];
rz(-1.0458535) q[1];
sx q[1];
rz(-0.24137361) q[1];
rz(-0.89208608) q[3];
sx q[3];
rz(-1.5329554) q[3];
sx q[3];
rz(2.2298262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.0084373077) q[2];
sx q[2];
rz(-1.6494273) q[2];
sx q[2];
rz(-0.56224242) q[2];
rz(1.0605313) q[3];
sx q[3];
rz(-0.74354592) q[3];
sx q[3];
rz(-2.8806768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19787191) q[0];
sx q[0];
rz(-1.3523538) q[0];
sx q[0];
rz(1.6725756) q[0];
rz(2.127227) q[1];
sx q[1];
rz(-1.0275774) q[1];
sx q[1];
rz(-1.8168824) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6727407) q[0];
sx q[0];
rz(-1.4434837) q[0];
sx q[0];
rz(1.5479814) q[0];
rz(-pi) q[1];
rz(-2.1340738) q[2];
sx q[2];
rz(-2.0118606) q[2];
sx q[2];
rz(-1.2457459) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8285671) q[1];
sx q[1];
rz(-1.9777858) q[1];
sx q[1];
rz(1.9810956) q[1];
rz(-pi) q[2];
rz(0.91388254) q[3];
sx q[3];
rz(-1.2607288) q[3];
sx q[3];
rz(-0.1012181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5614732) q[2];
sx q[2];
rz(-1.7529528) q[2];
sx q[2];
rz(-0.44357792) q[2];
rz(0.94868547) q[3];
sx q[3];
rz(-1.1921927) q[3];
sx q[3];
rz(-0.77478066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.401944) q[0];
sx q[0];
rz(-2.1304603) q[0];
sx q[0];
rz(-0.46052128) q[0];
rz(-3.0415688) q[1];
sx q[1];
rz(-2.1477551) q[1];
sx q[1];
rz(1.8519648) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.078482) q[0];
sx q[0];
rz(-1.0149628) q[0];
sx q[0];
rz(1.0312992) q[0];
x q[1];
rz(0.81823924) q[2];
sx q[2];
rz(-1.0050251) q[2];
sx q[2];
rz(-0.73275369) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3813078) q[1];
sx q[1];
rz(-0.21796255) q[1];
sx q[1];
rz(-1.0337619) q[1];
rz(-pi) q[2];
rz(2.4786948) q[3];
sx q[3];
rz(-1.7015966) q[3];
sx q[3];
rz(-0.42636426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.148968) q[2];
sx q[2];
rz(-0.31105369) q[2];
sx q[2];
rz(-2.0588493) q[2];
rz(0.096171245) q[3];
sx q[3];
rz(-1.440719) q[3];
sx q[3];
rz(-1.910803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7643395) q[0];
sx q[0];
rz(-2.9057673) q[0];
sx q[0];
rz(2.8073231) q[0];
rz(1.9175247) q[1];
sx q[1];
rz(-1.5806438) q[1];
sx q[1];
rz(-0.28265488) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7455505) q[0];
sx q[0];
rz(-2.1746785) q[0];
sx q[0];
rz(-1.0066443) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8106868) q[2];
sx q[2];
rz(-1.7754284) q[2];
sx q[2];
rz(2.3757039) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4947195) q[1];
sx q[1];
rz(-0.44534007) q[1];
sx q[1];
rz(2.7237707) q[1];
rz(-pi) q[2];
x q[2];
rz(1.71404) q[3];
sx q[3];
rz(-0.69056615) q[3];
sx q[3];
rz(2.7018202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.21215542) q[2];
sx q[2];
rz(-1.1634469) q[2];
sx q[2];
rz(-0.7412509) q[2];
rz(-2.6397928) q[3];
sx q[3];
rz(-1.3391677) q[3];
sx q[3];
rz(2.5575976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0989477) q[0];
sx q[0];
rz(-0.89530033) q[0];
sx q[0];
rz(-0.50869554) q[0];
rz(3.0264061) q[1];
sx q[1];
rz(-0.70786628) q[1];
sx q[1];
rz(2.4597816) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0081351) q[0];
sx q[0];
rz(-1.7383766) q[0];
sx q[0];
rz(2.125678) q[0];
x q[1];
rz(-2.4585312) q[2];
sx q[2];
rz(-0.94711727) q[2];
sx q[2];
rz(1.9649486) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.21916418) q[1];
sx q[1];
rz(-2.104496) q[1];
sx q[1];
rz(-2.1283172) q[1];
x q[2];
rz(-2.3144249) q[3];
sx q[3];
rz(-2.4359772) q[3];
sx q[3];
rz(-3.0527715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8455785) q[2];
sx q[2];
rz(-1.6342376) q[2];
sx q[2];
rz(-2.8005023) q[2];
rz(1.0836481) q[3];
sx q[3];
rz(-0.79569474) q[3];
sx q[3];
rz(-0.17102374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(-1.1682128) q[0];
sx q[0];
rz(-1.4773049) q[0];
sx q[0];
rz(2.1422577) q[0];
rz(-0.60733168) q[1];
sx q[1];
rz(-2.2317531) q[1];
sx q[1];
rz(-2.9324525) q[1];
rz(2.091999) q[2];
sx q[2];
rz(-2.1672719) q[2];
sx q[2];
rz(1.5530829) q[2];
rz(-2.7568983) q[3];
sx q[3];
rz(-1.3833429) q[3];
sx q[3];
rz(-0.27237567) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
