OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5821563) q[0];
sx q[0];
rz(-0.59098935) q[0];
sx q[0];
rz(-2.5581869) q[0];
rz(-0.18435873) q[1];
sx q[1];
rz(-2.1579722) q[1];
sx q[1];
rz(0.89259994) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32514206) q[0];
sx q[0];
rz(-2.4976375) q[0];
sx q[0];
rz(-2.08026) q[0];
rz(-pi) q[1];
x q[1];
rz(0.61172723) q[2];
sx q[2];
rz(-2.3731542) q[2];
sx q[2];
rz(-2.7035463) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6751911) q[1];
sx q[1];
rz(-0.89414222) q[1];
sx q[1];
rz(-2.7387709) q[1];
rz(0.75099545) q[3];
sx q[3];
rz(-1.6023811) q[3];
sx q[3];
rz(-0.75814523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2261752) q[2];
sx q[2];
rz(-1.6025275) q[2];
sx q[2];
rz(1.9809451) q[2];
rz(-2.9246269) q[3];
sx q[3];
rz(-2.6187077) q[3];
sx q[3];
rz(-1.0552361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0579257) q[0];
sx q[0];
rz(-2.1766429) q[0];
sx q[0];
rz(-2.5657186) q[0];
rz(1.2469762) q[1];
sx q[1];
rz(-1.2966825) q[1];
sx q[1];
rz(-1.974568) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.089433) q[0];
sx q[0];
rz(-2.8773327) q[0];
sx q[0];
rz(1.7961851) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0300693) q[2];
sx q[2];
rz(-1.7833372) q[2];
sx q[2];
rz(-0.18268798) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.45707073) q[1];
sx q[1];
rz(-1.6324537) q[1];
sx q[1];
rz(-0.80656273) q[1];
rz(-pi) q[2];
x q[2];
rz(0.59252177) q[3];
sx q[3];
rz(-2.3791109) q[3];
sx q[3];
rz(0.034686397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1077659) q[2];
sx q[2];
rz(-1.1788538) q[2];
sx q[2];
rz(2.9272184) q[2];
rz(0.073444627) q[3];
sx q[3];
rz(-0.44973222) q[3];
sx q[3];
rz(-2.8607821) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6630994) q[0];
sx q[0];
rz(-0.61613023) q[0];
sx q[0];
rz(-1.2269155) q[0];
rz(-2.7413209) q[1];
sx q[1];
rz(-1.8881256) q[1];
sx q[1];
rz(2.1267557) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0815711) q[0];
sx q[0];
rz(-2.2131753) q[0];
sx q[0];
rz(-1.6409671) q[0];
rz(0.24200183) q[2];
sx q[2];
rz(-1.9654462) q[2];
sx q[2];
rz(1.637527) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.240757) q[1];
sx q[1];
rz(-1.490331) q[1];
sx q[1];
rz(-0.97297538) q[1];
x q[2];
rz(2.879911) q[3];
sx q[3];
rz(-0.90173429) q[3];
sx q[3];
rz(2.3467968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0456475) q[2];
sx q[2];
rz(-1.0568876) q[2];
sx q[2];
rz(-0.8992368) q[2];
rz(2.4441161) q[3];
sx q[3];
rz(-1.2858425) q[3];
sx q[3];
rz(2.5337059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65524453) q[0];
sx q[0];
rz(-1.0826033) q[0];
sx q[0];
rz(-2.7096601) q[0];
rz(-0.63255429) q[1];
sx q[1];
rz(-2.7245941) q[1];
sx q[1];
rz(-2.5057709) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0974554) q[0];
sx q[0];
rz(-0.53253981) q[0];
sx q[0];
rz(-1.9299279) q[0];
rz(1.3859149) q[2];
sx q[2];
rz(-2.4061678) q[2];
sx q[2];
rz(2.3787969) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7446049) q[1];
sx q[1];
rz(-2.3481391) q[1];
sx q[1];
rz(-2.8810487) q[1];
x q[2];
rz(-2.5084247) q[3];
sx q[3];
rz(-2.79106) q[3];
sx q[3];
rz(1.7115953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9277966) q[2];
sx q[2];
rz(-2.9979604) q[2];
sx q[2];
rz(0.68871838) q[2];
rz(2.8074746) q[3];
sx q[3];
rz(-1.9709316) q[3];
sx q[3];
rz(-2.9978602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14389811) q[0];
sx q[0];
rz(-0.69902885) q[0];
sx q[0];
rz(-2.0671663) q[0];
rz(0.74514666) q[1];
sx q[1];
rz(-1.4829758) q[1];
sx q[1];
rz(2.863046) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6644088) q[0];
sx q[0];
rz(-1.1468977) q[0];
sx q[0];
rz(-0.82188481) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9859221) q[2];
sx q[2];
rz(-0.61741932) q[2];
sx q[2];
rz(1.9215259) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7601732) q[1];
sx q[1];
rz(-0.96254327) q[1];
sx q[1];
rz(-0.45254032) q[1];
rz(-pi) q[2];
rz(-0.99028011) q[3];
sx q[3];
rz(-0.84458447) q[3];
sx q[3];
rz(-0.57675225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4429861) q[2];
sx q[2];
rz(-1.7437982) q[2];
sx q[2];
rz(-2.8473575) q[2];
rz(0.081929835) q[3];
sx q[3];
rz(-0.51920813) q[3];
sx q[3];
rz(0.036227139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0329523) q[0];
sx q[0];
rz(-0.8529129) q[0];
sx q[0];
rz(-0.0090573514) q[0];
rz(2.5065705) q[1];
sx q[1];
rz(-0.68990866) q[1];
sx q[1];
rz(-0.10805282) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2436284) q[0];
sx q[0];
rz(-2.0192696) q[0];
sx q[0];
rz(-1.9237706) q[0];
rz(-pi) q[1];
rz(-2.8073505) q[2];
sx q[2];
rz(-0.40194449) q[2];
sx q[2];
rz(2.3964756) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6335771) q[1];
sx q[1];
rz(-2.0532236) q[1];
sx q[1];
rz(-2.0111994) q[1];
rz(-pi) q[2];
rz(-2.8764358) q[3];
sx q[3];
rz(-1.0401871) q[3];
sx q[3];
rz(-1.6924072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.99284995) q[2];
sx q[2];
rz(-1.381258) q[2];
sx q[2];
rz(-1.2711058) q[2];
rz(-3.0631915) q[3];
sx q[3];
rz(-1.4784808) q[3];
sx q[3];
rz(-1.1318077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(2.8909797) q[0];
sx q[0];
rz(-1.5761292) q[0];
sx q[0];
rz(2.4839731) q[0];
rz(-1.744386) q[1];
sx q[1];
rz(-1.9669292) q[1];
sx q[1];
rz(-2.2479642) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4511787) q[0];
sx q[0];
rz(-2.4126841) q[0];
sx q[0];
rz(2.6364987) q[0];
rz(-pi) q[1];
rz(2.1003175) q[2];
sx q[2];
rz(-2.6011701) q[2];
sx q[2];
rz(-1.668001) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6210729) q[1];
sx q[1];
rz(-0.90172651) q[1];
sx q[1];
rz(1.3135765) q[1];
rz(-pi) q[2];
rz(-1.6679538) q[3];
sx q[3];
rz(-0.81340862) q[3];
sx q[3];
rz(-2.017445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.39067337) q[2];
sx q[2];
rz(-0.94736391) q[2];
sx q[2];
rz(2.612109) q[2];
rz(-0.47618619) q[3];
sx q[3];
rz(-1.3323077) q[3];
sx q[3];
rz(2.7990394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3787518) q[0];
sx q[0];
rz(-1.6289926) q[0];
sx q[0];
rz(2.0513127) q[0];
rz(3.0293363) q[1];
sx q[1];
rz(-1.1039762) q[1];
sx q[1];
rz(1.1539248) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2147804) q[0];
sx q[0];
rz(-1.5367537) q[0];
sx q[0];
rz(2.4916617) q[0];
x q[1];
rz(2.9032193) q[2];
sx q[2];
rz(-2.5775238) q[2];
sx q[2];
rz(0.43018815) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9824144) q[1];
sx q[1];
rz(-0.7632066) q[1];
sx q[1];
rz(-3.0522507) q[1];
x q[2];
rz(2.047394) q[3];
sx q[3];
rz(-0.72924858) q[3];
sx q[3];
rz(0.13343982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5899137) q[2];
sx q[2];
rz(-2.7605197) q[2];
sx q[2];
rz(-2.7611458) q[2];
rz(-2.0137265) q[3];
sx q[3];
rz(-1.3600072) q[3];
sx q[3];
rz(1.5765223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.8001051) q[0];
sx q[0];
rz(-1.8999758) q[0];
sx q[0];
rz(-0.63968101) q[0];
rz(-1.2387964) q[1];
sx q[1];
rz(-1.4150554) q[1];
sx q[1];
rz(-1.9715086) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.026924883) q[0];
sx q[0];
rz(-1.2171193) q[0];
sx q[0];
rz(-0.35004079) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3001928) q[2];
sx q[2];
rz(-2.3286208) q[2];
sx q[2];
rz(-1.9780618) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3635892) q[1];
sx q[1];
rz(-0.52007857) q[1];
sx q[1];
rz(2.7705454) q[1];
rz(-pi) q[2];
rz(-3.0036003) q[3];
sx q[3];
rz(-0.90413168) q[3];
sx q[3];
rz(0.43825144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3141979) q[2];
sx q[2];
rz(-0.65076995) q[2];
sx q[2];
rz(-2.7588552) q[2];
rz(-2.2132204) q[3];
sx q[3];
rz(-1.9675156) q[3];
sx q[3];
rz(2.476957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2255573) q[0];
sx q[0];
rz(-1.6397497) q[0];
sx q[0];
rz(2.8826707) q[0];
rz(2.4312773) q[1];
sx q[1];
rz(-2.0878891) q[1];
sx q[1];
rz(-0.47992596) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5878764) q[0];
sx q[0];
rz(-1.9037316) q[0];
sx q[0];
rz(0.90162189) q[0];
rz(-pi) q[1];
rz(-1.8690228) q[2];
sx q[2];
rz(-2.4546461) q[2];
sx q[2];
rz(2.5145614) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6595315) q[1];
sx q[1];
rz(-1.2188984) q[1];
sx q[1];
rz(-2.5955276) q[1];
x q[2];
rz(2.6142526) q[3];
sx q[3];
rz(-0.96685997) q[3];
sx q[3];
rz(0.37249836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.84247983) q[2];
sx q[2];
rz(-0.92876902) q[2];
sx q[2];
rz(-1.8245565) q[2];
rz(1.2420098) q[3];
sx q[3];
rz(-0.95359355) q[3];
sx q[3];
rz(2.4035113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15923545) q[0];
sx q[0];
rz(-3.1095105) q[0];
sx q[0];
rz(1.6788917) q[0];
rz(-2.1622529) q[1];
sx q[1];
rz(-1.0995438) q[1];
sx q[1];
rz(-0.88811036) q[1];
rz(0.48158823) q[2];
sx q[2];
rz(-1.4143741) q[2];
sx q[2];
rz(1.0798567) q[2];
rz(0.76673037) q[3];
sx q[3];
rz(-2.7646716) q[3];
sx q[3];
rz(-0.9203831) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
