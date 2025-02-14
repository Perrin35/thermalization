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
rz(-0.24902046) q[0];
sx q[0];
rz(5.3244642) q[0];
sx q[0];
rz(9.1992314) q[0];
rz(2.0685937) q[1];
sx q[1];
rz(-0.45773157) q[1];
sx q[1];
rz(0.20888858) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0627014) q[0];
sx q[0];
rz(-1.3008504) q[0];
sx q[0];
rz(-0.45227082) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7643808) q[2];
sx q[2];
rz(-2.4754731) q[2];
sx q[2];
rz(-1.3782901) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1465197) q[1];
sx q[1];
rz(-2.0834298) q[1];
sx q[1];
rz(1.9216838) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.82868242) q[3];
sx q[3];
rz(-2.6912315) q[3];
sx q[3];
rz(0.27638232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1170342) q[2];
sx q[2];
rz(-0.25091761) q[2];
sx q[2];
rz(0.92570242) q[2];
rz(-0.77445817) q[3];
sx q[3];
rz(-1.2239933) q[3];
sx q[3];
rz(-1.8584049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9722612) q[0];
sx q[0];
rz(-2.4915578) q[0];
sx q[0];
rz(1.6768804) q[0];
rz(0.88893923) q[1];
sx q[1];
rz(-0.89193946) q[1];
sx q[1];
rz(1.3533786) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4735955) q[0];
sx q[0];
rz(-1.2819074) q[0];
sx q[0];
rz(-0.6205361) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.057463138) q[2];
sx q[2];
rz(-2.7749535) q[2];
sx q[2];
rz(2.5861248) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.96129744) q[1];
sx q[1];
rz(-2.5112481) q[1];
sx q[1];
rz(-0.35758361) q[1];
rz(-pi) q[2];
rz(0.7781182) q[3];
sx q[3];
rz(-1.3740842) q[3];
sx q[3];
rz(2.9784378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9951524) q[2];
sx q[2];
rz(-2.1285987) q[2];
sx q[2];
rz(-0.97274441) q[2];
rz(1.815833) q[3];
sx q[3];
rz(-1.1498068) q[3];
sx q[3];
rz(-2.5241234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72508088) q[0];
sx q[0];
rz(-1.4840115) q[0];
sx q[0];
rz(-3.1120279) q[0];
rz(-2.120453) q[1];
sx q[1];
rz(-2.857326) q[1];
sx q[1];
rz(0.31612843) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97702867) q[0];
sx q[0];
rz(-3.1119149) q[0];
sx q[0];
rz(0.82334955) q[0];
rz(1.7720376) q[2];
sx q[2];
rz(-1.4788718) q[2];
sx q[2];
rz(-2.6215009) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.21507922) q[1];
sx q[1];
rz(-1.2383019) q[1];
sx q[1];
rz(2.957587) q[1];
rz(-pi) q[2];
rz(-3.0532794) q[3];
sx q[3];
rz(-2.4188015) q[3];
sx q[3];
rz(0.73107728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.11008392) q[2];
sx q[2];
rz(-1.8434593) q[2];
sx q[2];
rz(0.52424866) q[2];
rz(0.60018572) q[3];
sx q[3];
rz(-2.6781008) q[3];
sx q[3];
rz(2.0704827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6780739) q[0];
sx q[0];
rz(-1.9215895) q[0];
sx q[0];
rz(-2.779261) q[0];
rz(2.433297) q[1];
sx q[1];
rz(-2.7999122) q[1];
sx q[1];
rz(-1.1452311) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1168829) q[0];
sx q[0];
rz(-1.4636466) q[0];
sx q[0];
rz(3.0610132) q[0];
rz(0.87605642) q[2];
sx q[2];
rz(-1.9578906) q[2];
sx q[2];
rz(0.53384483) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.64001361) q[1];
sx q[1];
rz(-0.7487491) q[1];
sx q[1];
rz(-1.9489991) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8484861) q[3];
sx q[3];
rz(-0.17129843) q[3];
sx q[3];
rz(0.78106487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6142673) q[2];
sx q[2];
rz(-0.82794398) q[2];
sx q[2];
rz(0.0086199363) q[2];
rz(-0.65420592) q[3];
sx q[3];
rz(-2.4253186) q[3];
sx q[3];
rz(2.8692828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24097405) q[0];
sx q[0];
rz(-2.8130377) q[0];
sx q[0];
rz(0.7884489) q[0];
rz(-2.12766) q[1];
sx q[1];
rz(-1.8221816) q[1];
sx q[1];
rz(-0.088277146) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7006314) q[0];
sx q[0];
rz(-0.92871237) q[0];
sx q[0];
rz(-2.2034646) q[0];
rz(-pi) q[1];
rz(0.16096756) q[2];
sx q[2];
rz(-2.5720398) q[2];
sx q[2];
rz(-0.19191027) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.83206415) q[1];
sx q[1];
rz(-0.58575478) q[1];
sx q[1];
rz(1.0447211) q[1];
rz(-pi) q[2];
rz(0.027242335) q[3];
sx q[3];
rz(-1.0758531) q[3];
sx q[3];
rz(-2.3875007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2225515) q[2];
sx q[2];
rz(-1.4107979) q[2];
sx q[2];
rz(-0.50773531) q[2];
rz(0.12568411) q[3];
sx q[3];
rz(-1.1832184) q[3];
sx q[3];
rz(2.3465033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3313726) q[0];
sx q[0];
rz(-2.3851244) q[0];
sx q[0];
rz(0.08547011) q[0];
rz(-0.62689176) q[1];
sx q[1];
rz(-1.9848928) q[1];
sx q[1];
rz(1.289182) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3835433) q[0];
sx q[0];
rz(-1.4669167) q[0];
sx q[0];
rz(1.770875) q[0];
rz(-1.7250546) q[2];
sx q[2];
rz(-2.3084894) q[2];
sx q[2];
rz(-1.8869893) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.61633718) q[1];
sx q[1];
rz(-2.8477745) q[1];
sx q[1];
rz(-2.1928685) q[1];
x q[2];
rz(-1.5103666) q[3];
sx q[3];
rz(-2.9039318) q[3];
sx q[3];
rz(0.56277871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.36279303) q[2];
sx q[2];
rz(-0.94393602) q[2];
sx q[2];
rz(1.2550521) q[2];
rz(-3.0826027) q[3];
sx q[3];
rz(-1.2041644) q[3];
sx q[3];
rz(0.98439938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4208218) q[0];
sx q[0];
rz(-1.868792) q[0];
sx q[0];
rz(0.76501784) q[0];
rz(1.4014686) q[1];
sx q[1];
rz(-0.88686371) q[1];
sx q[1];
rz(2.9041362) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9119806) q[0];
sx q[0];
rz(-1.53541) q[0];
sx q[0];
rz(-2.1982212) q[0];
rz(-pi) q[1];
rz(-0.075421926) q[2];
sx q[2];
rz(-0.89176501) q[2];
sx q[2];
rz(2.0184269) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.85126696) q[1];
sx q[1];
rz(-1.8614166) q[1];
sx q[1];
rz(-1.6162916) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6066437) q[3];
sx q[3];
rz(-0.69679835) q[3];
sx q[3];
rz(0.076059503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.456363) q[2];
sx q[2];
rz(-0.80358973) q[2];
sx q[2];
rz(2.2881499) q[2];
rz(1.338909) q[3];
sx q[3];
rz(-1.5693376) q[3];
sx q[3];
rz(0.16551031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6859739) q[0];
sx q[0];
rz(-1.2116665) q[0];
sx q[0];
rz(2.9562505) q[0];
rz(2.7475884) q[1];
sx q[1];
rz(-1.3961671) q[1];
sx q[1];
rz(-2.0448763) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1117437) q[0];
sx q[0];
rz(-1.1858479) q[0];
sx q[0];
rz(-1.0551595) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.11285891) q[2];
sx q[2];
rz(-0.90870171) q[2];
sx q[2];
rz(-2.2619132) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1450069) q[1];
sx q[1];
rz(-1.2357986) q[1];
sx q[1];
rz(-0.62535357) q[1];
rz(-pi) q[2];
x q[2];
rz(0.59456749) q[3];
sx q[3];
rz(-1.2487186) q[3];
sx q[3];
rz(2.8753124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.77157053) q[2];
sx q[2];
rz(-0.084150704) q[2];
sx q[2];
rz(2.8328075) q[2];
rz(1.2540981) q[3];
sx q[3];
rz(-1.8697238) q[3];
sx q[3];
rz(1.1312283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9517188) q[0];
sx q[0];
rz(-1.8872486) q[0];
sx q[0];
rz(-0.21981123) q[0];
rz(-1.9728569) q[1];
sx q[1];
rz(-1.3558148) q[1];
sx q[1];
rz(1.8692325) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.397906) q[0];
sx q[0];
rz(-2.1748161) q[0];
sx q[0];
rz(0.029725909) q[0];
rz(-pi) q[1];
rz(-0.15064871) q[2];
sx q[2];
rz(-2.6259171) q[2];
sx q[2];
rz(0.0055238481) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2337964) q[1];
sx q[1];
rz(-1.4005473) q[1];
sx q[1];
rz(-2.7336804) q[1];
x q[2];
rz(0.50140372) q[3];
sx q[3];
rz(-2.025017) q[3];
sx q[3];
rz(2.8387031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9610338) q[2];
sx q[2];
rz(-2.2982633) q[2];
sx q[2];
rz(0.45832222) q[2];
rz(0.35442963) q[3];
sx q[3];
rz(-1.3684045) q[3];
sx q[3];
rz(1.9663158) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71427041) q[0];
sx q[0];
rz(-1.9388119) q[0];
sx q[0];
rz(0.74139968) q[0];
rz(-1.7085913) q[1];
sx q[1];
rz(-0.42867908) q[1];
sx q[1];
rz(0.0892078) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87338221) q[0];
sx q[0];
rz(-2.6307627) q[0];
sx q[0];
rz(-1.9477773) q[0];
rz(0.54069767) q[2];
sx q[2];
rz(-0.6650228) q[2];
sx q[2];
rz(2.2407209) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9392689) q[1];
sx q[1];
rz(-1.2319733) q[1];
sx q[1];
rz(-0.044087709) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9815638) q[3];
sx q[3];
rz(-0.8770408) q[3];
sx q[3];
rz(-2.957537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0006813) q[2];
sx q[2];
rz(-0.050364308) q[2];
sx q[2];
rz(-2.4939406) q[2];
rz(-2.7378313) q[3];
sx q[3];
rz(-1.8881366) q[3];
sx q[3];
rz(3.0103179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(1.1972926) q[0];
sx q[0];
rz(-0.97351749) q[0];
sx q[0];
rz(1.0607006) q[0];
rz(-2.3464959) q[1];
sx q[1];
rz(-0.92082321) q[1];
sx q[1];
rz(-0.77650741) q[1];
rz(-2.2262103) q[2];
sx q[2];
rz(-1.6554828) q[2];
sx q[2];
rz(0.2632904) q[2];
rz(0.85579167) q[3];
sx q[3];
rz(-3.0223911) q[3];
sx q[3];
rz(0.08798616) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
