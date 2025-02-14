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
rz(-0.95872107) q[0];
sx q[0];
rz(-0.22554654) q[0];
rz(-1.072999) q[1];
sx q[1];
rz(3.5993242) q[1];
sx q[1];
rz(9.2158894) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0627014) q[0];
sx q[0];
rz(-1.3008504) q[0];
sx q[0];
rz(0.45227082) q[0];
rz(-1.8525847) q[2];
sx q[2];
rz(-2.1827841) q[2];
sx q[2];
rz(2.2300697) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.99507299) q[1];
sx q[1];
rz(-2.0834298) q[1];
sx q[1];
rz(-1.2199089) q[1];
x q[2];
rz(0.82868242) q[3];
sx q[3];
rz(-2.6912315) q[3];
sx q[3];
rz(2.8652103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1170342) q[2];
sx q[2];
rz(-2.890675) q[2];
sx q[2];
rz(0.92570242) q[2];
rz(-2.3671345) q[3];
sx q[3];
rz(-1.2239933) q[3];
sx q[3];
rz(1.8584049) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9722612) q[0];
sx q[0];
rz(-0.65003482) q[0];
sx q[0];
rz(-1.6768804) q[0];
rz(-0.88893923) q[1];
sx q[1];
rz(-2.2496532) q[1];
sx q[1];
rz(-1.7882141) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.618093) q[0];
sx q[0];
rz(-0.67636469) q[0];
sx q[0];
rz(-2.6690865) q[0];
rz(-pi) q[1];
rz(-3.0841295) q[2];
sx q[2];
rz(-0.36663917) q[2];
sx q[2];
rz(2.5861248) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2389489) q[1];
sx q[1];
rz(-1.7785935) q[1];
sx q[1];
rz(0.59955876) q[1];
rz(-0.27663265) q[3];
sx q[3];
rz(-2.3440954) q[3];
sx q[3];
rz(-1.5381787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9951524) q[2];
sx q[2];
rz(-2.1285987) q[2];
sx q[2];
rz(-0.97274441) q[2];
rz(1.815833) q[3];
sx q[3];
rz(-1.9917859) q[3];
sx q[3];
rz(2.5241234) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72508088) q[0];
sx q[0];
rz(-1.6575811) q[0];
sx q[0];
rz(3.1120279) q[0];
rz(1.0211396) q[1];
sx q[1];
rz(-0.28426668) q[1];
sx q[1];
rz(2.8254642) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8005978) q[0];
sx q[0];
rz(-1.5506239) q[0];
sx q[0];
rz(1.5925657) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7720376) q[2];
sx q[2];
rz(-1.4788718) q[2];
sx q[2];
rz(-0.52009174) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7252032) q[1];
sx q[1];
rz(-1.3969743) q[1];
sx q[1];
rz(1.2330139) q[1];
rz(-pi) q[2];
rz(0.088313266) q[3];
sx q[3];
rz(-2.4188015) q[3];
sx q[3];
rz(0.73107728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.11008392) q[2];
sx q[2];
rz(-1.2981334) q[2];
sx q[2];
rz(0.52424866) q[2];
rz(-2.5414069) q[3];
sx q[3];
rz(-0.46349183) q[3];
sx q[3];
rz(-2.0704827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6780739) q[0];
sx q[0];
rz(-1.9215895) q[0];
sx q[0];
rz(0.36233166) q[0];
rz(-0.70829567) q[1];
sx q[1];
rz(-0.34168044) q[1];
sx q[1];
rz(-1.9963616) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53745054) q[0];
sx q[0];
rz(-1.49068) q[0];
sx q[0];
rz(-1.6782922) q[0];
rz(-2.2655362) q[2];
sx q[2];
rz(-1.9578906) q[2];
sx q[2];
rz(-2.6077478) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.64001361) q[1];
sx q[1];
rz(-2.3928436) q[1];
sx q[1];
rz(-1.9489991) q[1];
rz(-pi) q[2];
rz(2.8484861) q[3];
sx q[3];
rz(-0.17129843) q[3];
sx q[3];
rz(-2.3605278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.5273253) q[2];
sx q[2];
rz(-2.3136487) q[2];
sx q[2];
rz(0.0086199363) q[2];
rz(2.4873867) q[3];
sx q[3];
rz(-0.71627408) q[3];
sx q[3];
rz(0.27230984) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9006186) q[0];
sx q[0];
rz(-2.8130377) q[0];
sx q[0];
rz(0.7884489) q[0];
rz(2.12766) q[1];
sx q[1];
rz(-1.8221816) q[1];
sx q[1];
rz(0.088277146) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4409613) q[0];
sx q[0];
rz(-0.92871237) q[0];
sx q[0];
rz(2.2034646) q[0];
x q[1];
rz(-0.16096756) q[2];
sx q[2];
rz(-2.5720398) q[2];
sx q[2];
rz(-2.9496824) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8535117) q[1];
sx q[1];
rz(-1.8520903) q[1];
sx q[1];
rz(1.049925) q[1];
rz(-pi) q[2];
rz(3.1143503) q[3];
sx q[3];
rz(-2.0657396) q[3];
sx q[3];
rz(0.75409192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.91904116) q[2];
sx q[2];
rz(-1.4107979) q[2];
sx q[2];
rz(-2.6338573) q[2];
rz(-3.0159085) q[3];
sx q[3];
rz(-1.9583743) q[3];
sx q[3];
rz(-2.3465033) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3313726) q[0];
sx q[0];
rz(-2.3851244) q[0];
sx q[0];
rz(3.0561225) q[0];
rz(0.62689176) q[1];
sx q[1];
rz(-1.1566999) q[1];
sx q[1];
rz(-1.8524106) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7580494) q[0];
sx q[0];
rz(-1.4669167) q[0];
sx q[0];
rz(-1.770875) q[0];
rz(-1.7250546) q[2];
sx q[2];
rz(-2.3084894) q[2];
sx q[2];
rz(1.2546033) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.35297817) q[1];
sx q[1];
rz(-1.4012238) q[1];
sx q[1];
rz(-1.3296849) q[1];
rz(-pi) q[2];
rz(1.6312261) q[3];
sx q[3];
rz(-2.9039318) q[3];
sx q[3];
rz(-2.5788139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7787996) q[2];
sx q[2];
rz(-0.94393602) q[2];
sx q[2];
rz(-1.2550521) q[2];
rz(-3.0826027) q[3];
sx q[3];
rz(-1.9374282) q[3];
sx q[3];
rz(2.1571933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7207709) q[0];
sx q[0];
rz(-1.868792) q[0];
sx q[0];
rz(2.3765748) q[0];
rz(1.7401241) q[1];
sx q[1];
rz(-0.88686371) q[1];
sx q[1];
rz(0.23745647) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7516605) q[0];
sx q[0];
rz(-0.62828763) q[0];
sx q[0];
rz(1.5105672) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2512194) q[2];
sx q[2];
rz(-1.6294663) q[2];
sx q[2];
rz(-0.49505297) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.73257414) q[1];
sx q[1];
rz(-1.5272101) q[1];
sx q[1];
rz(0.2909046) q[1];
x q[2];
rz(-0.87431425) q[3];
sx q[3];
rz(-1.5937991) q[3];
sx q[3];
rz(1.5222331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.456363) q[2];
sx q[2];
rz(-0.80358973) q[2];
sx q[2];
rz(-0.85344273) q[2];
rz(-1.338909) q[3];
sx q[3];
rz(-1.572255) q[3];
sx q[3];
rz(0.16551031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.6859739) q[0];
sx q[0];
rz(-1.2116665) q[0];
sx q[0];
rz(2.9562505) q[0];
rz(-2.7475884) q[1];
sx q[1];
rz(-1.7454255) q[1];
sx q[1];
rz(-2.0448763) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1857796) q[0];
sx q[0];
rz(-2.5086864) q[0];
sx q[0];
rz(-2.2586185) q[0];
rz(-pi) q[1];
x q[1];
rz(0.11285891) q[2];
sx q[2];
rz(-2.2328909) q[2];
sx q[2];
rz(-2.2619132) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2885159) q[1];
sx q[1];
rz(-0.6986874) q[1];
sx q[1];
rz(-0.53650155) q[1];
rz(-pi) q[2];
rz(2.6043209) q[3];
sx q[3];
rz(-0.66681615) q[3];
sx q[3];
rz(-1.7423354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3700221) q[2];
sx q[2];
rz(-3.0574419) q[2];
sx q[2];
rz(0.30878511) q[2];
rz(-1.2540981) q[3];
sx q[3];
rz(-1.2718688) q[3];
sx q[3];
rz(-2.0103644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18987385) q[0];
sx q[0];
rz(-1.254344) q[0];
sx q[0];
rz(-2.9217814) q[0];
rz(-1.1687357) q[1];
sx q[1];
rz(-1.7857779) q[1];
sx q[1];
rz(1.8692325) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74368661) q[0];
sx q[0];
rz(-0.96677654) q[0];
sx q[0];
rz(-3.1118667) q[0];
rz(-0.51080334) q[2];
sx q[2];
rz(-1.6448717) q[2];
sx q[2];
rz(-1.6965716) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2337964) q[1];
sx q[1];
rz(-1.7410454) q[1];
sx q[1];
rz(2.7336804) q[1];
rz(-pi) q[2];
x q[2];
rz(0.79325647) q[3];
sx q[3];
rz(-2.4783387) q[3];
sx q[3];
rz(-2.548717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1805588) q[2];
sx q[2];
rz(-0.84332931) q[2];
sx q[2];
rz(2.6832704) q[2];
rz(0.35442963) q[3];
sx q[3];
rz(-1.7731881) q[3];
sx q[3];
rz(1.1752769) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71427041) q[0];
sx q[0];
rz(-1.2027807) q[0];
sx q[0];
rz(-2.400193) q[0];
rz(1.4330014) q[1];
sx q[1];
rz(-0.42867908) q[1];
sx q[1];
rz(0.0892078) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8421546) q[0];
sx q[0];
rz(-2.0426867) q[0];
sx q[0];
rz(0.20345511) q[0];
rz(-pi) q[1];
x q[1];
rz(0.59190665) q[2];
sx q[2];
rz(-1.8940261) q[2];
sx q[2];
rz(2.9130251) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7877823) q[1];
sx q[1];
rz(-1.5292166) q[1];
sx q[1];
rz(1.9099243) q[1];
x q[2];
rz(0.73674006) q[3];
sx q[3];
rz(-1.2587446) q[3];
sx q[3];
rz(-1.4832352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1409113) q[2];
sx q[2];
rz(-0.050364308) q[2];
sx q[2];
rz(0.64765206) q[2];
rz(2.7378313) q[3];
sx q[3];
rz(-1.253456) q[3];
sx q[3];
rz(3.0103179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1972926) q[0];
sx q[0];
rz(-2.1680752) q[0];
sx q[0];
rz(-2.0808921) q[0];
rz(2.3464959) q[1];
sx q[1];
rz(-2.2207694) q[1];
sx q[1];
rz(2.3650852) q[1];
rz(0.10666974) q[2];
sx q[2];
rz(-2.223458) q[2];
sx q[2];
rz(-1.3724422) q[2];
rz(-1.480605) q[3];
sx q[3];
rz(-1.492751) q[3];
sx q[3];
rz(0.94730151) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
