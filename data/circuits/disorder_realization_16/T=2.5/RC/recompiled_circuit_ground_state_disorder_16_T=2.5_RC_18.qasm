OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.032441) q[0];
sx q[0];
rz(-1.1530131) q[0];
sx q[0];
rz(3.0670526) q[0];
rz(-1.5300765) q[1];
sx q[1];
rz(-0.059377436) q[1];
sx q[1];
rz(0.52409726) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5016035) q[0];
sx q[0];
rz(-1.934938) q[0];
sx q[0];
rz(0.95567165) q[0];
rz(-2.8048326) q[2];
sx q[2];
rz(-1.3408061) q[2];
sx q[2];
rz(0.57631341) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1519231) q[1];
sx q[1];
rz(-2.1300585) q[1];
sx q[1];
rz(2.2733746) q[1];
x q[2];
rz(-3.1007441) q[3];
sx q[3];
rz(-1.4788879) q[3];
sx q[3];
rz(-0.81113863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3268299) q[2];
sx q[2];
rz(-0.49850285) q[2];
sx q[2];
rz(-0.89547431) q[2];
rz(-2.879066) q[3];
sx q[3];
rz(-0.34590507) q[3];
sx q[3];
rz(-3.053022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.9240016) q[0];
sx q[0];
rz(-1.9460678) q[0];
sx q[0];
rz(3.0254645) q[0];
rz(-2.5092292) q[1];
sx q[1];
rz(-0.26591161) q[1];
sx q[1];
rz(-1.4858474) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6054942) q[0];
sx q[0];
rz(-1.5935531) q[0];
sx q[0];
rz(-0.066734138) q[0];
x q[1];
rz(-0.2367649) q[2];
sx q[2];
rz(-2.1000266) q[2];
sx q[2];
rz(0.53073149) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.33148674) q[1];
sx q[1];
rz(-1.6983733) q[1];
sx q[1];
rz(-1.677631) q[1];
rz(-pi) q[2];
x q[2];
rz(0.13074517) q[3];
sx q[3];
rz(-0.8746038) q[3];
sx q[3];
rz(-1.2519022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.31132013) q[2];
sx q[2];
rz(-2.1612284) q[2];
sx q[2];
rz(-2.2118528) q[2];
rz(0.43944198) q[3];
sx q[3];
rz(-1.5403055) q[3];
sx q[3];
rz(0.75696993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.065652549) q[0];
sx q[0];
rz(-0.75242281) q[0];
sx q[0];
rz(-2.2546076) q[0];
rz(-1.2607964) q[1];
sx q[1];
rz(-1.1075243) q[1];
sx q[1];
rz(-1.4874123) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9943574) q[0];
sx q[0];
rz(-1.5702797) q[0];
sx q[0];
rz(-1.6883255) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9145221) q[2];
sx q[2];
rz(-2.7393574) q[2];
sx q[2];
rz(-2.9410597) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9650813) q[1];
sx q[1];
rz(-1.9002894) q[1];
sx q[1];
rz(1.552187) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2441176) q[3];
sx q[3];
rz(-1.0138113) q[3];
sx q[3];
rz(-0.1989593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4999353) q[2];
sx q[2];
rz(-2.1782918) q[2];
sx q[2];
rz(-0.43528834) q[2];
rz(-2.7564202) q[3];
sx q[3];
rz(-1.9305072) q[3];
sx q[3];
rz(2.1729573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4061072) q[0];
sx q[0];
rz(-3.0574953) q[0];
sx q[0];
rz(-0.054585833) q[0];
rz(1.5054043) q[1];
sx q[1];
rz(-1.2137698) q[1];
sx q[1];
rz(-2.7240567) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62359257) q[0];
sx q[0];
rz(-1.5125649) q[0];
sx q[0];
rz(-0.43727711) q[0];
x q[1];
rz(-0.63234858) q[2];
sx q[2];
rz(-1.8397619) q[2];
sx q[2];
rz(-0.40695813) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3247973) q[1];
sx q[1];
rz(-2.9831605) q[1];
sx q[1];
rz(-3.0976803) q[1];
x q[2];
rz(0.066922234) q[3];
sx q[3];
rz(-1.8324781) q[3];
sx q[3];
rz(2.8953538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2893082) q[2];
sx q[2];
rz(-2.4384629) q[2];
sx q[2];
rz(-0.0011778041) q[2];
rz(-1.3834472) q[3];
sx q[3];
rz(-0.26418424) q[3];
sx q[3];
rz(-2.0980825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1498639) q[0];
sx q[0];
rz(-2.9809451) q[0];
sx q[0];
rz(-2.7657261) q[0];
rz(2.1916892) q[1];
sx q[1];
rz(-1.4774731) q[1];
sx q[1];
rz(-2.4163767) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26869088) q[0];
sx q[0];
rz(-2.9753472) q[0];
sx q[0];
rz(-0.77456559) q[0];
rz(-pi) q[1];
rz(0.8704758) q[2];
sx q[2];
rz(-2.0904125) q[2];
sx q[2];
rz(-0.82110559) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8661853) q[1];
sx q[1];
rz(-1.5650041) q[1];
sx q[1];
rz(-0.097472982) q[1];
x q[2];
rz(2.83738) q[3];
sx q[3];
rz(-0.72644573) q[3];
sx q[3];
rz(0.026487984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3788562) q[2];
sx q[2];
rz(-1.4521658) q[2];
sx q[2];
rz(-2.6689996) q[2];
rz(0.11911123) q[3];
sx q[3];
rz(-0.62370682) q[3];
sx q[3];
rz(2.3314893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8238207) q[0];
sx q[0];
rz(-2.9404984) q[0];
sx q[0];
rz(0.066545181) q[0];
rz(-2.9464974) q[1];
sx q[1];
rz(-1.0761484) q[1];
sx q[1];
rz(1.9871575) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6267363) q[0];
sx q[0];
rz(-1.8033228) q[0];
sx q[0];
rz(-1.0794845) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4835266) q[2];
sx q[2];
rz(-2.0258528) q[2];
sx q[2];
rz(-2.1966962) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.30224702) q[1];
sx q[1];
rz(-1.463942) q[1];
sx q[1];
rz(1.2059421) q[1];
x q[2];
rz(-2.1514417) q[3];
sx q[3];
rz(-1.5599453) q[3];
sx q[3];
rz(-0.059338245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.92516148) q[2];
sx q[2];
rz(-1.5387646) q[2];
sx q[2];
rz(-2.4161941) q[2];
rz(-1.7321436) q[3];
sx q[3];
rz(-1.9603445) q[3];
sx q[3];
rz(-0.60666549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.055939097) q[0];
sx q[0];
rz(-0.5760718) q[0];
sx q[0];
rz(2.2357909) q[0];
rz(0.55465758) q[1];
sx q[1];
rz(-2.3034818) q[1];
sx q[1];
rz(2.99756) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5921123) q[0];
sx q[0];
rz(-2.0680111) q[0];
sx q[0];
rz(-0.29659941) q[0];
rz(-pi) q[1];
rz(-1.7212825) q[2];
sx q[2];
rz(-2.624199) q[2];
sx q[2];
rz(-0.69240332) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1293948) q[1];
sx q[1];
rz(-1.8499814) q[1];
sx q[1];
rz(1.1711981) q[1];
rz(-2.9908189) q[3];
sx q[3];
rz(-1.322804) q[3];
sx q[3];
rz(-1.4797321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4732699) q[2];
sx q[2];
rz(-2.747135) q[2];
sx q[2];
rz(-2.1486166) q[2];
rz(0.18443491) q[3];
sx q[3];
rz(-2.1858229) q[3];
sx q[3];
rz(0.9930281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1297146) q[0];
sx q[0];
rz(-0.60420245) q[0];
sx q[0];
rz(-2.7469444) q[0];
rz(2.6036085) q[1];
sx q[1];
rz(-2.4514276) q[1];
sx q[1];
rz(-1.1916377) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79121548) q[0];
sx q[0];
rz(-1.0320623) q[0];
sx q[0];
rz(1.8693493) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7685652) q[2];
sx q[2];
rz(-1.3907888) q[2];
sx q[2];
rz(1.4733835) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0816457) q[1];
sx q[1];
rz(-2.4103202) q[1];
sx q[1];
rz(2.251723) q[1];
rz(-pi) q[2];
rz(-2.0268782) q[3];
sx q[3];
rz(-2.1563765) q[3];
sx q[3];
rz(-2.8600313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3523606) q[2];
sx q[2];
rz(-1.5287986) q[2];
sx q[2];
rz(-2.1996876) q[2];
rz(2.3801129) q[3];
sx q[3];
rz(-2.0165636) q[3];
sx q[3];
rz(3.1010845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2725459) q[0];
sx q[0];
rz(-1.7750374) q[0];
sx q[0];
rz(1.0719365) q[0];
rz(2.5275224) q[1];
sx q[1];
rz(-2.2449988) q[1];
sx q[1];
rz(0.44642064) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0508795) q[0];
sx q[0];
rz(-1.2901559) q[0];
sx q[0];
rz(-2.2801823) q[0];
rz(-0.69009366) q[2];
sx q[2];
rz(-1.9999256) q[2];
sx q[2];
rz(-2.5370425) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3371967) q[1];
sx q[1];
rz(-2.2257518) q[1];
sx q[1];
rz(-1.4319929) q[1];
rz(-pi) q[2];
x q[2];
rz(0.29576755) q[3];
sx q[3];
rz(-2.6172574) q[3];
sx q[3];
rz(-1.7400896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8322231) q[2];
sx q[2];
rz(-2.1840405) q[2];
sx q[2];
rz(-1.2584125) q[2];
rz(1.9084515) q[3];
sx q[3];
rz(-0.62658739) q[3];
sx q[3];
rz(1.8241749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.719139) q[0];
sx q[0];
rz(-0.8154251) q[0];
sx q[0];
rz(1.423214) q[0];
rz(-2.8990959) q[1];
sx q[1];
rz(-0.47912326) q[1];
sx q[1];
rz(-1.8877782) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72287382) q[0];
sx q[0];
rz(-2.6986045) q[0];
sx q[0];
rz(1.1305869) q[0];
rz(-1.5780588) q[2];
sx q[2];
rz(-1.6374131) q[2];
sx q[2];
rz(0.74502258) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1088646) q[1];
sx q[1];
rz(-2.0833587) q[1];
sx q[1];
rz(-0.3442007) q[1];
rz(-pi) q[2];
rz(-1.0494136) q[3];
sx q[3];
rz(-0.42179042) q[3];
sx q[3];
rz(0.78027314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8157876) q[2];
sx q[2];
rz(-0.82745224) q[2];
sx q[2];
rz(-0.023836689) q[2];
rz(1.2787974) q[3];
sx q[3];
rz(-0.33134225) q[3];
sx q[3];
rz(2.3648025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9780289) q[0];
sx q[0];
rz(-1.4373056) q[0];
sx q[0];
rz(1.9594255) q[0];
rz(0.72721807) q[1];
sx q[1];
rz(-1.4677508) q[1];
sx q[1];
rz(-0.8263091) q[1];
rz(-3.1221409) q[2];
sx q[2];
rz(-1.7218334) q[2];
sx q[2];
rz(-1.3109372) q[2];
rz(-1.6926077) q[3];
sx q[3];
rz(-1.5234608) q[3];
sx q[3];
rz(0.85352637) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
