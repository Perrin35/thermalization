OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.4322296) q[0];
sx q[0];
rz(-0.95786434) q[0];
sx q[0];
rz(0.14444484) q[0];
rz(0.56675178) q[1];
sx q[1];
rz(-0.52539879) q[1];
sx q[1];
rz(-2.1638343) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6173965) q[0];
sx q[0];
rz(-0.27271909) q[0];
sx q[0];
rz(2.8196536) q[0];
rz(-pi) q[1];
rz(0.30774967) q[2];
sx q[2];
rz(-1.0767184) q[2];
sx q[2];
rz(0.57927629) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2797151) q[1];
sx q[1];
rz(-1.4289083) q[1];
sx q[1];
rz(3.0008297) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7498738) q[3];
sx q[3];
rz(-2.350051) q[3];
sx q[3];
rz(0.66034987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4686761) q[2];
sx q[2];
rz(-1.1489931) q[2];
sx q[2];
rz(-2.2093175) q[2];
rz(2.9428234) q[3];
sx q[3];
rz(-1.1105744) q[3];
sx q[3];
rz(0.96536243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7070049) q[0];
sx q[0];
rz(-2.2362237) q[0];
sx q[0];
rz(2.7804651) q[0];
rz(1.4350285) q[1];
sx q[1];
rz(-1.3577434) q[1];
sx q[1];
rz(2.3235869) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0571787) q[0];
sx q[0];
rz(-1.1857496) q[0];
sx q[0];
rz(-1.8684698) q[0];
rz(-1.2641764) q[2];
sx q[2];
rz(-0.38796705) q[2];
sx q[2];
rz(-2.1167133) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.86707838) q[1];
sx q[1];
rz(-0.78650219) q[1];
sx q[1];
rz(0.92814501) q[1];
rz(-pi) q[2];
rz(1.9345476) q[3];
sx q[3];
rz(-1.517429) q[3];
sx q[3];
rz(1.3747017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.25257418) q[2];
sx q[2];
rz(-0.44712862) q[2];
sx q[2];
rz(-1.4206295) q[2];
rz(-1.8255) q[3];
sx q[3];
rz(-0.75794739) q[3];
sx q[3];
rz(2.7533598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4017568) q[0];
sx q[0];
rz(-0.49210423) q[0];
sx q[0];
rz(-0.87093583) q[0];
rz(-2.8254106) q[1];
sx q[1];
rz(-0.28156391) q[1];
sx q[1];
rz(0.2972163) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0270099) q[0];
sx q[0];
rz(-1.5073338) q[0];
sx q[0];
rz(-2.7364536) q[0];
rz(-pi) q[1];
x q[1];
rz(0.70298123) q[2];
sx q[2];
rz(-0.90548041) q[2];
sx q[2];
rz(-0.92232982) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3618468) q[1];
sx q[1];
rz(-1.5231113) q[1];
sx q[1];
rz(1.7576799) q[1];
rz(-pi) q[2];
rz(0.66629569) q[3];
sx q[3];
rz(-2.1994281) q[3];
sx q[3];
rz(2.6578238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7685984) q[2];
sx q[2];
rz(-0.82295376) q[2];
sx q[2];
rz(2.7139943) q[2];
rz(-1.9528495) q[3];
sx q[3];
rz(-0.62994981) q[3];
sx q[3];
rz(-0.52743131) q[3];
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
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4363842) q[0];
sx q[0];
rz(-1.5642865) q[0];
sx q[0];
rz(2.3676681) q[0];
rz(-2.4286843) q[1];
sx q[1];
rz(-2.1122825) q[1];
sx q[1];
rz(2.4598222) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6503158) q[0];
sx q[0];
rz(-0.33704764) q[0];
sx q[0];
rz(0.89462535) q[0];
rz(0.21708023) q[2];
sx q[2];
rz(-0.97765572) q[2];
sx q[2];
rz(2.8341688) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3817953) q[1];
sx q[1];
rz(-1.1361546) q[1];
sx q[1];
rz(2.1684907) q[1];
rz(-pi) q[2];
rz(1.9346312) q[3];
sx q[3];
rz(-1.893265) q[3];
sx q[3];
rz(-0.22842562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.008808) q[2];
sx q[2];
rz(-2.4286353) q[2];
sx q[2];
rz(2.183389) q[2];
rz(2.0751674) q[3];
sx q[3];
rz(-1.2833779) q[3];
sx q[3];
rz(-0.3796033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5565857) q[0];
sx q[0];
rz(-1.7946694) q[0];
sx q[0];
rz(-2.5812896) q[0];
rz(-2.141748) q[1];
sx q[1];
rz(-2.9381349) q[1];
sx q[1];
rz(-1.6220185) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1513838) q[0];
sx q[0];
rz(-2.4575893) q[0];
sx q[0];
rz(2.3955976) q[0];
rz(-pi) q[1];
rz(2.1278473) q[2];
sx q[2];
rz(-1.0336913) q[2];
sx q[2];
rz(2.277166) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6474364) q[1];
sx q[1];
rz(-1.372822) q[1];
sx q[1];
rz(0.76311771) q[1];
x q[2];
rz(-1.394746) q[3];
sx q[3];
rz(-2.5213443) q[3];
sx q[3];
rz(0.32905096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4524298) q[2];
sx q[2];
rz(-0.55183691) q[2];
sx q[2];
rz(2.9186644) q[2];
rz(0.034742268) q[3];
sx q[3];
rz(-1.3870753) q[3];
sx q[3];
rz(-3.070014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8054304) q[0];
sx q[0];
rz(-2.8283089) q[0];
sx q[0];
rz(-1.0700595) q[0];
rz(-1.7806212) q[1];
sx q[1];
rz(-2.762251) q[1];
sx q[1];
rz(1.2840575) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5491935) q[0];
sx q[0];
rz(-2.1201029) q[0];
sx q[0];
rz(2.4248289) q[0];
rz(-pi) q[1];
rz(-2.1806296) q[2];
sx q[2];
rz(-0.46513882) q[2];
sx q[2];
rz(-2.7679408) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8821841) q[1];
sx q[1];
rz(-2.5849197) q[1];
sx q[1];
rz(0.56682079) q[1];
rz(-pi) q[2];
rz(-2.1070126) q[3];
sx q[3];
rz(-1.3048733) q[3];
sx q[3];
rz(-1.769161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.47508919) q[2];
sx q[2];
rz(-1.6978426) q[2];
sx q[2];
rz(-2.5040023) q[2];
rz(-2.3049138) q[3];
sx q[3];
rz(-1.1058608) q[3];
sx q[3];
rz(1.3440514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5752983) q[0];
sx q[0];
rz(-0.42995444) q[0];
sx q[0];
rz(0.56754011) q[0];
rz(-0.42770806) q[1];
sx q[1];
rz(-1.5274915) q[1];
sx q[1];
rz(-2.2033851) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.339401) q[0];
sx q[0];
rz(-1.8759449) q[0];
sx q[0];
rz(-1.1423654) q[0];
rz(-pi) q[1];
rz(2.5222048) q[2];
sx q[2];
rz(-0.833138) q[2];
sx q[2];
rz(0.39820652) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2141014) q[1];
sx q[1];
rz(-0.92988211) q[1];
sx q[1];
rz(-0.44968857) q[1];
rz(0.35258099) q[3];
sx q[3];
rz(-2.4028824) q[3];
sx q[3];
rz(1.1383566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.87970916) q[2];
sx q[2];
rz(-1.9677013) q[2];
sx q[2];
rz(1.7555457) q[2];
rz(1.322768) q[3];
sx q[3];
rz(-1.1498007) q[3];
sx q[3];
rz(0.02903207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8742074) q[0];
sx q[0];
rz(-2.8081942) q[0];
sx q[0];
rz(-1.4338795) q[0];
rz(1.2738312) q[1];
sx q[1];
rz(-1.1359943) q[1];
sx q[1];
rz(0.83126718) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6600663) q[0];
sx q[0];
rz(-0.43018451) q[0];
sx q[0];
rz(-0.35920401) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7449042) q[2];
sx q[2];
rz(-2.499352) q[2];
sx q[2];
rz(-2.1070534) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2483406) q[1];
sx q[1];
rz(-0.89683956) q[1];
sx q[1];
rz(2.2158951) q[1];
x q[2];
rz(1.7940815) q[3];
sx q[3];
rz(-1.3459473) q[3];
sx q[3];
rz(0.68883483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5780118) q[2];
sx q[2];
rz(-1.3728023) q[2];
sx q[2];
rz(2.1772299) q[2];
rz(1.9780805) q[3];
sx q[3];
rz(-0.61146277) q[3];
sx q[3];
rz(2.4826629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37373856) q[0];
sx q[0];
rz(-2.4066194) q[0];
sx q[0];
rz(2.1642165) q[0];
rz(-1.7550229) q[1];
sx q[1];
rz(-1.8354548) q[1];
sx q[1];
rz(-1.1057373) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7439197) q[0];
sx q[0];
rz(-2.0154675) q[0];
sx q[0];
rz(1.4414653) q[0];
rz(-pi) q[1];
rz(-3.0366304) q[2];
sx q[2];
rz(-2.2347921) q[2];
sx q[2];
rz(-0.86519372) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9165009) q[1];
sx q[1];
rz(-1.8019925) q[1];
sx q[1];
rz(0.1901615) q[1];
rz(-0.68007277) q[3];
sx q[3];
rz(-2.701093) q[3];
sx q[3];
rz(2.3624453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.608312) q[2];
sx q[2];
rz(-2.7567342) q[2];
sx q[2];
rz(0.14979714) q[2];
rz(1.7685361) q[3];
sx q[3];
rz(-1.7356197) q[3];
sx q[3];
rz(-1.3214553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6479284) q[0];
sx q[0];
rz(-2.6239008) q[0];
sx q[0];
rz(3.0143484) q[0];
rz(1.4808902) q[1];
sx q[1];
rz(-2.7791185) q[1];
sx q[1];
rz(2.9737934) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.084899336) q[0];
sx q[0];
rz(-0.035824422) q[0];
sx q[0];
rz(-0.56532677) q[0];
x q[1];
rz(-0.22044214) q[2];
sx q[2];
rz(-0.61969212) q[2];
sx q[2];
rz(2.9881791) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.92229453) q[1];
sx q[1];
rz(-2.1592327) q[1];
sx q[1];
rz(1.5625619) q[1];
rz(-pi) q[2];
rz(-1.4526669) q[3];
sx q[3];
rz(-1.6932634) q[3];
sx q[3];
rz(2.1446705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1148791) q[2];
sx q[2];
rz(-2.202704) q[2];
sx q[2];
rz(2.3804469) q[2];
rz(0.090027697) q[3];
sx q[3];
rz(-1.0031676) q[3];
sx q[3];
rz(-0.95054039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.5933843) q[0];
sx q[0];
rz(-1.9822639) q[0];
sx q[0];
rz(-0.32250861) q[0];
rz(-2.7535915) q[1];
sx q[1];
rz(-1.7419659) q[1];
sx q[1];
rz(2.3566125) q[1];
rz(2.3315196) q[2];
sx q[2];
rz(-2.0740866) q[2];
sx q[2];
rz(-1.496051) q[2];
rz(1.2896982) q[3];
sx q[3];
rz(-0.55681183) q[3];
sx q[3];
rz(-1.1435215) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];