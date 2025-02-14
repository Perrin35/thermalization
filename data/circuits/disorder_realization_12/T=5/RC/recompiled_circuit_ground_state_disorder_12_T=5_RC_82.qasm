OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0491068) q[0];
sx q[0];
rz(4.4409039) q[0];
sx q[0];
rz(9.470603) q[0];
rz(1.2996281) q[1];
sx q[1];
rz(-1.8027432) q[1];
sx q[1];
rz(-1.8884698) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0631615) q[0];
sx q[0];
rz(-0.51126152) q[0];
sx q[0];
rz(1.0529165) q[0];
rz(-pi) q[1];
rz(0.21109493) q[2];
sx q[2];
rz(-2.935754) q[2];
sx q[2];
rz(-1.4343651) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3435077) q[1];
sx q[1];
rz(-1.4516593) q[1];
sx q[1];
rz(2.3878674) q[1];
rz(-pi) q[2];
rz(1.5398938) q[3];
sx q[3];
rz(-3.0140844) q[3];
sx q[3];
rz(2.7650583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8200298) q[2];
sx q[2];
rz(-1.4048615) q[2];
sx q[2];
rz(-2.826214) q[2];
rz(3.0553715) q[3];
sx q[3];
rz(-0.092970522) q[3];
sx q[3];
rz(0.84403795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.740199) q[0];
sx q[0];
rz(-1.711015) q[0];
sx q[0];
rz(-0.33541086) q[0];
rz(-0.3849349) q[1];
sx q[1];
rz(-2.3454869) q[1];
sx q[1];
rz(1.8523432) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5407643) q[0];
sx q[0];
rz(-2.7473831) q[0];
sx q[0];
rz(-2.9790131) q[0];
x q[1];
rz(1.4440437) q[2];
sx q[2];
rz(-0.64536649) q[2];
sx q[2];
rz(-2.5996948) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.34270479) q[1];
sx q[1];
rz(-1.2884166) q[1];
sx q[1];
rz(-0.39816959) q[1];
rz(-3.0063674) q[3];
sx q[3];
rz(-1.3298508) q[3];
sx q[3];
rz(-1.8640445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5424767) q[2];
sx q[2];
rz(-1.7581538) q[2];
sx q[2];
rz(-1.4625589) q[2];
rz(-0.8820495) q[3];
sx q[3];
rz(-2.4817011) q[3];
sx q[3];
rz(-1.5006458) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3682692) q[0];
sx q[0];
rz(-0.48650807) q[0];
sx q[0];
rz(2.4988556) q[0];
rz(0.87359387) q[1];
sx q[1];
rz(-1.830955) q[1];
sx q[1];
rz(-2.1547623) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29016748) q[0];
sx q[0];
rz(-2.254524) q[0];
sx q[0];
rz(-0.19398035) q[0];
rz(-0.061927783) q[2];
sx q[2];
rz(-1.0871097) q[2];
sx q[2];
rz(2.7360423) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.75007406) q[1];
sx q[1];
rz(-1.0315391) q[1];
sx q[1];
rz(0.92309322) q[1];
x q[2];
rz(1.1961206) q[3];
sx q[3];
rz(-1.5244532) q[3];
sx q[3];
rz(2.1425193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4387536) q[2];
sx q[2];
rz(-0.025391014) q[2];
sx q[2];
rz(-1.0566443) q[2];
rz(1.7993401) q[3];
sx q[3];
rz(-1.6868351) q[3];
sx q[3];
rz(-2.2630283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.037828) q[0];
sx q[0];
rz(-2.8570638) q[0];
sx q[0];
rz(1.1430662) q[0];
rz(-2.4567538) q[1];
sx q[1];
rz(-1.3416483) q[1];
sx q[1];
rz(-0.24982223) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.095465) q[0];
sx q[0];
rz(-0.85554142) q[0];
sx q[0];
rz(-2.1939795) q[0];
x q[1];
rz(-0.5742214) q[2];
sx q[2];
rz(-1.5674233) q[2];
sx q[2];
rz(0.33004883) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.494809) q[1];
sx q[1];
rz(-2.3338291) q[1];
sx q[1];
rz(-0.69112063) q[1];
x q[2];
rz(1.5644844) q[3];
sx q[3];
rz(-2.8996116) q[3];
sx q[3];
rz(2.210287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7217241) q[2];
sx q[2];
rz(-1.7896174) q[2];
sx q[2];
rz(-1.295759) q[2];
rz(1.7660247) q[3];
sx q[3];
rz(-1.4294759) q[3];
sx q[3];
rz(0.099055722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7527723) q[0];
sx q[0];
rz(-1.7701912) q[0];
sx q[0];
rz(0.8465299) q[0];
rz(1.1835774) q[1];
sx q[1];
rz(-1.1776244) q[1];
sx q[1];
rz(1.2394946) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13279937) q[0];
sx q[0];
rz(-1.6911048) q[0];
sx q[0];
rz(-1.8254721) q[0];
rz(-0.34244142) q[2];
sx q[2];
rz(-1.3991809) q[2];
sx q[2];
rz(-1.1521074) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.2787362) q[1];
sx q[1];
rz(-1.6267611) q[1];
sx q[1];
rz(-0.011933283) q[1];
rz(-pi) q[2];
rz(0.059885177) q[3];
sx q[3];
rz(-0.97447534) q[3];
sx q[3];
rz(-2.8879762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6401297) q[2];
sx q[2];
rz(-1.5093466) q[2];
sx q[2];
rz(2.787369) q[2];
rz(-0.7192449) q[3];
sx q[3];
rz(-0.57643276) q[3];
sx q[3];
rz(0.80938068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(1.3767913) q[0];
sx q[0];
rz(-2.9053595) q[0];
sx q[0];
rz(-0.34673044) q[0];
rz(2.4193343) q[1];
sx q[1];
rz(-1.5303231) q[1];
sx q[1];
rz(2.9647656) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.051561616) q[0];
sx q[0];
rz(-1.430184) q[0];
sx q[0];
rz(-2.0363801) q[0];
x q[1];
rz(-0.093164154) q[2];
sx q[2];
rz(-1.0736862) q[2];
sx q[2];
rz(0.99592956) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.74048282) q[1];
sx q[1];
rz(-1.4663457) q[1];
sx q[1];
rz(-2.9354276) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0079871) q[3];
sx q[3];
rz(-1.4848108) q[3];
sx q[3];
rz(-0.028349625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5727545) q[2];
sx q[2];
rz(-2.4080031) q[2];
sx q[2];
rz(1.0432165) q[2];
rz(2.9853232) q[3];
sx q[3];
rz(-2.0782491) q[3];
sx q[3];
rz(3.0472896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8909376) q[0];
sx q[0];
rz(-1.0303048) q[0];
sx q[0];
rz(-1.2283196) q[0];
rz(-2.7078775) q[1];
sx q[1];
rz(-2.3408196) q[1];
sx q[1];
rz(1.0450276) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21711536) q[0];
sx q[0];
rz(-1.7733367) q[0];
sx q[0];
rz(-3.0826352) q[0];
rz(-2.074998) q[2];
sx q[2];
rz(-2.3269751) q[2];
sx q[2];
rz(-2.2058918) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8974287) q[1];
sx q[1];
rz(-0.62627316) q[1];
sx q[1];
rz(1.3757964) q[1];
x q[2];
rz(0.52900903) q[3];
sx q[3];
rz(-0.85011357) q[3];
sx q[3];
rz(-1.0320272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9016483) q[2];
sx q[2];
rz(-2.1176691) q[2];
sx q[2];
rz(0.18590064) q[2];
rz(1.2402041) q[3];
sx q[3];
rz(-1.5408206) q[3];
sx q[3];
rz(1.0143636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0094725322) q[0];
sx q[0];
rz(-1.8931696) q[0];
sx q[0];
rz(-0.15952071) q[0];
rz(-2.3347143) q[1];
sx q[1];
rz(-1.9069549) q[1];
sx q[1];
rz(0.0085011403) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59183622) q[0];
sx q[0];
rz(-1.7248132) q[0];
sx q[0];
rz(-1.1598253) q[0];
rz(-3.0372856) q[2];
sx q[2];
rz(-2.3416714) q[2];
sx q[2];
rz(1.6751573) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9125625) q[1];
sx q[1];
rz(-2.8797315) q[1];
sx q[1];
rz(-2.6109241) q[1];
rz(-pi) q[2];
rz(-2.9656677) q[3];
sx q[3];
rz(-1.0930382) q[3];
sx q[3];
rz(-2.0328409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8622417) q[2];
sx q[2];
rz(-1.778435) q[2];
sx q[2];
rz(2.9478493) q[2];
rz(1.5396317) q[3];
sx q[3];
rz(-1.1774457) q[3];
sx q[3];
rz(-2.0120473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2368471) q[0];
sx q[0];
rz(-1.756825) q[0];
sx q[0];
rz(2.4054476) q[0];
rz(2.2537117) q[1];
sx q[1];
rz(-2.6887951) q[1];
sx q[1];
rz(3.091541) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0808635) q[0];
sx q[0];
rz(-1.9115149) q[0];
sx q[0];
rz(1.7572035) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8456837) q[2];
sx q[2];
rz(-2.8136133) q[2];
sx q[2];
rz(0.28959549) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.012573012) q[1];
sx q[1];
rz(-2.7677849) q[1];
sx q[1];
rz(-0.5859979) q[1];
rz(-2.3042559) q[3];
sx q[3];
rz(-0.76238197) q[3];
sx q[3];
rz(-1.1632944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8686409) q[2];
sx q[2];
rz(-1.7014039) q[2];
sx q[2];
rz(3.1325373) q[2];
rz(2.791259) q[3];
sx q[3];
rz(-0.30600268) q[3];
sx q[3];
rz(-0.30456021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41651273) q[0];
sx q[0];
rz(-2.0820936) q[0];
sx q[0];
rz(-0.9730202) q[0];
rz(-1.5618207) q[1];
sx q[1];
rz(-2.2353954) q[1];
sx q[1];
rz(0.80950338) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0334629) q[0];
sx q[0];
rz(-1.2133043) q[0];
sx q[0];
rz(-2.3603159) q[0];
x q[1];
rz(1.8000748) q[2];
sx q[2];
rz(-1.3827168) q[2];
sx q[2];
rz(-0.089872472) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1228408) q[1];
sx q[1];
rz(-1.1114618) q[1];
sx q[1];
rz(0.18045119) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6736835) q[3];
sx q[3];
rz(-2.10026) q[3];
sx q[3];
rz(0.61509815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1805264) q[2];
sx q[2];
rz(-2.3193391) q[2];
sx q[2];
rz(2.6194438) q[2];
rz(-2.7850049) q[3];
sx q[3];
rz(-0.57727376) q[3];
sx q[3];
rz(0.063551158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6575573) q[0];
sx q[0];
rz(-0.35150305) q[0];
sx q[0];
rz(1.2765314) q[0];
rz(-0.36318489) q[1];
sx q[1];
rz(-2.6138432) q[1];
sx q[1];
rz(1.6539727) q[1];
rz(1.034201) q[2];
sx q[2];
rz(-2.2727439) q[2];
sx q[2];
rz(-1.5287413) q[2];
rz(-0.47388062) q[3];
sx q[3];
rz(-0.6345748) q[3];
sx q[3];
rz(2.0848772) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
