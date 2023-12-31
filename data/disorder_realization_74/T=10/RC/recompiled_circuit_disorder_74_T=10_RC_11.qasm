OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.4607853) q[0];
sx q[0];
rz(-2.1587125) q[0];
sx q[0];
rz(2.13184) q[0];
rz(-0.17619625) q[1];
sx q[1];
rz(-2.2390525) q[1];
sx q[1];
rz(-1.2892105) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5754817) q[0];
sx q[0];
rz(-2.5876343) q[0];
sx q[0];
rz(-2.1411116) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1172328) q[2];
sx q[2];
rz(-1.5775541) q[2];
sx q[2];
rz(-0.18505219) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8569021) q[1];
sx q[1];
rz(-2.8077217) q[1];
sx q[1];
rz(2.2536709) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3025769) q[3];
sx q[3];
rz(-1.5032288) q[3];
sx q[3];
rz(-0.64642954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5499128) q[2];
sx q[2];
rz(-0.16273558) q[2];
sx q[2];
rz(-0.13554779) q[2];
rz(-3.0013951) q[3];
sx q[3];
rz(-2.086816) q[3];
sx q[3];
rz(-0.31301096) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41483375) q[0];
sx q[0];
rz(-0.93678513) q[0];
sx q[0];
rz(0.78080368) q[0];
rz(-0.26023284) q[1];
sx q[1];
rz(-2.5550911) q[1];
sx q[1];
rz(-1.82812) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86620264) q[0];
sx q[0];
rz(-1.35944) q[0];
sx q[0];
rz(-2.9160121) q[0];
rz(2.1140695) q[2];
sx q[2];
rz(-2.5894905) q[2];
sx q[2];
rz(0.13763025) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.62994176) q[1];
sx q[1];
rz(-0.18254666) q[1];
sx q[1];
rz(2.6444142) q[1];
x q[2];
rz(1.6766657) q[3];
sx q[3];
rz(-1.4956116) q[3];
sx q[3];
rz(2.4014846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.446622) q[2];
sx q[2];
rz(-1.6080674) q[2];
sx q[2];
rz(-2.1874645) q[2];
rz(-1.4387087) q[3];
sx q[3];
rz(-1.3043159) q[3];
sx q[3];
rz(-2.4296956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0838098) q[0];
sx q[0];
rz(-2.5397781) q[0];
sx q[0];
rz(2.4457248) q[0];
rz(-2.7867735) q[1];
sx q[1];
rz(-1.4777007) q[1];
sx q[1];
rz(-2.9755039) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27667339) q[0];
sx q[0];
rz(-1.2827946) q[0];
sx q[0];
rz(-2.9898781) q[0];
rz(-pi) q[1];
rz(-1.0679929) q[2];
sx q[2];
rz(-1.0223801) q[2];
sx q[2];
rz(1.3082248) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7158311) q[1];
sx q[1];
rz(-2.466723) q[1];
sx q[1];
rz(2.3900044) q[1];
rz(-pi) q[2];
rz(-0.2563128) q[3];
sx q[3];
rz(-0.87402065) q[3];
sx q[3];
rz(-0.13845201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7819536) q[2];
sx q[2];
rz(-2.226318) q[2];
sx q[2];
rz(2.1598143) q[2];
rz(2.0180457) q[3];
sx q[3];
rz(-2.8765364) q[3];
sx q[3];
rz(1.3004998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0825901) q[0];
sx q[0];
rz(-2.1868717) q[0];
sx q[0];
rz(0.29378763) q[0];
rz(-2.7000973) q[1];
sx q[1];
rz(-1.0499294) q[1];
sx q[1];
rz(-2.3667483) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4962909) q[0];
sx q[0];
rz(-3.0662363) q[0];
sx q[0];
rz(-2.9216179) q[0];
x q[1];
rz(0.36802937) q[2];
sx q[2];
rz(-1.6124915) q[2];
sx q[2];
rz(-2.7421943) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.47961071) q[1];
sx q[1];
rz(-1.4927215) q[1];
sx q[1];
rz(1.1675203) q[1];
rz(-pi) q[2];
rz(1.8294228) q[3];
sx q[3];
rz(-1.4050409) q[3];
sx q[3];
rz(2.7639619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.1321156) q[2];
sx q[2];
rz(-1.0813794) q[2];
sx q[2];
rz(-2.5715128) q[2];
rz(1.8360957) q[3];
sx q[3];
rz(-0.92731849) q[3];
sx q[3];
rz(1.501804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.775979) q[0];
sx q[0];
rz(-1.403084) q[0];
sx q[0];
rz(0.18138012) q[0];
rz(0.60896215) q[1];
sx q[1];
rz(-0.47157559) q[1];
sx q[1];
rz(-0.10770527) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2337423) q[0];
sx q[0];
rz(-1.6630917) q[0];
sx q[0];
rz(1.924563) q[0];
x q[1];
rz(2.5905215) q[2];
sx q[2];
rz(-2.1398009) q[2];
sx q[2];
rz(2.0713194) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9039893) q[1];
sx q[1];
rz(-0.9346107) q[1];
sx q[1];
rz(2.643305) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1343469) q[3];
sx q[3];
rz(-1.2754903) q[3];
sx q[3];
rz(1.9293279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.057377664) q[2];
sx q[2];
rz(-1.6485018) q[2];
sx q[2];
rz(-2.8520612) q[2];
rz(0.036751898) q[3];
sx q[3];
rz(-0.74220243) q[3];
sx q[3];
rz(-0.010820476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0832131) q[0];
sx q[0];
rz(-2.932974) q[0];
sx q[0];
rz(-1.6311197) q[0];
rz(-1.1389114) q[1];
sx q[1];
rz(-1.4803671) q[1];
sx q[1];
rz(-1.0169792) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4879918) q[0];
sx q[0];
rz(-1.2383608) q[0];
sx q[0];
rz(1.0792653) q[0];
rz(-0.20118841) q[2];
sx q[2];
rz(-1.5941465) q[2];
sx q[2];
rz(0.77229283) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.79343866) q[1];
sx q[1];
rz(-2.4373756) q[1];
sx q[1];
rz(2.3401295) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0820886) q[3];
sx q[3];
rz(-2.0492616) q[3];
sx q[3];
rz(2.2389776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5974474) q[2];
sx q[2];
rz(-1.0015254) q[2];
sx q[2];
rz(2.5069359) q[2];
rz(2.0641816) q[3];
sx q[3];
rz(-1.0297188) q[3];
sx q[3];
rz(-0.44657782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(0.25062659) q[0];
sx q[0];
rz(-2.4996596) q[0];
sx q[0];
rz(2.2977258) q[0];
rz(1.6183052) q[1];
sx q[1];
rz(-2.4089676) q[1];
sx q[1];
rz(1.2197781) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49918136) q[0];
sx q[0];
rz(-1.5326323) q[0];
sx q[0];
rz(-0.092060815) q[0];
x q[1];
rz(2.2336823) q[2];
sx q[2];
rz(-0.50953509) q[2];
sx q[2];
rz(2.8223035) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0724832) q[1];
sx q[1];
rz(-1.3246037) q[1];
sx q[1];
rz(0.90087955) q[1];
rz(-pi) q[2];
rz(3.0516902) q[3];
sx q[3];
rz(-1.9780759) q[3];
sx q[3];
rz(-2.390993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.62961489) q[2];
sx q[2];
rz(-2.2045366) q[2];
sx q[2];
rz(-2.8249595) q[2];
rz(3.1043502) q[3];
sx q[3];
rz(-1.7946295) q[3];
sx q[3];
rz(-0.75604701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0379631) q[0];
sx q[0];
rz(-2.0541971) q[0];
sx q[0];
rz(-2.8102002) q[0];
rz(1.1720852) q[1];
sx q[1];
rz(-2.4786699) q[1];
sx q[1];
rz(2.7323639) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17657875) q[0];
sx q[0];
rz(-1.2897549) q[0];
sx q[0];
rz(-0.9126419) q[0];
rz(-pi) q[1];
rz(-1.2572844) q[2];
sx q[2];
rz(-1.0672617) q[2];
sx q[2];
rz(1.3877649) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1815391) q[1];
sx q[1];
rz(-1.0078733) q[1];
sx q[1];
rz(-1.3905418) q[1];
rz(-pi) q[2];
rz(-2.6243022) q[3];
sx q[3];
rz(-1.358277) q[3];
sx q[3];
rz(1.1503435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1428712) q[2];
sx q[2];
rz(-1.2850392) q[2];
sx q[2];
rz(-0.55465737) q[2];
rz(-0.64152843) q[3];
sx q[3];
rz(-2.9252164) q[3];
sx q[3];
rz(-3.0564814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.2239969) q[0];
sx q[0];
rz(-1.5107091) q[0];
sx q[0];
rz(-0.0099442033) q[0];
rz(2.0857281) q[1];
sx q[1];
rz(-1.5751585) q[1];
sx q[1];
rz(1.508629) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54393923) q[0];
sx q[0];
rz(-2.2542852) q[0];
sx q[0];
rz(-2.9336714) q[0];
x q[1];
rz(0.091392322) q[2];
sx q[2];
rz(-2.1835727) q[2];
sx q[2];
rz(-0.44832715) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0548045) q[1];
sx q[1];
rz(-1.0541704) q[1];
sx q[1];
rz(-2.7727142) q[1];
rz(-pi) q[2];
rz(-0.056592654) q[3];
sx q[3];
rz(-0.86846272) q[3];
sx q[3];
rz(-0.0069204023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.14279723) q[2];
sx q[2];
rz(-1.5441511) q[2];
sx q[2];
rz(-2.7049086) q[2];
rz(1.8113332) q[3];
sx q[3];
rz(-2.4618849) q[3];
sx q[3];
rz(1.0288303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(-0.10903877) q[0];
sx q[0];
rz(-0.80820307) q[0];
sx q[0];
rz(2.57634) q[0];
rz(-2.8887707) q[1];
sx q[1];
rz(-2.1222474) q[1];
sx q[1];
rz(1.7601097) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1917282) q[0];
sx q[0];
rz(-1.0690332) q[0];
sx q[0];
rz(1.6185332) q[0];
rz(-2.9276766) q[2];
sx q[2];
rz(-2.6579318) q[2];
sx q[2];
rz(0.96825251) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.30365147) q[1];
sx q[1];
rz(-0.90180574) q[1];
sx q[1];
rz(2.6020357) q[1];
x q[2];
rz(-0.28678203) q[3];
sx q[3];
rz(-0.87773318) q[3];
sx q[3];
rz(-1.9758488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6910203) q[2];
sx q[2];
rz(-0.82442966) q[2];
sx q[2];
rz(2.5449469) q[2];
rz(-0.48554844) q[3];
sx q[3];
rz(-2.2462626) q[3];
sx q[3];
rz(-3.1141282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52994603) q[0];
sx q[0];
rz(-1.9367138) q[0];
sx q[0];
rz(1.1116897) q[0];
rz(-1.8687517) q[1];
sx q[1];
rz(-1.0796937) q[1];
sx q[1];
rz(-0.83273522) q[1];
rz(1.3163153) q[2];
sx q[2];
rz(-1.5263867) q[2];
sx q[2];
rz(2.9771794) q[2];
rz(1.4894555) q[3];
sx q[3];
rz(-1.0621536) q[3];
sx q[3];
rz(-2.5928706) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
