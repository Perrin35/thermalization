OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6808074) q[0];
sx q[0];
rz(2.1587125) q[0];
sx q[0];
rz(8.4150253) q[0];
rz(2.9653964) q[1];
sx q[1];
rz(5.3806452) q[1];
sx q[1];
rz(7.5723958) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6362808) q[0];
sx q[0];
rz(-1.2828151) q[0];
sx q[0];
rz(-2.050839) q[0];
x q[1];
rz(-0.024359811) q[2];
sx q[2];
rz(-1.5775541) q[2];
sx q[2];
rz(-0.18505219) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7154555) q[1];
sx q[1];
rz(-1.8278367) q[1];
sx q[1];
rz(0.21547683) q[1];
x q[2];
rz(1.8208002) q[3];
sx q[3];
rz(-2.8651926) q[3];
sx q[3];
rz(-1.1652511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5916799) q[2];
sx q[2];
rz(-2.9788571) q[2];
sx q[2];
rz(-0.13554779) q[2];
rz(3.0013951) q[3];
sx q[3];
rz(-2.086816) q[3];
sx q[3];
rz(-2.8285817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7267589) q[0];
sx q[0];
rz(-2.2048075) q[0];
sx q[0];
rz(2.360789) q[0];
rz(-2.8813598) q[1];
sx q[1];
rz(-0.5865016) q[1];
sx q[1];
rz(1.3134726) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1056054) q[0];
sx q[0];
rz(-0.30788883) q[0];
sx q[0];
rz(2.3769828) q[0];
rz(-0.30828373) q[2];
sx q[2];
rz(-1.1051902) q[2];
sx q[2];
rz(-0.75454933) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5116509) q[1];
sx q[1];
rz(-2.959046) q[1];
sx q[1];
rz(-0.49717848) q[1];
x q[2];
rz(3.0659862) q[3];
sx q[3];
rz(-1.4652271) q[3];
sx q[3];
rz(-2.3188863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.446622) q[2];
sx q[2];
rz(-1.5335252) q[2];
sx q[2];
rz(0.95412811) q[2];
rz(-1.4387087) q[3];
sx q[3];
rz(-1.8372767) q[3];
sx q[3];
rz(-0.71189705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.057782877) q[0];
sx q[0];
rz(-2.5397781) q[0];
sx q[0];
rz(-2.4457248) q[0];
rz(-0.35481915) q[1];
sx q[1];
rz(-1.6638919) q[1];
sx q[1];
rz(0.16608873) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9244708) q[0];
sx q[0];
rz(-0.32453254) q[0];
sx q[0];
rz(-1.0990418) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.66785779) q[2];
sx q[2];
rz(-0.72615004) q[2];
sx q[2];
rz(0.49612507) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3661763) q[1];
sx q[1];
rz(-2.011538) q[1];
sx q[1];
rz(0.52904769) q[1];
rz(-pi) q[2];
rz(0.8576287) q[3];
sx q[3];
rz(-1.3751251) q[3];
sx q[3];
rz(1.5426202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.35963905) q[2];
sx q[2];
rz(-2.226318) q[2];
sx q[2];
rz(0.98177838) q[2];
rz(-2.0180457) q[3];
sx q[3];
rz(-0.26505622) q[3];
sx q[3];
rz(1.3004998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.05900255) q[0];
sx q[0];
rz(-2.1868717) q[0];
sx q[0];
rz(0.29378763) q[0];
rz(0.44149533) q[1];
sx q[1];
rz(-1.0499294) q[1];
sx q[1];
rz(-2.3667483) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42472096) q[0];
sx q[0];
rz(-1.6443335) q[0];
sx q[0];
rz(-1.5872692) q[0];
rz(-pi) q[1];
rz(-0.11544322) q[2];
sx q[2];
rz(-2.7713159) q[2];
sx q[2];
rz(-1.8625129) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6619819) q[1];
sx q[1];
rz(-1.4927215) q[1];
sx q[1];
rz(1.1675203) q[1];
rz(2.9702441) q[3];
sx q[3];
rz(-1.8257986) q[3];
sx q[3];
rz(1.9920497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1321156) q[2];
sx q[2];
rz(-2.0602132) q[2];
sx q[2];
rz(-0.57007989) q[2];
rz(1.305497) q[3];
sx q[3];
rz(-0.92731849) q[3];
sx q[3];
rz(-1.501804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.775979) q[0];
sx q[0];
rz(-1.403084) q[0];
sx q[0];
rz(2.9602125) q[0];
rz(-0.60896215) q[1];
sx q[1];
rz(-0.47157559) q[1];
sx q[1];
rz(0.10770527) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90751326) q[0];
sx q[0];
rz(-0.36511746) q[0];
sx q[0];
rz(-1.8318729) q[0];
rz(-pi) q[1];
rz(-2.2568251) q[2];
sx q[2];
rz(-0.77026412) q[2];
sx q[2];
rz(-1.9212854) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4957461) q[1];
sx q[1];
rz(-1.1762113) q[1];
sx q[1];
rz(-0.87162019) q[1];
rz(-pi) q[2];
rz(-1.0531353) q[3];
sx q[3];
rz(-0.62873757) q[3];
sx q[3];
rz(-2.351458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.084215) q[2];
sx q[2];
rz(-1.4930909) q[2];
sx q[2];
rz(0.28953141) q[2];
rz(-3.1048408) q[3];
sx q[3];
rz(-0.74220243) q[3];
sx q[3];
rz(-0.010820476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.058379563) q[0];
sx q[0];
rz(-0.20861861) q[0];
sx q[0];
rz(-1.5104729) q[0];
rz(2.0026813) q[1];
sx q[1];
rz(-1.6612256) q[1];
sx q[1];
rz(1.0169792) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6769584) q[0];
sx q[0];
rz(-2.5559253) q[0];
sx q[0];
rz(2.2023489) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.20118841) q[2];
sx q[2];
rz(-1.5474461) q[2];
sx q[2];
rz(-0.77229283) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.79343866) q[1];
sx q[1];
rz(-0.70421709) q[1];
sx q[1];
rz(-0.80146316) q[1];
rz(0.53652699) q[3];
sx q[3];
rz(-1.121472) q[3];
sx q[3];
rz(2.7262053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.54414526) q[2];
sx q[2];
rz(-2.1400673) q[2];
sx q[2];
rz(2.5069359) q[2];
rz(1.0774111) q[3];
sx q[3];
rz(-1.0297188) q[3];
sx q[3];
rz(-2.6950148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25062659) q[0];
sx q[0];
rz(-0.64193305) q[0];
sx q[0];
rz(-2.2977258) q[0];
rz(1.6183052) q[1];
sx q[1];
rz(-2.4089676) q[1];
sx q[1];
rz(1.2197781) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49918136) q[0];
sx q[0];
rz(-1.6089604) q[0];
sx q[0];
rz(-0.092060815) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2336823) q[2];
sx q[2];
rz(-2.6320576) q[2];
sx q[2];
rz(-2.8223035) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.68901686) q[1];
sx q[1];
rz(-0.92454443) q[1];
sx q[1];
rz(-0.31022443) q[1];
rz(-pi) q[2];
rz(1.3656093) q[3];
sx q[3];
rz(-2.725051) q[3];
sx q[3];
rz(-0.52683559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.62961489) q[2];
sx q[2];
rz(-0.93705606) q[2];
sx q[2];
rz(2.8249595) q[2];
rz(-3.1043502) q[3];
sx q[3];
rz(-1.3469632) q[3];
sx q[3];
rz(-0.75604701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(2.1036296) q[0];
sx q[0];
rz(-1.0873955) q[0];
sx q[0];
rz(0.33139247) q[0];
rz(-1.1720852) q[1];
sx q[1];
rz(-2.4786699) q[1];
sx q[1];
rz(0.40922871) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0497826) q[0];
sx q[0];
rz(-0.70735332) q[0];
sx q[0];
rz(1.1298256) q[0];
rz(-pi) q[1];
rz(-2.6166603) q[2];
sx q[2];
rz(-1.2972752) q[2];
sx q[2];
rz(0.33821019) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.85247773) q[1];
sx q[1];
rz(-2.5534938) q[1];
sx q[1];
rz(0.27681338) q[1];
x q[2];
rz(-2.6243022) q[3];
sx q[3];
rz(-1.358277) q[3];
sx q[3];
rz(1.1503435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.99872148) q[2];
sx q[2];
rz(-1.8565535) q[2];
sx q[2];
rz(0.55465737) q[2];
rz(-0.64152843) q[3];
sx q[3];
rz(-2.9252164) q[3];
sx q[3];
rz(0.085111246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2239969) q[0];
sx q[0];
rz(-1.6308835) q[0];
sx q[0];
rz(-3.1316485) q[0];
rz(2.0857281) q[1];
sx q[1];
rz(-1.5751585) q[1];
sx q[1];
rz(-1.6329637) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5976534) q[0];
sx q[0];
rz(-0.88730747) q[0];
sx q[0];
rz(0.20792122) q[0];
x q[1];
rz(-1.6998859) q[2];
sx q[2];
rz(-2.522905) q[2];
sx q[2];
rz(-0.29030756) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.29533169) q[1];
sx q[1];
rz(-1.2519072) q[1];
sx q[1];
rz(-2.1178513) q[1];
rz(-3.085) q[3];
sx q[3];
rz(-0.86846272) q[3];
sx q[3];
rz(-3.1346723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9987954) q[2];
sx q[2];
rz(-1.5441511) q[2];
sx q[2];
rz(2.7049086) q[2];
rz(-1.8113332) q[3];
sx q[3];
rz(-0.67970777) q[3];
sx q[3];
rz(-2.1127624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0325539) q[0];
sx q[0];
rz(-0.80820307) q[0];
sx q[0];
rz(2.57634) q[0];
rz(2.8887707) q[1];
sx q[1];
rz(-2.1222474) q[1];
sx q[1];
rz(1.3814829) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0927267) q[0];
sx q[0];
rz(-0.5038358) q[0];
sx q[0];
rz(-3.0548274) q[0];
x q[1];
rz(1.6818468) q[2];
sx q[2];
rz(-2.0425218) q[2];
sx q[2];
rz(1.9327088) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6773721) q[1];
sx q[1];
rz(-2.3090625) q[1];
sx q[1];
rz(2.1470451) q[1];
x q[2];
rz(-2.284427) q[3];
sx q[3];
rz(-1.790159) q[3];
sx q[3];
rz(2.9227748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.4505724) q[2];
sx q[2];
rz(-2.317163) q[2];
sx q[2];
rz(-0.59664574) q[2];
rz(-0.48554844) q[3];
sx q[3];
rz(-0.89533007) q[3];
sx q[3];
rz(-0.027464494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(2.6116466) q[0];
sx q[0];
rz(-1.2048789) q[0];
sx q[0];
rz(-2.0299029) q[0];
rz(-1.8687517) q[1];
sx q[1];
rz(-1.0796937) q[1];
sx q[1];
rz(-0.83273522) q[1];
rz(-1.7455208) q[2];
sx q[2];
rz(-0.25824418) q[2];
sx q[2];
rz(1.2373409) q[2];
rz(-2.9968895) q[3];
sx q[3];
rz(-0.51454138) q[3];
sx q[3];
rz(0.71458057) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
