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
rz(-0.98288012) q[0];
sx q[0];
rz(-2.13184) q[0];
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
rz(-1.9292944) q[0];
sx q[0];
rz(-2.0295143) q[0];
sx q[0];
rz(2.8192768) q[0];
x q[1];
rz(3.1172328) q[2];
sx q[2];
rz(-1.5775541) q[2];
sx q[2];
rz(2.9565405) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.28469052) q[1];
sx q[1];
rz(-2.8077217) q[1];
sx q[1];
rz(2.2536709) q[1];
rz(-pi) q[2];
rz(-1.3207924) q[3];
sx q[3];
rz(-2.8651926) q[3];
sx q[3];
rz(-1.1652511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5916799) q[2];
sx q[2];
rz(-2.9788571) q[2];
sx q[2];
rz(0.13554779) q[2];
rz(3.0013951) q[3];
sx q[3];
rz(-1.0547767) q[3];
sx q[3];
rz(2.8285817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41483375) q[0];
sx q[0];
rz(-0.93678513) q[0];
sx q[0];
rz(-0.78080368) q[0];
rz(0.26023284) q[1];
sx q[1];
rz(-2.5550911) q[1];
sx q[1];
rz(-1.3134726) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.03598729) q[0];
sx q[0];
rz(-0.30788883) q[0];
sx q[0];
rz(0.76460989) q[0];
rz(-0.30828373) q[2];
sx q[2];
rz(-2.0364025) q[2];
sx q[2];
rz(0.75454933) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.62994176) q[1];
sx q[1];
rz(-0.18254666) q[1];
sx q[1];
rz(0.49717848) q[1];
x q[2];
rz(-3.0659862) q[3];
sx q[3];
rz(-1.6763655) q[3];
sx q[3];
rz(0.82270634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.69497067) q[2];
sx q[2];
rz(-1.5335252) q[2];
sx q[2];
rz(-0.95412811) q[2];
rz(-1.4387087) q[3];
sx q[3];
rz(-1.8372767) q[3];
sx q[3];
rz(-0.71189705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.057782877) q[0];
sx q[0];
rz(-2.5397781) q[0];
sx q[0];
rz(-2.4457248) q[0];
rz(0.35481915) q[1];
sx q[1];
rz(-1.4777007) q[1];
sx q[1];
rz(0.16608873) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8649193) q[0];
sx q[0];
rz(-1.858798) q[0];
sx q[0];
rz(-2.9898781) q[0];
rz(2.5327352) q[2];
sx q[2];
rz(-1.1470084) q[2];
sx q[2];
rz(2.599803) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7158311) q[1];
sx q[1];
rz(-0.67486963) q[1];
sx q[1];
rz(0.75158822) q[1];
x q[2];
rz(2.8852799) q[3];
sx q[3];
rz(-0.87402065) q[3];
sx q[3];
rz(-0.13845201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.35963905) q[2];
sx q[2];
rz(-0.91527462) q[2];
sx q[2];
rz(2.1598143) q[2];
rz(2.0180457) q[3];
sx q[3];
rz(-0.26505622) q[3];
sx q[3];
rz(-1.3004998) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05900255) q[0];
sx q[0];
rz(-2.1868717) q[0];
sx q[0];
rz(-0.29378763) q[0];
rz(2.7000973) q[1];
sx q[1];
rz(-2.0916633) q[1];
sx q[1];
rz(-2.3667483) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9967277) q[0];
sx q[0];
rz(-1.5872247) q[0];
sx q[0];
rz(3.0680455) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0261494) q[2];
sx q[2];
rz(-0.37027678) q[2];
sx q[2];
rz(1.8625129) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.47961071) q[1];
sx q[1];
rz(-1.6488711) q[1];
sx q[1];
rz(1.9740723) q[1];
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
rz(0.0094770771) q[2];
sx q[2];
rz(-2.0602132) q[2];
sx q[2];
rz(0.57007989) q[2];
rz(1.305497) q[3];
sx q[3];
rz(-0.92731849) q[3];
sx q[3];
rz(-1.501804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.775979) q[0];
sx q[0];
rz(-1.403084) q[0];
sx q[0];
rz(-2.9602125) q[0];
rz(2.5326305) q[1];
sx q[1];
rz(-2.6700171) q[1];
sx q[1];
rz(-0.10770527) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90751326) q[0];
sx q[0];
rz(-2.7764752) q[0];
sx q[0];
rz(1.8318729) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2147419) q[2];
sx q[2];
rz(-2.0276208) q[2];
sx q[2];
rz(-2.9608375) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4957461) q[1];
sx q[1];
rz(-1.9653814) q[1];
sx q[1];
rz(-2.2699725) q[1];
rz(-pi) q[2];
rz(-2.7961736) q[3];
sx q[3];
rz(-1.0343699) q[3];
sx q[3];
rz(-0.1766583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.057377664) q[2];
sx q[2];
rz(-1.6485018) q[2];
sx q[2];
rz(-2.8520612) q[2];
rz(-0.036751898) q[3];
sx q[3];
rz(-0.74220243) q[3];
sx q[3];
rz(0.010820476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.058379563) q[0];
sx q[0];
rz(-2.932974) q[0];
sx q[0];
rz(-1.5104729) q[0];
rz(2.0026813) q[1];
sx q[1];
rz(-1.4803671) q[1];
sx q[1];
rz(-1.0169792) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74422979) q[0];
sx q[0];
rz(-2.0332391) q[0];
sx q[0];
rz(0.37325333) q[0];
x q[1];
rz(1.5469656) q[2];
sx q[2];
rz(-1.7719291) q[2];
sx q[2];
rz(0.80326524) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0000227) q[1];
sx q[1];
rz(-2.0380028) q[1];
sx q[1];
rz(2.1187374) q[1];
rz(-pi) q[2];
rz(0.75628144) q[3];
sx q[3];
rz(-2.45621) q[3];
sx q[3];
rz(-1.7862198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5974474) q[2];
sx q[2];
rz(-2.1400673) q[2];
sx q[2];
rz(-2.5069359) q[2];
rz(-1.0774111) q[3];
sx q[3];
rz(-1.0297188) q[3];
sx q[3];
rz(-0.44657782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8909661) q[0];
sx q[0];
rz(-0.64193305) q[0];
sx q[0];
rz(0.84386688) q[0];
rz(-1.5232874) q[1];
sx q[1];
rz(-2.4089676) q[1];
sx q[1];
rz(-1.9218146) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6424113) q[0];
sx q[0];
rz(-1.6089604) q[0];
sx q[0];
rz(3.0495318) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.90791038) q[2];
sx q[2];
rz(-0.50953509) q[2];
sx q[2];
rz(2.8223035) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.68901686) q[1];
sx q[1];
rz(-2.2170482) q[1];
sx q[1];
rz(-0.31022443) q[1];
rz(-1.1620429) q[3];
sx q[3];
rz(-1.4882653) q[3];
sx q[3];
rz(-2.2857034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5119778) q[2];
sx q[2];
rz(-2.2045366) q[2];
sx q[2];
rz(-0.31663319) q[2];
rz(0.037242446) q[3];
sx q[3];
rz(-1.7946295) q[3];
sx q[3];
rz(0.75604701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0379631) q[0];
sx q[0];
rz(-2.0541971) q[0];
sx q[0];
rz(-0.33139247) q[0];
rz(-1.9695075) q[1];
sx q[1];
rz(-0.66292271) q[1];
sx q[1];
rz(0.40922871) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5361355) q[0];
sx q[0];
rz(-2.1989609) q[0];
sx q[0];
rz(0.34988846) q[0];
rz(-pi) q[1];
rz(2.6312469) q[2];
sx q[2];
rz(-2.5556457) q[2];
sx q[2];
rz(-0.79615359) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1815391) q[1];
sx q[1];
rz(-2.1337193) q[1];
sx q[1];
rz(-1.7510508) q[1];
x q[2];
rz(-0.41142923) q[3];
sx q[3];
rz(-0.55555389) q[3];
sx q[3];
rz(-3.0761449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1428712) q[2];
sx q[2];
rz(-1.2850392) q[2];
sx q[2];
rz(-2.5869353) q[2];
rz(-2.5000642) q[3];
sx q[3];
rz(-2.9252164) q[3];
sx q[3];
rz(-0.085111246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91759578) q[0];
sx q[0];
rz(-1.5107091) q[0];
sx q[0];
rz(3.1316485) q[0];
rz(1.0558646) q[1];
sx q[1];
rz(-1.5664342) q[1];
sx q[1];
rz(1.508629) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2752339) q[0];
sx q[0];
rz(-0.70952053) q[0];
sx q[0];
rz(1.8190246) q[0];
rz(-pi) q[1];
rz(1.4417068) q[2];
sx q[2];
rz(-2.522905) q[2];
sx q[2];
rz(2.8512851) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.29533169) q[1];
sx q[1];
rz(-1.2519072) q[1];
sx q[1];
rz(-1.0237414) q[1];
x q[2];
rz(-1.6375332) q[3];
sx q[3];
rz(-0.70422322) q[3];
sx q[3];
rz(3.0610386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.14279723) q[2];
sx q[2];
rz(-1.5441511) q[2];
sx q[2];
rz(-2.7049086) q[2];
rz(-1.8113332) q[3];
sx q[3];
rz(-2.4618849) q[3];
sx q[3];
rz(-1.0288303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-3.0325539) q[0];
sx q[0];
rz(-2.3333896) q[0];
sx q[0];
rz(-0.56525266) q[0];
rz(0.25282192) q[1];
sx q[1];
rz(-2.1222474) q[1];
sx q[1];
rz(-1.3814829) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0488659) q[0];
sx q[0];
rz(-0.5038358) q[0];
sx q[0];
rz(-0.086765246) q[0];
x q[1];
rz(1.6818468) q[2];
sx q[2];
rz(-1.0990708) q[2];
sx q[2];
rz(1.2088838) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6773721) q[1];
sx q[1];
rz(-0.83253011) q[1];
sx q[1];
rz(-2.1470451) q[1];
rz(-pi) q[2];
rz(-1.2425209) q[3];
sx q[3];
rz(-2.4007113) q[3];
sx q[3];
rz(1.5433943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.4505724) q[2];
sx q[2];
rz(-2.317163) q[2];
sx q[2];
rz(0.59664574) q[2];
rz(0.48554844) q[3];
sx q[3];
rz(-2.2462626) q[3];
sx q[3];
rz(-0.027464494) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52994603) q[0];
sx q[0];
rz(-1.2048789) q[0];
sx q[0];
rz(-2.0299029) q[0];
rz(-1.272841) q[1];
sx q[1];
rz(-2.0618989) q[1];
sx q[1];
rz(2.3088574) q[1];
rz(-1.8252774) q[2];
sx q[2];
rz(-1.5263867) q[2];
sx q[2];
rz(2.9771794) q[2];
rz(-2.6315401) q[3];
sx q[3];
rz(-1.4997713) q[3];
sx q[3];
rz(-0.98239519) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
