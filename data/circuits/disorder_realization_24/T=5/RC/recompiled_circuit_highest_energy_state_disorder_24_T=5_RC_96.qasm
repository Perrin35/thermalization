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
rz(0.46718207) q[0];
sx q[0];
rz(-1.2527569) q[0];
sx q[0];
rz(-2.3350265) q[0];
rz(2.2944577) q[1];
sx q[1];
rz(-2.6557014) q[1];
sx q[1];
rz(2.2810305) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0057324) q[0];
sx q[0];
rz(-1.4306698) q[0];
sx q[0];
rz(2.1707373) q[0];
rz(-pi) q[1];
rz(1.8612618) q[2];
sx q[2];
rz(-2.4142401) q[2];
sx q[2];
rz(2.1216105) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.367975) q[1];
sx q[1];
rz(-0.27918511) q[1];
sx q[1];
rz(0.18191819) q[1];
rz(-pi) q[2];
x q[2];
rz(0.99765649) q[3];
sx q[3];
rz(-1.1449779) q[3];
sx q[3];
rz(-2.4906858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2412771) q[2];
sx q[2];
rz(-0.52854717) q[2];
sx q[2];
rz(-0.47145525) q[2];
rz(0.49247646) q[3];
sx q[3];
rz(-1.9086647) q[3];
sx q[3];
rz(-1.2537778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8660698) q[0];
sx q[0];
rz(-0.38551426) q[0];
sx q[0];
rz(-2.8634014) q[0];
rz(1.9506075) q[1];
sx q[1];
rz(-1.8089801) q[1];
sx q[1];
rz(-1.2954378) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14034056) q[0];
sx q[0];
rz(-0.80632001) q[0];
sx q[0];
rz(-0.6722404) q[0];
rz(-pi) q[1];
rz(2.5668199) q[2];
sx q[2];
rz(-3.0468371) q[2];
sx q[2];
rz(3.0916328) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5767314) q[1];
sx q[1];
rz(-2.7612491) q[1];
sx q[1];
rz(-1.5647792) q[1];
rz(-pi) q[2];
rz(-2.6023618) q[3];
sx q[3];
rz(-2.1279716) q[3];
sx q[3];
rz(2.2702366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.94464716) q[2];
sx q[2];
rz(-1.4713919) q[2];
sx q[2];
rz(-0.50986457) q[2];
rz(1.6950131) q[3];
sx q[3];
rz(-2.6810985) q[3];
sx q[3];
rz(-0.50486008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6079717) q[0];
sx q[0];
rz(-1.6225659) q[0];
sx q[0];
rz(-2.1571889) q[0];
rz(3.1254752) q[1];
sx q[1];
rz(-2.0033483) q[1];
sx q[1];
rz(1.791753) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70671073) q[0];
sx q[0];
rz(-2.2684386) q[0];
sx q[0];
rz(-1.7876704) q[0];
rz(2.0651814) q[2];
sx q[2];
rz(-2.5177023) q[2];
sx q[2];
rz(-0.15984331) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4011154) q[1];
sx q[1];
rz(-0.7001895) q[1];
sx q[1];
rz(-0.062559532) q[1];
rz(-pi) q[2];
rz(-0.11589072) q[3];
sx q[3];
rz(-2.3724764) q[3];
sx q[3];
rz(2.9562841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4237889) q[2];
sx q[2];
rz(-2.4393647) q[2];
sx q[2];
rz(0.88199893) q[2];
rz(-1.9574022) q[3];
sx q[3];
rz(-0.94112527) q[3];
sx q[3];
rz(0.73630303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0940359) q[0];
sx q[0];
rz(-1.9147669) q[0];
sx q[0];
rz(3.0645698) q[0];
rz(-2.9432964) q[1];
sx q[1];
rz(-1.6402596) q[1];
sx q[1];
rz(-0.1870627) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2614005) q[0];
sx q[0];
rz(-2.7679256) q[0];
sx q[0];
rz(0.4445266) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9603029) q[2];
sx q[2];
rz(-1.9167056) q[2];
sx q[2];
rz(0.51900253) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6641561) q[1];
sx q[1];
rz(-2.2786447) q[1];
sx q[1];
rz(-1.493046) q[1];
rz(-pi) q[2];
rz(-1.2310394) q[3];
sx q[3];
rz(-1.8700784) q[3];
sx q[3];
rz(-1.9055942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5483115) q[2];
sx q[2];
rz(-2.686794) q[2];
sx q[2];
rz(-1.456267) q[2];
rz(0.71825394) q[3];
sx q[3];
rz(-1.3879489) q[3];
sx q[3];
rz(-0.83435241) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2497571) q[0];
sx q[0];
rz(-1.3351853) q[0];
sx q[0];
rz(0.43858132) q[0];
rz(-2.7843685) q[1];
sx q[1];
rz(-1.8635187) q[1];
sx q[1];
rz(-1.8908148) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8698602) q[0];
sx q[0];
rz(-1.1446867) q[0];
sx q[0];
rz(-0.10978384) q[0];
x q[1];
rz(1.3322796) q[2];
sx q[2];
rz(-0.3198238) q[2];
sx q[2];
rz(-0.25898283) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.631613) q[1];
sx q[1];
rz(-2.5812979) q[1];
sx q[1];
rz(2.3255682) q[1];
rz(-0.26024466) q[3];
sx q[3];
rz(-1.6992927) q[3];
sx q[3];
rz(-2.5082626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.232051) q[2];
sx q[2];
rz(-0.48018685) q[2];
sx q[2];
rz(1.5117744) q[2];
rz(-0.016228598) q[3];
sx q[3];
rz(-1.0954233) q[3];
sx q[3];
rz(-0.79286638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4474354) q[0];
sx q[0];
rz(-1.2801535) q[0];
sx q[0];
rz(1.9919027) q[0];
rz(2.7206874) q[1];
sx q[1];
rz(-1.4525388) q[1];
sx q[1];
rz(-0.044205753) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4251148) q[0];
sx q[0];
rz(-1.8938066) q[0];
sx q[0];
rz(-2.4130505) q[0];
rz(2.1109796) q[2];
sx q[2];
rz(-1.9760796) q[2];
sx q[2];
rz(-1.1073529) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.845712) q[1];
sx q[1];
rz(-0.39870495) q[1];
sx q[1];
rz(-1.7804342) q[1];
x q[2];
rz(1.2941235) q[3];
sx q[3];
rz(-1.0545316) q[3];
sx q[3];
rz(2.6003305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9629024) q[2];
sx q[2];
rz(-2.2447605) q[2];
sx q[2];
rz(0.94513354) q[2];
rz(-1.6804228) q[3];
sx q[3];
rz(-1.9943359) q[3];
sx q[3];
rz(-2.6233961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1951676) q[0];
sx q[0];
rz(-2.8312046) q[0];
sx q[0];
rz(-2.2921966) q[0];
rz(1.3769582) q[1];
sx q[1];
rz(-0.57146776) q[1];
sx q[1];
rz(-1.2845385) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.021668108) q[0];
sx q[0];
rz(-0.64739043) q[0];
sx q[0];
rz(-0.68418087) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6780186) q[2];
sx q[2];
rz(-1.2649415) q[2];
sx q[2];
rz(2.0448409) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4139288) q[1];
sx q[1];
rz(-1.5352121) q[1];
sx q[1];
rz(0.52558454) q[1];
rz(-2.128781) q[3];
sx q[3];
rz(-1.9922755) q[3];
sx q[3];
rz(-2.1881605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.66594243) q[2];
sx q[2];
rz(-1.3694171) q[2];
sx q[2];
rz(1.5744038) q[2];
rz(0.33268467) q[3];
sx q[3];
rz(-2.0761469) q[3];
sx q[3];
rz(0.98193297) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57629267) q[0];
sx q[0];
rz(-1.1197634) q[0];
sx q[0];
rz(2.0530307) q[0];
rz(2.2125878) q[1];
sx q[1];
rz(-1.4938415) q[1];
sx q[1];
rz(0.30379024) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.740363) q[0];
sx q[0];
rz(-1.0107733) q[0];
sx q[0];
rz(2.014888) q[0];
rz(-pi) q[1];
rz(2.1985198) q[2];
sx q[2];
rz(-0.20096261) q[2];
sx q[2];
rz(-0.64370868) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.15081295) q[1];
sx q[1];
rz(-2.7013198) q[1];
sx q[1];
rz(1.7899465) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.522842) q[3];
sx q[3];
rz(-1.1693952) q[3];
sx q[3];
rz(-1.5681745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2527689) q[2];
sx q[2];
rz(-1.4759651) q[2];
sx q[2];
rz(-0.26407537) q[2];
rz(-1.1325599) q[3];
sx q[3];
rz(-0.3370291) q[3];
sx q[3];
rz(-2.0866709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3573414) q[0];
sx q[0];
rz(-2.7989474) q[0];
sx q[0];
rz(-1.0145048) q[0];
rz(2.2134589) q[1];
sx q[1];
rz(-1.721563) q[1];
sx q[1];
rz(0.49044213) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3630842) q[0];
sx q[0];
rz(-0.16101232) q[0];
sx q[0];
rz(-1.9912316) q[0];
x q[1];
rz(-0.84607203) q[2];
sx q[2];
rz(-0.49642379) q[2];
sx q[2];
rz(-1.7299394) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.38476506) q[1];
sx q[1];
rz(-2.8382128) q[1];
sx q[1];
rz(-1.7089173) q[1];
rz(-pi) q[2];
rz(-0.2410335) q[3];
sx q[3];
rz(-2.0969982) q[3];
sx q[3];
rz(-2.0693501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.32533112) q[2];
sx q[2];
rz(-1.0518495) q[2];
sx q[2];
rz(0.10247792) q[2];
rz(1.8898194) q[3];
sx q[3];
rz(-1.7779558) q[3];
sx q[3];
rz(1.4832835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74806279) q[0];
sx q[0];
rz(-3.0952125) q[0];
sx q[0];
rz(1.8512132) q[0];
rz(-0.45404592) q[1];
sx q[1];
rz(-1.8171277) q[1];
sx q[1];
rz(-0.18347278) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2170524) q[0];
sx q[0];
rz(-1.3879225) q[0];
sx q[0];
rz(0.019269983) q[0];
x q[1];
rz(-1.5857592) q[2];
sx q[2];
rz(-0.97285473) q[2];
sx q[2];
rz(0.19549616) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3349592) q[1];
sx q[1];
rz(-1.0786966) q[1];
sx q[1];
rz(-1.5978807) q[1];
rz(-pi) q[2];
rz(0.091138361) q[3];
sx q[3];
rz(-2.1347858) q[3];
sx q[3];
rz(-2.6567949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4219249) q[2];
sx q[2];
rz(-2.0823961) q[2];
sx q[2];
rz(-0.025040778) q[2];
rz(2.5981564) q[3];
sx q[3];
rz(-0.70356026) q[3];
sx q[3];
rz(0.27633015) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80171361) q[0];
sx q[0];
rz(-1.0152974) q[0];
sx q[0];
rz(0.80831084) q[0];
rz(0.33357757) q[1];
sx q[1];
rz(-1.7671276) q[1];
sx q[1];
rz(-0.50552013) q[1];
rz(-1.2110151) q[2];
sx q[2];
rz(-2.790749) q[2];
sx q[2];
rz(-0.1505827) q[2];
rz(0.80584851) q[3];
sx q[3];
rz(-2.3060006) q[3];
sx q[3];
rz(-1.5455442) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
