OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.55943638) q[0];
sx q[0];
rz(-2.5506033) q[0];
sx q[0];
rz(-0.58340573) q[0];
rz(2.9572339) q[1];
sx q[1];
rz(-0.98362041) q[1];
sx q[1];
rz(2.2489927) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8164506) q[0];
sx q[0];
rz(-2.4976375) q[0];
sx q[0];
rz(-2.08026) q[0];
rz(-pi) q[1];
rz(-1.0640261) q[2];
sx q[2];
rz(-0.96553409) q[2];
sx q[2];
rz(-1.93047) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2979558) q[1];
sx q[1];
rz(-1.8814109) q[1];
sx q[1];
rz(-0.85308869) q[1];
x q[2];
rz(-0.046269429) q[3];
sx q[3];
rz(-2.390063) q[3];
sx q[3];
rz(2.3627594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9154174) q[2];
sx q[2];
rz(-1.5390652) q[2];
sx q[2];
rz(-1.9809451) q[2];
rz(-0.21696572) q[3];
sx q[3];
rz(-0.522885) q[3];
sx q[3];
rz(-1.0552361) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0579257) q[0];
sx q[0];
rz(-2.1766429) q[0];
sx q[0];
rz(2.5657186) q[0];
rz(1.2469762) q[1];
sx q[1];
rz(-1.8449102) q[1];
sx q[1];
rz(1.974568) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.052159667) q[0];
sx q[0];
rz(-2.8773327) q[0];
sx q[0];
rz(1.7961851) q[0];
rz(-pi) q[1];
rz(-2.8950047) q[2];
sx q[2];
rz(-1.0435259) q[2];
sx q[2];
rz(1.8794683) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0868623) q[1];
sx q[1];
rz(-2.3332101) q[1];
sx q[1];
rz(0.085309172) q[1];
rz(-pi) q[2];
rz(-0.59252177) q[3];
sx q[3];
rz(-0.7624818) q[3];
sx q[3];
rz(-3.1069063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.033826753) q[2];
sx q[2];
rz(-1.1788538) q[2];
sx q[2];
rz(-0.21437422) q[2];
rz(0.073444627) q[3];
sx q[3];
rz(-2.6918604) q[3];
sx q[3];
rz(-0.28081056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4784933) q[0];
sx q[0];
rz(-0.61613023) q[0];
sx q[0];
rz(-1.2269155) q[0];
rz(0.40027174) q[1];
sx q[1];
rz(-1.2534671) q[1];
sx q[1];
rz(-2.1267557) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0815711) q[0];
sx q[0];
rz(-0.9284174) q[0];
sx q[0];
rz(-1.6409671) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0929167) q[2];
sx q[2];
rz(-0.45959696) q[2];
sx q[2];
rz(-2.2082579) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.44753578) q[1];
sx q[1];
rz(-0.60255614) q[1];
sx q[1];
rz(1.4284929) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8869927) q[3];
sx q[3];
rz(-0.71101515) q[3];
sx q[3];
rz(2.7544114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0456475) q[2];
sx q[2];
rz(-1.0568876) q[2];
sx q[2];
rz(0.8992368) q[2];
rz(-0.69747654) q[3];
sx q[3];
rz(-1.8557502) q[3];
sx q[3];
rz(0.60788679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4863481) q[0];
sx q[0];
rz(-2.0589893) q[0];
sx q[0];
rz(-2.7096601) q[0];
rz(-2.5090384) q[1];
sx q[1];
rz(-0.41699854) q[1];
sx q[1];
rz(0.63582173) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83946562) q[0];
sx q[0];
rz(-1.3913904) q[0];
sx q[0];
rz(-2.0749712) q[0];
rz(-pi) q[1];
rz(-0.16480883) q[2];
sx q[2];
rz(-2.2909082) q[2];
sx q[2];
rz(-2.1317496) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3969877) q[1];
sx q[1];
rz(-2.3481391) q[1];
sx q[1];
rz(-0.26054392) q[1];
rz(-pi) q[2];
rz(-0.28663978) q[3];
sx q[3];
rz(-1.3661824) q[3];
sx q[3];
rz(0.74433792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9277966) q[2];
sx q[2];
rz(-0.14363229) q[2];
sx q[2];
rz(-0.68871838) q[2];
rz(-0.33411807) q[3];
sx q[3];
rz(-1.9709316) q[3];
sx q[3];
rz(-2.9978602) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9976945) q[0];
sx q[0];
rz(-0.69902885) q[0];
sx q[0];
rz(-1.0744263) q[0];
rz(-0.74514666) q[1];
sx q[1];
rz(-1.6586168) q[1];
sx q[1];
rz(2.863046) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2715831) q[0];
sx q[0];
rz(-2.2402813) q[0];
sx q[0];
rz(-0.55218009) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9859221) q[2];
sx q[2];
rz(-0.61741932) q[2];
sx q[2];
rz(1.9215259) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7601732) q[1];
sx q[1];
rz(-2.1790494) q[1];
sx q[1];
rz(-2.6890523) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1513125) q[3];
sx q[3];
rz(-2.2970082) q[3];
sx q[3];
rz(2.5648404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6986065) q[2];
sx q[2];
rz(-1.3977945) q[2];
sx q[2];
rz(-2.8473575) q[2];
rz(-3.0596628) q[3];
sx q[3];
rz(-0.51920813) q[3];
sx q[3];
rz(0.036227139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0329523) q[0];
sx q[0];
rz(-0.8529129) q[0];
sx q[0];
rz(-0.0090573514) q[0];
rz(-0.63502216) q[1];
sx q[1];
rz(-2.451684) q[1];
sx q[1];
rz(-3.0335398) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8979643) q[0];
sx q[0];
rz(-2.0192696) q[0];
sx q[0];
rz(-1.2178221) q[0];
rz(-pi) q[1];
rz(-1.4322386) q[2];
sx q[2];
rz(-1.1922622) q[2];
sx q[2];
rz(-2.0356503) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9891226) q[1];
sx q[1];
rz(-1.9580541) q[1];
sx q[1];
rz(0.5247922) q[1];
x q[2];
rz(-1.024527) q[3];
sx q[3];
rz(-1.3427991) q[3];
sx q[3];
rz(0.25817623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.99284995) q[2];
sx q[2];
rz(-1.7603346) q[2];
sx q[2];
rz(1.8704869) q[2];
rz(3.0631915) q[3];
sx q[3];
rz(-1.6631118) q[3];
sx q[3];
rz(-1.1318077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-0.25061297) q[0];
sx q[0];
rz(-1.5761292) q[0];
sx q[0];
rz(-0.65761956) q[0];
rz(-1.3972067) q[1];
sx q[1];
rz(-1.1746635) q[1];
sx q[1];
rz(-2.2479642) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0890869) q[0];
sx q[0];
rz(-2.1930709) q[0];
sx q[0];
rz(1.9786579) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.29427476) q[2];
sx q[2];
rz(-2.0308959) q[2];
sx q[2];
rz(2.0725046) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0221755) q[1];
sx q[1];
rz(-0.70964538) q[1];
sx q[1];
rz(2.8303353) q[1];
rz(-pi) q[2];
rz(0.10223933) q[3];
sx q[3];
rz(-0.7623626) q[3];
sx q[3];
rz(-1.8765212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7509193) q[2];
sx q[2];
rz(-2.1942287) q[2];
sx q[2];
rz(-2.612109) q[2];
rz(2.6654065) q[3];
sx q[3];
rz(-1.3323077) q[3];
sx q[3];
rz(2.7990394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7628409) q[0];
sx q[0];
rz(-1.5126001) q[0];
sx q[0];
rz(2.0513127) q[0];
rz(-0.11225637) q[1];
sx q[1];
rz(-2.0376164) q[1];
sx q[1];
rz(-1.1539248) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8114416) q[0];
sx q[0];
rz(-0.92130565) q[0];
sx q[0];
rz(1.5280456) q[0];
rz(0.23837337) q[2];
sx q[2];
rz(-2.5775238) q[2];
sx q[2];
rz(-0.43018815) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6653319) q[1];
sx q[1];
rz(-1.5090824) q[1];
sx q[1];
rz(-0.76121059) q[1];
rz(-pi) q[2];
x q[2];
rz(0.38903799) q[3];
sx q[3];
rz(-0.93718796) q[3];
sx q[3];
rz(0.47215677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.551679) q[2];
sx q[2];
rz(-2.7605197) q[2];
sx q[2];
rz(-2.7611458) q[2];
rz(2.0137265) q[3];
sx q[3];
rz(-1.3600072) q[3];
sx q[3];
rz(-1.5765223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34148759) q[0];
sx q[0];
rz(-1.8999758) q[0];
sx q[0];
rz(-0.63968101) q[0];
rz(-1.9027963) q[1];
sx q[1];
rz(-1.4150554) q[1];
sx q[1];
rz(1.9715086) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6696475) q[0];
sx q[0];
rz(-1.8983316) q[0];
sx q[0];
rz(1.9452842) q[0];
x q[1];
rz(-0.27530382) q[2];
sx q[2];
rz(-2.346056) q[2];
sx q[2];
rz(-0.77992935) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1994836) q[1];
sx q[1];
rz(-1.0892727) q[1];
sx q[1];
rz(1.7755309) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.13799237) q[3];
sx q[3];
rz(-0.90413168) q[3];
sx q[3];
rz(-0.43825144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3141979) q[2];
sx q[2];
rz(-0.65076995) q[2];
sx q[2];
rz(2.7588552) q[2];
rz(-0.9283723) q[3];
sx q[3];
rz(-1.9675156) q[3];
sx q[3];
rz(-2.476957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2255573) q[0];
sx q[0];
rz(-1.501843) q[0];
sx q[0];
rz(-0.25892192) q[0];
rz(-0.71031538) q[1];
sx q[1];
rz(-1.0537035) q[1];
sx q[1];
rz(-2.6616667) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5537162) q[0];
sx q[0];
rz(-1.9037316) q[0];
sx q[0];
rz(2.2399708) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9051022) q[2];
sx q[2];
rz(-2.2220526) q[2];
sx q[2];
rz(2.1361534) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.48206115) q[1];
sx q[1];
rz(-1.2188984) q[1];
sx q[1];
rz(0.54606502) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2009833) q[3];
sx q[3];
rz(-0.77946957) q[3];
sx q[3];
rz(1.1704695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2991128) q[2];
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
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15923545) q[0];
sx q[0];
rz(-3.1095105) q[0];
sx q[0];
rz(1.6788917) q[0];
rz(2.1622529) q[1];
sx q[1];
rz(-2.0420488) q[1];
sx q[1];
rz(2.2534823) q[1];
rz(0.32817763) q[2];
sx q[2];
rz(-0.50445088) q[2];
sx q[2];
rz(2.9403461) q[2];
rz(0.27771523) q[3];
sx q[3];
rz(-1.82901) q[3];
sx q[3];
rz(1.3808586) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
