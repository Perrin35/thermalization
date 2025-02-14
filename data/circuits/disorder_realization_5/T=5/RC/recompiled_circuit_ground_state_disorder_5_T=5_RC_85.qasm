OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.65052819) q[0];
sx q[0];
rz(-1.0125546) q[0];
sx q[0];
rz(0.92232409) q[0];
rz(2.2946279) q[1];
sx q[1];
rz(-1.4743409) q[1];
sx q[1];
rz(2.9234731) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0357649) q[0];
sx q[0];
rz(-2.0313861) q[0];
sx q[0];
rz(1.4683819) q[0];
x q[1];
rz(2.7635771) q[2];
sx q[2];
rz(-1.5326008) q[2];
sx q[2];
rz(-2.5390944) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.328124) q[1];
sx q[1];
rz(-1.3513997) q[1];
sx q[1];
rz(0.043299999) q[1];
x q[2];
rz(-3.0633801) q[3];
sx q[3];
rz(-1.999475) q[3];
sx q[3];
rz(-0.058905963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.073079022) q[2];
sx q[2];
rz(-1.2129236) q[2];
sx q[2];
rz(2.6050513) q[2];
rz(-2.1327298) q[3];
sx q[3];
rz(-0.77102414) q[3];
sx q[3];
rz(-2.3276276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5040078) q[0];
sx q[0];
rz(-1.3115839) q[0];
sx q[0];
rz(0.017039321) q[0];
rz(-2.0783966) q[1];
sx q[1];
rz(-0.76779643) q[1];
sx q[1];
rz(-0.95471901) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1621891) q[0];
sx q[0];
rz(-1.6260176) q[0];
sx q[0];
rz(2.9121141) q[0];
rz(1.0771712) q[2];
sx q[2];
rz(-2.8045678) q[2];
sx q[2];
rz(-0.21060196) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.20978759) q[1];
sx q[1];
rz(-1.796319) q[1];
sx q[1];
rz(2.8978845) q[1];
rz(-pi) q[2];
rz(1.6240261) q[3];
sx q[3];
rz(-1.5900776) q[3];
sx q[3];
rz(2.791491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.22975989) q[2];
sx q[2];
rz(-2.2288897) q[2];
sx q[2];
rz(-2.0274053) q[2];
rz(-2.0071425) q[3];
sx q[3];
rz(-2.9454234) q[3];
sx q[3];
rz(0.1964868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.590362) q[0];
sx q[0];
rz(-2.1901665) q[0];
sx q[0];
rz(0.1828585) q[0];
rz(-3.0602449) q[1];
sx q[1];
rz(-0.75129879) q[1];
sx q[1];
rz(2.523211) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8326022) q[0];
sx q[0];
rz(-1.497735) q[0];
sx q[0];
rz(-2.354524) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7612533) q[2];
sx q[2];
rz(-2.3229685) q[2];
sx q[2];
rz(-0.10701767) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.32633852) q[1];
sx q[1];
rz(-2.7732121) q[1];
sx q[1];
rz(-1.1756363) q[1];
x q[2];
rz(0.86520536) q[3];
sx q[3];
rz(-0.83704797) q[3];
sx q[3];
rz(-0.037038304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1387834) q[2];
sx q[2];
rz(-1.323779) q[2];
sx q[2];
rz(-0.36661026) q[2];
rz(1.9256516) q[3];
sx q[3];
rz(-1.1953019) q[3];
sx q[3];
rz(0.6634357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65130305) q[0];
sx q[0];
rz(-1.9090575) q[0];
sx q[0];
rz(2.41462) q[0];
rz(-2.2333249) q[1];
sx q[1];
rz(-2.5636702) q[1];
sx q[1];
rz(1.04331) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3573049) q[0];
sx q[0];
rz(-2.4803964) q[0];
sx q[0];
rz(-3.1327598) q[0];
rz(-pi) q[1];
rz(-1.7264446) q[2];
sx q[2];
rz(-0.37964941) q[2];
sx q[2];
rz(-1.2641035) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.88731447) q[1];
sx q[1];
rz(-1.320854) q[1];
sx q[1];
rz(-2.0292086) q[1];
rz(0.38049145) q[3];
sx q[3];
rz(-1.0416608) q[3];
sx q[3];
rz(-0.38705119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4116481) q[2];
sx q[2];
rz(-2.8204155) q[2];
sx q[2];
rz(1.3058861) q[2];
rz(-3.0236687) q[3];
sx q[3];
rz(-2.6730461) q[3];
sx q[3];
rz(-2.0809295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1082728) q[0];
sx q[0];
rz(-0.98618996) q[0];
sx q[0];
rz(1.5340075) q[0];
rz(-0.29574212) q[1];
sx q[1];
rz(-1.60138) q[1];
sx q[1];
rz(-1.4588413) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5285501) q[0];
sx q[0];
rz(-0.82697059) q[0];
sx q[0];
rz(-3.0463329) q[0];
x q[1];
rz(2.2590911) q[2];
sx q[2];
rz(-0.46445981) q[2];
sx q[2];
rz(-1.2683753) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4319358) q[1];
sx q[1];
rz(-1.767358) q[1];
sx q[1];
rz(0.86039575) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9146054) q[3];
sx q[3];
rz(-0.60887486) q[3];
sx q[3];
rz(0.46907779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.38833388) q[2];
sx q[2];
rz(-2.7497141) q[2];
sx q[2];
rz(-2.4957116) q[2];
rz(-1.215747) q[3];
sx q[3];
rz(-1.682621) q[3];
sx q[3];
rz(-1.7997883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.011768613) q[0];
sx q[0];
rz(-0.83368603) q[0];
sx q[0];
rz(-0.60761333) q[0];
rz(-1.1786849) q[1];
sx q[1];
rz(-1.2858398) q[1];
sx q[1];
rz(-0.36373055) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54485142) q[0];
sx q[0];
rz(-0.57692301) q[0];
sx q[0];
rz(1.1240918) q[0];
rz(-pi) q[1];
rz(1.9035788) q[2];
sx q[2];
rz(-1.6146891) q[2];
sx q[2];
rz(-1.8157168) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1235001) q[1];
sx q[1];
rz(-2.315633) q[1];
sx q[1];
rz(2.9781746) q[1];
x q[2];
rz(2.167114) q[3];
sx q[3];
rz(-1.5502366) q[3];
sx q[3];
rz(-2.1316776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5522449) q[2];
sx q[2];
rz(-3.037368) q[2];
sx q[2];
rz(-1.1927401) q[2];
rz(-0.73181152) q[3];
sx q[3];
rz(-1.4251499) q[3];
sx q[3];
rz(1.4881136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.8530387) q[0];
sx q[0];
rz(-0.1703425) q[0];
sx q[0];
rz(-1.5227675) q[0];
rz(1.6743926) q[1];
sx q[1];
rz(-1.3102945) q[1];
sx q[1];
rz(2.4640962) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.753669) q[0];
sx q[0];
rz(-0.69783995) q[0];
sx q[0];
rz(-1.4927255) q[0];
rz(-pi) q[1];
rz(-1.7588781) q[2];
sx q[2];
rz(-1.3408337) q[2];
sx q[2];
rz(2.8158958) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5657357) q[1];
sx q[1];
rz(-1.9565554) q[1];
sx q[1];
rz(-1.1123453) q[1];
x q[2];
rz(1.7301637) q[3];
sx q[3];
rz(-2.2969118) q[3];
sx q[3];
rz(2.1350525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.39685321) q[2];
sx q[2];
rz(-2.2262959) q[2];
sx q[2];
rz(1.8358561) q[2];
rz(-0.33852494) q[3];
sx q[3];
rz(-1.1208231) q[3];
sx q[3];
rz(0.6849851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5386388) q[0];
sx q[0];
rz(-2.9260577) q[0];
sx q[0];
rz(2.9576874) q[0];
rz(-1.454608) q[1];
sx q[1];
rz(-0.55211663) q[1];
sx q[1];
rz(-0.49585453) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3199098) q[0];
sx q[0];
rz(-1.3047072) q[0];
sx q[0];
rz(2.0771107) q[0];
x q[1];
rz(0.77347254) q[2];
sx q[2];
rz(-2.2272553) q[2];
sx q[2];
rz(0.17420775) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7396122) q[1];
sx q[1];
rz(-2.6079) q[1];
sx q[1];
rz(0.81783612) q[1];
rz(-1.0184278) q[3];
sx q[3];
rz(-1.6721069) q[3];
sx q[3];
rz(-0.61063572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0869861) q[2];
sx q[2];
rz(-0.84422529) q[2];
sx q[2];
rz(-1.1620713) q[2];
rz(1.453513) q[3];
sx q[3];
rz(-2.1789357) q[3];
sx q[3];
rz(-1.2801722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7206955) q[0];
sx q[0];
rz(-1.9891885) q[0];
sx q[0];
rz(0.56588093) q[0];
rz(-0.3903009) q[1];
sx q[1];
rz(-1.5696328) q[1];
sx q[1];
rz(-1.3077259) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54057656) q[0];
sx q[0];
rz(-1.8800288) q[0];
sx q[0];
rz(2.087167) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1928608) q[2];
sx q[2];
rz(-1.7687174) q[2];
sx q[2];
rz(1.4858311) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.81586397) q[1];
sx q[1];
rz(-0.96885175) q[1];
sx q[1];
rz(1.25156) q[1];
rz(-0.89687895) q[3];
sx q[3];
rz(-0.19052902) q[3];
sx q[3];
rz(0.49303699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2719443) q[2];
sx q[2];
rz(-0.41676909) q[2];
sx q[2];
rz(-2.1040253) q[2];
rz(-1.3197673) q[3];
sx q[3];
rz(-2.3128553) q[3];
sx q[3];
rz(-2.493609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2607525) q[0];
sx q[0];
rz(-2.1639731) q[0];
sx q[0];
rz(-2.4993437) q[0];
rz(-1.3308659) q[1];
sx q[1];
rz(-2.2763177) q[1];
sx q[1];
rz(2.9790402) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35643043) q[0];
sx q[0];
rz(-0.9097865) q[0];
sx q[0];
rz(-0.50161538) q[0];
rz(-pi) q[1];
rz(1.2106718) q[2];
sx q[2];
rz(-1.2136493) q[2];
sx q[2];
rz(1.2993297) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4288669) q[1];
sx q[1];
rz(-1.7460632) q[1];
sx q[1];
rz(-1.0840985) q[1];
rz(-pi) q[2];
rz(-0.58329669) q[3];
sx q[3];
rz(-1.5044893) q[3];
sx q[3];
rz(-0.60140677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6315397) q[2];
sx q[2];
rz(-1.8659464) q[2];
sx q[2];
rz(2.1007382) q[2];
rz(-0.96796525) q[3];
sx q[3];
rz(-2.0716045) q[3];
sx q[3];
rz(-0.66106838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80617245) q[0];
sx q[0];
rz(-1.5652884) q[0];
sx q[0];
rz(-0.096927222) q[0];
rz(2.9087635) q[1];
sx q[1];
rz(-2.1131344) q[1];
sx q[1];
rz(0.0075385787) q[1];
rz(2.3940968) q[2];
sx q[2];
rz(-1.9608254) q[2];
sx q[2];
rz(-1.335782) q[2];
rz(1.9540174) q[3];
sx q[3];
rz(-0.90771994) q[3];
sx q[3];
rz(-0.049605443) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
