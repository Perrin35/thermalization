OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.4133889) q[0];
sx q[0];
rz(-1.1336741) q[0];
sx q[0];
rz(1.5925621) q[0];
rz(1.6917317) q[1];
sx q[1];
rz(5.6258968) q[1];
sx q[1];
rz(13.110553) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.005851) q[0];
sx q[0];
rz(-0.07677456) q[0];
sx q[0];
rz(-2.6430921) q[0];
rz(-pi) q[1];
rz(1.5266225) q[2];
sx q[2];
rz(-1.6720811) q[2];
sx q[2];
rz(-3.0169808) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.7218329) q[1];
sx q[1];
rz(-2.3082323) q[1];
sx q[1];
rz(-2.1268197) q[1];
rz(-pi) q[2];
rz(-0.55239001) q[3];
sx q[3];
rz(-0.033621764) q[3];
sx q[3];
rz(0.53510964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.071775285) q[2];
sx q[2];
rz(-1.8775619) q[2];
sx q[2];
rz(1.7791746) q[2];
rz(0.028256265) q[3];
sx q[3];
rz(-1.3794206) q[3];
sx q[3];
rz(2.3513667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
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
rz(2.4966999) q[0];
sx q[0];
rz(-2.4364478) q[0];
sx q[0];
rz(-1.9702966) q[0];
rz(0.21121875) q[1];
sx q[1];
rz(-0.44208458) q[1];
sx q[1];
rz(1.404095) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2962759) q[0];
sx q[0];
rz(-0.8527841) q[0];
sx q[0];
rz(-0.57384558) q[0];
rz(2.8318769) q[2];
sx q[2];
rz(-2.6385912) q[2];
sx q[2];
rz(2.3878857) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5215056) q[1];
sx q[1];
rz(-1.5477991) q[1];
sx q[1];
rz(2.8584245) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8849765) q[3];
sx q[3];
rz(-0.49391541) q[3];
sx q[3];
rz(-1.4790505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.30963787) q[2];
sx q[2];
rz(-1.4761304) q[2];
sx q[2];
rz(2.1976166) q[2];
rz(-0.55654636) q[3];
sx q[3];
rz(-0.49749938) q[3];
sx q[3];
rz(1.336162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-2.8495162) q[0];
sx q[0];
rz(-2.3777666) q[0];
sx q[0];
rz(1.4235494) q[0];
rz(-0.81958333) q[1];
sx q[1];
rz(-2.3312566) q[1];
sx q[1];
rz(-2.5779285) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99479988) q[0];
sx q[0];
rz(-2.1437862) q[0];
sx q[0];
rz(-1.4865727) q[0];
rz(-0.59994772) q[2];
sx q[2];
rz(-1.7581345) q[2];
sx q[2];
rz(-1.3647321) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.015805294) q[1];
sx q[1];
rz(-1.2450706) q[1];
sx q[1];
rz(2.159352) q[1];
rz(-1.9824969) q[3];
sx q[3];
rz(-2.1512096) q[3];
sx q[3];
rz(-1.2702277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.50239572) q[2];
sx q[2];
rz(-2.0194619) q[2];
sx q[2];
rz(1.0602661) q[2];
rz(-0.075573102) q[3];
sx q[3];
rz(-0.85791701) q[3];
sx q[3];
rz(-0.11463541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8183427) q[0];
sx q[0];
rz(-2.8338354) q[0];
sx q[0];
rz(0.19317214) q[0];
rz(1.5974143) q[1];
sx q[1];
rz(-1.1491821) q[1];
sx q[1];
rz(-0.65778041) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0021129) q[0];
sx q[0];
rz(-1.6008899) q[0];
sx q[0];
rz(-2.9667903) q[0];
rz(-pi) q[1];
x q[1];
rz(0.31099702) q[2];
sx q[2];
rz(-1.1592835) q[2];
sx q[2];
rz(2.3222773) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0187877) q[1];
sx q[1];
rz(-2.0969166) q[1];
sx q[1];
rz(0.29463525) q[1];
x q[2];
rz(-2.8834881) q[3];
sx q[3];
rz(-1.2873642) q[3];
sx q[3];
rz(2.5131445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0458935) q[2];
sx q[2];
rz(-2.7313488) q[2];
sx q[2];
rz(-2.0945385) q[2];
rz(1.0632769) q[3];
sx q[3];
rz(-1.1626817) q[3];
sx q[3];
rz(-1.8803546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3445774) q[0];
sx q[0];
rz(-0.319096) q[0];
sx q[0];
rz(-1.0239333) q[0];
rz(-0.78760415) q[1];
sx q[1];
rz(-1.5658295) q[1];
sx q[1];
rz(-0.62686282) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89676566) q[0];
sx q[0];
rz(-1.5875495) q[0];
sx q[0];
rz(-1.6003952) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2704029) q[2];
sx q[2];
rz(-1.2671748) q[2];
sx q[2];
rz(0.52263573) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6080731) q[1];
sx q[1];
rz(-0.43577172) q[1];
sx q[1];
rz(-0.79969745) q[1];
rz(-pi) q[2];
rz(1.2792148) q[3];
sx q[3];
rz(-2.7334308) q[3];
sx q[3];
rz(-2.6274031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2287067) q[2];
sx q[2];
rz(-1.9084946) q[2];
sx q[2];
rz(2.1234925) q[2];
rz(-0.61156887) q[3];
sx q[3];
rz(-1.3221778) q[3];
sx q[3];
rz(-0.95156804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.6927032) q[0];
sx q[0];
rz(-1.792181) q[0];
sx q[0];
rz(-1.1992136) q[0];
rz(-1.1692283) q[1];
sx q[1];
rz(-1.0511845) q[1];
sx q[1];
rz(0.32454023) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6231461) q[0];
sx q[0];
rz(-1.070797) q[0];
sx q[0];
rz(-0.48796939) q[0];
rz(-pi) q[1];
rz(0.80771031) q[2];
sx q[2];
rz(-1.0683904) q[2];
sx q[2];
rz(-3.0903357) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4536344) q[1];
sx q[1];
rz(-1.0369685) q[1];
sx q[1];
rz(0.028793528) q[1];
rz(1.2732029) q[3];
sx q[3];
rz(-1.1614359) q[3];
sx q[3];
rz(0.62664947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.67363182) q[2];
sx q[2];
rz(-1.8362074) q[2];
sx q[2];
rz(-1.8661873) q[2];
rz(2.0868789) q[3];
sx q[3];
rz(-2.349699) q[3];
sx q[3];
rz(-1.595343) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6773029) q[0];
sx q[0];
rz(-1.3339366) q[0];
sx q[0];
rz(-1.7328847) q[0];
rz(-2.6858221) q[1];
sx q[1];
rz(-2.9401638) q[1];
sx q[1];
rz(1.2021525) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81284886) q[0];
sx q[0];
rz(-2.3658731) q[0];
sx q[0];
rz(-2.1562955) q[0];
rz(-pi) q[1];
rz(2.7289594) q[2];
sx q[2];
rz(-0.79421439) q[2];
sx q[2];
rz(2.7318294) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.251293) q[1];
sx q[1];
rz(-1.6541462) q[1];
sx q[1];
rz(-3.0492196) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.60336242) q[3];
sx q[3];
rz(-2.8238736) q[3];
sx q[3];
rz(0.41894223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1414286) q[2];
sx q[2];
rz(-1.2951853) q[2];
sx q[2];
rz(-0.84623519) q[2];
rz(-2.6464461) q[3];
sx q[3];
rz(-0.64619243) q[3];
sx q[3];
rz(-1.600986) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.071035944) q[0];
sx q[0];
rz(-1.2905916) q[0];
sx q[0];
rz(2.0284247) q[0];
rz(-0.69560266) q[1];
sx q[1];
rz(-2.7516987) q[1];
sx q[1];
rz(-1.6092469) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5672011) q[0];
sx q[0];
rz(-2.6041457) q[0];
sx q[0];
rz(0.48899942) q[0];
rz(0.74560994) q[2];
sx q[2];
rz(-0.62586212) q[2];
sx q[2];
rz(-1.2127884) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0678565) q[1];
sx q[1];
rz(-1.5508176) q[1];
sx q[1];
rz(2.0210135) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9701482) q[3];
sx q[3];
rz(-0.60895863) q[3];
sx q[3];
rz(2.3122548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9939076) q[2];
sx q[2];
rz(-2.9116178) q[2];
sx q[2];
rz(-2.2873986) q[2];
rz(1.2285852) q[3];
sx q[3];
rz(-1.5766141) q[3];
sx q[3];
rz(1.4860229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5937186) q[0];
sx q[0];
rz(-0.49867189) q[0];
sx q[0];
rz(-2.4771931) q[0];
rz(1.1876855) q[1];
sx q[1];
rz(-2.4408051) q[1];
sx q[1];
rz(1.857035) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5212095) q[0];
sx q[0];
rz(-1.8561583) q[0];
sx q[0];
rz(2.3032805) q[0];
x q[1];
rz(1.1906491) q[2];
sx q[2];
rz(-1.8188307) q[2];
sx q[2];
rz(1.5284577) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.43209546) q[1];
sx q[1];
rz(-2.2120683) q[1];
sx q[1];
rz(1.0048559) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0751238) q[3];
sx q[3];
rz(-2.7507938) q[3];
sx q[3];
rz(-1.9399411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5294042) q[2];
sx q[2];
rz(-0.89451423) q[2];
sx q[2];
rz(2.5908296) q[2];
rz(0.70358706) q[3];
sx q[3];
rz(-1.0102605) q[3];
sx q[3];
rz(1.6350869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4187014) q[0];
sx q[0];
rz(-2.7846865) q[0];
sx q[0];
rz(-3.0974467) q[0];
rz(-1.6126397) q[1];
sx q[1];
rz(-1.2395369) q[1];
sx q[1];
rz(-2.3619161) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44256193) q[0];
sx q[0];
rz(-2.583722) q[0];
sx q[0];
rz(2.872422) q[0];
rz(-pi) q[1];
x q[1];
rz(0.81515628) q[2];
sx q[2];
rz(-1.9197459) q[2];
sx q[2];
rz(2.7100035) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.838678) q[1];
sx q[1];
rz(-2.089114) q[1];
sx q[1];
rz(1.7811437) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1053106) q[3];
sx q[3];
rz(-0.84565425) q[3];
sx q[3];
rz(-2.5940965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6471275) q[2];
sx q[2];
rz(-0.72138849) q[2];
sx q[2];
rz(-1.54281) q[2];
rz(-2.4370082) q[3];
sx q[3];
rz(-1.5141809) q[3];
sx q[3];
rz(-3.1295479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71173944) q[0];
sx q[0];
rz(-1.5562417) q[0];
sx q[0];
rz(-0.65162311) q[0];
rz(1.343887) q[1];
sx q[1];
rz(-1.6603036) q[1];
sx q[1];
rz(2.4659326) q[1];
rz(2.007953) q[2];
sx q[2];
rz(-1.8532248) q[2];
sx q[2];
rz(2.1291477) q[2];
rz(-2.963831) q[3];
sx q[3];
rz(-0.64571417) q[3];
sx q[3];
rz(-0.4971102) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
