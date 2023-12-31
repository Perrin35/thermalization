OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.62587005) q[0];
sx q[0];
rz(-2.5928901) q[0];
sx q[0];
rz(-2.2572416) q[0];
rz(1.4305152) q[1];
sx q[1];
rz(-2.1880452) q[1];
sx q[1];
rz(1.5024827) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2384773) q[0];
sx q[0];
rz(-1.3912541) q[0];
sx q[0];
rz(-1.205501) q[0];
x q[1];
rz(-3.1332364) q[2];
sx q[2];
rz(-0.5651606) q[2];
sx q[2];
rz(-1.9985808) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1391746) q[1];
sx q[1];
rz(-2.8867509) q[1];
sx q[1];
rz(-2.8789218) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7847071) q[3];
sx q[3];
rz(-2.2483453) q[3];
sx q[3];
rz(-2.4524636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.78645906) q[2];
sx q[2];
rz(-2.3278475) q[2];
sx q[2];
rz(2.4856429) q[2];
rz(-1.2077228) q[3];
sx q[3];
rz(-1.9658807) q[3];
sx q[3];
rz(0.99457994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2475964) q[0];
sx q[0];
rz(-0.42034724) q[0];
sx q[0];
rz(2.7080652) q[0];
rz(-0.22878376) q[1];
sx q[1];
rz(-2.7119535) q[1];
sx q[1];
rz(0.0072335009) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1996961) q[0];
sx q[0];
rz(-1.4775839) q[0];
sx q[0];
rz(1.1774363) q[0];
rz(-pi) q[1];
rz(-0.29351182) q[2];
sx q[2];
rz(-1.8732757) q[2];
sx q[2];
rz(-2.7345865) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0821973) q[1];
sx q[1];
rz(-0.42947436) q[1];
sx q[1];
rz(2.5897964) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.90356566) q[3];
sx q[3];
rz(-1.9575319) q[3];
sx q[3];
rz(2.8500593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6146415) q[2];
sx q[2];
rz(-2.3336637) q[2];
sx q[2];
rz(-0.69765222) q[2];
rz(-3.0200322) q[3];
sx q[3];
rz(-1.9024885) q[3];
sx q[3];
rz(-0.30383032) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6737297) q[0];
sx q[0];
rz(-0.91600743) q[0];
sx q[0];
rz(1.3695705) q[0];
rz(1.2415775) q[1];
sx q[1];
rz(-1.4135655) q[1];
sx q[1];
rz(2.8799768) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75152552) q[0];
sx q[0];
rz(-1.5887504) q[0];
sx q[0];
rz(3.1319588) q[0];
rz(-pi) q[1];
rz(1.4357655) q[2];
sx q[2];
rz(-0.74880744) q[2];
sx q[2];
rz(2.7111862) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.47255998) q[1];
sx q[1];
rz(-1.4497888) q[1];
sx q[1];
rz(-2.3585412) q[1];
rz(-2.6286078) q[3];
sx q[3];
rz(-1.5981711) q[3];
sx q[3];
rz(-2.5740636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6802784) q[2];
sx q[2];
rz(-1.6363982) q[2];
sx q[2];
rz(0.20351163) q[2];
rz(2.2198548) q[3];
sx q[3];
rz(-1.8739871) q[3];
sx q[3];
rz(0.27954277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1638284) q[0];
sx q[0];
rz(-1.5777359) q[0];
sx q[0];
rz(1.5699566) q[0];
rz(-2.1381901) q[1];
sx q[1];
rz(-1.3137716) q[1];
sx q[1];
rz(1.8932231) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0715863) q[0];
sx q[0];
rz(-1.4741815) q[0];
sx q[0];
rz(3.0759401) q[0];
rz(-1.3893045) q[2];
sx q[2];
rz(-1.776812) q[2];
sx q[2];
rz(0.69603053) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6971671) q[1];
sx q[1];
rz(-2.4788692) q[1];
sx q[1];
rz(0.18130937) q[1];
rz(-pi) q[2];
rz(-0.29420935) q[3];
sx q[3];
rz(-0.74511408) q[3];
sx q[3];
rz(1.7780768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0288329) q[2];
sx q[2];
rz(-1.8303822) q[2];
sx q[2];
rz(0.83703414) q[2];
rz(-1.2083496) q[3];
sx q[3];
rz(-1.2669867) q[3];
sx q[3];
rz(2.3560431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0060624881) q[0];
sx q[0];
rz(-1.0536138) q[0];
sx q[0];
rz(2.3663882) q[0];
rz(-0.40183055) q[1];
sx q[1];
rz(-0.95087516) q[1];
sx q[1];
rz(2.2391589) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8474903) q[0];
sx q[0];
rz(-0.93203629) q[0];
sx q[0];
rz(1.5139447) q[0];
rz(-pi) q[1];
rz(0.47307737) q[2];
sx q[2];
rz(-0.50191754) q[2];
sx q[2];
rz(2.9830473) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8946998) q[1];
sx q[1];
rz(-2.1293853) q[1];
sx q[1];
rz(2.488399) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.0028354672) q[3];
sx q[3];
rz(-1.8875202) q[3];
sx q[3];
rz(-2.1383274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9465785) q[2];
sx q[2];
rz(-0.48823753) q[2];
sx q[2];
rz(-1.9449332) q[2];
rz(-1.442391) q[3];
sx q[3];
rz(-1.8701575) q[3];
sx q[3];
rz(1.4590013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8114132) q[0];
sx q[0];
rz(-2.9646962) q[0];
sx q[0];
rz(-2.6384171) q[0];
rz(1.4563837) q[1];
sx q[1];
rz(-1.0738942) q[1];
sx q[1];
rz(-0.20176372) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.266298) q[0];
sx q[0];
rz(-1.6365956) q[0];
sx q[0];
rz(-1.468303) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3927275) q[2];
sx q[2];
rz(-0.61486926) q[2];
sx q[2];
rz(-0.90235898) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.7355431) q[1];
sx q[1];
rz(-0.65945259) q[1];
sx q[1];
rz(2.4156648) q[1];
rz(2.4089036) q[3];
sx q[3];
rz(-1.2142039) q[3];
sx q[3];
rz(-2.4458812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9266944) q[2];
sx q[2];
rz(-1.639067) q[2];
sx q[2];
rz(-0.26724896) q[2];
rz(-2.3184508) q[3];
sx q[3];
rz(-0.079113364) q[3];
sx q[3];
rz(0.87583035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.293752) q[0];
sx q[0];
rz(-2.9822615) q[0];
sx q[0];
rz(0.057549495) q[0];
rz(1.6607704) q[1];
sx q[1];
rz(-1.7819504) q[1];
sx q[1];
rz(2.1988791) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7563815) q[0];
sx q[0];
rz(-2.4009973) q[0];
sx q[0];
rz(-2.5368607) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1812181) q[2];
sx q[2];
rz(-0.27563169) q[2];
sx q[2];
rz(-0.20197091) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.645694) q[1];
sx q[1];
rz(-1.2518479) q[1];
sx q[1];
rz(-2.0046528) q[1];
x q[2];
rz(-2.0833756) q[3];
sx q[3];
rz(-1.7985385) q[3];
sx q[3];
rz(0.77373576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1192347) q[2];
sx q[2];
rz(-1.0417577) q[2];
sx q[2];
rz(2.7589202) q[2];
rz(2.102397) q[3];
sx q[3];
rz(-1.8363876) q[3];
sx q[3];
rz(2.1634845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4246178) q[0];
sx q[0];
rz(-3.0452403) q[0];
sx q[0];
rz(0.27012816) q[0];
rz(-0.62942901) q[1];
sx q[1];
rz(-0.71297485) q[1];
sx q[1];
rz(-0.28392917) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0547202) q[0];
sx q[0];
rz(-0.76612872) q[0];
sx q[0];
rz(-2.2148569) q[0];
rz(-1.7173041) q[2];
sx q[2];
rz(-2.1657145) q[2];
sx q[2];
rz(-2.811424) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.64360196) q[1];
sx q[1];
rz(-2.0703348) q[1];
sx q[1];
rz(1.040578) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4817804) q[3];
sx q[3];
rz(-1.5471349) q[3];
sx q[3];
rz(-0.82994474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5876864) q[2];
sx q[2];
rz(-0.17051414) q[2];
sx q[2];
rz(-1.2109057) q[2];
rz(0.25990137) q[3];
sx q[3];
rz(-2.5172958) q[3];
sx q[3];
rz(0.62817812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-0.07638409) q[0];
sx q[0];
rz(-0.56088352) q[0];
sx q[0];
rz(-2.912345) q[0];
rz(2.8385838) q[1];
sx q[1];
rz(-1.7508933) q[1];
sx q[1];
rz(-1.680826) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6630336) q[0];
sx q[0];
rz(-1.8543108) q[0];
sx q[0];
rz(1.5525596) q[0];
rz(2.7621208) q[2];
sx q[2];
rz(-1.4037637) q[2];
sx q[2];
rz(2.6118979) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2158828) q[1];
sx q[1];
rz(-2.6162418) q[1];
sx q[1];
rz(-0.044563091) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7542354) q[3];
sx q[3];
rz(-0.82327561) q[3];
sx q[3];
rz(-0.92126095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.71172697) q[2];
sx q[2];
rz(-1.9248328) q[2];
sx q[2];
rz(0.6955859) q[2];
rz(-2.7097278) q[3];
sx q[3];
rz(-2.6769107) q[3];
sx q[3];
rz(2.4263884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5678976) q[0];
sx q[0];
rz(-1.9298113) q[0];
sx q[0];
rz(2.8046872) q[0];
rz(-0.20740549) q[1];
sx q[1];
rz(-2.1131056) q[1];
sx q[1];
rz(-2.7609603) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3363627) q[0];
sx q[0];
rz(-0.21348937) q[0];
sx q[0];
rz(1.9070894) q[0];
x q[1];
rz(2.400488) q[2];
sx q[2];
rz(-2.0217102) q[2];
sx q[2];
rz(-1.8795183) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.29585782) q[1];
sx q[1];
rz(-1.0746135) q[1];
sx q[1];
rz(-0.22175281) q[1];
rz(2.5467039) q[3];
sx q[3];
rz(-2.4747362) q[3];
sx q[3];
rz(-0.24645933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.0011065817) q[2];
sx q[2];
rz(-1.9753186) q[2];
sx q[2];
rz(2.005119) q[2];
rz(-3.100637) q[3];
sx q[3];
rz(-2.332873) q[3];
sx q[3];
rz(-1.3142746) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3363591) q[0];
sx q[0];
rz(-2.2362066) q[0];
sx q[0];
rz(2.8295828) q[0];
rz(-2.1144755) q[1];
sx q[1];
rz(-1.849091) q[1];
sx q[1];
rz(-1.0277933) q[1];
rz(0.74577352) q[2];
sx q[2];
rz(-0.33200982) q[2];
sx q[2];
rz(0.40340323) q[2];
rz(2.0102262) q[3];
sx q[3];
rz(-1.3888748) q[3];
sx q[3];
rz(0.5007762) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
