OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.80469552) q[0];
sx q[0];
rz(-1.0372294) q[0];
sx q[0];
rz(-2.7859935) q[0];
rz(2.8388677) q[1];
sx q[1];
rz(5.2390715) q[1];
sx q[1];
rz(10.704438) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5582433) q[0];
sx q[0];
rz(-1.8177176) q[0];
sx q[0];
rz(-0.14351828) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7945292) q[2];
sx q[2];
rz(-1.1587843) q[2];
sx q[2];
rz(1.7681233) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.31561786) q[1];
sx q[1];
rz(-0.74233913) q[1];
sx q[1];
rz(0.62645341) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7765331) q[3];
sx q[3];
rz(-1.7997026) q[3];
sx q[3];
rz(2.8347901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0063643) q[2];
sx q[2];
rz(-2.2029115) q[2];
sx q[2];
rz(-0.65650666) q[2];
rz(-0.73900977) q[3];
sx q[3];
rz(-2.6813172) q[3];
sx q[3];
rz(-2.7242993) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1332557) q[0];
sx q[0];
rz(-0.79611859) q[0];
sx q[0];
rz(2.546229) q[0];
rz(3.0796675) q[1];
sx q[1];
rz(-1.2281111) q[1];
sx q[1];
rz(0.48746902) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2330403) q[0];
sx q[0];
rz(-1.6516343) q[0];
sx q[0];
rz(3.0768865) q[0];
rz(0.42178085) q[2];
sx q[2];
rz(-1.1676719) q[2];
sx q[2];
rz(-0.031907206) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4844955) q[1];
sx q[1];
rz(-0.24689455) q[1];
sx q[1];
rz(0.8685649) q[1];
x q[2];
rz(0.08667605) q[3];
sx q[3];
rz(-1.8868586) q[3];
sx q[3];
rz(-2.812127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9618824) q[2];
sx q[2];
rz(-1.2625182) q[2];
sx q[2];
rz(2.7462192) q[2];
rz(-1.0428492) q[3];
sx q[3];
rz(-2.6104749) q[3];
sx q[3];
rz(3.1029207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9449126) q[0];
sx q[0];
rz(-2.0356464) q[0];
sx q[0];
rz(-0.19038598) q[0];
rz(3.0186675) q[1];
sx q[1];
rz(-0.38750896) q[1];
sx q[1];
rz(-0.22274676) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5771772) q[0];
sx q[0];
rz(-2.5860997) q[0];
sx q[0];
rz(0.76340686) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.22521714) q[2];
sx q[2];
rz(-1.9334963) q[2];
sx q[2];
rz(-0.97937102) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.60274117) q[1];
sx q[1];
rz(-0.95898421) q[1];
sx q[1];
rz(1.4879423) q[1];
rz(-pi) q[2];
rz(-2.1915216) q[3];
sx q[3];
rz(-2.7347703) q[3];
sx q[3];
rz(2.1086958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2591851) q[2];
sx q[2];
rz(-1.2732482) q[2];
sx q[2];
rz(-1.1941236) q[2];
rz(1.0549226) q[3];
sx q[3];
rz(-1.6566365) q[3];
sx q[3];
rz(0.96364337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7524183) q[0];
sx q[0];
rz(-0.27802935) q[0];
sx q[0];
rz(1.4032723) q[0];
rz(-1.4933043) q[1];
sx q[1];
rz(-1.5416668) q[1];
sx q[1];
rz(0.9202252) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.614914) q[0];
sx q[0];
rz(-0.81256142) q[0];
sx q[0];
rz(2.148669) q[0];
rz(-pi) q[1];
rz(2.9158981) q[2];
sx q[2];
rz(-2.7719471) q[2];
sx q[2];
rz(-0.29633488) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8060382) q[1];
sx q[1];
rz(-1.1127377) q[1];
sx q[1];
rz(2.4469417) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3654877) q[3];
sx q[3];
rz(-1.9234386) q[3];
sx q[3];
rz(0.40311381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9437287) q[2];
sx q[2];
rz(-1.7284164) q[2];
sx q[2];
rz(1.8614004) q[2];
rz(-2.3222893) q[3];
sx q[3];
rz(-1.8236482) q[3];
sx q[3];
rz(-1.7839446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.458805) q[0];
sx q[0];
rz(-1.3764494) q[0];
sx q[0];
rz(1.4047594) q[0];
rz(-0.76830307) q[1];
sx q[1];
rz(-2.4914425) q[1];
sx q[1];
rz(2.6884902) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4015409) q[0];
sx q[0];
rz(-1.6443559) q[0];
sx q[0];
rz(-0.037365035) q[0];
x q[1];
rz(0.034899072) q[2];
sx q[2];
rz(-1.7136095) q[2];
sx q[2];
rz(2.3226483) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6140155) q[1];
sx q[1];
rz(-1.7909389) q[1];
sx q[1];
rz(-3.0575391) q[1];
rz(-pi) q[2];
rz(-2.1702607) q[3];
sx q[3];
rz(-1.1820275) q[3];
sx q[3];
rz(2.142981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.12864628) q[2];
sx q[2];
rz(-0.83798989) q[2];
sx q[2];
rz(1.5117234) q[2];
rz(2.417918) q[3];
sx q[3];
rz(-1.2440163) q[3];
sx q[3];
rz(-0.85030142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1081651) q[0];
sx q[0];
rz(-1.7280248) q[0];
sx q[0];
rz(-2.9034555) q[0];
rz(0.37995964) q[1];
sx q[1];
rz(-2.0894876) q[1];
sx q[1];
rz(-1.628081) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0244004) q[0];
sx q[0];
rz(-1.8937366) q[0];
sx q[0];
rz(-0.49844235) q[0];
rz(0.061821826) q[2];
sx q[2];
rz(-2.0197868) q[2];
sx q[2];
rz(2.2668554) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7442419) q[1];
sx q[1];
rz(-0.89351082) q[1];
sx q[1];
rz(-2.1703297) q[1];
rz(-pi) q[2];
rz(-0.87168872) q[3];
sx q[3];
rz(-0.5629493) q[3];
sx q[3];
rz(-2.9339919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0825519) q[2];
sx q[2];
rz(-2.5566792) q[2];
sx q[2];
rz(-1.9980105) q[2];
rz(-3.011076) q[3];
sx q[3];
rz(-1.427622) q[3];
sx q[3];
rz(-2.9706764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1259574) q[0];
sx q[0];
rz(-0.97706777) q[0];
sx q[0];
rz(-2.7440199) q[0];
rz(1.827084) q[1];
sx q[1];
rz(-2.59771) q[1];
sx q[1];
rz(1.1869173) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0163527) q[0];
sx q[0];
rz(-1.851474) q[0];
sx q[0];
rz(-1.0513845) q[0];
rz(-pi) q[1];
x q[1];
rz(0.45849623) q[2];
sx q[2];
rz(-3.0756604) q[2];
sx q[2];
rz(2.3022848) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6656875) q[1];
sx q[1];
rz(-0.74531065) q[1];
sx q[1];
rz(1.3728632) q[1];
x q[2];
rz(0.62992923) q[3];
sx q[3];
rz(-1.3327206) q[3];
sx q[3];
rz(-1.5743953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.26178965) q[2];
sx q[2];
rz(-1.3061085) q[2];
sx q[2];
rz(-1.3558033) q[2];
rz(-1.6342182) q[3];
sx q[3];
rz(-2.5528788) q[3];
sx q[3];
rz(0.9986977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3223406) q[0];
sx q[0];
rz(-1.1871908) q[0];
sx q[0];
rz(-0.85574714) q[0];
rz(-3.1198655) q[1];
sx q[1];
rz(-1.0236434) q[1];
sx q[1];
rz(-2.1112679) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4521342) q[0];
sx q[0];
rz(-1.9118714) q[0];
sx q[0];
rz(2.336691) q[0];
rz(-pi) q[1];
rz(-3.0929504) q[2];
sx q[2];
rz(-2.0970793) q[2];
sx q[2];
rz(-0.7359879) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.1134539) q[1];
sx q[1];
rz(-2.8061562) q[1];
sx q[1];
rz(-1.0949275) q[1];
rz(-pi) q[2];
rz(0.52328531) q[3];
sx q[3];
rz(-2.0998294) q[3];
sx q[3];
rz(-2.7494591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8490303) q[2];
sx q[2];
rz(-1.2529255) q[2];
sx q[2];
rz(-3.0714152) q[2];
rz(2.3146546) q[3];
sx q[3];
rz(-1.4708054) q[3];
sx q[3];
rz(-0.58644811) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6475911) q[0];
sx q[0];
rz(-1.7140472) q[0];
sx q[0];
rz(-2.8253187) q[0];
rz(2.0896185) q[1];
sx q[1];
rz(-0.72967044) q[1];
sx q[1];
rz(1.1351599) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90056706) q[0];
sx q[0];
rz(-1.7591488) q[0];
sx q[0];
rz(0.39649773) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6584381) q[2];
sx q[2];
rz(-0.32155514) q[2];
sx q[2];
rz(-1.1631539) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.70442048) q[1];
sx q[1];
rz(-1.6295027) q[1];
sx q[1];
rz(1.0874332) q[1];
rz(-1.6039861) q[3];
sx q[3];
rz(-2.5182708) q[3];
sx q[3];
rz(-2.8896796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6241374) q[2];
sx q[2];
rz(-0.74666658) q[2];
sx q[2];
rz(-1.5448145) q[2];
rz(0.67772135) q[3];
sx q[3];
rz(-0.8422519) q[3];
sx q[3];
rz(0.28963447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.233376) q[0];
sx q[0];
rz(-0.69013086) q[0];
sx q[0];
rz(2.7897575) q[0];
rz(0.31967638) q[1];
sx q[1];
rz(-0.37477481) q[1];
sx q[1];
rz(0.19616729) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10712121) q[0];
sx q[0];
rz(-2.2959318) q[0];
sx q[0];
rz(-2.3496773) q[0];
x q[1];
rz(1.7061383) q[2];
sx q[2];
rz(-1.5208828) q[2];
sx q[2];
rz(2.4382255) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6876467) q[1];
sx q[1];
rz(-1.5554264) q[1];
sx q[1];
rz(-2.3049298) q[1];
rz(-pi) q[2];
rz(-2.4862643) q[3];
sx q[3];
rz(-2.1851551) q[3];
sx q[3];
rz(-2.4650246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7508042) q[2];
sx q[2];
rz(-1.6395586) q[2];
sx q[2];
rz(3.0604559) q[2];
rz(-1.9598512) q[3];
sx q[3];
rz(-0.90098444) q[3];
sx q[3];
rz(2.0666163) q[3];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1186196) q[0];
sx q[0];
rz(-1.5587627) q[0];
sx q[0];
rz(-1.8297304) q[0];
rz(1.3148057) q[1];
sx q[1];
rz(-2.343315) q[1];
sx q[1];
rz(1.5409484) q[1];
rz(2.9076004) q[2];
sx q[2];
rz(-1.5237612) q[2];
sx q[2];
rz(0.16618726) q[2];
rz(1.8004988) q[3];
sx q[3];
rz(-1.5256186) q[3];
sx q[3];
rz(-1.516173) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
