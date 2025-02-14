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
rz(0.12240527) q[0];
sx q[0];
rz(-0.91949099) q[0];
sx q[0];
rz(-1.1991731) q[0];
rz(-2.9543258) q[1];
sx q[1];
rz(-0.54228243) q[1];
sx q[1];
rz(-1.6220925) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47177256) q[0];
sx q[0];
rz(-0.61994055) q[0];
sx q[0];
rz(-0.13373904) q[0];
rz(-2.7678732) q[2];
sx q[2];
rz(-1.9961832) q[2];
sx q[2];
rz(-0.56217867) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4368545) q[1];
sx q[1];
rz(-0.66009843) q[1];
sx q[1];
rz(0.69003244) q[1];
rz(-pi) q[2];
rz(1.5222129) q[3];
sx q[3];
rz(-1.3268927) q[3];
sx q[3];
rz(0.78256202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2970994) q[2];
sx q[2];
rz(-0.62701925) q[2];
sx q[2];
rz(1.7681047) q[2];
rz(2.3376236) q[3];
sx q[3];
rz(-1.6425902) q[3];
sx q[3];
rz(1.9544301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74541575) q[0];
sx q[0];
rz(-0.71393037) q[0];
sx q[0];
rz(-2.7667238) q[0];
rz(0.51775852) q[1];
sx q[1];
rz(-1.9381783) q[1];
sx q[1];
rz(1.3688603) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8086209) q[0];
sx q[0];
rz(-0.00095168984) q[0];
sx q[0];
rz(0.10097058) q[0];
rz(-1.5977741) q[2];
sx q[2];
rz(-1.7671874) q[2];
sx q[2];
rz(2.2324004) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6102099) q[1];
sx q[1];
rz(-2.068559) q[1];
sx q[1];
rz(2.0396292) q[1];
x q[2];
rz(-1.3563223) q[3];
sx q[3];
rz(-1.703959) q[3];
sx q[3];
rz(1.565763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.123473) q[2];
sx q[2];
rz(-0.70088434) q[2];
sx q[2];
rz(0.71462053) q[2];
rz(-2.9361652) q[3];
sx q[3];
rz(-1.8193865) q[3];
sx q[3];
rz(-1.8429168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6767204) q[0];
sx q[0];
rz(-1.0370075) q[0];
sx q[0];
rz(2.6439457) q[0];
rz(2.5045555) q[1];
sx q[1];
rz(-0.62890816) q[1];
sx q[1];
rz(-0.10674891) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3321132) q[0];
sx q[0];
rz(-2.7710072) q[0];
sx q[0];
rz(-2.1912712) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2405143) q[2];
sx q[2];
rz(-1.5759829) q[2];
sx q[2];
rz(-0.55904311) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.52912583) q[1];
sx q[1];
rz(-2.4262316) q[1];
sx q[1];
rz(-1.6523408) q[1];
rz(0.62415267) q[3];
sx q[3];
rz(-2.660203) q[3];
sx q[3];
rz(2.2311051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.60483021) q[2];
sx q[2];
rz(-0.75211516) q[2];
sx q[2];
rz(-0.19482782) q[2];
rz(3.1365862) q[3];
sx q[3];
rz(-2.5275793) q[3];
sx q[3];
rz(-0.79202882) q[3];
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
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0358589) q[0];
sx q[0];
rz(-0.53426131) q[0];
sx q[0];
rz(2.1015097) q[0];
rz(-0.19373521) q[1];
sx q[1];
rz(-1.5015142) q[1];
sx q[1];
rz(-0.74660444) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0950985) q[0];
sx q[0];
rz(-1.455716) q[0];
sx q[0];
rz(-2.0128573) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7393635) q[2];
sx q[2];
rz(-2.4344846) q[2];
sx q[2];
rz(-0.76729028) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.50987096) q[1];
sx q[1];
rz(-1.6225909) q[1];
sx q[1];
rz(2.447261) q[1];
x q[2];
rz(-2.901731) q[3];
sx q[3];
rz(-1.6768394) q[3];
sx q[3];
rz(-0.21940809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.99183434) q[2];
sx q[2];
rz(-1.4703625) q[2];
sx q[2];
rz(2.9370918) q[2];
rz(-1.1931984) q[3];
sx q[3];
rz(-2.2843993) q[3];
sx q[3];
rz(-2.1804555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0123154) q[0];
sx q[0];
rz(-2.9670872) q[0];
sx q[0];
rz(2.0611064) q[0];
rz(-0.37733817) q[1];
sx q[1];
rz(-1.5107379) q[1];
sx q[1];
rz(-2.2468755) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52496929) q[0];
sx q[0];
rz(-1.7979413) q[0];
sx q[0];
rz(0.36360111) q[0];
x q[1];
rz(-0.65308833) q[2];
sx q[2];
rz(-0.90735596) q[2];
sx q[2];
rz(1.9651679) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1929568) q[1];
sx q[1];
rz(-2.0564449) q[1];
sx q[1];
rz(2.5705277) q[1];
rz(-pi) q[2];
rz(1.0977488) q[3];
sx q[3];
rz(-0.69907197) q[3];
sx q[3];
rz(3.0300273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6696024) q[2];
sx q[2];
rz(-2.695638) q[2];
sx q[2];
rz(-0.10776821) q[2];
rz(-0.17555155) q[3];
sx q[3];
rz(-1.8864417) q[3];
sx q[3];
rz(2.5811783) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.046722978) q[0];
sx q[0];
rz(-0.51114285) q[0];
sx q[0];
rz(1.0119337) q[0];
rz(1.5049505) q[1];
sx q[1];
rz(-0.71957809) q[1];
sx q[1];
rz(1.0171657) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5115248) q[0];
sx q[0];
rz(-1.2753133) q[0];
sx q[0];
rz(-2.9620671) q[0];
rz(-pi) q[1];
rz(0.58508137) q[2];
sx q[2];
rz(-0.85137109) q[2];
sx q[2];
rz(0.78851267) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6222657) q[1];
sx q[1];
rz(-2.9510088) q[1];
sx q[1];
rz(-3.0960073) q[1];
rz(-pi) q[2];
rz(-1.8886376) q[3];
sx q[3];
rz(-0.46964619) q[3];
sx q[3];
rz(1.8152678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.73198685) q[2];
sx q[2];
rz(-2.4797532) q[2];
sx q[2];
rz(0.99676639) q[2];
rz(-0.064920001) q[3];
sx q[3];
rz(-1.4109979) q[3];
sx q[3];
rz(-1.1499278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0885334) q[0];
sx q[0];
rz(-2.2011338) q[0];
sx q[0];
rz(2.233182) q[0];
rz(-2.9216595) q[1];
sx q[1];
rz(-1.5426153) q[1];
sx q[1];
rz(1.8904846) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0742038) q[0];
sx q[0];
rz(-0.54710623) q[0];
sx q[0];
rz(2.0873468) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8562137) q[2];
sx q[2];
rz(-1.6661465) q[2];
sx q[2];
rz(0.91172632) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0798222) q[1];
sx q[1];
rz(-2.1822189) q[1];
sx q[1];
rz(-2.3794214) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3014017) q[3];
sx q[3];
rz(-1.6641656) q[3];
sx q[3];
rz(1.8871117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0327586) q[2];
sx q[2];
rz(-2.3075576) q[2];
sx q[2];
rz(-1.533482) q[2];
rz(-1.1533302) q[3];
sx q[3];
rz(-1.0554375) q[3];
sx q[3];
rz(1.6495033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9913427) q[0];
sx q[0];
rz(-0.38337502) q[0];
sx q[0];
rz(0.94069329) q[0];
rz(-3.0062145) q[1];
sx q[1];
rz(-1.4043413) q[1];
sx q[1];
rz(1.0248331) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.301046) q[0];
sx q[0];
rz(-2.4468166) q[0];
sx q[0];
rz(-1.077681) q[0];
rz(-pi) q[1];
rz(2.5338855) q[2];
sx q[2];
rz(-1.6204699) q[2];
sx q[2];
rz(-2.7756754) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.35541818) q[1];
sx q[1];
rz(-0.23705951) q[1];
sx q[1];
rz(-0.60225822) q[1];
rz(-pi) q[2];
rz(0.41127326) q[3];
sx q[3];
rz(-1.6182634) q[3];
sx q[3];
rz(2.2293343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3465053) q[2];
sx q[2];
rz(-2.2515209) q[2];
sx q[2];
rz(-0.65004641) q[2];
rz(2.5551689) q[3];
sx q[3];
rz(-1.4805877) q[3];
sx q[3];
rz(-0.75064269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2670249) q[0];
sx q[0];
rz(-0.50849193) q[0];
sx q[0];
rz(1.1453999) q[0];
rz(-1.3151431) q[1];
sx q[1];
rz(-1.9086842) q[1];
sx q[1];
rz(-0.68793908) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4015539) q[0];
sx q[0];
rz(-1.7595707) q[0];
sx q[0];
rz(-1.1925634) q[0];
rz(-pi) q[1];
x q[1];
rz(0.39040651) q[2];
sx q[2];
rz(-2.1320754) q[2];
sx q[2];
rz(-2.3048603) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.44837828) q[1];
sx q[1];
rz(-0.82469392) q[1];
sx q[1];
rz(1.7918827) q[1];
x q[2];
rz(0.93291544) q[3];
sx q[3];
rz(-2.5543) q[3];
sx q[3];
rz(1.6577394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.73100662) q[2];
sx q[2];
rz(-1.1062016) q[2];
sx q[2];
rz(2.2753687) q[2];
rz(-0.90803641) q[3];
sx q[3];
rz(-2.3767545) q[3];
sx q[3];
rz(-2.0456555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1753801) q[0];
sx q[0];
rz(-1.1840273) q[0];
sx q[0];
rz(-0.41561919) q[0];
rz(0.79187727) q[1];
sx q[1];
rz(-1.2245347) q[1];
sx q[1];
rz(-2.9143639) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4675605) q[0];
sx q[0];
rz(-1.5825543) q[0];
sx q[0];
rz(3.1255162) q[0];
rz(-pi) q[1];
rz(3.0640494) q[2];
sx q[2];
rz(-2.3334586) q[2];
sx q[2];
rz(-1.6412954) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7586582) q[1];
sx q[1];
rz(-0.99731748) q[1];
sx q[1];
rz(-1.2632779) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0437743) q[3];
sx q[3];
rz(-1.1036345) q[3];
sx q[3];
rz(-2.710619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.26433429) q[2];
sx q[2];
rz(-0.72089583) q[2];
sx q[2];
rz(-2.7016675) q[2];
rz(2.544493) q[3];
sx q[3];
rz(-2.4720981) q[3];
sx q[3];
rz(-1.7841024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6482342) q[0];
sx q[0];
rz(-1.5879205) q[0];
sx q[0];
rz(1.2863202) q[0];
rz(-1.413912) q[1];
sx q[1];
rz(-1.021011) q[1];
sx q[1];
rz(0.11722142) q[1];
rz(-3.0138409) q[2];
sx q[2];
rz(-2.7314039) q[2];
sx q[2];
rz(1.9990767) q[2];
rz(2.8942378) q[3];
sx q[3];
rz(-2.3079899) q[3];
sx q[3];
rz(0.079903825) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
