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
rz(-0.51198045) q[0];
sx q[0];
rz(-2.5706302) q[0];
sx q[0];
rz(-2.9538739) q[0];
rz(0.81387782) q[1];
sx q[1];
rz(4.6133572) q[1];
sx q[1];
rz(7.5724966) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58285125) q[0];
sx q[0];
rz(-1.8076573) q[0];
sx q[0];
rz(1.357097) q[0];
rz(-pi) q[1];
rz(0.57065771) q[2];
sx q[2];
rz(-1.5957812) q[2];
sx q[2];
rz(-0.61793426) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6179168) q[1];
sx q[1];
rz(-2.6111228) q[1];
sx q[1];
rz(2.6350849) q[1];
rz(-pi) q[2];
rz(1.9240407) q[3];
sx q[3];
rz(-1.788967) q[3];
sx q[3];
rz(-2.6835203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8744371) q[2];
sx q[2];
rz(-1.1031373) q[2];
sx q[2];
rz(-0.77679408) q[2];
rz(-1.8393501) q[3];
sx q[3];
rz(-1.4812508) q[3];
sx q[3];
rz(-1.9384025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60748196) q[0];
sx q[0];
rz(-0.48180875) q[0];
sx q[0];
rz(2.1433461) q[0];
rz(-0.36969319) q[1];
sx q[1];
rz(-2.1001215) q[1];
sx q[1];
rz(2.0236156) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8000172) q[0];
sx q[0];
rz(-1.1693926) q[0];
sx q[0];
rz(-0.97824162) q[0];
x q[1];
rz(-0.68685617) q[2];
sx q[2];
rz(-1.1026898) q[2];
sx q[2];
rz(1.8091701) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4395968) q[1];
sx q[1];
rz(-2.0148104) q[1];
sx q[1];
rz(2.9786112) q[1];
x q[2];
rz(1.2106154) q[3];
sx q[3];
rz(-2.4781215) q[3];
sx q[3];
rz(-1.0960032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2758241) q[2];
sx q[2];
rz(-1.3752702) q[2];
sx q[2];
rz(2.7562874) q[2];
rz(-0.038711874) q[3];
sx q[3];
rz(-2.7888515) q[3];
sx q[3];
rz(-0.64046162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63967079) q[0];
sx q[0];
rz(-2.6646035) q[0];
sx q[0];
rz(1.3013526) q[0];
rz(0.057706984) q[1];
sx q[1];
rz(-0.4464018) q[1];
sx q[1];
rz(1.3267964) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2952174) q[0];
sx q[0];
rz(-1.0344939) q[0];
sx q[0];
rz(-0.55732507) q[0];
x q[1];
rz(1.7102555) q[2];
sx q[2];
rz(-1.7296975) q[2];
sx q[2];
rz(0.1316084) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2542299) q[1];
sx q[1];
rz(-2.0944203) q[1];
sx q[1];
rz(1.9033405) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1744099) q[3];
sx q[3];
rz(-0.76226888) q[3];
sx q[3];
rz(1.9844733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.79746276) q[2];
sx q[2];
rz(-2.3064488) q[2];
sx q[2];
rz(0.75378913) q[2];
rz(1.7720743) q[3];
sx q[3];
rz(-1.4621719) q[3];
sx q[3];
rz(-0.65521017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21514431) q[0];
sx q[0];
rz(-2.0560052) q[0];
sx q[0];
rz(1.3947067) q[0];
rz(-1.1766379) q[1];
sx q[1];
rz(-1.5086915) q[1];
sx q[1];
rz(-2.7746157) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4876539) q[0];
sx q[0];
rz(-1.5040845) q[0];
sx q[0];
rz(3.0848461) q[0];
rz(0.66195935) q[2];
sx q[2];
rz(-2.7267704) q[2];
sx q[2];
rz(2.7560459) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.48880348) q[1];
sx q[1];
rz(-1.7646043) q[1];
sx q[1];
rz(-0.36313063) q[1];
rz(2.9212679) q[3];
sx q[3];
rz(-0.37543618) q[3];
sx q[3];
rz(0.32849778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7299812) q[2];
sx q[2];
rz(-1.9133762) q[2];
sx q[2];
rz(-1.8709987) q[2];
rz(2.1805084) q[3];
sx q[3];
rz(-1.7226487) q[3];
sx q[3];
rz(-3.0911176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56080317) q[0];
sx q[0];
rz(-0.92629782) q[0];
sx q[0];
rz(1.3128989) q[0];
rz(-1.987223) q[1];
sx q[1];
rz(-1.1416953) q[1];
sx q[1];
rz(-2.5837574) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9168222) q[0];
sx q[0];
rz(-2.0799412) q[0];
sx q[0];
rz(-2.1780685) q[0];
rz(0.16061546) q[2];
sx q[2];
rz(-1.6470419) q[2];
sx q[2];
rz(2.3996224) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.5231042) q[1];
sx q[1];
rz(-0.84977154) q[1];
sx q[1];
rz(0.40722653) q[1];
rz(-2.6799503) q[3];
sx q[3];
rz(-1.4351234) q[3];
sx q[3];
rz(1.5380579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4981726) q[2];
sx q[2];
rz(-1.8489445) q[2];
sx q[2];
rz(0.84929973) q[2];
rz(2.9863206) q[3];
sx q[3];
rz(-2.1657491) q[3];
sx q[3];
rz(0.37080216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3575386) q[0];
sx q[0];
rz(-2.5548866) q[0];
sx q[0];
rz(-0.46911711) q[0];
rz(-0.47438374) q[1];
sx q[1];
rz(-0.7862888) q[1];
sx q[1];
rz(2.3862086) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3123388) q[0];
sx q[0];
rz(-0.27129284) q[0];
sx q[0];
rz(-1.5543227) q[0];
x q[1];
rz(1.9017436) q[2];
sx q[2];
rz(-0.96815434) q[2];
sx q[2];
rz(-1.861426) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1052221) q[1];
sx q[1];
rz(-1.3579988) q[1];
sx q[1];
rz(-0.12963055) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7128025) q[3];
sx q[3];
rz(-0.8800104) q[3];
sx q[3];
rz(0.59861983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0834121) q[2];
sx q[2];
rz(-1.8184793) q[2];
sx q[2];
rz(0.18784909) q[2];
rz(1.8958873) q[3];
sx q[3];
rz(-1.3789504) q[3];
sx q[3];
rz(-2.0511621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2096527) q[0];
sx q[0];
rz(-1.8745475) q[0];
sx q[0];
rz(2.7506822) q[0];
rz(-0.93005013) q[1];
sx q[1];
rz(-1.4314194) q[1];
sx q[1];
rz(1.8375058) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.025607312) q[0];
sx q[0];
rz(-1.5918333) q[0];
sx q[0];
rz(-1.9217092) q[0];
rz(-pi) q[1];
rz(-1.485059) q[2];
sx q[2];
rz(-1.9112327) q[2];
sx q[2];
rz(0.96735937) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.89829373) q[1];
sx q[1];
rz(-2.0212681) q[1];
sx q[1];
rz(-2.9551278) q[1];
x q[2];
rz(-2.845302) q[3];
sx q[3];
rz(-1.5888831) q[3];
sx q[3];
rz(2.2294105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8730674) q[2];
sx q[2];
rz(-2.8748942) q[2];
sx q[2];
rz(3.0282057) q[2];
rz(0.015965613) q[3];
sx q[3];
rz(-1.4957875) q[3];
sx q[3];
rz(-2.9769843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3575344) q[0];
sx q[0];
rz(-1.4736195) q[0];
sx q[0];
rz(-3.0175324) q[0];
rz(-0.1217753) q[1];
sx q[1];
rz(-0.75141326) q[1];
sx q[1];
rz(1.442499) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0589813) q[0];
sx q[0];
rz(-1.9537203) q[0];
sx q[0];
rz(3.0331066) q[0];
rz(-0.29192544) q[2];
sx q[2];
rz(-1.5291844) q[2];
sx q[2];
rz(2.6440563) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9297819) q[1];
sx q[1];
rz(-2.7719797) q[1];
sx q[1];
rz(-1.1951642) q[1];
rz(-pi) q[2];
rz(1.8201039) q[3];
sx q[3];
rz(-2.2852906) q[3];
sx q[3];
rz(-0.51333237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.96559912) q[2];
sx q[2];
rz(-1.7815353) q[2];
sx q[2];
rz(2.3967801) q[2];
rz(2.2143769) q[3];
sx q[3];
rz(-1.5085446) q[3];
sx q[3];
rz(1.0923227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48831478) q[0];
sx q[0];
rz(-1.4396242) q[0];
sx q[0];
rz(2.9162245) q[0];
rz(2.0291746) q[1];
sx q[1];
rz(-2.367159) q[1];
sx q[1];
rz(0.89881277) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32466896) q[0];
sx q[0];
rz(-0.82624871) q[0];
sx q[0];
rz(-0.81314317) q[0];
rz(-pi) q[1];
rz(-2.4821539) q[2];
sx q[2];
rz(-0.62022479) q[2];
sx q[2];
rz(-0.31977113) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.73279335) q[1];
sx q[1];
rz(-0.91161722) q[1];
sx q[1];
rz(-0.2562457) q[1];
rz(-pi) q[2];
rz(0.01047666) q[3];
sx q[3];
rz(-1.8636522) q[3];
sx q[3];
rz(-1.7550857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0972458) q[2];
sx q[2];
rz(-1.4298507) q[2];
sx q[2];
rz(-0.35935768) q[2];
rz(-0.25092956) q[3];
sx q[3];
rz(-2.0693306) q[3];
sx q[3];
rz(-0.76771626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74442416) q[0];
sx q[0];
rz(-3.011062) q[0];
sx q[0];
rz(0.8771483) q[0];
rz(-1.5769222) q[1];
sx q[1];
rz(-1.501187) q[1];
sx q[1];
rz(0.78561479) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79646677) q[0];
sx q[0];
rz(-0.32497999) q[0];
sx q[0];
rz(-2.3152405) q[0];
rz(0.53391407) q[2];
sx q[2];
rz(-1.9451687) q[2];
sx q[2];
rz(1.716734) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.12745107) q[1];
sx q[1];
rz(-1.2268865) q[1];
sx q[1];
rz(-1.4959072) q[1];
rz(-0.50650017) q[3];
sx q[3];
rz(-1.439038) q[3];
sx q[3];
rz(1.3598639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.91763532) q[2];
sx q[2];
rz(-2.3515297) q[2];
sx q[2];
rz(2.4033974) q[2];
rz(-2.2186642) q[3];
sx q[3];
rz(-1.3083369) q[3];
sx q[3];
rz(0.2963399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34687635) q[0];
sx q[0];
rz(-1.7807757) q[0];
sx q[0];
rz(-2.1697252) q[0];
rz(2.234266) q[1];
sx q[1];
rz(-2.0569888) q[1];
sx q[1];
rz(2.8203698) q[1];
rz(-0.99030607) q[2];
sx q[2];
rz(-1.2377501) q[2];
sx q[2];
rz(-0.35828423) q[2];
rz(-0.54013822) q[3];
sx q[3];
rz(-2.028699) q[3];
sx q[3];
rz(1.5069458) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
