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
rz(-1.2404233) q[0];
sx q[0];
rz(-1.197553) q[0];
sx q[0];
rz(2.9252606) q[0];
rz(-2.1842015) q[1];
sx q[1];
rz(-0.67617813) q[1];
sx q[1];
rz(1.6300936) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0759685) q[0];
sx q[0];
rz(-1.0145717) q[0];
sx q[0];
rz(1.5594622) q[0];
rz(-pi) q[1];
x q[1];
rz(1.897282) q[2];
sx q[2];
rz(-2.0225581) q[2];
sx q[2];
rz(-1.2690074) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.87147698) q[1];
sx q[1];
rz(-1.2387215) q[1];
sx q[1];
rz(-2.6144652) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7483398) q[3];
sx q[3];
rz(-0.90688469) q[3];
sx q[3];
rz(2.8877986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1067074) q[2];
sx q[2];
rz(-2.7896176) q[2];
sx q[2];
rz(2.0931639) q[2];
rz(0.18167051) q[3];
sx q[3];
rz(-0.964966) q[3];
sx q[3];
rz(-2.488193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(1.6492017) q[0];
sx q[0];
rz(-2.1972456) q[0];
sx q[0];
rz(2.6948068) q[0];
rz(0.88042879) q[1];
sx q[1];
rz(-1.7767521) q[1];
sx q[1];
rz(-0.78278881) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0271757) q[0];
sx q[0];
rz(-1.3265298) q[0];
sx q[0];
rz(3.0784803) q[0];
rz(-pi) q[1];
rz(-0.57116429) q[2];
sx q[2];
rz(-2.4413902) q[2];
sx q[2];
rz(-0.68298662) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0228717) q[1];
sx q[1];
rz(-1.1069655) q[1];
sx q[1];
rz(1.2341586) q[1];
rz(-pi) q[2];
rz(2.1332333) q[3];
sx q[3];
rz(-2.376412) q[3];
sx q[3];
rz(-2.6443114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.11144) q[2];
sx q[2];
rz(-1.4806662) q[2];
sx q[2];
rz(-2.1622369) q[2];
rz(2.7566946) q[3];
sx q[3];
rz(-1.220547) q[3];
sx q[3];
rz(0.16429193) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48591831) q[0];
sx q[0];
rz(-0.05412183) q[0];
sx q[0];
rz(2.3517877) q[0];
rz(-2.9549331) q[1];
sx q[1];
rz(-1.4013545) q[1];
sx q[1];
rz(0.99367118) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45897608) q[0];
sx q[0];
rz(-1.0280711) q[0];
sx q[0];
rz(1.2786464) q[0];
rz(2.8254059) q[2];
sx q[2];
rz(-2.5044166) q[2];
sx q[2];
rz(3.1115465) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.34383306) q[1];
sx q[1];
rz(-1.1145089) q[1];
sx q[1];
rz(-0.47618687) q[1];
rz(-pi) q[2];
rz(-1.2993811) q[3];
sx q[3];
rz(-2.2732353) q[3];
sx q[3];
rz(2.0347119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8197202) q[2];
sx q[2];
rz(-1.1178144) q[2];
sx q[2];
rz(-0.37128386) q[2];
rz(-0.34879455) q[3];
sx q[3];
rz(-2.0472417) q[3];
sx q[3];
rz(0.5955407) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45469859) q[0];
sx q[0];
rz(-2.1214387) q[0];
sx q[0];
rz(-1.4042847) q[0];
rz(-2.4644201) q[1];
sx q[1];
rz(-1.155747) q[1];
sx q[1];
rz(1.4926532) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9335404) q[0];
sx q[0];
rz(-1.3553936) q[0];
sx q[0];
rz(-0.24788863) q[0];
rz(-pi) q[1];
rz(-1.1392966) q[2];
sx q[2];
rz(-1.2988663) q[2];
sx q[2];
rz(-1.6587342) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3918152) q[1];
sx q[1];
rz(-1.7931374) q[1];
sx q[1];
rz(-0.27395497) q[1];
rz(2.2954313) q[3];
sx q[3];
rz(-1.807888) q[3];
sx q[3];
rz(2.2711636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6877785) q[2];
sx q[2];
rz(-2.1731589) q[2];
sx q[2];
rz(-0.34823927) q[2];
rz(1.4549152) q[3];
sx q[3];
rz(-1.7125407) q[3];
sx q[3];
rz(1.7961563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1763879) q[0];
sx q[0];
rz(-2.5768953) q[0];
sx q[0];
rz(0.89214605) q[0];
rz(0.46420321) q[1];
sx q[1];
rz(-1.2449539) q[1];
sx q[1];
rz(-1.2947327) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2703122) q[0];
sx q[0];
rz(-2.0593606) q[0];
sx q[0];
rz(1.0315007) q[0];
x q[1];
rz(1.4424397) q[2];
sx q[2];
rz(-1.0367298) q[2];
sx q[2];
rz(0.72557025) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.059493493) q[1];
sx q[1];
rz(-2.0096742) q[1];
sx q[1];
rz(2.782269) q[1];
x q[2];
rz(2.5312349) q[3];
sx q[3];
rz(-0.98516432) q[3];
sx q[3];
rz(-2.741339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8190454) q[2];
sx q[2];
rz(-2.5307541) q[2];
sx q[2];
rz(-0.56582212) q[2];
rz(-0.081341751) q[3];
sx q[3];
rz(-0.9809202) q[3];
sx q[3];
rz(2.5210023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6767947) q[0];
sx q[0];
rz(-0.066635266) q[0];
sx q[0];
rz(1.5860522) q[0];
rz(-2.0809035) q[1];
sx q[1];
rz(-1.565275) q[1];
sx q[1];
rz(-0.63180822) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6752967) q[0];
sx q[0];
rz(-1.491437) q[0];
sx q[0];
rz(0.23510374) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1366227) q[2];
sx q[2];
rz(-0.7938876) q[2];
sx q[2];
rz(-0.53295202) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.374208) q[1];
sx q[1];
rz(-0.65113089) q[1];
sx q[1];
rz(2.1506449) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.17978823) q[3];
sx q[3];
rz(-0.69151141) q[3];
sx q[3];
rz(-0.19516842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1598728) q[2];
sx q[2];
rz(-1.1184511) q[2];
sx q[2];
rz(2.6427606) q[2];
rz(-1.3129129) q[3];
sx q[3];
rz(-2.4098318) q[3];
sx q[3];
rz(1.4341199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6071534) q[0];
sx q[0];
rz(-0.3648912) q[0];
sx q[0];
rz(-0.69951192) q[0];
rz(0.46547678) q[1];
sx q[1];
rz(-0.87124467) q[1];
sx q[1];
rz(-0.93719283) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58317157) q[0];
sx q[0];
rz(-1.4025549) q[0];
sx q[0];
rz(2.9267163) q[0];
rz(-pi) q[1];
rz(1.466339) q[2];
sx q[2];
rz(-1.8091996) q[2];
sx q[2];
rz(-1.6676211) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4027462) q[1];
sx q[1];
rz(-0.35772309) q[1];
sx q[1];
rz(2.2132232) q[1];
rz(-pi) q[2];
rz(0.77009691) q[3];
sx q[3];
rz(-0.40989629) q[3];
sx q[3];
rz(2.485825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.297544) q[2];
sx q[2];
rz(-1.8397477) q[2];
sx q[2];
rz(2.76827) q[2];
rz(1.1897872) q[3];
sx q[3];
rz(-0.50896421) q[3];
sx q[3];
rz(0.51026195) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3968286) q[0];
sx q[0];
rz(-2.0592392) q[0];
sx q[0];
rz(2.3936791) q[0];
rz(0.76639908) q[1];
sx q[1];
rz(-0.26873573) q[1];
sx q[1];
rz(-0.0029729923) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9342614) q[0];
sx q[0];
rz(-0.89326477) q[0];
sx q[0];
rz(-0.24768655) q[0];
rz(0.019549088) q[2];
sx q[2];
rz(-1.4998933) q[2];
sx q[2];
rz(-2.8343458) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0743979) q[1];
sx q[1];
rz(-1.41681) q[1];
sx q[1];
rz(-0.043937307) q[1];
rz(0.79487309) q[3];
sx q[3];
rz(-0.27471457) q[3];
sx q[3];
rz(-0.87546722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.11091867) q[2];
sx q[2];
rz(-1.7577533) q[2];
sx q[2];
rz(-1.0670916) q[2];
rz(-3.057632) q[3];
sx q[3];
rz(-0.48560086) q[3];
sx q[3];
rz(2.3626204) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.344051) q[0];
sx q[0];
rz(-2.2028956) q[0];
sx q[0];
rz(3.0352266) q[0];
rz(2.1616409) q[1];
sx q[1];
rz(-1.6500902) q[1];
sx q[1];
rz(-0.76593691) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87493103) q[0];
sx q[0];
rz(-0.70737544) q[0];
sx q[0];
rz(-0.69702638) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7226549) q[2];
sx q[2];
rz(-0.21809245) q[2];
sx q[2];
rz(0.60148394) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.76141833) q[1];
sx q[1];
rz(-1.5360218) q[1];
sx q[1];
rz(2.4654287) q[1];
rz(-pi) q[2];
rz(-1.5670789) q[3];
sx q[3];
rz(-1.5848985) q[3];
sx q[3];
rz(-0.64025793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6672259) q[2];
sx q[2];
rz(-0.52046481) q[2];
sx q[2];
rz(2.4412947) q[2];
rz(2.4750366) q[3];
sx q[3];
rz(-1.3349814) q[3];
sx q[3];
rz(1.0880281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4661082) q[0];
sx q[0];
rz(-1.0270783) q[0];
sx q[0];
rz(1.745537) q[0];
rz(1.667977) q[1];
sx q[1];
rz(-1.8831848) q[1];
sx q[1];
rz(2.4748763) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5294801) q[0];
sx q[0];
rz(-0.94769883) q[0];
sx q[0];
rz(2.3701131) q[0];
rz(1.2856977) q[2];
sx q[2];
rz(-2.210493) q[2];
sx q[2];
rz(-2.6486047) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4426431) q[1];
sx q[1];
rz(-1.2546228) q[1];
sx q[1];
rz(1.7085646) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.75717775) q[3];
sx q[3];
rz(-0.72523967) q[3];
sx q[3];
rz(-2.4717769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.066976808) q[2];
sx q[2];
rz(-0.65698996) q[2];
sx q[2];
rz(2.2508049) q[2];
rz(-0.43205076) q[3];
sx q[3];
rz(-2.0712349) q[3];
sx q[3];
rz(-0.74705684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4074832) q[0];
sx q[0];
rz(-1.4984087) q[0];
sx q[0];
rz(1.8468504) q[0];
rz(1.2474077) q[1];
sx q[1];
rz(-0.71949646) q[1];
sx q[1];
rz(1.234642) q[1];
rz(-0.16130372) q[2];
sx q[2];
rz(-1.2304753) q[2];
sx q[2];
rz(-0.58438042) q[2];
rz(1.6040989) q[3];
sx q[3];
rz(-2.0596444) q[3];
sx q[3];
rz(2.8465908) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
