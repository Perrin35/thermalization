OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.11256448) q[0];
sx q[0];
rz(-1.636314) q[0];
sx q[0];
rz(-2.3583052) q[0];
rz(1.0822436) q[1];
sx q[1];
rz(-1.487027) q[1];
sx q[1];
rz(0.79190555) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5740252) q[0];
sx q[0];
rz(-1.8301395) q[0];
sx q[0];
rz(1.7631084) q[0];
rz(-pi) q[1];
rz(-1.804084) q[2];
sx q[2];
rz(-0.15809862) q[2];
sx q[2];
rz(0.563941) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.79719964) q[1];
sx q[1];
rz(-0.54423344) q[1];
sx q[1];
rz(0.060677008) q[1];
rz(-pi) q[2];
rz(0.79444076) q[3];
sx q[3];
rz(-1.1367961) q[3];
sx q[3];
rz(-1.7931012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5504494) q[2];
sx q[2];
rz(-0.14269665) q[2];
sx q[2];
rz(-3.0421416) q[2];
rz(2.8948696) q[3];
sx q[3];
rz(-1.5244502) q[3];
sx q[3];
rz(2.7439086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0196911) q[0];
sx q[0];
rz(-0.97625232) q[0];
sx q[0];
rz(-0.9285399) q[0];
rz(1.4506725) q[1];
sx q[1];
rz(-1.848315) q[1];
sx q[1];
rz(-2.4931152) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6715393) q[0];
sx q[0];
rz(-0.82255615) q[0];
sx q[0];
rz(0.7732735) q[0];
x q[1];
rz(1.8869867) q[2];
sx q[2];
rz(-1.0266227) q[2];
sx q[2];
rz(0.7989102) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8295171) q[1];
sx q[1];
rz(-1.5551928) q[1];
sx q[1];
rz(-2.9996458) q[1];
rz(-1.9514328) q[3];
sx q[3];
rz(-2.8320304) q[3];
sx q[3];
rz(2.7187687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4497946) q[2];
sx q[2];
rz(-0.75642502) q[2];
sx q[2];
rz(-0.85255426) q[2];
rz(0.84613386) q[3];
sx q[3];
rz(-1.5087912) q[3];
sx q[3];
rz(-1.030863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5071252) q[0];
sx q[0];
rz(-1.9159303) q[0];
sx q[0];
rz(2.0563828) q[0];
rz(-2.197544) q[1];
sx q[1];
rz(-1.6429106) q[1];
sx q[1];
rz(1.5012213) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0414934) q[0];
sx q[0];
rz(-1.3736808) q[0];
sx q[0];
rz(-1.5310578) q[0];
rz(-pi) q[1];
rz(-0.88884647) q[2];
sx q[2];
rz(-2.8434128) q[2];
sx q[2];
rz(-2.2928638) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2861917) q[1];
sx q[1];
rz(-0.58127379) q[1];
sx q[1];
rz(-1.3894677) q[1];
rz(0.15764938) q[3];
sx q[3];
rz(-1.356989) q[3];
sx q[3];
rz(-0.42663867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.40470716) q[2];
sx q[2];
rz(-0.53386226) q[2];
sx q[2];
rz(-0.85477465) q[2];
rz(2.2331623) q[3];
sx q[3];
rz(-1.5646489) q[3];
sx q[3];
rz(2.0250208) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92723769) q[0];
sx q[0];
rz(-0.4011918) q[0];
sx q[0];
rz(-1.1965363) q[0];
rz(-1.6664956) q[1];
sx q[1];
rz(-2.7207082) q[1];
sx q[1];
rz(1.1845142) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5329342) q[0];
sx q[0];
rz(-2.2581165) q[0];
sx q[0];
rz(1.88009) q[0];
rz(-0.20303161) q[2];
sx q[2];
rz(-2.3782999) q[2];
sx q[2];
rz(0.53084669) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0851497) q[1];
sx q[1];
rz(-2.5825325) q[1];
sx q[1];
rz(-2.0075625) q[1];
rz(-pi) q[2];
x q[2];
rz(0.27359815) q[3];
sx q[3];
rz(-1.3944542) q[3];
sx q[3];
rz(0.30418744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.92675942) q[2];
sx q[2];
rz(-2.6439809) q[2];
sx q[2];
rz(-1.8801749) q[2];
rz(0.43241209) q[3];
sx q[3];
rz(-1.132248) q[3];
sx q[3];
rz(0.47237083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0686907) q[0];
sx q[0];
rz(-1.2304767) q[0];
sx q[0];
rz(3.1121837) q[0];
rz(2.3853761) q[1];
sx q[1];
rz(-0.58719802) q[1];
sx q[1];
rz(-0.75278935) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29460337) q[0];
sx q[0];
rz(-1.5764109) q[0];
sx q[0];
rz(-0.0003729781) q[0];
rz(-pi) q[1];
rz(1.8886861) q[2];
sx q[2];
rz(-2.860002) q[2];
sx q[2];
rz(1.3791305) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.009234) q[1];
sx q[1];
rz(-2.1387324) q[1];
sx q[1];
rz(1.4465989) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0206971) q[3];
sx q[3];
rz(-2.0664762) q[3];
sx q[3];
rz(2.1695987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0011562) q[2];
sx q[2];
rz(-1.6480646) q[2];
sx q[2];
rz(-2.746554) q[2];
rz(-2.6664074) q[3];
sx q[3];
rz(-1.3522215) q[3];
sx q[3];
rz(-2.7648259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3756556) q[0];
sx q[0];
rz(-0.7203311) q[0];
sx q[0];
rz(-2.1790867) q[0];
rz(-0.51482254) q[1];
sx q[1];
rz(-1.6074901) q[1];
sx q[1];
rz(0.33822507) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7027959) q[0];
sx q[0];
rz(-1.3745752) q[0];
sx q[0];
rz(-0.81721925) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3338396) q[2];
sx q[2];
rz(-2.8050426) q[2];
sx q[2];
rz(1.5010264) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5023313) q[1];
sx q[1];
rz(-3.1384472) q[1];
sx q[1];
rz(1.5622713) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.148726) q[3];
sx q[3];
rz(-1.633612) q[3];
sx q[3];
rz(-2.738225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6414791) q[2];
sx q[2];
rz(-1.9794455) q[2];
sx q[2];
rz(-2.5872453) q[2];
rz(0.85136271) q[3];
sx q[3];
rz(-2.7832289) q[3];
sx q[3];
rz(-2.5135777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8165269) q[0];
sx q[0];
rz(-0.61282235) q[0];
sx q[0];
rz(-2.391173) q[0];
rz(0.57506192) q[1];
sx q[1];
rz(-1.776418) q[1];
sx q[1];
rz(0.65779984) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80245012) q[0];
sx q[0];
rz(-2.2070812) q[0];
sx q[0];
rz(-2.0872981) q[0];
rz(-pi) q[1];
rz(-1.5267685) q[2];
sx q[2];
rz(-2.2775473) q[2];
sx q[2];
rz(0.51191521) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7791351) q[1];
sx q[1];
rz(-1.5541847) q[1];
sx q[1];
rz(1.2846867) q[1];
rz(-pi) q[2];
rz(-1.9522454) q[3];
sx q[3];
rz(-1.7717517) q[3];
sx q[3];
rz(0.23924669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.49360069) q[2];
sx q[2];
rz(-2.1330264) q[2];
sx q[2];
rz(1.4823401) q[2];
rz(0.12065398) q[3];
sx q[3];
rz(-1.3841265) q[3];
sx q[3];
rz(2.306126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7928612) q[0];
sx q[0];
rz(-2.272235) q[0];
sx q[0];
rz(0.29801512) q[0];
rz(-1.0385849) q[1];
sx q[1];
rz(-2.5972001) q[1];
sx q[1];
rz(2.9878152) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9382291) q[0];
sx q[0];
rz(-1.8244317) q[0];
sx q[0];
rz(1.4934191) q[0];
rz(2.3390017) q[2];
sx q[2];
rz(-1.0804515) q[2];
sx q[2];
rz(-2.5732694) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7281935) q[1];
sx q[1];
rz(-1.4811885) q[1];
sx q[1];
rz(0.34743584) q[1];
x q[2];
rz(1.9977536) q[3];
sx q[3];
rz(-1.4597963) q[3];
sx q[3];
rz(-0.75253651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0457354) q[2];
sx q[2];
rz(-0.46892527) q[2];
sx q[2];
rz(1.9528961) q[2];
rz(1.8111604) q[3];
sx q[3];
rz(-1.9537787) q[3];
sx q[3];
rz(2.4082898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34340149) q[0];
sx q[0];
rz(-2.5372086) q[0];
sx q[0];
rz(-1.4991624) q[0];
rz(-0.54939735) q[1];
sx q[1];
rz(-1.5769985) q[1];
sx q[1];
rz(-2.8909491) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1090013) q[0];
sx q[0];
rz(-1.2622941) q[0];
sx q[0];
rz(-2.2453813) q[0];
x q[1];
rz(-0.05392404) q[2];
sx q[2];
rz(-0.26130518) q[2];
sx q[2];
rz(-0.037029412) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.10560606) q[1];
sx q[1];
rz(-0.79933724) q[1];
sx q[1];
rz(-1.4738333) q[1];
rz(-pi) q[2];
x q[2];
rz(0.19480494) q[3];
sx q[3];
rz(-1.6945632) q[3];
sx q[3];
rz(-0.74461246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9833019) q[2];
sx q[2];
rz(-0.13272186) q[2];
sx q[2];
rz(2.1251202) q[2];
rz(-0.086183444) q[3];
sx q[3];
rz(-2.1250171) q[3];
sx q[3];
rz(-0.76519722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8330399) q[0];
sx q[0];
rz(-2.2873531) q[0];
sx q[0];
rz(-0.41900751) q[0];
rz(2.6047756) q[1];
sx q[1];
rz(-0.62428004) q[1];
sx q[1];
rz(0.73582617) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2013071) q[0];
sx q[0];
rz(-2.2129411) q[0];
sx q[0];
rz(-1.8117732) q[0];
rz(0.64401099) q[2];
sx q[2];
rz(-2.0154872) q[2];
sx q[2];
rz(1.4294238) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9365471) q[1];
sx q[1];
rz(-0.94561316) q[1];
sx q[1];
rz(-1.8380828) q[1];
rz(-1.222615) q[3];
sx q[3];
rz(-2.0280119) q[3];
sx q[3];
rz(-1.3570076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.95911038) q[2];
sx q[2];
rz(-2.3190494) q[2];
sx q[2];
rz(-2.518892) q[2];
rz(-2.4032118) q[3];
sx q[3];
rz(-1.1850971) q[3];
sx q[3];
rz(0.27899376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5644792) q[0];
sx q[0];
rz(-2.8688685) q[0];
sx q[0];
rz(-2.175749) q[0];
rz(0.64361698) q[1];
sx q[1];
rz(-1.8543961) q[1];
sx q[1];
rz(-2.054945) q[1];
rz(2.9109091) q[2];
sx q[2];
rz(-1.9311957) q[2];
sx q[2];
rz(0.25372505) q[2];
rz(1.4293115) q[3];
sx q[3];
rz(-1.8093997) q[3];
sx q[3];
rz(2.3499478) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
