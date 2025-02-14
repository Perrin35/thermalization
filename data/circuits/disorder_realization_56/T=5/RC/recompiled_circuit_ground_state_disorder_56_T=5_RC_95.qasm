OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.0290282) q[0];
sx q[0];
rz(-1.5052786) q[0];
sx q[0];
rz(2.3583052) q[0];
rz(1.0822436) q[1];
sx q[1];
rz(-1.487027) q[1];
sx q[1];
rz(-2.3496871) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5675674) q[0];
sx q[0];
rz(-1.8301395) q[0];
sx q[0];
rz(-1.7631084) q[0];
x q[1];
rz(-1.3375086) q[2];
sx q[2];
rz(-0.15809862) q[2];
sx q[2];
rz(2.5776517) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.72630771) q[1];
sx q[1];
rz(-2.1139164) q[1];
sx q[1];
rz(-1.5341137) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5654985) q[3];
sx q[3];
rz(-0.8818501) q[3];
sx q[3];
rz(-2.9722633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.59114328) q[2];
sx q[2];
rz(-2.998896) q[2];
sx q[2];
rz(3.0421416) q[2];
rz(2.8948696) q[3];
sx q[3];
rz(-1.6171425) q[3];
sx q[3];
rz(-2.7439086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1219015) q[0];
sx q[0];
rz(-0.97625232) q[0];
sx q[0];
rz(0.9285399) q[0];
rz(-1.4506725) q[1];
sx q[1];
rz(-1.2932777) q[1];
sx q[1];
rz(0.64847747) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47005338) q[0];
sx q[0];
rz(-2.3190365) q[0];
sx q[0];
rz(-2.3683192) q[0];
x q[1];
rz(0.47467741) q[2];
sx q[2];
rz(-2.5203278) q[2];
sx q[2];
rz(0.23532993) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9916223) q[1];
sx q[1];
rz(-0.14279616) q[1];
sx q[1];
rz(-3.0317329) q[1];
x q[2];
rz(-0.11827151) q[3];
sx q[3];
rz(-1.2840446) q[3];
sx q[3];
rz(-0.82055446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6917981) q[2];
sx q[2];
rz(-2.3851676) q[2];
sx q[2];
rz(2.2890384) q[2];
rz(-2.2954588) q[3];
sx q[3];
rz(-1.6328014) q[3];
sx q[3];
rz(1.030863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5071252) q[0];
sx q[0];
rz(-1.2256624) q[0];
sx q[0];
rz(-2.0563828) q[0];
rz(0.94404864) q[1];
sx q[1];
rz(-1.6429106) q[1];
sx q[1];
rz(1.5012213) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2417898) q[0];
sx q[0];
rz(-0.2010303) q[0];
sx q[0];
rz(0.19636671) q[0];
rz(-pi) q[1];
rz(-0.19135059) q[2];
sx q[2];
rz(-1.3406959) q[2];
sx q[2];
rz(1.5887345) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0090364) q[1];
sx q[1];
rz(-1.6699797) q[1];
sx q[1];
rz(2.1445091) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.15764938) q[3];
sx q[3];
rz(-1.356989) q[3];
sx q[3];
rz(-2.714954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.40470716) q[2];
sx q[2];
rz(-2.6077304) q[2];
sx q[2];
rz(2.286818) q[2];
rz(0.9084304) q[3];
sx q[3];
rz(-1.5646489) q[3];
sx q[3];
rz(-2.0250208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92723769) q[0];
sx q[0];
rz(-2.7404009) q[0];
sx q[0];
rz(1.9450564) q[0];
rz(-1.475097) q[1];
sx q[1];
rz(-0.42088446) q[1];
sx q[1];
rz(-1.9570785) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23789433) q[0];
sx q[0];
rz(-1.8082976) q[0];
sx q[0];
rz(2.4303275) q[0];
rz(-pi) q[1];
rz(-2.938561) q[2];
sx q[2];
rz(-0.76329279) q[2];
sx q[2];
rz(-2.610746) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8912011) q[1];
sx q[1];
rz(-1.7970835) q[1];
sx q[1];
rz(1.0550642) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7538025) q[3];
sx q[3];
rz(-1.840045) q[3];
sx q[3];
rz(-1.8257917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2148332) q[2];
sx q[2];
rz(-0.49761179) q[2];
sx q[2];
rz(-1.8801749) q[2];
rz(-0.43241209) q[3];
sx q[3];
rz(-2.0093446) q[3];
sx q[3];
rz(-2.6692218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0686907) q[0];
sx q[0];
rz(-1.2304767) q[0];
sx q[0];
rz(-0.029408971) q[0];
rz(0.75621653) q[1];
sx q[1];
rz(-0.58719802) q[1];
sx q[1];
rz(0.75278935) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8653976) q[0];
sx q[0];
rz(-1.5704234) q[0];
sx q[0];
rz(-1.5651817) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2529066) q[2];
sx q[2];
rz(-0.28159062) q[2];
sx q[2];
rz(-1.7624621) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.009234) q[1];
sx q[1];
rz(-2.1387324) q[1];
sx q[1];
rz(-1.6949937) q[1];
x q[2];
rz(1.1208956) q[3];
sx q[3];
rz(-2.0664762) q[3];
sx q[3];
rz(-0.97199398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.14043643) q[2];
sx q[2];
rz(-1.6480646) q[2];
sx q[2];
rz(-2.746554) q[2];
rz(0.47518528) q[3];
sx q[3];
rz(-1.7893712) q[3];
sx q[3];
rz(2.7648259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3756556) q[0];
sx q[0];
rz(-0.7203311) q[0];
sx q[0];
rz(2.1790867) q[0];
rz(-0.51482254) q[1];
sx q[1];
rz(-1.5341026) q[1];
sx q[1];
rz(2.8033676) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0926577) q[0];
sx q[0];
rz(-0.83507628) q[0];
sx q[0];
rz(2.8754763) q[0];
rz(-pi) q[1];
x q[1];
rz(0.2372431) q[2];
sx q[2];
rz(-1.3297992) q[2];
sx q[2];
rz(2.4324696) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(0.40336762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6414791) q[2];
sx q[2];
rz(-1.1621472) q[2];
sx q[2];
rz(-2.5872453) q[2];
rz(0.85136271) q[3];
sx q[3];
rz(-0.35836372) q[3];
sx q[3];
rz(2.5135777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3250658) q[0];
sx q[0];
rz(-2.5287703) q[0];
sx q[0];
rz(0.75041962) q[0];
rz(-2.5665307) q[1];
sx q[1];
rz(-1.776418) q[1];
sx q[1];
rz(0.65779984) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.039666273) q[0];
sx q[0];
rz(-0.79625477) q[0];
sx q[0];
rz(2.5523561) q[0];
x q[1];
rz(0.70722981) q[2];
sx q[2];
rz(-1.5373188) q[2];
sx q[2];
rz(2.1113124) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7791351) q[1];
sx q[1];
rz(-1.5541847) q[1];
sx q[1];
rz(-1.2846867) q[1];
rz(-pi) q[2];
rz(-1.9522454) q[3];
sx q[3];
rz(-1.7717517) q[3];
sx q[3];
rz(-2.902346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.647992) q[2];
sx q[2];
rz(-1.0085663) q[2];
sx q[2];
rz(-1.4823401) q[2];
rz(-3.0209387) q[3];
sx q[3];
rz(-1.3841265) q[3];
sx q[3];
rz(2.306126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3487314) q[0];
sx q[0];
rz(-2.272235) q[0];
sx q[0];
rz(-0.29801512) q[0];
rz(-2.1030078) q[1];
sx q[1];
rz(-0.54439259) q[1];
sx q[1];
rz(2.9878152) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6385495) q[0];
sx q[0];
rz(-0.26493236) q[0];
sx q[0];
rz(-2.8517836) q[0];
rz(-0.80259097) q[2];
sx q[2];
rz(-2.0611412) q[2];
sx q[2];
rz(-0.56832321) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7419305) q[1];
sx q[1];
rz(-2.7832418) q[1];
sx q[1];
rz(-2.8835924) q[1];
x q[2];
rz(3.0197487) q[3];
sx q[3];
rz(-1.1466371) q[3];
sx q[3];
rz(-2.3736853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0457354) q[2];
sx q[2];
rz(-0.46892527) q[2];
sx q[2];
rz(-1.9528961) q[2];
rz(1.8111604) q[3];
sx q[3];
rz(-1.1878139) q[3];
sx q[3];
rz(-2.4082898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7981912) q[0];
sx q[0];
rz(-0.60438406) q[0];
sx q[0];
rz(1.6424302) q[0];
rz(-2.5921953) q[1];
sx q[1];
rz(-1.5645942) q[1];
sx q[1];
rz(-2.8909491) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3651705) q[0];
sx q[0];
rz(-0.93342268) q[0];
sx q[0];
rz(0.38743069) q[0];
rz(-pi) q[1];
rz(-1.5852087) q[2];
sx q[2];
rz(-1.8317128) q[2];
sx q[2];
rz(0.018785611) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1085775) q[1];
sx q[1];
rz(-2.3653154) q[1];
sx q[1];
rz(-0.099221512) q[1];
x q[2];
rz(-0.19480494) q[3];
sx q[3];
rz(-1.4470295) q[3];
sx q[3];
rz(-0.74461246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9833019) q[2];
sx q[2];
rz(-3.0088708) q[2];
sx q[2];
rz(1.0164725) q[2];
rz(3.0554092) q[3];
sx q[3];
rz(-1.0165756) q[3];
sx q[3];
rz(0.76519722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30855274) q[0];
sx q[0];
rz(-2.2873531) q[0];
sx q[0];
rz(2.7225851) q[0];
rz(0.53681701) q[1];
sx q[1];
rz(-0.62428004) q[1];
sx q[1];
rz(2.4057665) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2013071) q[0];
sx q[0];
rz(-2.2129411) q[0];
sx q[0];
rz(1.8117732) q[0];
rz(-pi) q[1];
rz(0.64401099) q[2];
sx q[2];
rz(-2.0154872) q[2];
sx q[2];
rz(1.4294238) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9089953) q[1];
sx q[1];
rz(-0.67282644) q[1];
sx q[1];
rz(-0.35079591) q[1];
x q[2];
rz(-0.48218547) q[3];
sx q[3];
rz(-1.8819359) q[3];
sx q[3];
rz(0.37267123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1824823) q[2];
sx q[2];
rz(-2.3190494) q[2];
sx q[2];
rz(2.518892) q[2];
rz(0.73838082) q[3];
sx q[3];
rz(-1.9564956) q[3];
sx q[3];
rz(-0.27899376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57711346) q[0];
sx q[0];
rz(-2.8688685) q[0];
sx q[0];
rz(-2.175749) q[0];
rz(0.64361698) q[1];
sx q[1];
rz(-1.8543961) q[1];
sx q[1];
rz(-2.054945) q[1];
rz(2.116133) q[2];
sx q[2];
rz(-2.7164216) q[2];
sx q[2];
rz(2.807775) q[2];
rz(0.24091992) q[3];
sx q[3];
rz(-1.708247) q[3];
sx q[3];
rz(0.74549992) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
