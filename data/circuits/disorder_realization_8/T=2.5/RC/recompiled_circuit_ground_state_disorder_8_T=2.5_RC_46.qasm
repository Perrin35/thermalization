OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.72579223) q[0];
sx q[0];
rz(-1.1644726) q[0];
sx q[0];
rz(-1.2512943) q[0];
rz(2.8407821) q[1];
sx q[1];
rz(-1.7164813) q[1];
sx q[1];
rz(0.11597522) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61455314) q[0];
sx q[0];
rz(-1.8924638) q[0];
sx q[0];
rz(3.0099203) q[0];
rz(-pi) q[1];
rz(1.9931727) q[2];
sx q[2];
rz(-2.6184888) q[2];
sx q[2];
rz(-0.98382271) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.310461) q[1];
sx q[1];
rz(-2.3678105) q[1];
sx q[1];
rz(-1.6602064) q[1];
rz(-pi) q[2];
x q[2];
rz(2.913353) q[3];
sx q[3];
rz(-0.55437088) q[3];
sx q[3];
rz(-1.3794464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7684324) q[2];
sx q[2];
rz(-3.0765522) q[2];
sx q[2];
rz(-1.4434641) q[2];
rz(2.4453435) q[3];
sx q[3];
rz(-1.5158451) q[3];
sx q[3];
rz(2.7878917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8732805) q[0];
sx q[0];
rz(-3.066635) q[0];
sx q[0];
rz(2.7798376) q[0];
rz(1.7573645) q[1];
sx q[1];
rz(-2.7499966) q[1];
sx q[1];
rz(1.0272383) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38798328) q[0];
sx q[0];
rz(-1.3902724) q[0];
sx q[0];
rz(-2.2645044) q[0];
rz(-pi) q[1];
rz(-1.7387668) q[2];
sx q[2];
rz(-0.78064519) q[2];
sx q[2];
rz(-1.7602242) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.69481315) q[1];
sx q[1];
rz(-0.6966335) q[1];
sx q[1];
rz(1.1778866) q[1];
rz(-pi) q[2];
rz(-0.35137041) q[3];
sx q[3];
rz(-0.23956579) q[3];
sx q[3];
rz(2.3860161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0072173) q[2];
sx q[2];
rz(-2.5773039) q[2];
sx q[2];
rz(1.0312274) q[2];
rz(2.4550896) q[3];
sx q[3];
rz(-1.7347696) q[3];
sx q[3];
rz(-0.321872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4903851) q[0];
sx q[0];
rz(-1.729916) q[0];
sx q[0];
rz(-1.7899293) q[0];
rz(1.6199813) q[1];
sx q[1];
rz(-2.9053575) q[1];
sx q[1];
rz(1.9550386) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8299415) q[0];
sx q[0];
rz(-1.1444939) q[0];
sx q[0];
rz(1.1521798) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6858052) q[2];
sx q[2];
rz(-1.6384587) q[2];
sx q[2];
rz(2.5316141) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4476984) q[1];
sx q[1];
rz(-1.5798301) q[1];
sx q[1];
rz(-1.5644521) q[1];
x q[2];
rz(-0.88079466) q[3];
sx q[3];
rz(-1.339151) q[3];
sx q[3];
rz(-0.85154545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1658907) q[2];
sx q[2];
rz(-1.2624792) q[2];
sx q[2];
rz(2.8603437) q[2];
rz(2.7665372) q[3];
sx q[3];
rz(-0.91383362) q[3];
sx q[3];
rz(2.8121172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81958812) q[0];
sx q[0];
rz(-0.81398886) q[0];
sx q[0];
rz(-2.8170012) q[0];
rz(1.5931256) q[1];
sx q[1];
rz(-0.32061583) q[1];
sx q[1];
rz(-2.2494907) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6911592) q[0];
sx q[0];
rz(-1.5884134) q[0];
sx q[0];
rz(1.5595647) q[0];
rz(-pi) q[1];
rz(0.11406819) q[2];
sx q[2];
rz(-2.4270396) q[2];
sx q[2];
rz(0.29334208) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2617787) q[1];
sx q[1];
rz(-0.38912762) q[1];
sx q[1];
rz(2.3107982) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6206498) q[3];
sx q[3];
rz(-2.5475603) q[3];
sx q[3];
rz(2.1332873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0804245) q[2];
sx q[2];
rz(-1.3865546) q[2];
sx q[2];
rz(2.5648153) q[2];
rz(-0.12442496) q[3];
sx q[3];
rz(-2.3915229) q[3];
sx q[3];
rz(-0.09593825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6782003) q[0];
sx q[0];
rz(-2.7233349) q[0];
sx q[0];
rz(2.8611355) q[0];
rz(0.29817835) q[1];
sx q[1];
rz(-3.0715946) q[1];
sx q[1];
rz(-0.15204522) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9579118) q[0];
sx q[0];
rz(-0.059558161) q[0];
sx q[0];
rz(-1.0899215) q[0];
x q[1];
rz(1.399081) q[2];
sx q[2];
rz(-2.6125312) q[2];
sx q[2];
rz(2.3429639) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3144532) q[1];
sx q[1];
rz(-1.9681853) q[1];
sx q[1];
rz(-2.023815) q[1];
rz(-pi) q[2];
rz(1.1880615) q[3];
sx q[3];
rz(-1.5788585) q[3];
sx q[3];
rz(-2.1915367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2966447) q[2];
sx q[2];
rz(-1.351202) q[2];
sx q[2];
rz(-0.039999261) q[2];
rz(-3.0039039) q[3];
sx q[3];
rz(-0.44860336) q[3];
sx q[3];
rz(-0.78534809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8765151) q[0];
sx q[0];
rz(-3.0410933) q[0];
sx q[0];
rz(2.5608089) q[0];
rz(0.68897092) q[1];
sx q[1];
rz(-0.13801485) q[1];
sx q[1];
rz(-1.8438011) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9823124) q[0];
sx q[0];
rz(-1.9907239) q[0];
sx q[0];
rz(-0.4585271) q[0];
x q[1];
rz(-1.6766729) q[2];
sx q[2];
rz(-2.5953802) q[2];
sx q[2];
rz(0.27413163) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2078198) q[1];
sx q[1];
rz(-1.1807654) q[1];
sx q[1];
rz(1.3323426) q[1];
x q[2];
rz(2.5291555) q[3];
sx q[3];
rz(-0.68842697) q[3];
sx q[3];
rz(-1.6226617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7350404) q[2];
sx q[2];
rz(-1.8031969) q[2];
sx q[2];
rz(1.5945386) q[2];
rz(2.7100587) q[3];
sx q[3];
rz(-1.1003234) q[3];
sx q[3];
rz(-0.6507473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65289098) q[0];
sx q[0];
rz(-3.1259013) q[0];
sx q[0];
rz(-0.59590644) q[0];
rz(1.3504922) q[1];
sx q[1];
rz(-0.48521438) q[1];
sx q[1];
rz(-2.5686666) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.030339) q[0];
sx q[0];
rz(-1.5736394) q[0];
sx q[0];
rz(0.009441998) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8469883) q[2];
sx q[2];
rz(-1.9825299) q[2];
sx q[2];
rz(0.48897435) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0966037) q[1];
sx q[1];
rz(-1.8879688) q[1];
sx q[1];
rz(2.5326292) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.75667824) q[3];
sx q[3];
rz(-1.8493493) q[3];
sx q[3];
rz(-0.070764489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.69221812) q[2];
sx q[2];
rz(-0.91392475) q[2];
sx q[2];
rz(-2.2268028) q[2];
rz(1.5335013) q[3];
sx q[3];
rz(-2.0017109) q[3];
sx q[3];
rz(-2.1515414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51946259) q[0];
sx q[0];
rz(-1.4254445) q[0];
sx q[0];
rz(0.064432681) q[0];
rz(2.875115) q[1];
sx q[1];
rz(-3.0173512) q[1];
sx q[1];
rz(-1.8775833) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2867133) q[0];
sx q[0];
rz(-0.42511308) q[0];
sx q[0];
rz(-0.69010587) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9585712) q[2];
sx q[2];
rz(-0.90803781) q[2];
sx q[2];
rz(-0.52802625) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.19169584) q[1];
sx q[1];
rz(-0.37572161) q[1];
sx q[1];
rz(2.5671178) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5573729) q[3];
sx q[3];
rz(-1.0845601) q[3];
sx q[3];
rz(-2.4013533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5233351) q[2];
sx q[2];
rz(-1.9674415) q[2];
sx q[2];
rz(-1.0443643) q[2];
rz(0.70762819) q[3];
sx q[3];
rz(-2.7510567) q[3];
sx q[3];
rz(-2.0948758) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3656965) q[0];
sx q[0];
rz(-1.5713659) q[0];
sx q[0];
rz(-2.7112992) q[0];
rz(-2.4039092) q[1];
sx q[1];
rz(-0.084641181) q[1];
sx q[1];
rz(2.784909) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.637054) q[0];
sx q[0];
rz(-1.86212) q[0];
sx q[0];
rz(-1.4715521) q[0];
rz(-pi) q[1];
rz(1.4595152) q[2];
sx q[2];
rz(-2.8751237) q[2];
sx q[2];
rz(0.9916208) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.99706471) q[1];
sx q[1];
rz(-2.4127919) q[1];
sx q[1];
rz(-2.4841734) q[1];
rz(1.4962219) q[3];
sx q[3];
rz(-1.6901724) q[3];
sx q[3];
rz(1.0159017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3714527) q[2];
sx q[2];
rz(-2.3944103) q[2];
sx q[2];
rz(2.746197) q[2];
rz(-1.6793518) q[3];
sx q[3];
rz(-1.3926927) q[3];
sx q[3];
rz(3.1126378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28703654) q[0];
sx q[0];
rz(-2.365878) q[0];
sx q[0];
rz(0.70242822) q[0];
rz(2.9632945) q[1];
sx q[1];
rz(-0.65419227) q[1];
sx q[1];
rz(-1.6184695) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.010945646) q[0];
sx q[0];
rz(-2.6551882) q[0];
sx q[0];
rz(2.0964699) q[0];
x q[1];
rz(0.014556292) q[2];
sx q[2];
rz(-1.6686898) q[2];
sx q[2];
rz(-2.720969) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5209841) q[1];
sx q[1];
rz(-1.0080326) q[1];
sx q[1];
rz(-1.4042793) q[1];
rz(0.63796343) q[3];
sx q[3];
rz(-2.5416982) q[3];
sx q[3];
rz(1.1092345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6452597) q[2];
sx q[2];
rz(-2.1767904) q[2];
sx q[2];
rz(-2.1920835) q[2];
rz(-2.005714) q[3];
sx q[3];
rz(-2.1954326) q[3];
sx q[3];
rz(-2.5778263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62354325) q[0];
sx q[0];
rz(-1.2508871) q[0];
sx q[0];
rz(2.293806) q[0];
rz(-1.7060948) q[1];
sx q[1];
rz(-1.2789627) q[1];
sx q[1];
rz(0.36437558) q[1];
rz(-0.87068229) q[2];
sx q[2];
rz(-1.489594) q[2];
sx q[2];
rz(1.7832047) q[2];
rz(-1.7617211) q[3];
sx q[3];
rz(-0.31287258) q[3];
sx q[3];
rz(2.3967299) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
