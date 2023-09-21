OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.52580994) q[0];
sx q[0];
rz(4.5594112) q[0];
sx q[0];
rz(8.863908) q[0];
rz(1.1129192) q[1];
sx q[1];
rz(-1.7634044) q[1];
sx q[1];
rz(1.2150432) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5117447) q[0];
sx q[0];
rz(-1.6748322) q[0];
sx q[0];
rz(-1.758979) q[0];
rz(1.2826074) q[2];
sx q[2];
rz(-2.211314) q[2];
sx q[2];
rz(0.033601947) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.558555) q[1];
sx q[1];
rz(-1.2416632) q[1];
sx q[1];
rz(-1.0298883) q[1];
rz(-2.9925572) q[3];
sx q[3];
rz(-1.3081074) q[3];
sx q[3];
rz(1.5460154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.78757301) q[2];
sx q[2];
rz(-0.95280567) q[2];
sx q[2];
rz(-2.9585178) q[2];
rz(2.7637774) q[3];
sx q[3];
rz(-1.0487882) q[3];
sx q[3];
rz(-2.8474076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8437682) q[0];
sx q[0];
rz(-2.4968708) q[0];
sx q[0];
rz(0.077117292) q[0];
rz(2.8027957) q[1];
sx q[1];
rz(-2.0270551) q[1];
sx q[1];
rz(1.5391301) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30277006) q[0];
sx q[0];
rz(-1.4937703) q[0];
sx q[0];
rz(-2.1588615) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3537172) q[2];
sx q[2];
rz(-2.2976544) q[2];
sx q[2];
rz(-0.94490563) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.96670818) q[1];
sx q[1];
rz(-1.0918573) q[1];
sx q[1];
rz(0.19660463) q[1];
rz(-pi) q[2];
rz(-1.8061403) q[3];
sx q[3];
rz(-1.9841521) q[3];
sx q[3];
rz(1.2873161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.845528) q[2];
sx q[2];
rz(-1.8639996) q[2];
sx q[2];
rz(-0.65845931) q[2];
rz(-0.15130875) q[3];
sx q[3];
rz(-2.1189809) q[3];
sx q[3];
rz(2.4466799) q[3];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028458683) q[0];
sx q[0];
rz(-0.78335339) q[0];
sx q[0];
rz(0.43310305) q[0];
rz(1.1921047) q[1];
sx q[1];
rz(-1.9299709) q[1];
sx q[1];
rz(0.55535299) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0387602) q[0];
sx q[0];
rz(-1.0881249) q[0];
sx q[0];
rz(2.32248) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0054587) q[2];
sx q[2];
rz(-1.8381422) q[2];
sx q[2];
rz(-2.8440059) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.024228) q[1];
sx q[1];
rz(-1.7324565) q[1];
sx q[1];
rz(2.2092186) q[1];
x q[2];
rz(0.20136307) q[3];
sx q[3];
rz(-2.106973) q[3];
sx q[3];
rz(-0.45504967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.73734036) q[2];
sx q[2];
rz(-2.3534687) q[2];
sx q[2];
rz(1.8910485) q[2];
rz(0.2441497) q[3];
sx q[3];
rz(-1.8593676) q[3];
sx q[3];
rz(1.6916493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8811532) q[0];
sx q[0];
rz(-2.6840211) q[0];
sx q[0];
rz(-0.81480169) q[0];
rz(1.762215) q[1];
sx q[1];
rz(-0.35019362) q[1];
sx q[1];
rz(0.25517685) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53084757) q[0];
sx q[0];
rz(-1.8662211) q[0];
sx q[0];
rz(0.3770963) q[0];
rz(-0.43222506) q[2];
sx q[2];
rz(-0.6859633) q[2];
sx q[2];
rz(-1.9908817) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.42923388) q[1];
sx q[1];
rz(-2.7731967) q[1];
sx q[1];
rz(-0.952094) q[1];
x q[2];
rz(-0.051024036) q[3];
sx q[3];
rz(-2.0912366) q[3];
sx q[3];
rz(2.737962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.8884376) q[2];
sx q[2];
rz(-1.5506813) q[2];
sx q[2];
rz(0.17318428) q[2];
rz(2.611768) q[3];
sx q[3];
rz(-0.14557043) q[3];
sx q[3];
rz(0.1023275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.859905) q[0];
sx q[0];
rz(-1.6485933) q[0];
sx q[0];
rz(1.7657071) q[0];
rz(-1.2777404) q[1];
sx q[1];
rz(-2.3294096) q[1];
sx q[1];
rz(3.0854991) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2984021) q[0];
sx q[0];
rz(-2.139233) q[0];
sx q[0];
rz(1.1355023) q[0];
rz(-pi) q[1];
rz(2.202583) q[2];
sx q[2];
rz(-2.5917705) q[2];
sx q[2];
rz(0.75013559) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5406815) q[1];
sx q[1];
rz(-1.8735421) q[1];
sx q[1];
rz(1.5555698) q[1];
rz(3.0523473) q[3];
sx q[3];
rz(-1.0122932) q[3];
sx q[3];
rz(0.23111471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.00099480199) q[2];
sx q[2];
rz(-0.91789118) q[2];
sx q[2];
rz(-2.7094005) q[2];
rz(-2.2473992) q[3];
sx q[3];
rz(-2.0420572) q[3];
sx q[3];
rz(-1.4661219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0284001) q[0];
sx q[0];
rz(-0.88554651) q[0];
sx q[0];
rz(2.4940441) q[0];
rz(1.2619069) q[1];
sx q[1];
rz(-1.6779265) q[1];
sx q[1];
rz(-0.9544968) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4650824) q[0];
sx q[0];
rz(-1.0445147) q[0];
sx q[0];
rz(2.9617873) q[0];
rz(-1.5231832) q[2];
sx q[2];
rz(-0.63112586) q[2];
sx q[2];
rz(-0.39436755) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2973605) q[1];
sx q[1];
rz(-0.73207049) q[1];
sx q[1];
rz(-0.75433235) q[1];
rz(-2.0733842) q[3];
sx q[3];
rz(-2.8121901) q[3];
sx q[3];
rz(2.9941032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.548617) q[2];
sx q[2];
rz(-1.9081215) q[2];
sx q[2];
rz(2.0992289) q[2];
rz(0.43867612) q[3];
sx q[3];
rz(-1.0500267) q[3];
sx q[3];
rz(1.3180805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8577268) q[0];
sx q[0];
rz(-0.23290578) q[0];
sx q[0];
rz(-2.3983811) q[0];
rz(1.5076393) q[1];
sx q[1];
rz(-2.4217024) q[1];
sx q[1];
rz(2.5315703) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9305206) q[0];
sx q[0];
rz(-2.0606344) q[0];
sx q[0];
rz(1.2852438) q[0];
x q[1];
rz(-2.2248613) q[2];
sx q[2];
rz(-1.6628633) q[2];
sx q[2];
rz(-2.9031861) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.66079084) q[1];
sx q[1];
rz(-1.5457488) q[1];
sx q[1];
rz(0.90344306) q[1];
rz(-pi) q[2];
rz(1.1915019) q[3];
sx q[3];
rz(-1.3790352) q[3];
sx q[3];
rz(0.8134884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.33621776) q[2];
sx q[2];
rz(-1.6992133) q[2];
sx q[2];
rz(-2.2231893) q[2];
rz(-1.5504799) q[3];
sx q[3];
rz(-2.1912626) q[3];
sx q[3];
rz(-0.38890719) q[3];
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
rz(2.780705) q[0];
sx q[0];
rz(-2.4724859) q[0];
sx q[0];
rz(-1.6280744) q[0];
rz(-2.6121415) q[1];
sx q[1];
rz(-2.0748731) q[1];
sx q[1];
rz(0.73658529) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.026222762) q[0];
sx q[0];
rz(-1.5407469) q[0];
sx q[0];
rz(3.1307334) q[0];
rz(-3.0326764) q[2];
sx q[2];
rz(-1.8913336) q[2];
sx q[2];
rz(-0.045217302) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8229586) q[1];
sx q[1];
rz(-1.6726603) q[1];
sx q[1];
rz(2.676079) q[1];
x q[2];
rz(-1.1235808) q[3];
sx q[3];
rz(-1.2774602) q[3];
sx q[3];
rz(2.6597027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0344051) q[2];
sx q[2];
rz(-1.9341058) q[2];
sx q[2];
rz(2.4576808) q[2];
rz(1.2290139) q[3];
sx q[3];
rz(-1.3701655) q[3];
sx q[3];
rz(1.7470523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1972315) q[0];
sx q[0];
rz(-1.5427538) q[0];
sx q[0];
rz(2.9558682) q[0];
rz(2.1445403) q[1];
sx q[1];
rz(-1.2652218) q[1];
sx q[1];
rz(0.7448147) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54430994) q[0];
sx q[0];
rz(-0.83680698) q[0];
sx q[0];
rz(0.44041667) q[0];
rz(1.1698193) q[2];
sx q[2];
rz(-2.3224761) q[2];
sx q[2];
rz(-0.70105201) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8150755) q[1];
sx q[1];
rz(-0.41985598) q[1];
sx q[1];
rz(-0.99680568) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3585988) q[3];
sx q[3];
rz(-1.2851614) q[3];
sx q[3];
rz(-1.0234969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1054489) q[2];
sx q[2];
rz(-0.75573409) q[2];
sx q[2];
rz(-1.194681) q[2];
rz(0.99669325) q[3];
sx q[3];
rz(-1.9255305) q[3];
sx q[3];
rz(2.1452346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3982518) q[0];
sx q[0];
rz(-2.353459) q[0];
sx q[0];
rz(-2.7375896) q[0];
rz(-3.1104654) q[1];
sx q[1];
rz(-1.6571836) q[1];
sx q[1];
rz(1.9706479) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6452713) q[0];
sx q[0];
rz(-1.5491345) q[0];
sx q[0];
rz(0.42692703) q[0];
rz(-2.3659336) q[2];
sx q[2];
rz(-1.3752898) q[2];
sx q[2];
rz(-2.069371) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.78632894) q[1];
sx q[1];
rz(-1.5621645) q[1];
sx q[1];
rz(-1.5481871) q[1];
rz(-pi) q[2];
rz(-1.8627432) q[3];
sx q[3];
rz(-1.1598831) q[3];
sx q[3];
rz(-1.7850072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4460454) q[2];
sx q[2];
rz(-1.7798767) q[2];
sx q[2];
rz(2.5496303) q[2];
rz(2.5752318) q[3];
sx q[3];
rz(-0.16470328) q[3];
sx q[3];
rz(-1.5238354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3175209) q[0];
sx q[0];
rz(-0.98012797) q[0];
sx q[0];
rz(-1.160887) q[0];
rz(-3.042165) q[1];
sx q[1];
rz(-1.8933404) q[1];
sx q[1];
rz(1.0642687) q[1];
rz(-2.4026985) q[2];
sx q[2];
rz(-1.0514435) q[2];
sx q[2];
rz(-2.4098868) q[2];
rz(3.1254461) q[3];
sx q[3];
rz(-1.2373677) q[3];
sx q[3];
rz(2.2131372) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
