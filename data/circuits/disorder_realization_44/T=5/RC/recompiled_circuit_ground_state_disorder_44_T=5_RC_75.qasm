OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3888336) q[0];
sx q[0];
rz(-1.692481) q[0];
sx q[0];
rz(-1.4504855) q[0];
rz(0.9736355) q[1];
sx q[1];
rz(-1.7042301) q[1];
sx q[1];
rz(-0.91926423) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96199233) q[0];
sx q[0];
rz(-1.5536947) q[0];
sx q[0];
rz(-0.03058612) q[0];
rz(0.12030258) q[2];
sx q[2];
rz(-1.4805111) q[2];
sx q[2];
rz(1.1995969) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.011292) q[1];
sx q[1];
rz(-2.4280972) q[1];
sx q[1];
rz(1.43047) q[1];
rz(-pi) q[2];
rz(-2.0415061) q[3];
sx q[3];
rz(-1.19095) q[3];
sx q[3];
rz(-2.0947411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3322525) q[2];
sx q[2];
rz(-1.6925749) q[2];
sx q[2];
rz(-1.4130886) q[2];
rz(0.20279065) q[3];
sx q[3];
rz(-1.7594124) q[3];
sx q[3];
rz(0.10281674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8973812) q[0];
sx q[0];
rz(-2.057071) q[0];
sx q[0];
rz(-0.71075034) q[0];
rz(-2.3731025) q[1];
sx q[1];
rz(-1.0667421) q[1];
sx q[1];
rz(1.01952) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9907889) q[0];
sx q[0];
rz(-1.5312563) q[0];
sx q[0];
rz(-0.44724748) q[0];
rz(-pi) q[1];
rz(1.7350082) q[2];
sx q[2];
rz(-2.7052042) q[2];
sx q[2];
rz(-1.1498888) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.92495698) q[1];
sx q[1];
rz(-1.5397738) q[1];
sx q[1];
rz(-1.9418632) q[1];
x q[2];
rz(-1.4367661) q[3];
sx q[3];
rz(-0.55826) q[3];
sx q[3];
rz(-2.9001065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3780313) q[2];
sx q[2];
rz(-1.2445933) q[2];
sx q[2];
rz(-1.9192609) q[2];
rz(1.9289121) q[3];
sx q[3];
rz(-2.1247532) q[3];
sx q[3];
rz(-0.71162629) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0843622) q[0];
sx q[0];
rz(-2.1570692) q[0];
sx q[0];
rz(-1.2179751) q[0];
rz(-0.51586622) q[1];
sx q[1];
rz(-0.54793826) q[1];
sx q[1];
rz(0.9224433) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24127125) q[0];
sx q[0];
rz(-1.6113971) q[0];
sx q[0];
rz(0.057828219) q[0];
rz(-2.416196) q[2];
sx q[2];
rz(-1.9624406) q[2];
sx q[2];
rz(1.894941) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6906185) q[1];
sx q[1];
rz(-0.96267525) q[1];
sx q[1];
rz(-2.4324424) q[1];
x q[2];
rz(0.32914583) q[3];
sx q[3];
rz(-1.7308047) q[3];
sx q[3];
rz(3.126006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1232274) q[2];
sx q[2];
rz(-2.0570698) q[2];
sx q[2];
rz(1.8017192) q[2];
rz(-2.5943622) q[3];
sx q[3];
rz(-2.7180143) q[3];
sx q[3];
rz(-2.4303998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0860586) q[0];
sx q[0];
rz(-2.9965897) q[0];
sx q[0];
rz(-2.9823629) q[0];
rz(3.1314462) q[1];
sx q[1];
rz(-1.0228913) q[1];
sx q[1];
rz(-2.8841282) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87436324) q[0];
sx q[0];
rz(-2.3152475) q[0];
sx q[0];
rz(2.6494725) q[0];
rz(-0.065071062) q[2];
sx q[2];
rz(-0.79539585) q[2];
sx q[2];
rz(0.8000904) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.60002335) q[1];
sx q[1];
rz(-0.63508717) q[1];
sx q[1];
rz(-2.9734008) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.61473989) q[3];
sx q[3];
rz(-0.90723824) q[3];
sx q[3];
rz(-3.1082982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.7633535) q[2];
sx q[2];
rz(-1.8469609) q[2];
sx q[2];
rz(-2.6687458) q[2];
rz(0.7044479) q[3];
sx q[3];
rz(-1.7770146) q[3];
sx q[3];
rz(0.64594597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6351629) q[0];
sx q[0];
rz(-1.1744873) q[0];
sx q[0];
rz(-0.065486431) q[0];
rz(-0.40924117) q[1];
sx q[1];
rz(-2.0114653) q[1];
sx q[1];
rz(1.7154891) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1714892) q[0];
sx q[0];
rz(-1.3521776) q[0];
sx q[0];
rz(-1.4755558) q[0];
rz(1.481856) q[2];
sx q[2];
rz(-1.1834025) q[2];
sx q[2];
rz(2.873444) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3196484) q[1];
sx q[1];
rz(-1.8488171) q[1];
sx q[1];
rz(-1.8026428) q[1];
rz(-pi) q[2];
rz(1.8614011) q[3];
sx q[3];
rz(-1.7365321) q[3];
sx q[3];
rz(-0.7208212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9466729) q[2];
sx q[2];
rz(-2.0544923) q[2];
sx q[2];
rz(-1.8348414) q[2];
rz(1.9715747) q[3];
sx q[3];
rz(-2.7368059) q[3];
sx q[3];
rz(-3.034333) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65866798) q[0];
sx q[0];
rz(-1.1558477) q[0];
sx q[0];
rz(-2.7401127) q[0];
rz(2.4688156) q[1];
sx q[1];
rz(-0.79877001) q[1];
sx q[1];
rz(-0.85404095) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14388785) q[0];
sx q[0];
rz(-1.524462) q[0];
sx q[0];
rz(1.7884939) q[0];
x q[1];
rz(-2.9561192) q[2];
sx q[2];
rz(-1.4234241) q[2];
sx q[2];
rz(2.3285248) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2931765) q[1];
sx q[1];
rz(-0.59976116) q[1];
sx q[1];
rz(2.9648215) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5402732) q[3];
sx q[3];
rz(-0.70402217) q[3];
sx q[3];
rz(-2.6773334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.65471571) q[2];
sx q[2];
rz(-2.2701023) q[2];
sx q[2];
rz(1.1775449) q[2];
rz(0.55142895) q[3];
sx q[3];
rz(-1.5827554) q[3];
sx q[3];
rz(-2.3059755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.7198782) q[0];
sx q[0];
rz(-0.28110176) q[0];
sx q[0];
rz(1.9792492) q[0];
rz(1.989919) q[1];
sx q[1];
rz(-1.3538227) q[1];
sx q[1];
rz(-0.91845671) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29354039) q[0];
sx q[0];
rz(-1.9013202) q[0];
sx q[0];
rz(-0.28384112) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2662042) q[2];
sx q[2];
rz(-1.3686485) q[2];
sx q[2];
rz(-0.018176807) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9388401) q[1];
sx q[1];
rz(-1.5781286) q[1];
sx q[1];
rz(1.5640902) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0189459) q[3];
sx q[3];
rz(-2.7635241) q[3];
sx q[3];
rz(-2.700875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.034417001) q[2];
sx q[2];
rz(-1.7566661) q[2];
sx q[2];
rz(2.6386063) q[2];
rz(-0.77477396) q[3];
sx q[3];
rz(-2.6796902) q[3];
sx q[3];
rz(-1.7983961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69304943) q[0];
sx q[0];
rz(-2.9149084) q[0];
sx q[0];
rz(-1.6402798) q[0];
rz(-0.34128183) q[1];
sx q[1];
rz(-1.7106067) q[1];
sx q[1];
rz(-2.0223845) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8728031) q[0];
sx q[0];
rz(-2.4478292) q[0];
sx q[0];
rz(-0.60141464) q[0];
rz(-pi) q[1];
x q[1];
rz(2.881024) q[2];
sx q[2];
rz(-2.0379279) q[2];
sx q[2];
rz(-1.1819983) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.70485605) q[1];
sx q[1];
rz(-2.4697127) q[1];
sx q[1];
rz(-0.19499548) q[1];
rz(1.7057034) q[3];
sx q[3];
rz(-1.9622318) q[3];
sx q[3];
rz(-1.5086482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1559653) q[2];
sx q[2];
rz(-2.1817744) q[2];
sx q[2];
rz(-2.7723374) q[2];
rz(0.30820942) q[3];
sx q[3];
rz(-0.9674558) q[3];
sx q[3];
rz(-1.6109899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2823328) q[0];
sx q[0];
rz(-1.8349324) q[0];
sx q[0];
rz(0.34307137) q[0];
rz(-1.9725017) q[1];
sx q[1];
rz(-1.7770551) q[1];
sx q[1];
rz(-1.7128568) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.13641) q[0];
sx q[0];
rz(-0.47629582) q[0];
sx q[0];
rz(-2.5695557) q[0];
x q[1];
rz(-2.3071204) q[2];
sx q[2];
rz(-1.021046) q[2];
sx q[2];
rz(-2.2817734) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0020121) q[1];
sx q[1];
rz(-1.5095469) q[1];
sx q[1];
rz(1.9467926) q[1];
rz(-1.021528) q[3];
sx q[3];
rz(-2.3695951) q[3];
sx q[3];
rz(0.87024161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2009361) q[2];
sx q[2];
rz(-0.20874615) q[2];
sx q[2];
rz(-1.3915871) q[2];
rz(-0.40677795) q[3];
sx q[3];
rz(-1.4429561) q[3];
sx q[3];
rz(0.87882915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10130356) q[0];
sx q[0];
rz(-0.60281301) q[0];
sx q[0];
rz(-1.3579177) q[0];
rz(1.2376002) q[1];
sx q[1];
rz(-2.1289181) q[1];
sx q[1];
rz(0.79992574) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42521503) q[0];
sx q[0];
rz(-2.3122462) q[0];
sx q[0];
rz(-0.31405507) q[0];
x q[1];
rz(-1.8736035) q[2];
sx q[2];
rz(-1.0410415) q[2];
sx q[2];
rz(-2.985266) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7042735) q[1];
sx q[1];
rz(-1.4156439) q[1];
sx q[1];
rz(-0.67351933) q[1];
rz(-pi) q[2];
x q[2];
rz(0.25000817) q[3];
sx q[3];
rz(-2.0431314) q[3];
sx q[3];
rz(1.3046622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.37821975) q[2];
sx q[2];
rz(-1.7174481) q[2];
sx q[2];
rz(0.073089449) q[2];
rz(2.8857005) q[3];
sx q[3];
rz(-0.81232324) q[3];
sx q[3];
rz(1.6528486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-1.273461) q[0];
sx q[0];
rz(-1.4556226) q[0];
sx q[0];
rz(1.8709394) q[0];
rz(1.3399667) q[1];
sx q[1];
rz(-1.3506964) q[1];
sx q[1];
rz(0.66257308) q[1];
rz(-1.5727829) q[2];
sx q[2];
rz(-2.28021) q[2];
sx q[2];
rz(1.8328666) q[2];
rz(-2.4347435) q[3];
sx q[3];
rz(-2.3656188) q[3];
sx q[3];
rz(-1.594365) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
