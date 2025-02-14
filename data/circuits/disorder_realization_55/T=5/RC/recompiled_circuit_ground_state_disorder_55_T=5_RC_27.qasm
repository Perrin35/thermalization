OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2321229) q[0];
sx q[0];
rz(-0.98804086) q[0];
sx q[0];
rz(0.35257941) q[0];
rz(-0.70631385) q[1];
sx q[1];
rz(-0.7847473) q[1];
sx q[1];
rz(0.027332505) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77307149) q[0];
sx q[0];
rz(-1.7681244) q[0];
sx q[0];
rz(0.34223603) q[0];
x q[1];
rz(2.4800221) q[2];
sx q[2];
rz(-2.9230766) q[2];
sx q[2];
rz(1.0929581) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.057291659) q[1];
sx q[1];
rz(-2.9708869) q[1];
sx q[1];
rz(0.18526669) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7165511) q[3];
sx q[3];
rz(-1.3363239) q[3];
sx q[3];
rz(-2.2343079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0693543) q[2];
sx q[2];
rz(-1.8202929) q[2];
sx q[2];
rz(-1.5864774) q[2];
rz(-2.41411) q[3];
sx q[3];
rz(-0.95557094) q[3];
sx q[3];
rz(2.6684842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.469406) q[0];
sx q[0];
rz(-1.3453901) q[0];
sx q[0];
rz(2.192002) q[0];
rz(-2.8166215) q[1];
sx q[1];
rz(-1.800671) q[1];
sx q[1];
rz(0.51345888) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9531813) q[0];
sx q[0];
rz(-1.9741749) q[0];
sx q[0];
rz(0.90103407) q[0];
rz(-0.84324381) q[2];
sx q[2];
rz(-1.8157037) q[2];
sx q[2];
rz(-1.5050735) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4411575) q[1];
sx q[1];
rz(-1.0714021) q[1];
sx q[1];
rz(-2.5446289) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1517215) q[3];
sx q[3];
rz(-1.9912915) q[3];
sx q[3];
rz(-2.6591501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.1309588) q[2];
sx q[2];
rz(-1.3291239) q[2];
sx q[2];
rz(2.0531674) q[2];
rz(2.4736577) q[3];
sx q[3];
rz(-0.1440983) q[3];
sx q[3];
rz(1.4151423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82069355) q[0];
sx q[0];
rz(-0.90105337) q[0];
sx q[0];
rz(1.8841085) q[0];
rz(3.056774) q[1];
sx q[1];
rz(-0.091400472) q[1];
sx q[1];
rz(2.4841097) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.062825052) q[0];
sx q[0];
rz(-1.4402534) q[0];
sx q[0];
rz(0.45147942) q[0];
x q[1];
rz(2.5490513) q[2];
sx q[2];
rz(-2.0503902) q[2];
sx q[2];
rz(-1.3035989) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.3448101) q[1];
sx q[1];
rz(-1.0940897) q[1];
sx q[1];
rz(-1.3784798) q[1];
x q[2];
rz(-0.90963962) q[3];
sx q[3];
rz(-1.8170905) q[3];
sx q[3];
rz(-1.9734188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.76454863) q[2];
sx q[2];
rz(-2.720764) q[2];
sx q[2];
rz(1.8910889) q[2];
rz(1.4416384) q[3];
sx q[3];
rz(-1.7504033) q[3];
sx q[3];
rz(2.3120248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4568951) q[0];
sx q[0];
rz(-2.1853515) q[0];
sx q[0];
rz(-2.5392505) q[0];
rz(-2.1777228) q[1];
sx q[1];
rz(-2.3606221) q[1];
sx q[1];
rz(-2.9197555) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5103127) q[0];
sx q[0];
rz(-1.5595734) q[0];
sx q[0];
rz(1.5785286) q[0];
x q[1];
rz(1.0413136) q[2];
sx q[2];
rz(-1.3544677) q[2];
sx q[2];
rz(-1.554503) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.79961046) q[1];
sx q[1];
rz(-2.1253808) q[1];
sx q[1];
rz(2.0444524) q[1];
rz(-pi) q[2];
rz(-0.45000123) q[3];
sx q[3];
rz(-1.6353938) q[3];
sx q[3];
rz(-2.3774176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7092789) q[2];
sx q[2];
rz(-1.1608492) q[2];
sx q[2];
rz(3.0529037) q[2];
rz(2.6426219) q[3];
sx q[3];
rz(-1.5191398) q[3];
sx q[3];
rz(-3.0626845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13224193) q[0];
sx q[0];
rz(-3.0920588) q[0];
sx q[0];
rz(0.90231878) q[0];
rz(0.29014507) q[1];
sx q[1];
rz(-2.2369308) q[1];
sx q[1];
rz(-1.1660928) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6571044) q[0];
sx q[0];
rz(-1.00923) q[0];
sx q[0];
rz(-1.3795555) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8639216) q[2];
sx q[2];
rz(-1.8050965) q[2];
sx q[2];
rz(-1.9512343) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2578967) q[1];
sx q[1];
rz(-1.2295112) q[1];
sx q[1];
rz(-2.8410556) q[1];
rz(-pi) q[2];
rz(-2.8899447) q[3];
sx q[3];
rz(-1.5960428) q[3];
sx q[3];
rz(-1.7612875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.12639283) q[2];
sx q[2];
rz(-2.3455399) q[2];
sx q[2];
rz(-1.8533206) q[2];
rz(1.0334233) q[3];
sx q[3];
rz(-1.1181592) q[3];
sx q[3];
rz(-1.4748632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7159336) q[0];
sx q[0];
rz(-3.0379744) q[0];
sx q[0];
rz(2.9820251) q[0];
rz(0.60699925) q[1];
sx q[1];
rz(-1.9885352) q[1];
sx q[1];
rz(2.0807696) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88638611) q[0];
sx q[0];
rz(-1.7632503) q[0];
sx q[0];
rz(-0.9631971) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.873704) q[2];
sx q[2];
rz(-1.3527591) q[2];
sx q[2];
rz(-1.5341369) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.80172864) q[1];
sx q[1];
rz(-0.9802981) q[1];
sx q[1];
rz(-2.9782359) q[1];
rz(-pi) q[2];
rz(0.14935519) q[3];
sx q[3];
rz(-2.5575482) q[3];
sx q[3];
rz(-2.4660939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1242421) q[2];
sx q[2];
rz(-2.3257747) q[2];
sx q[2];
rz(-2.7537277) q[2];
rz(1.806949) q[3];
sx q[3];
rz(-2.4692061) q[3];
sx q[3];
rz(1.0268802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89011985) q[0];
sx q[0];
rz(-2.1595182) q[0];
sx q[0];
rz(-0.89333308) q[0];
rz(-0.54126254) q[1];
sx q[1];
rz(-0.90507871) q[1];
sx q[1];
rz(1.062324) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8996595) q[0];
sx q[0];
rz(-0.77882732) q[0];
sx q[0];
rz(-0.3758442) q[0];
x q[1];
rz(-1.6387894) q[2];
sx q[2];
rz(-0.74882245) q[2];
sx q[2];
rz(2.6523726) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.17328611) q[1];
sx q[1];
rz(-1.5472551) q[1];
sx q[1];
rz(0.55577718) q[1];
rz(-pi) q[2];
x q[2];
rz(0.19936187) q[3];
sx q[3];
rz(-2.3219675) q[3];
sx q[3];
rz(0.79341753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5541151) q[2];
sx q[2];
rz(-2.1671961) q[2];
sx q[2];
rz(1.1678196) q[2];
rz(0.81985146) q[3];
sx q[3];
rz(-0.76698118) q[3];
sx q[3];
rz(-0.64275536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.576936) q[0];
sx q[0];
rz(-2.7695203) q[0];
sx q[0];
rz(1.8336953) q[0];
rz(-1.4564184) q[1];
sx q[1];
rz(-1.6273727) q[1];
sx q[1];
rz(-3.0110722) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3878138) q[0];
sx q[0];
rz(-0.31720933) q[0];
sx q[0];
rz(1.8075726) q[0];
rz(-pi) q[1];
rz(-2.0575664) q[2];
sx q[2];
rz(-2.3064838) q[2];
sx q[2];
rz(1.4722919) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.83090342) q[1];
sx q[1];
rz(-2.2210418) q[1];
sx q[1];
rz(1.5216773) q[1];
x q[2];
rz(2.42584) q[3];
sx q[3];
rz(-2.0945315) q[3];
sx q[3];
rz(-2.2534086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.1859583) q[2];
sx q[2];
rz(-2.3027577) q[2];
sx q[2];
rz(-1.4107417) q[2];
rz(-3.020982) q[3];
sx q[3];
rz(-1.4863622) q[3];
sx q[3];
rz(1.6284778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2738344) q[0];
sx q[0];
rz(-0.61536106) q[0];
sx q[0];
rz(-0.14061418) q[0];
rz(-0.73453844) q[1];
sx q[1];
rz(-1.7081552) q[1];
sx q[1];
rz(-0.75070423) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0854617) q[0];
sx q[0];
rz(-1.805725) q[0];
sx q[0];
rz(-0.30130193) q[0];
rz(-pi) q[1];
rz(-1.2241988) q[2];
sx q[2];
rz(-1.3470413) q[2];
sx q[2];
rz(0.72347298) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0929326) q[1];
sx q[1];
rz(-1.5092998) q[1];
sx q[1];
rz(1.7503121) q[1];
rz(-pi) q[2];
rz(3.0680411) q[3];
sx q[3];
rz(-0.86553177) q[3];
sx q[3];
rz(1.2770231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.004403) q[2];
sx q[2];
rz(-0.69978324) q[2];
sx q[2];
rz(2.3889551) q[2];
rz(-2.802012) q[3];
sx q[3];
rz(-1.7759674) q[3];
sx q[3];
rz(3.1201709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4362157) q[0];
sx q[0];
rz(-2.3973873) q[0];
sx q[0];
rz(-0.63802737) q[0];
rz(-0.72884196) q[1];
sx q[1];
rz(-1.809027) q[1];
sx q[1];
rz(-0.80150882) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34018286) q[0];
sx q[0];
rz(-1.7728285) q[0];
sx q[0];
rz(-0.8043181) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3036868) q[2];
sx q[2];
rz(-1.3896754) q[2];
sx q[2];
rz(1.6194026) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2406225) q[1];
sx q[1];
rz(-0.83604807) q[1];
sx q[1];
rz(3.088984) q[1];
rz(-pi) q[2];
rz(-0.92592327) q[3];
sx q[3];
rz(-1.7002673) q[3];
sx q[3];
rz(-2.7162939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.30404299) q[2];
sx q[2];
rz(-1.7767228) q[2];
sx q[2];
rz(-1.1207646) q[2];
rz(-2.4325727) q[3];
sx q[3];
rz(-2.6778383) q[3];
sx q[3];
rz(1.2254865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7365702) q[0];
sx q[0];
rz(-0.46157349) q[0];
sx q[0];
rz(2.9622958) q[0];
rz(0.90514056) q[1];
sx q[1];
rz(-2.1642579) q[1];
sx q[1];
rz(-0.49107818) q[1];
rz(2.7837743) q[2];
sx q[2];
rz(-0.87643788) q[2];
sx q[2];
rz(-1.4446148) q[2];
rz(0.2855339) q[3];
sx q[3];
rz(-1.0665421) q[3];
sx q[3];
rz(1.4530696) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
