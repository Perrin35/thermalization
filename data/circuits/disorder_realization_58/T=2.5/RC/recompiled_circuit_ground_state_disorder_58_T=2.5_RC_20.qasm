OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6405606) q[0];
sx q[0];
rz(-0.80067331) q[0];
sx q[0];
rz(-0.14525695) q[0];
rz(2.2486806) q[1];
sx q[1];
rz(-0.4141663) q[1];
sx q[1];
rz(2.034634) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.993385) q[0];
sx q[0];
rz(-2.275106) q[0];
sx q[0];
rz(-0.765245) q[0];
rz(-pi) q[1];
rz(0.28319226) q[2];
sx q[2];
rz(-1.257466) q[2];
sx q[2];
rz(-1.4669795) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2295215) q[1];
sx q[1];
rz(-2.1451277) q[1];
sx q[1];
rz(-0.34194754) q[1];
rz(-0.21671076) q[3];
sx q[3];
rz(-1.1997077) q[3];
sx q[3];
rz(0.75066523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.60078159) q[2];
sx q[2];
rz(-1.3593707) q[2];
sx q[2];
rz(-3.0493128) q[2];
rz(1.7021092) q[3];
sx q[3];
rz(-1.9974134) q[3];
sx q[3];
rz(-2.2030742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0094902078) q[0];
sx q[0];
rz(-1.9785896) q[0];
sx q[0];
rz(-2.5376885) q[0];
rz(-1.5314792) q[1];
sx q[1];
rz(-1.8096626) q[1];
sx q[1];
rz(-1.6139222) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22900362) q[0];
sx q[0];
rz(-0.67978501) q[0];
sx q[0];
rz(2.0849821) q[0];
rz(-pi) q[1];
rz(-2.3983922) q[2];
sx q[2];
rz(-1.1172406) q[2];
sx q[2];
rz(1.1519866) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.81948796) q[1];
sx q[1];
rz(-2.8675962) q[1];
sx q[1];
rz(-3.0986087) q[1];
rz(-1.5714721) q[3];
sx q[3];
rz(-1.6689894) q[3];
sx q[3];
rz(-0.46969603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4332726) q[2];
sx q[2];
rz(-1.4488139) q[2];
sx q[2];
rz(-2.0937008) q[2];
rz(0.6692872) q[3];
sx q[3];
rz(-0.71988121) q[3];
sx q[3];
rz(-1.1650813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79353756) q[0];
sx q[0];
rz(-0.78472835) q[0];
sx q[0];
rz(1.0070356) q[0];
rz(-2.8496565) q[1];
sx q[1];
rz(-2.597229) q[1];
sx q[1];
rz(2.4709591) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1517253) q[0];
sx q[0];
rz(-1.6581931) q[0];
sx q[0];
rz(-0.94957994) q[0];
rz(-pi) q[1];
x q[1];
rz(0.67523662) q[2];
sx q[2];
rz(-0.91885447) q[2];
sx q[2];
rz(2.3777131) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2082767) q[1];
sx q[1];
rz(-1.9066992) q[1];
sx q[1];
rz(0.64888727) q[1];
rz(-3.0257053) q[3];
sx q[3];
rz(-1.2440344) q[3];
sx q[3];
rz(0.22921697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2900419) q[2];
sx q[2];
rz(-0.85323492) q[2];
sx q[2];
rz(0.89149371) q[2];
rz(-2.0391035) q[3];
sx q[3];
rz(-0.99415556) q[3];
sx q[3];
rz(1.4055143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3201228) q[0];
sx q[0];
rz(-1.7550884) q[0];
sx q[0];
rz(2.084305) q[0];
rz(0.0013466324) q[1];
sx q[1];
rz(-0.82095447) q[1];
sx q[1];
rz(-2.3473306) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1201354) q[0];
sx q[0];
rz(-1.5891783) q[0];
sx q[0];
rz(-2.4796499) q[0];
rz(-pi) q[1];
rz(-2.4636872) q[2];
sx q[2];
rz(-2.6982582) q[2];
sx q[2];
rz(1.5405129) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9201607) q[1];
sx q[1];
rz(-2.2756612) q[1];
sx q[1];
rz(0.95734289) q[1];
rz(-pi) q[2];
rz(1.6861071) q[3];
sx q[3];
rz(-1.8666576) q[3];
sx q[3];
rz(1.4747696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0346251) q[2];
sx q[2];
rz(-0.62109533) q[2];
sx q[2];
rz(2.1194439) q[2];
rz(-1.8978097) q[3];
sx q[3];
rz(-0.94977489) q[3];
sx q[3];
rz(-1.7826084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6720402) q[0];
sx q[0];
rz(-1.7615027) q[0];
sx q[0];
rz(-0.45355466) q[0];
rz(1.0376616) q[1];
sx q[1];
rz(-2.0062168) q[1];
sx q[1];
rz(-2.3496148) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4290473) q[0];
sx q[0];
rz(-0.28158108) q[0];
sx q[0];
rz(1.7292119) q[0];
x q[1];
rz(-0.99292314) q[2];
sx q[2];
rz(-1.9067592) q[2];
sx q[2];
rz(-3.1022037) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.32159785) q[1];
sx q[1];
rz(-1.454108) q[1];
sx q[1];
rz(0.34081809) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4326454) q[3];
sx q[3];
rz(-1.589121) q[3];
sx q[3];
rz(1.5442755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3114634) q[2];
sx q[2];
rz(-2.4757803) q[2];
sx q[2];
rz(2.8361481) q[2];
rz(-2.9597802) q[3];
sx q[3];
rz(-1.612855) q[3];
sx q[3];
rz(0.73604933) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5859454) q[0];
sx q[0];
rz(-2.3190627) q[0];
sx q[0];
rz(0.12538759) q[0];
rz(1.5749982) q[1];
sx q[1];
rz(-1.4488723) q[1];
sx q[1];
rz(3.1094508) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18460546) q[0];
sx q[0];
rz(-1.3520387) q[0];
sx q[0];
rz(2.4167528) q[0];
x q[1];
rz(0.06216269) q[2];
sx q[2];
rz(-0.84033687) q[2];
sx q[2];
rz(-2.3671012) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.39592) q[1];
sx q[1];
rz(-0.88282871) q[1];
sx q[1];
rz(0.93979551) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4549082) q[3];
sx q[3];
rz(-2.3130696) q[3];
sx q[3];
rz(2.7477086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0142168) q[2];
sx q[2];
rz(-1.9275503) q[2];
sx q[2];
rz(1.1227013) q[2];
rz(3.0681916) q[3];
sx q[3];
rz(-0.95728907) q[3];
sx q[3];
rz(2.5286123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4323394) q[0];
sx q[0];
rz(-1.657635) q[0];
sx q[0];
rz(0.51026979) q[0];
rz(-2.7773652) q[1];
sx q[1];
rz(-0.42114708) q[1];
sx q[1];
rz(1.6417004) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.071935805) q[0];
sx q[0];
rz(-0.47635117) q[0];
sx q[0];
rz(1.3385962) q[0];
x q[1];
rz(-1.4128311) q[2];
sx q[2];
rz(-1.4539495) q[2];
sx q[2];
rz(2.771488) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3245974) q[1];
sx q[1];
rz(-2.5909068) q[1];
sx q[1];
rz(-0.43795069) q[1];
rz(-pi) q[2];
rz(2.4046201) q[3];
sx q[3];
rz(-1.342257) q[3];
sx q[3];
rz(1.6933954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4437272) q[2];
sx q[2];
rz(-1.5865734) q[2];
sx q[2];
rz(0.86722428) q[2];
rz(1.5602268) q[3];
sx q[3];
rz(-2.9428704) q[3];
sx q[3];
rz(1.8038512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7241868) q[0];
sx q[0];
rz(-3.0751808) q[0];
sx q[0];
rz(-0.41626406) q[0];
rz(1.9789713) q[1];
sx q[1];
rz(-1.9527718) q[1];
sx q[1];
rz(2.3847041) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9577741) q[0];
sx q[0];
rz(-1.2831492) q[0];
sx q[0];
rz(2.0145061) q[0];
x q[1];
rz(2.1483118) q[2];
sx q[2];
rz(-1.9406291) q[2];
sx q[2];
rz(-3.0209857) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6357506) q[1];
sx q[1];
rz(-1.381784) q[1];
sx q[1];
rz(2.0057949) q[1];
rz(-2.5439436) q[3];
sx q[3];
rz(-0.72859287) q[3];
sx q[3];
rz(0.33580175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.68652144) q[2];
sx q[2];
rz(-1.4342118) q[2];
sx q[2];
rz(-1.7835468) q[2];
rz(0.81651917) q[3];
sx q[3];
rz(-1.9101382) q[3];
sx q[3];
rz(0.16286287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4002976) q[0];
sx q[0];
rz(-1.4485899) q[0];
sx q[0];
rz(-1.4073538) q[0];
rz(-0.74527144) q[1];
sx q[1];
rz(-2.407357) q[1];
sx q[1];
rz(-1.5819246) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9055515) q[0];
sx q[0];
rz(-2.2187382) q[0];
sx q[0];
rz(-2.7929162) q[0];
rz(3.09868) q[2];
sx q[2];
rz(-2.4933698) q[2];
sx q[2];
rz(0.15951482) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0631536) q[1];
sx q[1];
rz(-0.4635103) q[1];
sx q[1];
rz(-1.2786675) q[1];
rz(-pi) q[2];
rz(0.99509434) q[3];
sx q[3];
rz(-1.5202209) q[3];
sx q[3];
rz(-2.3114396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.72426307) q[2];
sx q[2];
rz(-1.686692) q[2];
sx q[2];
rz(-1.1154741) q[2];
rz(-0.25350246) q[3];
sx q[3];
rz(-1.8799672) q[3];
sx q[3];
rz(-3.1326262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-2.1215006) q[0];
sx q[0];
rz(-1.8778863) q[0];
sx q[0];
rz(0.76250917) q[0];
rz(-2.1992042) q[1];
sx q[1];
rz(-1.0647048) q[1];
sx q[1];
rz(1.2368088) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.002703) q[0];
sx q[0];
rz(-1.5093439) q[0];
sx q[0];
rz(-0.49646722) q[0];
x q[1];
rz(0.89434172) q[2];
sx q[2];
rz(-0.80781762) q[2];
sx q[2];
rz(-1.3685014) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.33442228) q[1];
sx q[1];
rz(-1.5457193) q[1];
sx q[1];
rz(-1.5958171) q[1];
rz(-pi) q[2];
rz(-0.071466515) q[3];
sx q[3];
rz(-1.1671007) q[3];
sx q[3];
rz(1.57406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8952055) q[2];
sx q[2];
rz(-3.0463986) q[2];
sx q[2];
rz(0.0072366317) q[2];
rz(-0.17035189) q[3];
sx q[3];
rz(-1.6536313) q[3];
sx q[3];
rz(0.65792221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-1.2687179) q[0];
sx q[0];
rz(-1.7418516) q[0];
sx q[0];
rz(-2.0280784) q[0];
rz(0.089182236) q[1];
sx q[1];
rz(-1.5166278) q[1];
sx q[1];
rz(1.2022432) q[1];
rz(-2.9847445) q[2];
sx q[2];
rz(-1.4365755) q[2];
sx q[2];
rz(2.8683521) q[2];
rz(-2.8409957) q[3];
sx q[3];
rz(-1.7491805) q[3];
sx q[3];
rz(0.94284369) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
