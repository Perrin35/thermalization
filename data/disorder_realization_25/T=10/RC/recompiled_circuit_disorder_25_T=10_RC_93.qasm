OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.0619573) q[0];
sx q[0];
rz(-0.24793967) q[0];
sx q[0];
rz(2.3214582) q[0];
rz(-3.1007383) q[1];
sx q[1];
rz(-0.78894579) q[1];
sx q[1];
rz(3.0537002) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9633006) q[0];
sx q[0];
rz(-2.0110197) q[0];
sx q[0];
rz(-2.630169) q[0];
rz(-0.42638875) q[2];
sx q[2];
rz(-1.083056) q[2];
sx q[2];
rz(-2.3418155) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3808448) q[1];
sx q[1];
rz(-1.5729841) q[1];
sx q[1];
rz(-2.6373177) q[1];
rz(-pi) q[2];
rz(-1.8294931) q[3];
sx q[3];
rz(-1.5651349) q[3];
sx q[3];
rz(-2.2744846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1408046) q[2];
sx q[2];
rz(-1.452383) q[2];
sx q[2];
rz(-2.5722356) q[2];
rz(-1.6128444) q[3];
sx q[3];
rz(-2.5053535) q[3];
sx q[3];
rz(-1.3585842) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5434718) q[0];
sx q[0];
rz(-2.9765029) q[0];
sx q[0];
rz(-2.5879481) q[0];
rz(1.9042632) q[1];
sx q[1];
rz(-1.3668704) q[1];
sx q[1];
rz(1.233261) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.935826) q[0];
sx q[0];
rz(-0.92139771) q[0];
sx q[0];
rz(-2.4539024) q[0];
x q[1];
rz(-2.9821175) q[2];
sx q[2];
rz(-0.6808241) q[2];
sx q[2];
rz(-0.28123873) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.711449) q[1];
sx q[1];
rz(-1.5866382) q[1];
sx q[1];
rz(1.364691) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4741461) q[3];
sx q[3];
rz(-1.3818041) q[3];
sx q[3];
rz(-2.8611956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.1374986) q[2];
sx q[2];
rz(-1.6823021) q[2];
sx q[2];
rz(-0.0022350524) q[2];
rz(-0.8301174) q[3];
sx q[3];
rz(-0.79289645) q[3];
sx q[3];
rz(-1.8921651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24901351) q[0];
sx q[0];
rz(-0.61616388) q[0];
sx q[0];
rz(-2.2900443) q[0];
rz(-0.7710723) q[1];
sx q[1];
rz(-1.675019) q[1];
sx q[1];
rz(0.071203701) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.035633798) q[0];
sx q[0];
rz(-1.3741612) q[0];
sx q[0];
rz(1.3282302) q[0];
rz(-pi) q[1];
rz(-1.0686915) q[2];
sx q[2];
rz(-0.50902589) q[2];
sx q[2];
rz(-0.44391649) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.33176955) q[1];
sx q[1];
rz(-0.66546813) q[1];
sx q[1];
rz(-0.21824093) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2749388) q[3];
sx q[3];
rz(-0.94368499) q[3];
sx q[3];
rz(-1.31124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3123793) q[2];
sx q[2];
rz(-2.2177314) q[2];
sx q[2];
rz(-0.52345792) q[2];
rz(-0.33189014) q[3];
sx q[3];
rz(-0.060398014) q[3];
sx q[3];
rz(-2.6383242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4098542) q[0];
sx q[0];
rz(-2.8926352) q[0];
sx q[0];
rz(2.3160034) q[0];
rz(1.1666974) q[1];
sx q[1];
rz(-2.2749133) q[1];
sx q[1];
rz(1.694214) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71890812) q[0];
sx q[0];
rz(-1.3410765) q[0];
sx q[0];
rz(-2.9740888) q[0];
rz(3.1324603) q[2];
sx q[2];
rz(-1.7491241) q[2];
sx q[2];
rz(-0.4724617) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1737766) q[1];
sx q[1];
rz(-2.2602091) q[1];
sx q[1];
rz(-2.9348228) q[1];
x q[2];
rz(-0.36376603) q[3];
sx q[3];
rz(-1.3437265) q[3];
sx q[3];
rz(1.0432537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.46431413) q[2];
sx q[2];
rz(-1.7769287) q[2];
sx q[2];
rz(0.14492598) q[2];
rz(2.1302917) q[3];
sx q[3];
rz(-0.62711182) q[3];
sx q[3];
rz(3.1269126) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99575627) q[0];
sx q[0];
rz(-3.0881112) q[0];
sx q[0];
rz(-2.4575535) q[0];
rz(1.9513291) q[1];
sx q[1];
rz(-2.4763156) q[1];
sx q[1];
rz(-1.9285944) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1319259) q[0];
sx q[0];
rz(-1.3103232) q[0];
sx q[0];
rz(0.38513222) q[0];
rz(-2.2239387) q[2];
sx q[2];
rz(-2.0954164) q[2];
sx q[2];
rz(1.4337002) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.1274384) q[1];
sx q[1];
rz(-1.6886097) q[1];
sx q[1];
rz(-0.25728667) q[1];
rz(2.9641987) q[3];
sx q[3];
rz(-1.3131485) q[3];
sx q[3];
rz(0.28211668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.50903901) q[2];
sx q[2];
rz(-1.8703439) q[2];
sx q[2];
rz(0.63009134) q[2];
rz(0.57224327) q[3];
sx q[3];
rz(-2.4961491) q[3];
sx q[3];
rz(-2.2563289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1155788) q[0];
sx q[0];
rz(-2.2981839) q[0];
sx q[0];
rz(-2.8552326) q[0];
rz(2.5419366) q[1];
sx q[1];
rz(-2.0136132) q[1];
sx q[1];
rz(-0.58475959) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11767865) q[0];
sx q[0];
rz(-2.4524134) q[0];
sx q[0];
rz(-0.50424772) q[0];
rz(-pi) q[1];
x q[1];
rz(0.037126304) q[2];
sx q[2];
rz(-0.6147487) q[2];
sx q[2];
rz(0.45072134) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0271921) q[1];
sx q[1];
rz(-0.90585867) q[1];
sx q[1];
rz(2.7726638) q[1];
rz(-pi) q[2];
rz(2.9752475) q[3];
sx q[3];
rz(-1.5208941) q[3];
sx q[3];
rz(-0.14453105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9770603) q[2];
sx q[2];
rz(-0.88447276) q[2];
sx q[2];
rz(-1.5131081) q[2];
rz(-2.1785054) q[3];
sx q[3];
rz(-1.77308) q[3];
sx q[3];
rz(0.62455463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-3.0528089) q[0];
sx q[0];
rz(-1.7099021) q[0];
sx q[0];
rz(-0.54876304) q[0];
rz(0.051963003) q[1];
sx q[1];
rz(-0.49182645) q[1];
sx q[1];
rz(-2.5351977) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.711395) q[0];
sx q[0];
rz(-2.093962) q[0];
sx q[0];
rz(0.27336143) q[0];
rz(-1.1029926) q[2];
sx q[2];
rz(-0.71873795) q[2];
sx q[2];
rz(-1.4008092) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4013306) q[1];
sx q[1];
rz(-1.6832502) q[1];
sx q[1];
rz(-2.3349891) q[1];
rz(-pi) q[2];
rz(2.2783979) q[3];
sx q[3];
rz(-2.6226225) q[3];
sx q[3];
rz(-0.88245813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2018044) q[2];
sx q[2];
rz(-2.1494892) q[2];
sx q[2];
rz(0.50841224) q[2];
rz(2.0243747) q[3];
sx q[3];
rz(-0.17337392) q[3];
sx q[3];
rz(-1.6569051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8023119) q[0];
sx q[0];
rz(-0.34582129) q[0];
sx q[0];
rz(0.75396496) q[0];
rz(-0.94999653) q[1];
sx q[1];
rz(-1.4818622) q[1];
sx q[1];
rz(-0.22769134) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84275093) q[0];
sx q[0];
rz(-1.9831053) q[0];
sx q[0];
rz(-0.077667872) q[0];
rz(-pi) q[1];
rz(2.6954306) q[2];
sx q[2];
rz(-1.9036306) q[2];
sx q[2];
rz(-1.7762426) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.55351174) q[1];
sx q[1];
rz(-1.3462102) q[1];
sx q[1];
rz(-2.5976546) q[1];
rz(3.0759096) q[3];
sx q[3];
rz(-0.55995299) q[3];
sx q[3];
rz(-0.63510676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.078538744) q[2];
sx q[2];
rz(-0.20919122) q[2];
sx q[2];
rz(-0.83855808) q[2];
rz(2.7770384) q[3];
sx q[3];
rz(-1.3804599) q[3];
sx q[3];
rz(-0.84028876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7446328) q[0];
sx q[0];
rz(-1.7663706) q[0];
sx q[0];
rz(2.3378085) q[0];
rz(1.0546168) q[1];
sx q[1];
rz(-0.63916731) q[1];
sx q[1];
rz(-1.9931591) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6156857) q[0];
sx q[0];
rz(-1.4988168) q[0];
sx q[0];
rz(-3.1176438) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.9903332) q[2];
sx q[2];
rz(-0.60539421) q[2];
sx q[2];
rz(-2.6661172) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0828404) q[1];
sx q[1];
rz(-2.0420923) q[1];
sx q[1];
rz(2.6575762) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6059444) q[3];
sx q[3];
rz(-1.2190378) q[3];
sx q[3];
rz(2.7621321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0917197) q[2];
sx q[2];
rz(-2.1240978) q[2];
sx q[2];
rz(-0.77825528) q[2];
rz(-2.9459279) q[3];
sx q[3];
rz(-1.0396495) q[3];
sx q[3];
rz(1.0218609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8675999) q[0];
sx q[0];
rz(-2.1430528) q[0];
sx q[0];
rz(-0.25319779) q[0];
rz(-0.7397488) q[1];
sx q[1];
rz(-0.7363798) q[1];
sx q[1];
rz(-0.69828066) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90726501) q[0];
sx q[0];
rz(-0.31120473) q[0];
sx q[0];
rz(2.9051203) q[0];
x q[1];
rz(-1.6449498) q[2];
sx q[2];
rz(-1.5955131) q[2];
sx q[2];
rz(-2.2029684) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.89430289) q[1];
sx q[1];
rz(-1.786788) q[1];
sx q[1];
rz(2.0386019) q[1];
x q[2];
rz(-2.9000557) q[3];
sx q[3];
rz(-2.1250484) q[3];
sx q[3];
rz(-2.9709771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0239821) q[2];
sx q[2];
rz(-0.88279804) q[2];
sx q[2];
rz(-2.3013766) q[2];
rz(-0.64030567) q[3];
sx q[3];
rz(-2.2655723) q[3];
sx q[3];
rz(-1.1673814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8828076) q[0];
sx q[0];
rz(-0.85492815) q[0];
sx q[0];
rz(1.7229765) q[0];
rz(3.1124658) q[1];
sx q[1];
rz(-3.0283785) q[1];
sx q[1];
rz(1.821847) q[1];
rz(-0.48568934) q[2];
sx q[2];
rz(-2.4808241) q[2];
sx q[2];
rz(-0.11447699) q[2];
rz(1.7520262) q[3];
sx q[3];
rz(-0.69835865) q[3];
sx q[3];
rz(-2.9327304) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
