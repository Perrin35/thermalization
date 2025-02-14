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
rz(2.2223284) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60828078) q[0];
sx q[0];
rz(-1.601378) q[0];
sx q[0];
rz(-1.5879059) q[0];
rz(-pi) q[1];
rz(-1.4798574) q[2];
sx q[2];
rz(-1.6906066) q[2];
sx q[2];
rz(2.7594942) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8265958) q[1];
sx q[1];
rz(-0.8657786) q[1];
sx q[1];
rz(3.0211074) q[1];
rz(-pi) q[2];
rz(2.0415061) q[3];
sx q[3];
rz(-1.19095) q[3];
sx q[3];
rz(2.0947411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8093402) q[2];
sx q[2];
rz(-1.6925749) q[2];
sx q[2];
rz(-1.7285041) q[2];
rz(2.938802) q[3];
sx q[3];
rz(-1.3821802) q[3];
sx q[3];
rz(-3.0387759) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8973812) q[0];
sx q[0];
rz(-1.0845217) q[0];
sx q[0];
rz(-0.71075034) q[0];
rz(0.76849014) q[1];
sx q[1];
rz(-1.0667421) q[1];
sx q[1];
rz(1.01952) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9907889) q[0];
sx q[0];
rz(-1.6103364) q[0];
sx q[0];
rz(0.44724748) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4065845) q[2];
sx q[2];
rz(-2.7052042) q[2];
sx q[2];
rz(-1.1498888) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.92495698) q[1];
sx q[1];
rz(-1.5397738) q[1];
sx q[1];
rz(1.9418632) q[1];
x q[2];
rz(1.7048265) q[3];
sx q[3];
rz(-2.5833327) q[3];
sx q[3];
rz(2.9001065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.76356137) q[2];
sx q[2];
rz(-1.8969994) q[2];
sx q[2];
rz(-1.9192609) q[2];
rz(-1.9289121) q[3];
sx q[3];
rz(-2.1247532) q[3];
sx q[3];
rz(0.71162629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0843622) q[0];
sx q[0];
rz(-2.1570692) q[0];
sx q[0];
rz(1.9236176) q[0];
rz(2.6257264) q[1];
sx q[1];
rz(-2.5936544) q[1];
sx q[1];
rz(2.2191494) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9410011) q[0];
sx q[0];
rz(-3.0709478) q[0];
sx q[0];
rz(-0.61265041) q[0];
rz(-pi) q[1];
rz(2.5848021) q[2];
sx q[2];
rz(-2.3344731) q[2];
sx q[2];
rz(0.082376235) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.53239142) q[1];
sx q[1];
rz(-2.2431886) q[1];
sx q[1];
rz(2.3228541) q[1];
rz(-0.46307989) q[3];
sx q[3];
rz(-2.7768917) q[3];
sx q[3];
rz(-1.1187584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.1232274) q[2];
sx q[2];
rz(-1.0845228) q[2];
sx q[2];
rz(1.8017192) q[2];
rz(-0.54723048) q[3];
sx q[3];
rz(-0.42357835) q[3];
sx q[3];
rz(-2.4303998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0860586) q[0];
sx q[0];
rz(-0.14500293) q[0];
sx q[0];
rz(0.15922971) q[0];
rz(0.010146443) q[1];
sx q[1];
rz(-1.0228913) q[1];
sx q[1];
rz(-0.25746447) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5978569) q[0];
sx q[0];
rz(-0.86559767) q[0];
sx q[0];
rz(-1.0969093) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0765216) q[2];
sx q[2];
rz(-0.79539585) q[2];
sx q[2];
rz(0.8000904) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3066669) q[1];
sx q[1];
rz(-1.67027) q[1];
sx q[1];
rz(2.5132781) q[1];
rz(-0.93521714) q[3];
sx q[3];
rz(-2.2699589) q[3];
sx q[3];
rz(-0.88672598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3782392) q[2];
sx q[2];
rz(-1.2946318) q[2];
sx q[2];
rz(2.6687458) q[2];
rz(2.4371448) q[3];
sx q[3];
rz(-1.7770146) q[3];
sx q[3];
rz(-0.64594597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6351629) q[0];
sx q[0];
rz(-1.1744873) q[0];
sx q[0];
rz(3.0761062) q[0];
rz(0.40924117) q[1];
sx q[1];
rz(-1.1301273) q[1];
sx q[1];
rz(1.7154891) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.586389) q[0];
sx q[0];
rz(-0.23815933) q[0];
sx q[0];
rz(-0.40443964) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.481856) q[2];
sx q[2];
rz(-1.9581902) q[2];
sx q[2];
rz(2.873444) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9551505) q[1];
sx q[1];
rz(-1.7935866) q[1];
sx q[1];
rz(-0.28526116) q[1];
x q[2];
rz(2.0992005) q[3];
sx q[3];
rz(-0.33337731) q[3];
sx q[3];
rz(2.7957145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.19491974) q[2];
sx q[2];
rz(-2.0544923) q[2];
sx q[2];
rz(1.3067513) q[2];
rz(-1.9715747) q[3];
sx q[3];
rz(-2.7368059) q[3];
sx q[3];
rz(3.034333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4829247) q[0];
sx q[0];
rz(-1.985745) q[0];
sx q[0];
rz(-2.7401127) q[0];
rz(-2.4688156) q[1];
sx q[1];
rz(-2.3428226) q[1];
sx q[1];
rz(2.2875517) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6333165) q[0];
sx q[0];
rz(-2.9190953) q[0];
sx q[0];
rz(-1.7822687) q[0];
rz(-pi) q[1];
rz(0.18547345) q[2];
sx q[2];
rz(-1.7181686) q[2];
sx q[2];
rz(-2.3285248) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0103763) q[1];
sx q[1];
rz(-1.4713738) q[1];
sx q[1];
rz(2.5491284) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.61088125) q[3];
sx q[3];
rz(-1.1958836) q[3];
sx q[3];
rz(1.5883816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4868769) q[2];
sx q[2];
rz(-0.87149039) q[2];
sx q[2];
rz(-1.1775449) q[2];
rz(2.5901637) q[3];
sx q[3];
rz(-1.5827554) q[3];
sx q[3];
rz(2.3059755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7198782) q[0];
sx q[0];
rz(-2.8604909) q[0];
sx q[0];
rz(1.1623435) q[0];
rz(1.989919) q[1];
sx q[1];
rz(-1.3538227) q[1];
sx q[1];
rz(-0.91845671) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7699444) q[0];
sx q[0];
rz(-1.3027096) q[0];
sx q[0];
rz(-1.2275342) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2662042) q[2];
sx q[2];
rz(-1.3686485) q[2];
sx q[2];
rz(0.018176807) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.368093) q[1];
sx q[1];
rz(-1.5640904) q[1];
sx q[1];
rz(0.0073324629) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.17042589) q[3];
sx q[3];
rz(-1.9099351) q[3];
sx q[3];
rz(-2.2234774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1071757) q[2];
sx q[2];
rz(-1.7566661) q[2];
sx q[2];
rz(-2.6386063) q[2];
rz(-0.77477396) q[3];
sx q[3];
rz(-0.46190244) q[3];
sx q[3];
rz(1.7983961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69304943) q[0];
sx q[0];
rz(-0.22668426) q[0];
sx q[0];
rz(-1.6402798) q[0];
rz(0.34128183) q[1];
sx q[1];
rz(-1.7106067) q[1];
sx q[1];
rz(-1.1192082) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2687896) q[0];
sx q[0];
rz(-0.69376341) q[0];
sx q[0];
rz(-0.60141464) q[0];
rz(2.881024) q[2];
sx q[2];
rz(-2.0379279) q[2];
sx q[2];
rz(1.9595944) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.71257617) q[1];
sx q[1];
rz(-1.4498931) q[1];
sx q[1];
rz(-2.4790133) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.39466484) q[3];
sx q[3];
rz(-1.695444) q[3];
sx q[3];
rz(-3.0277071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.98562733) q[2];
sx q[2];
rz(-2.1817744) q[2];
sx q[2];
rz(-0.36925527) q[2];
rz(-2.8333832) q[3];
sx q[3];
rz(-0.9674558) q[3];
sx q[3];
rz(1.5306028) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8592598) q[0];
sx q[0];
rz(-1.8349324) q[0];
sx q[0];
rz(-0.34307137) q[0];
rz(1.169091) q[1];
sx q[1];
rz(-1.3645376) q[1];
sx q[1];
rz(-1.4287359) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.13641) q[0];
sx q[0];
rz(-2.6652968) q[0];
sx q[0];
rz(-0.57203697) q[0];
rz(-0.83115432) q[2];
sx q[2];
rz(-2.2546283) q[2];
sx q[2];
rz(1.9076951) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8642042) q[1];
sx q[1];
rz(-2.7608747) q[1];
sx q[1];
rz(1.4053132) q[1];
rz(-pi) q[2];
x q[2];
rz(0.47023021) q[3];
sx q[3];
rz(-0.93358002) q[3];
sx q[3];
rz(1.5772082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2009361) q[2];
sx q[2];
rz(-2.9328465) q[2];
sx q[2];
rz(1.7500056) q[2];
rz(2.7348147) q[3];
sx q[3];
rz(-1.6986366) q[3];
sx q[3];
rz(2.2627635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10130356) q[0];
sx q[0];
rz(-2.5387796) q[0];
sx q[0];
rz(-1.3579177) q[0];
rz(1.2376002) q[1];
sx q[1];
rz(-2.1289181) q[1];
sx q[1];
rz(-2.3416669) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3615204) q[0];
sx q[0];
rz(-1.3409541) q[0];
sx q[0];
rz(-2.337268) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2679891) q[2];
sx q[2];
rz(-1.0410415) q[2];
sx q[2];
rz(-2.985266) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4373191) q[1];
sx q[1];
rz(-1.4156439) q[1];
sx q[1];
rz(2.4680733) q[1];
rz(-pi) q[2];
rz(-1.1198197) q[3];
sx q[3];
rz(-0.52996892) q[3];
sx q[3];
rz(-1.3254904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.37821975) q[2];
sx q[2];
rz(-1.7174481) q[2];
sx q[2];
rz(3.0685032) q[2];
rz(-0.25589219) q[3];
sx q[3];
rz(-2.3292694) q[3];
sx q[3];
rz(1.488744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.273461) q[0];
sx q[0];
rz(-1.68597) q[0];
sx q[0];
rz(-1.2706533) q[0];
rz(-1.3399667) q[1];
sx q[1];
rz(-1.7908962) q[1];
sx q[1];
rz(-2.4790196) q[1];
rz(-1.5688098) q[2];
sx q[2];
rz(-0.86138267) q[2];
sx q[2];
rz(-1.3087261) q[2];
rz(2.1382016) q[3];
sx q[3];
rz(-1.0091253) q[3];
sx q[3];
rz(-0.71970018) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
