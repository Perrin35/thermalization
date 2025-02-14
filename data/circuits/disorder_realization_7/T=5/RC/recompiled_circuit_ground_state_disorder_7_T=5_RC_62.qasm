OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.4829798) q[0];
sx q[0];
rz(-2.7375672) q[0];
sx q[0];
rz(-0.39028302) q[0];
rz(-0.033493869) q[1];
sx q[1];
rz(3.6727603) q[1];
sx q[1];
rz(9.6277278) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4056978) q[0];
sx q[0];
rz(-2.4740088) q[0];
sx q[0];
rz(2.1632458) q[0];
rz(-1.6597117) q[2];
sx q[2];
rz(-1.5651859) q[2];
sx q[2];
rz(2.5054431) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6773078) q[1];
sx q[1];
rz(-0.57972902) q[1];
sx q[1];
rz(-2.7955758) q[1];
rz(-0.84439028) q[3];
sx q[3];
rz(-1.5914306) q[3];
sx q[3];
rz(0.72935361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1927294) q[2];
sx q[2];
rz(-0.84550965) q[2];
sx q[2];
rz(1.753099) q[2];
rz(3.0697611) q[3];
sx q[3];
rz(-2.6303232) q[3];
sx q[3];
rz(-0.82740074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0440867) q[0];
sx q[0];
rz(-0.16479099) q[0];
sx q[0];
rz(-0.49758115) q[0];
rz(1.7454106) q[1];
sx q[1];
rz(-1.0284245) q[1];
sx q[1];
rz(-2.5713249) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6570396) q[0];
sx q[0];
rz(-1.4385106) q[0];
sx q[0];
rz(-0.45251493) q[0];
x q[1];
rz(-2.3421418) q[2];
sx q[2];
rz(-2.852147) q[2];
sx q[2];
rz(-1.4674526) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.1062938) q[1];
sx q[1];
rz(-0.72162823) q[1];
sx q[1];
rz(-0.053573805) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2482613) q[3];
sx q[3];
rz(-0.81557298) q[3];
sx q[3];
rz(-1.4295242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0294864) q[2];
sx q[2];
rz(-1.4542397) q[2];
sx q[2];
rz(-2.4278329) q[2];
rz(1.7589689) q[3];
sx q[3];
rz(-2.6191923) q[3];
sx q[3];
rz(-2.6944323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56938982) q[0];
sx q[0];
rz(-2.1694006) q[0];
sx q[0];
rz(2.8045281) q[0];
rz(-2.6481248) q[1];
sx q[1];
rz(-2.4512873) q[1];
sx q[1];
rz(0.92672551) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2703646) q[0];
sx q[0];
rz(-2.7634826) q[0];
sx q[0];
rz(2.3367995) q[0];
x q[1];
rz(2.3078467) q[2];
sx q[2];
rz(-2.0485224) q[2];
sx q[2];
rz(1.1857978) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2973916) q[1];
sx q[1];
rz(-1.1549885) q[1];
sx q[1];
rz(-0.17900055) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.95229228) q[3];
sx q[3];
rz(-0.77889788) q[3];
sx q[3];
rz(-0.84918815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9530764) q[2];
sx q[2];
rz(-0.86543721) q[2];
sx q[2];
rz(-2.4600929) q[2];
rz(-1.513688) q[3];
sx q[3];
rz(-1.8047921) q[3];
sx q[3];
rz(-0.17076913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82591581) q[0];
sx q[0];
rz(-3.0382394) q[0];
sx q[0];
rz(-2.090825) q[0];
rz(0.61327618) q[1];
sx q[1];
rz(-2.3502217) q[1];
sx q[1];
rz(1.223986) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.104804) q[0];
sx q[0];
rz(-2.4147644) q[0];
sx q[0];
rz(2.6340371) q[0];
x q[1];
rz(0.87665571) q[2];
sx q[2];
rz(-2.3842065) q[2];
sx q[2];
rz(0.11981431) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5526841) q[1];
sx q[1];
rz(-0.97637227) q[1];
sx q[1];
rz(1.5556704) q[1];
rz(-pi) q[2];
rz(1.7324034) q[3];
sx q[3];
rz(-0.9847042) q[3];
sx q[3];
rz(-1.8047892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0497389) q[2];
sx q[2];
rz(-2.3564796) q[2];
sx q[2];
rz(2.4741057) q[2];
rz(1.5664172) q[3];
sx q[3];
rz(-2.9508041) q[3];
sx q[3];
rz(3.0070087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.62395537) q[0];
sx q[0];
rz(-1.0409545) q[0];
sx q[0];
rz(-2.5868296) q[0];
rz(1.293921) q[1];
sx q[1];
rz(-2.7298253) q[1];
sx q[1];
rz(2.256934) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3356757) q[0];
sx q[0];
rz(-1.5461304) q[0];
sx q[0];
rz(-1.0403344) q[0];
rz(2.6045962) q[2];
sx q[2];
rz(-2.3648713) q[2];
sx q[2];
rz(-3.1077488) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6044449) q[1];
sx q[1];
rz(-1.7824518) q[1];
sx q[1];
rz(2.0404633) q[1];
rz(2.7001722) q[3];
sx q[3];
rz(-1.316464) q[3];
sx q[3];
rz(2.1770432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.027792949) q[2];
sx q[2];
rz(-0.5641368) q[2];
sx q[2];
rz(-1.0998868) q[2];
rz(0.22282985) q[3];
sx q[3];
rz(-1.204071) q[3];
sx q[3];
rz(3.0075464) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9851538) q[0];
sx q[0];
rz(-0.30010656) q[0];
sx q[0];
rz(1.1837748) q[0];
rz(1.3166332) q[1];
sx q[1];
rz(-2.1778409) q[1];
sx q[1];
rz(-2.1060941) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98052227) q[0];
sx q[0];
rz(-2.5455406) q[0];
sx q[0];
rz(-2.1658394) q[0];
rz(-pi) q[1];
rz(0.42607985) q[2];
sx q[2];
rz(-1.0023062) q[2];
sx q[2];
rz(2.7531433) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.07044) q[1];
sx q[1];
rz(-2.8912918) q[1];
sx q[1];
rz(-2.3790199) q[1];
x q[2];
rz(-1.6430969) q[3];
sx q[3];
rz(-0.80567718) q[3];
sx q[3];
rz(-1.8206545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.7684795) q[2];
sx q[2];
rz(-0.36102411) q[2];
sx q[2];
rz(-2.7601472) q[2];
rz(-1.2162195) q[3];
sx q[3];
rz(-0.65665025) q[3];
sx q[3];
rz(-0.60429627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4878047) q[0];
sx q[0];
rz(-2.8233546) q[0];
sx q[0];
rz(2.8741264) q[0];
rz(-1.5221315) q[1];
sx q[1];
rz(-0.61264241) q[1];
sx q[1];
rz(0.47346514) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7077443) q[0];
sx q[0];
rz(-1.9421541) q[0];
sx q[0];
rz(1.9385563) q[0];
x q[1];
rz(0.47961819) q[2];
sx q[2];
rz(-1.7349744) q[2];
sx q[2];
rz(2.2043383) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5468041) q[1];
sx q[1];
rz(-1.4546397) q[1];
sx q[1];
rz(-2.919477) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7438873) q[3];
sx q[3];
rz(-0.4022214) q[3];
sx q[3];
rz(-0.88170748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.21753103) q[2];
sx q[2];
rz(-1.4950098) q[2];
sx q[2];
rz(1.7688497) q[2];
rz(-2.8352906) q[3];
sx q[3];
rz(-2.114571) q[3];
sx q[3];
rz(2.8021804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31996763) q[0];
sx q[0];
rz(-0.29682934) q[0];
sx q[0];
rz(2.6676275) q[0];
rz(-0.78556806) q[1];
sx q[1];
rz(-0.57890099) q[1];
sx q[1];
rz(0.89326352) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5441204) q[0];
sx q[0];
rz(-1.5878146) q[0];
sx q[0];
rz(1.3025137) q[0];
rz(-1.947587) q[2];
sx q[2];
rz(-0.8424785) q[2];
sx q[2];
rz(-1.467799) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.230669) q[1];
sx q[1];
rz(-1.603447) q[1];
sx q[1];
rz(-2.9327675) q[1];
x q[2];
rz(-1.1749218) q[3];
sx q[3];
rz(-1.4787889) q[3];
sx q[3];
rz(-2.8407872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.47024176) q[2];
sx q[2];
rz(-1.525815) q[2];
sx q[2];
rz(-0.5980171) q[2];
rz(-2.9371069) q[3];
sx q[3];
rz(-3.0308767) q[3];
sx q[3];
rz(-2.428875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.04190271) q[0];
sx q[0];
rz(-2.1429017) q[0];
sx q[0];
rz(2.4424851) q[0];
rz(0.39012575) q[1];
sx q[1];
rz(-0.68123078) q[1];
sx q[1];
rz(2.1652538) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6481144) q[0];
sx q[0];
rz(-2.8472487) q[0];
sx q[0];
rz(2.074138) q[0];
x q[1];
rz(1.2430698) q[2];
sx q[2];
rz(-1.9877745) q[2];
sx q[2];
rz(-2.8627739) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.38473746) q[1];
sx q[1];
rz(-2.6275674) q[1];
sx q[1];
rz(1.4207178) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5849708) q[3];
sx q[3];
rz(-1.8020013) q[3];
sx q[3];
rz(-1.7116261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1583027) q[2];
sx q[2];
rz(-0.96902865) q[2];
sx q[2];
rz(2.9950673) q[2];
rz(2.8807785) q[3];
sx q[3];
rz(-1.1365889) q[3];
sx q[3];
rz(2.8543616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(0.89676595) q[0];
sx q[0];
rz(-0.34938669) q[0];
sx q[0];
rz(2.8283327) q[0];
rz(0.80097711) q[1];
sx q[1];
rz(-1.6311092) q[1];
sx q[1];
rz(-0.40447485) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6666732) q[0];
sx q[0];
rz(-1.7099147) q[0];
sx q[0];
rz(2.9551201) q[0];
x q[1];
rz(0.54623917) q[2];
sx q[2];
rz(-1.7405542) q[2];
sx q[2];
rz(1.7129829) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.69533379) q[1];
sx q[1];
rz(-1.0647758) q[1];
sx q[1];
rz(-0.74600469) q[1];
rz(-pi) q[2];
rz(1.9451109) q[3];
sx q[3];
rz(-2.7817543) q[3];
sx q[3];
rz(1.8051749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6910088) q[2];
sx q[2];
rz(-0.49387026) q[2];
sx q[2];
rz(-2.3502926) q[2];
rz(-2.6386236) q[3];
sx q[3];
rz(-1.0478323) q[3];
sx q[3];
rz(-2.3414229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41225152) q[0];
sx q[0];
rz(-1.8251735) q[0];
sx q[0];
rz(1.770021) q[0];
rz(-0.62660632) q[1];
sx q[1];
rz(-1.4817487) q[1];
sx q[1];
rz(2.1201835) q[1];
rz(0.37143636) q[2];
sx q[2];
rz(-1.7219661) q[2];
sx q[2];
rz(0.013602376) q[2];
rz(3.088763) q[3];
sx q[3];
rz(-2.9078447) q[3];
sx q[3];
rz(-1.71753) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
