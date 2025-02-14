OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.74066585) q[0];
sx q[0];
rz(-1.292115) q[0];
sx q[0];
rz(2.3722755) q[0];
rz(1.5953335) q[1];
sx q[1];
rz(-2.9268664) q[1];
sx q[1];
rz(-0.62155849) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9233166) q[0];
sx q[0];
rz(-2.0686604) q[0];
sx q[0];
rz(0.68572361) q[0];
x q[1];
rz(3.123986) q[2];
sx q[2];
rz(-1.1531534) q[2];
sx q[2];
rz(2.524802) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4159704) q[1];
sx q[1];
rz(-0.86639222) q[1];
sx q[1];
rz(-2.7654057) q[1];
rz(3.0075843) q[3];
sx q[3];
rz(-1.4913627) q[3];
sx q[3];
rz(1.1955737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8605211) q[2];
sx q[2];
rz(-2.8513384) q[2];
sx q[2];
rz(2.2104635) q[2];
rz(0.18533254) q[3];
sx q[3];
rz(-1.0328181) q[3];
sx q[3];
rz(1.7271656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0463878) q[0];
sx q[0];
rz(-2.8075908) q[0];
sx q[0];
rz(-0.97292501) q[0];
rz(-2.3880549) q[1];
sx q[1];
rz(-2.2851508) q[1];
sx q[1];
rz(0.10261745) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30841973) q[0];
sx q[0];
rz(-3.0121332) q[0];
sx q[0];
rz(0.85516103) q[0];
rz(-pi) q[1];
rz(-0.13394103) q[2];
sx q[2];
rz(-1.6689577) q[2];
sx q[2];
rz(-2.5636473) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.953578) q[1];
sx q[1];
rz(-0.17718592) q[1];
sx q[1];
rz(-2.0582951) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9512734) q[3];
sx q[3];
rz(-2.7822436) q[3];
sx q[3];
rz(0.41216601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.40833452) q[2];
sx q[2];
rz(-1.4380941) q[2];
sx q[2];
rz(0.56780887) q[2];
rz(-0.37219498) q[3];
sx q[3];
rz(-1.8122383) q[3];
sx q[3];
rz(1.6224434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(2.5837625) q[0];
sx q[0];
rz(-0.74302858) q[0];
sx q[0];
rz(-1.8999735) q[0];
rz(2.325233) q[1];
sx q[1];
rz(-2.7528449) q[1];
sx q[1];
rz(-0.98966086) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69520607) q[0];
sx q[0];
rz(-1.1940845) q[0];
sx q[0];
rz(-1.8003617) q[0];
rz(-1.4098566) q[2];
sx q[2];
rz(-1.4795627) q[2];
sx q[2];
rz(-0.021878069) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6523) q[1];
sx q[1];
rz(-1.1649141) q[1];
sx q[1];
rz(1.7651221) q[1];
x q[2];
rz(-3.039547) q[3];
sx q[3];
rz(-1.1983173) q[3];
sx q[3];
rz(0.88119394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7648387) q[2];
sx q[2];
rz(-1.2412485) q[2];
sx q[2];
rz(1.6468916) q[2];
rz(0.74690789) q[3];
sx q[3];
rz(-2.5266095) q[3];
sx q[3];
rz(-0.99228215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84585369) q[0];
sx q[0];
rz(-0.60765147) q[0];
sx q[0];
rz(2.1920152) q[0];
rz(-1.9751366) q[1];
sx q[1];
rz(-1.8828853) q[1];
sx q[1];
rz(-2.9393401) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40914772) q[0];
sx q[0];
rz(-1.5209001) q[0];
sx q[0];
rz(1.6014535) q[0];
x q[1];
rz(-2.4522547) q[2];
sx q[2];
rz(-1.4408913) q[2];
sx q[2];
rz(-0.65016937) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5626853) q[1];
sx q[1];
rz(-2.1386801) q[1];
sx q[1];
rz(-0.61526362) q[1];
x q[2];
rz(-2.5519484) q[3];
sx q[3];
rz(-1.0760376) q[3];
sx q[3];
rz(2.445666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.64968455) q[2];
sx q[2];
rz(-1.4338355) q[2];
sx q[2];
rz(-1.8053619) q[2];
rz(2.6797471) q[3];
sx q[3];
rz(-2.0838085) q[3];
sx q[3];
rz(-1.0285146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.860054) q[0];
sx q[0];
rz(-3.0335732) q[0];
sx q[0];
rz(2.7463013) q[0];
rz(2.1832502) q[1];
sx q[1];
rz(-1.3805362) q[1];
sx q[1];
rz(2.7470632) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1720264) q[0];
sx q[0];
rz(-0.2663258) q[0];
sx q[0];
rz(-1.5035433) q[0];
rz(-pi) q[1];
rz(-1.2917323) q[2];
sx q[2];
rz(-1.9968281) q[2];
sx q[2];
rz(-2.5468777) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3580618) q[1];
sx q[1];
rz(-0.79800843) q[1];
sx q[1];
rz(-2.0302057) q[1];
rz(-1.35704) q[3];
sx q[3];
rz(-0.96352623) q[3];
sx q[3];
rz(0.34363817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.40256527) q[2];
sx q[2];
rz(-1.5274916) q[2];
sx q[2];
rz(2.1220477) q[2];
rz(-1.839365) q[3];
sx q[3];
rz(-3.0890833) q[3];
sx q[3];
rz(2.4359865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5810982) q[0];
sx q[0];
rz(-0.56466931) q[0];
sx q[0];
rz(-2.1288921) q[0];
rz(0.46214354) q[1];
sx q[1];
rz(-1.7299165) q[1];
sx q[1];
rz(0.046646811) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29829866) q[0];
sx q[0];
rz(-3.0180535) q[0];
sx q[0];
rz(-1.6575097) q[0];
x q[1];
rz(1.8671124) q[2];
sx q[2];
rz(-1.3256529) q[2];
sx q[2];
rz(-0.60173881) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0477284) q[1];
sx q[1];
rz(-2.1360435) q[1];
sx q[1];
rz(0.39611343) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6316333) q[3];
sx q[3];
rz(-0.3777856) q[3];
sx q[3];
rz(1.2041262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3095653) q[2];
sx q[2];
rz(-3.053061) q[2];
sx q[2];
rz(2.1868165) q[2];
rz(2.471762) q[3];
sx q[3];
rz(-1.9570743) q[3];
sx q[3];
rz(-0.36439103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1466115) q[0];
sx q[0];
rz(-0.28047383) q[0];
sx q[0];
rz(1.9914419) q[0];
rz(-0.096788302) q[1];
sx q[1];
rz(-2.0014747) q[1];
sx q[1];
rz(1.4801625) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32616781) q[0];
sx q[0];
rz(-2.302044) q[0];
sx q[0];
rz(-2.6214548) q[0];
rz(-pi) q[1];
rz(-1.3315013) q[2];
sx q[2];
rz(-1.2579489) q[2];
sx q[2];
rz(-2.8205736) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4188618) q[1];
sx q[1];
rz(-1.7081385) q[1];
sx q[1];
rz(1.2773499) q[1];
x q[2];
rz(1.1365165) q[3];
sx q[3];
rz(-1.4308813) q[3];
sx q[3];
rz(2.851448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4074771) q[2];
sx q[2];
rz(-2.263767) q[2];
sx q[2];
rz(-2.5131098) q[2];
rz(2.8001522) q[3];
sx q[3];
rz(-2.1554558) q[3];
sx q[3];
rz(-0.83610523) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2899365) q[0];
sx q[0];
rz(-2.4735232) q[0];
sx q[0];
rz(-0.063902721) q[0];
rz(-2.5996161) q[1];
sx q[1];
rz(-2.0109476) q[1];
sx q[1];
rz(-2.7625387) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5139007) q[0];
sx q[0];
rz(-1.3435044) q[0];
sx q[0];
rz(-0.17666575) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7070243) q[2];
sx q[2];
rz(-2.3751525) q[2];
sx q[2];
rz(-2.5640783) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8833911) q[1];
sx q[1];
rz(-2.5978226) q[1];
sx q[1];
rz(0.030042458) q[1];
rz(-pi) q[2];
rz(0.98446357) q[3];
sx q[3];
rz(-1.9317614) q[3];
sx q[3];
rz(0.20853768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9717676) q[2];
sx q[2];
rz(-2.0242033) q[2];
sx q[2];
rz(0.28905147) q[2];
rz(-1.0158094) q[3];
sx q[3];
rz(-0.93520516) q[3];
sx q[3];
rz(-2.9142006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0357901) q[0];
sx q[0];
rz(-0.38689026) q[0];
sx q[0];
rz(-2.8111358) q[0];
rz(-0.48078787) q[1];
sx q[1];
rz(-1.5590706) q[1];
sx q[1];
rz(-0.50723433) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.071166) q[0];
sx q[0];
rz(-2.2615763) q[0];
sx q[0];
rz(2.4450177) q[0];
rz(-pi) q[1];
rz(2.4749996) q[2];
sx q[2];
rz(-2.5799516) q[2];
sx q[2];
rz(3.0009342) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6751409) q[1];
sx q[1];
rz(-1.351739) q[1];
sx q[1];
rz(-1.1417927) q[1];
rz(2.5413569) q[3];
sx q[3];
rz(-1.5215989) q[3];
sx q[3];
rz(0.95358301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.5759739) q[2];
sx q[2];
rz(-0.86056346) q[2];
sx q[2];
rz(2.7536075) q[2];
rz(3.0787789) q[3];
sx q[3];
rz(-0.7395491) q[3];
sx q[3];
rz(1.1886103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0906618) q[0];
sx q[0];
rz(-3.1212786) q[0];
sx q[0];
rz(0.055572979) q[0];
rz(-2.9203501) q[1];
sx q[1];
rz(-2.1317) q[1];
sx q[1];
rz(1.5004858) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21723233) q[0];
sx q[0];
rz(-1.8460197) q[0];
sx q[0];
rz(-1.0479529) q[0];
rz(-3.0048641) q[2];
sx q[2];
rz(-1.5216344) q[2];
sx q[2];
rz(1.8135742) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5855687) q[1];
sx q[1];
rz(-0.50630664) q[1];
sx q[1];
rz(-0.8697747) q[1];
rz(-pi) q[2];
rz(2.5440823) q[3];
sx q[3];
rz(-2.198285) q[3];
sx q[3];
rz(1.6399872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9556433) q[2];
sx q[2];
rz(-2.9538437) q[2];
sx q[2];
rz(-1.1090247) q[2];
rz(0.016131314) q[3];
sx q[3];
rz(-1.4343836) q[3];
sx q[3];
rz(0.46512887) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3817417) q[0];
sx q[0];
rz(-2.1325337) q[0];
sx q[0];
rz(2.7251563) q[0];
rz(0.79553678) q[1];
sx q[1];
rz(-1.3858613) q[1];
sx q[1];
rz(-0.081079986) q[1];
rz(1.9151081) q[2];
sx q[2];
rz(-1.7034265) q[2];
sx q[2];
rz(-1.8272057) q[2];
rz(-0.78043191) q[3];
sx q[3];
rz(-1.9698236) q[3];
sx q[3];
rz(1.3079328) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
