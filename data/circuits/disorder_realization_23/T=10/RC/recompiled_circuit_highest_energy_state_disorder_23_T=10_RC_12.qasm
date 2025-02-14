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
rz(-2.5459557) q[0];
sx q[0];
rz(-0.41002265) q[0];
sx q[0];
rz(-2.3469143) q[0];
rz(-1.0533286) q[1];
sx q[1];
rz(5.3310634) q[1];
sx q[1];
rz(7.4330243) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4276221) q[0];
sx q[0];
rz(-2.8454878) q[0];
sx q[0];
rz(1.4480736) q[0];
x q[1];
rz(-1.2501405) q[2];
sx q[2];
rz(-0.76092623) q[2];
sx q[2];
rz(-2.2087532) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9147891) q[1];
sx q[1];
rz(-2.7046596) q[1];
sx q[1];
rz(2.7257257) q[1];
rz(-0.47969476) q[3];
sx q[3];
rz(-2.624676) q[3];
sx q[3];
rz(1.2933047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6940234) q[2];
sx q[2];
rz(-2.1489216) q[2];
sx q[2];
rz(0.64219323) q[2];
rz(1.0688952) q[3];
sx q[3];
rz(-1.5725719) q[3];
sx q[3];
rz(1.6525035) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.883413) q[0];
sx q[0];
rz(-3.1386107) q[0];
sx q[0];
rz(2.6302443) q[0];
rz(1.5072352) q[1];
sx q[1];
rz(-1.8664482) q[1];
sx q[1];
rz(1.3828329) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94039932) q[0];
sx q[0];
rz(-2.7015238) q[0];
sx q[0];
rz(1.0690734) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2425425) q[2];
sx q[2];
rz(-1.8457638) q[2];
sx q[2];
rz(-3.0877047) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1630711) q[1];
sx q[1];
rz(-0.51632612) q[1];
sx q[1];
rz(-2.2341245) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6750349) q[3];
sx q[3];
rz(-2.0917419) q[3];
sx q[3];
rz(-1.5261861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4649268) q[2];
sx q[2];
rz(-1.9366465) q[2];
sx q[2];
rz(-0.10291544) q[2];
rz(2.0788705) q[3];
sx q[3];
rz(-1.1954185) q[3];
sx q[3];
rz(-0.11500558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9677143) q[0];
sx q[0];
rz(-1.0391087) q[0];
sx q[0];
rz(0.1097196) q[0];
rz(0.32490718) q[1];
sx q[1];
rz(-1.1512681) q[1];
sx q[1];
rz(-2.6085764) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95105458) q[0];
sx q[0];
rz(-1.8673273) q[0];
sx q[0];
rz(1.426479) q[0];
rz(-1.5948086) q[2];
sx q[2];
rz(-1.5188076) q[2];
sx q[2];
rz(2.1980224) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.17688454) q[1];
sx q[1];
rz(-1.0057276) q[1];
sx q[1];
rz(2.3749897) q[1];
rz(2.332052) q[3];
sx q[3];
rz(-1.24461) q[3];
sx q[3];
rz(2.1128138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.816232) q[2];
sx q[2];
rz(-1.3282789) q[2];
sx q[2];
rz(1.6858961) q[2];
rz(0.085518941) q[3];
sx q[3];
rz(-1.0695846) q[3];
sx q[3];
rz(2.1948309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.7607255) q[0];
sx q[0];
rz(-0.65422288) q[0];
sx q[0];
rz(0.72232676) q[0];
rz(1.5707312) q[1];
sx q[1];
rz(-1.9941565) q[1];
sx q[1];
rz(-3.0991203) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0642424) q[0];
sx q[0];
rz(-2.0127065) q[0];
sx q[0];
rz(-1.7741706) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1257079) q[2];
sx q[2];
rz(-1.6611036) q[2];
sx q[2];
rz(1.7866194) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.153923) q[1];
sx q[1];
rz(-2.3312285) q[1];
sx q[1];
rz(-2.8778879) q[1];
rz(-pi) q[2];
rz(0.23353429) q[3];
sx q[3];
rz(-0.81538768) q[3];
sx q[3];
rz(2.5366207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8954358) q[2];
sx q[2];
rz(-1.7652721) q[2];
sx q[2];
rz(-0.24406544) q[2];
rz(-0.79143381) q[3];
sx q[3];
rz(-1.8297198) q[3];
sx q[3];
rz(0.70146504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1112082) q[0];
sx q[0];
rz(-3.0272439) q[0];
sx q[0];
rz(2.9682888) q[0];
rz(-0.93470848) q[1];
sx q[1];
rz(-2.2917031) q[1];
sx q[1];
rz(-0.25397837) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5565709) q[0];
sx q[0];
rz(-0.75977548) q[0];
sx q[0];
rz(1.4069446) q[0];
rz(-1.9333657) q[2];
sx q[2];
rz(-2.1340279) q[2];
sx q[2];
rz(-0.24224396) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7813999) q[1];
sx q[1];
rz(-1.5500229) q[1];
sx q[1];
rz(0.69333496) q[1];
rz(3.0766904) q[3];
sx q[3];
rz(-1.1838473) q[3];
sx q[3];
rz(-0.98777366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0737334) q[2];
sx q[2];
rz(-0.24418712) q[2];
sx q[2];
rz(1.4402639) q[2];
rz(-0.67617792) q[3];
sx q[3];
rz(-1.9192326) q[3];
sx q[3];
rz(0.079782709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9192231) q[0];
sx q[0];
rz(-1.3105404) q[0];
sx q[0];
rz(2.4134912) q[0];
rz(1.9758196) q[1];
sx q[1];
rz(-2.246942) q[1];
sx q[1];
rz(0.56070915) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8090208) q[0];
sx q[0];
rz(-0.14223465) q[0];
sx q[0];
rz(0.47551544) q[0];
x q[1];
rz(-0.54975551) q[2];
sx q[2];
rz(-1.6719596) q[2];
sx q[2];
rz(1.2899952) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5594031) q[1];
sx q[1];
rz(-1.3546556) q[1];
sx q[1];
rz(1.1690421) q[1];
rz(-1.6339763) q[3];
sx q[3];
rz(-1.6681021) q[3];
sx q[3];
rz(0.41690505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5460983) q[2];
sx q[2];
rz(-1.4608773) q[2];
sx q[2];
rz(-1.072849) q[2];
rz(0.32235518) q[3];
sx q[3];
rz(-0.41283804) q[3];
sx q[3];
rz(2.7505007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3748465) q[0];
sx q[0];
rz(-1.640919) q[0];
sx q[0];
rz(-2.7333976) q[0];
rz(2.8152668) q[1];
sx q[1];
rz(-2.7944481) q[1];
sx q[1];
rz(1.4761285) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8360739) q[0];
sx q[0];
rz(-1.7461494) q[0];
sx q[0];
rz(-0.039964635) q[0];
rz(-2.1340685) q[2];
sx q[2];
rz(-1.3671698) q[2];
sx q[2];
rz(-1.3290392) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1154815) q[1];
sx q[1];
rz(-0.12899765) q[1];
sx q[1];
rz(-1.1327101) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9093698) q[3];
sx q[3];
rz(-1.7360064) q[3];
sx q[3];
rz(0.31295199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8866715) q[2];
sx q[2];
rz(-1.485606) q[2];
sx q[2];
rz(0.86611789) q[2];
rz(0.73650375) q[3];
sx q[3];
rz(-2.0564506) q[3];
sx q[3];
rz(-2.2542663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12164584) q[0];
sx q[0];
rz(-3.096088) q[0];
sx q[0];
rz(1.8628927) q[0];
rz(1.0026898) q[1];
sx q[1];
rz(-1.2607144) q[1];
sx q[1];
rz(-0.444828) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3747975) q[0];
sx q[0];
rz(-2.594769) q[0];
sx q[0];
rz(1.8575791) q[0];
rz(-pi) q[1];
rz(-2.41909) q[2];
sx q[2];
rz(-1.926427) q[2];
sx q[2];
rz(-2.1788545) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.63193601) q[1];
sx q[1];
rz(-2.3705814) q[1];
sx q[1];
rz(-1.3994751) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.37737198) q[3];
sx q[3];
rz(-1.8131953) q[3];
sx q[3];
rz(2.9494772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.43850809) q[2];
sx q[2];
rz(-1.4139516) q[2];
sx q[2];
rz(-0.50602305) q[2];
rz(-0.94432008) q[3];
sx q[3];
rz(-3.1157065) q[3];
sx q[3];
rz(-2.4054312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29189062) q[0];
sx q[0];
rz(-2.1493752) q[0];
sx q[0];
rz(-0.010183656) q[0];
rz(-2.5306375) q[1];
sx q[1];
rz(-0.98038951) q[1];
sx q[1];
rz(-0.082286509) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5389494) q[0];
sx q[0];
rz(-1.9398488) q[0];
sx q[0];
rz(2.6239388) q[0];
rz(2.9111229) q[2];
sx q[2];
rz(-0.90417143) q[2];
sx q[2];
rz(0.33483349) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4520784) q[1];
sx q[1];
rz(-1.6321215) q[1];
sx q[1];
rz(-0.21558) q[1];
rz(-2.4769267) q[3];
sx q[3];
rz(-1.3997404) q[3];
sx q[3];
rz(1.4881575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8102707) q[2];
sx q[2];
rz(-1.6931345) q[2];
sx q[2];
rz(-0.70915478) q[2];
rz(1.6592615) q[3];
sx q[3];
rz(-2.5100561) q[3];
sx q[3];
rz(-0.46019301) q[3];
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
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.822478) q[0];
sx q[0];
rz(-2.3111486) q[0];
sx q[0];
rz(-0.6293695) q[0];
rz(1.7259701) q[1];
sx q[1];
rz(-1.6547838) q[1];
sx q[1];
rz(1.1782882) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8346658) q[0];
sx q[0];
rz(-0.7531868) q[0];
sx q[0];
rz(0.96145328) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0731931) q[2];
sx q[2];
rz(-1.8821239) q[2];
sx q[2];
rz(-2.5466998) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.1078892) q[1];
sx q[1];
rz(-0.7830355) q[1];
sx q[1];
rz(1.7262001) q[1];
rz(0.62507665) q[3];
sx q[3];
rz(-0.82603329) q[3];
sx q[3];
rz(2.1476098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2668931) q[2];
sx q[2];
rz(-1.5723672) q[2];
sx q[2];
rz(2.3229522) q[2];
rz(1.6308174) q[3];
sx q[3];
rz(-1.2997593) q[3];
sx q[3];
rz(-0.40217933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4928987) q[0];
sx q[0];
rz(-1.10981) q[0];
sx q[0];
rz(1.4640402) q[0];
rz(-1.2661487) q[1];
sx q[1];
rz(-1.9910973) q[1];
sx q[1];
rz(1.533351) q[1];
rz(0.89546236) q[2];
sx q[2];
rz(-1.7747028) q[2];
sx q[2];
rz(1.5595421) q[2];
rz(-2.9569382) q[3];
sx q[3];
rz(-1.6600556) q[3];
sx q[3];
rz(-0.32515812) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
