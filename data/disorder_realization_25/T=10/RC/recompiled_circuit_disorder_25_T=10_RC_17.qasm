OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0796354) q[0];
sx q[0];
rz(-2.893653) q[0];
sx q[0];
rz(-2.3214582) q[0];
rz(-3.1007383) q[1];
sx q[1];
rz(-0.78894579) q[1];
sx q[1];
rz(3.0537002) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9633006) q[0];
sx q[0];
rz(-1.1305729) q[0];
sx q[0];
rz(0.51142366) q[0];
x q[1];
rz(-2.0983661) q[2];
sx q[2];
rz(-1.9448148) q[2];
sx q[2];
rz(2.5803215) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3275798) q[1];
sx q[1];
rz(-0.50427932) q[1];
sx q[1];
rz(-3.1370647) q[1];
x q[2];
rz(-0.0058562972) q[3];
sx q[3];
rz(-1.8294888) q[3];
sx q[3];
rz(-0.70219016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1408046) q[2];
sx q[2];
rz(-1.452383) q[2];
sx q[2];
rz(-0.56935707) q[2];
rz(-1.6128444) q[3];
sx q[3];
rz(-2.5053535) q[3];
sx q[3];
rz(-1.3585842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59812087) q[0];
sx q[0];
rz(-0.16508979) q[0];
sx q[0];
rz(-2.5879481) q[0];
rz(1.9042632) q[1];
sx q[1];
rz(-1.3668704) q[1];
sx q[1];
rz(1.233261) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90399088) q[0];
sx q[0];
rz(-2.1008137) q[0];
sx q[0];
rz(2.3474098) q[0];
rz(-pi) q[1];
rz(-2.9821175) q[2];
sx q[2];
rz(-0.6808241) q[2];
sx q[2];
rz(2.8603539) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.065029) q[1];
sx q[1];
rz(-0.20670465) q[1];
sx q[1];
rz(-1.6480584) q[1];
rz(-pi) q[2];
x q[2];
rz(0.29970701) q[3];
sx q[3];
rz(-2.4518659) q[3];
sx q[3];
rz(-2.0852065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.0040940293) q[2];
sx q[2];
rz(-1.4592905) q[2];
sx q[2];
rz(3.1393576) q[2];
rz(-2.3114752) q[3];
sx q[3];
rz(-0.79289645) q[3];
sx q[3];
rz(1.8921651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8925791) q[0];
sx q[0];
rz(-0.61616388) q[0];
sx q[0];
rz(0.85154831) q[0];
rz(-0.7710723) q[1];
sx q[1];
rz(-1.4665736) q[1];
sx q[1];
rz(3.070389) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.035633798) q[0];
sx q[0];
rz(-1.3741612) q[0];
sx q[0];
rz(-1.3282302) q[0];
x q[1];
rz(1.0686915) q[2];
sx q[2];
rz(-2.6325668) q[2];
sx q[2];
rz(2.6976762) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8098231) q[1];
sx q[1];
rz(-0.66546813) q[1];
sx q[1];
rz(-0.21824093) q[1];
rz(-pi) q[2];
rz(2.4124868) q[3];
sx q[3];
rz(-2.2359072) q[3];
sx q[3];
rz(-0.8641181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3123793) q[2];
sx q[2];
rz(-0.92386121) q[2];
sx q[2];
rz(-2.6181347) q[2];
rz(-2.8097025) q[3];
sx q[3];
rz(-3.0811946) q[3];
sx q[3];
rz(0.50326842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7317384) q[0];
sx q[0];
rz(-2.8926352) q[0];
sx q[0];
rz(2.3160034) q[0];
rz(1.1666974) q[1];
sx q[1];
rz(-2.2749133) q[1];
sx q[1];
rz(-1.4473787) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4226845) q[0];
sx q[0];
rz(-1.8005162) q[0];
sx q[0];
rz(2.9740888) q[0];
rz(-pi) q[1];
rz(-3.1324603) q[2];
sx q[2];
rz(-1.3924686) q[2];
sx q[2];
rz(-0.4724617) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4059658) q[1];
sx q[1];
rz(-1.7298797) q[1];
sx q[1];
rz(-2.2707978) q[1];
x q[2];
rz(-0.36376603) q[3];
sx q[3];
rz(-1.7978661) q[3];
sx q[3];
rz(-1.0432537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6772785) q[2];
sx q[2];
rz(-1.364664) q[2];
sx q[2];
rz(0.14492598) q[2];
rz(-2.1302917) q[3];
sx q[3];
rz(-0.62711182) q[3];
sx q[3];
rz(-3.1269126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99575627) q[0];
sx q[0];
rz(-0.053481426) q[0];
sx q[0];
rz(-0.68403912) q[0];
rz(1.9513291) q[1];
sx q[1];
rz(-2.4763156) q[1];
sx q[1];
rz(1.2129983) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0096668) q[0];
sx q[0];
rz(-1.8312694) q[0];
sx q[0];
rz(-2.7564604) q[0];
rz(-pi) q[1];
x q[1];
rz(0.80981363) q[2];
sx q[2];
rz(-0.81293101) q[2];
sx q[2];
rz(2.6990858) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0141543) q[1];
sx q[1];
rz(-1.6886097) q[1];
sx q[1];
rz(-0.25728667) q[1];
rz(-pi) q[2];
rz(2.1608859) q[3];
sx q[3];
rz(-0.31168918) q[3];
sx q[3];
rz(-2.2463472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.50903901) q[2];
sx q[2];
rz(-1.2712487) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0260139) q[0];
sx q[0];
rz(-0.84340874) q[0];
sx q[0];
rz(-0.28636006) q[0];
rz(-0.59965602) q[1];
sx q[1];
rz(-1.1279794) q[1];
sx q[1];
rz(0.58475959) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8557381) q[0];
sx q[0];
rz(-1.2585088) q[0];
sx q[0];
rz(0.62494846) q[0];
rz(2.5271687) q[2];
sx q[2];
rz(-1.5493869) q[2];
sx q[2];
rz(-1.9911839) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6905745) q[1];
sx q[1];
rz(-1.283053) q[1];
sx q[1];
rz(-2.2698127) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9752475) q[3];
sx q[3];
rz(-1.6206985) q[3];
sx q[3];
rz(0.14453105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.16453234) q[2];
sx q[2];
rz(-0.88447276) q[2];
sx q[2];
rz(-1.6284846) q[2];
rz(-0.96308723) q[3];
sx q[3];
rz(-1.77308) q[3];
sx q[3];
rz(-0.62455463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-3.0528089) q[0];
sx q[0];
rz(-1.4316906) q[0];
sx q[0];
rz(0.54876304) q[0];
rz(0.051963003) q[1];
sx q[1];
rz(-0.49182645) q[1];
sx q[1];
rz(-2.5351977) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.711395) q[0];
sx q[0];
rz(-1.0476307) q[0];
sx q[0];
rz(0.27336143) q[0];
x q[1];
rz(-0.90784448) q[2];
sx q[2];
rz(-1.2693451) q[2];
sx q[2];
rz(2.9482258) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4013306) q[1];
sx q[1];
rz(-1.4583424) q[1];
sx q[1];
rz(-2.3349891) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2783979) q[3];
sx q[3];
rz(-2.6226225) q[3];
sx q[3];
rz(-2.2591345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9397883) q[2];
sx q[2];
rz(-0.99210343) q[2];
sx q[2];
rz(-0.50841224) q[2];
rz(1.1172179) q[3];
sx q[3];
rz(-0.17337392) q[3];
sx q[3];
rz(1.6569051) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8023119) q[0];
sx q[0];
rz(-0.34582129) q[0];
sx q[0];
rz(-0.75396496) q[0];
rz(-2.1915961) q[1];
sx q[1];
rz(-1.6597304) q[1];
sx q[1];
rz(2.9139013) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4906625) q[0];
sx q[0];
rz(-2.7224446) q[0];
sx q[0];
rz(-1.7463669) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.44616206) q[2];
sx q[2];
rz(-1.237962) q[2];
sx q[2];
rz(1.7762426) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5880809) q[1];
sx q[1];
rz(-1.3462102) q[1];
sx q[1];
rz(0.54393804) q[1];
rz(-1.6119192) q[3];
sx q[3];
rz(-2.1293981) q[3];
sx q[3];
rz(2.5839644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0630539) q[2];
sx q[2];
rz(-2.9324014) q[2];
sx q[2];
rz(-0.83855808) q[2];
rz(2.7770384) q[3];
sx q[3];
rz(-1.3804599) q[3];
sx q[3];
rz(2.3013039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39695981) q[0];
sx q[0];
rz(-1.3752221) q[0];
sx q[0];
rz(-0.80378419) q[0];
rz(1.0546168) q[1];
sx q[1];
rz(-2.5024253) q[1];
sx q[1];
rz(1.9931591) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.046612) q[0];
sx q[0];
rz(-1.5469095) q[0];
sx q[0];
rz(1.6427965) q[0];
x q[1];
rz(0.9903332) q[2];
sx q[2];
rz(-2.5361984) q[2];
sx q[2];
rz(0.47547542) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0587522) q[1];
sx q[1];
rz(-1.0995004) q[1];
sx q[1];
rz(2.6575762) q[1];
rz(2.7896342) q[3];
sx q[3];
rz(-1.5378012) q[3];
sx q[3];
rz(1.9381423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0917197) q[2];
sx q[2];
rz(-2.1240978) q[2];
sx q[2];
rz(-0.77825528) q[2];
rz(-0.19566472) q[3];
sx q[3];
rz(-2.1019432) q[3];
sx q[3];
rz(-2.1197317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2739928) q[0];
sx q[0];
rz(-0.99853981) q[0];
sx q[0];
rz(2.8883949) q[0];
rz(-0.7397488) q[1];
sx q[1];
rz(-2.4052129) q[1];
sx q[1];
rz(0.69828066) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90726501) q[0];
sx q[0];
rz(-2.8303879) q[0];
sx q[0];
rz(-0.23647232) q[0];
x q[1];
rz(1.2487222) q[2];
sx q[2];
rz(-3.0634355) q[2];
sx q[2];
rz(-0.31101481) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0638949) q[1];
sx q[1];
rz(-2.6297036) q[1];
sx q[1];
rz(1.1179395) q[1];
rz(1.0032759) q[3];
sx q[3];
rz(-1.7756117) q[3];
sx q[3];
rz(1.8703465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0239821) q[2];
sx q[2];
rz(-2.2587946) q[2];
sx q[2];
rz(-2.3013766) q[2];
rz(-0.64030567) q[3];
sx q[3];
rz(-0.87602031) q[3];
sx q[3];
rz(-1.9742112) q[3];
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
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
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
rz(2.5393456) q[2];
sx q[2];
rz(-1.8613653) q[2];
sx q[2];
rz(-2.0801434) q[2];
rz(-1.7520262) q[3];
sx q[3];
rz(-2.443234) q[3];
sx q[3];
rz(0.20886226) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];