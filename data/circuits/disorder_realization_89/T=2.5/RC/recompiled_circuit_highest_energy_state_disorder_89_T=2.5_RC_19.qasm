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
rz(0.97857082) q[0];
sx q[0];
rz(-0.11712722) q[0];
sx q[0];
rz(0.22476619) q[0];
rz(-0.4966785) q[1];
sx q[1];
rz(-1.574312) q[1];
sx q[1];
rz(1.1800676) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71043541) q[0];
sx q[0];
rz(-1.4340766) q[0];
sx q[0];
rz(-1.4275292) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1279951) q[2];
sx q[2];
rz(-0.44348479) q[2];
sx q[2];
rz(2.9421303) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5109053) q[1];
sx q[1];
rz(-1.0461591) q[1];
sx q[1];
rz(-1.252763) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1389569) q[3];
sx q[3];
rz(-1.4980928) q[3];
sx q[3];
rz(1.3951226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4687389) q[2];
sx q[2];
rz(-1.8127952) q[2];
sx q[2];
rz(-1.9358181) q[2];
rz(0.85637158) q[3];
sx q[3];
rz(-0.93256336) q[3];
sx q[3];
rz(2.2642853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67742753) q[0];
sx q[0];
rz(-2.7441315) q[0];
sx q[0];
rz(2.933266) q[0];
rz(-1.9147929) q[1];
sx q[1];
rz(-1.0700763) q[1];
sx q[1];
rz(-2.5770381) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3608332) q[0];
sx q[0];
rz(-1.6134904) q[0];
sx q[0];
rz(-1.7303995) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.54148285) q[2];
sx q[2];
rz(-0.80831203) q[2];
sx q[2];
rz(-1.186651) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3422505) q[1];
sx q[1];
rz(-2.2671351) q[1];
sx q[1];
rz(-1.8042248) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1222778) q[3];
sx q[3];
rz(-1.0310415) q[3];
sx q[3];
rz(-0.12899736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.63689256) q[2];
sx q[2];
rz(-1.2496639) q[2];
sx q[2];
rz(1.4862109) q[2];
rz(0.49556035) q[3];
sx q[3];
rz(-1.4080181) q[3];
sx q[3];
rz(2.1347617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14690742) q[0];
sx q[0];
rz(-0.72145307) q[0];
sx q[0];
rz(-1.7578693) q[0];
rz(-1.9231298) q[1];
sx q[1];
rz(-1.0885025) q[1];
sx q[1];
rz(-0.4247492) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61439862) q[0];
sx q[0];
rz(-1.6334115) q[0];
sx q[0];
rz(-1.539143) q[0];
rz(-0.23992691) q[2];
sx q[2];
rz(-1.8169962) q[2];
sx q[2];
rz(2.4407516) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.16097611) q[1];
sx q[1];
rz(-2.7572639) q[1];
sx q[1];
rz(-2.9381782) q[1];
x q[2];
rz(0.61664403) q[3];
sx q[3];
rz(-1.6034521) q[3];
sx q[3];
rz(-0.045696229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4212627) q[2];
sx q[2];
rz(-2.0739372) q[2];
sx q[2];
rz(3.1171411) q[2];
rz(1.9931741) q[3];
sx q[3];
rz(-1.0992071) q[3];
sx q[3];
rz(1.0816157) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9352683) q[0];
sx q[0];
rz(-2.9117888) q[0];
sx q[0];
rz(0.21388737) q[0];
rz(2.3230486) q[1];
sx q[1];
rz(-1.7844424) q[1];
sx q[1];
rz(-0.68665409) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.370433) q[0];
sx q[0];
rz(-2.4338182) q[0];
sx q[0];
rz(0.64162713) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0222203) q[2];
sx q[2];
rz(-1.8214263) q[2];
sx q[2];
rz(0.20285367) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.61942809) q[1];
sx q[1];
rz(-0.81541895) q[1];
sx q[1];
rz(-1.700089) q[1];
rz(-pi) q[2];
x q[2];
rz(0.43709932) q[3];
sx q[3];
rz(-1.0954787) q[3];
sx q[3];
rz(3.0393485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.116918) q[2];
sx q[2];
rz(-1.7322098) q[2];
sx q[2];
rz(-0.68946687) q[2];
rz(1.3245448) q[3];
sx q[3];
rz(-1.2169415) q[3];
sx q[3];
rz(2.569258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4631184) q[0];
sx q[0];
rz(-3.0480223) q[0];
sx q[0];
rz(1.0372739) q[0];
rz(-2.4689238) q[1];
sx q[1];
rz(-1.9692407) q[1];
sx q[1];
rz(0.79576463) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8032629) q[0];
sx q[0];
rz(-2.0764377) q[0];
sx q[0];
rz(0.36257669) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1188385) q[2];
sx q[2];
rz(-0.6217494) q[2];
sx q[2];
rz(-2.5202519) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9059674) q[1];
sx q[1];
rz(-0.86938953) q[1];
sx q[1];
rz(2.7011306) q[1];
rz(-pi) q[2];
rz(0.37973399) q[3];
sx q[3];
rz(-2.437161) q[3];
sx q[3];
rz(-1.3317127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.98046389) q[2];
sx q[2];
rz(-2.0420065) q[2];
sx q[2];
rz(0.83578342) q[2];
rz(-2.8849844) q[3];
sx q[3];
rz(-1.1107239) q[3];
sx q[3];
rz(-0.49032828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68818727) q[0];
sx q[0];
rz(-1.9904933) q[0];
sx q[0];
rz(1.9292462) q[0];
rz(-2.1064099) q[1];
sx q[1];
rz(-1.6500708) q[1];
sx q[1];
rz(2.129668) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8674171) q[0];
sx q[0];
rz(-1.3998446) q[0];
sx q[0];
rz(0.11049185) q[0];
rz(1.2669358) q[2];
sx q[2];
rz(-1.5603608) q[2];
sx q[2];
rz(-0.70994678) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3287836) q[1];
sx q[1];
rz(-1.145383) q[1];
sx q[1];
rz(-2.0219457) q[1];
rz(-0.38326855) q[3];
sx q[3];
rz(-1.3196527) q[3];
sx q[3];
rz(-1.0896089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6637491) q[2];
sx q[2];
rz(-0.19890824) q[2];
sx q[2];
rz(-2.193023) q[2];
rz(-0.47032022) q[3];
sx q[3];
rz(-1.0985273) q[3];
sx q[3];
rz(1.8868014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5941641) q[0];
sx q[0];
rz(-2.1305888) q[0];
sx q[0];
rz(-1.8773361) q[0];
rz(-2.3612379) q[1];
sx q[1];
rz(-2.5006313) q[1];
sx q[1];
rz(-2.9497214) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.980968) q[0];
sx q[0];
rz(-2.0349398) q[0];
sx q[0];
rz(1.1726517) q[0];
rz(-2.7816747) q[2];
sx q[2];
rz(-1.9267757) q[2];
sx q[2];
rz(2.6448665) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3650771) q[1];
sx q[1];
rz(-1.1264701) q[1];
sx q[1];
rz(1.278484) q[1];
x q[2];
rz(2.0597919) q[3];
sx q[3];
rz(-2.4823175) q[3];
sx q[3];
rz(2.6687572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.65036217) q[2];
sx q[2];
rz(-2.0872842) q[2];
sx q[2];
rz(2.3336156) q[2];
rz(0.84336495) q[3];
sx q[3];
rz(-1.9873514) q[3];
sx q[3];
rz(1.5850867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4845487) q[0];
sx q[0];
rz(-3.1175685) q[0];
sx q[0];
rz(-2.3727544) q[0];
rz(2.8688042) q[1];
sx q[1];
rz(-1.6285248) q[1];
sx q[1];
rz(1.3974894) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13279877) q[0];
sx q[0];
rz(-1.2075612) q[0];
sx q[0];
rz(-1.2433267) q[0];
x q[1];
rz(2.8238472) q[2];
sx q[2];
rz(-1.0771128) q[2];
sx q[2];
rz(-1.8733446) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8500415) q[1];
sx q[1];
rz(-1.631389) q[1];
sx q[1];
rz(1.4539976) q[1];
x q[2];
rz(-0.52941771) q[3];
sx q[3];
rz(-2.4000958) q[3];
sx q[3];
rz(1.3138258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.93691319) q[2];
sx q[2];
rz(-1.0656554) q[2];
sx q[2];
rz(1.6395456) q[2];
rz(1.3192734) q[3];
sx q[3];
rz(-0.84984142) q[3];
sx q[3];
rz(-1.5523065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4424292) q[0];
sx q[0];
rz(-0.84446877) q[0];
sx q[0];
rz(0.97802877) q[0];
rz(-0.49081048) q[1];
sx q[1];
rz(-1.699479) q[1];
sx q[1];
rz(-1.4960272) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.467062) q[0];
sx q[0];
rz(-0.26755324) q[0];
sx q[0];
rz(2.6567099) q[0];
rz(-1.7424351) q[2];
sx q[2];
rz(-1.9783859) q[2];
sx q[2];
rz(-1.7625347) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5514976) q[1];
sx q[1];
rz(-1.2480772) q[1];
sx q[1];
rz(-1.8264052) q[1];
x q[2];
rz(1.6469025) q[3];
sx q[3];
rz(-1.0400912) q[3];
sx q[3];
rz(0.87065856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.28828037) q[2];
sx q[2];
rz(-0.86709443) q[2];
sx q[2];
rz(1.8193998) q[2];
rz(0.36618048) q[3];
sx q[3];
rz(-1.7145232) q[3];
sx q[3];
rz(1.1314932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2243097) q[0];
sx q[0];
rz(-1.6095105) q[0];
sx q[0];
rz(-0.073534615) q[0];
rz(-1.3167943) q[1];
sx q[1];
rz(-1.2578745) q[1];
sx q[1];
rz(1.2988466) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.038485) q[0];
sx q[0];
rz(-1.2943876) q[0];
sx q[0];
rz(-1.3871219) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3234947) q[2];
sx q[2];
rz(-2.6404233) q[2];
sx q[2];
rz(0.71401087) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.83734456) q[1];
sx q[1];
rz(-2.6927444) q[1];
sx q[1];
rz(2.138109) q[1];
x q[2];
rz(-2.4019626) q[3];
sx q[3];
rz(-1.4907903) q[3];
sx q[3];
rz(-0.39749872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1765882) q[2];
sx q[2];
rz(-0.95866385) q[2];
sx q[2];
rz(-2.5625572) q[2];
rz(-1.5163007) q[3];
sx q[3];
rz(-2.5987891) q[3];
sx q[3];
rz(1.496284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(1.2753006) q[0];
sx q[0];
rz(-1.0466812) q[0];
sx q[0];
rz(-0.064591797) q[0];
rz(0.91118377) q[1];
sx q[1];
rz(-1.8235527) q[1];
sx q[1];
rz(1.8240737) q[1];
rz(-1.5742792) q[2];
sx q[2];
rz(-1.9949253) q[2];
sx q[2];
rz(-0.16535769) q[2];
rz(-1.961962) q[3];
sx q[3];
rz(-0.65476553) q[3];
sx q[3];
rz(-2.7162566) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
