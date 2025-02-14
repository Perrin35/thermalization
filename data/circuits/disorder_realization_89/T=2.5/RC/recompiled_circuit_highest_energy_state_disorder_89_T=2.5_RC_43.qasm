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
rz(3.0244654) q[0];
sx q[0];
rz(9.2000118) q[0];
rz(-0.4966785) q[1];
sx q[1];
rz(-1.574312) q[1];
sx q[1];
rz(1.1800676) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71043541) q[0];
sx q[0];
rz(-1.707516) q[0];
sx q[0];
rz(1.7140634) q[0];
x q[1];
rz(-3.1279951) q[2];
sx q[2];
rz(-2.6981079) q[2];
sx q[2];
rz(-2.9421303) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6306873) q[1];
sx q[1];
rz(-2.0954335) q[1];
sx q[1];
rz(1.8888297) q[1];
rz(-pi) q[2];
rz(-1.1389569) q[3];
sx q[3];
rz(-1.4980928) q[3];
sx q[3];
rz(1.3951226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4687389) q[2];
sx q[2];
rz(-1.3287975) q[2];
sx q[2];
rz(1.9358181) q[2];
rz(0.85637158) q[3];
sx q[3];
rz(-2.2090293) q[3];
sx q[3];
rz(-2.2642853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.67742753) q[0];
sx q[0];
rz(-2.7441315) q[0];
sx q[0];
rz(-2.933266) q[0];
rz(-1.9147929) q[1];
sx q[1];
rz(-1.0700763) q[1];
sx q[1];
rz(-2.5770381) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0923742) q[0];
sx q[0];
rz(-2.9764247) q[0];
sx q[0];
rz(-1.3081999) q[0];
rz(-2.0656086) q[2];
sx q[2];
rz(-0.90246614) q[2];
sx q[2];
rz(0.47023103) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1543252) q[1];
sx q[1];
rz(-2.4134153) q[1];
sx q[1];
rz(-2.8716692) q[1];
rz(-3.1222778) q[3];
sx q[3];
rz(-1.0310415) q[3];
sx q[3];
rz(3.0125953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5047001) q[2];
sx q[2];
rz(-1.2496639) q[2];
sx q[2];
rz(1.4862109) q[2];
rz(-2.6460323) q[3];
sx q[3];
rz(-1.7335745) q[3];
sx q[3];
rz(-2.1347617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9946852) q[0];
sx q[0];
rz(-2.4201396) q[0];
sx q[0];
rz(-1.7578693) q[0];
rz(1.2184628) q[1];
sx q[1];
rz(-1.0885025) q[1];
sx q[1];
rz(-0.4247492) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95441636) q[0];
sx q[0];
rz(-1.6023876) q[0];
sx q[0];
rz(-0.0626465) q[0];
rz(-pi) q[1];
rz(0.23992691) q[2];
sx q[2];
rz(-1.8169962) q[2];
sx q[2];
rz(0.70084106) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.220881) q[1];
sx q[1];
rz(-1.6466116) q[1];
sx q[1];
rz(-0.37714173) q[1];
rz(-1.6108171) q[3];
sx q[3];
rz(-0.95453018) q[3];
sx q[3];
rz(-1.54824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.72033) q[2];
sx q[2];
rz(-2.0739372) q[2];
sx q[2];
rz(0.024451582) q[2];
rz(-1.9931741) q[3];
sx q[3];
rz(-2.0423856) q[3];
sx q[3];
rz(-2.059977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(2.9352683) q[0];
sx q[0];
rz(-2.9117888) q[0];
sx q[0];
rz(2.9277053) q[0];
rz(-0.81854406) q[1];
sx q[1];
rz(-1.3571502) q[1];
sx q[1];
rz(-2.4549386) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59348561) q[0];
sx q[0];
rz(-1.0229551) q[0];
sx q[0];
rz(-1.0975029) q[0];
rz(2.0271717) q[2];
sx q[2];
rz(-2.5438512) q[2];
sx q[2];
rz(-2.1591612) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.86246762) q[1];
sx q[1];
rz(-1.4767932) q[1];
sx q[1];
rz(-0.75956168) q[1];
rz(-pi) q[2];
rz(-0.43709932) q[3];
sx q[3];
rz(-1.0954787) q[3];
sx q[3];
rz(-3.0393485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.116918) q[2];
sx q[2];
rz(-1.7322098) q[2];
sx q[2];
rz(0.68946687) q[2];
rz(1.8170478) q[3];
sx q[3];
rz(-1.2169415) q[3];
sx q[3];
rz(-2.569258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6784742) q[0];
sx q[0];
rz(-3.0480223) q[0];
sx q[0];
rz(1.0372739) q[0];
rz(2.4689238) q[1];
sx q[1];
rz(-1.172352) q[1];
sx q[1];
rz(-2.345828) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0507495) q[0];
sx q[0];
rz(-1.8863057) q[0];
sx q[0];
rz(1.0361703) q[0];
x q[1];
rz(-2.8383083) q[2];
sx q[2];
rz(-1.0192843) q[2];
sx q[2];
rz(1.1597275) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5111566) q[1];
sx q[1];
rz(-1.9025584) q[1];
sx q[1];
rz(0.81960441) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4733801) q[3];
sx q[3];
rz(-1.3283806) q[3];
sx q[3];
rz(-2.6072864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1611288) q[2];
sx q[2];
rz(-1.0995862) q[2];
sx q[2];
rz(-2.3058092) q[2];
rz(-2.8849844) q[3];
sx q[3];
rz(-2.0308688) q[3];
sx q[3];
rz(-2.6512644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68818727) q[0];
sx q[0];
rz(-1.1510993) q[0];
sx q[0];
rz(1.9292462) q[0];
rz(-2.1064099) q[1];
sx q[1];
rz(-1.4915219) q[1];
sx q[1];
rz(1.0119247) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31549227) q[0];
sx q[0];
rz(-1.6796711) q[0];
sx q[0];
rz(-1.7427765) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2669358) q[2];
sx q[2];
rz(-1.5812318) q[2];
sx q[2];
rz(2.4316459) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.95532596) q[1];
sx q[1];
rz(-1.1623993) q[1];
sx q[1];
rz(-2.6751861) q[1];
rz(-pi) q[2];
rz(1.3009146) q[3];
sx q[3];
rz(-1.2001451) q[3];
sx q[3];
rz(2.5605367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6637491) q[2];
sx q[2];
rz(-2.9426844) q[2];
sx q[2];
rz(0.94856962) q[2];
rz(-2.6712724) q[3];
sx q[3];
rz(-1.0985273) q[3];
sx q[3];
rz(1.2547913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54742852) q[0];
sx q[0];
rz(-2.1305888) q[0];
sx q[0];
rz(1.2642566) q[0];
rz(-0.7803548) q[1];
sx q[1];
rz(-0.64096132) q[1];
sx q[1];
rz(-2.9497214) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.980968) q[0];
sx q[0];
rz(-1.1066529) q[0];
sx q[0];
rz(-1.1726517) q[0];
x q[1];
rz(-1.1926417) q[2];
sx q[2];
rz(-1.2343692) q[2];
sx q[2];
rz(1.2044729) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.33436066) q[1];
sx q[1];
rz(-1.8340115) q[1];
sx q[1];
rz(-2.6802147) q[1];
rz(-pi) q[2];
rz(0.97080462) q[3];
sx q[3];
rz(-1.2789342) q[3];
sx q[3];
rz(1.4960904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.65036217) q[2];
sx q[2];
rz(-2.0872842) q[2];
sx q[2];
rz(-2.3336156) q[2];
rz(0.84336495) q[3];
sx q[3];
rz(-1.1542412) q[3];
sx q[3];
rz(-1.5850867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4845487) q[0];
sx q[0];
rz(-3.1175685) q[0];
sx q[0];
rz(2.3727544) q[0];
rz(2.8688042) q[1];
sx q[1];
rz(-1.6285248) q[1];
sx q[1];
rz(-1.7441033) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0087939) q[0];
sx q[0];
rz(-1.9340314) q[0];
sx q[0];
rz(-1.2433267) q[0];
rz(-pi) q[1];
rz(-1.0553841) q[2];
sx q[2];
rz(-1.8495108) q[2];
sx q[2];
rz(-2.9936522) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3859017) q[1];
sx q[1];
rz(-0.13151691) q[1];
sx q[1];
rz(2.0507858) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.66889735) q[3];
sx q[3];
rz(-1.2227158) q[3];
sx q[3];
rz(-2.477248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.93691319) q[2];
sx q[2];
rz(-2.0759373) q[2];
sx q[2];
rz(1.5020471) q[2];
rz(-1.3192734) q[3];
sx q[3];
rz(-2.2917512) q[3];
sx q[3];
rz(1.5892861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(-2.4424292) q[0];
sx q[0];
rz(-0.84446877) q[0];
sx q[0];
rz(2.1635639) q[0];
rz(-2.6507822) q[1];
sx q[1];
rz(-1.699479) q[1];
sx q[1];
rz(1.4960272) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.467062) q[0];
sx q[0];
rz(-0.26755324) q[0];
sx q[0];
rz(2.6567099) q[0];
rz(-pi) q[1];
rz(-0.37668682) q[2];
sx q[2];
rz(-0.44038195) q[2];
sx q[2];
rz(0.96681606) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.89801187) q[1];
sx q[1];
rz(-1.812938) q[1];
sx q[1];
rz(-2.8088074) q[1];
rz(-pi) q[2];
rz(1.4946901) q[3];
sx q[3];
rz(-2.1015014) q[3];
sx q[3];
rz(0.87065856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8533123) q[2];
sx q[2];
rz(-0.86709443) q[2];
sx q[2];
rz(-1.8193998) q[2];
rz(0.36618048) q[3];
sx q[3];
rz(-1.4270695) q[3];
sx q[3];
rz(-1.1314932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91728297) q[0];
sx q[0];
rz(-1.6095105) q[0];
sx q[0];
rz(0.073534615) q[0];
rz(-1.8247983) q[1];
sx q[1];
rz(-1.8837181) q[1];
sx q[1];
rz(-1.8427461) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6361489) q[0];
sx q[0];
rz(-2.8110286) q[0];
sx q[0];
rz(-0.57204582) q[0];
rz(-pi) q[1];
rz(2.0590563) q[2];
sx q[2];
rz(-1.4529145) q[2];
sx q[2];
rz(-1.0746909) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.21231743) q[1];
sx q[1];
rz(-1.3354509) q[1];
sx q[1];
rz(-1.1849684) q[1];
rz(-pi) q[2];
rz(2.4019626) q[3];
sx q[3];
rz(-1.6508023) q[3];
sx q[3];
rz(2.7440939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9650044) q[2];
sx q[2];
rz(-2.1829288) q[2];
sx q[2];
rz(-2.5625572) q[2];
rz(-1.5163007) q[3];
sx q[3];
rz(-0.54280353) q[3];
sx q[3];
rz(-1.496284) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8662921) q[0];
sx q[0];
rz(-1.0466812) q[0];
sx q[0];
rz(-0.064591797) q[0];
rz(2.2304089) q[1];
sx q[1];
rz(-1.3180399) q[1];
sx q[1];
rz(-1.3175189) q[1];
rz(-0.42413128) q[2];
sx q[2];
rz(-1.567622) q[2];
sx q[2];
rz(1.4040053) q[2];
rz(-0.95355036) q[3];
sx q[3];
rz(-1.8051157) q[3];
sx q[3];
rz(2.3122862) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
