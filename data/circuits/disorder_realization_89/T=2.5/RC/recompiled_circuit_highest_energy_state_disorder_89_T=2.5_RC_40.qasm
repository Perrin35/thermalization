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
rz(-2.9168265) q[0];
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
rz(0.84070228) q[0];
sx q[0];
rz(-1.4288752) q[0];
sx q[0];
rz(-0.13811708) q[0];
rz(-pi) q[1];
rz(-3.1279951) q[2];
sx q[2];
rz(-0.44348479) q[2];
sx q[2];
rz(2.9421303) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5109053) q[1];
sx q[1];
rz(-1.0461591) q[1];
sx q[1];
rz(1.252763) q[1];
rz(0.080022607) q[3];
sx q[3];
rz(-1.140174) q[3];
sx q[3];
rz(-2.9324556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.67285377) q[2];
sx q[2];
rz(-1.3287975) q[2];
sx q[2];
rz(-1.2057745) q[2];
rz(-0.85637158) q[3];
sx q[3];
rz(-2.2090293) q[3];
sx q[3];
rz(-0.87730733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4641651) q[0];
sx q[0];
rz(-2.7441315) q[0];
sx q[0];
rz(-2.933266) q[0];
rz(1.2267998) q[1];
sx q[1];
rz(-2.0715163) q[1];
sx q[1];
rz(-0.56455451) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3584261) q[0];
sx q[0];
rz(-1.7302528) q[0];
sx q[0];
rz(0.043242976) q[0];
rz(-pi) q[1];
rz(-2.0656086) q[2];
sx q[2];
rz(-2.2391265) q[2];
sx q[2];
rz(2.6713616) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.98726749) q[1];
sx q[1];
rz(-0.72817737) q[1];
sx q[1];
rz(-0.2699234) q[1];
x q[2];
rz(3.1222778) q[3];
sx q[3];
rz(-2.1105511) q[3];
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
rz(-2.6460323) q[3];
sx q[3];
rz(-1.4080181) q[3];
sx q[3];
rz(2.1347617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14690742) q[0];
sx q[0];
rz(-0.72145307) q[0];
sx q[0];
rz(1.3837234) q[0];
rz(-1.9231298) q[1];
sx q[1];
rz(-1.0885025) q[1];
sx q[1];
rz(-0.4247492) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0587416) q[0];
sx q[0];
rz(-0.070151873) q[0];
sx q[0];
rz(0.46746107) q[0];
rz(-pi) q[1];
x q[1];
rz(0.23992691) q[2];
sx q[2];
rz(-1.8169962) q[2];
sx q[2];
rz(0.70084106) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9207117) q[1];
sx q[1];
rz(-1.6466116) q[1];
sx q[1];
rz(0.37714173) q[1];
rz(0.61664403) q[3];
sx q[3];
rz(-1.5381406) q[3];
sx q[3];
rz(-3.0958964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.72033) q[2];
sx q[2];
rz(-2.0739372) q[2];
sx q[2];
rz(-3.1171411) q[2];
rz(1.1484185) q[3];
sx q[3];
rz(-2.0423856) q[3];
sx q[3];
rz(1.0816157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9352683) q[0];
sx q[0];
rz(-2.9117888) q[0];
sx q[0];
rz(0.21388737) q[0];
rz(-0.81854406) q[1];
sx q[1];
rz(-1.3571502) q[1];
sx q[1];
rz(0.68665409) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4249466) q[0];
sx q[0];
rz(-1.970463) q[0];
sx q[0];
rz(2.5406688) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0271717) q[2];
sx q[2];
rz(-0.59774146) q[2];
sx q[2];
rz(-0.98243143) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.80685341) q[1];
sx q[1];
rz(-0.76419965) q[1];
sx q[1];
rz(-0.13607009) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0874168) q[3];
sx q[3];
rz(-1.9566907) q[3];
sx q[3];
rz(1.8836881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0246747) q[2];
sx q[2];
rz(-1.7322098) q[2];
sx q[2];
rz(-2.4521258) q[2];
rz(1.8170478) q[3];
sx q[3];
rz(-1.9246512) q[3];
sx q[3];
rz(-0.57233468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-1.6784742) q[0];
sx q[0];
rz(-3.0480223) q[0];
sx q[0];
rz(-2.1043188) q[0];
rz(0.67266881) q[1];
sx q[1];
rz(-1.172352) q[1];
sx q[1];
rz(-0.79576463) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33832975) q[0];
sx q[0];
rz(-1.065155) q[0];
sx q[0];
rz(-2.779016) q[0];
x q[1];
rz(2.0227541) q[2];
sx q[2];
rz(-2.5198433) q[2];
sx q[2];
rz(-2.5202519) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5111566) q[1];
sx q[1];
rz(-1.9025584) q[1];
sx q[1];
rz(2.3219882) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8759868) q[3];
sx q[3];
rz(-2.2161336) q[3];
sx q[3];
rz(-0.84922817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1611288) q[2];
sx q[2];
rz(-1.0995862) q[2];
sx q[2];
rz(-2.3058092) q[2];
rz(-2.8849844) q[3];
sx q[3];
rz(-1.1107239) q[3];
sx q[3];
rz(2.6512644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4534054) q[0];
sx q[0];
rz(-1.9904933) q[0];
sx q[0];
rz(1.9292462) q[0];
rz(2.1064099) q[1];
sx q[1];
rz(-1.6500708) q[1];
sx q[1];
rz(1.0119247) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8674171) q[0];
sx q[0];
rz(-1.7417481) q[0];
sx q[0];
rz(0.11049185) q[0];
rz(-0.010936485) q[2];
sx q[2];
rz(-1.8746398) q[2];
sx q[2];
rz(-2.2774709) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.95532596) q[1];
sx q[1];
rz(-1.1623993) q[1];
sx q[1];
rz(-0.46640654) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5402734) q[3];
sx q[3];
rz(-0.45479247) q[3];
sx q[3];
rz(-3.0704344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6637491) q[2];
sx q[2];
rz(-2.9426844) q[2];
sx q[2];
rz(2.193023) q[2];
rz(-2.6712724) q[3];
sx q[3];
rz(-1.0985273) q[3];
sx q[3];
rz(1.2547913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5941641) q[0];
sx q[0];
rz(-2.1305888) q[0];
sx q[0];
rz(-1.2642566) q[0];
rz(-2.3612379) q[1];
sx q[1];
rz(-0.64096132) q[1];
sx q[1];
rz(2.9497214) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.980968) q[0];
sx q[0];
rz(-2.0349398) q[0];
sx q[0];
rz(1.1726517) q[0];
rz(-pi) q[1];
rz(-0.81249313) q[2];
sx q[2];
rz(-0.50069649) q[2];
sx q[2];
rz(-2.814584) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.75429186) q[1];
sx q[1];
rz(-2.6151492) q[1];
sx q[1];
rz(-2.597288) q[1];
rz(-pi) q[2];
rz(-2.7924813) q[3];
sx q[3];
rz(-2.1421332) q[3];
sx q[3];
rz(3.0219363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4912305) q[2];
sx q[2];
rz(-1.0543084) q[2];
sx q[2];
rz(0.80797705) q[2];
rz(-0.84336495) q[3];
sx q[3];
rz(-1.1542412) q[3];
sx q[3];
rz(-1.556506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.657044) q[0];
sx q[0];
rz(-3.1175685) q[0];
sx q[0];
rz(-0.76883823) q[0];
rz(-0.27278849) q[1];
sx q[1];
rz(-1.6285248) q[1];
sx q[1];
rz(-1.7441033) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5834806) q[0];
sx q[0];
rz(-1.2654103) q[0];
sx q[0];
rz(-2.7598513) q[0];
rz(-pi) q[1];
rz(-0.31774546) q[2];
sx q[2];
rz(-2.0644798) q[2];
sx q[2];
rz(1.8733446) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8694526) q[1];
sx q[1];
rz(-1.4542129) q[1];
sx q[1];
rz(-0.061007331) q[1];
x q[2];
rz(0.66889735) q[3];
sx q[3];
rz(-1.9188768) q[3];
sx q[3];
rz(-2.477248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.93691319) q[2];
sx q[2];
rz(-1.0656554) q[2];
sx q[2];
rz(1.6395456) q[2];
rz(-1.3192734) q[3];
sx q[3];
rz(-2.2917512) q[3];
sx q[3];
rz(-1.5523065) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69916344) q[0];
sx q[0];
rz(-0.84446877) q[0];
sx q[0];
rz(-0.97802877) q[0];
rz(0.49081048) q[1];
sx q[1];
rz(-1.699479) q[1];
sx q[1];
rz(1.4960272) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.967088) q[0];
sx q[0];
rz(-1.3347111) q[0];
sx q[0];
rz(-1.4437136) q[0];
rz(1.7424351) q[2];
sx q[2];
rz(-1.9783859) q[2];
sx q[2];
rz(1.7625347) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2792932) q[1];
sx q[1];
rz(-0.40888834) q[1];
sx q[1];
rz(-2.4942231) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4946901) q[3];
sx q[3];
rz(-1.0400912) q[3];
sx q[3];
rz(-0.87065856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8533123) q[2];
sx q[2];
rz(-2.2744982) q[2];
sx q[2];
rz(1.3221928) q[2];
rz(-0.36618048) q[3];
sx q[3];
rz(-1.7145232) q[3];
sx q[3];
rz(-1.1314932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91728297) q[0];
sx q[0];
rz(-1.5320822) q[0];
sx q[0];
rz(0.073534615) q[0];
rz(-1.8247983) q[1];
sx q[1];
rz(-1.2578745) q[1];
sx q[1];
rz(-1.2988466) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6361489) q[0];
sx q[0];
rz(-2.8110286) q[0];
sx q[0];
rz(-2.5695468) q[0];
rz(3.0082874) q[2];
sx q[2];
rz(-2.0553737) q[2];
sx q[2];
rz(2.7078748) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6886814) q[1];
sx q[1];
rz(-1.1961403) q[1];
sx q[1];
rz(2.8883347) q[1];
x q[2];
rz(-0.73963005) q[3];
sx q[3];
rz(-1.4907903) q[3];
sx q[3];
rz(0.39749872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1765882) q[2];
sx q[2];
rz(-2.1829288) q[2];
sx q[2];
rz(-0.57903543) q[2];
rz(1.5163007) q[3];
sx q[3];
rz(-0.54280353) q[3];
sx q[3];
rz(-1.6453086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-2.7174614) q[2];
sx q[2];
rz(-1.5739706) q[2];
sx q[2];
rz(-1.7375873) q[2];
rz(-1.1796307) q[3];
sx q[3];
rz(-2.4868271) q[3];
sx q[3];
rz(0.42533608) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
