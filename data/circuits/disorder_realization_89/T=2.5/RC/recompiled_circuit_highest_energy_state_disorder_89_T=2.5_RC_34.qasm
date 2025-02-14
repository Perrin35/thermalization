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
rz(-2.1630218) q[0];
sx q[0];
rz(-3.0244654) q[0];
sx q[0];
rz(2.9168265) q[0];
rz(2.6449142) q[1];
sx q[1];
rz(4.7159046) q[1];
sx q[1];
rz(11.386303) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84070228) q[0];
sx q[0];
rz(-1.4288752) q[0];
sx q[0];
rz(3.0034756) q[0];
rz(-pi) q[1];
rz(-2.6981437) q[2];
sx q[2];
rz(-1.5649619) q[2];
sx q[2];
rz(1.3590517) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9182723) q[1];
sx q[1];
rz(-1.2967355) q[1];
sx q[1];
rz(2.5943701) q[1];
x q[2];
rz(1.3985083) q[3];
sx q[3];
rz(-0.43753657) q[3];
sx q[3];
rz(-3.1222536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4687389) q[2];
sx q[2];
rz(-1.3287975) q[2];
sx q[2];
rz(1.2057745) q[2];
rz(2.2852211) q[3];
sx q[3];
rz(-2.2090293) q[3];
sx q[3];
rz(-0.87730733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.4641651) q[0];
sx q[0];
rz(-2.7441315) q[0];
sx q[0];
rz(-2.933266) q[0];
rz(-1.2267998) q[1];
sx q[1];
rz(-1.0700763) q[1];
sx q[1];
rz(-0.56455451) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0923742) q[0];
sx q[0];
rz(-0.16516797) q[0];
sx q[0];
rz(-1.3081999) q[0];
rz(-pi) q[1];
x q[1];
rz(0.73123572) q[2];
sx q[2];
rz(-1.952716) q[2];
sx q[2];
rz(-2.3637091) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1543252) q[1];
sx q[1];
rz(-2.4134153) q[1];
sx q[1];
rz(2.8716692) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6030231) q[3];
sx q[3];
rz(-2.6015266) q[3];
sx q[3];
rz(-0.16656729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.63689256) q[2];
sx q[2];
rz(-1.8919287) q[2];
sx q[2];
rz(-1.4862109) q[2];
rz(-2.6460323) q[3];
sx q[3];
rz(-1.4080181) q[3];
sx q[3];
rz(2.1347617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14690742) q[0];
sx q[0];
rz(-2.4201396) q[0];
sx q[0];
rz(-1.3837234) q[0];
rz(1.9231298) q[1];
sx q[1];
rz(-2.0530901) q[1];
sx q[1];
rz(-0.4247492) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1871763) q[0];
sx q[0];
rz(-1.6023876) q[0];
sx q[0];
rz(0.0626465) q[0];
rz(2.9016657) q[2];
sx q[2];
rz(-1.3245965) q[2];
sx q[2];
rz(0.70084106) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.16097611) q[1];
sx q[1];
rz(-0.38432877) q[1];
sx q[1];
rz(-2.9381782) q[1];
rz(-pi) q[2];
rz(-3.0851641) q[3];
sx q[3];
rz(-2.5241969) q[3];
sx q[3];
rz(1.4790725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.72033) q[2];
sx q[2];
rz(-1.0676554) q[2];
sx q[2];
rz(-3.1171411) q[2];
rz(1.9931741) q[3];
sx q[3];
rz(-2.0423856) q[3];
sx q[3];
rz(2.059977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20632437) q[0];
sx q[0];
rz(-2.9117888) q[0];
sx q[0];
rz(-2.9277053) q[0];
rz(2.3230486) q[1];
sx q[1];
rz(-1.7844424) q[1];
sx q[1];
rz(-0.68665409) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7711597) q[0];
sx q[0];
rz(-2.4338182) q[0];
sx q[0];
rz(0.64162713) q[0];
rz(-pi) q[1];
rz(1.0222203) q[2];
sx q[2];
rz(-1.3201664) q[2];
sx q[2];
rz(-0.20285367) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.80685341) q[1];
sx q[1];
rz(-0.76419965) q[1];
sx q[1];
rz(0.13607009) q[1];
x q[2];
rz(-2.7044933) q[3];
sx q[3];
rz(-1.0954787) q[3];
sx q[3];
rz(3.0393485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0246747) q[2];
sx q[2];
rz(-1.7322098) q[2];
sx q[2];
rz(0.68946687) q[2];
rz(-1.8170478) q[3];
sx q[3];
rz(-1.9246512) q[3];
sx q[3];
rz(-2.569258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4631184) q[0];
sx q[0];
rz(-0.093570396) q[0];
sx q[0];
rz(-1.0372739) q[0];
rz(-2.4689238) q[1];
sx q[1];
rz(-1.172352) q[1];
sx q[1];
rz(-0.79576463) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0908431) q[0];
sx q[0];
rz(-1.2552869) q[0];
sx q[0];
rz(-1.0361703) q[0];
rz(-pi) q[1];
rz(-2.1433709) q[2];
sx q[2];
rz(-1.8279982) q[2];
sx q[2];
rz(0.57359475) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.63043609) q[1];
sx q[1];
rz(-1.2390343) q[1];
sx q[1];
rz(0.81960441) q[1];
rz(-pi) q[2];
rz(2.4733801) q[3];
sx q[3];
rz(-1.3283806) q[3];
sx q[3];
rz(-2.6072864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1611288) q[2];
sx q[2];
rz(-2.0420065) q[2];
sx q[2];
rz(2.3058092) q[2];
rz(-2.8849844) q[3];
sx q[3];
rz(-1.1107239) q[3];
sx q[3];
rz(-0.49032828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.68818727) q[0];
sx q[0];
rz(-1.1510993) q[0];
sx q[0];
rz(1.9292462) q[0];
rz(-2.1064099) q[1];
sx q[1];
rz(-1.6500708) q[1];
sx q[1];
rz(-1.0119247) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2741756) q[0];
sx q[0];
rz(-1.3998446) q[0];
sx q[0];
rz(-3.0311008) q[0];
rz(3.1306562) q[2];
sx q[2];
rz(-1.8746398) q[2];
sx q[2];
rz(0.86412175) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1862667) q[1];
sx q[1];
rz(-1.1623993) q[1];
sx q[1];
rz(-2.6751861) q[1];
rz(-pi) q[2];
x q[2];
rz(0.38326855) q[3];
sx q[3];
rz(-1.82194) q[3];
sx q[3];
rz(2.0519837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.47784352) q[2];
sx q[2];
rz(-0.19890824) q[2];
sx q[2];
rz(2.193023) q[2];
rz(2.6712724) q[3];
sx q[3];
rz(-2.0430653) q[3];
sx q[3];
rz(-1.8868014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5941641) q[0];
sx q[0];
rz(-1.0110039) q[0];
sx q[0];
rz(1.8773361) q[0];
rz(-0.7803548) q[1];
sx q[1];
rz(-0.64096132) q[1];
sx q[1];
rz(-2.9497214) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77593747) q[0];
sx q[0];
rz(-1.2167551) q[0];
sx q[0];
rz(-2.6440622) q[0];
rz(-pi) q[1];
rz(-1.948951) q[2];
sx q[2];
rz(-1.9072235) q[2];
sx q[2];
rz(-1.9371197) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.807232) q[1];
sx q[1];
rz(-1.8340115) q[1];
sx q[1];
rz(-0.46137793) q[1];
x q[2];
rz(2.170788) q[3];
sx q[3];
rz(-1.8626584) q[3];
sx q[3];
rz(-1.6455022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4912305) q[2];
sx q[2];
rz(-1.0543084) q[2];
sx q[2];
rz(0.80797705) q[2];
rz(-2.2982277) q[3];
sx q[3];
rz(-1.1542412) q[3];
sx q[3];
rz(-1.5850867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4845487) q[0];
sx q[0];
rz(-3.1175685) q[0];
sx q[0];
rz(0.76883823) q[0];
rz(2.8688042) q[1];
sx q[1];
rz(-1.6285248) q[1];
sx q[1];
rz(1.3974894) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5114355) q[0];
sx q[0];
rz(-2.6574597) q[0];
sx q[0];
rz(0.70229437) q[0];
x q[1];
rz(2.8238472) q[2];
sx q[2];
rz(-2.0644798) q[2];
sx q[2];
rz(-1.2682481) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8694526) q[1];
sx q[1];
rz(-1.4542129) q[1];
sx q[1];
rz(0.061007331) q[1];
rz(-pi) q[2];
rz(0.66889735) q[3];
sx q[3];
rz(-1.2227158) q[3];
sx q[3];
rz(2.477248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.93691319) q[2];
sx q[2];
rz(-2.0759373) q[2];
sx q[2];
rz(1.5020471) q[2];
rz(-1.8223193) q[3];
sx q[3];
rz(-0.84984142) q[3];
sx q[3];
rz(-1.5523065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69916344) q[0];
sx q[0];
rz(-0.84446877) q[0];
sx q[0];
rz(2.1635639) q[0];
rz(-0.49081048) q[1];
sx q[1];
rz(-1.4421137) q[1];
sx q[1];
rz(1.4960272) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1745047) q[0];
sx q[0];
rz(-1.3347111) q[0];
sx q[0];
rz(1.697879) q[0];
rz(1.3991576) q[2];
sx q[2];
rz(-1.1632068) q[2];
sx q[2];
rz(1.7625347) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5514976) q[1];
sx q[1];
rz(-1.2480772) q[1];
sx q[1];
rz(1.3151874) q[1];
rz(1.4946901) q[3];
sx q[3];
rz(-1.0400912) q[3];
sx q[3];
rz(2.2709341) q[3];
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
rz(-0.36618048) q[3];
sx q[3];
rz(-1.4270695) q[3];
sx q[3];
rz(-2.0100994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91728297) q[0];
sx q[0];
rz(-1.5320822) q[0];
sx q[0];
rz(-3.068058) q[0];
rz(-1.3167943) q[1];
sx q[1];
rz(-1.2578745) q[1];
sx q[1];
rz(-1.8427461) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5054437) q[0];
sx q[0];
rz(-0.33056405) q[0];
sx q[0];
rz(-2.5695468) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0082874) q[2];
sx q[2];
rz(-2.0553737) q[2];
sx q[2];
rz(2.7078748) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4529113) q[1];
sx q[1];
rz(-1.9454524) q[1];
sx q[1];
rz(-2.8883347) q[1];
rz(-0.73963005) q[3];
sx q[3];
rz(-1.4907903) q[3];
sx q[3];
rz(-2.7440939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1765882) q[2];
sx q[2];
rz(-0.95866385) q[2];
sx q[2];
rz(-0.57903543) q[2];
rz(1.5163007) q[3];
sx q[3];
rz(-2.5987891) q[3];
sx q[3];
rz(1.6453086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8662921) q[0];
sx q[0];
rz(-1.0466812) q[0];
sx q[0];
rz(-0.064591797) q[0];
rz(-0.91118377) q[1];
sx q[1];
rz(-1.3180399) q[1];
sx q[1];
rz(-1.3175189) q[1];
rz(-3.1338794) q[2];
sx q[2];
rz(-2.7174502) q[2];
sx q[2];
rz(2.9677718) q[2];
rz(-2.1880423) q[3];
sx q[3];
rz(-1.336477) q[3];
sx q[3];
rz(-0.82930641) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
