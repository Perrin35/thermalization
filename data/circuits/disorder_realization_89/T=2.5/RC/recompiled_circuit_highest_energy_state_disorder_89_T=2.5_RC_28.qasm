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
rz(-0.4966785) q[1];
sx q[1];
rz(-1.574312) q[1];
sx q[1];
rz(-1.9615251) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4311572) q[0];
sx q[0];
rz(-1.707516) q[0];
sx q[0];
rz(1.7140634) q[0];
rz(-1.5643372) q[2];
sx q[2];
rz(-2.0142372) q[2];
sx q[2];
rz(2.9270767) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6306873) q[1];
sx q[1];
rz(-2.0954335) q[1];
sx q[1];
rz(1.8888297) q[1];
x q[2];
rz(2.0026358) q[3];
sx q[3];
rz(-1.6434998) q[3];
sx q[3];
rz(1.7464701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.67285377) q[2];
sx q[2];
rz(-1.8127952) q[2];
sx q[2];
rz(-1.2057745) q[2];
rz(-0.85637158) q[3];
sx q[3];
rz(-0.93256336) q[3];
sx q[3];
rz(0.87730733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67742753) q[0];
sx q[0];
rz(-2.7441315) q[0];
sx q[0];
rz(0.20832668) q[0];
rz(1.2267998) q[1];
sx q[1];
rz(-2.0715163) q[1];
sx q[1];
rz(-0.56455451) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3584261) q[0];
sx q[0];
rz(-1.4113398) q[0];
sx q[0];
rz(-0.043242976) q[0];
rz(0.73123572) q[2];
sx q[2];
rz(-1.952716) q[2];
sx q[2];
rz(-2.3637091) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.98726749) q[1];
sx q[1];
rz(-0.72817737) q[1];
sx q[1];
rz(2.8716692) q[1];
rz(-1.0309593) q[3];
sx q[3];
rz(-1.5542277) q[3];
sx q[3];
rz(-1.4318717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5047001) q[2];
sx q[2];
rz(-1.2496639) q[2];
sx q[2];
rz(-1.6553817) q[2];
rz(2.6460323) q[3];
sx q[3];
rz(-1.7335745) q[3];
sx q[3];
rz(-1.0068309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9946852) q[0];
sx q[0];
rz(-0.72145307) q[0];
sx q[0];
rz(1.7578693) q[0];
rz(1.2184628) q[1];
sx q[1];
rz(-1.0885025) q[1];
sx q[1];
rz(2.7168435) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.527194) q[0];
sx q[0];
rz(-1.5081811) q[0];
sx q[0];
rz(1.6024496) q[0];
rz(-pi) q[1];
rz(2.9016657) q[2];
sx q[2];
rz(-1.8169962) q[2];
sx q[2];
rz(2.4407516) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.37990824) q[1];
sx q[1];
rz(-1.9468005) q[1];
sx q[1];
rz(1.4892745) q[1];
x q[2];
rz(0.056428595) q[3];
sx q[3];
rz(-0.61739576) q[3];
sx q[3];
rz(-1.4790725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.72033) q[2];
sx q[2];
rz(-1.0676554) q[2];
sx q[2];
rz(-0.024451582) q[2];
rz(1.1484185) q[3];
sx q[3];
rz(-1.0992071) q[3];
sx q[3];
rz(-1.0816157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20632437) q[0];
sx q[0];
rz(-0.22980389) q[0];
sx q[0];
rz(2.9277053) q[0];
rz(2.3230486) q[1];
sx q[1];
rz(-1.3571502) q[1];
sx q[1];
rz(0.68665409) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.370433) q[0];
sx q[0];
rz(-0.70777445) q[0];
sx q[0];
rz(2.4999655) q[0];
x q[1];
rz(-1.0222203) q[2];
sx q[2];
rz(-1.3201664) q[2];
sx q[2];
rz(-2.938739) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5221646) q[1];
sx q[1];
rz(-0.81541895) q[1];
sx q[1];
rz(-1.4415037) q[1];
rz(0.88249256) q[3];
sx q[3];
rz(-2.5074041) q[3];
sx q[3];
rz(0.89804441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0246747) q[2];
sx q[2];
rz(-1.7322098) q[2];
sx q[2];
rz(2.4521258) q[2];
rz(-1.3245448) q[3];
sx q[3];
rz(-1.2169415) q[3];
sx q[3];
rz(-2.569258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4631184) q[0];
sx q[0];
rz(-3.0480223) q[0];
sx q[0];
rz(-2.1043188) q[0];
rz(-0.67266881) q[1];
sx q[1];
rz(-1.9692407) q[1];
sx q[1];
rz(-0.79576463) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.138863) q[0];
sx q[0];
rz(-0.61289591) q[0];
sx q[0];
rz(1.0010368) q[0];
rz(-0.99822171) q[2];
sx q[2];
rz(-1.3135944) q[2];
sx q[2];
rz(0.57359475) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.60476323) q[1];
sx q[1];
rz(-0.80789551) q[1];
sx q[1];
rz(1.1033586) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8759868) q[3];
sx q[3];
rz(-0.92545907) q[3];
sx q[3];
rz(-0.84922817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.98046389) q[2];
sx q[2];
rz(-2.0420065) q[2];
sx q[2];
rz(-2.3058092) q[2];
rz(0.25660822) q[3];
sx q[3];
rz(-2.0308688) q[3];
sx q[3];
rz(-2.6512644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68818727) q[0];
sx q[0];
rz(-1.1510993) q[0];
sx q[0];
rz(-1.9292462) q[0];
rz(-1.0351828) q[1];
sx q[1];
rz(-1.6500708) q[1];
sx q[1];
rz(-2.129668) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8674171) q[0];
sx q[0];
rz(-1.3998446) q[0];
sx q[0];
rz(0.11049185) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8746568) q[2];
sx q[2];
rz(-1.5603608) q[2];
sx q[2];
rz(-2.4316459) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3287836) q[1];
sx q[1];
rz(-1.9962096) q[1];
sx q[1];
rz(-2.0219457) q[1];
rz(-pi) q[2];
rz(2.7583241) q[3];
sx q[3];
rz(-1.82194) q[3];
sx q[3];
rz(1.0896089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.47784352) q[2];
sx q[2];
rz(-0.19890824) q[2];
sx q[2];
rz(-2.193023) q[2];
rz(-2.6712724) q[3];
sx q[3];
rz(-2.0430653) q[3];
sx q[3];
rz(-1.2547913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.54742852) q[0];
sx q[0];
rz(-1.0110039) q[0];
sx q[0];
rz(-1.8773361) q[0];
rz(2.3612379) q[1];
sx q[1];
rz(-0.64096132) q[1];
sx q[1];
rz(-2.9497214) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22671285) q[0];
sx q[0];
rz(-2.5396945) q[0];
sx q[0];
rz(-2.4826217) q[0];
rz(-pi) q[1];
rz(2.3290995) q[2];
sx q[2];
rz(-2.6408962) q[2];
sx q[2];
rz(-0.3270087) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7765155) q[1];
sx q[1];
rz(-2.0151225) q[1];
sx q[1];
rz(1.278484) q[1];
rz(-pi) q[2];
x q[2];
rz(0.34911134) q[3];
sx q[3];
rz(-0.99945949) q[3];
sx q[3];
rz(-3.0219363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4912305) q[2];
sx q[2];
rz(-1.0543084) q[2];
sx q[2];
rz(-0.80797705) q[2];
rz(-2.2982277) q[3];
sx q[3];
rz(-1.9873514) q[3];
sx q[3];
rz(1.5850867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.657044) q[0];
sx q[0];
rz(-3.1175685) q[0];
sx q[0];
rz(-2.3727544) q[0];
rz(-2.8688042) q[1];
sx q[1];
rz(-1.5130679) q[1];
sx q[1];
rz(1.3974894) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63015712) q[0];
sx q[0];
rz(-2.6574597) q[0];
sx q[0];
rz(0.70229437) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0447803) q[2];
sx q[2];
rz(-2.5616841) q[2];
sx q[2];
rz(1.2666262) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2721401) q[1];
sx q[1];
rz(-1.4542129) q[1];
sx q[1];
rz(-3.0805853) q[1];
x q[2];
rz(-2.0040183) q[3];
sx q[3];
rz(-0.94846359) q[3];
sx q[3];
rz(-2.4985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2046795) q[2];
sx q[2];
rz(-2.0759373) q[2];
sx q[2];
rz(1.6395456) q[2];
rz(-1.8223193) q[3];
sx q[3];
rz(-2.2917512) q[3];
sx q[3];
rz(1.5523065) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69916344) q[0];
sx q[0];
rz(-0.84446877) q[0];
sx q[0];
rz(-2.1635639) q[0];
rz(0.49081048) q[1];
sx q[1];
rz(-1.699479) q[1];
sx q[1];
rz(-1.6455654) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42616823) q[0];
sx q[0];
rz(-1.4472571) q[0];
sx q[0];
rz(-2.9036595) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3991576) q[2];
sx q[2];
rz(-1.1632068) q[2];
sx q[2];
rz(1.379058) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5514976) q[1];
sx q[1];
rz(-1.2480772) q[1];
sx q[1];
rz(-1.8264052) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0127527) q[3];
sx q[3];
rz(-2.6059756) q[3];
sx q[3];
rz(-2.420466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8533123) q[2];
sx q[2];
rz(-0.86709443) q[2];
sx q[2];
rz(1.8193998) q[2];
rz(-2.7754122) q[3];
sx q[3];
rz(-1.4270695) q[3];
sx q[3];
rz(-1.1314932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2243097) q[0];
sx q[0];
rz(-1.5320822) q[0];
sx q[0];
rz(-0.073534615) q[0];
rz(-1.8247983) q[1];
sx q[1];
rz(-1.8837181) q[1];
sx q[1];
rz(1.2988466) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6599346) q[0];
sx q[0];
rz(-1.7474239) q[0];
sx q[0];
rz(-2.8606979) q[0];
rz(-pi) q[1];
rz(2.0590563) q[2];
sx q[2];
rz(-1.4529145) q[2];
sx q[2];
rz(2.0669017) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6886814) q[1];
sx q[1];
rz(-1.1961403) q[1];
sx q[1];
rz(-0.25325798) q[1];
x q[2];
rz(0.73963005) q[3];
sx q[3];
rz(-1.4907903) q[3];
sx q[3];
rz(-0.39749872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9650044) q[2];
sx q[2];
rz(-2.1829288) q[2];
sx q[2];
rz(-0.57903543) q[2];
rz(-1.5163007) q[3];
sx q[3];
rz(-2.5987891) q[3];
sx q[3];
rz(1.496284) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8662921) q[0];
sx q[0];
rz(-1.0466812) q[0];
sx q[0];
rz(-0.064591797) q[0];
rz(-2.2304089) q[1];
sx q[1];
rz(-1.8235527) q[1];
sx q[1];
rz(1.8240737) q[1];
rz(0.0077132465) q[2];
sx q[2];
rz(-2.7174502) q[2];
sx q[2];
rz(2.9677718) q[2];
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
