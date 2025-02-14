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
rz(2.6449142) q[1];
sx q[1];
rz(4.7159046) q[1];
sx q[1];
rz(11.386303) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3008904) q[0];
sx q[0];
rz(-1.4288752) q[0];
sx q[0];
rz(0.13811708) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6981437) q[2];
sx q[2];
rz(-1.5766307) q[2];
sx q[2];
rz(1.3590517) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5109053) q[1];
sx q[1];
rz(-1.0461591) q[1];
sx q[1];
rz(-1.8888297) q[1];
rz(2.0026358) q[3];
sx q[3];
rz(-1.4980928) q[3];
sx q[3];
rz(1.3951226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.67285377) q[2];
sx q[2];
rz(-1.3287975) q[2];
sx q[2];
rz(1.2057745) q[2];
rz(-2.2852211) q[3];
sx q[3];
rz(-0.93256336) q[3];
sx q[3];
rz(-0.87730733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4641651) q[0];
sx q[0];
rz(-2.7441315) q[0];
sx q[0];
rz(-2.933266) q[0];
rz(-1.9147929) q[1];
sx q[1];
rz(-1.0700763) q[1];
sx q[1];
rz(-2.5770381) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0923742) q[0];
sx q[0];
rz(-2.9764247) q[0];
sx q[0];
rz(1.8333927) q[0];
rz(-pi) q[1];
rz(-0.73123572) q[2];
sx q[2];
rz(-1.1888767) q[2];
sx q[2];
rz(0.77788355) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7993421) q[1];
sx q[1];
rz(-0.87445757) q[1];
sx q[1];
rz(1.3373678) q[1];
x q[2];
rz(1.0309593) q[3];
sx q[3];
rz(-1.587365) q[3];
sx q[3];
rz(1.709721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5047001) q[2];
sx q[2];
rz(-1.2496639) q[2];
sx q[2];
rz(1.6553817) q[2];
rz(-0.49556035) q[3];
sx q[3];
rz(-1.7335745) q[3];
sx q[3];
rz(-1.0068309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9946852) q[0];
sx q[0];
rz(-2.4201396) q[0];
sx q[0];
rz(-1.3837234) q[0];
rz(-1.9231298) q[1];
sx q[1];
rz(-2.0530901) q[1];
sx q[1];
rz(0.4247492) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95441636) q[0];
sx q[0];
rz(-1.6023876) q[0];
sx q[0];
rz(-0.0626465) q[0];
rz(2.9016657) q[2];
sx q[2];
rz(-1.3245965) q[2];
sx q[2];
rz(0.70084106) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.37990824) q[1];
sx q[1];
rz(-1.9468005) q[1];
sx q[1];
rz(-1.6523182) q[1];
rz(0.61664403) q[3];
sx q[3];
rz(-1.5381406) q[3];
sx q[3];
rz(0.045696229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4212627) q[2];
sx q[2];
rz(-1.0676554) q[2];
sx q[2];
rz(3.1171411) q[2];
rz(1.1484185) q[3];
sx q[3];
rz(-1.0992071) q[3];
sx q[3];
rz(-1.0816157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
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
rz(-2.3230486) q[1];
sx q[1];
rz(-1.7844424) q[1];
sx q[1];
rz(0.68665409) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4249466) q[0];
sx q[0];
rz(-1.1711297) q[0];
sx q[0];
rz(-0.60092385) q[0];
rz(-pi) q[1];
rz(2.8501007) q[2];
sx q[2];
rz(-1.0412058) q[2];
sx q[2];
rz(-1.518371) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.80685341) q[1];
sx q[1];
rz(-0.76419965) q[1];
sx q[1];
rz(0.13607009) q[1];
rz(0.43709932) q[3];
sx q[3];
rz(-2.046114) q[3];
sx q[3];
rz(-3.0393485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0246747) q[2];
sx q[2];
rz(-1.7322098) q[2];
sx q[2];
rz(-0.68946687) q[2];
rz(1.8170478) q[3];
sx q[3];
rz(-1.9246512) q[3];
sx q[3];
rz(-0.57233468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6784742) q[0];
sx q[0];
rz(-0.093570396) q[0];
sx q[0];
rz(1.0372739) q[0];
rz(-2.4689238) q[1];
sx q[1];
rz(-1.172352) q[1];
sx q[1];
rz(2.345828) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.138863) q[0];
sx q[0];
rz(-0.61289591) q[0];
sx q[0];
rz(1.0010368) q[0];
rz(-2.8383083) q[2];
sx q[2];
rz(-1.0192843) q[2];
sx q[2];
rz(1.1597275) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2356253) q[1];
sx q[1];
rz(-0.86938953) q[1];
sx q[1];
rz(2.7011306) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4733801) q[3];
sx q[3];
rz(-1.8132121) q[3];
sx q[3];
rz(-2.6072864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.98046389) q[2];
sx q[2];
rz(-2.0420065) q[2];
sx q[2];
rz(2.3058092) q[2];
rz(0.25660822) q[3];
sx q[3];
rz(-2.0308688) q[3];
sx q[3];
rz(-2.6512644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68818727) q[0];
sx q[0];
rz(-1.9904933) q[0];
sx q[0];
rz(-1.2123464) q[0];
rz(2.1064099) q[1];
sx q[1];
rz(-1.6500708) q[1];
sx q[1];
rz(-2.129668) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4453011) q[0];
sx q[0];
rz(-2.9383349) q[0];
sx q[0];
rz(1.0023893) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2669358) q[2];
sx q[2];
rz(-1.5603608) q[2];
sx q[2];
rz(-2.4316459) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3287836) q[1];
sx q[1];
rz(-1.145383) q[1];
sx q[1];
rz(-1.119647) q[1];
rz(-2.7583241) q[3];
sx q[3];
rz(-1.3196527) q[3];
sx q[3];
rz(1.0896089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.47784352) q[2];
sx q[2];
rz(-0.19890824) q[2];
sx q[2];
rz(0.94856962) q[2];
rz(0.47032022) q[3];
sx q[3];
rz(-2.0430653) q[3];
sx q[3];
rz(-1.2547913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54742852) q[0];
sx q[0];
rz(-2.1305888) q[0];
sx q[0];
rz(1.8773361) q[0];
rz(0.7803548) q[1];
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
rz(-2.3656552) q[0];
sx q[0];
rz(-1.2167551) q[0];
sx q[0];
rz(0.49753041) q[0];
rz(-pi) q[1];
rz(0.35991798) q[2];
sx q[2];
rz(-1.9267757) q[2];
sx q[2];
rz(2.6448665) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3873008) q[1];
sx q[1];
rz(-0.52644345) q[1];
sx q[1];
rz(0.5443046) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0818008) q[3];
sx q[3];
rz(-0.65927514) q[3];
sx q[3];
rz(2.6687572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4912305) q[2];
sx q[2];
rz(-2.0872842) q[2];
sx q[2];
rz(2.3336156) q[2];
rz(-2.2982277) q[3];
sx q[3];
rz(-1.1542412) q[3];
sx q[3];
rz(-1.5850867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4845487) q[0];
sx q[0];
rz(-3.1175685) q[0];
sx q[0];
rz(2.3727544) q[0];
rz(-2.8688042) q[1];
sx q[1];
rz(-1.5130679) q[1];
sx q[1];
rz(-1.7441033) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63015712) q[0];
sx q[0];
rz(-2.6574597) q[0];
sx q[0];
rz(-0.70229437) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0553841) q[2];
sx q[2];
rz(-1.8495108) q[2];
sx q[2];
rz(-2.9936522) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3859017) q[1];
sx q[1];
rz(-0.13151691) q[1];
sx q[1];
rz(-2.0507858) q[1];
rz(-pi) q[2];
rz(-2.4726953) q[3];
sx q[3];
rz(-1.9188768) q[3];
sx q[3];
rz(-2.477248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2046795) q[2];
sx q[2];
rz(-1.0656554) q[2];
sx q[2];
rz(1.5020471) q[2];
rz(1.8223193) q[3];
sx q[3];
rz(-0.84984142) q[3];
sx q[3];
rz(-1.5892861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
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
rz(2.6507822) q[1];
sx q[1];
rz(-1.4421137) q[1];
sx q[1];
rz(1.4960272) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42616823) q[0];
sx q[0];
rz(-1.6943356) q[0];
sx q[0];
rz(-0.23793313) q[0];
rz(-pi) q[1];
rz(0.37668682) q[2];
sx q[2];
rz(-2.7012107) q[2];
sx q[2];
rz(0.96681606) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.89801187) q[1];
sx q[1];
rz(-1.3286546) q[1];
sx q[1];
rz(2.8088074) q[1];
rz(-pi) q[2];
rz(-1.6469025) q[3];
sx q[3];
rz(-1.0400912) q[3];
sx q[3];
rz(2.2709341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.28828037) q[2];
sx q[2];
rz(-0.86709443) q[2];
sx q[2];
rz(1.8193998) q[2];
rz(-0.36618048) q[3];
sx q[3];
rz(-1.4270695) q[3];
sx q[3];
rz(-2.0100994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
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
rz(-1.8837181) q[1];
sx q[1];
rz(1.2988466) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6599346) q[0];
sx q[0];
rz(-1.7474239) q[0];
sx q[0];
rz(2.8606979) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0082874) q[2];
sx q[2];
rz(-2.0553737) q[2];
sx q[2];
rz(-0.43371782) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9292752) q[1];
sx q[1];
rz(-1.8061418) q[1];
sx q[1];
rz(1.1849684) q[1];
rz(-pi) q[2];
rz(-1.6789087) q[3];
sx q[3];
rz(-0.834081) q[3];
sx q[3];
rz(-2.0410868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-1.6453086) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
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
rz(-0.0077132465) q[2];
sx q[2];
rz(-0.42414244) q[2];
sx q[2];
rz(-0.17382081) q[2];
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
