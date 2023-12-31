OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.6289829) q[0];
sx q[0];
rz(-2.1142168) q[0];
sx q[0];
rz(2.7789814) q[0];
rz(0.91712046) q[1];
sx q[1];
rz(2.6511104) q[1];
sx q[1];
rz(9.766415) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61603123) q[0];
sx q[0];
rz(-1.4481059) q[0];
sx q[0];
rz(0.37567715) q[0];
rz(-pi) q[1];
rz(-2.6356959) q[2];
sx q[2];
rz(-0.6542754) q[2];
sx q[2];
rz(1.8088532) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8475854) q[1];
sx q[1];
rz(-0.60677401) q[1];
sx q[1];
rz(-2.1104382) q[1];
rz(-pi) q[2];
rz(2.0658675) q[3];
sx q[3];
rz(-0.44701156) q[3];
sx q[3];
rz(1.8398374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5477156) q[2];
sx q[2];
rz(-1.9888069) q[2];
sx q[2];
rz(-2.1851052) q[2];
rz(0.18125136) q[3];
sx q[3];
rz(-0.58877188) q[3];
sx q[3];
rz(1.6199002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0797121) q[0];
sx q[0];
rz(-0.065608874) q[0];
sx q[0];
rz(-1.99615) q[0];
rz(2.0727797) q[1];
sx q[1];
rz(-0.83199465) q[1];
sx q[1];
rz(2.4157445) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94978588) q[0];
sx q[0];
rz(-1.0254142) q[0];
sx q[0];
rz(-1.1508862) q[0];
rz(0.59402324) q[2];
sx q[2];
rz(-1.960264) q[2];
sx q[2];
rz(2.4246852) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.010163807) q[1];
sx q[1];
rz(-1.5958852) q[1];
sx q[1];
rz(-1.2618834) q[1];
rz(-pi) q[2];
rz(1.782136) q[3];
sx q[3];
rz(-2.1826715) q[3];
sx q[3];
rz(2.8733727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6590865) q[2];
sx q[2];
rz(-1.3799474) q[2];
sx q[2];
rz(2.144311) q[2];
rz(-1.7287792) q[3];
sx q[3];
rz(-2.0480859) q[3];
sx q[3];
rz(2.2028082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41855758) q[0];
sx q[0];
rz(-0.96590531) q[0];
sx q[0];
rz(-2.0130656) q[0];
rz(-1.8380802) q[1];
sx q[1];
rz(-0.729527) q[1];
sx q[1];
rz(2.2361141) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14074621) q[0];
sx q[0];
rz(-1.8119537) q[0];
sx q[0];
rz(-1.3748037) q[0];
rz(-pi) q[1];
rz(2.0037903) q[2];
sx q[2];
rz(-2.275122) q[2];
sx q[2];
rz(-1.1906884) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.803373) q[1];
sx q[1];
rz(-1.4969016) q[1];
sx q[1];
rz(-2.282663) q[1];
rz(0.043025322) q[3];
sx q[3];
rz(-0.85775162) q[3];
sx q[3];
rz(-2.713664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.7604312) q[2];
sx q[2];
rz(-2.9734549) q[2];
sx q[2];
rz(-2.1543489) q[2];
rz(3.0958214) q[3];
sx q[3];
rz(-1.2922623) q[3];
sx q[3];
rz(0.18019095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-3.0320597) q[0];
sx q[0];
rz(-0.68234545) q[0];
sx q[0];
rz(-0.12606829) q[0];
rz(-2.9571422) q[1];
sx q[1];
rz(-1.0686864) q[1];
sx q[1];
rz(-1.1674081) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2950738) q[0];
sx q[0];
rz(-0.82058883) q[0];
sx q[0];
rz(1.2942765) q[0];
x q[1];
rz(-1.3912958) q[2];
sx q[2];
rz(-1.0738392) q[2];
sx q[2];
rz(-0.77979445) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.57707225) q[1];
sx q[1];
rz(-2.2153691) q[1];
sx q[1];
rz(1.2381899) q[1];
rz(-pi) q[2];
rz(0.91058369) q[3];
sx q[3];
rz(-1.0331717) q[3];
sx q[3];
rz(-2.0349353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.8644774) q[2];
sx q[2];
rz(-2.7272759) q[2];
sx q[2];
rz(0.63151044) q[2];
rz(-0.19550066) q[3];
sx q[3];
rz(-1.2771527) q[3];
sx q[3];
rz(-0.41803944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9773848) q[0];
sx q[0];
rz(-3.003484) q[0];
sx q[0];
rz(0.50267977) q[0];
rz(1.1133105) q[1];
sx q[1];
rz(-1.0031507) q[1];
sx q[1];
rz(-2.6745093) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0665849) q[0];
sx q[0];
rz(-1.9420615) q[0];
sx q[0];
rz(-2.9270372) q[0];
rz(-2.5555243) q[2];
sx q[2];
rz(-1.35891) q[2];
sx q[2];
rz(2.6168602) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0847229) q[1];
sx q[1];
rz(-1.9922171) q[1];
sx q[1];
rz(-2.9028068) q[1];
x q[2];
rz(2.983106) q[3];
sx q[3];
rz(-0.89156686) q[3];
sx q[3];
rz(2.2331626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6749394) q[2];
sx q[2];
rz(-1.249524) q[2];
sx q[2];
rz(-2.6494027) q[2];
rz(-2.8594033) q[3];
sx q[3];
rz(-1.7084833) q[3];
sx q[3];
rz(-1.5757489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2714587) q[0];
sx q[0];
rz(-0.54015714) q[0];
sx q[0];
rz(-0.085993275) q[0];
rz(-1.7135235) q[1];
sx q[1];
rz(-2.1336887) q[1];
sx q[1];
rz(-1.496398) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8498358) q[0];
sx q[0];
rz(-1.0514326) q[0];
sx q[0];
rz(-2.915328) q[0];
rz(2.1093844) q[2];
sx q[2];
rz(-1.7787873) q[2];
sx q[2];
rz(-1.5185192) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8346356) q[1];
sx q[1];
rz(-0.35528696) q[1];
sx q[1];
rz(-1.6389695) q[1];
rz(-pi) q[2];
rz(0.014724894) q[3];
sx q[3];
rz(-0.71299362) q[3];
sx q[3];
rz(-1.7409489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4795586) q[2];
sx q[2];
rz(-2.5223314) q[2];
sx q[2];
rz(2.7453444) q[2];
rz(0.59404343) q[3];
sx q[3];
rz(-1.0915353) q[3];
sx q[3];
rz(1.6408287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8783766) q[0];
sx q[0];
rz(-0.57290572) q[0];
sx q[0];
rz(2.0492045) q[0];
rz(-1.3573525) q[1];
sx q[1];
rz(-1.7204294) q[1];
sx q[1];
rz(-0.011627442) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25225885) q[0];
sx q[0];
rz(-0.47124915) q[0];
sx q[0];
rz(2.528119) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1161357) q[2];
sx q[2];
rz(-1.5989132) q[2];
sx q[2];
rz(2.1643929) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7203622) q[1];
sx q[1];
rz(-2.1167397) q[1];
sx q[1];
rz(-3.1139042) q[1];
x q[2];
rz(2.4680448) q[3];
sx q[3];
rz(-1.827497) q[3];
sx q[3];
rz(1.7208769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6860883) q[2];
sx q[2];
rz(-0.65767613) q[2];
sx q[2];
rz(-2.8536076) q[2];
rz(1.1503495) q[3];
sx q[3];
rz(-1.3214313) q[3];
sx q[3];
rz(-2.5434125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.8312254) q[0];
sx q[0];
rz(-2.6103525) q[0];
sx q[0];
rz(-1.9294552) q[0];
rz(-0.37777004) q[1];
sx q[1];
rz(-1.9394082) q[1];
sx q[1];
rz(-1.4935965) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50266788) q[0];
sx q[0];
rz(-1.272861) q[0];
sx q[0];
rz(1.7940679) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1475032) q[2];
sx q[2];
rz(-0.14483843) q[2];
sx q[2];
rz(-2.7763979) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4920602) q[1];
sx q[1];
rz(-0.82693716) q[1];
sx q[1];
rz(2.5745113) q[1];
x q[2];
rz(-1.0687374) q[3];
sx q[3];
rz(-1.6747253) q[3];
sx q[3];
rz(0.41307377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3528072) q[2];
sx q[2];
rz(-1.1151423) q[2];
sx q[2];
rz(-1.2517694) q[2];
rz(-0.88636032) q[3];
sx q[3];
rz(-2.3754407) q[3];
sx q[3];
rz(3.0322976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4733474) q[0];
sx q[0];
rz(-2.1003523) q[0];
sx q[0];
rz(-2.9633203) q[0];
rz(-3.0687304) q[1];
sx q[1];
rz(-2.5477414) q[1];
sx q[1];
rz(1.1791621) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7992135) q[0];
sx q[0];
rz(-1.8633435) q[0];
sx q[0];
rz(-0.53036687) q[0];
rz(-0.945325) q[2];
sx q[2];
rz(-1.9918459) q[2];
sx q[2];
rz(2.9025214) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0852016) q[1];
sx q[1];
rz(-1.3136275) q[1];
sx q[1];
rz(-2.609054) q[1];
x q[2];
rz(-2.7618802) q[3];
sx q[3];
rz(-2.0413997) q[3];
sx q[3];
rz(3.0203366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.69958413) q[2];
sx q[2];
rz(-0.99094355) q[2];
sx q[2];
rz(-2.1098095) q[2];
rz(-1.7992841) q[3];
sx q[3];
rz(-1.7248025) q[3];
sx q[3];
rz(-0.22527307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28829065) q[0];
sx q[0];
rz(-1.0422215) q[0];
sx q[0];
rz(3.0798966) q[0];
rz(-1.8348947) q[1];
sx q[1];
rz(-1.6625762) q[1];
sx q[1];
rz(2.0526989) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93931224) q[0];
sx q[0];
rz(-1.5081524) q[0];
sx q[0];
rz(-1.63562) q[0];
rz(2.225869) q[2];
sx q[2];
rz(-1.7172126) q[2];
sx q[2];
rz(-0.54619782) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.006146487) q[1];
sx q[1];
rz(-2.1278312) q[1];
sx q[1];
rz(1.0165434) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9519134) q[3];
sx q[3];
rz(-1.9072201) q[3];
sx q[3];
rz(-1.3996901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6178199) q[2];
sx q[2];
rz(-2.4536295) q[2];
sx q[2];
rz(-0.069996746) q[2];
rz(-2.2648515) q[3];
sx q[3];
rz(-1.8165959) q[3];
sx q[3];
rz(-0.35818067) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4927647) q[0];
sx q[0];
rz(-1.6236826) q[0];
sx q[0];
rz(0.59326011) q[0];
rz(1.6246673) q[1];
sx q[1];
rz(-1.0472052) q[1];
sx q[1];
rz(-0.44890961) q[1];
rz(-0.4060612) q[2];
sx q[2];
rz(-0.95352298) q[2];
sx q[2];
rz(2.7609115) q[2];
rz(-3.1003351) q[3];
sx q[3];
rz(-1.2357124) q[3];
sx q[3];
rz(0.99832051) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
