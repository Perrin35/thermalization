OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4797526) q[0];
sx q[0];
rz(-2.2979484) q[0];
sx q[0];
rz(2.9736829) q[0];
rz(-1.9703938) q[1];
sx q[1];
rz(-0.29532239) q[1];
sx q[1];
rz(3.0854316) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9553298) q[0];
sx q[0];
rz(-1.3583399) q[0];
sx q[0];
rz(2.0583378) q[0];
rz(-pi) q[1];
rz(0.92476966) q[2];
sx q[2];
rz(-1.3999108) q[2];
sx q[2];
rz(0.31121635) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9819298) q[1];
sx q[1];
rz(-0.61239457) q[1];
sx q[1];
rz(-1.0931404) q[1];
x q[2];
rz(-0.25986259) q[3];
sx q[3];
rz(-1.6136323) q[3];
sx q[3];
rz(-2.6916137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.37796676) q[2];
sx q[2];
rz(-2.8597735) q[2];
sx q[2];
rz(-2.7089233) q[2];
rz(-1.1928605) q[3];
sx q[3];
rz(-1.9038707) q[3];
sx q[3];
rz(0.38309923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52779657) q[0];
sx q[0];
rz(-2.6531117) q[0];
sx q[0];
rz(1.8288076) q[0];
rz(-0.20547543) q[1];
sx q[1];
rz(-2.165129) q[1];
sx q[1];
rz(-1.1516494) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0002999) q[0];
sx q[0];
rz(-0.89501689) q[0];
sx q[0];
rz(-2.7594901) q[0];
rz(-pi) q[1];
rz(-1.0863016) q[2];
sx q[2];
rz(-2.0995579) q[2];
sx q[2];
rz(2.8690086) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.68576751) q[1];
sx q[1];
rz(-1.1357422) q[1];
sx q[1];
rz(-2.0352092) q[1];
x q[2];
rz(0.63668164) q[3];
sx q[3];
rz(-0.59906206) q[3];
sx q[3];
rz(0.77081313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.1318704) q[2];
sx q[2];
rz(-1.4902318) q[2];
sx q[2];
rz(-0.22182626) q[2];
rz(-0.37718537) q[3];
sx q[3];
rz(-2.714034) q[3];
sx q[3];
rz(0.77243531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31056988) q[0];
sx q[0];
rz(-0.092386827) q[0];
sx q[0];
rz(-0.036852766) q[0];
rz(-2.316078) q[1];
sx q[1];
rz(-1.8258391) q[1];
sx q[1];
rz(0.056578606) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3765592) q[0];
sx q[0];
rz(-0.74900904) q[0];
sx q[0];
rz(0.88699938) q[0];
rz(0.47358863) q[2];
sx q[2];
rz(-2.8550672) q[2];
sx q[2];
rz(2.5055656) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.50443447) q[1];
sx q[1];
rz(-1.27099) q[1];
sx q[1];
rz(-1.7791041) q[1];
rz(-pi) q[2];
rz(-0.49719663) q[3];
sx q[3];
rz(-0.70780863) q[3];
sx q[3];
rz(-1.7266112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.27292192) q[2];
sx q[2];
rz(-1.6438831) q[2];
sx q[2];
rz(0.92612129) q[2];
rz(0.55666322) q[3];
sx q[3];
rz(-2.848048) q[3];
sx q[3];
rz(-2.0986957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91519231) q[0];
sx q[0];
rz(-0.72015786) q[0];
sx q[0];
rz(-0.74209374) q[0];
rz(2.0023951) q[1];
sx q[1];
rz(-0.4793872) q[1];
sx q[1];
rz(-0.46359584) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47953654) q[0];
sx q[0];
rz(-2.2502796) q[0];
sx q[0];
rz(2.7583073) q[0];
rz(-pi) q[1];
rz(0.33275231) q[2];
sx q[2];
rz(-1.5126265) q[2];
sx q[2];
rz(1.0156877) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.78471781) q[1];
sx q[1];
rz(-2.3349635) q[1];
sx q[1];
rz(-0.23826092) q[1];
rz(-pi) q[2];
x q[2];
rz(0.87426825) q[3];
sx q[3];
rz(-0.6797176) q[3];
sx q[3];
rz(0.10248871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3670369) q[2];
sx q[2];
rz(-2.98428) q[2];
sx q[2];
rz(-0.049499361) q[2];
rz(0.1285304) q[3];
sx q[3];
rz(-1.5481719) q[3];
sx q[3];
rz(3.0226504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0054935) q[0];
sx q[0];
rz(-0.43731421) q[0];
sx q[0];
rz(2.8438925) q[0];
rz(-2.659335) q[1];
sx q[1];
rz(-0.75459701) q[1];
sx q[1];
rz(0.94435) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1221065) q[0];
sx q[0];
rz(-1.5257611) q[0];
sx q[0];
rz(-1.7800063) q[0];
rz(-pi) q[1];
rz(2.3659533) q[2];
sx q[2];
rz(-1.8619985) q[2];
sx q[2];
rz(-1.2736543) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.797232) q[1];
sx q[1];
rz(-1.5096944) q[1];
sx q[1];
rz(-0.0021211591) q[1];
rz(0.51575757) q[3];
sx q[3];
rz(-0.48832794) q[3];
sx q[3];
rz(-2.8749136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8828316) q[2];
sx q[2];
rz(-1.0927039) q[2];
sx q[2];
rz(-3.0333701) q[2];
rz(-0.0023068874) q[3];
sx q[3];
rz(-1.5294411) q[3];
sx q[3];
rz(0.32430696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7047983) q[0];
sx q[0];
rz(-2.695485) q[0];
sx q[0];
rz(2.5571402) q[0];
rz(0.8862409) q[1];
sx q[1];
rz(-2.5247572) q[1];
sx q[1];
rz(0.054919682) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7391101) q[0];
sx q[0];
rz(-1.6533028) q[0];
sx q[0];
rz(-1.301618) q[0];
rz(-pi) q[1];
rz(3.1228742) q[2];
sx q[2];
rz(-1.8939549) q[2];
sx q[2];
rz(-1.2013555) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.37327787) q[1];
sx q[1];
rz(-1.0265961) q[1];
sx q[1];
rz(-1.1276223) q[1];
rz(0.98723282) q[3];
sx q[3];
rz(-1.9658486) q[3];
sx q[3];
rz(2.6810255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.75446689) q[2];
sx q[2];
rz(-3.0209164) q[2];
sx q[2];
rz(-1.0167271) q[2];
rz(2.5975442) q[3];
sx q[3];
rz(-2.7816911) q[3];
sx q[3];
rz(1.8090766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.465437) q[0];
sx q[0];
rz(-0.98452079) q[0];
sx q[0];
rz(-2.8570535) q[0];
rz(-0.94447213) q[1];
sx q[1];
rz(-1.1962793) q[1];
sx q[1];
rz(-2.231266) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0383354) q[0];
sx q[0];
rz(-0.084780134) q[0];
sx q[0];
rz(2.2708562) q[0];
rz(-pi) q[1];
rz(1.1548642) q[2];
sx q[2];
rz(-1.8850733) q[2];
sx q[2];
rz(-1.2736125) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7200304) q[1];
sx q[1];
rz(-1.7517462) q[1];
sx q[1];
rz(2.6471495) q[1];
rz(-0.99366412) q[3];
sx q[3];
rz(-2.2054407) q[3];
sx q[3];
rz(1.5821379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7515144) q[2];
sx q[2];
rz(-0.063515924) q[2];
sx q[2];
rz(-2.2195623) q[2];
rz(0.56728029) q[3];
sx q[3];
rz(-1.7276948) q[3];
sx q[3];
rz(-2.1218307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89408016) q[0];
sx q[0];
rz(-0.67665726) q[0];
sx q[0];
rz(-3.0122053) q[0];
rz(0.63240504) q[1];
sx q[1];
rz(-2.1148732) q[1];
sx q[1];
rz(0.30050373) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7104615) q[0];
sx q[0];
rz(-2.3346402) q[0];
sx q[0];
rz(-2.0956844) q[0];
x q[1];
rz(0.6408765) q[2];
sx q[2];
rz(-2.2705728) q[2];
sx q[2];
rz(-1.7388294) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.44725542) q[1];
sx q[1];
rz(-1.2303423) q[1];
sx q[1];
rz(-1.3854331) q[1];
rz(-1.9037876) q[3];
sx q[3];
rz(-2.3559642) q[3];
sx q[3];
rz(-1.9612519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.58632103) q[2];
sx q[2];
rz(-0.99493146) q[2];
sx q[2];
rz(-2.3596181) q[2];
rz(0.55109763) q[3];
sx q[3];
rz(-1.3827773) q[3];
sx q[3];
rz(-2.9836695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5732116) q[0];
sx q[0];
rz(-1.1567572) q[0];
sx q[0];
rz(3.0138299) q[0];
rz(2.5993775) q[1];
sx q[1];
rz(-0.95710373) q[1];
sx q[1];
rz(0.75884563) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.471506) q[0];
sx q[0];
rz(-2.717088) q[0];
sx q[0];
rz(-1.612081) q[0];
rz(-pi) q[1];
rz(2.5326469) q[2];
sx q[2];
rz(-1.4119838) q[2];
sx q[2];
rz(1.9132523) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.54770494) q[1];
sx q[1];
rz(-2.2624359) q[1];
sx q[1];
rz(1.7734852) q[1];
rz(-pi) q[2];
rz(1.5893448) q[3];
sx q[3];
rz(-0.66569257) q[3];
sx q[3];
rz(1.9513643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1252497) q[2];
sx q[2];
rz(-1.3602076) q[2];
sx q[2];
rz(0.49003595) q[2];
rz(-1.4222493) q[3];
sx q[3];
rz(-1.2079206) q[3];
sx q[3];
rz(1.0629883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35995099) q[0];
sx q[0];
rz(-2.521823) q[0];
sx q[0];
rz(3.066257) q[0];
rz(-0.8967337) q[1];
sx q[1];
rz(-1.1743841) q[1];
sx q[1];
rz(-2.5316701) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41994914) q[0];
sx q[0];
rz(-0.83166612) q[0];
sx q[0];
rz(1.9589817) q[0];
rz(-pi) q[1];
rz(1.8462734) q[2];
sx q[2];
rz(-1.2071949) q[2];
sx q[2];
rz(-1.3678577) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3366821) q[1];
sx q[1];
rz(-1.8098127) q[1];
sx q[1];
rz(0.75093021) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6033953) q[3];
sx q[3];
rz(-2.3231069) q[3];
sx q[3];
rz(-1.8173816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.23218368) q[2];
sx q[2];
rz(-2.3164618) q[2];
sx q[2];
rz(0.71371901) q[2];
rz(0.37832007) q[3];
sx q[3];
rz(-0.49351966) q[3];
sx q[3];
rz(-0.87987125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.338035) q[0];
sx q[0];
rz(-1.150158) q[0];
sx q[0];
rz(-1.5858571) q[0];
rz(-0.65080416) q[1];
sx q[1];
rz(-1.4918068) q[1];
sx q[1];
rz(3.0202958) q[1];
rz(0.030608721) q[2];
sx q[2];
rz(-1.7666713) q[2];
sx q[2];
rz(-0.91790744) q[2];
rz(0.25272947) q[3];
sx q[3];
rz(-1.1128294) q[3];
sx q[3];
rz(2.2313234) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
