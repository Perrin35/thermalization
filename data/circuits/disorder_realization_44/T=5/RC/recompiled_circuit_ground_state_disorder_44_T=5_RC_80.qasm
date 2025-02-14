OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7527591) q[0];
sx q[0];
rz(1.692481) q[0];
sx q[0];
rz(11.115885) q[0];
rz(0.9736355) q[1];
sx q[1];
rz(-1.7042301) q[1];
sx q[1];
rz(-0.91926423) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96199233) q[0];
sx q[0];
rz(-1.5878979) q[0];
sx q[0];
rz(0.03058612) q[0];
x q[1];
rz(-3.0212901) q[2];
sx q[2];
rz(-1.4805111) q[2];
sx q[2];
rz(-1.9419958) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8074933) q[1];
sx q[1];
rz(-1.4791282) q[1];
sx q[1];
rz(0.86218545) q[1];
rz(0.84897016) q[3];
sx q[3];
rz(-2.5458286) q[3];
sx q[3];
rz(1.1536382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3322525) q[2];
sx q[2];
rz(-1.6925749) q[2];
sx q[2];
rz(1.4130886) q[2];
rz(2.938802) q[3];
sx q[3];
rz(-1.7594124) q[3];
sx q[3];
rz(3.0387759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8973812) q[0];
sx q[0];
rz(-1.0845217) q[0];
sx q[0];
rz(0.71075034) q[0];
rz(0.76849014) q[1];
sx q[1];
rz(-2.0748506) q[1];
sx q[1];
rz(-1.01952) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33776721) q[0];
sx q[0];
rz(-0.44887421) q[0];
sx q[0];
rz(3.0503737) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1395733) q[2];
sx q[2];
rz(-1.6399472) q[2];
sx q[2];
rz(2.5716242) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4162035) q[1];
sx q[1];
rz(-0.37230154) q[1];
sx q[1];
rz(-1.4854234) q[1];
rz(-1.4367661) q[3];
sx q[3];
rz(-2.5833327) q[3];
sx q[3];
rz(2.9001065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3780313) q[2];
sx q[2];
rz(-1.8969994) q[2];
sx q[2];
rz(1.2223318) q[2];
rz(1.2126806) q[3];
sx q[3];
rz(-2.1247532) q[3];
sx q[3];
rz(0.71162629) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.057230495) q[0];
sx q[0];
rz(-0.98452345) q[0];
sx q[0];
rz(-1.9236176) q[0];
rz(-0.51586622) q[1];
sx q[1];
rz(-0.54793826) q[1];
sx q[1];
rz(0.9224433) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9003214) q[0];
sx q[0];
rz(-1.5301955) q[0];
sx q[0];
rz(-3.0837644) q[0];
rz(-2.0751245) q[2];
sx q[2];
rz(-2.2309003) q[2];
sx q[2];
rz(-2.4911027) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.565971) q[1];
sx q[1];
rz(-2.1346483) q[1];
sx q[1];
rz(-0.82852461) q[1];
rz(-pi) q[2];
rz(-0.46307989) q[3];
sx q[3];
rz(-2.7768917) q[3];
sx q[3];
rz(2.0228342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.018365232) q[2];
sx q[2];
rz(-1.0845228) q[2];
sx q[2];
rz(1.3398735) q[2];
rz(-0.54723048) q[3];
sx q[3];
rz(-2.7180143) q[3];
sx q[3];
rz(2.4303998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.055534) q[0];
sx q[0];
rz(-0.14500293) q[0];
sx q[0];
rz(2.9823629) q[0];
rz(-0.010146443) q[1];
sx q[1];
rz(-1.0228913) q[1];
sx q[1];
rz(-2.8841282) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87436324) q[0];
sx q[0];
rz(-0.82634514) q[0];
sx q[0];
rz(2.6494725) q[0];
x q[1];
rz(-1.5045549) q[2];
sx q[2];
rz(-2.3640354) q[2];
sx q[2];
rz(2.24868) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3336948) q[1];
sx q[1];
rz(-0.94606864) q[1];
sx q[1];
rz(1.4480535) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2063755) q[3];
sx q[3];
rz(-2.2699589) q[3];
sx q[3];
rz(2.2548667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3782392) q[2];
sx q[2];
rz(-1.2946318) q[2];
sx q[2];
rz(2.6687458) q[2];
rz(2.4371448) q[3];
sx q[3];
rz(-1.364578) q[3];
sx q[3];
rz(0.64594597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5064297) q[0];
sx q[0];
rz(-1.9671054) q[0];
sx q[0];
rz(-0.065486431) q[0];
rz(0.40924117) q[1];
sx q[1];
rz(-2.0114653) q[1];
sx q[1];
rz(1.4261036) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5201841) q[0];
sx q[0];
rz(-1.4778293) q[0];
sx q[0];
rz(0.21958242) q[0];
rz(-0.21435301) q[2];
sx q[2];
rz(-0.39696908) q[2];
sx q[2];
rz(-3.1052542) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8219442) q[1];
sx q[1];
rz(-1.8488171) q[1];
sx q[1];
rz(1.3389498) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9687443) q[3];
sx q[3];
rz(-1.2842872) q[3];
sx q[3];
rz(2.2423173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9466729) q[2];
sx q[2];
rz(-2.0544923) q[2];
sx q[2];
rz(1.8348414) q[2];
rz(1.170018) q[3];
sx q[3];
rz(-2.7368059) q[3];
sx q[3];
rz(-0.10725966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65866798) q[0];
sx q[0];
rz(-1.1558477) q[0];
sx q[0];
rz(0.40147993) q[0];
rz(-0.67277706) q[1];
sx q[1];
rz(-0.79877001) q[1];
sx q[1];
rz(-0.85404095) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9977048) q[0];
sx q[0];
rz(-1.524462) q[0];
sx q[0];
rz(-1.3530988) q[0];
rz(-0.18547345) q[2];
sx q[2];
rz(-1.7181686) q[2];
sx q[2];
rz(-0.81306785) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8484162) q[1];
sx q[1];
rz(-2.5418315) q[1];
sx q[1];
rz(2.9648215) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.61088125) q[3];
sx q[3];
rz(-1.9457091) q[3];
sx q[3];
rz(-1.5883816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.65471571) q[2];
sx q[2];
rz(-0.87149039) q[2];
sx q[2];
rz(-1.1775449) q[2];
rz(0.55142895) q[3];
sx q[3];
rz(-1.5588372) q[3];
sx q[3];
rz(2.3059755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7198782) q[0];
sx q[0];
rz(-2.8604909) q[0];
sx q[0];
rz(-1.1623435) q[0];
rz(-1.1516736) q[1];
sx q[1];
rz(-1.3538227) q[1];
sx q[1];
rz(-0.91845671) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7699444) q[0];
sx q[0];
rz(-1.3027096) q[0];
sx q[0];
rz(1.2275342) q[0];
rz(0.26084857) q[2];
sx q[2];
rz(-0.89226228) q[2];
sx q[2];
rz(1.7186173) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9388401) q[1];
sx q[1];
rz(-1.5781286) q[1];
sx q[1];
rz(1.5775024) q[1];
rz(1.1226467) q[3];
sx q[3];
rz(-0.3780685) q[3];
sx q[3];
rz(-2.700875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.1071757) q[2];
sx q[2];
rz(-1.7566661) q[2];
sx q[2];
rz(-2.6386063) q[2];
rz(-0.77477396) q[3];
sx q[3];
rz(-0.46190244) q[3];
sx q[3];
rz(-1.3431965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69304943) q[0];
sx q[0];
rz(-0.22668426) q[0];
sx q[0];
rz(-1.5013129) q[0];
rz(2.8003108) q[1];
sx q[1];
rz(-1.4309859) q[1];
sx q[1];
rz(-1.1192082) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.354108) q[0];
sx q[0];
rz(-1.2006) q[0];
sx q[0];
rz(0.60109971) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0518752) q[2];
sx q[2];
rz(-1.8029034) q[2];
sx q[2];
rz(-0.26929917) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1895442) q[1];
sx q[1];
rz(-2.2276926) q[1];
sx q[1];
rz(-1.4178965) q[1];
x q[2];
rz(1.7057034) q[3];
sx q[3];
rz(-1.9622318) q[3];
sx q[3];
rz(1.6329444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.98562733) q[2];
sx q[2];
rz(-0.95981821) q[2];
sx q[2];
rz(-0.36925527) q[2];
rz(-2.8333832) q[3];
sx q[3];
rz(-0.9674558) q[3];
sx q[3];
rz(1.5306028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8592598) q[0];
sx q[0];
rz(-1.8349324) q[0];
sx q[0];
rz(2.7985213) q[0];
rz(-1.169091) q[1];
sx q[1];
rz(-1.7770551) q[1];
sx q[1];
rz(1.7128568) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37821445) q[0];
sx q[0];
rz(-1.175048) q[0];
sx q[0];
rz(-1.8431435) q[0];
rz(-pi) q[1];
rz(-0.69099094) q[2];
sx q[2];
rz(-2.1804902) q[2];
sx q[2];
rz(2.8728849) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1395806) q[1];
sx q[1];
rz(-1.6320458) q[1];
sx q[1];
rz(-1.9467926) q[1];
rz(-pi) q[2];
x q[2];
rz(1.021528) q[3];
sx q[3];
rz(-2.3695951) q[3];
sx q[3];
rz(2.271351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2009361) q[2];
sx q[2];
rz(-2.9328465) q[2];
sx q[2];
rz(-1.7500056) q[2];
rz(-2.7348147) q[3];
sx q[3];
rz(-1.6986366) q[3];
sx q[3];
rz(0.87882915) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0402891) q[0];
sx q[0];
rz(-2.5387796) q[0];
sx q[0];
rz(1.783675) q[0];
rz(1.9039924) q[1];
sx q[1];
rz(-2.1289181) q[1];
sx q[1];
rz(-0.79992574) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7800723) q[0];
sx q[0];
rz(-1.8006386) q[0];
sx q[0];
rz(0.80432463) q[0];
rz(-0.47100701) q[2];
sx q[2];
rz(-0.60294294) q[2];
sx q[2];
rz(2.431536) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8854196) q[1];
sx q[1];
rz(-0.90682632) q[1];
sx q[1];
rz(1.3732984) q[1];
x q[2];
rz(-0.25000817) q[3];
sx q[3];
rz(-2.0431314) q[3];
sx q[3];
rz(1.8369305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.37821975) q[2];
sx q[2];
rz(-1.7174481) q[2];
sx q[2];
rz(-3.0685032) q[2];
rz(-2.8857005) q[3];
sx q[3];
rz(-2.3292694) q[3];
sx q[3];
rz(1.6528486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.273461) q[0];
sx q[0];
rz(-1.68597) q[0];
sx q[0];
rz(-1.2706533) q[0];
rz(-1.801626) q[1];
sx q[1];
rz(-1.3506964) q[1];
sx q[1];
rz(0.66257308) q[1];
rz(0.70941464) q[2];
sx q[2];
rz(-1.5723036) q[2];
sx q[2];
rz(0.2607762) q[2];
rz(-0.64107278) q[3];
sx q[3];
rz(-1.0985634) q[3];
sx q[3];
rz(-2.6177277) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
