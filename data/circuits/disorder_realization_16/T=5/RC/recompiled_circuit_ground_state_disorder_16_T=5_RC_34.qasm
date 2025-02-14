OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.0082173) q[0];
sx q[0];
rz(-1.2674588) q[0];
sx q[0];
rz(-0.01292364) q[0];
rz(0.68459964) q[1];
sx q[1];
rz(3.9404865) q[1];
sx q[1];
rz(10.482492) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36127451) q[0];
sx q[0];
rz(-2.4343505) q[0];
sx q[0];
rz(0.35298423) q[0];
rz(-pi) q[1];
x q[1];
rz(0.77785141) q[2];
sx q[2];
rz(-1.1925282) q[2];
sx q[2];
rz(-3.0562378) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2014241) q[1];
sx q[1];
rz(-0.88646171) q[1];
sx q[1];
rz(-1.0971054) q[1];
x q[2];
rz(-3.0989981) q[3];
sx q[3];
rz(-1.9622246) q[3];
sx q[3];
rz(3.1393876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.58618033) q[2];
sx q[2];
rz(-1.5427898) q[2];
sx q[2];
rz(3.0482698) q[2];
rz(-3.1210476) q[3];
sx q[3];
rz(-2.8829657) q[3];
sx q[3];
rz(1.7790022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.071863197) q[0];
sx q[0];
rz(-1.3988031) q[0];
sx q[0];
rz(-0.82759696) q[0];
rz(2.1780275) q[1];
sx q[1];
rz(-1.5255442) q[1];
sx q[1];
rz(-0.73659426) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6927836) q[0];
sx q[0];
rz(-2.8570456) q[0];
sx q[0];
rz(1.8440767) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3308649) q[2];
sx q[2];
rz(-2.3849769) q[2];
sx q[2];
rz(-0.87079853) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.735966) q[1];
sx q[1];
rz(-0.36639226) q[1];
sx q[1];
rz(2.165731) q[1];
x q[2];
rz(3.0071665) q[3];
sx q[3];
rz(-1.5930015) q[3];
sx q[3];
rz(-1.7355222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5654512) q[2];
sx q[2];
rz(-0.57500035) q[2];
sx q[2];
rz(-2.3973993) q[2];
rz(-0.60892504) q[3];
sx q[3];
rz(-2.3600793) q[3];
sx q[3];
rz(1.6412546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3610483) q[0];
sx q[0];
rz(-0.86549509) q[0];
sx q[0];
rz(-1.7678827) q[0];
rz(-0.79958493) q[1];
sx q[1];
rz(-2.1345963) q[1];
sx q[1];
rz(-1.132157) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0810185) q[0];
sx q[0];
rz(-1.3558398) q[0];
sx q[0];
rz(-3.0490387) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5930873) q[2];
sx q[2];
rz(-1.6509075) q[2];
sx q[2];
rz(-1.3673283) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.76946041) q[1];
sx q[1];
rz(-0.67382183) q[1];
sx q[1];
rz(0.12427434) q[1];
x q[2];
rz(3.0789496) q[3];
sx q[3];
rz(-2.5317305) q[3];
sx q[3];
rz(0.26882879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.75480294) q[2];
sx q[2];
rz(-0.58480442) q[2];
sx q[2];
rz(1.7542138) q[2];
rz(2.7925708) q[3];
sx q[3];
rz(-1.6882221) q[3];
sx q[3];
rz(2.6121228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.741852) q[0];
sx q[0];
rz(-2.1735503) q[0];
sx q[0];
rz(-2.7834748) q[0];
rz(-1.0707567) q[1];
sx q[1];
rz(-2.4929969) q[1];
sx q[1];
rz(-1.7020285) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0107797) q[0];
sx q[0];
rz(-0.81310105) q[0];
sx q[0];
rz(-1.2050425) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1522572) q[2];
sx q[2];
rz(-0.78063595) q[2];
sx q[2];
rz(2.6065741) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4434699) q[1];
sx q[1];
rz(-1.021153) q[1];
sx q[1];
rz(0.89068954) q[1];
x q[2];
rz(0.76684769) q[3];
sx q[3];
rz(-1.7550751) q[3];
sx q[3];
rz(0.16367463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7122571) q[2];
sx q[2];
rz(-1.9108994) q[2];
sx q[2];
rz(-2.5168391) q[2];
rz(1.5077) q[3];
sx q[3];
rz(-2.3816536) q[3];
sx q[3];
rz(-2.0190575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.7655012) q[0];
sx q[0];
rz(-0.88075817) q[0];
sx q[0];
rz(1.1464024) q[0];
rz(1.7064077) q[1];
sx q[1];
rz(-2.621666) q[1];
sx q[1];
rz(-0.82040876) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7026742) q[0];
sx q[0];
rz(-1.6587388) q[0];
sx q[0];
rz(-1.2984718) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7960288) q[2];
sx q[2];
rz(-0.67906717) q[2];
sx q[2];
rz(-1.0240384) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0972406) q[1];
sx q[1];
rz(-0.28438452) q[1];
sx q[1];
rz(2.3769955) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5375529) q[3];
sx q[3];
rz(-2.1602732) q[3];
sx q[3];
rz(1.6575898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4121805) q[2];
sx q[2];
rz(-2.086144) q[2];
sx q[2];
rz(2.6030276) q[2];
rz(1.4696848) q[3];
sx q[3];
rz(-1.4476176) q[3];
sx q[3];
rz(-0.058549747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3013714) q[0];
sx q[0];
rz(-0.138962) q[0];
sx q[0];
rz(3.050991) q[0];
rz(-0.36965707) q[1];
sx q[1];
rz(-1.5934817) q[1];
sx q[1];
rz(0.37818092) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2205296) q[0];
sx q[0];
rz(-0.53595966) q[0];
sx q[0];
rz(0.45951636) q[0];
rz(-pi) q[1];
rz(0.090094968) q[2];
sx q[2];
rz(-1.5290773) q[2];
sx q[2];
rz(-0.99064529) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.82434618) q[1];
sx q[1];
rz(-2.6326944) q[1];
sx q[1];
rz(-2.8864278) q[1];
rz(-pi) q[2];
x q[2];
rz(0.079433283) q[3];
sx q[3];
rz(-1.380668) q[3];
sx q[3];
rz(-3.0832689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6211264) q[2];
sx q[2];
rz(-1.4755321) q[2];
sx q[2];
rz(-3.0401518) q[2];
rz(-1.8488047) q[3];
sx q[3];
rz(-0.83270508) q[3];
sx q[3];
rz(-1.7267936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8349649) q[0];
sx q[0];
rz(-1.8649626) q[0];
sx q[0];
rz(-2.2391338) q[0];
rz(-2.9023671) q[1];
sx q[1];
rz(-1.7996412) q[1];
sx q[1];
rz(-0.36144027) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1162565) q[0];
sx q[0];
rz(-0.5098719) q[0];
sx q[0];
rz(1.340853) q[0];
x q[1];
rz(2.7290384) q[2];
sx q[2];
rz(-0.13880402) q[2];
sx q[2];
rz(-0.96497646) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.70009184) q[1];
sx q[1];
rz(-0.48466408) q[1];
sx q[1];
rz(1.086471) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.33935762) q[3];
sx q[3];
rz(-1.4213143) q[3];
sx q[3];
rz(0.71163346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0685737) q[2];
sx q[2];
rz(-1.34015) q[2];
sx q[2];
rz(2.2646591) q[2];
rz(-0.98313156) q[3];
sx q[3];
rz(-1.5777595) q[3];
sx q[3];
rz(-0.94299281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8125732) q[0];
sx q[0];
rz(-0.1952157) q[0];
sx q[0];
rz(-0.83874291) q[0];
rz(0.70873952) q[1];
sx q[1];
rz(-2.510431) q[1];
sx q[1];
rz(0.90604025) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6853537) q[0];
sx q[0];
rz(-2.2954012) q[0];
sx q[0];
rz(1.5989717) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9211414) q[2];
sx q[2];
rz(-1.3679149) q[2];
sx q[2];
rz(1.4500993) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4455394) q[1];
sx q[1];
rz(-2.1043244) q[1];
sx q[1];
rz(-0.75977709) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2165478) q[3];
sx q[3];
rz(-2.5028526) q[3];
sx q[3];
rz(2.6415792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.16443843) q[2];
sx q[2];
rz(-2.437037) q[2];
sx q[2];
rz(-2.0257115) q[2];
rz(-0.62853938) q[3];
sx q[3];
rz(-2.0597337) q[3];
sx q[3];
rz(2.5270497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81166613) q[0];
sx q[0];
rz(-2.3217432) q[0];
sx q[0];
rz(-2.9803168) q[0];
rz(2.2992086) q[1];
sx q[1];
rz(-1.7120275) q[1];
sx q[1];
rz(-2.9715723) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1415213) q[0];
sx q[0];
rz(-1.5789487) q[0];
sx q[0];
rz(-0.018281451) q[0];
rz(-0.44012897) q[2];
sx q[2];
rz(-1.7227731) q[2];
sx q[2];
rz(-2.6226247) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.61115757) q[1];
sx q[1];
rz(-2.6621205) q[1];
sx q[1];
rz(1.1392639) q[1];
x q[2];
rz(1.6415855) q[3];
sx q[3];
rz(-2.7578691) q[3];
sx q[3];
rz(-0.078670382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9876447) q[2];
sx q[2];
rz(-1.5445856) q[2];
sx q[2];
rz(-1.6551931) q[2];
rz(-0.060128309) q[3];
sx q[3];
rz(-0.82258737) q[3];
sx q[3];
rz(1.4415584) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63115591) q[0];
sx q[0];
rz(-3.1102409) q[0];
sx q[0];
rz(1.7106868) q[0];
rz(-2.778964) q[1];
sx q[1];
rz(-1.6094094) q[1];
sx q[1];
rz(-1.8261725) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5595488) q[0];
sx q[0];
rz(-1.2799731) q[0];
sx q[0];
rz(-1.8946339) q[0];
rz(-3.0193826) q[2];
sx q[2];
rz(-1.9463681) q[2];
sx q[2];
rz(-1.5236601) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3676995) q[1];
sx q[1];
rz(-2.4871792) q[1];
sx q[1];
rz(0.44632895) q[1];
x q[2];
rz(-0.4591367) q[3];
sx q[3];
rz(-0.92184421) q[3];
sx q[3];
rz(2.5970363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0843087) q[2];
sx q[2];
rz(-1.5466651) q[2];
sx q[2];
rz(-2.9810737) q[2];
rz(2.6594243) q[3];
sx q[3];
rz(-0.67626685) q[3];
sx q[3];
rz(-2.4619861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1228444) q[0];
sx q[0];
rz(-1.3958805) q[0];
sx q[0];
rz(-0.91176283) q[0];
rz(1.979076) q[1];
sx q[1];
rz(-1.5717506) q[1];
sx q[1];
rz(2.7350978) q[1];
rz(-1.4760426) q[2];
sx q[2];
rz(-1.4838525) q[2];
sx q[2];
rz(-2.9529689) q[2];
rz(0.58070498) q[3];
sx q[3];
rz(-2.28149) q[3];
sx q[3];
rz(-2.1547439) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
