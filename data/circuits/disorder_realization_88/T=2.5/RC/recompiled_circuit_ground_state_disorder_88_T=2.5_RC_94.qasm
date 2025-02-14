OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.2849046) q[0];
sx q[0];
rz(-1.8104799) q[0];
sx q[0];
rz(-0.80944219) q[0];
rz(-2.5419905) q[1];
sx q[1];
rz(-1.124958) q[1];
sx q[1];
rz(-2.6198299) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4057465) q[0];
sx q[0];
rz(-3.0015115) q[0];
sx q[0];
rz(1.5790042) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5810165) q[2];
sx q[2];
rz(-1.3180705) q[2];
sx q[2];
rz(1.2863024) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6859036) q[1];
sx q[1];
rz(-2.0200811) q[1];
sx q[1];
rz(-2.8390769) q[1];
rz(0.14484804) q[3];
sx q[3];
rz(-2.600935) q[3];
sx q[3];
rz(0.40504211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1521505) q[2];
sx q[2];
rz(-2.4336954) q[2];
sx q[2];
rz(-1.36261) q[2];
rz(1.6705492) q[3];
sx q[3];
rz(-1.5971284) q[3];
sx q[3];
rz(-2.7267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7489557) q[0];
sx q[0];
rz(-1.3480659) q[0];
sx q[0];
rz(-1.0697399) q[0];
rz(-2.3906129) q[1];
sx q[1];
rz(-0.90194482) q[1];
sx q[1];
rz(-0.78838563) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5716924) q[0];
sx q[0];
rz(-0.65061749) q[0];
sx q[0];
rz(2.8293737) q[0];
rz(-0.090615409) q[2];
sx q[2];
rz(-0.66747016) q[2];
sx q[2];
rz(-2.8634976) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.8395811) q[1];
sx q[1];
rz(-0.55802781) q[1];
sx q[1];
rz(1.1990947) q[1];
rz(2.0072924) q[3];
sx q[3];
rz(-1.1009163) q[3];
sx q[3];
rz(1.7738455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1145733) q[2];
sx q[2];
rz(-0.35832778) q[2];
sx q[2];
rz(0.9066073) q[2];
rz(1.00057) q[3];
sx q[3];
rz(-1.118719) q[3];
sx q[3];
rz(-0.51893836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24666102) q[0];
sx q[0];
rz(-2.7592359) q[0];
sx q[0];
rz(1.3038127) q[0];
rz(1.9815824) q[1];
sx q[1];
rz(-1.8142533) q[1];
sx q[1];
rz(-1.4132285) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50158635) q[0];
sx q[0];
rz(-1.4297856) q[0];
sx q[0];
rz(1.6580576) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.63615362) q[2];
sx q[2];
rz(-1.9463123) q[2];
sx q[2];
rz(-2.2454303) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1294413) q[1];
sx q[1];
rz(-2.0291162) q[1];
sx q[1];
rz(-1.3740416) q[1];
rz(-pi) q[2];
rz(-1.71536) q[3];
sx q[3];
rz(-2.4213397) q[3];
sx q[3];
rz(-1.6584058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7580737) q[2];
sx q[2];
rz(-1.5531837) q[2];
sx q[2];
rz(-0.12913945) q[2];
rz(-2.6376851) q[3];
sx q[3];
rz(-1.4174856) q[3];
sx q[3];
rz(2.5655139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1103766) q[0];
sx q[0];
rz(-1.2126558) q[0];
sx q[0];
rz(-0.84621286) q[0];
rz(0.57202488) q[1];
sx q[1];
rz(-1.6366448) q[1];
sx q[1];
rz(-1.2342854) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18674984) q[0];
sx q[0];
rz(-0.34663793) q[0];
sx q[0];
rz(2.8524272) q[0];
rz(2.8568966) q[2];
sx q[2];
rz(-2.0139593) q[2];
sx q[2];
rz(-0.62847947) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9124536) q[1];
sx q[1];
rz(-0.22446147) q[1];
sx q[1];
rz(2.771586) q[1];
rz(0.5192786) q[3];
sx q[3];
rz(-2.2745345) q[3];
sx q[3];
rz(-0.26070538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.3884376) q[2];
sx q[2];
rz(-2.2074102) q[2];
sx q[2];
rz(-1.4446806) q[2];
rz(1.2706903) q[3];
sx q[3];
rz(-1.570574) q[3];
sx q[3];
rz(-1.1893893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63950771) q[0];
sx q[0];
rz(-1.6725699) q[0];
sx q[0];
rz(-2.0462659) q[0];
rz(-3.0785839) q[1];
sx q[1];
rz(-1.6687702) q[1];
sx q[1];
rz(0.94620401) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99710714) q[0];
sx q[0];
rz(-1.5988264) q[0];
sx q[0];
rz(-1.0260242) q[0];
x q[1];
rz(1.191372) q[2];
sx q[2];
rz(-0.96948955) q[2];
sx q[2];
rz(1.7446818) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7530147) q[1];
sx q[1];
rz(-1.565897) q[1];
sx q[1];
rz(0.51112643) q[1];
x q[2];
rz(-1.9619476) q[3];
sx q[3];
rz(-1.3538651) q[3];
sx q[3];
rz(-2.2892078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5064064) q[2];
sx q[2];
rz(-2.0621641) q[2];
sx q[2];
rz(1.6335454) q[2];
rz(0.15277282) q[3];
sx q[3];
rz(-1.2341876) q[3];
sx q[3];
rz(2.2085786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9043115) q[0];
sx q[0];
rz(-0.58552423) q[0];
sx q[0];
rz(3.0329419) q[0];
rz(-0.12734224) q[1];
sx q[1];
rz(-1.8904949) q[1];
sx q[1];
rz(0.88643518) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21289794) q[0];
sx q[0];
rz(-1.587364) q[0];
sx q[0];
rz(0.51942708) q[0];
x q[1];
rz(1.7753052) q[2];
sx q[2];
rz(-1.5123774) q[2];
sx q[2];
rz(-1.2893334) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8005714) q[1];
sx q[1];
rz(-1.6280975) q[1];
sx q[1];
rz(-1.351746) q[1];
x q[2];
rz(0.16951321) q[3];
sx q[3];
rz(-1.680114) q[3];
sx q[3];
rz(-1.5284571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2966557) q[2];
sx q[2];
rz(-2.3462494) q[2];
sx q[2];
rz(-2.8080158) q[2];
rz(-0.9345471) q[3];
sx q[3];
rz(-0.75342527) q[3];
sx q[3];
rz(0.88753382) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6880671) q[0];
sx q[0];
rz(-1.4997361) q[0];
sx q[0];
rz(-2.8330579) q[0];
rz(2.535179) q[1];
sx q[1];
rz(-2.2589222) q[1];
sx q[1];
rz(0.44953129) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36469421) q[0];
sx q[0];
rz(-0.5185655) q[0];
sx q[0];
rz(3.0125951) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3541127) q[2];
sx q[2];
rz(-1.228294) q[2];
sx q[2];
rz(2.1975348) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7779904) q[1];
sx q[1];
rz(-2.5153179) q[1];
sx q[1];
rz(1.3961149) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1949439) q[3];
sx q[3];
rz(-0.80393857) q[3];
sx q[3];
rz(0.80292668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0799847) q[2];
sx q[2];
rz(-1.6500762) q[2];
sx q[2];
rz(-1.9645346) q[2];
rz(2.4401149) q[3];
sx q[3];
rz(-0.25291118) q[3];
sx q[3];
rz(-2.285932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6845282) q[0];
sx q[0];
rz(-0.51790154) q[0];
sx q[0];
rz(-1.6494226) q[0];
rz(-2.0860784) q[1];
sx q[1];
rz(-1.5558473) q[1];
sx q[1];
rz(0.42748705) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1221177) q[0];
sx q[0];
rz(-0.97163288) q[0];
sx q[0];
rz(2.504306) q[0];
rz(-pi) q[1];
rz(1.4094818) q[2];
sx q[2];
rz(-1.1448793) q[2];
sx q[2];
rz(0.49954712) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6969949) q[1];
sx q[1];
rz(-0.31365221) q[1];
sx q[1];
rz(3.0693377) q[1];
rz(-pi) q[2];
rz(0.31510809) q[3];
sx q[3];
rz(-1.3118532) q[3];
sx q[3];
rz(-1.068598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.18275729) q[2];
sx q[2];
rz(-1.3311102) q[2];
sx q[2];
rz(1.9367564) q[2];
rz(3.0086573) q[3];
sx q[3];
rz(-3.0429621) q[3];
sx q[3];
rz(1.0814063) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8590915) q[0];
sx q[0];
rz(-1.4097255) q[0];
sx q[0];
rz(2.6891563) q[0];
rz(0.85894194) q[1];
sx q[1];
rz(-0.57933885) q[1];
sx q[1];
rz(-1.4525684) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70266028) q[0];
sx q[0];
rz(-2.0468674) q[0];
sx q[0];
rz(0.66113801) q[0];
rz(0.33272393) q[2];
sx q[2];
rz(-2.1004268) q[2];
sx q[2];
rz(-2.4650448) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7048129) q[1];
sx q[1];
rz(-1.5680934) q[1];
sx q[1];
rz(2.0767434) q[1];
rz(-2.7730016) q[3];
sx q[3];
rz(-1.907371) q[3];
sx q[3];
rz(-2.9934538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.23898827) q[2];
sx q[2];
rz(-1.3024104) q[2];
sx q[2];
rz(-2.7461309) q[2];
rz(-2.0579193) q[3];
sx q[3];
rz(-0.62267059) q[3];
sx q[3];
rz(-1.7284988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.285242) q[0];
sx q[0];
rz(-0.1150035) q[0];
sx q[0];
rz(1.2913936) q[0];
rz(2.3639288) q[1];
sx q[1];
rz(-1.5592557) q[1];
sx q[1];
rz(-1.2819598) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7972231) q[0];
sx q[0];
rz(-0.31569052) q[0];
sx q[0];
rz(0.80202054) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7927112) q[2];
sx q[2];
rz(-2.2063428) q[2];
sx q[2];
rz(-1.9161621) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.79445542) q[1];
sx q[1];
rz(-2.0973118) q[1];
sx q[1];
rz(1.305425) q[1];
rz(1.3814209) q[3];
sx q[3];
rz(-0.52703372) q[3];
sx q[3];
rz(0.98637146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.164244) q[2];
sx q[2];
rz(-0.35173309) q[2];
sx q[2];
rz(2.742761) q[2];
rz(-2.2073958) q[3];
sx q[3];
rz(-2.4103006) q[3];
sx q[3];
rz(-3.0493693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4451404) q[0];
sx q[0];
rz(-2.2284989) q[0];
sx q[0];
rz(-3.1365119) q[0];
rz(2.483881) q[1];
sx q[1];
rz(-1.7291768) q[1];
sx q[1];
rz(-1.9531858) q[1];
rz(-1.2451386) q[2];
sx q[2];
rz(-0.79380582) q[2];
sx q[2];
rz(0.028887916) q[2];
rz(-1.4767968) q[3];
sx q[3];
rz(-0.58782676) q[3];
sx q[3];
rz(1.6997433) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
